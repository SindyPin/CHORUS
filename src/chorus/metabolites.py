"""Metabolite identifiers shared across cohorts.

Cohorts deposit metabolites under their own names ("L-Glutamine", "Glutamine", "Gln", "PC aa
C36:6"). ``metabolite_ids`` maps each deposited name to one shared identifier:

1. the name (with a platform prefix ``<platform>::`` removed) is standardised to a RefMet name by
   the Metabolomics Workbench RefMet name conversion service (batch POST; answers are cached in a
   CSV so that a rebuild does not query again);
2. the RefMet name is looked up in the RefMet table (``refmet.csv``, the full download), which
   gives the RefMet id, HMDB id, InChIKey, PubChem and KEGG ids;
3. the shared identifier is ``HMDB:<id>`` when RefMet gives an HMDB id; otherwise an HMDB id
   deposited by the authors for that name (``deposited``); otherwise ``RefMet:<refmet name>``;
   otherwise the deposited InChIKey (``InChIKey:<key>``); otherwise the name is unmapped.

The identifier is a function of the RefMet name when RefMet matches, so two cohorts whose names
standardise to the same RefMet name get the same identifier. Names that look like unannotated
features (``<m/z>_<RT>``, ``Unknown ...``) are not sent to the match service.
"""

from __future__ import annotations

import re
import time
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Callable, Dict, Iterable, List, Optional

import pandas as pd

BATCH_URL = "https://www.metabolomicsworkbench.org/databases/refmet/name_to_refmet_new_minID.php"
UNANNOTATED = re.compile(r"^(\d+(\.\d+)?_\d+(\.\d+)?|unknown\b.*|unk\b.*|x\s*-\s*\d+.*)$", re.IGNORECASE)


def strip_prefix(name: str) -> str:
    """Name without a ``<platform>::`` prefix (``GCMS::Glutamine`` -> ``Glutamine``). A double colon
    marks the prefix because single colons occur in lipid names (``PC 36:4``)."""
    s = str(name).strip()
    return s.split("::", 1)[1].strip() if "::" in s else s


def load_refmet(path: str | Path) -> pd.DataFrame:
    """RefMet table indexed by refmet_name (first row kept for a repeated name)."""
    t = pd.read_csv(path, dtype=str)
    return t.drop_duplicates("refmet_name").set_index("refmet_name")


def fetch_match(names: List[str], timeout: float = 120.0) -> Dict[str, str]:
    """Query name -> RefMet name ('' when there is no match) for a batch of names, from the RefMet
    name conversion service (``BATCH_URL``, one POST with the names separated by new lines; it
    accepts names with "/", which the REST match service answers with 404)."""
    import io

    data = urllib.parse.urlencode({"metabolite_name": "\n".join(names)}).encode()
    with urllib.request.urlopen(urllib.request.Request(BATCH_URL, data=data), timeout=timeout) as r:
        body = r.read().decode("utf-8", "replace")
    t = pd.read_csv(io.StringIO(body), sep="\t", dtype=str, keep_default_na=False) if body.strip() else pd.DataFrame(columns=["Input name", "Standardized name"])
    out = {}
    for q, rn in zip(t["Input name"], t["Standardized name"]):
        out[q] = "" if rn in ("-", "") else rn
    return {q: out.get(q, "") for q in names}


def match_names(names: Iterable[str], cache: str | Path, fetch: Optional[Callable[[List[str]], Dict[str, str]]] = fetch_match, chunk: int = 300, pause: float = 1.0) -> pd.Series:
    """Query name -> RefMet name ('' when unmatched), from ``cache`` (CSV: query, refmet_name,
    queried) and, for names not in it, from ``fetch`` in batches of ``chunk`` (None: no network,
    missing stays NaN). New answers are appended to the cache after every batch."""
    cache = Path(cache)
    known = pd.read_csv(cache, dtype=str, keep_default_na=False) if cache.exists() else pd.DataFrame(columns=["query", "refmet_name", "queried"])
    lut = dict(zip(known["query"], known["refmet_name"]))
    todo = [q for q in dict.fromkeys(names) if q not in lut]
    if fetch is not None:
        for i in range(0, len(todo), chunk):
            part = todo[i : i + chunk]
            try:
                ans = fetch(part)
            except Exception as e:  # network failure: leave the batch unresolved, record nothing
                print(f"[refmet] batch {i // chunk}: {e}")
                continue
            new = pd.DataFrame({"query": part, "refmet_name": [ans.get(q, "") for q in part], "queried": time.strftime("%Y-%m-%d")})
            lut.update(zip(new["query"], new["refmet_name"]))
            known = pd.concat([known, new], ignore_index=True)
            cache.parent.mkdir(parents=True, exist_ok=True)
            known.to_csv(cache, index=False)
            time.sleep(pause)
    return pd.Series({q: lut.get(q) for q in dict.fromkeys(names)}, dtype=object)


def metabolite_ids(names: Iterable[str], refmet: pd.DataFrame, matches: pd.Series, deposited: Optional[pd.DataFrame] = None) -> pd.DataFrame:
    """One row per deposited name: query (name sent to RefMet), refmet_name, refmet_id, hmdb_id,
    inchi_key, kegg_id, shared_id and id_source (hmdb_refmet, hmdb_deposited, refmet,
    inchikey_deposited, unannotated, unmatched). ``deposited`` (indexed by name) may hold columns
    hmdb_id and inchi_key given by the authors."""
    rows = []
    for n in dict.fromkeys(names):
        q = strip_prefix(n)
        r = {"name": n, "query": q, "refmet_name": None, "refmet_id": None, "hmdb_id": None, "inchi_key": None, "kegg_id": None}
        dep = deposited.loc[n] if deposited is not None and n in deposited.index else None
        dep_hmdb = str(dep.get("hmdb_id")) if dep is not None and pd.notna(dep.get("hmdb_id")) and str(dep.get("hmdb_id")).startswith("HMDB") else None
        dep_ik = str(dep.get("inchi_key")) if dep is not None and pd.notna(dep.get("inchi_key")) else None
        rn = matches.get(q)
        if UNANNOTATED.match(q):
            src, sid = "unannotated", None
        elif isinstance(rn, str) and rn and rn in refmet.index:
            m = refmet.loc[rn]
            r.update(refmet_name=rn, refmet_id=m.get("refmet_id"), hmdb_id=m.get("hmdb_id") if pd.notna(m.get("hmdb_id")) else None, inchi_key=m.get("inchi_key") if pd.notna(m.get("inchi_key")) else None, kegg_id=m.get("kegg_id") if pd.notna(m.get("kegg_id")) else None)
            if r["hmdb_id"]:
                src, sid = "hmdb_refmet", f"HMDB:{r['hmdb_id']}"
            elif dep_hmdb:
                src, sid = "hmdb_deposited", f"HMDB:{dep_hmdb}"
            else:
                src, sid = "refmet", f"RefMet:{rn}"
        elif isinstance(rn, str) and rn:  # standardised name absent from the table (for example "A or B")
            r.update(refmet_name=rn)
            src, sid = ("hmdb_deposited", f"HMDB:{dep_hmdb}") if dep_hmdb else ("refmet", f"RefMet:{rn}")
        elif dep_hmdb:
            src, sid = "hmdb_deposited", f"HMDB:{dep_hmdb}"
        elif dep_ik:
            src, sid = "inchikey_deposited", f"InChIKey:{dep_ik}"
        else:
            src, sid = "unmatched", None
        r.update(shared_id=sid, id_source=src)
        rows.append(r)
    return pd.DataFrame(rows)


def pairwise_shared(ids: Dict[str, Iterable[str]]) -> pd.DataFrame:
    """Cohort x cohort number of shared identifiers (diagonal: identifiers in the cohort)."""
    sets = {k: set(v) for k, v in ids.items()}
    keys = list(sets)
    return pd.DataFrame([[len(sets[a] & sets[b]) for b in keys] for a in keys], index=keys, columns=keys)
