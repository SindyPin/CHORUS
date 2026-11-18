"""Acute COVID cohorts: find, download, build and verify (Part B of the acute COVID validation).

The rule and bars of the validation are fixed in results/README.md ("Acute COVID validation: rule and bars"); this
script only finds and prepares the cohorts and runs no benchmark or noise ceiling analysis on them.

Subcommands, each resumable (finished work is cached or skipped) and logged to logs/acute_covid_sweep.log; the script
writes only under data/, docs/ and logs/, except that report writes its own new configs configs/acute_covid_*.yaml
(enabled: false; an existing config it did not write is never overwritten):

  search   query GEO (E-utilities), ArrayExpress and BioStudies, SRA BioProjects, PRIDE and ProteomeXchange,
           Metabolomics Workbench, Zenodo, Mendeley Data (DataCite) and Europe PMC (data availability statements of
           the papers found); every response is cached under data/acute_covid/cache/. Output: one row per candidate in
           docs/acute_covid_candidates.csv.
  screen   metadata only. Automatic exclusions (organism, single cell, tissue, sample count, assay) and the curated
           decisions (CURATED) for the candidates read by hand; decision columns added to the candidates file.
  fetch    downloads of the included cohorts (deposited matrices preferred), with resume, three retries, repository
           checksums where given, and size and SHA-256 of every file (data/acute_covid/files.csv).
  build    per cohort metadata in the CHORUS schema (data/meta/<cohort>_meta.csv) with the contrast classes and the
           label mapping table (docs/acute_covid_labels/<cohort>.csv), and the feature matrix with mapped ids.
  verify   checks a to k per cohort (data/acute_covid/verify/<cohort>.csv).
  report   docs/acute_covid_cohorts.md, data/meta/acute_covid_manifest.csv and configs (enabled: false).

  python scripts/acute_covid_sweep.py search|screen|fetch|build|verify|report [--cohort ID]
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import io
import json
import logging
import re
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import requests

ROOT = Path(__file__).resolve().parents[1]
WORK = ROOT / "data" / "acute_covid"
CACHE = WORK / "cache"
DOCS = ROOT / "docs"
CANDIDATES = DOCS / "acute_covid_candidates.csv"
LOG = ROOT / "logs" / "acute_covid_sweep.log"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
YEARS = "2020:2026"
COVID = "(COVID-19 OR SARS-CoV-2)"

log = logging.getLogger("acute_covid_sweep")


def setup_logging() -> None:
    LOG.parent.mkdir(parents=True, exist_ok=True)
    if not log.handlers:
        log.setLevel(logging.INFO)
        for h in (logging.FileHandler(LOG), logging.StreamHandler(sys.stdout)):
            h.setFormatter(logging.Formatter("%(asctime)s %(levelname)s %(message)s"))
            log.addHandler(h)


# ------------------------------------------------------------------------------------------------ http with cache
def _key(url: str, params: dict | None) -> str:
    return hashlib.sha256((url + json.dumps(params or {}, sort_keys=True)).encode()).hexdigest()[:24]


def get(url: str, params: dict | None = None, kind: str = "json", tries: int = 3, sleep: float = 0.4):
    """GET with a disk cache (data/acute_covid/cache) and three retries; returns parsed JSON or text."""
    CACHE.mkdir(parents=True, exist_ok=True)
    f = CACHE / f"{_key(url, params)}.{'json' if kind == 'json' else 'txt'}"
    if f.exists():
        return json.loads(f.read_text()) if kind == "json" else f.read_text()
    err = None
    for i in range(tries):
        try:
            r = requests.get(url, params=params, timeout=120, headers={"Accept": "application/json" if kind == "json" else "*/*"})
            if r.status_code == 404:
                f.write_text("null" if kind == "json" else "")
                return None if kind == "json" else ""
            r.raise_for_status()
            if kind == "json":
                data = r.json()
                f.write_text(json.dumps(data))
                time.sleep(sleep)
                return data
            f.write_text(r.text)
            time.sleep(sleep)
            return r.text
        except Exception as e:  # noqa: BLE001
            err = e
            time.sleep(2 * (i + 1))
    log.warning(f"GET failed after {tries} tries: {url} {params} ({err})")
    return None if kind == "json" else ""


# ------------------------------------------------------------------------------------------------ search
SEARCH_QUERIES = {
    "geo_rnaseq": ("gds", f'{COVID} AND (blood OR "whole blood" OR PBMC OR leukocytes OR plasma OR serum) AND "Homo sapiens"[Organism] '
                          f'AND "expression profiling by high throughput sequencing"[DataSet Type] AND gse[ETYP] AND {YEARS}[PDAT]'),
    "geo_array": ("gds", f'{COVID} AND (blood OR "whole blood" OR PBMC OR leukocytes) AND "Homo sapiens"[Organism] '
                         f'AND "expression profiling by array"[DataSet Type] AND gse[ETYP] AND {YEARS}[PDAT]'),
    "geo_protein": ("gds", f'{COVID} AND (plasma OR serum OR blood) AND "Homo sapiens"[Organism] AND (proteomics OR Olink OR SomaScan OR "protein profiling by protein array"[DataSet Type]) '
                           f'AND gse[ETYP] AND {YEARS}[PDAT]'),
    "bioproject": ("bioproject", f'{COVID} AND (blood OR PBMC OR plasma OR serum) AND human[orgn]'),
}
EPMC_BASE = ('ABSTRACT:("COVID-19" OR "SARS-CoV-2") AND OPEN_ACCESS:y AND PUB_YEAR:[2020 TO 2026] AND (ACCESSION_TYPE:geo OR ACCESSION_TYPE:pxd OR '
             'ACCESSION_TYPE:arrayexpress OR ACCESSION_TYPE:bioproject OR ACCESSION_TYPE:metabolights OR "zenodo" OR "Mendeley Data" OR "Metabolomics Workbench" OR "MassIVE")')
EPMC_QUERIES = [  # papers that deposit data (accession in the text, or a data repository named), with the layer and tissue in the abstract
    EPMC_BASE + ' AND ABSTRACT:(plasma OR serum) AND ABSTRACT:(Olink OR SomaScan OR proteomics OR proteome OR proteomic)',
    EPMC_BASE + ' AND ABSTRACT:("whole blood" OR PBMC OR blood OR leukocytes) AND ABSTRACT:(RNA-seq OR "RNA sequencing" OR transcriptome OR transcriptomic OR transcriptional)',
    EPMC_BASE + ' AND ABSTRACT:(plasma OR serum) AND ABSTRACT:(metabolomics OR metabolome OR metabolomic OR lipidomics)',
]
FROM_MEMORY = ["GSE157103", "GSE152418", "GSE161731", "GSE152641", "GSE171110", "GSE172114", "GSE198449", "GSE155454", "GSE163151", "MGH_Olink_Filbin2021"]
ON_DISK_ACUTE = {"GSE215865": "bulk_rna", "Su2022": "proteomics", "Cervia2024_SomaScan": "proteomics", "Zoodsma2022_Hannover": "proteomics", "Woodruff2023": "proteomics"}
ACC_RE = re.compile(r"\b(GSE\d{4,7}|PXD\d{6}|MSV\d{9}|ST\d{6}|E-MTAB-\d+|PRJNA\d+|EGAS\d+|EGAD\d+|syn\d{6,9})\b")
DOI_RE = re.compile(r"\b(10\.5281/zenodo\.\d+|10\.17632/[a-z0-9]+(?:\.\d+)?)", re.I)


def _esearch(db: str, term: str) -> list:
    ids, start = [], 0
    while True:
        r = get(EUTILS + "esearch.fcgi", {"db": db, "term": term, "retmax": 500, "retstart": start, "retmode": "json"})
        if not r:
            break
        got = r["esearchresult"].get("idlist", [])
        ids += got
        start += len(got)
        if not got or start >= int(r["esearchresult"]["count"]):
            break
    return ids


def _esummary(db: str, ids: list) -> list:
    out = []
    for i in range(0, len(ids), 200):
        r = get(EUTILS + "esummary.fcgi", {"db": db, "id": ",".join(ids[i:i + 200]), "retmode": "json"})
        if r and "result" in r:
            out += [r["result"][u] for u in r["result"].get("uids", [])]
    return out


def search_geo() -> list:
    rows = []
    for name in ("geo_rnaseq", "geo_array", "geo_protein"):
        db, term = SEARCH_QUERIES[name]
        ids = _esearch(db, term)
        log.info(f"search {name}: {len(ids)} series")
        for s in _esummary(db, ids):
            rows.append({"source": "GEO", "query": name, "accession": s.get("accession"), "title": s.get("title"), "description": s.get("summary", "")[:2000],
                         "organism": s.get("taxon"), "n_samples": s.get("n_samples"), "assay": s.get("gdstype"), "platform": "; ".join("GPL" + p for p in str(s.get("gpl", "")).split(";") if p),
                         "pmids": ";".join(str(p) for p in s.get("pubmedids", [])), "date": s.get("pdat"), "supplementary": s.get("suppfile", ""), "bioproject": s.get("bioproject", ""),
                         "url": f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={s.get('accession')}"})
    return rows


def search_arrayexpress() -> list:
    rows = []
    for q in ['"COVID-19" AND (blood OR PBMC OR plasma OR serum)', '"SARS-CoV-2" AND (blood OR PBMC OR plasma OR serum)']:
        page = 1
        while True:
            r = get("https://www.ebi.ac.uk/biostudies/api/v1/arrayexpress/search", {"query": q, "pageSize": 100, "page": page})
            hits = (r or {}).get("hits", [])
            rows += [{"source": "ArrayExpress", "query": q, "accession": h.get("accession"), "title": h.get("title"), "description": str(h.get("content", ""))[:2000],
                      "n_samples": None, "date": h.get("release_date"), "url": f"https://www.ebi.ac.uk/biostudies/arrayexpress/studies/{h.get('accession')}"} for h in hits]
            if len(hits) < 100:
                break
            page += 1
    log.info(f"search arrayexpress: {len(rows)} hits")
    return rows


def search_pride() -> list:
    rows = []
    for kw in ["COVID-19 plasma", "COVID-19 serum", "SARS-CoV-2 plasma", "SARS-CoV-2 serum"]:
        page = 0
        while True:
            r = get("https://www.ebi.ac.uk/pride/ws/archive/v2/search/projects", {"keyword": kw, "pageSize": 100, "page": page}) or []
            rows += [{"source": "PRIDE", "query": kw, "accession": p.get("accession"), "title": p.get("title"), "description": str(p.get("projectDescription", ""))[:2000],
                      "organism": "; ".join(p.get("organisms", []) or []), "date": p.get("publicationDate"), "url": f"https://www.ebi.ac.uk/pride/archive/projects/{p.get('accession')}"} for p in r]
            if len(r) < 100:
                break
            page += 1
    r = get("https://proteomecentral.proteomexchange.org/api/proxi/v0.1/datasets", {"resultType": "compact", "pageSize": 1000, "keywords": "COVID-19"})
    for d in (r or {}).get("datasets", []):
        rows.append({"source": "ProteomeXchange", "query": "COVID-19", "accession": d[0], "title": d[1], "organism": d[3], "description": "", "url": f"https://proteomecentral.proteomexchange.org/cgi/GetDataset?ID={d[0]}"})
    log.info(f"search pride and proteomexchange: {len(rows)} hits")
    return rows


def search_metabolomics_workbench() -> list:
    rows = []
    for w in ["COVID", "SARS-CoV-2"]:
        r = get(f"https://www.metabolomicsworkbench.org/rest/study/study_title/{w}/summary") or {}
        r = [r] if "study_id" in r else list(r.values())
        rows += [{"source": "Metabolomics Workbench", "query": w, "accession": s.get("study_id"), "title": s.get("study_title"), "description": str(s.get("study_summary", ""))[:2000],
                  "organism": s.get("species"), "n_samples": s.get("number_of_samples"), "assay": s.get("analysis_type"), "date": s.get("submit_date"),
                  "url": f"https://www.metabolomicsworkbench.org/data/DRCCMetadata.php?Mode=Study&StudyID={s.get('study_id')}"} for s in r]
    log.info(f"search metabolomics workbench: {len(rows)} hits")
    return rows


def search_zenodo() -> list:
    rows = []
    for q in [f'{COVID} AND (plasma OR serum OR blood) AND (Olink OR SomaScan OR proteomics OR metabolomics OR RNA-seq OR transcriptome)']:
        page = 1
        while True:
            r = get("https://zenodo.org/api/records", {"q": q, "type": "dataset", "size": 25, "page": page}) or {}
            hits = r.get("hits", {}).get("hits", [])
            rows += [{"source": "Zenodo", "query": q, "accession": f"zenodo.{h.get('id')}", "title": h.get("metadata", {}).get("title"),
                      "description": re.sub("<[^>]+>", " ", str(h.get("metadata", {}).get("description", "")))[:2000], "date": h.get("metadata", {}).get("publication_date"),
                      "url": h.get("links", {}).get("self_html")} for h in hits]
            if len(hits) < 25 or page >= 20:
                break
            page += 1
    log.info(f"search zenodo: {len(rows)} hits")
    return rows


def search_mendeley() -> list:
    rows = []
    q = f'{COVID} AND (plasma OR serum OR blood) AND (Olink OR SomaScan OR proteomics OR proteome OR metabolomics OR RNA-seq OR transcriptome)'
    page = 1
    while True:
        r = get("https://api.datacite.org/dois", {"query": q, "prefix": "10.17632", "page[size]": 100, "page[number]": page}) or {}
        data = r.get("data", [])
        for d in data:
            a = d.get("attributes", {})
            rows.append({"source": "Mendeley Data", "query": q, "accession": a.get("doi"), "title": "; ".join(t.get("title", "") for t in a.get("titles", [])),
                         "description": " ".join(x.get("description", "") or "" for x in a.get("descriptions", []))[:2000], "date": a.get("publicationYear"), "url": f"https://doi.org/{a.get('doi')}"})
        if len(data) < 100 or page >= 10:
            break
        page += 1
    log.info(f"search mendeley (datacite): {len(rows)} hits")
    return rows


def search_europepmc() -> list:
    """Papers on acute COVID blood omics and the accessions in their full text (data availability statements)."""
    rows = []
    for q in EPMC_QUERIES:
        cursor, n = "*", 0
        while True:
            r = get("https://www.ebi.ac.uk/europepmc/webservices/rest/search", {"query": q, "format": "json", "pageSize": 1000, "cursorMark": cursor, "resultType": "lite"}) or {}
            res = r.get("resultList", {}).get("result", [])
            for p in res:
                rows.append({"source": "Europe PMC", "query": q, "pmid": p.get("pmid"), "pmcid": p.get("pmcid"), "doi": p.get("doi"), "title": p.get("title"), "date": p.get("pubYear")})
            n += len(res)
            if not res or r.get("nextCursorMark") in (None, cursor) or n >= 3000:
                break
            cursor = r["nextCursorMark"]
    papers = pd.DataFrame(rows).drop_duplicates("pmcid")
    log.info(f"search europe pmc: {len(papers)} open access papers")
    from concurrent.futures import ThreadPoolExecutor

    todo = list(papers.dropna(subset=["pmcid"]).itertuples())
    with ThreadPoolExecutor(8) as ex:  # full texts in parallel (cached)
        texts = list(ex.map(lambda p: get(f"https://www.ebi.ac.uk/europepmc/webservices/rest/{p.pmcid}/fullTextXML", kind="text", sleep=0.1) or "", todo))
    out = []
    for p, txt in zip(todo, texts):
        accs = sorted(set(ACC_RE.findall(txt)) | {d.lower() for d in DOI_RE.findall(txt)})
        for a in accs:
            out.append({"source": "paper", "query": "Europe PMC full text", "accession": a, "title": p.title, "pmids": p.pmid, "doi_paper": p.doi, "pmcid": p.pmcid, "date": p.date})
    log.info(f"accessions in the papers' full text: {len(out)} (paper, accession) pairs")
    return out


def _bioproject_geo_links(uids: list) -> dict:
    """BioProject uid -> GEO series accessions linked to it (elink), so that a BioProject already in GEO is not a second
    candidate."""
    out = {}
    for i in range(0, len(uids), 20):
        r = get(EUTILS + "elink.fcgi", {"dbfrom": "bioproject", "db": "gds", "id": uids[i:i + 20], "retmode": "json"})
        for ls in (r or {}).get("linksets", []):
            gds = [x for db in ls.get("linksetdbs", []) for x in db.get("links", [])]
            for u in ls.get("ids", []):
                out[str(u)] = gds
    gds_all = sorted({g for v in out.values() for g in v})
    acc = {str(s["uid"]): s.get("accession") for s in _esummary("gds", gds_all)} if gds_all else {}
    return {u: sorted({acc.get(str(g)) for g in v if str(acc.get(str(g), "")).startswith("GSE")}) for u, v in out.items()}


def search_bioproject() -> list:
    db, term = SEARCH_QUERIES["bioproject"]
    ids = _esearch(db, term)
    log.info(f"search bioproject: {len(ids)} projects")
    links: dict = {}  # GEO links come from the GEO summaries (build_candidates); elink only for the ones left (screen)
    rows = []
    for s in _esummary(db, ids):
        rows.append({"source": "BioProject", "query": "bioproject", "accession": s.get("project_acc"), "title": s.get("project_title"), "description": str(s.get("project_description", ""))[:2000],
                     "organism": s.get("organism_name"), "n_samples": None, "assay": s.get("project_data_type"), "date": s.get("registration_date"),
                     "linked_geo": ";".join(links.get(str(s.get("uid")), [])), "url": f"https://www.ncbi.nlm.nih.gov/bioproject/{s.get('project_acc')}"})
    return rows


LAYER_WORDS = {"proteins": r"olink|somascan|somalogic|proteom|protein|aptamer|proximity extension|mass spectrometry|lc-ms/ms|dia-ms",
               "metabolites": r"metabolom|metabolit|lipidom",
               "RNA": r"rna-seq|rna seq|transcriptom|transcriptional|gene expression|mrna|expression profiling"}
TISSUE_WORDS = [("whole blood", r"whole blood|paxgene|tempus"), ("PBMC", r"\bpbmc|peripheral blood mononuclear"), ("leukocytes", r"leukocyte|leucocyte|buffy coat"),
                ("plasma", r"\bplasma\b"), ("serum", r"\bserum\b|\bsera\b"), ("blood", r"\bblood\b")]


def _layer(r) -> str:
    if r["source"] in ("PRIDE", "ProteomeXchange") or str(r["accession"]).startswith(("PXD", "MSV")):
        return "proteins"
    if r["source"] == "Metabolomics Workbench" or str(r["accession"]).startswith("ST"):
        return "metabolites"
    if r["source"] == "GEO" and r.get("query") in ("geo_rnaseq", "geo_array"):
        return "RNA"
    txt = f"{r.get('title', '')} {r.get('description', '')}".lower()
    for layer, pat in LAYER_WORDS.items():
        if re.search(pat, txt):
            return layer
    return "unknown"


def _tissue(txt: str) -> str:
    txt = txt.lower()
    return "; ".join(t for t, pat in TISSUE_WORDS if re.search(pat, txt)) or "not stated"


def _access(acc: str, source: str) -> str:
    if acc.startswith(("EGAS", "EGAD", "phs", "syn")):
        return "controlled (application)"
    return "open"


def _deposited(r) -> str:
    if r["source"] == "GEO":
        supp = str(r.get("supplementary") or "")
        return f"processed in GEO supplement ({supp})" if supp and supp != "nan" else "raw reads only (SRA)"
    if r["source"] in ("PRIDE", "ProteomeXchange"):
        return "MS raw files and search results (matrix to check)"
    if r["source"] == "Metabolomics Workbench":
        return "processed (Metabolomics Workbench data tables)"
    if r["source"] in ("Zenodo", "Mendeley Data"):
        return "files of the record (to check)"
    return "to check"


def _labels(r) -> str:
    if r["source"] == "GEO":
        return "GEO sample characteristics"
    if r["source"] == "Metabolomics Workbench":
        return "Metabolomics Workbench sample factors"
    if r["source"] in ("PRIDE", "ProteomeXchange"):
        return "SDRF or paper supplement (to check)"
    return "record files or paper supplement (to check)"


def enrich_paper_accessions(H: pd.DataFrame) -> pd.DataFrame:
    """Repository records for accessions seen only in papers (GEO, PRIDE, BioProject, Zenodo, Mendeley Data,
    Metabolomics Workbench), so that they can be screened like the search hits."""
    have = set(H.loc[H["source"] != "paper", "accession"].astype(str))
    accs = sorted(set(H.loc[H["source"] == "paper", "accession"].astype(str)) - have)
    rows = []
    gse = [a for a in accs if a.startswith("GSE")]
    for i in range(0, len(gse), 100):
        ids = _esearch("gds", " OR ".join(f"{a}[ACCN]" for a in gse[i:i + 100]) + " AND gse[ETYP]")
        for s in _esummary("gds", ids):
            rows.append({"source": "GEO", "query": "paper accession", "accession": s.get("accession"), "title": s.get("title"), "description": s.get("summary", "")[:2000],
                         "organism": s.get("taxon"), "n_samples": s.get("n_samples"), "assay": s.get("gdstype"), "platform": "; ".join("GPL" + p for p in str(s.get("gpl", "")).split(";") if p),
                         "pmids": ";".join(str(p) for p in s.get("pubmedids", [])), "date": s.get("pdat"), "supplementary": s.get("suppfile", ""), "bioproject": s.get("bioproject", ""),
                         "url": f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={s.get('accession')}"})
    for a in [a for a in accs if a.startswith("PXD")]:
        p = get(f"https://www.ebi.ac.uk/pride/ws/archive/v2/projects/{a}")
        if p:
            rows.append({"source": "PRIDE", "query": "paper accession", "accession": a, "title": p.get("title"), "description": str(p.get("projectDescription", ""))[:2000],
                         "organism": "; ".join(o.get("name", "") for o in p.get("organisms", []) or []), "date": p.get("publicationDate"), "url": f"https://www.ebi.ac.uk/pride/archive/projects/{a}"})
    bp = [a for a in accs if a.startswith("PRJ")]
    if bp:
        ids = []
        for i in range(0, len(bp), 100):
            ids += _esearch("bioproject", " OR ".join(f"{a}[PRJA]" for a in bp[i:i + 100]))
        for s in _esummary("bioproject", ids):
            rows.append({"source": "BioProject", "query": "paper accession", "accession": s.get("project_acc"), "title": s.get("project_title"),
                         "description": str(s.get("project_description", ""))[:2000], "organism": s.get("organism_name"), "assay": s.get("project_data_type"),
                         "url": f"https://www.ncbi.nlm.nih.gov/bioproject/{s.get('project_acc')}"})
    for a in [a for a in accs if a.startswith(("10.5281", "10.17632"))]:
        r = get(f"https://api.datacite.org/dois/{a}")
        at = ((r or {}).get("data") or {}).get("attributes", {})
        if at:
            rows.append({"source": "Zenodo" if a.startswith("10.5281") else "Mendeley Data", "query": "paper accession", "accession": a,
                         "title": "; ".join(t.get("title", "") for t in at.get("titles", [])), "description": " ".join(x.get("description", "") or "" for x in at.get("descriptions", []))[:2000],
                         "date": at.get("publicationYear"), "url": f"https://doi.org/{a}"})
    for a in [a for a in accs if re.fullmatch(r"ST\d{6}", a)]:
        s = get(f"https://www.metabolomicsworkbench.org/rest/study/study_id/{a}/summary") or {}
        if s:
            rows.append({"source": "Metabolomics Workbench", "query": "paper accession", "accession": a, "title": s.get("study_title"), "description": str(s.get("study_summary", ""))[:2000],
                         "organism": s.get("species"), "n_samples": s.get("number_of_samples"), "assay": s.get("analysis_type"), "url": f"https://www.metabolomicsworkbench.org/data/DRCCMetadata.php?Mode=Study&StudyID={a}"})
    log.info(f"repository records for accessions seen only in papers: {len(rows)} of {len(accs)}")
    return pd.concat([H, pd.DataFrame(rows)], ignore_index=True)


def enrich_proteomexchange(H: pd.DataFrame) -> pd.DataFrame:
    """Description, species and instruments of ProteomeXchange and MassIVE datasets from their PROXI record (the
    search listing gives only the title)."""
    H = H.copy()
    m = H["source"].isin(["ProteomeXchange", "MassIVE"]) | (H["accession"].astype(str).str.startswith(("PXD", "MSV")) & H["description"].fillna("").eq(""))
    for i in H.index[m]:
        acc = str(H.at[i, "accession"])
        r = get(f"https://proteomecentral.proteomexchange.org/api/proxi/v0.1/datasets/{acc}")
        if not r:
            continue
        H.at[i, "description"] = str(r.get("description", ""))[:2000]
        sp = [t.get("value", t.get("name", "")) for g in r.get("species", []) for t in g.get("terms", []) if t.get("name") in ("taxonomy: scientific name", "NEWT")] if r.get("species") else []
        if sp:
            H.at[i, "organism"] = "; ".join(sp)
        ins = [t.get("name", "") for g in r.get("instruments", []) for t in ([g] if "name" in g else g.get("terms", []))] if r.get("instruments") else []
        if ins:
            H.at[i, "platform"] = "; ".join(ins)[:200]
    log.info(f"ProteomeXchange and MassIVE records read: {int(m.sum())}")
    return H


def canonical(acc: str) -> str:
    """One key per record: 10.5281/zenodo.N -> zenodo.N; a Mendeley Data DOI without its version."""
    acc = str(acc)
    m = re.fullmatch(r"10\.5281/zenodo\.(\d+)", acc, re.I)
    if m:
        return f"zenodo.{m.group(1)}"
    m = re.fullmatch(r"(10\.17632/[a-z0-9]+)(?:\.\d+)?", acc, re.I)
    return m.group(1).lower() if m else acc


def build_candidates(H: pd.DataFrame) -> pd.DataFrame:
    """One row per data accession. Hits found only through a paper add that paper; an ArrayExpress copy of a GEO
    series (E-GEOD) and a BioProject with a GEO series are kept as duplicates of the GEO series in the screen."""
    H = H.copy()
    H["accession"] = H["accession"].astype(str).map(canonical)
    data = H[H["source"] != "paper"].copy()
    papers = H[H["source"] == "paper"].copy()
    rows = []
    for acc, d in data.groupby("accession", sort=False):
        r = d.iloc[0].to_dict()
        r["sources"] = "; ".join(sorted(set(d["source"]))) + ("; paper" if acc in set(papers["accession"]) else "")
        r["queries"] = "; ".join(sorted(set(d["query"].astype(str))))
        rows.append(r)
    have = {r["accession"] for r in rows}
    for acc, d in papers.groupby("accession", sort=False):  # accessions seen only in papers
        if acc in have:
            continue
        r = d.iloc[0].to_dict()
        r["sources"], r["queries"] = "paper", "Europe PMC full text"
        if acc.startswith("GSE"):
            r["source"] = "GEO"
        elif acc.startswith(("PXD", "MSV")):
            r["source"] = "PRIDE" if acc.startswith("PXD") else "MassIVE"
        elif acc.startswith("ST"):
            r["source"] = "Metabolomics Workbench"
        elif acc.startswith("10.5281/zenodo"):
            r["source"] = "Zenodo"
        elif acc.startswith("10.17632"):
            r["source"] = "Mendeley Data"
        elif acc.startswith(("EGAS", "EGAD", "syn")):
            r["source"] = "controlled access archive"
        else:
            r["source"] = "other"
        rows.append(r)
    C = pd.DataFrame(rows)
    if "bioproject" in data:
        bp2geo = data.dropna(subset=["bioproject"]).groupby("bioproject")["accession"].apply(lambda v: ";".join(sorted(set(v))))
        C["linked_geo"] = [bp2geo.get(a, r) if s == "BioProject" else r for a, s, r in zip(C["accession"], C["source"], C.get("linked_geo", pd.Series([None] * len(C))))]
    pm = papers.groupby("accession").agg(paper_pmids=("pmids", lambda v: ";".join(sorted({str(x) for x in v if pd.notna(x)}))),
                                          paper_dois=("doi_paper", lambda v: ";".join(sorted({str(x) for x in v if pd.notna(x)}))),
                                          paper_titles=("title", lambda v: " | ".join(sorted({str(x) for x in v if pd.notna(x)}))[:500]))
    C = C.merge(pm, left_on="accession", right_index=True, how="left")
    C["paper_pmid"] = [";".join(sorted({x for x in f"{a};{b}".split(";") if x and x != "nan"})) for a, b in zip(C.get("pmids", ""), C["paper_pmids"].fillna(""))]
    C["paper_doi"] = C["paper_dois"].fillna("")
    txt = (C["title"].fillna("") + " " + C["description"].fillna("") + " " + C["paper_titles"].fillna(""))
    C["layer"] = [_layer(r) for r in C.to_dict("records")]
    C["tissue"] = [_tissue(t) for t in txt]
    C["access_type"] = [_access(a, s) for a, s in zip(C["accession"], C["source"])]
    C["deposited"] = [_deposited(r) for r in C.to_dict("records")]
    C["labels_where"] = [_labels(r) for r in C.to_dict("records")]
    cols = ["accession", "source", "sources", "title", "paper_pmid", "paper_doi", "layer", "tissue", "platform", "n_samples", "access_type", "deposited", "labels_where",
            "organism", "assay", "date", "linked_geo", "url", "paper_titles", "description"]
    return C.reindex(columns=cols)


def search() -> None:
    hits = []
    for f in (search_geo, search_bioproject, search_arrayexpress, search_pride, search_metabolomics_workbench, search_zenodo, search_mendeley, search_europepmc):
        try:
            hits += f()
        except Exception as e:  # noqa: BLE001
            log.error(f"{f.__name__} failed: {e}")
    H = pd.DataFrame(hits)
    for acc in FROM_MEMORY:  # candidates from memory, to verify
        if acc not in set(H.get("accession", pd.Series(dtype=str)).astype(str)):
            H = pd.concat([H, pd.DataFrame([{"source": "memory", "query": "from memory", "accession": acc, "title": "candidate from memory (to verify)"}])], ignore_index=True)
    H = enrich_paper_accessions(H)
    H = enrich_proteomexchange(H)
    WORK.mkdir(parents=True, exist_ok=True)
    H.to_csv(WORK / "search_hits.csv", index=False)
    log.info(f"search: {len(H)} hits written to {WORK / 'search_hits.csv'}")
    C = build_candidates(H)
    DOCS.mkdir(exist_ok=True)
    C.drop(columns=["description"]).to_csv(CANDIDATES, index=False)
    C.to_csv(WORK / "candidates_full.csv", index=False)
    log.info(f"search: {len(C)} candidates in {CANDIDATES} ({C['layer'].value_counts().to_dict()})")


# ------------------------------------------------------------------------------------------------ screen (automatic part)
AUTO_EXCLUDE = [  # (pattern on title, description and paper titles, reason); checked in this order
    (r"single[- ]cell|scrna|sc-rna|\bscrnaseq|cite-seq|10x genomics|snrna|single-nucle", "single cell data"),
    (r"\bmice\b|\bmouse\b|murine|hamster|macaque|ferret|non-human primate|rhesus", "not human"),
    (r"cell line|organoid|in vitro|infected cells|calu-3|vero e6|caco-2|hek293|ipsc|a549|ex vivo stimulat|stimulated with", "in vitro or cell culture"),
    (r"\bvaccin", "vaccination study"),
    (r"\bchild|children|pediatric|paediatric|mis-c|kawasaki|infant|neonat", "not adults"),
    (r"methylation|methylome|atac-seq|chip-seq|t cell receptor|tcr repertoire|bcr repertoire|\bmirna|microrna|small rna|circrna|exosom|extracellular vesicle", "other assay or analyte"),
]
NON_BLOOD = r"nasal|nasopharyn|swab|bronchoalveolar|\bbal\b|lung tissue|autopsy|post-mortem|postmortem|saliva|urine|cerebrospinal|placenta|tissue biops|kidney|brain|heart|airway epitheli"
POST_ACUTE = r"long covid|post[- ]covid|\bpasc\b|post-acute|postacute|sequelae"
ACUTE = r"acute|hospitali|intensive care|\bicu\b|severity|severe|critical|admission|admitted|emergency"


def auto_screen(C: pd.DataFrame) -> pd.DataFrame:
    """Automatic decisions from the candidate metadata; a candidate that passes gets "needs decision" until it is
    curated (CURATED)."""
    out = []
    geo = set(C.loc[C["source"] == "GEO", "accession"])
    for r in C.to_dict("records"):
        txt = f"{r.get('title', '')} {r.get('description', '')} {r.get('paper_titles', '')}".lower()
        acc = str(r["accession"])
        dec, why = "needs decision", "passes the automatic screen; not read by hand"
        org = str(r.get("organism") or "")
        if acc.startswith("E-GEOD-") and f"GSE{acc.split('-')[-1]}" in geo:
            dec, why = "exclude", f"duplicate: ArrayExpress copy of GSE{acc.split('-')[-1]}"
        elif r["source"] == "BioProject" and str(r.get("linked_geo") or "") not in ("", "nan"):
            dec, why = "exclude", f"duplicate: BioProject of {r['linked_geo']}"
        elif r["access_type"] != "open":
            dec, why = "exclude", "not open access (application needed)"
        elif org and org != "nan" and "homo sapiens" not in org.lower() and "human" not in org.lower():
            dec, why = "exclude", f"not human ({org})"
        elif r["layer"] == "unknown":
            dec, why = "exclude", "not an RNA, protein or metabolite dataset (by its record)"
        elif not re.search(r"covid|sars-cov|ncov|coronavirus", txt):
            dec, why = "exclude", "not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record)"
        else:
            for pat, reason in AUTO_EXCLUDE:
                if re.search(pat, txt):
                    dec, why = "exclude", reason
                    break
            else:
                blood = r["tissue"] != "not stated"
                if not blood and re.search(NON_BLOOD, txt):
                    dec, why = "exclude", "not blood"
                elif re.search(POST_ACUTE, txt) and not re.search(ACUTE, txt):
                    dec, why = "exclude", "post-acute only (no acute samples by its record)"
                elif pd.notna(r.get("n_samples")) and str(r.get("n_samples")).isdigit() and int(r["n_samples"]) < 20:
                    dec, why = "exclude", f"{r['n_samples']} samples: cannot reach 10 participants per class"
                elif r["layer"] == "RNA" and r["tissue"] in ("plasma", "serum", "plasma; serum"):
                    dec, why = "exclude", "RNA from plasma or serum (cell-free RNA), not whole blood or PBMC"
                elif r["layer"] in ("proteins", "metabolites") and not re.search(r"plasma|serum|blood", r["tissue"]):
                    dec, why = ("exclude", "not blood") if re.search(NON_BLOOD, txt) else ("needs decision", "tissue not stated in the record; not read by hand")
        out.append({"accession": acc, "decision": dec, "reason": why})
    return pd.DataFrame(out)


PROTEIN_NOT_PLASMA = (r"interactom|protease|\bspike\b|nsp\d|orf\d|replication|antibod|epitope|immunopeptid|glycan|glycosyl|virion|viral protein|crosslink|"
                      r"semen|breast ?milk|breastmilk|cerebrospinal|\bcsf\b|saliva|urine|swab|sputum|tissue|biops|autops|lung|liver|kidney|brain|ovar|placenta|adipose|"
                      r"platelet|neutrophil|monocyte|macrophage|lymphocyte|b cells|t cells|\bpbmc|cell line|organoid|in vitro|infected cells|hamster|mouse|exosom|vesicle|"
                      r"\bgut\b|stool|feces|faecal|intestin")
COVID_VALUE = r"covid|sars-cov-2|sars cov 2|sarscov2|ncov|positive"
OTHER_VALUE = r"healthy|control|donor|uninfected|negative|non-covid|noncovid|normal|baseline|pre-pandemic"
SEVERITY_VALUE = r"icu|intensive|critical|severe|ventilat|intubat|who|mild|moderate|fatal|death|died"


def overlaps(included: list) -> dict:
    """GEO series read in the screen that share GSM accessions or BioSamples with an included GEO cohort; returns
    {other accession: (included cohort, number shared)} where the included cohort is the larger of the two."""
    out = {}
    files = sorted((WORK / "geo").glob("*_samples.csv"))
    ids = {}
    for f in files:
        try:
            d = pd.read_csv(f, dtype=str)
        except Exception:  # noqa: BLE001 (empty file: no series matrix)
            continue
        s = set(d.get("geo_accession", pd.Series(dtype=str)).dropna()) | set(d.get("biosample", pd.Series(dtype=str)).dropna())
        ids[f.name.split("_")[0]] = s
    for cid in included:
        acc = COHORTS[cid]["accession"]
        mine = ids.get(acc, set())
        for other, s in ids.items():
            if other == acc or not mine:
                continue
            n = len(mine & s)
            if n and len(s) <= len(mine):
                out[other] = (cid, n)
    return out


def characteristics_estimate(d: pd.DataFrame) -> dict:
    """Automatic estimate from GEO sample characteristics, samples (not participants): the largest number of samples
    with a COVID-like value and with a control-like value in one column, and whether a column holds severity words.
    Used only to rank and to exclude candidates that cannot reach 10 per class."""
    best = {"est_covid": 0, "est_other": 0, "est_column": "", "severity_column": ""}
    cols = [c for c in d.columns if c.startswith("ch_") or c in ("source_name_ch1", "title")]
    for c in cols:
        v = d[c].astype(str).str.lower()
        cov = v.str.contains(COVID_VALUE) & ~v.str.contains(r"non-covid|noncovid|negative|non covid")
        oth = v.str.contains(OTHER_VALUE) & ~cov
        if min(cov.sum(), oth.sum()) > min(best["est_covid"], best["est_other"]):
            best.update(est_covid=int(cov.sum()), est_other=int(oth.sum()), est_column=c)
        if c != "title" and v.str.contains(SEVERITY_VALUE).sum() >= 10 and v.nunique() <= 20:
            best["severity_column"] = c
    return best


def screen() -> None:
    C = pd.read_csv(WORK / "candidates_full.csv", dtype=str)
    A = auto_screen(C).set_index("accession")
    # GEO candidates still open: read their sample characteristics (metadata only)
    open_geo = [a for a in A.index[A["decision"] == "needs decision"] if a.startswith("GSE")]
    from concurrent.futures import ThreadPoolExecutor

    def _read(a):
        try:
            geo_samples(a)
        except Exception as e:  # noqa: BLE001
            log.warning(f"{a}: series matrix not read ({e})")

    with ThreadPoolExecutor(8) as ex:
        list(ex.map(_read, open_geo))
    est = {}
    for a in open_geo:
        try:
            d = geo_samples(a)
        except Exception as e:  # noqa: BLE001
            log.warning(f"{a}: series matrix not read ({e})")
            continue
        e = characteristics_estimate(d)
        e["n_samples_matrix"] = len(d)
        e["tissues"] = "; ".join(sorted(set(d.get("source_name_ch1", pd.Series(dtype=str)).astype(str))))[:200]
        est[a] = e
        if min(e["est_covid"], e["est_other"]) < 10 and not e["severity_column"]:
            A.loc[a, ["decision", "reason"]] = ["exclude", f"no characteristic splits 10 COVID samples from 10 others and no severity field (automatic; best {e['est_column']}: {e['est_covid']} against {e['est_other']})"]
    E = pd.DataFrame(est).T
    E.index.name = "accession"
    E.to_csv(WORK / "screen_estimates.csv")
    # proteins: records about cells, tissues, other fluids or viral proteins are not plasma or serum proteomes
    for a in A.index[A["decision"] == "needs decision"]:
        r = C.set_index("accession").loc[a] if a in set(C["accession"]) else None
        if r is None or (r["layer"] if isinstance(r, pd.Series) else r.iloc[0]["layer"]) != "proteins":
            continue
        r = r if isinstance(r, pd.Series) else r.iloc[0]
        txt = f"{r.get('title', '')} {r.get('description', '')}".lower()
        if re.search(PROTEIN_NOT_PLASMA, txt):
            A.loc[a, ["decision", "reason"]] = ["exclude", "not a plasma or serum proteome of patients (cells, tissue, other fluid or viral proteins)"]
        elif re.search(POST_ACUTE, txt) and not re.search(r"\bacute\b|hospitali|severity|severe", txt):
            A.loc[a, ["decision", "reason"]] = ["exclude", "post-acute only (no acute samples by its record)"]
    # mass spectrometry projects still open: without result tables there is no abundance matrix to use
    for a in [a for a in A.index[A["decision"] == "needs decision"] if a.startswith("PXD")]:
        fs = get(f"https://www.ebi.ac.uk/pride/ws/archive/v2/projects/{a}/files", {"pageSize": 1000}) or []
        names = [f.get("fileName", "") for f in fs]
        tables = [n for n in names if re.search(r"\.(txt|tsv|csv|xlsx?)$", n, re.I) and not re.search(r"\.raw\.|mgf|mzml", n, re.I)]
        sdrf = [n for n in names if "sdrf" in n.lower()]
        if not names:
            A.loc[a, ["decision", "reason"]] = ["exclude", "no files listed in PRIDE (not public in PRIDE, or hosted elsewhere)"]
        elif not tables:
            A.loc[a, ["decision", "reason"]] = ["exclude", f"MS raw files only ({len(names)} files); no abundance matrix deposited"]
        else:
            A.loc[a, "reason"] = f"abundance tables deposited ({len(tables)}{', SDRF' if sdrf else ''}); per sample labels to take from the paper; not read by hand"
    # RNA candidates still open once the cap of 10 RNA cohorts is reached are left unread
    n_rna = sum(1 for cur in CURATED.values() if cur["decision"] == "include" for c in cur.get("cohorts", []) if COHORTS[c]["layer"] == "RNA")
    if n_rna >= 10:
        rna_open = [a for a in A.index[A["decision"] == "needs decision"] if a in set(C.loc[C["layer"] == "RNA", "accession"]) and a not in CURATED]
        A.loc[rna_open, "reason"] = A.loc[rna_open, "reason"] + "; the cap of 10 RNA cohorts was reached, so it was not read further"
    met_open = [a for a in A.index[A["decision"] == "needs decision"] if a in set(C.loc[C["layer"] == "metabolites", "accession"]) and a not in CURATED]
    A.loc[met_open, "reason"] = A.loc[met_open, "reason"] + "; metabolite cohorts were not read further (the targets of this sweep were 10 RNA and 6 protein cohorts)"
    prot_open = [a for a in A.index[A["decision"] == "needs decision"] if a in set(C.loc[C["layer"] == "proteins", "accession"]) and a not in CURATED]
    A.loc[prot_open, "reason"] = A.loc[prot_open, "reason"] + "; not read further after the protein candidates with open matrices and labels were read"
    # curated decisions (read by hand) override the automatic ones
    for acc, cur in CURATED.items():
        A.loc[acc, ["decision", "reason"]] = [cur["decision"], cur["reason"]]
    rows = []
    for acc, cur in CURATED.items():
        for cid in cur.get("cohorts", []):
            spec = COHORTS[cid]
            L = cohort_labels(cid)
            c = class_counts(L)
            rows.append({"accession": acc, "cohort": cid, "layer": spec["layer"], **c, "contrasts": ";".join(contrasts_available(c)), "deposited_matrix": spec.get("matrix_kind", "")})
    K = pd.DataFrame(rows)
    # the include cap: 10 RNA and 6 protein cohorts, by participants and by availability of a processed matrix
    if len(K):
        K["n_contrast"] = K.apply(lambda r: max(min(r["infection_acute_covid"], r["infection_never_infected"]), min(r["severity_severe_critical"], r["severity_mild_moderate"])), axis=1)
        K["rank_key"] = K["deposited_matrix"].ne("raw").astype(int) * 10**6 + K["n_participants"]
        for layer, cap in (("RNA", 10), ("proteins", 6), ("metabolites", 10**6)):
            k = K[(K["layer"] == layer) & (K["contrasts"] != "")].sort_values("rank_key", ascending=False)
            over = k.index[cap:]
            for i in over:
                acc = K.loc[i, "accession"]
                A.loc[acc, ["decision", "reason"]] = ["exclude", f"not needed: the cap of {cap} {layer} cohorts was reached by larger cohorts or ones with a processed matrix"]
        K.to_csv(WORK / "screen_counts.csv", index=False)
    # overlap: GEO series sharing samples (GSM) or BioSamples with an included cohort; the larger one is kept
    ov = overlaps([c for cur in CURATED.values() if cur["decision"] == "include" for c in cur.get("cohorts", []) if COHORTS[c].get("meta", "geo") == "geo"])
    for other, (kept, n) in ov.items():
        if A.loc[other, "decision"] != "include" if other in A.index else False:
            A.loc[other, ["decision", "reason"]] = ["exclude", f"duplicate: shares {n} samples or BioSamples with {kept} (kept, larger)"]
    pd.DataFrame([{"accession": k, "duplicate_of": v[0], "n_shared": v[1]} for k, v in ov.items()]).to_csv(WORK / "overlaps.csv", index=False)
    D = C.drop(columns=["description"]).merge(A.reset_index(), on="accession", how="left").merge(E.reset_index(), on="accession", how="left")
    if len(K):
        agg = K.groupby("accession").agg(cohorts=("cohort", ";".join), contrasts=("contrasts", lambda v: " | ".join(f"{c}" for c in v)),
                                          counts=("cohort", lambda v: ""))
        cnt = K.assign(txt=lambda t: t["cohort"] + ": infection " + t["infection_acute_covid"].astype(str) + " acute / " + t["infection_never_infected"].astype(str)
                       + " never infected, severity " + t["severity_severe_critical"].astype(str) + " severe or critical / " + t["severity_mild_moderate"].astype(str) + " mild or moderate (participants)")
        agg["counts"] = cnt.groupby("accession")["txt"].apply(" | ".join)
        D = D.merge(agg.reset_index(), on="accession", how="left")
    for acc, cur in CURATED.items():
        for k_ in ("duplicate_of", "notes"):
            if cur.get(k_):
                D.loc[D["accession"] == acc, k_] = cur[k_]
    D.to_csv(CANDIDATES, index=False)
    s = D["decision"].value_counts().to_dict()
    log.info(f"screen: {s}; included cohorts: {K.loc[K['accession'].isin(D.loc[D['decision'] == 'include', 'accession']), 'cohort'].tolist() if len(K) else []}")


# ------------------------------------------------------------------------------------------------ GEO metadata
def geo_dir(gse: str) -> str:
    return f"https://ftp.ncbi.nlm.nih.gov/geo/series/{gse[:-3]}nnn/{gse}"


def geo_samples(gse: str) -> pd.DataFrame:
    """Per sample metadata of a GEO series from its series matrix file(s) (metadata lines only); cached as
    data/acute_covid/geo/<gse>_samples.csv. Characteristics are split into one column per key (``ch_<key>``)."""
    out = WORK / "geo" / f"{gse}_samples.csv"
    if out.exists():
        return pd.read_csv(out, dtype=str)
    out.parent.mkdir(parents=True, exist_ok=True)
    listing = get(geo_dir(gse) + "/matrix/", kind="text") or ""
    files = sorted(set(re.findall(rf'href="({gse}[^"]*_series_matrix\.txt\.gz)"', listing)))
    frames = []
    for fn in files:
        lines = []
        with requests.get(f"{geo_dir(gse)}/matrix/{fn}", stream=True, timeout=300) as r:  # metadata lines only
            r.raise_for_status()
            with gzip.GzipFile(fileobj=r.raw) as gz:
                for raw_line in gz:
                    ln = raw_line.decode("utf-8", "replace").rstrip("\n").rstrip("\r")
                    lines.append(ln)
                    if ln.startswith("!series_matrix_table_begin"):
                        break
        recs: dict = {}
        for ln in lines:
            if ln.startswith("!series_matrix_table_begin"):
                break
            if not ln.startswith("!Sample_"):
                continue
            key, *vals = ln.split("\t")
            vals = [v.strip('"') for v in vals]
            key = key[len("!Sample_"):]
            if key == "characteristics_ch1":
                # one line per characteristic; key from "key: value"
                for i, v in enumerate(vals):
                    k, _, val = v.partition(": ")
                    recs.setdefault(f"ch_{k.strip().lower().replace(' ', '_')}" if val else f"ch_{len([x for x in recs if x.startswith('ch_')])}", {})[i] = val if val else v
            elif key == "relation":
                for i, v in enumerate(vals):
                    if v.startswith("BioSample"):
                        recs.setdefault("biosample", {})[i] = v.split("/")[-1]
                    elif v.startswith("SRA"):
                        recs.setdefault("sra", {})[i] = v.split("=")[-1]
            elif key in ("title", "geo_accession", "source_name_ch1", "description", "platform_id", "molecule_ch1", "organism_ch1", "submission_date"):
                recs[key] = dict(enumerate(vals))
        frames.append(pd.DataFrame(recs))
    df = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    df.to_csv(out, index=False)
    return df


# ------------------------------------------------------------------------------------------------ labels
# Classes. infection: acute_covid, never_infected, other_disease (other infections or critical illness: their own class,
# never a control), excluded (for example convalescent). severity: severe_critical (intensive care, mechanical
# ventilation, WHO scale 6 or more, or death), mild_moderate, or not assigned.
INFECTION = ("acute_covid", "never_infected", "other_disease", "excluded")
SEVERITY = ("severe_critical", "mild_moderate")
LABEL_COLUMNS = ["sample_id", "participant_id", "order", "label_original", "infection_class", "severity_original", "severity_class",
                 "days", "time_anchor", "in_acute_window", "sex", "age", "batch_fields"]
LABELERS: dict = {}


def labeler(cohort: str):
    def deco(f):
        LABELERS[cohort] = f
        return f
    return deco


def _norm_sex(v) -> str:
    v = str(v).strip().lower()
    return {"m": "male", "male": "male", "f": "female", "female": "female"}.get(v, "unknown")


def _num(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return np.nan


def one_per_participant(L: pd.DataFrame) -> pd.DataFrame:
    """One sample per participant: the first acute COVID sample when the participant has one, otherwise the first
    sample (by ``order``)."""
    L = L.copy()
    L["_acute"] = (L["infection_class"] == "acute_covid").astype(int)
    L = L.sort_values(["participant_id", "_acute", "order"], ascending=[True, False, True], kind="mergesort")
    return L.drop_duplicates("participant_id").drop(columns="_acute").sort_index()


def class_counts(L: pd.DataFrame) -> dict:
    """Participants per class and contrast after one sample per participant."""
    K = one_per_participant(L)
    inf = K["infection_class"].value_counts()
    sev = K.loc[K["infection_class"] == "acute_covid", "severity_class"].value_counts()
    return {"n_samples": len(L), "n_participants": int(K["participant_id"].nunique()),
            "infection_acute_covid": int(inf.get("acute_covid", 0)), "infection_never_infected": int(inf.get("never_infected", 0)),
            "other_disease": int(inf.get("other_disease", 0)), "excluded": int(inf.get("excluded", 0)),
            "severity_severe_critical": int(sev.get("severe_critical", 0)), "severity_mild_moderate": int(sev.get("mild_moderate", 0))}


def contrasts_available(c: dict, minimum: int = 10) -> list:
    out = []
    if min(c["infection_acute_covid"], c["infection_never_infected"]) >= minimum:
        out.append("infection")
    if min(c["severity_severe_critical"], c["severity_mild_moderate"]) >= minimum:
        out.append("severity")
    return out


def mapping_table(L: pd.DataFrame) -> pd.DataFrame:
    """Original text to class, one row per distinct combination, with the number of samples."""
    return (L.groupby(["label_original", "infection_class", "severity_original", "severity_class"], dropna=False).size()
            .rename("n_samples").reset_index())


def _frame(**cols) -> pd.DataFrame:
    df = pd.DataFrame(cols)
    for c in LABEL_COLUMNS:
        if c not in df:
            df[c] = np.nan
    return df[LABEL_COLUMNS + [c for c in df.columns if c not in LABEL_COLUMNS]]


@labeler("GSE161731")
def _gse161731(d: pd.DataFrame) -> pd.DataFrame:
    # McClain 2021: whole blood; cohort COVID-19 / healthy / CoV other / Bacterial / Influenza; time since onset early
    # (up to 10 days), middle (11 to 21) or late (over 21) for the COVID-19 samples (paper Methods)
    inf = d["ch_cohort"].map({"COVID-19": "acute_covid", "healthy": "never_infected"}).fillna("other_disease")
    order = d["ch_time_since_onset"].map({"early": 0, "middle": 1, "late": 2}).fillna(0)
    window = d["ch_time_since_onset"].map({"early": True, "middle": True, "late": np.nan})
    return _frame(sample_id=d["geo_accession"], participant_id="GSE161731:" + d["ch_subject_id"], order=order, label_original="cohort=" + d["ch_cohort"] + "; time_since_onset=" + d["ch_time_since_onset"],
                  infection_class=inf, severity_original="hospitalized=" + d["ch_hospitalized"], severity_class=np.nan, days=np.nan,
                  time_anchor=np.where(inf == "acute_covid", "symptom_onset", None), in_acute_window=np.where(inf == "acute_covid", window, np.nan),
                  sex=d["ch_gender"].map(_norm_sex), age=d["ch_age"].map(_num), batch_fields="batch=" + d["ch_batch"], batch=d["ch_batch"])


@labeler("GSE152641")
def _gse152641(d: pd.DataFrame) -> pd.DataFrame:
    # Thair 2021: whole blood of COVID-19 patients (hospitalised, Italy and the US; paper Methods) and healthy controls
    inf = d["ch_disease"].map({"COVID19": "acute_covid", "Healthy control": "never_infected"})
    return _frame(sample_id=d["geo_accession"], participant_id="GSE152641:" + d["geo_accession"], order=0, label_original="disease=" + d["ch_disease"],
                  infection_class=inf, severity_original=np.nan, severity_class=np.nan, days=np.nan, time_anchor=None,
                  in_acute_window=np.where(inf == "acute_covid", True, np.nan), sex=d["ch_sex"].map(_norm_sex), age=d["ch_age"].map(_num), batch_fields="")


@labeler("GSE171110")
def _gse171110(d: pd.DataFrame) -> pd.DataFrame:
    # Levy 2021: whole blood of patients with severe COVID-19 admitted to intensive care, and healthy donors
    inf = d["ch_status"].map({"Severe Covid-19": "acute_covid", "Healthy": "never_infected"})
    return _frame(sample_id=d["geo_accession"], participant_id="GSE171110:" + d["geo_accession"], order=0, label_original="status=" + d["ch_status"],
                  infection_class=inf, severity_original=np.where(inf == "acute_covid", "Severe Covid-19", None), severity_class=np.where(inf == "acute_covid", "severe_critical", None),
                  days=np.nan, in_acute_window=np.where(inf == "acute_covid", True, np.nan), sex="unknown", age=np.nan, batch_fields="")


@labeler("GSE152418")
def _gse152418(d: pd.DataFrame) -> pd.DataFrame:
    # Arunachalam 2020: PBMC; Healthy, COVID-19 with severity Moderate, Severe or ICU, and one Convalescent sample.
    # Severity: ICU -> severe_critical; Moderate -> mild_moderate; Severe (supplemental oxygen, WHO 5) is neither.
    inf = d["ch_disease_state"].map({"COVID-19": "acute_covid", "Healthy": "never_infected", "Convalescent": "excluded"})
    sev = d["ch_severity"].map({"ICU": "severe_critical", "Moderate": "mild_moderate"})
    days = d["ch_days_post_symptom_onset"].map(_num)
    return _frame(sample_id=d["geo_accession"], participant_id="GSE152418:" + d["geo_accession"], order=0, label_original="disease_state=" + d["ch_disease_state"],
                  infection_class=inf, severity_original="severity=" + d["ch_severity"], severity_class=np.where(inf == "acute_covid", sev, None), days=days,
                  time_anchor=np.where(inf == "acute_covid", "symptom_onset", None), in_acute_window=np.where(inf == "acute_covid", days <= 30, np.nan),
                  sex=d["ch_gender"].map(_norm_sex), age=np.nan, batch_fields="")


@labeler("GSE172114")
def _gse172114(d: pd.DataFrame) -> pd.DataFrame:
    # Carapito 2022: whole blood of hospitalised COVID-19 patients, critical (intensive care) against non-critical
    sev = d["ch_disease_state"].map({"Critical": "severe_critical", "Non-critical": "mild_moderate"})
    return _frame(sample_id=d["geo_accession"], participant_id="GSE172114:" + d["geo_accession"], order=0, label_original="disease=" + d["ch_disease"],
                  infection_class="acute_covid", severity_original="disease_state=" + d["ch_disease_state"], severity_class=sev, days=np.nan,
                  time_anchor="admission", in_acute_window=True, sex="unknown", age=np.nan, batch_fields="")


@labeler("GSE157103")
def _gse157103(d: pd.DataFrame) -> pd.DataFrame:
    # Overmyer 2021: leukocytes from whole blood of hospitalised patients with and without COVID-19 (non-COVID: other
    # hospitalised illness, own class); ICU admission -> severe_critical, otherwise mild_moderate
    inf = d["ch_disease_state"].map({"COVID-19": "acute_covid", "non-COVID-19": "other_disease"})
    sev = d["ch_icu"].map({"yes": "severe_critical", "no": "mild_moderate"})
    return _frame(sample_id=d["geo_accession"], participant_id="GSE157103:" + d["title"].str.replace(r"\W", "", regex=True), order=0, label_original="disease_state=" + d["ch_disease_state"],
                  infection_class=inf, severity_original="icu=" + d["ch_icu"] + "; mechanical_ventilation=" + d["ch_mechanical_ventilation"],
                  severity_class=np.where(inf == "acute_covid", sev, None), days=np.nan, time_anchor="admission", in_acute_window=True,
                  sex=d["ch_sex"].map(_norm_sex), age=d["ch_age_(years)"].map(_num), batch_fields="")


@labeler("GSE155454")
def _gse155454(d: pd.DataFrame) -> pd.DataFrame:
    # Fong 2022 (no PMID in GEO): whole blood, COVID-19 patients (several draws) and healthy controls
    inf = d["ch_diagnosis"].map({"COVID-19": "acute_covid", "Healthy": "never_infected"})
    return _frame(sample_id=d["geo_accession"], participant_id="GSE155454:" + d["ch_patient_id"], order=d["title"].str.extract(r"(\d+)$")[0].map(_num).fillna(0),
                  label_original="diagnosis=" + d["ch_diagnosis"], infection_class=inf, severity_original="severity=" + d["ch_severity"], severity_class=np.nan,
                  days=np.nan, in_acute_window=np.where(inf == "acute_covid", True, np.nan), sex="unknown", age=np.nan, batch_fields="")


@labeler("GSE163151")
def _gse163151(d: pd.DataFrame) -> pd.DataFrame:
    # Ng 2021: nasal swabs and whole blood; only the whole blood samples are considered
    d = d[d["source_name_ch1"] == "Whole Blood"]
    inf = d["ch_disease_state"].map({"COVID-19": "acute_covid", "Donor control": "never_infected"}).fillna("other_disease")
    return _frame(sample_id=d["geo_accession"], participant_id="GSE163151:" + d["title"], order=0, label_original="disease_state=" + d["ch_disease_state"],
                  infection_class=inf, severity_original="severity=" + d["ch_disease_severity"], severity_class=np.nan, days=np.nan, in_acute_window=np.nan,
                  sex="unknown", age=np.nan, batch_fields="")


@labeler("GSE198449")
def _gse198449(d: pd.DataFrame) -> pd.DataFrame:
    # CHARM (Marine recruits, PAXgene blood, repeated draws): acute_covid = a sample with a detected SARS-CoV-2 PCR;
    # never_infected = a participant whose every sample is group Control with PCR not detected; other samples excluded
    # (pre-infection, post-infection, exposed, immune, inconclusive). Days are days since t0 of the study.
    pid = "GSE198449:" + d["ch_participant_id"]
    pcr = d["ch_pcr_test_for_sars-cov-2"]
    never = d.groupby("ch_participant_id").apply(lambda g: bool((g["ch_group"] == "Control").all() and (g["ch_pcr_test_for_sars-cov-2"] == "Not").all()))
    inf = np.where(pcr == "Detected", "acute_covid", np.where(d["ch_participant_id"].map(never), "never_infected", "excluded"))
    t = d["ch_sample_collection_time_point_(days_since_t0)"].map(_num)
    return _frame(sample_id=d["geo_accession"], participant_id=pid, order=t, label_original="group=" + d["ch_group"] + "; pcr=" + pcr,
                  infection_class=inf, severity_original="group=" + d["ch_group"], severity_class=np.nan, days=np.where(inf == "acute_covid", 0, np.nan),
                  time_anchor=np.where(inf == "acute_covid", "positive_test", None), in_acute_window=np.where(inf == "acute_covid", True, np.nan),
                  sex=d["ch_sex"].map(_norm_sex), age=d["ch_age"].map(_num),
                  batch_fields="plate=" + d["ch_rna-seq_plate"] + "; enrollment_batch=" + d["ch_enrollment_batch"] + "; t0=" + d["ch_enrollment_date_t0"],
                  plate=d["ch_rna-seq_plate"], enrollment_batch=d["ch_enrollment_batch"], collection_date=d["ch_enrollment_date_t0"])


# ------------------------------------------------------------------------------------------------ downloads
RECORDS = WORK / "records"  # files of repository records (metadata and data), shared by screen and fetch


def sha256_of(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def md5_of(path: Path) -> str:
    h = hashlib.md5()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def download(url: str, dest: Path, md5: str = None, sha256: str = None, tries: int = 3) -> dict:
    """Resumable download (HTTP range) with three retries; a repository checksum is verified when given. A finished
    file is recorded in data/acute_covid/files.csv with its size and SHA-256 and is not downloaded again."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    reg = WORK / "files.csv"
    R = pd.read_csv(reg, dtype=str) if reg.exists() else pd.DataFrame(columns=["path", "url", "size", "sha256", "repository_checksum", "checksum_ok", "date"])
    rel = str(dest.relative_to(ROOT))
    if dest.exists() and rel in set(R["path"]):
        return R[R["path"] == rel].iloc[-1].to_dict()
    part = dest.with_name(dest.name + ".part")
    err = None
    for i in range(tries):
        try:
            have = part.stat().st_size if part.exists() else 0
            hdr = {"Range": f"bytes={have}-"} if have else {}
            with requests.get(url, stream=True, timeout=300, headers=hdr) as r:
                if r.status_code == 416:  # already complete
                    pass
                else:
                    r.raise_for_status()
                    mode = "ab" if have and r.status_code == 206 else "wb"
                    with open(part, mode) as f:
                        for b in r.iter_content(1 << 20):
                            f.write(b)
            got_md5 = md5_of(part) if md5 else None
            got_sha = sha256_of(part)
            ok = (md5 is None or got_md5 == md5) and (sha256 is None or got_sha == sha256)
            if not ok:
                part.unlink()
                raise ValueError(f"checksum mismatch ({'md5 ' + str(got_md5) if md5 else 'sha256 ' + got_sha})")
            part.rename(dest)
            row = {"path": rel, "url": url, "size": dest.stat().st_size, "sha256": got_sha, "repository_checksum": f"md5:{md5}" if md5 else (f"sha256:{sha256}" if sha256 else ""),
                   "checksum_ok": bool(md5 or sha256) or "", "date": time.strftime("%Y-%m-%d %H:%M:%S")}
            pd.concat([R, pd.DataFrame([row])], ignore_index=True).to_csv(reg, index=False)
            log.info(f"downloaded {rel} ({row['size']} bytes{', checksum verified' if row['checksum_ok'] else ''})")
            return row
        except Exception as e:  # noqa: BLE001
            err = e
            log.warning(f"download attempt {i + 1} of {tries} failed: {url} ({e})")
            time.sleep(5 * (i + 1))
    raise RuntimeError(f"download failed after {tries} tries: {url} ({err})")


def mendeley_files(dataset: str) -> dict:
    m = get(f"https://data.mendeley.com/public-api/datasets/{dataset}")
    files = get(f"https://data.mendeley.com/public-api/datasets/{dataset}/files", {"folder_id": "root", "version": m["version"]}) or []
    return {f["filename"]: {"url": f["content_details"]["download_url"], "sha256": f["content_details"].get("sha256_hash"), "size": f["content_details"].get("size"),
                            "version": m["version"]} for f in files}


def zenodo_files(record: str) -> dict:
    r = get(f"https://zenodo.org/api/records/{record}")
    return {f["key"]: {"url": f["links"]["self"], "md5": f["checksum"].split(":")[-1] if f.get("checksum", "").startswith("md5:") else None, "size": f["size"]} for f in r.get("files", [])}


def record_file(source: str, record: str, filename: str, dest_dir: Path = None) -> Path:
    """A file of a Mendeley Data or Zenodo record, downloaded once (with its checksum) to data/acute_covid/records/."""
    files = mendeley_files(record) if source == "mendeley" else zenodo_files(record)
    f = files[filename]
    dest = (dest_dir or RECORDS / f"{source}_{record}") / filename
    download(f["url"], dest, md5=f.get("md5"), sha256=f.get("sha256"))
    return dest


def geo_suppl_file(gse: str, filename: str, dest_dir: Path) -> Path:
    dest = dest_dir / filename
    download(f"{geo_dir(gse)}/suppl/{filename}", dest)
    return dest


# ------------------------------------------------------------------------------------------------ non-GEO metadata
def meta_filbin() -> pd.DataFrame:
    p = record_file("mendeley", "nf853r8xsj", "Clinical Metadata.xlsx")
    return pd.read_excel(p, sheet_name="Subject-level metadata")


def meta_gisby_plasma() -> pd.DataFrame:
    return pd.read_csv(record_file("zenodo", "4287557", "plasma_sample_level.csv"))


def meta_combat() -> pd.DataFrame:
    import tarfile

    p = record_file("zenodo", "6120249", "CBD-KEY-CLINVAR.tar.gz")
    with tarfile.open(p) as t:
        return pd.read_csv(t.extractfile("CBD-KEY-CLINVAR/COMBAT_CLINVAR_for_processed.txt"), sep="\t")


@labeler("Filbin2021_MGH")
def _filbin(d: pd.DataFrame) -> pd.DataFrame:
    # Filbin 2021: emergency department patients with acute respiratory distress, COVID-19 positive (COVID 1) or
    # negative (COVID 0: other illness, own class); one plasma draw per participant at day 0 (enrolment, in the
    # emergency department). Acuity max: 1 died, 2 intubated or ventilated (WHO 7 to 9) -> severe_critical; 4
    # hospitalised without oxygen, 5 discharged -> mild_moderate; 3 hospitalised with supplemental oxygen (any type,
    # WHO 5 or 6) -> not assigned. Participants without a day 0 draw are left out.
    d = d[d["D0_draw"] == 1]
    inf = d["COVID"].map({1: "acute_covid", 0: "other_disease"})
    sev = d["Acuity max"].map({1: "severe_critical", 2: "severe_critical", 4: "mild_moderate", 5: "mild_moderate"})
    return _frame(sample_id=d["Public ID"].astype(str) + "_D0", participant_id="Filbin2021_MGH:" + d["Public ID"].astype(str), order=0,
                  label_original="COVID=" + d["COVID"].astype(str), infection_class=inf, severity_original="Acuity max=" + d["Acuity max"].astype(str),
                  severity_class=np.where(inf == "acute_covid", sev, None), days=0, time_anchor="admission", in_acute_window=True,
                  sex="unknown", age=np.nan, batch_fields="", age_category=d["Age cat"])


@labeler("Gisby2021_plasma")
def _gisby(d: pd.DataFrame) -> pd.DataFrame:
    # Gisby 2021 (eLife): haemodialysis patients, SARS-CoV-2 positive (several draws) or negative (one draw), Olink
    # plasma. WHO peak severity (WHO COVID-19 clinical management): critical, or fatal disease -> severe_critical;
    # mild or moderate -> mild_moderate; severe (oxygen for severe pneumonia, WHO progression scale 5 or 6) -> not
    # assigned. Days from first symptoms where recorded; acute window: up to 30 days.
    inf = d["Case_Control"].map({"POSITIVE": "acute_covid", "NEGATIVE": "never_infected"})
    sev = np.where(d["Fatal_Disease"] == True, "severe_critical", d["WHO_Severity_Peak"].map({"critical": "severe_critical", "mild": "mild_moderate", "moderate": "mild_moderate"}))  # noqa: E712
    days = d["Time_From_First_Symptoms"]
    order = d["SampleID"].str.extract(r"_(\d+)$")[0].astype(float)
    return _frame(sample_id=d["SampleID"], participant_id="Gisby2021:" + d["Individual_ID"].astype(str), order=order,
                  label_original="Case_Control=" + d["Case_Control"], infection_class=inf,
                  severity_original="WHO_Severity_Peak=" + d["WHO_Severity_Peak"].astype(str) + "; Fatal_Disease=" + d["Fatal_Disease"].astype(str),
                  severity_class=np.where(inf == "acute_covid", sev, None), days=days, time_anchor=np.where(inf == "acute_covid", "symptom_onset", None),
                  in_acute_window=np.where(inf == "acute_covid", np.where(days.isna(), np.nan, days <= 30), np.nan),
                  sex=d["Sex"].map(_norm_sex), age=np.nan, batch_fields="plate=" + d["Plate_ID"], plate=d["Plate_ID"], age_band=d["Age"])


@labeler("COMBAT2022_WB")
def _combat(d: pd.DataFrame) -> pd.DataFrame:
    # COMBAT (Oxford): whole blood RNA-seq samples. COVID_MILD, COVID_SEV, COVID_CRIT, COVID_HCW_MILD,
    # COVID_HCW_SEVERE, COVID_HCW_ASYMPTOMATIC and COVID_LDN are acute COVID-19; HV healthy volunteers are never
    # infected; Sepsis, Flu and Sepsis_CONV are other disease; COVID_CONV and COVID_HCW_CONV (convalescent) are
    # excluded. Severity: COVID_CRIT (intensive care or ventilation) or death within 28 days -> severe_critical;
    # COVID_MILD and COVID_HCW_MILD -> mild_moderate; COVID_SEV (oxygen, WHO 5 or 6) and the rest -> not assigned.
    d = d[d["RNASeq_sample_ID"].notna()].copy()
    src = d["Source"]
    acute = src.isin(["COVID_MILD", "COVID_SEV", "COVID_CRIT", "COVID_HCW_MILD", "COVID_HCW_SEVERE", "COVID_HCW_ASYMPTOMATIC", "COVID_LDN"])
    inf = np.where(acute, "acute_covid", np.where(src == "HV", "never_infected", np.where(src.isin(["Sepsis", "Flu", "Sepsis_CONV"]), "other_disease", "excluded")))
    sev = np.where(d["Death28"] == 1, "severe_critical", src.map({"COVID_CRIT": "severe_critical", "COVID_MILD": "mild_moderate", "COVID_HCW_MILD": "mild_moderate"}))
    order = d["COMBAT_participant_timepoint_ID"].str.extract(r"(\d+)$")[0].astype(float)
    return _frame(sample_id=d["RNASeq_sample_ID"], participant_id="COMBAT:" + d["COMBAT_ID"].astype(str), order=order, label_original="Source=" + src,
                  infection_class=inf, severity_original="Source=" + src + "; Death28=" + d["Death28"].astype(str), severity_class=np.where(inf == "acute_covid", sev, None),
                  days=np.nan, time_anchor=np.where(inf == "acute_covid", "admission", None), in_acute_window=np.where(inf == "acute_covid", True, np.nan),
                  sex=d["Sex"].map({0: "male", 1: "female"}), age=np.nan, batch_fields="institute=" + d["Institute"], institute=d["Institute"], age_band=d["Age"])


@labeler("GSE300129")
def _gse300129(d: pd.DataFrame) -> pd.DataFrame:
    # Mechanically ventilated COVID-19 patients in intensive care (thrombotic or not) and healthy controls
    inf = d["ch_status"].map({"Healthy_Control": "never_infected", "Thrombotic": "acute_covid", "Non_Thrombotic": "acute_covid"})
    return _frame(sample_id=d["geo_accession"], participant_id="GSE300129:" + d["geo_accession"], order=0, label_original="status=" + d["ch_status"],
                  infection_class=inf, severity_original=np.where(inf == "acute_covid", "mechanically ventilated (series summary)", None),
                  severity_class=np.where(inf == "acute_covid", "severe_critical", None), days=np.nan, time_anchor=np.where(inf == "acute_covid", "admission", None),
                  in_acute_window=np.where(inf == "acute_covid", True, np.nan), sex=d["ch_gender"].map(_norm_sex), age=d["ch_age"].map(_num), batch_fields="")


@labeler("GSE243219")
def _gse243219(d: pd.DataFrame) -> pd.DataFrame:
    # PBMC of COVID-19 patients, sepsis patients (other disease) and healthy donors; the series also holds miRNA-seq
    # of the same people, so only the samples titled "(RNA-seq)" are used (platform GPL30209; the miRNA-seq samples are
    # on GPL24676)
    d = d[d["title"].str.contains(r"\(RNA-seq\)", regex=True)]
    inf = d["ch_disease_state"].map({"COVID-19": "acute_covid", "healthy donor": "never_infected", "Sepsis": "other_disease"})
    return _frame(sample_id=d["geo_accession"], participant_id="GSE243219:" + d["title"].str.split(" ").str[0], order=0, label_original="disease_state=" + d["ch_disease_state"],
                  infection_class=inf, severity_original=np.nan, severity_class=np.nan, days=np.nan, time_anchor=np.where(inf == "acute_covid", "admission", None),
                  in_acute_window=np.where(inf == "acute_covid", True, np.nan), sex="unknown", age=np.nan, batch_fields="platform=" + d["platform_id"], platform_id=d["platform_id"])


@labeler("GSE213313")
def _gse213313(d: pd.DataFrame) -> pd.DataFrame:
    # Serial whole blood of 19 critical and 15 severe but non-critical (not in intensive care) COVID-19 patients and
    # 11 healthy controls; T1 is the first sample. Critical -> severe_critical; Non-critical -> mild_moderate.
    inf = d["ch_disease_state"].map({"COVID-19": "acute_covid", "Healthy": "never_infected"})
    sev = d["ch_severity"].map({"Critical": "severe_critical", "Non-critical": "mild_moderate"})
    return _frame(sample_id=d["geo_accession"], participant_id="GSE213313:" + d["ch_individual"], order=d["ch_time"].str.extract(r"(\d+)")[0].astype(float),
                  label_original="disease_state=" + d["ch_disease_state"] + "; time=" + d["ch_time"], infection_class=inf, severity_original="severity=" + d["ch_severity"],
                  severity_class=np.where(inf == "acute_covid", sev, None), days=np.nan, time_anchor=np.where(inf == "acute_covid", "admission", None),
                  in_acute_window=np.where(inf == "acute_covid", True, np.nan), sex="unknown", age=np.nan, batch_fields="")


@labeler("GSE207015")
def _gse207015(d: pd.DataFrame) -> pd.DataFrame:
    # SomaScan 7k plasma of pregnant and non-pregnant adults with SARS-CoV-2 infection (Cases) and uninfected controls.
    # NIH clinical spectrum: Critical -> severe_critical; Asymptomatic, Mild, Moderate -> mild_moderate; Severe (oxygen
    # saturation below 94 percent, without intensive care) -> not assigned.
    inf = d["ch_group"].map({"Cases": "acute_covid", "Controls": "never_infected"})
    nih = d["ch_clinical_spectrum_nih_classification"]
    sev = nih.map({"Critical": "severe_critical", "Asymptomatic": "mild_moderate", "Mild": "mild_moderate", "Moderate": "mild_moderate"})
    pregnant = np.where(d["ch_ga"].astype(str) != "NA", "pregnant", "not pregnant")
    return _frame(sample_id=d["geo_accession"], participant_id="GSE207015:" + d["ch_individual"], order=0, label_original="group=" + d["ch_group"],
                  infection_class=inf, severity_original="NIH=" + nih, severity_class=np.where(inf == "acute_covid", sev, None), days=np.nan,
                  in_acute_window=np.where(inf == "acute_covid", np.nan, np.nan), sex=d["ch_sex"].map(_norm_sex), age=d["ch_age"].map(_num),
                  batch_fields="pregnancy=" + pd.Series(pregnant, index=d.index), pregnancy=pregnant)


# ------------------------------------------------------------------------------------------------ cohorts and decisions
# COHORTS: one entry per cohort that may be built. layer RNA or proteins; tissue in the schema vocabulary; platform;
# meta: "geo" or a function; matrix_kind: counts, normalised or raw; files: what fetch downloads; paper: PMID and DOI.
COHORTS: dict = {
    "GSE198449": dict(accession="GSE198449", layer="RNA", tissue="whole_blood", platform="Illumina NovaSeq 6000, PAXgene whole blood RNA-seq", matrix_kind="counts",
                      files=[("geo_suppl", "GSE198449_featureCounts.txt.gz")], paper="PMID 36323307; DOI 10.1016/j.cels.2022.10.005 (CHARM, Cell Systems 2022)"),
    "COMBAT2022_WB": dict(accession="zenodo.6120249", layer="RNA", tissue="whole_blood", platform="whole blood RNA-seq (COMBAT)", matrix_kind="counts", meta=meta_combat,
                          files=[("zenodo", "6120249", "CBD-KEY-RNASEQ-WB.tar.gz"), ("zenodo", "6120249", "CBD-KEY-CLINVAR.tar.gz")],
                          paper="PMID 35216673; DOI 10.1016/j.cell.2022.01.012 (COMBAT consortium)"),
    "GSE161731": dict(accession="GSE161731", layer="RNA", tissue="whole_blood", platform="Illumina NovaSeq 6000 whole blood RNA-seq", matrix_kind="counts",
                      files=[("geo_suppl", "GSE161731_counts.csv.gz"), ("geo_suppl", "GSE161731_counts_key.csv.gz"), ("geo_suppl", "GSE161731_key.csv.gz")],
                      paper="PMID 33597532; DOI 10.1038/s41467-021-21289-y (McClain et al.)"),
    "GSE152641": dict(accession="GSE152641", layer="RNA", tissue="whole_blood", platform="Illumina NovaSeq 6000 whole blood RNA-seq", matrix_kind="counts",
                      files=[("geo_suppl", "GSE152641_Inflammatix_COVID19_counts_entrez.csv.gz")], paper="PMID 33437935; DOI 10.1016/j.isci.2020.101947 (Thair et al.)"),
    "GSE300129": dict(accession="GSE300129", layer="RNA", tissue="whole_blood", platform="Affymetrix array (series matrix values)", matrix_kind="normalised",
                      files=[("geo_matrix", "GSE300129_series_matrix.txt.gz"), ("geo_gpl", "GPL34133")], paper="PMID 41695759"),
    "GSE243219": dict(accession="GSE243219", layer="RNA", tissue="PBMC", platform="Illumina RNA-seq of PBMC", matrix_kind="counts",
                      files=[("geo_suppl", "GSE243219_RAW.tar")], paper="PMID 39160575"),
    "GSE172114": dict(accession="GSE172114", layer="RNA", tissue="whole_blood", platform="Illumina NovaSeq 6000 whole blood RNA-seq (RSEM, TMM)", matrix_kind="normalised",
                      files=[("geo_suppl", "GSE172114_rsem_gene_count_matrix_TMM_69samples.csv.gz")], paper="PMID 34698500; DOI 10.1126/scitranslmed.abj7521 (Carapito et al.)"),
    "GSE171110": dict(accession="GSE171110", layer="RNA", tissue="whole_blood", platform="Illumina HiSeq 2500 whole blood RNA-seq", matrix_kind="counts",
                      files=[("geo_suppl", "GSE171110_Data_RNAseq_raw_counts_geo.txt.gz")], paper="PMID 34127958; DOI 10.1172/jci.insight.148895 (Levy et al.)"),
    "GSE213313": dict(accession="GSE213313", layer="RNA", tissue="whole_blood", platform="Affymetrix array (series matrix values)", matrix_kind="normalised",
                      files=[("geo_matrix", "GSE213313_series_matrix.txt.gz"), ("geo_gpl", "GPL21185")], paper="PMID 37365222"),
    "GSE152418": dict(accession="GSE152418", layer="RNA", tissue="PBMC", platform="Illumina NovaSeq 6000 PBMC RNA-seq", matrix_kind="counts",
                      files=[("geo_suppl", "GSE152418_p20047_Study1_RawCounts.txt.gz")], paper="PMID 32788292; DOI 10.1126/science.abc6261 (Arunachalam et al.)"),
    "Filbin2021_MGH": dict(accession="10.17632/nf853r8xsj", layer="proteins", tissue="plasma", platform="Olink Explore 1536 (NPX)", matrix_kind="normalised", meta=meta_filbin,
                           files=[("mendeley", "nf853r8xsj", "Olink Proteomics.xlsx"), ("mendeley", "nf853r8xsj", "Clinical Metadata.xlsx"), ("mendeley", "nf853r8xsj", "Supplemental-Table-2-Olink-Assays-NPX-values_v2.xlsx")],
                           paper="DOI 10.1016/j.xcrm.2021.100287 (Filbin et al. 2021, Cell Reports Medicine)",
                           unlabelled_expected="later draws (day 3, day 7, event) are not labelled: one draw per participant, day 0"),
    "Gisby2021_plasma": dict(accession="zenodo.4287557", layer="proteins", tissue="plasma", platform="Olink Target 96 (NPX)", matrix_kind="normalised", meta=meta_gisby_plasma,
                             files=[("zenodo", "4287557", "plasma_npx_level.csv"), ("zenodo", "4287557", "plasma_sample_level.csv"), ("zenodo", "4287557", "README.pdf")],
                             paper="PMID 33704068; DOI 10.7554/eLife.64827 (Gisby et al.)"),
    "GSE207015": dict(accession="GSE207015", layer="proteins", tissue="plasma", platform="SomaScan 7k v4.1 (series matrix values)", matrix_kind="normalised",
                      files=[("geo_matrix", "GSE207015_series_matrix.txt.gz"), ("geo_gpl", "GPL32354")], paper="see GEO record"),
}

_RNA_REASON = "included: open, per sample labels in the repository, adults, blood, acute; deposited matrix"
CURATED: dict = {
    **{acc: dict(decision="include", reason=f"{_RNA_REASON}; {why}", cohorts=[cid]) for acc, cid, why in [
        ("GSE198449", "GSE198449", "infection (first PCR positive sample against participants never positive)"),
        ("zenodo.6120249", "COMBAT2022_WB", "infection and severity (whole blood RNA-seq archive of the COMBAT record)"),
        ("GSE161731", "GSE161731", "infection"),
        ("GSE152641", "GSE152641", "infection"),
        ("GSE300129", "GSE300129", "infection (ventilated patients against healthy controls)"),
        ("GSE243219", "GSE243219", "infection (PBMC)"),
        ("GSE172114", "GSE172114", "severity (critical against non-critical)"),
        ("GSE171110", "GSE171110", "infection (10 healthy donors)"),
        ("GSE213313", "GSE213313", "infection and severity"),
        ("GSE152418", "GSE152418", "infection (PBMC)"),
        ("10.17632/nf853r8xsj", "Filbin2021_MGH", "severity (COVID-negative patients are other disease, so no infection contrast)"),
        ("zenodo.4287557", "Gisby2021_plasma", "infection and severity; all participants on haemodialysis"),
        ("GSE207015", "GSE207015", "infection and severity; pregnant and non-pregnant adults"),
    ]},
    "GSE157103": dict(decision="needs decision", reason="leukocytes isolated from whole blood (not whole blood or PBMC as the rule states); severity 50 against 50 (ICU); non-COVID patients are other disease"),
    "GSE217948": dict(decision="needs decision", reason="PREDICT-19: 'control' (71) is not defined in the record; per sample disease state missing for 193 samples"),
    "GSE167000": dict(decision="needs decision", reason="who the 30 SARS-CoV-2 negative subjects are (healthy or symptomatic patients) is not stated in the record"),
    "GSE227080": dict(decision="needs decision", reason="NanoString immunology panel (579 genes) on PBMC; the 59 COVID-negative subjects are not described; severity class MOD_SEV merges moderate and severe"),
    "GSE212866": dict(decision="needs decision", reason="SuperSeries of arrays; the negatives are hospitalised patients (other disease); severity only as ARDS against no ARDS"),
    "GSE183071": dict(decision="needs decision", reason="NanoString panel; blood subset of 54 samples, class counts per participant not read"),
    "GSE215865": dict(decision="exclude", reason="already in CHORUS (acute samples on disk); the bars use new cohorts only"),
    "10.17632/tzydswhhb5": dict(decision="exclude", reason="duplicate: Su et al. 2020 INCOV cohort, the participants of Su2022 already in CHORUS", duplicate_of="Su2022"),
    "10.17632/73j3nzjxd3": dict(decision="exclude", reason="COVIDome: COVID-negative participants are symptomatic patients (other disease) and no severity field: no contrast"),
    "10.17632/xss9pdtc4f": dict(decision="exclude", reason="no per sample protein values deposited (summary tables only)"),
    "10.17632/734559nz2k": dict(decision="exclude", reason="sample list only; no data deposited"),
    "10.17632/44p9ck395x": dict(decision="exclude", reason="fewer than 10 per class: 8 healthy; severity 5 fatal against 10 mild (severe without intensive care not assigned)"),
    "zenodo.3895886": dict(decision="exclude", reason="analysis scripts only; no data deposited"),
    "zenodo.8252959": dict(decision="exclude", reason="48-plex cytokine panel, not a proteomics platform; no healthy controls; 'severe' without intensive care not assigned"),
    "PXD037486": dict(decision="exclude", reason="7 healthy controls (fewer than 10 per class)"),
    "PXD029376": dict(decision="exclude", reason="glycopeptide identifications only; no per sample protein matrix"),
    "GSE225349": dict(decision="exclude", reason="COVID-negative participants are hospitalised patients (other disease) and no severity field: no contrast"),
    "GSE178399": dict(decision="exclude", reason="COVID-19 patients only; outcome (died or survived) has no mild or moderate class"),
    "GSE283501": dict(decision="exclude", reason="COVID-19 ARDS patients only: no contrast"),
    "GSE163151": dict(decision="exclude", reason="whole blood subset: 7 COVID-19 samples (fewer than 10 per class)"),
    "GSE155454": dict(decision="exclude", reason="6 healthy samples of one shared id (fewer than 10 per class); no severity scale"),
    "GSE221234": dict(decision="exclude", reason="COVID-19 only; severity at the BQC19 biobank (application needed)"),
    "GSE264443": dict(decision="exclude", reason="COVID-19 only; phase at sampling, not a per participant severity"),
    "GSE211378": dict(decision="exclude", reason="convalescent, not acute"),
    "10.17632/2mc6rrc5j3": dict(decision="exclude", reason="duplicate: COVIDome datasets, version 2 of 10.17632/73j3nzjxd3 (no contrast)", duplicate_of="10.17632/73j3nzjxd3"),
    "PXD022817": dict(decision="exclude", reason="COVIDome (Colorado) MS plasma: negatives are symptomatic patients (other disease) and no severity field"),
    "PXD046103": dict(decision="exclude", reason="22 patients across severity groups and few negative controls (fewer than 10 per class)"),
    "PXD023521": dict(decision="exclude", reason="SDRF only; no abundance matrix deposited"),
    "PXD029437": dict(decision="exclude", reason="acute COVID-19 patients only (no controls, no severity classes); peptide-level MRM table"),
    "PXD055973": dict(decision="exclude", reason="hospitalised patients only; acquisition setup file, no abundance matrix"),
    "E-MTAB-12777": dict(decision="exclude", reason="antibody array; severity class 'critical or severe' merges severe (oxygen) with critical, so it cannot be mapped"),
    "zenodo.15883238": dict(decision="exclude", reason="IVIG-treated COVID-19 patients only: no contrast"),
    "10.17632/2cbxgsn7vx": dict(decision="needs decision", reason="Osaka Olink NPX (49 acute and 34 later COVID-19 samples, 20 healthy): no participant ids, so one sample per participant cannot be checked; no paper linked"),
}


def included_cohorts() -> list:
    return [c for cur in CURATED.values() if cur["decision"] == "include" for c in cur.get("cohorts", [])]


def cohort_dir(cid: str) -> Path:
    layer = COHORTS[cid]["layer"]
    return ROOT / "data" / {"RNA": "counts", "proteins": "proteomics", "metabolites": "metabolomics"}[layer] / cid


def fetch(cohort: str = None) -> None:
    """Deposited files of the included cohorts into data/<layer>/<cohort>/raw/ (resume, three retries, checksums)."""
    import shutil

    free = shutil.disk_usage(ROOT).free / 1e9
    log.info(f"fetch: {free:.0f} GB free (at least 20 GB are kept free)")
    if free < 20:
        raise RuntimeError("fewer than 20 GB free: fetch stopped")
    for cid in [cohort] if cohort else included_cohorts():
        raw = cohort_dir(cid) / "raw"
        for spec in COHORTS[cid]["files"]:
            if shutil.disk_usage(ROOT).free / 1e9 < 20:
                raise RuntimeError("fewer than 20 GB free: fetch stopped")
            kind = spec[0]
            if kind == "geo_suppl":
                download(f"{geo_dir(COHORTS[cid]['accession'])}/suppl/{spec[1]}", raw / spec[1])
            elif kind == "geo_matrix":
                download(f"{geo_dir(COHORTS[cid]['accession'])}/matrix/{spec[1]}", raw / spec[1])
            elif kind == "geo_gpl":  # platform annotation (SOFT)
                gpl = spec[1]
                download(f"https://ftp.ncbi.nlm.nih.gov/geo/platforms/{gpl[:-3]}nnn/{gpl}/soft/{gpl}_family.soft.gz", raw / f"{gpl}_family.soft.gz")
            elif kind in ("zenodo", "mendeley"):
                files = zenodo_files(spec[1]) if kind == "zenodo" else mendeley_files(spec[1])
                f = files[spec[2]]
                download(f["url"], raw / Path(spec[2]).name, md5=f.get("md5"), sha256=f.get("sha256"))
        log.info(f"fetch {cid}: done")


# ------------------------------------------------------------------------------------------------ build: features
HGNC = ROOT / "data" / "ref" / "hgnc_complete_set.txt"


def hgnc_maps() -> dict:
    """Symbol, previous or alias symbol, Entrez and UniProt to Ensembl gene id (HGNC); an alias or previous symbol is
    used only when it points to a single gene and is not an approved symbol."""
    h = pd.read_csv(HGNC, sep="\t", dtype=str, usecols=["symbol", "prev_symbol", "alias_symbol", "entrez_id", "ensembl_gene_id", "uniprot_ids"]).dropna(subset=["ensembl_gene_id"])
    sym = dict(zip(h["symbol"], h["ensembl_gene_id"]))
    extra = {}
    for col in ("prev_symbol", "alias_symbol"):
        for s, e in zip(h[col], h["ensembl_gene_id"]):
            for a in str(s).split("|") if pd.notna(s) else []:
                extra.setdefault(a, set()).add(e)
    for a, es in extra.items():
        if a not in sym and len(es) == 1:
            sym[a] = next(iter(es))
    entrez = dict(zip(h["entrez_id"].dropna(), h.loc[h["entrez_id"].notna(), "ensembl_gene_id"]))
    return {"symbol": sym, "entrez": entrez}


def collapse(X: pd.DataFrame, ids: pd.Series, how: str) -> tuple:
    """Rename features by ``ids`` (NaN dropped) and resolve duplicates: counts are summed, other values keep the feature
    with the highest mean. Returns (matrix, report)."""
    n0 = len(X)
    ids = ids.reindex(X.index)
    keep = ids.notna()
    X, ids = X[keep], ids[keep]
    dup_targets = int(ids.duplicated(keep=False).groupby(ids).any().sum())
    if how == "sum":
        Y = X.groupby(ids.values).sum(min_count=1)
    else:
        m = X.mean(axis=1)
        order = m.sort_values(ascending=False).index
        first = ~ids.loc[order].duplicated()
        Y = X.loc[order[first.values]]
        Y.index = ids.loc[Y.index].values
    rep = {"n_features_in": n0, "n_mapped": int(keep.sum()), "mapping_rate": round(keep.sum() / n0, 4) if n0 else np.nan,
           "n_targets_with_duplicates": dup_targets, "duplicate_rule": "summed" if how == "sum" else "feature with the highest mean kept", "n_features_out": len(Y)}
    return Y.sort_index(), rep


def ensembl_ids(index: pd.Index, kind: str, maps: dict) -> pd.Series:
    s = pd.Series(index.astype(str), index=index)
    if kind == "ensembl":
        return s.str.extract(r"^(ENSG\d+)")[0]
    if kind == "symbol":
        return s.map(maps["symbol"])
    if kind == "entrez":
        return s.str.replace(r"_at$", "", regex=True).map(maps["entrez"])
    raise ValueError(kind)


def read_series_matrix(path: Path) -> pd.DataFrame:
    with gzip.open(path, "rt") as f:
        for ln in f:
            if ln.startswith("!series_matrix_table_begin"):
                break
        X = pd.read_csv(f, sep="\t", index_col=0, comment="!")
    return X.apply(pd.to_numeric, errors="coerce")


def read_gpl_table(path: Path) -> pd.DataFrame:
    """The platform table of a GEO SOFT family file (streamed; the sample sections are not read)."""
    rows, header, on = [], None, False
    with gzip.open(path, "rt", errors="replace") as f:
        for ln in f:
            if ln.startswith("!platform_table_begin"):
                on = True
                continue
            if ln.startswith("!platform_table_end"):
                break
            if on:
                parts = ln.rstrip("\n").split("\t")
                if header is None:
                    header = parts
                else:
                    rows.append(parts + [""] * (len(header) - len(parts)))
    return pd.DataFrame(rows, columns=header).replace("", np.nan)


def _geo_title_map(acc: str, key=lambda t: t) -> dict:
    d = geo_samples(acc)
    return {key(t): g for t, g in zip(d["title"], d["geo_accession"])}


def read_matrix(cid: str) -> tuple:
    """(matrix with features in rows and the sample ids of the labels in columns, feature id kind, value kind, notes)."""
    raw = cohort_dir(cid) / "raw"
    if cid == "GSE161731":
        X = pd.read_csv(raw / "GSE161731_counts.csv.gz", index_col=0)
        X.columns = X.columns.map(_geo_title_map(cid, lambda t: t.replace("whole blood, ", "")))
        return X, "ensembl", "counts", "columns are the rna_id of the GEO titles"
    if cid == "GSE152641":
        X = pd.read_csv(raw / "GSE152641_Inflammatix_COVID19_counts_entrez.csv.gz", index_col=0)
        X.index = X.index.astype(str)
        X.columns = X.columns.map(_geo_title_map(cid))
        return X, "entrez", "counts", "rows are Entrez gene ids"
    if cid == "GSE171110":
        X = pd.read_csv(raw / "GSE171110_Data_RNAseq_raw_counts_geo.txt.gz", sep="\t", index_col=0, skiprows=[1, 2])
        X = X.apply(pd.to_numeric, errors="coerce")
        X.columns = X.columns.map(_geo_title_map(cid, lambda t: t.split(" ")[0]))
        return X, "symbol", "counts (expected counts, not all integers)", "rows are gene symbols; columns are the Sample_sheet ids of the GEO titles"
    if cid == "GSE172114":
        X = pd.read_csv(raw / "GSE172114_rsem_gene_count_matrix_TMM_69samples.csv.gz", index_col=0)
        X.columns = X.columns.map(_geo_title_map(cid, lambda t: t.split("-")[0]))
        return X, "ensembl", "normalised (TMM, log scale as deposited)", "rows are Ensembl id with version and symbol"
    if cid == "GSE152418":
        X = pd.read_csv(raw / "GSE152418_p20047_Study1_RawCounts.txt.gz", sep="\t", index_col=0)
        X.columns = X.columns.map(_geo_title_map(cid))
        return X, "ensembl", "counts", ""
    if cid == "GSE198449":
        X = pd.read_csv(raw / "GSE198449_featureCounts.txt.gz", sep="\t", index_col=0)
        X = X.drop(columns=[c for c in ("Chr", "Start", "End", "Strand", "Length") if c in X.columns])
        X.columns = X.columns.map(_geo_title_map(cid))
        return X, "ensembl", "counts", "featureCounts"
    if cid == "GSE243219":
        import tarfile

        cols = {}
        with tarfile.open(raw / "GSE243219_RAW.tar") as t:
            for m in t.getmembers():
                if m.name.endswith("_raw_count.txt.gz"):
                    s = pd.read_csv(io.BytesIO(gzip.decompress(t.extractfile(m).read())), sep="\t", index_col=0).iloc[:, 0]
                    cols[m.name.split("_")[0]] = s
        X = pd.DataFrame(cols)
        return X, "symbol", "counts (expected counts, not all integers)", "one file per sample (GSM in the file name); rows are gene symbols"
    if cid == "COMBAT2022_WB":
        import tarfile

        with tarfile.open(raw / "CBD-KEY-RNASEQ-WB.tar.gz") as t:
            X = pd.read_csv(t.extractfile("CBD-KEY-RNASEQ-WB/Raw_count_data_143_60683.txt"), sep="\t", index_col=0)
        return X, "ensembl", "counts", "columns are RNASeq_sample_ID"
    if cid == "GSE300129":
        X = read_series_matrix(raw / "GSE300129_series_matrix.txt.gz")
        return X, "entrez", "normalised (log2 array intensities, series matrix)", "probe sets named <Entrez id>_at (GPL34133, ENTREZID column)"
    if cid == "GSE213313":
        X = read_series_matrix(raw / "GSE213313_series_matrix.txt.gz")
        gpl = read_gpl_table(raw / "GPL21185_family.soft.gz").set_index("ID")
        ens = gpl["ENSEMBL_ID"].where(gpl["ENSEMBL_ID"].astype(str).str.startswith("ENSG"))
        sym = gpl["GENE_SYMBOL"].map(hgnc_maps()["symbol"])
        X.attrs["ids"] = ens.fillna(sym)
        return X, "precomputed", "normalised (log2 array intensities, series matrix)", "Agilent probes (GPL21185): ENSEMBL_ID, else GENE_SYMBOL through HGNC"
    if cid == "Filbin2021_MGH":
        X = pd.read_excel(raw / "Olink Proteomics.xlsx", sheet_name="Olink Proteomics").set_index("Public ID").drop(columns=["Day"]).T
        X.columns = X.columns.astype(str)
        a = pd.read_excel(raw / "Supplemental-Table-2-Olink-Assays-NPX-values_v2.xlsx", sheet_name="2A-Olink-Assay", header=1).set_index("OlinkID")
        X.attrs["ids"] = a["UniProt"].reindex(X.index)
        return X, "precomputed", "NPX (log2 scale)", "Olink assay ids to UniProt by Supplemental Table 2A"
    if cid == "Gisby2021_plasma":
        d = pd.read_csv(raw / "plasma_npx_level.csv")
        X = d.pivot_table(index="UniProt", columns="SampleID", values="NPX", aggfunc="mean")
        X.attrs["ids"] = pd.Series(X.index, index=X.index)
        return X, "precomputed", "NPX (log2 scale)", "long table pivoted; an assay measured on two panels is averaged by the pivot"
    if cid == "GSE207015":
        X = read_series_matrix(raw / "GSE207015_series_matrix.txt.gz")
        gpl = read_gpl_table(raw / "GPL32354_family.soft.gz").set_index("ID")
        u = gpl["UniProt"].where(~gpl["UniProt"].astype(str).str.contains(r"[|,; ]"))  # aptamers to one protein only
        X.attrs["ids"] = u
        return X, "precomputed", "SomaScan RFU (log2 as deposited)", "aptamers to UniProt by GPL32354; aptamers to several proteins dropped"
    raise KeyError(cid)


def build(cohort: str = None) -> None:
    """Matrix with mapped features and metadata in the CHORUS schema for each included cohort, and its label mapping
    table (committed, docs/acute_covid_labels/<cohort>.csv)."""
    sys.path.insert(0, str(ROOT / "scripts"))
    from build_cohort_meta import COLUMNS  # noqa: E402
    from chorus.schema import ClinicalSchema, validate_metadata

    maps = hgnc_maps()
    (DOCS / "acute_covid_labels").mkdir(exist_ok=True)
    (WORK / "build").mkdir(exist_ok=True)
    for cid in [cohort] if cohort else included_cohorts():
        spec = COHORTS[cid]
        L = cohort_labels(cid).reset_index(drop=True)
        X, kind, value_kind, notes = read_matrix(cid)
        X = X.loc[:, ~X.columns.isna()]
        n_cols_in = X.shape[1]
        dup_cols = int(pd.Index(X.columns).duplicated().sum())
        X = X.loc[:, ~pd.Index(X.columns).duplicated()]
        ids = X.attrs.get("ids") if kind == "precomputed" else ensembl_ids(X.index, kind, maps)
        how = "sum" if value_kind.startswith("counts") else "max_mean"
        M, rep = collapse(X, ids, how)
        if spec["layer"] == "RNA":
            M.index.name = "gene_id"
            out = cohort_dir(cid) / f"{cid}_matrix.csv.gz"
            omics = "bulk_rna"
        else:
            M.index.name = "uniprot"
            out = cohort_dir(cid) / f"{cid}_matrix.csv"
            omics = "olink" if "Olink" in spec["platform"] else "somascan"
        in_both = [s for s in L["sample_id"] if s in set(M.columns)]
        M = M[in_both]
        M.to_csv(out)
        Lm = L[L["sample_id"].isin(in_both)].copy()
        sev_band = np.where(Lm["severity_class"] == "severe_critical", "icu_or_ventilated", "unknown")
        meta = pd.DataFrame({
            "sample_id": Lm["sample_id"], "subject_id": Lm["participant_id"], "study_id": cid,
            "lc_status": Lm["infection_class"].map({"acute_covid": "unknown", "never_infected": "never_infected", "other_disease": "disease_comparator"}).fillna("unknown"),
            "case_definition": "not_applicable", "sample_matrix": spec["tissue"], "omics_layer": omics,
            "days_since_infection": pd.to_numeric(Lm["days"], errors="coerce").where(lambda v: v >= 0), "time_anchor": Lm["time_anchor"].where(Lm["time_anchor"].isin(["symptom_onset", "positive_test", "admission"]), None),
            "acute_severity": sev_band, "age": Lm["age"], "sex": Lm["sex"].fillna("unknown"), "platform": spec["platform"],
            "control_type": np.where(Lm["infection_class"] == "never_infected", "never_infected_unspecified", "not_applicable"),
            "sample_phase": np.where(Lm["infection_class"] == "acute_covid", "acute", np.where(Lm["infection_class"] == "never_infected", "not_applicable", "unknown")),
            "timepoint": Lm["order"].astype(str), "notes": "label: " + Lm["label_original"].astype(str) + "; severity: " + Lm["severity_original"].astype(str),
        })
        meta = meta.reindex(columns=COLUMNS)
        validate_metadata(meta, ClinicalSchema(), strict=True)
        meta["days_from_anchor_original"] = Lm["days"].values  # negative when drawn before the anchor (schema field left empty)
        for c in ["infection_class", "severity_class", "severity_original", "label_original", "in_acute_window", "order", "batch_fields"] + [c for c in Lm.columns if c not in LABEL_COLUMNS]:
            meta[f"acute_{c}" if c in ("order",) else c] = Lm[c].values
        meta["one_sample_per_participant"] = meta["sample_id"].isin(one_per_participant(Lm)["sample_id"])
        meta.to_csv(ROOT / "data" / "meta" / f"{cid}_meta.csv", index=False)
        mapping_table(L).to_csv(DOCS / "acute_covid_labels" / f"{cid}.csv", index=False)
        info = {"cohort": cid, "matrix": str(out.relative_to(ROOT)), "meta": f"data/meta/{cid}_meta.csv", "feature_kind": kind, "value_kind": value_kind, "notes": notes,
                "n_columns_in_file": n_cols_in, "n_duplicated_columns": dup_cols, "n_labelled_samples": len(L), "n_matrix_samples": int(X.shape[1]),
                "n_in_both": len(in_both), "labelled_not_in_matrix": sorted(set(L["sample_id"]) - set(X.columns))[:50],
                "matrix_not_labelled": sorted(map(str, set(X.columns) - set(L["sample_id"])))[:50], **rep, "matrix_sha256": sha256_of(out)}
        (WORK / "build" / f"{cid}.json").write_text(json.dumps(info, indent=1, default=str))
        log.info(f"build {cid}: {M.shape[0]} features x {M.shape[1]} samples, mapping {rep['mapping_rate']}, {len(L) - len(in_both)} labelled samples not in the matrix")


# ------------------------------------------------------------------------------------------------ verify
# Participants per class in the paper (or the record text) for check b; "source" says where the numbers are from.
PAPER_COUNTS = {
    "GSE198449": dict(source="PMID 36323307 abstract (whole CHARM cohort, not the RNA-seq subset): 2,885 recruits, 1,170 infected", numbers={}),
    "COMBAT2022_WB": dict(source="record README (CBD-KEY-RNASEQ-WB): 143 samples from 123 individuals", numbers={"n_samples": 143, "n_participants": 123}),
    "GSE161731": dict(source="PMID 33597532 abstract: 46 subjects with COVID-19", numbers={"infection_acute_covid": 46}),
    "GSE152641": dict(source="PMID 33437935 abstract: 62 patients with COVID-19 and 24 healthy controls", numbers={"infection_acute_covid": 62, "infection_never_infected": 24}),
    "GSE300129": dict(source="PMID 41695759 abstract: 57 mechanically ventilated patients (36 with and 21 without thrombosis)", numbers={"infection_acute_covid": 57}),
    "GSE243219": dict(source="PMID 39160575 abstract: 35 COVID-19, 22 bacterial sepsis, 15 healthy", numbers={"infection_acute_covid": 35, "infection_never_infected": 15, "other_disease": 22}),
    "GSE172114": dict(source="PMID 34698500 abstract (whole cohort): 47 critical, 25 non-critical, 22 healthy", numbers={"severity_severe_critical": 47, "severity_mild_moderate": 25}),
    "GSE171110": dict(source="PMID 34127958 abstract (whole study): 61 severe hospitalised patients", numbers={"infection_acute_covid": 61}),
    "GSE213313": dict(source="PMID 37365222 abstract: 17 severe, 15 moderate, 11 healthy (GEO: 19 critical, 15 non-critical)",
                      numbers={"severity_severe_critical": 17, "severity_mild_moderate": 15, "infection_never_infected": 11}),
    "GSE152418": dict(source="PMID 32788292 abstract (all assays): 76 patients, 69 healthy", numbers={"infection_acute_covid": 76, "infection_never_infected": 69}),
    "Filbin2021_MGH": dict(source="Filbin 2021 abstract (DOI 10.1016/j.xcrm.2021.100287): 306 COVID-19 patients, 78 symptomatic controls",
                           numbers={"infection_acute_covid": 306, "other_disease": 78}),
    "Gisby2021_plasma": dict(source="PMID 33704068 abstract: 256 samples from 55 patients with COVID-19; 51 non-infected patients",
                             numbers={"infection_acute_covid": 55, "infection_never_infected": 51}),
    "GSE207015": dict(source="no paper linked in GEO; the GEO record lists 124 cases and 70 controls", numbers={}),
}
IFN_GENES = {"IFI27": "ENSG00000165949", "IFI44L": "ENSG00000137959", "IFI6": "ENSG00000126709", "ISG15": "ENSG00000187608", "RSAD2": "ENSG00000134321"}
CXCL10 = "P02778"
BATCH_COLUMNS = ["batch", "plate", "enrollment_batch", "collection_date", "institute", "platform_id", "pregnancy"]


def _row(check, status, value, details=""):
    return {"check": check, "status": status, "value": value, "details": details}


def load_built(cid: str) -> tuple:
    spec = COHORTS[cid]
    path = cohort_dir(cid) / (f"{cid}_matrix.csv.gz" if spec["layer"] == "RNA" else f"{cid}_matrix.csv")
    M = pd.read_csv(path, index_col=0)
    meta = pd.read_csv(ROOT / "data" / "meta" / f"{cid}_meta.csv", dtype={"sample_id": str}, low_memory=False)
    info = json.loads((WORK / "build" / f"{cid}.json").read_text())
    return M, meta, info


def log_values(M: pd.DataFrame, value_kind: str) -> pd.DataFrame:
    if value_kind.startswith("counts"):
        return np.log2(M.div(M.sum(axis=0), axis=1) * 1e6 + 1)
    return M


def hedges(x: np.ndarray, y: np.ndarray) -> float:
    x, y = x[np.isfinite(x)], y[np.isfinite(y)]
    if len(x) < 2 or len(y) < 2:
        return np.nan
    s = np.sqrt(((len(x) - 1) * x.var(ddof=1) + (len(y) - 1) * y.var(ddof=1)) / (len(x) + len(y) - 2))
    if s == 0:
        return np.nan
    j = 1 - 3 / (4 * (len(x) + len(y)) - 9)
    return float(j * (x.mean() - y.mean()) / s)


def cramers_v(a: pd.Series, b: pd.Series) -> float:
    t = pd.crosstab(a, b)
    if t.shape[0] < 2 or t.shape[1] < 2:
        return np.nan
    from scipy.stats import chi2_contingency

    chi2 = chi2_contingency(t, correction=False)[0]
    return float(np.sqrt(chi2 / (t.values.sum() * (min(t.shape) - 1))))


def mad_outliers(v: pd.Series, k: float = 3.0) -> list:
    v = v.dropna()
    med = v.median()
    mad = (v - med).abs().median()
    if not np.isfinite(mad):
        return []
    if mad == 0:  # all but a few samples equal: any sample that differs from the median is outside
        return list(v.index[v != med])
    return list(v.index[(v - med).abs() > k * mad])


def contrast_labels(meta: pd.DataFrame) -> dict:
    """Binary labels (one sample per participant) for each contrast with at least 10 participants per class."""
    K = meta[meta["one_sample_per_participant"]]
    out = {}
    inf = K[K["infection_class"].isin(["acute_covid", "never_infected"])]
    if min((inf["infection_class"] == "acute_covid").sum(), (inf["infection_class"] == "never_infected").sum()) >= 10:
        out["infection"] = pd.Series((inf["infection_class"] == "acute_covid").astype(int).values, index=inf["sample_id"])
    sev = K[(K["infection_class"] == "acute_covid") & K["severity_class"].isin(list(SEVERITY))]
    if min((sev["severity_class"] == "severe_critical").sum(), (sev["severity_class"] == "mild_moderate").sum()) >= 10:
        out["severity"] = pd.Series((sev["severity_class"] == "severe_critical").astype(int).values, index=sev["sample_id"])
    return out


def cv_auroc(X: pd.DataFrame, y: pd.Series, n_features: int = 500, repeats: int = 5) -> float:
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import roc_auc_score
    from sklearn.model_selection import StratifiedKFold
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler
    from sklearn.impute import SimpleImputer

    Z = X.loc[:, y.index].T
    Z = Z.loc[:, Z.notna().mean() >= 0.8]
    top = Z.var().sort_values(ascending=False).index[:n_features]  # unsupervised filter
    Z = Z[top].values
    aucs = []
    for r in range(repeats):
        p = np.zeros(len(y))
        for tr, te in StratifiedKFold(5, shuffle=True, random_state=r).split(Z, y.values):
            m = make_pipeline(SimpleImputer(strategy="median"), StandardScaler(), LogisticRegression(C=0.1, max_iter=2000, class_weight="balanced")).fit(Z[tr], y.values[tr])
            p[te] = m.predict_proba(Z[te])[:, 1]
        aucs.append(roc_auc_score(y.values, p))
    return float(np.mean(aucs))


def split_half(X: pd.DataFrame, y: pd.Series, n: int = 100, seed: int = 0) -> float:
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    import concordance_noise_ceiling as nc  # noqa: E402

    Z = X.loc[:, y.index].T
    Z = Z.loc[:, Z.notna().mean() >= 0.8].to_numpy(float)
    yy = y.values
    rng = np.random.default_rng(seed)
    r = [nc.spearman(nc.g_of(Z, yy, a), nc.g_of(Z, yy, b)) for a, b in (nc.half_split(yy, rng) for _ in range(n))]
    return float(nc.sb(np.nanmean(r))) if np.isfinite(r).any() else np.nan


def verify_cohort(cid: str) -> pd.DataFrame:
    spec = COHORTS[cid]
    M, meta, info = load_built(cid)
    rows = []
    # a integrity
    empty = [c for c in M.columns if M[c].isna().all() or (M[c].fillna(0) == 0).all()]
    dims = f"matrix {M.shape[0]} features x {M.shape[1]} samples; file had {info['n_matrix_samples']} sample columns; labelled samples in the record {info['n_labelled_samples']}; in both {info['n_in_both']}"
    extra_ok = spec.get("unlabelled_expected")
    clean = info["n_in_both"] == info["n_labelled_samples"] and (info["n_matrix_samples"] == info["n_labelled_samples"] or bool(extra_ok))
    st = "fail" if info["n_in_both"] == 0 else ("pass" if (clean and not empty and not info["n_duplicated_columns"]) else "flag")
    if extra_ok:
        dims += f" (unlabelled samples expected: {extra_ok})"
    rows.append(_row("a integrity", st, dims, f"empty columns {len(empty)}; duplicated columns {info['n_duplicated_columns']}; labelled not in matrix {info['labelled_not_in_matrix'][:10]}; matrix not labelled {info['matrix_not_labelled'][:10]}"))
    # b counts against the paper
    c = class_counts(meta.rename(columns={"subject_id": "participant_id", "acute_order": "order"}))
    pc = PAPER_COUNTS[cid]
    diffs = {k: (c.get(k), v) for k, v in pc["numbers"].items() if c.get(k) != v}
    st = "flag" if (diffs or not pc["numbers"]) else "pass"
    mine = ", ".join(f"{k} {v}" for k, v in c.items() if v)
    rows.append(_row("b counts", st, mine, f"{pc['source']}; differences (built, paper): {diffs}" if pc["numbers"] else f"{pc['source']}; no comparable number"))
    # c value type
    vals = M.to_numpy(float)
    fin = vals[np.isfinite(vals)]
    integer = float(np.mean(np.equal(np.mod(fin, 1), 0))) if fin.size else np.nan
    is_counts = info["value_kind"].startswith("counts")
    lib = M.sum(axis=0) if is_counts else None
    detected = (M > 0).sum(axis=0) if is_counts else M.notna().sum(axis=0)
    txt = (f"{info['value_kind']}; integer fraction {integer:.3f}; range {np.nanmin(vals):.3g} to {np.nanmax(vals):.3g}; zeros {np.mean(fin == 0):.3f}; missing {np.mean(~np.isfinite(vals)):.3f}; "
           f"mapping {info['mapping_rate']} ({info['n_mapped']} of {info['n_features_in']}; {info['n_targets_with_duplicates']} targets with duplicates, {info['duplicate_rule']})")
    if is_counts:
        txt += f"; library size median {lib.median():.3g} (range {lib.min():.3g} to {lib.max():.3g})"
    st = "flag" if (is_counts and integer < 0.99 and "expected" not in info["value_kind"]) else "pass"
    rows.append(_row("c value type", st, txt, info["notes"]))
    L = log_values(M, info["value_kind"])
    # d sex
    if spec["layer"] == "RNA":
        from chorus.sex import infer_sex

        try:
            S = infer_sex(L.T.reindex(columns=[g for g in ("ENSG00000229807", "ENSG00000129824") if g in L.index]))
            rep = meta.set_index("sample_id")["sex"].reindex(S.index)
            both = rep.isin(["female", "male"]) & S["inferred_sex"].isin(["female", "male"])
            agree = float((rep[both] == S.loc[both, "inferred_sex"]).mean()) if both.any() else np.nan
            bad = list(S.index[both & (rep != S["inferred_sex"])])
            p = rep[both].value_counts(normalize=True)
            chance = float((p ** 2).sum()) if both.any() else np.nan  # agreement expected if sex labels were assigned at random
            st = "n/a" if not both.any() else ("fail" if agree <= chance + 0.05 else ("flag" if agree < 0.9 else "pass"))
            rows.append(_row("d sex", st, (f"agreement {agree:.3f} on {int(both.sum())} samples (chance {chance:.3f})" + ("; no better than chance: the sample metadata may not match the expression columns" if st == "fail" else "")) if both.any() else "no reported sex to compare",
                             f"inferred {S['inferred_sex'].value_counts().to_dict()} ({S['inferred_sex_method'].iloc[0]}); disagreeing {bad[:30]}"))
            sex_inferred = S["inferred_sex"]
        except Exception as e:  # noqa: BLE001
            rows.append(_row("d sex", "n/a", "not inferred", str(e)[:200]))
            sex_inferred = None
    else:
        rows.append(_row("d sex", "n/a", "protein panel without sex chromosome markers", ""))
        sex_inferred = None
    # e duplicates within the cohort (cross cohort pairs are added by verify())
    feats = L.loc[L.notna().all(axis=1)].var(axis=1).sort_values(ascending=False).index[:2000]
    R = np.corrcoef(L.loc[feats].T.to_numpy(float))
    pid = meta.set_index("sample_id")["subject_id"].reindex(L.columns).to_numpy()
    iu = np.triu_indices(len(R), 1)
    hi = [(L.columns[i], L.columns[j], round(R[i, j], 4)) for i, j in zip(*iu) if R[i, j] > 0.99 and pid[i] != pid[j]]
    same = sum(1 for i, j in zip(*iu) if R[i, j] > 0.99 and pid[i] == pid[j])
    rows.append(_row("e duplicates (within)", "flag" if hi else "pass", f"{len(hi)} pairs of different participants with correlation above 0.99",
                     f"{hi[:20]}; pairs of the same participant above 0.99: {same}; correlation on the {len(feats)} most variable complete features"))
    # f outliers
    mc = pd.Series((R.sum(axis=1) - 1) / (len(R) - 1), index=L.columns)
    out = {"library size": mad_outliers(lib) if is_counts else [], "detected features": mad_outliers(detected), "mean correlation": mad_outliers(mc)}
    n_out = len(set(sum(out.values(), [])))
    rows.append(_row("f outliers", "flag" if n_out else "pass", f"{n_out} samples more than 3 MAD from the cohort median",
                     "; ".join(f"{k}: {len(v)} {v[:10]}" for k, v in out.items())))
    # g label direction, h signal, i confounding, k demographics on one sample per participant
    con = contrast_labels(meta)
    if "infection" in con:
        y = con["infection"]
        if spec["layer"] == "RNA":
            gs = {s: hedges(L.loc[e, y.index[y == 1]].to_numpy(float), L.loc[e, y.index[y == 0]].to_numpy(float)) for s, e in IFN_GENES.items() if e in L.index}
        else:
            gs = {"CXCL10": hedges(L.loc[CXCL10, y.index[y == 1]].to_numpy(float), L.loc[CXCL10, y.index[y == 0]].to_numpy(float))} if CXCL10 in L.index else {}
        neg = sum(1 for v in gs.values() if np.isfinite(v) and v < 0)
        st = "n/a" if not gs else ("flag" if neg > len(gs) / 2 else "pass")
        rows.append(_row("g label direction", st, f"{neg} of {len(gs)} markers lower in acute COVID", "; ".join(f"{k} g={v:.2f}" for k, v in gs.items())))
    else:
        rows.append(_row("g label direction", "n/a", "no infection contrast", ""))
    for k_, y in con.items():
        rows.append(_row(f"h signal ({k_})", "report", f"cross validated AUROC {cv_auroc(L, y):.3f}; split half reliability {split_half(L, y):.3f}",
                         f"{int(y.sum())} against {int((y == 0).sum())} participants; logistic regression on the 500 most variable features, 5 x 5-fold"))
    if not con:
        rows.append(_row("h signal", "n/a", "no contrast with 10 per class in the built data", ""))
    bcols = [b for b in BATCH_COLUMNS if b in meta.columns and meta[b].notna().any()]
    if bcols and con:
        vs, worst = [], 0
        for k_, y in con.items():
            mm = meta.set_index("sample_id").loc[y.index]
            for b in bcols:
                v = cramers_v(y, mm[b].astype(str))
                vs.append(f"{k_} x {b}: V={v:.2f}" if np.isfinite(v) else f"{k_} x {b}: one level")
                worst = max(worst, v if np.isfinite(v) else 0)
        rows.append(_row("i confounding", "flag" if worst >= 0.8 else "pass", f"largest Cramer's V {worst:.2f}", "; ".join(vs)))
    else:
        rows.append(_row("i confounding", "n/a", "no batch field recorded" if not bcols else "no contrast", ""))
    K = meta[meta["one_sample_per_participant"] & (meta["infection_class"] == "acute_covid")]
    d = pd.to_numeric(K["days_from_anchor_original"], errors="coerce")
    win = K["in_acute_window"].map({True: True, False: False, "True": True, "False": False, 1.0: True, 0.0: False})
    outside = int(((d > 30) | (d < 0) | (win == False)).sum())  # noqa: E712
    txt = f"{len(K)} acute participants; days known for {int(d.notna().sum())}"
    if d.notna().any():
        txt += f"; median {d.median():.0f} (IQR {d.quantile(0.25):.0f} to {d.quantile(0.75):.0f}, range {d.min():.0f} to {d.max():.0f})"
    rows.append(_row("j timing", "flag" if outside else "pass", txt, f"{outside} outside 0 to 30 days or marked outside the window; window unknown for {int(win.isna().sum())}; anchor {K['time_anchor'].dropna().unique().tolist()}"))
    Kall = meta[meta["one_sample_per_participant"]]
    dem = []
    for cls, g in Kall.groupby("infection_class"):
        a = pd.to_numeric(g["age"], errors="coerce")
        dem.append(f"{cls}: n {len(g)}, age median {a.median():.0f} (missing {int(a.isna().sum())}), sex {g['sex'].value_counts().to_dict()}" if a.notna().any()
                   else f"{cls}: n {len(g)}, age missing, sex {g['sex'].value_counts().to_dict()}")
    rows.append(_row("k demographics", "report", "; ".join(dem), f"age band columns: {[c for c in ('age_band', 'age_category') if c in meta.columns]}"))
    V = pd.DataFrame(rows)
    V.insert(0, "cohort", cid)
    if sex_inferred is not None:
        meta = meta.set_index("sample_id")
        meta["inferred_sex"] = sex_inferred.reindex(meta.index)
        meta.reset_index().to_csv(ROOT / "data" / "meta" / f"{cid}_meta.csv", index=False)
    return V


def cross_duplicates(cids: list) -> pd.DataFrame:
    """Pairs of samples from different cohorts of the same layer with correlation above 0.99 on the shared features
    (the 2,000 most variable complete shared features)."""
    data = {}
    for cid in cids:
        M, meta, info = load_built(cid)
        data[cid] = log_values(M, info["value_kind"])
    rows = []
    for i, a in enumerate(cids):
        for b in cids[i + 1:]:
            if COHORTS[a]["layer"] != COHORTS[b]["layer"]:
                continue
            shared = data[a].index.intersection(data[b].index)
            A, B = data[a].loc[shared], data[b].loc[shared]
            ok = A.notna().all(axis=1) & B.notna().all(axis=1)
            feats = (A[ok].var(axis=1) + B[ok].var(axis=1)).sort_values(ascending=False).index[:2000]
            if len(feats) < 50:
                rows.append({"cohort_a": a, "cohort_b": b, "n_shared_features": len(feats), "n_pairs_above_0.99": np.nan, "max_correlation": np.nan})
                continue
            Za = (A.loc[feats] - A.loc[feats].mean(axis=0)) / A.loc[feats].std(axis=0)
            Zb = (B.loc[feats] - B.loc[feats].mean(axis=0)) / B.loc[feats].std(axis=0)
            R = (Za.T.to_numpy() @ Zb.to_numpy()) / (len(feats) - 1)
            rows.append({"cohort_a": a, "cohort_b": b, "n_shared_features": len(feats), "n_pairs_above_0.99": int((R > 0.99).sum()), "max_correlation": float(np.nanmax(R))})
    return pd.DataFrame(rows)


def verify(cohort: str = None) -> None:
    out = WORK / "verify"
    out.mkdir(exist_ok=True)
    cids = [cohort] if cohort else included_cohorts()
    for cid in cids:
        V = verify_cohort(cid)
        V.to_csv(out / f"{cid}.csv", index=False)
        log.info(f"verify {cid}: " + "; ".join(f"{r.check} {r.status}" for r in V.itertuples()))
    if not cohort:
        X = cross_duplicates(cids)
        X.to_csv(out / "cross_cohort_duplicates.csv", index=False)
        log.info(f"verify cross cohort: {int(X['n_pairs_above_0.99'].fillna(0).sum())} sample pairs above 0.99 across cohorts")


def cohort_status(V: pd.DataFrame, X: pd.DataFrame, cid: str) -> tuple:
    flags = [f"{r.check} (FAIL): {r.value}" for r in V.itertuples() if r.status == "fail"] + [f"{r.check}: {r.value}" for r in V.itertuples() if r.status == "flag"]
    xs = X[((X["cohort_a"] == cid) | (X["cohort_b"] == cid)) & (X["n_pairs_above_0.99"].fillna(0) > 0)]
    flags += [f"e duplicates (across): {int(r['n_pairs_above_0.99'])} pairs with {r['cohort_b'] if r['cohort_a'] == cid else r['cohort_a']}" for _, r in xs.iterrows()]
    failed = (V["status"] == "fail").any() or not any(V["check"].str.startswith("h signal ("))
    return ("failed" if failed else ("verified with flags" if flags else "verified")), flags


# ------------------------------------------------------------------------------------------------ report
def _md_escape(s) -> str:
    return str(s).replace("|", "/").replace("\n", " ") if pd.notna(s) else ""


def write_configs(status: dict) -> list:
    """One config per layer and contrast (enabled: false), with the cohorts verified or verified with flags."""
    written = []
    for layer in ("RNA", "proteins"):
        for contrast in ("infection", "severity"):
            cids = [c for c in included_cohorts() if COHORTS[c]["layer"] == layer and status[c][0] != "failed" and contrast in status[c][2]]
            if not cids:
                continue
            name = f"acute_covid_{'rna' if layer == 'RNA' else 'proteomics'}_{contrast}"
            lines = [f"# Acute COVID validation (Part B, 18 November 2026): {layer}, {contrast} contrast. Built by scripts/acute_covid_sweep.py; cohorts",
                     "# verified or verified with flags (docs/acute_covid_cohorts.md). Every cohort is enabled: false; nothing has been run on them.",
                     "# One sample per participant: the meta column one_sample_per_participant marks the sample to keep (the first acute",
                     "# sample, otherwise the first sample). The classes are in infection_class and severity_class.",
                     f"layer: {'bulk_rna' if layer == 'RNA' else 'proteomics'}", "cohorts:"]
            for c in cids:
                s = COHORTS[c]
                m = f"data/counts/{c}/{c}_matrix.csv.gz" if layer == "RNA" else f"data/proteomics/{c}/{c}_matrix.csv"
                extra = f"is_counts: {'true' if s['matrix_kind'] == 'counts' else 'false'}, provenance: author_matrix" if layer == "RNA" else "feature_space: uniprot"
                lines.append(f"  - {{study_id: {c}, matrix: {m}, meta: data/meta/{c}_meta.csv, {extra}, accession: \"{s['accession']}\", enabled: false}}")
                lines.append(f"    # {status[c][0]}" + (f"; flags: {'; '.join(f.split(':')[0] for f in status[c][1])}" if status[c][1] else ""))
            if contrast == "infection":
                lines += ["preprocessing:", "  relabel_acute: true            # sample_phase acute -> lc_status acute", "  keep_lc_status: [acute, never_infected]",
                          "  deduplicate: true", "  features: per_cohort", "benchmark:", "  positive: acute", "  negative: never_infected"]
            else:
                lines += ["preprocessing:", "  deduplicate: true", "  features: per_cohort", "benchmark:",
                          "  # the severity label is the meta column severity_class; the chorus loader reads lc_status only, so Part C maps it",
                          "  label_column: severity_class", "  positive: severe_critical", "  negative: mild_moderate"]
            lines.append(f"out: results/{name}")
            p = ROOT / "configs" / f"{name}.yaml"
            if p.exists() and not p.read_text().startswith("# Acute COVID validation (Part B"):
                raise RuntimeError(f"{p} exists and was not written by this script: not overwritten")
            p.write_text("\n".join(lines) + "\n")
            written.append(str(p.relative_to(ROOT)))
    return written


def report() -> None:
    D = pd.read_csv(CANDIDATES, dtype=str)
    X = pd.read_csv(WORK / "verify" / "cross_cohort_duplicates.csv")
    reg = pd.read_csv(WORK / "files.csv", dtype=str)
    status, cards, man = {}, [], []
    for cid in included_cohorts():
        V = pd.read_csv(WORK / "verify" / f"{cid}.csv", dtype=str)
        st, flags = cohort_status(V, X, cid)
        M, meta, info = load_built(cid)
        c = class_counts(meta.rename(columns={"subject_id": "participant_id", "acute_order": "order"}))
        con = contrasts_available(c)
        status[cid] = (st, flags, con)
        spec = COHORTS[cid]
        acc = [a for a, cur in CURATED.items() if cid in cur.get("cohorts", [])][0]
        raw = reg[reg["path"].str.startswith(str(cohort_dir(cid).relative_to(ROOT)))]
        man.append({"cohort": cid, "accession": acc, "layer": spec["layer"], "tissue": spec["tissue"], "platform": spec["platform"], "paper": spec["paper"],
                    "matrix": info["matrix"], "matrix_sha256": info["matrix_sha256"], "meta": info["meta"], "meta_sha256": sha256_of(ROOT / info["meta"]),
                    "label_mapping": f"docs/acute_covid_labels/{cid}.csv", "raw_files": "; ".join(f"{r.path} ({r.size} bytes, sha256 {r.sha256}{', ' + r.repository_checksum + ' verified' if r.checksum_ok == 'True' else ''})" for r in raw.itertuples()),
                    "n_features": M.shape[0], "n_samples": M.shape[1], **{k: v for k, v in c.items()}, "contrasts": ";".join(con), "status": st, "flags": " | ".join(flags)})
        rows = "\n".join(f"| {_md_escape(r.check)} | {r.status} | {_md_escape(r.value)} | {_md_escape(r.details)[:400]} |" for r in V.itertuples())
        mp = pd.read_csv(DOCS / "acute_covid_labels" / f"{cid}.csv", dtype=str)
        cards.append(f"""### {cid}

- **Source:** {acc} ({spec['layer']}, {spec['tissue'].replace('_', ' ')}); files: {', '.join(Path(r.path).name for r in raw.itertuples())}
- **Paper:** {spec['paper']}
- **Platform:** {spec['platform']}; {info['value_kind']}; {M.shape[0]} features x {M.shape[1]} samples
- **Participants (one sample each):** infection {c['infection_acute_covid']} acute COVID against {c['infection_never_infected']} never infected; severity {c['severity_severe_critical']} severe or critical against {c['severity_mild_moderate']} mild or moderate; {c['other_disease']} other disease (own class); {c['excluded']} excluded
- **Contrasts with at least 10 per class:** {', '.join(con) or 'none'}
- **Status:** {st}{(' (' + '; '.join(flags) + ')') if flags else ''}
- **Label mapping:** `docs/acute_covid_labels/{cid}.csv` ({len(mp)} rows)

| check | status | result | details |
|---|---|---|---|
{rows}
""")
    pd.DataFrame(man).to_csv(ROOT / "data" / "meta" / "acute_covid_manifest.csv", index=False)
    cfgs = write_configs(status)
    order = {"include": 0, "needs decision": 1, "exclude": 2}
    T = D.assign(_o=D["decision"].map(order)).sort_values(["_o", "layer", "source", "accession"])
    table = "\n".join(f"| {_md_escape(r.accession)} | {_md_escape(r.source)} | {_md_escape(r.layer)} | {_md_escape(str(r.title)[:90])} | {r.decision} | {_md_escape(r.reason)[:220]} |" for r in T.itertuples())
    n = D["decision"].value_counts().to_dict()
    xs = X[X["n_pairs_above_0.99"].fillna(0) > 0]
    lines = [f"""# Acute COVID cohorts (Part B of the acute COVID validation)

Built 18 November 2026 by `scripts/acute_covid_sweep.py` (search, screen, fetch, build, verify, report). The rule and bars of
the validation are in `results/README.md` ("Acute COVID validation: rule and bars"). No benchmark or noise ceiling analysis
has been run on these cohorts. Data files are not committed; the label mapping tables (`docs/acute_covid_labels/`), this
file, `docs/acute_covid_candidates.csv` and the configs are.

**Candidates:** {len(D)} ({n.get('include', 0)} included, {n.get('exclude', 0)} excluded, {n.get('needs decision', 0)} needing a decision).
Sources: GEO (E-utilities), ArrayExpress and BioStudies, SRA BioProjects, PRIDE and ProteomeXchange, Metabolomics
Workbench, Zenodo, Mendeley Data (DataCite) and the accessions in the full text of {''}open access papers (Europe PMC).

**Severity mapping used for every cohort.** severe or critical: intensive care, mechanical ventilation, WHO progression
scale 6 or more, or death; mild or moderate: a class the study defines below that line (ambulatory, hospitalised without
oxygen, or explicitly not critical / not in intensive care). A class that straddles the line (for example "severe" meaning
supplemental oxygen of any type, or a merged "moderate/severe" class) is not assigned. Each cohort's table in
`docs/acute_covid_labels/` gives the original text and the class.

**Included cohorts.**

| cohort | layer | tissue | infection (acute / never infected) | severity (severe or critical / mild or moderate) | status |
|---|---|---|---|---|---|
""" + "\n".join(f"| {m['cohort']} | {m['layer']} | {m['tissue']} | {m['infection_acute_covid']} / {m['infection_never_infected']} | {m['severity_severe_critical']} / {m['severity_mild_moderate']} | {m['status']} |" for m in man) + f"""

Contrasts reach at least 10 participants per class where listed in the cards. Cross cohort duplicates (correlation above
0.99, same layer): {'none' if xs.empty else '; '.join(f"{r.cohort_a} and {r.cohort_b}: {int(r['n_pairs_above_0.99'])}" for _, r in xs.iterrows())}.

Configs (enabled: false): {', '.join(f'`{c}`' for c in cfgs)}.

## Included cohorts
""", *cards, f"""
## All candidates

{n.get('include', 0)} included, {n.get('needs decision', 0)} needing a decision, {n.get('exclude', 0)} excluded. The same table with every column is
`docs/acute_covid_candidates.csv`.

| accession | source | layer | title | decision | reason |
|---|---|---|---|---|---|
{table}
"""]
    (DOCS / "acute_covid_cohorts.md").write_text("\n".join(lines))
    log.info(f"report: docs/acute_covid_cohorts.md, data/meta/acute_covid_manifest.csv, configs {cfgs}")
    print(f"candidates: {len(D)}; included {n.get('include', 0)}; excluded {n.get('exclude', 0)}; needs decision {n.get('needs decision', 0)}")
    for layer in ("RNA", "proteins"):
        for contrast in ("infection", "severity"):
            cs = [m for m in man if m["layer"] == layer and contrast in m["contrasts"]]
            key = ("infection_acute_covid", "infection_never_infected") if contrast == "infection" else ("severity_severe_critical", "severity_mild_moderate")
            print(f"{layer} {contrast}: " + ", ".join(f"{m['cohort']} {m[key[0]]}/{m[key[1]]} ({m['status']})" for m in cs))
    for m in man:
        if m["flags"]:
            print(f"flags {m['cohort']}: {m['flags']}")
    n_inf = sum(1 for m in man if m["layer"] == "RNA" and "infection" in m["contrasts"] and m["status"] != "failed")
    print(f"RNA cohorts verified for the infection contrast: {n_inf}" + ("  -- FEWER THAN 4" if n_inf < 4 else ""))


def _sanity_curated() -> None:
    for acc, cur in CURATED.items():
        for cid in cur.get("cohorts", []):
            assert cid in COHORTS and cid in LABELERS, cid


def cohort_labels(cid: str) -> pd.DataFrame:
    spec = COHORTS[cid]
    if spec.get("meta", "geo") == "geo":
        return LABELERS[cid](geo_samples(spec["accession"]))
    return LABELERS[cid](spec["meta"]())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("command", choices=["search", "screen", "fetch", "build", "verify", "report"])
    ap.add_argument("--cohort")
    a = ap.parse_args()
    setup_logging()
    log.info(f"=== {a.command} {a.cohort or ''}")
    globals()[a.command]() if a.cohort is None else globals()[a.command](a.cohort)


if __name__ == "__main__":
    main()
