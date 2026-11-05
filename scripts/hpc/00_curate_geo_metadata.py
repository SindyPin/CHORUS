"""LLM assisted, human verified curation of GEO sample characteristics into the CHORUS schema.

Step 1 (this script, no LLM): parse a GEO family SOFT file into a tidy table of raw
characteristics per GSM, with the series summary and design text, and emit a curation
worksheet with the CHORUS columns empty.

Step 2 (LLM, optional): `--draft` sends each sample's raw characteristics plus the series
design text to a model and asks for a JSON draft of the schema fields with a supporting
quote per field. Drafts are written to a separate column set (draft_*) and are never used
until a human copies them into the final columns. Provenance is recorded as
curation_provenance = "LLM_assisted <model> <date>; verified by <initials>".

Usage:
  python 00_curate_geo_metadata.py GSE275334_family.soft.gz --out GSE275334_curation.csv
  python 00_curate_geo_metadata.py GSE275334_family.soft.gz --out ... --draft --model claude-opus-5-5
"""

from __future__ import annotations

import argparse
import gzip
import json
import re
from collections import defaultdict
from pathlib import Path

import pandas as pd

SCHEMA_FIELDS = [
    "subject_id", "lc_status", "case_definition", "sample_matrix", "days_since_infection", "time_anchor",
    "acute_severity", "hospitalised", "age", "sex", "bmi", "vaccinated_before_infection", "variant_era", "notes",
]


def parse_soft(path: Path):
    op = gzip.open if str(path).endswith(".gz") else open
    series, samples, cur = {}, {}, None
    with op(path, "rt", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("^SAMPLE"):
                cur = line.split("=", 1)[1].strip()
                samples[cur] = defaultdict(list)
            elif line.startswith("!Series_"):
                k, _, v = line[1:].partition(" = ")
                series.setdefault(k, []).append(v)
            elif cur and line.startswith("!Sample_"):
                k, _, v = line[1:].partition(" = ")
                samples[cur][k].append(v)
    return series, samples


def tidy(series, samples) -> pd.DataFrame:
    rows = []
    for gsm, d in samples.items():
        chars = {}
        for c in d.get("Sample_characteristics_ch1", []):
            k, _, v = c.partition(":")
            chars[k.strip().lower().replace(" ", "_")] = v.strip()
        rows.append(
            {
                "sample_id": gsm,
                "title": " ".join(d.get("Sample_title", [])),
                "source": " ".join(d.get("Sample_source_name_ch1", [])),
                "raw_characteristics": json.dumps(chars),
                "sra": ";".join(x for x in d.get("Sample_relation", []) if "SRA" in x),
                **{f"raw__{k}": v for k, v in chars.items()},
            }
        )
    df = pd.DataFrame(rows)
    for f in SCHEMA_FIELDS:
        df[f] = ""
    df["study_id"] = " ".join(series.get("Series_geo_accession", [""]))
    df["omics_layer"] = "bulk_rna"
    df["curation_provenance"] = ""
    return df


def draft_with_llm(df: pd.DataFrame, series, model: str) -> pd.DataFrame:
    """Draft schema fields with an LLM; requires the anthropic package and an API key."""
    import datetime
    import anthropic

    client = anthropic.Anthropic()
    design = "\n".join(series.get("Series_summary", []) + series.get("Series_overall_design", []))
    prompt_head = (
        "You curate Long COVID omics sample metadata into a fixed schema. Use only the text given. "
        "For each field return a value from the allowed vocabulary or null, plus a short supporting quote. "
        "Vocabularies: lc_status {LC, COVID_recovered, never_infected, disease_comparator, unknown}; "
        "case_definition {WHO_2021, NICE_NG188, CDC_4wk, NASEM_2024, RECOVER_PASC_score, self_report, study_specific}; "
        "sample_matrix {plasma, serum, whole_blood, PBMC, saliva, urine, BAL, CSF, other}; "
        "time_anchor {symptom_onset, positive_test, admission, discharge, unknown}; "
        "acute_severity {ambulatory, hospitalised_no_oxygen, hospitalised_oxygen, icu_or_ventilated, unknown}; "
        "sex {female, male, other, unknown}; days_since_infection integer days; age in years. "
        "Return JSON with keys: " + ", ".join(SCHEMA_FIELDS) + ", and quotes (object field -> quote).\n\nSERIES DESIGN:\n" + design[:6000]
    )
    drafts = []
    for _, r in df.iterrows():
        msg = client.messages.create(model=model, max_tokens=800, messages=[{"role": "user", "content": prompt_head + f"\n\nSAMPLE {r['sample_id']}: title={r['title']}; source={r['source']}; characteristics={r['raw_characteristics']}"}])
        txt = msg.content[0].text
        try:
            js = json.loads(re.search(r"\{.*\}", txt, re.S).group(0))
        except Exception:
            js = {}
        drafts.append({f"draft_{k}": js.get(k) for k in SCHEMA_FIELDS} | {"draft_quotes": json.dumps(js.get("quotes", {}))})
    out = pd.concat([df.reset_index(drop=True), pd.DataFrame(drafts)], axis=1)
    out["curation_provenance"] = f"LLM_assisted {model} {datetime.date.today()}; NOT YET VERIFIED"
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("soft")
    ap.add_argument("--out", required=True)
    ap.add_argument("--draft", action="store_true")
    ap.add_argument("--model", default="claude-opus-5-5")
    a = ap.parse_args()
    series, samples = parse_soft(Path(a.soft))
    df = tidy(series, samples)
    if a.draft:
        df = draft_with_llm(df, series, a.model)
    df.to_csv(a.out, index=False)
    print(f"{len(df)} samples -> {a.out}")
    print("series summary:", " ".join(series.get("Series_summary", []))[:500])


if __name__ == "__main__":
    main()
