"""Matrices and ClinicalSchema metadata of the transcriptomic cohorts added on 16 November 2026 (docs/count_matrices.md,
"Cohorts added 16 November 2026").

PRJNA1149196 (Sinclair 2024, Nat Microbiol, PMID 39478108, PMC11602718). Whole blood total RNA-seq, two cohorts by site,
written as two studies because site, library kit, sequencer and sampling time differ:
  PRJNA1149196_SA   cohort 1, South Australia (Central Adelaide Health Network), DNBSEQ-T7, Nugen Universal Plus total
                    RNA with ribosomal and globin depletion; PASC, PASC-CVS, Recovered at ~44 and ~68 weeks after a
                    March-April 2020 infection, and uninfected controls
  PRJNA1149196_QLD  cohort 2, Queensland (COVID OZGenetics), NovaSeq 6000, Illumina stranded total RNA with Ribo-Zero
                    plus; PASC-CVS and Recovered, one sample (BioSample "1 year post infection"; paper: more than 520 days)
Counts: salmon against the GENCODE v47 index of GSE224615 (docs/provenance.md). SA: the first 5,000,000 read pairs of
each run (scripts/salmon_quant_ena.py --max-pairs). QLD: all reads (--runs, data/quant/PRJNA1149196_QLD_full), because
its FASTQ files are sorted by alignment position and a prefix holds only chr1 and chr10 reads. Labels: BioSample attributes ``disease`` and
``Timepoint`` (data/meta/raw_PRJNA1149196/biosamples.csv, from data/meta/raw_PRJNA1149196/biosamples.xml).
The subject is the BioSample ``isolate`` with hyphens removed (SA_LC-068 and SA_LC068 are written both ways).
Matrix columns are the ENA sample_alias.

GSE297499 (Orlinick 2026, J Infect Dis, PMID 40920586, PMC12811860). PBMC bulk RNA-seq of women with neuropsychiatric
Long COVID and controls with documented infection and no persistent symptoms (Yale COVID Mind Study); the CSF samples
of the series are not used. Values: the authors' TPM (HISAT2, StringTie; GSE297499_final_processed_tpm_values_pbmc)
keyed by gene symbol, mapped to Ensembl with scripts/map_symbols_to_ensembl.symbol_map (rows of one Ensembl gene
summed on the linear scale). The loader takes them as author_normalised (log2(x + 1)). Labels: GEO characteristic
``treatment`` (NP-Long COVID, Control), joined on the sample title (underscores removed: "2057_V2_PBMC" is the
column "2057V2_PBMC").

  python scripts/build_rna_cohorts_2026_11.py
"""

from __future__ import annotations

import gzip
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from build_cohort_meta import finish  # noqa: E402
from map_symbols_to_ensembl import symbol_map  # noqa: E402

META = ROOT / "data" / "meta"
QUANT = ROOT / "data" / "quant" / "PRJNA1149196"
QUANT_FULL = {"PRJNA1149196_QLD": ROOT / "data" / "quant" / "PRJNA1149196_QLD_full"}  # quantified on all reads
GSE297499 = ROOT / "data" / "counts" / "GSE297499"

SINCLAIR_STATUS = {"PASC": "LC", "PASC-CVS": "LC", "PASC_CVS": "LC", "Recovered": "COVID_recovered", "Control": "never_infected"}
SINCLAIR_DAYS = {"44 weeks post infection": 308, "68 weeks post infection": 476}
SINCLAIR_NOTE = {
    "PASC": "PASC: persistent symptoms more than 12 weeks after infection with ongoing discomfort at sampling",
    "PASC-CVS": "PASC-CVS: chest pain or tightness more than 12 weeks after infection",
    "Recovered": "Recovered: infected, no persistent symptoms more than 12 weeks after infection",
    "Control": "Healthy: recruited before community transmission or seronegative (paper: contemporaneous, negative spike serology)",
}


def series_matrix_samples(path: Path) -> pd.DataFrame:
    """One row per GSM from a GEO series matrix: title, gsm and one ch_<key> column per characteristic."""
    fields: dict = {}
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("!Sample_"):
                k, *v = line.rstrip("\n").split("\t")
                fields.setdefault(k, []).append([x.strip('"') for x in v])
    df = pd.DataFrame({"title": fields["!Sample_title"][0], "gsm": fields["!Sample_geo_accession"][0]})
    for ch in fields.get("!Sample_characteristics_ch1", []):
        key = ch[0].split(":", 1)[0].strip().replace(" ", "_")
        df[f"ch_{key}"] = [x.split(":", 1)[1].strip() if ":" in x else x for x in ch]
    return df


def gse297499() -> pd.DataFrame:
    tpm = pd.read_csv(GSE297499 / "GSE297499_final_processed_tpm_values_pbmc.txt.gz", sep="\t", index_col=0)
    m = symbol_map(pd.read_csv(ROOT / "data" / "ref" / "hgnc_complete_set.txt", sep="\t", dtype=str))
    ens = tpm.index.to_series().map(m)
    Y = tpm[ens.notna()].groupby(ens[ens.notna()]).sum()
    Y.index.name = "gene_id"
    Y.to_csv(GSE297499 / "GSE297499_pbmc_tpm_ensembl.tsv.gz", sep="\t")
    print(f"GSE297499: {ens.notna().sum()} of {len(tpm)} symbols mapped; {len(Y)} Ensembl genes x {Y.shape[1]} PBMC samples")
    geo = series_matrix_samples(GSE297499 / "GSE297499_series_matrix.txt.gz")
    geo = geo[geo["ch_tissue"] == "PBMC"].copy()
    geo["column"] = geo["title"].str.replace("_", "", regex=False).str.replace("PBMC$", "_PBMC", regex=True)
    missing = set(tpm.columns) - set(geo["column"])
    if missing:
        raise ValueError(f"GSE297499 TPM columns without a GSM: {sorted(missing)}")
    geo = geo.set_index("column").loc[list(tpm.columns)]
    subj = geo.index.str.extract(r"^(\d+)")[0].values
    df = pd.DataFrame({
        "sample_id": geo.index,
        "subject_id": "GSE297499:" + pd.Series(subj),
        "lc_status": geo["ch_treatment"].map({"NP-Long COVID": "LC", "Control": "COVID_recovered"}).values,
        "case_definition": np.where(geo["ch_treatment"] == "NP-Long COVID", "WHO_2021", "not_applicable"),
        "control_type": np.where(geo["ch_treatment"] == "NP-Long COVID", "not_applicable", "COVID_recovered"),
        "sample_phase": "post_acute",
        "sample_matrix": "PBMC",
        "sex": "female",
        "platform": "Illumina NovaSeq X, polyA RNA; author TPM (HISAT2, StringTie)",
        "collection_country": "United States",
        "gsm": geo["gsm"].values,
        "timepoint": "single_visit",
        "notes": np.where(geo["ch_treatment"] == "NP-Long COVID",
                          "self-reported neuropsychiatric Long COVID, WHO definition, documented infection (paper Methods)",
                          "documented SARS-CoV-2 infection without persistent symptoms (paper Methods)"),
    })
    if df["lc_status"].isna().any():
        raise ValueError("GSE297499: unmapped treatment value")
    return finish(df, "GSE297499")


def prjna1149196() -> dict:
    bs = pd.read_csv(META / "raw_PRJNA1149196" / "biosamples.csv", dtype=str)
    runs = pd.read_csv(QUANT / "ena_filereport.tsv", sep="\t", dtype=str)
    bs = bs.merge(runs[["run_accession", "sample_alias"]], left_on="sample_name", right_on="sample_alias", how="left")
    if bs["run_accession"].isna().any():
        raise ValueError("PRJNA1149196: BioSample without a run")
    bs["site"] = np.where(bs["biomaterial_provider"] == "Kirsty Short", "QLD", "SA")
    out = {}
    for site in ["SA", "QLD"]:
        b = bs[bs["site"] == site]
        dis = b["disease"].str.replace("_", "-", regex=False)
        df = pd.DataFrame({
            "sample_id": b["sample_alias"].values,
            "subject_id": ("PRJNA1149196:" + b["isolate"].str.replace("-", "", regex=False)).values,
            "lc_status": dis.map(SINCLAIR_STATUS).values,
            "case_definition": np.where(dis.map(SINCLAIR_STATUS) == "LC", "study_specific", "not_applicable"),
            "control_type": dis.map({"Recovered": "COVID_recovered", "Control": "never_infected_contemporary"}).fillna("not_applicable").values,
            "sample_phase": np.where(dis == "Control", "not_applicable", "post_acute"),
            "sample_matrix": "whole_blood",
            "days_since_infection": b["Timepoint"].map(SINCLAIR_DAYS).values,
            "time_anchor": np.where(b["Timepoint"].isin(list(SINCLAIR_DAYS)), "positive_test", None),
            "sex": b["sex"].replace({"missing": "unknown"}).values,
            "platform": ("DNBSEQ-T7, Nugen Universal Plus total RNA (rRNA and globin depleted)" if site == "SA"
                         else "Illumina NovaSeq 6000, Illumina stranded total RNA (Ribo-Zero plus)"),
            "collection_country": "Australia",
            "overlaps_with_study": "GSE169687 (same South Australian participants at 12-24 weeks; ids not linkable)" if site == "SA" else None,
            "timepoint": b["Timepoint"].values,
            "notes": (dis.map(lambda d: SINCLAIR_NOTE[d.replace("PASC-CVS", "PASC-CVS")]) + "; BioSample " + b["biosample"] + ", run " + b["run_accession"]
                      + ("; days_since_infection is the nominal visit (44 or 68 weeks)" if site == "SA" else "; BioSample 1 year post infection, paper more than 520 days")).values,
        })
        if df["lc_status"].isna().any():
            raise ValueError(f"PRJNA1149196_{site}: unmapped disease value")
        sid = f"PRJNA1149196_{site}"
        out[sid] = finish(df, sid)
    return out


def split_counts(meta: dict) -> None:
    """One gene count matrix per site from the merged salmon output (columns = sample_alias)."""
    for sid, df in meta.items():
        merged = QUANT_FULL.get(sid, QUANT) / "salmon.merged.gene_counts.tsv"
        if not merged.exists():
            print(f"{sid}: {merged} not written yet; matrix not written")
            continue
        X = pd.read_csv(merged, sep="\t", index_col=0)
        cols = [c for c in df["sample_id"] if c in X.columns]
        missing = sorted(set(df["sample_id"]) - set(cols))
        X[cols].to_csv(QUANT / f"{sid}_gene_counts.tsv.gz", sep="\t", float_format="%.3f")
        print(f"{sid}: {len(cols)} samples in the count matrix from {merged.parent.name}" + (f"; not quantified: {missing}" if missing else ""))


def main():
    built = {"GSE297499": gse297499(), **prjna1149196()}
    for sid, df in built.items():
        df.to_csv(META / f"{sid}_meta.csv", index=False)
        print(f"{sid}: {len(df)} samples, {df['subject_id'].nunique()} subjects; lc_status {df['lc_status'].value_counts().to_dict()}")
    split_counts({k: v for k, v in built.items() if k.startswith("PRJNA1149196")})


if __name__ == "__main__":
    main()
