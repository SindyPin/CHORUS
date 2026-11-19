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

GSE228320 (Molinero 2022, Biomed Pharmacother 155:113617, PMID 36058144, PMC9424524; added 17 November 2026). Whole
blood RNA-seq (Tempus tubes) of survivors of COVID-19 ARDS (CIBERESUCICOVID, Spain) at the pulmonary evaluation 3 months
after hospital discharge. Values: the deposited raw counts (GSE228320_Raw_counts.txt.gz, Ensembl, columns are the GSM
descriptions MC1_n). Label: the GEO characteristic ``dlco (%_of_predicted)``: below 60 is a case (the title says
sequela), 60 or above a control (the title says control); the build stops if a title disagrees with its DLCO value. The
case is a lung function outcome, not a Long COVID definition: lc_status LC with case_definition pulmonary_sequelae,
controls COVID_recovered (ARDS survivors with mild or no diffusion impairment, other symptoms not recorded). A
sensitivity cohort: it is in no main benchmark config.

CORSIP2026 (Espin 2026, Front Cell Infect Microbiol, PMID 42328163, PMC13275656; added 17 November 2026). Whole blood
(PAXgene, final study visit, January 2023) of British Columbia paramedics of the CORSIP cohort, profiled with the
NanoString nCounter PanCancer Immune Profiling Panel (a targeted panel, not RNA-seq). Source: the GitHub repository
estefespin/CORSIP_omics (CC0) at the commit in docs/provenance.md. Values: RawNanoStringData.csv, raw nCounter counts;
the 14 positive and negative control probes are dropped (they are used for the QC flags), the endogenous and
housekeeping genes are kept and mapped from symbol to Ensembl as for GSE297499; samples: the 37 of SubjectsNoRep.csv
(the 4 other raw samples are replicates without a label). Labels: SubjectsNoRep.csv LCStatus (Case = LC, WHO 2021 from
the 24 month questionnaire; Control = infected without LC symptoms, COVID_recovered), checked against CORSIP_Full.csv.
A control without any record of infection in CORSIP_Full.csv (positive test, reactive nucleocapsid serology or symptom
date) gets control_type convalescent_unverified. QC (flags, not exclusions; NanoString nSolver defaults): positive
control scaling factor outside 0.3 to 3, housekeeping scaling factor outside 0.1 to 10, positive control linearity R2
below 0.95, POS_E not above the mean plus 2 SD of the negative controls.

  python scripts/build_rna_cohorts_2026_11.py            (PRJNA1149196 and GSE297499)
  python scripts/build_rna_cohorts_2026_11.py GSE228320  (one study)
PRJNA1184005 (Missailidis 2024, Front Immunol 15:1450853, PMID 39691709, PMC11649547; added 18 November 2026). PBMC
RNA-seq (NovaSeq 6000; 13 paired-end runs and one single-end run, vacSG112), La Trobe University, Australia. Counts:
salmon 1.10.3 on all reads with the GENCODE v47 index of GSE224615 (scripts/salmon_quant_ena.py, data/quant/PRJNA1184005).
Matrix columns are the ENA sample_alias, which equals the BioSample isolate. Labels: from the author (S. Annesley,
personal communication, 18 November 2026; PRJNA1184005_LABELS below); the deposit and the paper have no per sample label.
LC: WHO description (paper Methods), case_definition WHO_2021. COVID_recovered ("Post COVID"): age and sex matched
participants who recovered (paper Methods); the author's control definition is pending confirmation and is recorded in
control_definition_status. Sex and age: BioSample attributes (data/meta/PRJNA1184005_biosamples.csv). Per sample time
since infection is not deposited (paper Table 1: mean 6.1 months LC, 5.5 months recovered).

  python scripts/build_rna_cohorts_2026_11.py CORSIP2026 (one study)
  python scripts/build_rna_cohorts_2026_11.py PRJNA1184005 (one study)
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
GSE228320 = ROOT / "data" / "counts" / "GSE228320"
CORSIP = ROOT / "data" / "counts" / "CORSIP2026"
CORSIP_REPO = CORSIP / "raw" / "CORSIP_omics"
PRJNA1184005_BIOSAMPLES = META / "PRJNA1184005_biosamples.csv"
PRJNA1184005_COUNTS = ROOT / "data" / "quant" / "PRJNA1184005" / "salmon.merged.gene_counts.tsv"
PRJNA1184005_QC = ROOT / "data" / "quant" / "PRJNA1184005" / "salmon_qc.tsv"
# S. Annesley, personal communication, 18 November 2026 (LC; "Post COVID" = COVID_recovered)
PRJNA1184005_LABELS = {
    **{s: "LC" for s in ["SG22", "vacSG60", "vacSG116", "vacSG117", "vacSG47", "vacSG55", "vacSG61"]},
    **{s: "COVID_recovered" for s in ["vacSG102", "vacSG105", "vacSG111", "vacSG21", "vacSG34", "vacSG70", "vacSG112"]},
}
NANOSTRING_POS_FM = {"POS_A": 128, "POS_B": 32, "POS_C": 8, "POS_D": 2, "POS_E": 0.5, "POS_F": 0.125}

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


def gse228320() -> pd.DataFrame:
    geo = series_matrix_samples(GSE228320 / "GSE228320_series_matrix.txt.gz")
    with gzip.open(GSE228320 / "GSE228320_series_matrix.txt.gz", "rt") as fh:
        desc = [line.rstrip("\n").split("\t")[1:] for line in fh if line.startswith("!Sample_description")]
    geo["column"] = [x.strip('"') for x in desc[1]]  # second description line: MC1_n, the count matrix column
    X = pd.read_csv(GSE228320 / "GSE228320_Raw_counts.txt.gz", sep="\t", index_col=0, nrows=1)
    if set(X.columns) != set(geo["column"]):
        raise ValueError("GSE228320: count matrix columns and GSM descriptions differ")
    dlco = pd.to_numeric(geo["ch_dlco_(%_of_predicted)"], errors="raise")
    case = dlco < 60
    title_case = geo["title"].str.contains("sequela")
    if not (case == title_case).all() or not (title_case ^ geo["title"].str.contains("control")).all():
        raise ValueError(f"GSE228320: DLCO and title disagree for {list(geo.loc[case != title_case, 'gsm'])}")
    df = pd.DataFrame({
        "sample_id": geo["column"].values,
        "subject_id": ("GSE228320:" + geo["column"]).values,
        "lc_status": np.where(case, "LC", "COVID_recovered"),
        "case_definition": np.where(case, "pulmonary_sequelae", "not_applicable"),
        "control_type": np.where(case, "not_applicable", "COVID_recovered"),
        "sample_phase": "post_acute",
        "sample_matrix": "whole_blood",
        "time_anchor": "discharge",
        "acute_severity": "icu_or_ventilated",
        "hospitalised": True,
        "age": pd.to_numeric(geo["ch_age"]).values,
        "sex": geo["ch_Sex"].str.lower().values,
        "platform": "Illumina RNA-seq, whole blood (Tempus); deposited raw counts",
        "collection_country": "Spain",
        "gsm": geo["gsm"].values,
        "timepoint": "3_months_after_discharge",
        "dlco_pct_predicted": dlco.values,
        "notes": np.where(case, "survivor of COVID-19 ARDS with DLCO below 60% of predicted 3 months after discharge (pulmonary sequelae; not a Long COVID definition)",
                          "survivor of COVID-19 ARDS with DLCO 60% of predicted or above 3 months after discharge"),
    })
    out = finish(df, "GSE228320")  # keeps the schema columns, gsm and timepoint
    out["dlco_pct_predicted"] = dlco.values
    return out


def nanostring_qc(raw: pd.DataFrame, housekeeping: list) -> pd.DataFrame:
    """Per sample nSolver style QC of raw nCounter counts (samples x probes): scaling factors, linearity, detection limit."""
    pos = raw[list(NANOSTRING_POS_FM)].astype(float)
    neg = raw[[c for c in raw.columns if c.startswith("NEG_")]].astype(float)
    gm = lambda d: np.exp(np.log(d.clip(lower=1)).mean(axis=1))  # noqa: E731
    pos_gm, hk_gm = gm(pos), gm(raw[housekeeping].astype(float))
    x = np.log2(list(NANOSTRING_POS_FM.values()))
    r2 = pos.apply(lambda r: np.corrcoef(x, np.log2(r.clip(lower=1)))[0, 1] ** 2, axis=1)
    lod = neg.mean(axis=1) + 2 * neg.std(axis=1, ddof=1)
    q = pd.DataFrame({
        "pos_scaling_factor": pos_gm.mean() / pos_gm,
        "hk_scaling_factor": hk_gm.mean() / hk_gm,
        "pos_linearity_r2": r2,
        "pos_e": pos["POS_E"],
        "neg_mean_plus_2sd": lod,
    })
    reasons = pd.DataFrame({
        "positive control scaling factor outside 0.3 to 3": ~q["pos_scaling_factor"].between(0.3, 3),
        "housekeeping scaling factor outside 0.1 to 10": ~q["hk_scaling_factor"].between(0.1, 10),
        "positive control linearity R2 below 0.95": q["pos_linearity_r2"] < 0.95,
        "POS_E not above the negative control mean plus 2 SD": q["pos_e"] <= q["neg_mean_plus_2sd"],
    })
    q["flag"] = reasons.any(axis=1)
    q["flag_reason"] = reasons.apply(lambda r: "; ".join(r.index[r.values]), axis=1)
    return q


def corsip2026() -> pd.DataFrame:
    raw = pd.read_csv(CORSIP_REPO / "RawNanoStringData.csv").set_index("SampleID")
    subj = pd.read_csv(CORSIP_REPO / "SubjectsNoRep.csv")
    full = pd.read_csv(CORSIP_REPO / "CORSIP_Full.csv", low_memory=False)
    norm_genes = set(pd.read_csv(CORSIP_REPO / "NormNoRepNoHK.csv", nrows=0).columns) - {"SampleID"}
    controls = [c for c in raw.columns if c.startswith(("POS_", "NEG_"))]
    housekeeping = sorted(set(raw.columns) - set(controls) - norm_genes)  # removed from the normalised file as housekeeping
    if set(subj["SampleID"]) - set(raw.index):
        raise ValueError("CORSIP2026: labelled samples without raw counts")
    J = subj.merge(full, left_on="SampleID", right_on="SAMPLE_ID", how="left", suffixes=("", "_full"))
    if J["SAMPLE_ID"].isna().any() or not (J["LCStatus"] == J["LCStatus_full"]).all():
        raise ValueError("CORSIP2026: SubjectsNoRep and CORSIP_Full disagree on a sample or its LCStatus")
    raw_used = raw.loc[subj["SampleID"]]
    qc = nanostring_qc(raw_used, housekeeping)
    qc.rename_axis("sample_id").to_csv(CORSIP / "CORSIP2026_nanostring_qc.tsv", sep="\t", float_format="%.4f")
    m = symbol_map(pd.read_csv(ROOT / "data" / "ref" / "hgnc_complete_set.txt", sep="\t", dtype=str))
    genes = [c for c in raw.columns if c not in controls]
    ens = pd.Series(genes, index=genes).map(m)
    Y = raw_used[genes].T
    Y = Y[ens.notna().values].groupby(ens[ens.notna()].values).sum()
    Y.index.name = "gene_id"
    Y.to_csv(CORSIP / "CORSIP2026_nanostring_counts_ensembl.tsv.gz", sep="\t")
    print(f"CORSIP2026: {ens.notna().sum()} of {len(genes)} panel genes mapped ({len(housekeeping)} housekeeping); unmapped {sorted(ens.index[ens.isna()])}; "
          f"{len(Y)} Ensembl genes x {Y.shape[1]} samples; {int(qc['flag'].sum())} samples flagged at QC")
    tests = J[[f"CovPosTest{i}" for i in range(1, 9)]].notna().any(axis=1)
    n_react = J[[c for c in J.columns if c.startswith("CBS_Roche_N")]].apply(lambda r: (r.astype(str) == "Reactive").any(), axis=1)
    evidence = tests | n_react | J["CovPosSympDate"].notna()
    lc = J["LCStatus"] == "Case"
    first_test = pd.to_datetime(J[[f"CovPosTest{i}" for i in range(1, 9)]].stack(), errors="coerce").groupby(level=0).min().reindex(J.index)
    df = pd.DataFrame({
        "sample_id": J["SampleID"].values,
        "subject_id": ("CORSIP2026:" + J["StudyID"]).values,
        "lc_status": np.where(lc, "LC", "COVID_recovered"),
        "case_definition": np.where(lc, "WHO_2021", "not_applicable"),
        "control_type": np.where(lc, "not_applicable", np.where(evidence, "COVID_recovered", "convalescent_unverified")),
        "sample_phase": "post_acute",
        "sample_matrix": "whole_blood",
        "acute_severity": "ambulatory",
        "age": J["Age"].values,
        "sex": J["SexBirth"].str.lower().values,
        "platform": "NanoString nCounter PanCancer Immune Profiling Panel (770 genes), whole blood PAXgene; raw counts",
        "collection_country": "Canada",
        "timepoint": "final_visit_2023_01",
        "notes": np.where(lc, "LC: WHO clinical case definition from the 24 month questionnaire (paper Methods)",
                          np.where(evidence, "control: infected, no LC symptoms during the study",
                                   "control: no positive test, reactive nucleocapsid serology or symptom date in CORSIP_Full.csv; the participant answered never having had COVID-19"))
                 + "; PAXgene blood at the final visit (January 2023); first positive test " + first_test.dt.strftime("%Y-%m-%d").fillna("not recorded").values
                 + np.where(qc["flag"].values, "; NanoString QC flag: " + qc["flag_reason"].values, ""),
    })
    out = finish(df, "CORSIP2026")
    out["first_positive_test"] = first_test.dt.strftime("%Y-%m-%d").values
    out["nanostring_qc_flag"] = qc["flag"].values
    return out


def prjna1184005() -> pd.DataFrame:
    bs = pd.read_csv(PRJNA1184005_BIOSAMPLES, dtype=str)
    if set(bs["isolate"]) != set(PRJNA1184005_LABELS):
        raise ValueError(f"PRJNA1184005: isolates and author labels differ: {sorted(set(bs['isolate']) ^ set(PRJNA1184005_LABELS))}")
    if PRJNA1184005_COUNTS.exists():
        cols = set(pd.read_csv(PRJNA1184005_COUNTS, sep="\t", nrows=0).columns) - {"gene_id"}
        if cols != set(bs["isolate"]):
            raise ValueError(f"PRJNA1184005: count matrix columns and isolates differ: {sorted(cols ^ set(bs['isolate']))}")
    st = bs["isolate"].map(PRJNA1184005_LABELS)
    lc = (st == "LC").to_numpy()
    df = pd.DataFrame({
        "sample_id": bs["isolate"].values,
        "subject_id": ("PRJNA1184005:" + bs["isolate"]).values,
        "lc_status": st.values,
        "case_definition": np.where(lc, "WHO_2021", "not_applicable"),
        "control_type": np.where(lc, "not_applicable", "COVID_recovered"),
        "sample_phase": "post_acute",
        "sample_matrix": "PBMC",
        "age": pd.to_numeric(bs["age"]).values,
        "sex": bs["sex"].str.lower().values,
        "platform": "Illumina NovaSeq 6000, PBMC RNA-seq" + np.where(bs["isolate"] == "vacSG112", " (single-end run)", " (paired-end)"),
        "collection_country": "Australia",
        "timepoint": "single_visit",
        "notes": np.where(lc, "LC (WHO description, paper Methods); label from the author (S. Annesley, personal communication, 18 November 2026)",
                          "Post COVID (COVID_recovered): age and sex matched, recovered (paper Methods); label from the author (S. Annesley, personal "
                          "communication, 18 November 2026); control definition pending confirmation")
                 + "; BioSample " + bs["biosample"] + ", collected " + bs["collection_date"],
    })
    out = finish(df, "PRJNA1184005")
    out["label_source"] = "author: S. Annesley, personal communication, 2026-11-18"
    out["control_definition_status"] = np.where(lc, "not_applicable", "pending_confirmation")
    out["collection_date"] = bs["collection_date"].values
    out["biosample"] = bs["biosample"].values
    if PRJNA1184005_QC.exists():  # GSE224615 rule (below 40%, or more than 2.5 MAD below the median): a flag, not an exclusion
        qc = pd.read_csv(PRJNA1184005_QC, sep="\t", dtype={"sample": str}).set_index("sample").reindex(bs["isolate"])
        out["salmon_mapping_rate"] = qc["mapping_rate"].values
        out["salmon_qc_flag"] = qc["flag"].astype(str).str.lower().eq("true").values
        out["notes"] = out["notes"] + np.where(out["salmon_qc_flag"], "; salmon QC flag: " + qc["flag_reason"].fillna("").values, "")
    return out


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
    if sys.argv[1:] == ["GSE228320"]:
        built = {"GSE228320": gse228320()}
    elif sys.argv[1:] == ["CORSIP2026"]:
        built = {"CORSIP2026": corsip2026()}
    elif sys.argv[1:] == ["PRJNA1184005"]:
        built = {"PRJNA1184005": prjna1184005()}
    else:
        built = {"GSE297499": gse297499(), **prjna1149196()}
    for sid, df in built.items():
        df.to_csv(META / f"{sid}_meta.csv", index=False)
        print(f"{sid}: {len(df)} samples, {df['subject_id'].nunique()} subjects; lc_status {df['lc_status'].value_counts().to_dict()}")
    split_counts({k: v for k, v in built.items() if k.startswith("PRJNA1149196")})


if __name__ == "__main__":
    main()
