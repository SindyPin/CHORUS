"""Round2 step 0: labels, clinical variables, anchors and partners.

Reads the cohort metadata on disk (data/meta/<cohort>_meta.csv, written in Round1), the curated definitions below
(taken from the Round1 cohort cards: Round1/docs/*.md and data/meta/label_sources.md, curation_notes.md) and, where
they are open, the clinical or symptom tables downloaded into data/clinical/ (see CLINICAL). Writes:

  output/labels_audit.csv          one row per cohort on disk: layer, tissue, platform, LC definition and label
                                   source, comparison group, time since infection, sample timing, participants per
                                   class and phase, clinical variables, flags, which tasks it can serve
  output/clinical_variables.csv    one row per cohort and variable: kind, source, participants recorded
  output/anchors_and_partners.csv  per task and layer (horizontal): anchor, partner, class counts, shared features;
                                   and every case where the same participants have two layers (vertical)
  output/clinical_links.csv        downloaded clinical tables: rows, rows linked to our participants
  inputs.csv                       files read, with SHA-256

Rules (Round2/PROTOCOL.md): one sample per participant; anchor = largest cohort for the task and layer with at
least 10 participants per class; partner = next largest of the same layer (also at least 10 per class);
GSE270045, GSE228320 and CORSIP2026 are listed but not used as anchor or partner. A label is never guessed: a
participant without a stated label is "unknown".

  python Round2/0_Labels/build_labels.py            (from Code/)
"""

from __future__ import annotations

import glob
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

CODE = Path(__file__).resolve().parents[2]
STEP = Path(__file__).resolve().parent
OUT = STEP / "output"
META = CODE / "data" / "meta"
SEED = 20261120
MIN_PER_CLASS = 10
LABEL_FLAGS = {
    "GSE270045": "labels confounded with collection",
    "GSE228320": "non standard case definition (pulmonary sequelae after ARDS)",
    "CORSIP2026": "gene panel only (NanoString, 770 genes)",
}
LAYER = {  # broad layer used for "same layer" (anchor and partner)
    "bulk_rna": "transcriptomics", "olink": "proteomics", "somascan": "proteomics", "ms_proteomics": "proteomics",
    "metabolomics": "metabolomics", "methylation": "methylation",
}
# participants whose label sources disagree: unknown in Round2 (Round1 metadata unchanged)
UNKNOWN_IN_ROUND2 = {
    "GSE226260_WB": {"GSE226260:BI4036": "GEO titles and the paper's Medication sheet say Recovered; its Demographics and Symptom List sheets say PASC",
                     "GSE226260:BI4002": "GEO titles and the paper's Medication sheet say PASC; its Demographics and Symptom List sheets say Recovered"},
}
# phase for the metadata files written before sample_phase was a schema column
DEFAULT_PHASE = {"GSE215865": "acute", "GSE224615": "post_acute", "GSE226260_PBMC": "post_acute",
                 "GSE251849": "post_acute", "GSE270045": "post_acute", "GSE226260_WB": "post_acute"}

# Curated per cohort (sources in the Round1 docs named in "source_doc"). lc: how LC is defined; label: where the
# per-participant label comes from; comparison: how the comparison group is defined; time: time since infection.
C = {}


def cur(study, paper, lc, label, comparison, time, doc, doubt=""):
    C[study] = dict(paper=paper, lc_definition=lc, label_source=label, comparison_group=comparison,
                    time_since_infection=time, source_doc=doc, doubts=doubt)


R1 = "Round1/docs/"
# transcriptomics, LC versus recovered
cur("PRJNA1149196_SA", "Sinclair 2024 Nat Microbiol, PMID 39478108 (cohort 1, South Australia)",
    "PASC: persistent symptoms more than 12 weeks after infection with ongoing discomfort at sampling (study_specific)",
    "BioSample attribute 'disease' (PASC / Recovered / Control)", "Recovered: infected, no persistent symptoms more than 12 weeks after infection; Control: uninfected (never_infected)",
    "nominal visits 44 and 68 weeks after infection (BioSample 'Timepoint')", "data/meta/PRJNA1149196_SA_meta.csv notes; " + R1 + "dataset_sweep_2026-11.md",
    "repeated visits per participant; days are the nominal visit, not a per sample date")
cur("PRJNA1149196_QLD", "Sinclair 2024 Nat Microbiol, PMID 39478108 (cohort 2, Queensland)",
    "PASC-CVS: PASC with cardiovascular symptoms (study_specific)", "BioSample attribute 'disease'",
    "Recovered: infected, no persistent symptoms more than 12 weeks after infection", "1 year after infection (BioSample); paper more than 520 days",
    "data/meta/PRJNA1149196_QLD_meta.csv notes", "")
cur("GSE226260_WB", "GSE226260 (Aid 2025), whole blood subset",
    "GEO 'current long_covid_status' = PASC (study_specific)", "GEO characteristic 'current long_covid_status'",
    "GEO 'Recovered' (COVID_recovered)", "days from the sample title, 15 to 1678 (Acute 15-30, 3-6 M, 6-12 M, >12 M); anchor (onset or test) not stated",
    "data/meta/curation_notes.md", "the WB and PBMC subsets use different id systems, so people who gave both cannot be linked; 4 acute draws (all PASC)")
cur("GSE226260_PBMC", "GSE226260 (Aid 2025), PBMC subset",
    "GEO 'pasc status' = PASC (study_specific)", "GEO characteristic 'pasc status'", "GEO 'NOPASC' (COVID_recovered)",
    "not given: GEO has collection dates only", "data/meta/curation_notes.md",
    "sample timing taken as post-acute from the PASC / NOPASC status; no infection date per sample")
cur("GSE224615", "Yin 2024 (LIINC), GSE224615", "GEO 'lc status' = LC; WHO_2021 (paper)", "GEO characteristic 'lc status'",
    "GEO 'Non-LC' (COVID_recovered)", "8 months after infection (cohort visit; 243 days set for every sample)",
    "data/meta/curation_notes.md", "paper cohort is 27 / 16; GEO has 23 / 13")
cur("GSE251849", "Greene 2024", "GEO group 'Long COVID' or 'Long COVID brain fog'; symptoms over 12 weeks (paper)",
    "count matrix id prefix, checked against GEO 'group' counts", "GEO 'COVID convalescent' (COVID_recovered); 'Healthy' never_infected",
    "not given per sample", "data/meta/curation_notes.md", "matrix columns cannot be mapped to GSM, so sex is unknown per column")
cur("GSE270045", "Shahbaz 2025 Cell Rep Med, PMID 41205594", "LC with ME/CFS, at least 12 months after infection (paper)",
    "GEO title; paper says the 'Healthy Control' samples are recovered", "COVID-19 recovered without symptoms (paper)",
    "at least 12 months (paper); no per sample time", "data/meta/label_sources.md",
    "flagged: labels confounded with collection; 19 LC in GEO against 18 in the paper")
cur("GSE228320", "GSE228320 (post-ARDS cohort)", "pulmonary sequelae: DLCO below 60% of predicted 3 months after discharge (non standard)",
    "GEO characteristics (DLCO group)", "ARDS survivors with DLCO 60% of predicted or above", "3 months after hospital discharge",
    "data/meta/GSE228320_meta.csv notes", "flagged: non standard case definition; all participants had ARDS")
cur("GSE297499", "Orlinick 2026 J Infect Dis, PMID 40920586", "self-reported neuropsychiatric Long COVID, WHO definition, documented infection",
    "GEO characteristic 'treatment' (NP-Long COVID / Control)", "women with documented SARS-CoV-2 infection without persistent symptoms",
    "not given per sample", R1 + "dataset_sweep_2026-11.md", "")
cur("PRJNA1184005", "Missailidis 2024 Front Immunol, PMID 39691709", "LC (WHO description, paper Methods)",
    "author (S. Annesley, personal communication, 18 November 2026)", "post COVID (COVID_recovered); control definition pending confirmation",
    "not given per sample (collection 2021)", "data/meta/label_sources.md", "control definition pending confirmation from the author")
cur("CORSIP2026", "Espin 2026 Front Cell Infect Microbiol, PMID 42328163", "WHO clinical case definition from the 24 month questionnaire",
    "repository participant table (github.com/estefespin/CORSIP_omics)", "infected healthcare workers without LC at 24 months",
    "final visit January 2023 (about 11 months after a first positive test in February 2022)", "data/meta/CORSIP2026_meta.csv notes",
    "flagged: gene panel only")
cur("GSE215865", "Thompson 2023 (Mount Sinai), GSE215865", "TACO table column Long_COVID = 1 (definition not given in the TACO table)",
    "TACO clinical table (SindyPin/TACO, commit cac2afb) joined on GEO 'subject id'; original Synapse syn35874390 not checked",
    "Long_COVID = 0 and COVID positive in GEO (COVID_recovered); 0 and negative is never_infected_or_negative",
    "acute hospital draws (T0 to T13, days since the first study sample)", "data/meta/label_sources.md",
    "labels from a secondary table, not yet checked against Synapse; 1 subject LC in TACO but COVID negative in GEO (unknown); 20 subjects without a label")
for s, p in (("GSE169687", "Ryan 2022 BMC Med, PMID 35027067"), ("GSE267625", "Blankestijn 2024 JACI, PMID 38830512"),
             ("GSE222253", "Front Immunol 2023, PMID 37680625")):
    cur(s, p, "not in hand (requested from the authors)", "none on disk", "see data/meta/label_sources.md", "see data/meta/label_sources.md",
        "data/meta/label_sources.md", "no per participant LC label: every participant unknown")
cur("Fineschi2025_PBMC", "Fineschi 2025 Front Immunol, PMID 40519915 (PBMC RNA-seq, figshare 10.17044/scilifelab.28832492)",
    "PAT = patients meeting WHO post COVID-19 condition (as Fineschi2026)", "SamplesPC.txt 'Condition' (PAT / CTL)",
    "CTL = infected, recovered fully without persistent symptoms", "about 3 years (blood April-May 2023)", R1 + "proteomics_cohorts.md",
    "only the label file is on disk; the count matrix (salmon.merged.gene_counts.tsv) is not downloaded")
# proteomics
cur("Buggert2025", "Gao 2025 Nat Immunol (Buggert lab)", "UK: NICE NG188, symptoms at least 12 weeks; Sweden: no case definition given",
    "deposited metadata", "UK healthy convalescent individuals; the 30 Swedish controls were not profiled", "UK medians 416 days (cases) and 268 days (controls) from diagnosis; none per sample",
    R1 + "proteomics_cohorts.md", "the 95 Swedish participants are all LC, so country and LC status overlap; 2 plate samples without metadata are unknown")
cur("Fineschi2026", "Fineschi 2026 Front Immunol, PMID 41808838", "WHO post COVID-19 condition", "SampleID prefix PA / KO (checked against the PBMC deposit's Condition)",
    "infected in the same period, recovered fully without persistent symptoms", "mean symptom duration 34 months; no per sample date",
    R1 + "proteomics_cohorts.md", "")
cur("Chowdhury2026", "Chowdhury 2026 Sci Rep, PMID 41998026", "PCC: persistent symptoms consistent with WHO criteria", "PRIDE SDRF",
    "convalescent: PCR negative 1 month after infection, symptom free at 3 months", "3 months after a PCR confirmed infection",
    R1 + "proteomics_cohorts.md", "5 runs coded Lc in file names but C in the SDRF (SDRF used); overlap with Chowdhury2024 cannot be excluded")
cur("Chowdhury2024", "Chowdhury 2024 Front Physiol", "PCC: at least one candidate symptom 3 months after a positive PCR (study_specific)",
    "PRIDE SDRF", "9 healthy hospital co-workers called convalescent without a per person infection record (convalescent_unverified)",
    "PCC duration from 1 month (group level)", R1 + "proteomics_cohorts.md", "controls' infection not documented per person; exercise intervention between the two LC draws")
cur("Woodruff2023", "Woodruff 2023 Nat Commun, PMID 37452024", "Emory PASC clinic: new or worsening symptoms more than 14 days after onset (study_specific)",
    "deposit column group_hr", "CR: uncomplicated recoveries sampled at similar intervals", "days post symptom onset per sample (PASC median 111, CR median 99)",
    R1 + "proteomics_cohorts.md", "38 of 93 PASC samples are under 90 days after onset; sex recorded for 2 of 93 PASC")
cur("Cervia2024_SomaScan", "Cervia-Hasler 2024 Science", "one or more persisting COVID-19 related symptoms at 6 months, smell or taste only excluded (PACS_6M_woDys)",
    "deposited PACS_6M_woDys", "PACS_6M_woDys = 0: recovered within a month, 'recovered Long Covid' and smell or taste only, mixed",
    "acute samples and 6 months after infection (per sample days)", R1 + "proteomics_cohorts.md", "the recovered group mixes three groups the deposit cannot separate")
C["Cervia2024_MS"] = dict(C["Cervia2024_SomaScan"], time_since_infection="6 months after infection (per sample days)")
cur("Su2022", "Su 2022 Cell (INCOV)", "at least one PASC symptom at T3 by interview and chart review (study_specific)", "authors' PASC groups (Mendeley clinical tables)",
    "no PASC at T3", "T1 diagnosis, T2 about a week later, T3 about 2 to 3 months after onset", R1 + "proteomics_cohorts.md",
    "single assessment at about 2 to 3 months, no duration requirement; values are authors' age, sex and BMI adjusted Z scores")
cur("Bansal2025", "Bansal 2026 Commun Med, PMID 41975235", "WHO_2021: symptoms 3 months after infection lasting at least 2 months",
    "deposited 'Long_COVID status'", "Recovered: symptoms shorter than 90 days", "6 to 9 months after the first infection, before vaccination",
    R1 + "proteomics_cohorts.md", "")
cur("Kovarik2023_plasma", "Kovarik 2023 iScience", "chronic fatigue and/or lack of concentration with at least one more chronic symptom, at least 3 months after infection",
    "MaxQuant columns matched to Table S1 labels", "symptom free at the blood draw, at least 3 months after infection, anti-N and anti-S positive",
    "at least 3 months; no per sample time", R1 + "proteomics_cohorts.md", "")
cur("Captur2022", "Captur 2022 eBioMedicine (COVIDsortium)", "1 or more persistent symptom out to 12 months (paper); not deposited",
    "none deposited (data access portal)", "not deposited", "serial samples weeks -1 to 6 around the first positive test",
    R1 + "proteomics_cohorts.md", "no per participant label, timepoint or infection status: every sample unknown")
# metabolomics
C["Su2022_metabolomics"] = dict(C["Su2022"])
cur("Martinez2024", "Martinez 2024 J Proteome Res, PMID 38566450", "Workbench factor 'Long COVID' (paper behind a subscription; case definition not read)",
    "Metabolomics Workbench ST003103 sample factors", "'Recovered from acute COVID-19' (fully recovered from mild COVID-19)", "12+ weeks (catalogue)",
    R1 + "metabolomics_cohorts.md", "case definition not read")
cur("Lopez2023", "Lopez-Hernandez 2023 Sci Rep, PMID 37528111", "at least one persistent symptom (class A fewer than five, class B five or more)",
    "Mendeley deposit group", "recovered: no persistent symptoms; CONTROL: tested negative in 2020 (never_infected)",
    "about two years after infection; acute sample of 22 patients", R1 + "metabolomics_cohorts.md", "deposit counts differ from the paper by one per group")
cur("Saito2024", "Saito 2024 Front Immunol, PMID 38304426", "LC meeting ME/CFS criteria about 12 months after infection", "data file group letter (L)",
    "R: recovered without symptoms (group letter C taken as R by elimination); HC healthy", "about 12 months; acute ICU samples have no outcome",
    R1 + "metabolomics_cohorts.md", "group C taken as R by elimination; clinical table not linkable per sample")
cur("Holmes2021", "Holmes 2021 J Proteome Res, PMID 34009992", "not resolved: four symptom composites with undocumented coding",
    "none usable", "no controls deposited", "about 3 months", R1 + "metabolomics_cohorts.md", "binary composites disagree (23 against 18): every participant unknown")
# methylation
cur("GSE197152", "Balnis 2022 Clin Epigenetics, PMID 35871090", "all 15 patients seen one year after discharge expressed multiple PASC symptoms (paper, group statement)",
    "paper text (no per sample field)", "none (no recovered group)", "one year after infection", R1 + "metabolomics_cohorts.md",
    "the LC label rests on a statement about the whole group, not a per sample field")
cur("GSE247869", "Calzari 2024 Clin Epigenetics, PMID 39164752", "PASC subgroup of 28 (at least one persistent symptom), not per sample",
    "none public", "not public", "6 months", R1 + "metabolomics_cohorts.md", "every participant unknown")
cur("GSE174818", "Balnis 2021 Clin Epigenetics, PMID 34034806", "none (acute)", "none", "non-COVID-19 patients", "during the acute illness",
    R1 + "metabolomics_cohorts.md", "no LC outcome")
NO_LC = ("GSE198449 COMBAT2022_WB GSE161731 GSE152641 GSE300129 GSE243219 GSE172114 GSE171110 GSE213313 GSE152418 "
         "Filbin2021_MGH Gisby2021_plasma GSE207015 Zoodsma2022_Breda Zoodsma2022_Hannover Zoodsma2022_Nijmegen").split()
for s in NO_LC:
    cur(s, "see Round1/docs/acute_covid_cohorts.md" if not s.startswith("Zoodsma") else "Zoodsma 2022 Front Immunol",
        "none: no LC outcome recorded", "none", "acute COVID-19 or never infected (acute COVID validation cohorts)" if not s.startswith("Zoodsma")
        else "post-COVID-19 convalescents (Hannover) were largely symptom free; no LC split", "acute phase (Hannover convalescents 13 to 75 days after symptom resolution)",
        R1 + ("acute_covid_cohorts.md" if not s.startswith("Zoodsma") else "proteomics_cohorts.md"), "no LC outcome: cannot serve the three tasks")
EXTRA = {  # cohorts with labels on disk but no ClinicalSchema metadata file
    "GSE197152": dict(layer="methylation", tissue="whole_blood", platform="Illumina HumanMethylationEPIC (GPL21145)", LC_post=15),
    "GSE247869": dict(layer="methylation", tissue="whole_blood", platform="Illumina HumanMethylationEPIC (GPL21145)", unknown=96),
    "GSE174818": dict(layer="methylation", tissue="whole_blood", platform="Illumina HumanMethylationEPIC (GPL21145)", unknown=128),
    "GSE169687": dict(layer="bulk_rna", tissue="PBMC", platform="Illumina NovaSeq 6000", unknown=69, never=14),
    "GSE267625": dict(layer="bulk_rna", tissue="whole_blood", platform="Illumina NovaSeq 6000", unknown=95),
    "GSE222253": dict(layer="bulk_rna", tissue="whole_blood", platform="Illumina NovaSeq 6000", unknown=24),
}

# clinical and symptom variables: kind; "defines_label" = the variable is part of the LC definition
KIND = {
    "age": "demographic", "sex": "demographic", "bmi": "demographic", "race": "demographic", "age_band": "demographic",
    "age_category": "demographic",
    "acute_severity": "severity", "hospitalised": "severity", "covid19_severity_max": "severity",
    "covid19_severity_grade_max": "severity", "peak_acute_wos": "severity", "wos_at_draw": "severity",
    "who_score_authors": "severity", "disease_severity_score_authors": "severity", "disease_status_subgroup": "severity",
    "severity_original": "severity", "encounter_location": "severity",
    "comorbidity_count": "comorbidity", "comorbidity": "comorbidity",
    "vaccinated_before_infection": "vaccination", "vaccinated_before_sampling": "vaccination",
    "days_since_infection": "timing", "days_since_symptom_resolution": "timing",
    "symptom_domains": "symptoms", "n_pasc_symptoms": "symptoms", "pasc_symptoms": "symptoms", "pasc_groups": "symptoms",
    "pasc_count_group": "symptoms", "lc_subgroup": "symptoms", "composite_binary": "symptoms",
    "composite_binary_rec": "symptoms", "pacs_6m_wodys": "symptoms", "pacs_12m": "symptoms",
    "borg_cr10": "score", "dyspnoea12": "score", "nijmegen": "score", "dlco_pct_predicted": "score",
    "dexamethasone": "treatment", "intervention_arm": "treatment", "variant_era": "infection", "pregnancy": "other",
    "antibody_S1_ratio": "serology", "antibody_NCP_ratio": "serology",
}
DEFINES_LABEL = {"symptom_domains", "n_pasc_symptoms", "pasc_symptoms", "pasc_groups", "pasc_count_group", "lc_subgroup",
                 "composite_binary", "composite_binary_rec", "pacs_6m_wodys", "pacs_12m", "dlco_pct_predicted",
                 "borg_cr10", "dyspnoea12", "nijmegen"}  # symptom scores count as defining the label

# open per participant clinical or symptom tables not in the Round1 metadata: (cohort, file, url, licence). Files
# with a url are downloaded into data/clinical/<cohort>/ when missing; the others were already on disk from Round1.
SPRINGER = "https://static-content.springer.com/esm/art%3A10.1038%2F"
CLINICAL = [
    ("GSE224615", "data/clinical/GSE224615/41590_2023_1724_MOESM1_ESM.pdf",
     SPRINGER + "s41590-023-01724-6/MediaObjects/41590_2023_1724_MOESM1_ESM.pdf", "CC BY 4.0 (Yin 2024 Nat Immunol, Supplementary Table 3)"),
    ("GSE226260_WB", "data/clinical/GSE226260/41590_2025_2353_MOESM4_ESM.xlsx",
     SPRINGER + "s41590-025-02353-x/MediaObjects/41590_2025_2353_MOESM4_ESM.xlsx", "CC BY 4.0 (Nat Immunol 2025, PMID 41388153; 2023-2024 cohort)"),
    ("GSE226260_PBMC", "data/clinical/GSE226260/41590_2025_2353_MOESM2_ESM.xlsx",
     SPRINGER + "s41590-025-02353-x/MediaObjects/41590_2025_2353_MOESM2_ESM.xlsx", "CC BY 4.0 (Nat Immunol 2025, PMID 41388153; 2020-2021 cohort)"),
    ("GSE226260_PBMC", "data/clinical/GSE226260/41590_2025_2353_MOESM3_ESM.xlsx",
     SPRINGER + "s41590-025-02353-x/MediaObjects/41590_2025_2353_MOESM3_ESM.xlsx", "CC BY 4.0 (sample sheet; downloaded, not linked: see README)"),
    ("Lopez2023", "data/metabolomics/Lopez2023/raw/DATA AVAILABILITY.xls", "", "CC BY 4.0 (Mendeley 10.17632/8zfdjsypd8.1; on disk from Round1)"),
    ("Fineschi2026", "data/counts/Fineschi2025_check/SamplesPC.txt", "", "CC BY 4.0 (figshare 10.17044/scilifelab.28832492; on disk from Round1)"),
    ("CORSIP2026", "data/counts/CORSIP2026/raw/CORSIP_omics/CORSIP_Full.csv", "", "CC0 (github.com/estefespin/CORSIP_omics; on disk from Round1)"),
]

INPUTS: dict = {}


def fetch(path: str, url: str) -> Path:
    p = CODE / path
    if not p.exists() and url:
        import urllib.request

        p.parent.mkdir(parents=True, exist_ok=True)
        urllib.request.urlretrieve(url, p)
    INPUTS[path] = p
    return p


def yes(v) -> bool:
    return str(v).strip().lower() in ("1", "1.0", "yes", "true")


def clinical_links(parts: dict) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Link each open table to our participants. Returns (links, variables); the linked tables are written to
    data/clinical/<cohort>/<cohort>_clinical_linked.csv (not committed)."""
    links, var_rows = [], []

    def record(cohort, path, url, licence, T, key, label_col, label_map, variables, note=""):
        m, ph, P = parts[cohort]
        lab = m.drop_duplicates("subject_id").set_index("subject_id")["lc_status"].astype(str)
        T = T.copy()
        T.columns = [str(c).strip() for c in T.columns]
        T["subject_id"] = T[key]
        linked = T[T["subject_id"].isin(lab.index)]
        agree = disagree = 0
        if label_col:
            ours = linked["subject_id"].map(lab)
            theirs = linked[label_col].map(label_map)
            ok = theirs.notna() & ours.isin(["LC", "COVID_recovered"])
            agree, disagree = int((ours[ok] == theirs[ok]).sum()), int((ours[ok] != theirs[ok]).sum())
        out = CODE / "data/clinical" / cohort / f"{cohort}_clinical_linked.csv"
        out.parent.mkdir(parents=True, exist_ok=True)
        linked.to_csv(out, index=False)
        labelled = P["LC_post"] | P["REC_post"] | P["LC_acute"] | P["REC_acute"]
        links.append(dict(study_id=cohort, file=path, url=url, licence=licence, sha256=sha(CODE / path), rows=len(T),
                          rows_linked=len(linked), labelled_participants=len(labelled),
                          labelled_participants_linked=int(linked.loc[linked["subject_id"].isin(labelled), "subject_id"].nunique()),
                          label_agrees=agree, label_disagrees=disagree, linked_table=str(out.relative_to(CODE)), note=note))
        for v, (col, kind, defines) in variables.items():
            okv = linked[col].notna() & ~linked[col].astype(str).str.strip().isin(["", "nan", "NA"])
            var_rows.append(dict(study_id=cohort, variable=v, kind=kind, defines_label=defines, source=path,
                                 n_participants_recorded=int(linked.loc[okv, "subject_id"].nunique()),
                                 n_labelled_participants_recorded=int(linked.loc[okv & linked["subject_id"].isin(labelled), "subject_id"].nunique()),
                                 n_labelled_participants=len(labelled)))

    for cohort, path, url, licence in CLINICAL:
        p = fetch(path, url)
        if cohort == "GSE224615":  # Supplementary Table 3, pages 3-4 of the PDF: Person ID, age, race, status, sex, ...
            import re

            import pypdf

            text = "\n".join(pg.extract_text() for pg in pypdf.PdfReader(p).pages[2:4])
            rx = re.compile(r"^(\d+)[*&^]? (\d+) (\S+) (R|LC) (Female|Male) (Yes|No) (\S+)\s*$")
            T = pd.DataFrame([g.groups() for g in map(rx.match, text.splitlines()) if g],
                             columns=["person_id", "age", "race", "status", "sex", "hospitalized", "who_category_acute"])
            titles = read(META / "labels_geo_and_taco.csv", header=None, names=["study", "gsm", "title", "subset", "lc"])
            titles = titles[titles["study"] == "GSE224615"]
            gsm2person = dict(zip(titles["gsm"], titles["title"].str.split("_").str[0]))
            m = parts[cohort][0]
            person2subject = {gsm2person[g]: s for g, s in zip(m["sample_id"], m["subject_id"]) if g in gsm2person}
            T["subject_id_key"] = T["person_id"].map(person2subject)
            record(cohort, path, url, licence, T, "subject_id_key", "status", {"R": "COVID_recovered", "LC": "LC"},
                   {"age": ("age", "demographic", False), "race": ("race", "demographic", False),
                    "hospitalised": ("hospitalized", "severity", False), "acute_severity_who": ("who_category_acute", "severity", False)},
                   "PDF table parsed; joined by the person id in the GEO sample title (e.g. 942_Non-PASC_WT)")
        elif cohort == "GSE226260_WB":
            D = pd.read_excel(p, sheet_name="Demographics", header=2)
            S = pd.read_excel(p, sheet_name="Symptom List", header=2)
            S["SubjectID"] = S["Biorep ID"].astype(str).str.replace("-", "", regex=False)
            T = D.merge(S[["SubjectID", "Symptom_count", "symptom_list"]], on="SubjectID", how="left")
            T["key"] = "GSE226260:" + T["SubjectID"].astype(str)
            record(cohort, path, url, licence, T, "key", "PASC", {"PASC": "LC", "Recovered": "COVID_recovered"},
                   {"age": ("Age", "demographic", False), "sex": ("Sex at birth", "demographic", False), "race": ("Race", "demographic", False),
                    "ethnicity": ("Ethnicity", "demographic", False), "days_since_infection_at_sample": ("TPI", "timing", False),
                    "preexisting_conditions": ("preexisting_conditions_cc", "comorbidity", False),
                    "symptom_count": ("Symptom_count", "symptoms", True), "symptom_list": ("symptom_list", "symptoms", True)},
                   "sheets Demographics and Symptom List; SubjectID BI4001 (symptom sheet BI-4001); BI4036 and BI4002 disagree with GEO and are unknown in Round2")
        elif cohort == "GSE226260_PBMC" and "MOESM2" in path:
            T = pd.read_excel(p, header=2)
            T["key"] = "GSE226260:" + T["SN"].astype(str).str.replace(r"\.0$", "", regex=True)
            record(cohort, path, url, licence, T, "key", "Evidence of PASC", {"PASC": "LC", "Recovered": "COVID_recovered"},
                   {"age": ("Age (12/31/22)", "demographic", False), "sex": ("Gender", "demographic", False), "race": ("Race", "demographic", False),
                    "ethnicity": ("Ethnicity", "demographic", False), "vaccine_product": ("Vaccine product", "vaccination", False)},
                   "joined on SN (GEO subject id)")
        elif cohort == "Lopez2023":
            T = pd.read_excel(p, engine="xlrd", dtype={"Folio": str})
            T.columns = [str(c).strip() for c in T.columns]
            folio, group = T["Folio"].astype(str).str.strip(), T["GROUP"].astype(str).str.strip()
            # the Round1 builder's id rule (Round1/scripts/metabolomics/build_cohorts.py)
            T["sample_id"] = np.where(group == "CONTROL", folio, np.where(group == "COVID-19", "P" + folio.str.replace(".", "_", regex=False) + "_acute", "P" + folio + "_post"))
            m = parts[cohort][0]
            T["key"] = T["sample_id"].map(dict(zip(m["sample_id"].astype(str), m["subject_id"])))
            sym = ["FEVER", "ARTHRALGIA", "MYALGIA", "CEPHALEA", "RHINORRHEA", "COUGH", "DYSPNEA", "THORACIC PAIN", "ANOSMYA",
                   "DYSGEUSIA", "DIARRHEA", "LOSS OF CONCENTRATION", "FATIGUE", "VISUAL ALTERATION", "ANXIETY", "CONJUNCTIVITIS",
                   "DEPRESSION", "ABDOMINAL PAIN", "ODYNOPHAGIA"]
            v = {"age": ("AGE", "demographic", False), "sex_code": ("SEX", "demographic", False), "smoking": ("SMOKING", "other", False),
                 "t2dm": ("T2DM", "comorbidity", False), "hypertension": ("ARTERIAL HYPERTENSION", "comorbidity", False),
                 "obesity": ("OBESITY", "comorbidity", False)}
            v.update({f"symptom_{s.lower().replace(' ', '_')}": (s, "symptoms", True) for s in sym})
            record(cohort, path, url, licence, T, "key", "", {}, v,
                   "same file as the metabolite values (Folio); symptom and sex coding not documented (every never infected control has FEVER = 1, so 1 does not mean present); acute or persistent not stated")
        elif cohort == "Fineschi2026":
            T = pd.read_csv(p, sep="\t")
            T["key"] = "Fineschi2026_" + T["Sample"].astype(str).str.replace(".", "-", regex=False)
            record(cohort, path, url, licence, T, "key", "Condition", {"PAT": "LC", "CTL": "COVID_recovered"},
                   {"infection_wave": ("Wave", "infection", False), "post_covid_symptom_severity_score": ("Score", "symptoms", True),
                    "dyspnoea_item": ("Dyspnoea", "symptoms", True)},
                   "PBMC study label file (Fineschi 2025); the same people as the Olink samples; same draw not stated")
        elif cohort == "CORSIP2026":
            T = pd.read_csv(p, low_memory=False)
            m = parts[cohort][0]
            T["key"] = T["SAMPLE_ID"].map(dict(zip(m["sample_id"], m["subject_id"])))
            symp = [c for c in T.columns if c.startswith("Symp")]
            record(cohort, path, url, licence, T, "key", "LCStatus", {"Case": "LC", "Control": "COVID_recovered"},
                   {"vaccinated": ("Vaccinated", "vaccination", False), "hospitalised_ever_18m": ("EverCovHosp18", "severity", False),
                    "weight": ("Weight", "demographic", False), "height": ("Height", "demographic", False),
                    f"symptoms ({len(symp)} columns, baseline to 18 months)": (symp[0], "symptoms", True)},
                   "flagged cohort (gene panel only); joined on SAMPLE_ID")
    return pd.DataFrame(links), pd.DataFrame(var_rows)


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def read(path: Path, **kw) -> pd.DataFrame:
    INPUTS[str(path.relative_to(CODE))] = path
    return pd.read_csv(path, low_memory=False, **kw)


def phase_of(study: str, m: pd.DataFrame) -> pd.Series:
    if "sample_phase" in m and m["sample_phase"].notna().any():
        ph = m["sample_phase"].fillna("unknown").astype(str)
    else:
        ph = pd.Series(DEFAULT_PHASE.get(study, "unknown"), index=m.index)
    if study == "GSE226260_WB":  # the four 'Acute' draws (days 15 to 30)
        ph = ph.where(m["timepoint"].astype(str) != "Acute", "acute")
    lab = m["lc_status"].astype(str)
    return ph.where(~lab.isin(["never_infected"]), "not_applicable")


def participants(m: pd.DataFrame, phase: pd.Series, ph: str, lab: str) -> set:
    return set(m.loc[(phase == ph) & (m["lc_status"].astype(str) == lab), "subject_id"].astype(str))


def config_entries() -> dict:
    seen = {}
    for f in sorted(glob.glob(str(CODE / "Round1/configs/*.yaml"))):
        for c in yaml.safe_load(open(f))["cohorts"]:
            if (CODE / c["matrix"]).exists():
                seen.setdefault(c["study_id"], (f, c))
    return seen


def features(study: str, entries: dict, layer: str) -> set:
    """Feature keys of the cohort as the Round1 loader reads them (genes as Ensembl ids, proteins as UniProt keys,
    metabolites as shared identifiers); before any expression filter."""
    from chorus.io import Cohort

    f, c = entries[study]
    INPUTS[str(Path(c["matrix"]))] = CODE / c["matrix"]
    co = Cohort.from_files(c["study_id"], layer, str(CODE / c["matrix"]), str(CODE / c["meta"]),
                           features_in_rows=c.get("features_in_rows", True), is_counts=c.get("is_counts", False),
                           feature_space=c.get("feature_space", "ensembl_gene"))
    return set(map(str, co.X.columns[co.X.notna().any(axis=0)]))


def audit() -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    rows, clin, parts = [], [], {}
    for f in sorted(META.glob("*_meta.csv")):
        study = f.name[:-9]
        if study not in C:
            continue
        m = read(f)
        for subj in UNKNOWN_IN_ROUND2.get(study, {}):
            m.loc[m["subject_id"].astype(str) == subj, "lc_status"] = "unknown"
        ph = phase_of(study, m)
        P = {k: participants(m, ph, p, lab) for k, (p, lab) in {
            "LC_post": ("post_acute", "LC"), "REC_post": ("post_acute", "COVID_recovered"),
            "LC_acute": ("acute", "LC"), "REC_acute": ("acute", "COVID_recovered")}.items()}
        P["never"] = set(m.loc[m["lc_status"].astype(str) == "never_infected", "subject_id"].astype(str))
        known = set().union(*P.values())
        P["unknown"] = set(m["subject_id"].astype(str)) - known
        parts[study] = (m, ph, P)
        lab_ids = P["LC_post"] | P["REC_post"] | P["LC_acute"] | P["REC_acute"]
        vars_recorded = []
        for v, kind in KIND.items():
            if v in m and m[v].notna().any() and not (m[v].astype(str).str.strip().isin(["", "unknown", "nan"]).all()):
                ok = m[v].notna() & ~m[v].astype(str).str.strip().isin(["", "unknown", "nan"])
                n_all = m.loc[ok, "subject_id"].nunique()
                n_lab = m.loc[ok & m["subject_id"].astype(str).isin(lab_ids), "subject_id"].nunique()
                clin.append(dict(study_id=study, variable=v, kind=kind, defines_label=v in DEFINES_LABEL,
                                 source=f"data/meta/{f.name}", n_participants_recorded=n_all,
                                 n_labelled_participants_recorded=n_lab, n_labelled_participants=len(lab_ids)))
                vars_recorded.append(v)
        layer_raw = str(m["omics_layer"].dropna().iloc[0])
        rows.append(dict(study_id=study, layer=LAYER.get(layer_raw, layer_raw), layer_detail=layer_raw,
                         tissue="; ".join(sorted(m["sample_matrix"].dropna().astype(str).unique())),
                         platform="; ".join(sorted(m["platform"].dropna().astype(str).unique()))[:300] if "platform" in m else "",
                         **{k: C[study][k] for k in ("paper", "lc_definition", "label_source", "comparison_group", "time_since_infection")},
                         sample_timing="; ".join(f"{k} {v}" for k, v in ph.value_counts().items()),
                         n_LC_post_acute=len(P["LC_post"]), n_recovered_post_acute=len(P["REC_post"]),
                         n_LC_outcome_acute=len(P["LC_acute"]), n_recovered_outcome_acute=len(P["REC_acute"]),
                         n_never_infected=len(P["never"]), n_unknown=len(P["unknown"]),
                         n_samples=len(m), n_participants=m["subject_id"].nunique(),
                         clinical_variables="; ".join(vars_recorded),
                         label_flag=LABEL_FLAGS.get(study, ""),
                         doubts="; ".join(x for x in [C[study]["doubts"]] + [f"{k} set unknown: {v}" for k, v in UNKNOWN_IN_ROUND2.get(study, {}).items()] if x), source_doc=C[study]["source_doc"]))
    for study, e in EXTRA.items():
        rows.append(dict(study_id=study, layer=LAYER[e["layer"]], layer_detail=e["layer"], tissue=e["tissue"], platform=e["platform"],
                         **{k: C[study][k] for k in ("paper", "lc_definition", "label_source", "comparison_group", "time_since_infection")},
                         sample_timing="", n_LC_post_acute=e.get("LC_post", 0), n_recovered_post_acute=0, n_LC_outcome_acute=0,
                         n_recovered_outcome_acute=0, n_never_infected=e.get("never", 0), n_unknown=e.get("unknown", 0),
                         n_samples=np.nan, n_participants=e.get("LC_post", 0) + e.get("unknown", 0) + e.get("never", 0),
                         clinical_variables="", label_flag="", doubts=C[study]["doubts"], source_doc=C[study]["source_doc"]))
    A = pd.DataFrame(rows)
    A["serves_3_Detection"] = (A["n_LC_post_acute"] > 0) & (A["n_recovered_post_acute"] > 0)
    A["serves_2_EarlyDetection"] = (A["n_LC_outcome_acute"] > 0) & (A["n_recovered_outcome_acute"] > 0)
    A["serves_1_Subtyping"] = A["n_LC_post_acute"] >= MIN_PER_CLASS
    return A, pd.DataFrame(clin), parts


TASKS = {  # task: (class 1 count column, class 0 count column)
    "3_Detection": ("n_LC_post_acute", "n_recovered_post_acute"),
    "2_EarlyDetection": ("n_LC_outcome_acute", "n_recovered_outcome_acute"),
    "1_Subtyping": ("n_LC_post_acute", "n_recovered_post_acute"),
}


def anchors(A: pd.DataFrame, entries: dict) -> pd.DataFrame:
    rows = []
    feats = {}
    for task, (c1, c0) in TASKS.items():
        for layer, d in A.groupby("layer", sort=True):
            d = d.assign(size=d[c1] if task == "1_Subtyping" else d[c1] + d[c0])
            ok = d[(d[c1] >= MIN_PER_CLASS) & (d[c0] >= MIN_PER_CLASS) & (d["label_flag"] == "")].sort_values(["size", "study_id"], ascending=[False, True])
            note = ""
            if task == "1_Subtyping" and ok.empty:  # no recovered group: LC participants only, 200 most variable features
                ok = d[(d[c1] >= MIN_PER_CLASS) & (d["label_flag"] == "")].sort_values(["size", "study_id"], ascending=[False, True])
                note = "no cohort with 10 recovered: anchor chosen on LC participants; features = 200 most variable" if len(ok) else ""
            listed = d[(d[c1] > 0) | (d[c0] > 0)].sort_values("size", ascending=False)
            if listed.empty:
                continue
            anchor = ok["study_id"].iloc[0] if len(ok) else ""
            partner = ok["study_id"].iloc[1] if len(ok) > 1 else ""
            shared = np.nan
            if anchor and partner:
                for s in (anchor, partner):
                    if s not in feats and s in entries:
                        feats[s] = features(s, entries, A.set_index("study_id").at[s, "layer_detail"])
                if anchor in feats and partner in feats:
                    shared = len(feats[anchor] & feats[partner])
            def counts(s):
                r = d.set_index("study_id").loc[s]
                return f"{int(r[c1])} / {int(r[c0])}"
            rows.append(dict(kind="horizontal", task=task, layer=layer, anchor=anchor, anchor_LC_vs_recovered=counts(anchor) if anchor else "",
                             anchor_n_features=len(feats.get(anchor, ())) or np.nan, partner=partner,
                             partner_LC_vs_recovered=counts(partner) if partner else "", partner_n_features=len(feats.get(partner, ())) or np.nan,
                             n_shared_features=shared,
                             candidates="; ".join(f"{r.study_id} {int(r[c1])}/{int(r[c0])}" + (" (flagged)" if r.label_flag else "")
                                                  + ("" if (r[c1] >= MIN_PER_CLASS and r[c0] >= MIN_PER_CLASS) else " (<10 per class)")
                                                  for _, r in listed.iterrows()),
                             note=note or ("no partner: only one cohort of this layer has 10 per class" if anchor and not partner else
                                           ("no anchor: no cohort of this layer has 10 per class" if not anchor else ""))))
    return pd.DataFrame(rows), feats


def vertical(parts: dict, A: pd.DataFrame, clin: pd.DataFrame) -> pd.DataFrame:
    rows = []

    def pair(a, b, cohort, link):
        ma, pa, Pa = parts[a]
        mb, pb, Pb = parts[b]
        for task, keys in (("3_Detection", ("LC_post", "REC_post")), ("2_EarlyDetection", ("LC_acute", "REC_acute"))):
            n1 = len(Pa[keys[0]] & Pb[keys[0]])
            n0 = len(Pa[keys[1]] & Pb[keys[1]])
            if n1 or n0 or len(Pa[keys[0]]) or len(Pb[keys[0]]):
                rows.append(dict(kind="vertical", task=task, cohort=cohort, layer_1=a, layer_2=b, link=link,
                                 participants_both_LC=n1, participants_both_recovered=n0,
                                 meets_10_per_class=bool(n1 >= MIN_PER_CLASS and n0 >= MIN_PER_CLASS), status="on disk, linked by subject_id"))
    pair("Cervia2024_SomaScan", "Cervia2024_MS", "Cervia2024", "same serum sample_id and subject_id (Round1/docs/paired_cohorts.md)")
    pair("Su2022", "Su2022_metabolomics", "Su2022", "same blood draw ids (Round1/docs/paired_cohorts.md)")
    for cohort, l1, l2, status in (
        ("Fineschi2026 / Fineschi2025", "Fineschi2026 (Olink)", "Fineschi2025 PBMC RNA-seq", "RNA count matrix open (CC BY 4.0) but not on disk; 110 people overlap (60 PAT, 50 CTL); same draw not stated"),
        ("CORSIP2026", "CORSIP2026 (NanoString)", "serum metabolomics and lipidomics (repository)", "flagged cohort (gene panel only); metabolomics not on disk"),
        ("PRJNA1149196", "PRJNA1149196 (RNA)", "plasma cytokines; PRIDE PXD050202 cited", "not checked or not on disk"),
        ("GSE297499", "PBMC RNA-seq", "CSF cell RNA-seq (same series)", "CSF not built; 17 / 9 PBMC is below 10 per class"),
        ("GSE226260", "whole blood", "PBMC", "different id systems: people who gave both cannot be linked"),
        ("Kovarik2023", "plasma proteomics", "eicosanoids, Biocrates, triglycerides", "supplement column labels do not link the tables per sample (not built)")):
        rows.append(dict(kind="vertical", task="", cohort=cohort, layer_1=l1, layer_2=l2, link="", status=status))
    # clinical variables as a layer (non-label variables on labelled participants)
    for study, d in clin.groupby("study_id"):
        d = d[~d["defines_label"] & (d["n_labelled_participants_recorded"] > 0) & (d["kind"] != "timing")]
        if d.empty:
            continue
        a = A.set_index("study_id").loc[study]
        for task, (c1, c0) in (("3_Detection", ("n_LC_post_acute", "n_recovered_post_acute")),
                               ("2_EarlyDetection", ("n_LC_outcome_acute", "n_recovered_outcome_acute"))):
            if a[c1] and a[c0]:
                rows.append(dict(kind="vertical", task=task, cohort=study, layer_1=f"{study} ({a['layer']})", layer_2="clinical",
                                 link="same participants (metadata)", participants_both_LC=int(a[c1]), participants_both_recovered=int(a[c0]),
                                 meets_10_per_class=bool(a[c1] >= MIN_PER_CLASS and a[c0] >= MIN_PER_CLASS),
                                 clinical_variables="; ".join(f"{r.variable} ({r.n_labelled_participants_recorded}/{r.n_labelled_participants})" for r in d.itertuples()),
                                 status="on disk; variables that define the label (symptoms) excluded"))
    return pd.DataFrame(rows)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    A, clin, parts = audit()
    L, linked_vars = clinical_links(parts)
    L.to_csv(OUT / "clinical_links.csv", index=False)
    clin = pd.concat([clin, linked_vars], ignore_index=True)
    extra = linked_vars.groupby("study_id")["variable"].apply(lambda s: "; ".join(s + " (linked table)"))
    A["clinical_variables"] = [("; ".join(x for x in (a, extra.get(s, "")) if x)) for s, a in zip(A["study_id"], A["clinical_variables"].fillna(""))]
    entries = config_entries()
    H, feats = anchors(A, entries)
    V = vertical(parts, A, clin)
    A.to_csv(OUT / "labels_audit.csv", index=False)
    lab = []
    for study, (m, ph, P) in parts.items():
        d = m.assign(sample_phase_round2=ph.values)[["study_id", "subject_id", "sample_id", "lc_status", "sample_phase_round2"]]
        lab.append(d.assign(label_flag=LABEL_FLAGS.get(study, "")))
    pd.concat(lab, ignore_index=True).to_csv(OUT / "sample_labels.csv", index=False)
    clin.to_csv(OUT / "clinical_variables.csv", index=False)
    pd.concat([H, V], ignore_index=True).to_csv(OUT / "anchors_and_partners.csv", index=False)
    pd.DataFrame([dict(path=k, sha256=sha(p), bytes=p.stat().st_size) for k, p in sorted(INPUTS.items())]).to_csv(STEP / "inputs.csv", index=False)
    print(H[["task", "layer", "anchor", "anchor_LC_vs_recovered", "partner", "partner_LC_vs_recovered", "n_shared_features"]].to_string(index=False))
    print(V[["task", "cohort", "layer_1", "layer_2", "participants_both_LC", "participants_both_recovered"]].to_string(index=False))


if __name__ == "__main__":
    main()
