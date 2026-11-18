# Author count matrices: verdicts

Downloaded on 2026-11-09 from the GEO supplementary folders into `data/counts/<GSE>/` (git-ignored). The logCPM and
voom weight matrices of GSE215865 were not downloaded because they are derived from the raw count file. Per-GSM sample
tables (title, source, characteristics, BioSample, SRX) were written to `data/meta/<GSE>_geo_samples.csv` and used
for the id checks. No sequencing reads were downloaded.

## One line per file

| file | verdict |
|---|---|
| GSE215865_rnaseq_raw_count_matrix.csv.gz | integer counts; genes x samples; Ensembl versioned (ENSG…\.5); 58,929 genes x 1,392 samples; ids are not GSM titles but match `blood_sample_id` + `_` + `library_prep_plate` 1,392/1,392 (e.g. `Subj_16ae238fT0_Plate_2`), unique, so the mapping to GSM is exact. |
| GSE226260_AdditionalSamples.rawCounts.csv.gz | integer counts; genes x samples; Ensembl unversioned; 78,932 genes x 103 samples; ids equal GSM titles 103/103, all of them the whole blood subset. |
| GSE226260_AdditionalSamples.logNormalizedCounts.csv.gz | normalised (log scale, negative values); genes x samples; Ensembl; 78,932 x 103; ids equal GSM titles 103/103 (whole blood). Not needed: use the raw counts. |
| GSE226260_combined.ExprData.updated.csv.gz | integer counts; genes x samples; Ensembl unversioned; 24,203 genes x 228 samples (86 PBMC + 142 'Blood' acute); ids equal GSM titles 228/228. **Malformed tail**: 17,345 extra rows with an empty gene id, each holding a single value in column `p22179-s049_601096-14-Dec-2020`; dropped. Treat sample s049 (GSM of that title) with suspicion until checked. |
| GSE267625_P4O2_COVID_blood_transcriptome.csv.gz | not integer (24.7% of values fractional, e.g. 4.008; non-negative, column sums ~31 M): estimated counts (Salmon/RSEM style), not normalised; genes x samples; Ensembl unversioned; 34,515 genes x 111 samples; ids are `<subject>_<index>` (e.g. `AA002_1`) and match no GSM title, but the subject sequence of the 111 columns is identical to the subject sequence of the GSM titles in GSM order, so a positional mapping is plausible. 25 subjects have two samples (3 and 12 months), so the timepoint can only be recovered by position. Prefer re-quantification (reads are public). |
| GSE251849_Counts.txt.gz | not integer (10.2% fractional, non-negative, column sums ~15 M): estimated counts; genes x samples; Ensembl unversioned; 62,710 genes x 23 samples; ids `Control001-007`, `Convalescent001-005`, `LongCOVID001-011` match no GSM title. Groups map by count to Healthy 7 / COVID convalescent 5 / Long COVID 6 + brain fog 5, but which LongCOVID columns are brain fog is unknown. LC vs recovered labels are unaffected by that. |
| GSE270045_LC_counts.tsv.gz | not integer (48.5% fractional, column sums ~3 M): normalised or estimated, not raw counts; genes x samples; **gene symbols** (incl. `5S_rRNA`); 28,889 genes x 36 samples; ids match the `sample id` characteristic 36/36 (not the GSM titles). Re-quantify from SRA. |
| GSE169687_RAW.tar (152 featureCounts files) | integer counts; one file per GSM, merged to genes x samples `GSE169687_merged_counts.csv.gz`; Ensembl unversioned (GRCh38.101, HISAT2 + featureCounts, stranded); 60,671 genes x 152 samples; file names carry the GSM accession, 152/152 matched. |
| GSE222253_RAW.tar (35 count files) | integer counts (HTSeq style, five `__` summary rows dropped); one file per GSM, merged to `GSE222253_merged_counts.csv.gz`; Ensembl unversioned; 60,664 genes x 35 samples; file names carry the GSM accession, 35/35 matched. |

## GSE226260 split

GEO puts three cohorts into this series, distinguished by `source_name` and the `tissue` and `cell type`
characteristics. They were written as separate matrices (genes x GSM, integer counts, Ensembl):

| subset | rule | file | genes x samples | labels in GEO |
|---|---|---|---|---|
| whole blood | tissue = Whole blood | `GSE226260_wholeblood_counts.csv.gz` | 78,932 x 103 | current long_covid_status: 52 PASC / 51 Recovered; time of blood draw 3–6 M, 6–12 M, >12 M, 4 Acute |
| PBMC | cell type = PBMC | `GSE226260_PBMC_counts.csv.gz` | 24,203 x 86 | pasc status: 72 PASC / 14 NOPASC |
| acute blood | source = Blood (no cell type) | `GSE226260_blood_acute_counts.csv.gz` | 24,203 x 142 | pasc severity only (acute), no LC field |

The whole blood and PBMC subsets share 23,864 genes. The whole blood subset comes from a different count file
with a different annotation, so the subsets must be harmonised (or re-quantified together) before they are compared.

## PRJNA1184005 (Missailidis 2024, Front Immunol 15:1450853, PMID 39691709, PMC11649547)

- 14 BioSamples (SAMN44656608–SAMN44656621) and 14 SRA runs, all PBMC, NovaSeq 6000, 13 paired-end plus 1
  single-end (vacSG112). Attributes are age, sex, tissue, collection date, geo_loc_name and isolate (`SG22`,
  `vacSG21` … `vacSG117`); written to `data/meta/PRJNA1184005_biosamples.csv`. **There is no LC or recovered attribute.**
- No processed matrix: the BioProject is raw reads only and has no GEO series. The paper's data statement points
  only to PRJNA1184005, and the article has no supplementary files. Table 1 gives only group summaries (7 LC aged
  23–53, 7 recovered aged 23–63). The paper's tables are a group summary and a 70-gene DE list.
- Consequence: re-quantify from SRA (`requantify`). Labels must come from the authors (see `data/meta/label_sources.md`).
- **Built 18 November 2026** (`scripts/build_rna_cohorts_2026_11.py PRJNA1184005`).
  - **Labels from the author** (S. Annesley, personal communication, 18 November 2026):
    - LC: SG22, vacSG60, vacSG116, vacSG117, vacSG47, vacSG55 and vacSG61.
    - COVID_recovered ("Post COVID"): vacSG102, vacSG105, vacSG111, vacSG21, vacSG34, vacSG70 and vacSG112.
    - The labels agree with the deduction from age: vacSG34, vacSG102 and vacSG112 are recovered.
  - **Control definition pending confirmation.** The paper describes the recovered group as age and sex matched and
    recovered. Until the author confirms the definition, the metadata carries `control_definition_status`
    `pending_confirmation` for the 7 recovered samples, and `label_source` names the communication.
  - **Case definition.** LC is `WHO_2021`: the paper's Methods use the WHO description.
  - **Metadata.** `sample_phase` is post_acute, `sample_matrix` PBMC and `timepoint` single_visit. Age and sex come
    from the BioSample attributes (`data/meta/PRJNA1184005_biosamples.csv`), and the meta also has `collection_date`
    and `biosample`. The output is `data/meta/PRJNA1184005_meta.csv`.
  - **Groups.**
    - LC: 7 samples, mean age 41.4 (25 to 53), 85.7% female.
    - Recovered: 7 samples, mean age 45.4 (23 to 63), 71.4% female.
    - Paper Table 1: LC mean 41 (23 to 53), recovered mean 45 (23 to 63), 85.7% and 71.4% female. **The lower end of
      the LC age range differs:** 25 from the BioSample ages with the author labels, 23 in Table 1. The youngest
      sample, vacSG21 (23), is labelled recovered.
  - **Matrix.** Salmon on all reads, with the GSE224615 index and options, aggregated to genes:
    `data/quant/PRJNA1184005/salmon.merged.gene_counts.tsv`, 78,277 genes x 14 samples. Matrix columns are the
    isolates, which equal the ENA `sample_alias`. Procedure, QC and checksums are in `docs/provenance.md`.
  - **QC.** Mapping rate median 86.2% (79.0 to 90.0). Applying the GSE224615 rule within this aggregation flags one
    sample: vacSG21, 79.0%, more than 2.5 MAD below the median, against a lower bound of 80.0%. This is a flag, not an
    exclusion: it is recorded in `salmon_qc_flag` and `salmon_mapping_rate`.
  - **Configs.** The entries in the transcriptomic configs are unchanged (`train_only: true`, `enabled: false`), and
    no benchmark has been rerun with this cohort.

## Cohorts added 16 November 2026

Builder: `scripts/build_rna_cohorts_2026_11.py`. Provenance and checksums: `docs/provenance.md` ("Cohorts added 16
November 2026"). Configs: `configs/lc_transcriptomics_blood_v4.yaml` and its sensitivity variant
`configs/lc_transcriptomics_blood_v4_noQCflag.yaml`.

### PRJNA1149196 (Sinclair 2024, Nat Microbiol, PMID 39478108, PMC11602718)

- **Data:** 101 BioSamples and 101 paired-end runs of whole-blood total RNA (1.11 TB of FASTQ at ENA, 16.2 to 90.7 M
  read pairs per run, median 49.2 M). No processed matrix is deposited.
- **Two studies by site.** BioSample `biomaterial_provider` "Kirsty Short" marks cohort 2; all others are cohort 1.
  Site, library kit, sequencer and sampling time differ, so they are written as two studies:

| study | site | library, sequencer | samples (subjects) | LC | COVID_recovered | never_infected |
|---|---|---|---|---|---|---|
| PRJNA1149196_SA | South Australia (Central Adelaide Health Network) | Nugen Universal Plus total RNA (rRNA and globin depleted), DNBSEQ-T7 | 92 (54) | 55 (28 at 44 weeks, 27 at 68 weeks) | 22 (11 and 11) | 15 |
| PRJNA1149196_QLD | Queensland (COVID OZGenetics) | Illumina stranded total RNA (Ribo-Zero plus), NovaSeq 6000 | 9 (9) | 5 | 4 | 0 |

- **Labels:** BioSample `disease`.
  - `PASC` (SA 48 samples): persistent symptoms more than 12 weeks after infection, with ongoing discomfort at
    sampling.
  - `PASC-CVS` (SA 7, QLD 5): chest pain or tightness more than 12 weeks after infection.
  - Both PASC groups are LC, `case_definition` study_specific.
  - `Recovered` (SA 22, QLD 4): infected, no persistent symptoms more than 12 weeks after infection; COVID_recovered,
    control_type COVID_recovered.
  - `Control` (SA 15): uninfected (recruited before community transmission or negative spike serology);
    never_infected, control_type never_infected_contemporary. Not in the LC versus COVID_recovered contrast.
- **Time:** BioSample `Timepoint`. SA: "44 weeks post infection" (`days_since_infection` 308) and "68 weeks post
  infection" (476), nominal visit times anchored at the positive test (March to April 2020 infections). QLD: "1 year
  post infection" (the paper says more than 520 days); `days_since_infection` is left empty.
- **Subjects:** BioSample `isolate` with hyphens removed (`SA_LC-068` and `SA_LC068` are both written). 38 SA subjects
  have samples at both visits; with `one_sample_per_subject` the earliest, the 44 week sample, is kept.
- **Sex:** BioSample `sex` (SA 49 male, 43 female); QLD "missing", written as unknown.
- **Overlap:** cohort 1 is the same South Australian participants as GSE169687 (Ryan 2022, 12 to 24 weeks), which is
  disabled for lack of labels. The ids of the two deposits cannot be linked. `overlaps_with_study` records this.
- **Counts:** salmon 1.10.3 against the GENCODE v47 index of GSE224615.
  - **SA:** on the **first 5,000,000 read pairs of each run** (`scripts/salmon_quant_ena.py --max-pairs 5000000`). The
    full reads of the 101 runs (1.11 TB) were not downloaded.
  - **QLD:** on **all reads** of its 9 runs (`--runs`, 27.4 GB, `data/quant/PRJNA1149196_QLD_full`). Its FASTQ files are
    ordered by alignment position, so the 5 M pair prefix held only chr1 and chr10 reads (75.8% to 95.7% of the
    counts on chr1); the prefix matrix is not used.
  - Procedure and QC: `docs/provenance.md`. Matrix columns are the ENA `sample_alias`. The site matrices are
    `data/quant/PRJNA1149196/PRJNA1149196_{SA,QLD}_gene_counts.tsv.gz`.
- **Deviation from GSE224615:** GSE224615 was quantified on all reads; SA uses only the first 5 M pairs, which is fewer
  reads per sample than the full runs (16 M to 91 M pairs).
- **QC:** all runs quantified (three SA runs after one retry). Mapping rate median 39.56% in SA (34.72 to 47.05,
  library type ISF, 5 M pairs) and 59.37% in QLD (50.81 to 67.33, ISR, all reads). With the GSE224615 rule (below 40%,
  or more than 2.5 MAD below the median of the aggregation) 56 SA samples and no QLD sample are flagged; none is
  excluded (`docs/provenance.md`).

### GSE297499 (Orlinick 2026, J Infect Dis, PMID 40920586, PMC12811860)

- **Data:** PBMC bulk RNA-seq of women from the Yale COVID Mind Study (Illumina NovaSeq X, polyA). Only the PBMC samples
  are used; the CSF samples of the series are not.
- **Values:** the authors' TPM (HISAT2, StringTie; `GSE297499_final_processed_tpm_values_pbmc.txt.gz`), keyed by gene
  symbol, mapped to Ensembl (36,599 of 56,609 symbols; rows of one gene summed), 36,499 genes x 26 samples. Not counts:
  the config sets `is_counts: false`, `provenance: author_normalised`, and the loader applies log2(x + 1).
- **Labels:** GEO characteristic `treatment`. `NP-Long COVID` (17): self-reported neuropsychiatric Long COVID meeting
  the WHO definition, with documented infection; LC, `case_definition` WHO_2021. `Control` (9): documented SARS-CoV-2
  infection without persistent symptoms; COVID_recovered, control_type COVID_recovered.
- **Join:** the TPM columns equal the GEO sample titles with underscores removed (`2057_V2_PBMC` is column
  `2057V2_PBMC`); all 26 columns matched. One sample per subject. Sex female for all (study design).
- **Time:** no per-sample date or days since infection in GEO; `timepoint` single_visit.

## GSE228320 (Molinero 2022, Biomed Pharmacother 155:113617, PMID 36058144, PMC9424524; added 17 November 2026)

- **Data.** Whole blood RNA-seq (Tempus tubes) of 50 survivors of COVID-19 ARDS (CIBERESUCICOVID, Spain), drawn at the
  pulmonary evaluation 3 months after hospital discharge. `GSE228320_Raw_counts.txt.gz` holds integer counts
  (Ensembl gene ids, 60,683 genes × 50 samples, column sums 4.9 M to 12.4 M); its columns are the GSM descriptions
  (`MC1_n`). The normalised count file of the series was not downloaded.
- **Counts before building.** GEO characteristic `dlco (%_of_predicted)`: 18 samples below 60 (24.4 to 59.0) and 32
  at 60 or above (60.0 to 92.5). The sample titles say sequela or control, and every title agrees with its DLCO value
  (the builder stops otherwise). Sex: 35 male, 15 female.
- **Labels.** The case is a lung function outcome, not a Long COVID definition:
  - DLCO below 60% of predicted: `lc_status` LC with `case_definition` pulmonary_sequelae (a vocabulary term added for
    this cohort).
  - DLCO 60% or above: `COVID_recovered`. These are ARDS survivors with mild or no diffusion impairment; other symptoms
    are not recorded.
  - `dlco_pct_predicted` keeps the value.
  - Every sample has `acute_severity` icu_or_ventilated, `time_anchor` discharge and `timepoint`
    3_months_after_discharge.
- **Role.** A sensitivity cohort. It is in no benchmark config and must not enter the main benchmark.
- **Builder.** `python scripts/build_rna_cohorts_2026_11.py GSE228320` writes `data/meta/GSE228320_meta.csv`; the
  matrix is read as deposited.

## CORSIP2026 (Espin 2026, Front Cell Infect Microbiol, PMID 42328163, PMC13275656; added 17 November 2026)

- **Source.** GitHub repository estefespin/CORSIP_omics (CC0 1.0), cloned at commit
  8bbb8421de7c60afff328d1675b978f00f73914b into `data/counts/CORSIP2026/raw/CORSIP_omics` (104 files, 81 MB). Its files:
  - NanoString: `RawNanoStringData.csv` (raw counts), `NormNoRepNoHK.csv` (nSolver normalised, replicates and
    housekeeping genes removed) and `PanCancer.xlsx` (the panel gene list).
  - Metadata: `SubjectsNoRep.csv` (SampleID, StudyID, SexBirth, Age, Race, LCStatus), `CORSIP_Full.csv` (47
    participants × 336 questionnaire, test and serology columns) and `Joined_BMI.csv`.
  - Metabolomics: `Metabolomics_CORSIP.csv` and the HILIC, RP and lipid files, positive and negative mode.
  - Analysis output: differential expression tables, validation tables, figures, R scripts and `.RData`. None of these
    was run or loaded.
- **Platform.** The expression data are **NanoString nCounter PanCancer Immune Profiling** counts of whole blood
  (PAXgene), a targeted panel, not RNA-seq. There are no reads and no RNA-seq matrix. The raw file has 41 samples × 804
  probes: 750 endogenous genes, 40 housekeeping genes, 6 positive and 8 negative controls. The endogenous probes include
  SARS-CoV-2 M, N and S, CMV UL83 and EBV LMP2.
- **Counts before building.**
  - `SubjectsNoRep.csv`: 37 samples, `LCStatus` Case 8 and Control 29, one per participant.
  - `CORSIP_Full.csv`: 47 participants, Case 12 and Control 35 (the paper's 12 and 35). 37 of them have a NanoString
    sample, and their `LCStatus`, sex and age agree with `SubjectsNoRep.csv`.
  - The 4 other raw samples (PX-02, PX-19, PX-30, PX-43) are not in either metadata file. They correlate at 0.98 with
    retained samples, fitting the README's "replicate samples removed". They are not used.
- **Definitions** (paper Methods): CORSIP British Columbia paramedics, not hospitalised.
  - LC (Case): WHO clinical case definition from the 24-month questionnaire: new symptoms lasting at least 2 months,
    starting within 3 months of a confirmed infection. `case_definition` WHO_2021.
  - Control: confirmed SARS-CoV-2 infection (positive PCR or antigen test, nucleocapsid antibodies, or a symptom-based
    diagnosis) and no LC symptoms during the study: `COVID_recovered`.
  - Infection evidence per control in `CORSIP_Full.csv`: 27 of the 29 have a positive test or a reactive Roche
    nucleocapsid result at some collection.
  - PX-34 and PX-47 have neither, and no symptom date, and answered "I have never had COVID-19". They keep the authors'
    label but get `control_type` convalescent_unverified (as Chowdhury2024).
- **Time.**
  - PAXgene blood was drawn only at the final study visit (January 2023, paper Methods). The repository has no
    per-sample PAXgene date, so `days_since_infection` is empty.
  - The first positive test date is kept in `first_positive_test`, for every LC sample and 21 of the 29 controls.
  - `CORSIP_Full.csv` also has `MonthsSinceInfection1` to `6`, which are months between infection and the serum
    collections.
- **Other layers of the same people.** Untargeted serum metabolomics before and after infection, for 12 participants
  (6 LC, 6 control; subjects S1 to S12, timepoints 0 and 1). The paper says two LC cases are in both assays. The
  repository does not link the S ids to the PX sample or participant ids, so the layers cannot be paired. The
  metabolomics is not built.
- **Matrix and identifiers.** `scripts/build_rna_cohorts_2026_11.py CORSIP2026`:
  - The raw counts of the 790 endogenous and housekeeping genes are mapped from symbol to Ensembl with the HGNC table,
    as for GSE297499. 779 map, onto 778 Ensembl genes.
  - Not mapped: BAGE, TARP, the four KIR group probes and the five viral probes.
  - Output: `data/counts/CORSIP2026/CORSIP2026_nanostring_counts_ensembl.tsv.gz`, 778 genes × 37 samples, integer.
  - The loader treats it as counts (log CPM over the panel). All 778 genes pass the expression filter.
- **QC** (flags, not exclusions; nSolver defaults; `data/counts/CORSIP2026/CORSIP2026_nanostring_qc.tsv`):
  - Rules: positive control scaling factor 0.3 to 3, housekeeping scaling factor 0.1 to 10, positive control linearity
    R² at least 0.95, POS_E above the negative control mean plus 2 SD.
  - Results: positive control factors 0.81 to 1.31, housekeeping factors 0.69 to 1.33, R² 0.974 to 0.993, POS_E 81 to
    172 against 16 to 32. No sample is flagged.
  - The RIN ≥ 7 rule of the paper was applied by the authors before profiling.
- **Gene space.** 495 of the 778 genes are among the 10,053 genes shared by the `lc_transcriptomics_blood_v4` cohorts.
  A config that intersects genes across cohorts and includes CORSIP2026 therefore keeps at most 495 genes. How it
  enters the version 5 rerun is a config decision; no config was changed and no benchmark was rerun.
- **Participant overlap.** None of the cohorts in the benchmark recruited from CORSIP or from British Columbia
  paramedics:
  - The Canadian cohorts are Chowdhury2024 and Chowdhury2026 (Sherbrooke, Quebec, BQC19).
  - The repository's `Validation/` folder uses GSE224615 and GSE275334 as external validation sets, not as the same
    people.
  - Ids cannot be linked across deposits, so the check rests on recruitment site and study.
