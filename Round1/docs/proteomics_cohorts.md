# Proteomics cohorts

The plasma and serum proteomics cohorts in `docs/proteomics_inventory.csv` that need no login. This file records what was
built and how. It does not interpret anything. The benchmarks run on these cohorts are described in `results/README.md`
("Proteomics").

- **Downloads:** `scripts/proteomics/download.sh`; source, licence, checksums and date are in `docs/provenance.md`
  ("Proteomics cohorts").
- **Builders:** one per study, `scripts/proteomics/build_<study>.py`. All of them write through
  `scripts/proteomics/common.py:write_cohort`:
  - `data/proteomics/<cohort>/<cohort>_matrix.csv`: UniProt key (rows) by sample (columns), one column per measured
    sample or run;
  - `data/proteomics/<cohort>/<cohort>_features.csv`: which source assays, aptamers or protein groups make each key;
  - `data/meta/<cohort>_meta.csv`: ClinicalSchema metadata, validated strictly.

  `data/` is not committed.
- **Shared feature space:** the raw values reach it only through the layer adapters of `chorus.adapters` (`olink`,
  `somascan`, `ms_proteomics`). The key is the canonical UniProt accession. Multi-protein assays are keyed by their
  sorted accessions joined with "_".
- **Processing:** the builders apply only vendor QC flags. There is no normalisation, imputation, bridging or cross-plate
  correction.
- **Tables:** `python scripts/proteomics/build_all.py` reruns every builder and the two generated tables below;
  `--tables-only` redraws the tables only.

**Schema additions** (`chorus.schema`):
- layer `somascan`;
- `control_type`: COVID_recovered, convalescent_unverified (added 15 November 2026: a control the authors call
  convalescent without a per person record of infection), never_infected_prepandemic, never_infected_contemporary,
  never_infected_unspecified, not_applicable;
- layer `proteomics` (added 15 November 2026): Olink, SomaScan and MS cohorts pooled in one benchmark config;
- `sample_phase`: pre_infection, acute, post_acute, not_applicable, unknown;
- `timepoint` (visit label).

Acute samples carry the patient's later outcome in `lc_status` when the authors give it.

## Curation decisions (15 November 2026)

Applied in the builders, so the files in `data/` and the tables below already reflect them.

| cohort | decision | effect |
|---|---|---|
| Cervia2024 | the 5 healthy controls that match later COVID-19 subjects are linked to those subjects, and their healthy samples are dropped | Cervia2024_SomaScan: 39 to 34 healthy samples, 265 to 260 samples. The COVID-19 subjects' rows (SomaScan and MS) carry `linked_healthy_subject_id` and a note. `later_covid_subject_id` is removed. |
| Chowdhury2024 | the authors' labels are kept; the controls get `control_type` convalescent_unverified; the unlabelled runs are excluded | 55 to 51 runs (`DIA_Pre_13`, `DIA_Post_13`, `DIA_C59_60a`, `DIA_Sheela` excluded); 405 proteins, none dropped with them. `lc_status` of the controls stays COVID_recovered. |
| Zoodsma2022 | the 16 Nijmegen patients not in the authors' "allowed to use" file stay excluded | unchanged |
| Captur2022 | not used | the cohort is still built (`unknown` labels), but no benchmark config lists it |

<!-- generated:start -->
### Cohorts

Samples are matrix columns (runs); subjects in parentheses. n LC and n COVID recovered count post-acute and unknown-phase samples only (acute samples are in n acute whatever their later outcome).

| cohort | layer | samples (subjects) | n LC | n COVID recovered | n never infected | n acute | n other or unknown | sample matrix | platform | proteins | other layer on the same people |
|---|---|---|---|---|---|---|---|---|---|---|---|
| Bansal2025 | olink | 56 (56) | 12 (12) | 20 (20) | 24 | 0 | 0 | plasma | Olink Target 96 Inflammation and Neurology (NPX) | 182 | not checked |
| Buggert2025 | olink | 221 (199) | 159 (146) | 58 (51) | 0 | 0 | 4 | plasma | Olink Explore 3072 (Explore 3.0.5); Olink Explore 3072 (Explore 3.10.0); Olink Explore 3072 (Explore 3.6.0) | 2903 | no (flow cytometry, serology on request only) |
| Captur2022 | ms_proteomics | 320 (157) | 0 (0) | 0 (0) | 0 | 0 | 320 | plasma | Targeted LC-MRM (custom ~91 protein inflammation/neuroinflammation panel, two MRM methods; Skyline); peptide area / yeast enolase internal standard | 100 | none public |
| Cervia2024_MS | ms_proteomics | 119 (115) | 40 (40) | 73 (73) | 0 | 0 | 6 | serum | LC-MS/MS DIA (Orbitrap Eclipse; DIA-NN library-free, maxLFQ) | 308 | yes: SomaScan (Cervia2024_SomaScan, same sera) |
| Cervia2024_SomaScan | somascan | 260 (147) | 40 (40) | 73 (73) | 34 | 113 | 0 | serum | SomaScan 7k (v4; 7596 SOMAmers deposited) | 6402 | yes: MS proteomics (Cervia2024_MS, same sera at 6 months); monocyte scRNA-seq, PBMC flow, VirScan (Mendeley 10.17632/dvf6yvrg4x.1, subsets); Olink Inflammation (Ceglarek 2025) |
| Chowdhury2024 | ms_proteomics | 51 (30) | 42 (21) | 9 (9) | 0 | 0 | 0 | plasma | LC-MS/MS diaPASEF timsTOF Pro (DIA-NN 1.8.1 library-free, MaxLFQ pg_matrix) | 405 | no (glutathione, 8-OHdG, ELISA, flow cytometry; PCC participants from the BQC19 biobank, no omics link given) |
| Chowdhury2026 | ms_proteomics | 150 (150) | 53 (53) | 62 (62) | 35 | 0 | 0 | plasma | LC-MS/MS diaPASEF timsTOF Pro (DIA-NN 1.8.1 library-free, MaxLFQ pg_matrix) | 473 | not checked |
| Fineschi2026 | olink | 165 (165) | 92 (92) | 73 (73) | 0 | 0 | 0 | plasma | Olink Explore 384 Inflammation (Explore 3.11.1) | 358 | not checked |
| Kovarik2023_plasma | ms_proteomics | 39 (39) | 13 (13) | 13 (13) | 13 | 0 | 0 | plasma | LC-MS/MS DDA-PASEF timsTOF Pro (MaxQuant 1.6.17.0 LFQ, match between runs) | 323 | eicosanoids, Biocrates metabolomics, triglycerides in the paper's supplementary tables, but their column labels do not link the tables per sample (docs/paired_cohorts.md); no repository accession |
| Su2022 | olink | 817 (493) | 75 (75) | 50 (50) | 289 | 390 | 13 | plasma | Olink Target 96 (Cardiovascular II, Inflammation, Immune Response, Metabolism, Organ Damage); authors' Z scores | 441 | yes: plasma metabolomics (Table S2.2), PBMC scRNA-seq/CITE-seq (ArrayExpress E-MTAB-10129); TCR-seq without accession |
| Woodruff2023 | olink | 149 (149) | 93 (93) | 27 (27) | 0 | 29 | 0 | plasma | Olink Explore 3072 (NPX) | 2924 | not checked |
| Zoodsma2022_Breda | olink | 576 (293) | 0 (0) | 0 (0) | 5 | 542 | 29 | plasma | Olink Target 96 Inflammation (v.3022); NPX | 92 | none public |
| Zoodsma2022_Hannover | olink | 261 (260) | 0 (0) | 0 (0) | 62 | 13 | 186 | plasma | Olink Target 96 Cardiometabolic (v.3603), Cardiovascular II (v.5006), Inflammation (v.3022), Neurology (v.8013); NPX | 359 | no omics (flow counts, serology in the annotation; Bonifacius 2021) |
| Zoodsma2022_Nijmegen | olink | 550 (144) | 0 (0) | 0 (0) | 0 | 550 | 0 | plasma | Olink Target 96 Cardiometabolic (v.3603), Cardiovascular II (v.5006), Inflammation (v.3021); NPX; Olink Target 96 Cardiometabolic (v.3603), Cardiovascular II (v.5006), Inflammation (v.3022); NPX; Olink Target 96 Cardiometabolic (v.3603), Inflammation (v.3021); NPX | 269 | no omics (flow cytometry; Janssen 2021) |

### Proteins shared between cohorts

Number of UniProt keys in common (diagonal: proteins in the cohort).

| | Bansal2025 | Buggert2025 | Captur2022 | Cervia2024_MS | Cervia2024_SomaScan | Chowdhury2024 | Chowdhury2026 | Fineschi2026 | Kovarik2023_plasma | Su2022 | Woodruff2023 | Zoodsma2022_Breda | Zoodsma2022_Hannover | Zoodsma2022_Nijmegen |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Bansal2025 | 182 | 179 | 5 | 0 | 174 | 1 | 0 | 76 | 0 | 89 | 179 | 92 | 182 | 92 |
| Buggert2025 | 179 | 2903 | 50 | 126 | 2147 | 159 | 198 | 358 | 133 | 430 | 2903 | 92 | 355 | 268 |
| Captur2022 | 5 | 50 | 100 | 26 | 74 | 31 | 40 | 6 | 29 | 9 | 51 | 5 | 11 | 11 |
| Cervia2024_MS | 0 | 126 | 26 | 308 | 177 | 218 | 224 | 5 | 187 | 4 | 127 | 0 | 23 | 23 |
| Cervia2024_SomaScan | 174 | 2147 | 74 | 177 | 6402 | 216 | 287 | 307 | 186 | 386 | 2163 | 85 | 344 | 255 |
| Chowdhury2024 | 1 | 159 | 31 | 218 | 216 | 405 | 334 | 4 | 237 | 5 | 160 | 0 | 34 | 33 |
| Chowdhury2026 | 0 | 198 | 40 | 224 | 287 | 334 | 473 | 5 | 263 | 11 | 199 | 0 | 45 | 45 |
| Fineschi2026 | 76 | 358 | 6 | 5 | 307 | 4 | 5 | 358 | 2 | 186 | 358 | 63 | 107 | 94 |
| Kovarik2023_plasma | 0 | 133 | 29 | 187 | 186 | 237 | 263 | 2 | 323 | 4 | 134 | 0 | 27 | 27 |
| Su2022 | 89 | 430 | 9 | 4 | 386 | 5 | 11 | 186 | 4 | 441 | 433 | 89 | 174 | 174 |
| Woodruff2023 | 179 | 2903 | 51 | 127 | 2163 | 160 | 199 | 358 | 134 | 433 | 2924 | 92 | 356 | 269 |
| Zoodsma2022_Breda | 92 | 92 | 5 | 0 | 85 | 0 | 0 | 63 | 0 | 89 | 92 | 92 | 92 | 92 |
| Zoodsma2022_Hannover | 182 | 355 | 11 | 23 | 344 | 34 | 45 | 107 | 27 | 174 | 356 | 92 | 359 | 269 |
| Zoodsma2022_Nijmegen | 92 | 268 | 11 | 23 | 255 | 33 | 45 | 94 | 27 | 174 | 269 | 92 | 269 | 269 |
<!-- generated:end -->

## Per cohort: definitions and what could not be recovered

Quotes are from the papers' Methods or Results. `case_definition` is the ClinicalSchema value assigned.

### Buggert2025 (Gao Y, ..., Buggert M 2025 Nat Immunol; Olink Explore 3072)
- **LC (UK, University Hospital Llandough).** "Cases were diagnosed according to the National Institute for Health and
  Care Excellence guideline NG188 ... symptoms that persisted for at least 12 weeks after the initial diagnosis of acute
  COVID-19. One persistent symptom was sufficient". All participants had mild, molecularly confirmed acute COVID-19, and
  none was hospitalised. `case_definition` NICE_NG188.
- **Recovered.** "healthy convalescent individuals (controls)". The Olink subset is 51 LC and 51 convalescent, selected
  "data-driven" from 70 and 70.
- **Sweden (Karolinska).** 95 long COVID cases only; the 30 Swedish controls were not profiled. No case definition is
  given (`study_specific`). KLIMP_KS0010 and KS0280 are on the plates but not in the deposited metadata, so they are
  `unknown`.
- **Not recovered:**
  - days since infection (the UK medians are 416 days for cases and 268 for controls, from diagnosis; none for Sweden);
  - vaccination and variant, which are reported by group only;
  - Swedish severity.
- **Values and runs:**
  - NPX, plate-control normalised, from three Olink projects with Explore software 3.0.5, 3.6.0 and 3.10.0.
  - 16 UK samples are on plates 1 and 2, and 6 Swedish samples on plates 2 and 3 (the authors' bridges). Each run is its
    own column with a shared `subject_id`, and the runs are not bridge-normalised.
  - Values whose QC_Warning or Assay_Warning is not PASS are missing (12,860).
  - NT-proBNP carries no UniProt accession in the export and is dropped.

### Zoodsma2022 (Front Immunol 2022, "Targeted proteomics identifies circulating biomarkers associated with active COVID-19 and post-COVID-19"; Olink Target 96)
Three sub-cohorts with different people, written as three cohorts from the raw NPX exports of the repository. The
authors' bridged and processed files are not used.
- **Hannover:** 186 post-COVID-19 convalescents sampled 13 to 75 days after symptom resolution, and 62 healthy blood
  donors (seronegative, recruited April to October 2020, so `never_infected_contemporary`).
- **Breda:** hospitalised COVID-19, longitudinal.
- **Nijmegen:** hospitalised COVID-19, longitudinal, in three separately normalised runs.
- **No LC or recovered split exists.** The Discussion says: "We are not able to relate the protein concentrations to
  chronic COVID-19 symptoms since the enrolled post-COVID-19 patients were largely symptom-free". So the Hannover
  convalescents are `unknown` with `post_acute` phase, and the hospitalised samples are `acute` with `unknown` outcome.
- **Excluded:** 16 Nijmegen patients that the authors' `patient_information_Aline.xlsx` does not list as usable (63
  samples, plus bridge copies in the other runs), and one control well.
- **Bridges:** samples of one sub-cohort measured in another run are kept in that run with `bridge_sample` set.
- **Not recovered:**
  - days since infection for Hannover (only days since symptom resolution, in its own column) and for Nijmegen (week
    labels only);
  - vaccination;
  - symptoms;
  - oxygen status of non-ICU patients.
- **Panels:** Inflammation, Cardiometabolic, Cardiovascular II and Neurology (Hannover); Inflammation only (Breda);
  Inflammation, Cardiometabolic and Cardiovascular II (Nijmegen). Inflammation is v.3021 in Nijmegen runs 1 and 2.

### Cervia2024_SomaScan and Cervia2024_MS (Cervia-Hasler 2024 Science)
The same serum samples measured by SomaScan (acute, 6 months and healthy) and by DIA mass spectrometry (6 months).
`sample_id` and `subject_id` are shared.
- **LC.** "Long Covid was defined as one or more persisting COVID-19–related symptoms, which could not be explained by an
  alternative diagnosis ... 6-month Long Covid ... persisting Long Covid symptoms at 6-month ... follow-up". Patients with
  only smell or taste changes (8 of 48) were excluded from LC; that is the deposited label PACS_6M_woDys.
  `case_definition` is study_specific because of that exclusion and no stated duration criterion.
- **Recovered.** PACS_6M_woDys = 0 (73 patients). This group mixes three groups that the deposited labels cannot
  separate:
  - patients recovered within one month (49);
  - "recovered Long Covid" (symptoms for more than 4 weeks, none at 6 months; 16);
  - the 8 smell- or taste-only patients.
- **Healthy.** 39 seronegative people recruited at the same time (April 2020 to April 2021), so
  `never_infected_contemporary`; 34 after curation (below).
- **Acute samples.** 113, with the 6-month outcome in `lc_status`.
- **Flags:**
  - 5 healthy controls (COV_ZH_90/181/185/197/198_ID) match later COVID-19 subjects (`..._2_ID`) by number, sex and age.
    Curation decision: they are taken to be those people before infection. Their healthy samples are dropped, and
    the later subject's rows record the link in `linked_healthy_subject_id`.
  - 6 MS samples have no SomaScan or label row and are `unknown`.
- **Values:**
  - SomaScan: SomaScan v4 with 7596 aptamers deposited (7k); the 7289 human aptamers are kept. The deposited values are
    log10 RFU, converted to RFU before the adapter's log2. No plate columns are deposited.
  - MS: DIA-NN maxLFQ with the authors' drift correction, 308 proteins.
- **Not recovered:** variant era, symptom lists, and the sub-groups of the recovered label.

### Su2022 (Su Y 2022 Cell, INCOV; Olink Target 96, authors' Z scores)
- **Values.** The deposited values are not NPX. They are "batch-corrected ... further adjusted for age, sex and BMI ...
  converted into Z-scores using the means and the standard deviations estimated for the residuals in the matched control
  samples". No plate, batch or QC flags are deposited.
- **Panels and mapping.** Five panels: Cardiovascular II, Immune Response, Inflammation, Metabolism and Organ Damage.
  Symbols are mapped to UniProt through HGNC (453 of 454 assays; CALCA is ambiguous and dropped).
- **LC.** PASC is assessed at the T3 draw (about 2 to 3 months after onset) by interview and chart review: "PASC
  symptoms were deemed to be present when reported symptoms from the interview and EHR are consistent". LC means at least
  one PASC symptom (the authors' 1–3 and ≥4 PASC groups, 75 patients). Recovered is "no PASC" (50). `case_definition` is
  study_specific (a single time point, no duration requirement).
- **Healthy.** 289 participants of the Arivale wellness programme, `never_infected_unspecified` (pre-pandemic collection
  is not stated).
- **Acute samples.** T1 (diagnosis) and T2 (about a week later) are acute and carry the patient's T3 outcome. 11 draws
  labelled PD1 or PD2 are not described in the paper and are `unknown`.
- **HAARVI.** These 100 patients are not in the deposited proteomics.
- **Not recovered:** vaccination, variant, collection dates, and healthy-control data beyond age, sex and BMI.

### Kovarik2023_plasma (Kovarik 2023 iScience; MaxQuant LFQ, timsTOF DDA-PASEF)
- **Groups.** 13 long COVID, 13 recovered and 13 healthy, recruited May to July 2021 in Vienna, one plasma run each. The
  sample codes (`LC`, `pc`, `i`) are resolved to groups by matching the MaxQuant LFQ columns to the per-sample columns of
  the paper's Table S1.
- **LC.** At least 3 months after PCR-positive infection, "chronic fatigue and/or severe chronic lack of concentrations
  ... combined with at least one more chronic symptom". `case_definition` study_specific.
- **Recovered.** At least 3 months after infection, "symptom free in anamnesis at the time of blood draw", anti-N and
  anti-S positive.
- **Healthy.** Fully vaccinated, anti-N negative (`never_infected_contemporary`).
- **Not recovered (group summaries only):** time since infection, sex, age, severity, vaccination before infection,
  variant.
- **Values.** MaxQuant LFQ intensity. 48 groups are dropped: those flagged only identified by site, reverse, or
  potential contaminant.
- **Not a cohort.** PXD036970 (macrophage supernatants) is a U937 cell-line experiment with no participant material.

### Captur2022 (Captur 2022 eBioMedicine; targeted MRM, COVIDsortium healthcare workers)
- **No labels deposited.** The deposited files carry no timepoint, infection status, symptom outcome, age or sex.
  `Supplementary_Data_1.xlsx` is pathway enrichment output only. The clinical data are "available following a data
  access request through the COVIDsortium data access portal". Every sample is therefore `unknown` in `lc_status` and
  `sample_phase`.
- **Paper definitions.**
  - Cases: PCR-confirmed infection between weeks 0 and 6, sampled serially from 1 week before to 6 weeks after the first
    positive test.
  - Uninfected: seronegative throughout and sampled at baseline only.
  - Persistent symptoms: "1 or more persistent symptom out to 12 months".
- **Values.** Per peptide, the quantifier transition area divided by the spiked yeast enolase standard (light/heavy for
  ALDOA), then the log2 mean over the protein's peptides. Two injections per sample are two different MRM panels and are
  joined into one column.
- **Proteins.** 100 proteins, against 91 in the paper.
- **Instrument.** The PRIDE record says Agilent 6410; the paper says Waters Xevo TQ-XS.
- **IDs.** 10 barcodes have conflicting participant numbers between injections (the Injection1 number is used), and one
  conflicting well is dropped.

### Chowdhury2024 (Chowdhury 2024 Front Physiol; DIA-NN)
- **Design.** An 8-week supervised exercise intervention in people with post-COVID condition (PCC). Participants were
  "allocated sequentially" (not randomised) to exercise or usual care, and plasma was drawn at recruitment and a week
  after the intervention. Recruitment ran from March 2022 to April 2023 in Sherbrooke, Canada.
- **LC (PCC).** "PCC was diagnosed when a patient presents at least one candidate symptom 3 months following a positive
  COVID-19 PCR test", with WHO-listed symptoms. PCC duration ranges from 1 month (Table 1), so `case_definition` is
  study_specific.
- **LC samples.** 21 people with runs before and after the intervention: 12 exercise, 9 no exercise. Both draws are
  `post_acute` and share `subject_id`. `timepoint` and `intervention_arm` record the visit and the arm.
- **Controls.** "No-PCC, convalescent individuals without PCC", 9 healthy hospital co-workers. `lc_status` is the
  authors' label, `COVID_recovered`, and `control_type` is `convalescent_unverified` (curation decision):
  - No infection test or date is reported for them.
  - The Limitations section asks for future "inclusion of infected individuals who do not develop PCC".
  - The recovered label rests on the Table 1 footnote only.
- **Run mapping.** Runs are matched to `metadata.sdrf.tsv` by file name.
  - Two SDRF rows (pre-covex26 and pre-covex27) point to the post files; they are mapped to `DIA_Pre_26` and
    `DIA_Pre_27`, which is flagged.
  - Runs absent from the SDRF have no label: `DIA_Pre_13` and `DIA_Post_13` (probably the 10th no-exercise
    participant, not deposited), `DIA_C59_60a` and `DIA_Sheela`. Curation decision: they are excluded from the
    written cohort.
- **Values.** DIA-NN 1.8.1 `report.pg_matrix.tsv` (MaxLFQ, library-free, UniProt UP000005640 including TrEMBL).
  - Id column `Protein.Group`.
  - DIA-NN orders a group's accessions alphabetically. The adapter therefore takes the first accession that is in the
    HGNC Swiss-Prot set when there is one (`preferred_accessions`). This leaves 23 of 405 keys outside that set, against
    24 without the preference (Kovarik2023_plasma: 3 against 5; Cervia2024_MS: 0).
  - The authors' own analysis used `unique_genes_matrix.tsv` and removed "non-mover" samples; neither is done here.
- **Not recovered:** days since infection (only group PCC duration), severity, vaccination, BMI and comorbidities (group
  level), variant era.

### Vijayakumar2022 (Vijayakumar 2022 Immunity; Olink, Dryad)
Not downloaded. Dryad refuses scripted downloads (see `docs/provenance.md`). No matrix was built.

## Cohorts added 16 November 2026

Leads of `docs/dataset_sweep_2026-11.md`. Downloads, checksums and licences are in `docs/provenance.md` ("Cohorts added
16 November 2026"). They enter the version 4 configs (`configs/lc_proteomics_affinity_v4.yaml`,
`configs/lc_proteomics_ms_v4.yaml`, `configs/concordance_proteomics_all_v4.yaml`); the earlier configs are unchanged.

### Fineschi2026 (Fineschi 2026 Front Immunol; Olink Explore 384 Inflammation)
- **Source.** SciLifeLab figshare 10.17044/scilifelab.31316779 (CC BY 4.0): `XL-4180_raw_data.csv` (long NPX export,
  two plates, Explore 3.11.1), README, MANIFEST and the Olink analysis report. Paper PMID 41808838, PMC12968220.
- **Labels.** There is no separate sample sheet. The group is the prefix of `SampleID`: PA-nn (patient) and KO-nn
  (kontroll). The export has 92 PA and 73 KO ids, the counts the README and the paper give (92 patients, 73 controls).
  In the authors' PBMC RNA-seq deposit (below) the same prefixes carry an explicit `Condition`: every PA id is PAT and
  every KO id is CTL.
- **LC.** Met the WHO post COVID-19 condition criteria (paper Methods, Participants); `case_definition` WHO_2021.
- **Controls.** "Contracted SARS-CoV-2 during a similar timeframe but recovered fully without any persistent symptoms or
  complications", matched by age and sex: `COVID_recovered`.
- **Timing.** Infection between February 2020 and June 2022, none hospitalised; blood drawn in April-May 2023 or
  November-December 2024; mean symptom duration 34 months (paper Table 1). Per person dates are not deposited, so
  `days_since_infection` is empty. Uppsala, Sweden.
- **Values.** NPX; values whose `QC_Warning` or `Assay_Warning` is not PASS are set missing (1,650 values); values below
  LOD are kept; Olink control wells dropped. 368 assays, 358 UniProt keys.
- **Same people as the PBMC RNA-seq set 10.17044/scilifelab.28832492 (Fineschi 2025 Front Immunol, PMID 40519915).**
  - That deposit holds `salmon.merged.gene_counts.tsv` (nf-core/rnaseq 3.15.1) and `SamplesPC.txt` (110 labelled
    samples: 60 PAT, 50 CTL; columns SampleName, Sample, Condition, Wave, Score, Dyspnoea), CC BY 4.0. Only the label
    file was downloaded (`data/counts/Fineschi2025_check/`).
  - All 110 labelled PBMC sample ids (KO.01 written KO-01) are Olink sample ids: 60 of the 92 patients and 50 of the
    73 controls. The 2026 paper says the PASC cohort of the earlier transcriptome study is "also included in this
    study". The PBMC blood draws were in April-May 2023, the first of the two Olink sampling periods; whether the
    plasma of the overlapping people is from the same draw is not stated.
  - `overlaps_with_study` is set for the 110 people. The PBMC set is not built here.

### Chowdhury2026 (Chowdhury 2026 Sci Rep; DIA-NN)
- **Source.** PRIDE PXD066724 (released 2026-03-24): `diann.7z` (DIA-NN 1.8.1 output; only the summary tables are
  extracted), `metadata.sdrf.tsv` and `PXD066724_community_annotated.sdrf.tsv`. Paper PMID 41998026, PMC13260905.
- **People.** BQC19 biobank, Université de Sherbrooke hospital research centre; samples March 2020 to October 2021;
  plasma 3 months after a PCR confirmed infection (paper Table 1: no infection 35, convalescent 62, PCC 53). One SDRF
  row and one run per person.
- **LC.** PCC: "persistent symptoms consistent with WHO criteria in effect during the sample collection period";
  `case_definition` WHO_2021.
- **Controls.** Convalescent: PCR negative 1 month after infection and "free of infection-related symptoms at 3 months
  post-infection": `COVID_recovered`. Uninfected: RT-PCR negative at sampling: `never_infected`,
  `never_infected_contemporary` (not in the LC versus recovered contrast).
- **Run mapping.** The DIA-NN search has 173 runs; the 150 the SDRF names are used. Not used: 20 reruns of first batch
  samples (`sample_NNNb`), one pool (`mix_of_8`), and the `_NE_` runs of samples 213 and 215 (the SDRF names their
  `_Nc_` runs).
- **Label conflict.** Second batch file names carry a group code (Lc, c, cR, N, NE, NT, Nc). For samples 226, 227,
  231, 232 and 233 the code is Lc while both SDRFs give C. The SDRF is used: with it the totals equal the paper's,
  with the file codes they would not (58 LC, 57 convalescent). The code is kept in `file_group_code` and the 5 runs are
  named in `notes`.
- **Batches.** `ms_batch`: first (`Sheela_DIA_sample_NNN`; 43 LC, 47 convalescent, 10 uninfected) and second
  (`..._DIA_TRAP_400nl_flow`; 10 LC, 15 convalescent, 25 uninfected).
- **Values.** As Chowdhury2024: `report.pg_matrix.tsv`, `Protein.Group` through the ms_proteomics adapter with the HGNC
  Swiss-Prot accessions preferred (11 of 473 keys outside that set), zeros missing, log2, no normalisation.
  573 protein groups, 473 keys. The paper reports 235 proteins after its own filtering, which is not applied here.
- **Overlap with Chowdhury2024.** Same laboratory and biobank. Chowdhury2024 recruited March 2022 to April 2023
  (exercise study) and Chowdhury2026 sampled March 2020 to October 2021; ids are not linkable, so overlap cannot be
  excluded. `overlaps_with_study` records this.
- **Not recovered:** age (SDRF "not available"), severity, vaccination, per person dates.

### Huapaya2024 (NHLBI figshare 10.25444/nhlbi.26383798; SomaScan): checked, not added
- `V2 SOMA COVID_project_data.xlsx` (public domain) has sheets Acute_phase (110 rows), Acute_recovery_phases (64),
  Convalescent_phase (34), Omicron_data (40) and Protein_list.
- The non-protein columns are ID, cluster, age, sex, race, variant, vaccination, ordinal scale, days since onset,
  vaccination group, time or phase and sample date. No sheet has a PASC, symptom or lung disease column. The paper's
  18 PASC at 10 months therefore cannot be assigned per sample, and the cohort is not built.

## Cohorts added 17 November 2026

Downloads, checksums and licences are in `docs/provenance.md` ("Cohorts added 17 November 2026"). Both cohorts are
built (`scripts/proteomics/build_Woodruff2023.py`, `scripts/proteomics/build_Bansal2025.py`). No config uses them yet,
and no benchmark was rerun.

### Woodruff2023 (Woodruff 2023 Nat Commun; Olink Explore 3072)
- **Source.** Zenodo 10.5281/zenodo.8092298 (CC BY 4.0), the analysis code release `Covid-v0.4.0.tar.gz`. Only
  `data/complete.csv` is taken out of the archive. It is a long table with one row per sample and assay: `sample`,
  `group_hr`, `olinkID`, `uniprot`, `gene`, `NPX`, `group`, `sex`, `race`, `disease_phase`,
  `initial_disease_severity` and `dpso` (days post symptom onset). Paper PMID 37452024, PMC10349085.
- **Counts before building.** 158 samples × 2,943 assays:

  | `group_hr` | `disease_phase` | samples |
  |---|---|---|
  | PASC | recovery | 93 |
  | CR | recovery | 27 |
  | ICU Pre-Dex | acute | 16 |
  | ICU Post-Dex | acute | 13 |
  | SLE | sle | 9 (dropped) |

  Each sample id occurs once and no subject repeats across samples.
- **Labels.** `group_hr` per sample:
  - PASC = LC, `case_definition` study_specific. These are Emory PASC clinic patients with documented infection or
    meeting the CDC clinical case definition, who had new or worsening symptoms more than 14 days after onset (paper
    Methods).
  - CR = `COVID_recovered`: donors with uncomplicated recoveries, sampled at similar intervals after onset.
  - ICU Pre-Dex and ICU Post-Dex: `sample_phase` acute and `lc_status` unknown, since no later outcome is given. They
    are outside the LC contrast; the column `dexamethasone` (before or after) records the group.
- **CR count.** The paper gives 26 CR donors; the table has 27. The deposit's `data/README.md` notes that sample
  1433-COV-3140 was listed in both PASC and CR in the Olink table and as CR in the metadata. `complete.csv` lists it
  once, as CR, and it is kept as CR.
- **Timing.** `days_since_infection` = `dpso`, `time_anchor` symptom_onset:
  - PASC: median 111 days; 2 below 28 days and 38 below 90.
  - CR: median 99 days; 3 below 28 days and 13 below 90.

  Every PASC and CR sample has `sample_phase` post_acute, as the deposit's `disease_phase` (recovery) gives.
- **Other fields.**
  - Severity from `initial_disease_severity`: out = ambulatory, icu = icu_or_ventilated, hosp = hospitalised with the
    band unknown (oxygen not recorded).
  - Sex is recorded for 2 of the 93 PASC samples and for every CR and ICU sample.
- **Values.** NPX as deposited. The table has no QC column, so nothing is masked. 2,943 assays give 2,924 UniProt keys
  (assays of one protein on several panels are averaged).

### Bansal2025 (Bansal 2026 Commun Med; Olink Target 96 Inflammation and Neurology)
- **Source.** Zenodo 10.5281/zenodo.15237087 (latest version 15237088, CC BY 4.0), `baseline57.xlsx`: the
  cross-sectional baseline set, one row per participant and assay with the participant's metadata. The longitudinal
  workbook is not used. Paper PMID 41975235, PMC13442878 (COVID PROFILE cohort, Victoria, Australia).
- **Group columns and counts before building** (57 participants × 182 assays, `Time` Pre-vax for every row):

  | `study_group` | `Long_COVID status` | participants |
  |---|---|---|
  | COVID-19 recovered participant | Long COVID | 12 |
  | COVID-19 recovered participant | Recovered | 21 |
  | Non-COVID-19 healthy participant | Healthy | 24 |

  `disease_status_subgroup`: Moderate/Mild 30, Severe 2, Asymptomatic 1, empty for the healthy participants.
- **Labels.** `Long_COVID status` per participant:
  - Long COVID = LC, WHO_2021: symptoms persisting or new 3 months after infection and lasting at least 2 months
    (paper Methods).
  - Recovered = `COVID_recovered`: symptoms shorter than 90 days.
  - Healthy = `never_infected` (`never_infected_contemporary`; unvaccinated and SARS-CoV-2 naive, sampled at
    enrolment).
  - The COVID-19 samples were drawn 6 to 9 months after the first infection, before vaccination.
- **Values.** NPX. Values whose `QC_Warning` is not Pass are set missing (1,370 of 10,374), as for Buggert2025 and
  Fineschi2026. The warning is set per participant and panel:
  - 13 participants lose one panel (90 or 92 values).
  - Participant 3 (Recovered) has the warning on both panels, so all its values are missing and it is not in the
    matrix.
  - Built: 56 participants, 12 LC, 20 COVID_recovered and 24 never_infected; 182 assays and 182 UniProt keys.
  - The `Normalization` column is not used.

### Checked, not added (17 November 2026)
- **PRIDE PXD051740** (SWATH-MS of plasma supernatant and pellet, TripleTOF 5600; J Infect Public Health 2024,
  doi:10.1016/j.jiph.2024.102571).
  - The project lists 265 raw files, named by a number and a fraction (`45_P.wiff`, `77_S.wiff`), and one spectral
    library (`post_covid_swath_library__FDR.xlsx`). It has no SDRF or sample table.
  - The description says 65 people with post-COVID condition and 9 healthy controls. There is no recovered group, and
    no per-sample label is deposited.
- **PRIDE PXD041762** (MRM of 269 peptides, Agilent 6490; Mohammed 2024 J Proteome Res, PMID 38048423).
  - 74 hospitalised patients were sampled at admission and at 3 and 6 months.
  - The project has 195 mzML files named by plate well (`REC-2349_H4.mzML`) and two Skyline documents. It has no SDRF.
  - The Skyline document's replicates are well positions (A1, B1, …), and its first 3 MB hold no annotation
    definitions. No per-sample LC label is deposited.
- **Berezhnoy 2023 Front Immunol** (PMID 37228606, PMC10203989, NMR metabolomics and cytokines). The supplement
  (Europe PMC) has two files:
  - `Table_1.xlsx`: Supplementary Tables 1 to 10, which are group sizes, feature lists, test results and group means
    with SDs or fold changes (Acute, HC, LTCS, Recov).
  - `DataSheet_1.pdf`: seven pages of supplementary figures.

  Neither has a per-sample table.
