# Proteomics cohorts

The plasma and serum proteomics cohorts in `docs/proteomics_inventory.csv` that need no login. This file records what was
built and how. It does not interpret anything, and no benchmark has been run on these cohorts.

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
- `control_type`: COVID_recovered, never_infected_prepandemic, never_infected_contemporary, never_infected_unspecified,
  not_applicable;
- `sample_phase`: pre_infection, acute, post_acute, not_applicable, unknown;
- `timepoint` (visit label).

Acute samples carry the patient's later outcome in `lc_status` when the authors give it.

<!-- generated:start -->
### Cohorts

Samples are matrix columns (runs); subjects in parentheses. n LC and n COVID recovered count post-acute and unknown-phase samples only (acute samples are in n acute whatever their later outcome).

| cohort | layer | samples (subjects) | n LC | n COVID recovered | n never infected | n acute | n other or unknown | sample matrix | platform | proteins | other layer on the same people |
|---|---|---|---|---|---|---|---|---|---|---|---|
| Buggert2025 | olink | 221 (199) | 159 (146) | 58 (51) | 0 | 0 | 4 | plasma | Olink Explore 3072 (Explore 3.0.5); Olink Explore 3072 (Explore 3.10.0); Olink Explore 3072 (Explore 3.6.0) | 2903 | no (flow cytometry, serology on request only) |
| Captur2022 | ms_proteomics | 320 (157) | 0 (0) | 0 (0) | 0 | 0 | 320 | plasma | Targeted LC-MRM (custom ~91 protein inflammation/neuroinflammation panel, two MRM methods; Skyline); peptide area / yeast enolase internal standard | 100 | none public |
| Cervia2024_MS | ms_proteomics | 119 (115) | 40 (40) | 73 (73) | 0 | 0 | 6 | serum | LC-MS/MS DIA (Orbitrap Eclipse; DIA-NN library-free, maxLFQ) | 308 | yes: SomaScan (Cervia2024_SomaScan, same sera) |
| Cervia2024_SomaScan | somascan | 265 (152) | 40 (40) | 73 (73) | 39 | 113 | 0 | serum | SomaScan 7k (v4; 7596 SOMAmers deposited) | 6402 | yes: MS proteomics (Cervia2024_MS, same sera at 6 months); monocyte scRNA-seq, PBMC flow, VirScan (Mendeley 10.17632/dvf6yvrg4x.1, subsets); Olink Inflammation (Ceglarek 2025) |
| Chowdhury2024 | ms_proteomics | 55 (33) | 42 (21) | 9 (9) | 0 | 0 | 4 | plasma | LC-MS/MS diaPASEF timsTOF Pro (DIA-NN 1.8.1 library-free, MaxLFQ pg_matrix) | 405 | no (glutathione, 8-OHdG, ELISA, flow cytometry; PCC participants from the BQC19 biobank, no omics link given) |
| Kovarik2023_plasma | ms_proteomics | 39 (39) | 13 (13) | 13 (13) | 13 | 0 | 0 | plasma | LC-MS/MS DDA-PASEF timsTOF Pro (MaxQuant 1.6.17.0 LFQ, match between runs) | 323 | yes, per sample in the paper's supplementary tables: eicosanoids, Biocrates metabolomics, triglycerides; no repository accession |
| Su2022 | olink | 817 (493) | 75 (75) | 50 (50) | 289 | 390 | 13 | plasma | Olink Target 96 (Cardiovascular II, Inflammation, Immune Response, Metabolism, Organ Damage); authors' Z scores | 441 | yes: plasma metabolomics (Table S2.2), PBMC scRNA-seq/CITE-seq (ArrayExpress E-MTAB-10129); TCR-seq without accession |
| Zoodsma2022_Breda | olink | 576 (293) | 0 (0) | 0 (0) | 5 | 542 | 29 | plasma | Olink Target 96 Inflammation (v.3022); NPX | 92 | none public |
| Zoodsma2022_Hannover | olink | 261 (260) | 0 (0) | 0 (0) | 62 | 13 | 186 | plasma | Olink Target 96 Cardiometabolic (v.3603), Cardiovascular II (v.5006), Inflammation (v.3022), Neurology (v.8013); NPX | 359 | no omics (flow counts, serology in the annotation; Bonifacius 2021) |
| Zoodsma2022_Nijmegen | olink | 550 (144) | 0 (0) | 0 (0) | 0 | 550 | 0 | plasma | Olink Target 96 Cardiometabolic (v.3603), Cardiovascular II (v.5006), Inflammation (v.3021); NPX; Olink Target 96 Cardiometabolic (v.3603), Cardiovascular II (v.5006), Inflammation (v.3022); NPX; Olink Target 96 Cardiometabolic (v.3603), Inflammation (v.3021); NPX | 269 | no omics (flow cytometry; Janssen 2021) |

### Proteins shared between cohorts

Number of UniProt keys in common (diagonal: proteins in the cohort).

| | Buggert2025 | Captur2022 | Cervia2024_MS | Cervia2024_SomaScan | Chowdhury2024 | Kovarik2023_plasma | Su2022 | Zoodsma2022_Breda | Zoodsma2022_Hannover | Zoodsma2022_Nijmegen |
|---|---|---|---|---|---|---|---|---|---|---|
| Buggert2025 | 2903 | 50 | 126 | 2147 | 159 | 133 | 430 | 92 | 355 | 268 |
| Captur2022 | 50 | 100 | 26 | 74 | 31 | 29 | 9 | 5 | 11 | 11 |
| Cervia2024_MS | 126 | 26 | 308 | 177 | 218 | 187 | 4 | 0 | 23 | 23 |
| Cervia2024_SomaScan | 2147 | 74 | 177 | 6402 | 216 | 186 | 386 | 85 | 344 | 255 |
| Chowdhury2024 | 159 | 31 | 218 | 216 | 405 | 237 | 5 | 0 | 34 | 33 |
| Kovarik2023_plasma | 133 | 29 | 187 | 186 | 237 | 323 | 4 | 0 | 27 | 27 |
| Su2022 | 430 | 9 | 4 | 386 | 5 | 4 | 441 | 89 | 174 | 174 |
| Zoodsma2022_Breda | 92 | 5 | 0 | 85 | 0 | 0 | 89 | 92 | 92 | 92 |
| Zoodsma2022_Hannover | 355 | 11 | 23 | 344 | 34 | 27 | 174 | 92 | 359 | 269 |
| Zoodsma2022_Nijmegen | 268 | 11 | 23 | 255 | 33 | 27 | 174 | 92 | 269 | 269 |
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
  `never_infected_contemporary`.
- **Acute samples.** 113, with the 6-month outcome in `lc_status`.
- **Flags:**
  - 5 healthy controls (COV_ZH_90/181/185/197/198_ID) match later COVID-19 subjects (`..._2_ID`) by number, sex and age.
    The authors' labels are kept, and the link is in `later_covid_subject_id`.
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
- **Controls.** "No-PCC, convalescent individuals without PCC", 9 healthy hospital co-workers, so `COVID_recovered`.
  - No infection test or date is reported for them.
  - The Limitations section asks for future "inclusion of infected individuals who do not develop PCC".
  - The recovered label rests on the Table 1 footnote only.
- **Run mapping.** Runs are matched to `metadata.sdrf.tsv` by file name.
  - Two SDRF rows (pre-covex26 and pre-covex27) point to the post files; they are mapped to `DIA_Pre_26` and
    `DIA_Pre_27`, which is flagged.
  - Runs absent from the SDRF are `unknown`: `DIA_Pre_13` and `DIA_Post_13` (probably the 10th no-exercise participant,
    not deposited), `DIA_C59_60a` and `DIA_Sheela`.
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
