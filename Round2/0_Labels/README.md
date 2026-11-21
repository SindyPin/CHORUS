# 0_Labels: labels, clinical variables, anchors and partners

## Purpose

Record, for every cohort on disk:

- how its Long COVID (LC) label and comparison group are defined;
- where each label comes from;
- when the samples were taken;
- how many participants fall in each class;
- which clinical and symptom variables are recorded.

From these, pick the anchor and partner of each task and layer by the rules of `Round2/PROTOCOL.md`, and list every
case where the same participants have two layers. Every Round2 task reads its labels from here.

## Input

- **Cohort metadata:** `data/meta/<cohort>_meta.csv`, the ClinicalSchema tables written by the Round1 builders.
  44 files are read; the 3 methylation and 3 unlabelled RNA series have no metadata file, so their counts come from
  the Round1 cards.
- **Definitions:** taken from the Round1 cohort cards:
  - `Round1/docs/proteomics_cohorts.md`
  - `Round1/docs/metabolomics_cohorts.md`
  - `Round1/docs/paired_cohorts.md`
  - `Round1/docs/acute_covid_cohorts.md`
  - `Round1/docs/dataset_sweep_2026-11.md`
  - `data/meta/label_sources.md`
  - `data/meta/curation_notes.md`
- **Matrices:** the Round1 matrices, through the Round1 configs and `chorus.io.Cohort`, for the feature counts.
- **Open clinical tables:** listed in the table below.
- **Checksums:** every file read, with its SHA-256, is in `inputs.csv`.

| cohort | file | source | licence | linked |
|---|---|---|---|---|
| GSE224615 | `data/clinical/GSE224615/41590_2023_1724_MOESM1_ESM.pdf` (Supplementary Table 3) | Yin 2024 Nat Immunol, Springer | CC BY 4.0 | 36 of 36 participants, through the person id in the GEO sample title (`942_Non-PASC_WT`) |
| GSE226260_WB | `data/clinical/GSE226260/41590_2025_2353_MOESM4_ESM.xlsx` (Demographics, Symptom List) | Nat Immunol 2025, PMID 41388153 | CC BY 4.0 | 38 of 38 subjects (`BI4001`) |
| GSE226260_PBMC | `data/clinical/GSE226260/41590_2025_2353_MOESM2_ESM.xlsx` | same paper | CC BY 4.0 | 35 of 36 labelled participants (SN = GEO subject id) |
| GSE226260 | `data/clinical/GSE226260/41590_2025_2353_MOESM3_ESM.xlsx` (sample sheet) | same paper | CC BY 4.0 | downloaded, not linked: its patient ids match only 35 of the 135 non-BI GEO subject ids |
| Lopez2023 | `data/metabolomics/Lopez2023/raw/DATA AVAILABILITY.xls` (on disk from Round1) | Mendeley 10.17632/8zfdjsypd8.1 | CC BY 4.0 | 48 of 48 labelled participants (Folio, by the Round1 id rule) |
| Fineschi2026 | `data/counts/Fineschi2025_check/SamplesPC.txt` (on disk from Round1) | figshare 10.17044/scilifelab.28832492 | CC BY 4.0 | 110 of 165 (the people also in the PBMC study) |
| CORSIP2026 | `data/counts/CORSIP2026/raw/CORSIP_omics/CORSIP_Full.csv` (on disk from Round1) | github.com/estefespin/CORSIP_omics | CC0 | 37 of 37 |

The URLs, checksums, row counts and label agreement are in `output/clinical_links.csv`. The linked tables are written
to `data/clinical/<cohort>/<cohort>_clinical_linked.csv`; they are not committed.

**Checked, nothing per participant open:** GSE297499, GSE251849, GSE228320, Kovarik2023, Fineschi2026 (Olink deposit),
and Sinclair 2024 (PRJNA1149196).
- PRIDE PXD050202 gives a condition per proteomics run, but has no key to the SRA samples.
- GSE270045 (Shahbaz 2025, S Table 1) has age, sex and hospitalisation for LC1–LC78 and R1–R62, but no mapping to the
  GEO ids. It is not downloaded.
- Saito2024 Supplementary Table 1 uses LC1–LC32 with gaps, while the data use L1–L30. It is not linked (the Round1
  finding).

**Already in the Round1 metadata:**
- Chowdhury2024 age (SDRF);
- the Cervia2024 clinical variables;
- the Bansal2025 baseline variables. The serology columns of `baseline57.xlsx` were not added.

## Method

Script: `build_labels.py`, run from `Code/` as `python Round2/0_Labels/build_labels.py`.

1. **Phase.** Each sample's phase is the metadata `sample_phase`. Four older tables have no such column:
   - GSE215865 is set to acute;
   - GSE224615, GSE226260_PBMC, GSE251849 and GSE270045 are set to post-acute;
   - GSE226260_WB is post-acute except its four `Acute` draws.
2. **Counting.** Participants are counted per class and phase, one per participant:
   - post-acute LC and recovered (task 3);
   - acute samples with a later LC or recovered outcome (task 2);
   - never infected;
   - unknown.

   A label is never inferred. A participant without a stated label is unknown, and so is a participant whose sources
   disagree (GSE226260_WB, below).
3. **Eligibility per task:**
   - 3_Detection: post-acute LC against recovered.
   - 2_EarlyDetection: acute samples with the later outcome.
   - 1_Subtyping: post-acute LC participants.
     - The anchor is ranked by LC count, and features come from LC against recovered.
     - Since Amendment 1 (item 2), the anchor and partner also need at least 10 recovered. Cohorts with cases only
       are listed, not used.
     - An anchor with fewer than 50 LC is marked "too small to interpret".
4. **Anchor and partner.**
   - **Layer.** Layers are grouped as transcriptomics, proteomics (Olink, SomaScan and MS together), metabolomics and
     methylation.
   - **Anchor.** The largest cohort of the layer with at least 10 participants per class. Excluded:
     - GSE270045, GSE228320 and CORSIP2026;
     - since Amendment 1 (item 1), the full Buggert2025. It is split by `collection_country` into Buggert2025_UK
       (51 / 51) and Buggert2025_SE (95 LC, cases only), which share the Buggert2025 matrix.
   - **Partner.** The next largest cohort under the same rule. **Interpretation:** the partner must also have 10 per
     class.
   - **Anchor sensitivity (Amendment 1, item 3).** The three largest eligible cohorts (two if only two are eligible) are
     listed in `sensitivity_anchors` when their sizes span fewer than five participants. **Interpretation:** "differ by
     fewer than five" is read as the range of the three sizes.
   - **Different platforms or scales (Amendment 1, item 5).** `different_platform_or_scale` is true when the anchor's
     and partner's platform descriptions differ.
     - **Reading of 20 November 2026:** two cohorts quantified with the same pipeline and the same reference count as
       the same scale, so B_naive is interpretable for them.
     - Here these are the Round1 re-quantifications: salmon 1.10.3 against `data/ref/gencode_v47/salmon_index_v47_k31`,
       per each run's `cmd_info.json`. They cover GSE224615, PRJNA1149196_SA, PRJNA1149196_QLD and PRJNA1184005.
     - So PRJNA1149196_SA with GSE224615 (RNA detection and subtyping) is the same scale. Every other pair is
       flagged.
   - **Shared features.** Feature keys present in both matrices as the Round1 loader reads them:
     - transcriptomics: Ensembl genes, before any expression filter;
     - proteomics: UniProt keys;
     - metabolomics: shared metabolite identifiers.
5. **Vertical cases.** Two kinds are listed:
   - Two omics layers on the same participants: Cervia2024 and Su2022, counted from the shared `subject_id`. The
     others are listed with the reason they are not usable.
   - Omics plus clinical variables, for every cohort with recorded variables.
     - **Rule (Amendment 1, item 4).** `clinical_layer` is true when at least three variables are recorded for at
       least 80 percent of the task's participants.
     - **Excluded variables:** those that define the label (symptom lists, counts and scores), timing variables and
       treatment variables.
     - **Early detection:** only variables known at the acute stage count. These are demographics, comorbidity, acute
       severity (peak severity of the acute illness included), infection wave and vaccination before infection.

## Output

| file | content |
|---|---|
| `output/labels_audit.csv` | one row per cohort (52, with Buggert2025_UK and Buggert2025_SE): layer, tissue, platform, paper, LC definition, label source, comparison group, time since infection, sample timing, participants per class and phase, clinical variables, label flag, doubts, which tasks it can serve |
| `output/clinical_variables.csv` | one row per cohort and variable (241 rows, 40 cohorts): kind, whether it defines the label, source file, participants recorded (all and labelled) |
| `output/anchors_and_partners.csv` | horizontal rows: per task and layer, the anchor, partner, class counts, feature counts, shared features and every candidate with its counts. Vertical rows: two-layer cases and omics plus clinical cases, with participants per class |
| `output/clinical_links.csv` | each open clinical table: URL, licence, SHA-256, rows, rows linked, label agreement |
| `output/sample_labels.csv` | per sample: study, subject, sample, label used in Round2, symptoms only label (Amendment 7; the same as the label except GSE215865), phase, flag. Not committed (`.gitignore`), since it lists every sample id |
| `inputs.csv` | files read, with SHA-256 |

## Result

### Horizontal (Amendment 1)

Participants are LC / recovered; "shared" is the number of shared features. Scale:
- **Same scale:** the two RNA pairs (PRJNA1149196_SA with GSE224615) share pipeline and reference, so B_naive is
  interpretable.
- **Different platforms or scales:** every other pair, so only B_scaled is interpretable there.

| task | layer | anchor | partner | shared | note |
|---|---|---|---|---|---|
| 3_Detection | transcriptomics | PRJNA1149196_SA 28 / 11 | GSE224615 23 / 13 | 78,277 genes (before filtering) | sensitivity anchors: PRJNA1149196_SA, GSE224615 and GSE226260_WB (39, 36 and 35 participants) |
| 3_Detection | proteomics | Fineschi2026 92 / 73 | Su2022 75 / 50 | 186 | Buggert2025 flagged; Buggert2025_UK 51 / 51 ranks below |
| 3_Detection | metabolomics | Su2022_metabolomics 68 / 45 | Martinez2024 40 / 40 | 158 | |
| 3_Detection | methylation | none | | | GSE197152 15 LC, cases only |
| 2_EarlyDetection | transcriptomics | GSE215865 159 / 29 (Amendment 7 labels) | none | | |
| 2_EarlyDetection | proteomics | Su2022 75 / 50 | Cervia2024_SomaScan 40 / 73 | 386 | |
| 2_EarlyDetection | metabolomics | Su2022_metabolomics 74 / 50 | none | | Lopez2023 18 / 4 |
| 1_Subtyping | transcriptomics | PRJNA1149196_SA (28 LC) | GSE224615 (23 LC) | 78,277 | too small to interpret (fewer than 50 LC) |
| 1_Subtyping | proteomics | Woodruff2023 (93 LC, 27 recovered) | Fineschi2026 (92 LC) | 358 | Buggert2025_SE (95 LC, cases only) listed, not used |
| 1_Subtyping | metabolomics | Su2022_metabolomics (68 LC) | Martinez2024 (40 LC) | 158 | |
| 1_Subtyping | methylation | none | | | GSE197152 cases only |

Changes from the table before Amendment 1:
- **Proteomics detection:** Buggert2025 to Fineschi2026 as anchor, and Fineschi2026 to Su2022 as partner.
- **Proteomics subtyping:** Buggert2025 to Woodruff2023 as anchor, and Woodruff2023 to Fineschi2026 as partner.
- **Methylation subtyping:** now has no anchor.

Two choices are close:
- **Proteomics subtyping:** Woodruff2023 (93 LC) and Fineschi2026 (92) differ by one, but the third, Su2022 (75), is
  further away, so item 3 does not apply under the range reading.
- **RNA detection:** GSE226260_WB would have had 37 participants (17 / 20) before BI4036 and BI4002 became unknown.

### Vertical

- **Cervia2024:** SomaScan plus MS, 40 LC / 73 recovered at 6 months. There is no acute MS sample.
- **Su2022:** Olink plus metabolomics.
  - Post-acute: 68 / 45.
  - Acute with outcome: 74 / 50.
- **Omics plus clinical:** 26 rows. 17 meet the clinical layer rule, and 14 of those also have 10 per class:
  - 3_Detection: Bansal2025, Buggert2025 (flagged), Buggert2025_UK, Cervia2024_MS, Cervia2024_SomaScan,
    GSE224615, GSE226260_WB, GSE228320 (flagged), Lopez2023, Su2022 and Su2022_metabolomics.
  - 2_EarlyDetection: Cervia2024_SomaScan, Su2022 and Su2022_metabolomics.

  The variables recorded for at least 80 percent are in `clinical_variables_80pct`.
  - **Fails the rule:** GSE215865 (sex only), Woodruff2023 (hospitalisation only; sex missing for 91 of 93 PASC),
    Chowdhury2026 (sex only), Fineschi2026 and PRJNA1149196_SA.
  - **Lopez2023:** its sex and comorbidity coding is not documented (below).
- **Listed, not usable now:**
  - Fineschi2026 / Fineschi2025: 110 shared people, but the PBMC RNA matrix is open and not downloaded.
  - CORSIP2026: flagged, and its metabolomics is not on disk and not keyed to the PX ids.
  - PRJNA1149196: no key to the proteomics runs.
  - GSE297499: the CSF layer is not built, and PBMC has 17 / 9.
  - GSE226260 WB and PBMC: the id systems are unlinked.
  - Kovarik2023: supplement labels do not link per sample.

## Doubts about labels

Each doubt is also in the `doubts` column of `labels_audit.csv`.

**Flagged cohorts (listed, not used as anchor or partner):**
- **GSE270045:** labels confounded with collection; GEO has 19 LC against 18 in the paper.
- **GSE228320:** a non standard case definition (DLCO after ARDS).
- **Buggert2025 (full cohort, Amendment 1):** country coincides with LC status. Buggert2025_UK and Buggert2025_SE are
  used instead.
- **CORSIP2026:** a gene panel only.

**Changed in Round2:**
- **GSE226260_WB, BI4036 and BI4002: set to unknown in Round2** (Round1 metadata unchanged). Within the paper's own
  supplement:
  - the Demographics and Symptom List sheets call BI4036 PASC and BI4002 Recovered;
  - the Medication sheet and the GEO titles say the opposite.

**Labels that rest on a single or secondary source:**
- **GSE215865 (Amendment 7, 21 November 2026):** the TACO labels are withdrawn. Their rule kept each subject's first
  row, so they reflect the PCR based recovery status at admission, not post COVID symptoms
  (`GSE215865_checklist/taco_rule.py` reproduces all 489).
  - Labels now come from the post COVID checklist of the Synapse table (syn35876735), rebuilt in
    `GSE215865_checklist/`: main, the paper's definition (159 LC, 29 recovered of the 508 GEO subjects); sensitivity,
    symptoms only (156 / 32), in the column `lc_status_symptoms_only` of `sample_labels.csv`.
  - COVID negative in GEO: never_infected (71). Never answered the checklist: unknown (249), excluded.
  - The same labels are used for GSE225349 (SomaScan, `GSE225349/`).
  - This cohort is still the RNA anchor for early detection (159 / 29); it has no partner.
- **GSE197152:** "LC" rests on the paper's statement that all 15 patients had PASC symptoms. There is no per sample
  field.
- **PRJNA1184005:** the labels come from the author, and the control definition is pending confirmation.

**Comparison groups that mix or overlap:**
- **Buggert2025:** the 95 Swedish participants are all LC, and the Swedish controls were not profiled. Since
  Amendment 1 the cohort is split by country; only Buggert2025_UK has both classes.
- **Cervia2024:** "recovered" mixes:
  - patients recovered within a month;
  - "recovered Long Covid";
  - smell- or taste-only patients.
- **Chowdhury2024:** the 9 controls are convalescent without a per person infection record.
- **Chowdhury2026:** 5 runs are coded Lc in the file names and C in the SDRF, and the SDRF is used. Overlap with
  Chowdhury2024 cannot be excluded.
- **Saito2024:** group C is taken as recovered by elimination.

**Definitions that differ from WHO:**
- **Su2022:** PASC is assessed once, at about 2 to 3 months, with no duration requirement.
- **Woodruff2023:** PASC is new or worsening symptoms more than 14 days after onset, and 38 of 93 PASC samples are
  under 90 days after onset.
- **Martinez2024:** the case definition was not read (the paper is behind a subscription).

**Counts that differ from the paper:**
- **GSE224615:** GEO has 23 / 13; the paper cohort is 27 / 16.
- **Lopez2023:** the deposit differs from the paper by one per group.

**Labels or timing not assignable:**
- **GSE226260_PBMC:** post-acute is assumed from the PASC / NOPASC status. GEO gives collection dates only.
- **GSE251849:** matrix columns cannot be mapped to GSM, so sex per column is unknown.
- **Every participant unknown:** Holmes2021 (undocumented symptom coding), Captur2022 (no labels deposited),
  GSE247869 (PASC subgroup not public), GSE169687, GSE267625 and GSE222253 (labels requested from the authors).

**Clinical variables:**
- **Woodruff2023:** sex is recorded for 2 of 93 PASC against all 27 recovered, so missing sex alone separates the
  classes.
- **Lopez2023:** the sex and symptom coding of the deposit is not documented. Every never infected control has
  FEVER = 1, so 1 does not mean present. The 19 symptom columns are recorded but not interpreted, and they count as
  label defining.
- **Fineschi2026:** the PBMC study's Score (Post-COVID Symptom Severity Score) and Dyspnoea are linked for 110 people.
  They count as label defining. Whether the PBMC and plasma draws are the same is not stated.
- **Buggert2025:** the Borg CR10, dyspnoea-12 and Nijmegen scores are treated as label defining and left out of the
  clinical layer.
