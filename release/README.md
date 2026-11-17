# CHORUS harmonised datasets, release v1

**Frozen:** 17 November 2026, tag `release-v1`. The files were built by `scripts/build_release.py` from the commit given
in `manifest.csv` (column `commit`), and the tag holds that commit's code and configs. Rebuild with
`python scripts/build_release.py`; it needs the local `data/` folder, which is not distributed.
`tests/test_release.py` checks that every matrix reproduces the values the CHORUS benchmark loads for the same cohort,
samples and features.

## Files

- `manifest.csv` has one row per file with these columns:
  - `path` and `sha256`.
  - `kind`: matrix, samples, splits, signature or baselines.
  - `layer` and `cohorts`.
  - `roles`: the role of each cohort.
  - `n_rows` and `n_features`.
  - `harmoniser`: none, center_scale or fsqn.
  - `n_LC` and `n_COVID_recovered`.
  - `config` and `commit`: what the file was built from.
  - `flags`.
  - `source`: how the test reloads it.
- Matrices are gzip CSV files with samples in rows. The first column `row_id` is the benchmark's pooled id
  (`<cohort>::<sample>`), or the sample id for the paired Su2022 files. The features follow, with 7 significant digits.
- Feature ids are Ensembl gene ids (RNA), UniProt accessions (proteins; a complex is `P1_P2`) and `HMDB:<id>` or
  `RefMet:<name>` (metabolites).
- Values:
  - RNA: log2 CPM, with the library size over all genes.
  - Proteins: as deposited (NPX, log2 RFU, log2 intensity), with missing values filled per cohort as in the benchmark
    configs.
  - Metabolites: as in the benchmark configs.
- Each `*_samples.csv` holds the metadata of the rows of the matrices next to it: cohort, sample and subject id,
  `y_LC` (1 = LC, 0 = COVID_recovered), sex, age, days_since_infection and time anchor where recorded (empty
  otherwise), sample matrix and timepoint.

## Roles and flags

Roles (`roles` in the manifest):
- `training_and_held_out`: a cohort used for training and held out in the CHORUS leave-one-study-out benchmarks.
- `training` / `training_only`: used for training, never held out. **PRJNA1149196_QLD is training only** (5 LC / 4
  COVID_recovered).
- `test_only`: never used to train or to build a pooled estimate. **CORSIP2026** (NanoString panel) and **GSE228320**
  (case = DLCO below 60% of predicted, a lung function outcome, not a Long COVID definition) are test_only. Here they
  appear only as their own Hedges g columns in `signature/transcriptomics_signature.csv`, outside the pooled estimates.
  The early detection test sets (Su2022 and Cervia2024_SomaScan acute samples) are also marked test_only.
- `concordance_only`: Saito2024 (metabolomics), in no CHORUS benchmark; it shares 15 to 28 metabolites with the others.

Flags (`flags` in the manifest):
- **GSE270045: `label_confounded_with_collection`.** LC and COVID_recovered were collected separately, so group is
  completely confounded with sample collection (`docs/gse270045_check.md`). The transcriptomic signature therefore also
  gives the pooled estimate without it (`pooled_unflagged_*`).
- **Cervia2024_MS: `same_participants_as_Cervia2024_SomaScan`** (the same sera on two platforms). It has its own g
  column in the proteomic signature but is not pooled, so those participants count once.

**Woodruff2023 and Bansal2025** (version 5) are in the proteomics files of `clustering/` (affinity v5; Woodruff2023 also
in Explore v5) and `signature/`.

## early_detection/ (for a prediction method developed elsewhere)

The task is the later outcome, LC versus COVID_recovered, predicted from a sample drawn in the acute phase.

- **Development set:** `gse215865_development_log2cpm.csv.gz`, `_samples.csv`, `_splits.csv`.
  - GSE215865 (Thompson 2023) hospital whole blood: the earliest library of each subject, 424 subjects (302 LC, 122
    COVID_recovered), on the 15,116 genes expressed at timepoint 0.
  - `is_timepoint0` marks the 345 subjects whose earliest library is the timepoint 0 library (241 LC, 104 recovered).
    The other 79 have a later first library.
  - Labels are the later Long COVID outcome from the TACO clinical table.
  - The split file gives, for the timepoint 0 samples, the folds of the CHORUS reference baseline: 10 repeats of
    stratified 5-fold CV (`fold_repeat0` to `fold_repeat9`; repeat r uses StratifiedKFold random_state r). It is empty
    for the other 79.
- **Frozen test sets:** `acute_proteins_386_none.csv.gz`, `acute_proteins_386_center_scale.csv.gz`, `_samples.csv`,
  `_splits.csv`.
  - Acute samples of Su2022 (Olink, 75 LC / 50 recovered) and Cervia2024_SomaScan (40 / 73), one per subject (the
    earliest acute draw), on their 386 shared proteins.
  - `none` is without harmonisation. `center_scale` is a per-cohort z score of each protein, fitted on that cohort's own
    samples without labels.
  - The split file gives the two directions of the baseline: train on one cohort, test on the other.
- **Multi-layer:** `su2022_paired_olink.csv.gz` (441 proteins) and `su2022_paired_metabolomics.csv.gz` (891
  metabolites).
  - The acute draws of the same 124 Su2022 participants (74 LC / 50 recovered); missing values are kept as missing.
  - `su2022_paired_splits.csv` holds the 10 × 5-fold CV folds of the CHORUS baseline.
- **Baselines** (`baselines.csv`), the numbers to beat:

  | set | evaluation | best CHORUS baseline (AUROC) |
  |---|---|---|
  | acute_proteins_386 | train on one cohort, test on the other | Cervia2024_SomaScan: 0.644 (identity / differential / logreg, trained on Su2022); Su2022: 0.608 (center_scale / differential / rf, trained on Cervia2024_SomaScan) |
  | su2022_paired | 10 × 5-fold CV, mean over folds | 0.516 (late_fusion) |
  | gse215865_development (timepoint 0) | 10 × 5-fold CV, mean over folds | 0.569 (differential / rf); published TACO 0.675 (not recomputed) |

  Every baseline configuration, with intervals and the label-permutation p values of the horizontal analysis, is in
  `baselines.csv` (from `results/early_detection`).

**The test sets are to be scored once.** Develop and tune on GSE215865, or within the given folds. Then score the frozen
acute protein test sets one time, in the given direction, and report that result against the baselines above.

## clustering/ (for a subtyping method developed elsewhere)

LC samples only, per layer.

| layer | files | cohorts | LC samples | features | config |
|---|---|---|---|---|---|
| transcriptomics_v4 | `transcriptomics_v4_{none,center_scale,fsqn}.csv.gz` | 8 (PRJNA1149196_QLD training only; GSE270045 flagged) | 148 | 10,053 shared genes | `configs/lc_transcriptomics_blood_v4.yaml` |
| proteomics_affinity_v5 | `proteomics_affinity_v5_*` | Buggert2025, Cervia2024_SomaScan, Su2022, Fineschi2026, Woodruff2023, Bansal2025 | 458 | 55 shared proteins | `configs/lc_proteomics_affinity_v5.yaml` |
| proteomics_explore_v5 | `proteomics_explore_v5_*` | Buggert2025, Fineschi2026, Woodruff2023 | 331 | 358 | `configs/lc_proteomics_explore_v5.yaml` |
| proteomics_ms_v4 | `proteomics_ms_v4_*` | Cervia2024_MS (flagged), Chowdhury2024, Kovarik2023_plasma, Chowdhury2026 | 127 | 179 | `configs/lc_proteomics_ms_v4.yaml` |
| metabolomics | `metabolomics_*` | Martinez2024, Lopez2023, Su2022_metabolomics | 138 | 38 | `configs/lc_metabolomics.yaml` |

- **Versions:**
  - `none`: without harmonisation.
  - `center_scale`: per-cohort z score of each feature.
  - `fsqn`: feature-specific quantile normalisation of each cohort to the pooled distribution of that feature.
  - Both harmonisers are fitted on the LC samples of all cohorts of the layer, without labels. They are the two that
    removed the cohort structure in `results/clustering_readiness_v5`.
- **Samples file:** `<layer>_samples.csv` gives cohort, sex, age and days since infection where recorded. Cohort is
  kept so that the adjusted Rand index between clusters and cohort can be checked.
- **Caution:** clusters that follow the cohort column are technical, not subtypes.

## signature/ (for drug repurposing)

One file per layer: `transcriptomics_signature.csv`, `proteomics_signature.csv`, `metabolomics_signature.csv`. One row
per feature, with `symbol` (HGNC symbol, or RefMet name for metabolites).

- **`g_<cohort>`:** Hedges g of LC minus COVID_recovered within each cohort, on that cohort's own features (post-acute
  samples, one per subject; configs `concordance_rna_v5`, `concordance_proteomics_all_v5`,
  `concordance_metabolomics_all`). It is empty where the cohort does not measure the feature.
- **`pooled_*`:** a DerSimonian-Laird random-effects estimate over the pooled cohorts, for features measured in at least
  2 of them. Columns: g, se, 95% interval, p, Q, Q_p, I2, tau2, `n_cohorts` and `n_agree_sign` (cohorts whose g has the
  sign of the pooled g).
- **`pooled_direction_inconsistent`:** True when at least one cohort's g has the opposite sign to the pooled g, so this
  is a strict rule. `n_agree_sign` and `n_cohorts` allow looser thresholds.
- **Not pooled:** CORSIP2026 and GSE228320 (test_only) and Cervia2024_MS (same participants) keep their own g columns
  but are not pooled.
- **Transcriptomics also:** `pooled_unflagged_*` gives the same estimate without GSE270045.

Counts per file:

| file | features | pooled (2 or more cohorts) | pooled cohorts (max) | direction consistent | direction inconsistent | pooled p < 0.05 and consistent |
|---|---|---|---|---|---|---|
| `transcriptomics_signature.csv` | 42781 | 20083 | 8 | 2499 | 17584 | 525 |
| `proteomics_signature.csv` | 7455 | 3199 | 7 | 528 | 2671 | 161 |
| `metabolomics_signature.csv` | 1317 | 226 | 4 | 109 | 117 | 13 |

