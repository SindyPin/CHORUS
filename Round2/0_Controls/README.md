# 0_Controls

Known answer tasks (PROTOCOL.md with Amendments 1, 2 and 3). They show what adding one cohort, or one layer, does to
the prediction of two variables that every blood sample carries:
- **sex:** classification, female = 1 and male = 0;
- **age:** regression, in years.

**All results rest on 15 folds:** 5 folds by 3 repeats, the first 3 of the 10 seeded repeats (seed 20261120;
Amendment 3).
- Step 01 ran up to 50 folds for most jobs before Amendment 3. Those raw files are kept, and every summary reads folds 0
  to 14 only.
- The first 15 folds of the 10-repeat split are the same as the 3-repeat split; this was checked.

## Answer per layer

"Improves" means the paired difference between arms, on the same folds and oriented so that positive is better, has a
95% interval (corrected resampled t test, Nadeau and Bengio) that excludes zero. Primary metric: AUROC for sex, mean
absolute error for age. Counts are numbers of models.

### Does adding one cohort improve the prediction? (B_scaled vs A_shared)

Every pair here is measured on different platforms or with different pipelines. B_scaled is the interpretable arm.
B_naive is reported only (Amendment 1, item 5).

| control | layer | anchor + partner | improves | no change | worse | answer |
|---|---|---|---|---|---|---|
| sex | transcriptomics | GSE215865 + GSE198449 | 0 | 19 | 0 | no |
| sex | proteomics | Su2022 + Zoodsma2022_Hannover | 0 | 7 | 11 | no |
| sex | metabolomics | Su2022_metabolomics + Holmes2021 | 0 | 18 | 0 | no |
| age | transcriptomics | GSE198449 + GSE161731 | 0 | 9 | 8 | no |
| age | proteomics | Su2022 + Zoodsma2022_Hannover | 0 | 12 | 4 | no |
| age | metabolomics | Su2022_metabolomics + Lopez2023 | 0 | 11 | 5 | no |

B_naive vs A_shared, reported only:
- sex: no model improves in any layer;
- age: one model improves (metabolomics).

### Does adding one layer improve the prediction? (joined and late fusion vs each layer alone)

The same participants are used in every arm:
- Su2022: layer 1 is Olink (441 proteins) and layer 2 is plasma metabolomics (679); 225 participants.
- Cervia2024: layer 1 is SomaScan (6,402) and layer 2 is MS proteomics (308); 113 participants.

| control | pair | joined vs layer 1 | joined vs layer 2 | late fusion vs layer 1 | late fusion vs layer 2 | answer |
|---|---|---|---|---|---|---|
| sex | Su2022 | 5 / 13 / 0 | 7 / 11 / 0 | 3 / 15 / 0 | 9 / 9 / 0 | yes in some models, over either layer alone; none gets worse |
| sex | Cervia2024 | 0 / 18 / 0 | 16 / 2 / 0 | 0 / 18 / 0 | 17 / 1 / 0 | over MS proteomics alone: yes in most models; over SomaScan alone: no |
| age | Su2022 | 7 / 9 / 0 | 2 / 14 / 0 | 11 / 5 / 0 | 1 / 15 / 0 | over Olink alone: yes in about half or more of the models; over metabolomics alone: in 1 or 2 models |
| age | Cervia2024 | 0 / 16 / 0 | 14 / 2 / 0 | 0 / 13 / 3 | 16 / 0 / 0 | over MS proteomics alone: yes in most or all models; over SomaScan alone: no (late fusion is worse in 3) |

Cells are improves / no change / worse.

### Primary metric by arm: median over models (min to max)

sex, AUROC

| layer | A_all | A_shared | B_naive | B_scaled |
|---|---|---|---|---|
| transcriptomics | 0.999 (0.648 to 1.000) | 0.999 (0.648 to 1.000) | 1.000 (0.613 to 1.000) | 1.000 (0.642 to 1.000) |
| proteomics | 0.569 (0.520 to 0.738) | 0.590 (0.530 to 0.714) | 0.569 (0.494 to 0.666) | 0.563 (0.499 to 0.612) |
| metabolomics | 0.649 (0.613 to 0.700) | 0.632 (0.594 to 0.693) | 0.636 (0.591 to 0.689) | 0.639 (0.577 to 0.701) |

| pair | layer 1 | layer 2 | joined | late fusion |
|---|---|---|---|---|
| Su2022 | 0.800 (0.674 to 0.855) | 0.794 (0.706 to 0.817) | 0.832 (0.691 to 0.888) | 0.841 (0.690 to 0.891) |
| Cervia2024 | 0.998 (0.974 to 1.000) | 0.882 (0.774 to 0.928) | 0.994 (0.927 to 1.000) | 0.989 (0.953 to 0.997) |

age, mean absolute error in years (lower is better)

| layer | A_all | A_shared | B_naive | B_scaled |
|---|---|---|---|---|
| transcriptomics | 1.545 (1.405 to 3.239) | 1.535 (1.405 to 3.214) | 1.705 (1.495 to 3.541) | 1.706 (1.483 to 8.246) |
| proteomics | 11.008 (10.237 to 15.311) | 11.206 (10.238 to 15.029) | 11.280 (10.359 to 13.717) | 11.427 (10.380 to 14.955) |
| metabolomics | 10.181 (9.038 to 14.547) | 10.965 (10.064 to 15.503) | 10.953 (10.218 to 11.931) | 11.077 (10.291 to 14.891) |

| pair | layer 1 | layer 2 | joined | late fusion |
|---|---|---|---|---|
| Su2022 | 10.904 (9.750 to 13.845) | 10.311 (8.560 to 14.719) | 10.105 (8.260 to 13.324) | 10.029 (8.378 to 12.032) |
| Cervia2024 | 6.896 (5.236 to 12.404) | 10.072 (7.799 to 12.634) | 7.030 (5.240 to 12.572) | 7.816 (5.730 to 10.173) |

Facts needed to read these tables:
- **Su2022 values (Olink and metabolomics) are the authors' Z scores.** They are already adjusted for age, sex and BMI
  (Round1 cohort card).
- **The age RNA anchor GSE198449 records ages 18 to 34 only.**

## Models not run, and why (Amendment 2, item 3)

Model grid:
- 13 classical models;
- the TACO ensemble;
- TabPFN v2, TabPFN 2.5, TabICL v2, TabDPT and Mitra;
- on RNA only, BulkFormer and BulkRNABert.

| model | where | reason |
|---|---|---|
| Mitra | every task and arm (40 jobs) | Amendment 2: time budget exceeded on this hardware. It took 510 to 769 s per fit on the CPU, with no usable GPU, measured on step 01 fold 0. The timing rows are kept in each step's `timing.csv` and the fold 0 fits in `01_one_cohort/output/mitra_timed/`. |
| TabPFN 2.5 | every age (regression) job (18 jobs), and the 2 age late fusion jobs that depend on them | licence or token. The 2.5 regressor weights are not on this machine, and with no TabPFN API key here (`TABPFN_TOKEN` or `~/.cache/tabpfn/auth_token`), tabpfn cannot download them. It was retried under Amendment 3, item 3, and failed again. The sex (classification) jobs ran. |
| Gaussian naive Bayes | every age job (20 jobs) | no regression form |
| BulkRNABert | every RNA job (8 jobs) | time. It took 10.2 s per sample on the CPU and ran out of memory on MPS at 2 samples, so embedding 475 to 508 samples would take 81 to 86 minutes, over the 15 minute limit per fit. |

Everything is listed in `04_summary/output/failures.csv` and `models_not_run.csv`. No other job failed.

## Cohorts

The anchor and partner are chosen by the protocol rule among the cohorts that record the variable, whatever their LC
status:
- the anchor is the largest eligible cohort of the layer, and the partner the next largest;
- for sex, at least 10 per class; for age, at least 20 participants;
- one sample per participant: the earliest by days since infection, then by time point, then by sample id.

| control | layer | anchor (n) | partner (n) | shared features |
|---|---|---|---|---|
| sex | transcriptomics | GSE215865 (508) | GSE198449 (475) | 58,884 |
| sex | proteomics | Su2022 (493) | Zoodsma2022_Hannover (260) | 174 |
| sex | metabolomics | Su2022_metabolomics (377) | Holmes2021 (27) | 35 |
| age | transcriptomics | GSE198449 (475) | GSE161731 (151) | 56,939 |
| age | proteomics | Su2022 (493) | Zoodsma2022_Hannover (260) | 174 |
| age | metabolomics | Su2022_metabolomics (377) | Lopez2023 (85) | 61 |

- **Excluded from the choice:** the flagged cohorts, Buggert2025 (full cohort) and GSE300129 (Round1 sex check
  failed).
- **Added cohorts:** Holmes2021 and Zoodsma2022_Nijmegen are on disk but in no Round1 config. They are added in
  `common.py` (`ON_DISK`). Nijmegen has 144 participants and changes no choice.
- **Sex as recorded.** Every chosen cohort records sex in its metadata; none has only sex inferred from expression. So
  no chrX or chrY genes were removed.
- **No anchor sensitivity.** In no layer do the three largest cohorts differ by fewer than five participants.
- **No pair shares a pipeline and reference** (Amendment 1, item 5, as read on 20 November 2026), so B_naive is
  reported only.

## Methods in short

**Preprocessing.**
- RNA counts are converted to log2 CPM.
- Other layers are used as the Round1 loader returns them, on the log scale.
- Median imputation and standardisation use training fold statistics.

**Feature selection, on the training fold.**
- sex: Welch |t|, top 200;
- age: |Spearman|, top 200;
- step 03: 100 features per layer.

**Inner search.** 3-fold grid search on the fixed grids of PROTOCOL.md.
- Regression forms: logistic regression becomes Ridge (alpha 0.1, 1, 10, 100). Elastic net becomes ElasticNet (alpha
  0.01, 0.1, 1; l1_ratio 0.2, 0.5, 0.8). The linear and RBF SVMs become SVR, and kNN, the tree models, the boosting
  models and the MLP become their regressors, on the same grids.
- CatBoost (300 iterations) and histogram gradient boosting (default 100 iterations) each search 4 settings. That
  already met Amendment 3, item 2, so nothing changed.

**Other models.**
- **TACO ensemble:** the Round1 `make_classifier("taco")`, without TabPFN. For age, the same four members as
  regressors.
- **BulkFormer-37M:** sample-level embeddings from raw counts, then tuned logistic or ridge regression.
- **Mitra:** was run in a separate environment, `Code/.venv-mitra`.

**Arms.**
- **B_scaled:** each cohort z-scored on its own statistics; the anchor's held-out fold is scaled with its training
  fold statistics.
- **Late fusion:** the mean of the layer 1 and layer 2 predictions of the same model.

## Folders

| folder | content |
|---|---|
| `01_one_cohort` | candidates, anchors and partners, A_all |
| `02_horizontal_one_plus_one` | A_shared, B_naive, B_scaled |
| `03_vertical_one_plus_one` | layer_1, layer_2, joined, late_fusion |
| `04_summary` | metrics by arm, paired differences, counts, failures, answers, models not run |
| `common.py` | shared code: cohorts, folds, preprocessing, models, runner, timing |
| `mitra_worker.py` | Mitra in its own environment |
| `write_inputs.py` | writes `inputs.csv` (files read, with SHA-256) in every step folder |

Each step can resume. Before the full run, each step times fold 0 of every job and prints the estimate; the step runs
only when the estimate is within 12 h (Amendments 1 to 3). The estimated wall clock for 15 folds, models run, was:

| step | estimate |
|---|---|
| 01 | 1.0 h |
| 02 | 2.6 h |
| 03 | 0.4 h |

Not committed: `cache/` (arrays and participant ids) and `logs/`.

Run from `Code/`:

```
python Round2/0_Controls/01_one_cohort/run.py
python Round2/0_Controls/02_horizontal_one_plus_one/run.py
python Round2/0_Controls/03_vertical_one_plus_one/run.py
python Round2/0_Controls/04_summary/run.py
python Round2/0_Controls/write_inputs.py
```
