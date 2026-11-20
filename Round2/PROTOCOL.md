# Round2 protocol

Written and committed on 20 November 2026, before any Round2 experiment. Every Round2 task follows it. A task
that cannot follow a rule says so in its README, and does not change the rule.

## Tasks

- `1_Subtyping`: clustering of Long COVID (LC) participants.
- `2_EarlyDetection`: LC against recovered, from samples taken in the acute phase (outcome known later).
- `3_Detection`: LC against recovered, from samples taken after recovery from the acute infection.

`0_Labels` holds the labels, clinical variables, anchors and partners that every task reads.

## Randomness and resampling

- Seed 20261120 for everything (folds, inner searches, models, clustering, bootstrap).
- Repeated stratified 5 fold cross validation, 10 repeats (50 folds), on the anchor cohort only. The folds are
  drawn once per task and layer and reused by every arm and model.

## Anchor and partner

- Anchor: the largest cohort for the task and layer with at least 10 participants per class.
- Partner: the next largest cohort of the same layer.
- Cohorts with flagged labels are listed but not used as anchor or partner:
  - GSE270045: labels confounded with collection.
  - GSE228320: non standard case definition.
  - CORSIP2026: gene panel only.
- One sample per participant.

## Horizontal arms

Always tested on the anchor's held out fold.

| arm | training data | features |
|---|---|---|
| A_all | anchor training fold | all anchor features |
| A_shared | anchor training fold | features shared with the partner |
| B_naive | anchor training fold plus all partner samples | shared features, no adjustment |
| B_scaled | as B_naive | shared features, each feature z scored within each cohort without labels |

In B_scaled:
- The partner is scaled with its own mean and standard deviation.
- The anchor's training fold is scaled with its own mean and standard deviation.
- The anchor's held out fold is scaled with the training fold's mean and standard deviation.

## Vertical arms

Use the same participants and the same folds. The folds are drawn on the participants who have both layers.

- Layer 1 alone.
- Layer 2 alone.
- Both joined (early fusion).
- The average of the two single layer models' probabilities (late fusion).

Clinical variables count as a layer where they are recorded.

## Preprocessing

All statistics come from the training fold.

- RNA: log2 CPM.
- Proteins and metabolites: on the log scale as deposited.
- Microbiome: centred log ratio.
- Missing values: filled with the training fold median.
- Features: standardised with training fold statistics.

## Feature selection

Feature selection is the same in every task and model.

- Welch t test between the two classes, on the training fold only.
- Keep the top 200 features by absolute statistic, or all features if there are fewer.
- Vertical arms keep 100 per layer. The joined arm therefore has up to 200 features.
- Clustering:
  - Features come from LC against recovered in the cohorts used.
  - If there is no recovered group, use the 200 most variable features, and say so.

## Classification models

Hyperparameters are tuned by an inner stratified 3 fold search on the training fold, scored by AUROC. The grids
below are fixed here.

| model | fixed settings | grid |
|---|---|---|
| logistic regression | L2, lbfgs, max_iter 5000 | C 0.01, 0.1, 1, 10 |
| elastic net | logistic, saga, max_iter 5000 | C 0.01, 0.1, 1; l1_ratio 0.2, 0.5, 0.8 |
| linear SVM | Platt probabilities | C 0.01, 0.1, 1, 10 |
| RBF SVM | Platt probabilities | C 0.1, 1, 10; gamma scale, 0.001, 0.01 |
| k nearest neighbours | | k 3, 5, 11, 21; weights uniform, distance |
| Gaussian naive Bayes | | var_smoothing 1e-9, 1e-6, 1e-3 |
| random forest | 500 trees | max_features sqrt, 0.3; min_samples_leaf 1, 5 |
| extra trees | 500 trees | max_features sqrt, 0.3; min_samples_leaf 1, 5 |
| histogram gradient boosting | | learning_rate 0.05, 0.1; max_depth 3, None |
| XGBoost | 300 trees | max_depth 2, 4; learning_rate 0.05, 0.1 |
| LightGBM | 300 trees, min_child_samples 5 | num_leaves 7, 15; learning_rate 0.05, 0.1 |
| CatBoost | 300 iterations | depth 4, 6; learning_rate 0.05, 0.1 |
| multilayer perceptron | max_iter 1000 | hidden layers (64), (64, 32); alpha 1e-4, 1e-2 |

A grid point that fails inside the inner search is skipped and recorded. If every point fails, the model fails.

## Foundation models

Foundation models are used as released, with default settings and no tuning:

- TabPFN v2
- TabPFN 2.5
- TabICLv2
- TabDPT
- Mitra

For gene expression there are two more. Each gives embeddings, followed by logistic regression tuned as above:

- BulkFormer (smallest model first)
- BulkRNABert

Before use, check that each one exists and installs. Use the version found and record it in the task's README and
in `Round2/environment.txt`. Tokens and licences are read from the environment and never printed or written.

## Failures

Every model is wrapped. The model fails, for that task and arm, when it:

- cannot be installed;
- needs a licence or token that is missing;
- runs out of memory;
- has no usable GPU;
- exceeds 15 minutes per fit.

For each failure, write a row to `failures.csv` with columns:

| column | content |
|---|---|
| task | |
| arm | |
| model | |
| reason | one of: install, licence or token, memory, GPU, time, error |
| message | |
| device | |
| memory | |

Then continue. Never drop a model silently.

## Classification metrics

Metrics are computed on each held out fold:

- AUROC (primary).
- AUPRC.
- Balanced accuracy, F1, MCC, sensitivity and specificity, at threshold 0.5.
- Brier score.
- Log loss, with probabilities clipped to 1e-15.
- Expected calibration error, with 10 equal width bins.

## Comparison

- Arm pairs: for every model and metric, the paired difference between arms on the same folds. This covers
  every pair of horizontal arms and every pair of vertical arms.
- Sign: the difference is oriented so that a positive value is better. Lower is better for Brier score, log
  loss and expected calibration error.
- Interval: a 95 percent interval from the corrected resampled t test (Nadeau and Bengio), over the 50 folds:
  - variance = (1/50 + n_test/n_train) × the variance of the fold differences;
  - t distribution with 49 degrees of freedom.
- Wording:
  - "Improves" means the interval is above zero.
  - "Gets worse" means the interval is below zero.
  - "Does not change" means the interval includes zero.
- Report: per metric, how many models improve, do not change, or get worse.

## Folder layout of every step

Each step has its own folder containing:

- `README.md`: purpose, input, method, output, and the result in plain words.
- `inputs.csv`: the files read, with SHA-256 checksums.
- `output/`.
- One script.

Data files are never committed.

## Amendment 1 (20 November 2026)

Added on 20 November 2026, before any Round2 experiment. The text above is unchanged; where they differ, this
amendment applies.

1. Buggert2025 is split by country into Buggert2025_UK and Buggert2025_SE, because country coincides with LC status
   in the full cohort. The full cohort is flagged and not used.
2. Subtyping: the anchor and the partner need at least 10 recovered participants, so that the same feature selection
   applies. Cohorts with cases only are listed, not used. An anchor with fewer than 50 LC cases is run but marked
   "too small to interpret".
3. Anchor sensitivity: when the three largest cohorts of a layer differ by fewer than five participants, steps 01
   and 02 are repeated with each of them as anchor.
4. Clinical variables count as a layer only when:
   - at least three variables are recorded for at least 80 percent of participants;
   - none of them defines the label;
   - for early detection, only variables known at the acute stage are used.
5. For pairs measured on different platforms or scales, B_naive is run and reported, and the README states that only
   B_scaled is interpretable.
6. One more model: the Round1 TACO ensemble (random forest, extra trees, SVC and MLP with soft voting).
7. Regression tasks (used for the age control):
   - models: the regression form of every model;
   - feature selection: absolute Spearman correlation with the outcome on the training fold, top 200;
   - metrics: mean absolute error (primary), root mean squared error, R squared and Spearman correlation;
   - "improves" means the paired difference in mean absolute error favours the combined arm, with its interval
     excluding zero.
8. Run control:
   - Every step can resume after interruption.
   - Before a full run, time one fold of every model and print the estimated total.
   - If the estimate exceeds 12 hours, stop and report the estimate per model.
9. New folders: Round2/0_Controls (known answer tasks) and Round2/4_Signature (signature discovery).

## Amendment 2 (20 November 2026)

Added on 20 November 2026, after the step 01 timing pass of Round2/0_Controls and before any full run. The text above
and Amendment 1 are unchanged; where they differ, this amendment applies.

1. Mitra is removed from the model grid on this hardware.
   - Reason: about 10 to 12 minutes per fit on CPU with no usable GPU. That is about 28 of the 32 estimated hours of
     step 01 alone.
   - It is recorded in failures.csv for every task and arm with:
     - the reason "time budget exceeded on this hardware";
     - the measured seconds per fit;
     - the device;
     - the estimated total.
   - The Mitra timing rows are kept in timing.csv.
2. After this removal, the 12 hour limit applies per step. If a step still exceeds it, the step stops and the
   estimate is reported per model, as in Amendment 1, item 8.
3. Every summary README lists the models that were not run, and why, next to the results.

## Amendment 3 (20 November 2026)

Added on 20 November 2026, to keep each step within the time limit. The text above and Amendments 1 and 2 are
unchanged; where they differ, this amendment applies.

1. Cross validation is 5 folds and 3 repeats for every task, arm and model. These are the first 3 of the 10 seeded
   repeats. For steps already run with 10 repeats, the raw files are kept and every summary is computed on the first
   3 repeats, so that every step uses the same folds.
2. CatBoost and histogram gradient boosting use a reduced search: at most 4 settings each and at most 300
   iterations.
3. The TabPFN 2.5 jobs that failed with a licence error before the licence was accepted are rerun. Their rows are
   removed from failures.csv once they succeed.
4. Before each step, the estimated time is printed as before. The 12 hour limit per step stays.

Every summary README states that its results rest on 15 folds.

## Amendment 5 (21 November 2026)

Added on 21 November 2026, before any run under it. The text above and the earlier amendments are unchanged; where
they differ, this amendment applies.

### A. New foundation models

TabPFN 3, TabPFN 3.5 and TabPFN 3.5 Fast are added for every task.
- Each is pinned by its explicit version, never the package default.
- The package version and the checkpoint name are recorded.
- The versions that the existing tabpfn_v2 and tabpfn_2_5 jobs actually loaded are checked and stated.
- If a licence has to be accepted, the work stops and reports exactly which one.
- The time limits still apply. A version too slow on this hardware is written to failures.csv with the measured time.

### B. Clinical variables as predictors of Long COVID

This replaces the clinical layer rule in Amendment 1, item 4.

1. Eligible variables: anything recorded about the participant that does not define the label. This includes:
   - demographics (age, sex, BMI, ethnicity);
   - acute illness (severity, hospitalisation, intensive care, WHO scale);
   - comorbidities;
   - vaccination before infection;
   - infection wave or variant;
   - smoking.
2. Excluded, because they define the label or come after it:
   - symptoms, symptom counts and scores;
   - questionnaire items used to define LC;
   - LC subgroups or trajectories;
   - treatment given for LC;
   - outcome measures used as the case definition (for example DLCO in GSE228320).

   Time since infection is not used as a predictor in the main analysis, because recruitment design often sets it.
   It is run as a separate sensitivity analysis and marked as such.
3. Early detection: only variables known at the acute stage.
4. Missing values and encoding:
   - A variable is kept in a cohort if it is recorded for at least 50 percent of participants.
   - Missing values are filled with the training fold median (or most frequent category), and a missing indicator is
     added.
   - Categorical variables are one hot encoded inside the training fold.
   - Severity is mapped to the WHO ordinal scale where possible, and the mapping is recorded.
5. Feature selection does not apply to the clinical layer, which keeps all eligible variables. Omics keep the protocol
   rule.
6. Every variable decision is written to Round2/5_Clinical/variables.csv: kept, or excluded and why, and the missing
   rate per cohort.

### C. New folder: Round2/5_Clinical

Run for detection and early detection, with the same seed, 15 folds, models and metrics as the protocol, plus the new
TabPFN versions.

- **01_univariate.** For each cohort and each eligible variable: the AUROC of the variable alone, and the odds ratio
  with its 95 percent interval. Where a variable is recorded in at least three cohorts, a random effects pooled
  estimate and heterogeneity (I squared) are added. It answers: which clinical variables separate LC from recovered,
  and is the effect consistent across cohorts?
- **02_clinical_alone.** Every model on all eligible clinical variables, per cohort with at least 10 per class. Also
  run with age and sex only.
- **03_horizontal_one_plus_one.** The clinical model on the anchor (the largest cohort with clinical data) plus a
  partner, on shared variables, with arms A_shared, B_naive and B_scaled. Each person is counted once: Cervia SomaScan
  and mass spectrometry are the same people, and so are the Su2022 layers.
- **04_vertical.** In every cohort with both omics and clinical data, on the same folds: omics alone, clinical alone,
  omics plus clinical joined, and late fusion. The key differences are:
  - omics plus clinical minus clinical alone (does omics add anything beyond clinical data?);
  - omics plus clinical minus omics alone.
- **05_summary.** Tables for every model and metric, the number of models that improve per comparison, failures.csv,
  and a README that answers in plain words:
  - which clinical variables predict Long COVID;
  - does a clinical model beat the omics models;
  - does omics add value beyond clinical data;
  - does combining cohorts help the clinical model?

### D. GSE215865 (early detection)

If the full clinical table from Synapse syn35874390 is in data/clinical/GSE215865, it is used after removing any
variable that defines the label. If it is not there, the report states that the clinical layer was tested with sex
only.

## Amendment 6 (21 November 2026)

Added on 21 November 2026, before any model of Round2/5_Clinical was run. Amendment 4 is reserved for the size
matched controls, which are still to come. The text above and the earlier amendments are unchanged; where they
differ, this amendment applies.

In 5_Clinical step 03 (horizontal, one plus one), the partner is chosen as follows:
- the unit that shares the most kept variables with the anchor;
- it must share at least age and sex;
- ties go to the larger unit.

If no unit qualifies, this is reported and step 03 is skipped for that task.
