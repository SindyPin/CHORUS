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
