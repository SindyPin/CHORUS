# 01_one_cohort

A_all: the anchor alone with all its features, for every model, by control and layer. Script: `run.py`.
**Results rest on 15 folds:** 5 folds by the first 3 of the 10 seeded repeats (Amendment 3).

- This step ran with 10 repeats before Amendment 3. Its raw files keep up to 50 folds, but every summary
  (`../04_summary`) reads folds 0 to 14 only.
- Estimate before the run, for 15 folds of the models run: 1.0 h (budget 12 h per step).

## Outputs (`output/`)

| file | content |
|---|---|
| `candidates.csv` | every cohort recording sex or age: participants (one sample each), women and men or age range, exclusion, eligibility |
| `control_anchors.csv` | anchor and partner per control and layer (the protocol rule), sensitivity anchors (none), same scale (none) |
| `cohorts_used.csv` | anchors and partners: participants, features, shared features |
| `raw/<job>.csv`, `raw/<job>_pred.csv` | metrics and held out predictions per fold |
| `timing.csv` | fold 0 seconds and the estimate per job; Mitra rows kept, marked not run (Amendment 2) |
| `failures.csv` | failures of this step |
| `mitra_timed/` | the six Mitra fold 0 fits that measured its time, kept out of `raw/` so no summary counts them |

`inputs.csv` lists every file read, with its SHA-256.

## Models not run, and why

| model | reason |
|---|---|
| Mitra | Amendment 2: time budget exceeded on this hardware. It took 510 to 769 s per fit on the CPU (fold 0 here). |
| TabPFN 2.5, age | licence or token: no TabPFN API key on this machine to download the regressor weights; retried, and it failed again |
| Gaussian naive Bayes, age | no regression form |
| BulkRNABert | time: 10.2 s per sample on the CPU, so 81 to 86 min per embedding, over the 15 minute limit |
