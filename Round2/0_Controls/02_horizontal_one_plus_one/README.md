# 02_horizontal_one_plus_one

One cohort plus one, on the anchors and partners of step 01. Script: `run.py`. Every arm is tested on the anchor's
held out fold (the folds of step 01).

| arm | training set |
|---|---|
| A_shared | anchor training fold, features shared with the partner |
| B_naive | the same plus every partner participant, without adjustment |
| B_scaled | the same, after z scoring each feature within each cohort without labels |

**Results rest on 15 folds:** 5 folds by the first 3 of the 10 seeded repeats (Amendment 3).

- **B_naive is reported only.** No pair is measured with the same pipeline and reference, so B_scaled is the
  interpretable arm (Amendment 1, item 5).
- **Estimate before the run:** 2.6 h for the models run (budget 12 h per step); before the Amendment 5 run: 3.7 h.

## Outputs (`output/`)

| file | content |
|---|---|
| `pairs.csv` | per control and layer: participants, features, shared features, same scale, interpretable arms |
| `raw/` | metrics and held out predictions per job and fold |
| `timing.csv` | fold 0 seconds and the estimate per job; Mitra rows use the step 01 measurement, marked not run |
| `failures.csv` | failures of this step |

`inputs.csv` lists every file read, with its SHA-256.

## Models not run, and why

| model | reason |
|---|---|
| Mitra | Amendment 2: time budget exceeded on this hardware, measured in step 01 |
| Gaussian naive Bayes, age | no regression form |
| BulkRNABert | time: over the 15 minute limit per embedding |

TabPFN 2.5 age jobs: they failed on the licence before (no TabPFN API key on this machine); they ran under Amendment 5 once
the key was saved, and their failure rows were removed (Amendment 3, item 3).
TabPFN 3, 3.5 and 3.5 Fast (Amendment 5, A) ran on every job of this step, on the same 15 folds; no new failure.
