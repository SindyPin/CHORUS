# 03_vertical_one_plus_one

One layer plus one, on the same participants. Script: `run.py`.

| pair | layer 1 | layer 2 | participants |
|---|---|---|---|
| Su2022 | Olink | plasma metabolomics | 225 |
| Cervia2024 | SomaScan | MS proteomics | 113 |

Participants: those with the control recorded and a sample in both layers (the same sample id), one sample each (the
earliest).

| arm | features |
|---|---|
| layer_1 | layer 1 alone, 100 features |
| layer_2 | layer 2 alone, 100 features |
| joined | both layers, 100 features per layer selected separately (early fusion) |
| late_fusion | mean of the layer_1 and layer_2 predictions of the same model, fold by fold (no extra fit) |

**Results rest on 15 folds:** 5 folds by the first 3 of the 10 seeded repeats (Amendment 3).
- Su2022 values in both layers are the authors' Z scores, adjusted for age, sex and BMI.
- Estimate before the run: 0.4 h for the models run (budget 12 h per step).

## Outputs (`output/`)

| file | content |
|---|---|
| `cohorts_used.csv` | per control and pair: participants, features per layer, women and men or age range |
| `raw/` | metrics and held out predictions per job and fold |
| `timing.csv` | fold 0 timings and estimates |
| `failures.csv` | failures of this step |

`inputs.csv` lists every file read, with its SHA-256.

## Models not run, and why

| model | reason |
|---|---|
| Mitra | Amendment 2: time budget exceeded on this hardware, measured in step 01. It has no late fusion either. |
| TabPFN 2.5, age | licence or token: no TabPFN API key on this machine to download the regressor weights; retried, and it failed again. It has no late fusion either. |
| Gaussian naive Bayes, age | no regression form |

The expression models are not part of this step, because neither pair has RNA.
