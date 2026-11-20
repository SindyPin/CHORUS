# 04_summary

Summary of steps 01 to 03. Script: `run.py`. The plain answers per control and layer are in `../README.md`.

**Results rest on 15 folds:** 5 folds by the first 3 of the 10 seeded repeats (Amendment 3), for every task, arm and
model. Step 01 jobs that ran more folds before Amendment 3 are read on folds 0 to 14 only.

## Outputs (`output/`)

| file | content |
|---|---|
| `metrics_by_arm.csv` | every control, layer, arm and model: mean and SD of every metric over the 15 folds |
| `paired_differences.csv` | every model, comparison and metric: paired difference on the same folds (positive = better), 95% interval from the corrected resampled t test (15 folds, test-to-train ratio of the anchor's split), verdict |
| `counts.csv` | per comparison and metric: models that improve, do not change, get worse, are not evaluable |
| `answers.csv` | the primary metric (AUROC for sex, MAE for age) for the comparisons that answer the question: B_naive and B_scaled vs A_shared, and joined and late fusion vs each layer |
| `failures.csv` | every failure of steps 01 to 03 |
| `models_not_run.csv` | per model and reason: number of jobs and where |

Metrics:
- **sex:** AUROC, AUPRC, balanced accuracy, F1, MCC, sensitivity, specificity, Brier score, log loss and ECE (10
  bins);
- **age:** MAE, RMSE, R squared and Spearman correlation.

Comparisons:
- **horizontal:** A_shared vs A_all; B_naive, B_scaled vs A_shared; B_scaled vs B_naive; B_naive, B_scaled vs A_all;
- **vertical:** joined and late fusion vs each layer; late fusion vs joined; layer 2 vs layer 1.

`inputs.csv` lists the step outputs read, with their SHA-256.

## Models not run, and why (Amendment 2, item 3)

| model | jobs | reason |
|---|---|---|
| Mitra | 40: every task and arm | Amendment 2: time budget exceeded on this hardware. It took 510 to 769 s per fit on the CPU, with no usable GPU. |
| TabPFN 2.5 | 18 age jobs, plus 2 age late fusion jobs | licence or token. There is no TabPFN API key on this machine to download the regressor weights; retried under Amendment 3, and it failed again. |
| Gaussian naive Bayes | 20 age jobs | no regression form |
| BulkRNABert | 8 RNA jobs | time: 10.2 s per sample on the CPU, so 81 to 86 min per embedding, over the 15 minute limit |

Models per comparison:
- sex: 18 (19 on RNA, with BulkFormer);
- age: 16 (17 on RNA).
