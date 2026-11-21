# GSE215865: labels from the post COVID checklist, compared with TACO

**Amendment 7 (21 November 2026): the TACO labels are withdrawn.** Round2 now uses the checklist labels built here for
GSE215865 and GSE225349:
- **main:** the paper's definition below (`round2_main`);
- **sensitivity:** symptoms only, the `_Ever` symptom items without health, quality of life or mood
  (`round2_symptoms_only`), reported separately on its own folds;
- subjects who never answered are unknown and excluded under both.

| GEO subjects (508) | main | symptoms only |
|---|---|---|
| LC (COVID positive, PASC) | 159 | 156 |
| recovered (COVID positive, no PASC) | 29 | 32 |
| never infected (COVID negative in GEO) | 71 | 71 |
| unknown (never answered) | 249 | 249 |

The three subjects that differ are PASC under the main definition only through health, quality of life or mood.
Round1 metadata and results are unchanged; the Round1 results that used the TACO labels are listed in
`affected_round1.csv` (not rerun).

## Source and definition

**Source:** the Mount Sinai COVID-19 Biobank clinical table (Synapse syn35876735; Thompson et al., Nat Med
2023;29:236-246). Provenance is in `data/meta/label_sources.md`. The table stays in `data/clinical/GSE215865/` and is
not committed.

**Definition**, from the paper's Methods ("Clinical data generation") and Table 1:
- Checklist items are of two kinds:
  - symptoms (true = present);
  - clinical deterioration: health, quality of life and daily mood (true = "Worse").
- Answers are collapsed over every checklist a person completed: an item is true if it was true at any time. The
  table's `*_Ever` columns are these collapsed answers.
- PASC means at least one item is true. No PASC means every answered item is false. No answer means unknown.

**Variant.** The paper removed items reported by 5 or fewer RNA-seq participants before modelling. No item is that
rare here, so the variant changes no label.

## Results (`output/summary.csv`, `output/comparison_counts.csv`)

### Who answered

| group | answered | PASC | no PASC |
|---|---|---|---|
| Synapse table (767 subjects) | 251 | 212 | 39 |
| paper, Table 1 | 232 | 195 | 37 |
| paper, Table 1, RNA-seq | 165 | 140 | 25 |
| GEO, RNA-seq (508 subjects) | 190 | 161 | 29 |

The Synapse release has 19 more respondents than the paper's Table 1.

### GEO subjects: TACO against the checklist

| category | subjects |
|---|---|
| agree: LC (TACO 1, PASC) | 131 |
| agree: no LC (TACO 0, no PASC) | 9 |
| disagree: TACO 1, checklist no PASC | 18 |
| disagree: TACO 0, checklist PASC | 25 |
| disagree: TACO 0, checklist PASC, COVID negative in GEO | 1 |
| TACO 1, no checklist answer | 154 |
| TACO 0, no checklist answer | 151 |
| no TACO label | 19 |

TACO labels 305 GEO subjects who never answered the checklist. So TACO's `Long_COVID` is not the checklist definition
of the paper.

### Early detection participants under the withdrawn TACO labels (424: acute sample, LC or recovered)

| TACO label | checklist PASC | checklist no PASC | no checklist answer |
|---|---|---|---|
| LC (302) | 130 | 18 | 154 |
| recovered (122) | 25 | 9 | 88 |

The 43 disagreements and every subject's labels are in `data/clinical/GSE215865/GSE215865_checklist_labels.csv`
(participant level, not committed).

## How TACO built its labels (`taco_rule.py`, `output/taco_rule_comparison.csv`, `output/taco_never_answered.csv`)

The TACO notebook `ML_Models/DCE_CV_Abation.ipynb` (github.com/SindyPin/TACO, commit cac2afb) reads this same table
(identical content) and labels each row 1 if any daily `Post_COVID19_Symptom_*` is True, or
`SARSCoV2_Recovery_Status` is "Not_Recovered", or `Post_COVID19_Health_Ever` is "Worse". It then keeps each subject's
first row, which is the earliest day. That rule reproduces all 489 TACO labels; without the daily symptoms it still
does, so in practice the label is the recovery status (most recent PCR positive) on the earliest day, plus 30
subjects with health "Worse". The 154 TACO LC and 151 TACO non LC subjects who never answered the checklist were
labelled by that status alone.

| rule (489 TACO subjects) | agree | disagree | TACO labels, rule cannot |
|---|---|---|---|
| a: `_Ever` symptom items | 137 | 47 | 305 |
| b: a + health, quality of life or mood "Worse" (paper) | 140 | 44 | 305 |
| c: daily `Post_COVID19_*`, all rows, missing as false | 291 | 198 | 0 |
| d: a, b or c + symptoms at presentation | 285 to 374 | 98 to 115 | 0 to 106 |
| e: the TACO rule (first row) | 489 | 0 | 0 |

## Scripts

`compare_labels.py` (labels and comparison) and `taco_rule.py` (the TACO rule and the variants), run from `Code/`.
`taco_rule.py` takes the path of `Data/411_DCEG_Paths_average.csv` at commit cac2afb (its first two columns).
