# GSE215865: labels from the post COVID checklist, compared with TACO

**For review. The Round2 labels are not changed:** GSE215865 still uses the TACO labels.

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

### Round2 early detection participants (424: acute sample, LC or recovered)

| Round2 (TACO) label | checklist PASC | checklist no PASC | no checklist answer |
|---|---|---|---|
| LC (302) | 130 | 18 | 154 |
| recovered (122) | 25 | 9 | 88 |

The 43 disagreements and every subject's labels are in `data/clinical/GSE215865/GSE215865_checklist_labels.csv`
(participant level, not committed).

## Script

`compare_labels.py`, run from `Code/`.
