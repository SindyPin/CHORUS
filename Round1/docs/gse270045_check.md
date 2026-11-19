# GSE270045 check: why the held-out AUROC is inverted

**Verdict:** the labels and sample ids are correct, so the metadata is unchanged. GEO "Long Covid" maps to LC and GEO
"Healthy Control" to COVID_recovered, as the paper states. Sex inferred from expression matches the paper's sex counts in
both groups. Of the 64 genes whose direction the paper names, 63 have the same sign in our matrix. The inversion is
therefore not a labelling error. Two properties of the cohort cause it:
1. The genes that separate LC from recovered in the training cohorts run the opposite way in GSE270045.
2. Group is completely confounded with sample collection, so a technical difference cannot be separated from biology.

The cohort stays in the primary config with this caveat. The final results should include a sensitivity run that
leaves it out.

Checked on 2026-11-10 against `data/counts/GSE270045/GSE270045_LC_counts.tsv.gz` (author matrix, gene symbols) and
Shahbaz et al. 2025, Cell Rep Med (PMID 41205594, PMC12711683). The computation is `scripts/check_gse270045.py`.

## 1. GEO title to lc_status, with sex from expression

GEO has no sex field for this series, and XIST is not in the author matrix. Sex was inferred from RPS4Y1, which is
39–166 CPM in males and at most 3.3 CPM in females. KDM5D and EIF1AY were not used because their X-chromosome
counterparts give background in females.

| GSM | GEO title | sample id (matrix column) | lc_status | RPS4Y1 CPM | sex (RPS4Y1) | library (M) |
|---|---|---|---|---|---|---|
| GSM8333270 | Healthy Control 1 | CCI031 | COVID_recovered | 0.0 | female | 2.38 |
| GSM8333271 | Healthy Control 2 | CCI041 | COVID_recovered | 63.5 | male | 4.25 |
| GSM8333272 | Healthy Control 3 | CCI042 | COVID_recovered | 84.8 | male | 2.30 |
| GSM8333273 | Healthy Control 4 | CCI050 | COVID_recovered | 0.0 | female | 4.13 |
| GSM8333274 | Healthy Control 5 | CCI064 | COVID_recovered | 0.0 | female | 2.72 |
| GSM8333275 | Healthy Control 6 | CCI069 | COVID_recovered | 49.3 | male | 3.42 |
| GSM8333276 | Healthy Control 7 | CCI077 | COVID_recovered | 0.0 | female | 3.39 |
| GSM8333277 | Healthy Control 8 | CCI088 | COVID_recovered | 0.0 | female | 2.82 |
| GSM8333278 | Healthy Control 9 | CCI12 | COVID_recovered | 1.2 | female | 2.43 |
| GSM8333279 | Healthy Control 10 | CCI136 | COVID_recovered | 70.0 | male | 2.23 |
| GSM8333280 | Healthy Control 11 | CCI138 | COVID_recovered | 0.0 | female | 2.99 |
| GSM8333281 | Healthy Control 12 | CCI58 | COVID_recovered | 0.3 | female | 3.25 |
| GSM8333282 | Healthy Control 13 | CCI93 | COVID_recovered | 71.1 | male | 2.87 |
| GSM8333283 | Healthy Control 14 | HP1193 | COVID_recovered | 0.3 | female | 3.48 |
| GSM8333284 | Healthy Control 15 | HP1194 | COVID_recovered | 0.3 | female | 3.64 |
| GSM8333285 | Healthy Control 16 | HP1301 | COVID_recovered | 0.0 | female | 3.79 |
| GSM8333286 | Healthy Control 17 | HP851 | COVID_recovered | 0.0 | female | 4.91 |
| GSM8333287 | Long Covid 1 | S013_11 | LC | 38.8 | male | 7.18 |
| GSM8333288 | Long Covid 2 | S014_11 | LC | 0.0 | female | 4.42 |
| GSM8333289 | Long Covid 3 | S016_11 | LC | 0.1 | female | 7.13 |
| GSM8333290 | Long Covid 4 | S020_11 | LC | 0.3 | female | 3.99 |
| GSM8333291 | Long Covid 5 | S023_11 | LC | 0.0 | female | 3.69 |
| GSM8333292 | Long Covid 6 | S024_11 | LC | 0.0 | female | 7.23 |
| GSM8333293 | Long Covid 7 | S025_11 | LC | 0.0 | female | 2.60 |
| GSM8333294 | Long Covid 8 | S026_11 | LC | 0.0 | female | 4.65 |
| GSM8333295 | Long Covid 9 | S027_11 | LC | 0.2 | female | 4.92 |
| GSM8333296 | Long Covid 10 | S029_11 | LC | 0.0 | female | 2.16 |
| GSM8333297 | Long Covid 11 | S030_11 | LC | 3.3 | female | 0.91 |
| GSM8333298 | Long Covid 12 | S032_11 | LC | 0.0 | female | 0.90 |
| GSM8333299 | Long Covid 13 | S033_11 | LC | 166.2 | male | 2.02 |
| GSM8333300 | Long Covid 14 | S034_11 | LC | 0.0 | female | 0.65 |
| GSM8333301 | Long Covid 15 | S035_11 | LC | 1.5 | female | 0.68 |
| GSM8333302 | Long Covid 16 | S038_11 | LC | 0.0 | female | 1.42 |
| GSM8333303 | Long Covid 17 | S039_11 | LC | 1.4 | female | 1.43 |
| GSM8333304 | Long Covid 18 | S041_11 | LC | 0.0 | female | 1.96 |
| GSM8333305 | Long Covid 19 | S042_11 | LC | 161.7 | male | 8.68 |

| group | female | male | paper |
|---|---|---|---|
| COVID_recovered (17) | 12 | 5 | 12 RF + 5 RM: matches |
| LC (19) | 16 | 3 | 15 LCF + 3 LCM: matches. The one sample GEO has beyond the paper's 18 is female |

Every title maps to the expected label. The predicted sex split matches the paper's per-group counts, so the matrix
columns are attached to the right groups. Per-sample sex is still not recorded in GEO, so the schema `sex` field is
left empty: the inferred value is expression-based, not reported sex.

## 2. Direction against the paper

There are no DE tables in the paper. The supplement has two files: mmc1 holds supplementary figures and mmc2 is the
article. The paper reports directions only for named genes, and separately for female (LCF vs RF) and male (LCM vs RM)
patients. Hedges g (LC minus recovered, log2 CPM, 15,574 genes with CPM ≥ 1 in half the samples) was computed within
each sex using the RPS4Y1 calls.

- **LCF vs RF:** all 35 named genes agree in sign. Up: ZNF469, BRINP2, FEZF2, HOXC12, SKOR2, ESR1, ESR2, ESRRG, AR,
  NR3C1, CD69, CD302, CD226, CD200R1, MEIS2. Down: WDR62, MUC12, SYCE1, SYCE1L, CD3E, CD14, CD7, CD320, CD101, CD82,
  CD81, CD6, CD68, CXCR3, IL32, IL15RA, IL27, IL2RB, IL12RB1, IL11RA.
- **LCM vs RM:** 28 of 29 agree. The exception is ARTN, reported up and here g = −0.78, with only 3 vs 5 samples.

The top 50 genes by |g| over both sexes are 27 up and 23 down, and none is named in the paper. Within every one of
the 50, the female and male strata agree in sign:

CIC (-4.32), TAF1D (+4.11), USP33 (+4.09), MOB1A (+3.94), SIPA1 (-3.92), PHF20L1 (+3.91), HDAC5 (-3.85), MAP1S (-3.81), PRRC2A (-3.80), NUFIP2 (+3.78), PPP6R1 (-3.76), BLOC1S6 (+3.72), SCAF1 (-3.71), MAP7D1 (-3.69), RBM5 (+3.67), CYSLTR1 (+3.63), MIPOL1 (+3.59), GDAP2 (+3.56), OR4F17 (+3.56), NMNAT3 (+3.56), RAP1A (+3.56), SLC38A10 (-3.55), PLEKHM2 (-3.53), PTPRC (+3.52), FAM83A (+3.51), SH2D3C (-3.51), ABHD13 (+3.51), EPN1 (-3.46), TBC1D8B (+3.40), PGAP1 (+3.40), GRN (-3.39), STX7 (+3.35), CTSD (-3.34), MLLT1 (-3.34), MARCHF7 (+3.33), MKLN1 (+3.33), ARAF (-3.33), CHERP (-3.32), GAK (-3.30), LMBRD1 (+3.29), PIGN (+3.28), DDX6 (+3.27), CFAP54 (+3.26), ACTN4 (-3.26), ZNF710 (-3.25), ARHGAP45 (-3.25), ZNF33A (+3.21), PTGES3 (+3.21), ARHGEF2 (-3.19), ZFP36L2 (-3.19)


## 3. Why the transfer AUROC is inverted

- **The signal is unusually large and global.** 4,306 of 15,574 genes have |g| > 1 and 737 have |g| > 2. The top
  genes are broadly expressed housekeeping and regulatory genes (CIC, PRRC2A, HDAC5, PTPRC, CTSD), not a focused
  immune programme. PC2, which explains 12.5% of the variance, separates the groups: mean +24 for recovered and
  −22 for LC.
- **Group is confounded with collection.** All 19 LC samples belong to one id series (`S0xx_11`). All 17 controls
  belong to others (`CCI…`, `HP…`), so this looks like two separate collections or sequencing runs. All 5 LC samples
  with libraries under 1 M reads are in the S series. Batch and phenotype cannot be separated within this cohort, and
  per-study centring (center_scale) cannot remove a batch that coincides with the label.
- **The opposite direction is specific to the selected genes.** Over all 10,795 genes shared with GSE226260 whole
  blood, the GSE270045 effects correlate weakly and positively (r = 0.12). But among the 200 genes the training
  studies select when GSE270045 is held out, the sign agrees in only 30% (differential selector) to 34% (invariant)
  of genes, against every training study alike. With the variance selector, which ignores the label, agreement is
  50%. The genes that separate LC from recovered elsewhere move the other way here, so a correctly labelled cohort
  is scored upside down.
- **The phenotype differs.** LC here is post-COVID ME/CFS at least 12 months after infection (Canadian Consensus
  Criteria), against broader LC definitions in the training cohorts. This may add to the reversal but cannot be
  separated from the batch effect.

## Actions

- Metadata unchanged: labels and ids verified.
- Keep GSE270045 in the primary run and report a leave-it-out sensitivity run alongside it. When GSE270045 is held
  out, read its AUROC as a test of transfer to a confounded ME/CFS cohort, not of label quality.
- Ask the lead contact (S. Elahi, elahi@ualberta.ca) two questions: whether the LC and control groups were collected
  and sequenced in separate batches, and which LC sample is extra relative to the paper's 18. Per-sample batch would
  allow a batch covariate. Raw reads are public, so re-quantifying from SRA would at least remove differences in
  pipeline.
