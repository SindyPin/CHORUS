# Paired cohorts: several layers on the same people (Test 2)

This file records which cohorts have two or more omics layers measured on the same samples, how the samples were
matched, and what could not be matched. It does not interpret anything. The analyses are run by
`scripts/paired_layers.py` (module `chorus.eval.paired`) and described in `results/README.md` ("Test 2").

**Matching rule.** A sample is paired when the same `sample_id` is in every layer's matrix and its labels agree
across the layers' metadata. The runner checks the agreement and stops if the labels disagree. For each contrast, one
sample per subject is kept: the earliest by timepoint.

## Summary

| cohort | layers | subjects with both layers | LC vs COVID_recovered (paired) | acute vs never_infected (paired) | status |
|---|---|---|---|---|---|
| Cervia2024 | SomaScan 7k (serum) + DIA-MS (serum) | 113 | 40 / 73 | none: the MS layer has no acute or healthy sample | run (LC only) |
| Su2022 | Olink Target 96 (plasma) + untargeted metabolomics (plasma) | 225 | 68 / 45 | 198 / 27 | run (both contrasts) |
| Kovarik2023 | MS proteomics + eicosanoids + Biocrates + triglycerides | none matched | none | none (no acute samples) | not built: no per-sample key (below) |

## Cervia2024: SomaScan plus MS on the same sera

- **Files:** `Cervia2024_SomaScan` and `Cervia2024_MS`, built from the same Mendeley deposit. The builder gives the
  same serum sample the same `sample_id` in both layers (`docs/proteomics_cohorts.md`).
- **Paired:** 113 six-month samples of 113 subjects (40 LC, 73 COVID_recovered). Every SomaScan six-month sample
  has an MS run.
- **Not paired:**
  - 6 MS samples have no SomaScan sample and no label: COV_ZH_181.3, COV_ZH_185.5, COV_ZH_197.3, COV_ZH_198.3,
    COV_ZH_28.2 and COV_ZH_30.2.
  - The 113 acute and 34 healthy SomaScan samples have no MS run (the MS deposit covers the six-month timepoint
    only), so the acute versus never_infected control cannot be run in this cohort.
- **Shared feature space:** both layers are keyed by UniProt, so the LC effect can be compared protein by protein
  (`effect_correlation.csv`).

## Su2022: Olink plus plasma metabolomics on the same draws

- **Files:** `Su2022` (Table S2.1) and `Su2022_metabolomics` (Table S2.2, `scripts/proteomics/build_Su2022_metabolomics.py`).
  - **Sample ids:** both sheets use the same rule: `HC-<id>` for healthy donors, the INCOV blood draw id for patients.
  - **Metadata:** both layers take their metadata from the same function and clinical tables
    (`build_Su2022.sample_metadata`).
  - **Values:** the metabolomics values are the authors' batch corrected, age, sex and BMI adjusted Z scores, like the
    Olink values. They are keyed by metabolite name as deposited, 891 metabolites, with no missing values.
- **Draws in each layer:**

| draws | Olink | metabolomics | both | Olink only | metabolomics only |
|---|---|---|---|---|---|
| acute (T1, T2) | 390 | 380 | 378 | 12 | 2 |
| healthy | 289 | 178 | 27 | 262 | 151 |
| post-acute LC | 75 | 68 | 68 | 7 | 0 |
| post-acute COVID_recovered | 50 | 45 | 45 | 5 | 0 |
| post-acute, no label | 2 | 1 | 1 | 1 | 0 |
| unknown phase (PD1, PD2) | 11 | 11 | 10 | 1 | 1 |

- **Paired subjects:** 225 subjects have at least one draw in both layers.
  - The LC contrast uses the 68 + 45 T3 draws.
  - The acute contrast uses each patient's earliest paired acute draw (198 patients) against the 27 healthy donors in
    both sheets.
- **Not matched:**
  - 151 healthy donors are in the metabolomics sheet only, and 262 in the Olink sheet only. The deposit does not say
    why the two healthy sets differ.
  - The metabolomics-only patient draws are INCOV203-T1, INCOV203-T2 and INCOV139-PD2.
- **No shared feature space:** proteins and metabolites, so there is no effect correlation between layers. Only the
  correlation of the out-of-fold probabilities is reported.

## Kovarik2023: proteomics plus eicosanoids, Biocrates and triglycerides (not built)

**Downloaded.** The supplement (Tables S1 to S4, `mmc2.xlsx` to `mmc5.xlsx`) was downloaded and checked
(`docs/provenance.md`). Each table has one column per sample, named `healthy_<n>`, `recovered_<n>` and
`long-covid_<n>`, n = 1 to 13:

| table | content | analytes | values |
|---|---|---|---|
| S1 | proteomics | 278 proteins | log2 LFQ |
| S2 | eicosanoids | 97 | nAUC |
| S3 | Biocrates MxP Quant 500 | 497 metabolites, 229 of them triglycerides | log2 µM |
| S4 | triglyceride heatmap | 126 triglycerides | "normalized values" |

**Linking S1 to the proteomics runs.** The proteomics builder links each S1 column to a MaxQuant run by exact equality
of the LFQ values (`table_s1_label` in `data/meta/Kovarik2023_plasma_meta.csv`). For S2, S3 and S4 the only per-sample
key is the column label. The paper does not say whether `healthy_6` is the same person in every table.

**Check of the labels.** S3 and S4 share 126 triglyceride species, so the labels can be tested on them:
- Under the shared labels, the per-lipid Spearman correlation across the 39 samples has a median of 0.10.
- When each S4 column is compared with every S3 column (values standardised per lipid), the best match is the same
  label for 1 of 39 columns, and the best correlations lie between 0.12 and 0.29. No column order of S3 reproduces S4.
- S4 is standardised per lipid (mean 0, SD 1 for every lipid), not "divided by the average concentration" as its
  caption says.

So the labels do not link two of the supplement's own tables, and they cannot be relied on to link S2 and S3 to the
proteomics samples. The Kovarik2023 paired set is not built. The analysis would need the per-sample key of the
eicosanoid and Biocrates measurements from the authors.

## Other candidates not built here

From `scripts/proteomics/other_layers.py`:
- **Su2022:** PBMC scRNA-seq and CITE-seq (ArrayExpress E-MTAB-10129) are a cell-level layer of other blood material and
  were not downloaded.
- **Cervia2024:** monocyte scRNA-seq, PBMC flow cytometry and VirScan cover subsets and were not taken. An Olink
  Inflammation layer (Ceglarek 2025) is not in the deposit.
- **Buggert2025, Captur2022, Chowdhury2024, Zoodsma2022:** no other public omics layer on the same people (flow
  cytometry, serology, ELISA or none).
