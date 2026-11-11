# Acute transfer experiment

**Question.** Does a Long COVID signature learned from post-acute blood transcriptomes also separate future LC from
recovered patients in blood drawn during the acute hospital stay?

**Config.** `configs/lc_acute_transfer.yaml`.

**Run.**

```
chorus run configs/lc_acute_transfer.yaml --within
```

**Design.**
- **Training.** The post-acute cohorts that are enabled in `configs/lc_transcriptomics_blood.yaml`: currently
  GSE226260_WB, GSE226260_PBMC, GSE224615, GSE251849 and GSE270045.
  - Each is marked `train_only: true`, so `run_loso` never holds them out and GSE215865 is the only held-out study.
  - When a cohort is enabled in the primary config, add it here as `train_only` as well.
- **Test.** GSE215865 (Thompson 2023, Mount Sinai): whole blood drawn during acute COVID-19 admissions.
  - Labels follow the step 4 rule in `scripts/build_cohort_meta.py`: 302 LC versus 122 COVID_recovered subjects.
  - The 64 COVID-negative subjects (`never_infected`) fall outside the contrast.
  - The 20 subjects without a usable label are `unknown` and are dropped.
- **One sample per subject.** `one_sample_per_subject` keeps each subject's earliest library.
  - GEO labels only 21 libraries `T0`, so "T0" here means the first available draw.
  - That is day 0 for 418 of 508 subjects and day 1 to 31 for the other 90, whose earlier libraries were not deposited.
  - Repeat libraries of the same draw (several plates) are reduced to the first one in file order.
- **Fitting.** Harmonisers, selectors and classifiers are fitted on the training cohorts only.
  - Harmonisers are then applied to GSE215865 in add-on mode, exactly as in the leave-one-study-out benchmark.
  - With `nested_select_harmoniser`, the harmoniser is chosen by an inner leave-one-study-out over the training cohorts.
- **Reference.** `--within` adds a within-GSE215865 cross-validated AUROC as the optimism reference. The gap between it
  and the transfer AUROC is the cost of crossing from post-acute to acute blood.

**Caveats to report with the results.**
- **Mixed tissue.** Training mixes whole blood and PBMC while the test set is whole blood. Rerun with
  `matrices: [whole_blood]`, which currently leaves GSE226260_WB and GSE270045 as training cohorts.
- **Different label timing.** Training labels describe the state at sampling; GSE215865 labels describe an outcome about
  a year after discharge (Thompson 2023). The LC definitions also differ by cohort.
- **Severity.** Every GSE215865 subject was hospitalised, while most training cohorts are mainly mild infections.
  Severity cannot be adjusted for because it is not in the GEO metadata.
- **No prior selectors.** The TACO 411 DCE genes and the MRCONTROL 32 causal genes were both derived from GSE215865 with its
  Long COVID labels: TACO trained on the Mount Sinai blood RNA-seq, and MRCONTROL used it for its RNA-seq step. Selecting
  features with them would carry information from the test labels into training, so this config uses only the variance,
  differential and invariant selectors. The prior selectors stay in the primary config, which does not include GSE215865.
- **Covariates.** The only covariate is `sample_matrix`, because age and sex are recorded only for GSE226260_PBMC.
