# Results

This file records how each result folder was produced. It does not interpret the results.

## Environment

- Every run below used the arm64 virtual environment `~/.venvs/chorus-arm64` (Homebrew Python 3.12.16), called as
  `arch -arm64 ~/.venvs/chorus-arm64/bin/...`.
- Package versions are pinned in `docs/requirements-chorus-arm64.txt`. Key versions: numpy 2.5.4, pandas 3.1.1,
  scikit-learn 1.9.1, torch 2.15.0, tabpfn 9.2.0, inmoose 0.9.1, xgboost 3.4.1.
- `docs/provenance.md` describes the environments and the GSE224615 re-quantification.

## TabPFN status

The TabPFN licence token file `~/.config/chorus/tabpfn_token` was absent when the real-data runs started. They ran with
`CHORUS_DISABLE_TABPFN=1`, so the `taco` classifier is the soft voting ensemble of **rf, extratrees, svc and mlp only,
without the TabPFN base learner**. The simulation runs do not use `taco`.

## Runs

Runtimes are wall-clock seconds from `results/runtimes.tsv`, all on the same laptop. The runs of 11 and 12 November were
made one after another. On 13 November `lc_transcriptomics_blood_within_v2` ran at the same time as the two morning
control runs (the label permutation and the sample_matrix control removed since), so those runtimes overlap. The later
13 November runs (from `infer_sex` on) ran one after another. `infer_sex` was run twice with identical output; the
runtime is that of the second run. The version 2 reruns of 13 November (from 11:19) ran one after another; see "Version 2
changes" for what they replaced. The version 3 reruns of 13 November (from 15:12) ran one after another; see "Version 3
changes". The runs of 15 November (from 06:49) ran one after another; see "Proteomics (15 November 2026)". The two
Test 2 runs (from 07:43) ran one after another; see "Test 2: layers on the same people (15 November 2026)". The Test 3
run (17:24) ran alone; see "Test 3: mosaic programmes (15 November 2026)". The metabolomics runs of 16 November ran one
after another (the all-cohort concordance was rerun into its final folder at 00:19:51); see "Metabolomics (16 November
2026)". The version 4 proteomics runs of 16 November (from 07:27) ran one after another while the PRJNA1149196 salmon
quantification ran in the background, so their runtimes overlap with it. The reporting reruns of 16 November (from
10:51) ran one after another, also during the quantification; see "Version 4 cohort set and reporting fixes (16 November
2026)". The version 4 transcriptomic runs of 17 November (from 00:36) ran one after another; the runtime of
`lc_transcriptomics_blood_v4` includes time the laptop was asleep. The six follow-up runs of 17 November (from
09:19) ran at the same time, so their runtimes overlap; see "Follow-up analyses (17 November 2026)". The version 5 runs of
17 November (from 10:17) ran one after another; see "Version 5 (17 November 2026)". The CORSIP2026 checks of
17 November (from 11:36) overlapped: the panel495 benchmark, then the permutation, random gene and concordance chunks in
parallel; see "CORSIP2026 panel transfer: verification". The further analyses of 17 and 18 November ran as listed (the
guided training units in 10 parallel workers); see "Further analyses: noise ceiling, guided training, full range".

| folder | command | config | runtime (s) |
|---|---|---|---|
| `lc_transcriptomics_blood_v1` | `chorus run configs/lc_transcriptomics_blood.yaml --within` (written to `lc_transcriptomics_blood`, copied) | configs/lc_transcriptomics_blood.yaml | 1053 (started 2026-11-12 09:19:36) |
| `lc_transcriptomics_blood` | the v1 run, then `chorus run configs/lc_transcriptomics_blood.yaml --only-harmonisers combat_addon,chorus,nested` (version 2), then `--only-harmonisers combat_addon,chorus_protected,nested` (version 3) | configs/lc_transcriptomics_blood.yaml | 397 for the version 2 rerun (started 2026-11-13 11:37:23), 404 for the version 3 rerun (started 15:12:49) |
| `lc_transcriptomics_blood_noGSE270045_v1` | `chorus run configs/lc_transcriptomics_blood_noGSE270045.yaml --within` (copied) | configs/lc_transcriptomics_blood_noGSE270045.yaml | 805 (started 2026-11-12 09:37:09) |
| `lc_transcriptomics_blood_noGSE270045` | the v1 run, then `chorus run configs/lc_transcriptomics_blood_noGSE270045.yaml --only-harmonisers combat_addon,chorus,nested` (version 2), then `--only-harmonisers combat_addon,chorus_protected,nested` (version 3) | configs/lc_transcriptomics_blood_noGSE270045.yaml | 300 for the version 2 rerun (started 2026-11-13 11:44:00), 312 for the version 3 rerun (started 15:19:33) |
| `lc_acute_transfer_v1` | `chorus run configs/lc_acute_transfer.yaml --within` (copied) | configs/lc_acute_transfer.yaml | 231 (started 2026-11-12 09:50:34) |
| `lc_acute_transfer` | the v1 run, then `chorus run configs/lc_acute_transfer.yaml --only-harmonisers combat_addon,chorus,nested` (version 2), then `--only-harmonisers combat_addon,chorus_protected,nested` (version 3) | configs/lc_acute_transfer.yaml | 101 for the version 2 rerun (started 2026-11-13 11:49:00), 105 for the version 3 rerun (started 15:24:45) |
| `simulation_factors_v1` | `scripts/simulation_benchmark.py --out results/simulation_factors` (folder renamed) | built into the script | 231 (started 2026-11-11 09:09:10) |
| `simulation_factors` | `scripts/simulation_benchmark.py --out results/simulation_factors` | built into the script | 475 (version 3, started 2026-11-13 15:30:07) |
| `simulation_shifts_only_v1` | `scripts/simulation_benchmark.py --out results/simulation_shifts_only --n-factors 0` (folder renamed) | built into the script | 346 (started 2026-11-11 09:13:14) |
| `simulation_shifts_only` | `scripts/simulation_benchmark.py --out results/simulation_shifts_only --n-factors 0` | built into the script | 466 (version 3, started 2026-11-13 15:38:02) |
| `lc_transcriptomics_blood_within_v2` | `chorus run configs/lc_transcriptomics_blood.yaml --within-only --out results/lc_transcriptomics_blood_within_v2` | configs/lc_transcriptomics_blood.yaml | 1320 (started 2026-11-13 07:05:33) |
| `control_label_permutation` | `scripts/control_runs.py permute configs/lc_transcriptomics_blood.yaml --out results/control_label_permutation --n-perm 20` | configs/lc_transcriptomics_blood.yaml, changed by the script | 178 (version 3, started 2026-11-13 15:46:30) |
| `control_sex_label/inferred_sex.csv` | `scripts/infer_sex.py configs/lc_transcriptomics_blood.yaml --copy results/control_sex_label/inferred_sex.csv` | configs/lc_transcriptomics_blood.yaml (raw matrices) | 16 (started 2026-11-13 09:08:07) |
| `control_sex_label_v1` | `scripts/control_runs.py sex configs/lc_transcriptomics_blood.yaml --out results/control_sex_label` (folder renamed) | configs/lc_transcriptomics_blood.yaml, changed by the script | 84 (started 2026-11-13 08:56:20) |
| `control_sex_label` | `scripts/control_runs.py sex configs/lc_transcriptomics_blood.yaml --out results/control_sex_label` | configs/lc_transcriptomics_blood.yaml, changed by the script | 112 (rerun with `forced_features`, started 2026-11-15 06:57:24; version 3: 113) |
| `lc_transcriptomics_blood_harmonisation_v1` | `chorus harmonisation configs/lc_transcriptomics_blood.yaml --out results/lc_transcriptomics_blood_harmonisation` (folder renamed) | configs/lc_transcriptomics_blood.yaml | 100 (started 2026-11-13 09:01:09) |
| `lc_transcriptomics_blood_harmonisation` | `chorus harmonisation configs/lc_transcriptomics_blood.yaml --out results/lc_transcriptomics_blood_harmonisation` | configs/lc_transcriptomics_blood.yaml | 104 (version 3, started 2026-11-13 15:26:30) |
| `concordance` | `scripts/concordance.py configs/lc_transcriptomics_blood.yaml --out results/concordance` | configs/lc_transcriptomics_blood.yaml | 48 (started 2026-11-13 09:03:45); reporting rerun 49 (started 2026-11-16 10:51:30) |
| `lc_proteomics_affinity` | `chorus run configs/lc_proteomics_affinity.yaml --within` | configs/lc_proteomics_affinity.yaml | 236 (started 2026-11-15 06:49:03) |
| `lc_proteomics_ms` | `chorus run configs/lc_proteomics_ms.yaml --within` | configs/lc_proteomics_ms.yaml | 226 (started 2026-11-15 06:52:59) |
| `control_acute_proteomics` | `chorus run configs/control_acute_proteomics.yaml --within` | configs/control_acute_proteomics.yaml | 23 (started 2026-11-15 06:56:45) |
| `lc_proteomics_affinity/harmonisation` | `chorus harmonisation configs/lc_proteomics_affinity.yaml --out results/lc_proteomics_affinity/harmonisation` | configs/lc_proteomics_affinity.yaml | 4 (started 2026-11-15 06:57:08) |
| `lc_proteomics_ms/harmonisation` | `chorus harmonisation configs/lc_proteomics_ms.yaml --out results/lc_proteomics_ms/harmonisation` | configs/lc_proteomics_ms.yaml | 1 (started 2026-11-15 06:57:12) |
| `lc_proteomics_affinity/concordance` | `scripts/concordance.py configs/lc_proteomics_affinity.yaml --out results/lc_proteomics_affinity/concordance --top 50` | configs/lc_proteomics_affinity.yaml | 2 (started 2026-11-15 06:57:13); reporting rerun 3 (started 2026-11-16 10:52:19) |
| `lc_proteomics_ms/concordance` | `scripts/concordance.py configs/lc_proteomics_ms.yaml --out results/lc_proteomics_ms/concordance --top 50` | configs/lc_proteomics_ms.yaml | 1 (started 2026-11-15 06:57:16); reporting rerun 1 (started 2026-11-16 10:52:22) |
| `concordance_proteomics_all` | `scripts/concordance.py configs/concordance_proteomics_all.yaml --out results/concordance_proteomics_all --top 50` | configs/concordance_proteomics_all.yaml | 7 (started 2026-11-15 06:57:17); reporting rerun 10 (started 2026-11-16 10:52:23) |
| `control_acute_transcriptomics` | no run: fewer than three cohorts qualify (`cohorts.csv`) | none | none |
| `paired_Cervia2024` | `scripts/paired_layers.py --cohort Cervia2024 --out results/paired_Cervia2024` | built into the script | 2059 (started 2026-11-15 07:43:05); `effect_correlation.csv` rewritten by `--pairs-only`, 2 (started 2026-11-16 10:53:01) |
| `paired_Su2022` | `scripts/paired_layers.py --cohort Su2022 --out results/paired_Su2022` | built into the script | 1029 (started 2026-11-15 08:17:24) |
| `mosaic_programmes` | `scripts/mosaic_programmes.py --out results/mosaic_programmes` | configs/mosaic_lc_rna.yaml, configs/concordance_proteomics_all.yaml, configs/mosaic_acute_rna.yaml, configs/mosaic_acute_protein.yaml | 239 (started 2026-11-15 17:24:32; acute 28, lc 209); `pairs_<set>.csv` written by `--contrasts acute,lc --pairs-only`, 1 (started 2026-11-16 10:52:59) |
| `lc_metabolomics` | `chorus run configs/lc_metabolomics.yaml --within` | configs/lc_metabolomics.yaml | 205 (started 2026-11-16 00:16:10) |
| `lc_metabolomics/harmonisation` | `chorus harmonisation configs/lc_metabolomics.yaml --out results/lc_metabolomics/harmonisation --n-pcs 10` | configs/lc_metabolomics.yaml | 1 (started 2026-11-16 00:19:35) |
| `lc_metabolomics/concordance` | `scripts/concordance.py configs/lc_metabolomics.yaml --out results/lc_metabolomics/concordance --top 20` | configs/lc_metabolomics.yaml | 1 (started 2026-11-16 00:19:36); reporting rerun 2 (started 2026-11-16 10:52:54) |
| `lc_metabolomics/concordance_all` | `scripts/concordance.py configs/concordance_metabolomics_all.yaml --out results/lc_metabolomics/concordance_all --top 50` | configs/concordance_metabolomics_all.yaml | 2 (started 2026-11-16 00:19:51); reporting rerun 3 (started 2026-11-16 10:52:56) |
| `lc_proteomics_affinity_v4` | `chorus run configs/lc_proteomics_affinity_v4.yaml --within` | configs/lc_proteomics_affinity_v4.yaml | 421 (started 2026-11-16 07:27:40) |
| `lc_proteomics_affinity_v4/harmonisation` | `chorus harmonisation configs/lc_proteomics_affinity_v4.yaml --out results/lc_proteomics_affinity_v4/harmonisation` | configs/lc_proteomics_affinity_v4.yaml | 3 (started 2026-11-16 07:34:58) |
| `lc_proteomics_affinity_v4/concordance` | `scripts/concordance.py configs/lc_proteomics_affinity_v4.yaml --out results/lc_proteomics_affinity_v4/concordance --top 50` | configs/lc_proteomics_affinity_v4.yaml | 2 (started 2026-11-16 07:35:01); reporting rerun 4 (started 10:52:33) |
| `lc_proteomics_ms_v4` | `chorus run configs/lc_proteomics_ms_v4.yaml --within` | configs/lc_proteomics_ms_v4.yaml | 403 (started 2026-11-16 07:40:23) |
| `lc_proteomics_ms_v4/harmonisation` | `chorus harmonisation configs/lc_proteomics_ms_v4.yaml --out results/lc_proteomics_ms_v4/harmonisation` | configs/lc_proteomics_ms_v4.yaml | 3 (started 2026-11-16 07:47:06) |
| `lc_proteomics_ms_v4/concordance` | `scripts/concordance.py configs/lc_proteomics_ms_v4.yaml --out results/lc_proteomics_ms_v4/concordance --top 50` | configs/lc_proteomics_ms_v4.yaml | 1 (started 2026-11-16 07:47:09); reporting rerun 1 (started 10:52:37) |
| `concordance_proteomics_all_v4` | `scripts/concordance.py configs/concordance_proteomics_all_v4.yaml --out results/concordance_proteomics_all_v4 --top 50` | configs/concordance_proteomics_all_v4.yaml | 14 (started 2026-11-16 07:47:10); reporting rerun 14 (started 10:52:39) |
| `lc_transcriptomics_blood_v4` | `chorus run configs/lc_transcriptomics_blood_v4.yaml --within` | configs/lc_transcriptomics_blood_v4.yaml | 21455, including laptop sleep (started 2026-11-17 00:36:55) |
| `lc_transcriptomics_blood_harmonisation_v4` | `chorus harmonisation configs/lc_transcriptomics_blood_v4.yaml --out results/lc_transcriptomics_blood_harmonisation_v4` | configs/lc_transcriptomics_blood_v4.yaml | 162 (started 2026-11-17 06:34:30) |
| `concordance_v4` | `scripts/concordance.py configs/lc_transcriptomics_blood_v4.yaml --out results/concordance_v4` | configs/lc_transcriptomics_blood_v4.yaml | 81 (started 2026-11-17 06:37:12) |
| `mosaic_programmes_v4` | `scripts/mosaic_programmes.py --out results/mosaic_programmes_v4 --contrasts lc_v4` | configs/mosaic_lc_rna_v4.yaml, configs/concordance_proteomics_all_v4.yaml | 167 (started 2026-11-17 06:38:33) |
| `lc_transcriptomics_blood_v4_noQCflag` | `chorus run configs/lc_transcriptomics_blood_v4_noQCflag.yaml --within` | configs/lc_transcriptomics_blood_v4_noQCflag.yaml | 3547 (started 2026-11-17 06:41:20) |
| `early_detection/vertical` | `scripts/paired_layers.py --cohort Su2022 --contrasts early --out results/early_detection/vertical` | built into the script | 436 (started 2026-11-17 09:19:31) |
| `early_detection/horizontal` | `scripts/early_detection.py horizontal configs/early_detection_horizontal.yaml --out results/early_detection/horizontal --n-perm 100` | configs/early_detection_horizontal.yaml | 742 (started 2026-11-17 09:21:55) |
| `early_detection/reference` | `scripts/early_detection.py reference configs/early_detection_reference.yaml --out results/early_detection/reference` | configs/early_detection_reference.yaml | 39 (started 2026-11-17 09:21:55) |
| `clustering_readiness` | `scripts/clustering_readiness.py --out results/clustering_readiness` | configs/lc_transcriptomics_blood_v4.yaml, configs/lc_proteomics_affinity_v4.yaml | 273 (started 2026-11-17 09:23:32) |
| `acute_persistence` | `scripts/acute_persistence.py --out results/acute_persistence` | configs/control_acute_proteomics.yaml, configs/concordance_proteomics_all_v4.yaml, configs/mosaic_acute_rna.yaml, configs/mosaic_acute_protein.yaml, configs/mosaic_lc_rna_v4.yaml | 187 (started 2026-11-17 09:25:23) |
| `concordance_vs_transfer` | `scripts/concordance_vs_transfer.py --out results/concordance_vs_transfer` | the loso_results.csv of lc_transcriptomics_blood_v4, lc_proteomics_affinity_v4, lc_proteomics_ms_v4, lc_metabolomics | 16 (started 2026-11-17 09:26:09) |
| `lc_proteomics_affinity_v5` | `chorus run configs/lc_proteomics_affinity_v5.yaml --within` | configs/lc_proteomics_affinity_v5.yaml | 466 (started 2026-11-17 10:17:01) |
| `lc_proteomics_explore_v5` | `chorus run configs/lc_proteomics_explore_v5.yaml --within` | configs/lc_proteomics_explore_v5.yaml | 236 (started 2026-11-17 10:24:47) |
| `concordance_proteomics_all_v5` | `scripts/concordance.py configs/concordance_proteomics_all_v5.yaml --out results/concordance_proteomics_all_v5 --top 50` | configs/concordance_proteomics_all_v5.yaml | 12 (started 2026-11-17 10:28:43) |
| `concordance_v5` | `scripts/concordance.py configs/concordance_rna_v5.yaml --out results/concordance_v5` | configs/concordance_rna_v5.yaml | 145 (started 2026-11-17 10:28:55) |
| `mosaic_programmes_v5` | `scripts/mosaic_programmes.py --out results/mosaic_programmes_v5 --contrasts lc_v5` | configs/mosaic_lc_rna_v4.yaml, configs/concordance_proteomics_all_v5.yaml | 174 (started 2026-11-17 10:31:20) |
| `panel_transfer_CORSIP2026` | `chorus run configs/panel_transfer_CORSIP2026.yaml` | configs/panel_transfer_CORSIP2026.yaml | 104 (started 2026-11-17 10:34:14) |
| `panel_transfer_CORSIP2026_noUnverified` | `chorus run configs/panel_transfer_CORSIP2026_noUnverified.yaml` | configs/panel_transfer_CORSIP2026_noUnverified.yaml | 104 (started 2026-11-17 10:35:58) |
| `lc_transcriptomics_blood_v4_GSE228320` | `chorus run configs/lc_transcriptomics_blood_v4_GSE228320.yaml` | configs/lc_transcriptomics_blood_v4_GSE228320.yaml | 141 (started 2026-11-17 10:37:42) |
| `concordance_vs_transfer_v5` | `scripts/concordance_vs_transfer.py --out results/concordance_vs_transfer_v5 --blocks transcriptomics=results/lc_transcriptomics_blood_v4/loso_results.csv,proteomics_affinity=results/lc_proteomics_affinity_v5/loso_results.csv,proteomics_ms=results/lc_proteomics_ms_v4/loso_results.csv,metabolomics=results/lc_metabolomics/loso_results.csv,proteomics_explore=results/lc_proteomics_explore_v5/loso_results.csv --not-pooled proteomics_explore` | the loso_results.csv of lc_transcriptomics_blood_v4, lc_proteomics_affinity_v5, lc_proteomics_ms_v4, lc_metabolomics, lc_proteomics_explore_v5 | 15 (started 2026-11-17 10:40:03) |
| `clustering_readiness_v5` | `scripts/clustering_readiness.py --out results/clustering_readiness_v5 --blocks rna=configs/lc_transcriptomics_blood_v4.yaml,protein_affinity=configs/lc_proteomics_affinity_v5.yaml,protein_explore=configs/lc_proteomics_explore_v5.yaml` | configs/lc_transcriptomics_blood_v4.yaml, configs/lc_proteomics_affinity_v5.yaml, configs/lc_proteomics_explore_v5.yaml | 328 (started 2026-11-17 10:40:18) |
| `acute_persistence_v5` | `scripts/acute_persistence.py --out results/acute_persistence_v5 --target-protein configs/concordance_proteomics_all_v5.yaml --acute-check Woodruff2023` | configs/control_acute_proteomics.yaml, configs/concordance_proteomics_all_v5.yaml, configs/mosaic_acute_rna.yaml, configs/mosaic_acute_protein.yaml, configs/mosaic_lc_rna_v4.yaml | 180 (started 2026-11-17 10:45:46) |
| `lc_transcriptomics_blood_v4_panel495` | `chorus run configs/lc_transcriptomics_blood_v4_panel495.yaml` | configs/lc_transcriptomics_blood_v4_panel495.yaml | 301 (started 2026-11-17 11:36:16) |
| `panel_transfer_CORSIP2026_checks` | `scripts/panel_transfer_checks.py confounding`; `permute --n-perm 200` in four chunks of 50; `concordance`; `permutation_summary` | configs/panel_transfer_CORSIP2026.yaml | confounding 87; permutation chunks 311 to 314 and concordance 146, in parallel (started 2026-11-17 11:39:41); summary 1 |
| `lc_transcriptomics_blood_v4_panel495/random_genes` | `scripts/panel_transfer_checks.py random --n-sets 20` in four chunks of 5 | configs/lc_transcriptomics_blood_v4.yaml, configs/lc_transcriptomics_blood_v4_panel495.yaml (grid) | 1157 to 1167 per chunk, in parallel (started 2026-11-17 11:47:01; a first attempt at 11:39:41 failed in every fit because the gene ids were numpy strings, and was deleted) |
| `concordance_noise_ceiling` | `scripts/concordance_noise_ceiling.py --out results/concordance_noise_ceiling` | configs/concordance_rna_v5.yaml, configs/concordance_proteomics_all_v5.yaml, configs/concordance_metabolomics_all.yaml | 232 (started 2026-11-17 22:21:00) |
| `concordance_guided` | `scripts/concordance_guided.py all --out results/concordance_guided --workers 10` | the configs of the six v4 and v5 blocks (grid set in the script) | 990 (started 2026-11-18 07:07:55); summaries rewritten with `summarise` after a fix of the bar check |
| `concordance_vs_transfer_full` | `scripts/concordance_vs_transfer_full.py --out results/concordance_vs_transfer_full` | the loso_results.csv of the LC blocks and of control_acute_proteomics and control_sex_label | 14 (started 2026-11-18 07:08:31) |
| `enrichment.csv` | `scripts/enrichment.py configs/lc_transcriptomics_blood.yaml results/lc_transcriptomics_blood/selected_features.json --out results/enrichment.csv` | configs/lc_transcriptomics_blood.yaml (expressed universe) | 42 (version 3, started 2026-11-13 15:45:48) |

- `lc_transcriptomics_blood`: primary leave-one-study-out benchmark, full grid (9 harmonisers, among them `chorus_protected` since version 3, plus nested selection,
  6 selector configurations, 4 classifiers), with the within-study cross-validation reference.
- `lc_transcriptomics_blood_noGSE270045`: identical except GSE270045 is disabled (sensitivity run;
  `docs/gse270045_check.md`).
- `lc_acute_transfer`: post-acute cohorts as training only, GSE215865 (acute admissions) as the single held-out study
  (`docs/acute_transfer.md`).
- `simulation_factors`: `scripts/simulation_benchmark.py` with its defaults (6 simulated studies, 3 sample-level
  technical factors).
- `simulation_shifts_only`: the same with `--n-factors 0` (study location and scale shifts only).
- `lc_transcriptomics_blood_within_v2`: within-study cross-validation only (`--within-only`) for the primary config, with
  the revised method: stratified 5-fold CV repeated 10 times, metrics per fold averaged over the 50 folds, and the pooled
  value kept for comparison. The leave-one-study-out benchmark was not rerun.
- `control_label_permutation`: label permutation control, 20 permutations. See "Control runs" below.
- `control_sex_label`: positive control with inferred sex as the label, on all filtered genes, without chrX and chrY
  genes, and with the sex genes forced in. It replaces the `sample_matrix` control (`control_matrix_label`, removed; it is in git history at commit
  aad458e). See "Inferred sex" and "Control runs" below.
- `lc_transcriptomics_blood_harmonisation`: mixing and conservation metrics of every harmoniser on every LOSO training
  pool of the primary config. See "Harmonisation diagnostics" below.
- `concordance`: within-cohort LC effects and their pairwise Spearman correlation across cohorts. See "Concordance" below.
- `enrichment.csv`: enrichment of each fold's invariant genes (from `lc_transcriptomics_blood/selected_features.json`) in
  the prior gene lists. See "Enrichment" below. Rerun in version 3 on the current `selected_features.json`.
- Folders ending in `_v1`: the results before the version 2 changes, kept for comparison. See "Version 2 changes" below.
  The version 2 results replaced in version 3 were not copied to new folders; they are in git at commit 44dffee.

The cohort table each real-data run loaded (samples and labels per study) is printed at the top of its log,
`results/<folder>.log`. The logs stay local (not committed). Simulation runs write it to `cohorts.csv`;
`control_sex_label/cohorts.csv` gives the sex counts per cohort.

**Two versions of within_study_cv.csv.** The `within_study_cv.csv` files in `lc_transcriptomics_blood`,
`lc_transcriptomics_blood_noGSE270045`, `lc_acute_transfer`, the simulation and smoke folders come from the earlier method:
one stratified 5-fold split with the out-of-fold predictions pooled into one AUROC per study. Their `auroc` equals
`auroc_pooled_repeat0` in `lc_transcriptomics_blood_within_v2` (checked: largest difference 0).

`smoke_transcriptomics` and `smoke_acute_transfer` are smoke tests from 12 November 2026 on the same cohorts with a reduced
grid: harmonisers identity and center_scale, selectors variance, differential and invariant, classifier logreg, no nested
selection, 200 bootstrap resamples.

## Cohorts in the real-data runs

| run | held out | training only |
|---|---|---|
| lc_transcriptomics_blood | GSE226260_WB, GSE226260_PBMC, GSE224615, GSE251849, GSE270045 | none |
| lc_transcriptomics_blood_noGSE270045 | GSE226260_WB, GSE226260_PBMC, GSE224615, GSE251849 | none |
| lc_acute_transfer | GSE215865 | GSE226260_WB, GSE226260_PBMC, GSE224615, GSE251849, GSE270045 |
| control_label_permutation | as lc_transcriptomics_blood | none |
| control_sex_label (and _v1) | GSE226260_WB, GSE226260_PBMC, GSE224615, GSE270045 (GSE251849 skipped: 4 males, below 5 per class) | none |
| lc_transcriptomics_blood_harmonisation | as lc_transcriptomics_blood (training pools of its folds) | none |
| lc_proteomics_affinity (and its harmonisation and concordance) | Buggert2025, Cervia2024_SomaScan, Su2022 | none |
| lc_proteomics_ms (and its harmonisation and concordance) | Cervia2024_MS, Chowdhury2024, Kovarik2023_plasma | none |
| control_acute_proteomics | Cervia2024_SomaScan, Su2022, Zoodsma2022_Hannover | none |
| concordance_proteomics_all | no folds; the six cohorts of the two blocks | none |
| lc_proteomics_affinity_v4 (and its harmonisation and concordance) | Buggert2025, Cervia2024_SomaScan, Su2022, Fineschi2026 | none |
| lc_proteomics_ms_v4 (and its harmonisation and concordance) | Cervia2024_MS, Chowdhury2024, Kovarik2023_plasma, Chowdhury2026 | none |
| concordance_proteomics_all_v4 | no folds; the eight cohorts of the two v4 blocks | none |
| lc_transcriptomics_blood_v4 (and `lc_transcriptomics_blood_harmonisation_v4`) | GSE226260_WB, GSE226260_PBMC, GSE224615, GSE251849, GSE270045, PRJNA1149196_SA, GSE297499 | PRJNA1149196_QLD (5 / 4 samples) |
| lc_transcriptomics_blood_v4_noQCflag | as lc_transcriptomics_blood_v4, GSE224615 without GSM7027478, GSM7027490 and GSM7027501 | PRJNA1149196_QLD |
| concordance_v4 | no folds; the eight cohorts of lc_transcriptomics_blood_v4 | none |
| mosaic_programmes_v4/lc_v4 | test cohorts: the RNA and protein cohorts except PRJNA1149196_QLD (5 / 4); Cervia2024_MS as sensitivity | PRJNA1149196_QLD |
| early_detection/horizontal | Su2022 and Cervia2024_SomaScan (acute samples), each trained on the other | none |
| early_detection/vertical, early_detection/reference | no folds across cohorts; within Su2022 (paired) and within GSE215865 | none |
| clustering_readiness | no labels; LC cases of the lc_transcriptomics_blood_v4 and lc_proteomics_affinity_v4 cohorts, each cohort left out in turn | none |
| lc_proteomics_affinity_v5 | Buggert2025, Cervia2024_SomaScan, Su2022, Fineschi2026, Woodruff2023, Bansal2025 | none |
| lc_proteomics_explore_v5 | Buggert2025, Fineschi2026, Woodruff2023 | none |
| concordance_proteomics_all_v5, concordance_v5 | no folds; the ten protein cohorts, and the eight v4 RNA cohorts with CORSIP2026 and GSE228320 | none |
| mosaic_programmes_v5/lc_v5 | as mosaic_programmes_v4/lc_v4, with Woodruff2023 and Bansal2025 | PRJNA1149196_QLD |
| panel_transfer_CORSIP2026 (and _noUnverified) | CORSIP2026 | the eight lc_transcriptomics_blood_v4 cohorts |
| lc_transcriptomics_blood_v4_GSE228320 | GSE228320 | the eight lc_transcriptomics_blood_v4 cohorts |
| acute_persistence_v5 | scored: Buggert2025, Fineschi2026, Woodruff2023, Bansal2025 (proteins); as acute_persistence (programmes) | as acute_persistence |
| acute_persistence | scored: Buggert2025, Fineschi2026 (proteins); the RNA cohorts of mosaic_lc_rna_v4 except PRJNA1149196_QLD (programmes) | the acute cohorts of control_acute_proteomics (and GSE215865 for programmes) |

## GSE224615 samples flagged at QC

GSE224615 was re-quantified with salmon (`docs/provenance.md`). A sample is flagged when its mapping rate is below 40% or
more than 2.5 MAD (raw) below the cohort median: median 81.57%, MAD 1.77, lower bound 77.13%.

| sample | mapping rate (%) | mitochondrial fraction | reason |
|---|---|---|---|
| GSM7027478 | 52.53 | 0.580 | more than 2.5 MAD below median |
| GSM7027490 | 59.07 | 0.457 | more than 2.5 MAD below median |
| GSM7027501 | 56.84 | 0.052 | more than 2.5 MAD below median |

**All three are included in every run above except `lc_transcriptomics_blood_v4_noQCflag`**, the sensitivity run of
`lc_transcriptomics_blood_v4` without them (`exclude_samples` in `configs/lc_transcriptomics_blood_v4_noQCflag.yaml`;
see "Version 4 cohort set and reporting fixes", section 3). The full per-sample QC
table is `data/quant/GSE224615/salmon_qc.tsv` (local, not committed).

## Files in each folder

Only the CSV files, this README and `runtimes.tsv` are committed. `selected_features.json`, the figures and the logs stay
local.

- `loso_results.csv`: one row per held-out study, harmoniser, selector and classifier.
- `loso_summary.csv`: AUROC across held-out studies per harmoniser, selector and classifier.
- `loso_predictions.csv`: per-sample predicted probabilities.
- `within_study_cv.csv`: 5-fold stratified cross-validation inside each study (written with `--within`).
- `selected_features.json`: genes selected per held-out study, harmoniser and selector (real-data runs).
- `within_study_cv_folds.csv`: one row per study, selector, classifier, repeat and fold (`lc_transcriptomics_blood_within_v2`).
- `permutation_results.csv`, `permutation_summary.csv`: `control_label_permutation` only.
- `cohorts.csv`: simulation runs and `control_sex_label`. `figure_simulation.png`: simulation runs only.
- `control_sex_label`: `inferred_sex.csv`, `cohorts.csv`, `feature_sets.csv`, `forced_genes_by_cohort.csv`,
  `forced_genes_by_fold.csv`, `auroc_by_held_out.csv`, and the loso files in `all_genes/`, `no_xy/` and
  `forced_sex_genes/`. `control_sex_label_v1` has no forced set and no forced gene tables.
- `lc_transcriptomics_blood_harmonisation`: `<held_out>/<harmoniser>.csv` (one table per fold and harmoniser) and
  `all_metrics.csv` (all of them stacked).
- `concordance`: `effects.csv`, `spearman_all_genes.csv`, `spearman_all_genes_p.csv`, `spearman_top500.csv`,
  `spearman_top500_p.csv`, `top500_genes.csv`, `cohort_annotation.csv`, and since the reporting rerun of 16 November
  `n_shared.csv`, `n_shared_top500.csv` and `pairs.csv` (see "Version 4 cohort set and reporting fixes"). Every
  concordance folder has these files, with its own `top<N>`.
- Proteomics folders (`lc_proteomics_affinity`, `lc_proteomics_ms`, `control_acute_proteomics`): the loso and
  within-study files above, plus `cohorts.csv` (samples, class counts and `n_values_imputed` per cohort) and
  `features_observed_fraction.csv` (per protein and cohort, the fraction of the cohort's samples with a value before
  imputation; one row per protein used, the proteins shared by all cohorts of the config).
  `lc_proteomics_affinity` and `lc_proteomics_ms` also hold `harmonisation/` (as `lc_transcriptomics_blood_harmonisation`)
  and `concordance/` (as `concordance`, with `top50` in place of `top500`, and `n_shared.csv`).
- `concordance_proteomics_all`: as `concordance`, with `top50` files and `n_shared.csv` (proteins with an effect in both
  cohorts of each pair).
- `control_acute_transcriptomics`: `cohorts.csv` only (see "Proteomics (15 November 2026)").
- `paired_<cohort>/<contrast>/` (`lc`, `acute`): `samples.csv`, `features.csv`, `fold_metrics.csv`, `summary.csv`,
  `predictions.csv`, `probability_correlation.csv`, and for Cervia2024 `effects.csv` and `effect_correlation.csv`;
  `paired_<cohort>/run_info.json` gives the class counts and why a contrast was skipped (see "Test 2").
- `mosaic_programmes/`: `genesets.csv`, `run_info.json` and, per contrast (`acute`, `lc`), `cohorts.csv`,
  `coverage.csv`, `effects.csv`, `effects_var.csv`, `concordance_<set>.csv` with `_p` and `n_shared_<set>.csv`
  (`all`, `hallmark`, `reactome`), `partition.csv`, `loso.csv`, `loso_predictions.csv`, and for `lc`
  `bridge_cervia.csv` (see "Test 3"), and `pairs_<set>.csv` (pairs with `n_shared` and `fewer_than_10_shared`).
  `scores.csv.gz` is not committed.
- `lc_metabolomics/`: as `lc_proteomics_affinity` (with `harmonisation/` and `concordance/` with `top20` files), plus
  `concordance_all/` (four cohorts, pairwise shared, `top50` files) and `mapping/` (`mapping_summary.csv`,
  `n_shared.csv`, `shared_by_k.csv`); see "Metabolomics".
- `lc_proteomics_affinity_v4/`, `lc_proteomics_ms_v4/`: as `lc_proteomics_affinity`; `concordance_proteomics_all_v4/`:
  as `concordance_proteomics_all`.
- `concordance_epigenomics/`: `cohorts.csv` only (groups and why no effect was computed; see "Metabolomics").
- `lc_transcriptomics_blood_v4/`, `lc_transcriptomics_blood_v4_noQCflag/`: as `lc_transcriptomics_blood`, with
  `within_study_cv_folds.csv` (the current within-study method); `lc_transcriptomics_blood_harmonisation_v4/`: as
  `lc_transcriptomics_blood_harmonisation`; `concordance_v4/`: as `concordance`; `mosaic_programmes_v4/`:
  `genesets.csv`, `run_info.json` and `lc_v4/` with the files of `mosaic_programmes/lc` (`scores.csv.gz` not
  committed).

- `early_detection/`: `horizontal/` (`cohorts.csv`, the loso files, `permutation_results.csv`,
  `permutation_summary.csv`, `run_info.json`), `vertical/` (as `paired_<cohort>`, contrast `early`) and `reference/`
  (`cohorts.csv`, `within_study_cv.csv` with `published_taco_auroc`, `within_study_cv_folds.csv`, `run_info.json`).
- `concordance_vs_transfer/`: `rows.csv`, `correlation.csv`, `variance.csv`, `run_info.json`.
- `clustering_readiness/`: `cohorts.csv`, `clusters.csv`, `stability.csv`, `stability_summary.csv`, `composition.csv`,
  `run_info.json`.
- `acute_persistence/protein/`, `acute_persistence/rna_programmes/`: `per_cohort.csv`, `meta.csv`, `model.csv`,
  `sharing.csv`, `scores.csv` (acute score per target sample), `run_info.json`.

- Version 5: `lc_proteomics_affinity_v5/`, `lc_proteomics_explore_v5/` as `lc_proteomics_affinity_v4`
  (without `harmonisation/` and `concordance/`); `concordance_proteomics_all_v5/`, `concordance_v5/` as
  `concordance_proteomics_all`; `mosaic_programmes_v5/` with `lc_v5/` as `mosaic_programmes_v4/lc_v4`;
  `panel_transfer_CORSIP2026/`, `panel_transfer_CORSIP2026_noUnverified/`, `lc_transcriptomics_blood_v4_GSE228320/`:
  `loso_results.csv`, `loso_summary.csv`, `loso_predictions.csv`; `concordance_vs_transfer_v5/`,
  `clustering_readiness_v5/`, `acute_persistence_v5/` as their version 4 folders, plus `acute_check.csv` in
  `acute_persistence_v5/protein/`.

- `panel_transfer_CORSIP2026_checks/`: `corsip_demographics.csv`, `corsip_single_predictors.csv`,
  `within_sex_auroc.csv`, `adjusted_association.csv`, `refit_check.csv`, `top_weights.csv`,
  `corsip_sex_difference.csv`, `permutation_null.csv` (and its four chunks `permutation_null_NNN.csv`),
  `permutation_summary.csv`, `effects_495.csv`, `concordance_with_corsip.csv`, `concordance_summary.json`,
  `leave_one_training_cohort_out.csv`. `lc_transcriptomics_blood_v4_panel495/`: the loso files, `cohorts.csv` is not
  written (RNA), `comparison_by_held_out.csv`, `median_of_mean_auroc.csv`, and `random_genes/` (`setNN/genes.txt`,
  `setNN/loso_results.csv`, `summary_by_set.csv`).

- `concordance_noise_ceiling/`: `cohorts.csv`, `pairs.csv`, `summary.csv`, `run_info.json`.
  `concordance_guided/`: `arms.csv`, `differences.csv`, `bars.csv`, `empty.csv`, `skipped.csv`, `per_cohort.csv`,
  `per_split.csv`, `agreement.csv`, `run_info.json`, `units/`, `blocks/<block>/`. `concordance_vs_transfer_full/`:
  `points.csv`, `correlation.csv`, `monotone.csv`, `curve.csv`, `run_info.json`.

## Column glossary

**loso_results.csv**

| column | meaning |
|---|---|
| held_out | study left out and used as the test set |
| harmoniser | harmoniser in the grid; `nested` means it was chosen by an inner leave-one-study-out over the training studies |
| harmoniser_chosen | harmoniser actually applied (differs from `harmoniser` only for `nested`) |
| leak_risk | true when the harmoniser is transductive (fitted on all studies, without labels) rather than on training studies only |
| n_factors | number of technical factors the `chorus` harmoniser removed in this fit (also for `nested` rows that chose it); empty for other harmonisers and in the `_v1` folders |
| selector | selector configuration: `variance`, `differential`, `invariant`, or `prior:<gene list>:<mode>[:k<k>]` |
| n_features | number of model inputs: selected genes plus covariate columns |
| n_features_dropped_incomplete | features left out of the fold because a training or held-out sample had no value (only `control_sex_label`, set `forced_sex_genes`; 0 otherwise) |
| n_strict_invariant | genes that passed the strict invariance screen (invariant selector; empty otherwise) |
| sign_agreement | over the selected genes, the fraction whose Hedges g between LC and COVID_recovered in the held-out study has the same sign as the pooled training effect; uses the held-out labels and is a diagnostic only. Empty when the selector has no effect table (variance, prior_only) or the held-out study has fewer than 3 samples per class |
| effect_spearman | Spearman correlation between the held-out Hedges g and the pooled training effect over the same genes; empty under the same conditions |
| classifier | `logreg` (L2 logistic regression, C 0.5, balanced class weights), `rf` (random forest), `taco` (soft voting ensemble; in every real-data run **without TabPFN**, `CHORUS_DISABLE_TABPFN=1`, so rf, extratrees, svc and mlp only), `cross_study` (one logreg per training study with at least 3 samples per class, the members' probabilities averaged with weights proportional to each member's mean AUROC on the other training studies minus 0.5, floored at 0; Patil and Parmigiani 2018) |
| n_train | training samples |
| n_train_studies | training studies |
| auroc, auprc | area under the ROC and precision-recall curves on the held-out study |
| balanced_accuracy, f1, mcc | at a probability threshold of 0.5; mcc is 0 when only one class is predicted |
| brier | Brier score |
| ece | expected calibration error |
| n | held-out samples |
| prevalence | fraction of LC among held-out samples |
| auroc_lo, auroc_hi | 95% percentile interval of the AUROC from a class-stratified bootstrap over held-out samples (`n_boot` resamples) |
| inner_scores | for `nested`: mean inner AUROC per candidate harmoniser |
| seconds | time for this harmoniser and selector fit plus this classifier |
| error | present only on rows where a fit failed; the message |

**loso_summary.csv**: per harmoniser, selector and classifier, the `mean`, `median`, `std`, `min` and `max` of AUROC
across held-out studies, `n_studies`, and the mean AUPRC (`auprc_mean`) and Brier score (`brier_mean`). Failed rows are
excluded from these values. Since 16 November 2026 (commit 4ca16b1) the summary also has `n_errors` (failed held-out
rows) and `error` (their distinct messages), and a configuration that failed in every held-out study keeps a row
with `n_studies` 0 and empty metrics (classifier `not_fitted` when the failure came before any classifier). Summaries
written before then (other folders) drop such configurations; `lc_metabolomics` and `control_acute_proteomics` were
rewritten with the new rule.

**loso_predictions.csv**: `sample`, `held_out`, `harmoniser`, `selector`, `classifier`, true label `y` (1 = LC) and
predicted probability `p`.

**within_study_cv.csv** (`lc_transcriptomics_blood_within_v2`):

| column | meaning |
|---|---|
| study, selector, classifier | as above |
| auroc, auprc, balanced_accuracy, f1, mcc, brier, ece | computed on each test fold, then averaged over the `n_folds` folds |
| auroc_sd | standard deviation of the fold AUROCs |
| auroc_pooled | out-of-fold predictions pooled across the 5 folds of a repeat into one AUROC, averaged over repeats |
| auroc_pooled_repeat0 | the pooled AUROC of the first repeat, which uses the split of the earlier method |
| n, prevalence | samples in the study and fraction LC |
| n_folds, n_failed_folds | folds with metrics, and folds whose fit failed |
| n_splits, n_repeats | 5 and 10. Repeat r uses split seed `random_state + r` |

`within_study_cv_folds.csv` has the per-fold values with `repeat` and `fold`; there, `n` and `prevalence` refer to the
test fold.

**within_study_cv.csv** (earlier method, other folders): `study`, `selector`, `classifier` and the same metrics as in
loso_results.csv, from one pooled 5-fold stratified cross-validation inside the study. The invariant selector falls back to differential, because one study gives one
environment. `cross_study` is skipped. Studies with fewer than 5 samples in a class are skipped. Training-only studies
are included. This applies to both versions.

**permutation_results.csv**: the loso_results.csv columns plus `permutation` (0 to 19).

**permutation_summary.csv**: per `held_out`, `harmoniser`, `selector` and `classifier`, the `mean`, `std`, `min` and `max`
of AUROC over permutations, and `n_permutations`.

## Control runs

Both use `scripts/control_runs.py` and start from the cohorts, preprocessing and covariates of
`configs/lc_transcriptomics_blood.yaml`. TabPFN is not involved, because neither control uses `taco`.

**control_sex_label (positive control).**
- **Label:** `inferred_sex` (female positive, male negative; see "Inferred sex"), on the samples of the primary contrast
  (LC or COVID_recovered) after the primary preprocessing. Ambiguous samples are left out (one, GSM7069947).
- **Counts:** `cohorts.csv`. GSE251849 has 4 males in the contrast and is skipped as a held-out study
  (`min_test_per_class` 5); it is still used for training.
- **Grid:** harmonisers identity, center_scale, combat_addon, chorus and chorus_protected (since version 3); selectors variance (k 500) and invariant (k 200);
  classifier logreg; no nested selection; 1000 bootstrap resamples. Covariates and `protect` are `sample_matrix` as in
  the primary config.
- **Three feature sets** (`feature_sets.csv`; `_v1` has the first two):
  - `all_genes`: the 10,442 genes that pass the primary expression filter. These include 320 chrX genes and no chrY
    gene; XIST and RPS4Y1 are not among them, because the filter requires CPM of at least 1 in half the samples of every
    cohort.
  - `no_xy`: the same without chrX and chrY genes (10,122 genes). Chromosomes come from the GENCODE v47 tx2gene table
    (`data/quant/GSE224615/tx2gene_gencode_v47.tsv`); genes on both chrX and chrY (pseudoautosomal) count as both.
    Every filtered gene had a chromosome.
  - `forced_sex_genes` (version 2): `all_genes` plus XIST, RPS4Y1 and the chrY genes (chrY only, not
    pseudoautosomal) present in at least one cohort's matrix, 681 candidates in all. In each cohort a candidate is
    `used` when it has a value for every sample and, for chrY genes other than RPS4Y1, when its log2 CPM (prior count 2)
    is at least log2(1) in at least half of the cohort's inferred male samples (the config's expression thresholds
    applied to males); otherwise it is set to missing in that cohort (`forced_genes_by_cohort.csv`: `in_matrix`,
    `frac_males_expressed`, `n_males`, `used`). In each fold a candidate is a feature only when it is used in every
    training cohort and in the held-out cohort, because the held-out samples must be scored
    (`forced_genes_by_fold.csv`: `in_every_training_cohort`, `in_held_out`, `used`). The other candidates are dropped
    in that fold (`n_features_dropped_incomplete`, 673 in every fold).
  - In every fold the 8 forced genes used are RPS4Y1, KDM5D, DDX3Y, ZFY, USP9Y, TTTY14, UTY and EIF1AY. XIST is used in
    no fold: GSE270045 has no XIST row and is in every fold, as a training or the held-out cohort. In the GSE270045
    fold, XIST, TMSB4Y and LINC00278 are used in every training cohort but not in GSE270045.
  - In `selected_features.json` of `forced_sex_genes`, all 8 forced genes are among the selected genes for the variance
    selector with identity, center_scale and combat_addon in every fold, none for chorus with the variance selector,
    and none for the invariant selector. In the GSE224615 fold (identity), the 8 genes have I squared between 92 and 97
    and Q_p 0, so they fail the strict invariance screen (I squared at most 50), and 254 genes pass it; the invariant
    rows of `forced_sex_genes` therefore equal those of `all_genes`. This paragraph describes the version 3 files; since
    the 15 November rerun the forced genes are appended to every selector's output ("Proteomics (15 November 2026)", item 6).
- **auroc_by_held_out.csv**: `feature_set`, `held_out`, `harmoniser`, `selector`, `classifier`, `n_features`, `auroc`,
  `auroc_lo`, `auroc_hi`, `error`. The full rows are in `all_genes/loso_results.csv` and `no_xy/loso_results.csv`.

**control_label_permutation.**
- **Permutation:** within each cohort, `lc_status` is permuted among the samples labelled LC or COVID_recovered, so class
  counts per cohort are unchanged.
- **Benchmark:** the primary leave-one-study-out benchmark with harmonisers identity and chorus, the invariant selector
  (k 200) and logreg. Covariates and `protect` are `sample_matrix` as in the primary config. There is no nested selection.
- **Repeats:** 20 permutations (seed 0), 1000 bootstrap resamples per fit (the interval columns are not summarised).

## Inferred sex

`scripts/infer_sex.py` (module `chorus.sex`) infers sex inside each cohort of the primary config, on all samples of the
raw matrix (before one sample per subject and before the expression filter).

- **Input:** log2 CPM of XIST (ENSG00000229807) and RPS4Y1 (ENSG00000129824), library size over all genes, prior count 2.
- **Split:** a two-component Gaussian mixture (full covariance, `reg_covar` 0.25, 5 initialisations, seed 0). The
  component with the higher mean XIST minus RPS4Y1 is female.
- **Ambiguous:** posterior of the assigned component below 0.95 (`low_posterior`), or squared Mahalanobis distance to it
  above the chi-square 0.999 quantile (`outlier`).
- **GSE270045:** the author matrix has no XIST row, so the split uses RPS4Y1 alone (lower component female). The method
  string says so.
- **Stored:** columns `inferred_sex` and `inferred_sex_method` (method, thresholds, script and date) in
  `data/meta/<study_id>_meta.csv`, and the per-sample table `data/meta/inferred_sex.csv` (local), copied to
  `control_sex_label/inferred_sex.csv`. `scripts/build_cohort_meta.py` merges it back when it rebuilds the tables.

| cohort | samples | female | male | ambiguous | genes | agreement with recorded sex |
|---|---|---|---|---|---|---|
| GSE226260_WB | 103 | 67 | 36 | 0 | XIST, RPS4Y1 | no recorded sex |
| GSE226260_PBMC | 86 | 69 | 16 | 1 | XIST, RPS4Y1 | 85 of 85 non-ambiguous |
| GSE224615 | 36 | 20 | 16 | 0 | XIST, RPS4Y1 | 36 of 36 |
| GSE251849 | 23 | 18 | 5 | 0 | XIST, RPS4Y1 | no recorded sex |
| GSE270045 | 36 | 28 | 8 | 0 | RPS4Y1 only | no recorded sex |

- The ambiguous sample is GSM7069947 (GSE226260_PBMC, recorded female): XIST -0.22 and RPS4Y1 -0.36 log2 CPM, squared
  Mahalanobis distance 35.3. It is the sample with the malformed matrix tail (`docs/gse226260_check.md`).
- No subject with several samples has different inferred sex across its samples.

**inferred_sex.csv** columns: `study_id`, `sample_id`, `inferred_sex`, `posterior` (of the assigned component),
`mahalanobis2`, `reason` (empty, `low_posterior` or `outlier`), `xist_logcpm`, `rps4y1_logcpm`, `provenance`,
`recorded_sex` (the curated `sex` column; empty where not recorded).

## Harmonisation diagnostics

`chorus harmonisation` (module `chorus.eval.harmonisation`) uses the folds of the primary benchmark. For each held-out
study, the training pool is the labelled (LC or COVID_recovered) samples of the other studies. Each of the 9 harmonisers
(8 in `_v1`; `chorus_protected` added in version 3) is fitted as in the benchmark: on the training pool with its labels where the operator uses them, and for
`combat_pooled` on the pool plus the held-out samples without labels. The pool is compared before (log CPM as loaded)
and after.

- **Genes.** Mixing and LC effect metrics use the 10,442 filtered genes. XIST, RPS4Y1 and the marker genes are added to
  the harmonised matrix (none of the 10 markers needed adding; XIST and RPS4Y1 did). A gene is added in a fold only
  when every sample the harmoniser is fitted on has a value. GSE270045 has no XIST row, so XIST enters only the fold
  where GSE270045 is held out, and not for `combat_pooled` (which also fits on the held-out samples). Column
  `extra_genes` lists what was added.
- **Markers.** Neutrophil: FCGR3B, CSF3R, CXCR2, FPR1, MNDA. Lymphocyte: CD3E, CD3D, CD79A, MS4A1, CD2 (the lists in
  `chorus.harmonise.transcriptomics.composition_covariates`).
- **Sex labels** are `inferred_sex`; ambiguous samples are left out of the sex metrics.

| metric | family | definition |
|---|---|---|
| silhouette_study | technical_mixing | mean silhouette width with study as the label, Euclidean, on the top 50 principal components of the pool after each gene is standardised across the pool |
| knn_other_study | technical_mixing | mean over samples of the fraction of the 20 nearest neighbours (same space) from another study |
| knn_other_study_expected_random | technical_mixing | the same fraction under random mixing, mean of (N - n_study) / (N - 1) |
| lc_effect_spearman | biological_conservation | per training study, Spearman correlation over genes of the within-study Hedges g (LC minus COVID_recovered) before versus after; `mean_over_studies` is their mean |
| sex_auroc_XIST | biological_conservation | AUROC of XIST for female versus male over the pool |
| sex_auroc_RPS4Y1 | biological_conservation | AUROC of RPS4Y1 for male versus female over the pool |
| silhouette_sex | biological_conservation | mean silhouette width by inferred sex on XIST and RPS4Y1 (each standardised across the pool), on those of the two available |
| marker_score_spearman_neutrophil, marker_score_spearman_lymphocyte | biological_conservation | score = mean over the marker genes of the values standardised across the pool; Spearman correlation of the score before versus after, over the pool and within each study (`mean_over_studies` is the mean of the within-study values) |

**Table columns** (`<held_out>/<harmoniser>.csv` and `all_metrics.csv`): `held_out`, `harmoniser`, `family`, `metric`,
`scope` (`pool`, a study id or `mean_over_studies`), `before` and `after` (mixing and sex separation), `before_vs_after`
(LC effect and marker comparisons), `n` (samples, or studies for `mean_over_studies` of lc_effect_spearman), `n_genes`,
`n_pcs`, `k`, `leak_risk`, `n_factors` (`chorus` and `chorus_protected` only; not in `_v1`), `extra_genes`. A failed fit would give one row with `family` = `error` and the message in
`scope`; there were none.

**combat_addon with protect.** The `_v1` folder was made with the version 1 `AddonComBat`, in which a protected
covariate constant within a study moved the study location (see "Version 2 changes"). The current folder was made
with the version 3 `AddonComBat` (see "Version 3 changes").

## Concordance

`scripts/concordance.py` on the primary config, after its preprocessing (one sample per subject, deduplication).

- **Effect:** inside each cohort, Hedges g of LC minus COVID_recovered per gene on log CPM, over the 10,442 common
  expressed genes (the genes that pass the expression filter in every cohort).
- **Correlation:** pairwise Spearman correlation of g between cohorts, over all common genes and over the 500 genes
  with the largest mean |g| across the five cohorts.
- **Annotation** (`cohort_annotation.csv`), from the metadata: `case_definition` (of LC samples), `control_type`
  (`lc_status` of the negative class), `control_case_definition`, `time_anchor` (`missing` when the column is empty),
  `sample_matrix`, `n_LC`, `n_control`, `other_lc_status` (other statuses present in the cohort and not used), and
  median, minimum and maximum `days_since_infection` per class (empty when not recorded).

**Files:** `effects.csv` (`gene_id`, `symbol` from the HGNC complete set, one column of g per cohort);
`spearman_all_genes.csv` and `spearman_top500.csv` (cohort by cohort correlation matrices); `*_p.csv` (two-sided p
values of the same pairs); `top500_genes.csv` (`gene_id`, `symbol`, `mean_abs_g`).

## Enrichment

`scripts/enrichment.py` takes every `<held_out>|<harmoniser>|invariant` entry of
`results/lc_transcriptomics_blood/selected_features.json`: 5 folds by 10 harmonisers (including `nested` and
`chorus_protected`), 200 genes each.

- **Test:** one-sided hypergeometric (overlap at least as large as observed) over the expressed universe, the 10,442
  genes that pass the primary expression filter, which is the set the selector chose from. List genes outside the
  universe are dropped first.
- **Lists:**
  - `TACO_411`: `priors/taco_411_dce_genes.txt`, 411 genes, 333 in the universe.
  - `MRCONTROL_32`: `priors/mrcontrol_32_genes.txt`, 32 genes, 25 in the universe.
  - `review_complement`: `priors/review_complement_9.txt`. The file holds 11 UniProt accessions (C4A, C4BPB, C5, CFHR1,
    CFHR5, MBL2, FCN2, MASP2, C1QA, C1QB, C1QC). They map to 11 Ensembl genes through the `uniprot_ids` column of the
    HGNC complete set, and 2 are in the universe.
- **q_bh:** Benjamini-Hochberg over all 150 rows (135 before version 3).

**enrichment.csv** columns: `held_out`, `harmoniser`, `selector`, `gene_list`, `list_size` (genes after mapping),
`universe`, `list_in_universe`, `n_selected`, `overlap`, `expected` (n_selected x list_in_universe / universe),
`fold_enrichment` (overlap / expected), `p_hypergeom`, `q_bh`, `overlap_genes` (symbols).

## Version 2 changes (13 November 2026)

Three code changes, after which the affected results were rerun. The results before them are kept in the `_v1` folders.
The `_v1` results therefore differ from the current ones by both the `combat_addon` and the `chorus` change.

**1. AddonComBat (commit ee8248a).** When a protected covariate is constant within a training study, it is now dropped
from that study's design, so it no longer moves the study location. In version 1 the design had an intercept plus the
protected columns plus the label; a protected column constant within the study (as `sample_matrix` is in every cohort)
made the least-squares intercept non-unique, and the minimum-norm intercept of `numpy.linalg.lstsq` was used as the
study location. The label column is unchanged, so the location is still the intercept at label 0 when the label is in
the design. Test: two studies at means 0 and 5 with a one-hot protect by study; after the fix both study means are within
0.3 of 2.5 and equal to the result without protect (version 1 gave a difference of 6.8 between the study means in this test).

**2. Sex control (commit 564f9e5, with fixes f268a12 and c0cd0bd).** The third feature set
`forced_sex_genes` (see "Control runs"). `LosoConfig.drop_incomplete_features` drops, per fold, features with a missing
value in a training or held-out sample.

**3. chorus harmoniser (commit 34c59fb).**
- **Number of factors.** Chosen by parallel analysis on the null gene block of the training pool: squared singular
  values of the centred null block against the 0.95 quantile of those from 20 copies with each column permuted
  independently (seed 0); k is the number of leading components above it, and can be zero (no factor removed). Version 1
  used the number of training studies plus 4 (`k_method: studies_plus`, still available).
- **Protected marker scores** (option `protect_scores`, off by default and not used in the version 2 reruns; used by
  `chorus_protected` since version 3): named scores
  (`neutrophil`, `lymphocyte`, with the marker genes of the harmonisation diagnostics, or a dict of gene lists) enter
  the joint loading fit alongside the outcome and are dropped at transform time. A score is the mean of its genes
  standardised across the training pool. Its genes are left out of the null set, and the factors and the parallel
  analysis use the null block after the scores are projected out of it. Set it in a config with
  `benchmark: harmoniser_params: {chorus: {protect_scores: [neutrophil, lymphocyte]}}`.

**Reruns.**

| folder | what was rerun | kept as |
|---|---|---|
| `lc_transcriptomics_blood` | rows `combat_addon`, `chorus` and `nested` (`chorus run --only-harmonisers`); other rows, `within_study_cv.csv` and `within_study_cv_folds.csv` unchanged | `lc_transcriptomics_blood_v1` (full copy) |
| `lc_transcriptomics_blood_noGSE270045` | the same | `lc_transcriptomics_blood_noGSE270045_v1` |
| `lc_acute_transfer` | the same | `lc_acute_transfer_v1` |
| `simulation_factors`, `simulation_shifts_only` | whole run (no protect there, so only `chorus` and `nested` rows can change) | `simulation_factors_v1`, `simulation_shifts_only_v1` |
| `control_sex_label` | whole run, now with three feature sets | `control_sex_label_v1` |
| `lc_transcriptomics_blood_harmonisation` | whole run | `lc_transcriptomics_blood_harmonisation_v1` |

- `--only-harmonisers` reruns the named rows with the same folds and seeds; `nested` still chooses among all harmonisers
  of the config. The rerun rows replace the old ones in `loso_results.csv`, `loso_predictions.csv` and
  `selected_features.json`, and `loso_summary.csv` is recomputed. Rows of the other harmonisers were checked to be
  identical to `_v1` in all five rerun folders and in `control_sex_label` (identity and center_scale).
- Factors chosen by parallel analysis (`n_factors`): primary benchmark 5 to 8 per fold (version 1: 8); noGSE270045 4 to
  8 (version 1: 7); acute transfer 7 (version 1: 9); harmonisation diagnostics 5 to 8; sex control 4 to 8;
  `simulation_factors` 18 and 19 (the simulation has 3 sample-level factors and 6 studies); `simulation_shifts_only` 4.
- No fit failed in any rerun.
- **Not rerun:** `enrichment.csv` (it reads the version 1 `selected_features.json`, whose `combat_addon`, `chorus` and
  `nested` entries have since been replaced), `control_label_permutation` (uses `chorus`), `concordance` and
  `lc_transcriptomics_blood_within_v2` (no harmoniser), and the smoke runs. `enrichment.csv` and
  `control_label_permutation` were rerun in version 3.

## Version 3 changes (13 November 2026)

Two code changes, after which the affected results were rerun. The version 2 results they replace are in git at commit
44dffee. The `_v1` folders are unchanged.

**1. AddonComBat (commit 34b2eec).**
- **Diagnosis.** After version 2, the `combat_addon` training pool kept study structure in the harmonisation
  diagnostics: pool silhouette by study 0.21 to 0.26 and other-study neighbour fraction 0.07 to 0.12, against -0.03 to
  -0.02 and 0.57 to 0.63 for `center_scale`. On the fold with GSE226260_WB held out, the version 2 operator was refitted
  with parts of it switched off (scratch script, not kept). The columns are:
  - the maximum over training studies of |study mean minus pool mean|, median over genes, after fit and transform;
  - the pool silhouette by study;
  - the other-study neighbour fraction.

  | variant (version 2 code) | max abs study mean deviation | silhouette_study | knn_other_study |
  |---|---|---|---|
  | as run: protect `sample_matrix` and label in the design, `shrink` 0.5 | 1.378 | 0.229 | 0.102 |
  | label only in the design, `shrink` 0.5 | 1.378 | 0.229 | 0.102 |
  | protect and label, `shrink` 0 | 0.637 | 0.058 | 0.400 |
  | label only, `shrink` 0 | 0.637 | 0.058 | 0.400 |
  | no design, `shrink` 0.5 | 1.326 | 0.208 | 0.123 |
  | no design, `shrink` 0 | 0.000 | -0.033 | 0.605 |
  | `center_scale` | 0.000 | -0.033 | 0.605 |

  - **Standardisation:** the training studies were shifted and scaled to the grand mean and pooled variance per gene of
    the training pool.
  - **Shrinkage:** each study's per gene location deviation from the grand mean was replaced by 0.5 times itself plus
    0.5 times its average over genes. The log variance ratio was shrunk the same way. So half of each gene's study
    offset, less its across-gene average, was kept.
  - **Covariates:** `sample_matrix` is constant within every cohort and was already dropped from the design by the
    version 2 fix (the rows with and without it are identical). The label stayed in the design uncentred, so the
    location was the intercept at label 0 (the control mean), and studies with different case fractions kept different
    means.
- **Fix.**
  - The label and protected covariates are centred within each training study before the least-squares fit, so the
    intercept is the study mean. They still remove their variation from the variance estimate.
  - `shrink` defaults to 0 (it was 0.5). The option is kept.
  - The held-out study is unchanged: location and scale come from that study alone, without covariates.
  - On the same fold the fixed operator gives max deviation 0.000, silhouette -0.043 and other-study fraction 0.610.
- **Test** (`test_combat_addon_training_study_means_equal_after_transform`):
  - Setup: three studies with per gene offsets and scales, case fractions 0.2, 0.5 and 0.8, a label effect on half the
    genes, and a protected covariate that varies within studies.
  - After fit and transform, the per study gene means equal the pool means within 1e-8, with and without protect.
  - With `shrink` 0.5 they differ by more than 0.1.
  - It fails on the version 2 code (largest deviation 6.2).

**2. chorus_protected (commit bcf509d).** A second chorus variant, registered as harmoniser `chorus_protected`:
`InvarianceAnchored(protect_scores=["neutrophil", "lymphocyte"])`.
- Its other settings are those of `chorus` (factor count by parallel analysis).
- Its marker genes are those of the harmonisation diagnostics. All 10 are among the filtered genes of the three
  benchmarks.
- It was added to `harmonisers` in the three benchmark configs and in the sex control grid. The harmonisation
  diagnostics use the config list, so they include it.
- New config field `nested_candidates` (`LosoConfig.nested_candidates`): the harmonisers that `nested` chooses among,
  default `harmonisers`. The three configs set it to the version 2 list, so `chorus_protected` is not a nested
  candidate.
- It is not in the simulations. The simulated genes (`G00000` to `G01999`) contain no marker genes, so the protected
  scores cannot be formed (the operator raises an error when a named score has no genes). The simulation grid is
  unchanged.

**3. --only-harmonisers merge (commit 2b58d1b).** A harmoniser new to the results is placed in the `run_loso` loop order,
here `chorus_protected` between `chorus` and `nested`.

**Reruns.**

| folder | what was rerun |
|---|---|
| `lc_transcriptomics_blood`, `lc_transcriptomics_blood_noGSE270045`, `lc_acute_transfer` | rows `combat_addon`, `chorus_protected` (new) and `nested` (`--only-harmonisers combat_addon,chorus_protected,nested`); other rows and the within-study files unchanged |
| `simulation_factors`, `simulation_shifts_only` | whole run (only `combat_addon` and `nested` rows can change) |
| `control_sex_label` | whole run, now with `chorus_protected` |
| `lc_transcriptomics_blood_harmonisation` | whole run, now with `chorus_protected` |
| `enrichment.csv` | whole run on the current primary `selected_features.json` |
| `control_label_permutation` | whole run (identity and the current `chorus`) |

- **Unchanged rows checked against version 2:**
  - In the three benchmarks, both simulations and `control_sex_label`, the rows of the harmonisers not rerun are
    identical.
  - In the harmonisation diagnostics, every harmoniser other than `combat_addon` gives the same values.
  - In `control_label_permutation` the identity rows are identical.
- **enrichment.csv:** compared with the file computed before version 2, the `chorus` and `nested` rows changed (4 each),
  the `combat_addon` rows are unchanged, and the 15 `chorus_protected` rows are new.
- No fit failed in any rerun.
- **Factors chosen by parallel analysis for `chorus_protected` (`n_factors`):** primary benchmark 5 to 9 per fold;
  noGSE270045 4 to 10; acute transfer 5; harmonisation diagnostics 5 to 9.
- **Not rerun:** `concordance` and `lc_transcriptomics_blood_within_v2` (no harmoniser), and the smoke runs.

## Code changes after version 3 without reruns (14 November 2026)

- **Sex control, forced sex genes (commit efcf77b).**
  - **Why the rows were equal.** In `control_sex_label`, the `forced_sex_genes` rows equal the `all_genes` rows for the
    invariant selector only. The forced genes reached the candidate features, but the invariant selector dropped them
    (I squared 92 to 97 against its limit of 50), so the selected features were the same.
  - **Variance selector.** With identity, center_scale and combat_addon, the variance selector already kept all 8 forced
    genes. With chorus and chorus_protected it kept none.
  - **Fix.** The forced genes present in a fold are now appended to every selector's output (`LosoConfig.forced_features`,
    result column `n_forced_added`).
  - **Status.** `control_sex_label` was not rerun on 14 November. It was rerun on 15 November; see "Proteomics (15
    November 2026)", item 6.
- **Layer adapters (commit b5083c4).** The transcriptomics loader now goes through the `bulk_rna` adapter. The loaded
  matrices of `configs/lc_transcriptomics_blood.yaml` and `configs/lc_acute_transfer.yaml` were checked to be identical
  (hash of every cohort matrix) before and after the change.

## Proteomics (15 November 2026)

**1. Curation applied first** (in the builders; `docs/proteomics_cohorts.md`, "Curation decisions"):
- **Cervia2024:** the 5 healthy controls that match later COVID-19 subjects are linked to them
  (`linked_healthy_subject_id`), and their healthy samples are dropped (SomaScan: 39 to 34 healthy).
- **Chowdhury2024:** the authors' labels are kept, the controls have `control_type` convalescent_unverified, and the 4
  unlabelled runs are excluded.
- **Zoodsma2022:** the 16 Nijmegen patients stay excluded.
- **Captur2022:** not used.

The schema gains `control_type` convalescent_unverified and the pooled layer `proteomics`.

**Loader options added for these configs** (`chorus.cli.load_config`, `preprocessing`). Transcriptomics configs do not
use them.

| option | effect |
|---|---|
| `drop_sample_phases` | drops samples whose `sample_phase` is listed (Test 1: `acute`, because acute draws carry the later outcome in `lc_status`) |
| `keep_lc_status` | keeps only the listed `lc_status` values, before one sample per subject |
| `timepoint_order` | timepoint labels from earliest to latest, for one sample per subject. Chowdhury2024 keeps the pre-intervention sample, and Su2022 and Zoodsma keep the earliest acute draw. |
| `relabel_acute` | `lc_status` becomes `acute` for every sample with `sample_phase` acute (the curated value is kept in `lc_status_original`) |
| `impute` | per cohort and protein, from that cohort's own retained values, without labels. `study_median` for affinity (values masked at QC); `study_min` for MS (missing mostly below detection). A cohort entry can set its own `impute`. |
| `features: union` | every cohort keeps its own proteins (pairwise concordance); the default is the proteins shared by all cohorts |

- **Shared proteins.** For a cohort that is not counts, a protein counts as present when it has at least one value
  among the retained samples.
- **Other choices.** No covariates and no `protect` are used. Repeated runs of one Buggert2025 sample (20 subjects) are
  reduced to the first in file order.

**2. Test 1, LC versus COVID_recovered** (`lc_proteomics_affinity`, `lc_proteomics_ms`).
- **Samples:** one post-acute sample per subject. Leave-one-study-out within each block, on the proteins shared by all
  three cohorts of the block.
- **Grid:**
  - harmonisers identity, center_scale, rank, fsqn, combat_addon, chorus, and nested selection over these six;
  - selectors variance (k 100), differential (k 50), invariant (k 50);
  - classifiers logreg and rf;
  - 1000 bootstrap resamples;
  - within-study CV with the per-fold method (5 folds x 10 repeats, `--within`).
- **Selector sizes:** k was chosen for these blocks; the transcriptomics grid uses 500/200/200 on 10,442 genes.
- **Failed fits:** none, in 126 rows per block.

| block | cohort | samples (LC / COVID_recovered) | proteins used | values imputed |
|---|---|---|---|---|
| affinity | Buggert2025 | 197 (146 / 51) | 376 | 642 (median) |
| affinity | Cervia2024_SomaScan | 113 (40 / 73) | 376 | 0 |
| affinity | Su2022 | 125 (75 / 50) | 376 | 0 |
| MS | Cervia2024_MS | 113 (40 / 73) | 182 | 907 (minimum) |
| MS | Chowdhury2024 | 30 (21 / 9) | 182 | 56 (minimum) |
| MS | Kovarik2023_plasma | 26 (13 / 13) | 182 | 646 (minimum) |

The per-protein observed fractions are in `features_observed_fraction.csv`. In the MS block, 49 (Cervia2024_MS), 11
(Chowdhury2024) and 69 (Kovarik2023_plasma) of the 182 proteins have at least one missing value; the lowest observed
fractions are 0.44, 0.20 and 0.08.

**Su2022 and the harmonisers.** The Su2022 values are the authors' Z scores: batch corrected, adjusted for age, sex and
BMI, and scaled by the mean and standard deviation of the residuals of the matched healthy controls. The other two
affinity cohorts are NPX (Buggert2025) and log2 RFU (Cervia2024_SomaScan). What each harmoniser does to Su2022:

| harmoniser | on Su2022 | meaningful for Su2022 |
|---|---|---|
| identity | pools the Z scores with NPX and log2 RFU on their own scales | no: the units differ by construction |
| center_scale | per protein z-score within Su2022's retained samples, a linear map of the authors' values; the same operation as on the other cohorts | yes. The authors' age, sex and BMI adjustment stays, and is not applied to the other cohorts. |
| combat_addon | per protein location and scale of the study mapped to the training pool's grand mean and pooled variance | yes, as center_scale. The pool's grand mean and variance mix the three units. |
| fsqn | maps each protein's within-study ranks onto the training pool's quantiles of that protein | the within-study ranks are unaffected by the authors' per-protein transform. The reference quantiles are taken from the raw training values, so they mix units whenever cohorts with different units are in training (any affinity fold). |
| rank | within-sample rank across proteins | no: in Su2022 it orders proteins by deviation from healthy controls; in NPX and RFU cohorts, by measured level |
| chorus | centring by the training grand mean, factors from the null proteins of the pooled training values, then a per-study z-score | partly. The final per-study z-score is as center_scale. The factors are estimated on the pool before any per-study scaling, so they are fitted on mixed units. |
| nested | chooses one of the above per fold | as the harmoniser chosen (`harmoniser_chosen`) |

**3. Diagnostics and concordance.**
- **Harmonisation diagnostics:** `lc_proteomics_affinity/harmonisation` and `lc_proteomics_ms/harmonisation`, same
  metrics as for transcriptomics.
  - There are no sex rows (the cohorts have no `inferred_sex`).
  - The marker rows have `n_genes` 0 and empty values, because the marker genes are Ensembl ids, not UniProt.
  - n_pcs is 50.
  - No failed fit.
- **Concordance per block:** `lc_proteomics_affinity/concordance` and `lc_proteomics_ms/concordance`, on the 376 and
  182 shared proteins after the same preprocessing and imputation as the benchmark. `--top 50`.
  - Symbols come from the HGNC `uniprot_ids` column; a multi-protein key gets its symbols joined by "_", and an
    unmapped key has no symbol.
- **Concordance across all six cohorts:** `concordance_proteomics_all`, from `configs/concordance_proteomics_all.yaml`.
  - Same sample selection, with imputation per cohort (median for affinity, minimum for MS).
  - Each cohort is kept on its own proteins (`features: union`, 7,419 in all; 38 without a symbol), and each pair is
    correlated over the proteins both have.
  - Deduplication is off, because Cervia2024_SomaScan and Cervia2024_MS are the same people.
  - `n_shared.csv` gives the number of proteins per pair. Su2022 shares 4 or 5 proteins with each MS cohort, so those
    three correlations rest on 4 or 5 values; they are flagged `fewer_than_10_shared` in `pairs.csv`.
  - `spearman_top50.csv`: the 50 proteins with the largest mean |g| over the union of proteins; no pair shares 10 of
    them (0 or 1, `n_shared_top50.csv`), so every off-diagonal value is NaN.
  - Chowdhury2024 has 404 proteins here: one of its 405 has no value in the 30 retained samples.

**4. Positive control in proteomics, acute versus never_infected** (`control_acute_proteomics`).
- **Label:** acute = `sample_phase` acute, whatever the later outcome; never_infected = curated `lc_status`. The earliest
  acute draw per subject is kept.
- **Grid:** the six harmonisers plus nested; selectors variance (k 100) and invariant (k 50); logreg; within-study CV.
- **Cohorts:** those with both classes after one sample per subject and removal of subjects measured in two cohorts.
  - **Zoodsma2022_Breda is disabled.** Its 5 never infected samples are bridge samples of Hannover healthy donors
    already in Zoodsma2022_Hannover (TK_10, TK_33, TK_55, TK_63, TK_65), so it keeps acute samples only.
  - **Zoodsma2022_Nijmegen** has no never infected sample.
  - **Zoodsma2022_Hannover's acute samples** are 13 bridge samples (12 patients) of the Nijmegen hospitalised cohort.

| cohort | acute | never_infected | proteins used | values imputed |
|---|---|---|---|---|
| Cervia2024_SomaScan | 113 | 34 | 163 | 0 |
| Su2022 | 204 | 289 | 163 | 0 |
| Zoodsma2022_Hannover | 12 | 62 | 163 | 82 |

- **chorus failed in every fold** (6 error rows, both selectors): "only 0 null genes passed the screen" (Cervia2024_SomaScan
  held out), "only 6" (Su2022) and "only 4" (Zoodsma2022_Hannover). The default null screen of `InvarianceAnchored` was
  not changed.
- **Why the chorus rows are absent (checked 15 November).** They are present in `loso_results.csv` as error rows
  (column `error`, no metrics). They were absent from `loso_summary.csv` and `loso_predictions.csv`, which kept only
  fits that succeeded. Since 16 November `loso_summary.csv` keeps them: one row per selector with classifier
  `not_fitted`, `n_studies` 0, `n_errors` 3 and the messages in `error` (the summary was rewritten from
  `loso_results.csv`; the other rows are unchanged). `loso_predictions.csv` still has no chorus rows, as no
  prediction was made.
  - **The screen.** `InvarianceAnchored` needs at least 20 null proteins. A null protein has a meta-analytic p above
    0.2, an absolute Hedges g below 0.35 in every training cohort, I squared at most 50, and training variance at or
    above the 20th percentile.
  - **Per fold,** on the 163 proteins of the acute contrast:

    | held out | p above 0.2 | max abs g below 0.35 | I squared at most 50 | variance ok | all four | needed |
    |---|---|---|---|---|---|---|
    | Cervia2024_SomaScan | 40 | 5 | 40 | 130 | 0 | 20 |
    | Su2022 | 85 | 17 | 47 | 130 | 6 | 20 |
    | Zoodsma2022_Hannover | 73 | 14 | 29 | 130 | 4 | 20 |

  - **Decision.** The failure is the operator refusing to fit, not a bug, so it is documented and the screen
    thresholds were not relaxed. Relaxing them would make "null" proteins that are associated with the label in the
    training cohorts.
- `nested` chooses among all six candidates. A candidate that fails in an inner fold scores missing there.

**5. Positive control in transcriptomics, acute versus never_infected** (`control_acute_transcriptomics/cohorts.csv`).
- **Survey:** every transcriptomics cohort on disk was checked for acute COVID-19 samples and never infected or healthy
  samples. Only GSE215865 has both: 1,218 libraries (424 subjects) of COVID-19 positive admissions, and 151 libraries
  (64 subjects) of COVID-19 negative admissions curated as never_infected. Those are acutely ill patients without
  COVID-19, not healthy donors.
- **GSE226260 acute set (p20169):** 104 acute samples with a severity. The other 38 samples (36 subjects) have
  "pasc severity: N/A" and no other field. The paper reports 35 uninfected participants, but GEO labels no sample as
  uninfected.
- **Other cohorts:** GSE169687 (14 never infected) and GSE251849 (7) have no acute samples. GSE222253 (11 in hospital)
  has no uninfected samples.
- **Not run:** at most one cohort qualifies (two if the GSE226260 N/A samples were taken as uninfected), against the
  three required.

**6. Sex control rerun** (`control_sex_label`, with `LosoConfig.forced_features`).
- **Forced genes added:** `n_forced_added` is 8 in every `forced_sex_genes` row except the variance selector with
  identity, center_scale and combat_addon, where it is 0 because the selector already chose all 8. The forced genes
  therefore reach the model in every row.
- **Rows equal to all_genes:** 2 of the 40 `forced_sex_genes` rows have the same AUROC as `all_genes` (GSE224615 held
  out, invariant selector, center_scale and combat_addon).
- **all_genes and no_xy:** identical to the version 3 files except `seconds` and the new `n_forced_added` column (0).
- **Earlier statement superseded:** the statement under "Control runs" that the invariant rows of `forced_sex_genes`
  equal those of `all_genes` describes the version 3 files.

**Files not committed:** the logs (`results/<run>.log`, the diagnostics and concordance logs named with `_` in place of
`/`) and `selected_features.json`.

## Test 2: layers on the same people (15 November 2026)

**Paired sets** (`docs/paired_cohorts.md`). A sample is paired when its `sample_id` is in every layer and its labels
agree across the layers' metadata. One sample per subject is kept: the earliest by timepoint.

| cohort | layers (features used) | LC vs COVID_recovered | acute vs never_infected |
|---|---|---|---|
| Cervia2024 | SomaScan (6402), MS (301 of 308; 7 dropped as observed in fewer than half the samples) | 40 / 73 | not run: the MS layer has no acute or healthy sample |
| Su2022 | Olink (441), metabolomics (891, `data/metabolomics/Su2022_metabolomics`, Table S2.2) | 68 / 45 (T3) | 198 / 27 (earliest paired acute draw; healthy donors in both sheets) |
| Kovarik2023 | not built | | |

**Kovarik2023 not built.** The supplement was downloaded and checked, but its tables share only column labels
(`healthy_<n>` and the like). Those labels do not link the triglyceride table (S4) to the Biocrates table (S3) on the
126 triglycerides they share, so they cannot be relied on to link the eicosanoid and Biocrates tables to the proteomics
runs (`docs/paired_cohorts.md`).

**Models** (`chorus.eval.paired.PairedCV`). Repeated stratified 5-fold CV, 10 repeats, split seed `random_state + r`;
every model is fitted and scored on the same folds.
- **Missing values:** filled inside each fold from that fold's training samples (MS: per protein minimum; the other
  layers: median). Features observed in fewer than half of the cohort's paired samples are dropped first (label free).
- **Classifier:** logreg (the CHORUS `logreg`: standardised, C 0.5, balanced classes) for every model except DIABLO.

| model | per fold |
|---|---|
| `layer:<name>` | features standardised on the training fold, top 50 by absolute training effect (differential selector), logreg |
| `early_fusion` | the per-layer selections concatenated, each layer divided by the square root of its number of selected features |
| `late_fusion` | `chorus.integrate.mosaic.LateFusion` fitted on inner 5-fold out-of-fold layer probabilities of the training fold, applied to the layer models refitted on the whole training fold |
| `mofa` | MOFA+ (mofapy2 0.7.5), 10 factors, one group, layers as views (feature names prefixed by the layer), all features, fitted on the training fold; test samples projected onto the training loadings (`mofa_factors`); logreg on the factors |
| `diablo` | mixOmics 6.32.0 `block.splsda` in R 4.5.3, 2 components, keepX 20 per block and component, design 0.1 off the diagonal, `scale = TRUE`; the score is the positive-class column of `predict(...)$WeightedPredict` at component 2, clipped to [0, 1]. It is not a probability, so its Brier score and ECE are those of this score. |

**Agreement between layers.**
- `effect_correlation.csv` (Cervia2024 only, the one pair with a shared feature space): Spearman correlation of the
  Hedges g (LC minus COVID_recovered, all 113 paired samples) between SomaScan and MS over the UniProt keys in both
  layers (173).
- `probability_correlation.csv`: per repeat and pair of models, Spearman and Pearson correlation across subjects of
  the out-of-fold probabilities. Every subject is predicted once per repeat. The pairs include the layer models with
  each other and with the fusion models.

**Positive control:** acute versus never_infected, the same models and tables, in Su2022 only (the only paired cohort
with both classes).

**Failures:** none. Every model has 50 folds in each run (`n_failed_folds` 0).

**Column glossary.**
- `fold_metrics.csv`: `model`, `repeat`, `fold`, `n` (test samples), `prevalence`, and the metrics of
  loso_results.csv on that fold; `error` is filled on a failed fold.
- `summary.csv`: per model, the mean over folds of each metric, `auroc_sd`, `n_folds` and `n_failed_folds`.
- `predictions.csv`: `sample`, `repeat`, `fold`, `model`, `y`, `p`, `error`.
- `samples.csv`: `sample_id`, `subject_id`, `y`, `timepoint`, `label`.
- `features.csv`: `layer`, `n_features_cohort`, `n_dropped_missing`, `n_used`, `impute`.
- `effects.csv`: Hedges g per feature (rows) and layer (columns).
- `effect_correlation.csv`: `layer_a`, `layer_b`, `n_shared`, `spearman`, `p`.

**Not committed:** the logs (`results/paired_<cohort>.log`).

## Test 3: mosaic programmes (15 November 2026)

**Method:** `docs/mosaic.md`.
- **Programmes:** MSigDB 2026.2 Hallmark (50) and Reactome (1,867), CC BY 4.0, mapped by Entrez id through HGNC to Ensembl and UniProt.
- **Scores:** in each cohort, the mean of the members standardised within that cohort, with at least 5 measured members.
- **Environment:** cohort and layer.
- **Harmonisers and classifier:** identity and center_scale, with logreg.
- **Run order:** the positive control (acute versus never_infected) ran first.

**Gene set mapping** (`genesets.csv`). Hallmark: 7,322 gene memberships, 7,321 found in HGNC with an Ensembl id, 7,317 with a UniProt accession. Reactome: 105,477 memberships, 105,472 in HGNC, 105,464 with an Ensembl id, 105,307 with a UniProt accession.

### Positive control: acute versus never_infected (`acute/`)

**Cohorts and coverage** (`cohorts.csv`, `coverage.csv`):

| environment | role | positive / negative | features measured | Hallmark scored (of 50) | Reactome scored (of 1,867) | median fraction of members measured (scored programmes) |
|---|---|---|---|---|---|---|
| rna:GSE215865 | primary | 424 / 64 | 15173 | 50 | 1628 | 1.00 |
| protein:Cervia2024_SomaScan | primary | 113 / 34 | 6402 | 50 | 1540 | 1.00 |
| protein:Su2022 | primary | 204 / 289 | 441 | 33 | 177 | 0.11 |
| protein:Zoodsma2022_Hannover | primary | 12 / 62 | 359 | 26 | 122 | 0.12 |

**Concordance of programme effects** (`concordance_all.csv`; Spearman, with the number of programmes scored in both in parentheses; Hallmark and Reactome alone in `concordance_hallmark.csv` and `concordance_reactome.csv`):

| | rna:GSE215865 | protein:Cervia2024_SomaScan | protein:Su2022 | protein:Zoodsma2022_Hannover |
|---|---|---|---|---|
| rna:GSE215865 |  | -0.02 (1486) | 0.15 (210) | 0.13 (148) |
| protein:Cervia2024_SomaScan | -0.02 (1486) |  | 0.27 (210) | 0.31 (148) |
| protein:Su2022 | 0.15 (210) | 0.27 (210) |  | 0.37 (123) |
| protein:Zoodsma2022_Hannover | 0.13 (148) | 0.31 (148) | 0.37 (123) |  |

**Invariance partition** (`partition.csv`):

| analysis | environments | programmes partitioned (scored in every environment) | invariant | null | unstable | none | not_partitioned |
|---|---|---|---|---|---|---|---|
| primary | 4 | 123 | 0 | 0 | 123 | 0 | 1794 |

**Cross-layer leave one cohort out** (`loso.csv`, logreg; bootstrap 95% interval, 1000 resamples):

| scheme | test | role | positive / negative | programmes | AUROC identity (95% CI) | AUROC center_scale (95% CI) |
|---|---|---|---|---|---|---|
| protein_to_rna | rna:GSE215865 | primary | 424 / 64 | 123 | 0.637 (0.579 to 0.695) | 0.681 (0.626 to 0.733) |
| rna_plus_other_protein | protein:Cervia2024_SomaScan | primary | 113 / 34 | 123 | 0.472 (0.364 to 0.571) | 0.576 (0.471 to 0.683) |
| rna_plus_other_protein | protein:Su2022 | primary | 204 / 289 | 123 | 0.910 (0.880 to 0.937) | 0.829 (0.786 to 0.865) |
| rna_plus_other_protein | protein:Zoodsma2022_Hannover | primary | 12 / 62 | 123 | 0.996 (0.984 to 1.000) | 0.960 (0.909 to 0.993) |

### LC versus COVID_recovered (`lc/`)

**Cohorts and coverage:**

| environment | role | positive / negative | features measured | Hallmark scored (of 50) | Reactome scored (of 1,867) | median fraction of members measured (scored programmes) |
|---|---|---|---|---|---|---|
| rna:GSE226260_WB | primary | 17 / 20 | 12801 | 50 | 1617 | 0.89 |
| rna:GSE226260_PBMC | primary | 28 / 8 | 13866 | 50 | 1648 | 0.92 |
| rna:GSE224615 | primary | 23 / 13 | 15734 | 50 | 1660 | 0.93 |
| rna:GSE251849 | primary | 11 / 5 | 13872 | 50 | 1639 | 0.92 |
| rna:GSE270045 | primary | 19 / 17 | 15561 | 50 | 1739 | 0.97 |
| protein:Buggert2025 | primary | 146 / 51 | 2903 | 50 | 967 | 0.39 |
| protein:Cervia2024_SomaScan | primary | 40 / 73 | 6402 | 50 | 1540 | 0.91 |
| protein:Su2022 | primary | 75 / 50 | 441 | 33 | 177 | 0.09 |
| protein:Chowdhury2024 | primary | 21 / 9 | 404 | 22 | 134 | 0.14 |
| protein:Kovarik2023_plasma | primary | 13 / 13 | 314 | 19 | 144 | 0.12 |
| protein:Cervia2024_MS | sensitivity | 40 / 73 | 308 | 19 | 113 | 0.16 |

**Concordance of programme effects** (`concordance_all.csv`, as above; Cervia2024_MS is the sensitivity column, the same sera as Cervia2024_SomaScan):

| | rna:GSE226260_WB | rna:GSE226260_PBMC | rna:GSE224615 | rna:GSE251849 | rna:GSE270045 | protein:Buggert2025 | protein:Cervia2024_SomaScan | protein:Su2022 | protein:Chowdhury2024 | protein:Kovarik2023_plasma | protein:Cervia2024_MS |
|---|---|---|---|---|---|---|---|---|---|---|---|
| rna:GSE226260_WB |  | -0.06 (1662) | 0.25 (1665) | 0.21 (1657) | 0.19 (1661) | -0.04 (966) | 0.23 (1481) | 0.11 (210) | 0.38 (145) | -0.01 (153) | 0.35 (121) |
| rna:GSE226260_PBMC | -0.06 (1662) |  | -0.03 (1685) | -0.17 (1677) | -0.13 (1691) | -0.23 (976) | -0.08 (1501) | -0.13 (210) | 0.33 (148) | -0.26 (156) | 0.35 (124) |
| rna:GSE224615 | 0.25 (1665) | -0.03 (1685) |  | 0.24 (1681) | -0.08 (1701) | -0.07 (977) | 0.21 (1506) | -0.04 (210) | 0.50 (147) | -0.07 (156) | 0.60 (124) |
| rna:GSE251849 | 0.21 (1657) | -0.17 (1677) | 0.24 (1681) |  | -0.34 (1680) | -0.15 (974) | 0.16 (1496) | -0.07 (210) | -0.18 (149) | 0.04 (158) | -0.20 (126) |
| rna:GSE270045 | 0.19 (1661) | -0.13 (1691) | -0.08 (1701) | -0.34 (1680) |  | 0.07 (991) | 0.05 (1544) | 0.24 (210) | 0.30 (151) | 0.05 (159) | 0.31 (127) |
| protein:Buggert2025 | -0.04 (966) | -0.23 (976) | -0.07 (977) | -0.15 (974) | 0.07 (991) |  | -0.04 (1011) | -0.08 (210) | 0.16 (152) | -0.27 (160) | 0.17 (129) |
| protein:Cervia2024_SomaScan | 0.23 (1481) | -0.08 (1501) | 0.21 (1506) | 0.16 (1496) | 0.05 (1544) | -0.04 (1011) |  | 0.08 (210) | -0.01 (156) | 0.25 (163) | 0.21 (132) |
| protein:Su2022 | 0.11 (210) | -0.13 (210) | -0.04 (210) | -0.07 (210) | 0.24 (210) | -0.08 (210) | 0.08 (210) |  | 0.01 (86) | 0.03 (88) | -0.10 (73) |
| protein:Chowdhury2024 | 0.38 (145) | 0.33 (148) | 0.50 (147) | -0.18 (149) | 0.30 (151) | 0.16 (152) | -0.01 (156) | 0.01 (86) |  | -0.09 (139) | 0.72 (125) |
| protein:Kovarik2023_plasma | -0.01 (153) | -0.26 (156) | -0.07 (156) | 0.04 (158) | 0.05 (159) | -0.27 (160) | 0.25 (163) | 0.03 (88) | -0.09 (139) |  | -0.19 (124) |
| protein:Cervia2024_MS | 0.35 (121) | 0.35 (124) | 0.60 (124) | -0.20 (126) | 0.31 (127) | 0.17 (129) | 0.21 (132) | -0.10 (73) | 0.72 (125) | -0.19 (124) |  |

**Invariance partition** (`partition.csv`; `sensitivity_Cervia2024_MS` replaces Cervia2024_SomaScan with Cervia2024_MS):

| analysis | environments | programmes partitioned (scored in every environment) | invariant | null | unstable | none | not_partitioned |
|---|---|---|---|---|---|---|---|
| primary | 10 | 77 | 0 | 0 | 77 | 0 | 1840 |
| sensitivity_Cervia2024_MS | 10 | 67 | 0 | 0 | 67 | 0 | 1850 |

**Cross-layer leave one cohort out** (`loso.csv`). Cervia2024_MS is a test cohort only, under `rna_to_protein`:

| scheme | test | role | positive / negative | programmes | AUROC identity (95% CI) | AUROC center_scale (95% CI) |
|---|---|---|---|---|---|---|
| rna_to_protein | protein:Buggert2025 | primary | 146 / 51 | 961 | 0.568 (0.487 to 0.654) | 0.428 (0.348 to 0.513) |
| rna_to_protein | protein:Cervia2024_SomaScan | primary | 40 / 73 | 1470 | 0.536 (0.429 to 0.649) | 0.497 (0.382 to 0.608) |
| rna_to_protein | protein:Su2022 | primary | 75 / 50 | 210 | 0.659 (0.565 to 0.756) | 0.651 (0.561 to 0.747) |
| rna_to_protein | protein:Chowdhury2024 | primary | 21 / 9 | 144 | 0.614 (0.344 to 0.847) | 0.471 (0.238 to 0.704) |
| rna_to_protein | protein:Kovarik2023_plasma | primary | 13 / 13 | 153 | 0.562 (0.325 to 0.799) | 0.491 (0.260 to 0.728) |
| rna_to_protein | protein:Cervia2024_MS | sensitivity | 40 / 73 | 121 | 0.523 (0.410 to 0.638) | 0.468 (0.352 to 0.583) |
| protein_to_rna | rna:GSE226260_WB | primary | 17 / 20 | 77 | 0.682 (0.476 to 0.844) | 0.638 (0.453 to 0.806) |
| protein_to_rna | rna:GSE226260_PBMC | primary | 28 / 8 | 77 | 0.518 (0.295 to 0.746) | 0.473 (0.254 to 0.705) |
| protein_to_rna | rna:GSE224615 | primary | 23 / 13 | 77 | 0.548 (0.334 to 0.763) | 0.515 (0.294 to 0.739) |
| protein_to_rna | rna:GSE251849 | primary | 11 / 5 | 77 | 0.364 (0.109 to 0.655) | 0.309 (0.054 to 0.600) |
| protein_to_rna | rna:GSE270045 | primary | 19 / 17 | 77 | 0.851 (0.703 to 0.963) | 0.845 (0.709 to 0.954) |
| leave_one_cohort_out | rna:GSE226260_WB | primary | 17 / 20 | 77 | 0.688 (0.494 to 0.841) | 0.647 (0.459 to 0.812) |
| leave_one_cohort_out | rna:GSE226260_PBMC | primary | 28 / 8 | 77 | 0.464 (0.268 to 0.679) | 0.469 (0.254 to 0.705) |
| leave_one_cohort_out | rna:GSE224615 | primary | 23 / 13 | 77 | 0.552 (0.324 to 0.763) | 0.532 (0.311 to 0.746) |
| leave_one_cohort_out | rna:GSE251849 | primary | 11 / 5 | 77 | 0.400 (0.091 to 0.710) | 0.345 (0.073 to 0.655) |
| leave_one_cohort_out | rna:GSE270045 | primary | 19 / 17 | 77 | 0.789 (0.625 to 0.926) | 0.715 (0.542 to 0.870) |
| leave_one_cohort_out | protein:Buggert2025 | primary | 146 / 51 | 77 | 0.578 (0.485 to 0.667) | 0.551 (0.456 to 0.631) |
| leave_one_cohort_out | protein:Cervia2024_SomaScan | primary | 40 / 73 | 77 | 0.537 (0.427 to 0.642) | 0.554 (0.447 to 0.661) |
| leave_one_cohort_out | protein:Su2022 | primary | 75 / 50 | 77 | 0.512 (0.407 to 0.618) | 0.478 (0.377 to 0.583) |
| leave_one_cohort_out | protein:Chowdhury2024 | primary | 21 / 9 | 77 | 0.545 (0.338 to 0.730) | 0.529 (0.323 to 0.730) |
| leave_one_cohort_out | protein:Kovarik2023_plasma | primary | 13 / 13 | 77 | 0.704 (0.485 to 0.893) | 0.781 (0.568 to 0.953) |

**Bridge check in Cervia2024** (`bridge_cervia.csv`). SomaScan and MS programme scores of the same 113 sera, Spearman correlation per programme across the sera; `own`: each layer's own measured members, `shared`: only members measured in both layers.

| member set | collection | programmes scored in both layers | median Spearman | minimum | maximum | programmes with Spearman p < 0.05 |
|---|---|---|---|---|---|---|
| own | hallmark | 19 | 0.09 | -0.15 | 0.63 | 4 |
| own | reactome | 113 | 0.14 | -0.15 | 0.73 | 47 |
| shared | hallmark | 16 | 0.34 | 0.05 | 0.58 | 13 |
| shared | reactome | 84 | 0.34 | -0.02 | 0.71 | 62 |

**Column glossary.**
- `cohorts.csv`: `env`, `layer`, `cohort`, `role` (primary or sensitivity), `n_positive`, `n_negative`, `n_features_measured`, `n_scored_<collection>`, `n_programmes_<collection>`.
- `coverage.csv`: `layer`, `env`, `programme`, `collection`, `n_members` (members in the layer's feature space), `n_measured`, `fraction_measured`, `scored`.
- `effects.csv`, `effects_var.csv`: per programme, `collection` and Hedges g (or its variance) per environment; missing where the programme is not scored.
- `partition.csv`: `analysis`, `programme`, `collection`, the meta-analysis columns of `chorus.features.meta_analyse` (`effect`, `se`, `z`, `p`, `Q`, `Q_p`, `I2`, `tau2`, `k`, `sign_consistency`), `set`, `max_abs_g`, `n_env_scored`, `n_env`, `environments`; the meta-analysis columns are empty for `not_partitioned`.
- `loso.csv`: `scheme`, `test_role`, `train`, `test`, `n_train`, `n_test_pos`, `n_test_neg`, `n_programmes`, `harmoniser`, `classifier`, the metrics of loso_results.csv, `auroc_lo`, `auroc_hi`, `skipped` (empty: nothing was skipped).
- `loso_predictions.csv`: `scheme`, `test`, `harmoniser`, `classifier`, `sample` (`<cohort>::<sample_id>`), `y`, `p`.
- `bridge_cervia.csv`: `member_set`, `programme`, `collection`, `n_samples`, `spearman`, `spearman_p`, `pearson`, `pearson_p`, `n_members_somascan`, `n_members_ms`, `n_members_shared`.

**Code and config changes:**
- `chorus.integrate.programmes`, new.
- `preprocessing: features: per_cohort`: every cohort keeps the features that pass its own filter.
- `acute_lc_status`: per cohort, for `relabel_acute`.
- New configs: `configs/mosaic_lc_rna.yaml`, `configs/mosaic_acute_rna.yaml` and `configs/mosaic_acute_protein.yaml`.

**Not committed:** `scores.csv.gz` and the log (`results/mosaic_programmes.log`).

## Metabolomics (16 November 2026)

**Cohorts and identifiers** (`docs/metabolomics_cohorts.md`).
- **Built:** Martinez2024 (Metabolomics Workbench ST003103), Lopez2023 (Mendeley 8zfdjsypd8), Saito2024
  (Front Immunol supplementary Table 2) and Holmes2021 (Zenodo 4713391).
- **Not built:** Mikuteit2023: its supplement holds group summaries only.
- **Labels:** Holmes2021 is built but `lc_status` unknown, because the deposit does not explain its symptom composites.
- **Identifiers:** metabolite names are mapped to a shared identifier through RefMet: `HMDB:<id>`, else
  `RefMet:<name>` (`chorus.metabolites`, adapter `metabolomics_id`).
- **Mapping tables:** `lc_metabolomics/mapping/`:
  - `mapping_summary.csv`: per cohort, names mapped by each rule;
  - `n_shared.csv`: identifiers shared by every pair, Su2022_metabolomics included;
  - `shared_by_k.csv`: identifiers in at least k cohorts, and shared by each triple of labelled cohorts.

**Why these three cohorts for Test 1.** Of the labelled cohorts (Martinez2024, Lopez2023, Saito2024,
Su2022_metabolomics), one triple shares 30 or more identifiers: Martinez2024, Lopez2023 and Su2022_metabolomics (38).
- **Other triples:** Martinez2024 + Lopez2023 + Saito2024 share 15; Martinez2024 + Saito2024 + Su2022_metabolomics
  share 17; Lopez2023 + Saito2024 + Su2022_metabolomics share 25.
- **All four:** share 15.
- **So:** Test 1 runs on that triple and its 38 identifiers (`configs/lc_metabolomics.yaml`). Saito2024 enters the
  concordance only.

**Test 1 setup** (`lc_metabolomics`). Leave one study out, as for proteomics:
- **Harmonisers:** identity, center_scale, rank, fsqn, combat_addon and chorus, plus nested.
- **Selectors:** variance, differential and invariant. Their sizes are the proteomics fractions of the shared features
  (100 and 50 of 376) applied to 38: variance 10, differential and invariant 5.
- **Classifiers:** logreg and rf; nested selection on.
- **Within-study CV:** 5 folds × 10 repeats (per-fold method).
- **Samples:** post-acute, one per subject: Martinez2024 40 LC / 40 COVID_recovered, Lopez2023 30 / 18,
  Su2022_metabolomics 68 / 45 (T3).
- **Imputation:** `study_min`. 29 values were filled in Martinez2024 (one sample has no GC-MS values) and 1 in
  Lopez2023.
- **Su2022_metabolomics** holds the authors' adjusted Z scores, so the harmoniser table of "Su2022 and the
  harmonisers" (Proteomics, 2.) applies to it in the same way. The other two cohorts are log2 areas and log2 µM.

**chorus.** All 9 chorus rows (3 held-out studies × 3 selectors) carry an error: only 12 or 14 null metabolites
passed the screen, against the 20 required. As for the acute proteomics control, the thresholds were not relaxed.
They are error rows in `loso_results.csv` (no metrics, no predictions). Since 16 November `loso_summary.csv` keeps them
as three rows (one per selector, classifier `not_fitted`, `n_studies` 0, `n_errors` 3, messages in `error`); the other
rows of the rewritten summary are unchanged.
The harmonisation diagnostics carry the same message for chorus. `nested` chooses among the other five candidates
when chorus fails in an inner fold. Choices over the 18 nested rows: rank 8, center_scale 4, combat_addon 2, identity 2, fsqn 2.

**LOSO, mean AUROC across the three held-out studies (minimum to maximum)** (`loso_summary.csv`; per held-out
rows with bootstrap intervals in `loso_results.csv`):

| harmoniser | differential logreg | differential rf | invariant logreg | invariant rf | variance logreg | variance rf |
|---|---|---|---|---|---|---|
| identity | 0.594 (0.504 to 0.696) | 0.536 (0.425 to 0.626) | 0.586 (0.492 to 0.7) | 0.519 (0.5 to 0.534) | 0.528 (0.489 to 0.58) | 0.439 (0.369 to 0.547) |
| center_scale | 0.621 (0.577 to 0.667) | 0.596 (0.502 to 0.649) | 0.607 (0.575 to 0.644) | 0.501 (0.428 to 0.551) | 0.569 (0.538 to 0.593) | 0.586 (0.495 to 0.71) |
| rank | 0.554 (0.528 to 0.574) | 0.553 (0.479 to 0.639) | 0.539 (0.476 to 0.582) | 0.486 (0.479 to 0.491) | 0.567 (0.485 to 0.622) | 0.558 (0.444 to 0.71) |
| fsqn | 0.569 (0.498 to 0.622) | 0.544 (0.519 to 0.57) | 0.547 (0.51 to 0.594) | 0.475 (0.41 to 0.514) | 0.563 (0.469 to 0.668) | 0.505 (0.445 to 0.56) |
| combat_addon | 0.619 (0.581 to 0.666) | 0.593 (0.475 to 0.653) | 0.607 (0.575 to 0.643) | 0.501 (0.419 to 0.564) | 0.557 (0.455 to 0.721) | 0.519 (0.473 to 0.6) |
| nested | 0.606 (0.574 to 0.666) | 0.598 (0.502 to 0.652) | 0.574 (0.566 to 0.582) | 0.48 (0.428 to 0.534) | 0.525 (0.469 to 0.622) | 0.533 (0.444 to 0.71) |

**Within-study CV, mean AUROC over folds** (`within_study_cv.csv`; per fold in `within_study_cv_folds.csv`). With
k = 5, differential and invariant select the same metabolites within one study, so their rows are equal.

| study | differential logreg | differential rf | invariant logreg | invariant rf | variance logreg | variance rf |
|---|---|---|---|---|---|---|
| Lopez2023 | 0.661 | 0.559 | 0.661 | 0.559 | 0.621 | 0.59 |
| Martinez2024 | 0.804 | 0.765 | 0.804 | 0.765 | 0.779 | 0.821 |
| Su2022_metabolomics | 0.511 | 0.508 | 0.511 | 0.508 | 0.519 | 0.513 |

**Harmonisation diagnostics** (`lc_metabolomics/harmonisation`). These are the transcriptomics metrics, run with
`--n-pcs 10` because there are 38 features. `all_metrics.csv` has 213 rows. Chorus has one row per held-out study, with family `error` and the null screen message in `scope`.

**Concordance** (`lc_metabolomics/concordance`: the three Test 1 cohorts on their 38 identifiers, with `top20` files;
`lc_metabolomics/concordance_all`: the four labelled cohorts each on its own identifiers, pairwise over those shared,
`n_shared.csv`, with `top50` files). Spearman correlation of Hedges g (LC minus COVID_recovered):

Test 1 cohorts, 38 identifiers:

| | Martinez2024 | Lopez2023 | Su2022_metabolomics |
|---|---|---|---|
| Martinez2024 |  | 0.112 (38) | 0.532 (38) |
| Lopez2023 | 0.112 (38) |  | 0.345 (38) |
| Su2022_metabolomics | 0.532 (38) | 0.345 (38) |  |

All labelled cohorts, pairwise over shared identifiers (number shared in parentheses):

| | Martinez2024 | Lopez2023 | Saito2024 | Su2022_metabolomics |
|---|---|---|---|---|
| Martinez2024 |  | 0.156 (41) | -0.287 (17) | 0.058 (158) |
| Lopez2023 | 0.156 (41) |  | 0.007 (28) | 0.213 (61) |
| Saito2024 | -0.287 (17) | 0.007 (28) |  | 0.104 (66) |
| Su2022_metabolomics | 0.058 (158) | 0.213 (61) | 0.104 (66) |  |

- **`top50` files of `concordance_all`:** the 50 identifiers with the largest mean |g| are taken over the union of
  the four cohorts. Each pair of cohorts has 1 to 4 of them in common (`n_shared_top50.csv`). Before the reporting
  rerun of 16 November the file gave correlations of plus or minus 1 on 2 to 4 values; it now gives NaN for every
  off-diagonal pair (fewer than 10 shared).
- In `concordance_all/pairs.csv` no pair over all identifiers has fewer than 10 shared (the smallest is 17,
  Martinez2024 and Saito2024).

**Epigenomics** (`concordance_epigenomics/cohorts.csv`, `docs/metabolomics_cohorts.md`). GSE174818, GSE197152 and
GSE247869 were downloaded (EPIC, whole blood); the groups are documented, but no per-CpG effect was computed.
- **GSE174818:** acute samples only.
- **GSE197152:** 15 one-year samples, all described with PASC symptoms, so there is no recovered group.
- **GSE247869:** 96 six-month samples. The paper's PASC subgroup (28) is not identifiable per sample in GEO or in
  the supplement.
- **So:** no series has per-sample LC and COVID_recovered labels.

**Code and config changes:**
- `chorus.metabolites` (RefMet matching with a cache, identifier rules, pairwise shared counts);
- the `metabolomics_id` adapter;
- `scripts/metabolomics/build_cohorts.py` and `scripts/metabolomics/map_ids.py`;
- `scripts/concordance.py`, which labels `metabolite_id` features with the RefMet name;
- `configs/lc_metabolomics.yaml` and `configs/concordance_metabolomics_all.yaml`;
- xlrd 2.0.2 added to the venv.

**Not committed:** the logs and `selected_features.json`.

## Version 4 cohort set and reporting fixes (16 November 2026)

The results of earlier versions are kept in their folders; the version 4 runs write to new folders with `_v4`.

**1. Cohorts added** (definitions: `docs/proteomics_cohorts.md` and `docs/count_matrices.md`, sections "Cohorts added
16 November 2026"; provenance and checksums: `docs/provenance.md`).

| study | layer | source | samples used (LC / COVID_recovered) | label |
|---|---|---|---|---|
| Fineschi2026 | plasma, Olink Explore 384 Inflammation | SciLifeLab 10.17044/scilifelab.31316779 | 165 (92 / 73) | sample id prefix PA (LC) or KO (recovered) |
| Chowdhury2026 | plasma, DIA-MS (DIA-NN 1.8.1) | PRIDE PXD066724 | 115 (53 / 62); 35 never_infected not in the contrast | SDRF `characteristics[disease]` |
| PRJNA1149196_SA | whole blood total RNA | ENA PRJNA1149196, cohort 1 | see section 3 | BioSample `disease` |
| PRJNA1149196_QLD | whole blood total RNA | ENA PRJNA1149196, cohort 2 | see section 3 (training only) | BioSample `disease` |
| GSE297499 | PBMC, author TPM | GEO GSE297499, PBMC only | see section 3 | GEO `treatment` |

- **Not added:** NHLBI 10.25444/nhlbi.26383798 (Huapaya2024, SomaScan): the workbook has no per-sample PASC or
  symptom column.
- **Participant overlap, SciLifeLab:** all 110 labelled samples of the PBMC RNA-seq set 10.17044/scilifelab.28832492
  (`SamplesPC.txt`) are Fineschi2026 Olink participants (60 PA, 50 KO); `overlaps_with_study` marks them. The PBMC
  RNA-seq set is not a cohort of the transcriptomic benchmark.
- **Other overlaps recorded, not resolvable from the deposits:** Chowdhury2026 and Chowdhury2024 (same laboratory and
  BQC19 Sherbrooke biobank); PRJNA1149196_SA and GSE169687 (same South Australian participants; GSE169687 is disabled).
- **Chowdhury2026 label check:** five second-batch runs (samples 226, 227, 231, 232, 233) have file code `Lc` and SDRF
  disease C; the SDRF is used (`file_group_code` and the notes keep the code).

**2. Test 1 in proteomics, version 4** (`lc_proteomics_affinity_v4`, `lc_proteomics_ms_v4`).
- **Configs:** `configs/lc_proteomics_affinity_v4.yaml` (Fineschi2026 added) and `configs/lc_proteomics_ms_v4.yaml`
  (Chowdhury2026 added). Grid, preprocessing, imputation and selector sizes as in "Proteomics (15 November 2026)".
- **Failed fits:** none, in 168 rows per block.
- **chorus factors** (`n_factors`, the same for every selector and classifier of a fold): affinity Buggert2025 1,
  Cervia2024_SomaScan 5, Fineschi2026 1, Su2022 1; MS Cervia2024_MS 1, Chowdhury2024 2, Chowdhury2026 1,
  Kovarik2023_plasma 1.
- **Nested choices** (24 rows per block): affinity identity 10, chorus 6, rank 2, combat_addon 2, center_scale 2, fsqn 2;
  MS rank 8, identity 6, chorus 4, combat_addon 2, center_scale 2, fsqn 2.

| block | cohort | samples (LC / COVID_recovered) | proteins used | values imputed |
|---|---|---|---|---|
| affinity | Buggert2025 | 197 (146 / 51) | 166 | 88 (median) |
| affinity | Cervia2024_SomaScan | 113 (40 / 73) | 166 | 0 |
| affinity | Su2022 | 125 (75 / 50) | 166 | 0 |
| affinity | Fineschi2026 | 165 (92 / 73) | 166 | 0 |
| MS | Cervia2024_MS | 113 (40 / 73) | 179 | 842 (minimum) |
| MS | Chowdhury2024 | 30 (21 / 9) | 179 | 56 (minimum) |
| MS | Kovarik2023_plasma | 26 (13 / 13) | 179 | 593 (minimum) |
| MS | Chowdhury2026 | 115 (53 / 62) | 179 | 198 (minimum) |

- The proteins shared by all cohorts fall from 376 to 166 (affinity; Fineschi2026 measures 358) and from 182 to 179 (MS).
- In the MS block 47 (Cervia2024_MS), 11 (Chowdhury2024), 66 (Kovarik2023_plasma) and 21 (Chowdhury2026) of the 179
  proteins have at least one missing value; the lowest observed fractions are 0.44, 0.20, 0.08 and 0.23. In the affinity
  block one protein of Buggert2025 has missing values (lowest observed fraction 0.55).

**Mean AUROC across the four held-out studies (minimum to maximum)** (`loso_summary.csv`; per held-out rows with
bootstrap intervals in `loso_results.csv`).

Affinity block (`lc_proteomics_affinity_v4`):

| harmoniser | differential logreg | differential rf | invariant logreg | invariant rf | variance logreg | variance rf |
|---|---|---|---|---|---|---|
| identity | 0.504 (0.416 to 0.577) | 0.486 (0.391 to 0.558) | 0.555 (0.504 to 0.593) | 0.475 (0.267 to 0.574) | 0.523 (0.503 to 0.538) | 0.441 (0.309 to 0.5) |
| center_scale | 0.467 (0.444 to 0.488) | 0.449 (0.389 to 0.518) | 0.524 (0.466 to 0.584) | 0.461 (0.39 to 0.531) | 0.488 (0.44 to 0.578) | 0.471 (0.45 to 0.493) |
| rank | 0.540 (0.48 to 0.63) | 0.606 (0.511 to 0.709) | 0.558 (0.527 to 0.593) | 0.582 (0.497 to 0.702) | 0.484 (0.427 to 0.602) | 0.559 (0.474 to 0.696) |
| fsqn | 0.482 (0.413 to 0.559) | 0.494 (0.401 to 0.591) | 0.452 (0.374 to 0.556) | 0.472 (0.408 to 0.577) | 0.460 (0.437 to 0.512) | 0.516 (0.444 to 0.569) |
| combat_addon | 0.470 (0.448 to 0.493) | 0.456 (0.386 to 0.514) | 0.524 (0.471 to 0.584) | 0.467 (0.408 to 0.524) | 0.515 (0.465 to 0.583) | 0.506 (0.451 to 0.564) |
| chorus | 0.509 (0.438 to 0.603) | 0.548 (0.488 to 0.626) | 0.548 (0.472 to 0.629) | 0.551 (0.478 to 0.64) | 0.513 (0.47 to 0.577) | 0.550 (0.507 to 0.615) |
| nested | 0.494 (0.438 to 0.577) | 0.478 (0.386 to 0.561) | 0.536 (0.472 to 0.593) | 0.483 (0.267 to 0.574) | 0.477 (0.437 to 0.538) | 0.555 (0.459 to 0.696) |

MS block (`lc_proteomics_ms_v4`):

| harmoniser | differential logreg | differential rf | invariant logreg | invariant rf | variance logreg | variance rf |
|---|---|---|---|---|---|---|
| identity | 0.498 (0.272 to 0.683) | 0.507 (0.257 to 0.765) | 0.472 (0.213 to 0.624) | 0.561 (0.499 to 0.659) | 0.537 (0.479 to 0.56) | 0.514 (0.423 to 0.635) |
| center_scale | 0.494 (0.385 to 0.635) | 0.423 (0.246 to 0.548) | 0.509 (0.396 to 0.619) | 0.396 (0.192 to 0.576) | 0.495 (0.365 to 0.582) | 0.454 (0.34 to 0.574) |
| rank | 0.546 (0.433 to 0.746) | 0.538 (0.459 to 0.624) | 0.566 (0.385 to 0.751) | 0.546 (0.467 to 0.728) | 0.490 (0.361 to 0.539) | 0.401 (0.272 to 0.487) |
| fsqn | 0.416 (0.278 to 0.519) | 0.452 (0.287 to 0.586) | 0.414 (0.112 to 0.565) | 0.369 (0.0769 to 0.595) | 0.406 (0.302 to 0.507) | 0.463 (0.302 to 0.566) |
| combat_addon | 0.485 (0.391 to 0.608) | 0.402 (0.186 to 0.566) | 0.511 (0.396 to 0.624) | 0.394 (0.172 to 0.561) | 0.482 (0.397 to 0.592) | 0.459 (0.257 to 0.571) |
| chorus | 0.480 (0.401 to 0.561) | 0.412 (0.272 to 0.517) | 0.504 (0.379 to 0.64) | 0.457 (0.337 to 0.524) | 0.403 (0.26 to 0.524) | 0.420 (0.225 to 0.529) |
| nested | 0.468 (0.278 to 0.608) | 0.447 (0.287 to 0.597) | 0.532 (0.379 to 0.619) | 0.439 (0.337 to 0.495) | 0.513 (0.479 to 0.55) | 0.464 (0.423 to 0.545) |

**Diagnostics and concordance.**
- `lc_proteomics_affinity_v4/harmonisation`, `lc_proteomics_ms_v4/harmonisation`: as for the three-cohort blocks
  (n_pcs 50, no sex rows, marker rows empty); 432 rows each, no failed fit.
- `lc_proteomics_affinity_v4/concordance`, `lc_proteomics_ms_v4/concordance`: on the 166 and 179 shared proteins,
  `--top 50`.
- `concordance_proteomics_all_v4`: the eight cohorts of both blocks, each on its own proteins (`features: union`), from
  `configs/concordance_proteomics_all_v4.yaml`. Proteins shared by the new cohorts (`n_shared.csv`): Fineschi2026 with
  Buggert2025 358, Cervia2024_SomaScan 307, Su2022 186, and with the four MS cohorts 2 to 5; Chowdhury2026 with
  Chowdhury2024 334, Kovarik2023_plasma 260, Cervia2024_MS 224, Cervia2024_SomaScan 286, Buggert2025 197, Su2022 10,
  Fineschi2026 5. Chowdhury2026 has 472 proteins here: one of its 473 has no value in the 115 retained samples.

**3. Test 1 in transcriptomics, version 4** (`lc_transcriptomics_blood_v4`).
- **Config:** `configs/lc_transcriptomics_blood_v4.yaml`: the five cohorts of `lc_transcriptomics_blood` plus
  PRJNA1149196_SA, PRJNA1149196_QLD and GSE297499. Grid, preprocessing and selectors as in `lc_transcriptomics_blood`
  (9 harmonisers plus nested selection, 6 selector configurations, classifiers logreg, rf, taco and cross_study).
- **PRJNA1149196_QLD** is quantified on all reads; SA on the first 5,000,000 read pairs of each run (`docs/provenance.md`,
  "PRJNA1149196"). A first version 4 run made with a QLD read prefix (16 November) was not kept: the QLD FASTQ files are
  ordered by alignment position, so the prefix held only chr1 and chr10 reads and QLD passed the expression filter
  with 5,914 genes, which cut the genes shared by all cohorts to 2,574.
- **Samples used** (one sample per subject, earliest post-acute; the `never_infected` samples are not in the contrast):

| cohort | matrix | samples (LC / COVID_recovered) | never_infected (not used) | role |
|---|---|---|---|---|
| GSE226260_WB | whole blood | 37 (17 / 20) | 0 | held out and training |
| GSE226260_PBMC | PBMC | 36 (28 / 8) | 0 | held out and training |
| GSE224615 | PBMC | 36 (23 / 13) | 0 | held out and training |
| GSE251849 | PBMC | 16 (11 / 5) | 7 | held out and training |
| GSE270045 | whole blood | 36 (19 / 17) | 0 | held out and training |
| PRJNA1149196_SA | whole blood | 39 (28 / 11) | 15 | held out and training |
| PRJNA1149196_QLD | whole blood | 9 (5 / 4) | 0 | training only (fewer than 5 recovered; `min_test_per_class` 5) |
| GSE297499 | PBMC | 26 (17 / 9) | 0 | held out and training |

- **Genes:** 10,053 pass the expression filter in every cohort (CPM at least 1 in at least half of the samples; 10,442
  in `lc_transcriptomics_blood`). PRJNA1149196_SA keeps its 44-week sample per subject (38 subjects have two visits).
  GSE297499 enters as author TPM (`author_normalised`, log2(x + 1)).
- **Training pools:** 7 training studies per fold, 196 to 219 samples.
- **Failed fits** (`loso_results.csv`, column `error`): `chorus` and `chorus_protected` in two folds, every selector:
  "only 18 null genes passed the screen" (GSE224615 held out) and "only 19" (PRJNA1149196_SA held out). Their summary
  rows have `n_studies` 5. Nested selection chose `chorus` in two fits of the GSE224615 fold (selectors differential and
  TACO 411 only), which failed with the same message, so those nested rows have `n_studies` 6. `loso_summary.csv` keeps
  one `not_fitted` row per failed configuration (14). The default null screen was not changed.
- **chorus factors** (`n_factors`, the folds that fitted): `chorus` 4 to 6, `chorus_protected` 3 to 5.
- **Nested choices** (162 rows): fsqn 52, quantile_ref 32, identity 28, center_scale 16, combat_pooled 16, chorus 10
  (2 of them the failed fits), combat_addon 8.

**Mean AUROC across the held-out studies (minimum to maximum)** (`loso_summary.csv`; per held-out rows with bootstrap
intervals in `loso_results.csv`). Seven held-out studies unless a count in brackets says otherwise; "not fitted": no
fold fitted. Selectors: variance (500 most variable genes), differential (200), invariant (200), "TACO 411 only"
(`priors/taco_411_dce_genes.txt`, the genes in the matrix, `prior_only`), "TACO 411 + 200" (those genes together with
the 200 invariant genes, `union`), "MR 32 + 200" (`priors/mrcontrol_32_genes.txt` together with the 200 invariant
genes). `n_features` in `loso_results.csv`: 502, 202, 202, 318, 503 to 512 and 225 to 227.

`logreg`:

| harmoniser | variance | differential | invariant | TACO 411 only | TACO 411 + 200 | MR 32 + 200 |
|---|---|---|---|---|---|---|
| identity | 0.413 (0.2 to 0.616) | 0.490 (0.266 to 0.689) | 0.376 (0.145 to 0.594) | 0.464 (0.155 to 0.806) | 0.410 (0.226 to 0.791) | 0.364 (0.109 to 0.579) |
| center_scale | 0.493 (0.309 to 0.662) | 0.439 (0.218 to 0.606) | 0.402 (0.182 to 0.532) | 0.482 (0.0464 to 0.721) | 0.427 (0.142 to 0.676) | 0.390 (0.218 to 0.504) |
| quantile_ref | 0.465 (0.273 to 0.675) | 0.441 (0.327 to 0.542) | 0.398 (0.102 to 0.602) | 0.470 (0.152 to 0.826) | 0.412 (0.124 to 0.609) | 0.396 (0.167 to 0.582) |
| fsqn | 0.466 (0.4 to 0.606) | 0.400 (0.0182 to 0.579) | 0.496 (0.309 to 0.632) | 0.504 (0.291 to 0.68) | 0.462 (0.236 to 0.644) | 0.468 (0.255 to 0.603) |
| rank | 0.435 (0.218 to 0.601) | 0.499 (0.396 to 0.605) | 0.428 (0.0341 to 0.614) | 0.408 (0.105 to 0.668) | 0.411 (0.0279 to 0.608) | 0.412 (0.0217 to 0.575) |
| combat_addon | 0.455 (0.272 to 0.682) | 0.441 (0.218 to 0.632) | 0.398 (0.182 to 0.526) | 0.474 (0.0495 to 0.747) | 0.422 (0.13 to 0.679) | 0.381 (0.218 to 0.498) |
| combat_pooled | 0.502 (0.341 to 0.716) | 0.456 (0.218 to 0.608) | 0.385 (0.109 to 0.547) | 0.486 (0.096 to 0.726) | 0.421 (0.183 to 0.694) | 0.384 (0.164 to 0.518) |
| chorus | 0.488 (0.182 to 0.683) [5] | 0.293 (0.121 to 0.4) [5] | 0.385 (0.0545 to 0.559) [5] | 0.492 (0.0867 to 0.712) [5] | 0.421 (0.176 to 0.671) [5] | 0.361 (0.0545 to 0.518) [5] |
| chorus_protected | 0.395 (0.105 to 0.682) [5] | 0.355 (0.238 to 0.453) [5] | 0.341 (0.109 to 0.571) [5] | 0.471 (0.0681 to 0.719) [5] | 0.410 (0.155 to 0.662) [5] | 0.338 (0.0909 to 0.547) [5] |
| nested | 0.426 (0.272 to 0.565) | 0.494 (0.345 to 0.621) [6] | 0.437 (0.26 to 0.594) | 0.444 (0.0867 to 0.68) [6] | 0.406 (0.13 to 0.609) | 0.424 (0.255 to 0.594) |

`rf`:

| harmoniser | variance | differential | invariant | TACO 411 only | TACO 411 + 200 | MR 32 + 200 |
|---|---|---|---|---|---|---|
| identity | 0.459 (0.325 to 0.618) | 0.456 (0.299 to 0.725) | 0.466 (0.29 to 0.686) | 0.415 (0.276 to 0.582) | 0.421 (0.282 to 0.559) | 0.471 (0.257 to 0.626) |
| center_scale | 0.465 (0.276 to 0.737) | 0.389 (0.0944 to 0.631) | 0.424 (0.17 to 0.691) | 0.458 (0.135 to 0.813) | 0.417 (0.138 to 0.716) | 0.394 (0.133 to 0.612) |
| quantile_ref | 0.483 (0.364 to 0.606) | 0.549 (0.315 to 0.769) | 0.459 (0.275 to 0.715) | 0.498 (0.339 to 0.694) | 0.433 (0.209 to 0.616) | 0.445 (0.297 to 0.597) |
| fsqn | 0.495 (0.325 to 0.694) | 0.364 (0.0759 to 0.626) | 0.372 (0.138 to 0.572) | 0.468 (0.241 to 0.831) | 0.377 (0.147 to 0.671) | 0.370 (0.0619 to 0.599) |
| rank | 0.523 (0.393 to 0.65) | 0.496 (0.186 to 0.679) | 0.448 (0.234 to 0.699) | 0.492 (0.4 to 0.614) | 0.398 (0.113 to 0.65) | 0.442 (0.189 to 0.686) |
| combat_addon | 0.415 (0.228 to 0.624) | 0.390 (0.0464 to 0.629) | 0.424 (0.173 to 0.679) | 0.474 (0.229 to 0.809) | 0.413 (0.0681 to 0.719) | 0.404 (0.135 to 0.632) |
| combat_pooled | 0.482 (0.246 to 0.766) | 0.396 (0.0759 to 0.626) | 0.410 (0.128 to 0.66) | 0.447 (0.166 to 0.74) | 0.417 (0.105 to 0.725) | 0.411 (0.138 to 0.644) |
| chorus | 0.392 (0.113 to 0.604) [5] | 0.351 (0.0743 to 0.542) [5] | 0.398 (0.195 to 0.624) [5] | 0.454 (0.226 to 0.662) [5] | 0.386 (0.118 to 0.65) [5] | 0.398 (0.153 to 0.596) [5] |
| chorus_protected | 0.379 (0.031 to 0.712) [5] | 0.391 (0.0588 to 0.641) [5] | 0.391 (0.0743 to 0.585) [5] | 0.399 (0.0325 to 0.8) [5] | 0.396 (0.0464 to 0.635) [5] | 0.391 (0.0418 to 0.601) [5] |
| nested | 0.427 (0.228 to 0.737) | 0.383 (0.0759 to 0.725) [6] | 0.416 (0.17 to 0.715) | 0.458 (0.226 to 0.831) [6] | 0.390 (0.0681 to 0.616) | 0.393 (0.133 to 0.597) |

`taco` (soft-voting ensemble of rf, extratrees, svc and mlp; TabPFN disabled, `CHORUS_DISABLE_TABPFN=1`):

| harmoniser | variance | differential | invariant | TACO 411 only | TACO 411 + 200 | MR 32 + 200 |
|---|---|---|---|---|---|---|
| identity | 0.458 (0.218 to 0.678) | 0.497 (0.263 to 0.709) | 0.399 (0.109 to 0.612) | 0.448 (0.3 to 0.588) | 0.462 (0.273 to 0.568) | 0.416 (0.109 to 0.559) |
| center_scale | 0.465 (0.299 to 0.659) | 0.389 (0.161 to 0.609) | 0.378 (0.155 to 0.562) | 0.458 (0.0929 to 0.724) | 0.390 (0.0681 to 0.671) | 0.367 (0.111 to 0.591) |
| quantile_ref | 0.507 (0.327 to 0.718) | 0.492 (0.304 to 0.665) | 0.429 (0.198 to 0.609) | 0.528 (0.364 to 0.709) | 0.428 (0.158 to 0.621) | 0.403 (0.266 to 0.565) |
| fsqn | 0.447 (0.315 to 0.665) | 0.347 (0.0364 to 0.594) | 0.363 (0.18 to 0.568) | 0.455 (0.254 to 0.682) | 0.387 (0.173 to 0.624) | 0.376 (0.111 to 0.556) |
| rank | 0.458 (0.36 to 0.659) | 0.490 (0.309 to 0.686) | 0.421 (0.158 to 0.662) | 0.465 (0.18 to 0.712) | 0.391 (0.096 to 0.656) | 0.419 (0.0743 to 0.65) |
| combat_addon | 0.415 (0.142 to 0.676) | 0.388 (0.136 to 0.605) | 0.369 (0.133 to 0.562) | 0.445 (0.0836 to 0.729) | 0.385 (0.0433 to 0.676) | 0.361 (0.108 to 0.571) |
| combat_pooled | 0.444 (0.279 to 0.671) | 0.400 (0.152 to 0.619) | 0.378 (0.146 to 0.582) | 0.451 (0.136 to 0.715) | 0.386 (0.0557 to 0.676) | 0.367 (0.0898 to 0.603) |
| chorus | 0.405 (0.183 to 0.662) [5] | 0.334 (0.0279 to 0.529) [5] | 0.354 (0.183 to 0.606) [5] | 0.470 (0.142 to 0.741) [5] | 0.394 (0.115 to 0.644) [5] | 0.354 (0.155 to 0.585) [5] |
| chorus_protected | 0.389 (0.0495 to 0.744) [5] | 0.366 (0.0867 to 0.585) [5] | 0.359 (0.102 to 0.624) [5] | 0.464 (0.0712 to 0.744) [5] | 0.408 (0.0774 to 0.691) [5] | 0.342 (0.0681 to 0.603) [5] |
| nested | 0.414 (0.142 to 0.659) | 0.436 (0.282 to 0.709) [6] | 0.393 (0.155 to 0.609) | 0.428 (0.142 to 0.682) [6] | 0.424 (0.0433 to 0.621) | 0.382 (0.111 to 0.536) |

`cross_study` (one logreg per training study, weighted by its mean AUROC on the other training studies minus 0.5):

| harmoniser | variance | differential | invariant | TACO 411 only | TACO 411 + 200 | MR 32 + 200 |
|---|---|---|---|---|---|---|
| identity | 0.485 (0.255 to 0.694) | 0.561 (0.196 to 0.908) | 0.487 (0.29 to 0.609) | 0.518 (0 to 0.727) | 0.434 (0.294 to 0.7) | 0.491 (0.272 to 0.671) |
| center_scale | 0.420 (0.182 to 0.602) | 0.458 (0.238 to 0.627) | 0.371 (0.0867 to 0.571) | 0.586 (0.435 to 0.691) | 0.413 (0.0991 to 0.653) | 0.363 (0.0743 to 0.571) |
| quantile_ref | 0.479 (0.345 to 0.592) | 0.406 (0.228 to 0.603) | 0.431 (0.156 to 0.612) | 0.449 (0.0217 to 0.694) | 0.415 (0.134 to 0.612) | 0.438 (0.143 to 0.603) |
| fsqn | 0.501 (0.412 to 0.672) | 0.326 (0.182 to 0.592) | 0.376 (0.121 to 0.659) | 0.504 (0.433 to 0.575) | 0.394 (0.139 to 0.644) | 0.364 (0.096 to 0.618) |
| rank | 0.544 (0.375 to 0.732) | 0.401 (0 to 0.778) | 0.439 (0.0526 to 0.697) | 0.448 (0.277 to 0.662) | 0.406 (0.0526 to 0.696) | 0.448 (0.0712 to 0.691) |
| combat_addon | 0.483 (0.255 to 0.683) | 0.458 (0.211 to 0.641) | 0.368 (0.0743 to 0.571) | 0.590 (0.438 to 0.697) | 0.413 (0.102 to 0.653) | 0.363 (0.0774 to 0.576) |
| combat_pooled | 0.518 (0.291 to 0.696) | 0.451 (0.214 to 0.634) | 0.370 (0.0836 to 0.568) | 0.579 (0.451 to 0.691) | 0.405 (0.0929 to 0.65) | 0.370 (0.0743 to 0.576) |
| chorus | 0.406 (0.223 to 0.529) [5] | 0.328 (0.00929 to 0.523) [5] | 0.350 (0.096 to 0.594) [5] | 0.331 (0.0557 to 0.482) [5] | 0.338 (0.0557 to 0.568) [5] | 0.340 (0.096 to 0.574) [5] |
| chorus_protected | 0.431 (0.255 to 0.624) [5] | 0.357 (0.0155 to 0.51) [5] | 0.333 (0.0186 to 0.579) [5] | 0.356 (0.139 to 0.497) [5] | 0.363 (0.0433 to 0.594) [5] | 0.332 (0.0433 to 0.579) [5] |
| nested | 0.442 (0.182 to 0.672) | 0.493 (0.195 to 0.908) [6] | 0.375 (0.0867 to 0.612) | 0.385 (0.0557 to 0.569) [6] | 0.384 (0.102 to 0.591) | 0.383 (0.0743 to 0.603) |

**Ten highest means** (configurations with at least one fitted fold). The `taco` row ran without TabPFN
(`CHORUS_DISABLE_TABPFN=1`).

| harmoniser | selector | classifier | mean | median | min | max | n_studies |
|---|---|---|---|---|---|---|---|
| combat_addon | TACO 411 only | cross_study | 0.590 | 0.582 | 0.438 | 0.697 | 7 |
| center_scale | TACO 411 only | cross_study | 0.586 | 0.601 | 0.435 | 0.691 | 7 |
| combat_pooled | TACO 411 only | cross_study | 0.579 | 0.557 | 0.451 | 0.691 | 7 |
| identity | differential | cross_study | 0.561 | 0.522 | 0.196 | 0.908 | 7 |
| quantile_ref | differential | rf | 0.549 | 0.562 | 0.315 | 0.769 | 7 |
| rank | variance | cross_study | 0.544 | 0.536 | 0.375 | 0.732 | 7 |
| quantile_ref | TACO 411 only | taco | 0.528 | 0.512 | 0.364 | 0.709 | 7 | (no TabPFN)
| rank | variance | rf | 0.523 | 0.527 | 0.393 | 0.650 | 7 |
| identity | TACO 411 only | cross_study | 0.518 | 0.601 | 0.000 | 0.727 | 7 |
| combat_pooled | variance | cross_study | 0.518 | 0.496 | 0.291 | 0.696 | 7 |

**Within-study cross-validation** (`within_study_cv.csv`, `within_study_cv_folds.csv`; stratified 5-fold CV repeated 10
times, mean over the 50 folds, no failed fold; no harmoniser, the classifiers other than `cross_study`). AUROC
with the variance selector, and the range over the six selectors:

| study | n (prevalence LC) | logreg | rf | taco (no TabPFN) | range over selectors, all three classifiers |
|---|---|---|---|---|---|
| GSE226260_WB | 37 (0.46) | 0.722 | 0.721 | 0.722 | 0.693 to 0.814 |
| GSE226260_PBMC | 36 (0.78) | 0.886 | 0.854 | 0.870 | 0.613 to 0.908 |
| GSE224615 | 36 (0.64) | 0.663 | 0.521 | 0.572 | 0.271 to 0.663 |
| GSE251849 | 16 (0.69) | 0.783 | 0.810 | 0.743 | 0.710 to 0.893 |
| GSE270045 | 36 (0.53) | 1.000 | 0.997 | 1.000 | 0.997 to 1.000 |
| PRJNA1149196_SA | 39 (0.72) | 0.416 | 0.476 | 0.428 | 0.247 to 0.476 |
| GSE297499 | 26 (0.65) | 0.805 | 0.833 | 0.829 | 0.802 to 0.847 |

**Diagnostics** (`lc_transcriptomics_blood_harmonisation_v4`): as `lc_transcriptomics_blood_harmonisation` (n_pcs 50,
k 20), on the 10,053 genes; 3,131 rows. RPS4Y1 is added in every fold and XIST where every fitted sample has a value
(`extra_genes`). Four `error` rows: `chorus` and `chorus_protected` with GSE224615 and with PRJNA1149196_SA held out,
with the messages above.

**Concordance** (`concordance_v4`, `scripts/concordance.py configs/lc_transcriptomics_blood_v4.yaml`, `--top 500`): the
eight cohorts on the 10,053 genes (PRJNA1149196_QLD included; 5 / 4 samples). Every pair shares the 10,053 genes (`n_shared.csv`) and
500 of the top 500 (`n_shared_top500.csv`); no pair is flagged in `pairs.csv`. `cohort_annotation.csv`: `control_type` is `COVID_recovered` for PRJNA1149196_SA, PRJNA1149196_QLD and
GSE297499 and `missing` for the five earlier cohorts (metadata built before the column existed).

**Sensitivity run without the three GSE224615 samples flagged at QC** (`lc_transcriptomics_blood_v4_noQCflag`,
`configs/lc_transcriptomics_blood_v4_noQCflag.yaml`: `lc_transcriptomics_blood_v4` with `exclude_samples:
[GSM7027478, GSM7027490, GSM7027501]` for GSE224615, applied before any preprocessing).
- GSE224615 has 33 samples (21 LC / 12 COVID_recovered); the 10,053 filtered genes are the same as in the main run.
- **Failed fits:** `chorus` and `chorus_protected` only with GSE224615 held out ("only 18 null genes passed the
  screen"), and nested selection in two fits of that fold (selectors differential and TACO 411 only), as in the main
  run. With PRJNA1149196_SA held out both chorus variants fitted here. 1,638 rows, 14 with an error (main run 1,602 and
  26).
- **Change from the main run** (fits present in both, 1,576): AUROC difference noQCflag minus main, median per held-out
  study from -0.020 (GSE297499) to 0.004 (GSE226260_WB); over all fits from -0.356 to 0.395. Over the 240 summary
  configurations fitted in both runs, the mean AUROC changes by -0.134 to 0.096 (median -0.004).
- **Ten highest means** (`loso_summary.csv`; no `taco` row among them):

| harmoniser | selector | classifier | mean | median | min | max | n_studies |
|---|---|---|---|---|---|---|---|
| center_scale | TACO 411 only | cross_study | 0.585 | 0.569 | 0.448 | 0.706 | 7 |
| combat_addon | TACO 411 only | cross_study | 0.584 | 0.564 | 0.448 | 0.702 | 7 |
| combat_pooled | TACO 411 only | cross_study | 0.577 | 0.564 | 0.455 | 0.694 | 7 |
| rank | variance | cross_study | 0.562 | 0.539 | 0.375 | 0.690 | 7 |
| quantile_ref | variance | rf | 0.541 | 0.509 | 0.480 | 0.732 | 7 |
| rank | TACO 411 only | rf | 0.525 | 0.534 | 0.257 | 0.759 | 7 |
| quantile_ref | TACO 411 only | rf | 0.518 | 0.492 | 0.306 | 0.771 | 7 |
| chorus | TACO 411 only | logreg | 0.517 | 0.612 | 0.133 | 0.706 | 6 |
| rank | variance | rf | 0.516 | 0.497 | 0.400 | 0.752 | 7 |
| quantile_ref | differential | rf | 0.515 | 0.464 | 0.321 | 0.710 | 7 |

- **Within-study CV of GSE224615** (33 samples): variance selector logreg 0.589, rf 0.539, taco (no TabPFN) 0.515;
  0.317 to 0.589 over selectors and classifiers (main run 0.271 to 0.663). The other studies' within-study values are
  identical to the main run.

**4. Test 3, mosaic programmes, version 4** (`mosaic_programmes_v4/lc_v4`; `scripts/mosaic_programmes.py --out
results/mosaic_programmes_v4 --contrasts lc_v4`, configs `configs/mosaic_lc_rna_v4.yaml` and
`configs/concordance_proteomics_all_v4.yaml`). Method, programmes, harmonisers and classifier as in "Test 3"; LC versus
COVID_recovered on the eight RNA cohorts of `lc_transcriptomics_blood_v4` (PRJNA1149196_QLD on all reads) and the
eight proteomic cohorts of the version 4 blocks (Cervia2024_MS as the sensitivity environment).

**Cohorts and coverage** (`cohorts.csv`, `coverage.csv`). The fraction of members measured is over the members in the
layer's feature space, which is the union of the layer's cohorts: GSE297499 (36,499 genes) adds members to the RNA
space, so the fractions of the earlier RNA cohorts are lower than in `mosaic_programmes/lc` while their measured
members and scored programmes are the same.

| environment | role | positive / negative | features measured | Hallmark scored (of 50) | Reactome scored (of 1,867) | median fraction of members measured (scored programmes) |
|---|---|---|---|---|---|---|
| rna:GSE226260_WB | primary | 17 / 20 | 12801 | 50 | 1617 | 0.81 |
| rna:GSE226260_PBMC | primary | 28 / 8 | 13866 | 50 | 1648 | 0.83 |
| rna:GSE224615 | primary | 23 / 13 | 15734 | 50 | 1660 | 0.85 |
| rna:GSE251849 | primary | 11 / 5 | 13872 | 50 | 1639 | 0.84 |
| rna:GSE270045 | primary | 19 / 17 | 15561 | 50 | 1739 | 0.86 |
| rna:PRJNA1149196_SA | primary | 28 / 11 | 17751 | 50 | 1667 | 0.87 |
| rna:PRJNA1149196_QLD | primary | 5 / 4 | 16055 | 50 | 1656 | 0.86 |
| rna:GSE297499 | primary | 17 / 9 | 36499 | 50 | 1799 | 0.96 |
| protein:Buggert2025 | primary | 146 / 51 | 2903 | 50 | 967 | 0.38 |
| protein:Cervia2024_SomaScan | primary | 40 / 73 | 6402 | 50 | 1540 | 0.90 |
| protein:Su2022 | primary | 75 / 50 | 441 | 33 | 177 | 0.09 |
| protein:Chowdhury2024 | primary | 21 / 9 | 404 | 22 | 134 | 0.14 |
| protein:Kovarik2023_plasma | primary | 13 / 13 | 314 | 19 | 144 | 0.12 |
| protein:Fineschi2026 | primary | 92 / 73 | 358 | 25 | 116 | 0.08 |
| protein:Chowdhury2026 | primary | 53 / 62 | 472 | 30 | 202 | 0.15 |
| protein:Cervia2024_MS | sensitivity | 40 / 73 | 308 | 19 | 113 | 0.16 |

**Concordance of programme effects** (`concordance_<set>.csv`, `n_shared_<set>.csv`, `pairs_<set>.csv`). Smallest
number of shared programmes: 63 (`all`, Fineschi2026 with Cervia2024_MS), 15 (`hallmark`, Fineschi2026 with
Kovarik2023_plasma and with Cervia2024_MS), 48 (`reactome`, Fineschi2026 with Cervia2024_MS); no pair is flagged.

**Invariance partition** (`partition.csv`):

| analysis | environments | programmes partitioned (scored in every environment) | invariant | null | unstable | none | not_partitioned |
|---|---|---|---|---|---|---|---|
| primary | 15 | 63 | 0 | 0 | 63 | 0 | 1854 |
| sensitivity_Cervia2024_MS | 15 | 55 | 0 | 0 | 55 | 0 | 1862 |

**Cross-layer leave one cohort out** (`loso.csv`). PRJNA1149196_QLD is skipped as a test cohort (`skipped`: fewer than
5 test samples in a class) and stays in training:

| scheme | test | role | positive / negative | programmes | AUROC identity (95% CI) | AUROC center_scale (95% CI) |
|---|---|---|---|---|---|---|
| rna_to_protein | protein:Buggert2025 | primary | 146 / 51 | 959 | 0.460 (0.377 to 0.544) | 0.434 (0.350 to 0.514) |
| rna_to_protein | protein:Cervia2024_SomaScan | primary | 40 / 73 | 1466 | 0.485 (0.367 to 0.602) | 0.460 (0.346 to 0.577) |
| rna_to_protein | protein:Su2022 | primary | 75 / 50 | 210 | 0.639 (0.541 to 0.733) | 0.610 (0.516 to 0.704) |
| rna_to_protein | protein:Chowdhury2024 | primary | 21 / 9 | 144 | 0.630 (0.360 to 0.868) | 0.587 (0.365 to 0.794) |
| rna_to_protein | protein:Kovarik2023_plasma | primary | 13 / 13 | 153 | 0.456 (0.225 to 0.692) | 0.556 (0.314 to 0.781) |
| rna_to_protein | protein:Fineschi2026 | primary | 92 / 73 | 141 | 0.493 (0.405 to 0.589) | 0.514 (0.426 to 0.602) |
| rna_to_protein | protein:Chowdhury2026 | primary | 53 / 62 | 216 | 0.527 (0.422 to 0.631) | 0.520 (0.413 to 0.627) |
| rna_to_protein | protein:Cervia2024_MS | sensitivity | 40 / 73 | 121 | 0.461 (0.345 to 0.575) | 0.416 (0.300 to 0.528) |
| protein_to_rna | rna:GSE226260_WB | primary | 17 / 20 | 63 | 0.750 (0.565 to 0.909) | 0.724 (0.524 to 0.868) |
| protein_to_rna | rna:GSE226260_PBMC | primary | 28 / 8 | 63 | 0.438 (0.232 to 0.661) | 0.446 (0.236 to 0.674) |
| protein_to_rna | rna:GSE224615 | primary | 23 / 13 | 63 | 0.545 (0.328 to 0.779) | 0.559 (0.334 to 0.783) |
| protein_to_rna | rna:GSE251849 | primary | 11 / 5 | 63 | 0.364 (0.055 to 0.691) | 0.418 (0.145 to 0.709) |
| protein_to_rna | rna:GSE270045 | primary | 19 / 17 | 63 | 0.885 (0.762 to 0.975) | 0.882 (0.765 to 0.978) |
| protein_to_rna | rna:PRJNA1149196_SA | primary | 28 / 11 | 63 | 0.568 (0.344 to 0.776) | 0.526 (0.289 to 0.747) |
| protein_to_rna | rna:PRJNA1149196_QLD | primary | 5 / 4 | 63 | fewer than 5 test samples in a class | fewer than 5 test samples in a class |
| protein_to_rna | rna:GSE297499 | primary | 17 / 9 | 63 | 0.392 (0.163 to 0.647) | 0.438 (0.196 to 0.660) |
| leave_one_cohort_out | rna:GSE226260_WB | primary | 17 / 20 | 63 | 0.776 (0.594 to 0.915) | 0.765 (0.588 to 0.897) |
| leave_one_cohort_out | rna:GSE226260_PBMC | primary | 28 / 8 | 63 | 0.429 (0.241 to 0.638) | 0.438 (0.245 to 0.652) |
| leave_one_cohort_out | rna:GSE224615 | primary | 23 / 13 | 63 | 0.602 (0.375 to 0.809) | 0.532 (0.308 to 0.746) |
| leave_one_cohort_out | rna:GSE251849 | primary | 11 / 5 | 63 | 0.382 (0.091 to 0.691) | 0.455 (0.145 to 0.782) |
| leave_one_cohort_out | rna:GSE270045 | primary | 19 / 17 | 63 | 0.876 (0.749 to 0.966) | 0.845 (0.703 to 0.954) |
| leave_one_cohort_out | rna:PRJNA1149196_SA | primary | 28 / 11 | 63 | 0.562 (0.325 to 0.779) | 0.555 (0.325 to 0.773) |
| leave_one_cohort_out | rna:PRJNA1149196_QLD | primary | 5 / 4 | 63 | fewer than 5 test samples in a class | fewer than 5 test samples in a class |
| leave_one_cohort_out | rna:GSE297499 | primary | 17 / 9 | 63 | 0.346 (0.144 to 0.562) | 0.503 (0.248 to 0.719) |
| leave_one_cohort_out | protein:Buggert2025 | primary | 146 / 51 | 63 | 0.638 (0.547 to 0.721) | 0.513 (0.429 to 0.598) |
| leave_one_cohort_out | protein:Cervia2024_SomaScan | primary | 40 / 73 | 63 | 0.573 (0.472 to 0.680) | 0.567 (0.465 to 0.671) |
| leave_one_cohort_out | protein:Su2022 | primary | 75 / 50 | 63 | 0.489 (0.386 to 0.586) | 0.448 (0.347 to 0.550) |
| leave_one_cohort_out | protein:Chowdhury2024 | primary | 21 / 9 | 63 | 0.661 (0.450 to 0.841) | 0.471 (0.259 to 0.677) |
| leave_one_cohort_out | protein:Kovarik2023_plasma | primary | 13 / 13 | 63 | 0.456 (0.219 to 0.687) | 0.580 (0.337 to 0.805) |
| leave_one_cohort_out | protein:Fineschi2026 | primary | 92 / 73 | 63 | 0.592 (0.504 to 0.676) | 0.584 (0.498 to 0.666) |
| leave_one_cohort_out | protein:Chowdhury2026 | primary | 53 / 62 | 63 | 0.564 (0.458 to 0.667) | 0.535 (0.440 to 0.636) |

**Bridge check in Cervia2024** (`bridge_cervia.csv`), identical to `mosaic_programmes/lc` (the same two Cervia2024
layers):

| member set | collection | programmes scored in both layers | median Spearman | minimum | maximum | programmes with Spearman p < 0.05 |
|---|---|---|---|---|---|---|
| own | hallmark | 19 | 0.09 | -0.15 | 0.63 | 4 |
| own | reactome | 113 | 0.14 | -0.15 | 0.73 | 47 |
| shared | hallmark | 16 | 0.34 | 0.05 | 0.58 | 13 |
| shared | reactome | 84 | 0.34 | -0.02 | 0.71 | 62 |

**Not committed:** `scores.csv.gz` and the log.

**5. Reporting fixes.** All concordance outputs listed in "Runs" were recomputed with the new code ("reporting rerun");
the effects and the all-feature matrices are identical to the previous files in every folder.
- **`spearman_top<N>.csv` minimum.** A pair of cohorts that shares fewer than 10 of the `top<N>` features is now NaN
  (correlation and p value). `n_shared_top<N>.csv` gives those counts. Changed cells: only in
  `lc_metabolomics/concordance_all` (the six pairs that had 1 to 4 shared identifiers, before plus or minus 1 or a value
  on 2 to 4 points). In `concordance_proteomics_all` and `concordance_proteomics_all_v4` the off-diagonal top-50 values
  were already NaN (0 or 1 shared) and stay NaN. The other concordance folders have no pair below 10.
- **Pairs with fewer than 10 shared features flagged.** Every consumer of the shared counts writes a per-pair table
  with the column `fewer_than_10_shared`: `pairs.csv` in every concordance folder (feature sets `all` and `top<N>`),
  `pairs_<set>.csv` in `mosaic_programmes/<contrast>/` and `mosaic_programmes_v4/lc_v4/`, and the column in
  `paired_Cervia2024/lc/effect_correlation.csv`. The all-feature matrices keep every pair. Flagged pairs over all
  features:
  - `concordance_proteomics_all`: Su2022 with Cervia2024_MS (4), Chowdhury2024 (5), Kovarik2023_plasma (4).
  - `concordance_proteomics_all_v4`: the same three, and Fineschi2026 with Cervia2024_MS (5), Chowdhury2024 (4),
    Kovarik2023_plasma (2) and Chowdhury2026 (5).
  - None in the other concordance folders (`concordance_v4`: 10,053 for every pair), in `mosaic_programmes` (smallest 16
    programmes, Hallmark, `lc`), in `mosaic_programmes_v4/lc_v4` (smallest 15, Hallmark) or in `paired_Cervia2024` (173).
- **`cohort_annotation.csv`.** The column `control_type` (which held the `lc_status` of the negative class) is renamed
  `negative_class_lc_status`. A new column `control_type` gives the schema `control_type` of the negative class from
  the metadata: `convalescent_unverified` for Chowdhury2024 and `COVID_recovered` for the other proteomic and
  metabolomic cohorts, and for PRJNA1149196_SA, PRJNA1149196_QLD and GSE297499 in `concordance_v4`. The five earlier
  transcriptomic cohorts (`concordance`, `concordance_v4`) show `missing`: their metadata files were built before the
  schema had `control_type` and were not rebuilt.
- **Classifiers.** The `classifier` row of the loso_results.csv glossary now defines all four. `taco` ran without
  TabPFN (`CHORUS_DISABLE_TABPFN=1`) in every real-data run; this is stated wherever a `taco` value is reported (section 3).

## Follow-up analyses (17 November 2026)

Four analyses on data already in the repository; no existing benchmark was rerun. Each wrote to a new folder. The runs
started between 09:19 and 09:26 and ran at the same time, so their runtimes overlap. "Label free harmonisers" here means
identity, center_scale, quantile_ref, fsqn and rank: they are fitted without labels and applied to a new cohort in
addon mode. combat_addon and chorus use the training labels, and combat_pooled is transductive.

### A. Early detection baseline (`early_detection/`)

The later outcome, LC versus COVID_recovered, predicted from a sample drawn in the acute phase. The samples have
`sample_phase` acute (GSE215865 has no `sample_phase`, and every library there was drawn during the hospital stay)
and a later `lc_status` of LC or COVID_recovered. Each subject contributes one sample, the earliest. The script writes
the subjects per group to `cohorts.csv` before fitting anything, and stops when a group has fewer than 10. No group
did:

| cohort | analysis | LC subjects | COVID_recovered subjects |
|---|---|---|---|
| Su2022 (Olink) | horizontal | 75 | 50 |
| Cervia2024_SomaScan | horizontal | 40 | 73 |
| Su2022 (Olink and metabolomics, paired) | vertical | 74 | 50 |
| GSE215865 (timepoint 0) | reference | 241 | 104 |

**1. Horizontal** (`early_detection/horizontal/`, `scripts/early_detection.py horizontal`,
`configs/early_detection_horizontal.yaml`).
- **Data:** Su2022 and Cervia2024_SomaScan on the 386 proteins they share, imputed with each cohort's median.
- **Design:** trained on one cohort and tested on the other, in both directions (`run_loso` with
  `min_train_studies: 1`). Harmonisers: identity, center_scale, rank, fsqn, combat_addon and chorus. Selectors:
  differential and invariant, with k = 50. Classifiers: logreg and rf. Intervals come from 1,000 stratified
  bootstrap resamples.
- **chorus:** fails in both directions ("invariance anchored harmonisation needs at least two labelled training
  studies"); its rows are error rows.
- **invariant selector:** with one training cohort, no protein passes the heterogeneity screen
  (`n_strict_invariant` 0). The selector then falls back to ranking proteins by effect. The two selectors give the
  same AUROC in 6 of the 20 harmoniser, classifier and direction combinations.
- **Label permutation control:** `lc_status` is shuffled within each cohort (class counts kept) and the grid is rerun
  100 times (`permutation_results.csv`). p = (1 + permutations with AUROC ≥ observed) / 101
  (`permutation_summary.csv`).
- **identity with rf, tested on Cervia2024_SomaScan:** the model gives every test sample the same probability, so the
  AUROC is 0.500 and the interval is 0.50 to 0.50.

| test (train on the other) | harmoniser | selector | classifier | AUROC (95% interval) | permutation mean | p |
|---|---|---|---|---|---|---|
| Su2022 | identity | differential | logreg | 0.554 (0.45 to 0.66) | 0.501 | 0.119 |
| Su2022 | identity | differential | rf | 0.503 (0.44 to 0.57) | 0.494 | 0.426 |
| Su2022 | center_scale | differential | logreg | 0.564 (0.47 to 0.67) | 0.505 | 0.139 |
| Su2022 | center_scale | differential | rf | 0.608 (0.51 to 0.71) | 0.504 | 0.040 |
| Su2022 | rank | differential | logreg | 0.432 (0.33 to 0.54) | 0.506 | 0.891 |
| Su2022 | rank | differential | rf | 0.439 (0.34 to 0.53) | 0.509 | 0.851 |
| Su2022 | fsqn | differential | logreg | 0.530 (0.43 to 0.63) | 0.500 | 0.287 |
| Su2022 | fsqn | differential | rf | 0.580 (0.48 to 0.68) | 0.498 | 0.030 |
| Su2022 | combat_addon | differential | logreg | 0.566 (0.47 to 0.67) | 0.505 | 0.139 |
| Su2022 | combat_addon | differential | rf | 0.607 (0.51 to 0.71) | 0.504 | 0.040 |
| Su2022 | identity | invariant | logreg | 0.529 (0.42 to 0.63) | 0.504 | 0.277 |
| Su2022 | identity | invariant | rf | 0.491 (0.42 to 0.57) | 0.494 | 0.554 |
| Su2022 | center_scale | invariant | logreg | 0.574 (0.47 to 0.68) | 0.506 | 0.109 |
| Su2022 | center_scale | invariant | rf | 0.591 (0.49 to 0.69) | 0.505 | 0.069 |
| Su2022 | rank | invariant | logreg | 0.430 (0.33 to 0.54) | 0.508 | 0.931 |
| Su2022 | rank | invariant | rf | 0.471 (0.37 to 0.57) | 0.504 | 0.723 |
| Su2022 | fsqn | invariant | logreg | 0.507 (0.40 to 0.61) | 0.502 | 0.505 |
| Su2022 | fsqn | invariant | rf | 0.603 (0.51 to 0.71) | 0.504 | 0.020 |
| Su2022 | combat_addon | invariant | logreg | 0.578 (0.48 to 0.68) | 0.508 | 0.059 |
| Su2022 | combat_addon | invariant | rf | 0.603 (0.50 to 0.70) | 0.503 | 0.030 |
| Cervia2024_SomaScan | identity | differential | logreg | 0.644 (0.52 to 0.75) | 0.492 | 0.010 |
| Cervia2024_SomaScan | identity | differential | rf | 0.500 (0.50 to 0.50) | 0.499 | 0.723 |
| Cervia2024_SomaScan | center_scale | differential | logreg | 0.555 (0.43 to 0.67) | 0.495 | 0.119 |
| Cervia2024_SomaScan | center_scale | differential | rf | 0.448 (0.34 to 0.56) | 0.492 | 0.842 |
| Cervia2024_SomaScan | rank | differential | logreg | 0.564 (0.46 to 0.67) | 0.500 | 0.119 |
| Cervia2024_SomaScan | rank | differential | rf | 0.544 (0.43 to 0.65) | 0.493 | 0.178 |
| Cervia2024_SomaScan | fsqn | differential | logreg | 0.561 (0.45 to 0.67) | 0.493 | 0.079 |
| Cervia2024_SomaScan | fsqn | differential | rf | 0.484 (0.37 to 0.60) | 0.492 | 0.564 |
| Cervia2024_SomaScan | combat_addon | differential | logreg | 0.554 (0.43 to 0.66) | 0.495 | 0.129 |
| Cervia2024_SomaScan | combat_addon | differential | rf | 0.446 (0.33 to 0.56) | 0.491 | 0.832 |
| Cervia2024_SomaScan | identity | invariant | logreg | 0.644 (0.52 to 0.75) | 0.493 | 0.010 |
| Cervia2024_SomaScan | identity | invariant | rf | 0.500 (0.50 to 0.50) | 0.495 | 0.723 |
| Cervia2024_SomaScan | center_scale | invariant | logreg | 0.555 (0.43 to 0.67) | 0.498 | 0.109 |
| Cervia2024_SomaScan | center_scale | invariant | rf | 0.448 (0.34 to 0.56) | 0.494 | 0.812 |
| Cervia2024_SomaScan | rank | invariant | logreg | 0.535 (0.43 to 0.64) | 0.499 | 0.257 |
| Cervia2024_SomaScan | rank | invariant | rf | 0.546 (0.44 to 0.64) | 0.498 | 0.257 |
| Cervia2024_SomaScan | fsqn | invariant | logreg | 0.537 (0.42 to 0.65) | 0.493 | 0.238 |
| Cervia2024_SomaScan | fsqn | invariant | rf | 0.501 (0.39 to 0.61) | 0.498 | 0.485 |
| Cervia2024_SomaScan | combat_addon | invariant | logreg | 0.554 (0.43 to 0.66) | 0.497 | 0.119 |
| Cervia2024_SomaScan | combat_addon | invariant | rf | 0.446 (0.33 to 0.56) | 0.493 | 0.802 |

**2. Vertical** (`early_detection/vertical/early/`, `scripts/paired_layers.py --cohort Su2022 --contrasts early`).
- **Contrast:** the new `early` contrast of `paired_layers.py`, LC versus COVID_recovered on acute samples. Each
  subject contributes its earliest acute draw, T1 before T2.
- **Data:** 124 subjects measured in both layers. All 441 Olink proteins and all 891 metabolites are used; none was
  dropped for missingness.
- **Models and folds:** the Test 2 models (single layers, early_fusion, late_fusion, mofa, diablo) on the same
  repeated stratified 5-fold CV, 10 repeats, with k = 50 features per layer. Means are over the 50 folds.

| model | AUROC (mean over 50 folds) | SD | AUPRC | Brier | failed folds |
|---|---|---|---|---|---|
| layer:olink | 0.481 | 0.107 | 0.616 | 0.359 | 0 |
| layer:metabolomics | 0.490 | 0.128 | 0.653 | 0.359 | 0 |
| early_fusion | 0.449 | 0.118 | 0.614 | 0.401 | 0 |
| late_fusion | 0.516 | 0.109 | 0.649 | 0.244 | 0 |
| mofa | 0.439 | 0.112 | 0.620 | 0.276 | 0 |
| diablo | 0.442 | 0.118 | 0.610 | 0.318 | 0 |

**3. Reference** (`early_detection/reference/`, `scripts/early_detection.py reference`,
`configs/early_detection_reference.yaml`).
- **Samples:** GSE215865, the earliest library of each subject, kept only when it is the timepoint 0 library. This uses
  the new `keep_timepoints` preprocessing option.
- **Genes:** 15,116 genes pass the expression filter on these samples.
- **Models:** within-cohort repeated stratified 5-fold CV, 10 repeats, with the differential selector (k = 200) and
  logreg and rf. The invariant selector is not run because within one cohort it is the differential selector.
- **Published value:** the TACO AUROC (0.675) is in the `published_taco_auroc` column. It is quoted, not recomputed.

| classifier | selector | AUROC (mean over 50 folds) | SD | pooled AUROC (mean over repeats) | published TACO |
|---|---|---|---|---|---|
| logreg | differential | 0.512 | 0.068 | 0.514 | 0.675 |
| rf | differential | 0.569 | 0.060 | 0.567 | 0.675 |

### B. Concordance and transfer (`concordance_vs_transfer/`)

`scripts/concordance_vs_transfer.py` reads the `loso_results.csv` of `lc_transcriptomics_blood_v4`,
`lc_proteomics_affinity_v4`, `lc_proteomics_ms_v4` and `lc_metabolomics`.

- **Rows:** every held-out cohort with the invariant selector, logreg and a label free harmoniser (no error rows).
  quantile_ref was run only in the transcriptomic block. This gives 79 rows over 18 held-out cohorts (`rows.csv`).
- **What is related:** each row's AUROC to its `effect_spearman` and `sign_agreement`. Both are computed from the
  held-out labels after harmonisation and are diagnostics, not model inputs.
- **Statistics:** Spearman rho with two 2,000-resample bootstrap intervals. The cohort bootstrap resamples held-out
  cohorts together with all their rows; the row bootstrap resamples rows. In the pooled analysis both resample within
  block. The permutation p value comes from 10,000 shuffles of the metric across rows (within block for pooled),
  two-sided on |rho|, with a smallest possible value of 1/10,001.

| analysis | metric | rows | cohorts | Spearman | cohort bootstrap | row bootstrap | p (permutation) |
|---|---|---|---|---|---|---|---|
| transcriptomics | effect_spearman | 35 | 7 | 0.356 | -0.25 to 0.68 | 0.03 to 0.62 | 0.0390 |
| transcriptomics | sign_agreement | 35 | 7 | 0.266 | -0.27 to 0.50 | -0.06 to 0.55 | 0.1250 |
| proteomics_affinity | effect_spearman | 16 | 4 | 0.494 | -0.24 to 0.98 | -0.12 to 0.90 | 0.0570 |
| proteomics_affinity | sign_agreement | 16 | 4 | 0.189 | 0.06 to 0.87 | -0.47 to 0.73 | 0.4823 |
| proteomics_ms | effect_spearman | 16 | 4 | 0.841 | -0.17 to 0.98 | 0.45 to 0.97 | 0.0002 |
| proteomics_ms | sign_agreement | 16 | 4 | 0.892 | 0.25 to 0.96 | 0.62 to 0.99 | 0.0001 |
| metabolomics | effect_spearman | 12 | 3 | -0.015 | -0.45 to 0.77 | -0.64 to 0.71 | 0.9680 |
| metabolomics | sign_agreement | 12 | 3 | 0.000 | -0.38 to 0.74 | -0.63 to 0.72 | 1.0000 |
| pooled | effect_spearman | 79 | 18 | 0.511 | 0.25 to 0.68 | 0.35 to 0.66 | 0.0001 |
| pooled | sign_agreement | 79 | 18 | 0.426 | 0.23 to 0.56 | 0.24 to 0.60 | 0.0004 |

- **Variance shares** (`variance.csv`): a two-way ANOVA of AUROC on held-out cohort and harmoniser, with type II sums
  of squares divided by the total sum of squares. The pooled design is unbalanced because quantile_ref appears in one
  block only. The F tests and p values are in the file.

| analysis | held-out cohort | harmoniser | residual |
|---|---|---|---|
| transcriptomics | 0.539 | 0.082 | 0.379 |
| proteomics_affinity | 0.291 | 0.462 | 0.247 |
| proteomics_ms | 0.766 | 0.125 | 0.109 |
| metabolomics | 0.495 | 0.213 | 0.292 |
| pooled | 0.635 | 0.011 | 0.323 |

### C. Clustering readiness (`clustering_readiness/`)

`scripts/clustering_readiness.py` uses LC cases only, and the labels are not used otherwise. No gene list is used,
including the MRCONTROL list.

- **Blocks:**
  - RNA: the cohorts of `lc_transcriptomics_blood_v4` on their 10,053 shared genes. 148 LC samples: GSE226260_WB 17,
    GSE226260_PBMC 28, GSE224615 23, GSE251849 11, GSE270045 19, PRJNA1149196_SA 28, PRJNA1149196_QLD 5,
    GSE297499 17.
  - Affinity proteins: the cohorts of `lc_proteomics_affinity_v4` on their 166 shared proteins. 353 LC samples:
    Buggert2025 146, Cervia2024_SomaScan 40, Su2022 75, Fineschi2026 92.
- **Steps, per label free harmoniser:**
  1. The harmoniser is fitted on all LC samples and applied to them.
  2. RNA keeps the 2,000 genes with the largest variance after harmonisation; proteins keep all features. Features are
     standardised and reduced to 20 principal components.
  3. Clustering for k = 2 to 6 with three methods: k-means (50 starts); hierarchical (Ward); and consensus (100 k-means
     runs on 80% subsamples, then average linkage on 1 − consensus).
- **Outputs:**
  - Per configuration (`clusters.csv`): the ARI between cluster and cohort, normalised mutual information, silhouette
    and the smallest cluster. For consensus clustering it also gives the PAC: the fraction of sample pairs with a
    consensus value between 0.1 and 0.9.
  - Stability (`stability.csv`, `stability_summary.csv`): each cohort is left out in turn and harmoniser, features,
    components and clusters are refitted. Stability is the ARI between these clusters and the full-run clusters on the
    remaining samples.
  - Cluster-by-cohort counts: `composition.csv`.
- **Tables below:** k-means only; the other methods are in the files.
- **Ranges over all three methods and k:**
  - RNA, ARI with cohort: center_scale 0.00 to 0.04, fsqn 0.00 to 0.02, identity, quantile_ref and rank 0.24 to 0.89.
  - RNA, mean stability: identity 0.80 to 0.96, quantile_ref 0.84 to 0.98, rank 0.82 to 0.98, center_scale 0.46 to
    0.79, fsqn 0.27 to 0.84.
  - Proteins, ARI with cohort: center_scale −0.02 to 0.07, fsqn −0.02 to 0.03, identity 0.19 to 0.75, quantile_ref 0.17
    to 0.79, rank 0.08 to 0.87.
  - Proteins, mean stability: 0.13 to 0.96, with the lowest values from hierarchical clustering after center_scale and
    fsqn.

RNA, ARI between k-means clusters and cohort:

| harmoniser | k=2 | k=3 | k=4 | k=5 | k=6 |
|---|---|---|---|---|---|
| identity | 0.27 | 0.44 | 0.53 | 0.71 | 0.89 |
| center_scale | 0.01 | 0.00 | 0.01 | 0.01 | 0.03 |
| quantile_ref | 0.24 | 0.48 | 0.59 | 0.77 | 0.89 |
| fsqn | 0.01 | -0.00 | 0.02 | 0.01 | 0.02 |
| rank | 0.24 | 0.48 | 0.59 | 0.77 | 0.89 |

RNA, stability (mean ARI with the full run over the eight left-out cohorts), k-means:

| harmoniser | k=2 | k=3 | k=4 | k=5 | k=6 |
|---|---|---|---|---|---|
| identity | 0.83 | 0.96 | 0.83 | 0.82 | 0.94 |
| center_scale | 0.51 | 0.60 | 0.79 | 0.67 | 0.47 |
| quantile_ref | 0.98 | 0.92 | 0.84 | 0.91 | 0.94 |
| fsqn | 0.83 | 0.65 | 0.56 | 0.59 | 0.53 |
| rank | 0.89 | 0.95 | 0.87 | 0.95 | 0.94 |

RNA, PAC of consensus clustering:

| harmoniser | k=2 | k=3 | k=4 | k=5 | k=6 |
|---|---|---|---|---|---|
| identity | 0.25 | 0.00 | 0.21 | 0.22 | 0.14 |
| center_scale | 0.70 | 0.47 | 0.36 | 0.41 | 0.47 |
| quantile_ref | 0.00 | 0.00 | 0.00 | 0.07 | 0.09 |
| fsqn | 0.07 | 0.58 | 0.62 | 0.60 | 0.59 |
| rank | 0.00 | 0.00 | 0.17 | 0.00 | 0.05 |

Affinity proteins, ARI between k-means clusters and cohort:

| harmoniser | k=2 | k=3 | k=4 | k=5 | k=6 |
|---|---|---|---|---|---|
| identity | 0.19 | 0.30 | 0.32 | 0.53 | 0.66 |
| center_scale | 0.01 | 0.04 | 0.01 | 0.04 | 0.04 |
| quantile_ref | 0.33 | 0.37 | 0.79 | 0.58 | 0.63 |
| fsqn | -0.02 | 0.02 | 0.03 | 0.01 | 0.01 |
| rank | 0.14 | 0.35 | 0.78 | 0.60 | 0.64 |

Affinity proteins, stability (mean over the four left-out cohorts), k-means:

| harmoniser | k=2 | k=3 | k=4 | k=5 | k=6 |
|---|---|---|---|---|---|
| identity | 0.75 | 0.65 | 0.63 | 0.69 | 0.78 |
| center_scale | 0.92 | 0.71 | 0.68 | 0.56 | 0.72 |
| quantile_ref | 0.59 | 0.64 | 0.57 | 0.80 | 0.77 |
| fsqn | 0.57 | 0.71 | 0.70 | 0.63 | 0.61 |
| rank | 0.51 | 0.62 | 0.58 | 0.82 | 0.80 |

Affinity proteins, PAC of consensus clustering:

| harmoniser | k=2 | k=3 | k=4 | k=5 | k=6 |
|---|---|---|---|---|---|
| identity | 0.00 | 0.13 | 0.15 | 0.10 | 0.03 |
| center_scale | 0.05 | 0.26 | 0.29 | 0.34 | 0.24 |
| quantile_ref | 0.31 | 0.09 | 0.02 | 0.02 | 0.07 |
| fsqn | 0.06 | 0.32 | 0.15 | 0.20 | 0.32 |
| rank | 0.24 | 0.04 | 0.04 | 0.06 | 0.03 |

### D. Persistence of the acute signature (`acute_persistence/`)

`scripts/acute_persistence.py`. The acute score is the logreg log odds of acute. In each cohort it is compared between
LC and COVID_recovered:
- Hedges g (LC minus recovered) with a 1,000-resample bootstrap interval, and AUROC with a stratified bootstrap
  interval.
- Where `days_since_infection` is recorded with more than one value in a cohort, the OLS coefficient of LC in
  score ~ LC + days, divided by the pooled within-group SD (`g_adjusted`).
- A DerSimonian–Laird random effects meta-analysis across the scored cohorts: of g, and of `g_days`, which uses
  `g_adjusted` where available and g elsewhere.

Subjects in the training samples (same `subject_id`) are removed from the targets. Files: `per_cohort.csv`,
`meta.csv`, `model.csv`, `sharing.csv`, `scores.csv` and `run_info.json`.

**Proteins** (`acute_persistence/protein/`).
- **Training:** the samples and 163 shared proteins of `control_acute_proteomics` (Cervia2024_SomaScan, Su2022,
  Zoodsma2022_Hannover; acute versus never_infected). The model is fitted once per label free harmoniser and never
  refitted.
- **Targets:** the post-acute LC and COVID_recovered samples of `configs/concordance_proteomics_all_v4.yaml`.
- **Removed as training subjects:**
  - every post-acute subject of Cervia2024_SomaScan (40 LC / 73 recovered) and Su2022 (75 / 50);
  - Cervia2024_MS (40 / 73), whose subjects are those of Cervia2024_SomaScan.
- **Sharing the proteins** (measuring at least half of the 163):
  - Buggert2025 measures 162 and Fineschi2026 84.
  - The MS cohorts measure 2 to 5: Cervia2024_MS 2, Chowdhury2024 3, Kovarik2023_plasma 2, Chowdhury2026 5.
- **Model proteins:** the training proteins were restricted to the 84 measured in both Buggert2025 and Fineschi2026.
  For each harmoniser, the invariant selector (k = 50) chose proteins after harmonisation of those 84. The harmoniser
  was then refitted on the 50 chosen proteins and logreg was fitted on them.
- **Scored cohorts:** Buggert2025 (146 LC / 51 recovered) and Fineschi2026 (92 / 73). Neither records
  `days_since_infection`, so `g_days` equals g.

| harmoniser | cohort | LC / recovered | g (95% interval) | AUROC (95% interval) | g adjusted for days |
|---|---|---|---|---|---|
| identity | Buggert2025 | 146 / 51 | -0.82 (-1.10 to -0.56) | 0.282 (0.22 to 0.35) |  |
| identity | Fineschi2026 | 92 / 73 | 0.36 (0.08 to 0.63) | 0.588 (0.50 to 0.67) |  |
| center_scale | Buggert2025 | 146 / 51 | -0.77 (-1.07 to -0.47) | 0.296 (0.22 to 0.38) |  |
| center_scale | Fineschi2026 | 92 / 73 | 0.41 (0.15 to 0.71) | 0.608 (0.53 to 0.69) |  |
| quantile_ref | Buggert2025 | 146 / 51 | -0.84 (-1.12 to -0.59) | 0.266 (0.20 to 0.34) |  |
| quantile_ref | Fineschi2026 | 92 / 73 | -0.23 (-0.53 to 0.05) | 0.435 (0.35 to 0.52) |  |
| fsqn | Buggert2025 | 146 / 51 | -0.67 (-1.04 to -0.33) | 0.321 (0.24 to 0.41) |  |
| fsqn | Fineschi2026 | 92 / 73 | 0.40 (0.11 to 0.70) | 0.606 (0.52 to 0.69) |  |
| rank | Buggert2025 | 146 / 51 | -0.61 (-0.96 to -0.30) | 0.337 (0.24 to 0.43) |  |
| rank | Fineschi2026 | 92 / 73 | -0.07 (-0.35 to 0.23) | 0.478 (0.39 to 0.57) |  |

Meta-analysis (k = 2):

| harmoniser | estimate | k | adjusted | pooled g (95% interval) | p | I2 (%) | tau2 |
|---|---|---|---|---|---|---|---|
| identity | g | 2 | 0 | -0.23 (-1.39 to 0.93) | 0.701 | 96.2 | 0.675 |
| identity | g_days | 2 | 0 | -0.23 (-1.39 to 0.93) | 0.701 | 96.2 | 0.675 |
| center_scale | g | 2 | 0 | -0.18 (-1.33 to 0.98) | 0.763 | 96.2 | 0.666 |
| center_scale | g_days | 2 | 0 | -0.18 (-1.33 to 0.98) | 0.763 | 96.2 | 0.666 |
| quantile_ref | g | 2 | 0 | -0.53 (-1.13 to 0.06) | 0.078 | 85.6 | 0.157 |
| quantile_ref | g_days | 2 | 0 | -0.53 (-1.13 to 0.06) | 0.078 | 85.6 | 0.157 |
| fsqn | g | 2 | 0 | -0.14 (-1.19 to 0.91) | 0.798 | 95.4 | 0.550 |
| fsqn | g_days | 2 | 0 | -0.14 (-1.19 to 0.91) | 0.798 | 95.4 | 0.550 |
| rank | g | 2 | 0 | -0.34 (-0.87 to 0.19) | 0.212 | 82.2 | 0.120 |
| rank | g_days | 2 | 0 | -0.34 (-0.87 to 0.19) | 0.212 | 82.2 | 0.120 |

**RNA programmes** (`acute_persistence/rna_programmes/`).
- **Model:** the mosaic acute model of `scripts/mosaic_programmes.py`. Programme scores of the acute environments are
  rna:GSE215865, protein:Cervia2024_SomaScan, protein:Su2022 and protein:Zoodsma2022_Hannover (acute versus
  never_infected). Harmonisers identity and center_scale, logreg on the 123 programmes scored in every training
  environment and every target cohort, fitted once.
- **Targets:** the RNA cohorts of `configs/mosaic_lc_rna_v4.yaml`. No subject was removed as a training subject.
- **Not scored:** PRJNA1149196_QLD (5 / 4, below 5 per group).
- **Days adjustment:** possible only in GSE226260_WB. GSE224615 records a single value (243) and PRJNA1149196_SA, after
  one sample per subject, a single value (308).

| harmoniser | cohort | LC / recovered | g (95% interval) | AUROC (95% interval) | g adjusted for days |
|---|---|---|---|---|---|
| identity | GSE226260_WB | 17 / 20 | -0.55 (-1.30 to 0.11) | 0.312 (0.15 to 0.50) | -0.68 |
| identity | GSE226260_PBMC | 28 / 8 | 0.67 (0.06 to 1.48) | 0.714 (0.54 to 0.87) |  |
| identity | GSE224615 | 23 / 13 | 0.14 (-0.50 to 0.75) | 0.478 (0.27 to 0.69) |  |
| identity | GSE251849 | 11 / 5 | -0.72 (-1.94 to 0.48) | 0.309 (0.07 to 0.64) |  |
| identity | GSE270045 | 19 / 17 | -0.13 (-0.75 to 0.59) | 0.480 (0.29 to 0.67) |  |
| identity | PRJNA1149196_SA | 28 / 11 | -0.05 (-0.91 to 0.84) | 0.494 (0.26 to 0.72) |  |
| identity | GSE297499 | 17 / 9 | -0.73 (-1.49 to -0.04) | 0.294 (0.12 to 0.52) |  |
| center_scale | GSE226260_WB | 17 / 20 | -0.17 (-0.83 to 0.42) | 0.444 (0.27 to 0.62) | -0.32 |
| center_scale | GSE226260_PBMC | 28 / 8 | 0.38 (-0.58 to 1.30) | 0.540 (0.30 to 0.75) |  |
| center_scale | GSE224615 | 23 / 13 | 0.13 (-0.64 to 0.69) | 0.435 (0.25 to 0.63) |  |
| center_scale | GSE251849 | 11 / 5 | -0.92 (-2.33 to 0.10) | 0.255 (0.04 to 0.56) |  |
| center_scale | GSE270045 | 19 / 17 | -0.21 (-0.82 to 0.46) | 0.467 (0.28 to 0.67) |  |
| center_scale | PRJNA1149196_SA | 28 / 11 | -0.08 (-0.88 to 0.71) | 0.497 (0.27 to 0.73) |  |
| center_scale | GSE297499 | 17 / 9 | -0.48 (-1.42 to 0.31) | 0.333 (0.14 to 0.57) |  |

Meta-analysis (k = 7):

| harmoniser | estimate | k | adjusted | pooled g (95% interval) | p | I2 (%) | tau2 |
|---|---|---|---|---|---|---|---|
| identity | g | 7 | 0 | -0.16 (-0.52 to 0.19) | 0.365 | 35.1 | 0.079 |
| identity | g_days | 7 | 1 | -0.18 (-0.56 to 0.19) | 0.329 | 41.2 | 0.102 |
| center_scale | g | 7 | 0 | -0.13 (-0.41 to 0.15) | 0.369 | 0.0 | 0.000 |
| center_scale | g_days | 7 | 1 | -0.16 (-0.43 to 0.12) | 0.275 | 0.0 | 0.000 |


## Version 5 (17 November 2026)

The cohorts added on 17 November (`docs/proteomics_cohorts.md`, `docs/count_matrices.md`) enter version 5.
- **Proteins:** Woodruff2023 (Olink Explore 3072) and Bansal2025 (Olink Target 96). Bansal2025 keeps the QC rule as
  built, 56 participants.
- **RNA:** CORSIP2026 (NanoString panel) is a test cohort only. GSE228320 (DLCO below 60% of predicted) is a
  sensitivity cohort only.
- **Grid:** all the v5 runs use a reduced grid.
  - Harmonisers: identity, center_scale, rank, fsqn, combat_addon, chorus, chorus_protected and nested (nested chooses
    among the first six).
  - Selectors: variance (k 100 for proteins, 500 for RNA), differential and invariant (k 50 for proteins, 200 for RNA).
  - Classifiers: logreg and rf; 1,000 bootstrap resamples.
- **Not rerun:** the main transcriptomic (`lc_transcriptomics_blood_v4`), MS (`lc_proteomics_ms_v4`) and metabolomics
  (`lc_metabolomics`) benchmarks have no new main cohort, so their version 4 results stand as version 5.
- **Run order:** the runs went one after another from 10:17 (`results/runtimes.tsv`).

### 1. Comparison of version 4 and version 5 per block

Leave-one-study-out results restricted to the reduced grid:
- The v4 affinity block has no chorus_protected, so it has 42 configurations rather than 48.
- "Median of mean AUROC" is the median over the configurations that have a result in every fold of their mean AUROC
  across held-out cohorts.
- The best configuration is harmoniser / selector / classifier.

| block | version | cohorts | held-out cohorts | shared features | configurations (reduced grid) | with a result in every fold | median of mean AUROC | best configuration | best mean AUROC (min to max) |
|---|---|---|---|---|---|---|---|---|---|
| affinity | v4 | 4 cohorts | 4 | 166 | 42 | 42 | 0.499 | rank / differential / rf | 0.606 (0.51 to 0.71) |
| affinity | v5 | + Woodruff2023, Bansal2025 | 6 | 55 | 48 | 36 | 0.479 | rank / variance / rf | 0.563 (0.50 to 0.70) |
| Olink Explore | v5 | new block | 3 | 358 | 48 | 42 | 0.429 | identity / variance / rf | 0.593 (0.47 to 0.67) |
| MS | v4 = v5 | no new cohort; not rerun | 4 | 179 | 42 | 42 | 0.470 | rank / invariant / logreg | 0.566 (0.38 to 0.75) |
| transcriptomics | v4 = v5 | no new main cohort; not rerun | 7 | 10053 | 48 | 34 | 0.436 | rank / variance / rf | 0.523 (0.39 to 0.65) |
| metabolomics | v4 = v5 | no new cohort; not rerun | 3 | 38 | 42 | 36 | 0.555 | center_scale / differential / logreg | 0.621 (0.58 to 0.67) |

Affinity block: median AUROC over the reduced-grid configurations, per held-out cohort.

| held_out | v4 | v5 |
|---|---|---|
| Bansal2025 |  | 0.566 |
| Buggert2025 | 0.484 | 0.328 |
| Cervia2024_SomaScan | 0.483 | 0.529 |
| Fineschi2026 | 0.556 | 0.501 |
| Su2022 | 0.477 | 0.460 |
| Woodruff2023 |  | 0.502 |

- **Affinity v5** (`lc_proteomics_affinity_v5`):
  - Samples: Buggert2025 146 LC / 51 recovered, Cervia2024_SomaScan 40 / 73, Su2022 75 / 50, Fineschi2026 92 / 73,
    Woodruff2023 93 / 27, Bansal2025 12 / 20.
  - Proteins: the six cohorts share 55, against 166 in v4, because Bansal2025 measures 182 (Target 96). The variance
    selector (k 100) therefore keeps all 55.
  - Imputation: Bansal2025's 165 vendor-masked values are filled with its median.
  - chorus fails in every fold ("only 4 to 10 null genes passed the screen").
  - chorus_protected fails in every fold ("no genes of the protected scores"): its neutrophil and lymphocyte scores are
    gene lists that do not exist in UniProt space. The same holds in the Explore block.
  - Nested choices: fsqn in the Su2022 and Woodruff2023 folds; rank in the Cervia2024_SomaScan and Fineschi2026 folds;
    center_scale or fsqn in the Buggert2025 fold; identity or fsqn in the Bansal2025 fold.
  - Within-study CV is in `within_study_cv.csv`.
- **Olink Explore v5** (`lc_proteomics_explore_v5`): Buggert2025, Fineschi2026 and Woodruff2023 on their 358 shared
  proteins, with Buggert2025's 197 vendor-masked values filled with its median.
  - chorus runs in every fold; chorus_protected fails in every fold.
  - Nested chooses fsqn or identity in the Buggert2025 fold, rank in the Fineschi2026 fold, and rank or fsqn in the
    Woodruff2023 fold.

### 2. Concordance and mosaic programmes

- **`concordance_proteomics_all_v5`** (`configs/concordance_proteomics_all_v5.yaml`): the eight v4 protein cohorts plus
  Woodruff2023 and Bansal2025, each on its own proteins.
  - 45 pairs. 11 share fewer than 10 proteins over all features: the 7 of v4, plus Bansal2025 with Cervia2024_MS (0),
    Chowdhury2024 (1), Kovarik2023_plasma (0) and Chowdhury2026 (0).
  - Among the top 50, 44 of the 45 pairs are below 10 and are NaN.
  - The new pairs over all shared proteins (Spearman of the per-protein Hedges g, n_shared):

    | pair | n_shared | Spearman |
    |---|---|---|
    | Buggert2025, Woodruff2023 | 2,903 | −0.295 |
    | Cervia2024_SomaScan, Woodruff2023 | 2,163 | 0.171 |
    | Su2022, Woodruff2023 | 433 | 0.115 |
    | Fineschi2026, Woodruff2023 | 358 | 0.103 |
    | Chowdhury2026, Woodruff2023 | 198 | 0.154 |
    | Chowdhury2024, Woodruff2023 | 160 | 0.198 |
    | Kovarik2023_plasma, Woodruff2023 | 129 | −0.156 |
    | Cervia2024_MS, Woodruff2023 | 127 | 0.300 |
    | Buggert2025, Bansal2025 | 179 | −0.181 |
    | Woodruff2023, Bansal2025 | 179 | 0.017 |
    | Cervia2024_SomaScan, Bansal2025 | 174 | 0.182 |
    | Su2022, Bansal2025 | 89 | 0.292 |
    | Fineschi2026, Bansal2025 | 76 | −0.053 |

    All pairs and their p values are in `pairs.csv`.
- **`concordance_v5`** (`configs/concordance_rna_v5.yaml`, `features: per_cohort`): the eight v4 RNA cohorts plus
  CORSIP2026 and GSE228320, each pair on the genes both cohorts express.
  - No pair has fewer than 10 shared genes over all features. CORSIP2026 shares 553 to 704 genes per pair; GSE228320
    shares 11,187 to 12,723.
  - Among the top 500, the 9 CORSIP2026 pairs share 0 or 1 gene and are NaN.
  - CORSIP2026 with the v4 cohorts: Spearman −0.237 (GSE251849) to 0.502 (GSE270045). GSE228320 with the v4 cohorts:
    −0.171 to 0.154. CORSIP2026 with GSE228320: 0.230.
- **`mosaic_programmes_v5/lc_v5`** (`--contrasts lc_v5`): the protein side of lc_v4 plus Woodruff2023 (50 Hallmark, 983
  Reactome programmes scored) and Bansal2025 (16 and 53). The RNA side is unchanged.
  - Two Hallmark pairs share fewer than 10 programmes: Bansal2025 with Kovarik2023_plasma (8) and with Cervia2024_MS (9).
  - Partition: 36 programmes are scored in every environment (63 in lc_v4), and all 36 are unstable.
  - Test AUROC (identity / center_scale):
    - Woodruff2023: rna_to_protein 0.745 / 0.771 on 972 programmes; leave-one-cohort-out 0.449 / 0.513 on 36.
    - Bansal2025: rna_to_protein 0.700 / 0.675 on 69 programmes; leave-one-cohort-out 0.462 / 0.508 on 36.
  - Medians over test cohorts (identity / center_scale): leave_one_cohort_out 0.545 / 0.554, protein_to_rna 0.616 /
    0.549, rna_to_protein 0.510 / 0.538.

### 3. Panel transfer to CORSIP2026 (`panel_transfer_CORSIP2026`, `panel_transfer_CORSIP2026_noUnverified`)

- **Design:** training on the eight v4 transcriptomic cohorts, all `train_only`, on the 495 genes that pass the
  expression filter in every cohort including CORSIP2026 (the genes CORSIP2026 measures). Tested on CORSIP2026 only,
  which is never used for training. Covariate and protected variable: sample_matrix, as in v4. The 202 and 497 features
  include the two sample_matrix columns.
- **chorus and chorus_protected** fail ("only 0 null genes passed the screen").
- **Nested choices:** combat_addon (differential), identity (invariant), identity (variance).
- **Second run:** without the two `convalescent_unverified` controls PX-34 and PX-47.

| harmoniser | selector | classifier | CORSIP2026, 8 LC / 29 recovered | CORSIP2026 without PX-34 and PX-47, 8 / 27 |
|---|---|---|---|---|
| identity | variance | logreg | 0.724 (0.51 to 0.89) | 0.736 (0.54 to 0.91) |
| identity | variance | rf | 0.631 (0.44 to 0.80) | 0.632 (0.46 to 0.80) |
| center_scale | variance | logreg | 0.828 (0.65 to 0.96) | 0.819 (0.63 to 0.95) |
| center_scale | variance | rf | 0.903 (0.77 to 0.99) | 0.928 (0.82 to 1.00) |
| rank | variance | logreg | 0.655 (0.45 to 0.84) | 0.667 (0.47 to 0.85) |
| rank | variance | rf | 0.534 (0.23 to 0.80) | 0.530 (0.24 to 0.79) |
| fsqn | variance | logreg | 0.724 (0.51 to 0.91) | 0.713 (0.51 to 0.88) |
| fsqn | variance | rf | 0.776 (0.55 to 0.94) | 0.785 (0.58 to 0.95) |
| combat_addon | variance | logreg | 0.810 (0.62 to 0.94) | 0.810 (0.64 to 0.94) |
| combat_addon | variance | rf | 0.838 (0.63 to 0.97) | 0.870 (0.70 to 0.99) |
| nested | variance | logreg | 0.724 (0.51 to 0.89) | 0.736 (0.54 to 0.91) |
| nested | variance | rf | 0.631 (0.44 to 0.80) | 0.632 (0.46 to 0.80) |
| identity | differential | logreg | 0.789 (0.58 to 0.94) | 0.801 (0.58 to 0.96) |
| identity | differential | rf | 0.769 (0.62 to 0.90) | 0.752 (0.61 to 0.89) |
| center_scale | differential | logreg | 0.802 (0.63 to 0.94) | 0.796 (0.62 to 0.94) |
| center_scale | differential | rf | 0.843 (0.65 to 0.97) | 0.863 (0.69 to 0.99) |
| rank | differential | logreg | 0.793 (0.62 to 0.94) | 0.787 (0.61 to 0.94) |
| rank | differential | rf | 0.300 (0.13 to 0.50) | 0.289 (0.12 to 0.48) |
| fsqn | differential | logreg | 0.595 (0.38 to 0.81) | 0.588 (0.37 to 0.81) |
| fsqn | differential | rf | 0.832 (0.65 to 0.97) | 0.866 (0.71 to 0.98) |
| combat_addon | differential | logreg | 0.797 (0.62 to 0.93) | 0.792 (0.61 to 0.94) |
| combat_addon | differential | rf | 0.884 (0.71 to 0.99) | 0.884 (0.71 to 1.00) |
| nested | differential | logreg | 0.797 (0.62 to 0.93) | 0.792 (0.61 to 0.94) |
| nested | differential | rf | 0.884 (0.71 to 0.99) | 0.884 (0.71 to 1.00) |
| identity | invariant | logreg | 0.690 (0.46 to 0.88) | 0.690 (0.47 to 0.89) |
| identity | invariant | rf | 0.625 (0.34 to 0.87) | 0.637 (0.36 to 0.89) |
| center_scale | invariant | logreg | 0.621 (0.39 to 0.83) | 0.616 (0.40 to 0.82) |
| center_scale | invariant | rf | 0.774 (0.57 to 0.94) | 0.794 (0.61 to 0.94) |
| rank | invariant | logreg | 0.547 (0.31 to 0.78) | 0.542 (0.31 to 0.77) |
| rank | invariant | rf | 0.543 (0.37 to 0.72) | 0.549 (0.37 to 0.75) |
| fsqn | invariant | logreg | 0.759 (0.53 to 0.96) | 0.759 (0.53 to 0.94) |
| fsqn | invariant | rf | 0.853 (0.70 to 0.97) | 0.882 (0.75 to 0.99) |
| combat_addon | invariant | logreg | 0.642 (0.41 to 0.83) | 0.639 (0.43 to 0.83) |
| combat_addon | invariant | rf | 0.797 (0.60 to 0.95) | 0.819 (0.64 to 0.97) |
| nested | invariant | logreg | 0.690 (0.46 to 0.88) | 0.690 (0.47 to 0.89) |
| nested | invariant | rf | 0.625 (0.34 to 0.87) | 0.637 (0.36 to 0.89) |

### 4. Sensitivity: GSE228320 as an extra held-out cohort (`lc_transcriptomics_blood_v4_GSE228320`)

- **Design:** training on the eight v4 transcriptomic cohorts, all `train_only`, on the 9,556 genes they share with
  GSE228320. Tested on GSE228320 only: LC means DLCO below 60% of predicted (`case_definition` pulmonary_sequelae), not
  a Long COVID definition.
- **Scope:** the v4 folds are not changed, and GSE228320 never trains a model. This is reported apart from the main
  benchmark.
- **chorus and chorus_protected** fail ("only 12 null genes passed the screen").
- **Nested choices:** identity (differential), fsqn (invariant), center_scale (variance).

| harmoniser | selector | classifier | GSE228320 AUROC (95% interval), 18 DLCO below 60 / 32 at 60 or above |
|---|---|---|---|
| identity | variance | logreg | 0.595 (0.42 to 0.77) |
| identity | variance | rf | 0.351 (0.19 to 0.52) |
| center_scale | variance | logreg | 0.564 (0.40 to 0.72) |
| center_scale | variance | rf | 0.470 (0.30 to 0.64) |
| rank | variance | logreg | 0.595 (0.44 to 0.74) |
| rank | variance | rf | 0.486 (0.30 to 0.66) |
| fsqn | variance | logreg | 0.458 (0.29 to 0.64) |
| fsqn | variance | rf | 0.576 (0.40 to 0.74) |
| combat_addon | variance | logreg | 0.488 (0.32 to 0.65) |
| combat_addon | variance | rf | 0.575 (0.40 to 0.74) |
| nested | variance | logreg | 0.564 (0.40 to 0.72) |
| nested | variance | rf | 0.470 (0.30 to 0.64) |
| identity | differential | logreg | 0.523 (0.37 to 0.68) |
| identity | differential | rf | 0.253 (0.12 to 0.40) |
| center_scale | differential | logreg | 0.564 (0.39 to 0.74) |
| center_scale | differential | rf | 0.550 (0.39 to 0.72) |
| rank | differential | logreg | 0.479 (0.29 to 0.67) |
| rank | differential | rf | 0.453 (0.28 to 0.61) |
| fsqn | differential | logreg | 0.517 (0.36 to 0.67) |
| fsqn | differential | rf | 0.543 (0.37 to 0.71) |
| combat_addon | differential | logreg | 0.573 (0.39 to 0.74) |
| combat_addon | differential | rf | 0.514 (0.36 to 0.67) |
| nested | differential | logreg | 0.523 (0.37 to 0.68) |
| nested | differential | rf | 0.253 (0.12 to 0.40) |
| identity | invariant | logreg | 0.495 (0.32 to 0.68) |
| identity | invariant | rf | 0.598 (0.42 to 0.78) |
| center_scale | invariant | logreg | 0.514 (0.36 to 0.68) |
| center_scale | invariant | rf | 0.562 (0.39 to 0.72) |
| rank | invariant | logreg | 0.465 (0.29 to 0.64) |
| rank | invariant | rf | 0.576 (0.42 to 0.73) |
| fsqn | invariant | logreg | 0.549 (0.38 to 0.71) |
| fsqn | invariant | rf | 0.648 (0.50 to 0.79) |
| combat_addon | invariant | logreg | 0.523 (0.37 to 0.69) |
| combat_addon | invariant | rf | 0.604 (0.42 to 0.76) |
| nested | invariant | logreg | 0.549 (0.38 to 0.71) |
| nested | invariant | rf | 0.648 (0.50 to 0.79) |

### 5. Follow-up analyses on the version 5 cohorts

- **B** (`concordance_vs_transfer_v5`): the blocks are transcriptomics v4, affinity v5, MS v4, metabolomics, and the
  new Explore v5 block.
  - The Explore block is analysed alone and left out of the pooled analysis (`--not-pooled`), because its three
    held-out cohorts are also held out in affinity v5.
  - Affinity v5 (24 rows, 6 cohorts): Spearman with effect_spearman 0.577 (cohort bootstrap −0.12 to 0.87, permutation
    p 0.005), with sign_agreement 0.669 (0.13 to 0.90, p 0.0006). In v4 these were 0.494 and 0.189.
  - Explore v5 (12 rows, 3 cohorts): 0.664 (0.21 to 0.95, p 0.022) and 0.931 (0.32 to 0.95, p 0.0002).
  - Pooled (87 rows, 20 cohorts): 0.503 (0.27 to 0.66) and 0.512 (0.33 to 0.64), both p 0.0001. In v4 these were 0.511
    and 0.426.
  - Variance share of the held-out cohort (and of the harmoniser):

    | analysis | held-out cohort | harmoniser |
    |---|---|---|
    | affinity v5 | 0.615 | 0.064 |
    | Explore v5 | 0.380 | 0.192 |
    | pooled | 0.645 | 0.012 |

    Transcriptomics, MS and metabolomics are unchanged from v4.
- **C** (`clustering_readiness_v5`): RNA (the v4 cohorts, as before), affinity v5 and Explore v5.
  - Why the RNA block is v4: CORSIP2026 and GSE228320 are test-only or sensitivity cohorts, and CORSIP2026 would cut
    the shared genes to 495.
  - The RNA results are identical to `clustering_readiness`.
  - LC samples: affinity v5 has 458 (Buggert2025 146, Cervia2024_SomaScan 40, Su2022 75, Fineschi2026 92, Woodruff2023
    93, Bansal2025 12); Explore v5 has 331.
  - ARI between cluster and cohort, range over the three methods and k = 2 to 6:

    | block | center_scale | fsqn | identity | quantile_ref | rank |
    |---|---|---|---|---|---|
    | affinity v5 | 0.00 to 0.02 | 0.00 to 0.01 | 0.13 to 0.96 | 0.34 to 0.91 | 0.33 to 0.87 |
    | Explore v5 | 0.01 to 0.07 | 0.00 to 0.02 | 0.01 to 0.79 | 0.13 to 0.90 | 0.12 to 0.90 |

  - Mean leave-one-cohort-out stability: 0.30 to 0.94 in affinity v5 and 0.31 to 0.91 in Explore v5.
- **D** (`acute_persistence_v5`): the protein targets are `configs/concordance_proteomics_all_v5.yaml`.
  - Cohorts that share the proteins and keep subjects: Buggert2025, Fineschi2026, Woodruff2023 (163 of the 163 training
    proteins) and Bansal2025 (82 of 163).
  - The model proteins were therefore restricted to the 55 measured in all four, against 84 in v4. The invariant
    selector picks 50 of them.
  - Every post-acute subject of Cervia2024 and Su2022 is still removed as a training subject. The MS cohorts measure 2
    to 5 of the proteins.
  - Woodruff2023 records `days_since_infection` (days after symptom onset) for all 120 scored samples, so `g_days` uses
    its adjusted g. Its g adjusted for days (identity / center_scale / quantile_ref / fsqn / rank) is 0.86 / 0.62 / 0.76
    / 0.80 / 0.29.
  - The RNA programme part is unchanged; its targets are the v4 RNA cohorts, and its output equals
    `acute_persistence/rna_programmes`.

| harmoniser | cohort | LC / recovered | g (95% interval) | AUROC (95% interval) | g adjusted for days |
|---|---|---|---|---|---|
| identity | Buggert2025 | 146 / 51 | -0.21 (-0.51 to 0.06) | 0.440 (0.35 to 0.52) |  |
| identity | Fineschi2026 | 92 / 73 | 0.10 (-0.18 to 0.41) | 0.517 (0.43 to 0.61) |  |
| identity | Woodruff2023 | 93 / 27 | 0.78 (0.50 to 1.11) | 0.712 (0.62 to 0.81) | 0.86 |
| identity | Bansal2025 | 12 / 20 | 0.15 (-0.62 to 0.80) | 0.500 (0.28 to 0.71) |  |
| center_scale | Buggert2025 | 146 / 51 | -0.62 (-0.92 to -0.35) | 0.318 (0.24 to 0.40) |  |
| center_scale | Fineschi2026 | 92 / 73 | 0.33 (0.06 to 0.61) | 0.593 (0.51 to 0.68) |  |
| center_scale | Woodruff2023 | 93 / 27 | 0.55 (0.20 to 0.94) | 0.661 (0.55 to 0.76) | 0.62 |
| center_scale | Bansal2025 | 12 / 20 | 0.39 (-0.28 to 1.06) | 0.583 (0.38 to 0.78) |  |
| quantile_ref | Buggert2025 | 146 / 51 | -0.37 (-0.61 to -0.14) | 0.406 (0.33 to 0.48) |  |
| quantile_ref | Fineschi2026 | 92 / 73 | -0.09 (-0.38 to 0.21) | 0.481 (0.39 to 0.57) |  |
| quantile_ref | Woodruff2023 | 93 / 27 | 0.69 (0.26 to 1.16) | 0.681 (0.56 to 0.79) | 0.76 |
| quantile_ref | Bansal2025 | 12 / 20 | -0.26 (-0.99 to 0.62) | 0.475 (0.27 to 0.70) |  |
| fsqn | Buggert2025 | 146 / 51 | -0.49 (-0.84 to -0.18) | 0.361 (0.28 to 0.44) |  |
| fsqn | Fineschi2026 | 92 / 73 | 0.41 (0.14 to 0.69) | 0.597 (0.52 to 0.68) |  |
| fsqn | Woodruff2023 | 93 / 27 | 0.72 (0.38 to 1.09) | 0.700 (0.59 to 0.80) | 0.79 |
| fsqn | Bansal2025 | 12 / 20 | 0.67 (-0.05 to 1.50) | 0.629 (0.40 to 0.83) |  |
| rank | Buggert2025 | 146 / 51 | -0.25 (-0.52 to 0.02) | 0.419 (0.34 to 0.50) |  |
| rank | Fineschi2026 | 92 / 73 | 0.01 (-0.29 to 0.31) | 0.504 (0.42 to 0.59) |  |
| rank | Woodruff2023 | 93 / 27 | 0.24 (-0.21 to 0.71) | 0.562 (0.43 to 0.70) | 0.28 |
| rank | Bansal2025 | 12 / 20 | -0.15 (-0.90 to 0.73) | 0.508 (0.30 to 0.75) |  |

Meta-analysis (k = 4):

| harmoniser | estimate | k | adjusted | pooled g (95% interval) | p | I2 (%) | tau2 |
|---|---|---|---|---|---|---|---|
| identity | g | 4 | 0 | 0.19 (-0.23 to 0.61) | 0.382 | 76.6 | 0.136 |
| identity | g_days | 4 | 1 | 0.22 (-0.25 to 0.69) | 0.365 | 81.5 | 0.178 |
| center_scale | g | 4 | 0 | 0.14 (-0.45 to 0.74) | 0.636 | 88.1 | 0.311 |
| center_scale | g_days | 4 | 1 | 0.16 (-0.45 to 0.77) | 0.600 | 88.9 | 0.332 |
| quantile_ref | g | 4 | 0 | -0.00 (-0.47 to 0.46) | 0.984 | 80.3 | 0.169 |
| quantile_ref | g_days | 4 | 1 | 0.02 (-0.48 to 0.51) | 0.952 | 83.5 | 0.207 |
| fsqn | g | 4 | 0 | 0.30 (-0.31 to 0.90) | 0.336 | 88.6 | 0.327 |
| fsqn | g_days | 4 | 1 | 0.32 (-0.31 to 0.94) | 0.318 | 89.5 | 0.351 |
| rank | g | 4 | 0 | -0.05 (-0.25 to 0.16) | 0.666 | 11.2 | 0.005 |
| rank | g_days | 4 | 1 | -0.03 (-0.26 to 0.20) | 0.785 | 24.5 | 0.013 |

**Acute check** (`acute_check.csv`):
- The protein model scores Woodruff2023's 29 ICU acute samples together with its 120 post-acute LC and COVID_recovered
  samples, harmonised as one cohort. The scores therefore differ from those in `per_cohort.csv`, where the 120 were
  harmonised alone.
- The acute samples are used only here.

| harmoniser | acute / post-acute | mean score acute | post-acute LC | post-acute recovered | g acute minus post-acute | AUROC acute vs post-acute (95% interval) |
|---|---|---|---|---|---|---|
| identity | 29 / 120 | 2.13 | -0.89 | -1.97 | 2.15 | 0.917 (0.84 to 0.97) |
| center_scale | 29 / 120 | 8.51 | 0.45 | -1.25 | 2.24 | 0.915 (0.82 to 0.97) |
| quantile_ref | 29 / 120 | 3.08 | -2.51 | -4.40 | 2.18 | 0.945 (0.90 to 0.98) |
| fsqn | 29 / 120 | 4.77 | 0.69 | -1.43 | 1.36 | 0.798 (0.68 to 0.90) |
| rank | 29 / 120 | 2.48 | -3.53 | -4.25 | 2.13 | 0.943 (0.89 to 0.98) |


## CORSIP2026 panel transfer: verification (17 November 2026)

Checks of `panel_transfer_CORSIP2026`; the existing results are read, not changed.
- **Script:** `scripts/panel_transfer_checks.py`. Outputs: `panel_transfer_CORSIP2026_checks/` and
  `lc_transcriptomics_blood_v4_panel495/`.
- **Configurations checked:** the five with the highest CORSIP2026 AUROC (top1 to top5) and the median configuration
  (the lower middle of the 36 that ran).
- **Refits:** each configuration is refitted exactly as the benchmark fits it. The refits reproduce the stored AUROCs
  (0.903 / 0.884 / 0.884 / 0.853 / 0.843 / 0.759; `refit_check.csv`).

**Does the result survive?**
- **Sex and age do not explain it.**
  - Neither differs significantly between LC and recovered in CORSIP2026.
  - Within each sex the five best configurations keep AUROC 0.76 to 0.96.
  - Adjusted for sex and age, the score of each of them stays associated with LC: odds ratio per SD 11.7 to 27.8, all
    intervals above 1, but very wide.
  - Their largest weights are autosomal immune genes, none sex-differential at FDR 0.05.
- **It beats the label permutation control:** p 0.010 to 0.025 for the five best and 0.030 for the median
  configuration. These configurations were chosen after the fact out of 36, and the p values are not corrected for
  that.
- **It is not a general property of the panel genes.** On the same 495 genes, the version 4 leave-one-study-out
  benchmark is within the range of 20 random sets of 495 expressed genes: median of mean AUROC 0.452, against
  0.382 to 0.496 for the random sets (17 of
  20 below the panel) and 0.424 on all genes.
- **It rests on one training cohort.**
  - Without GSE270045, the training cohort whose LC effects agree most with CORSIP2026 (Spearman 0.534 over the 495
    genes), the CORSIP2026 AUROC of the six configurations falls to 0.39 to 0.55.
  - Without any other single cohort it stays at 0.55 to 0.93.
- **It rests on 8 LC cases:** the interval of the best configuration is 0.77 to 0.99.

The transfer is real in these checks, but it should be reported as driven by GSE270045 and based on 8 cases, not as
evidence that the panel genes carry a signal that generalises across cohorts.

### 1. Confounding in CORSIP2026

| group | n | female / male | age mean (SD) | age range |
|---|---|---|---|---|
| LC | 8 | 5 / 3 | 43.1 (7.5) | 31 to 52 |
| COVID_recovered | 29 | 11 / 18 | 37.7 (8.4) | 24 to 59 |

Sex alone and age alone, as scores for LC (`corsip_single_predictors.csv`):

| predictor | AUROC for LC | test | p |
|---|---|---|---|
| sex (female = 1) | 0.623 | Fisher exact | 0.254 |
| age | 0.692 | Mann-Whitney | 0.104 |

Transfer AUROC within each sex (both classes have at least 3 samples in each sex; every configuration in
`within_sex_auroc.csv`):

| role | configuration | AUROC all (8 / 29) | female (5 LC / 11) | male (3 LC / 18) |
|---|---|---|---|---|
| top1 | center_scale / variance / rf | 0.903 | 0.882 | 0.917 |
| top2 | combat_addon / differential / rf | 0.884 | 0.818 | 0.963 |
| top3 | nested / differential / rf | 0.884 | 0.818 | 0.963 |
| top4 | fsqn / invariant / rf | 0.853 | 0.764 | 0.944 |
| top5 | center_scale / differential / rf | 0.843 | 0.800 | 0.907 |
| median | fsqn / invariant / logreg | 0.759 | 0.855 | 0.574 |

Logistic regression of LC on the standardised model score, sex and standardised age (`adjusted_association.csv`, every
configuration; Wald interval). With 8 events the estimates are unstable; 11 of the 36 configurations have p below 0.05.

| role | configuration | odds ratio per SD of score (95% interval) | p |
|---|---|---|---|
| top1 | center_scale / variance / rf | 18.38 (1.90 to 177.73) | 0.012 |
| top2 | combat_addon / differential / rf | 27.78 (1.84 to 420.41) | 0.016 |
| top3 | nested / differential / rf | 27.78 (1.84 to 420.41) | 0.016 |
| top4 | fsqn / invariant / rf | 11.66 (1.76 to 77.15) | 0.011 |
| top5 | center_scale / differential / rf | 15.67 (1.63 to 150.54) | 0.017 |
| median | fsqn / invariant / logreg | 3.38 (0.86 to 13.30) | 0.081 |

**Largest model weights** (`top_weights.csv`, 20 per configuration):
- No chrY gene is on the panel; 18 panel genes are on chrX.
- In the top 20 of the five best configurations there is no chrX or chrY gene. The median configuration has one, IL2RG
  (chrX).
- Sex difference in CORSIP2026: a Welch t test of log CPM, female minus male, Benjamini-Hochberg over the 495 genes
  (`corsip_sex_difference.csv`).
  - Two panel genes are sex-differential at FDR 0.05, NFKB2 and RELB. Neither is in any top 20.
  - Nominally (p below 0.05, FDR 0.39 to 0.61) the top 20 include PTGS2, CD274, NFKBIA, IFI35, IL1B, OAS3 and TNFAIP3
    in the best configurations, and SERPING1, PTGS2, IL2RG, LYN and IFIH1 in the median one.

| role | configuration | weight | 10 largest (symbol, chromosome) |
|---|---|---|---|
| top1 | center_scale / variance / rf | impurity importance | REL (chr2), PTGS2 (chr1), CD274 (chr9), CXCL8 (chr4), CCL4 (chr17), PSMB10 (chr16), ITGB2 (chr21), NUP107 (chr12), LAMP1 (chr13), TANK (chr2) |
| top2 | combat_addon / differential / rf | impurity importance | ITGB2 (chr21), REL (chr2), PTGS2 (chr1), LAMP1 (chr13), CCL4 (chr17), CD274 (chr9), USP8 (chr15), RIPK2 (chr8), GPI (chr19), IFI35 (chr17) |
| top3 | nested / differential / rf | impurity importance | ITGB2 (chr21), REL (chr2), PTGS2 (chr1), LAMP1 (chr13), CCL4 (chr17), CD274 (chr9), USP8 (chr15), RIPK2 (chr8), GPI (chr19), IFI35 (chr17) |
| top4 | fsqn / invariant / rf | impurity importance | NFKBIA (chr14), CXCL8 (chr4), CCL4 (chr17), PTGS2 (chr1), OAS3 (chr12), RIPK2 (chr8), IFI35 (chr17), CST3 (chr20), PSMB10 (chr16), MICB (chr6) |
| top5 | center_scale / differential / rf | impurity importance | CCL4 (chr17), ITGB2 (chr21), REL (chr2), IFI35 (chr17), RIPK2 (chr8), CD274 (chr9), GPI (chr19), LAMP1 (chr13), USP8 (chr15), PSMB10 (chr16) |
| median | fsqn / invariant / logreg | coefficient | CD96 (chr3), C3 (chr19), NKG7 (chr19), SLAMF1 (chr1), SERPING1 (chr11), PTGS2 (chr1), IL2RG (chrX), CCL5 (chr17), NFATC2 (chr20), SELL (chr1) |

### 2. Permutation control

`lc_status` is permuted within each training cohort (class counts kept; CORSIP2026 labels unchanged) 200 times, and the
checked configurations are refitted and score CORSIP2026 (`permutation_null.csv`, `permutation_summary.csv`).
- p = (1 + permutations with AUROC at or above the observed) / (1 + permutations).
- Permutation i uses the seed (0, i); the 200 ran as four parallel chunks.
- In 26 permutations nested chose chorus in its inner leave-one-study-out, and the outer chorus fit then failed ("only 4
  null genes passed the screen"). Those permutations have no nested value, so top3 has 174.

| role | configuration | observed AUROC | permutations with a result | null mean (SD) | null 5% / 50% / 95% | null max | p |
|---|---|---|---|---|---|---|---|
| top1 | center_scale / variance / rf | 0.903 | 200 | 0.497 (0.182) | 0.23 / 0.47 / 0.78 | 0.922 | 0.010 |
| top2 | combat_addon / differential / rf | 0.884 | 200 | 0.500 (0.192) | 0.22 / 0.48 / 0.80 | 0.927 | 0.015 |
| top3 | nested / differential / rf | 0.884 | 174 | 0.493 (0.167) | 0.22 / 0.50 / 0.75 | 0.948 | 0.011 |
| top4 | fsqn / invariant / rf | 0.853 | 200 | 0.495 (0.177) | 0.23 / 0.47 / 0.77 | 0.886 | 0.020 |
| top5 | center_scale / differential / rf | 0.843 | 200 | 0.497 (0.189) | 0.22 / 0.48 / 0.80 | 0.948 | 0.025 |
| median | fsqn / invariant / logreg | 0.759 | 200 | 0.501 (0.138) | 0.27 / 0.51 / 0.72 | 0.815 | 0.030 |

### 3. Panel-restricted benchmark (`lc_transcriptomics_blood_v4_panel495`)

- **Design:** the version 4 transcriptomic leave-one-study-out (`configs/lc_transcriptomics_blood_v4_panel495.yaml`) on
  the 495 genes, reduced version 5 grid, without within-study CV.
  - The 495 genes are listed in `configs/gene_lists/corsip2026_panel495.txt`, applied with the new `restrict_features`
    preprocessing option.
  - Log CPM uses each cohort's library size over all genes, so the values equal those of `lc_transcriptomics_blood_v4`.
- **chorus and chorus_protected** fail in every fold ("only 0 to 3 null genes passed the screen").
- **Random gene control:** 20 random sets of 495 genes drawn from the 10,053 genes expressed in every v4 cohort (set i
  with the seed (0, i); `random_genes/setNN/genes.txt` and `loso_results.csv`), same grid.
- **Comparison:** with `lc_transcriptomics_blood_v4` restricted to the same reduced grid. Values are medians over the
  configurations of each held-out cohort's AUROC (`comparison_by_held_out.csv`; median of mean AUROC per set in
  `random_genes/summary_by_set.csv`).

| held out | all genes: median (max) | panel 495: median (max) | 20 random sets: median, mean (range) | random sets below the panel |
|---|---|---|---|---|
| GSE224615 | 0.569 (0.689) | 0.603 (0.706) | 0.556 (0.491 to 0.617) | 18 |
| GSE226260_PBMC | 0.373 (0.683) | 0.485 (0.683) | 0.332 (0.267 to 0.410) | 20 |
| GSE226260_WB | 0.612 (0.737) | 0.605 (0.776) | 0.580 (0.452 to 0.712) | 12 |
| GSE251849 | 0.327 (0.536) | 0.436 (0.727) | 0.341 (0.245 to 0.436) | 18 |
| GSE270045 | 0.272 (0.686) | 0.217 (0.564) | 0.256 (0.057 to 0.397) | 6 |
| GSE297499 | 0.471 (0.650) | 0.343 (0.611) | 0.484 (0.382 to 0.582) | 0 |
| PRJNA1149196_SA | 0.385 (0.575) | 0.352 (0.497) | 0.425 (0.382 to 0.555) | 0 |

Median over configurations of the mean AUROC across held-out cohorts: all genes 0.424, panel 0.452, random sets
0.430 on average (0.382 to 0.496).

### 4. Concordance and leave one training cohort out

- **Concordance:** Hedges g of LC versus COVID_recovered per cohort on the 495 genes (`effects_495.csv`); Spearman
  correlation and sign agreement of CORSIP2026 with each training cohort (`concordance_with_corsip.csv`).
- **Leave one training cohort out:** the six checked configurations refitted without that cohort and CORSIP2026 scored
  (`leave_one_training_cohort_out.csv`).
- **Full-data AUROC** (top1 to top5, median): 0.903 / 0.884 / 0.884 / 0.853 / 0.843 / 0.759.
- Across the eight training cohorts, the Spearman correlation between concordance with CORSIP2026 and the mean AUROC
  change when the cohort is left out is -0.619 (p 0.102).

| training cohort left out | Spearman with CORSIP2026 (p) | sign agreement | AUROC without, top1 | AUROC without, top2 | AUROC without, top3 | AUROC without, top4 | AUROC without, top5 | AUROC without, median | mean change |
|---|---|---|---|---|---|---|---|---|---|
| GSE226260_WB | 0.090 (0.046) | 0.525 | 0.795 | 0.890 | 0.817 | 0.806 | 0.806 | 0.789 | -0.037 |
| GSE226260_PBMC | -0.086 (0.057) | 0.483 | 0.871 | 0.897 | 0.838 | 0.808 | 0.838 | 0.642 | -0.038 |
| GSE224615 | -0.011 (0.813) | 0.459 | 0.808 | 0.838 | 0.554 | 0.832 | 0.838 | 0.772 | -0.080 |
| GSE251849 | -0.249 (0.000) | 0.406 | 0.853 | 0.836 | 0.765 | 0.836 | 0.849 | 0.853 | -0.022 |
| GSE270045 | 0.534 (0.000) | 0.701 | 0.552 | 0.494 | 0.489 | 0.386 | 0.489 | 0.483 | -0.372 |
| PRJNA1149196_SA | 0.122 (0.007) | 0.545 | 0.849 | 0.845 | 0.845 | 0.756 | 0.862 | 0.711 | -0.043 |
| PRJNA1149196_QLD | -0.073 (0.106) | 0.493 | 0.888 | 0.931 | 0.815 | 0.843 | 0.892 | 0.728 | -0.005 |
| GSE297499 | 0.032 (0.477) | 0.511 | 0.806 | 0.821 | 0.823 | 0.761 | 0.823 | 0.659 | -0.072 |

## Cohort export and manifest (17 November 2026)

`scripts/export_cohorts.py` writes the analysis-ready LC versus COVID_recovered cohorts to `data/exports/` (not
committed) and the manifest to `exports/manifest.csv` (committed).
- **Files:** each cohort is written on its own features, as `<layer>/<cohort>_matrix.csv.gz` (features × samples) with
  `<layer>/<cohort>_meta.csv`.
- **Layers and configs:**
  - RNA: `configs/concordance_rna_v5.yaml`, the eight v4 cohorts plus CORSIP2026 and GSE228320. log2 CPM on each
    cohort's expressed genes; GSE297499 is log2(TPM + 1).
  - Proteomics: `configs/concordance_proteomics_all_v5.yaml`, the ten version 5 cohorts including Woodruff2023 and
    Bansal2025, imputed per cohort.
  - Metabolomics: `configs/concordance_metabolomics_all.yaml`.
- **Manifest:** one row per cohort with role, flags, sample and feature counts, value type, platform, accession,
  overlaps, source config and the SHA-256 of both files.
- **Roles:**
  - `train_test`: 20 cohorts.
  - `train_only`: PRJNA1149196_QLD.
  - `test_only`: CORSIP2026 and GSE228320.
  - `concordance_only`: Saito2024, which is in no leave-one-study-out benchmark.
- **Flags:** GSE270045 is `label_confounded_with_collection`: group is completely confounded with sample collection
  (`docs/gse270045_check.md`). No other cohort is flagged.
- **Run:** 2026-11-17, about 170 s. It was run twice; the second run only added the `concordance_only` role for Saito2024
  and wrote the same files.

## Further analyses: noise ceiling, guided training, full range (17 November 2026)

**Bars, stated before the analyses were run** (copied from the request; the results below are reported against them):

- **A. Noise ceiling of concordance** (`concordance_noise_ceiling/`). Disagreement between cohorts counts as real in a
  layer if more than half of the cohort pairs of that layer have an observed correlation below the 2.5th percentile
  of the correlation expected if both cohorts came from one population.
- **B. Concordance-guided training** (`concordance_guided/`).
  - Guided training is useful if guided selection (arm 2) beats both all training cohorts (arm 1) and random selection
    of the same number of cohorts (arm 4), with the interval of the paired difference above zero, in at least two
    blocks.
  - It adds value over local data if calibration plus guided selection (arm 7) beats calibration alone (arm 5) in at
    least two blocks.
- **C. Full range of concordance against transfer** (`concordance_vs_transfer_full/`). No numerical bar was given.
  Reported: the pooled Spearman with and without the positive-control points, and whether the LC points lie on the
  same monotone curve as the controls.

### Results against the bars

**A. Noise ceiling: the disagreement is not real by the bar in any layer.**

- **Bar:** more than half of a layer's pairs below the 2.5th percentile of the one-population expectation.
- **Fraction of pairs below it:** 18 of 45 in transcriptomics (0.40), 18 of 40 evaluable in proteomics (0.45), 1 of 6
  in metabolomics (0.17).
- **Proteomics:** five pairs share fewer than 3 proteins and are not evaluable.
- **Observed against expected:** in every layer the median observed correlation is below the median expected one.
- **Reliability:** several cohorts barely reproduce their own LC effect across random halves (split-half r at or below
  0.05: GSE224615, PRJNA1149196_SA, PRJNA1149196_QLD, Su2022_metabolomics; GSE228320 0.051). Pairs with such cohorts cannot
  agree beyond chance.

`scripts/concordance_noise_ceiling.py`, 200 splits, seed 0. Files: `cohorts.csv`, `pairs.csv` (observed, split-half
correlations and reliabilities on the shared features, disattenuated, expected mean, 2.5th and 97.5th percentiles,
below_2p5), `summary.csv`.
- **Disattenuation:** a disattenuated value is missing when a reliability on the shared features is not positive (22
  RNA and 9 protein pairs).
- **Expected correlation:** computed from the larger cohort, a stratified subsample of the smaller cohort's class counts
  (capped at half of each class) against the rest. A negative split-half r is possible when a cohort has almost no LC
  signal: the two halves partition the same samples.

| layer | pairs (evaluable) | below the 2.5th percentile | fraction | bar met (more than half) | median observed | median expected (one population) | median disattenuated |
|---|---|---|---|---|---|---|---|
| transcriptomics | 45 (45) | 18 | 0.400 | no | -0.032 | 0.171 | -0.108 |
| proteomics | 45 (40) | 18 | 0.450 | no | 0.102 | 0.373 | 0.164 |
| metabolomics | 6 (6) | 1 | 0.167 | no | 0.081 | 0.227 | 0.238 |

| layer | cohort | LC / recovered | features | split-half r (SD) | Spearman-Brown reliability |
|---|---|---|---|---|---|
| transcriptomics | GSE226260_WB | 17 / 20 | 12801 | 0.520 (0.146) | 0.684 |
| transcriptomics | GSE226260_PBMC | 28 / 8 | 13866 | 0.154 (0.089) | 0.267 |
| transcriptomics | GSE224615 | 23 / 13 | 15734 | -0.179 (0.165) | -0.435 |
| transcriptomics | GSE251849 | 11 / 5 | 13872 | 0.210 (0.104) | 0.347 |
| transcriptomics | GSE270045 | 19 / 17 | 15561 | 0.765 (0.039) | 0.867 |
| transcriptomics | PRJNA1149196_SA | 28 / 11 | 17751 | -0.049 (0.084) | -0.104 |
| transcriptomics | PRJNA1149196_QLD | 5 / 4 | 16055 | 0.001 (0.109) | 0.002 |
| transcriptomics | GSE297499 | 17 / 9 | 36499 | 0.254 (0.083) | 0.405 |
| transcriptomics | CORSIP2026 | 8 / 29 | 778 | 0.397 (0.078) | 0.568 |
| transcriptomics | GSE228320 | 18 / 32 | 14466 | 0.051 (0.136) | 0.097 |
| proteomics | Buggert2025 | 146 / 51 | 2903 | 0.832 (0.029) | 0.908 |
| proteomics | Cervia2024_SomaScan | 40 / 73 | 6402 | 0.133 (0.125) | 0.234 |
| proteomics | Su2022 | 75 / 50 | 441 | 0.124 (0.133) | 0.221 |
| proteomics | Cervia2024_MS | 40 / 73 | 308 | 0.397 (0.095) | 0.569 |
| proteomics | Chowdhury2024 | 21 / 9 | 404 | 0.277 (0.076) | 0.434 |
| proteomics | Kovarik2023_plasma | 13 / 13 | 314 | 0.566 (0.069) | 0.723 |
| proteomics | Fineschi2026 | 92 / 73 | 358 | 0.409 (0.124) | 0.581 |
| proteomics | Chowdhury2026 | 53 / 62 | 472 | 0.082 (0.097) | 0.151 |
| proteomics | Woodruff2023 | 93 / 27 | 2924 | 0.449 (0.057) | 0.620 |
| proteomics | Bansal2025 | 12 / 20 | 182 | 0.148 (0.091) | 0.257 |
| metabolomics | Martinez2024 | 40 / 40 | 511 | 0.539 (0.098) | 0.701 |
| metabolomics | Lopez2023 | 30 / 18 | 106 | 0.315 (0.120) | 0.480 |
| metabolomics | Saito2024 | 30 / 15 | 312 | 0.480 (0.098) | 0.648 |
| metabolomics | Su2022_metabolomics | 68 / 45 | 679 | 0.003 (0.090) | 0.007 |

**B. Concordance-guided training: the bar for guided training is not met (one block of the two required), and the bar
for value over local data is not met (no block).**

- **Arm 2 against 1 and 4:** guided selection beats all cohorts (arm 1) and random selection (arm 4), with both paired
  intervals above zero, only in proteomics_affinity_v5: +0.056 (0.008 to 0.108) and +0.040 (0.009 to 0.088).
  - In proteomics_ms_v4 the interval against random is above zero but the one against all cohorts is not.
  - In transcriptomics_v4 the differences are the largest (+0.127 and +0.114), but their intervals include zero.
- **Arm 7 against 5:** calibration plus guided cohorts is below the calibration part alone in every block (−0.022 to
  −0.086).
- **Calibration alone:** the calibration part alone (arm 5) has the highest mean AUROC of all arms in every block.

**Deviation from the specification** (decided on 18 November 2026 for run time; the bars are unchanged):
- **Grid:** harmonisers identity, center_scale and combat_addon only (no nested, rank, fsqn, chorus or
  chorus_protected); selectors invariant and differential (no variance); logreg and rf.
- **rf:** 200 trees on one core per fit (the benchmark uses 500).
- **Arm 4:** 5 random draws per split instead of 20.
- **Splits:** 20 calibration/test splits as specified.

Design (`scripts/concordance_guided.py`):
- **Per held-out cohort:** 20 stratified 50/50 calibration/test splits.
- **Agreement:** the Spearman correlation, over the block's features, of the calibration part's Hedges g with each
  training cohort's g.
- **Harmonisation:** the harmoniser is fitted on the arm's training studies and applied to the test part as a new
  cohort, without labels. The calibration part is its own study.
- **AUROC:** per split, the mean over the 12 configurations (arm 4 also over its draws); per cohort, the mean over
  splits; per block, the mean over cohorts.
- **Intervals:** 2,000-resample cohort bootstraps. With 3 to 6 cohorts per block these intervals are coarse.
- **Paired differences:** use the splits where both arms are defined.

Held-out cohorts:
- **Skipped** (calibration part with fewer than 4 samples per class): GSE251849 (6 LC / 2 recovered) and
  PRJNA1149196_QLD (2 / 2; QLD is not held out in the benchmark either).
- **Empty guided selection** (no training cohort with positive agreement): Lopez2023 (metabolomics) 1 of 20; Su2022 (proteomics_affinity_v4) 5 of 20; Buggert2025 (proteomics_affinity_v5) 6 of 20; Woodruff2023 (proteomics_explore_v5) 8 of 20; Chowdhury2026 (proteomics_ms_v4) 5 of 20; Kovarik2023_plasma (proteomics_ms_v4) 20 of 20; GSE226260_PBMC (transcriptomics_v4) 2 of 20. In every other split at least one
  cohort agrees.
- **Size of the guided selection:** mean 1.35 (Explore v5) to 2.88 (transcriptomics) training cohorts.
- **Errors:** no fit failed (64,128 configuration results).

Files: `arms.csv`, `differences.csv`, `bars.csv`, `empty.csv`, `skipped.csv`, `per_cohort.csv`, `per_split.csv`,
`agreement.csv`, `units/` (every fit), and `blocks/<block>/` (the same summaries per block, written as each block
finished).

Mean AUROC on the test parts (95% cohort bootstrap):

| block | cohorts | 1 all | 2 guided | 3 weighted | 4 random | 5 calib | 6 calib + all | 7 calib + guided |
|---|---|---|---|---|---|---|---|---|
| metabolomics | 3 | 0.578 (0.52 to 0.61) | 0.565 (0.52 to 0.60) | 0.565 (0.52 to 0.59) | 0.569 (0.52 to 0.60) | 0.615 (0.51 to 0.76) | 0.573 (0.50 to 0.62) | 0.573 (0.50 to 0.63) |
| proteomics_affinity_v4 | 4 | 0.483 (0.45 to 0.53) | 0.489 (0.46 to 0.52) | 0.499 (0.46 to 0.53) | 0.489 (0.44 to 0.53) | 0.661 (0.56 to 0.77) | 0.617 (0.54 to 0.70) | 0.625 (0.55 to 0.70) |
| proteomics_affinity_v5 | 6 | 0.492 (0.41 to 0.56) | 0.545 (0.48 to 0.62) | 0.554 (0.48 to 0.63) | 0.506 (0.45 to 0.54) | 0.632 (0.55 to 0.70) | 0.557 (0.51 to 0.59) | 0.611 (0.54 to 0.67) |
| proteomics_explore_v5 | 3 | 0.424 (0.35 to 0.50) | 0.483 (0.41 to 0.54) | 0.515 (0.42 to 0.64) | 0.431 (0.35 to 0.49) | 0.777 (0.71 to 0.83) | 0.643 (0.62 to 0.69) | 0.690 (0.61 to 0.74) |
| proteomics_ms_v4 | 4 | 0.459 (0.34 to 0.57) | 0.541 (0.49 to 0.59) | 0.545 (0.48 to 0.62) | 0.498 (0.44 to 0.53) | 0.716 (0.57 to 0.87) | 0.569 (0.53 to 0.60) | 0.565 (0.51 to 0.61) |
| transcriptomics_v4 | 6 | 0.461 (0.38 to 0.55) | 0.587 (0.47 to 0.68) | 0.621 (0.49 to 0.75) | 0.473 (0.41 to 0.53) | 0.687 (0.52 to 0.86) | 0.603 (0.49 to 0.75) | 0.659 (0.50 to 0.82) |

Paired differences:

| block | comparison | cohorts | mean paired difference (95% cohort bootstrap) | interval above zero |
|---|---|---|---|---|
| metabolomics | 2 − 1 | 3 | -0.012 (-0.029 to -0.001) | no |
| proteomics_affinity_v4 | 2 − 1 | 4 | +0.010 (-0.009 to +0.032) | no |
| proteomics_affinity_v5 | 2 − 1 | 6 | +0.056 (+0.008 to +0.108) | yes |
| proteomics_explore_v5 | 2 − 1 | 3 | +0.061 (-0.007 to +0.123) | no |
| proteomics_ms_v4 | 2 − 1 | 3 | +0.023 (-0.035 to +0.066) | no |
| transcriptomics_v4 | 2 − 1 | 6 | +0.127 (-0.018 to +0.292) | no |
| metabolomics | 2 − 4 | 3 | -0.005 (-0.008 to -0.001) | no |
| proteomics_affinity_v4 | 2 − 4 | 4 | -0.000 (-0.013 to +0.013) | no |
| proteomics_affinity_v5 | 2 − 4 | 6 | +0.040 (+0.009 to +0.088) | yes |
| proteomics_explore_v5 | 2 − 4 | 3 | +0.051 (-0.000 to +0.095) | no |
| proteomics_ms_v4 | 2 − 4 | 3 | +0.044 (+0.015 to +0.069) | yes |
| transcriptomics_v4 | 2 − 4 | 6 | +0.114 (-0.025 to +0.265) | no |
| metabolomics | 7 − 5 | 3 | -0.042 (-0.128 to +0.013) | no |
| proteomics_affinity_v4 | 7 − 5 | 4 | -0.037 (-0.071 to -0.004) | no |
| proteomics_affinity_v5 | 7 − 5 | 6 | -0.022 (-0.044 to -0.004) | no |
| proteomics_explore_v5 | 7 − 5 | 3 | -0.086 (-0.097 to -0.077) | no |
| proteomics_ms_v4 | 7 − 5 | 3 | -0.070 (-0.149 to -0.015) | no |
| transcriptomics_v4 | 7 − 5 | 6 | -0.028 (-0.051 to -0.005) | no |

**C. Full range: concordance tracks transfer over LC and controls together, but the LC points do not lie on the
controls' curve.**

- **Data:** the positive controls run through the same pipeline add 20 points from 7 held-out cohorts:
  - acute versus never infected, proteomics (Cervia2024_SomaScan, Su2022, Zoodsma2022_Hannover; 12 points);
  - inferred sex, transcriptomics, all genes (GSE224615, GSE226260_WB, GSE226260_PBMC, GSE270045; 8 points).

  The acute versus never-infected transcriptomic control has no point: it was never run, because only one cohort
  (GSE215865) has both classes.
- **Pooled Spearman with effect_spearman:** 0.503 for the 87 LC points; 0.641 (0.46 to 0.75) for LC and controls
  together. With sign_agreement: 0.512 and 0.613.
- **One increasing isotonic curve over all points:** R2 0.52 with effect_spearman.
  - The LC points lie below it (mean residual −0.034) and the controls above it (+0.146).
  - Separate curves for LC and controls reduce the squared error from 2.06 to 1.16. No permutation of the group labels
    did as well (p 0.0005, the smallest possible with 2,000 permutations).
  - The same holds with sign_agreement, and with the sex control on genes outside chrX and chrY (no_xy).
  - At the same agreement, the LC points reach a lower AUROC than the control points. Within the range of agreement the
    controls cover, LC points sit 0.035 to 0.064 below the joint curve on average, depending on the
    metric and sex gene set.

`scripts/concordance_vs_transfer_full.py`. Files: `points.csv`, `correlation.csv`, `monotone.csv`, `curve.csv`.

| points | metric | rows | cohorts | Spearman (cohort bootstrap) | p |
|---|---|---|---|---|---|
| LC | effect_spearman | 87 | 20 | 0.503 (0.27 to 0.66) | 0.0001 |
| LC | sign_agreement | 87 | 20 | 0.512 (0.33 to 0.64) | 0.0001 |
| controls | effect_spearman | 20 | 7 | 0.470 (-0.34 to 0.94) | 0.0143 |
| controls | sign_agreement | 20 | 7 | 0.727 (0.10 to 0.82) | 0.0001 |
| LC and controls | effect_spearman | 107 | 27 | 0.641 (0.46 to 0.75) | 0.0001 |
| LC and controls | sign_agreement | 107 | 27 | 0.613 (0.42 to 0.73) | 0.0001 |
| LC and controls no xy | effect_spearman | 107 | 27 | 0.635 (0.46 to 0.74) | 0.0001 |
| LC and controls no xy | sign_agreement | 107 | 27 | 0.624 (0.45 to 0.73) | 0.0001 |

| points | metric | R2 joint isotonic fit | mean residual LC / control | squared error joint → separate | p (separate curves) | LC points in control range (mean residual) |
|---|---|---|---|---|---|---|
| LC and controls | effect_spearman | 0.520 | -0.034 / +0.146 | 2.06 → 1.16 | 0.0005 | 51 of 87 (-0.064) |
| LC and controls | sign_agreement | 0.473 | -0.040 / +0.175 | 2.26 → 1.13 | 0.0005 | 70 of 87 (-0.043) |
| LC and controls no xy | effect_spearman | 0.512 | -0.030 / +0.130 | 1.85 → 1.14 | 0.0005 | 58 of 87 (-0.046) |
| LC and controls no xy | sign_agreement | 0.482 | -0.035 / +0.151 | 1.96 → 1.11 | 0.0005 | 77 of 87 (-0.035) |
