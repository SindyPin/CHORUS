# Literature stream 3: harmonising plasma proteomics and metabolomics across studies and platforms (searched 4 Nov 2026)

## Papers and tools
Olink PEA
1. OlinkAnalyze bridging vignette (CRAN). Median of paired NPX differences on shared bridge samples; 8 to 64 bridge samples depending on product; datasets without shared samples should not be combined by bridging.
2. OlinkAnalyze cross product bridging (Explore 3072 to HT): about 2,900 overlapping assays; bridgeable if R2 greater than 0.8 or 10 to 90 percent NPX range difference below 1 or median counts at least 150; KS statistic chooses median centring versus quantile smoothing; failing assays flagged NotBridgeable.
3. Olink Explore HT validation document 2024: 116 samples on HT and 3072, 2,839 overlapping assays, median Pearson r 0.88.
4. Rooney MR et al. 2025, Clin Chem 71(6):677, doi 10.1093/clinchem/hvaf030. SomaScan 11k versus Olink Explore HT; HT median CV 35.7 percent; between platform correlations bimodal, about one tenth of 4,443 shared proteins with r at least 0.8.
5. UK Biobank PPP normalisation note (Olink 3072 B0 to B6). Plate control, within batch intensity, across batch medians relative to batch 1; relies on randomisation.
6. Heylen D et al. 2024/2025, Brief Bioinform 26(1):bbae657, doi 10.1093/bib/bbae657. SPOC correction using synthetic plasma controls as surrogate bridge (Target 96); rank based inverse normal transformation unreliable for cross cohort pooling; ComBat and CONSTANd unsuitable for PEA.
7. SCALLOP INF1 analysis plan. Rank based inverse normal transformation per cohort, per cohort regression, METAL meta-analysis; effect sizes pooled, not NPX.
8. Olink FAQ on overlapping analytes: absolute NPX differs across panels; ratios between samples preserved.

Cross platform affinity proteomics
9. Pietzner M et al. 2021, Nat Commun 12:6822, doi 10.1038/s41467-021-27164-0. Fenland, 937 SomaScan to Olink pairs; median rho 0.38.
10. Katz DH et al. 2022, Sci Adv, doi 10.1126/sciadv.abm5164.
11. Rooney MR et al. 2022/2023, Clin Chem 69(1):68. ARIC; mean Spearman 0.46.
12. Eldjarn GH et al. 2023, Nature 622, doi 10.1038/s41586-023-06563-x. 1,514 Icelanders; median Spearman 0.33 over 1,848 shared proteins.
13. Wang B et al. 2025, Nat Commun, doi 10.1038/s41467-025-56935-2. China Kadoorie Biobank; median rho 0.29.
14. Barente A, Zhao Z, Zheng Z, Wang T, Sun BB et al. 2026, Nat Commun 18, doi 10.1038/s41467-026-78443-7 (details from coverage). 11 cohorts; 757 consistent probe pairs; disease models lose performance across platforms; optimum correlation cutoff 0.3 to 0.4.
15. cpiVAE (Bader lab), bioRxiv 2025, doi 10.1101/2025.11.20.689433. Joint VAE imputing Olink from SomaScan; median per protein r 0.609.

MS proteomics
16. Dai C, Pfeuffer J et al. 2024, Nat Methods 21:1603, doi 10.1038/s41592-024-02343-1. quantms: 83 ProteomeXchange datasets, 29,354 files; 4,993 plasma proteins pooled; no batch correction step.
17. Choi M et al. 2020, Nat Methods, MassIVE.quant (PubMed 32929271).
18. Voss H et al. 2022, Nat Commun 13:3523, doi 10.1038/s41467-022-31007-x. HarmonizR: ComBat or limma on sub matrices by missingness pattern; only 28 to 58 percent of proteins present in every batch.
19. Cuklina J et al. 2021, Mol Syst Biol 17:e10240, doi 10.15252/msb.202110240. proBatch: median centring under heavy missingness; ComBat needs complete features.
20. Chen Q et al. 2025, Nat Commun 16:9735, doi 10.1038/s41467-025-64718-y. Seven correction methods; protein level correction most robust; ratio to reference best for prediction in 1,431 plasma samples.
21. Yu Y et al. 2023, Genome Biol 24, doi 10.1186/s13059-023-03047-z. Quartet reference materials; ratio based scaling tops ranking; ComBat fails in confounded designs.
22. Harris L et al. 2023, J Proteome Res 22:3427, doi 10.1021/acs.jproteome.3c00205. Imputation benchmark.
23. Kardell O et al. 2025, J Proteome Res 24(3):1017, doi 10.1021/acs.jproteome.4c00644. Eight laboratory round robin: 300 to 400 identifications per run; 71 proteins reproducibly detected everywhere.
24. Pre analytical drivers of bias in plasma proteomics, bioRxiv 2025, doi 10.1101/2025.05.07.652495. Contamination indices; do not mix EDTA, heparin and serum.

MS plus affinity
25. Sissala N et al. 2025, Commun Chem 8:327, doi 10.1038/s42004-025-01753-2. Olink Explore 3072 versus HiRIEF LC MS; 1,129 shared; median Spearman 0.59.

Metabolomics
26. RefMet (Metabolomics Workbench; Nat Methods 2020).
27. Patt A et al. 2025, J Proteome Res 24(5), doi 10.1021/acs.jproteome.4c01051. metLinkR.
28. Villalba H et al. 2023, Metabolites 13(12):1167, doi 10.3390/metabo13121167. MMS.
29. Metanorm, bioRxiv 2025, doi 10.1101/2025.09.30.679445. QC drift normalisers.
30. Pang Z et al. 2024, MetaboAnalyst 6.0, NAR.

## Pooling Olink without bridge samples
Vendor: do not place NPX from different projects on one scale without shared samples; only about 2,900 of 3072 assays map to HT and some fail bridgeability. Consortium practice (SCALLOP): never pool NPX; transform within cohort, estimate per cohort, meta-analyse. For ML the only safely comparable quantities are within study ranks or standardised values and study level effect sizes. Within study z scoring assumes comparable disease to control composition; rank INT is unreliable across cohorts (Heylen). Practice for CHORUS: restrict to assays common to all panel generations, drop assays with high below LOD fraction, z score or rank within study, treat pooled models as rank level signatures, note that panel version is confounded with year and site.

## Pooling public MS plasma proteomics
Re-search with one FASTA, one FDR procedure and one quantification (quantms, DIA-NN, FragPipe). Usable intersection is small: 28 to 58 percent of proteins in every batch (HarmonizR), 71 reproducible proteins across eight laboratories, 300 to 400 identifications per neat plasma run. Compute contamination indices per sample. Correct at protein level by median centring or ratio to reference; use missingness aware tools, do not impute then ComBat. Expect overlap in the low hundreds and MNAR missingness; restrict to consistently quantified proteins or encode missingness.

## Gaps
No method evaluates pooling Olink across panel generations without bridges or which transformation preserves comparability. Cross platform work is Olink versus SomaScan; no joint Olink plus MS model at UniProt level. quantms harmonises identification, not batch. No Long COVID study pools proteomic cohorts at sample level; no framework jointly harmonises proteomics, metabolomics and clinical covariates for ML.
