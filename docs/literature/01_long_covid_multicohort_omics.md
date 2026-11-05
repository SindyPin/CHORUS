# Literature stream 1: multi cohort, pooled and externally validated Long COVID omics (searched 4 Nov 2026)

## Relevant papers

Multi cohort pooling and meta-analysis
1. Maison DP, 2025, Mitochondrion 85:102072, doi 10.1016/j.mito.2025.102072. 21 GEO BioProjects, 1,272 individuals (673 acute, 444 no known COVID, 132 recovered, 23 LC). Uniform reprocessing (STAR GRCh38) then DESeq2 with study as covariate. No ML, no external validation. All 23 LC samples from one study: study and condition confounded. Most important comparator.
2. Varga B, 2026, Front Syst Biol 6:1873899, doi 10.3389/fsysb.2026.1873899. Five datasets, 394 samples: discovery GSE275334 (15 LC vs 18 HC), validation GSE270045 (19 LC/ME-CFS vs 17 HC) and GSE226260 (124 PASC vs 63 recovered), plus GSE157103 and E-ENAD-46. Mean imputation and z scoring within CV folds; no batch correction. MCC-REFS feature selection to 7 genes (PPP2CB, SOCS3, ARG1, IL6R, ECHS1, FZD2, TRGV3/5), MLP; cross dataset AUC 0.85, 0.68, 0.78, 0.77. Nearest competitor; no leave one study out, no harmonisation.
3. Lopez-Rincon A, 2025, medRxiv 10.1101/2025.02.07.25321856 (predecessor of 2). GSE275334 discovery, GSE270045 and GSE157103 test. REFS inside nested 10 fold CV; no batch correction method.
4. Eroglu E and Firat Y, 2026, Zenodo 10.5281/zenodo.21398550. Six GEO cohorts (GSE267625, GSE169687, GSE152418, GSE222253, GSE265753, GSE235938), 343 individuals. ComBat plus log2 CPM, XGBoost on 149 MitoCarta genes, nested CV AUC 0.999, external AUC 0.984. Train versus test cohorts not stated; implausible for cross cohort transfer; ComBat applied across all cohorts (leakage pattern). Cautionary citation only.
5. Avdiu B, 2025, master's thesis, FH Wiener Neustadt, github.com/besAvdiu303/LONG-COVID-001. Ten GEO series, 515 PBMC samples; kallisto, DESeq2 with study covariate, limma removeBatchEffect; GSVA, WGCNA. No ML; README states all LC samples from one study.
6. Khoo LY, 2025, Sci Rep 15:19644, doi 10.1038/s41598-025-03640-1. Ten GEO RNA-seq datasets, acute COVID versus controls; DESeq2 with batches, SVA, RankProd meta-analysis. Template not yet applied to LC.
7. Tasoula A, 2026, Front Immunol 17:1776555, doi 10.3389/fimmu.2026.1776555. Cross study multi omics synthesis (hamster and human), OXPHOS suppression; no harmonisation framework or model.
8. Maestri E and Kwon WJ, 2025, medRxiv 10.64898/2025.12.04.25341650 (STOP-PASC). Public Olink meta-analysis across nine cohorts, 590 samples: discovery six cohorts (374 samples; 52 HC, 76 infected recovered, 246 LC), validation three cohorts (216; 22 HC, 61 IR, 133 LC). 60 protein Long COVID signature, summary AUROC 0.783 discovery, 0.759 validation. No classifier; authors state a locked classifier still needs development. Key proteomics comparator.
9. Farztdinov V, Ralser M, 2024, medRxiv 10.1101/2024.03.14.24304224. Zurich PACS proteomics cohort reanalysed and tested against two Charite cohorts (SWATH-MS 200 proteins; DIA-MS 233 proteins). After age, sex and BMI balancing the complement signal did not replicate. Best motivation for leave one study out.

Independent cohort validation
10. Aid M, Barouch DH, 2025/2026, Nat Immunol 27:61-71, doi 10.1038/s41590-025-02353-x, GSE226260. Discovery MassCPR/BIDMC (n 142) and independent 2023-2024 BIDMC/RECOVER validation (n 38); PBMC RNA-seq plus SomaScan v4.1; random forest on ssGSEA with 70/30 split; signatures replicated in the independent cohort.
11. Gao Y, Cai C, Adamo S, 2025, Nat Immunol 26:692-705, doi 10.1038/s41590-025-02135-5. UK (70 LC, 70 convalescent) and Swedish (95 LC, 30 controls), Olink Explore 3072; no protein survived correction in the second cohort.
12. Gabernet G, 2025/2026, JCI 193698 (preprint 10.1101/2025.02.12.637926). IMPACC, 513 participants, 20 hospitals, six omics layers; SPEAR supervised Bayesian factor model versus MOFA; internal 80/20 split (AUROC 0.69); heme metabolism signature validated in an independent 102 participant cohort.
13. Klein J et al. (MY-LC), 2023, Nature, doi 10.1038/s41586-023-06651-y. 268 analysed; LASSO; cortisol replicated externally, galectin-1 and EBV gp42 did not; full model never externally validated.
14. Cai E, Kouznetsova VL, Tsigelny IF, 2025, Metabolites 15(12):801, doi 10.3390/metabo15120801. Metabolite lists pooled from publications; MLP test AUC 0.90; could not separate PASC from fibromyalgia. Pools at list level, not sample level.
15. Lammi V et al., 2025, Nat Genet, doi 10.1038/s41588-025-02100-w. 24 studies, 16 countries, 6,450 cases; FOXP4. Multi cohort harmonisation achieved in LC genetics.
16. Fricke F, Muller-Hilke B, 2025, Front Immunol 16:1710783, GSE286325; 44 LC, 44 convalescent, scRNA; ad hoc validation against two public healthy datasets.
17. Espin E, Tebbutt SJ, 2026, Front Cell Infect Microbiol 16:1808564. Discovery cohort (12 LC, 35 controls) with external replication in GSE224615 (27 LC, 16 recovered) and GSE275334; fold change direction checks, no harmonisation, no ML.

Single cohort ML without external validation (baselines being improved on)
18. Ten Cate V, Wild PS, 2025, Research Square 10.21203/rs.3.rs-6024310/v1. Gutenberg cohorts, 495 survivors, Olink Explore 3072; LASSO with LOOCV; overall AUC 0.62; no replication possible.
19. Liew F et al. (PHOSP-COVID), 2024, Nat Immunol, doi 10.1038/s41590-024-01778-0. 657 adults (424 LC, 233 recovered), Olink Explore 384 Inflammation; ridge regression, nested CV; internal validation only.
20. Doni Jayavelu N, Altman MC, 2026, Commun Med 6:1. IMPACC clinical and immunologic features; AUROC 0.64 to 0.66; no external cohort.
21. Meyer F, Traidl S, 2026, PLOS ONE 21(2):e0338451. Single centre n 76; random forest 80/20; "validation" is internal.
22. Garcia-Hidalgo MC, Torres G, 2026, J Transl Med 24:237. Single centre n 42; multilayer within one cohort.
23. Hou Z, 2025, Immun Inflamm Dis 13, doi 10.1002/iid3.70137. Hub gene panel (CEP55, CDCA2, MELK, DEPDC1B), ROC 0.876; no validation cohort.
24. Sayyab S et al. (CoVUm), 2025, medRxiv 10.1101/2025.02.28.25323075. EPIC methylation, 22 vs 22; ComBat for slide and array; no independent cohort.
25. Missailidis D, Annesley SJ, 2024, Front Immunol 15:1450853, PRJNA1184005; 7 vs 7 PBMC; LILRB1 and LILRB2.
26. Reviews: Pinero S et al., 2025, Crit Rev Clin Lab Sci, doi 10.1080/10408363.2025.2583083 (101 studies; gap statement). Bansal A, 2026, Front Immunol 17:1783192 (narrative review; no validation in independent cohorts; batch effects unadjusted).
27. Pinero S et al., 2025, PLoS Comput Biol 21:e1013725 (MRCONTROL; GSE215865; no external cohort).
28. EGA EGAS50000000142, 31 LC patients, paired snRNA and snATAC, unpublished.

## State of the field
No systematic harmonised multi cohort Long COVID omics resource or benchmark exists. Pooled differential expression: yes, but confounded. Cross dataset ML transfer: one small lineage without harmonisation. Cross cohort proteomic signature: once, without a locked classifier. Replication evidence is negative in two independent attempts. Large single consortium resources stop at internal CV. Nobody has published a harmonised sample level multi cohort resource with a leave one study out benchmark and early detection models.

## Introduction references
Pinero 2025 review; Maestri 2025; Maison 2025; Varga 2026; Farztdinov 2024; Gao 2025; Liew 2024; Lammi 2025; Gabernet 2025.

Caveats: Eroglu and Firat (Zenodo) unverified; Hou 2025 accessions unconfirmed; Tasoula 2026 cohort sizes not accessible.
