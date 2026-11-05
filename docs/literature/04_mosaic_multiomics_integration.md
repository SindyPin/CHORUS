# Literature stream 4: mosaic multi omics integration across cohorts with partially overlapping layers (searched 4 Nov 2026)

## Papers and tools
Factor models tolerating missing views and groups
1. Argelaguet R et al. 2020, Genome Biol 21:111, doi 10.1186/s13059-020-02015-1. MOFA+: groups and views, missing samples per view allowed; no native query projection. MOFA2.
2. Velten B et al. 2022, Nat Methods. MEFISTO (temporal and spatial covariates across groups).
3. Hirst DP, Terezol M, Cantini L, Villoutreix P, Vignes M, Baudot A 2025, Genome Biol 26, doi 10.1186/s13059-025-03675-7. MOTL: transfer MOFA factors from a large learning set to a small target.
4. Rohart F et al. 2017, BMC Bioinformatics 18:128. MINT (single layer, multi study); DIABLO requires matched samples.
5. nipalsMCIA (Bioconductor; Bioinformatics 41(1) btaf015, 2025).
6. Rappoport N, Shamir R 2019, Bioinformatics. NEMO: partial multi omic clustering without imputation.
7. Baena Miret S, Reverter F, Vegas E 2024, PLOS ONE 19(7), doi 10.1371/journal.pone.0307482. bwm: supervised model for block wise missing omics.
8. MosaiClusteR (CRAN), unverified.
9. JIVE, intNMF, SNF, iCluster: matched samples in standard forms.

Mosaic and diagonal single cell methods usable at sample level
10. Ghazanfar S, Guibentif C, Marioni JC 2024, Nat Biotechnol 42:284, doi 10.1038/s41587-023-01766-z. StabMap: mosaic data topology, shortest path projection, supervised reference reduction, query mapping. Most direct conceptual match.
11. He Z et al. 2024, Nat Biotechnol, doi 10.1038/s41587-023-02040-y. MIDAS product of experts VAE.
12. Zhang Z et al. 2023, Nat Commun 14, doi 10.1038/s41467-023-36066-2. scMoMaT.
13. Kriebel AR, Welch JD 2022, Nat Commun 13:780, doi 10.1038/s41467-022-28431-4. UINMF.
14. Litinetskaya A, Lotfollahi M et al., Multigrate, bioRxiv 10.1101/2022.03.16.484643.
15. Arriola M et al. 2024, arXiv 2405.11280. SC5: cohort as domain, omic as modality; imputes a modality absent from a whole cohort. Only method found that explicitly treats cohort as domain.
16. SCMMIB benchmark, bioRxiv 2025, doi 10.1101/2025.04.01.646578. 40 algorithms; MIDAS, MultiVI, GLUE, StabMap leaders by setting.
17. Baiao AR et al. 2025, Brief Bioinform 26(4) bbaf355, doi 10.1093/bib/bbaf355. Review of missing modality handling.

Deep learning with missing modalities at patient level
18. Wang T et al. 2021, Nat Commun 12, doi 10.1038/s41467-021-23774-w. MOGONET (complete modalities).
19. Tabakhi S, Chen C, Lu H 2025, arXiv 2506.22901. MAGNET.
20. Park S et al. 2025, ICML workshop. MOIRA.
21. Nambiar A, Melendez C, Noble WS 2026, bioRxiv 10.64898/2026.02.04.703630. MIMIR.
22. MOGEDN (PMC12753307).
23. Yun S et al. 2024, NeurIPS, arXiv 2410.08245. Flex MoE.
24. Lepe Soltero D et al. 2025, arXiv 2503.18856. MODIS.
25. Boussena M et al. 2026, arXiv 2605.25050. MSB: late fusion survival stacking with blockwise missing values, missingness indicators, meta learner; PIONeeR 17 centres and 21 TCGA cohorts.
26. Mallick H et al. 2024, Stat Med 43:983, doi 10.1002/sim.9953. IntegratedLearner (Bayesian stacking).
27. Choi JM, Chae H 2026, Sci Rep 16, doi 10.1038/s41598-026-42355-9. MoDAmix.
28. Zhou H et al. 2024, arXiv 2406.19611. Survey of missing modality methods.

Transfer to single omics cohorts and cross omics translation
29. Liu L et al. 2026, Nat Biomed Eng 10:1672, doi 10.1038/s41551-025-01528-z. scTranslator.
30. Schluter J, Schonhuth A 2025, arXiv 2609.24004. ReTrans (mRNA to protein, about half of proteins reach PCC 0.5).
31. Nandy S, Ma Z, arXiv 2407.19030. OrchAMP.

Federated harmonisation
32. Matschinske J et al. 2023, JMIR 25:e42621. FeatureCloud.
33. Burankova Y et al. 2024/2026, arXiv 2412.05894. fedRBE (federated limma removeBatchEffect).
34. Silva S et al. 2026, arXiv 2601.14314. Fed ComBat.
35. Hoang B et al. 2024, KDD, arXiv 2405.15081. Distributed Cluster ComBat.
36. Swarm Learning (Warnat Herresthal 2021 Nature).

Cross cohort benchmarks
37. Wissel D et al. 2025, Brief Bioinform 26(5) bbaf521, doi 10.1093/bib/bbaf521. SurvBoard.
38. Liang Y et al. 2026, arXiv 2604.05478. Immunotherapy predictors near chance on unseen cohorts.
39. Scicluna BP et al. 2025, Nat Med, doi 10.1038/s41591-025-03964-5. Consensus sepsis subtypes with ComBat co normalisation validated in two external cohorts.
40. Gabernet G et al. 2025, JCI, doi 10.1172/JCI193698. IMPACC multi omics recovery factor.

## Best fit for 8 to 15 small blood cohorts
MOFA+ with cohorts as groups is the most mature; no native query projection; factors dominated by the layer each cohort has. StabMap matches the topology (cohorts as nodes, shared layers as edges, bilayer cohorts as bridges) and supports query projection. SC5 is the only cohort as domain method but single cell only. MIDAS and Multigrate are tested far above 20 to 400 samples per cohort. For prediction, late fusion with missingness aware stacking (MSB style, or MINT style per layer models) is safest. Federated tools only if pooling is impossible.

## Gaps
No framework couples per layer cross cohort harmonisation, mosaic integration and leave one cohort out predictive validation with clinical covariates. Every mosaic method was validated on single cell data or single cohort TCGA splits. SurvBoard uses TCGA cancer types as cohorts without harmonisation. Include StabMap, MOFA+ groups, SC5 and MSB as baselines.
