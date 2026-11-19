# Literature stream 2: harmonising bulk RNA-seq across studies for classifiers that generalise (searched 4 Nov 2026)

## Papers
Batch and study effect removal
1. Johnson, Li, Rabinovic 2007, Biostatistics, doi 10.1093/biostatistics/kxj037. ComBat (sva).
2. Leek and Storey 2007, PLoS Genet, doi 10.1371/journal.pgen.0030161. SVA.
3. Zhang, Parmigiani, Johnson 2020, NAR Genom Bioinform 2(3):lqaa078, doi 10.1093/nargab/lqaa078. ComBat-seq; benchmarked for DE, not for cross study prediction.
4. Zhang X 2025, Comput Struct Biotechnol J 27:58, doi 10.1016/j.csbj.2024.12.010. ComBat-ref (reference batch); no classification benchmark.
5. Risso et al. 2014, Nat Biotechnol, doi 10.1038/nbt.2931. RUVSeq.
6. Gagnon-Bartsch and Speed 2012, Biostatistics, doi 10.1093/biostatistics/kxr034. RUV-2; RUV-4 (2013 report; R package ruv).
7. Molania et al. 2023, Nat Biotechnol 41:82, doi 10.1038/s41587-022-01440-w. RUV-III-PRPS.
8. Parker, Corrada Bravo, Leek 2014, PeerJ 2:e561, doi 10.7717/peerj.561. Frozen SVA for prediction; fails when batch and outcome are nearly confounded.
9. Hornung, Causeur, Bernau, Boulesteix 2016/2017, Bioinformatics 33:397 (verify DOI 10.1093/bioinformatics/btw650). Addon batch adjustment; cross study validation over 34 datasets; ComBat, mean centring, ratio-A and ratio-G addon all helped; frozen SVA and FAbatch often hurt; addon quantile normalisation never hurt. Package bapred.
10. Nygaard, Rodland, Hovig 2016, Biostatistics 17:29, doi 10.1093/biostatistics/kxv027. ComBat with group covariate on unbalanced designs inflates confidence.
11. Zindler et al. 2020, BMC Bioinformatics 21:271, doi 10.1186/s12859-020-03559-6. Simulating ComBat: false positives.
12. Soneson, Gerster, Delorenzi 2014, PLoS ONE 9:e100335, doi 10.1371/journal.pone.0100335. Batch confounded with class biases CV; use study wise resampling.
13. Hui, Kong, Goh 2024, Brief Bioinform 25(6):bbae515, doi 10.1093/bib/bbae515. Review: imbalance causes mis correction; new data is a new batch; external data is the true test.
14. Van, Alvarez, Mize, Gannavarapu, Chintham Reddy, Nasoz, Han 2024, BMC Bioinformatics 25:181, doi 10.1186/s12859-024-05801-x. 16 preprocessing pipelines, SVM tissue classifiers TCGA to GTEx and ICGC/GEO; limma plus QN plus scaling best on GTEx, every batch correction hurt on ICGC/GEO; no pipeline generalised to both. Closest existing work.
15. Sun, Zhang, Liu, Zheng, Zou 2024/2026, bioRxiv 10.1101/2024.09.30.615938. Ten cross platform methods; supervised methods leak labels; limma best; QN maps new samples to a fixed reference.
16. BatchQC 2.x (Bioconductor, wejlab/BatchQC). Diagnostics.
17. Harmony and MNN on bulk: no dedicated benchmark found.

Cross platform and cross study normalisation for ML transfer
18. Thompson, Tan, Greene 2016, PeerJ 4:e1621, doi 10.7717/peerj.1621. TDM.
19. Franks, Cai, Whitfield 2018, Bioinformatics 34:1868, doi 10.1093/bioinformatics/bty026. FSQN.
20. Skubleny et al. 2024, BMC Bioinformatics 25:136, doi 10.1186/s12859-024-05759-w. FSQN versus FSMVN equivalent.
21. Wu et al. 2025, Mol Biomed 6:22, doi 10.1186/s43556-025-00263-w (review of TSP, k-TSP, GERs; switchBox, multiclassPairs); Lu et al. 2024, J Transl Med, doi 10.1186/s12967-024-04940-2 (rank based trees).
22. Jin et al. 2025, Brief Bioinform 26(1):bbae661, doi 10.1093/bib/bbae661. scPAGE2 gene pair differences, 7 external bulk cohorts, mean AUROC 0.947.

Multi study integration and cross study validation
23. Rohart, Eslami, Matigian, Bougeard, Le Cao 2017, BMC Bioinformatics 18:128, doi 10.1186/s12859-017-1553-8. MINT (mixOmics).
24. Choi and Chae 2026, Sci Rep 16:12341, doi 10.1038/s41598-026-42355-9. MoDAmix adversarial domain adaptation; within dataset CV only.
25. Bernau et al. 2014, Bioinformatics 30:i105, doi 10.1093/bioinformatics/btu279. Cross study validation; within study CV overstates and reorders.
26. Riester et al. 2014, JNCI 106:dju048, doi 10.1093/jnci/dju048. Leave one dataset out on curatedOvarianData.
27. Patil and Parmigiani 2018, PNAS 115:2578, doi 10.1073/pnas.1708283115. Cross study learners beat merging as heterogeneity grows; Guan, Parmigiani, Patil arXiv 1905.07382.
28. Haynes et al. 2017, PSB 22:144. MetaIntegrator (leave one dataset out meta-analysis).
29. Korkmaz 2026, arXiv 2604.10965. bioLeak: leaky pipeline AUC 0.81 versus guarded 0.61 on 14 ovarian studies.

Domain adaptation and foundation models
30. Dradjat, Hamidi, Hanczar 2026, arXiv 2603.08062. Adversarial domain adaptation across TCGA, ARCHS4, GTEx.
31. Padron-Manrique et al. 2025, arXiv 2504.10343. DANN; 5 fold CV only.
32. Kang et al. 2025, bioRxiv 10.1101/2025.06.11.659222. BulkFormer (522,769 GEO/ARCHS4 profiles); no study held out evaluation.
33. Zhou, Agarwal, Gopinath, Kassis 2025, bioRxiv 10.1101/2025.08.15.670537. TabPFN on 58,000 RNA-seq samples: 500 feature limit; below classical methods on IBD; no cross cohort test.

Uniform reprocessing
34. Wang, Armenia, Zhang et al. 2018, Sci Data 5:180061, doi 10.1038/sdata.2018.61. After uniform STAR/RSEM, PC2 still data source; ComBat per tissue still needed. Re-quantification necessary but not sufficient.
35. Arora et al. 2020, Sci Rep 10:2734, doi 10.1038/s41598-020-59516-z. Pipeline choice confounds integration.
36. Wilks et al. 2021, Genome Biol 22:323 (recount3); Lachmann et al. 2018, Nat Commun 9:1366 (ARCHS4); refine.bio; nf-core/rnaseq (Ewels 2020).

Blood composition
37. O'Connell 2023, Sci Rep 13:15514, doi 10.1038/s41598-023-41443-4. NLR dominates whole blood expression; TMM, median ratio and QN obscure it.
38. Gautam et al. 2019, PLoS ONE 14:e0225137, doi 10.1371/journal.pone.0225137. Whole blood and PBMC separate cleanly.
39. Deconvolution (CIBERSORTx, MuSiC, BayesPrism; Brief Bioinform 2025 bbaf031 benchmark); no paper tested composition adjustment as a harmonisation step for a pooled classifier.

## Best practice and failure modes
Re-quantify from raw reads (necessary, not sufficient). Learn correction on training studies only; apply addon or frozen to the held out study. Supervised correction leaks labels. Report under leave one study out; within study CV overstates by about 0.2 AUC in bioLeak's case study. Simple location scale methods do as well as elaborate ones; QN and FSQN map new samples without refitting. Failure modes: study confounded with case ratio (ComBat with covariate removes or fabricates biology; without it removes disease signal); over correction with small studies; no universal pipeline (validate the choice out of study); composition dominates blood; a new study is a new batch (rank based representations trade power for robustness); ensembles beat merging under heterogeneity.

## Existing end to end tools and the gap
bapred (addon adjustments, no study resampling or covariates), MINT (own learner only), MetaIntegrator (signature discovery), bioLeak (split plans and audits, no ComBat), BatchQC (diagnostics). Nothing couples a menu of addon harmonisers, leave one study out validation that selects among them, and harmonised clinical covariates with composition handling.
