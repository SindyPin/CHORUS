# Literature stream 5: cohort as environment, causal and invariant learning for cross cohort transfer (searched 4 Nov 2026)

## Papers
Foundations
1. Peters, Buhlmann, Meinshausen 2016, JRSS B, doi 10.1111/rssb.12167. Invariant causal prediction. R InvariantCausalPrediction.
2. Meinshausen, Hauser, Mooij, Peters, Versteeg, Buhlmann 2016, PNAS, doi 10.1073/pnas.1510493113. ICP on yeast expression with knockouts as environments. Only omics ICP application located.
3. Rojas-Carulla, Scholkopf, Turner, Peters 2018, JMLR 19. Invariant models for causal transfer learning.
4. Arjovsky, Bottou, Gulrajani, Lopez-Paz 2019, arXiv 1907.02893. IRM.
5. Rothenhausler, Meinshausen, Buhlmann, Peters 2021, JRSS B 83:215, doi 10.1111/rssb.12398. Anchor regression.
6. Pfister, Williams, Peters, Aebersold, Buhlmann 2021, Ann Appl Stat 15:1220, doi 10.1214/21-AOAS1487. Stabilised regression (proteomics application). R StabilizedRegression.
7. Sagawa, Koh, Hashimoto, Liang 2020, ICLR, arXiv 1911.08731. Group DRO.
8. Wu, Yin, Wang, Cunningham, Blei 2025, arXiv 2506.22675. Bayesian invariance modeling (BIP), 6169 features on Kemmeren yeast.
9. Londschien, Burger, Ratsch, Buhlmann 2026, arXiv 2507.21783. Anchor regression with ICU database as anchor (nine databases, about 400k patients); 1 to 3 percent OOD gains; leave one environment out tuning under selects. Closest analogue to cohort as environment in clinical data.
10. Nastl and Hardt 2024, NeurIPS, arXiv 2402.09891. On 16 tabular tasks with natural shifts, all feature predictors beat causal subsets in and out of domain; IRM and REx no advantage. Key negative result to answer.

Omics and biomedical applications
11. Perez, Hong, Zweig, Azizi 2025, bioRxiv 10.1101/2025.09.22.677881. Domain invariant feature learning (CDANN plus CenterLoss) for patient level prediction from single cell; held out city or donor; roughly halves the ERM drop.
12. Ng, Lyu, Jiang, Li, Chen 2025, bioRxiv 10.64898/2025.11.30.691365. InvGraphCox: cohorts as environments, leave one domain out survival on multi cohort cancer transcriptomics.
13. O'Donoghue et al. 2021, arXiv 2111.07348. IRM across radiation RNA-seq experiments.
14. Yu, Liu, Li, Ding, Le 2019, IEEE TPAMI, doi 10.1109/TPAMI.2019.2908373. MCFS multi source causal feature selection (invariant Markov blanket style features across sources). Thuc Le coauthor.
15. Raab, Johnson, Piccolo 2026, bioRxiv 10.64898/2026.06.24.734352. Batch effects on ML classifiers across six TB cohorts; supervised ComBat can invert class signal; internal CV optimistic.
16. Patil and Parmigiani 2018, PNAS; Ramchandran, Patil, Parmigiani 2020, PSB; Guan, Parmigiani, Patil arXiv 1905.07382. Cross study learning.

Causal views of batch correction
17. Bridgeford et al. 2025, Imaging Neuroscience 3, doi 10.1162/imag_a_00458. Batch as exposure; positivity and ignorability; causalBatch (CRAN).
18. Wang, Chaudhari, Davatzikos 2021, MICCAI, doi 10.1007/978-3-030-87199-4_17. Flow based causal harmonisation.
19. Dong et al. 2023, Nat Methods, doi 10.1038/s41592-023-02040-5. CINEMA-OT.
20. Dong, Li, Wang 2025, IJCAI. scConCRL.
21. Gao et al. 2025, Nat Commun 16:6775, doi 10.1038/s41467-025-62008-1. CausCell.
22. scCausalVI 2025, bioRxiv 10.1101/2025.02.02.636136.
23. scDisInFact 2024, Nat Commun, doi 10.1038/s41467-024-45227-w.

Causal discovery across datasets
24. Mooij, Magliacane, Claassen 2020, JMLR 21(99). Joint causal inference.
25. Huang, Zhang et al. 2020, JMLR 21(89). CD-NOD.
26. Zhou, He, Ni 2022, arXiv 2201.12392. CHOD on TCGA breast.
27. Jablonski, Pirkl, Cevid, Buhlmann, Beerenwinkel 2022, Bioinformatics 38:1550, doi 10.1093/bioinformatics/btab847. Differential causal effects (dce).
28. Bang and Didelez 2025, arXiv 2503.21526. Tiered FCI for overlapping cohorts.

Hybrid causal selection then classifier
29. Sens et al. 2024, Genome Res 34:1276, doi 10.1101/gr.279252.124. PRiMeR.
30. CausalDRIFT 2025, medRxiv 10.1101/2025.07.10.25331298.
31. Pinero et al. 2025/2026, TACO (own work): DCE plus TabPFN, single cohort, no harmonisation.

Tabular foundation models under shift
32. Cheng, Jia, Zhou, Li, Guo 2025, arXiv 2505.16226. TabPFN v2 in open environments.
33. Loza et al. 2026, arXiv 2607.26000. TabPFN best ID and OOD AUC but shift gap equals other models.

Not verified: StableMate (Le Cao group); group DRO omics applications; arXiv 2508.12798.

## Has "study as environment, select invariant features, then harmonise" been done?
Not as a combined pipeline in omics. Pieces exist separately: anchor on database ID in ICU data; adversarial invariance for single cell patient prediction; cohorts as environments in cancer survival with mean alignment; formal ICP or stabilised regression only with experimental regimes as environments; causal batch correction stops at correction. No paper selects features by invariance of Y given X across observational cohorts, harmonises, and shows a leave one cohort out ML benefit over variance or DE selection. The nearest evidence is negative and non omics (Nastl and Hardt) and Raab et al. show the adjustment protocol dominates transfer.

## Most defensible novelty claim
First study to treat Long COVID cohorts as causal environments and to test, under strict leave one cohort out validation with harmonisation fitted inside the training fold, whether invariance or causally selected features (DCE, ICP style screen, stabilised regression) yield an early detection classifier that transfers better than variance selected features harmonised the same way, with per layer harmonisation and a tabular foundation model. Frame as a test, include all features and DE baselines, report the internal versus cross cohort gap. Do not claim novelty for causal feature selection, causal batch correction or TabPFN under shift.
