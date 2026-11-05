# Simulation benchmark (scripts/simulation_benchmark.py, seed 0, 4 November 2026)

Six simulated cohorts (40 to 200 samples each, 2,000 features, 30 invariant causal genes, 60 cohort specific shortcut genes, per feature location and scale study effects, case fractions 0.2 to 0.8). Grid: 6 harmonisers plus nested selection, 3 selectors, 3 classifiers, leave one study out, plus within study five fold CV.

Mean held out AUROC (L2 logistic regression):

| harmoniser | variance (300) | differential (60) | invariant (60) |
|---|---|---|---|
| identity | 0.670 | 0.823 | 0.876 |
| center_scale | 0.594 | 0.824 | 0.882 |
| quantile_ref | 0.668 | 0.823 | 0.855 |
| fsqn | 0.536 | 0.812 | 0.897 |
| rank | 0.573 | 0.785 | 0.858 |
| combat_addon | 0.549 | 0.827 | 0.877 |
| nested (chosen by inner LOSO) | 0.636 | 0.808 | 0.888 |

Within study five fold CV with the same selectors: variance 0.94, differential 1.00, invariant 1.00. The optimism gap is therefore 0.1 to 0.4 AUROC depending on the selector.

Mean Brier score (lower is better), invariant selector: identity 0.174, center_scale 0.164, fsqn 0.149, combat_addon 0.214, rank 0.180. Harmonisation improves calibration more than it improves AUROC, because AUROC within a held out cohort is invariant to an additive shift of all scores, whereas probabilities and thresholded metrics are not. This is the reason the package reports Brier, expected calibration error and balanced accuracy alongside AUROC.

Reading. (1) Feature selection with cohort as environment dominates: the invariant screen selected 15 of 20 causal genes and no shortcut genes in the unit test, differential expression selected 12 causal and 8 shortcut genes, and variance selection (the TACO baseline) selected none of the causal genes. (2) Harmonisation helps on top of invariant selection by a small AUROC margin (fsqn 0.897 versus identity 0.876) and a larger Brier margin. (3) Nested selection of the harmoniser recovers most of the gain of the best fixed harmoniser without using the held out cohort (0.888 versus 0.897). (4) Transductive ComBat was not run here; on real data it is included as a flagged comparator.

Caveat. The simulation has only location and scale study effects and linear signal; real cohorts add composition, platform specific features and label heterogeneity, so these numbers calibrate the pipeline rather than predict the Long COVID result.
