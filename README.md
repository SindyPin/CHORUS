# CHORUS: Cross cohort Harmonisation of Omics for Robust and Unified early detection Signatures

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python 3.10+](https://img.shields.io/badge/python-3.10%2B-blue.svg)](https://www.python.org/)

CHORUS is a Python framework and benchmark for pooling independent Long COVID omics cohorts
and testing, honestly, whether pooling improves machine learning detection. It couples

1. a harmonised clinical metadata schema for public Long COVID cohorts (case definition,
   control type, time since infection and its anchor, severity, sample matrix),
2. per layer technical harmonisation operators (bulk RNA-seq, Olink, MS proteomics,
   metabolomics) with fit/transform semantics that respect leave-one-study-out evaluation,
3. feature selection that treats each cohort as an environment and keeps features whose
   association with Long COVID is invariant across cohorts, optionally combined with causal
   priors from MRCONTROL, TACO and SPLIT,
4. mosaic multi-omics integration for cohorts that measured different layers, and
5. a leave-one-study-out benchmark in which the harmoniser, the selector and the classifier
   are all fitted inside the training fold, with a nested option that chooses the
   harmonisation strategy by inner cross study validation.

The package was built for the paper *CHORUS: cross cohort harmonisation of Long COVID omics
improves and externally validates machine learning detection* (ISMB/ECCB 2027 submission).

## Why

Across 101 Long COVID omics studies, machine learning models report high within cohort
performance and are almost never externally validated (Pinero et al. 2025, Crit Rev Clin
Lab Sci). Published signatures have failed to replicate across cohorts (Farztdinov et al.
2024; Gao et al. 2025). Pooled re-analyses exist but either confound study with condition
or skip batch harmonisation. CHORUS makes the harmonisation choice part of what is
validated, and never lets a held-out cohort leak into the model that is scored on it.

## Install

```bash
pip install -e .                 # core: numpy, pandas, scipy, scikit-learn, statsmodels, pyyaml
pip install -e ".[all]"          # adds inmoose (ComBat), tabpfn, mofapy2, xgboost, plotting
```

## Quick start

```bash
chorus simulate --out sim --n-studies 6        # six simulated cohorts with study effects
chorus run sim/config.yaml --within            # LOSO benchmark plus within-study CV reference
python scripts/simulation_benchmark.py         # full grid and summary figure
```

Python:

```python
from chorus import Cohort, CohortSet
from chorus.harmonise.transcriptomics import filter_low_counts, to_continuous, select_tissue
from chorus.eval import LosoConfig, run_loso, summarise

cohorts = [Cohort.from_files("GSE275334", "bulk_rna", "data/GSE275334_counts.tsv", "data/GSE275334_meta.csv", is_counts=True, accession="GSE275334"),
           Cohort.from_files("GSE267625", "bulk_rna", "data/GSE267625_counts.tsv", "data/GSE267625_meta.csv", is_counts=True, accession="GSE267625")]
cs = select_tissue(CohortSet(cohorts), ["whole_blood"]).deduplicate_subjects()
cs = to_continuous(cs, filter_low_counts(cs))
cfg = LosoConfig(harmonisers=["identity", "center_scale", "fsqn", "rank", "combat_addon"],
                 selectors=[{"name": "variance", "k": 500}, {"name": "invariant", "k": 200}],
                 classifiers=["logreg_en", "rf", "taco", "cross_study"],
                 covariates=["age", "sex"], nested_select_harmoniser=True)
results, predictions, selected = run_loso(cs, cfg)
print(summarise(results))
```

## The evaluation rule

For every held-out study: the harmoniser is fitted on the training studies and applied to
the held-out study in addon mode (label free, using only that study's own values); the
selector uses the training studies as environments; the classifier never sees the held-out
study. The one transductive comparator (`combat_pooled`, ComBat on all studies without
labels) is flagged `leak_risk` in every results row. Within-study cross validation is
reported alongside so that the optimism gap is visible.

## Harmonisers

| name | what it does | addon safe |
|---|---|---|
| identity | pool as is | yes |
| center / center_scale | per study, per feature centring (and scaling) | yes |
| quantile_ref | sample wise quantile normalisation to a training reference | yes |
| fsqn | feature specific quantile mapping of each study to the training distribution | yes |
| rank | within sample rank transform | yes (no fit) |
| gene_pair | top scoring pair indicators selected on training labels | yes |
| combat_addon | location/scale ComBat with addon estimation for a new study | yes |
| combat_pooled | ComBat (inmoose) on all studies, transductive comparator | no, flagged |

## Selectors

`variance` (most variable genes, the TACO baseline), `differential` (meta-analytic effect),
`invariant` (random effects meta-analysis with heterogeneity screen: pooled p, I squared,
Cochran Q, sign consistency), `prior` (external causal lists, alone or intersected or unioned
with the invariant screen), `all`.

## Classifiers

`logreg_en`, `rf`, `extratrees`, `svc`, `mlp`, `xgb`, `tabpfn`, `taco` (TACO soft voting
ensemble), `cross_study` (per study learners weighted by cross study AUROC, the merging
versus ensembling comparator of Patil and Parmigiani 2018).

## Metadata schema

`chorus.schema.ClinicalSchema` defines the harmonised columns and vocabularies. Required:
`sample_id, subject_id, study_id, lc_status {LC, COVID_recovered, never_infected,
disease_comparator, unknown}, case_definition {WHO_2021, NICE_NG188, CDC_4wk, NASEM_2024,
RECOVER_PASC_score, self_report, study_specific}, sample_matrix, omics_layer`. Optional:
`days_since_infection, time_anchor, acute_severity, hospitalised, age, sex, bmi,
vaccinated_before_infection, variant_era, comorbidity_count, symptom_domains, platform,
overlaps_with_study, notes`. `chorus.schema.strict_who_subset` gives the sensitivity analysis
mask for WHO compatible cases.

## HPC recipes

`scripts/hpc/` holds Slurm templates for uniform re-quantification of GEO/SRA bulk RNA-seq
with nf-core/rnaseq (Salmon, GENCODE), for DIA-NN re-search of PRIDE plasma proteomics, and
for assembling CHORUS inputs from their outputs.

## Citation

Pinero S, et al. CHORUS: cross cohort harmonisation of Long COVID omics improves and
externally validates machine learning detection. In preparation (ISMB/ECCB 2027).

Related work by the group: MRCONTROL (Pinero et al. 2025, PLoS Comput Biol,
doi 10.1371/journal.pcbi.1013725), TACO (Pinero et al. 2026, Mach Learn Sci Technol,
doi 10.1088/2632-2153/ae62ca), SPLIT (Pinero et al. 2026, under review), and the critical
review of 101 Long COVID omics studies (Pinero et al. 2025, Crit Rev Clin Lab Sci,
doi 10.1080/10408363.2025.2583083).

## License

MIT
