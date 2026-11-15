# Mosaic integration through programme scores (Test 3)

This file describes the method of Test 3. In Test 3, different omics layers are measured in different cohorts and are
bridged at the level of biological programmes. It does not interpret anything. The code is `chorus.integrate.programmes`
(scoring) and `scripts/mosaic_programmes.py` (the analyses). The results are in `results/mosaic_programmes/` and are
described in `results/README.md` ("Test 3").

## Gene sets

- **Source:** MSigDB release 2026.2 (human), Hallmark (`h.all.v2026.2.Hs`, 50 sets) and Reactome
  (`c2.cp.reactome.v2026.2.Hs`, 1,867 sets). Both were downloaded on 2026-11-15 from
  `https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2026.2.Hs/` into `data/ref/msigdb/` (local). The Entrez
  and symbol GMT files were downloaded; the analysis uses the Entrez files. MD5s are in `docs/provenance.md`.
- **Licence:** MSigDB v2022.1 and later is distributed under the Creative Commons Attribution 4.0 International
  licence (© Broad Institute, MIT and the Regents of the University of California). Additional terms on the MSigDB
  licence page apply only to the KEGG, BioCarta and AAAS/STKE derived sets, none of which are used here.
- **Mapping:** each Entrez id is looked up in the HGNC complete set (`data/ref/hgnc_complete_set.txt`, columns
  `entrez_id`, `ensembl_gene_id`, `uniprot_ids`). This gives the Ensembl gene id (RNA) and every UniProt accession of
  the gene (protein). An Entrez id that is not in HGNC stays in the set unmapped. `results/mosaic_programmes/genesets.csv`
  gives per programme the genes, those found in HGNC, those with an Ensembl id and those with a UniProt accession, and
  the number of accessions.
- **Proteins:** a protein feature is a member of a programme when its UniProt accession belongs to a member gene. A
  multi-protein key (`P1_P2`, Olink and SomaScan assays of complexes) is a member when any of its accessions is, and it
  counts as one member.

## Programme scores

`chorus.integrate.programmes.score_programmes`, separately in each cohort:

1. **Measured features.** A feature is measured in the cohort when it has at least one value there and a non-zero
   standard deviation over the cohort's samples.
2. **Standardisation.** Each measured feature is standardised over the cohort's samples (mean 0, SD 1). No labels are
   used.
3. **Score.** The score of a programme in a sample is the mean of the standardised values of the programme's measured
   members.
4. **Minimum members.** A programme is scored in a cohort only when at least 5 of its members are measured there.
   Otherwise its score is missing for every sample of that cohort.

The score of a cohort depends only on that cohort's own samples, so it is the same whether or not the cohort is held
out.

**Coverage** (`coverage.csv`): one row per cohort and programme, with these columns:
- the members in the layer's feature space;
- the members measured;
- the fraction measured;
- whether the programme is scored.

**Input values per layer.** These are the values after each config's preprocessing:
- **RNA:** log2 CPM. Every cohort keeps the genes that pass its own expression filter: CPM ≥ 1 in at least half of its
  samples (`preprocessing: features: per_cohort`, new). The other genes are missing in that cohort.
- **Protein:** each cohort keeps its own proteins (`features: union`), with the within-cohort imputation of the config
  (study_median for affinity platforms, study_min for MS).
- **Samples:** only the samples of the contrast are kept and scored. The standardisation is over those samples.

## Cohorts and samples

**Environment.** Each cohort is one environment, named by layer and cohort (`rna:GSE224615`, `protein:Su2022`).

**LC versus COVID_recovered:**
- **RNA** (`configs/mosaic_lc_rna.yaml`): GSE226260_WB, GSE226260_PBMC, GSE224615, GSE251849 and GSE270045. These have
  the sample selection of `configs/lc_transcriptomics_blood.yaml`: whole blood and PBMC, the earliest post-acute
  sample per subject, and subjects in two cohorts removed. Only LC and COVID_recovered samples are kept.
- **Protein** (`configs/concordance_proteomics_all.yaml`, the Test 1 sample selection): Buggert2025,
  Cervia2024_SomaScan, Su2022, Chowdhury2024 and Kovarik2023_plasma. Cervia2024_MS is the same sera as
  Cervia2024_SomaScan, so it is a sensitivity environment:
  - it is reported in the effects and concordance tables;
  - it is a test cohort only (never trained on);
  - it replaces Cervia2024_SomaScan in a second partition (`sensitivity_Cervia2024_MS`).

**Acute versus never_infected** (positive control):
- **RNA** (`configs/mosaic_acute_rna.yaml`): GSE215865, the only RNA cohort on disk with both classes. Its metadata
  has no `sample_phase`, and every library was drawn during a hospital admission. So every library of a COVID-19
  positive subject (curated LC or COVID_recovered by later outcome) is relabelled acute through the new per-cohort
  config key `acute_lc_status`. The COVID-19 negative admissions are never_infected, and lc_status unknown is
  dropped. One library per subject, the earliest.
- **Protein** (`configs/mosaic_acute_protein.yaml`): Cervia2024_SomaScan, Su2022 and Zoodsma2022_Hannover, with the
  sample selection of `configs/control_acute_proteomics.yaml`, each cohort on its own proteins.

## Analyses

**1. Programme effects** (`effects.csv`, `effects_var.csv`). In each environment, Hedges g of every scored
programme, positive minus negative class, with its variance (`chorus.features.per_study_effects`).

**2. Concordance** (`concordance_<set>.csv`, `_p.csv`, `n_shared_<set>.csv`).
- Spearman correlation of the programme effects between every pair of environments, over the programmes scored in
  both, with the two-sided p value and the number of programmes.
- `<set>`: `all`, `hallmark` or `reactome`.

**3. Invariance partition** (`partition.csv`).
- **Method:** `InvarianceAnchored.partition` of the CHORUS harmoniser applied to the programme scores, with the
  environments as the studies. The thresholds are the operator defaults.
  - **Meta-analysis:** per programme, a random-effects (DerSimonian and Laird) meta-analysis of the per-environment g.
  - **invariant:** p ≤ 0.05, I² ≤ 50, Cochran Q p ≥ 0.05, and the same sign in every environment.
  - **null:** p > 0.2, max |g| < 0.35, I² ≤ 50, and score variance at or above the 20th percentile.
  - **unstable:** max |g| ≥ 0.5 and not invariant.
  - **none:** the remaining programmes.
- **No cap:** the operator's cap on the number of null features (800) is lifted (`n_controls` 10⁹).
- **Programmes partitioned:** the operator assumes every feature is present in every environment. So only the
  programmes scored in every environment of the analysis are partitioned. The others are listed with set
  `not_partitioned` and the number of environments that score them (`n_env_scored`).
- **Analyses:** `primary` (the primary environments) and, for LC, `sensitivity_Cervia2024_MS`.

**4. Cross-layer leave one cohort out in programme space** (`loso.csv`, `loso_predictions.csv`).
- **Features:** the programmes scored in every training cohort and in the test cohort.
- **Harmonisers:** `identity` and `center_scale`, fitted on the training cohorts. The test cohort is transformed with
  its own label-free statistics (addon).
- **Classifier:** the CHORUS `logreg` (standardised, C 0.5, balanced classes, seed 0).
- **Metrics:** on the test cohort, with a stratified bootstrap 95% interval for AUROC (1000 resamples, seed 0).
- **Minimum test size:** a test cohort with fewer than 5 samples in a class would be skipped; none was.
- **Schemes:**

| contrast | scheme | training | test |
|---|---|---|---|
| acute | `protein_to_rna` | the three protein cohorts | GSE215865 |
| acute | `rna_plus_other_protein` | GSE215865 and two protein cohorts | the third protein cohort |
| lc | `rna_to_protein` | the five RNA cohorts | each protein cohort, and Cervia2024_MS (`test_role` sensitivity) |
| lc | `protein_to_rna` | the five protein cohorts (not Cervia2024_MS) | each RNA cohort |
| lc | `leave_one_cohort_out` | the other nine of the ten primary cohorts | each of the ten |

**5. Bridge check in Cervia2024** (`bridge_cervia.csv`).
- **Samples:** the SomaScan and MS layers of the 113 sera with both layers (LC and COVID_recovered).
- **Scoring:** each layer is scored on these samples (standardised within the layer).
- **Correlations:** per programme scored in both layers, Spearman and Pearson correlation across the sera of the two
  layers' scores.
- **Member sets:**
  - `own`: each layer's own measured members;
  - `shared`: only the members measured in both layers, so the scores rest on the same proteins measured on two
    platforms.
- **Counts:** the member numbers of each layer and the shared members are given per programme.

## Outputs not committed

- **Scores:** `scores.csv.gz` (samples × programmes per contrast). They are derived from the local data.
- **Log:** `results/mosaic_programmes.log`.
