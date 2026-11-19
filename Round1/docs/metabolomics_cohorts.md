# Metabolomics cohorts (and the epigenomics check)

This file records the public metabolomics cohorts with a Long COVID (LC) versus COVID_recovered contrast, how each was
built, how LC and recovered were defined, and how the metabolite names were mapped to a shared identifier. It does not
interpret anything.

- **Code:**
  - builders: `scripts/metabolomics/build_cohorts.py`;
  - identifier mapping: `scripts/metabolomics/map_ids.py`, using `chorus.metabolites` and the `metabolomics_id` adapter.
- **Analyses:** `results/README.md` ("Metabolomics").
- **Provenance and checksums:** `docs/provenance.md` ("Metabolomics cohorts").

## Sources

The candidates come from the catalogue (`Omics_Data.xlsx`, sheet Main_Sources, rows with input data type
Metabolomics).

| catalogue row | cohort | source | used |
|---|---|---|---|
| 75 | Martinez2024 | Metabolomics Workbench ST003103 | built |
| 72 | Lopez2023 | Mendeley Data 10.17632/8zfdjsypd8.1 | built |
| 73 | Holmes2021 | Zenodo 10.5281/zenodo.4713391 (concept record 4713390) | built, labels not resolved (below) |
| 70 | Saito2024 | Front Immunol supplementary Table 2 (PMC10830702) | built |
| 69 | Mikuteit2023 | Amino Acids supplementary file (PMC10564820) | not built: QC results and group summaries only, no per sample values |
| (Su2022) | Su2022_metabolomics | Mendeley Data 10.17632/96v329bg7g.1, Table S2.2 (built for Test 2) | mapped with the others |

**Other catalogue rows** (71 Taenzer 2023 urine, 74 Machado 2024 saliva, 76 Anda 2024, 77 Mendes 2024 saliva,
78 Ansone 2024, 79 Sorokina 2023) give "Under request" as their data source and were not pursued.

## How each cohort was built and labelled

Values are log2 of the deposited values; zero and negative values are missing. Su2022 is kept as deposited (the
authors' adjusted Z scores).

**Martinez2024**
- **Source:** Martinez S et al. 2024, J Proteome Res, PMID 38566450.
- **Samples:** plasma, 40 LC and 40 recovered.
- **Platforms:** untargeted CE-ESI(+/-)-TOF-MS and GC-Q-MS, plus LC-ESI(+/-)-QTOF-MS lipidomics. The five Workbench
  analyses AN005077 to AN005081 are joined, with names prefixed `LCMSPOS::`, `LCMSNEG::`, `GCMS::`, `CEMSPOS::` and
  `CEMSNEG::`. The 12 QC samples are dropped.
- **Missing data:**
  - The CE features are named by m/z and migration time (unannotated).
  - One sample has no GC-MS values.
- **Labels:** the Workbench factor of each sample is "Long COVID" or "Recovered from acute COVID-19". The abstract
  compares "Long COVID patients" with "individuals fully recovered from acute mild COVID-19". The paper is behind a
  subscription, so no further case definition was read.
- **Timing:** the catalogue gives 12+ weeks after infection.

**Lopez2023**
- **Source:** López-Hernández Y et al. 2023, Sci Rep, PMID 37528111, PMC10394026.
- **Samples:** plasma, Mexico.
- **Platform:** targeted TMIC PRIME LC-MS/MS and FIA-MS/MS, µM.
- **Definitions** (paper, Methods):
  - "Recovered patients were classified as those who did not report persistent symptoms. Long COVID was considered if
    patients reported at least one persistent neurologic, psychiatric, gastrointestinal, cardiologic, respiratory, or
    systemic symptom."
  - Class A is fewer than five persistent symptoms; class B is five or more. Sampling was about two years after
    infection.
  - Negative controls are "stored plasma samples from 37 individuals who tested negative for SARS-CoV-2 in 2020".
- **Built:**
  - LC = POST-COVID class A or class B; COVID_recovered = POST-COVID RECOVERED.
  - never_infected (contemporary) = CONTROL.
  - The COVID-19 group is the acute sample of 22 of the post-COVID patients: Folio `<n>.1` pairs with post-COVID Folio
    `<n>`, so these subjects share a subject id. These samples are `sample_phase` acute and carry the later outcome.
- **Counts:** the deposit has 18 class A, 12 class B, 18 recovered and 38 controls. The paper reports 17 class A,
  13 class B and 37 controls.
- **Sex:** coded 1/2 in the deposit without a key. It is kept as `sex_code_deposited`, and `sex` is not filled.
- **Metabolites:** 107 columns, from Creatinine to Octadecenoylcarnitine. Ratios, IL-17, leptin and clinical
  chemistry are not taken.

**Saito2024**
- **Source:** Saito S et al. 2024, Front Immunol, PMID 38304426.
- **Samples:** plasma, Edmonton, Canada.
- **Platform:** dansyl chemical isotope labelling LC-MS. The values are the peak ratio of each sample to a
  13C-labelled pool.
- **Features:** the peak pairs with a Tier 1 (171) or Tier 2 (299) identification, named by the identified compound.
  The authors' external identifier (HMDB or KEGG) is kept in `Saito2024_deposited_ids.csv`, and an HMDB id there is
  used when RefMet gives none.
- **Definitions** (paper): LC are Long COVID patients meeting ME/CFS criteria about 12 months after infection.
  Recovered (R) are previously infected people who recovered without symptoms. HC are healthy controls (negative
  SARS-CoV-2 immunoassays, mostly recruited before the pandemic). A are ICU patients with acute COVID-19.
- **Sample groups:**
  - The data file names its groups A, C, HC and L (15, 15, 15 and 30 samples).
  - The clinical table (Supplementary Table 1) names them A, R, HC and LC with the same sizes.
  - **C is taken as R**, by elimination: the letter C is not explained in the paper.
  - The sample numbers differ between the two tables (L1 to L30 against LC1 to LC32 with gaps), so age, sex and
    sampling dates are not linked per sample.
  - The acute samples are `lc_status` unknown (outcome not given), `sample_phase` acute.

**Holmes2021**
- **Source:** Holmes E et al. 2021, J Proteome Res, PMID 34009992.
- **Samples:** plasma of the 27 non-hospitalised follow-up patients, about 3 months after infection.
- **Not deposited:** the paper's 41 controls and 18 hospitalised patients, and its NMR data.
- **Platform:** targeted LC-MS:
  - amino acids, with `Quantity` × 2 for the dilution step, as the README says;
  - tryptophan pathway, `conc`. The per-sample file for COV01172 is not read: the README says its values must be
    divided by the molecular weight, which the file does not give, so its tryptophan values are missing.
- **Labels: not resolved.**
  - The sample list gives four symptom composites (`composite`, `compositeBinary`, `compositeRec`,
    `compositeBinaryRec`) and per-symptom codes (11, 12, 21, 22, ...), but the deposit does not explain the coding.
  - The two binary composites disagree: 23 versus 18 patients with "at least 1 sympt". The abstract reports 57% with
    persistent symptoms.
  - So `lc_status` is unknown and both binary composites are kept as columns. Holmes2021 is mapped but not used in
    the LC analyses.

**Su2022_metabolomics**
- **Source:** Test 2 (`docs/paired_cohorts.md`).
- **Values:** untargeted plasma metabolomics, the authors' batch-corrected and age-, sex- and BMI-adjusted Z scores.
- **Labels:** LC and COVID_recovered at the post-acute draw T3, by the INCOV PASC definition used for the Olink layer
  (`docs/proteomics_cohorts.md`).

## Shared identifiers

`chorus.metabolites.metabolite_ids` maps each deposited name in four steps:

1. **Prefix.** The platform prefix `<platform>::` is removed.
2. **RefMet match.** The name is standardised by the RefMet match service. Names of unannotated features (`<m/z>_<RT>`,
   `Unknown ...`) are not sent.
3. **RefMet table.** The RefMet name is looked up in the RefMet table, which gives the HMDB id, InChIKey and RefMet id.
4. **Identifier.** The identifier is the first of:
   - `HMDB:<id>` from RefMet;
   - an HMDB id deposited by the authors (Saito2024);
   - `RefMet:<name>`;
   - a deposited InChIKey.

**Duplicates.** When several names of one cohort map to one identifier (the same metabolite on two platforms of
Martinez2024, or two spellings), the name with the most observed values is kept. Values are not averaged.

<!-- generated:start -->
### Cohorts

Samples per class after the builders (post-acute samples for LC and COVID_recovered; acute samples are counted under acute whatever their later outcome).

| cohort | samples (subjects) | LC | COVID_recovered | never_infected | acute | label unknown | names |
|---|---|---|---|---|---|---|---|
| Martinez2024 | 80 (80) | 40 | 40 | 0 | 0 | 0 | 820 |
| Lopez2023 | 108 (86) | 30 | 18 | 38 | 22 | 0 | 107 |
| Saito2024 | 75 (75) | 30 | 15 | 15 | 15 | 0 | 470 |
| Holmes2021 | 27 (27) | 0 | 0 | 0 | 0 | 27 | 51 |
| Su2022_metabolomics | 683 (377) | 68 | 45 | 178 | 380 | 12 | 891 |

### Mapping per cohort

Names: metabolite columns as deposited. unannotated: m/z and RT features, not sent to RefMet. matched_refmet: names standardised by RefMet. hmdb_refmet, hmdb_deposited, refmet, inchikey_deposited: names mapped by each rule. identifiers: after duplicate names of one identifier are dropped.

| study_id | n_names | n_unannotated | n_matched_refmet | n_with_hmdb | n_hmdb_refmet | n_hmdb_deposited | n_refmet | n_inchikey_deposited | n_unmatched | n_names_mapped | n_identifiers | n_names_dropped_as_duplicate |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Martinez2024 | 820 | 223 | 584 | 146 | 146 | 0 | 438 | 0 | 13 | 584 | 511 | 73 |
| Lopez2023 | 107 | 0 | 107 | 72 | 72 | 0 | 35 | 0 | 0 | 107 | 106 | 1 |
| Saito2024 | 470 | 0 | 312 | 288 | 286 | 2 | 26 | 0 | 156 | 314 | 312 | 2 |
| Holmes2021 | 51 | 0 | 50 | 48 | 48 | 0 | 2 | 0 | 1 | 50 | 48 | 2 |
| Su2022_metabolomics | 891 | 109 | 691 | 489 | 489 | 0 | 202 | 0 | 91 | 691 | 679 | 12 |

### Identifiers shared between cohorts

Diagonal: identifiers in the cohort.

| | Martinez2024 | Lopez2023 | Saito2024 | Holmes2021 | Su2022_metabolomics |
|---|---|---|---|---|---|
| Martinez2024 | 511 | 41 | 17 | 20 | 158 |
| Lopez2023 | 41 | 106 | 28 | 27 | 61 |
| Saito2024 | 17 | 28 | 312 | 32 | 66 |
| Holmes2021 | 20 | 27 | 32 | 48 | 35 |
| Su2022_metabolomics | 158 | 61 | 66 | 35 | 679 |

### Identifiers shared by several cohorts

| what | cohorts | n |
|---|---|---|
| identifiers in at least 1 cohorts | all | 1323 |
| identifiers in at least 2 cohorts | all | 236 |
| identifiers in at least 3 cohorts | all | 56 |
| identifiers in at least 4 cohorts | all | 27 |
| identifiers in at least 5 cohorts | all | 14 |
| identifiers shared by all 3 cohorts | Martinez2024;Lopez2023;Saito2024 | 15 |
| identifiers shared by all 3 cohorts | Martinez2024;Lopez2023;Su2022_metabolomics | 38 |
| identifiers shared by all 3 cohorts | Martinez2024;Saito2024;Su2022_metabolomics | 17 |
| identifiers shared by all 3 cohorts | Lopez2023;Saito2024;Su2022_metabolomics | 25 |
| identifiers shared by all 4 cohorts | Martinez2024;Lopez2023;Saito2024;Su2022_metabolomics | 15 |
<!-- generated:end -->

## Epigenomics (concordance only): not computable

**Series.** The three series of the request were downloaded into `data/epigenomics/<series>/` (local): the processed
matrix (beta values with detection p values) and the series matrix of each. GEO gives no checksums; the SHA-256 are in
`docs/provenance.md`. All three are Illumina HumanMethylationEPIC (GPL21145) on whole blood. The groups are in
`results/concordance_epigenomics/cohorts.csv`.

| series | paper | samples | groups in GEO | LC / COVID_recovered per sample |
|---|---|---|---|---|
| GSE174818 | Balnis 2021 Clin Epigenetics (PMID 34034806) | 128 | covid status 1 (102, hospitalised COVID-19, during infection), 0 (26, non-COVID-19 patients); death 0/1 | none: drawn during the acute illness |
| GSE197152 | Balnis 2022 Clin Epigenetics (PMID 35871090) | 15 | post COVID-19, one year after infection; the individual ids are those of GSE174818 acute samples | no recovered group: the paper says the 15 patients seen one year after discharge expressed multiple PASC symptoms |
| GSE247869 | Calzari 2024 Clin Epigenetics (PMID 39164752) | 96 | Post_Covid_6Months; sex and age only | not public: the paper's Additional file 1 reports a PASC subgroup of 28 (at least one persistent symptom), but no per-sample status appears in GEO, the matrix columns or the nine supplementary files |

**Not computed.** The per-CpG LC effect needs LC and COVID_recovered samples in at least two series. None of the three
has both with per-sample labels, so no effect and no correlation were computed. GSE247869 would qualify if the per-sample
PASC status of its 96 patients were obtained from the authors (Calzari L, Gentilini D, Parati G).
