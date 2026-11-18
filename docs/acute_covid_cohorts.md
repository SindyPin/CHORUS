# Acute COVID cohorts (Part B of the acute COVID validation)

Built 18 November 2026 by `scripts/acute_covid_sweep.py` (search, screen, fetch, build, verify, report). The rule and bars of
the validation are in `results/README.md` ("Acute COVID validation: rule and bars"). No benchmark or noise ceiling analysis
has been run on these cohorts. Data files are not committed; the label mapping tables (`docs/acute_covid_labels/`), this
file, `docs/acute_covid_candidates.csv` and the configs are.

**Candidates:** 3095 (13 included, 2825 excluded, 257 needing a decision).
Sources: GEO (E-utilities), ArrayExpress and BioStudies, SRA BioProjects, PRIDE and ProteomeXchange, Metabolomics
Workbench, Zenodo, Mendeley Data (DataCite) and the accessions in the full text of open access papers (Europe PMC).

**Severity mapping used for every cohort.** severe or critical: intensive care, mechanical ventilation, WHO progression
scale 6 or more, or death; mild or moderate: a class the study defines below that line (ambulatory, hospitalised without
oxygen, or explicitly not critical / not in intensive care). A class that straddles the line (for example "severe" meaning
supplemental oxygen of any type, or a merged "moderate/severe" class) is not assigned. Each cohort's table in
`docs/acute_covid_labels/` gives the original text and the class.

**Included cohorts.**

| cohort | layer | tissue | infection (acute / never infected) | severity (severe or critical / mild or moderate) | status |
|---|---|---|---|---|---|
| GSE198449 | RNA | whole_blood | 255 / 60 | 0 / 0 | verified with flags |
| COMBAT2022_WB | RNA | whole_blood | 77 / 10 | 19 / 30 | verified with flags |
| GSE161731 | RNA | whole_blood | 46 / 15 | 0 / 0 | verified with flags |
| GSE152641 | RNA | whole_blood | 62 / 24 | 0 / 0 | verified with flags |
| GSE300129 | RNA | whole_blood | 57 / 55 | 57 / 0 | failed |
| GSE243219 | RNA | PBMC | 35 / 15 | 0 / 0 | verified with flags |
| GSE172114 | RNA | whole_blood | 69 / 0 | 46 / 23 | verified with flags |
| GSE171110 | RNA | whole_blood | 44 / 10 | 44 / 0 | verified with flags |
| GSE213313 | RNA | whole_blood | 34 / 11 | 19 / 15 | verified with flags |
| GSE152418 | RNA | PBMC | 16 / 17 | 4 / 4 | verified with flags |
| Filbin2021_MGH | proteins | plasma | 305 / 0 | 109 / 63 | verified with flags |
| Gisby2021_plasma | proteins | plasma | 55 / 50 | 13 / 28 | verified with flags |
| GSE207015 | proteins | plasma | 124 / 70 | 56 / 44 | verified with flags |

Contrasts reach at least 10 participants per class where listed in the cards. Cross cohort duplicates (correlation above
0.99, same layer): none.

Configs (enabled: false): `configs/acute_covid_rna_infection.yaml`, `configs/acute_covid_rna_severity.yaml`, `configs/acute_covid_proteomics_infection.yaml`, `configs/acute_covid_proteomics_severity.yaml`.

## Included cohorts

### GSE198449

- **Source:** GSE198449 (RNA, whole blood); files: GSE198449_featureCounts.txt.gz
- **Paper:** PMID 36323307; DOI 10.1016/j.cels.2022.10.005 (CHARM, Cell Systems 2022)
- **Platform:** Illumina NovaSeq 6000, PAXgene whole blood RNA-seq; counts; 58884 features x 1858 samples
- **Participants (one sample each):** infection 255 acute COVID against 60 never infected; severity 0 severe or critical against 0 mild or moderate; 0 other disease (own class); 160 excluded
- **Contrasts with at least 10 per class:** infection
- **Status:** verified with flags (b counts: n_samples 1858, n_participants 475, infection_acute_covid 255, infection_never_infected 60, excluded 160; c value type: counts; integer fraction 0.882; range 0 to 1.09e+07; zeros 0.607; missing 0.000; mapping 1.0 (58929 of 58929; 45 targets with duplicates, summed); library size median 1.75e+07 (range 2.14e+06 to 1.19e+08); f outliers: 583 samples more than 3 MAD from the cohort median)
- **Label mapping:** `docs/acute_covid_labels/GSE198449.csv` (34 rows)

| check | status | result | details |
|---|---|---|---|
| a integrity | pass | matrix 58884 features x 1858 samples; file had 1858 sample columns; labelled samples in the record 1858; in both 1858 | empty columns 0; duplicated columns 0; labelled not in matrix []; matrix not labelled [] |
| b counts | flag | n_samples 1858, n_participants 475, infection_acute_covid 255, infection_never_infected 60, excluded 160 | PMID 36323307 abstract (whole CHARM cohort, not the RNA-seq subset): 2,885 recruits, 1,170 infected; no comparable number |
| c value type | flag | counts; integer fraction 0.882; range 0 to 1.09e+07; zeros 0.607; missing 0.000; mapping 1.0 (58929 of 58929; 45 targets with duplicates, summed); library size median 1.75e+07 (range 2.14e+06 to 1.19e+08) | featureCounts |
| d sex | pass | agreement 1.000 on 1852 samples (chance 0.672) | inferred {'male': 1469, 'female': 383, 'ambiguous': 6} (gmm2:XIST+RPS4Y1); disagreeing [] |
| e duplicates (within) | pass | 0 pairs of different participants with correlation above 0.99 | []; pairs of the same participant above 0.99: 0; correlation on the 2000 most variable complete features |
| f outliers | flag | 583 samples more than 3 MAD from the cohort median | library size: 292 ['GSM5947048', 'GSM5947070', 'GSM5947203', 'GSM5947233', 'GSM5947307', 'GSM5947440', 'GSM5947607', 'GSM5947615', 'GSM5947622', 'GSM5947821']; detected features: 207 ['GSM5947166', 'GSM5947179', 'GSM5947184', 'GSM5947186', 'GSM5947198', 'GSM5947202', 'GSM5947203', 'GSM5947207', 'GSM5947210', 'GSM5947219']; mean correlation: 235 ['GSM5947077', 'GSM5947079', 'GSM5947131', 'GSM594714 |
| g label direction | pass | 0 of 5 markers lower in acute COVID | IFI27 g=2.16; IFI44L g=1.65; IFI6 g=1.34; ISG15 g=1.33; RSAD2 g=1.59 |
| h signal (infection) | report | cross validated AUROC 0.971; split half reliability 0.782 | 255 against 60 participants; logistic regression on the 500 most variable features, 5 x 5-fold |
| i confounding | pass | largest Cramer's V 0.68 | infection x plate: V=0.68; infection x enrollment_batch: V=0.03; infection x collection_date: V=0.56 |
| j timing | pass | 255 acute participants; days known for 255; median 0 (IQR 0 to 0, range 0 to 0) | 0 outside 0 to 30 days or marked outside the window; window unknown for 0; anchor ['positive_test'] |
| k demographics | report | acute_covid: n 255, age median 19 (missing 0), sex {'male': 215, 'female': 40}; excluded: n 160, age median 20 (missing 0), sex {'male': 116, 'female': 44}; never_infected: n 60, age median 19 (missing 0), sex {'male': 43, 'female': 17} | age band columns: [] |

### COMBAT2022_WB

- **Source:** zenodo.6120249 (RNA, whole blood); files: CBD-KEY-RNASEQ-WB.tar.gz, CBD-KEY-CLINVAR.tar.gz
- **Paper:** PMID 35216673; DOI 10.1016/j.cell.2022.01.012 (COMBAT consortium)
- **Platform:** whole blood RNA-seq (COMBAT); counts; 60683 features x 143 samples
- **Participants (one sample each):** infection 77 acute COVID against 10 never infected; severity 19 severe or critical against 30 mild or moderate; 36 other disease (own class); 0 excluded
- **Contrasts with at least 10 per class:** infection, severity
- **Status:** verified with flags (a integrity: matrix 60683 features x 143 samples; file had 143 sample columns; labelled samples in the record 144; in both 143; f outliers: 29 samples more than 3 MAD from the cohort median)
- **Label mapping:** `docs/acute_covid_labels/COMBAT2022_WB.csv` (9 rows)

| check | status | result | details |
|---|---|---|---|
| a integrity | flag | matrix 60683 features x 143 samples; file had 143 sample columns; labelled samples in the record 144; in both 143 | empty columns 0; duplicated columns 0; labelled not in matrix ['N00049-Ja005T-TRGa']; matrix not labelled [] |
| b counts | pass | n_samples 143, n_participants 123, infection_acute_covid 77, infection_never_infected 10, other_disease 36, severity_severe_critical 19, severity_mild_moderate 30 | record README (CBD-KEY-RNASEQ-WB): 143 samples from 123 individuals; differences (built, paper): {} |
| c value type | pass | counts; integer fraction 1.000; range 0 to 1.98e+06; zeros 0.541; missing 0.000; mapping 1.0 (60683 of 60683; 0 targets with duplicates, summed); library size median 3.96e+07 (range 2.3e+07 to 5.75e+07) | columns are RNASeq_sample_ID |
| d sex | pass | agreement 1.000 on 143 samples (chance 0.504) | inferred {'male': 78, 'female': 65} (gmm2:XIST+RPS4Y1); disagreeing [] |
| e duplicates (within) | pass | 0 pairs of different participants with correlation above 0.99 | []; pairs of the same participant above 0.99: 1; correlation on the 2000 most variable complete features |
| f outliers | flag | 29 samples more than 3 MAD from the cohort median | library size: 8 ['G05064-Ja005T-TRGa', 'G05112-Ja005T-TRGa', 'G05145-Ja005T-TRGa', 'N00024-Ja001T-TRGa', 'N00047-Ja005T-TRGa', 'S00028-Ja001T-TRGa', 'S00033-Ja001T-TRGa', 'S00109-Ja003T-TRGa']; detected features: 18 ['N00004-Ja001T-TRGa', 'N00006-Ja001T-TRGa', 'N00008-Ja001T-TRGa', 'N00008-Ja003T-TRGa', 'N00009-Ja001T-TRGa', 'N00010-Ja001T-TRGa', 'N00012-Ja003T-TRGa', 'N00014-Ja001T-TRGa', 'N00023 |
| g label direction | pass | 0 of 5 markers lower in acute COVID | IFI27 g=1.74; IFI44L g=1.24; IFI6 g=1.14; ISG15 g=0.89; RSAD2 g=1.13 |
| h signal (infection) | report | cross validated AUROC 0.994; split half reliability 0.781 | 77 against 10 participants; logistic regression on the 500 most variable features, 5 x 5-fold |
| h signal (severity) | report | cross validated AUROC 0.976; split half reliability 0.899 | 19 against 30 participants; logistic regression on the 500 most variable features, 5 x 5-fold |
| i confounding | pass | largest Cramer's V 0.00 | infection x institute: one level; severity x institute: one level |
| j timing | pass | 77 acute participants; days known for 0 | 0 outside 0 to 30 days or marked outside the window; window unknown for 0; anchor ['admission'] |
| k demographics | report | acute_covid: n 77, age missing, sex {'male': 43, 'female': 34}; never_infected: n 10, age missing, sex {'male': 6, 'female': 4}; other_disease: n 36, age missing, sex {'male': 21, 'female': 15} | age band columns: ['age_band'] |

### GSE161731

- **Source:** GSE161731 (RNA, whole blood); files: GSE161731_counts.csv.gz, GSE161731_counts_key.csv.gz, GSE161731_key.csv.gz
- **Paper:** PMID 33597532; DOI 10.1038/s41467-021-21289-y (McClain et al.)
- **Platform:** Illumina NovaSeq 6000 whole blood RNA-seq; counts; 60675 features x 195 samples
- **Participants (one sample each):** infection 46 acute COVID against 15 never infected; severity 0 severe or critical against 0 mild or moderate; 91 other disease (own class); 0 excluded
- **Contrasts with at least 10 per class:** infection
- **Status:** verified with flags (a integrity: matrix 60675 features x 195 samples; file had 195 sample columns; labelled samples in the record 198; in both 195; e duplicates (within): 2 pairs of different participants with correlation above 0.99; f outliers: 46 samples more than 3 MAD from the cohort median)
- **Label mapping:** `docs/acute_covid_labels/GSE161731.csv` (7 rows)

| check | status | result | details |
|---|---|---|---|
| a integrity | flag | matrix 60675 features x 195 samples; file had 195 sample columns; labelled samples in the record 198; in both 195 | empty columns 0; duplicated columns 0; labelled not in matrix ['GSM4913634', 'GSM4913644', 'GSM4913645']; matrix not labelled [] |
| b counts | pass | n_samples 195, n_participants 152, infection_acute_covid 46, infection_never_infected 15, other_disease 91 | PMID 33597532 abstract: 46 subjects with COVID-19; differences (built, paper): {} |
| c value type | pass | counts; integer fraction 1.000; range 0 to 1.84e+06; zeros 0.667; missing 0.000; mapping 1.0 (60675 of 60675; 0 targets with duplicates, summed); library size median 5.11e+07 (range 1.63e+06 to 1.3e+08) | columns are the rna_id of the GEO titles |
| d sex | pass | agreement 0.982 on 169 samples (chance 0.502) | inferred {'female': 97, 'male': 93, 'ambiguous': 5} (gmm2:XIST+RPS4Y1); disagreeing ['GSM4913529', 'GSM4913567', 'GSM4913622'] |
| e duplicates (within) | flag | 2 pairs of different participants with correlation above 0.99 | [('GSM4913615', 'GSM4913636', np.float64(0.9912)), ('GSM4913639', 'GSM4913642', np.float64(0.995))]; pairs of the same participant above 0.99: 0; correlation on the 2000 most variable complete features |
| f outliers | flag | 46 samples more than 3 MAD from the cohort median | library size: 20 ['GSM4913501', 'GSM4913510', 'GSM4913520', 'GSM4913538', 'GSM4913566', 'GSM4913571', 'GSM4913604', 'GSM4913613', 'GSM4913620', 'GSM4913623']; detected features: 3 ['GSM4913571', 'GSM4913623', 'GSM4913674']; mean correlation: 28 ['GSM4913535', 'GSM4913560', 'GSM4913564', 'GSM4913568', 'GSM4913569', 'GSM4913570', 'GSM4913571', 'GSM4913572', 'GSM4913573', 'GSM4913574'] |
| g label direction | pass | 0 of 5 markers lower in acute COVID | IFI27 g=1.44; IFI44L g=0.45; IFI6 g=0.37; ISG15 g=0.45; RSAD2 g=0.43 |
| h signal (infection) | report | cross validated AUROC 0.950; split half reliability 0.420 | 46 against 15 participants; logistic regression on the 500 most variable features, 5 x 5-fold |
| i confounding | pass | largest Cramer's V 0.20 | infection x batch: V=0.20 |
| j timing | pass | 46 acute participants; days known for 0 | 0 outside 0 to 30 days or marked outside the window; window unknown for 6; anchor ['symptom_onset'] |
| k demographics | report | acute_covid: n 46, age median 33 (missing 1), sex {'male': 24, 'female': 22}; never_infected: n 15, age median 18 (missing 0), sex {'male': 9, 'female': 6}; other_disease: n 91, age median 18 (missing 0), sex {'female': 42, 'male': 35, 'unknown': 14} | age band columns: [] |

### GSE152641

- **Source:** GSE152641 (RNA, whole blood); files: GSE152641_Inflammatix_COVID19_counts_entrez.csv.gz
- **Paper:** PMID 33437935; DOI 10.1016/j.isci.2020.101947 (Thair et al.)
- **Platform:** Illumina NovaSeq 6000 whole blood RNA-seq; counts; 20155 features x 86 samples
- **Participants (one sample each):** infection 62 acute COVID against 24 never infected; severity 0 severe or critical against 0 mild or moderate; 0 other disease (own class); 0 excluded
- **Contrasts with at least 10 per class:** infection
- **Status:** verified with flags (f outliers: 19 samples more than 3 MAD from the cohort median)
- **Label mapping:** `docs/acute_covid_labels/GSE152641.csv` (2 rows)

| check | status | result | details |
|---|---|---|---|
| a integrity | pass | matrix 20155 features x 86 samples; file had 86 sample columns; labelled samples in the record 86; in both 86 | empty columns 0; duplicated columns 0; labelled not in matrix []; matrix not labelled [] |
| b counts | pass | n_samples 86, n_participants 86, infection_acute_covid 62, infection_never_infected 24 | PMID 33437935 abstract: 62 patients with COVID-19 and 24 healthy controls; differences (built, paper): {} |
| c value type | pass | counts; integer fraction 1.000; range 0 to 1.3e+06; zeros 0.137; missing 0.000; mapping 0.9851 (20155 of 20460; 0 targets with duplicates, summed); library size median 1.28e+07 (range 6.74e+06 to 2.32e+07) | rows are Entrez gene ids |
| d sex | pass | agreement 0.953 on 86 samples (chance 0.553) | inferred {'male': 59, 'female': 27} (gmm2:XIST+RPS4Y1); disagreeing ['GSM4622656', 'GSM4622672', 'GSM4622699', 'GSM4622704'] |
| e duplicates (within) | pass | 0 pairs of different participants with correlation above 0.99 | []; pairs of the same participant above 0.99: 0; correlation on the 2000 most variable complete features |
| f outliers | flag | 19 samples more than 3 MAD from the cohort median | library size: 12 ['GSM4622628', 'GSM4622633', 'GSM4622642', 'GSM4622645', 'GSM4622647', 'GSM4622650', 'GSM4622651', 'GSM4622652', 'GSM4622654', 'GSM4622656']; detected features: 2 ['GSM4622642', 'GSM4622680']; mean correlation: 9 ['GSM4622630', 'GSM4622645', 'GSM4622657', 'GSM4622660', 'GSM4622669', 'GSM4622678', 'GSM4622680', 'GSM4622687', 'GSM4622701'] |
| g label direction | pass | 0 of 5 markers lower in acute COVID | IFI27 g=2.22; IFI44L g=1.45; IFI6 g=1.18; ISG15 g=1.04; RSAD2 g=1.25 |
| h signal (infection) | report | cross validated AUROC 0.992; split half reliability 0.905 | 62 against 24 participants; logistic regression on the 500 most variable features, 5 x 5-fold |
| i confounding | nan | no batch field recorded |  |
| j timing | pass | 62 acute participants; days known for 0 | 0 outside 0 to 30 days or marked outside the window; window unknown for 0; anchor [] |
| k demographics | report | acute_covid: n 62, age median 62 (missing 1), sex {'male': 40, 'female': 22}; never_infected: n 24, age median 36 (missing 0), sex {'male': 17, 'female': 7} | age band columns: [] |

### GSE300129

- **Source:** GSE300129 (RNA, whole blood); files: GSE300129_series_matrix.txt.gz, GPL34133_family.soft.gz
- **Paper:** PMID 41695759
- **Platform:** Affymetrix array (series matrix values); normalised (log2 array intensities, series matrix); 18387 features x 112 samples
- **Participants (one sample each):** infection 57 acute COVID against 55 never infected; severity 57 severe or critical against 0 mild or moderate; 0 other disease (own class); 0 excluded
- **Contrasts with at least 10 per class:** infection
- **Status:** failed (d sex (FAIL): agreement 0.741 on 112 samples (chance 0.755); no better than chance: the sample metadata may not match the expression columns; f outliers: 8 samples more than 3 MAD from the cohort median; g label direction: 5 of 5 markers lower in acute COVID)
- **Label mapping:** `docs/acute_covid_labels/GSE300129.csv` (3 rows)

| check | status | result | details |
|---|---|---|---|
| a integrity | pass | matrix 18387 features x 112 samples; file had 112 sample columns; labelled samples in the record 112; in both 112 | empty columns 0; duplicated columns 0; labelled not in matrix []; matrix not labelled [] |
| b counts | pass | n_samples 112, n_participants 112, infection_acute_covid 57, infection_never_infected 55, severity_severe_critical 57 | PMID 41695759 abstract: 57 mechanically ventilated patients (36 with and 21 without thrombosis); differences (built, paper): {} |
| c value type | pass | normalised (log2 array intensities, series matrix); integer fraction 0.000; range 0.993 to 13.1; zeros 0.000; missing 0.000; mapping 0.9898 (18387 of 18576; 0 targets with duplicates, feature with the highest mean kept) | probe sets named <Entrez id>_at (GPL34133, ENTREZID column) |
| d sex | fail | agreement 0.741 on 112 samples (chance 0.755); no better than chance: the sample metadata may not match the expression columns | inferred {'male': 95, 'female': 17} (gmm2:RPS4Y1 (no XIST row in the matrix)); disagreeing ['GSM9054079', 'GSM9054080', 'GSM9054082', 'GSM9054085', 'GSM9054089', 'GSM9054095', 'GSM9054102', 'GSM9054110', 'GSM9054115', 'GSM9054116', 'GSM9054119', 'GSM9054124', 'GSM9054128', 'GSM9054134', 'GSM9054138', 'GSM9054144', 'GSM9054148', 'GSM9054152', 'GSM9054154', 'GSM9054164', 'GSM9054165', 'GSM9054167',  |
| e duplicates (within) | pass | 0 pairs of different participants with correlation above 0.99 | []; pairs of the same participant above 0.99: 0; correlation on the 2000 most variable complete features |
| f outliers | flag | 8 samples more than 3 MAD from the cohort median | library size: 0 []; detected features: 0 []; mean correlation: 8 ['GSM9054087', 'GSM9054111', 'GSM9054120', 'GSM9054123', 'GSM9054139', 'GSM9054143', 'GSM9054148', 'GSM9054150'] |
| g label direction | flag | 5 of 5 markers lower in acute COVID | IFI27 g=-0.34; IFI44L g=-0.02; IFI6 g=-0.09; ISG15 g=-0.26; RSAD2 g=-0.17 |
| h signal (infection) | report | cross validated AUROC 0.625; split half reliability 0.677 | 57 against 55 participants; logistic regression on the 500 most variable features, 5 x 5-fold |
| i confounding | nan | no batch field recorded |  |
| j timing | pass | 57 acute participants; days known for 0 | 0 outside 0 to 30 days or marked outside the window; window unknown for 0; anchor ['admission'] |
| k demographics | report | acute_covid: n 57, age median 65 (missing 0), sex {'male': 49, 'female': 8}; never_infected: n 55, age median 66 (missing 0), sex {'male': 47, 'female': 8} | age band columns: [] |

### GSE243219

- **Source:** GSE243219 (RNA, PBMC); files: GSE243219_RAW.tar
- **Paper:** PMID 39160575
- **Platform:** Illumina RNA-seq of PBMC; counts (expected counts, not all integers); 24392 features x 72 samples
- **Participants (one sample each):** infection 35 acute COVID against 15 never infected; severity 0 severe or critical against 0 mild or moderate; 22 other disease (own class); 0 excluded
- **Contrasts with at least 10 per class:** infection
- **Status:** verified with flags (f outliers: 14 samples more than 3 MAD from the cohort median)
- **Label mapping:** `docs/acute_covid_labels/GSE243219.csv` (3 rows)

| check | status | result | details |
|---|---|---|---|
| a integrity | pass | matrix 24392 features x 72 samples; file had 72 sample columns; labelled samples in the record 72; in both 72 | empty columns 0; duplicated columns 0; labelled not in matrix []; matrix not labelled [] |
| b counts | pass | n_samples 72, n_participants 72, infection_acute_covid 35, infection_never_infected 15, other_disease 22 | PMID 39160575 abstract: 35 COVID-19, 22 bacterial sepsis, 15 healthy; differences (built, paper): {} |
| c value type | pass | counts (expected counts, not all integers); integer fraction 0.811; range 0 to 1.33e+07; zeros 0.297; missing 0.000; mapping 0.9279 (24464 of 26364; 70 targets with duplicates, summed); library size median 2.22e+07 (range 1.37e+07 to 2.9e+07) | one file per sample (GSM in the file name); rows are gene symbols |
| d sex | nan | no reported sex to compare | inferred {'male': 47, 'female': 24, 'ambiguous': 1} (gmm2:XIST+RPS4Y1); disagreeing [] |
| e duplicates (within) | pass | 0 pairs of different participants with correlation above 0.99 | []; pairs of the same participant above 0.99: 0; correlation on the 2000 most variable complete features |
| f outliers | flag | 14 samples more than 3 MAD from the cohort median | library size: 4 ['GSM7780772', 'GSM7780779', 'GSM7780803', 'GSM7780807']; detected features: 8 ['GSM7780772', 'GSM7780773', 'GSM7780776', 'GSM7780781', 'GSM7780824', 'GSM7780841', 'GSM7780842', 'GSM7780843']; mean correlation: 6 ['GSM7780776', 'GSM7780797', 'GSM7780828', 'GSM7780830', 'GSM7780842', 'GSM7780843'] |
| g label direction | pass | 0 of 5 markers lower in acute COVID | IFI27 g=2.50; IFI44L g=1.33; IFI6 g=1.18; ISG15 g=1.10; RSAD2 g=1.43 |
| h signal (infection) | report | cross validated AUROC 0.995; split half reliability 0.813 | 35 against 15 participants; logistic regression on the 500 most variable features, 5 x 5-fold |
| i confounding | pass | largest Cramer's V 0.00 | infection x platform_id: one level |
| j timing | pass | 35 acute participants; days known for 0 | 0 outside 0 to 30 days or marked outside the window; window unknown for 0; anchor ['admission'] |
| k demographics | report | acute_covid: n 35, age missing, sex {'unknown': 35}; never_infected: n 15, age missing, sex {'unknown': 15}; other_disease: n 22, age missing, sex {'unknown': 22} | age band columns: [] |

### GSE172114

- **Source:** GSE172114 (RNA, whole blood); files: GSE172114_rsem_gene_count_matrix_TMM_69samples.csv.gz
- **Paper:** PMID 34698500; DOI 10.1126/scitranslmed.abj7521 (Carapito et al.)
- **Platform:** Illumina NovaSeq 6000 whole blood RNA-seq (RSEM, TMM); normalised (TMM, log scale as deposited); 15957 features x 69 samples
- **Participants (one sample each):** infection 69 acute COVID against 0 never infected; severity 46 severe or critical against 23 mild or moderate; 0 other disease (own class); 0 excluded
- **Contrasts with at least 10 per class:** severity
- **Status:** verified with flags (b counts: n_samples 69, n_participants 69, infection_acute_covid 69, severity_severe_critical 46, severity_mild_moderate 23; f outliers: 6 samples more than 3 MAD from the cohort median)
- **Label mapping:** `docs/acute_covid_labels/GSE172114.csv` (2 rows)

| check | status | result | details |
|---|---|---|---|
| a integrity | pass | matrix 15957 features x 69 samples; file had 69 sample columns; labelled samples in the record 69; in both 69 | empty columns 0; duplicated columns 0; labelled not in matrix []; matrix not labelled [] |
| b counts | flag | n_samples 69, n_participants 69, infection_acute_covid 69, severity_severe_critical 46, severity_mild_moderate 23 | PMID 34698500 abstract (whole cohort): 47 critical, 25 non-critical, 22 healthy; differences (built, paper): {'severity_severe_critical': (46, 47), 'severity_mild_moderate': (23, 25)} |
| c value type | pass | normalised (TMM, log scale as deposited); integer fraction 0.001; range -4.28 to 18.3; zeros 0.000; missing 0.000; mapping 1.0 (15957 of 15957; 0 targets with duplicates, feature with the highest mean kept) | rows are Ensembl id with version and symbol |
| d sex | nan | no reported sex to compare | inferred {'male': 51, 'female': 18} (gmm2:XIST+RPS4Y1); disagreeing [] |
| e duplicates (within) | pass | 0 pairs of different participants with correlation above 0.99 | []; pairs of the same participant above 0.99: 0; correlation on the 2000 most variable complete features |
| f outliers | flag | 6 samples more than 3 MAD from the cohort median | library size: 0 []; detected features: 0 []; mean correlation: 6 ['GSM5242130', 'GSM5242149', 'GSM5242180', 'GSM5242181', 'GSM5242189', 'GSM5242191'] |
| g label direction | nan | no infection contrast |  |
| h signal (severity) | report | cross validated AUROC 0.951; split half reliability 0.922 | 46 against 23 participants; logistic regression on the 500 most variable features, 5 x 5-fold |
| i confounding | nan | no batch field recorded |  |
| j timing | pass | 69 acute participants; days known for 0 | 0 outside 0 to 30 days or marked outside the window; window unknown for 0; anchor ['admission'] |
| k demographics | report | acute_covid: n 69, age missing, sex {'unknown': 69} | age band columns: [] |

### GSE171110

- **Source:** GSE171110 (RNA, whole blood); files: GSE171110_Data_RNAseq_raw_counts_geo.txt.gz
- **Paper:** PMID 34127958; DOI 10.1172/jci.insight.148895 (Levy et al.)
- **Platform:** Illumina HiSeq 2500 whole blood RNA-seq; counts (expected counts, not all integers); 25146 features x 54 samples
- **Participants (one sample each):** infection 44 acute COVID against 10 never infected; severity 44 severe or critical against 0 mild or moderate; 0 other disease (own class); 0 excluded
- **Contrasts with at least 10 per class:** infection
- **Status:** verified with flags (b counts: n_samples 54, n_participants 54, infection_acute_covid 44, infection_never_infected 10, severity_severe_critical 44; f outliers: 7 samples more than 3 MAD from the cohort median)
- **Label mapping:** `docs/acute_covid_labels/GSE171110.csv` (2 rows)

| check | status | result | details |
|---|---|---|---|
| a integrity | pass | matrix 25146 features x 54 samples; file had 54 sample columns; labelled samples in the record 54; in both 54 | empty columns 0; duplicated columns 0; labelled not in matrix []; matrix not labelled [] |
| b counts | flag | n_samples 54, n_participants 54, infection_acute_covid 44, infection_never_infected 10, severity_severe_critical 44 | PMID 34127958 abstract (whole study): 61 severe hospitalised patients; differences (built, paper): {'infection_acute_covid': (44, 61)} |
| c value type | pass | counts (expected counts, not all integers); integer fraction 0.396; range 0 to 1.15e+06; zeros 0.249; missing 0.000; mapping 0.8331 (25148 of 30185; 2 targets with duplicates, summed); library size median 3.61e+07 (range 2.73e+07 to 5.71e+07) | rows are gene symbols; columns are the Sample_sheet ids of the GEO titles |
| d sex | nan | no reported sex to compare | inferred {'male': 47, 'female': 7} (gmm2:RPS4Y1 (no XIST row in the matrix)); disagreeing [] |
| e duplicates (within) | pass | 0 pairs of different participants with correlation above 0.99 | []; pairs of the same participant above 0.99: 0; correlation on the 2000 most variable complete features |
| f outliers | flag | 7 samples more than 3 MAD from the cohort median | library size: 5 ['GSM5219000', 'GSM5219001', 'GSM5219037', 'GSM5219038', 'GSM5219039']; detected features: 3 ['GSM5218994', 'GSM5219001', 'GSM5219039']; mean correlation: 1 ['GSM5219021'] |
| g label direction | pass | 0 of 5 markers lower in acute COVID | IFI27 g=2.07; IFI44L g=0.84; IFI6 g=1.19; ISG15 g=0.53; RSAD2 g=0.85 |
| h signal (infection) | report | cross validated AUROC 0.998; split half reliability 0.896 | 44 against 10 participants; logistic regression on the 500 most variable features, 5 x 5-fold |
| i confounding | nan | no batch field recorded |  |
| j timing | pass | 44 acute participants; days known for 0 | 0 outside 0 to 30 days or marked outside the window; window unknown for 0; anchor [] |
| k demographics | report | acute_covid: n 44, age missing, sex {'unknown': 44}; never_infected: n 10, age missing, sex {'unknown': 10} | age band columns: [] |

### GSE213313

- **Source:** GSE213313 (RNA, whole blood); files: GSE213313_series_matrix.txt.gz, GPL21185_family.soft.gz
- **Paper:** PMID 37365222
- **Platform:** Affymetrix array (series matrix values); normalised (log2 array intensities, series matrix); 14866 features x 94 samples
- **Participants (one sample each):** infection 34 acute COVID against 11 never infected; severity 19 severe or critical against 15 mild or moderate; 0 other disease (own class); 0 excluded
- **Contrasts with at least 10 per class:** infection, severity
- **Status:** verified with flags (b counts: n_samples 94, n_participants 45, infection_acute_covid 34, infection_never_infected 11, severity_severe_critical 19, severity_mild_moderate 15; f outliers: 5 samples more than 3 MAD from the cohort median)
- **Label mapping:** `docs/acute_covid_labels/GSE213313.csv` (7 rows)

| check | status | result | details |
|---|---|---|---|
| a integrity | pass | matrix 14866 features x 94 samples; file had 94 sample columns; labelled samples in the record 94; in both 94 | empty columns 0; duplicated columns 0; labelled not in matrix []; matrix not labelled [] |
| b counts | flag | n_samples 94, n_participants 45, infection_acute_covid 34, infection_never_infected 11, severity_severe_critical 19, severity_mild_moderate 15 | PMID 37365222 abstract: 17 severe, 15 moderate, 11 healthy (GEO: 19 critical, 15 non-critical); differences (built, paper): {'severity_severe_critical': (19, 17)} |
| c value type | pass | normalised (log2 array intensities, series matrix); integer fraction 0.000; range -2.22 to 18.5; zeros 0.000; missing 0.000; mapping 0.6882 (17527 of 25469; 2280 targets with duplicates, feature with the highest mean kept) | Agilent probes (GPL21185): ENSEMBL_ID, else GENE_SYMBOL through HGNC |
| d sex | nan | no reported sex to compare | inferred {'male': 59, 'female': 35} (gmm2:XIST+RPS4Y1); disagreeing [] |
| e duplicates (within) | pass | 0 pairs of different participants with correlation above 0.99 | []; pairs of the same participant above 0.99: 0; correlation on the 2000 most variable complete features |
| f outliers | flag | 5 samples more than 3 MAD from the cohort median | library size: 0 []; detected features: 0 []; mean correlation: 5 ['GSM6578729', 'GSM6578743', 'GSM6578758', 'GSM6578778', 'GSM6578811'] |
| g label direction | pass | 0 of 5 markers lower in acute COVID | IFI27 g=4.90; IFI44L g=2.66; IFI6 g=2.64; ISG15 g=2.01; RSAD2 g=2.60 |
| h signal (infection) | report | cross validated AUROC 1.000; split half reliability 0.885 | 34 against 11 participants; logistic regression on the 500 most variable features, 5 x 5-fold |
| h signal (severity) | report | cross validated AUROC 0.710; split half reliability 0.470 | 19 against 15 participants; logistic regression on the 500 most variable features, 5 x 5-fold |
| i confounding | nan | no batch field recorded |  |
| j timing | pass | 34 acute participants; days known for 0 | 0 outside 0 to 30 days or marked outside the window; window unknown for 0; anchor ['admission'] |
| k demographics | report | acute_covid: n 34, age missing, sex {'unknown': 34}; never_infected: n 11, age missing, sex {'unknown': 11} | age band columns: [] |

### GSE152418

- **Source:** GSE152418 (RNA, PBMC); files: GSE152418_p20047_Study1_RawCounts.txt.gz
- **Paper:** PMID 32788292; DOI 10.1126/science.abc6261 (Arunachalam et al.)
- **Platform:** Illumina NovaSeq 6000 PBMC RNA-seq; counts; 60683 features x 34 samples
- **Participants (one sample each):** infection 16 acute COVID against 17 never infected; severity 4 severe or critical against 4 mild or moderate; 0 other disease (own class); 1 excluded
- **Contrasts with at least 10 per class:** infection
- **Status:** verified with flags (b counts: n_samples 34, n_participants 34, infection_acute_covid 16, infection_never_infected 17, excluded 1, severity_severe_critical 4, severity_mild_moderate 4; f outliers: 5 samples more than 3 MAD from the cohort median)
- **Label mapping:** `docs/acute_covid_labels/GSE152418.csv` (5 rows)

| check | status | result | details |
|---|---|---|---|
| a integrity | pass | matrix 60683 features x 34 samples; file had 34 sample columns; labelled samples in the record 34; in both 34 | empty columns 0; duplicated columns 0; labelled not in matrix []; matrix not labelled [] |
| b counts | flag | n_samples 34, n_participants 34, infection_acute_covid 16, infection_never_infected 17, excluded 1, severity_severe_critical 4, severity_mild_moderate 4 | PMID 32788292 abstract (all assays): 76 patients, 69 healthy; differences (built, paper): {'infection_acute_covid': (16, 76), 'infection_never_infected': (17, 69)} |
| c value type | pass | counts; integer fraction 1.000; range 0 to 1.17e+06; zeros 0.520; missing 0.000; mapping 1.0 (60683 of 60683; 0 targets with duplicates, summed); library size median 1.59e+07 (range 1.3e+07 to 2.05e+07) |  |
| d sex | pass | agreement 1.000 on 34 samples (chance 0.502) | inferred {'female': 18, 'male': 16} (gmm2:XIST+RPS4Y1); disagreeing [] |
| e duplicates (within) | pass | 0 pairs of different participants with correlation above 0.99 | []; pairs of the same participant above 0.99: 0; correlation on the 2000 most variable complete features |
| f outliers | flag | 5 samples more than 3 MAD from the cohort median | library size: 0 []; detected features: 3 ['GSM4614993', 'GSM4614995', 'GSM4614996']; mean correlation: 5 ['GSM4614992', 'GSM4614993', 'GSM4614994', 'GSM4614995', 'GSM4614996'] |
| g label direction | pass | 0 of 5 markers lower in acute COVID | IFI27 g=3.51; IFI44L g=0.55; IFI6 g=0.77; ISG15 g=0.57; RSAD2 g=0.40 |
| h signal (infection) | report | cross validated AUROC 1.000; split half reliability 0.744 | 16 against 17 participants; logistic regression on the 500 most variable features, 5 x 5-fold |
| i confounding | nan | no batch field recorded |  |
| j timing | pass | 16 acute participants; days known for 16; median 14 (IQR 9 to 16, range 2 to 23) | 0 outside 0 to 30 days or marked outside the window; window unknown for 0; anchor ['symptom_onset'] |
| k demographics | report | acute_covid: n 16, age missing, sex {'female': 9, 'male': 7}; excluded: n 1, age missing, sex {'male': 1}; never_infected: n 17, age missing, sex {'female': 9, 'male': 8} | age band columns: [] |

### Filbin2021_MGH

- **Source:** 10.17632/nf853r8xsj (proteins, plasma); files: Olink Proteomics.xlsx, Clinical Metadata.xlsx, Supplemental-Table-2-Olink-Assays-NPX-values_v2.xlsx
- **Paper:** DOI 10.1016/j.xcrm.2021.100287 (Filbin et al. 2021, Cell Reports Medicine)
- **Platform:** Olink Explore 1536 (NPX); NPX (log2 scale); 1420 features x 382 samples
- **Participants (one sample each):** infection 305 acute COVID against 0 never infected; severity 109 severe or critical against 63 mild or moderate; 77 other disease (own class); 0 excluded
- **Contrasts with at least 10 per class:** severity
- **Status:** verified with flags (a integrity: matrix 1420 features x 382 samples; file had 784 sample columns; labelled samples in the record 383; in both 382 (unlabelled samples expected: later draws (day 3, day 7, event) are not labelled: one draw per participant, day 0); b counts: n_samples 382, n_participants 382, infection_acute_covid 305, other_disease 77, severity_severe_critical 109, severity_mild_moderate 63; f outliers: 53 samples more than 3 MAD from the cohort median)
- **Label mapping:** `docs/acute_covid_labels/Filbin2021_MGH.csv` (10 rows)

| check | status | result | details |
|---|---|---|---|
| a integrity | flag | matrix 1420 features x 382 samples; file had 784 sample columns; labelled samples in the record 383; in both 382 (unlabelled samples expected: later draws (day 3, day 7, event) are not labelled: one draw per participant, day 0) | empty columns 0; duplicated columns 0; labelled not in matrix ['180_D0']; matrix not labelled ['100_D3', '101_D3', '101_D7', '102_D3', '102_D7', '102_DE', '106_D3', '106_D7', '108_D3', '114_D3'] |
| b counts | flag | n_samples 382, n_participants 382, infection_acute_covid 305, other_disease 77, severity_severe_critical 109, severity_mild_moderate 63 | Filbin 2021 abstract (DOI 10.1016/j.xcrm.2021.100287): 306 COVID-19 patients, 78 symptomatic controls; differences (built, paper): {'infection_acute_covid': (305, 306), 'other_disease': (77, 78)} |
| c value type | pass | NPX (log2 scale); integer fraction 0.000; range -7.51 to 15.7; zeros 0.000; missing 0.000; mapping 1.0 (1429 of 1429; 3 targets with duplicates, feature with the highest mean kept) | Olink assay ids to UniProt by Supplemental Table 2A |
| d sex | nan | protein panel without sex chromosome markers |  |
| e duplicates (within) | pass | 0 pairs of different participants with correlation above 0.99 | []; pairs of the same participant above 0.99: 0; correlation on the 1417 most variable complete features |
| f outliers | flag | 53 samples more than 3 MAD from the cohort median | library size: 0 []; detected features: 11 ['57_D0', '58_D0', '131_D0', '140_D0', '151_D0', '245_D0', '259_D0', '300_D0', '350_D0', '365_D0']; mean correlation: 44 ['2_D0', '16_D0', '21_D0', '23_D0', '29_D0', '59_D0', '78_D0', '110_D0', '116_D0', '143_D0'] |
| g label direction | nan | no infection contrast |  |
| h signal (severity) | report | cross validated AUROC 0.923; split half reliability 0.888 | 109 against 63 participants; logistic regression on the 500 most variable features, 5 x 5-fold |
| i confounding | nan | no batch field recorded |  |
| j timing | pass | 305 acute participants; days known for 305; median 0 (IQR 0 to 0, range 0 to 0) | 0 outside 0 to 30 days or marked outside the window; window unknown for 0; anchor ['admission'] |
| k demographics | report | acute_covid: n 305, age missing, sex {'unknown': 305}; other_disease: n 77, age missing, sex {'unknown': 77} | age band columns: ['age_category'] |

### Gisby2021_plasma

- **Source:** zenodo.4287557 (proteins, plasma); files: plasma_npx_level.csv, plasma_sample_level.csv, README.pdf
- **Paper:** PMID 33704068; DOI 10.7554/eLife.64827 (Gisby et al.)
- **Platform:** Olink Target 96 (NPX); NPX (log2 scale); 436 features x 307 samples
- **Participants (one sample each):** infection 55 acute COVID against 50 never infected; severity 13 severe or critical against 28 mild or moderate; 0 other disease (own class); 0 excluded
- **Contrasts with at least 10 per class:** infection, severity
- **Status:** verified with flags (b counts: n_samples 307, n_participants 105, infection_acute_covid 55, infection_never_infected 50, severity_severe_critical 13, severity_mild_moderate 28; f outliers: 29 samples more than 3 MAD from the cohort median; j timing: 55 acute participants; days known for 52; median 12 (IQR 8 to 14, range 0 to 37))
- **Label mapping:** `docs/acute_covid_labels/Gisby2021_plasma.csv` (6 rows)

| check | status | result | details |
|---|---|---|---|
| a integrity | pass | matrix 436 features x 307 samples; file had 307 sample columns; labelled samples in the record 307; in both 307 | empty columns 0; duplicated columns 0; labelled not in matrix []; matrix not labelled [] |
| b counts | flag | n_samples 307, n_participants 105, infection_acute_covid 55, infection_never_infected 50, severity_severe_critical 13, severity_mild_moderate 28 | PMID 33704068 abstract: 256 samples from 55 patients with COVID-19; 51 non-infected patients; differences (built, paper): {'infection_never_infected': (50, 51)} |
| c value type | pass | NPX (log2 scale); integer fraction 0.000; range -2.38 to 17.5; zeros 0.000; missing 0.002; mapping 1.0 (436 of 436; 0 targets with duplicates, feature with the highest mean kept) | long table pivoted; an assay measured on two panels is averaged by the pivot |
| d sex | nan | protein panel without sex chromosome markers |  |
| e duplicates (within) | pass | 0 pairs of different participants with correlation above 0.99 | []; pairs of the same participant above 0.99: 162; correlation on the 261 most variable complete features |
| f outliers | flag | 29 samples more than 3 MAD from the cohort median | library size: 0 []; detected features: 15 ['C10_6', 'C40_26', 'C167_106', 'C23_172', 'C7_178', 'C9_191', 'C20_226', 'C27_249', 'C36_287', 'C18_302']; mean correlation: 14 ['C65_11', 'C10_40', 'C10_47', 'C32_188', 'C9_213', 'C11_231', 'C11_237', 'C24_261', 'C99_299', 'C49_305'] |
| g label direction | pass | 0 of 1 markers lower in acute COVID | CXCL10 g=1.00 |
| h signal (infection) | report | cross validated AUROC 0.906; split half reliability 0.862 | 55 against 50 participants; logistic regression on the 500 most variable features, 5 x 5-fold |
| h signal (severity) | report | cross validated AUROC 0.898; split half reliability 0.740 | 13 against 28 participants; logistic regression on the 500 most variable features, 5 x 5-fold |
| i confounding | pass | largest Cramer's V 0.40 | infection x plate: V=0.25; severity x plate: V=0.40 |
| j timing | flag | 55 acute participants; days known for 52; median 12 (IQR 8 to 14, range 0 to 37) | 2 outside 0 to 30 days or marked outside the window; window unknown for 3; anchor ['symptom_onset'] |
| k demographics | report | acute_covid: n 55, age missing, sex {'male': 39, 'female': 16}; never_infected: n 50, age missing, sex {'male': 36, 'female': 14} | age band columns: ['age_band'] |

### GSE207015

- **Source:** GSE207015 (proteins, plasma); files: GSE207015_series_matrix.txt.gz, GPL32354_family.soft.gz
- **Paper:** see GEO record
- **Platform:** SomaScan 7k v4.1 (series matrix values); SomaScan RFU (log2 as deposited); 6333 features x 194 samples
- **Participants (one sample each):** infection 124 acute COVID against 70 never infected; severity 56 severe or critical against 44 mild or moderate; 0 other disease (own class); 0 excluded
- **Contrasts with at least 10 per class:** infection, severity
- **Status:** verified with flags (b counts: n_samples 194, n_participants 194, infection_acute_covid 124, infection_never_infected 70, severity_severe_critical 56, severity_mild_moderate 44; e duplicates (within): 1 pairs of different participants with correlation above 0.99; f outliers: 25 samples more than 3 MAD from the cohort median)
- **Label mapping:** `docs/acute_covid_labels/GSE207015.csv` (6 rows)

| check | status | result | details |
|---|---|---|---|
| a integrity | pass | matrix 6333 features x 194 samples; file had 194 sample columns; labelled samples in the record 194; in both 194 | empty columns 0; duplicated columns 0; labelled not in matrix []; matrix not labelled [] |
| b counts | flag | n_samples 194, n_participants 194, infection_acute_covid 124, infection_never_infected 70, severity_severe_critical 56, severity_mild_moderate 44 | no paper linked in GEO; the GEO record lists 124 cases and 70 controls; no comparable number |
| c value type | pass | SomaScan RFU (log2 as deposited); integer fraction 0.000; range 3.07 to 18.6; zeros 0.000; missing 0.000; mapping 0.989 (7208 of 7288; 771 targets with duplicates, feature with the highest mean kept) | aptamers to UniProt by GPL32354; aptamers to several proteins dropped |
| d sex | nan | protein panel without sex chromosome markers |  |
| e duplicates (within) | flag | 1 pairs of different participants with correlation above 0.99 | [('GSM6268178', 'GSM6268180', np.float64(0.9926))]; pairs of the same participant above 0.99: 0; correlation on the 2000 most variable complete features |
| f outliers | flag | 25 samples more than 3 MAD from the cohort median | library size: 0 []; detected features: 0 []; mean correlation: 25 ['GSM6268055', 'GSM6268065', 'GSM6268084', 'GSM6268087', 'GSM6268102', 'GSM6268106', 'GSM6268150', 'GSM6268155', 'GSM6268160', 'GSM6268164'] |
| g label direction | pass | 0 of 1 markers lower in acute COVID | CXCL10 g=1.42 |
| h signal (infection) | report | cross validated AUROC 0.972; split half reliability 0.915 | 124 against 70 participants; logistic regression on the 500 most variable features, 5 x 5-fold |
| h signal (severity) | report | cross validated AUROC 0.915; split half reliability 0.867 | 56 against 44 participants; logistic regression on the 500 most variable features, 5 x 5-fold |
| i confounding | pass | largest Cramer's V 0.00 | infection x pregnancy: one level; severity x pregnancy: one level |
| j timing | pass | 124 acute participants; days known for 0 | 0 outside 0 to 30 days or marked outside the window; window unknown for 124; anchor [] |
| k demographics | report | acute_covid: n 124, age median 35 (missing 0), sex {'female': 102, 'male': 22}; never_infected: n 70, age median 38 (missing 0), sex {'female': 48, 'male': 22} | age band columns: [] |


## All candidates

13 included, 257 needing a decision, 2825 excluded. The same table with every column is
`docs/acute_covid_candidates.csv`.

| accession | source | layer | title | decision | reason |
|---|---|---|---|---|---|
| GSE152418 | GEO | RNA | Systems biological assessment of immunity to severe and mild COVID-19 infections | include | included: open, per sample labels in the repository, adults, blood, acute; deposited matrix; infection (PBMC) |
| GSE152641 | GEO | RNA | Transcriptomic Similarities and Differences in Host Response between SARS-CoV-2 and Other  | include | included: open, per sample labels in the repository, adults, blood, acute; deposited matrix; infection |
| GSE161731 | GEO | RNA | Dysregulated transcriptional responses to SARS-CoV-2 in the periphery support novel diagno | include | included: open, per sample labels in the repository, adults, blood, acute; deposited matrix; infection |
| GSE171110 | GEO | RNA | Whole blood transciptome analysis in severe Covid-19  patients and healthy donors | include | included: open, per sample labels in the repository, adults, blood, acute; deposited matrix; infection (10 healthy donors) |
| GSE172114 | GEO | RNA | Identification of driver genes for severe forms of COVID-19 in a deeply phenotyped young p | include | included: open, per sample labels in the repository, adults, blood, acute; deposited matrix; severity (critical against non-critical) |
| GSE198449 | GEO | RNA | RNA-seq data to study immune response to SARS-CoV-2 virus infection in COVID-19 Health Act | include | included: open, per sample labels in the repository, adults, blood, acute; deposited matrix; infection (first PCR positive sample against participants never positive) |
| GSE213313 | GEO | RNA | Serial whole blood transcriptomic analysis demonstrates activation of neutrophils, platele | include | included: open, per sample labels in the repository, adults, blood, acute; deposited matrix; infection and severity |
| GSE243219 | GEO | RNA | Comparative Analysis of Host Responses Between Sepsis and COVID-19: A Prospective Observat | include | included: open, per sample labels in the repository, adults, blood, acute; deposited matrix; infection (PBMC) |
| GSE300129 | GEO | RNA | Whole blood gene expression data of thrombotic and non-thrombotic critically ill COVID-19  | include | included: open, per sample labels in the repository, adults, blood, acute; deposited matrix; infection (ventilated patients against healthy controls) |
| GSE207015 | GEO | proteins | COVID-19 drives a distinct plasma proteome in pregnant and non-pregnant individuals | include | included: open, per sample labels in the repository, adults, blood, acute; deposited matrix; infection and severity; pregnant and non-pregnant adults |
| 10.17632/nf853r8xsj | Mendeley Data | proteins | Plasma proteomics reveals tissue-specific cell death and mediators of cell-cell interactio | include | included: open, per sample labels in the repository, adults, blood, acute; deposited matrix; severity (COVID-negative patients are other disease, so no infection contrast) |
| zenodo.4287557 | Zenodo | proteins | Longitudinal proteomic profiling of high-risk patients with COVID-19 reveals markers of se | include | included: open, per sample labels in the repository, adults, blood, acute; deposited matrix; infection and severity; all participants on haemodialysis |
| zenodo.6120249 | Zenodo | proteins | A blood atlas of COVID-19 defines hallmarks of disease severity and specificity: Associate | include | included: open, per sample labels in the repository, adults, blood, acute; deposited matrix; infection and severity (whole blood RNA-seq archive of the COMBAT record) |
| E-MTAB-10022 | ArrayExpress | RNA | Blood RNA-seq data from individuals during the UK Covid-19 lockdown, spanning sequential w | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| E-MTAB-10298 | ArrayExpress | RNA | ncRNA-Seq of human blood with different disease state of COVID-19 | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| E-MTAB-10431 | ArrayExpress | RNA | High titers and low fucosylation of early anti-SARS-CoV-2 IgG promote hyper-inflammation b | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| E-MTAB-10527 | ArrayExpress | RNA | Blood RNA sequencing was performed on a cohort of adults attending the Emergency Departmen | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| E-MTAB-10926 | ArrayExpress | RNA | Characterisation of the blood RNA host response underpinning severity in COVID-19 patients | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| E-MTAB-11579 | ArrayExpress | RNA | Transcriptomic profile of human mesenchymal stem cells after statins treatment | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| E-MTAB-12791 | ArrayExpress | RNA | Analysis of blood and nasal epithelial transcriptomes to identify mechanisms associated wi | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| E-MTAB-12993 | ArrayExpress | RNA | RNA-seq of longitudinal human blood and nose swab samples from a controlled human infectio | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| E-MTAB-13307 | ArrayExpress | RNA | Blood RNA sequencing was performed on a cohort of adults attending the Emergency Departmen | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| E-MTAB-14545 | ArrayExpress | RNA | RNA-seq of PBMCs isolated from patients with various cancers and COVID-19 infection status | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| E-MTAB-15068 | ArrayExpress | RNA | Plasmablasts in COVID-19: Dual Role as Immune Effectors or Amplifiers of Inflammation | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| E-MTAB-15584 | ArrayExpress | RNA | RNA-seq of human Peripheral blood mononuclear cells (PBMCs) from COVID-19 patients requiri | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| E-MTAB-8871 | ArrayExpress | RNA | Transcriptomic analysis of immune response in healthy controls and COVID-19 cases using th | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| E-MTAB-9719 | ArrayExpress | RNA | Daily transcriptomic profiling of immune response in COVID-19 cases and healthy controls u | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| E-MTAB-9721 | ArrayExpress | RNA | Temporal dynamics of the host molecular responses underlying severe COVID-19 progression a | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| PRJEB47548 | BioProject | RNA | Characterisation of the blood RNA host response underpinning severity in COVID-19 patients | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| PRJNA1090478 | BioProject | RNA | Bone marrow megakaryocyte's transcriptome in SARS-CoV-2 infection | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| PRJNA1111567 | BioProject | RNA | Cytokine-induced memory-like responses in endothelial cells link chronic inflammation to v | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| PRJNA1261932 | BioProject | RNA | COVID19 human transcriptome | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| PRJNA201433 | BioProject | RNA | Identification of differentially expressed transcripts and pathways one week and six month | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| PRJNA230906 | BioProject | RNA | Global gene expression changes in peripheral blood during H7N9 infection | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| PRJNA232593 | BioProject | RNA | Variation in RNA-Seq transcriptome profiles of peripheral whole blood from healthy individ | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| PRJNA305001 | BioProject | RNA | IL-10 dysregulation in acute mountain sickness revealed by transcriptome analysis | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| PRJNA329148 | BioProject | RNA | Whole blood RNA sequencing in idiopathic pneumonia syndrome reveals a unique transcriptomi | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| PRJNA380820 | BioProject | RNA | Transcriptional responses induced by controlled human malaria infection (CHMI) | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| PRJNA430406 | BioProject | RNA | ROP: Dumpster Diving in RNA-sequencing to find the source of 1 trillion reads across diver | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| PRJNA434133 | BioProject | RNA | Microbiome and Inflammatory Interactions in Obese and Severe Asthmatic Adults | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| PRJNA454445 | BioProject | RNA | TimeSignature: A Universal Method for Robust Detection of Circadian State from Gene Expres | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| PRJNA454694 | BioProject | RNA | Sequencing based maternal whole blood expression changes with gestational age and labor in | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| PRJNA476781 | BioProject | RNA | Discovering in vivo cytokine eQTL interactions from a lupus clinical trial | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| PRJNA493832 | BioProject | RNA | The Regulation of IFN Type I Pathway Related Genes RSAD2 and ETV7 Specifically Indicate An | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| PRJNA494155 | BioProject | RNA | Differential analysis of gene expression profiles in peripheral blood cells of patients wi | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| PRJNA511891 | BioProject | RNA | Acne flare-up in treatment process: investigating the difference of gene expression in per | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| PRJNA562305 | BioProject | RNA | F508del homozygous cystic fibrosis patients and healthy controls [RNA-seq] | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| PRJNA607120 | BioProject | RNA | Evaluation of transcriptomic regulations behind metabolic syndrome in obese and lean subje | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| PRJNA613909 | BioProject | RNA | Inflammasome Signaling and Impaired Vascular Health in Psoriasis | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| PRJNA630674 | BioProject | RNA | Gene expression classifiers identify bacterial and viral etiologies of acute respiratory i | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| PRJNA632871 | BioProject | RNA | Comparative evaluation for the globin gene depletion methods  for mRNA sequencing using th | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| PRJNA638753 | BioProject | RNA | COVID-19 transcriptional response in blood | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| PRJNA680771 | BioProject | RNA | A two-gene signature for diagnosis of tuberculosis in persons with advanced HIV | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| PRJNA683803 | BioProject | RNA | Transcriptomic biomarker pathways associated with death in HIV-infected patients with cryp | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| PRJNA686397 | BioProject | RNA | Blood RNA signature RISK4LEP predicts leprosy years before clinical onset | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| PRJNA691231 | BioProject | RNA | Multi-omic approach identifies a transcriptional network coupling innate immune response t | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| PRJNA691419 | BioProject | RNA | A diagnostic host response biosignature for COVID-19 from RNA profiling of nasal swabs and | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| PRJNA691933 | BioProject | RNA | COVID-19 High-Depth Blood PAXgene RNAseq | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| PRJNA692253 | BioProject | RNA | A blood RNA transcriptome profile for COVID-19 | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| PRJNA703029 | BioProject | RNA | PolyA RNA-seq from whole blood of subjects with Trisomy 21 ± clonal hematopoiesis | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| PRJNA735444 | BioProject | RNA | Transcriptomic analysis of treated MSCs | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| PRJNA777938 | BioProject | RNA | Transcriptomic analysis of critically ill COVID-19 patients and other patients with AHRF i | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| PRJNA779249 | BioProject | RNA | Homo sapiens Transcriptome or Gene expression | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| PRJNA852889 | BioProject | RNA | Blood Transcriptome Profiles of COVID-19 patients admitted to the Intensive Care Unit | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| PRJNA861390 | BioProject | RNA | Profiling of lung SARS-CoV-2 and influenza virus infection dissects virus-specific host re | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| PRJNA862565 | BioProject | RNA | Transcriptome of leukocytes from acute COVID-19 patients | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| PRJNA895325 | BioProject | RNA | Transcriptomic analysis of PBMC in long covid | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| PRJNA912696 | BioProject | RNA | Dynamic Activity in cis-Regulatory Elements of Leukocytes Identifies Transcription Factor  | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| PRJNA912895 | BioProject | RNA | Homo sapiens Raw sequence reads | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| PRJNA912903 | BioProject | RNA | Immune dysfunction in COVID-19 Associated Mucormycosis | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| PRJNA924681 | BioProject | RNA | Homo sapiens Transcriptome or Gene expression | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| PRJNA953515 | BioProject | RNA | The specialized pro-resolving lipid mediator Protectin D1 affects macrophages differentiat | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| PRJNA966114 | BioProject | RNA | Transcriptomic profile of peripheral blood neutrophils in COVID-19 and IPF (Homo sapiens R | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| GSE101705 | GEO | RNA | Existing blood transcriptional classifiers accurately discriminate active tuberculosis fro | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| GSE107994 | GEO | RNA | A modular transcriptional signature identifies phenotypic heterogeneity of human tuberculo | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| GSE108497 | GEO | RNA | Blood immunomonitoring uncovers pregnancy-induced immune programs and their dysregulation  | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| GSE111368 | GEO | RNA | Progression of whole blood transcriptional signatures from interferon-induced to neutrophi | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| GSE128078 | GEO | RNA | Whole blood human transcriptome and virome analysis of ME/CFS patients experiencing post-e | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| GSE13015 | GEO | RNA | Genomic Transcriptional Profiling Identifies a Blood Biomarker Signature for the Diagnosis | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| GSE137342 | GEO | RNA | Blood transcriptome of human sepsis in an Indian cohort | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| GSE149437 | GEO | RNA | Whole blood transcriptomics data in preterm birth and early preeclampsia | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| GSE154272 | GEO | RNA | Long-term longitudinal analysis of peripheral blood transcriptome in COVID-19 patients | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| GSE157103 | GEO | RNA | Large-scale Multi-omic Analysis of COVID-19 Severity | needs decision | leukocytes isolated from whole blood (not whole blood or PBMC as the rule states); severity 50 against 50 (ICU); non-COVID patients are other disease |
| GSE158038 | GEO | RNA | In severe COVID-19, SARS-CoV-2 induces a chronic, TGF-β-dominated adaptive immune response | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| GSE161778 | GEO | RNA | Longitudinal multi-omics identifies responses of megakaryocytes, erythroid cells and plasm | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| GSE162086 | GEO | RNA | Low avidity T cell responses to SARS-CoV-2 in unexposed individuals and severe COVID-19 | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| GSE162562 | GEO | RNA | Immune transcriptomes of highly exposed SARS-CoV-2 asymptomatic seropositive versus serone | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| GSE162835 | GEO | RNA | Host transcriptomic profiling of COVID-19 patients with mild, moderate, and severe clinica | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| GSE164381 | GEO | RNA | Distinct B cell repertoires characterize patients with mild and severe COVID-19 | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| GSE165461 | GEO | RNA | The expansion of unconventional CD56dim CD16neg cells is a hallmark of the natural killer  | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| GSE166424 | GEO | RNA | Asymptomatic COVID-19: disease tolerance with efficient anti-viral immunity against SARS-C | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| GSE167000 | GEO | RNA | PolyA RNA-seq from whole blood of Sars-COV2-negative and -positive subjects | needs decision | who the 30 SARS-CoV-2 negative subjects are (healthy or symptomatic patients) is not stated in the record |
| GSE167118 | GEO | RNA | Clonal expansion and activation of tissue-resident memory-like Th17 cells expressing GM-CS | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| GSE168286 | GEO | RNA | Reduced A-to-I editing of endogenous Alu RNAs in severe COVID-19 disease | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| GSE169687 | GEO | RNA | Long-term perturbation of the blood transcriptome months after SARS-CoV-2 infection. | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| GSE172450 | GEO | RNA | Longitudinal profiling of the whole blood transcriptional response of human volunteers to  | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| GSE174083 | GEO | RNA | Large-scale genomic study reveals robust activation of the immune system following advance | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| GSE177025 | GEO | RNA | Effects of steroids on peripheral blood gene expression in COVID-19 patients according to  | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| GSE177477 | GEO | RNA | Type I interferon pathways in SARS-CoV-2-infected individuals from Pakistan | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| GSE178757 | GEO | RNA | Rapid, scalable assessment of SARS-CoV-2 cellular immunity by whole-blood PCR | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| GSE178967 | GEO | RNA | Baseline signatures associated with clinical, virologic, and immunologic outcomes in patie | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| GSE179627 | GEO | RNA | Transcriptome analysis of PBMCs reveals distinct immune response in the asymptomatic and r | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| GSE179850 | GEO | RNA | Gene expression analysis of whole blood obtained from COVID-19 in the Osaka cohort | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| GSE181958 | GEO | RNA | Transcriptomic profiling of host response genes in convalescent COVID-19 donors and health | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| GSE182616 | GEO | RNA | Early Transcriptomic Response to Burn injury: Severe Burns are Associated with Immune Path | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| GSE183071 | GEO | RNA | A multi-tissue study of immune gene expression profiling highlights the key role of the na | needs decision | NanoString panel; blood subset of 54 samples, class counts per participant not read |
| GSE183533 | GEO | RNA | Whole transcriptome sequencing from covid-19 and control lung biopsies | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| GSE184401 | GEO | RNA | Identifying Risk Factors for Secondary Infection Post-SARS-CoV-2 Infection | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| GSE185263 | GEO | RNA | Predicting sepsis severity at first clinical presentation: the role of endotypes and mecha | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| GSE185557 | GEO | RNA | Maternal-Fetal Immune Responses in Pregnant Women Infected with SARS-CoV-2 | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| GSE188847 | GEO | RNA | Severe COVID-19 is associated with molecular signatures of aging in the human brain | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| GSE192902 | GEO | RNA | Maternal plasma cell-free RNA sequencing for preeclamptic and normotensive mothers | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| GSE193022 | GEO | RNA | Blood leukocyte transcriptional modules and differentially expressed genes associated with | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| GSE193708 | GEO | RNA | Protective immune trajectories in early viral containment of non-pneumonic SARS-CoV-2 infe | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| GSE195796 | GEO | RNA | Blood transcriptomes of SARS-CoV-2 infected kidney transplant recipients demonstrate immun | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| GSE196117 | GEO | RNA | Transcriptional profiling of whole blood leukocytes obtained from critically ill patients  | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| GSE196399 | GEO | RNA | Comparative gene and lncRNA expression profiling in patients with severe community acquire | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| GSE196728 | GEO | RNA | The effect of high altitudes and time-of-day on human whole-blood transcriptome | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| GSE196822 | GEO | RNA | Whole blood transcriptome from Covid-19 patients in India | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| GSE196988 | GEO | RNA | Durable Alterations of Hematopoiesis and Chromatin post-COVID-19 [MultiomeRNA] | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| GSE198256 | GEO | RNA | Functional reprogramming of monocytes in acute and convalescent severe COVID-19 patients [ | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| GSE198257 | GEO | RNA | Functional reprogramming of monocytes in acute and convalescent severe COVID-19 patients | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| GSE198864 | GEO | RNA | Human lungs show limited permissiveness for SARS-CoV-2 due to scarce ACE2 levels but virus | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| GSE208587 | GEO | RNA | Risk assessment with gene expression markers in sepsis development | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| GSE211979 | GEO | RNA | Redox Imbalance in COVID-19 Pathophysiology | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| GSE212865 | GEO | RNA | Dynamics of gene expression profiling by microarrays and identification of high-risk patie | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| GSE212866 | GEO | RNA | Dynamics of gene expression profiling by microarrays and identification of high-risk patie | needs decision | SuperSeries of arrays; the negatives are hospitalised patients (other disease); severity only as ARDS against no ARDS |
| GSE214647 | GEO | RNA | Abnormal thrombosis and neutrophil activation increase hospital-acquired sacral pressure i | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| GSE217799 | GEO | RNA | Cardiac-related long non-coding RNA profiling in the peripheral blood mononuclear cells of | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| GSE217948 | GEO | RNA | COVID-19 patient samples collected by the PREDICT-19 consortium | needs decision | PREDICT-19: 'control' (71) is not defined in the record; per sample disease state missing for 193 samples |
| GSE221066 | GEO | RNA | Dynamic Activity in cis-Regulatory Elements of Leukocytes Identifies Transcription Factor  | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| GSE221067 | GEO | RNA | Dynamic Activity in cis-Regulatory Elements of Leukocytes Identifies Transcription Factor  | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| GSE222253 | GEO | RNA | Post-COVID symptoms and endotypes are associated with the inability to modulate the trajec | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| GSE222393 | GEO | RNA | Unresolved Immune Dysfunction is Lethal in Both COVID-19 and non-COVID-19 Sepsis | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| GSE223885 | GEO | RNA | Deficient radiation transcription response in COVID-19 patients | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| GSE226260 | GEO | RNA | Long COVID involves activation of proinflammatory and immune exhaustion pathways | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| GSE226677 | GEO | RNA | Immune dysfunction in COVID-19 Associated Mucormycosis | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| GSE227080 | GEO | RNA | Early differentially expressed immunological genes   in mild and severe COVID-19 | needs decision | NanoString immunology panel (579 genes) on PBMC; the 59 COVID-negative subjects are not described; severity class MOD_SEV merges moderate and severe |
| GSE227116 | GEO | RNA | Gene expression profiles of 65 blood samples recovered from 46 COVID-19 convalescent patie | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| GSE228320 | GEO | RNA | Whole blood transcriptional profiling of pulmonary functional sequelae in ARDS- secondary  | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| GSE231660 | GEO | RNA | Repair of Acute Respiratory Distress Syndrome in COVID-19 by Stromal Cells (REALIST-COVID  | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| GSE232628 | GEO | RNA | Expansion of profibrotic monocyte-derived alveolar macrophages in patients with persistent | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| GSE233626 | GEO | RNA | Reanalysis of peripheral blood gene expression of IFIH1 rs1990760 variants in critically-i | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| GSE233627 | GEO | RNA | Reanalysis of peripheral blood gene expression of IFIH1 rs1990760 variants in critically-i | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| GSE234939 | GEO | RNA | Durable Alterations of Hematopoiesis and Chromatin post-COVID-19 [bulkRNA-Seq] | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| GSE236318 | GEO | RNA | TI caused by SARS-CoV-2 | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| GSE236562 | GEO | RNA | Recapitulation of pathophysiological features of AD in SARS-CoV-2 infected subjects | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| GSE237180 | GEO | RNA | Distinct pulmonary and systemic effects of dexamethasone in severe COVID-19. | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| GSE239450 | GEO | RNA | SARS-CoV-2 infection elucidates unique features of pregnancy-specific immunity [GEX] | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| GSE239452 | GEO | RNA | SARS-CoV-2 infection elucidates unique features of pregnancy-specific immunity | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| GSE239588 | GEO | RNA | Pharmacological inhibition of bromodomain and extra-terminal proteins induces NRF-2-mediat | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| GSE249139 | GEO | RNA | Placental Privilege: Evidence of Organ Resilience in severe COVID-19 | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| GSE270045 | GEO | RNA | Upregulation of olfactory receptors and neuronal-associated genes highlights complex immun | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| GSE27383 | GEO | RNA | Marked Reduction of AKT1 Expression and Deregulation of AKT1-associated Pathways in Periph | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| GSE275334 | GEO | RNA | Immune Exhaustion in ME/CFS and long COVID | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| GSE282464 | GEO | RNA | Molecular landscape of respiratory infection: A large-scale, multi-centre blood transcript | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| GSE28991 | GEO | RNA | Genome-wide gene expression analysis of human whole-blood samples during acute dengue dise | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| GSE29385 | GEO | RNA | Influenzavirus serotype association to global whole blood transcriptional changes | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| GSE297499 | GEO | RNA | CSF immune cell alterations in women with neuropsychiatric Long COVID | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| GSE322560 | GEO | RNA | Durable interferon-linked blood signatures in COVID-19 convalescence: a longitudinal pre/p | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| GSE32707 | GEO | RNA | Inflammasome-regulated cytokines are critical mediators of lung injury in humans. | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| GSE346693 | GEO | RNA | Early enhanced TLR3/MDA5 responsiveness in individuals displaying resistance to SARS-CoV-2 | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| GSE37250 | GEO | RNA | Genome-wide transcriptional profiling of HIV positive and negative adults with active tube | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| GSE40012 | GEO | RNA | Blood transcriptome of human bacterial and influenza A pneumonia | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| GSE46449 | GEO | RNA | Expression data from Patients with Bipolar (BP) Disorder and Matched Control Subjects | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| GSE48000 | GEO | RNA | Whole blood gene expression profiles distinguish clinical phenotypes of venous thromboembo | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| GSE54514 | GEO | RNA | Whole blood transcriptome of survivors and nonsurvivors of sepsis | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| GSE54629 | GEO | RNA | Whole-blood transcriptomic profiling highlights rituximab response in rheumatoid arthritis | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| GSE60244 | GEO | RNA | Transcriptional Profiling is Superior to Procalcitonin to Discriminate Bacterial vs. Viral | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| GSE61821 | GEO | RNA | Patient-based transcriptome-wide analysis identify interferon and ubiquination pathways as | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| GSE6613 | GEO | RNA | Parkinson's disease vs. controls, whole blood | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| GSE68310 | GEO | RNA | Host transcriptional response to influenza and other acute respiratory viral infections –  | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| GSE76262 | GEO | RNA | Expression data of induced sputum from Unbiased BIOmarkers in Prediction of REspiratory Di | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| GSE84076 | GEO | RNA | Host blood trancriptional profiles during Mycobacterium tuberculosis infection. | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| GSE85307 | GEO | RNA | Early Pregnancy Gene Expression Profiling of Women with Preeclampsia: The “Vitamin D Anten | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| GSE86200 | GEO | RNA | The role of vitamin D in the transcriptional program of human pregnancy. | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| GSE90081 | GEO | RNA | Discovery of new candidate genes for rheumatoid arthritis through integration of genetic a | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| GSE93606 | GEO | RNA | Host-Microbial interactions in Idiopathic Pulmonary Fibrosis | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| 10.17632/8wxhhykfnh | Mendeley Data | RNA | Transcriptomic signatures of SARS-CoV-2 and Influenza virus infections | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| 10.17632/m4zdf3mg8c | Mendeley Data | RNA | Transcriptomics of SARS-CoV-2’s coinfection potential: implications for the host targetome | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| 10.17632/t4twwtvv7r | Mendeley Data | RNA | Transcriptomic Similarities and Differences in Host Response between SARS-CoV-2 and Other  | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| 10.17632/tg37fc2s3y | Mendeley Data | RNA | Human coronaviruses in idiopathic Parkinson's disease: Implications of SARS-CoV-2's modula | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| zenodo.15587564 | Zenodo | RNA | AI-driven multimodal algorithm predicts immunotherapy and targeted therapy outcomes in cle | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| zenodo.15722617 | Zenodo | RNA | Tracing the Evolutionary Pathway of SARS-CoV-2 Through RNA Sequencing Analysis | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| zenodo.4283753 | Zenodo | RNA | A diagnostic host response biosignature for COVID-19 from RNA profiling of nasal swabs and | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| zenodo.4757764 | Zenodo | RNA | Preprocessing of Public RNA-sequencing Datasets to Facilitate Downstream Analyses of Human | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| E-MTAB-1548 | other | RNA | A 6-mRNA host response classifier in whole blood predicts outcomes in COVID-19 and other a | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| E-MTAB-5921 | other | RNA | Transcriptomic signatures and repurposing drugs for COVID-19 patients: findings of bioinfo | needs decision | passes the automatic screen; not read by hand; the cap of 10 RNA cohorts was reached, so it was not read further |
| PRJNA789983 | BioProject | metabolites | Rapid Acceleration of Diagnostics - Radical (RADx-rad): COVID-19 Detection through Scent A | needs decision | passes the automatic screen; not read by hand; metabolite cohorts were not read further (the targets of this sweep were 10 RNA and 6 protein cohorts) |
| PRJNA789984 | BioProject | metabolites | Rapid Acceleration of Diagnostics - Radical (RADx-rad): COVID-19 Detection through Scent A | needs decision | passes the automatic screen; not read by hand; metabolite cohorts were not read further (the targets of this sweep were 10 RNA and 6 protein cohorts) |
| 10.17632/3h96n97xrb | Mendeley Data | metabolites | SARS-CoV-2 infection dysregulates the metabolomic and lipidomic profiles of serum. Bruzzon | needs decision | passes the automatic screen; not read by hand; metabolite cohorts were not read further (the targets of this sweep were 10 RNA and 6 protein cohorts) |
| 10.17632/7d63nm2nrm | Mendeley Data | metabolites | IMMUNOMETABOLISM COVID-19 | needs decision | passes the automatic screen; not read by hand; metabolite cohorts were not read further (the targets of this sweep were 10 RNA and 6 protein cohorts) |
| 10.17632/7fnt3nfhdv | Mendeley Data | metabolites | Targeted metabolomics identifies high performing diagnostic and prognostic biomarkers for  | needs decision | tissue not stated in the record; not read by hand; metabolite cohorts were not read further (the targets of this sweep were 10 RNA and 6 protein cohorts) |
| 10.17632/bg45mx8rxd | Mendeley Data | metabolites | A Comparative NMR-based Metabolomics Study of Lung Parenchyma of Severe COVID‑19 Patients | needs decision | passes the automatic screen; not read by hand; metabolite cohorts were not read further (the targets of this sweep were 10 RNA and 6 protein cohorts) |
| 10.17632/x9tw3knwsd | Mendeley Data | metabolites | Immunometabolic Signatures Predict Risk of Progression to Sepsis in COVID-19 | needs decision | passes the automatic screen; not read by hand; metabolite cohorts were not read further (the targets of this sweep were 10 RNA and 6 protein cohorts) |
| ST001709 | Metabolomics Workbench | metabolites | SARS-CoV-2 infection rewires host cell metabolism and is potentially susceptible to mTORC1 | needs decision | tissue not stated in the record; not read by hand; metabolite cohorts were not read further (the targets of this sweep were 10 RNA and 6 protein cohorts) |
| ST001735 | Metabolomics Workbench | metabolites | The COVIDome Explorer Researcher Portal (Red Blood Cells) | needs decision | passes the automatic screen; not read by hand; metabolite cohorts were not read further (the targets of this sweep were 10 RNA and 6 protein cohorts) |
| ST001736 | Metabolomics Workbench | metabolites | The COVIDome Explorer Researcher Portal (blood plasma) | needs decision | passes the automatic screen; not read by hand; metabolite cohorts were not read further (the targets of this sweep were 10 RNA and 6 protein cohorts) |
| ST001849 | Metabolomics Workbench | metabolites | Longitudinal Metabolomics of Human Plasma Reveals Robust Prognostic Markers of COVID-19 Di | needs decision | passes the automatic screen; not read by hand; metabolite cohorts were not read further (the targets of this sweep were 10 RNA and 6 protein cohorts) |
| ST001933 | Metabolomics Workbench | metabolites | Absolute quantification of plasma cytokines and metabolome reveals the glycylproline regul | needs decision | passes the automatic screen; not read by hand; metabolite cohorts were not read further (the targets of this sweep were 10 RNA and 6 protein cohorts) |
| ST001993 | Metabolomics Workbench | metabolites | CE-MS based metabolomics to study plasma samples that reveal new pathways implicated in SA | needs decision | passes the automatic screen; not read by hand; metabolite cohorts were not read further (the targets of this sweep were 10 RNA and 6 protein cohorts) |
| ST002016 | Metabolomics Workbench | metabolites | Metabolomics of COVID patients | needs decision | tissue not stated in the record; not read by hand; metabolite cohorts were not read further (the targets of this sweep were 10 RNA and 6 protein cohorts) |
| ST002087 | Metabolomics Workbench | metabolites | Profiling metabolites and lipoproteins in COMETA, an Italian cohort of COVID-19 patients | needs decision | tissue not stated in the record; not read by hand; metabolite cohorts were not read further (the targets of this sweep were 10 RNA and 6 protein cohorts) |
| ST002194 | Metabolomics Workbench | metabolites | Metabolic profiling at COVID-19 onset shows disease severity and sex-specific dysregulatio | needs decision | tissue not stated in the record; not read by hand; metabolite cohorts were not read further (the targets of this sweep were 10 RNA and 6 protein cohorts) |
| ST002291 | Metabolomics Workbench | metabolites | Integrated metabolic and inflammatory signatures associated with severity, fatality, and r | needs decision | tissue not stated in the record; not read by hand; metabolite cohorts were not read further (the targets of this sweep were 10 RNA and 6 protein cohorts) |
| ST002301 | Metabolomics Workbench | metabolites | Serum metabolomics profiling identifies new predictive biomarkers for disease severity in  | needs decision | passes the automatic screen; not read by hand; metabolite cohorts were not read further (the targets of this sweep were 10 RNA and 6 protein cohorts) |
| ST002404 | Metabolomics Workbench | metabolites | Profiling Metabolites and Lipoproteins in COMETA, an Italian Cohort of COVID-19 patients-p | needs decision | tissue not stated in the record; not read by hand; metabolite cohorts were not read further (the targets of this sweep were 10 RNA and 6 protein cohorts) |
| ST002477 | Metabolomics Workbench | metabolites | Neutrophil metabolomics in COVID-19 | needs decision | tissue not stated in the record; not read by hand; metabolite cohorts were not read further (the targets of this sweep were 10 RNA and 6 protein cohorts) |
| ST002829 | Metabolomics Workbench | metabolites | Nucleotide, phospholipid, and kynurenine metabolites are robustly associated with COVID-19 | needs decision | passes the automatic screen; not read by hand; metabolite cohorts were not read further (the targets of this sweep were 10 RNA and 6 protein cohorts) |
| ST002845 | Metabolomics Workbench | metabolites | Methylprednisolone therapy induces differential metabolic trajectories in severe COVID-19  | needs decision | tissue not stated in the record; not read by hand; metabolite cohorts were not read further (the targets of this sweep were 10 RNA and 6 protein cohorts) |
| ST002984 | Metabolomics Workbench | metabolites | Amino acid catabolite markers for early prognostication of pneumonia in patients with COVI | needs decision | tissue not stated in the record; not read by hand; metabolite cohorts were not read further (the targets of this sweep were 10 RNA and 6 protein cohorts) |
| ST003674 | Metabolomics Workbench | metabolites | Metabolomic and lipoproteomic differences and similarities among different types of pneumo | needs decision | tissue not stated in the record; not read by hand; metabolite cohorts were not read further (the targets of this sweep were 10 RNA and 6 protein cohorts) |
| ST003922 | Metabolomics Workbench | metabolites | Time-series metabolomic profiling of SARS-CoV-2 infection: Possible prognostic biomarkers  | needs decision | tissue not stated in the record; not read by hand; metabolite cohorts were not read further (the targets of this sweep were 10 RNA and 6 protein cohorts) |
| ST004195 | Metabolomics Workbench | metabolites | Untargeted metabolomics profiling of HBEpCs infected with wild-type or deISGylation-defici | needs decision | tissue not stated in the record; not read by hand; metabolite cohorts were not read further (the targets of this sweep were 10 RNA and 6 protein cohorts) |
| zenodo.12737754 | Zenodo | metabolites | Plasma metabolomic signatures for copy number variants and COVID-19 risk loci in Northern  | needs decision | passes the automatic screen; not read by hand; metabolite cohorts were not read further (the targets of this sweep were 10 RNA and 6 protein cohorts) |
| PRJDB18101 | BioProject | proteins | Japan COVID-19 Task Force | needs decision | passes the automatic screen; not read by hand; not read further after the protein candidates with open matrices and labels were read |
| PRJEB43380 | BioProject | proteins | The Trans-omics Landscape of COVID-19 | needs decision | passes the automatic screen; not read by hand; not read further after the protein candidates with open matrices and labels were read |
| PRJNA1291229 | BioProject | proteins | IVIG-treated COVID-19 patients | needs decision | passes the automatic screen; not read by hand; not read further after the protein candidates with open matrices and labels were read |
| PRJNA1303036 | BioProject | proteins | Pentraxin 3 (PTX-3) and severness of COVID | needs decision | passes the automatic screen; not read by hand; not read further after the protein candidates with open matrices and labels were read |
| PRJNA1465860 | BioProject | proteins | Characterizing the Immune Response to COVID-19 Infection | needs decision | passes the automatic screen; not read by hand; not read further after the protein candidates with open matrices and labels were read |
| PRJNA1465861 | BioProject | proteins | Characterizing the Immune Response to COVID-19 Infection | needs decision | passes the automatic screen; not read by hand; not read further after the protein candidates with open matrices and labels were read |
| PRJNA1476232 | BioProject | proteins | Truway Health Long COVID Cohort Host Genomics Study | needs decision | passes the automatic screen; not read by hand; not read further after the protein candidates with open matrices and labels were read |
| PRJNA251404 | BioProject | proteins | Integrative DNA, RNA and protein evidence connects TREML4 to coronary artery calcification | needs decision | passes the automatic screen; not read by hand; not read further after the protein candidates with open matrices and labels were read |
| PRJNA357628 | BioProject | proteins | Expression of long non-coding RNAs in autoimmunity and linkage to enhancer function and au | needs decision | passes the automatic screen; not read by hand; not read further after the protein candidates with open matrices and labels were read |
| PRJNA664368 | BioProject | proteins | Antiretroviral treatment induced decrease in immune activation contributes to reduced susc | needs decision | passes the automatic screen; not read by hand; not read further after the protein candidates with open matrices and labels were read |
| PRJNA882800 | BioProject | proteins | The longitudinal landscaping of inflammatory plasma factors reveals heterogeneous but core | needs decision | passes the automatic screen; not read by hand; not read further after the protein candidates with open matrices and labels were read |
| PRJNA950540 | BioProject | proteins | Multimodal immune profiling to determine mechanisms of COVID-19 clinical trajectory in Uga | needs decision | passes the automatic screen; not read by hand; not read further after the protein candidates with open matrices and labels were read |
| PRJNA950542 | BioProject | proteins | Multimodal immune profiling to determine mechanisms of COVID-19 clinical trajectory in Uga | needs decision | passes the automatic screen; not read by hand; not read further after the protein candidates with open matrices and labels were read |
| GSE107104 | GEO | proteins | Increased expression of immune regulatory mediators independently detects tuberculosis in  | needs decision | passes the automatic screen; not read by hand; not read further after the protein candidates with open matrices and labels were read |
| GSE151513 | GEO | proteins | Experimental and natural evidence of SARS-CoV-2 infection-induced activation of type I int | needs decision | tissue not stated in the record; not read by hand; not read further after the protein candidates with open matrices and labels were read |
| GSE36177 | GEO | proteins | Genome-wide gene expression analysis of the whole blood leukocyte transcriptional response | needs decision | passes the automatic screen; not read by hand; not read further after the protein candidates with open matrices and labels were read |
| GSE65682 | GEO | proteins | Genome-wide blood transcriptional profiling in critically ill patients - MARS consortium | needs decision | passes the automatic screen; not read by hand; not read further after the protein candidates with open matrices and labels were read |
| MSV000085703 | MassIVE | proteins | Decreased plasma cartilage acidic protein 1 in COVID-19. | needs decision | passes the automatic screen; not read by hand; not read further after the protein candidates with open matrices and labels were read |
| MSV000087099 | MassIVE | proteins | Divergent and self-reactive immune responses in the CNS of COVID-19 patients with neurolog | needs decision | tissue not stated in the record; not read by hand; not read further after the protein candidates with open matrices and labels were read |
| MSV000087115 | MassIVE | proteins | Divergent and self-reactive immune responses in the CNS of COVID-19 patients with neurolog | needs decision | tissue not stated in the record; not read by hand; not read further after the protein candidates with open matrices and labels were read |
| MSV000088382 | MassIVE | proteins | Protein Posttranslational Signatures Identified in COVID-19 Patient Plasma. | needs decision | passes the automatic screen; not read by hand; not read further after the protein candidates with open matrices and labels were read |
| MSV000089568 | MassIVE | proteins | Integrated multiomics analysis to infer COVID-19 biological insights. | needs decision | tissue not stated in the record; not read by hand; not read further after the protein candidates with open matrices and labels were read |
| MSV000089614 | MassIVE | proteins | Viral afterlife: SARS-CoV-2 as a reservoir of immunomimetic peptides that reassemble into  | needs decision | tissue not stated in the record; not read by hand; not read further after the protein candidates with open matrices and labels were read |
| MSV000091179 | MassIVE | proteins | Into the Dark Serum Proteome: Personalized Features of IgG1 and IgA1 Repertoires in Severe | needs decision | passes the automatic screen; not read by hand; not read further after the protein candidates with open matrices and labels were read |
| 10.17632/2cbxgsn7vx | Mendeley Data | proteins | Olink proteomics in COVID-19 | needs decision | Osaka Olink NPX (49 acute and 34 later COVID-19 samples, 20 healthy): no participant ids, so one sample per participant cannot be checked; no paper linked |
| 10.17632/5prg4tsz7n | Mendeley Data | proteins | Multiomics personalized network analyses highlight progressive disruption of central metab | needs decision | tissue not stated in the record; not read by hand; not read further after the protein candidates with open matrices and labels were read |
| 10.17632/8dyncfg3y8 | Mendeley Data | proteins | Blood alternative splicing biomarker for infectious disease. Zhang et al. | needs decision | passes the automatic screen; not read by hand; not read further after the protein candidates with open matrices and labels were read |
| 10.17632/k33vwbdfsn | Mendeley Data | proteins | Broad Auto-Reactive IgM Responses are Common in Critically Ill Patients Including Those wi | needs decision | tissue not stated in the record; not read by hand; not read further after the protein candidates with open matrices and labels were read |
| 10.17632/vmr9zrrgky | Mendeley Data | proteins | Supplementary information for "Plasma Proteomics identify Biomarkers and Pathogenesis of C | needs decision | passes the automatic screen; not read by hand; not read further after the protein candidates with open matrices and labels were read |
| PXD018874 | PRIDE | proteins | High-throughput proteomics platform for clinical applications | needs decision | abundance tables deposited (1); per sample labels to take from the paper; not read by hand; not read further after the protein candidates with open matrices and labels were read |
| PXD037682 | PRIDE | proteins | ALSPAC human inflammation plasma proteomics | needs decision | abundance tables deposited (1); per sample labels to take from the paper; not read by hand; not read further after the protein candidates with open matrices and labels were read |
| PXD070765 | PRIDE | proteins | Cross-platform Clinical Proteomics using the Charite Open Standard for Plasma Proteomics ( | needs decision | abundance tables deposited (1); per sample labels to take from the paper; not read by hand; not read further after the protein candidates with open matrices and labels were read |
| PXD018682 | ProteomeXchange | proteins | MS analysis of SARS-CoV2 proteins from patient samples | needs decision | abundance tables deposited (2, SDRF); per sample labels to take from the paper; not read by hand; not read further after the protein candidates with open matrices and labels were read |
| PXD019423 | ProteomeXchange | proteins | MS analysis of SARS-CoV2 proteins from patient samples | needs decision | abundance tables deposited (1, SDRF); per sample labels to take from the paper; not read by hand; not read further after the protein candidates with open matrices and labels were read |
| PXD021581 | ProteomeXchange | proteins | Prognostic accuracy of Mass Spectrometric Analysis of Plasma in COVID-19 | needs decision | abundance tables deposited (1, SDRF); per sample labels to take from the paper; not read by hand; not read further after the protein candidates with open matrices and labels were read |
| PXD025138 | ProteomeXchange | proteins | Prognostic accuracy of MALDI-TOF mass spectrometric analysis of plasma in COVID-19 | needs decision | abundance tables deposited (2, SDRF); per sample labels to take from the paper; not read by hand; not read further after the protein candidates with open matrices and labels were read |
| PXD030860 | ProteomeXchange | proteins | The serum of COVID-19 asymptomatic patients  in circulating angiogenic cells ex-vivo | needs decision | abundance tables deposited (1, SDRF); per sample labels to take from the paper; not read by hand; not read further after the protein candidates with open matrices and labels were read |
| PXD034620 | ProteomeXchange | proteins | Cardiovascular-related proteomic changes in ECFCs exposed to the serum of COVID-19 patient | needs decision | abundance tables deposited (1); per sample labels to take from the paper; not read by hand; not read further after the protein candidates with open matrices and labels were read |
| PXD036020 | ProteomeXchange | proteins | Circulating peptidome is strongly altered in COVID-19 patients | needs decision | abundance tables deposited (1, SDRF); per sample labels to take from the paper; not read by hand; not read further after the protein candidates with open matrices and labels were read |
| PXD036609 | ProteomeXchange | proteins | Human serum LC-MS/MS for long-term persistent COVID-19 patients | needs decision | abundance tables deposited (5, SDRF); per sample labels to take from the paper; not read by hand; not read further after the protein candidates with open matrices and labels were read |
| PXD037400 | ProteomeXchange | proteins | Quantitative Proteomics Analysis of COVID-19 Patients | needs decision | abundance tables deposited (1); per sample labels to take from the paper; not read by hand; not read further after the protein candidates with open matrices and labels were read |
| PXD040205 | ProteomeXchange | proteins | Fast proteomics with dia-PASEF and analytical flow-rate chromatography  | needs decision | abundance tables deposited (7, SDRF); per sample labels to take from the paper; not read by hand; not read further after the protein candidates with open matrices and labels were read |
| PXD047901 | ProteomeXchange | proteins | Deep plasma proteomics with data-independent acquisition: A fastlane towards biomarkers id | needs decision | abundance tables deposited (13); per sample labels to take from the paper; not read by hand; not read further after the protein candidates with open matrices and labels were read |
| PXD051341 | ProteomeXchange | proteins | Proteomic profiling of severe pneumonia fluids | needs decision | abundance tables deposited (4, SDRF); per sample labels to take from the paper; not read by hand; not read further after the protein candidates with open matrices and labels were read |
| PXD066427 | ProteomeXchange | proteins | Proteins profiling of ICU COVID-19 with versus without comorbidities using Proteomics: a c | needs decision | abundance tables deposited (2); per sample labels to take from the paper; not read by hand; not read further after the protein candidates with open matrices and labels were read |
| zenodo.10142777 | Zenodo | proteins | Longitudinal plasma proteomics reveals biomarkers of alveolar-capillary barrier disruption | needs decision | passes the automatic screen; not read by hand; not read further after the protein candidates with open matrices and labels were read |
| zenodo.10469692 | Zenodo | proteins | Longitudinal plasma proteomics reveals biomarkers of alveolar-capillary barrier disruption | needs decision | passes the automatic screen; not read by hand; not read further after the protein candidates with open matrices and labels were read |
| zenodo.13363117 | Zenodo | proteins | Long COVID IRIS Study Olink Proteomics Dataset | needs decision | passes the automatic screen; not read by hand; not read further after the protein candidates with open matrices and labels were read |
| zenodo.13942581 | Zenodo | proteins | COVID-19: Longitudinal immune correlates studies icddr,b (GH-VAP-SI-55) | needs decision | passes the automatic screen; not read by hand; not read further after the protein candidates with open matrices and labels were read |
| zenodo.14687731 | Zenodo | proteins | Metabolomic profile of severe COVID-19 and a signature predictive of progression towards s | needs decision | passes the automatic screen; not read by hand; not read further after the protein candidates with open matrices and labels were read |
| zenodo.14772494 | Zenodo | proteins | Olink data for 'Identification of soluble biomarkers that associate with distinct manifest | needs decision | passes the automatic screen; not read by hand; not read further after the protein candidates with open matrices and labels were read |
| zenodo.17609054 | Zenodo | proteins | Afrah-Khairallah/somascan-host-proteins-sarscov2: Release of analysis scripts | needs decision | tissue not stated in the record; not read by hand; not read further after the protein candidates with open matrices and labels were read |
| zenodo.21640246 | Zenodo | proteins | Olink proteomic dataset from a cohort of patients with PASC | needs decision | passes the automatic screen; not read by hand; not read further after the protein candidates with open matrices and labels were read |
| zenodo.4329382 | Zenodo | proteins | Covid-19 automated diagnosis and risk assessment through Metabolomics and Machine Learning | needs decision | passes the automatic screen; not read by hand; not read further after the protein candidates with open matrices and labels were read |
| zenodo.6659752 | Zenodo | proteins | pietznerm/elf5_covid_19: ELF5_proteomics | needs decision | tissue not stated in the record; not read by hand; not read further after the protein candidates with open matrices and labels were read |
| E-MTAB-10169 | ArrayExpress | RNA | A single-cell atlas of lymphocyte adaptive immune repertoires and transcriptomes reveals a | exclude | single cell data |
| E-MTAB-10299 | ArrayExpress | RNA | miRNA-Seq of human blood with different disease state of COVID-19 | exclude | other assay or analyte |
| E-MTAB-10562 | ArrayExpress | RNA | small RNA sequencing from plasma samples of SARS-CoV-2 infected patients | exclude | other assay or analyte |
| E-MTAB-10636 | ArrayExpress | RNA | Single-cell V(D)J-seq of CD3+ and CD19+ lymphocytes from PBMCs of COVID-19 patients | exclude | single cell data |
| E-MTAB-11315 | ArrayExpress | RNA | NanoString profiling of host response in healthy volunteers vaccinated with ARCT-021 self- | exclude | vaccination study |
| E-MTAB-11845 | ArrayExpress | RNA | Single-cell RNA-sequencing of spike-specific CD4 and CD8 T cells after COVID-19 vaccinatio | exclude | single cell data |
| E-MTAB-11973 | ArrayExpress | RNA | RNAseq analyses in HEK293T-ACE2 cells infected with SARS-CoV-2 (VOC Delta, 24hours) vs unt | exclude | in vitro or cell culture |
| E-MTAB-12651 | ArrayExpress | RNA | Longitudinal scRNAseq profiling of memory B cells in mRNA vaccinated individuals post BA.1 | exclude | single cell data |
| E-MTAB-12716 | ArrayExpress | RNA | scRNA-seq data related to Booster mRNA vaccination post-SARS-CoV-2 infection enhances func | exclude | single cell data |
| E-MTAB-12793 | ArrayExpress | RNA | Diagnosis of multisystem inflammatory syndrome in children (MIS-C) by a whole-blood transc | exclude | not adults |
| E-MTAB-12829 | ArrayExpress | RNA | Microarray analysis of healthcare workers that received the COVID-19 mRNA vaccination | exclude | vaccination study |
| E-MTAB-13544 | ArrayExpress | RNA | Single cell RNA-Seq of post mortem COVID19 cases from Malawi (lung, nasal & blood) | exclude | single cell data |
| E-MTAB-13593 | ArrayExpress | RNA | TCR Rep-Seq of COVID-19 affected individuals | exclude | other assay or analyte |
| E-MTAB-13721 | ArrayExpress | RNA | scRNA-seq, V(D)J seq, CITE-seq, LIBRA-seq of human CD4+, CD8+ and B lymphocytes PBMCs of 5 | exclude | single cell data |
| E-MTAB-13919 | ArrayExpress | RNA | RNA-seq of human Caco-2 cells treated with a nutraceutical formula (solution-3), or vehicl | exclude | in vitro or cell culture |
| E-MTAB-14405 | ArrayExpress | RNA | PROS1 released by lung basal cells limits inflammation in epithelial and monocytes during  | exclude | single cell data |
| E-MTAB-14606 | ArrayExpress | RNA | Peripheral Blood Mononuclear Cell Electrophysiology as a Marker for Immune Response to Exp | exclude | single cell data |
| E-MTAB-15066 | ArrayExpress | RNA | Coverage-Immuno | exclude | vaccination study |
| E-MTAB-15195 | ArrayExpress | RNA | Early gene expression after COVID-19 mRNA (BNT162b2) vaccination that is associated with s | exclude | vaccination study |
| E-MTAB-16531 | ArrayExpress | RNA | Single-cell RNA and B cell receptor profiling of B cell response to SARS-CoV-2 vaccination | exclude | single cell data |
| E-MTAB-16728 | ArrayExpress | RNA | Single-cell RNA-seq of human primary pDCs before and after 24 h stimulation with TLR9 agon | exclude | single cell data |
| E-MTAB-9221 | ArrayExpress | RNA | RNA-seq of RBC depleted whole blood from COVID-19 patients and controls | exclude | single cell data |
| E-MTAB-9652 | ArrayExpress | RNA | Integrated longitudinal immunophenotyping, transcriptional and repertoire analyses depict  | exclude | single cell data |
| PRJDB12990 | BioProject | RNA | Four cases of death from systemic immune response syndrome after COVID-19 vaccination | exclude | not human |
| PRJDB14078 | BioProject | RNA | Human immune and gut microbial parameters associated with inter-individual variations in B | exclude | vaccination study |
| PRJDB14085 | BioProject | RNA | Human immune and gut microbial parameters associated with inter-individual variations in B | exclude | vaccination study |
| PRJDB14382 | BioProject | RNA | single cell RNA-seq of COVID-19 patients (with ARDS, n=8; without ARDS, n=8) and healthy c | exclude | single cell data |
| PRJDB37770 | BioProject | RNA | Single-Cell Transcriptomic Analysis of Acute Encephalopathy/Encephalitis Associated with S | exclude | single cell data |
| PRJEB38935 | BioProject | RNA | RNA-seq of RBC depleted whole blood from COVID-19 patients and controls | exclude | single cell data |
| PRJEB44229 | BioProject | RNA | Single cell RNA-seq and TCR-seq for PBMC from COVID-19 patients with persistent SARS-CoV-2 | exclude | single cell data |
| PRJEB45698 | BioProject | RNA | small RNA sequencing from plasma samples of SARS-CoV-2 infected patients | exclude | other assay or analyte |
| PRJEB57285 | BioProject | RNA | Droplet microfluidic analysis of the human B cell immune repertoire over the course of SAR | exclude | RNA from plasma or serum (cell-free RNA), not whole blood or PBMC |
| PRJEB60435 | BioProject | RNA | Longitudinal scRNAseq profiling of memory B cells in mRNA vaccinated individuals post BA.1 | exclude | single cell data |
| PRJEB60784 | BioProject | RNA | Analysis of blood and nasal epithelial transcriptomes to identify mechanisms associated wi | exclude | vaccination study |
| PRJEB60962 | BioProject | RNA | A Gene-expression Module Identifies Circulating Immune Cells with Enhanced Recruitment to  | exclude | single cell data |
| PRJEB65846 | BioProject | RNA | single cell RNAseq of peripheral CD4+ T cells in the context SARS-CoV-2 infection | exclude | single cell data |
| PRJEB71809 | BioProject | RNA | Single cell RNA-Seq of post mortem COVID19 cases from Malawi (lung, nasal & blood) | exclude | single cell data |
| PRJEB74103 | BioProject | RNA | Transcriptomic profiling of mRNAs, lncRNAs, circRNAs and alternative splicing patterns of  | exclude | other assay or analyte |
| PRJEB74742 | BioProject | RNA | scRNA-seq, V(D)J seq, CITE-seq, LIBRA-seq of human CD4+, CD8+ and B lymphocytes PBMCs of 5 | exclude | single cell data |
| PRJEB84380 | BioProject | RNA | Deciphering COVID-19 Host-Pathogen Interactions: Epitranscriptomic Insights | exclude | other assay or analyte |
| PRJEB90224 | BioProject | RNA | Immune Repertoire Dynamics in Vaccinated Kidney Transplant Recipients | exclude | vaccination study |
| PRJNA1000970 | BioProject | RNA | Systems biological assessment of the temporal dynamics of immunity to a viral infection in | exclude | duplicate: BioProject of GSE239787 |
| PRJNA1006282 | BioProject | RNA | Developing enhanced immunotherapy using NKG2A knockout human pluripotent stem cell-derived | exclude | duplicate: BioProject of GSE241039 |
| PRJNA1014418 | BioProject | RNA | Differential gene expressions of various interferon, antiviral response genes as well as t | exclude | duplicate: BioProject of GSE242722 |
| PRJNA1017409 | BioProject | RNA | Integrative biology of patients acutely infected with yellow fever virus | exclude | duplicate: BioProject of GSE243442 |
| PRJNA1017417 | BioProject | RNA | Comparative Analysis of Host Responses Between Sepsis and COVID-19: A Prospective Observat | exclude | duplicate: BioProject of GSE243217 |
| PRJNA1017418 | BioProject | RNA | Comparative Analysis of Host Responses Between Sepsis and COVID-19: A Prospective Observat | exclude | other assay or analyte |
| PRJNA1017864 | BioProject | RNA | Longitudinal gene expression profiling of self-collected blood samples in COVID-19+ and he | exclude | duplicate: BioProject of GSE243348 |
| PRJNA1022135 | BioProject | RNA | Memory B cells dominate the early antibody-secreting cell response to  SARS-CoV-2 mRNA vac | exclude | duplicate: BioProject of GSE244297 |
| PRJNA1029158 | BioProject | RNA | Immunomodulatory Potential of Dapphitkhai and Chanthaleela Thai Traditional Recipes in the | exclude | duplicate: BioProject of GSE245616 |
| PRJNA1033573 | BioProject | RNA | Systems biology analysis reveals distinct molecular signatures promoting immune responsive | exclude | duplicate: BioProject of GSE246525 |
| PRJNA1035131 | BioProject | RNA | Systems biology analysis reveals distinct molecular signatures promoting immune responsive | exclude | duplicate: BioProject of GSE246937 |
| PRJNA1037486 | BioProject | RNA | Blood transcriptome profiling reveals distinct gene networks induced by mRNA vaccination a | exclude | duplicate: BioProject of GSE247401 |
| PRJNA1040889 | BioProject | RNA | mRNA COVID-19 vaccine elicits potent adaptive immune response without the acute inflammati | exclude | duplicate: BioProject of GSE247901 |
| PRJNA1040890 | BioProject | RNA | mRNA COVID-19 vaccine elicits potent adaptive immune response without the acute inflammati | exclude | duplicate: BioProject of GSE247902 |
| PRJNA1040891 | BioProject | RNA | mRNA COVID-19 vaccine elicits potent adaptive immune response without the acute inflammati | exclude | duplicate: BioProject of GSE247904 |
| PRJNA1040892 | BioProject | RNA | mRNA COVID-19 vaccine elicits potent adaptive immune response without the acute inflammati | exclude | duplicate: BioProject of GSE247903 |
| PRJNA1040893 | BioProject | RNA | mRNA COVID-19 vaccine elicits potent adaptive immune response without the acute inflammati | exclude | duplicate: BioProject of GSE247906 |
| PRJNA1040894 | BioProject | RNA | mRNA COVID-19 vaccine elicits potent adaptive immune response without the acute inflammati | exclude | duplicate: BioProject of GSE247905 |
| PRJNA1040895 | BioProject | RNA | mRNA COVID-19 vaccine elicits potent adaptive immune response without the acute inflammati | exclude | duplicate: BioProject of GSE247909 |
| PRJNA1040896 | BioProject | RNA | mRNA COVID-19 vaccine elicits potent adaptive immune response without the acute inflammati | exclude | duplicate: BioProject of GSE247908 |
| PRJNA1040897 | BioProject | RNA | mRNA COVID-19 vaccine elicits potent adaptive immune response without the acute inflammati | exclude | duplicate: BioProject of GSE247907 |
| PRJNA1040898 | BioProject | RNA | mRNA COVID-19 vaccine elicits potent adaptive immune response without the acute inflammati | exclude | duplicate: BioProject of GSE247913 |
| PRJNA1040899 | BioProject | RNA | mRNA COVID-19 vaccine elicits potent adaptive immune response without the acute inflammati | exclude | duplicate: BioProject of GSE247912 |
| PRJNA1040900 | BioProject | RNA | mRNA COVID-19 vaccine elicits potent adaptive immune response without the acute inflammati | exclude | duplicate: BioProject of GSE247911 |
| PRJNA1040901 | BioProject | RNA | mRNA COVID-19 vaccine elicits potent adaptive immune response without the acute inflammati | exclude | duplicate: BioProject of GSE247910 |
| PRJNA1040902 | BioProject | RNA | mRNA COVID-19 vaccine elicits potent adaptive immune response without the acute inflammati | exclude | duplicate: BioProject of GSE247916 |
| PRJNA1040903 | BioProject | RNA | mRNA COVID-19 vaccine elicits potent adaptive immune response without the acute inflammati | exclude | duplicate: BioProject of GSE247915 |
| PRJNA1040904 | BioProject | RNA | mRNA COVID-19 vaccine elicits potent adaptive immune response without the acute inflammati | exclude | duplicate: BioProject of GSE247914 |
| PRJNA1044775 | BioProject | RNA | Single-cell transcriptome of PBMC response to Omicron breakthrough infection | exclude | duplicate: BioProject of GSE248556 |
| PRJNA1047537 | BioProject | RNA | Placental Privilege: Evidence of Organ Resilience in severe COVID-19 | exclude | duplicate: BioProject of GSE249139 |
| PRJNA1047643 | BioProject | RNA | A new humanized mouse mounting mature class-switched, hypermutated and neutralizing antibo | exclude | not human |
| PRJNA1047676 | BioProject | RNA | Antigen-specific CD4+ T cells exhibit distinct transcriptional phenotypes in the lymph nod | exclude | duplicate: BioProject of GSE249313 |
| PRJNA1049350 | BioProject | RNA | Identification of a human hematopoietic stem cell subset that retains memory of inflammato | exclude | duplicate: BioProject of GSE249474 |
| PRJNA1049356 | BioProject | RNA | Identification of a human hematopoietic stem cell subset that retains memory of inflammato | exclude | duplicate: BioProject of GSE249476 |
| PRJNA1049358 | BioProject | RNA | Identification of a human hematopoietic stem cell subset that retains memory of inflammato | exclude | duplicate: BioProject of GSE249478 |
| PRJNA1051620 | BioProject | RNA | Impaired innate and adaptive immune responses to BNT162b2 SARS-CoV-2 vaccination in system | exclude | duplicate: BioProject of GSE250023 |
| PRJNA1051621 | BioProject | RNA | Impaired innate and adaptive immune responses to BNT162b2 SARS-CoV-2 vaccination in system | exclude | duplicate: BioProject of GSE250024 |
| PRJNA1055552 | BioProject | RNA | Blood–brain barrier disruption and sustained systemic inflammation in individuals with lon | exclude | duplicate: BioProject of GSE251849 |
| PRJNA1057237 | BioProject | RNA | SARS-CoV-2 infection evokes transcriptional changes in Calu-3 cells and cultured primary m | exclude | duplicate: BioProject of GSE252056 |
| PRJNA1061045 | BioProject | RNA | Frequency, kinetics and determinants of viable SARS-CoV-2 in bioaerosols from ambulatory C | exclude | duplicate: BioProject of GSE252508 |
| PRJNA1067903 | BioProject | RNA | Comparable CD8 + T cell responses to SARS-CoV-2 vaccination in single-cell transcriptomics | exclude | single cell data |
| PRJNA1068943 | BioProject | RNA | TGFβ links EBV to multisystem inflammatory syndrome in children | exclude | duplicate: BioProject of GSE254179 |
| PRJNA1069682 | BioProject | RNA | SARS-CoV-2 causes dysfunction in human iPSC-derived brain microvascular endothelial like c | exclude | duplicate: BioProject of GSE254309 |
| PRJNA1074394 | BioProject | RNA | AS03 adjuvant enhances the magnitude, persistence, and clonal breadth of memory B cell res | exclude | vaccination study |
| PRJNA1080455 | BioProject | RNA | The Transcriptomic Analysis of Profibrinolytic (PLAUR) and Fibrinolytic Inhibitor (SERPINE | exclude | duplicate: BioProject of GSE259214 |
| PRJNA1080569 | BioProject | RNA | Prolonged exposure to lung-derived cytokines is associated with activation of microglia in | exclude | duplicate: BioProject of GSE259276 |
| PRJNA1097598 | BioProject | RNA | Interferon-responsive neutrophils and macrophages extricate SARS-CoV-2 Omicron critical pa | exclude | duplicate: BioProject of GSE263469 |
| PRJNA1099287 | BioProject | RNA | Comparative single-cell analysis reveals IFN-g as a driver of respiratory sequelae post CO | exclude | duplicate: BioProject of GSE263817 |
| PRJNA1104143 | BioProject | RNA | Single cell gene expression profile of peripheral whole blood cells from long COVID patien | exclude | duplicate: BioProject of GSE265753 |
| PRJNA1105077 | BioProject | RNA | Transcriptomic analysis after SARS-CoV-2 mRNA vaccination reveals a specific gene signatur | exclude | duplicate: BioProject of GSE265975 |
| PRJNA1107496 | BioProject | RNA | Hematopoietic Stem Cell Transplant and Healthy Donor Single cell RNA-seq (scRNA-seq) combi | exclude | duplicate: BioProject of GSE266519 |
| PRJNA1111463 | BioProject | RNA | Cytokine-induced memory-like responses in endothelial cells link chronic inflammation to v | exclude | duplicate: BioProject of GSE267418 |
| PRJNA1112027 | BioProject | RNA | Whole blood transcriptome in long-COVID patients reveals association with lung function an | exclude | duplicate: BioProject of GSE267625 |
| PRJNA1115213 | BioProject | RNA | Extracellular vesicles alter trophoblast function in pregnancies complicated by COVID-19 | exclude | duplicate: BioProject of GSE268230 |
| PRJNA1117966 | BioProject | RNA | Blood DNA methylation in post-acute sequelae of COVID-19 (PASC): a prospective cohort stud | exclude | other assay or analyte |
| PRJNA1117967 | BioProject | RNA | Blood DNA methylation in post-acute sequelae of COVID-19 (PASC): a prospective cohort stud | exclude | other assay or analyte |
| PRJNA1123677 | BioProject | RNA | Differential blood transcriptional signatures induced by SARS-CoV-2 Variants Alpha and Del | exclude | duplicate: BioProject of GSE269780 |
| PRJNA1124166 | BioProject | RNA | Long COVID cohort single-cell RNA sequencing | exclude | duplicate: BioProject of GSE270479 |
| PRJNA1124905 | BioProject | RNA | Upregulation of olfactory receptors and neuronal-associated genes highlights complex immun | exclude | duplicate: BioProject of GSE270045 |
| PRJNA1127560 | BioProject | RNA | Maturing Neutrophils of Lower Density Associate with Thrombocytopenia in Puumala Orthohant | exclude | duplicate: BioProject of GSE270609 |
| PRJNA1128725 | BioProject | RNA | The life-saving benefit of dexamethasone in severe COVID-19 is linked to a reversal of mon | exclude | single cell data |
| PRJNA1136514 | BioProject | RNA | Long-term efficacy of carbon dioxide laser in patients with Hailey-Hailey disease. Assessm | exclude | duplicate: BioProject of GSE272359 |
| PRJNA1136544 | BioProject | RNA | Transcriptomic Profiling of Neutrophils and Low-Density Granulocytes in COVID-19 Patients | exclude | duplicate: BioProject of GSE272381 |
| PRJNA1137881 | BioProject | RNA | Features of chronic urticaria after COVID-19 mRNA vaccine, a real-life cohort study | exclude | duplicate: BioProject of GSE272645 |
| PRJNA1140233 | BioProject | RNA | Identification of immune signatures associated with both SARS-CoV-2 infection and lung tra | exclude | duplicate: BioProject of GSE273137 |
| PRJNA1140250 | BioProject | RNA | Longitudinal transcriptomic analysis reveals persistent enrichment of iron homeostasis and | exclude | duplicate: BioProject of GSE273149 |
| PRJNA1150265 | BioProject | RNA | Immune Exhaustion in ME/CFS and long COVID | exclude | duplicate: BioProject of GSE275334 |
| PRJNA1152718 | BioProject | RNA | Adaptive reprogramming of NK cells during recovery from severe COVID-19 | exclude | single cell data |
| PRJNA1164162 | BioProject | RNA | Remodeling of the chromatin landscape in peripheral  blood cells in patients with severe D | exclude | other assay or analyte |
| PRJNA1175747 | BioProject | RNA | Monocytes in Macrophage Activation Syndrome | exclude | single cell data |
| PRJNA1175748 | BioProject | RNA | Monocytes in Macrophage Activation Syndrome | exclude | single cell data |
| PRJNA1182788 | BioProject | RNA | Widespread gene-environment interactions shape the immune response to SARS-CoV-2 infection | exclude | duplicate: BioProject of GSE281237 |
| PRJNA1186068 | BioProject | RNA | Detailed delineation of the immune response to COVID-19 mRNA vaccination via high-temporal | exclude | duplicate: BioProject of GSE281864 |
| PRJNA1190389 | BioProject | RNA | Remodeling of the chromatin landscape in peripheral blood cells in patients with severe De | exclude | other assay or analyte |
| PRJNA1195123 | BioProject | RNA | Comparative Single-Cell Analyses in Infants Revealed RSV-Specific Declines in NK Cell Resp | exclude | duplicate: BioProject of GSE283746 |
| PRJNA1195124 | BioProject | RNA | Comparative Single-Cell Analyses in Infants Revealed RSV-Specific Declines in NK Cell Resp | exclude | single cell data |
| PRJNA1199065 | BioProject | RNA | Single Cell Profiling of Transcriptional and Epigenetic Differences Between Virus and Seve | exclude | single cell data |
| PRJNA1199066 | BioProject | RNA | Single Cell Profiling of Transcriptional and Epigenetic Differences Between Virus and Seve | exclude | single cell data |
| PRJNA1208100 | BioProject | RNA | Transcriptomic alterations in PBMCs of Long COVID patients | exclude | single cell data |
| PRJNA1208929 | BioProject | RNA | Transcriptomes of classical monocytes in post-COVID patients suggest predisposing conditio | exclude | duplicate: BioProject of GSE286325 |
| PRJNA1209384 | BioProject | RNA | Haplotype-resolved full-length transcriptome analysis in single cells by scCirRL-seq | exclude | single cell data |
| PRJNA1214738 | BioProject | RNA | Temporal profiling of human lymphoid tissues reveals coordinated defence against viral cha | exclude | duplicate: BioProject of GSE287808 |
| PRJNA1214844 | BioProject | RNA | Human naïve-like T follicular regulatory cells are primed for differentiation into mature  | exclude | duplicate: BioProject of GSE287876 |
| PRJNA1214845 | BioProject | RNA | Human naïve-like T follicular regulatory cells are primed for differentiation into mature  | exclude | duplicate: BioProject of GSE287877 |
| PRJNA1216985 | BioProject | RNA | Multi-Omics Analysis of the Immunological State at Baseline and Following SARS-CoV-2 mRNA  | exclude | duplicate: BioProject of GSE288346 |
| PRJNA1218445 | BioProject | RNA | Multi-Omics Analysis of the Immunological State at Baseline and Following SARS-CoV-2 mRNA  | exclude | duplicate: BioProject of GSE288562 |
| PRJNA1218672 | BioProject | RNA | Atherosclerosis licenses for an exceeding immune response in COVID-19 disease by interfero | exclude | single cell data |
| PRJNA1223889 | BioProject | RNA | Maternal acute SARS-CoV-2 infection impairs preimplantation embryo development and reprogr | exclude | duplicate: BioProject of GSE289648 |
| PRJNA1226538 | BioProject | RNA | Minimal impact of ivermectin on immune response and transcriptional profiles in naïve adul | exclude | duplicate: BioProject of GSE290188 |
| PRJNA1251981 | BioProject | RNA | Differential response of human plasmacytoid pre-dendritic cells to SARS-CoV-2 variants | exclude | duplicate: BioProject of GSE294888 |
| PRJNA1256396 | BioProject | RNA | Immune factors produced by PBMCs upon stimulation with Lactobacillus delbrueckii ssp. bulg | exclude | duplicate: BioProject of GSE295798 |
| PRJNA1256938 | BioProject | RNA | High-resolution single cell profiling of SARS-CoV-2 vaccine responsiveness in healthy adul | exclude | single cell data |
| PRJNA1256939 | BioProject | RNA | High-resolution single cell profiling of SARS-CoV-2 vaccine responsiveness in healthy adul | exclude | single cell data |
| PRJNA1261088 | BioProject | RNA | Transcriptome and repertoire analysis of PBMC and CSF cells of a patient with cerebral dem | exclude | duplicate: BioProject of GSE296814 |
| PRJNA1262383 | BioProject | RNA | Biological modifications of the immune response to COVID-19 vaccine in patients treated wi | exclude | vaccination study |
| PRJNA1265093 | BioProject | RNA | CSF immune cell alterations in women with neuropsychiatric Long COVID | exclude | duplicate: BioProject of GSE297499 |
| PRJNA1265151 | BioProject | RNA | Discovery of Key Regulators in Classical Monocyte Phenotypes Linked to COVID-19 Severity u | exclude | duplicate: BioProject of GSE297529 |
| PRJNA1270988 | BioProject | RNA | BAFF-R Expression as a Predictive Biomarker for COVID-19 Vaccine Non-Responsiveness in Ant | exclude | duplicate: BioProject of GSE298709 |
| PRJNA1277703 | BioProject | RNA | Single-cell omics analysis of age-related immune signatures induced by Omicron sublineage  | exclude | duplicate: BioProject of GSE299953 |
| PRJNA1279990 | BioProject | RNA | Single-cell Transcriptomics Reveal Signatures of Altered B cell Response and T cell Senesc | exclude | duplicate: BioProject of GSE300265 |
| PRJNA1289130 | BioProject | RNA | Human brain organoids: a new model to study Cryptococcus neoformans neurotropism | exclude | duplicate: BioProject of GSE302166 |
| PRJNA1291601 | BioProject | RNA | Plasma of COVID-19 patients and healthy subjects | exclude | other assay or analyte |
| PRJNA1293800 | BioProject | RNA | Single cell RNA-Seq of Healthy, COVID-19 Infected, recovered, and longtudinal recovered sa | exclude | duplicate: BioProject of GSE303208 |
| PRJNA1301831 | BioProject | RNA | MAIT and other innate-like T cells integrate adaptive immune responses to modulate interva | exclude | duplicate: BioProject of GSE304583 |
| PRJNA1301834 | BioProject | RNA | MAIT and other innate-like T cells integrate adaptive immune responses to modulate interva | exclude | duplicate: BioProject of GSE304584 |
| PRJNA1303183 | BioProject | RNA | Inflammasome activation and oxidative stress in SARS-CoV-2 breakthrough infection and symp | exclude | duplicate: BioProject of GSE304891 |
| PRJNA1334201 | BioProject | RNA | Enhanced interferon-stimulated gene signatures in B cells following influenza B virus infe | exclude | duplicate: BioProject of GSE309140 |
| PRJNA1364758 | BioProject | RNA | Interrogating antiviral antibody responses with multiplexed, high-throughput serum assays | exclude | in vitro or cell culture |
| PRJNA1368911 | BioProject | RNA | UBA1-CDK16: A Female-Specific Chimeric RNA Emerging Through Evolution and Involved in Immu | exclude | duplicate: BioProject of GSE311311 |
| PRJNA1372037 | BioProject | RNA | Symptomatic SARS-CoV-2 breakthrough infections broaden the repertoire of Spike-reactive CD | exclude | duplicate: BioProject of GSE311993 |
| PRJNA1375552 | BioProject | RNA | Post-COVID impairment of memory T cell responses to community-acquired pathogens can be re | exclude | duplicate: BioProject of GSE312633 |
| PRJNA1391649 | BioProject | RNA | Prothrombotic antibody origin in COVID-19 | exclude | single cell data |
| PRJNA1415219 | BioProject | RNA | Identification of a human hematopoietic stem cell subset that retains memory of inflammato | exclude | single cell data |
| PRJNA1416361 | BioProject | RNA | Integrative analysis reveals ciliary dysfunction and compartmentalized LCN2 expression as  | exclude | duplicate: BioProject of GSE318035 |
| PRJNA1416376 | BioProject | RNA | Identification of a human hematopoietic stem cell subset that retains memory of inflammato | exclude | duplicate: BioProject of GSE318050 |
| PRJNA1428738 | BioProject | RNA | Dysregulated NK-cell gene expression defines the enduring symptoms of long COVID-19 | exclude | duplicate: BioProject of GSE320507 |
| PRJNA1430423 | BioProject | RNA | Durable interferon-linked blood signatures in COVID-19 convalescence: a longitudinal pre/p | exclude | duplicate: BioProject of GSE322560 |
| PRJNA1445095 | BioProject | RNA | Integrative analysis reveals ciliary dysfunction and compartmentalized LCN2 expression as  | exclude | duplicate: BioProject of GSE326313 |
| PRJNA1454901 | BioProject | RNA | Neutrophil extracellular traps aggravate COVID-19 related acute macular neuroretinopathy b | exclude | duplicate: BioProject of GSE328399 |
| PRJNA1458120 | BioProject | RNA | Zanubrutinib monotherapy for IgG4-related head and neck disease | exclude | duplicate: BioProject of GSE329190 |
| PRJNA1462287 | BioProject | RNA | Pathogenic IgG sialylation drives macrophage immunometabolic dysfunction and lung immunopa | exclude | duplicate: BioProject of GSE329962 |
| PRJNA1465100 | BioProject | RNA | Single-cell transcriptomic and immune repertoire profiling of COVID-19 mRNA vaccination re | exclude | single cell data |
| PRJNA1470541 | BioProject | RNA | Baseline interferon signaling in monocytes and antibody-mediated innate activation at inje | exclude | duplicate: BioProject of GSE333235 |
| PRJNA1473281 | BioProject | RNA | Baseline interferon signaling in monocytes and antibody-mediated innate activation at inje | exclude | not human |
| PRJNA1474377 | BioProject | RNA | A pilot study for treatment of chronic cough after COVID-19 infection by nebulization ther | exclude | duplicate: BioProject of GSE334259 |
| PRJNA1526215 | BioProject | RNA | Early enhanced TLR3/MDA5 responsiveness in individuals displaying resistance to SARS-CoV-2 | exclude | duplicate: BioProject of GSE346693 |
| PRJNA315611 | BioProject | RNA | A blood RNA signature for tuberculosis disease risk: a prospective cohort study | exclude | duplicate: BioProject of GSE79362 |
| PRJNA327986 | BioProject | RNA | Host blood trancriptional profiles during Mycobacterium tuberculosis infection. | exclude | duplicate: BioProject of GSE84076 |
| PRJNA352396 | BioProject | RNA | Delineating survival from a fatal outcome in Ebola virus disease in humans | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| PRJNA354367 | BioProject | RNA | Discovery of new candidate genes for rheumatoid arthritis through integration of genetic a | exclude | duplicate: BioProject of GSE90081 |
| PRJNA358580 | BioProject | RNA | RNA sequencing data of whole blood cells of normal glucose tolerant (NGT) and gestational  | exclude | other assay or analyte |
| PRJNA378794 | BioProject | RNA | Genome Scale Analysis of miRNA and mRNA regulation during preterm labor [whole blood] | exclude | other assay or analyte |
| PRJNA384259 | BioProject | RNA | A third generation therapeutic vaccine for human visceral leishmaniasis and post kala azar | exclude | vaccination study |
| PRJNA390194 | BioProject | RNA | RNA-Seq of COPD patients' BAL | exclude | not blood |
| PRJNA390289 | BioProject | RNA | Whole blood RNA-seq of pediatric cases of natural Chikungunya infection | exclude | duplicate: BioProject of GSE99992 |
| PRJNA398240 | BioProject | RNA | Metformin alters human host responses to Mycobacterium tuberculosis in-vitro and in health | exclude | in vitro or cell culture |
| PRJNA401870 | BioProject | RNA | Transcriptional responses to RTS,S standard regimen (RRR) and Ad35.CS+R,R vaccination in c | exclude | vaccination study |
| PRJNA418179 | BioProject | RNA | Progenitor identification and SARS-CoV-2 infection in human distal lung organoids | exclude | single cell data |
| PRJNA427575 | BioProject | RNA | Isolated Iliac Cryptococcosis in an Immunocompetent Patient | exclude | not human |
| PRJNA437114 | BioProject | RNA | Tuberculous Meningitis in Pediatric Patients | exclude | not adults |
| PRJNA484801 | BioProject | RNA | ATAC-seq data | exclude | duplicate: BioProject of GSE118189 |
| PRJNA492956 | BioProject | RNA | A peripheral blood gene expression signature diagnoses subclinical acute rejection [discov | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| PRJNA504827 | BioProject | RNA | Tuberculous Meningitis in Pediatric Patients [CSF] | exclude | not adults |
| PRJNA507472 | BioProject | RNA | Subjects acutely Infected with Chikungunya Virus | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| PRJNA521868 | BioProject | RNA | Double-stranded Alu RNAs, type 1 IFN responses, and relapsing remitting multiple sclerosis | exclude | duplicate: BioProject of GSE126427 |
| PRJNA526259 | BioProject | RNA | Whole blood human transcriptome and virome analysis of ME/CFS patients experiencing post-e | exclude | duplicate: BioProject of GSE128078 |
| PRJNA533086 | BioProject | RNA | RNA-seq transcript and gene data on zika exposed and zika naïve samples | exclude | duplicate: BioProject of GSE129882 |
| PRJNA588242 | BioProject | RNA | Transcriptomic Profiles of Confirmed Pediatric Tuberculosis Patients and Household Contact | exclude | not adults |
| PRJNA591657 | BioProject | RNA | dual RNA-seq of human and Plasmodium falciparum from early childhood Malian blood samples | exclude | not adults |
| PRJNA597217 | BioProject | RNA | ScRNA-seq analysis and functional characterization of myeloid lincRNAs differentially expr | exclude | duplicate: BioProject of GSE142503 |
| PRJNA600846 | BioProject | RNA | RNA sequencing data of peripheral blood isolated from Korean patients with Crohn's disease | exclude | duplicate: BioProject of GSE143507 |
| PRJNA601661 | BioProject | RNA | Shear Stresses Associated with Cardiopulmonary Bypass Activate Expression Inflammatory Cyt | exclude | in vitro or cell culture |
| PRJNA605907 | BioProject | RNA | 2019 novel coronavirus | exclude | not human (Severe acute respiratory syndrome coronavirus 2) |
| PRJNA608900 | BioProject | RNA | Transcriptome profile of leukocyte RNA from children with steroid sensitive and resistant  | exclude | not adults |
| PRJNA622298 | BioProject | RNA | Single-cell RNA-seq of engineered healthy human kidney organoids | exclude | duplicate: BioProject of GSE147863 |
| PRJNA629752 | BioProject | RNA | Immunophenotyping of COVID-19 and Influenza Underscores the Association of Type I IFN Resp | exclude | duplicate: BioProject of GSE149689 |
| PRJNA630932 | BioProject | RNA | Single-cell gene expression sequencing of PBMCs from COVID-19-infected, influenza-infected | exclude | single cell data |
| PRJNA633393 | BioProject | RNA | A single-cell atlas of the peripheral immune response to severe COVID-19 | exclude | duplicate: BioProject of GSE150728 |
| PRJNA634090 | BioProject | RNA | Single-cell RNA-seq of PBMCs from severe COVID-19 patients | exclude | duplicate: BioProject of GSE150861 |
| PRJNA634938 | BioProject | RNA | Blocking of the CD80/86 axis as a therapeutic approach to prevent progression to more seve | exclude | duplicate: BioProject of GSE151161 |
| PRJNA635580 | BioProject | RNA | A shift towards an emergency myelopoiesis in critical COVID-19 patients | exclude | single cell data |
| PRJNA637911 | BioProject | RNA | Multiplexed single-cell transcriptomic analysis of normal and impaired lung development in | exclude | not human (Mus musculus) |
| PRJNA638653 | BioProject | RNA | Tuberculosis in Malnourished Individuals | exclude | not adults |
| PRJNA639275 | BioProject | RNA | Systems biological assessment of immunity to severe and mild COVID-19 infections | exclude | duplicate: BioProject of GSE152418 |
| PRJNA639278 | BioProject | RNA | Episignatures stratifying ADNP syndrome show modest correlation with phenotype [RNA-seq] | exclude | other assay or analyte |
| PRJNA639893 | BioProject | RNA | Transcriptomic Similarities and Differences in Host Response between SARS-CoV-2 and Other  | exclude | duplicate: BioProject of GSE152641 |
| PRJNA643420 | BioProject | RNA | Immunologic characterization of COVID-19 patients with hematologic cancer: biologic and cl | exclude | duplicate: BioProject of GSE153610 |
| PRJNA644579 | BioProject | RNA | Severely ill COVID-19 patients display impaired exhaustion features in SARS-CoV-2-reactive | exclude | duplicate: BioProject of GSE153931 |
| PRJNA645700 | BioProject | RNA | Long-term longitudinal analysis of peripheral blood transcriptome in COVID-19 patients | exclude | duplicate: BioProject of GSE154272 |
| PRJNA646704 | BioProject | RNA | Single cell RNAseq reveal covid19 patient immune response | exclude | duplicate: BioProject of GSE154567 |
| PRJNA647880 | BioProject | RNA | RNAseq analysis of blood from sepsis patients and healthy controls | exclude | duplicate: BioProject of GSE154918 |
| PRJNA648425 | BioProject | RNA | Human dendritic cell and macrophage response to SARS-CoV-2 infection | exclude | duplicate: BioProject of GSE155106 |
| PRJNA649786 | BioProject | RNA | Whole blood transcriptome analysis reveals SARS-CoV-2 ORF8 as a potential therapeutic and  | exclude | duplicate: BioProject of GSE155454 |
| PRJNA661032 | BioProject | RNA | Deciphering the state of immune silence in fatal COVID-19 patients | exclude | duplicate: BioProject of GSE157344 |
| PRJNA661222 | BioProject | RNA | Next Generation Sequencing Quantitative Analysis on COVID-19 Kidney tissue | exclude | duplicate: BioProject of GSE157403 |
| PRJNA662985 | BioProject | RNA | Whole transcriptome sequencing of peripheral blood mononuclear cells from patients with CO | exclude | duplicate: BioProject of GSE157859 |
| PRJNA663837 | BioProject | RNA | In severe COVID-19, SARS-CoV-2 induces a chronic, TGF-β-dominated adaptive immune response | exclude | duplicate: BioProject of GSE158030 |
| PRJNA663838 | BioProject | RNA | In severe COVID-19, SARS-CoV-2 induces a chronic, TGF-β-dominated adaptive immune response | exclude | duplicate: BioProject of GSE158034 |
| PRJNA663839 | BioProject | RNA | In severe COVID-19, SARS-CoV-2 induces a chronic, TGF-β-dominated adaptive immune response | exclude | duplicate: BioProject of GSE158036 |
| PRJNA663865 | BioProject | RNA | Large-scale single-cell analysis reveals critical immune characteristics of COVID-19 patie | exclude | duplicate: BioProject of GSE158055 |
| PRJNA663975 | BioProject | RNA | Human primary brain microvascular endothelial cells transcriptome with/without SARS-CoV-2  | exclude | not blood |
| PRJNA666884 | BioProject | RNA | Early Predicting COVID-19 Severity Value with Extracellular Vesicles and Extracellular RNA | exclude | other assay or analyte |
| PRJNA670179 | BioProject | RNA | Hydroxychloroquine inhibits trained immunity | exclude | duplicate: BioProject of GSE159678 |
| PRJNA672981 | BioProject | RNA | Macrophage expression and prognostic significance of the long pentraxin PTX3 in COVID-19 | exclude | duplicate: BioProject of GSE160351 |
| PRJNA679264 | BioProject | RNA | Dysregulated transcriptional responses to SARS-CoV-2 in the periphery support novel diagno | exclude | duplicate: BioProject of GSE161731 |
| PRJNA679331 | BioProject | RNA | Longitudinal multi-omics identifies responses of megakaryocytes, erythroid cells and plasm | exclude | duplicate: BioProject of GSE161777 |
| PRJNA682211 | BioProject | RNA | Immune transcriptomes of highly exposed SARS-CoV-2 asymptomatic seropositive versus serone | exclude | duplicate: BioProject of GSE162562 |
| PRJNA684327 | BioProject | RNA | Single-cell RNA sequencing of Tocilizumab-treated peripheral blood mononuclear cells as an | exclude | duplicate: BioProject of GSE163014 |
| PRJNA685117 | BioProject | RNA | A diagnostic host response biosignature for COVID-19 from RNA profiling of nasal swabs and | exclude | duplicate: BioProject of GSE163151 |
| PRJNA685749 | BioProject | RNA | Spontaneous NLRP3 inflammasome-driven IL-1-β secretion is induced in severe COVID-19 patie | exclude | duplicate: BioProject of GSE163317 |
| PRJNA687221 | BioProject | RNA | Self-organized stem cell-derived human lung buds with proximo-distal patterning and novel  | exclude | duplicate: BioProject of GSE163698 |
| PRJNA690251 | BioProject | RNA | Integrated analysis of multimodal single-cell data | exclude | duplicate: BioProject of GSE164378 |
| PRJNA690255 | BioProject | RNA | Distinct B cell repertoires characterize patients with mild and severe COVID-19 [GEX] | exclude | duplicate: BioProject of GSE164379 |
| PRJNA690256 | BioProject | RNA | Distinct B cell repertoires characterize patients with mild and severe COVID-19 [VDJ] | exclude | single cell data |
| PRJNA692106 | BioProject | RNA | The transcriptional profiles of severe COVID-19 | exclude | duplicate: BioProject of GSE164805 |
| PRJNA693167 | BioProject | RNA | A Single Cell Atlas of Alveolar Development | exclude | not human (Mus musculus) |
| PRJNA693201 | BioProject | RNA | Identification of distinct immune cell subsets associated with various clinical presentati | exclude | duplicate: BioProject of GSE165080 |
| PRJNA693202 | BioProject | RNA | DNA methylation and expression profiles of whole blood in Parkinson’s disease (RNA-Seq) | exclude | duplicate: BioProject of GSE165082 |
| PRJNA693504 | BioProject | RNA | Single-cell RNA sequencing of PBMCs sampled from COVID-19 patients | exclude | duplicate: BioProject of GSE165182 |
| PRJNA697969 | BioProject | RNA | Single-cell immunophenotyping of the fetal immune response to maternal SARS-CoV-2 infectio | exclude | duplicate: BioProject of GSE165193 |
| PRJNA699562 | BioProject | RNA | Understanding COVID-19 in children and adults | exclude | duplicate: BioProject of GSE166190 |
| PRJNA699856 | BioProject | RNA | Transcriptomic characteristics and impaired immune function of patients who retest positiv | exclude | duplicate: BioProject of GSE166253 |
| PRJNA700863 | BioProject | RNA | Asymptomatic COVID-19: disease tolerance with efficient anti-viral immunity against SARS-C | exclude | duplicate: BioProject of GSE166424 |
| PRJNA701260 | BioProject | RNA | A study of differential circRNA and lncRNA expressions in COVID-19-infected peripheral blo | exclude | other assay or analyte |
| PRJNA702537 | BioProject | RNA | Single Cell RNA-Seq of PBMC from COVID-19 patients and Control subjects | exclude | duplicate: BioProject of GSE166992 |
| PRJNA702558 | BioProject | RNA | PolyA RNA-seq from whole blood of Sars-COV2-negative and -positive subjects | exclude | duplicate: BioProject of GSE167000 |
| PRJNA702763 | BioProject | RNA | A monocyte/dendritic cell molecular signature of SARS-CoV2-related multisystem  inflammato | exclude | duplicate: BioProject of GSE167028 |
| PRJNA702765 | BioProject | RNA | A monocyte/dendritic cell molecular signature of SARS-CoV2-related multisystem  inflammato | exclude | duplicate: BioProject of GSE167029 |
| PRJNA703473 | BioProject | RNA | Transcriptomic profiling of PBMCs from healthy persons and persons recovered from COVID-19 | exclude | duplicate: BioProject of GSE167158 |
| PRJNA703474 | BioProject | RNA | Single-cell RNA sequencing of PBMCs from convalescent COVID-19 patients compared to health | exclude | duplicate: BioProject of GSE167159 |
| PRJNA704611 | BioProject | RNA | Histone H4 is increased in SARS-CoV-2-infected lung epithelial cells and linked to COVID-1 | exclude | duplicate: BioProject of GSE167410 |
| PRJNA705602 | BioProject | RNA | Genome-wide PBMCs Profiling Reveals Early Detection of Immunothrombotic Dysregulation in C | exclude | duplicate: BioProject of GSE167930 |
| PRJNA707223 | BioProject | RNA | Peripheral blood gene expression according to IFIH1 rs1990760 variants in critically-ill C | exclude | duplicate: BioProject of GSE168400 |
| PRJNA708208 | BioProject | RNA | Transcriptional profiling of human KIR+ and KIR- CD8 T cells by bulk RNA-seq | exclude | duplicate: BioProject of GSE168527 |
| PRJNA713777 | BioProject | RNA | IFN-γ and TNF-α drive a CXCL10+ CCL2+ macrophage phenotype expanded in severe COVID-19 lun | exclude | duplicate: BioProject of GSE168710 |
| PRJNA716306 | BioProject | RNA | scRNAseq of Antigen Presenting Cells in severe COVID-19 patients | exclude | duplicate: BioProject of GSE169346 |
| PRJNA716611 | BioProject | RNA | Immune responses to Coronavirus disease-19 with Thymosin α 1 intervention: an observationa | exclude | duplicate: BioProject of GSE169441 |
| PRJNA717662 | BioProject | RNA | Long-term perturbation of the blood transcriptome months after SARS-CoV-2 infection. | exclude | duplicate: BioProject of GSE169687 |
| PRJNA718199 | BioProject | RNA | Single-cell RNA sequencing reveals ex vivo signatures of SARS-CoV-2-reactive T cells throu | exclude | duplicate: BioProject of GSE171037 |
| PRJNA718249 | BioProject | RNA | Case Report: Deep Immune Profiling Links Neutrophil and Plasmablast Dysregulation to Micro | exclude | duplicate: BioProject of GSE171052 |
| PRJNA718349 | BioProject | RNA | Whole blood transciptome analysis in severe Covid-19  patients and healthy donors | exclude | duplicate: BioProject of GSE171110 |
| PRJNA720088 | BioProject | RNA | Single-cell RNA-seq of human peripheral blood mononuclear cells in COVID-19 | exclude | duplicate: BioProject of GSE171555 |
| PRJNA727024 | BioProject | RNA | Differential dynamics of peripheral immune responses to acute SARS-CoV-2 infection in the  | exclude | single cell data |
| PRJNA731508 | BioProject | RNA | Regulatory KIR+RA+ T cells accumulate with age and are highly activated during viral respi | exclude | duplicate: BioProject of GSE174779 |
| PRJNA731518 | BioProject | RNA | RNAseq analysis of PS+ vs PS- non-naive CD8 T cells sorted from peripheral blood of COVID- | exclude | duplicate: BioProject of GSE174786 |
| PRJNA731855 | BioProject | RNA | Molecular profiling reveals human biomarker signatures in patients after allogeneic stem c | exclude | duplicate: BioProject of GSE174825 |
| PRJNA732416 | BioProject | RNA | Complement activation induces excessive T cell cytotoxicity in severe COVID-19 | exclude | duplicate: BioProject of GSE175450 |
| PRJNA734178 | BioProject | RNA | Identification and characterization of human CD34+ and CD34dim/- neutrophil-committed prog | exclude | duplicate: BioProject of GSE175880 |
| PRJNA735258 | BioProject | RNA | Immune signatures underlying post-acute COVID-19 lung sequelae | exclude | duplicate: BioProject of GSE176201 |
| PRJNA735729 | BioProject | RNA | RNA sequencing of blood monocyte subsets in COVID patients | exclude | duplicate: BioProject of GSE176290 |
| PRJNA736300 | BioProject | RNA | Impaired antigen-specific memory B cell and plasma cell responses including lack of specif | exclude | duplicate: BioProject of GSE176442 |
| PRJNA736437 | BioProject | RNA | Cardiometabolic MicroRNAs are Associated with severeerity and 28-Day Intensive Care Unit M | exclude | duplicate: BioProject of GSE176498 |
| PRJNA736483 | BioProject | RNA | Genetic ancestry effects on the response to viral infection are pervasive but cell type sp | exclude | single cell data |
| PRJNA736848 | BioProject | RNA | Effects of steroids on peripheral blood gene expression in COVID-19 patients according to  | exclude | duplicate: BioProject of GSE177025 |
| PRJNA737074 | BioProject | RNA | Type I interferon pathways in SARS-CoV-2-infected individuals from Pakistan | exclude | duplicate: BioProject of GSE177477 |
| PRJNA737416 | BioProject | RNA | Single Cell RNAseq and VDJ sequence of ACE2 reactive and non-reactive B-cells isolated fro | exclude | duplicate: BioProject of GSE178123 |
| PRJNA739210 | BioProject | RNA | An AI-guided invariant signature places MIS-C with Kawasaki disease in a continuum of host | exclude | duplicate: BioProject of GSE178491 |
| PRJNA740323 | BioProject | RNA | Rapid, scalable assessment of SARS-CoV-2 cellular immunity by whole-blood PCR | exclude | duplicate: BioProject of GSE178757 |
| PRJNA741195 | BioProject | RNA | Severe COVID-19 is Characterized by an Impaired Type I Interferon Response and Elevated Le | exclude | duplicate: BioProject of GSE178824 |
| PRJNA741686 | BioProject | RNA | Baseline signatures associated with clinical, virologic, and immunologic outcomes in patie | exclude | duplicate: BioProject of GSE178967 |
| PRJNA743862 | BioProject | RNA | RNAseq analysis of human CD4+ regulatory Tregs and Tconvs in COVID-19 and healthy donors i | exclude | duplicate: BioProject of GSE179448 |
| PRJNA744408 | BioProject | RNA | Transcriptome analysis of PBMCs reveals distinct immune response in the asymptomatic and r | exclude | duplicate: BioProject of GSE179627 |
| PRJNA744851 | BioProject | RNA | SARS-CoV-2 antigen exposure history shapes phenotypes and specificity of memory CD8 T cell | exclude | single cell data |
| PRJNA745333 | BioProject | RNA | Gene expression analysis of whole blood obtained from COVID-19 in the Osaka cohort | exclude | duplicate: BioProject of GSE179850 |
| PRJNA745903 | BioProject | RNA | SARS-CoV-2 infects the brain vasculature in COVID-19 patients and a hiPSC-derived model of | exclude | duplicate: BioProject of GSE179923 |
| PRJNA748962 | BioProject | RNA | TWAS-based translational genomics approach identifies IL10RB as the top candidate gene for | exclude | duplicate: BioProject of GSE180622 |
| PRJNA750398 | BioProject | RNA | Next Generation Sequencing of sorted central memory CD4 T cells (TCM) from young and elder | exclude | duplicate: BioProject of GSE181005 |
| PRJNA750442 | BioProject | RNA | Convalescent COVID-19 patients without comorbidities display similar immunophenotypes over | exclude | duplicate: BioProject of GSE181032 |
| PRJNA752081 | BioProject | RNA | Upregulation of cytokine signalling in platelet increases risk of thrombophilia in severe  | exclude | duplicate: BioProject of GSE181439 |
| PRJNA754621 | BioProject | RNA | Clinical and immunological signatures of severe COVID-19 in previously healthy patients wi | exclude | duplicate: BioProject of GSE182123 |
| PRJNA756948 | BioProject | RNA | TLR7 sensing of SARS-CoV-2 and IFN-I response by plasmacytoid dendritic cells leads to the | exclude | duplicate: BioProject of GSE182572 |
| PRJNA761110 | BioProject | RNA | A multi-tissue study of immune gene expression profiling highlights the key role of the na | exclude | duplicate: BioProject of GSE183071 |
| PRJNA761726 | BioProject | RNA | The innate and adaptive immune landscape of SARS-CoV-2-associated multisystem inflammatory | exclude | duplicate: BioProject of GSE183716 |
| PRJNA764083 | BioProject | RNA | Untimely TGFβ responses in COVID-19 limit antiviral functions of NK cells | exclude | duplicate: BioProject of GSE184329 |
| PRJNA764515 | BioProject | RNA | Identifying Risk Factors for Secondary Infection Post-SARS-CoV-2 Infection | exclude | duplicate: BioProject of GSE184401 |
| PRJNA767499 | BioProject | RNA | Bioinformatic comparisons of SARS-CoV-2 infection in human lungs and in K18-hACE2 and mACE | exclude | not human |
| PRJNA769665 | BioProject | RNA | Maternal-Fetal Immune Responses in Pregnant Women Infected with SARS-CoV-2 | exclude | duplicate: BioProject of GSE185557 |
| PRJNA773094 | BioProject | RNA | Transcriptional reprogramming of infiltrating neutrophils drives lung disease in severe CO | exclude | duplicate: BioProject of GSE186267 |
| PRJNA773942 | BioProject | RNA | Baricitinib attenuates the proinflammatory phase1 of COVID-19 driven by lung-infiltrating  | exclude | duplicate: BioProject of GSE186460 |
| PRJNA775870 | BioProject | RNA | Induction and YF-17D stimulation of THP-1 derived dendritic cells. | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| PRJNA777877 | BioProject | RNA | PBMC single cell transcriptomics of COVID-19 patients | exclude | duplicate: BioProject of GSE188172 |
| PRJNA779439 | BioProject | RNA | Epigenetic and transcriptomic remodeling in monocytes of severe COVID-19 patients mirrors  | exclude | duplicate: BioProject of GSE188573 |
| PRJNA783299 | BioProject | RNA | Small RNA-seq analysis of serum miRNAs in patients with severe COVID-19 | exclude | duplicate: BioProject of GSE189506 |
| PRJNA785113 | BioProject | RNA | COVID-19 mRNA vaccination via high-temporal frequency blood transcriptomics | exclude | duplicate: BioProject of GSE190001 |
| PRJNA785539 | BioProject | RNA | RNA Sequencing in COVID-19 patients identifies neutrophil activation biomarkers as promisi | exclude | duplicate: BioProject of GSE189990 |
| PRJNA787951 | BioProject | RNA | Immune transcriptome analysis of COVID-19 patients infected with SARS-CoV-2 variants carry | exclude | duplicate: BioProject of GSE190680 |
| PRJNA788374 | BioProject | RNA | mRNA vaccination in octogenarian nuns 15 and 20 months after recovery from COVID-19 elicit | exclude | duplicate: BioProject of GSE190747 |
| PRJNA788889 | BioProject | RNA | Investigation of SARS-CoV-2 variants that escape from cellular immunity in the absence of  | exclude | duplicate: BioProject of GSE190895 |
| PRJNA789645 | BioProject | RNA | Comprehensive comparison of adaptive immune response by inactivated SARS-CoV-2 vaccine bet | exclude | duplicate: BioProject of GSE191088 |
| PRJNA792807 | BioProject | RNA | Resistin Forms a Network with Cytokines Associated with Endothelial Damage and Worse Outco | exclude | duplicate: BioProject of GSE192707 |
| PRJNA794143 | BioProject | RNA | Blood leukocyte transcriptional modules and differentially expressed genes associated with | exclude | duplicate: BioProject of GSE193022 |
| PRJNA795800 | BioProject | RNA | Elucidation of the liver pathophysiology of COVID-19 patients using liver-on-a-chips | exclude | duplicate: BioProject of GSE193330 |
| PRJNA797780 | BioProject | RNA | Single-cell transcriptomic analysis of blood T cells from healthy controls, MS patients an | exclude | duplicate: BioProject of GSE193770 |
| PRJNA797951 | BioProject | RNA | ScRNA-seq Expression of APOC2 and IFI27 Identifies Four Alveolar Macrophage Superclusters  | exclude | duplicate: BioProject of GSE193782 |
| PRJNA800565 | BioProject | RNA | IFN-γ and TNF-α drive a CXCL10+ CCL2+ macrophage phenotype expanded in severe COVID-19 lun | exclude | single cell data |
| PRJNA800566 | BioProject | RNA | IFN-γ and TNF-α drive a CXCL10+ CCL2+ macrophage phenotype expanded in severe COVID-19 lun | exclude | single cell data |
| PRJNA800739 | BioProject | RNA | RNA-Sequencing of whole blood from non-human primate-neonate exposed to SARS-CoV-2. | exclude | not human |
| PRJNA802194 | BioProject | RNA | Blood transcriptomes of SARS-CoV-2 infected kidney transplant recipients demonstrate immun | exclude | duplicate: BioProject of GSE195796 |
| PRJNA802669 | BioProject | RNA | Covid-19 Infection During Pregnancy Induces Differential Gene Expression in Human Cord Blo | exclude | duplicate: BioProject of GSE195938 |
| PRJNA807340 | BioProject | RNA | Transcriptome analysis of umbilical cord blood CD34+ cells of neonates delivered from SARS | exclude | duplicate: BioProject of GSE196802 |
| PRJNA807370 | BioProject | RNA | Whole blood transcriptome from Covid-19 patients in India | exclude | duplicate: BioProject of GSE196822 |
| PRJNA808248 | BioProject | RNA | Durable Alterations of Hematopoiesis and Chromatin post-COVID-19 [MultiomeRNA] | exclude | duplicate: BioProject of GSE196988 |
| PRJNA809356 | BioProject | RNA | Transcriptomic clustering of critically-ill COVID-19 patients [Day 1-RNA-Seq] | exclude | duplicate: BioProject of GSE197204 |
| PRJNA809566 | BioProject | RNA | Transcriptomic clustering of critically-ill COVID-19 patients [Day 1-miRNA-Seq] | exclude | other assay or analyte |
| PRJNA811519 | BioProject | RNA | Exposure of PBMCs to SARS-CoV-2, as opposed to SARS-CoV, induces cell-intrinsic type I inn | exclude | duplicate: BioProject of GSE197665 |
| PRJNA813159 | BioProject | RNA | Impact of an Immune Modulator Mycobacterium-w on Adaptive Natural Killer Cells and Protect | exclude | duplicate: BioProject of GSE197976 |
| PRJNA814290 | BioProject | RNA | Functional reprogramming of monocytes in acute and convalescent severe COVID-19 patients [ | exclude | duplicate: BioProject of GSE198256 |
| PRJNA815324 | BioProject | RNA | RNA-seq data to study immune response to SARS-CoV-2 virus infection in COVID-19 Health Act | exclude | duplicate: BioProject of GSE198449 |
| PRJNA815866 | BioProject | RNA | Immune transcriptome in adult-onset Still’s disease with mild flare following administrati | exclude | duplicate: BioProject of GSE198549 |
| PRJNA816679 | BioProject | RNA | IL-8 dynamics and clinical outcomes in COVID-19 | exclude | single cell data |
| PRJNA817356 | BioProject | RNA | Homo sapiens Raw sequence reads | exclude | other assay or analyte |
| PRJNA817521 | BioProject | RNA | Single cell RNA sequencing reveals immunological rewiring at the maternal fetal interface  | exclude | single cell data |
| PRJNA819082 | BioProject | RNA | Homo sapiens Targeted Locus (Loci) | exclude | single cell data |
| PRJNA819422 | BioProject | RNA | Single-cell transcriptomics reveal a unique memory-like NK cell subset that accumulates wi | exclude | duplicate: BioProject of GSE199337 |
| PRJNA821445 | BioProject | RNA | COVID-19 Vaccine Immune Responses Study (COVIRS) | exclude | duplicate: BioProject of GSE199750 |
| PRJNA821533 | BioProject | RNA | Transcriptional Differences Between COVID-19 and Sepsis | exclude | duplicate: BioProject of GSE199816 |
| PRJNA827153 | BioProject | RNA | Biologic correlates of COVID-19 resolution using single cell RNA sequecing | exclude | duplicate: BioProject of GSE200896 |
| PRJNA830341 | BioProject | RNA | Single-cell multi-omics analyses in the DILfrequency clinical trail samples | exclude | duplicate: BioProject of GSE201197 |
| PRJNA832132 | BioProject | RNA | Prior vaccination exceeds prior infection in eliciting innate and humoral immune responses | exclude | duplicate: BioProject of GSE201530 |
| PRJNA832138 | BioProject | RNA | Heterologous ChAdOx1-BNT162b2 vaccination in a Korean cohort induces a robust genetic immu | exclude | duplicate: BioProject of GSE201534 |
| PRJNA832139 | BioProject | RNA | Heterologous ChAdOx1-BNT162b2 vaccination in a Korean cohort induces a robust genetic immu | exclude | duplicate: BioProject of GSE201533 |
| PRJNA832479 | BioProject | RNA | Impact of two and three dose COVID-19 vaccination schedule on transcriptional, humoral and | exclude | duplicate: BioProject of GSE201642 |
| PRJNA833037 | BioProject | RNA | Circulating microRNA signatures associated with disease severity and outcome in COVID-19 p | exclude | in vitro or cell culture |
| PRJNA835413 | BioProject | RNA | Epigenetic activation of antiviral sensors and effectors of interferon response pathways d | exclude | other assay or analyte |
| PRJNA837285 | BioProject | RNA | Global Transcriptomic analysis of PBMC derived from COVID-19 acute-severe patients and con | exclude | duplicate: BioProject of GSE202805 |
| PRJNA844257 | BioProject | RNA | Limited cross-variant immune response from SARS-CoV-2 Omicron BA.2 in naïve but not previo | exclude | duplicate: BioProject of GSE205244 |
| PRJNA844974 | BioProject | RNA | SARS-CoV-2 mRNA Dual Immunization Induces Innate Transcriptional Signatures, Establishes T | exclude | duplicate: BioProject of GSE205402 |
| PRJNA847067 | BioProject | RNA | Impact of pregravid obesity on anti-microbial fetal monocyte response | exclude | not human |
| PRJNA849923 | BioProject | RNA | Transcriptomic clustering of critically-ill COVID-19 patients [Day 4 RNAseq] | exclude | duplicate: BioProject of GSE206264 |
| PRJNA849968 | BioProject | RNA | SARS-CoV-2 induced immune perturbations in infants vary with disease severity and differ f | exclude | duplicate: BioProject of GSE206289 |
| PRJNA856339 | BioProject | RNA | RNA Sequencing Analysis of Peripheral Blood Mononuclear Cells Associated with COVID-19 Det | exclude | duplicate: BioProject of GSE207625 |
| PRJNA857312 | BioProject | RNA | Transcriptional effects of SARS-CoV-2 on human  M-CSF-derived macrophages | exclude | duplicate: BioProject of GSE207840 |
| PRJNA858262 | BioProject | RNA | Transcriptome analysis to explore the shared and distinct immune signature between severe  | exclude | other assay or analyte |
| PRJNA859437 | BioProject | RNA | Elevated Myl9 reflects the Myl9-containing microthrombi in SARS-CoV-2-induced exudative va | exclude | duplicate: BioProject of GSE208337 |
| PRJNA860181 | BioProject | RNA | Risk assessment with gene expression markers in sepsis development | exclude | duplicate: BioProject of GSE208587 |
| PRJNA863278 | BioProject | RNA | Response to mRNA vaccine BNT162b2 in hemodialysis patients | exclude | duplicate: BioProject of GSE209985 |
| PRJNA865889 | BioProject | RNA | scRNA-sequencing on whole lung single-cell suspensions isolated from both NDUFS2 cKO and c | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| PRJNA867004 | BioProject | RNA | TFEB and TFE3 are required for efficient cellular response against beta-coronavirus [RNA-s | exclude | duplicate: BioProject of GSE210676 |
| PRJNA873078 | BioProject | RNA | Redox Imbalance in COVID-19 Pathophysiology | exclude | duplicate: BioProject of GSE211979 |
| PRJNA875643 | BioProject | RNA | Monocyte dysfunction in COVID-19-associated mucormycosis | exclude | RNA from plasma or serum (cell-free RNA), not whole blood or PBMC |
| PRJNA876121 | BioProject | RNA | T cell cellular stress and reticulocyte signatures, but not loss of naïve T lymphocytes, c | exclude | duplicate: BioProject of GSE212584 |
| PRJNA880519 | BioProject | RNA | Serial whole blood transcriptomic analysis demonstrates activation of neutrophils, platele | exclude | duplicate: BioProject of GSE213313 |
| PRJNA884500 | BioProject | RNA | A longitudinal molecular and cellular lung atlas of lethal SARS-CoV-2 infection in K18-hAC | exclude | single cell data |
| PRJNA885120 | BioProject | RNA | A rapid method for detection of antigen-specific B cells | exclude | vaccination study |
| PRJNA889444 | BioProject | RNA | Transcriptional patterns of peripheral blood mononuclear cells (PBMCs) in SARS-CoV-2 infec | exclude | duplicate: BioProject of GSE215262 |
| PRJNA890265 | BioProject | RNA | ChIP-seq and RNA-seq on COVID-19 PBMC | exclude | other assay or analyte |
| PRJNA890374 | BioProject | RNA | Single cell transcriptomic data from CD4+ T cells from children with MIS-C | exclude | single cell data |
| PRJNA890375 | BioProject | RNA | Single cell transcriptomic data from CD4+ T cells from children with MIS-C | exclude | single cell data |
| PRJNA891054 | BioProject | RNA | Molecular states during acute COVID-19 reveal distinct etiologies of long-term sequelae | exclude | duplicate: BioProject of GSE215865 |
| PRJNA891278 | BioProject | RNA | Multiomic chromatin accessibility and gene expression in peripheral blood predicts future  | exclude | duplicate: BioProject of GSE215915 |
| PRJNA891717 | BioProject | RNA | Heterogeneity of neutrophils and inflammatory responses in patients with COVID-19 and heal | exclude | duplicate: BioProject of GSE216020 |
| PRJNA893694 | BioProject | RNA | Human IRF1 governs phagocytic IFN-γ immunity to mycobacteria but not cell-intrinsic IFN-a/ | exclude | duplicate: BioProject of GSE216488 |
| PRJNA893695 | BioProject | RNA | Human IRF1 governs phagocytic IFN-γ immunity to mycobacteria but not cell-intrinsic IFN-a/ | exclude | duplicate: BioProject of GSE216486 |
| PRJNA894159 | BioProject | RNA | Transcriptomic profiling of peripheral immune cells reveals distinct innate immune signatu | exclude | duplicate: BioProject of GSE216529 |
| PRJNA895695 | BioProject | RNA | SARS-CoV-2-induced disruption of a vascular bed in a microphysiological system caused by t | exclude | duplicate: BioProject of GSE216851 |
| PRJNA898633 | BioProject | RNA | Neutrophil profiles of pediatric COVID-19 and multisystem inflammatory syndrome in childre | exclude | duplicate: BioProject of GSE217370 |
| PRJNA900574 | BioProject | RNA | Cardiac-related long non-coding RNA profiling in the peripheral blood mononuclear cells of | exclude | duplicate: BioProject of GSE217799 |
| PRJNA901461 | BioProject | RNA | COVID-19 patient samples collected by the PREDICT-19 consortium | exclude | duplicate: BioProject of GSE217948;GSE282464 |
| PRJNA910224 | BioProject | RNA | Transcriptomic analysis of double-negative B cell subsets in IgG4-related disease | exclude | duplicate: BioProject of GSE220582 |
| PRJNA910999 | BioProject | RNA | Time-series transcriptome analysis of peripheral blood mononuclear cells obtained from ind | exclude | duplicate: BioProject of GSE220682 |
| PRJNA912697 | BioProject | RNA | Dynamic Activity in cis-Regulatory Elements of Leukocytes Identifies Transcription Factor  | exclude | duplicate: BioProject of GSE221066 |
| PRJNA913576 | BioProject | RNA | Predicting severity in COVID-19 disease using sepsis blood gene expression signatures | exclude | duplicate: BioProject of GSE221234 |
| PRJNA918688 | BioProject | RNA | Post-COVID symptoms and endotypes are associated with the inability to modulate the trajec | exclude | duplicate: BioProject of GSE222253 |
| PRJNA921862 | BioProject | RNA | Unresolved Immune Dysfunction is Lethal in Both COVID-19 and non-COVID-19 Sepsis | exclude | duplicate: BioProject of GSE222393 |
| PRJNA922090 | BioProject | RNA | Demultiplexing Ig repertoires by parallel mRNA/DNA sequencing shows major differential alt | exclude | duplicate: BioProject of GSE222431 |
| PRJNA922506 | BioProject | RNA | Mapping disease-associated regulatory circuits by cell type from single-cell multiomics da | exclude | single cell data |
| PRJNA927243 | BioProject | RNA | Deficient radiation transcription response in COVID-19 patients | exclude | duplicate: BioProject of GSE223885 |
| PRJNA929329 | BioProject | RNA | Treatment of acute respiratory distress syndrome caused by COVID-19 with human umbilical c | exclude | other assay or analyte |
| PRJNA929502 | BioProject | RNA | COVD-19 infection during pregnancy induces a hematopoietic stem cell signature in cord blo | exclude | duplicate: BioProject of GSE224063 |
| PRJNA929954 | BioProject | RNA | Transcriptional effects of siRNA-mediated MAFB  knockdown on SARS-CoV-2-infected human M-C | exclude | duplicate: BioProject of GSE224131 |
| PRJNA930046 | BioProject | RNA | Single-cell multi-omics analysis revealing immune features of inactivated SARS-CoV-2 vacci | exclude | duplicate: BioProject of GSE224198 |
| PRJNA932783 | BioProject | RNA | Transcriptional effects of SARS-CoV-2 and their VLPs (virus-like particles) on human GM-CS | exclude | duplicate: BioProject of GSE224845 |
| PRJNA934283 | BioProject | RNA | Assessing chromatin accessibility in human classical monocytes following vaccination with  | exclude | duplicate: BioProject of GSE225165 |
| PRJNA934362 | BioProject | RNA | Effects of the circulating environment of COVID-19 on platelet and neutrophil behavior | exclude | duplicate: BioProject of GSE225217 |
| PRJNA940730 | BioProject | RNA | RNA sequencing analysis of lung epithelial cells in conditional Ndufs2 knockout mice with  | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| PRJNA940746 | BioProject | RNA | RNA sequencing analysis of lung epithelial cells in conditional SDHD knockout mice and con | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| PRJNA940973 | BioProject | RNA | RNA sequencing analysis of lung epithelial cells in conditional Ndufs2 knockout mice treat | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| PRJNA940986 | BioProject | RNA | RNA sequencing analysis of AT2 cells from conditional Ndufs2 knockout mice and control mic | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| PRJNA940992 | BioProject | RNA | RNA-seq of cultured lung epithelial cells from P6 wild-type mice with and without Piericid | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| PRJNA941219 | BioProject | RNA | Immune dysfunction in COVID-19 Associated Mucormycosis | exclude | duplicate: BioProject of GSE226677 |
| PRJNA942421 | BioProject | RNA | Whole-genome and single-cell RNA sequencing data for Korean COVID-19 patient samples | exclude | single cell data |
| PRJNA942422 | BioProject | RNA | Whole-genome and single-cell RNA sequencing data for Korean COVID-19 patient samples | exclude | single cell data |
| PRJNA943190 | BioProject | RNA | Early differentially expressed immunological genes   in mild and severe COVID-19 | exclude | duplicate: BioProject of GSE227080 |
| PRJNA943250 | BioProject | RNA | Gene expression profiles of 65 blood samples recovered from 46 COVID-19 convalescent patie | exclude | duplicate: BioProject of GSE227116 |
| PRJNA944676 | BioProject | RNA | Gene expression profile of peripheral blood mononuclear cells from mild and severe COVID-1 | exclude | duplicate: BioProject of GSE227341 |
| PRJNA945865 | BioProject | RNA | Effects of Dexamethasone in COVID-19 patients | exclude | duplicate: BioProject of GSE227585 |
| PRJNA947172 | BioProject | RNA | The protective role of the HSP90 inhibitor, STA-9090, in lungs of SARS-CoV-2-infected Syri | exclude | duplicate: BioProject of GSE227846 |
| PRJNA947544 | BioProject | RNA | A SINGLE-CELL RNA-SEQ ANALYSIS UNRAVELS THE HETEROGENEITY OF PRIMARY CULTURED HUMAN CORNEA | exclude | duplicate: BioProject of GSE227942 |
| PRJNA949252 | BioProject | RNA | Whole blood transcriptional profiling of pulmonary functional sequelae in ARDS- secondary  | exclude | duplicate: BioProject of GSE228320 |
| PRJNA951615 | BioProject | RNA | Morphologic and molecular analysis of liver injury after SARS-CoV-2 vaccination reveals di | exclude | duplicate: BioProject of GSE228779 |
| PRJNA951716 | BioProject | RNA | COVID-19 vaccination is associated with attenuation of inflammatory responses during break | exclude | duplicate: BioProject of GSE228839 |
| PRJNA951717 | BioProject | RNA | COVID-19 vaccination is associated with attenuation of inflammatory responses during break | exclude | duplicate: BioProject of GSE228840 |
| PRJNA951718 | BioProject | RNA | COVID-19 vaccination is associated with attenuation of inflammatory responses during break | exclude | vaccination study |
| PRJNA954407 | BioProject | RNA | Morphologic and molecular analysis of liver injury after SARS-CoV-2 vaccination reveals di | exclude | vaccination study |
| PRJNA954694 | BioProject | RNA | scRNA-seq of SARS-CoV-2 RBD specific B cells | exclude | single cell data |
| PRJNA958094 | BioProject | RNA | Cytokinopathy with aberrant cytotoxic lymphocytes and profibrotic myeloid response in SARS | exclude | duplicate: BioProject of GSE230227 |
| PRJNA964769 | BioProject | RNA | Upper respiratory and peripheral blood gene expression profiles of children, adolescents,  | exclude | duplicate: BioProject of GSE231409 |
| PRJNA966240 | BioProject | RNA | SARS-COV-2 putative mi-RNAs may target host immune system and cellular mechanisms to media | exclude | other assay or analyte |
| PRJNA966242 | BioProject | RNA | SARS-COV-2 putative mi-RNAs may target host immune system and cellular mechanisms to media | exclude | other assay or analyte |
| PRJNA967265 | BioProject | RNA | Repair of Acute Respiratory Distress Syndrome in COVID-19 by Stromal Cells (REALIST-COVID  | exclude | duplicate: BioProject of GSE231660 |
| PRJNA971244 | BioProject | RNA | SARS-CoV-2 viral RNA disseminates to hamster toes and associates with localized IFN-I prod | exclude | duplicate: BioProject of GSE232226 |
| PRJNA972668 | BioProject | RNA | Unique molecular signatures sustained in circulating monocytes and regulatory T cells inco | exclude | duplicate: BioProject of GSE232523 |
| PRJNA973610 | BioProject | RNA | Distinct and overlapping immunological responses to SARS-CoV-2 and Mycobacterium tuberculo | exclude | duplicate: BioProject of GSE232705 |
| PRJNA975506 | BioProject | RNA | Memory B cell development elicited by mRNA booster vaccinations in the elderly | exclude | duplicate: BioProject of GSE233230 |
| PRJNA975860 | BioProject | RNA | Immune landscape of PBMC based on single-cell rna sequancing in End Stage Renal Disease | exclude | duplicate: BioProject of GSE233315 |
| PRJNA976606 | BioProject | RNA | Clonal landscape of autoantibody-secreting plasmablasts in COVID-19 patients | exclude | duplicate: BioProject of GSE233522 |
| PRJNA983711 | BioProject | RNA | Association of Neutrophil Gene Expression with 28-Day Ventilator Weaning in COVID-19 ARDS  | exclude | duplicate: BioProject of GSE234902 |
| PRJNA983718 | BioProject | RNA | Association of Neutrophil Gene Expression with 28-Day Ventilator Weaning in COVID-19 ARDS  | exclude | duplicate: BioProject of GSE234903 |
| PRJNA983776 | BioProject | RNA | Durable Alterations of Hematopoiesis and Chromatin post-COVID-19 [bulkRNA-Seq] | exclude | duplicate: BioProject of GSE234939 |
| PRJNA985376 | BioProject | RNA | Integrated Analysis of Multimodal Single-Cell Datasets for SARS-CoV-2 Vaccination | exclude | single cell data |
| PRJNA985377 | BioProject | RNA | Integrated Analysis of Multimodal Single-Cell Datasets for SARS-CoV-2 Vaccination | exclude | single cell data |
| PRJNA988058 | BioProject | RNA | Single‐cell RNA sequencing reveals characteristics of myeloid cells in pulmonary post-acut | exclude | duplicate: BioProject of GSE235938 |
| PRJNA988459 | BioProject | RNA | Distinct pulmonary and systemic effects of dexamethasone in severe COVID-19 | exclude | duplicate: BioProject of GSE236030 |
| PRJNA989171 | BioProject | RNA | Key patient demographics shape innate immune topography in non-critical hypoxic COVID-19 p | exclude | duplicate: BioProject of GSE236177 |
| PRJNA991080 | BioProject | RNA | Serum from COVID-19 patients promotes endothelial cells activation through protease-activa | exclude | in vitro or cell culture |
| PRJNA991578 | BioProject | RNA | Recapitulation of pathophysiological features of AD in SARS-CoV-2 infected subjects | exclude | duplicate: BioProject of GSE236562 |
| PRJNA991973 | BioProject | RNA | Stimulation of PSTPIP1 to trigger proinflammatory responses in asymptomatic SARS-CoV-2 inf | exclude | duplicate: BioProject of GSE236651 |
| PRJNA993792 | BioProject | RNA | Distinct pulmonary and systemic effects of dexamethasone in severe COVID-19 | exclude | duplicate: BioProject of GSE237109 |
| PRJNA995237 | BioProject | RNA | Repertoire-based scRNA-Seq phenotyping delineates responding helper T cell subsets | exclude | single cell data |
| PRJNA997309 | BioProject | RNA | A Phase I, First-in-Human Study to Evaluate the Safety and Tolerability, Pharmacokinetics  | exclude | duplicate: BioProject of GSE237965 |
| GSE100150 | GEO | RNA | A Transcriptome Fingerprinting Assay for Clinical Immune Monitoring | exclude | vaccination study |
| GSE100504 | GEO | RNA | Primary human airway epithelial cell transcriptome response to a wild type Mers-CoV (icMER | exclude | single cell data |
| GSE101207 | GEO | RNA | High-throughput single cell transcriptome analysis and CRISPR screen identify key β cell-s | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE101709 | GEO | RNA | Gene expression analysis of Influenza vaccine response in Young and Old - Year 4 | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE101710 | GEO | RNA | Gene expression analysis of Influenza vaccine response in Young and Old - Year 5 | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE10211 | GEO | RNA | Airway Epithelial Cell Response to Sendai virus infection | exclude | not human (Mus musculus) |
| GSE103174 | GEO | RNA | Expression data from lung tissue in mild-moderate COPD | exclude | single cell data |
| GSE103232 | GEO | RNA | Deep sequencing of endogenous mRNA from Homo sapiens from four age categories in blood | exclude | single cell data |
| GSE10325 | GEO | RNA | Expression data from human peripheral blood subsets | exclude | no characteristic splits 10 COVID samples from 10 others and no severity field (automatic; best : 0 against 0) |
| GSE103354 | GEO | RNA | A revised airway epithelial hierarchy includes CFTR-expressing ionocytes | exclude | not human (Mus musculus) |
| GSE103544 | GEO | RNA | Single-cell RNA-seq analysis of human CD14+ monocytes | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE103825 | GEO | RNA | A conserved transcriptional response to intranasal Ebola virus exposure in nonhuman primat | exclude | not human (Macaca fascicularis) |
| GSE103918 | GEO | RNA | Single cell transcriptomic profiling of pluripotent stem cell-derived SCGB3A2+ airway epit | exclude | single cell data |
| GSE104714 | GEO | RNA | Clustering gene expression time series data using an infinite Gaussian process mixture mod | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE106757 | GEO | RNA | Analysis of monocyte subsets from pretreatment RCC patients compared to healthy donors. | exclude | vaccination study |
| GSE106840 | GEO | RNA | Analysis of human classical and non-classical monocytes compared to their dendritic cell c | exclude | vaccination study |
| GSE107011 | GEO | RNA | RNA-Seq profiling of 29 immune cell types and peripheral blood mononuclear cells | exclude | single cell data |
| GSE107105 | GEO | RNA | Microarray analysis of freshly isolated synovial fibroblast subsets in patients with rheum | exclude | no characteristic splits 10 COVID samples from 10 others and no severity field (automatic; best : 0 against 0) |
| GSE107746 | GEO | RNA | Transcriptome Landscape of Human Folliculogenesis Reveals Oocytes and Granulosa Cells Inte | exclude | single cell data |
| GSE107747 | GEO | RNA | Single-cell RNA sequencing of circulating tumor cells in hepatocellular carcinoma | exclude | single cell data |
| GSE108134 | GEO | RNA | Ambient Pollution Related Reprogramming of the Human Small Airway Epithelial Transcriptome | exclude | single cell data |
| GSE108363 | GEO | RNA | Mycobacterium tuberculosis Exploits a Molecular Off Switch of the Immune System for Intrac | exclude | in vitro or cell culture |
| GSE108571 | GEO | RNA | Single-cell RNA sequencing of vascularized brain organoids | exclude | single cell data |
| GSE10896 | GEO | RNA | Impact of curcumin on human monocytes (U937 cells) exposed to oxidative stress | exclude | no characteristic splits 10 COVID samples from 10 others and no severity field (automatic; best : 0 against 0) |
| GSE109037 | GEO | RNA | 10x Genomics Drop-seq single-cell RNA-seq of isolated Adult human spermatogonia, spermatoc | exclude | single cell data |
| GSE110236 | GEO | RNA | IFNb/LPS activates convergence of inflammatory gene expression patterns between adipocytes | exclude | not human (Mus musculus) |
| GSE110907 | GEO | RNA | Identification of tumor suppressor microRNAs by integrative miRNA and mRNA sequencing of m | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE110973 | GEO | RNA | Ex-vivo Human Hematopoietic Stem Cell Expansion Requires Coordination of Cellular Reprogra | exclude | single cell data |
| GSE111014 | GEO | RNA | Chromatin mapping and single-cell immune profiling defines the temporal dynamics of ibruti | exclude | single cell data |
| GSE111360 | GEO | RNA | Organoid modeling of the tumor immune microenvironment | exclude | single cell data |
| GSE112101 | GEO | RNA | RNA-seq of nine primary human cell types exposed in vitro to methylprednisolone | exclude | in vitro or cell culture |
| GSE112104 | GEO | RNA | Gene expression profiling identifies a gene signature that predicts progression to active  | exclude | no characteristic splits 10 COVID samples from 10 others and no severity field (automatic; best : 0 against 0) |
| GSE112356 | GEO | RNA | A genomic study of the contribution of DNA methylation to regulatory evolution in primates | exclude | not human |
| GSE112570 | GEO | RNA | Single Cell RNA-Seq profiling human embryonic kidney cortex cells | exclude | single cell data |
| GSE112845 | GEO | RNA | PBMC Fixation and Processing for Chromium Single-Cell RNA Sequencing | exclude | single cell data |
| GSE113036 | GEO | RNA | Pericyte-to-neuron reprogramming via unfolding of a neural stem cell-like program | exclude | single cell data |
| GSE113196 | GEO | RNA | Single Cell RNA sequencing of Adult Human Breast Epithelial Cells [Individual 4..7] | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE11343 | GEO | RNA | Rosiglitazone Treatment Reduces Diabetic Neuropathy in STZ treated DBA/2J mice | exclude | not human (Mus musculus) |
| GSE11352 | GEO | RNA | Timecourse of estradiol (10nM) exposure in MCF7 breast cancer cells. | exclude | 18 samples: cannot reach 10 participants per class |
| GSE11393 | GEO | RNA | Monocyte gene expression profiling in familial combined hyperlipidemia and its modificatio | exclude | 9 samples: cannot reach 10 participants per class |
| GSE114052 | GEO | RNA | AR regulated transcriptome in prostate cancer cells | exclude | 12 samples: cannot reach 10 participants per class |
| GSE114284 | GEO | RNA | IRF1 regulates IFN dependent and independent gene expression | exclude | in vitro or cell culture |
| GSE114466 | GEO | RNA | Transcriptional profile from peripheral blood in patients infected with avian influenza H7 | exclude | no characteristic splits 10 COVID samples from 10 others and no severity field (automatic; best : 0 against 0) |
| GSE114530 | GEO | RNA | Single cell RNA-sequencing of human fetal kidneys | exclude | single cell data |
| GSE114588 | GEO | RNA | Characterization of antibiotic resistance and host-microbiome interactions in the human up | exclude | no characteristic splits 10 COVID samples from 10 others and no severity field (automatic; best : 0 against 0) |
| GSE114724 | GEO | RNA | Single-cell Map of Diverse Immune Phenotypes in the Breast Tumor Microenvironment - 5' RNA | exclude | single cell data |
| GSE114802 | GEO | RNA | Single cell RNA-Seq of four human kidney organoids | exclude | single cell data |
| GSE115007 | GEO | RNA | Single-cell RNA-seq analysis of human tumor ascites dendritic cells and macrophages | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE115149 | GEO | RNA | Single cell sequencing of human mesenchymal stem/stromal cells (MSCs) upon inflammatory li | exclude | single cell data |
| GSE115189 | GEO | RNA | Single cell profiling of peripheral blood mononuclear cells from healthy human donor | exclude | single cell data |
| GSE115469 | GEO | RNA | Dissecting the human liver cellular landscape by single cell RNA-seq reveals novel intrahe | exclude | single cell data |
| GSE115571 | GEO | RNA | Single-cell transcriptomics reveals distinct microglia signatures under inflammation | exclude | not human (Mus musculus) |
| GSE115770 | GEO | RNA | A network of transcriptome modules demonstrates mechanistic pathways of both virus induced | exclude | single cell data |
| GSE115785 | GEO | RNA | Whole blood RNA-Sequencing from 12 Rhesus Macaques following  Ebola Makona Infection | exclude | not human (Macaca mulatta) |
| GSE115830 | GEO | RNA | Gene expression analysis of human bronchial airway epithelial BEAS-2B cells stimulated by  | exclude | not blood |
| GSE116113 | GEO | RNA | A single-cell RNA-Seq of human hepatocyte-derived liver progenitor-like cells | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE116250 | GEO | RNA | RNA-seq of heart failure in human left ventricles | exclude | not blood |
| GSE117156 | GEO | RNA | Single cell dissection of plasma cell heterogeneity in symptomatic and asymptomatic myelom | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE117211 | GEO | RNA | Kidney micro-organoids generated in suspension culture | exclude | single cell data |
| GSE117261 | GEO | RNA | Systems Analysis of the Human Pulmonary Arterial Hypertension Lung Transcriptome | exclude | not blood |
| GSE117403 | GEO | RNA | Single Cell RNA-sequencing of cell types isolated by FACS from normal human prostates | exclude | single cell data |
| GSE117415 | GEO | RNA | RNA-seq in SUNE-1 cells after downregulation of DANCR expression | exclude | not blood |
| GSE117824 | GEO | RNA | RNA-seq of Single-Cell Genotyping of Transcriptomes | exclude | single cell data |
| GSE117827 | GEO | RNA | Host Gene Expression in Nasal and Blood Samples for the Diagnosis of Viral Respiratory Inf | exclude | not adults |
| GSE11809 | GEO | RNA | Lung transcriptome from M. avium infected WT, IFNg-/- and IRF1-/- mice | exclude | not human (Mus musculus) |
| GSE118127 | GEO | RNA | Single-cell analysis of adult human ovary using 10X genomics | exclude | single cell data |
| GSE118184 | GEO | RNA | Comparative analysis of kidney organoid and adult human kidney single cell and single nucl | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE118189 | GEO | RNA | ATAC-seq data | exclude | other assay or analyte |
| GSE118257 | GEO | RNA | Altered oligodendrocyte heterogeneity in Multiple sclerosis. | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE118259 | GEO | RNA | Gene expression profile of Endothelial colony-forming cells (ECFCs) isolated from patients | exclude | 13 samples: cannot reach 10 participants per class |
| GSE118773 | GEO | RNA | Cell-to-cell variation in defective virus expression and effect on host response during in | exclude | single cell data |
| GSE11886 | GEO | RNA | Gene expression analysis of macrophages from ankylosing spondylitis patients reveals inter | exclude | no characteristic splits 10 COVID samples from 10 others and no severity field (automatic; best : 0 against 0) |
| GSE11906 | GEO | RNA | Quality Control in Microarray Assessment of Gene Expression in Human Airway Epithelium | exclude | single cell data |
| GSE119117 | GEO | RNA | Longitudinal transcriptomic characterization of the immune response to acute hepatitis C v | exclude | no characteristic splits 10 COVID samples from 10 others and no severity field (automatic; best : 0 against 0) |
| GSE119212 | GEO | RNA | Characterization of the transcriptional landscape of human developing hippocampus by scRNA | exclude | single cell data |
| GSE119373 | GEO | RNA | single cell RNA-sequencing data from CD4+CD25+ human cells | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE119506 | GEO | RNA | Single-cell RNA-seq analysis of human tonsil CD14+ cells | exclude | single cell data |
| GSE119507 | GEO | RNA | Single-cell RNA-seq analysis of human tonsil CD4 T cells | exclude | single cell data |
| GSE119561 | GEO | RNA | Fate-mapping within human iPSC-derived kidney organoids reveals conserved mammalian nephro | exclude | single cell data |
| GSE119594 | GEO | RNA | Precise Gene Editing Preserves Hematopoietic Stem Cell Function Following Transient p53-Me | exclude | single cell data |
| GSE119856 | GEO | RNA | Global transcriptional profiling unveils the interferon network in blood and tissues acros | exclude | not human (Mus musculus) |
| GSE119926 | GEO | RNA | Single cell RNA-seq analysis of medulloblastoma | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE119939 | GEO | RNA | IL-1 inhibition in systemic sclerosis skin. | exclude | no characteristic splits 10 COVID samples from 10 others and no severity field (automatic; best : 0 against 0) |
| GSE120221 | GEO | RNA | Human Bone Marrow Assessment by Single Cell RNA Sequencing, Mass Cytometry and Flow Cytome | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE120446 | GEO | RNA | Human Bone Marrow Assessment by Single Cell RNA Sequencing, Mass Cytometry and Flow Cytome | exclude | single cell data |
| GSE121267 | GEO | RNA | Molecular and functional heterogeneity of IL-10-producing CD4+ T cells [Human 10x] | exclude | single cell data |
| GSE121567 | GEO | RNA | Next Generation Sequencing Facilitates Quantitative Analysis of PSKγ1ox and control Transc | exclude | not human (Arabidopsis thaliana) |
| GSE121600 | GEO | RNA | Single-cell RNA-seq reveals novel cell differentiation dynamics during human airway epithe | exclude | single cell data |
| GSE121611 | GEO | RNA | Single-Cell Transcriptomic Analysis of Human Lung Reveals Complex Multicellular Changes Du | exclude | not human (Mus musculus) |
| GSE122342 | GEO | RNA | Fusion of Human Brain Organoids Assembles Bidirectional Connectivity Between Thalamus and  | exclude | single cell data |
| GSE122703 | GEO | RNA | Single-cell profiling of cutaneous T-cell lymphoma reveals underlying heterogeneity predic | exclude | single cell data |
| GSE122846 | GEO | RNA | Single cell profiling of CD4+ T cells from human fetal intestine | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE122960 | GEO | RNA | Single-Cell Transcriptomic Analysis of Human Lung Reveals Complex Multicellular Changes Du | exclude | single cell data |
| GSE123926 | GEO | RNA | Single cell RNA-seq of 25 thousand epithelial cells from xenografts of triple-negative bre | exclude | single cell data |
| GSE12446 | GEO | RNA | Endometrium of hormone-treated postmenopausal women | exclude | no characteristic splits 10 COVID samples from 10 others and no severity field (automatic; best : 0 against 0) |
| GSE124472 | GEO | RNA | Single Cell RNA-Seq profiling of human embryonic kidney outer and inner cortical cells and | exclude | single cell data |
| GSE124494 | GEO | RNA | Single-cell survey of human lymphatics unveils marked endothelial cell heterogeneity and m | exclude | single cell data |
| GSE124742 | GEO | RNA | Pancreas single cell patch-seq links physiologic dysfunction in diabetes to transcriptomic | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE124810 | GEO | RNA | Epigenetic modulation of β-cells by interferon-α via PNPT11-miR-26a-TET2 triggers autoimmu | exclude | not human |
| GSE125680 | GEO | RNA | Unleashing Type-2 Dendritic Cells | exclude | single cell data |
| GSE125970 | GEO | RNA | Single-cell transcriptome analysis reveals differential nutrient absorption functions in h | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE126030 | GEO | RNA | A single cell reference map for human blood and tissue T cell activation | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE126085 | GEO | RNA | Single-cell RNA-sequencing of CD14+ monocyte differentiation with M-CSF stimulus | exclude | single cell data |
| GSE126124 | GEO | RNA | Concordance between gene expression in peripheral whole blood and colonic tissue in childr | exclude | not adults |
| GSE126427 | GEO | RNA | Double-stranded Alu RNAs, type 1 IFN responses, and relapsing remitting multiple sclerosis | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE127465 | GEO | RNA | Single cell transcriptomics of human and mouse lung cancers reveals conserved myeloid popu | exclude | single cell data |
| GSE127792 | GEO | RNA | Inflammatory Mediator Expression Associated With Antibody Response Induced by Live Attenua | exclude | vaccination study |
| GSE128033 | GEO | RNA | Proliferating SPP1/MERTK-expressing macrophages in idiopathic pulmonary fibrosis | exclude | single cell data |
| GSE128066 | GEO | RNA | A Bayesian mixture model for clustering droplet-based single cell transcriptomic data from | exclude | single cell data |
| GSE128101 | GEO | RNA | Subcellular pathways shared by afflicted patients and mutant mice identify a new drug trea | exclude | not human |
| GSE128169 | GEO | RNA | Single-cell analysis reveals fibroblast heterogeneity and myofibroblasts in systemic scler | exclude | single cell data |
| GSE128266 | GEO | RNA | A Unique Epigenomic Landscape Defines Human Erythropoiesis (ATAC-seq) | exclude | single cell data |
| GSE128518 | GEO | RNA | Lipid-associated macrophages control metabolic homeostasis in a Trem2-dependent manner | exclude | single cell data |
| GSE129845 | GEO | RNA | Single-cell transcriptomic map of the human and mouse bladders | exclude | single cell data |
| GSE129882 | GEO | RNA | RNA-seq transcript and gene data on zika exposed and zika naïve samples | exclude | no characteristic splits 10 COVID samples from 10 others and no severity field (automatic; best : 0 against 0) |
| GSE130073 | GEO | RNA | Single cell RNA-sequencing of iPS cells derived multicellular human liver organoids | exclude | single cell data |
| GSE130151 | GEO | RNA | Single-cell sequencing of Klinefelter testicular biopsies | exclude | single cell data |
| GSE130238 | GEO | RNA | Nested oscillatory dynamics in cortical organoids model early human brain network developm | exclude | single cell data |
| GSE130318 | GEO | RNA | Single cell analysis of Monolayer and Organoid Primary Prostate Human Epithelial Cells Der | exclude | single cell data |
| GSE130430 | GEO | RNA | The Developmental Heterogeneity of Human Natural Killer Cells Defined by Single-cell Trans | exclude | single cell data |
| GSE130636 | GEO | RNA | Single-cell RNA-Seq Investigation of Foveal and Peripheral Expression in the Human Retina | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE130888 | GEO | RNA | Single cell transcriptome of peritoneal cells | exclude | single cell data |
| GSE131685 | GEO | RNA | Single-cell RNA sequencing of human kidney | exclude | single cell data |
| GSE131928 | GEO | RNA | single cell RNA-seq analysis of adult and paediatric IDH-wildtype Glioblastomas | exclude | single cell data |
| GSE132040 | GEO | RNA | Tabula Muris Senis: Bulk sequencing | exclude | not human (Mus musculus) |
| GSE132228 | GEO | RNA | Social network architecture of human immune cells in ZIKV-infected patients highlights a u | exclude | no characteristic splits 10 COVID samples from 10 others and no severity field (automatic; best : 0 against 0) |
| GSE132802 | GEO | RNA | Single-cell RNA-sequencing of skin, fresh blood, and cultured peripheral blood mononuclear | exclude | single cell data |
| GSE132951 | GEO | RNA | Small RNA sequencing of LA-N-2 and LAN-5 cells exposed to ciliary neurotrophic factor for  | exclude | other assay or analyte |
| GSE133028 | GEO | RNA | CSF B cells in relapsing multiple sclerosis are driven towards an antigen-experienced, inf | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE133378 | GEO | RNA | Diagnosing enterovirus meningitis via blood transcriptomics: an alternative for lumbar pun | exclude | not adults |
| GSE133704 | GEO | RNA | Epigenomics and Single-cell Sequencing Define a Developmental Hierarchy in Langerhans Cell | exclude | single cell data |
| GSE134174 | GEO | RNA | Dissecting the cellular specificity of smoking effects and reconstructing lineages in the  | exclude | single cell data |
| GSE134209 | GEO | RNA | Microarray data of human moDC treated with IFNβ and TNFα | exclude | no characteristic splits 10 COVID samples from 10 others and no severity field (automatic; best : 0 against 0) |
| GSE13428 | GEO | RNA | Gene Expression Profiling of Rat Hippocampus Following Exposure to the Acetylcholinesteras | exclude | not human (Rattus norvegicus) |
| GSE134432 | GEO | RNA | Single-cell analysis of gene expression variation and phenotype switching in melanoma | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE135069 | GEO | RNA | Antiviral potential of human IFN-α subtypes against influenza H3N2 infection in human lung | exclude | single cell data |
| GSE135188 | GEO | RNA | Intermittent exposure to whole cigarette smoke alters the differentiation of primary small | exclude | single cell data |
| GSE135618 | GEO | RNA | Transcriptional profiling of microglia; current state of the art and future perspectives | exclude | single cell data |
| GSE135779 | GEO | RNA | A single cell approach to map cellular subsets involved in Systemic Lupus Erythematosus (S | exclude | single cell data |
| GSE135893 | GEO | RNA | Single-cell RNA-sequencing reveals profibrotic roles of distinct epithelial and mesenchyma | exclude | single cell data |
| GSE136103 | GEO | RNA | Resolving the fibrotic niche of human liver cirrhosis using single-cell transcriptomics | exclude | single cell data |
| GSE136314 | GEO | RNA | Kidney organoid reproducibility across multiple human iPSC lines and diminished off target | exclude | single cell data |
| GSE136394 | GEO | RNA | Single-cell transcriptome analysis reveals gene signatures associated with T-cell persiste | exclude | single cell data |
| GSE137537 | GEO | RNA | Single-cell Transcriptomic Atlas of the Human Retina Identifies Cell Types Associated with | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE137540 | GEO | RNA | Single-cell transcriptome profiling reveals neutrophil heterogeneity in homeostasis and in | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE137864 | GEO | RNA | Single-cell transcriptome landscape of adult HSPC | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE138266 | GEO | RNA | Integrated single cell analysis of blood and cerebrospinal fluid leukocytes in multiple sc | exclude | single cell data |
| GSE138516 | GEO | RNA | Ets1 confers context-dependent activation of Notch signaling in T-cell leukemia [ChIP-Seq] | exclude | single cell data |
| GSE138734 | GEO | RNA | The RNA Atlas expands the catalog of human non-coding RNAs | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE138852 | GEO | RNA | A single-cell atlas of the human cortex reveals drivers of transcriptional changes in  Alz | exclude | single cell data |
| GSE13899 | GEO | RNA | Quercetin supplementation and CD14+ monocyte gene expression | exclude | 6 samples: cannot reach 10 participants per class |
| GSE139249 | GEO | RNA | Discovery of Specialized NK Cell Populations Infiltrating Human Melanoma Metastases | exclude | single cell data |
| GSE139324 | GEO | RNA | Immune landscape of viral- and carcinogen-drived head and neck cancer | exclude | single cell data |
| GSE13933 | GEO | RNA | Trachea Epithelium as a “ Canary” for Cigarette Smoking-induced Biologic Phenotype of Smal | exclude | single cell data |
| GSE139516 | GEO | RNA | Transcriptomic Analysis Of circRNAs/miRNAs/mRNAs upon Middle East Respiratory Syndrome Cor | exclude | in vitro or cell culture |
| GSE140203 | GEO | RNA | Integrative single-cell chromatin and transcriptome profiling uncovers cell-type specific  | exclude | single cell data |
| GSE140312 | GEO | RNA | Single-Cell RNA-Seq Analysis of Metastasis-Primary Tumor Dissimilarity in a Patient with G | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE140989 | GEO | RNA | ­Single Cell Transcriptomics  identifies Focal Segmental Glomerulosclerosis Remission  End | exclude | single cell data |
| GSE141259 | GEO | RNA | Longitudinal single-cell transcriptomics analysis of mouse lung upon bleomycin-induced inj | exclude | not human (Mus musculus) |
| GSE141634 | GEO | RNA | single-cell RNA-seq in alveolar organoids cultured in MTECplus medium | exclude | not human (Mus musculus) |
| GSE142503 | GEO | RNA | ScRNA-seq analysis and functional characterization of myeloid lincRNAs differentially expr | exclude | single cell data |
| GSE143507 | GEO | RNA | RNA sequencing data of peripheral blood isolated from Korean patients with Crohn's disease | exclude | no characteristic splits 10 COVID samples from 10 others and no severity field (automatic; best : 0 against 0) |
| GSE144106 | GEO | RNA | RNA-sequencing in MCL cell lines after silencing cyclin D1 | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE144136 | GEO | RNA | Single-nucleus RNA-seq in the post-mortem brain in major depressive disorder | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE144430 | GEO | RNA | Diversity of peripheral blood human NK cells identified by single cell RNA sequencing | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE145031 | GEO | RNA | Inflammatory Signals induce AT2 Cell-Derived Damage-Associated Transient Progenitors that  | exclude | not human (Mus musculus) |
| GSE145926 | GEO | RNA | Single-cell landscape of bronchoalveolar immune cells in COVID-19 patients | exclude | single cell data |
| GSE146507 | GEO | RNA | RNA-sequencing of L. monocytogenes under PrfA-inducing conditions refines the PrfA regulon | exclude | not human (Listeria monocytogenes) |
| GSE146943 | GEO | RNA | Single-cell sequencing of pleural effusion reveals immune response of COV2019 infection | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE146974 | GEO | RNA | Deep learning enables accurate clustering with batch effect removal in single-cell RNA-seq | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE147313 | GEO | RNA | IRF1 governs the differential Interferon-Stimulated Gene responses in human monocytes and  | exclude | in vitro or cell culture |
| GSE147505 | GEO | RNA | RNA-seq analysis in fibroblast IMR90 CTRL/MYC/MYC+SF3A3 KD | exclude | 9 samples: cannot reach 10 participants per class |
| GSE147507 | GEO | RNA | Transcriptional response to SARS-CoV-2 infection | exclude | no characteristic splits 10 COVID samples from 10 others and no severity field (automatic; best title: 35 against 2) |
| GSE147528 | GEO | RNA | Molecular characterization of selectively vulnerable neurons in Alzheimer’s Disease | exclude | single cell data |
| GSE147863 | GEO | RNA | Single-cell RNA-seq of engineered healthy human kidney organoids | exclude | single cell data |
| GSE147975 | GEO | RNA | Identification of A Drug Blocking SARS-CoV-2 Infection using Human Pluripotent Stem Cell-d | exclude | single cell data |
| GSE148071 | GEO | RNA | Comprehensive Profiling of Cancer Cells and Their Microenvironment in Advanced NSCLC | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE148077 | GEO | RNA | Cell atlas of the human fovea and peripheral retina | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE148113 | GEO | RNA | Identification of Candidate COVID-19 Therapeutics using hPSC-derived Lung Organoids | exclude | in vitro or cell culture |
| GSE1482 | GEO | RNA | Murine Neurofibromatosis | exclude | not human (Mus musculus) |
| GSE148697 | GEO | RNA | Identification of Lead Drug Candidates for COVID-19 based on Drug Screening Using Human Pl | exclude | in vitro or cell culture |
| GSE148729 | GEO | RNA | Gene expression profiling of SARS-CoV-1/2 infected human cell lines at bulk and single-cel | exclude | single cell data |
| GSE148829 | GEO | RNA | SARS-CoV-2 receptor ACE2 is an interferon-stimulated gene in human airway epithelial cells | exclude | single cell data |
| GSE148881 | GEO | RNA | Host-viral infection maps at single-cell resolution reveal signatures of severe COVID-19 p | exclude | single cell data |
| GSE149443 | GEO | RNA | Benchmarking of viruses detection from single-cell RNA sequencing data | exclude | not human (Mus musculus) |
| GSE149583 | GEO | RNA | Time-course analysis of transcriptomic changes of mesothelial cells after stress exposure | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE149687 | GEO | RNA | Attenuation of SARS-CoV-2 Infection by Losartan in Human Kidney Organoids | exclude | in vitro or cell culture |
| GSE149689 | GEO | RNA | Immunophenotyping of COVID-19 and Influenza Underscores the Association of Type I IFN Resp | exclude | single cell data |
| GSE150132 | GEO | RNA | GM-CSF-producing T helper cells: a distinct subset  of T helper cells that acquire Th1 phe | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE150311 | GEO | RNA | Single-cell RNA sequencing of intermediates during human primed and naïve reprogramming. | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE150316 | GEO | RNA | Spectrum of Viral Load and Host Response Seen in Autopsies of SARS-CoV-2 Infected Lungs | exclude | no characteristic splits 10 COVID samples from 10 others and no severity field (automatic; best : 0 against 0) |
| GSE150392 | GEO | RNA | RNA-Seq of Human iPSC-cardiomyocytes infected with SARS-CoV-2 | exclude | in vitro or cell culture |
| GSE150708 | GEO | RNA | Modeling COVID-19 with Human Pluripotent Stem Cell-Derived Cells Reveals Synergistic Effec | exclude | 5 samples: cannot reach 10 participants per class |
| GSE15072 | GEO | RNA | Mitochondrial dysregulation and oxidative stress in patients with chronic kidney disease | exclude | no characteristic splits 10 COVID samples from 10 others and no severity field (automatic; best : 0 against 0) |
| GSE150728 | GEO | RNA | A single-cell atlas of the peripheral immune response to severe COVID-19 | exclude | single cell data |
| GSE150861 | GEO | RNA | Single-cell RNA-seq of PBMCs from severe COVID-19 patients | exclude | single cell data |
| GSE151161 | GEO | RNA | Blocking of the CD80/86 axis as a therapeutic approach to prevent progression to more seve | exclude | single cell data |
| GSE151263 | GEO | RNA | Single Cell RNA Sequencing Identifies an Early Monocyte Gene Signature in Acute Respirator | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE151764 | GEO | RNA | Two distinct immunopathological profiles in lungs of lethal COVID-19 | exclude | vaccination study |
| GSE151803 | GEO | RNA | Identification of Lead Drug Candidates for COVID-19 based on Drug Screening Using Human Pl | exclude | in vitro or cell culture |
| GSE151878 | GEO | RNA | Cardiomyocytes Infected by SARS-CoV-2 Recruit Cardiotoxic Macrophages [scRNA-seq dataset] | exclude | single cell data |
| GSE151879 | GEO | RNA | Cardiomyocytes Infected by SARS-CoV-2 Recruit Cardiotoxic Macrophages [bulk RNA-seq datase | exclude | not blood |
| GSE151928 | GEO | RNA | Bronchoalveolar lavage of healthy human adults | exclude | single cell data |
| GSE152075 | GEO | RNA | In vivo antiviral host response to SARS-CoV-2 by viral load, sex, and age [dataset I] | exclude | no characteristic splits 10 COVID samples from 10 others and no severity field (automatic; best : 0 against 0) |
| GSE152522 | GEO | RNA | Imbalance of regulatory and cytotoxic SARS-CoV-2-reactive CD4+ T cells in COVID-19 | exclude | single cell data |
| GSE152586 | GEO | RNA | Transcriptomic analysis of human alveolar type II cell organoids post infection of SARS-Co | exclude | in vitro or cell culture |
| GSE153218 | GEO | RNA | Bulk RNA sequencing of SARS-CoV-2 infected alveolar type I- and type II like cells | exclude | in vitro or cell culture |
| GSE153277 | GEO | RNA | Global transcriptomic analyses of pluripotent stem cell-derived type II alveolar epithelia | exclude | single cell data |
| GSE153370 | GEO | RNA | Single cell transcriptome of fresh and UM171-expanded CD34+ cord blood cells | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE153428 | GEO | RNA | Expression data from infant nasal airway epithelial cells treated with IFN-lambda, poly(I: | exclude | not adults |
| GSE153610 | GEO | RNA | Immunologic characterization of COVID-19 patients with hematologic cancer: biologic and cl | exclude | no characteristic splits 10 COVID samples from 10 others and no severity field (automatic; best : 0 against 0) |
| GSE153835 | GEO | RNA | Single-cell analysis reveals divergent responses of human dendritic cells to the MVA vacci | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE153931 | GEO | RNA | Severely ill COVID-19 patients display impaired exhaustion features in SARS-CoV-2-reactive | exclude | single cell data |
| GSE153970 | GEO | RNA | Primary Human Airway Epithelial Cultures infected with SARS-CoV-2 | exclude | not blood |
| GSE154244 | GEO | RNA | Expression profiling of SARS-CoV-2 infected patients by RNA-seq | exclude | not blood |
| GSE154311 | GEO | RNA | RNA-seq on CD16Int versus CD16High in severe COVID-19 patients | exclude | not blood |
| GSE154567 | GEO | RNA | Single cell RNAseq reveal covid19 patient immune response | exclude | single cell data |
| GSE15457 | GEO | RNA | Gene Expression Profiling of Cellular Dose-Response to Hypochlorous Acid (HOCl) | exclude | not human (Mus musculus) |
| GSE154613 | GEO | RNA | Modulating the transcriptional landscape of SARS-CoV-2 as an effective method for developi | exclude | in vitro or cell culture |
| GSE154741 | GEO | RNA | Autocrine vitamin D signaling switches off pro-inflammatory programs of Th1 cells | exclude | not blood |
| GSE154918 | GEO | RNA | RNAseq analysis of blood from sepsis patients and healthy controls | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE154998 | GEO | RNA | Transcriptional Profiling of Leukocytes in Critically Ill COVID19 Patients | exclude | 14 samples: cannot reach 10 participants per class |
| GSE155106 | GEO | RNA | Human dendritic cell and macrophage response to SARS-CoV-2 infection | exclude | 4 samples: cannot reach 10 participants per class |
| GSE155198 | GEO | RNA | Multi-omic data integration allows baseline immune signatures to predict hepatitis B vacci | exclude | vaccination study |
| GSE155223 | GEO | RNA | Single-Cell Omics Reveals Dyssynchrony of the Innate and Adaptive Immune System in Progres | exclude | single cell data |
| GSE155224 | GEO | RNA | Single-Cell Omics Reveals Dyssynchrony of the Innate and Adaptive Immune System in Progres | exclude | single cell data |
| GSE155241 | GEO | RNA | Identification of Drugs Blocking SARS-CoV-2 Infection using hPSC-derived Lung and Colonic  | exclude | in vitro or cell culture |
| GSE155249 | GEO | RNA | Circuits between infected macrophages and T cells in SARS-CoV-2 pneumonia | exclude | single cell data |
| GSE155454 | GEO | RNA | Whole blood transcriptome analysis reveals SARS-CoV-2 ORF8 as a potential therapeutic and  | exclude | 6 healthy samples of one shared id (fewer than 10 per class); no severity scale |
| GSE155635 | GEO | RNA | RNAseq analysis of blood from influenza infected patients and healthy controls-SIG2019 | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE155673 | GEO | RNA | Systems biological assessment of immunity to severe and mild COVID-19 infections | exclude | single cell data |
| GSE15573 | GEO | RNA | Immunity and Defense Genes in Peripheral Blood Mononuclear Cells of Rheumatoid Arthritis p | exclude | no characteristic splits 10 COVID samples from 10 others and no severity field (automatic; best : 0 against 0) |
| GSE155897 | GEO | RNA | Immunological Biomarkers of Fatal COVID-19: A Study of 868 Patients | exclude | no characteristic splits 10 COVID samples from 10 others and no severity field (automatic; best : 0 against 0) |
| GSE155986 | GEO | RNA | Transcriptomic profiling of human corona virus (HCoV)-229E -infected human cells | exclude | 15 samples: cannot reach 10 participants per class |
| GSE156063 | GEO | RNA | Upper airway gene expression differentiates COVID-19 from other acute respiratory illnesse | exclude | single cell data |
| GSE1563 | GEO | RNA | Kidney Transplant Rejection and Tissue Injury by Gene Profiling of Biopsies and Peripheral | exclude | no characteristic splits 10 COVID samples from 10 others and no severity field (automatic; best : 0 against 0) |
| GSE157007 | GEO | RNA | Multidimensional single-cell analysis of human peripheral blood reveals characteristic fea | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE157055 | GEO | RNA | Human lung organoid for modeling infection and disease conditions (Lung I) | exclude | in vitro or cell culture |
| GSE157057 | GEO | RNA | Human lung organoid for modeling infection and disease conditions (Lung II) | exclude | in vitro or cell culture |
| GSE157059 | GEO | RNA | AI-guided Discovery of a Paradoxical Immune Response in Fatal COVID-19 (Ileum) | exclude | not blood |
| GSE157240 | GEO | RNA | RNAseq analysis of blood from respiratory viruses-infected patients and healthy controls | exclude | no characteristic splits 10 COVID samples from 10 others and no severity field (automatic; best : 0 against 0) |
| GSE157344 | GEO | RNA | Deciphering the state of immune silence in fatal COVID-19 patients | exclude | single cell data |
| GSE157403 | GEO | RNA | Next Generation Sequencing Quantitative Analysis on COVID-19 Kidney tissue | exclude | 1 samples: cannot reach 10 participants per class |
| GSE157490 | GEO | RNA | A high-resolution temporal atlas of the SARS-CoV-2 translatome and transcriptome | exclude | no characteristic splits 10 COVID samples from 10 others and no severity field (automatic; best : 0 against 0) |
| GSE157789 | GEO | RNA | Dexamethasone modulates immature neutrophils and interferon programming in severe COVID-19 | exclude | single cell data |
| GSE15781 | GEO | RNA | New specific molecular targets for radiochemotherapy in colorectal cancer | exclude | no characteristic splits 10 COVID samples from 10 others and no severity field (automatic; best : 0 against 0) |
| GSE157827 | GEO | RNA | Single-nucleus transcriptome analysis reveals dysregulation of angiogenic endothelial cell | exclude | single cell data |
| GSE157852 | GEO | RNA | Human Pluripotent Stem Cell-Derived Neural Cells and Brain Organoids Reveal SARS-CoV-2 Neu | exclude | in vitro or cell culture |
| GSE157859 | GEO | RNA | Whole transcriptome sequencing of peripheral blood mononuclear cells from patients with CO | exclude | other assay or analyte |
| GSE158030 | GEO | RNA | In severe COVID-19, SARS-CoV-2 induces a chronic, TGF-β-dominated adaptive immune response | exclude | single cell data |
| GSE158034 | GEO | RNA | In severe COVID-19, SARS-CoV-2 induces a chronic, TGF-β-dominated adaptive immune response | exclude | single cell data |
| GSE158036 | GEO | RNA | In severe COVID-19, SARS-CoV-2 induces a chronic, TGF-β-dominated adaptive immune response | exclude | single cell data |
| GSE158055 | GEO | RNA | Large-scale single-cell analysis reveals critical immune characteristics of COVID-19 patie | exclude | single cell data |
| GSE158127 | GEO | RNA | Lung transplantation for patients with severe COVID-19 | exclude | single cell data |
| GSE158237 | GEO | RNA | Enterocyte Superoxide Dismutase 2 Deletion Drives  Hyperinsulinemia and Obesity | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE158374 | GEO | RNA | Genes with 5′ terminal oligopyrimidine tracts preferentially escape global suppression of  | exclude | no characteristic splits 10 COVID samples from 10 others and no severity field (automatic; best : 0 against 0) |
| GSE158484 | GEO | RNA | Functions of SARS-CoV-2 ORF3 during Virus Infection | exclude | in vitro or cell culture |
| GSE158558 | GEO | RNA | Intracranial aneurysm sequencing data | exclude | 8 samples: cannot reach 10 participants per class |
| GSE158592 | GEO | RNA | RNAseq analysis of blood from influenza infected patients and healthy controls SIG-2018 | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE158835 | GEO | RNA | MAIT cell activation augments adenovirus vector vaccine immunogenicity | exclude | not human |
| GSE158930 | GEO | RNA | Translational landscape of SARS-CoV-2 and infected cells | exclude | in vitro or cell culture |
| GSE158948 | GEO | RNA | A customizable purification assay designed for small-RNA biomarker identification and func | exclude | not human |
| GSE159191 | GEO | RNA | SARS-CoV-2 infected cells present HLA-I peptides from canonical and out-of-frame ORFs | exclude | in vitro or cell culture |
| GSE159201 | GEO | RNA | iPSC-derived Human Intestinal Organoids are permissive to SARS-CoV-2 infection | exclude | in vitro or cell culture |
| GSE159546 | GEO | RNA | Glucocorticoid receptor axis in lung cancer | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE159554 | GEO | RNA | Human Tubular Epithelial Cells Activate a Coordinated Stress Response after Serum Exposure | exclude | not human |
| GSE159556 | GEO | RNA | scRNA-seq analysis of SARS-CoV-2 infected human islets | exclude | single cell data |
| GSE159585 | GEO | RNA | The pulmonary vasculature in lethal COVID-19 and idiopathic pulmonary fibrosis at single c | exclude | single cell data |
| GSE159586 | GEO | RNA | Human Tubular Epithelial Cells Activate a Coordinated Stress Response after Serum Exposure | exclude | not human |
| GSE159678 | GEO | RNA | Hydroxychloroquine inhibits trained immunity | exclude | post-acute only (no acute samples by its record) |
| GSE159717 | GEO | RNA | SARS-CoV-2 infection of human pancreatic islets | exclude | single cell data |
| GSE159812 | GEO | RNA | Dysregulation of brain and choroid plexus cell types in severe COVID-19 | exclude | single cell data |
| GSE16008 | GEO | RNA | Gene expression from bronchial and nasal epithelial cell samples of healthy current and ne | exclude | not blood |
| GSE160230 | GEO | RNA | hPSC-derived Airway Organoids-based Screen Reveals the Role of HIF1/Glycolysis Axis in SAR | exclude | in vitro or cell culture |
| GSE160231 | GEO | RNA | hPSC-derived Airway Organoids-based Screen Reveals the Role of HIF1/Glycolysis Axis in SAR | exclude | single cell data |
| GSE160232 | GEO | RNA | hPSC-derived Airway Organoids-based Screen Reveals the Role of HIF1/Glycolysis Axis in SAR | exclude | in vitro or cell culture |
| GSE160351 | GEO | RNA | Macrophage expression and prognostic significance of the long pentraxin PTX3 in COVID-19 | exclude | single cell data |
| GSE1606 | GEO | RNA | Whole brain gene expression in X-monosomic and normal female mice | exclude | not human (Mus musculus) |
| GSE160631 | GEO | RNA | RNA-seq profiles on hESC-derived macrophage | exclude | no characteristic splits 10 COVID samples from 10 others and no severity field (automatic; best : 0 against 0) |
| GSE160864 | GEO | RNA | Discovery of a druggable copper-signaling pathway that drives inflammation | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE161200 | GEO | RNA | Temporal and dosage-dependent transcriptional response to SARS-CoV-2 in lung, trachea, bra | exclude | not human (Mesocricetus auratus) |
| GSE161225 | GEO | RNA | Next Generation Sequencing Analysis of  Transcriptomes in maculopapular drug rash with COV | exclude | not human |
| GSE161765 | GEO | RNA | Transcriptional checkpoints control Neutrophil-driven inflammation | exclude | not human (Mus musculus) |
| GSE161777 | GEO | RNA | Longitudinal multi-omics identifies responses of megakaryocytes, erythroid cells and plasm | exclude | single cell data |
| GSE161916 | GEO | RNA | Interferon-regulated genetic programs and JAK/STAT pathway activate the intronic promoter  | exclude | other assay or analyte |
| GSE161917 | GEO | RNA | JAK inhibitors dampen activation of interferon-activated transcriptomes and the SARS-CoV-2 | exclude | other assay or analyte |
| GSE161918 | GEO | RNA | Time-resolved Systems Immunology Reveals a Late Juncture Linked to Fatal COVID-19 | exclude | single cell data |
| GSE162113 | GEO | RNA | Metabolic reprogramming and epigenetic changes of vital organs in SARS-CoV-2 induced syste | exclude | not human |
| GSE162131 | GEO | RNA | Infection of human Nasal Epithelial Cells with SARS-CoV-2 and a 382-nt deletion isolate la | exclude | not blood |
| GSE162316 | GEO | RNA | A virus-derived microRNA targets immune response genes during SARS-CoV-2 infection [RNA-se | exclude | in vitro or cell culture |
| GSE162317 | GEO | RNA | A virus-derived microRNA targets immune response genes during SARS-CoV-2 infection [sRNA-s | exclude | in vitro or cell culture |
| GSE162318 | GEO | RNA | A virus-derived microRNA targets immune response genes during SARS-CoV-2 infection | exclude | other assay or analyte |
| GSE162323 | GEO | RNA | SARS-CoV-2 utilizes a multipronged strategy to suppress host protein synthesis | exclude | in vitro or cell culture |
| GSE162498 | GEO | RNA | Contribution of resident and circulating precursors to tumor-infiltrating CD8+ T cell popu | exclude | single cell data |
| GSE162632 | GEO | RNA | Single-cell RNA-sequencing reveals pervasive but highly cell type-specific genetic ancestr | exclude | single cell data |
| GSE162736 | GEO | RNA | SARS-CoV-2 infects human pluripotent stem cell-derived cardiomyocytes, impairing electrica | exclude | in vitro or cell culture |
| GSE162834 | GEO | RNA | Next-Generation Standard RNA-Sequencing Analysis of human monocytes stimulated by HMGB1 wi | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE162911 | GEO | RNA | A single-cell and spatial atlas of autopsy tissues reveals pathology and cellular targets  | exclude | single cell data |
| GSE163005 | GEO | RNA | Characterizing cerebrospinal fluid immunity in neurological manifestations of COVID-19 | exclude | single cell data |
| GSE163014 | GEO | RNA | Single-cell RNA sequencing of Tocilizumab-treated peripheral blood mononuclear cells as an | exclude | single cell data |
| GSE163151 | GEO | RNA | A diagnostic host response biosignature for COVID-19 from RNA profiling of nasal swabs and | exclude | whole blood subset: 7 COVID-19 samples (fewer than 10 per class) |
| GSE163317 | GEO | RNA | Spontaneous NLRP3 inflammasome-driven IL-1-β secretion is induced in severe COVID-19 patie | exclude | 8 samples: cannot reach 10 participants per class |
| GSE163426 | GEO | RNA | Tracheal aspirate RNA sequencing identifies distinct immunological features of COVID-19 AR | exclude | no characteristic splits 10 COVID samples from 10 others and no severity field (automatic; best ch_disease_state: 53 against 5) |
| GSE163668 | GEO | RNA | Global Absence and Targeting of Protective Immune States in Severe COVID-19 | exclude | single cell data |
| GSE163698 | GEO | RNA | Self-organized stem cell-derived human lung buds with proximo-distal patterning and novel  | exclude | single cell data |
| GSE16385 | GEO | RNA | Expression data from human macrophages | exclude | in vitro or cell culture |
| GSE163959 | GEO | RNA | Transcriptome analysis of Severe acute respiratory syndrome coronavirus 2(SARS-CoV-2)-infe | exclude | not blood |
| GSE164013 | GEO | RNA | A single-cell and spatial atlas of autopsy tissues reveals pathology and cellular targets  | exclude | single cell data |
| GSE164073 | GEO | RNA | Response to SARS-CoV-2 infection in cornea, limbus and sclera from human donors | exclude | vaccination study |
| GSE1643 | GEO | RNA | Normal Variance of Gene Expression in the Human Lung | exclude | single cell data |
| GSE164332 | GEO | RNA | Frontal Cortex whole transcriptome sequencing in COVID-19 elderly patients | exclude | single cell data |
| GSE164378 | GEO | RNA | Integrated analysis of multimodal single-cell data | exclude | single cell data |
| GSE164379 | GEO | RNA | Distinct B cell repertoires characterize patients with mild and severe COVID-19 [GEX] | exclude | single cell data |
| GSE164386 | GEO | RNA | Outcome of SARS-CoV-2 infection is linked to MAIT cell activation and cytotoxicity | exclude | 8 samples: cannot reach 10 participants per class |
| GSE164485 | GEO | RNA | The landscape of human brain immune response in patients with severe COVID-19 | exclude | single cell data |
| GSE164595 | GEO | RNA | The role of XIST in X-inactivation maintenance in B cells (RNA-seq) | exclude | single cell data |
| GSE164784 | GEO | RNA | Antihypertensive drug treatment and susceptibility to SARS-CoV-2 infection in human PSC-de | exclude | in vitro or cell culture |
| GSE164805 | GEO | RNA | The transcriptional profiles of severe COVID-19 | exclude | 15 samples: cannot reach 10 participants per class |
| GSE164890 | GEO | RNA | RNA-seq of seven altitude-adapted female speed skaters before and after exercise | exclude | 14 samples: cannot reach 10 participants per class |
| GSE164948 | GEO | RNA | Single-cell analysis identifies shared and distinct immune features of COVID-19, Influenza | exclude | single cell data |
| GSE165025 | GEO | RNA | BRD2 inhibition blocks SARS-CoV-2 infection by reducing transcription of the host cell rec | exclude | in vitro or cell culture |
| GSE165080 | GEO | RNA | Identification of distinct immune cell subsets associated with various clinical presentati | exclude | single cell data |
| GSE165082 | GEO | RNA | DNA methylation and expression profiles of whole blood in Parkinson’s disease (RNA-Seq) | exclude | other assay or analyte |
| GSE165182 | GEO | RNA | Single-cell RNA sequencing of PBMCs sampled from COVID-19 patients | exclude | single cell data |
| GSE165193 | GEO | RNA | Single-cell immunophenotyping of the fetal immune response to maternal SARS-CoV-2 infectio | exclude | single cell data |
| GSE165200 | GEO | RNA | Adult Stem Cell-derived Complete Lung Organoid Models Emulate Lung Disease in COVID-19 | exclude | in vitro or cell culture |
| GSE165328 | GEO | RNA | 25-Hydroxycholesterol Inhibits KSHV and EBV Infections and Activates Inflammatory Cytokine | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE165373 | GEO | RNA | Transcriptional analysis of SARS-CoV-2 specific CD4 T cells | exclude | 10 samples: cannot reach 10 participants per class |
| GSE16538 | GEO | RNA | Genome-wide gene expression profile analysis in pulmonary sarcoidosis | exclude | single cell data |
| GSE165747 | GEO | RNA | Modified Vaccinia Ankara Vector-Based Vaccine Protects Macaques from SARS-CoV-2 Infection, | exclude | not human (Macaca mulatta) |
| GSE165906 | GEO | RNA | Single-cell analysis of the epigenomic and transcriptional landscape of innate immunity to | exclude | single cell data |
| GSE165907 | GEO | RNA | Single-cell analysis of the epigenomic and transcriptional landscape of innate immunity to | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE166190 | GEO | RNA | Understanding COVID-19 in children and adults | exclude | not adults |
| GSE166209 | GEO | RNA | Pharmacological Activation of STING blocks SARS-CoV-2 infection | exclude | not human |
| GSE166253 | GEO | RNA | Transcriptomic characteristics and impaired immune function of patients who retest positiv | exclude | single cell data |
| GSE166489 | GEO | RNA | Immune dysregulation and autoreactivity correlate with disease severity in SARS-CoV-2-asso | exclude | single cell data |
| GSE166530 | GEO | RNA | Deciphering the host transcriptional response to SARS-CoV-2 infection | exclude | no characteristic splits 10 COVID samples from 10 others and no severity field (automatic; best title: 36 against 5) |
| GSE166651 | GEO | RNA | SARS-CoV-2 escapes CD8 T cell surveillance via mutations in MHC-I restricted epitopes [10x | exclude | 8 samples: cannot reach 10 participants per class |
| GSE166992 | GEO | RNA | Single Cell RNA-Seq of PBMC from COVID-19 patients and Control subjects | exclude | single cell data |
| GSE167026 | GEO | RNA | RNA-seq analysis of canine samples | exclude | not human (Canis lupus familiaris) |
| GSE167028 | GEO | RNA | A monocyte/dendritic cell molecular signature of SARS-CoV2-related multisystem  inflammato | exclude | single cell data |
| GSE167029 | GEO | RNA | A monocyte/dendritic cell molecular signature of SARS-CoV2-related multisystem  inflammato | exclude | single cell data |
| GSE167075 | GEO | RNA | m6A modification in SARS-CoV-2 virus regulates host cell innate immune response | exclude | in vitro or cell culture |
| GSE167131 | GEO | RNA | Host-virus chimeric events in SARS-CoV2 infected cells are infrequent and artifactual | exclude | in vitro or cell culture |
| GSE167158 | GEO | RNA | Transcriptomic profiling of PBMCs from healthy persons and persons recovered from COVID-19 | exclude | 6 samples: cannot reach 10 participants per class |
| GSE167159 | GEO | RNA | Single-cell RNA sequencing of PBMCs from convalescent COVID-19 patients compared to health | exclude | single cell data |
| GSE167250 | GEO | RNA | Glucocorticoid signaling in pancreatic islets modulates gene regulatory programs and genet | exclude | in vitro or cell culture |
| GSE16728 | GEO | RNA | Characterization of whole blood gene expression profiles in sickle-cell disease patients u | exclude | no characteristic splits 10 COVID samples from 10 others and no severity field (automatic; best : 0 against 0) |
| GSE167363 | GEO | RNA | Dynamic changes in human single cell transcriptional signatures during fatal sepsis | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE167410 | GEO | RNA | Histone H4 is increased in SARS-CoV-2-infected lung epithelial cells and linked to COVID-1 | exclude | 16 samples: cannot reach 10 participants per class |
| GSE167480 | GEO | RNA | SARS-CoV-2 spike protein interacts with and activates TLR4 | exclude | 6 samples: cannot reach 10 participants per class |
| GSE1675 | GEO | RNA | SHR and WKY rat adrenal glands | exclude | not human (Rattus norvegicus) |
| GSE167747 | GEO | RNA | SARS-CoV-2 infects the human kidney and drives fibrosis in kidney organoids | exclude | single cell data |
| GSE167930 | GEO | RNA | Genome-wide PBMCs Profiling Reveals Early Detection of Immunothrombotic Dysregulation in C | exclude | no characteristic splits 10 COVID samples from 10 others and no severity field (automatic; best : 0 against 0) |
| GSE168017 | GEO | RNA | Impaired antibacterial immune signaling and changes in the lung microbiome precedes second | exclude | single cell data |
| GSE168215 | GEO | RNA | Single cell RNA-seq analysis of bronchial epithelium in patients with severe SARS-CoV-2 pn | exclude | single cell data |
| GSE168388 | GEO | RNA | Maturation signatures of conventional dendritic cell subtypes in COVID-19 suggest direct v | exclude | in vitro or cell culture |
| GSE168400 | GEO | RNA | Peripheral blood gene expression according to IFIH1 rs1990760 variants in critically-ill C | exclude | no characteristic splits 10 COVID samples from 10 others and no severity field (automatic; best : 0 against 0) |
| GSE168453 | GEO | RNA | Longitudinal single-cell epitope and RNA-sequencing reveals the immunological impact of ty | exclude | single cell data |
| GSE168527 | GEO | RNA | Transcriptional profiling of human KIR+ and KIR- CD8 T cells by bulk RNA-seq | exclude | not human |
| GSE168710 | GEO | RNA | IFN-γ and TNF-α drive a CXCL10+ CCL2+ macrophage phenotype expanded in severe COVID-19 lun | exclude | single cell data |
| GSE16879 | GEO | RNA | Mucosal expression profiling in patients with inflammatory bowel disease before and after  | exclude | no characteristic splits 10 COVID samples from 10 others and no severity field (automatic; best : 0 against 0) |
| GSE168998 | GEO | RNA | Transcriptomic analysis revealed increased expression of genes involved in keratinization  | exclude | 6 samples: cannot reach 10 participants per class |
| GSE169241 | GEO | RNA | hPSC-derived cells to model Macrophage-mediated inflammation in COVID19 Hearts  | exclude | not blood |
| GSE169346 | GEO | RNA | scRNAseq of Antigen Presenting Cells in severe COVID-19 patients | exclude | single cell data |
| GSE169441 | GEO | RNA | Immune responses to Coronavirus disease-19 with Thymosin α 1 intervention: an observationa | exclude | single cell data |
| GSE169471 | GEO | RNA | Transcriptional profiling of lung cell populations in idiopathic pulmonary arterial hypert | exclude | single cell data |
| GSE169503 | GEO | RNA | Identification and characterization of a SARS-CoV-2 specific CD8 T cell response with immu | exclude | single cell data |
| GSE169504 | GEO | RNA | Systemic Tissue and Cellular Disruption from SARS-CoV-2 Infection revealed in COVID-19 Aut | exclude | not blood |
| GSE17025 | GEO | RNA | Gene Expression Analysis of Stage I Endometrial Cancers | exclude | no characteristic splits 10 COVID samples from 10 others and no severity field (automatic; best : 0 against 0) |
| GSE171037 | GEO | RNA | Single-cell RNA sequencing reveals ex vivo signatures of SARS-CoV-2-reactive T cells throu | exclude | single cell data |
| GSE171052 | GEO | RNA | Case Report: Deep Immune Profiling Links Neutrophil and Plasmablast Dysregulation to Micro | exclude | not adults |
| GSE171080 | GEO | RNA | Comparative Transcriptomic and Molecular Pathway Analyses of HL-CZ Human Pro-Monocytic Cel | exclude | in vitro or cell culture |
| GSE171130 | GEO | RNA | SUMOylation restrains the degradation of ACE2 through TOLLIP-mediated selective autophagy  | exclude | no characteristic splits 10 COVID samples from 10 others and no severity field (automatic; best : 0 against 0) |
| GSE171381 | GEO | RNA | SARS-CoV-2 infection in pregnancy is associated with robust inflammatory response at the m | exclude | single cell data |
| GSE171524 | GEO | RNA | Columbia University/NYP COVID-19 Lung Atlas | exclude | single cell data |
| GSE171555 | GEO | RNA | Single-cell RNA-seq of human peripheral blood mononuclear cells in COVID-19 | exclude | single cell data |
| GSE17156 | GEO | RNA | Gene expression signatures of symptomatic respiratory viral infection in adults | exclude | no characteristic splits 10 COVID samples from 10 others and no severity field (automatic; best : 0 against 0) |
| GSE171668 | GEO | RNA | COVID-19 tissue atlases reveal SARS-COV-2 pathology and cellular targets | exclude | single cell data |
| GSE171703 | GEO | RNA | SARS-CoV-2 infection-induced human B cells | exclude | single cell data |
| GSE171964 | GEO | RNA | CITE-seq of 6 patients treated with the Pfizer BioNTech BNT162b2 mRNA vaccine | exclude | single cell data |
| GSE171995 | GEO | RNA | SARS-CoV-2 infection in pregnancy is associated with robust inflammatory response at the m | exclude | not blood |
| GSE172348 | GEO | RNA | Thyroid Hormone effects on hiPSC-derived cardiomyocytes (hiPSC-CM). Identification of targ | exclude | in vitro or cell culture |
| GSE173317 | GEO | RNA | Blood transcriptomes of anti-SARS-CoV2 antibody  positive healthy individuals with prior   | exclude | vaccination study |
| GSE173488 | GEO | RNA | Next generation sequencing transcriptome analysis of patient-derived macrophages stimulate | exclude | in vitro or cell culture |
| GSE173590 | GEO | RNA | Longitudinal trajectories of mild/moderate SARS-CoV2 infection reveals key precursors and  | exclude | single cell data |
| GSE174072 | GEO | RNA | Multi-omic profiling reveals widespread dysregulation of innate immunity and hematopoiesis | exclude | single cell data |
| GSE174398 | GEO | RNA | Switching of OAS1 splicing isoforms mitigates SARS-CoV-2 infection | exclude | 6 samples: cannot reach 10 participants per class |
| GSE174621 | GEO | RNA | BNT162b2 vaccination enhances interferon-JAK-STAT-regulated antiviral programs in COVID-19 | exclude | vaccination study |
| GSE174657 | GEO | RNA | Implicating effector genes at COVID-19 GWAS loci using promoter-focused Capture-C in disea | exclude | 8 samples: cannot reach 10 participants per class |
| GSE174658 | GEO | RNA | Implicating effector genes at COVID-19 GWAS loci using promoter-focused Capture-C in disea | exclude | no characteristic splits 10 COVID samples from 10 others and no severity field (automatic; best : 0 against 0) |
| GSE174668 | GEO | RNA | A multi-omics investigation on the composition and physiological function of extracellular | exclude | other assay or analyte |
| GSE174745 | GEO | RNA | SARS-CoV-2 Infection Causes Dopaminergic Neuron Senescence | exclude | in vitro or cell culture |
| GSE174779 | GEO | RNA | Regulatory KIR+RA+ T cells accumulate with age and are highly activated during viral respi | exclude | no characteristic splits 10 COVID samples from 10 others and no severity field (automatic; best : 0 against 0) |
| GSE174786 | GEO | RNA | RNAseq analysis of PS+ vs PS- non-naive CD8 T cells sorted from peripheral blood of COVID- | exclude | 8 samples: cannot reach 10 participants per class |
| GSE174825 | GEO | RNA | Molecular profiling reveals human biomarker signatures in patients after allogeneic stem c | exclude | no characteristic splits 10 COVID samples from 10 others and no severity field (automatic; best : 0 against 0) |
| GSE175450 | GEO | RNA | Complement activation induces excessive T cell cytotoxicity in severe COVID-19 | exclude | single cell data |
| GSE175880 | GEO | RNA | Identification and characterization of human CD34+ and CD34dim/- neutrophil-committed prog | exclude | single cell data |
| GSE176201 | GEO | RNA | Immune signatures underlying post-acute COVID-19 lung sequelae | exclude | single cell data |
| GSE17624 | GEO | RNA | Expression data from human Ishikawa cells treated with Bisphenol A | exclude | in vitro or cell culture |
| GSE176269 | GEO | RNA | Human nasal epithelial and immune cell responses to SARS-CoV-2 versus influenza A virus | exclude | single cell data |
| GSE176290 | GEO | RNA | RNA sequencing of blood monocyte subsets in COVID patients | exclude | other assay or analyte |
| GSE176442 | GEO | RNA | Impaired antigen-specific memory B cell and plasma cell responses including lack of specif | exclude | single cell data |
| GSE176479 | GEO | RNA | Microvascular Endothelial Cell RNA-Seq: Treatment with Platelet Releasate from COVID-19 &  | exclude | 14 samples: cannot reach 10 participants per class |
| GSE176480 | GEO | RNA | Assessment of the platelet transcriptome in patients with COVID-19 versus controls | exclude | 18 samples: cannot reach 10 participants per class |
| GSE176498 | GEO | RNA | Cardiometabolic MicroRNAs are Associated with severeerity and 28-Day Intensive Care Unit M | exclude | other assay or analyte |
| GSE178123 | GEO | RNA | Single Cell RNAseq and VDJ sequence of ACE2 reactive and non-reactive B-cells isolated fro | exclude | single cell data |
| GSE178331 | GEO | RNA | mRNA profiling in pooled human umbilical vein endothelial cells stimulated with Coronaviru | exclude | in vitro or cell culture |
| GSE178388 | GEO | RNA | Downregulation of exhausted cytotoxic T cells in gene expression networks of multisystem i | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE178404 | GEO | RNA | Identification of immune correlates of fatal outcomes in critically ill COVID-19 patients | exclude | single cell data |
| GSE178491 | GEO | RNA | An AI-guided invariant signature places MIS-C with Kawasaki disease in a continuum of host | exclude | not adults |
| GSE178824 | GEO | RNA | Severe COVID-19 is Characterized by an Impaired Type I Interferon Response and Elevated Le | exclude | other assay or analyte |
| GSE178942 | GEO | RNA | SARS-CoV-2 infected and uninfected VeroE6 and A549-ACE2: Bemcentinib and liposome treatmen | exclude | in vitro or cell culture |
| GSE179156 | GEO | RNA | Up-regulation of ACE2, the SARS-CoV-2 Receptor, in Asthmatics on Maintenance Inhaled Corti | exclude | RNA from plasma or serum (cell-free RNA), not whole blood or PBMC |
| GSE179184 | GEO | RNA | Mapping the 3D genome and epigenome in SARS-CoV-2 infected human cells | exclude | no characteristic splits 10 COVID samples from 10 others and no severity field (automatic; best : 0 against 0) |
| GSE179251 | GEO | RNA | SARS-CoV-2 Nsp14 mediates the effects of viral infection on the host cell transcriptome | exclude | no characteristic splits 10 COVID samples from 10 others and no severity field (automatic; best : 0 against 0) |
| GSE179277 | GEO | RNA | Upper airway transcriptional profiling reveals a more robust innate and adaptive immune re | exclude | not adults |
| GSE179347 | GEO | RNA | An integrated approach to identify environmental modulators of genetic risk factors for co | exclude | in vitro or cell culture |
| GSE179448 | GEO | RNA | RNAseq analysis of human CD4+ regulatory Tregs and Tconvs in COVID-19 and healthy donors i | exclude | post-acute only (no acute samples by its record) |
| GSE179548 | GEO | RNA | Activation of FXR by bile acids in hepatic stellate cells plays a critical role in liver c | exclude | single cell data |
| GSE179923 | GEO | RNA | SARS-CoV-2 infects the brain vasculature in COVID-19 patients and a hiPSC-derived model of | exclude | in vitro or cell culture |
| GSE180226 | GEO | RNA | Pulmonary epithelial and endothelial cell damage, loss of repair, inhibition of fibrinolys | exclude | RNA from plasma or serum (cell-free RNA), not whole blood or PBMC |
| GSE180268 | GEO | RNA | Functional HPV-specific PD-1+ stem-like CD8 T cells in head and neck cancer | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE180298 | GEO | RNA | Uncovering perturbations in human hematopoiesis associated with healthy aging and myeloid  | exclude | single cell data |
| GSE180389 | GEO | RNA | Probiotic consortia improve anti-viral immunity to SARS-CoV-2 in Ferrets [RNA-seq] | exclude | not human |
| GSE180391 | GEO | RNA | Probiotic consortia improve anti-viral immunity to SARS-CoV-2 in Ferrets | exclude | not human |
| GSE180578 | GEO | RNA | People critically ill with COVID-19 exhibit peripheral immune profiles predictive of morta | exclude | single cell data |
| GSE180594 | GEO | RNA | Alterations in circulating monocytes predict COVID-19 severity and include chromatin modif | exclude | no characteristic splits 10 COVID samples from 10 others and no severity field (automatic; best ch_infection_status: 26 against 7) |
| GSE180622 | GEO | RNA | TWAS-based translational genomics approach identifies IL10RB as the top candidate gene for | exclude | in vitro or cell culture |
| GSE18090 | GEO | RNA | Gene Expression Profiling During Early Acute Febrile Stage of Dengue Infection Can Predict | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE181005 | GEO | RNA | Next Generation Sequencing of sorted central memory CD4 T cells (TCM) from young and elder | exclude | 12 samples: cannot reach 10 participants per class |
| GSE181032 | GEO | RNA | Convalescent COVID-19 patients without comorbidities display similar immunophenotypes over | exclude | 15 samples: cannot reach 10 participants per class |
| GSE181046 | GEO | RNA | Protective heterologous T cell immunity in COVID-19 induced by MMR and Tdap vaccine antige | exclude | single cell data |
| GSE181238 | GEO | RNA | Effect of SARS-CoV-2 on the human placenta | exclude | not blood |
| GSE181439 | GEO | RNA | Upregulation of cytokine signalling in platelet increases risk of thrombophilia in severe  | exclude | RNA from plasma or serum (cell-free RNA), not whole blood or PBMC |
| GSE181787 | GEO | RNA | X-linked recessive TLR7 deficiency in 1% of men under 60 years  with life-threatening COVI | exclude | in vitro or cell culture |
| GSE181827 | GEO | RNA | Short duplex RNAs induce potent immunostimulation via  Hoogsteen G-G base pairing-mediated | exclude | in vitro or cell culture |
| GSE181897 | GEO | RNA | Optimal pooled single-cell RNA-sequencing elucidates gene by environment interactions in c | exclude | single cell data |
| GSE182123 | GEO | RNA | Clinical and immunological signatures of severe COVID-19 in previously healthy patients wi | exclude | single cell data |
| GSE182297 | GEO | RNA | Molecular Profiling of COVID-19 Autopsies Uncovers Novel Disease Mechanisms [RNA-Seq] | exclude | single cell data |
| GSE182572 | GEO | RNA | TLR7 sensing of SARS-CoV-2 and IFN-I response by plasmacytoid dendritic cells leads to the | exclude | no characteristic splits 10 COVID samples from 10 others and no severity field (automatic; best : 0 against 0) |
| GSE182917 | GEO | RNA | Postmortem high-dimensional immune profiling of severe COVID-19 patients reveals distinct  | exclude | not blood |
| GSE183276 | GEO | RNA | An atlas of healthy and injured cell states and niches in the human kidney [scCv3] | exclude | single cell data |
| GSE18344 | GEO | RNA | The Transcriptome of Nrf2-Deficient Mice in Cigarette Smoke-Induced Emphysematous Changes | exclude | not human (Mus musculus) |
| GSE183716 | GEO | RNA | The innate and adaptive immune landscape of SARS-CoV-2-associated multisystem inflammatory | exclude | not adults |
| GSE183974 | GEO | RNA | Single-cell RNAseq analysis (10X Genomics Chromium) of residual bronchoalveolar lavage flu | exclude | single cell data |
| GSE184329 | GEO | RNA | Untimely TGFβ responses in COVID-19 limit antiviral functions of NK cells | exclude | single cell data |
| GSE184330 | GEO | RNA | Single cell RNA sequencing and TCR repertoire analysis of MIS-C affected patients versus h | exclude | single cell data |
| GSE184536 | GEO | RNA | The NF-κB transcriptional footprint is essential for SARS-CoV-2 replication | exclude | single cell data |
| GSE18454 | GEO | RNA | Analysis of normal lung cells treated with 5-aza-dC to induce DNA demethylation | exclude | in vitro or cell culture |
| GSE184569 | GEO | RNA | Enteric coronavirus infection and treatment modeled with an immunocompetent human intestin | exclude | in vitro or cell culture |
| GSE184806 | GEO | RNA | Investigation of human peripheral T cells responsive to SARS-CoV-2 antigens | exclude | single cell data |
| GSE184932 | GEO | RNA | Profilling of viral and host transcriptome in SARS-CoV-2 infected A549-ACE2 cells with VIC | exclude | in vitro or cell culture |
| GSE184949 | GEO | RNA | Similarities and differences in the acute-phase response to SARS-CoV-2 in rhesus macaques  | exclude | not human (Macaca mulatta; Chlorocebus sabaeus) |
| GSE185432 | GEO | RNA | Manufacture of GMP-compliant functional SARS-COV-2 cytotoxic T lymphocytes (CTLs) utilizin | exclude | 8 samples: cannot reach 10 participants per class |
| GSE185863 | GEO | RNA | Transcriptomic profiles revealed down regulation of low-density lipoprotein particle recep | exclude | 16 samples: cannot reach 10 participants per class |
| GSE186020 | GEO | RNA | Validate whether the targeting of crRNAs affected the host cell gene expression | exclude | in vitro or cell culture |
| GSE186213 | GEO | RNA | Spatial transcriptomic characterization of the lung parenchyma during COVID-19 pneumonitis | exclude | not blood |
| GSE186267 | GEO | RNA | Transcriptional reprogramming of infiltrating neutrophils drives lung disease in severe CO | exclude | single cell data |
| GSE186460 | GEO | RNA | Baricitinib attenuates the proinflammatory phase1 of COVID-19 driven by lung-infiltrating  | exclude | no characteristic splits 10 COVID samples from 10 others and no severity field (automatic; best title: 10 against 8) |
| GSE188172 | GEO | RNA | PBMC single cell transcriptomics of COVID-19 patients | exclude | single cell data |
| GSE18842 | GEO | RNA | Gene expression analysis of human lung cancer and control samples | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE188429 | GEO | RNA | Allelic variation in Class I HLA determines pre-existing memory responses to SARS-CoV-2 th | exclude | single cell data |
| GSE188477 | GEO | RNA | Variant specific global loss of cellular m6A RNA methylation during SARS-CoV-2 infection | exclude | vaccination study |
| GSE188573 | GEO | RNA | Epigenetic and transcriptomic remodeling in monocytes of severe COVID-19 patients mirrors  | exclude | single cell data |
| GSE188678 | GEO | RNA | A 2-Gene Host Signature for Improved Accuracy of COVID-19 Diagnosis Agnostic to Viral Vari | exclude | no characteristic splits 10 COVID samples from 10 others and no severity field (automatic; best : 0 against 0) |
| GSE188945 | GEO | RNA | A new tractable method for generating Human Alveolar Macrophage Like cells in vitro to stu | exclude | in vitro or cell culture |
| GSE189039 | GEO | RNA | Prior BNT162b vaccination enhances interferon-JAK-STAT regulated antiviral programs in hos | exclude | vaccination study |
| GSE189199 | GEO | RNA | Vascular Damage, Thromboinflammation, Plasmablast Activation, T-Cell Dysregulation and Pat | exclude | not blood |
| GSE189263 | GEO | RNA | Transcriptomes of SARS-CoV-2 naïve and recovered individuals vaccinated with inactivated v | exclude | vaccination study |
| GSE18948 | GEO | RNA | Personalized medicine in psoriasis: developing a genomic classifier to predict histologica | exclude | 16 samples: cannot reach 10 participants per class |
| GSE189506 | GEO | RNA | Small RNA-seq analysis of serum miRNAs in patients with severe COVID-19 | exclude | other assay or analyte |
| GSE189900 | GEO | RNA | Unambiguous detection of SARS-CoV-2 subgenomic mRNAs with single cell RNA sequencing | exclude | single cell data |
| GSE189990 | GEO | RNA | RNA Sequencing in COVID-19 patients identifies neutrophil activation biomarkers as promisi | exclude | other assay or analyte |
| GSE190001 | GEO | RNA | Detailed delineation of the immune response to COVID-19 mRNA vaccination via high-temporal | exclude | vaccination study |
| GSE190067 | GEO | RNA | Transcriptomic characterisation of haematopoietic stem and progenitor cells from adult hum | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE190494 | GEO | RNA | Nanostring GeoMx digital spatial profliing of myeloid cells from COVID-19 and normal lung | exclude | not blood |
| GSE190496 | GEO | RNA | Gene expression profiling of FFPE normal and COVID19 lung tissues | exclude | not blood |
| GSE190680 | GEO | RNA | Immune transcriptome analysis of COVID-19 patients infected with SARS-CoV-2 variants carry | exclude | vaccination study |
| GSE190747 | GEO | RNA | mRNA vaccination in octogenarian nuns 15 and 20 months after recovery from COVID-19 elicit | exclude | vaccination study |
| GSE190839 | GEO | RNA | Recruitment of highly cytotoxic CD8+ T cell receptors in mild SARS-CoV-2 infection | exclude | in vitro or cell culture |
| GSE190895 | GEO | RNA | Investigation of SARS-CoV-2 variants that escape from cellular immunity in the absence of  | exclude | 6 samples: cannot reach 10 participants per class |
| GSE191088 | GEO | RNA | Comprehensive comparison of adaptive immune response by inactivated SARS-CoV-2 vaccine bet | exclude | vaccination study |
| GSE191089 | GEO | RNA | Comprehensive comparison of adaptive immune response by inactivated SARS-CoV-2 vaccine bet | exclude | single cell data |
| GSE191232 | GEO | RNA | Conjunctival epithelial cells resist productive SARS-CoV-2 infection | exclude | single cell data |
| GSE191328 | GEO | RNA | Longitudinal multi-omics analysis identifies early blood-based predictors of anti-TNF ther | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE19151 | GEO | RNA | Whole Blood Gene Expression Profiles Distinguish Patients with Single versus Recurrent Ven | exclude | no characteristic splits 10 COVID samples from 10 others and no severity field (automatic; best : 0 against 0) |
| GSE19190 | GEO | RNA | Distinct epithelial gene expression phenotypes in childhood respiratory allergy | exclude | single cell data |
| GSE192391 | GEO | RNA | Single cell transcriptomics of PBMCs during severe COVID-19 | exclude | single cell data |
| GSE192483 | GEO | RNA | Single cell RNA sequencing analysis of fresh resected lung tissues from tuberculosis patie | exclude | single cell data |
| GSE192656 | GEO | RNA | What goes around comes around: artificial circular RNAs bypass cellular antiviral response | exclude | not human |
| GSE192707 | GEO | RNA | Resistin Forms a Network with Cytokines Associated with Endothelial Damage and Worse Outco | exclude | 15 samples: cannot reach 10 participants per class |
| GSE19272 | GEO | RNA | Genome-wide comparison between IL-17 and combined TNF-alpha/ IL-17 induced genes in primar | exclude | not human (Mus musculus) |
| GSE193096 | GEO | RNA | IL-4 receptor blockade is a global repressor of naïve B cell development and responses in  | exclude | single cell data |
| GSE19314 | GEO | RNA | Sarcoidosis-specific markers from whole blood gene expression | exclude | no characteristic splits 10 COVID samples from 10 others and no severity field (automatic; best : 0 against 0) |
| GSE193169 | GEO | RNA | SARS-CoV-2 diverges from other betacoronaviruses in only partially activating the IRE1a/XB | exclude | not human |
| GSE193330 | GEO | RNA | Elucidation of the liver pathophysiology of COVID-19 patients using liver-on-a-chips | exclude | 12 samples: cannot reach 10 participants per class |
| GSE19339 | GEO | RNA | Differential gene expression in thrombus-derived white blood cells of patients with acute  | exclude | 8 samples: cannot reach 10 participants per class |
| GSE193722 | GEO | RNA | SARS-CoV-2 Infection Induces Ferroptosis of Human Pluripotent Stem Cell-Derived Sinoatrial | exclude | no characteristic splits 10 COVID samples from 10 others and no severity field (automatic; best ch_infection: 15 against 6) |
| GSE193732 | GEO | RNA | Cap-independent translation and a precisely localized RNA sequence enable SARS-CoV-2 to co | exclude | no characteristic splits 10 COVID samples from 10 others and no severity field (automatic; best : 0 against 0) |
| GSE193770 | GEO | RNA | Single-cell transcriptomic analysis of blood T cells from healthy controls, MS patients an | exclude | single cell data |
| GSE193782 | GEO | RNA | ScRNA-seq Expression of APOC2 and IFI27 Identifies Four Alveolar Macrophage Superclusters  | exclude | single cell data |
| GSE194061 | GEO | RNA | Limited extent and consequences of pancreatic SARS-CoV-2 infection | exclude | single cell data |
| GSE194331 | GEO | RNA | Gene Expression Profiling of Acute Pancreatitis | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE194378 | GEO | RNA | Influenza vaccination reveals sex dimorphic imprints of prior mild COVID-19 [RNA-Seq] | exclude | single cell data |
| GSE19444 | GEO | RNA | Blood Transcriptional Profiles of Active and Latent TB (UK Test Set) | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE19510 | GEO | RNA | Transcriptional response of normal human lung WI-38 fibroblasts to benzo[a]pyrene diol epo | exclude | in vitro or cell culture |
| GSE195541 | GEO | RNA | Type I Interferon Transcriptional Network Regulates Expression of Coinhibitory Receptors i | exclude | in vitro or cell culture |
| GSE195542 | GEO | RNA | Type I Interferon Transcriptional Network Regulates Expression of Coinhibitory Receptors i | exclude | in vitro or cell culture |
| GSE19567 | GEO | RNA | Imatinib and Nilotinib induced molecular changes in Philadelphia chromosome-positive Chron | exclude | 12 samples: cannot reach 10 participants per class |
| GSE195673 | GEO | RNA | Germinal centre-driven maturation of B cell response to SARS-CoV-2 mRNA vaccination | exclude | vaccination study |
| GSE195938 | GEO | RNA | Covid-19 Infection During Pregnancy Induces Differential Gene Expression in Human Cord Blo | exclude | not adults |
| GSE196113 | GEO | RNA | SARS-CoV-2 disrupts respiratory vascular barriers by suppressing Claudin-5 expression | exclude | 6 samples: cannot reach 10 participants per class |
| GSE196198 | GEO | RNA | Longitudinal single-cell multiomics data for different severity of COVID-19 | exclude | single cell data |
| GSE196350 | GEO | RNA | RNAseq analysis of blood from influenza infected patients and healthy controls | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE196455 | GEO | RNA | Viral mimicry of Interleukin-17A by SARS-COV-2 ORF8 | exclude | 12 samples: cannot reach 10 participants per class |
| GSE196456 | GEO | RNA | Clonal dynamics of SARS-CoV-2-specific T cells in children and adults with COVID-19 | exclude | not adults |
| GSE196802 | GEO | RNA | Transcriptome analysis of umbilical cord blood CD34+ cells of neonates delivered from SARS | exclude | vaccination study |
| GSE19681 | GEO | RNA | Global gene expression profile of primary human leukemic blasts | exclude | 10 samples: cannot reach 10 participants per class |
| GSE196820 | GEO | RNA | A public antibody class recognizes an S2 epitope exposed on open conformations of SARS-CoV | exclude | single cell data |
| GSE196990 | GEO | RNA | Durable Alterations of Hematopoiesis and Chromatin post-COVID-19 | exclude | post-acute only (no acute samples by its record) |
| GSE197204 | GEO | RNA | Transcriptomic clustering of critically-ill COVID-19 patients [Day 1-RNA-Seq] | exclude | no characteristic splits 10 COVID samples from 10 others and no severity field (automatic; best : 0 against 0) |
| GSE197259 | GEO | RNA | Transcriptomic clustering of critically-ill COVID-19 patients | exclude | no characteristic splits 10 COVID samples from 10 others and no severity field (automatic; best : 0 against 0) |
| GSE197266 | GEO | RNA | Single-cell sequencing reveals homogeneity and heterogeneity in the molecular cytopatholog | exclude | not human (Mus musculus) |
| GSE19743 | GEO | RNA | A large-scale clinical study of gene expression response to severe burn injury | exclude | not adults |
| GSE197440 | GEO | RNA | A Novel MIP-1 Expressing Macrophage Subtype in Bronchoalveolar Lavage Fluid from Healthy V | exclude | single cell data |
| GSE197665 | GEO | RNA | Exposure of PBMCs to SARS-CoV-2, as opposed to SARS-CoV, induces cell-intrinsic type I inn | exclude | 4 samples: cannot reach 10 participants per class |
| GSE197976 | GEO | RNA | Impact of an Immune Modulator Mycobacterium-w on Adaptive Natural Killer Cells and Protect | exclude | vaccination study |
| GSE198281 | GEO | RNA | Proliferative Receptor Binding Domain (RBD)-specific CD4+ T cells correlate with neutraliz | exclude | vaccination study |
| GSE198486 | GEO | RNA | Total RNA expression in extracellular vesicles (EVs) derived from untreated and lipid nano | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE198533 | GEO | RNA | mRNA expression profiling of PBMCs from Behcet’s disease | exclude | single cell data |
| GSE198549 | GEO | RNA | Immune transcriptome in adult-onset Still’s disease with mild flare following administrati | exclude | vaccination study |
| GSE199272 | GEO | RNA | Transcriptome Analaysis of humanized lungs in COVID-19 [scRNAseq] | exclude | single cell data |
| GSE199286 | GEO | RNA | Monocytic-Myeloid Derived Suppressor Cells Suppress T-Cell Responses in Recovered SARS CoV | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE199337 | GEO | RNA | Single-cell transcriptomics reveal a unique memory-like NK cell subset that accumulates wi | exclude | single cell data |
| GSE199750 | GEO | RNA | COVID-19 Vaccine Immune Responses Study (COVIRS) | exclude | not human |
| GSE19976 | GEO | RNA | Gene expression analysis of lung biopsies from patients with two different forms of pulmon | exclude | not blood |
| GSE199816 | GEO | RNA | Transcriptional Differences Between COVID-19 and Sepsis | exclude | in vitro or cell culture |
| GSE199911 | GEO | RNA | Gene dysregulation in acute HIV-1 infection – early transcriptomic analysis reveals the cr | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE200046 | GEO | RNA | SEACells: Inference of transcriptional and epigenomic cellular states from single-cell gen | exclude | single cell data |
| GSE200274 | GEO | RNA | Next generation sequencing transcriptome analysis of patient-derived macrophages after SAR | exclude | vaccination study |
| GSE200561 | GEO | RNA | Transcriptome Analaysis of humanized lungs in COVID-19 [Bulk_RNAseq] | exclude | not human (Mus musculus) |
| GSE200596 | GEO | RNA | A live attenuated vaccine confers superior mucosal and systemic immunity to SARS-CoV-2 var | exclude | not human (Mesocricetus auratus) |
| GSE200849 | GEO | RNA | Tracheal aspirate bulk and single cell RNA sequencing from patients with ARDS | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE200896 | GEO | RNA | Biologic correlates of COVID-19 resolution using single cell RNA sequecing | exclude | single cell data |
| GSE201088 | GEO | RNA | Single Cell Multiomics of Healthy, COVID-19 positive and Recovered Individuals | exclude | single cell data |
| GSE20114 | GEO | RNA | DHA supplementation suppresses LDLR and OLR1 gene expression in blood cells from hypertrig | exclude | 15 samples: cannot reach 10 participants per class |
| GSE201197 | GEO | RNA | Single-cell multi-omics analyses in the DILfrequency clinical trail samples | exclude | single cell data |
| GSE201325 | GEO | RNA | Whole transcriptome exploration revealed the dysregulated RNAs involved in SARS-CoV-2 spik | exclude | in vitro or cell culture |
| GSE201333 | GEO | RNA | Tabula Sapiens | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE201463 | GEO | RNA | Antigenic determinants of SARS-CoV-2-specific CD4+ T cell lines reveals M protein-driven d | exclude | single cell data |
| GSE201523 | GEO | RNA | Trisomy 21 human iPSC derived interneuron progenitors | exclude | single cell data |
| GSE201530 | GEO | RNA | Prior vaccination exceeds prior infection in eliciting innate and humoral immune responses | exclude | vaccination study |
| GSE201533 | GEO | RNA | Heterologous ChAdOx1-BNT162b2 vaccination in a Korean cohort induces a robust genetic immu | exclude | single cell data |
| GSE201534 | GEO | RNA | Heterologous ChAdOx1-BNT162b2 vaccination in a Korean cohort induces a robust genetic immu | exclude | single cell data |
| GSE201535 | GEO | RNA | Heterologous ChAdOx1-BNT162b2 vaccination in Korean cohort induces robust immune and antib | exclude | single cell data |
| GSE201642 | GEO | RNA | Impact of two and three dose COVID-19 vaccination schedule on transcriptional, humoral and | exclude | vaccination study |
| GSE20211 | GEO | RNA | LMP-420: a novel purine nucleoside analogue with potent cytotoxic effects for chronic lymp | exclude | in vitro or cell culture |
| GSE202182 | GEO | RNA | A specific molecular signature in SARS-CoV-2 infected kidney biopsies | exclude | not blood |
| GSE202262 | GEO | RNA | Antigen-specific CD8 T cells against SARS-CoV-2 | exclude | vaccination study |
| GSE202805 | GEO | RNA | Global Transcriptomic analysis of PBMC derived from COVID-19 acute-severe patients and con | exclude | single cell data |
| GSE203001 | GEO | RNA | Transcriptional Profiling of SARS-CoV-2 and IAV across Syrian golden hamster brain regions | exclude | not human (Mesocricetus auratus) |
| GSE203081 | GEO | RNA | Single-Cell RNA Sequencing to Decipher the Immunogenicity of ChAdOx1 nCoV-19/AZD1222 and m | exclude | single cell data |
| GSE203220 | GEO | RNA | Human IRF1 governs macrophagic IFN-γ immunity to mycobacteria | exclude | in vitro or cell culture |
| GSE20346 | GEO | RNA | Expression profiling of critically ill influenza and bacterial pneumonia patients, also in | exclude | vaccination study |
| GSE204962 | GEO | RNA | A metabolic associated fatty liver disease risk variant in MBOAT7 regulates toll like rece | exclude | 12 samples: cannot reach 10 participants per class |
| GSE205014 | GEO | RNA | Cytokines and lipid mediators of inflammation in lungs of SARS-CoV-2 infected mice | exclude | not human (Mus musculus) |
| GSE205244 | GEO | RNA | Limited cross-variant immune response from SARS-CoV-2 Omicron BA.2 in naïve but not previo | exclude | post-acute only (no acute samples by its record) |
| GSE205402 | GEO | RNA | SARS-CoV-2 mRNA Dual Immunization Induces Innate Transcriptional Signatures, Establishes T | exclude | vaccination study |
| GSE205487 | GEO | RNA | Pro-inflammatory activity of SARS-CoV-2 viral peptide-dsRNA complex on human dermal microv | exclude | 6 samples: cannot reach 10 participants per class |
| GSE205606 | GEO | RNA | Single-cell RNA-seq analysis identifies distinct myeloid cells in COVID19 vaccination- ass | exclude | single cell data |
| GSE205607 | GEO | RNA | Single-cell RNA-seq analysis identifies distinct myeloid cells in COVID19 vaccination- ass | exclude | single cell data |
| GSE206023 | GEO | RNA | Transcriptomes analysis of COVID-19 Inactivated Vaccine and Protein Subunit Vaccine as Thi | exclude | vaccination study |
| GSE206134 | GEO | RNA | Augmenting SARS-CoV-2 vaccine induced spike protein-specific CD4+ T cell responses through | exclude | single cell data |
| GSE206147 | GEO | RNA | Single-cell RNA sequencing analysis of a COVID-19– associated maculopapular rash in a pati | exclude | single cell data |
| GSE206263 | GEO | RNA | Differential chromatin accessibility in peripheral blood mononuclear cells underlies COVID | exclude | single cell data |
| GSE206264 | GEO | RNA | Transcriptomic clustering of critically-ill COVID-19 patients [Day 4 RNAseq] | exclude | no characteristic splits 10 COVID samples from 10 others and no severity field (automatic; best : 0 against 0) |
| GSE206265 | GEO | RNA | Influenza vaccination reveals sex dimorphic imprints of prior mild COVID-19 [single-cell C | exclude | single cell data |
| GSE206283 | GEO | RNA | Differential chromatin accessibility in peripheral blood mononuclear cells underlies COVID | exclude | single cell data |
| GSE206284 | GEO | RNA | Differential chromatin accessibility in peripheral blood mononuclear cells underlies COVID | exclude | single cell data |
| GSE206289 | GEO | RNA | SARS-CoV-2 induced immune perturbations in infants vary with disease severity and differ f | exclude | single cell data |
| GSE206510 | GEO | RNA | IL-13 stimulation of primary human bronchial epithelial cells in control individuals and p | exclude | in vitro or cell culture |
| GSE206635 | GEO | RNA | SARS-CoV-2 disrupt respiratory vascular barrier by suppressing Claudin-5 expression (addit | exclude | 5 samples: cannot reach 10 participants per class |
| GSE207077 | GEO | RNA | β-hydroxybutyrate restores T cell metabolism and function in severe COVID-19 | exclude | not human |
| GSE207151 | GEO | RNA | HCV infection increases the expression of ACE2 receptor, leading to enhanced entry of both | exclude | 6 samples: cannot reach 10 participants per class |
| GSE207270 | GEO | RNA | Switching of OAS1 splicing isoforms mitigates SARS-CoV-2 infection [Caco-2] | exclude | in vitro or cell culture |
| GSE207473 | GEO | RNA | Serological mechanisms driving cardiotoxicity in COVID-19 | exclude | RNA from plasma or serum (cell-free RNA), not whole blood or PBMC |
| GSE207625 | GEO | RNA | RNA Sequencing Analysis of Peripheral Blood Mononuclear Cells Associated with COVID-19 Det | exclude | 12 samples: cannot reach 10 participants per class |
| GSE207840 | GEO | RNA | Transcriptional effects of SARS-CoV-2 on human  M-CSF-derived macrophages | exclude | no characteristic splits 10 COVID samples from 10 others and no severity field (automatic; best : 0 against 0) |
| GSE207981 | GEO | RNA | Identification of host dependency factors shared by multiple SARS-CoV-2 variants of concer | exclude | vaccination study |
| GSE208076 | GEO | RNA | Genome wide expression from COVID-19 patients and controls' lung | exclude | not blood |
| GSE208116 | GEO | RNA | Evasion of antiviral innate immunity by SARS-CoV-2 NSP1 is exerted through H3K9 di-methyla | exclude | in vitro or cell culture |
| GSE208320 | GEO | RNA | Early suppression of antiviral host response and protocadherins by SARS-CoV-2 Spike protei | exclude | vaccination study |
| GSE208337 | GEO | RNA | Elevated Myl9 reflects the Myl9-containing microthrombi in SARS-CoV-2-induced exudative va | exclude | single cell data |
| GSE208342 | GEO | RNA | Transcriptome and microbiome profile of nasopharyngeal samples among COVID-19 patients wit | exclude | not blood |
| GSE208684 | GEO | RNA | SARS-COV2 in intestinal organoids | exclude | in vitro or cell culture |
| GSE209676 | GEO | RNA | Analysis on SARS-CoV-2-M198-206-specific CD8+ T cells | exclude | single cell data |
| GSE209985 | GEO | RNA | Response to mRNA vaccine BNT162b2 in hemodialysis patients | exclude | vaccination study |
| GSE210229 | GEO | RNA | CD8+ T-cell memory induced by successive SARS-CoV-2 mRNA vaccinations is characterized by  | exclude | single cell data |
| GSE210231 | GEO | RNA | CD8+ T cell memory induced by successive SARS-CoV-2 mRNA vaccinations is characterized by  | exclude | vaccination study |
| GSE210432 | GEO | RNA | TET2 Mutant Clonal Hematopoiesis Enhances Severity of COVID-19 Related Cytokine Release Sy | exclude | single cell data |
| GSE210433 | GEO | RNA | TET2 Mutant Clonal Hematopoiesis Enhances Severity of COVID-19 Related Cytokine Release Sy | exclude | single cell data |
| GSE210435 | GEO | RNA | TET2 Mutant Clonal Hematopoiesis Enhances Severity of COVID-19 Related Cytokine Release Sy | exclude | no characteristic splits 10 COVID samples from 10 others and no severity field (automatic; best : 0 against 0) |
| GSE210676 | GEO | RNA | TFEB and TFE3 are required for efficient cellular response against beta-coronavirus [RNA-s | exclude | 12 samples: cannot reach 10 participants per class |
| GSE211240 | GEO | RNA | MDA5 guards against infection by surveying cellular RNA homeostasis [I] | exclude | in vitro or cell culture |
| GSE211378 | GEO | RNA | Whole Blood profiling of COVID convalescent and Healthy donors with nCounter | exclude | convalescent, not acute |
| GSE211567 | GEO | RNA | Host-response Transcriptional Biomarkers Accurately Discriminate Bacterial and Viral Infec | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE212041 | GEO | RNA | Longitudinal characterization of circulating neutrophils uncovers phenotypes associated wi | exclude | in vitro or cell culture |
| GSE212584 | GEO | RNA | T cell cellular stress and reticulocyte signatures, but not loss of naïve T lymphocytes, c | exclude | single cell data |
| GSE212861 | GEO | RNA | Dynamics of gene expression profiling by RNA-Seq and identification of high-risk patients  | exclude | no characteristic splits 10 COVID samples from 10 others and no severity field (automatic; best : 0 against 0) |
| GSE213004 | GEO | RNA | Distinct responses of newly identified monocyte subsets to advanced gastrointestinal cance | exclude | single cell data |
| GSE213192 | GEO | RNA | SARS-CoV-2 infection and recovery in children: Distinct T cell responses in MIS-C compared | exclude | not adults |
| GSE21329 | GEO | RNA | Multi-tissue, selective PPARγ modulation of insulin sensitivity and metabolic pathways in  | exclude | not human (Rattus norvegicus) |
| GSE21351 | GEO | RNA | Ubiquinol-induced gene expression signatures are translated into reduced erythropoiesis an | exclude | not human |
| GSE213516 | GEO | RNA | Gene expression profile at single cell level of peripheral blood mononuclear cells (PBMCs) | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE214150 | GEO | RNA | The transcriptome of platelets from healthy donors treated with SARS-CoV-2 envelope (E) pr | exclude | 6 samples: cannot reach 10 participants per class |
| GSE214247 | GEO | RNA | SARS-CoV-2-associated ssRNAs activate human neutrophils in a TLR8-dependent fashion  | exclude | no characteristic splits 10 COVID samples from 10 others and no severity field (automatic; best : 0 against 0) |
| GSE214624 | GEO | RNA | PD-1highCXCR5–CD4+ Peripheral Helper T (Tph) cells Promote CXCR3+ Plasmablasts in the earl | exclude | in vitro or cell culture |
| GSE214663 | GEO | RNA | MDA5 guards against infection by surveying cellular RNA homeostasis [II] | exclude | in vitro or cell culture |
| GSE214664 | GEO | RNA | MDA5 guards against infection by surveying cellular RNA homeostasis | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE215262 | GEO | RNA | Transcriptional patterns of peripheral blood mononuclear cells (PBMCs) in SARS-CoV-2 infec | exclude | single cell data |
| GSE215433 | GEO | RNA | Impact of SARS-CoV-2 ORF6 and its variant polymorphisms on host responses and viral pathog | exclude | not human |
| GSE215802 | GEO | RNA | Adaptive immune responses to SARS-CoV-2 persist in the pharyngeal lymphoid tissue of child | exclude | single cell data |
| GSE215814 | GEO | RNA | IFN-λ treatment accelerates SARS-CoV-2 clearance despite age-related delays in the inducti | exclude | single cell data |
| GSE215865 | GEO | RNA | Molecular states during acute COVID-19 reveal distinct etiologies of long-term sequelae | exclude | already in CHORUS (acute samples on disk); the bars use new cohorts only |
| GSE215915 | GEO | RNA | Multiomic chromatin accessibility and gene expression in peripheral blood predicts future  | exclude | single cell data |
| GSE216020 | GEO | RNA | Heterogeneity of neutrophils and inflammatory responses in patients with COVID-19 and heal | exclude | single cell data |
| GSE216049 | GEO | RNA | Organoid modeling of lung-intrinsic immune responses to SARS-CoV-2 infection | exclude | in vitro or cell culture |
| GSE216486 | GEO | RNA | Human IRF1 governs phagocytic IFN-γ immunity to mycobacteria but not cell-intrinsic IFN-a/ | exclude | single cell data |
| GSE216488 | GEO | RNA | Human IRF1 governs phagocytic IFN-γ immunity to mycobacteria but not cell-intrinsic IFN-a/ | exclude | in vitro or cell culture |
| GSE216529 | GEO | RNA | Transcriptomic profiling of peripheral immune cells reveals distinct innate immune signatu | exclude | no characteristic splits 10 COVID samples from 10 others and no severity field (automatic; best : 0 against 0) |
| GSE216704 | GEO | RNA | Loss of GM-CSF-dependent instruction of alveolar macrophages in COVID-19 provides a ration | exclude | 19 samples: cannot reach 10 participants per class |
| GSE216811 | GEO | RNA | Augmenting vaccine efficacy against delta variant with ‘Mycobacterium-w’ mediated modulati | exclude | vaccination study |
| GSE216851 | GEO | RNA | SARS-CoV-2-induced disruption of a vascular bed in a microphysiological system caused by t | exclude | in vitro or cell culture |
| GSE217099 | GEO | RNA | Evaluation of broad anti-coronavirus activity of autophagy-related compounds using human a | exclude | in vitro or cell culture |
| GSE217181 | GEO | RNA | Selective translational control of cellular and viral mRNAs by RPS3 mRNA binding | exclude | 16 samples: cannot reach 10 participants per class |
| GSE217370 | GEO | RNA | Neutrophil profiles of pediatric COVID-19 and multisystem inflammatory syndrome in childre | exclude | not adults |
| GSE217918 | GEO | RNA | A TLR7-nanoparticle adjuvant promotes broad immune  responses against heterologous strains | exclude | not human |
| GSE21802 | GEO | RNA | Hosts responses in critical disease caused by pandemic H1N1 | exclude | other assay or analyte |
| GSE218033 | GEO | RNA | Human IRF1 governs phagocytic IFN-gamma immunity to mycobacteria | exclude | in vitro or cell culture |
| GSE218065 | GEO | RNA | Single cell RNA-seq for human and bovine PBMCs stimulated with M. tuberculosis | exclude | single cell data |
| GSE218576 | GEO | RNA | Single-cell multi-omics decodes regulatory programs during development of mouse secondary  | exclude | not human (Mus musculus) |
| GSE219086 | GEO | RNA | Increased Antibody Titers but Induced T-Cell AICD and Apoptosis Response in COVID-19 Conva | exclude | vaccination study |
| GSE219098 | GEO | RNA | Tracking B Cell Responses to the SARS-CoV-2 mRNA-1273 Vaccine | exclude | single cell data |
| GSE21927 | GEO | RNA | Expression data from bone marrow derived- and tumor induced- CD11b+ MDSC | exclude | not human (Mus musculus) |
| GSE21942 | GEO | RNA | Expression data from peripheral blood mononuclear cells in multiple sclerosis patients and | exclude | no characteristic splits 10 COVID samples from 10 others and no severity field (automatic; best : 0 against 0) |
| GSE220046 | GEO | RNA | Epigenetic regulation of SARS-Cov-2 spike protein modulates cellular viability and nascent | exclude | 9 samples: cannot reach 10 participants per class |
| GSE220076 | GEO | RNA | Associations of miRNAs with blood phenotypes and ICU admission in COVID-19 patients [mRNA- | exclude | other assay or analyte |
| GSE220078 | GEO | RNA | Associations of miRNAs with blood phenotypes and ICU admission in COVID-19 patients | exclude | other assay or analyte |
| GSE220190 | GEO | RNA | Mapping disease-associated regulatory circuits by cell type from single-cell multiomics da | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE22047 | GEO | RNA | Modulation of Cystatin A Expression in Human Airway Epithelium Related to Genotype, Smokin | exclude | single cell data |
| GSE220582 | GEO | RNA | Transcriptomic analysis of double-negative B cell subsets in IgG4-related disease | exclude | 15 samples: cannot reach 10 participants per class |
| GSE220659 | GEO | RNA | Activation of Coagulation and Proinflammatory Pathways in Thrombosis with Thrombocytopenia | exclude | vaccination study |
| GSE220682 | GEO | RNA | Time-series transcriptome analysis of peripheral blood mononuclear cells obtained from ind | exclude | vaccination study |
| GSE22098 | GEO | RNA | Whole blood transcriptional profiles of patients with active tuberculosis (TB) and other i | exclude | vaccination study |
| GSE221234 | GEO | RNA | Predicting severity in COVID-19 disease using sepsis blood gene expression signatures | exclude | COVID-19 only; severity at the BQC19 biobank (application needed) |
| GSE221316 | GEO | RNA | Class switch towards non-inflammatory, spike-specific IgG4 antibodies after repeated SARS- | exclude | single cell data |
| GSE221319 | GEO | RNA | Class switch towards non-inflammatory, spike-specific IgG4 antibodies after repeated SARS- | exclude | single cell data |
| GSE221320 | GEO | RNA | Class switch towards non-inflammatory, spike-specific IgG4 antibodies after repeated SARS- | exclude | vaccination study |
| GSE222128 | GEO | RNA | Influenza vaccination reveals sex dimorphic imprints of prior mild COVID-19 | exclude | vaccination study |
| GSE222431 | GEO | RNA | Demultiplexing Ig repertoires by parallel mRNA/DNA sequencing shows major differential alt | exclude | no characteristic splits 10 COVID samples from 10 others and no severity field (automatic; best : 0 against 0) |
| GSE22255 | GEO | RNA | Blood genomic expression profile for ischemic stroke (IS) | exclude | no characteristic splits 10 COVID samples from 10 others and no severity field (automatic; best : 0 against 0) |
| GSE222872 | GEO | RNA | Primary ChAdOx1 vaccination does not reactivate pre-existing, cross-reactive immunity | exclude | vaccination study |
| GSE223139 | GEO | RNA | Effect of culture/microexpansion in various conditions on virus-specific CD4+ T cells | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE223140 | GEO | RNA | Effect of culture/microexpansion in various conditions on virus-specific CD8+ T cells | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE223142 | GEO | RNA | Transcriptomic analysis reveals optimal cytokine combinations for SARS-CoV-2-specific T ce | exclude | no characteristic splits 10 COVID samples from 10 others and no severity field (automatic; best : 0 against 0) |
| GSE223236 | GEO | RNA | Cytotoxic T cells targeting spike glycoprotein are associated with hybrid immunity to SARS | exclude | single cell data |
| GSE223646 | GEO | RNA | Compartmentalization of human B cell responses in early life through induction of bronchus | exclude | not adults |
| GSE223793 | GEO | RNA | Transcriptomic evidence of lung repair in paediatric ARDS survival | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE224028 | GEO | RNA | Intradermal skin test with mRNA vaccines as a surrogate marker of T cell immunity in immun | exclude | single cell data |
| GSE224063 | GEO | RNA | COVD-19 infection during pregnancy induces a hematopoietic stem cell signature in cord blo | exclude | 17 samples: cannot reach 10 participants per class |
| GSE224131 | GEO | RNA | Transcriptional effects of siRNA-mediated MAFB  knockdown on SARS-CoV-2-infected human M-C | exclude | RNA from plasma or serum (cell-free RNA), not whole blood or PBMC |
| GSE224198 | GEO | RNA | Single-cell multi-omics analysis revealing immune features of inactivated SARS-CoV-2 vacci | exclude | single cell data |
| GSE224615 | GEO | RNA | Long COVID manifests with T cell dysregulation, inflammation, and an uncoordinated adaptiv | exclude | not adults |
| GSE224845 | GEO | RNA | Transcriptional effects of SARS-CoV-2 and their VLPs (virus-like particles) on human GM-CS | exclude | no characteristic splits 10 COVID samples from 10 others and no severity field (automatic; best : 0 against 0) |
| GSE22491 | GEO | RNA | Transcriptomic profiling of peripheral blood mononuclear cells from patients with Parkinso | exclude | 18 samples: cannot reach 10 participants per class |
| GSE225165 | GEO | RNA | Assessing chromatin accessibility in human classical monocytes following vaccination with  | exclude | vaccination study |
| GSE225217 | GEO | RNA | Effects of the circulating environment of COVID-19 on platelet and neutrophil behavior | exclude | 10 samples: cannot reach 10 participants per class |
| GSE225220 | GEO | RNA | Nucleic acid biomarkers of immune response and cell and tissue damage in children with COV | exclude | not adults |
| GSE225221 | GEO | RNA | Nucleic acid biomarkers of immune response and cell and tissue damage in children with COV | exclude | not adults |
| GSE225223 | GEO | RNA | Nucleic acid biomarkers of immune response and cell and tissue damage in children with COV | exclude | not adults |
| GSE226131 | GEO | RNA | Single-cell RNA sequencing of hematopoietic stem and progenitor cells from young and aged  | exclude | not human (Mus musculus) |
| GSE227341 | GEO | RNA | Gene expression profile of peripheral blood mononuclear cells from mild and severe COVID-1 | exclude | 16 samples: cannot reach 10 participants per class |
| GSE227562 | GEO | RNA | SARS-CoV-2 Omicron boosting induces de novo B cell response in humans | exclude | no characteristic splits 10 COVID samples from 10 others and no severity field (automatic; best : 0 against 0) |
| GSE227585 | GEO | RNA | Effects of Dexamethasone in COVID-19 patients | exclude | in vitro or cell culture |
| GSE227846 | GEO | RNA | The protective role of the HSP90 inhibitor, STA-9090, in lungs of SARS-CoV-2-infected Syri | exclude | not human |
| GSE227942 | GEO | RNA | A SINGLE-CELL RNA-SEQ ANALYSIS UNRAVELS THE HETEROGENEITY OF PRIMARY CULTURED HUMAN CORNEA | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE227971 | GEO | RNA | SARS-CoV-2 infection establishes an enhanced, stable, and age-independent CD8+ T cell resp | exclude | in vitro or cell culture |
| GSE228236 | GEO | RNA | Single-cell analysis of bronchoalveolar cells in inflammatory and fibrotic post-COVID lung | exclude | single cell data |
| GSE228779 | GEO | RNA | Morphologic and molecular analysis of liver injury after SARS-CoV-2 vaccination reveals di | exclude | vaccination study |
| GSE228839 | GEO | RNA | COVID-19 vaccination is associated with attenuation of inflammatory responses during break | exclude | vaccination study |
| GSE228840 | GEO | RNA | COVID-19 vaccination is associated with attenuation of inflammatory responses during break | exclude | vaccination study |
| GSE228842 | GEO | RNA | COVID-19 vaccination is associated with attenuation of inflammatory responses during break | exclude | vaccination study |
| GSE229169 | GEO | RNA | Comparative single cell epigenomic analysis of gene regulatory programs in the rodent and  | exclude | single cell data |
| GSE229187 | GEO | RNA | CITE-seq of peripheral immune cells after SARS-COV-2 vaccination in patients with MM | exclude | single cell data |
| GSE229459 | GEO | RNA | Morphologic and molecular analysis of liver injury after SARS-CoV-2 vaccination reveals di | exclude | vaccination study |
| GSE230022 | GEO | RNA | Transcriptome of lung endothelial cells in SARS-CoV-2 infection model | exclude | not human (Mus musculus) |
| GSE230227 | GEO | RNA | Cytokinopathy with aberrant cytotoxic lymphocytes and profibrotic myeloid response in SARS | exclude | single cell data |
| GSE230785 | GEO | RNA | Transcriptomic Profiling of SARS-CoV-2 Spike Protein S1 Subunit Mediated Changes in Human  | exclude | vaccination study |
| GSE23117 | GEO | RNA | Gene expression in minor salivary gland of patients with Sjogren's syndrome (SS) and contr | exclude | not blood |
| GSE231409 | GEO | RNA | Upper respiratory and peripheral blood gene expression profiles of children, adolescents,  | exclude | not adults |
| GSE231733 | GEO | RNA | Bulk RNA-sequencing of A549 cells transduced with Ad5 empty vector (GFP) or Ad5 expressing | exclude | in vitro or cell culture |
| GSE232226 | GEO | RNA | SARS-CoV-2 viral RNA disseminates to hamster toes and associates with localized IFN-I prod | exclude | not human |
| GSE232523 | GEO | RNA | Unique molecular signatures sustained in circulating monocytes and regulatory T cells inco | exclude | single cell data |
| GSE232616 | GEO | RNA | Expansion of profibrotic monocyte-derived alveolar macrophages in patients with persistent | exclude | single cell data |
| GSE232623 | GEO | RNA | Expansion of profibrotic monocyte-derived alveolar macrophages in patients with persistent | exclude | single cell data |
| GSE232627 | GEO | RNA | Expansion of profibrotic monocyte-derived alveolar macrophages in patients with persistent | exclude | single cell data |
| GSE232705 | GEO | RNA | Distinct and overlapping immunological responses to SARS-CoV-2 and Mycobacterium tuberculo | exclude | single cell data |
| GSE233046 | GEO | RNA | Inflammation durably imprints memory CD4+ T cells | exclude | in vitro or cell culture |
| GSE233230 | GEO | RNA | Memory B cell development elicited by mRNA booster vaccinations in the elderly | exclude | vaccination study |
| GSE233315 | GEO | RNA | Immune landscape of PBMC based on single-cell rna sequancing in End Stage Renal Disease | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE233506 | GEO | RNA | Microfluidic immuno-serolomic assay reveals systems level association with COVID-19 pathol | exclude | single cell data |
| GSE233522 | GEO | RNA | Clonal landscape of autoantibody-secreting plasmablasts in COVID-19 patients | exclude | 9 samples: cannot reach 10 participants per class |
| GSE234008 | GEO | RNA | Transcriptomes analysis of COVID-19 Inactivated Vaccine and Protein Subunit Vaccine as Thi | exclude | vaccination study |
| GSE234118 | GEO | RNA | Mpox infection protects against re-challenge in rhesus macaques | exclude | not human (Macaca mulatta) |
| GSE234207 | GEO | RNA | The difference between neutrophils of HC, COVID-19 patients admitted in ICU and wards. | exclude | no characteristic splits 10 COVID samples from 10 others and no severity field (automatic; best : 0 against 0) |
| GSE234902 | GEO | RNA | Association of Neutrophil Gene Expression with 28-Day Ventilator Weaning in COVID-19 ARDS  | exclude | single cell data |
| GSE234903 | GEO | RNA | Association of Neutrophil Gene Expression with 28-Day Ventilator Weaning in COVID-19 ARDS  | exclude | single cell data |
| GSE234904 | GEO | RNA | Association of Neutrophil Gene Expression with 28-Day Ventilator Weaning in COVID-19 ARDS  | exclude | 7 samples: cannot reach 10 participants per class |
| GSE234918 | GEO | RNA | The Transcriptional and Phenotypic Characteristics that Define Alveolar Macrophage Subsets | exclude | single cell data |
| GSE235050 | GEO | RNA | scRNA-seq characterization of immune cell types among individuals with or without Long COV | exclude | single cell data |
| GSE235094 | GEO | RNA | Human monocytes exposed to SARS-CoV-2 display features of trained immunity producing high  | exclude | vaccination study |
| GSE235136 | GEO | RNA | Trimannose-coupled antimiR-21 for macrophage-targeted inhalation treatment of acute inflam | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE235164 | GEO | RNA | Mouse genome rewriting and tailoring of three important disease loci | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE235232 | GEO | RNA | Transcriptome analysis of transplanted vascular organoid grafts under cranial window at th | exclude | not human |
| GSE235433 | GEO | RNA | SARS-CoV-2 infection triggers pro-atherogenic inflammatory responses in human coronary [MA | exclude | not blood |
| GSE235437 | GEO | RNA | SARS-CoV-2 infection triggers pro-atherogenic inflammatory responses in human coronary | exclude | no characteristic splits 10 COVID samples from 10 others and no severity field (automatic; best title: 30 against 8) |
| GSE235938 | GEO | RNA | Single‐cell RNA sequencing reveals characteristics of myeloid cells in pulmonary post-acut | exclude | single cell data |
| GSE236030 | GEO | RNA | Distinct pulmonary and systemic effects of dexamethasone in severe COVID-19 | exclude | single cell data |
| GSE236177 | GEO | RNA | Key patient demographics shape innate immune topography in non-critical hypoxic COVID-19 p | exclude | single cell data |
| GSE236294 | GEO | RNA | Immunometabolic effect of nitric oxide on human macrophages challenged with the SARS-CoV2- | exclude | not blood |
| GSE236651 | GEO | RNA | Stimulation of PSTPIP1 to trigger proinflammatory responses in asymptomatic SARS-CoV-2 inf | exclude | no characteristic splits 10 COVID samples from 10 others and no severity field (automatic; best : 0 against 0) |
| GSE23687 | GEO | RNA | Expression data from SPARKS CHARMS JIA cohort | exclude | not adults |
| GSE237092 | GEO | RNA | Impact of Age and Sex on neuroinflammation due to SARS-CoV2 infection in a murine model. | exclude | not human (Mus musculus) |
| GSE237109 | GEO | RNA | Distinct pulmonary and systemic effects of dexamethasone in severe COVID-19 | exclude | single cell data |
| GSE237132 | GEO | RNA | Identification of AHR as a pro-viral host factor of human betacoronaviruses SARS-CoV-2 and | exclude | no characteristic splits 10 COVID samples from 10 others and no severity field (automatic; best title: 1 against 7) |
| GSE237576 | GEO | RNA | Effects of influenza A and SARS-CoV-2 coinfection on interferon stimulated gene expression | exclude | in vitro or cell culture |
| GSE237965 | GEO | RNA | A Phase I, First-in-Human Study to Evaluate the Safety and Tolerability, Pharmacokinetics  | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE238079 | GEO | RNA | Characterizing immune signatures of seasonal coronaviruses HCoV-229E and -OC43 in human na | exclude | not blood |
| GSE239585 | GEO | RNA | Pharmacological inhibition of bromodomain and extra-terminal proteins induces NRF-2-mediat | exclude | not blood |
| GSE239595 | GEO | RNA | Bulk RNA sequencing in SARS-CoV-2 infected patients | exclude | not blood |
| GSE239787 | GEO | RNA | Systems biological assessment of the temporal dynamics of immunity to a viral infection in | exclude | single cell data |
| GSE239818 | GEO | RNA | Unveiling novel conserved HIV-1 open reading frames encoding T-cell antigens using ribosom | exclude | in vitro or cell culture |
| GSE241039 | GEO | RNA | Developing enhanced immunotherapy using NKG2A knockout human pluripotent stem cell-derived | exclude | not human |
| GSE2411 | GEO | RNA | Mechanical Ventilation Augments LPS-induced Lung Injury | exclude | not human (Mus musculus) |
| GSE241135 | GEO | RNA | Human RP105 monoclonal antibody enables the production of antigen-specific antibody | exclude | single cell data |
| GSE241385 | GEO | RNA | CSF single cell transcriptomic changes associated with post-acute sequelae of SARS CoV 2 i | exclude | single cell data |
| GSE241717 | GEO | RNA | NULISA: a novel proteomic liquid biopsy platform with attomolar sensitivity and high multi | exclude | 3 samples: cannot reach 10 participants per class |
| GSE241843 | GEO | RNA | SARS-CoV-2 spike protein induces the cytokine release syndrome by stimulating T cells to p | exclude | 12 samples: cannot reach 10 participants per class |
| GSE241849 | GEO | RNA | Absolute quantitative and base-resolution sequencing reveals comprehensive landscape of ps | exclude | single cell data |
| GSE242127 | GEO | RNA | A Hub Gene Signature as a Therapeutic Target and Biomarker for Sepsis and Geriatric Sepsis | exclude | single cell data |
| GSE24250 | GEO | RNA | The transcriptional modulator H2AFY marks Huntington's disease activity in men and mice | exclude | not human |
| GSE242722 | GEO | RNA | Differential gene expressions of various interferon, antiviral response genes as well as t | exclude | in vitro or cell culture |
| GSE242775 | GEO | RNA | Determinants of de novo B cell responses to drifted epitopes in post-vaccination SARS-CoV- | exclude | vaccination study |
| GSE243217 | GEO | RNA | Comparative Analysis of Host Responses Between Sepsis and COVID-19: A Prospective Observat | exclude | other assay or analyte |
| GSE243268 | GEO | RNA | Multi-target mode of action of Sulfodyne, a stabilized Sulforaphane, against pathogenic ef | exclude | not human |
| GSE243348 | GEO | RNA | Longitudinal gene expression profiling of self-collected blood samples in COVID-19+ and he | exclude | vaccination study |
| GSE243442 | GEO | RNA | Integrative biology of patients acutely infected with yellow fever virus | exclude | vaccination study |
| GSE244026 | GEO | RNA | Switching of OAS1 splicing isoforms mitigates SARS-CoV-2 infection | exclude | 6 samples: cannot reach 10 participants per class |
| GSE244173 | GEO | RNA | Single-cell sequencing of 12 PBMC samples to study the different immune response after SAR | exclude | single cell data |
| GSE244297 | GEO | RNA | Memory B cells dominate the early antibody-secreting cell response to  SARS-CoV-2 mRNA vac | exclude | single cell data |
| GSE244484 | GEO | RNA | Transcriptomic analysis of sorted lung cells revealed a proviral activity of the NF-κB pat | exclude | single cell data |
| GSE244485 | GEO | RNA | Transcriptomic analysis of sorted lung cells revealed a proviral activity of the NF-κB pat | exclude | single cell data |
| GSE244486 | GEO | RNA | Transcriptomic analysis of sorted lung cells revealed a proviral activity of the NF-κB pat | exclude | single cell data |
| GSE244488 | GEO | RNA | Transcriptomic analysis of sorted lung cells revealed a proviral activity of the NF-κB pat | exclude | no characteristic splits 10 COVID samples from 10 others and no severity field (automatic; best : 0 against 0) |
| GSE244714 | GEO | RNA | Nanobodies agains Sars Cov 2 Nsp9 protein protects cells from virus infection | exclude | in vitro or cell culture |
| GSE244745 | GEO | RNA | B cell-reactive neoantigens boost immune response in personalized cancer vaccination | exclude | single cell data |
| GSE245466 | GEO | RNA | Transcriptional effects of SARS-CoV-2 E protein and altered BET function | exclude | 8 samples: cannot reach 10 participants per class |
| GSE245616 | GEO | RNA | Immunomodulatory Potential of Dapphitkhai and Chanthaleela Thai Traditional Recipes in the | exclude | 9 samples: cannot reach 10 participants per class |
| GSE245922 | GEO | RNA | ACE2-independent SARS-CoV-2 virus entry through cell surface GRP78 on monocytes - evidence | exclude | 8 samples: cannot reach 10 participants per class |
| GSE246407 | GEO | RNA | Discovery of a novel inhibitor of macropinocytosis with antiviral activity | exclude | 6 samples: cannot reach 10 participants per class |
| GSE246452 | GEO | RNA | Comparing the transcriptomic profiles of anti-viral candidates knockout cells | exclude | single cell data |
| GSE246491 | GEO | RNA | Resistance Potential of the HLA-A2-restricted Immunodominant SARS-CoV-2 Specific CD8+ T Ce | exclude | vaccination study |
| GSE246525 | GEO | RNA | Systems biology analysis reveals distinct molecular signatures promoting immune responsive | exclude | single cell data |
| GSE246535 | GEO | RNA | T cell responses in SARS-CoV-2 vaccinees with sustained or declined antibody titer | exclude | single cell data |
| GSE246536 | GEO | RNA | Single-nucleus transcriptional profiling of the placenta reveals the syncytiotrophoblast s | exclude | single cell data |
| GSE246937 | GEO | RNA | Systems biology analysis reveals distinct molecular signatures promoting immune responsive | exclude | single cell data |
| GSE246938 | GEO | RNA | Systems biology analysis reveals distinct molecular signatures promoting immune responsive | exclude | vaccination study |
| GSE247186 | GEO | RNA | Monocyte-mediated thrombosis linked to circulating tissue factor and immune-paralysis in C | exclude | 12 samples: cannot reach 10 participants per class |
| GSE247401 | GEO | RNA | Blood transcriptome profiling reveals distinct gene networks induced by mRNA vaccination a | exclude | vaccination study |
| GSE247406 | GEO | RNA | Heat inactivation of Nipah virus for downstream single cell RNA sequencing | exclude | single cell data |
| GSE247488 | GEO | RNA | Non-classical CD45RBlo memory B-cells are the majority of circulating antigen-specific B-c | exclude | single cell data |
| GSE24752 | GEO | RNA | Microarray Analysis of differential gene expression in peripheral blood cells of patients  | exclude | 6 samples: cannot reach 10 participants per class |
| GSE247582 | GEO | RNA | Single-cell analysis of CD4+ cytotoxic T lymphocytes in human oral squamous cell carcinoma | exclude | single cell data |
| GSE247719 | GEO | RNA | A Panoramic View of Cell Population Dynamics in Mammalian Aging | exclude | not human (Mus musculus) |
| GSE247901 | GEO | RNA | mRNA COVID-19 vaccine elicits potent adaptive immune response without the acute inflammati | exclude | single cell data |
| GSE247902 | GEO | RNA | mRNA COVID-19 vaccine elicits potent adaptive immune response without the acute inflammati | exclude | single cell data |
| GSE247903 | GEO | RNA | mRNA COVID-19 vaccine elicits potent adaptive immune response without the acute inflammati | exclude | single cell data |
| GSE247904 | GEO | RNA | mRNA COVID-19 vaccine elicits potent adaptive immune response without the acute inflammati | exclude | single cell data |
| GSE247905 | GEO | RNA | mRNA COVID-19 vaccine elicits potent adaptive immune response without the acute inflammati | exclude | single cell data |
| GSE247906 | GEO | RNA | mRNA COVID-19 vaccine elicits potent adaptive immune response without the acute inflammati | exclude | single cell data |
| GSE247907 | GEO | RNA | mRNA COVID-19 vaccine elicits potent adaptive immune response without the acute inflammati | exclude | single cell data |
| GSE247908 | GEO | RNA | mRNA COVID-19 vaccine elicits potent adaptive immune response without the acute inflammati | exclude | single cell data |
| GSE247909 | GEO | RNA | mRNA COVID-19 vaccine elicits potent adaptive immune response without the acute inflammati | exclude | single cell data |
| GSE247910 | GEO | RNA | mRNA COVID-19 vaccine elicits potent adaptive immune response without the acute inflammati | exclude | single cell data |
| GSE247911 | GEO | RNA | mRNA COVID-19 vaccine elicits potent adaptive immune response without the acute inflammati | exclude | single cell data |
| GSE247912 | GEO | RNA | mRNA COVID-19 vaccine elicits potent adaptive immune response without the acute inflammati | exclude | single cell data |
| GSE247913 | GEO | RNA | mRNA COVID-19 vaccine elicits potent adaptive immune response without the acute inflammati | exclude | single cell data |
| GSE247914 | GEO | RNA | mRNA COVID-19 vaccine elicits potent adaptive immune response without the acute inflammati | exclude | single cell data |
| GSE247915 | GEO | RNA | mRNA COVID-19 vaccine elicits potent adaptive immune response without the acute inflammati | exclude | single cell data |
| GSE247916 | GEO | RNA | mRNA COVID-19 vaccine elicits potent adaptive immune response without the acute inflammati | exclude | single cell data |
| GSE24814 | GEO | RNA | Role and function of STAT5 in BCR-ABL1 driven pre-B cells | exclude | not human (Mus musculus) |
| GSE248556 | GEO | RNA | Single-cell transcriptome of PBMC response to Omicron breakthrough infection | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE248669 | GEO | RNA | Single-cell characteristics of longitudinal immune profiling in subjects with breakthrough | exclude | single cell data |
| GSE249313 | GEO | RNA | CD4+ T cells exhibit distinct transcriptional phenotypes in the lymph nodes and blood foll | exclude | single cell data |
| GSE249474 | GEO | RNA | Identification of a human hematopoietic stem cell subset that retains memory of inflammato | exclude | single cell data |
| GSE249476 | GEO | RNA | Identification of a human hematopoietic stem cell subset that retains memory of inflammato | exclude | single cell data |
| GSE249478 | GEO | RNA | Identification of a human hematopoietic stem cell subset that retains memory of inflammato | exclude | single cell data |
| GSE249500 | GEO | RNA | Transcriptomic Analysis of THP-1 Cells  Electroporated with GFP-EV or GFP-SARS-CoV-2-ORF6 | exclude | 4 samples: cannot reach 10 participants per class |
| GSE249513 | GEO | RNA | High-dimensional comparison of monocytes and T cells in post-COVID and idiopathic pulmonar | exclude | single cell data |
| GSE250023 | GEO | RNA | Impaired innate and adaptive immune responses to BNT162b2 SARS-CoV-2 vaccination in system | exclude | single cell data |
| GSE250024 | GEO | RNA | Impaired innate and adaptive immune responses to BNT162b2 SARS-CoV-2 vaccination in system | exclude | single cell data |
| GSE250025 | GEO | RNA | Impaired innate and adaptive immune responses to BNT162b2 SARS-CoV-2 vaccination in system | exclude | vaccination study |
| GSE250377 | GEO | RNA | SARS-CoV-2 infection response in iPSC-derived alveolar epithelial cell type-2 | exclude | in vitro or cell culture |
| GSE251849 | GEO | RNA | Blood–brain barrier disruption and sustained systemic inflammation in individuals with lon | exclude | in vitro or cell culture |
| GSE252056 | GEO | RNA | SARS-CoV-2 infection evokes transcriptional changes in Calu-3 cells and cultured primary m | exclude | not human |
| GSE252508 | GEO | RNA | Frequency, kinetics and determinants of viable SARS-CoV-2 in bioaerosols from ambulatory C | exclude | in vitro or cell culture |
| GSE252959 | GEO | RNA | Recruitment of plasma cells from follicular and extrafollicular immune reactions to the bo | exclude | single cell data |
| GSE253234 | GEO | RNA | Lung epithelium releases growth differentiation factor 15 in response to pathogen-mediated | exclude | not human (Mus musculus) |
| GSE253845 | GEO | RNA | Neural Network-Assisted Humanization of COVID-19 Hamster scRNAseq Data Reveals Matching Se | exclude | not human (Phodopus roborovskii) |
| GSE253857 | GEO | RNA | Recruitment of plasma cells from follicular and extrafollicular immune reactions to the bo | exclude | single cell data |
| GSE253862 | GEO | RNA | Recruitment of plasma cells from follicular and extrafollicular immune reactions to the bo | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE254179 | GEO | RNA | TGFβ links EBV to multisystem inflammatory syndrome in children | exclude | not adults |
| GSE254309 | GEO | RNA | SARS-CoV-2 causes dysfunction in human iPSC-derived brain microvascular endothelial like c | exclude | in vitro or cell culture |
| GSE254537 | GEO | RNA | Inhibition of CXCL10 and IFN-γ ameliorates myocarditis in preclinical models of SARS-CoV-  | exclude | not human |
| GSE254647 | GEO | RNA | RNA expression differences in innate immune response to COVID19 with advancing age. | exclude | no characteristic splits 10 COVID samples from 10 others and no severity field (automatic; best : 0 against 0) |
| GSE254699 | GEO | RNA | Transcriptional profiling of SARS-CoV-2-infected Calu-3 cells. | exclude | in vitro or cell culture |
| GSE254732 | GEO | RNA | Inhibition of PCSK9 with evolocumab modulates lipoproteins and monocyte activation in high | exclude | no characteristic splits 10 COVID samples from 10 others and no severity field (automatic; best : 0 against 0) |
| GSE25504 | GEO | RNA | Whole blood mRNA expression profiling of host molecular networks in neonatal sepsis | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE25507 | GEO | RNA | Autism and increased paternal age related changes in global levels of gene expression regu | exclude | not human |
| GSE255647 | GEO | RNA | Characterizing spatial and global host innate immune responses against SARS-CoV-1 and -2 i | exclude | in vitro or cell culture |
| GSE256068 | GEO | RNA | Identification of gene regulatory networks affected across drug-resistant epilepsies | exclude | not blood |
| GSE25640 | GEO | RNA | Expression data from wild type or FIZZ2 knockout murine lungs | exclude | not human (Mus musculus) |
| GSE2565 | GEO | RNA | Phosgene exposure in the mouse lung | exclude | not human (Mus musculus) |
| GSE25724 | GEO | RNA | Expression data from type 2 diabetic and non-diabetic isolated human islets | exclude | 13 samples: cannot reach 10 participants per class |
| GSE25837 | GEO | RNA | Differential gene expression in circulating leukocytes from chronically lonely individuals | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE259214 | GEO | RNA | The Transcriptomic Analysis of Profibrinolytic (PLAUR) and Fibrinolytic Inhibitor (SERPINE | exclude | RNA from plasma or serum (cell-free RNA), not whole blood or PBMC |
| GSE259276 | GEO | RNA | Prolonged exposure to lung-derived cytokines is associated with activation of microglia in | exclude | not human |
| GSE260763 | GEO | RNA | Expansion of memory Vδ2 T cells following SARS-CoV-2 vaccination revealed by temporal sing | exclude | single cell data |
| GSE26104 | GEO | RNA | Search for specific biomarkers of IFN-beta bioactivity in patients with MS | exclude | in vitro or cell culture |
| GSE261278 | GEO | RNA | Maintenance and functional regulation of immune memory to COVID-19 vaccines in tissues | exclude | vaccination study |
| GSE261419 | GEO | RNA | Transcriptomic analysis of Omicron Pseudovirus(PsV)-infected HSPCs(hematopoietic stem and  | exclude | 12 samples: cannot reach 10 participants per class |
| GSE261862 | GEO | RNA | Prior vaccination prevents overactivation of innate immune responses during COVID-19 break | exclude | single cell data |
| GSE261966 | GEO | RNA | Epitope-specific T cell responses after SARS-CoV-2 vaccination [scRNA-dextramer_s2] | exclude | single cell data |
| GSE261967 | GEO | RNA | Epitope-specific T cell responses after SARS-CoV-2 vaccination [scRNA_dextramer_mixed_time | exclude | single cell data |
| GSE262379 | GEO | RNA | Paired analysis of host and pathogen genomes identifies determinants of human tuberculosis | exclude | in vitro or cell culture |
| GSE262861 | GEO | RNA | Improvement of immune dysregulation and health-related quality of life in individuals with | exclude | single cell data |
| GSE263469 | GEO | RNA | Interferon-responsive neutrophils and macrophages extricate SARS-CoV-2 Omicron critical pa | exclude | RNA from plasma or serum (cell-free RNA), not whole blood or PBMC |
| GSE2638 | GEO | RNA | HMEC gene profile after TNF-stimulation | exclude | 6 samples: cannot reach 10 participants per class |
| GSE263817 | GEO | RNA | Comparative single-cell analysis reveals IFN-g as a driver of respiratory sequelae post CO | exclude | single cell data |
| GSE2639 | GEO | RNA | HUVEC gene profile after TNF-stimulation | exclude | 8 samples: cannot reach 10 participants per class |
| GSE2640 | GEO | RNA | strain dependent pulmonary bleomycin response | exclude | not human (Mus musculus) |
| GSE264196 | GEO | RNA | A 50-gene high-risk profile predictive of COVID-19 and Idiopathic Pulmonary Fibrosis morta | exclude | single cell data |
| GSE264443 | GEO | RNA | Dynamic gene expression analysis reveals distinct severity phases of immune and cellular d | exclude | COVID-19 only; phase at sampling, not a per participant severity |
| GSE265753 | GEO | RNA | Single cell gene expression profile of peripheral whole blood cells from long COVID patien | exclude | single cell data |
| GSE265975 | GEO | RNA | Transcriptomic analysis after SARS-CoV-2 mRNA vaccination reveals a specific gene signatur | exclude | vaccination study |
| GSE266519 | GEO | RNA | Hematopoietic Stem Cell Transplant and Healthy Donor Single cell RNA-seq (scRNA-seq) combi | exclude | single cell data |
| GSE267418 | GEO | RNA | Cytokine-induced memory-like responses in endothelial cells link chronic inflammation to v | exclude | no characteristic splits 10 COVID samples from 10 others and no severity field (automatic; best : 0 against 0) |
| GSE267625 | GEO | RNA | Whole blood transcriptome in long-COVID patients reveals association with lung function an | exclude | post-acute only (no acute samples by its record) |
| GSE26817 | GEO | RNA | Impact of IRF1 on transcriptome changes in Huh-7 human hepatoma cells and STAT1-deficient  | exclude | in vitro or cell culture |
| GSE268230 | GEO | RNA | Extracellular vesicles alter trophoblast function in pregnancies complicated by COVID-19 | exclude | in vitro or cell culture |
| GSE268369 | GEO | RNA | Effects of envelope- or membrane-protein segments of SARS-CoV-2 on gene expression of HUVE | exclude | in vitro or cell culture |
| GSE268516 | GEO | RNA | MicroRNA expression profling from 14 adult whole blood samples with flow cytometry cell fr | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE269423 | GEO | RNA | Longitudinal Multiomic Profiling and Corticosteroid Modulation of the Immediate Innate Imm | exclude | single cell data |
| GSE269780 | GEO | RNA | Differential blood transcriptional signatures induced by SARS-CoV-2 Variants Alpha and Del | exclude | no characteristic splits 10 COVID samples from 10 others and no severity field (automatic; best : 0 against 0) |
| GSE270479 | GEO | RNA | Sex differences and immune correlates of Long COVID development, persistence, and resoluti | exclude | 15 samples: cannot reach 10 participants per class |
| GSE270518 | GEO | RNA | Immunological memory diversity in the human upper airway | exclude | single cell data |
| GSE270609 | GEO | RNA | Maturing Neutrophils of Lower Density Associate with Thrombocytopenia in Puumala Orthohant | exclude | 16 samples: cannot reach 10 participants per class |
| GSE271284 | GEO | RNA | Rewired type I IFN signaling is linked to age-dependent differences in COVID-19 | exclude | in vitro or cell culture |
| GSE271378 | GEO | RNA | Exploring host cellular responses to SARS-CoV-2 variants using RNA sequencing of infected  | exclude | not blood |
| GSE271988 | GEO | RNA | IFN-ɣ derived from activated human CD4+ T cells inhibits the replication of SARS-CoV-2 dep | exclude | in vitro or cell culture |
| GSE272359 | GEO | RNA | Long-term efficacy of carbon dioxide laser in patients with Hailey-Hailey disease. Assessm | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE272381 | GEO | RNA | Transcriptomic Profiling of Neutrophils and Low-Density Granulocytes in COVID-19 Patients | exclude | no characteristic splits 10 COVID samples from 10 others and no severity field (automatic; best title: 29 against 7) |
| GSE272392 | GEO | RNA | Transcriptional response to the SARS-CoV-2 Delta variant in patients with severe COVID-19 | exclude | 12 samples: cannot reach 10 participants per class |
| GSE272645 | GEO | RNA | Features of chronic urticaria after COVID-19 mRNA vaccine, a real-life cohort study | exclude | vaccination study |
| GSE2729 | GEO | RNA | Rotavirus activates B but impairs T lymphocytes | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE273120 | GEO | RNA | Dasatinib Overcomes Acquired Resistance to BRAF/MEK Inhibitor Therapy | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE273137 | GEO | RNA | Identification of immune signatures associated with both SARS-CoV-2 infection and lung tra | exclude | single cell data |
| GSE273149 | GEO | RNA | Longitudinal transcriptomic analysis reveals persistent enrichment of iron homeostasis and | exclude | no characteristic splits 10 COVID samples from 10 others and no severity field (automatic; best : 0 against 0) |
| GSE273994 | GEO | RNA | Antigen-specific T helper cells predict intensity and longevity of immune responses to SAR | exclude | in vitro or cell culture |
| GSE274026 | GEO | RNA | Host transcriptome in human 239T cells infected with HPIV3 expressing wild-type virus vers | exclude | no characteristic splits 10 COVID samples from 10 others and no severity field (automatic; best : 0 against 0) |
| GSE274310 | GEO | RNA | Host transcriptome in human lung A549 cells infected with SARS-CoV-2 expressing wild-type  | exclude | not human |
| GSE274826 | GEO | RNA | Hyperglycemia-triggered lipid peroxidation destabilizes STAT4 and impairs anti-viral Th1 r | exclude | not human |
| GSE276140 | GEO | RNA | Unveiling Alternate Pathways for SARS-CoV-2 Infection via Extracellular Vesicle-Mediated T | exclude | other assay or analyte |
| GSE276544 | GEO | RNA | MVA-based SARS-CoV-2 vaccine candidates encoding different spike protein conformations ind | exclude | vaccination study |
| GSE276765 | GEO | RNA | Disruption of Immune Responses By Type I Diabetes Exacerbates SARS-CoV-2 Mediated Lung Inj | exclude | not human (Mus musculus) |
| GSE276775 | GEO | RNA | Disruption of Immune Responses By Type I Diabetes Exacerbates SARS-CoV-2 Mediated Lung Inj | exclude | not human (Mus musculus) |
| GSE277576 | GEO | RNA | Effects of Ancestral SARS-CoV-2 versus Variants of Concern Delta and Omicron on Megakaryoc | exclude | 16 samples: cannot reach 10 participants per class |
| GSE277947 | GEO | RNA | Investigating cellular transcriptomic alterations caused by VEEV replicase and SARS-CoV-2  | exclude | vaccination study |
| GSE279753 | GEO | RNA | The respective roles of TMPRSS2 and cathepsins for SARS-CoV-2 infection in human respirato | exclude | not human |
| GSE281135 | GEO | RNA | The identification of a SARs-CoV2 S1 protein derived peptide with super-antigen-like stimu | exclude | not human |
| GSE281237 | GEO | RNA | Widespread gene-environment interactions shape the immune response to SARS-CoV-2 infection | exclude | single cell data |
| GSE281864 | GEO | RNA | Detailed delineation of the immune response to COVID-19 mRNA vaccination via high-temporal | exclude | vaccination study |
| GSE282172 | GEO | RNA | SLAMF7 Defines Subsets of Human Effector CD8 T Cells. | exclude | no characteristic splits 10 COVID samples from 10 others and no severity field (automatic; best : 0 against 0) |
| GSE28221 | GEO | RNA | Peripheral Blood Mononuclear Cell Gene Expression Profiles May Predict Poor Outcome in Idi | exclude | no characteristic splits 10 COVID samples from 10 others and no severity field (automatic; best : 0 against 0) |
| GSE282879 | GEO | RNA | Natural product inhibition of SARS-CoV-2 guanine-N7 RNA cap methyltransferase shows potent | exclude | 8 samples: cannot reach 10 participants per class |
| GSE282892 | GEO | RNA | Unveiling Multi-Omics and Immune Landscapes in COVID-19 Patients' Ovaries Microenvironment | exclude | no characteristic splits 10 COVID samples from 10 others and no severity field (automatic; best title: 9 against 17) |
| GSE283746 | GEO | RNA | Comparative Single-Cell Analyses in Infants Revealed RSV-Specific Declines in NK Cell Resp | exclude | single cell data |
| GSE284430 | GEO | RNA | Regulation of CD45 isoforms during human effector and memory CD8 T cell differentiation: I | exclude | vaccination study |
| GSE284624 | GEO | RNA | Intranasal administration of a multiplex targeting SCLMs@FEN1-hpDNAs with synergistic baic | exclude | not human |
| GSE28489 | GEO | RNA | miRNA expression profiling of human immune cell subsets (HUG) | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE285173 | GEO | RNA | Single-cell RNA sequencing of PBMCs from COVID-19 recovered adolescents (n = 3) and adults | exclude | single cell data |
| GSE285756 | GEO | RNA | Early innate immune response and evolution of a SARS-CoV-2 furin cleavage site inactive va | exclude | no characteristic splits 10 COVID samples from 10 others and no severity field (automatic; best : 0 against 0) |
| GSE286325 | GEO | RNA | Transcriptomes of classical monocytes in post-COVID patients suggest predisposing conditio | exclude | single cell data |
| GSE286429 | GEO | RNA | Gene expression profiles of nasal epithelium in long COVID and post-COVID | exclude | 12 samples: cannot reach 10 participants per class |
| GSE287796 | GEO | RNA | Investigating the Risk of Arrhythmogenesis Associated with Fentanyl Abuse Using Mice and H | exclude | not human |
| GSE287808 | GEO | RNA | Temporal profiling of human lymphoid tissues reveals coordinated defence against viral cha | exclude | not blood |
| GSE287876 | GEO | RNA | Human naïve-like T follicular regulatory cells are primed for differentiation into mature  | exclude | in vitro or cell culture |
| GSE287877 | GEO | RNA | Human naïve-like T follicular regulatory cells are primed for differentiation into mature  | exclude | in vitro or cell culture |
| GSE288251 | GEO | RNA | Human long noncoding RNA, VILMIR, is induced by major respiratory viral infections and mod | exclude | single cell data |
| GSE288346 | GEO | RNA | Multi-Omics Analysis of the Immunological State at Baseline and Following SARS-CoV-2 mRNA  | exclude | single cell data |
| GSE288562 | GEO | RNA | Multi-Omics Analysis of the Immunological State at Baseline and Following SARS-CoV-2 mRNA  | exclude | single cell data |
| GSE288962 | GEO | RNA | Combination immunotherapy induces post-intervention control of HIV | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE289384 | GEO | RNA | Cell surface TMPRSS11E interacts with CD98 to mediate SARS-CoV2 efficient entry of M1 macr | exclude | 5 samples: cannot reach 10 participants per class |
| GSE289435 | GEO | RNA | Single-cell transcriptional atlas of hematopoiesis reveals genetic and hierarchy-based det | exclude | single cell data |
| GSE289648 | GEO | RNA | Maternal acute SARS-CoV-2 infection impairs preimplantation embryo development and reprogr | exclude | single cell data |
| GSE290006 | GEO | RNA | Reduced vaccine-induced germinal center outputs in inflammatory bowel disease patients tre | exclude | single cell data |
| GSE290079 | GEO | RNA | SARS-CoV-2 NSP14 Inhibitor Confers Potent Antiviral Efficacy and Reverses NSP14-triggered  | exclude | not human |
| GSE290188 | GEO | RNA | Minimal impact of ivermectin on immune response and transcriptional profiles in naïve adul | exclude | in vitro or cell culture |
| GSE291409 | GEO | RNA | BALF gene expression profiling reveals innate immune activation bias in persistent post-CO | exclude | not blood |
| GSE291745 | GEO | RNA | Macrophages direct location-dependent recall of B cell memory to vaccination [scRNA-seq an | exclude | single cell data |
| GSE292810 | GEO | RNA | Germinal center-mediated broadening of B cell responses to SARS-CoV-2 booster immunization | exclude | not human |
| GSE29366 | GEO | RNA | Whole Blood Transcriptional Response to Pediatric Influenza Infection | exclude | not adults |
| GSE293707 | GEO | RNA | CD55 upregulation in T cells of COVID-19 patients suppresses type-I interferon responses [ | exclude | single cell data |
| GSE293708 | GEO | RNA | CD55 upregulation in T cells of COVID-19 patients suppresses type-I interferon responses [ | exclude | single cell data |
| GSE293875 | GEO | RNA | Modeling of lung-liver interaction during infection in a human fluidic organ-on-a-chip | exclude | 9 samples: cannot reach 10 participants per class |
| GSE294031 | GEO | RNA | A next-generation inhaled aerosol COVID-19 vaccine boosts lung immunity in human | exclude | vaccination study |
| GSE294713 | GEO | RNA | Proton-secreting cells modulate mucosal immune surveillance in the male reproductive tract | exclude | not human (Mus musculus) |
| GSE294888 | GEO | RNA | Differential response of human plasmacytoid pre-dendritic cells to SARS-CoV-2 variants | exclude | no characteristic splits 10 COVID samples from 10 others and no severity field (automatic; best : 0 against 0) |
| GSE295532 | GEO | RNA | Mycophenolic acid treatment drives the emergence of novel SARS-CoV-2 variants | exclude | 12 samples: cannot reach 10 participants per class |
| GSE295798 | GEO | RNA | Immune factors produced by PBMCs upon stimulation with Lactobacillus delbrueckii ssp. bulg | exclude | in vitro or cell culture |
| GSE296228 | GEO | RNA | The robust inflammatory responses caused by MjHKU4r-S mediated cell fusion | exclude | not human |
| GSE2964 | GEO | RNA | Motexafin Gadolinium and Zinc Induce Oxidative Stress Responses and Apoptosis in B-Cell Ly | exclude | in vitro or cell culture |
| GSE296814 | GEO | RNA | Transcriptome and repertoire analysis of PBMC and CSF cells of a patient with cerebral dem | exclude | single cell data |
| GSE297329 | GEO | RNA | Evolution of SARS-CoV-2 T cell responses as a function of multiple COVID-19 boosters | exclude | single cell data |
| GSE297509 | GEO | RNA | LINC01857 regulates antigen-presenting cell viability and antimicrobial defense in humans | exclude | 19 samples: cannot reach 10 participants per class |
| GSE297529 | GEO | RNA | Discovery of Key Regulators in Classical Monocyte Phenotypes Linked to COVID-19 Severity u | exclude | single cell data |
| GSE298237 | GEO | RNA | Metabolic quiescence of naïve-like memory T cells precedes and maintains antigen-specific  | exclude | single cell data |
| GSE29828 | GEO | RNA | Expression Profiling of Mixed Lineage Leukemia Cells Treated with a Potent Small-Molecule  | exclude | in vitro or cell culture |
| GSE298709 | GEO | RNA | BAFF-R Expression as a Predictive Biomarker for COVID-19 Vaccine Non-Responsiveness in Ant | exclude | in vitro or cell culture |
| GSE29949 | GEO | RNA | Gene expression comparison among spleen dendritic cells, brain microglia, brain dendritic  | exclude | not human (Mus musculus) |
| GSE299953 | GEO | RNA | Single-cell omics analysis of age-related immune signatures induced by Omicron sublineage  | exclude | single cell data |
| GSE300265 | GEO | RNA | Single-cell Transcriptomics Reveal Signatures of Altered B cell Response and T cell Senesc | exclude | single cell data |
| GSE300696 | GEO | RNA | Identification of reactive CpGs and RNA expression in early COVID-19 through cis-eQTM anal | exclude | other assay or analyte |
| GSE301682 | GEO | RNA | An mRNA-based broad-spectrum antiviral inspired by ISG15 deficiency protects against viral | exclude | not human |
| GSE301914 | GEO | RNA | Chemically modified tRNA enhances the translation capacity of mRNA rich in cognate codons  | exclude | vaccination study |
| GSE302166 | GEO | RNA | Human brain organoids: a new model to study Cryptococcus neoformans neurotropism | exclude | not human |
| GSE302357 | GEO | RNA | Coronavirus protein interaction mapping in bat and human cells identifies molecular and ge | exclude | no characteristic splits 10 COVID samples from 10 others and no severity field (automatic; best : 0 against 0) |
| GSE302938 | GEO | RNA | Neutrophil Elastase-Mediated Spike Cleavage Promotes SARS-CoV-2 Entry and Aggravates Infla | exclude | not human |
| GSE303020 | GEO | RNA | Transcriptome analysis of cells after SARS-CoV-2 infection and treatment with  MEK1/2-inhi | exclude | in vitro or cell culture |
| GSE303208 | GEO | RNA | Single cell RNA seq of Healthy, covid-19 infected, recovered and long recovered individual | exclude | single cell data |
| GSE303324 | GEO | RNA | SARS-CoV-2 Variant Boosters and Breakthrough Infections Drive Functional Differences in th | exclude | vaccination study |
| GSE303335 | GEO | RNA | Disease-associated genetic variants can cause missense effects in tissue-specific protein  | exclude | 12 samples: cannot reach 10 participants per class |
| GSE303984 | GEO | RNA | Dosing interval is a major factor determining the quality of T cells induced by SARS-CoV-2 | exclude | single cell data |
| GSE304243 | GEO | RNA | Human B cells induced by primary SARS-CoV-2 infection | exclude | single cell data |
| GSE304583 | GEO | RNA | MAIT and other innate-like T cells integrate adaptive immune responses to modulate interva | exclude | vaccination study |
| GSE304584 | GEO | RNA | MAIT and other innate-like T cells integrate adaptive immune responses to modulate interva | exclude | vaccination study |
| GSE304643 | GEO | RNA | Characterization of SARS-CoV-2–Specific B Cell Responses in Breakthrough and Multiple Infe | exclude | vaccination study |
| GSE304891 | GEO | RNA | Inflammasome activation and oxidative stress in SARS-CoV-2 breakthrough infection and symp | exclude | single cell data |
| GSE304973 | GEO | RNA | Synthetic mRNA Vaccines and Transcriptomic Dysregulation: Evidence from New-Onset Adverse  | exclude | vaccination study |
| GSE305372 | GEO | RNA | Human lungs maintain tissue-resident memory T cells against a broad spectrum of pathogens. | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE305774 | GEO | RNA | STAT3 controls IL-6-induced epigenomic remodeling and restrains a latent STAT1-associated  | exclude | 12 samples: cannot reach 10 participants per class |
| GSE305775 | GEO | RNA | STAT3 controls IL-6-induced epigenomic remodeling and restrains a latent STAT1-associated  | exclude | 8 samples: cannot reach 10 participants per class |
| GSE307386 | GEO | RNA | Coronavirus NSP14 Drives Internal m7G Modification to Rewire Host Splicing and Promote Vir | exclude | 18 samples: cannot reach 10 participants per class |
| GSE308611 | GEO | RNA | Pulmonary fibrosis after COVID-19 is characterized by airway abnormalities and elevated cl | exclude | single cell data |
| GSE308951 | GEO | RNA | A platelet transcriptomic signature of thromboinflammation predicts cardiovascular risk | exclude | no characteristic splits 10 COVID samples from 10 others and no severity field (automatic; best : 0 against 0) |
| GSE309140 | GEO | RNA | Enhanced interferon-stimulated gene signatures in B cells following influenza B virus infe | exclude | single cell data |
| GSE309153 | GEO | RNA | RNA-seq profiling following a 5’ isomiR  transfection | exclude | other assay or analyte |
| GSE309271 | GEO | RNA | N1-Methylpseudouridine Directly Modulates Translation Dynamics | exclude | vaccination study |
| GSE309595 | GEO | RNA | Longitudinal Analysis of Multiple Myeloma: Therapeutic Response and Immune Microenvironmen | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE309865 | GEO | RNA | Effect of Sars-CoV-2 Nucleocapsid protein or Spike protein on human monocyte-derived macro | exclude | in vitro or cell culture |
| GSE3100 | GEO | RNA | Cystic Fibrosis Mouse Lung Profiles | exclude | not human (Mus musculus) |
| GSE310111 | GEO | RNA | Dual antiviral and anti-inflammatory effects of BAY 11-7082 in SARS-CoV-2 infection via in | exclude | 12 samples: cannot reach 10 participants per class |
| GSE310441 | GEO | RNA | Antigen-reactive CD4 T cells after SARS-CoV-2 vaccination show divergent phenotypic states | exclude | single cell data |
| GSE310442 | GEO | RNA | Antigen-reactive CD4 T cells after SARS-CoV-2 vaccination show divergent phenotypic states | exclude | single cell data |
| GSE311311 | GEO | RNA | UBA1-CDK16: A Female-Specific Chimeric RNA Emerging Through Evolution and Involved in Immu | exclude | 2 samples: cannot reach 10 participants per class |
| GSE311600 | GEO | RNA | SAMHD1 depletion restricts SARS-CoV-2 infection by suppressing HNF1-dependent ACE2 express | exclude | in vitro or cell culture |
| GSE311993 | GEO | RNA | Symptomatic SARS-CoV-2 breakthrough infections broaden the repertoire of Spike-reactive CD | exclude | single cell data |
| GSE312633 | GEO | RNA | Post-COVID impairment of memory T cell responses to community-acquired pathogens can be re | exclude | no characteristic splits 10 COVID samples from 10 others and no severity field (automatic; best : 0 against 0) |
| GSE312773 | GEO | RNA | Fibroblast orchestration of inflammaging via NF-kB activation | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE313868 | GEO | RNA | Transcriptomic evidence of placental dysfunction in gestational COVID associated preeclamp | exclude | not blood |
| GSE313871 | GEO | RNA | A reversible broad-spectrum antiviral targets the human V-ATPase VO domain | exclude | not human |
| GSE3140 | GEO | RNA | Gene expression profiles of adult peripheral and cord blood mononuclear cells altered by l | exclude | in vitro or cell culture |
| GSE314223 | GEO | RNA | T cell immunity to SARS-CoV-2 vaccination in inflammatory bowel disease patients treated w | exclude | single cell data |
| GSE314706 | GEO | RNA | Single-cell profiling of paired nasal brushing and tissue samples reveals distinct cellula | exclude | single cell data |
| GSE315495 | GEO | RNA | SARS-CoV-2 spike protein exerts an anti-cancer effect by regulating MEG3 and BCYRN1 in A54 | exclude | in vitro or cell culture |
| GSE315668 | GEO | RNA | Robust effector memory features of human T-bet-hi B cells  induced by repeated mRNA vaccin | exclude | vaccination study |
| GSE31773 | GEO | RNA | Comparison of mRNA expression in circulating T-cells from patients with severe asthma | exclude | other assay or analyte |
| GSE318035 | GEO | RNA | Integrative analysis reveals ciliary dysfunction and compartmentalized LCN2 expression as  | exclude | single cell data |
| GSE318050 | GEO | RNA | Identification of a human hematopoietic stem cell subset that retains memory of inflammato | exclude | single cell data |
| GSE318176 | GEO | RNA | Integrative Epigenomic and Transcriptomic Analysis of Peripheral Blood Across Different SA | exclude | single cell data |
| GSE318946 | GEO | RNA | A comparative multiomic atlas of human immunity to multiple COVID-19 vaccines | exclude | vaccination study |
| GSE320507 | GEO | RNA | Dysregulated NK-cell gene expression defines the enduring symptoms of long COVID-19 | exclude | single cell data |
| GSE320554 | GEO | RNA | A comparative multiomic atlas of human immunity to multiple COVID-19 vaccines | exclude | single cell data |
| GSE32280 | GEO | RNA | Expression data from peripheral blood lymphocytes of major depression, subsyndromal sympto | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE324489 | GEO | RNA | Whole-Protein Screening and Multi-Modal Profiling of Antigen-Specific CD4+ T Cells at Sing | exclude | single cell data |
| GSE32496 | GEO | RNA | GENOME-WIDE CpG ISLAND METHYLATION ANALYSIS IN NON-SMALL CELL LUNG CANCER PATIENTS [Affyme | exclude | in vitro or cell culture |
| GSE325054 | GEO | RNA | Time-Resolved Multi-Omics Identify Biomarkers of Immediate Reactions to mRNA Vaccination [ | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE325987 | GEO | RNA | Potent type-specific de novo antibodies complement broadly reactive imprinted antibodies i | exclude | vaccination study |
| GSE326313 | GEO | RNA | Integrative analysis reveals ciliary dysfunction and compartmentalized LCN2 expression as  | exclude | single cell data |
| GSE3268 | GEO | RNA | Squamous Lung Cancer, Paired Samples | exclude | 10 samples: cannot reach 10 participants per class |
| GSE327545 | GEO | RNA | SARS-CoV-2 nucleocapsid induces hyperinflammation and vascular leakage through the Toll-li | exclude | not blood |
| GSE328399 | GEO | RNA | Neutrophil extracellular traps aggravate COVID-19 related acute macular neuroretinopathy b | exclude | 6 samples: cannot reach 10 participants per class |
| GSE329190 | GEO | RNA | Zanubrutinib monotherapy for IgG4-related head and neck disease | exclude | single cell data |
| GSE329962 | GEO | RNA | Pathogenic IgG sialylation drives macrophage immunometabolic dysfunction and lung immunopa | exclude | single cell data |
| GSE330092 | GEO | RNA | 3D organoids containing endothelial and neural cells generation by serial inductions of di | exclude | single cell data |
| GSE330843 | GEO | RNA | The CRAC channel inhibitor Auxora interrupts inflammatory circuits between alveolar macrop | exclude | no characteristic splits 10 COVID samples from 10 others and no severity field (automatic; best : 0 against 0) |
| GSE33266 | GEO | RNA | SM001: SARS CoV MA15 infection of C57Bl/6 mouse model – Data from 4 viral doses at 1, 2, 4 | exclude | not human (Mus musculus) |
| GSE333235 | GEO | RNA | Baseline interferon signaling in monocytes and antibody-mediated innate activation at inje | exclude | not human |
| GSE333382 | GEO | RNA | Time-dependent processing of dengue virus polyprotein yields multiple capsid forms that di | exclude | not human |
| GSE334259 | GEO | RNA | A pilot study for treatment of chronic cough after COVID-19 infection by nebulization ther | exclude | other assay or analyte |
| GSE334857 | GEO | RNA | CD4⁺ T Cell Dysregulation in Long COVID: Insights from an Unvaccinated Cohort | exclude | single cell data |
| GSE336394 | GEO | RNA | CD4 T cell subset mRNA expression microarray | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE33647 | GEO | RNA | Early Peroxisome proliferator-activated receptor gamma regulated genes involved in expansi | exclude | not human (Mus musculus) |
| GSE336577 | GEO | RNA | Integrated Transcriptomic and Immune Analysis Reveals Distinct Mitochondrial and Immune Si | exclude | no characteristic splits 10 COVID samples from 10 others and no severity field (automatic; best : 0 against 0) |
| GSE338264 | GEO | RNA | Loss of cellular RNA homeostasis contributes to MDA5 activation during virus infection | exclude | single cell data |
| GSE341347 | GEO | RNA | Balancing lipid synthesis and oxidation promotes B cell response to vaccination in Solid O | exclude | single cell data |
| GSE34205 | GEO | RNA | Transcriptional profile of PBMCs in patients with acute RSV or Influenza infection | exclude | not adults |
| GSE342810 | GEO | RNA | Isoform-level analysis reveals conserved early changes in transcript usage during human va | exclude | vaccination study |
| GSE343379 | GEO | RNA | Narsoplimab (anti-MASP2) modulates long COVID patient plasma-induced transcriptional respo | exclude | 8 samples: cannot reach 10 participants per class |
| GSE34378 | GEO | RNA | Aging Experiment | exclude | not human (Mus musculus) |
| GSE344986 | GEO | RNA | Human lung organoid modeling of antiviral tissue-resident T cell responses | exclude | in vitro or cell culture |
| GSE34526 | GEO | RNA | Differential Gene Expression in Granulosa Cells from Polycystic Ovary Syndrome Patients wi | exclude | 10 samples: cannot reach 10 participants per class |
| GSE348245 | GEO | RNA | STAT3 controls IL-6-induced epigenomic remodeling and restrains a latent STAT1-associated  | exclude | 12 samples: cannot reach 10 participants per class |
| GSE348568 | GEO | RNA | Single-nucleus transcriptomic and T-cell-receptor profiling of SARS-CoV-2-infected humaniz | exclude | single cell data |
| GSE3512 | GEO | RNA | Fat-type specific effects of high fat diets | exclude | not human (Rattus norvegicus) |
| GSE35455 | GEO | RNA | Human peripheral blood cells (PB): Effect of Tocilizumab/MRA compared to methotrexate(MTX) | exclude | no characteristic splits 10 COVID samples from 10 others and no severity field (automatic; best : 0 against 0) |
| GSE36474 | GEO | RNA | Comparison of bone-marrow mesenchymal stromal cells from multiple myeloma patients and hea | exclude | 7 samples: cannot reach 10 participants per class |
| GSE36552 | GEO | RNA | Tracing pluripotency of human early embryos and embryonic stem cells by single cell RNA-se | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE36678 | GEO | RNA | Expression data SK-N-SH cells treated with risperidone | exclude | in vitro or cell culture |
| GSE36791 | GEO | RNA | Gene expression profiling of blood in ruptured intracranial aneurysms | exclude | no characteristic splits 10 COVID samples from 10 others and no severity field (automatic; best : 0 against 0) |
| GSE36807 | GEO | RNA | Genome-wide analysis of Crohn's disease and ulcerative colitis biopsy samples. | exclude | no characteristic splits 10 COVID samples from 10 others and no severity field (automatic; best : 0 against 0) |
| GSE36930 | GEO | RNA | Comparative transcriptomic profiling of peripheral blood mononuclear cells (PBMCs) isolate | exclude | no characteristic splits 10 COVID samples from 10 others and no severity field (automatic; best : 0 against 0) |
| GSE36969 | GEO | RNA | Global gene expression in response to SARS mouse adapted variants MA-15, MA-15epsilon, and | exclude | not human (Mus musculus) |
| GSE36980 | GEO | RNA | Expression data from post mortem Alzheimer's disease brains | exclude | not human |
| GSE37147 | GEO | RNA | Bronchial airway gene expression reflects a COPD-associated field of injury that changes w | exclude | single cell data |
| GSE37441 | GEO | RNA | Relief of feedback inhibition of HER3 transcription by RAF and MEK inhibitors attenuates t | exclude | in vitro or cell culture |
| GSE37827 | GEO | RNA | SCL006,icSARS CoV Urbani or icSARS Bat SRBD (spike receptor binding domain from the wild t | exclude | not human |
| GSE38136 | GEO | RNA | Hepatic gene expression in BALB/c mice fed a quercetin diet | exclude | not human (Mus musculus) |
| GSE38246 | GEO | RNA | Transcriptional response to dengue infection in peripheral blood mononuclear cells of Nica | exclude | not adults |
| GSE38376 | GEO | RNA | Acquired resistance to lapatinib | exclude | in vitro or cell culture |
| GSE38642 | GEO | RNA | Expression data from human pancreatic islets | exclude | no characteristic splits 10 COVID samples from 10 others and no severity field (automatic; best : 0 against 0) |
| GSE38713 | GEO | RNA | Expression data from intestinal mucosa of patients with UC | exclude | no characteristic splits 10 COVID samples from 10 others and no severity field (automatic; best : 0 against 0) |
| GSE38900 | GEO | RNA | Genome-wide analysis of whole blood transcriptional response to Respiratory Syncytial Viru | exclude | not adults |
| GSE39059 | GEO | RNA | Changes in microRNA and mRNA expression with differentiation of human bronchial epithelial | exclude | single cell data |
| GSE3915 | GEO | RNA | PKCa -/- vs wt small intestine | exclude | not human (Mus musculus) |
| GSE39324 | GEO | RNA | Gene expression profiling of the effects of Cyclophosphamide on peripheral blood mononucle | exclude | no characteristic splits 10 COVID samples from 10 others and no severity field (automatic; best : 0 against 0) |
| GSE3952 | GEO | RNA | Organ - specific expression profiles of rat treated with targretin, 9-cis retinoic acid, a | exclude | not human (Rattus norvegicus) |
| GSE39653 | GEO | RNA | Differential Gene Expression in Patients with Mood Disorders | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE40165 | GEO | RNA | Immune correlates of early prednisolone therapy in children and adolescents with dengue | exclude | not adults |
| GSE40184 | GEO | RNA | The blood transcriptional signature of chronic HCV [Affymetrix data] | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE40396 | GEO | RNA | Whole blood transcriptional signature distinguishes viral infection from bacterial infecti | exclude | not adults |
| GSE41789 | GEO | RNA | Senescence gene signature of radiation fibrosis | exclude | not human (Mus musculus) |
| GSE41831 | GEO | RNA | Remission in Juvenile Idiopathic Arthritis Achieved on Different Drug Regimens | exclude | not adults |
| GSE42057 | GEO | RNA | Peripheral blood mononuclear cell gene expression in chronic obstructive pulmonary disease | exclude | no characteristic splits 10 COVID samples from 10 others and no severity field (automatic; best : 0 against 0) |
| GSE42808 | GEO | RNA | Expression data from human umbilical vein endothelial cells (HUVEC) treated with DMSO, tel | exclude | 6 samples: cannot reach 10 participants per class |
| GSE43458 | GEO | RNA | Gene expression profiling of lung adenocarcinomas and normal lung tissue | exclude | single cell data |
| GSE43695 | GEO | RNA | Comparative analysis of gene expression in Fra-1+/+ and Fra-1-/- mice lung | exclude | not human (Mus musculus) |
| GSE4386 | GEO | RNA | Transcriptomics in cardiac surgery | exclude | no characteristic splits 10 COVID samples from 10 others and no severity field (automatic; best : 0 against 0) |
| GSE45764 | GEO | RNA | Effect of influenza vaccination on PBMC and B cell gene expression profiles in healthy hum | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE45919 | GEO | RNA | Peripheral blood gene expression in human experiencing primary EBV infection (HT12) | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE4607 | GEO | RNA | Systemic inflammatory response syndrome and septic shock | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE4619 | GEO | RNA | Gene expression profiling of CD34+ cells from MDS patients and normal controls | exclude | no characteristic splits 10 COVID samples from 10 others and no severity field (automatic; best : 0 against 0) |
| GSE4668 | GEO | RNA | Expression data from hormone-starved MCF7/BUS cell culture | exclude | no characteristic splits 10 COVID samples from 10 others and no severity field (automatic; best : 0 against 0) |
| GSE47189 | GEO | RNA | Transcriptome-based network analysis reveals a spectrum model of human macrophage activati | exclude | single cell data |
| GSE47460 | GEO | RNA | Gene expression profiling of chronic lung disease for the Lung Genomics Research Consortiu | exclude | single cell data |
| GSE47516 | GEO | RNA | Gene expression alterations in the cerebellum and granule neurons of Cstb-/- mouse are ass | exclude | not human (Mus musculus) |
| GSE48438 | GEO | RNA | Expression data from osteoblastic lineage cells isolated from normal and leukemic mice | exclude | not human (Mus musculus) |
| GSE48466 | GEO | RNA | Expression data from well-differentiated human bronchial epithelial cells infected with H1 | exclude | in vitro or cell culture |
| GSE49126 | GEO | RNA | Transcriptomic profiling of peripheral blood mononuclear cells from patients with Parkinso | exclude | no characteristic splits 10 COVID samples from 10 others and no severity field (automatic; best : 0 against 0) |
| GSE49414 | GEO | RNA | Gene expression profiling in TPC1 cell line treated with the RTK inhibitor RPI-1 | exclude | in vitro or cell culture |
| GSE50000 | GEO | RNA | SM003 - icSARS CoV, SARS MA15 wild type and SARS BatSRBD mutant virus infections of C57BL6 | exclude | not human (Mus musculus) |
| GSE5007 | GEO | RNA | Expression data for LMF | exclude | 6 samples: cannot reach 10 participants per class |
| GSE5090 | GEO | RNA | PCOS patients vs control subjects | exclude | 17 samples: cannot reach 10 participants per class |
| GSE51392 | GEO | RNA | Expression data from airway epithelial cells stimulated with Poly(I:C) from patients with  | exclude | single cell data |
| GSE52428 | GEO | RNA | Host gene expression signatures of influenza A H1N1 and H3N2 virus infection in adults | exclude | no characteristic splits 10 COVID samples from 10 others and no severity field (automatic; best : 0 against 0) |
| GSE52471 | GEO | RNA | Dominant Th1 and Minimal Th17 Skewing in Discoid Lupus Revealed by Transcriptomic Comparis | exclude | no characteristic splits 10 COVID samples from 10 others and no severity field (automatic; best : 0 against 0) |
| GSE52819 | GEO | RNA | Vitamin D treatment of M.tb. infected macrophages | exclude | 12 samples: cannot reach 10 participants per class |
| GSE53960 | GEO | RNA | A rat RNA-Seq transcriptomic Bodymap across eleven organs and four developmental stages | exclude | not human (Rattus norvegicus) |
| GSE53987 | GEO | RNA | Microarray profiling of PFC, HPC and STR from subjects with schizophrenia, bipolar, MDD or | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE54536 | GEO | RNA | Involvement of endocytosis and alternative splicing in the formation of the pathological p | exclude | 10 samples: cannot reach 10 participants per class |
| GSE54562 | GEO | RNA | Expression data from human brain anterior cingulate cortex - including control samples and | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE54563 | GEO | RNA | Expression data from human brain anterior cingulate cortex - including control samples and | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE54564 | GEO | RNA | Expression data from human brain amygdala - including control samples and samples with maj | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE54565 | GEO | RNA | Expression data from human brain anterior cingulate cortex - including control samples and | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE54566 | GEO | RNA | Expression data from human brain amygdala - including control samples and samples with maj | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE54567 | GEO | RNA | Expression data from human brain dorsolateral prefrontal cortex - including control sample | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE54568 | GEO | RNA | Expression data from human brain dorsolateral prefrontal cortex - including control sample | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE54570 | GEO | RNA | Expression data from human brain dorsolateral prefrontal cortex - including control sample | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE54571 | GEO | RNA | Expression data from human brain anterior cingulate cortex - including control samples and | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE54572 | GEO | RNA | Expression data from human brain anterior cingulate cortex - including control samples and | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE54575 | GEO | RNA | Expression data from human brain orbital ventral prefrontal cortex - including control sam | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE54593 | GEO | RNA | Expression data from Arabidopsis leaves | exclude | not human (Arabidopsis thaliana) |
| GSE5504 | GEO | RNA | human peripheral blood derived monocytes, LPS stimulation time-series | exclude | 10 samples: cannot reach 10 participants per class |
| GSE55851 | GEO | RNA | Whole genome gene expression profiling in CD4+ T-cells, FACS-sorted based on CADM1/CD7 exp | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE56192 | GEO | RNA | Transcriptomic analysis of the Novel Middle East Respiratory Syndrome Coronavirus (Human,  | exclude | in vitro or cell culture |
| GSE5679 | GEO | RNA | Comparative gene expression profile of PPARg and RARa ligand treated human dendritic cells | exclude | 11 samples: cannot reach 10 participants per class |
| GSE57065 | GEO | RNA | Early and dynamic changes in gene expression in septic shock patients: a genome wide appro | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE57148 | GEO | RNA | Characterizing gene expression in lung tissue of COPD subjects using RNA-seq | exclude | single cell data |
| GSE58208 | GEO | RNA | Gene expression profiling of PBMC from normal individuals, chronic hepatitis B carriers an | exclude | no characteristic splits 10 COVID samples from 10 others and no severity field (automatic; best : 0 against 0) |
| GSE58294 | GEO | RNA | Gene Expression Following Cardioembolic Stroke | exclude | no characteristic splits 10 COVID samples from 10 others and no severity field (automatic; best : 0 against 0) |
| GSE59114 | GEO | RNA | Single cell RNA_Seq in aging of hematopoietic stem and progenitor cells | exclude | not human (Mus musculus) |
| GSE59635 | GEO | RNA | Gene expression analysis of Influenza vaccine response in Young and Old - Year 1 - PBMC | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE59654 | GEO | RNA | Gene expression analysis of Influenza vaccine response in Young and Old - Year 2 - PBMC | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE5972 | GEO | RNA | Gene expression profiling of patients with severe acute respiratory syndrome (SARS) | exclude | no characteristic splits 10 COVID samples from 10 others and no severity field (automatic; best : 0 against 0) |
| GSE59743 | GEO | RNA | Gene expression analysis of Influenza vaccine response in Young and Old - Year 3 - PBMC | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE60993 | GEO | RNA | Assessment and Diagnostic Relevance of Novel Serum Biomarkers for Early Decision of ST-ele | exclude | no characteristic splits 10 COVID samples from 10 others and no severity field (automatic; best : 0 against 0) |
| GSE61608 | GEO | RNA | Transcriptional Landscape of B cell Early Activation | exclude | not human (Mus musculus) |
| GSE61754 | GEO | RNA | Transcriptomic profiling facilitates classification of response to influenza challenge | exclude | vaccination study |
| GSE62191 | GEO | RNA | Gene expression profiles of patients with schizophrenia, bipolar disorder and healthy cont | exclude | no characteristic splits 10 COVID samples from 10 others and no severity field (automatic; best : 0 against 0) |
| GSE6264 | GEO | RNA | Transcriptional Profiling of Subjects from the Iowa Adoption Studies | exclude | in vitro or cell culture |
| GSE62646 | GEO | RNA | Altered gene expression pattern in peripheral blood mononuclear cells in patients with acu | exclude | no characteristic splits 10 COVID samples from 10 others and no severity field (automatic; best : 0 against 0) |
| GSE6269 | GEO | RNA | Gene expression patterns in blood leukocytes discriminate patients with acute infections | exclude | not adults |
| GSE63473 | GEO | RNA | Highly parallel genome-wide expression profiling of individual cells using nanoliter dropl | exclude | no characteristic splits 10 COVID samples from 10 others and no severity field (automatic; best : 0 against 0) |
| GSE64016 | GEO | RNA | Single cell RNA-seq data of human hESCs to evaluate Oscope - a statistical pipeline for id | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE64456 | GEO | RNA | Defining RNA Transcriptional Biosignatures to Distinguish Febrile Infants 60 Days of Age a | exclude | not adults |
| GSE64614 | GEO | RNA | Distal-to-Proximal Re-patterning of Small Airway Epithelium in Smoking-associated Chronic  | exclude | single cell data |
| GSE6511 | GEO | RNA | Four weeks expression data of the antipsychotics Clozapine and Haloperidol in the mouse br | exclude | not human (Mus musculus) |
| GSE65213 | GEO | RNA | Leukocyte gene expression as a function of perceived social isolation (PSI): Year 5 | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE65298 | GEO | RNA | Leukocyte gene expression as a function of perceived social isolation (PSI): Year 7 | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE65317 | GEO | RNA | Leukocyte gene expression as a function of perceived social isolation (PSI): Year 8 | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE65341 | GEO | RNA | Leukocyte gene expression as a function of perceived social isolation (PSI): Year 9 | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE65403 | GEO | RNA | Leukocyte gene expression as a function of perceived social isolation (PSI): Year 10 | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE65517 | GEO | RNA | Comparative analysis of monocytes from healthy donors, patients with metastatic breast can | exclude | single cell data |
| GSE6591 | GEO | RNA | Comparative analysis of gene expression in C57BL/6J and DBA/2J aging lungs. | exclude | not human (Mus musculus) |
| GSE66360 | GEO | RNA | A Whole Blood Molecular Signature for the Identification of Acute Myocardial Infarction Wi | exclude | no characteristic splits 10 COVID samples from 10 others and no severity field (automatic; best : 0 against 0) |
| GSE66494 | GEO | RNA | Development of gene expression profiles in human chronic kidney disease | exclude | not blood |
| GSE66936 | GEO | RNA | Expression data from intermediate monocytes from healthy donors and autoimmune uveitis pat | exclude | no characteristic splits 10 COVID samples from 10 others and no severity field (automatic; best : 0 against 0) |
| GSE67059 | GEO | RNA | Whole Blood Transcriptional Profiling Differentiates Between Asymptomatic and Symptomatic  | exclude | not adults |
| GSE6731 | GEO | RNA | Genome-wide gene expression differences between Crohn’s and ulcerative colitis from endosc | exclude | no characteristic splits 10 COVID samples from 10 others and no severity field (automatic; best : 0 against 0) |
| GSE67368 | GEO | RNA | Changes in gene expression after epigentic targeted therapy in T-cell prolymphocytic leuke | exclude | 12 samples: cannot reach 10 participants per class |
| GSE67472 | GEO | RNA | Airway epithelial gene expression in asthma versus healthy controls | exclude | single cell data |
| GSE67561 | GEO | RNA | Gene expression profiling of T47D breast cancer cells following progestin administration i | exclude | in vitro or cell culture |
| GSE6802 | GEO | RNA | Gene expression analysis of bronchial epithelial cells stimulated by different airway path | exclude | not blood |
| GSE68526 | GEO | RNA | Peripheral blood transcriptome profiles from an RNA Pilot Study within the United States H | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE68820 | GEO | RNA | Mouse lung tissue transcriptome response to a mouse-adapted strain of SARS-CoV in wild typ | exclude | not human (Mus musculus) |
| GSE6887 | GEO | RNA | Gene expression profile of peripheral blood lymphocytes: comparison between melanoma patie | exclude | no characteristic splits 10 COVID samples from 10 others and no severity field (automatic; best : 0 against 0) |
| GSE69529 | GEO | RNA | Elucidating the etiology and molecular pathogenicity of infectious diarrhea by high throug | exclude | not adults |
| GSE70528 | GEO | RNA | Gene expression annalysis of peripheral blood cells in patients with chronic kidney diseas | exclude | 19 samples: cannot reach 10 participants per class |
| GSE70630 | GEO | RNA | Single cell RNA-seq analysis of oligodendroglioma | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE7148 | GEO | RNA | Differential gene expression in peripheral blood leukocytes from CHASRS particpants | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE72492 | GEO | RNA | Gene expression in the pancreas of healthy control, auto-antibody positive, and type 1 dia | exclude | no characteristic splits 10 COVID samples from 10 others and no severity field (automatic; best title: 6 against 7) |
| GSE7261 | GEO | RNA | Yeast mRNA decay analysis after addition of 3 µ/mL thiolutin | exclude | not human (Saccharomyces cerevisiae) |
| GSE72748 | GEO | RNA | Transcriptomes of peripheral blood mononuclear cells from a Guillain-Barre Syndrome patien | exclude | 6 samples: cannot reach 10 participants per class |
| GSE73072 | GEO | RNA | Host gene expression signatures of H1N1, H3N2, HRV, RSV virus infection in adults | exclude | no characteristic splits 10 COVID samples from 10 others and no severity field (automatic; best : 0 against 0) |
| GSE73461 | GEO | RNA | Diagnosis of Kawasaki Disease in children using host RNA expression [Discovery_Dataset] | exclude | not adults |
| GSE73685 | GEO | RNA | Unique inflammatory transcriptome profiles at the maternal fetal interface and onset of hu | exclude | not adults |
| GSE75214 | GEO | RNA | Mucosal gene expression profiling in patients with inflammatory bowel disease study | exclude | no characteristic splits 10 COVID samples from 10 others and no severity field (automatic; best : 0 against 0) |
| GSE7524 | GEO | RNA | Expression profiling in RA disease pre and post anti-TNF treatment | exclude | 6 samples: cannot reach 10 participants per class |
| GSE75715 | GEO | RNA | TAp73 is a master regulator of airway multiciliogenesis [RNA-Seq] | exclude | not human (Mus musculus) |
| GSE76293 | GEO | RNA | Neutrophils from ARDS patients and healthy volunteers, and healthy volunteer samples treat | exclude | no characteristic splits 10 COVID samples from 10 others and no severity field (automatic; best : 0 against 0) |
| GSE76492 | GEO | RNA | Mitochondrial and oxidative stress genes are differentially expressed in neutrophils of sJ | exclude | 14 samples: cannot reach 10 participants per class |
| GSE76826 | GEO | RNA | Blood transcriptomic markers in patients with late-onset major depressive disorder | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE77087 | GEO | RNA | Nasopharyngeal microbiota, host transcriptome and disease severity in children with respir | exclude | not adults |
| GSE77939 | GEO | RNA | Gene expression data from immunological non-responders and responders to antiretroviral th | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE77962 | GEO | RNA | Adipose tissue gene expression is differentially regulated with different rates of weight  | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE7868 | GEO | RNA | Expression data from LNCaP cell line | exclude | in vitro or cell culture |
| GSE79362 | GEO | RNA | A blood RNA signature for tuberculosis disease risk: a prospective cohort study | exclude | no characteristic splits 10 COVID samples from 10 others and no severity field (automatic; best ch_qft: 86 against 6) |
| GSE80680 | GEO | RNA | Global transcriptome analysis reveals aging-related genes in Atp8b1 mutant lungs | exclude | not human (Mus musculus) |
| GSE81246 | GEO | RNA | Efficient immune responses in immunocompetent individuals developing symptomatic CMV infec | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE81252 | GEO | RNA | Multilineage communication regulates human liver bud self-organization from pluripotency I | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE81682 | GEO | RNA | The transcriptional landscape of mouse blood stem/progenitor cell transitions at single ce | exclude | not human (Mus musculus) |
| GSE83233 | GEO | RNA | An inhaled dose of budesonide induces genes involved in transcription and signaling in hum | exclude | not blood |
| GSE84133 | GEO | RNA | A single-cell transcriptomic map of the human and mouse pancreas reveals inter- and intra- | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE8597 | GEO | RNA | Gene expression analysis of hormone treated MCF7 breast cancer cells in the presence or ab | exclude | 16 samples: cannot reach 10 participants per class |
| GSE86623 | GEO | RNA | Mouse strain-specific chromatin and transcriptional states modulate the DNA-damage respons | exclude | not human (Mus musculus) |
| GSE8762 | GEO | RNA | Lymphocyte gene expression data from moderate stage HD patients and controls | exclude | no characteristic splits 10 COVID samples from 10 others and no severity field (automatic; best : 0 against 0) |
| GSE8835 | GEO | RNA | Chronic lymphocytic leukemia cells induce changes in gene expression of CD4 and CD8 T cell | exclude | no characteristic splits 10 COVID samples from 10 others and no severity field (automatic; best : 0 against 0) |
| GSE88966 | GEO | RNA | Depot dependent effects of dexamethasone on gene expression in human omental and abdominal | exclude | no characteristic splits 10 COVID samples from 10 others and no severity field (automatic; best : 0 against 0) |
| GSE89232 | GEO | RNA | Single-cell RNA-seq reveals heterogeneity within human pre-cDCs | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE89970 | GEO | RNA | Genome-wide map of HAEC chromatin landscape under resting and TNFa, IL1b, and OxPAPC stimu | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE9006 | GEO | RNA | Gene expression in PBMCs from children with diabetes | exclude | not adults |
| GSE90074 | GEO | RNA | Gene expression data from Phase 2 of the SAMARA study (Supporting a Multi-disciplinary App | exclude | no characteristic splits 10 COVID samples from 10 others and no severity field (automatic; best : 0 against 0) |
| GSE90748 | GEO | RNA | Mycobacterium bovis Bacillus Calmette-Guerin (BCG) reprograms melanoma microenvironment to | exclude | in vitro or cell culture |
| GSE9128 | GEO | RNA | Expression data from heart failure vs control peripheral blood mononuclear cells. | exclude | 11 samples: cannot reach 10 participants per class |
| GSE92332 | GEO | RNA | A single-cell survey of the small intestinal epithelium | exclude | not human (Mus musculus) |
| GSE92415 | GEO | RNA | Characterization of molecular response to Golimumab in Ulcerative Colitis by mucosal biops | exclude | no characteristic splits 10 COVID samples from 10 others and no severity field (automatic; best : 0 against 0) |
| GSE93284 | GEO | RNA | Alterations in bronchial airway microRNA expression for lung cancer detection | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE94820 | GEO | RNA | Single-cell RNA-seq reveals new types of human blood dendritic cells, monocytes and progen | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE9509 | GEO | RNA | A transcriptional repressor and co-repressor induced by the STAT3-regulated anti-inflammat | exclude | not human (Mus musculus) |
| GSE95233 | GEO | RNA | Fractalkine receptor CX3CR1 and leukocyte Ig-like receptor B2 LILRB2 are prognostic biomar | exclude | no characteristic splits 10 COVID samples from 10 others and no severity field (automatic; best : 0 against 0) |
| GSE9686 | GEO | RNA | Human colon expression in healthy, CD, treated CD, and UC | exclude | not adults |
| GSE97104 | GEO | RNA | Molecular transitions in early progenitors during human cord blood hematopoiesis | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE97320 | GEO | RNA | Differential gene expression profiles in peripheral blood in Northeast Chinese Han people  | exclude | 6 samples: cannot reach 10 participants per class |
| GSE97590 | GEO | RNA | Next Generation Sequencing RNA Seq data from UKE Phase I rVSV ZEBOV vaccine clinical trial | exclude | vaccination study |
| GSE98201 | GEO | RNA | Medial Ganglionic Eminence and Cortical Organoids Model Human Brain Development and Intern | exclude | single cell data |
| GSE98622 | GEO | RNA | Transcriptome profile of a murine renal bilateral ischemia reperfusion model 2 hours to 12 | exclude | not human (Mus musculus) |
| GSE98793 | GEO | RNA | Gene expression analysis in whole blood samples obtained from donors diagnosed with MDD co | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE99254 | GEO | RNA | T cell landscape of non-small cell lung cancer revealed by deep single-cell RNA sequencing | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE99992 | GEO | RNA | Whole blood RNA-seq of pediatric cases of natural Chikungunya infection | exclude | not adults |
| 10.17632/3jdywv5jrv | Mendeley Data | RNA | Profiling B cell immunodominance after SARS-CoV-2 infection reveals antibody evolution to  | exclude | single cell data |
| 10.17632/8bhhk7c5n9 | Mendeley Data | RNA | Whole-cortex in situ sequencing reveals peripheral input-dependent cellular and area ident | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| 10.17632/9sw3trkv8p | Mendeley Data | RNA | Complement activation induces excessive T cell cytotoxicity in severe COVID-19. Georg et a | exclude | single cell data |
| 10.17632/9vmdtjnh4k | Mendeley Data | RNA | Cross-Neutralization of Emerging SARS-CoV-2 Variants of Concern by Antibodies Targeting Di | exclude | single cell data |
| 10.17632/bvd647g3mc | Mendeley Data | RNA | Early IFN-a signatures and persistent dysfunction are distinguishing features of NK cells  | exclude | single cell data |
| 10.17632/csc9p34d5t | Mendeley Data | RNA | Mucosal Associated Invariant T (MAIT) Cell Responses Differ by Sex in COVID-19 by Yu et al | exclude | single cell data |
| 10.17632/dvp4y5ttd5 | Mendeley Data | RNA | COVID-19 immune features revealed by a large-scale single cell transcriptome atlas. Ren et | exclude | single cell data |
| 10.17632/hwxhw2sxys | Mendeley Data | RNA | Severe COVID-19 is marked by a dysregulated myeloid cell compartment. Schulte-Schrepping e | exclude | single cell data |
| 10.17632/mkmv3mgz86 | Mendeley Data | RNA | Limited TCR repertoire and ENTPD1 dysregulation mark late-stage COVID-19, Wang et al. | exclude | other assay or analyte |
| 10.17632/vrxdg7mm6x | Mendeley Data | RNA | Implicating Gene and Cell Networks Responsible for Differential COVID-19 Host Responses vi | exclude | single cell data |
| zenodo.10246007 | Zenodo | RNA | scRNA-seq revealed the rules for CDR3 length pairing in TCR beta and alpha chains and BCR  | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| zenodo.10405106 | Zenodo | RNA | Data for "Profiling the transcriptomic age of single-cells in humans" | exclude | single cell data |
| zenodo.10795478 | Zenodo | RNA | Data from: Association between gut health and gut microbiota in a polluted environment | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| zenodo.10916993 | Zenodo | RNA | Gene Expression Signatures in Blood from a West African Sepsis Cohort Define Host Response | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| zenodo.11067349 | Zenodo | RNA | Comparative Network Analysis of Single-Cell RNA Sequencing Data Unveils Deep Signatures of | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| zenodo.11100300 | Zenodo | RNA | Longitudinal single-cell RNA sequencing of 55,260 capillary PBMCs from 2 donors with ME/CF | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| zenodo.11103555 | Zenodo | RNA | Benchmark datasets for testing AIRR-seq data processing with pyIgMap pipeline | exclude | single cell data |
| zenodo.11207128 | Zenodo | RNA | Data from: Dissecting gene activation and chromatin remodeling dynamics in single human ce | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| zenodo.11217984 | Zenodo | RNA | Gene expression of PFAS exposed human liver spheroids and quality control | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| zenodo.11238538 | Zenodo | RNA | Reconstitution of the multiple myeloma microenvironment following lymphodepletion with BCM | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| zenodo.12536701 | Zenodo | RNA | Development of human pancreatic cancer avatars as a model for dynamic immune landscape pro | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| zenodo.12626634 | Zenodo | RNA | Single cell multiomic analysis identifies key genes differentially expressed in innate lym | exclude | single cell data |
| zenodo.12643405 | Zenodo | RNA | Data from: Impacts of weathered microplastic ingestion on gastrointestinal microbial commu | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| zenodo.14117450 | Zenodo | RNA | Defective interferon signaling in the circulating monocytes of type 2 diabetic mice | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| zenodo.14619808 | Zenodo | RNA | Gene expression profiling of dexamethasone effects in multiple tissues | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| zenodo.14872873 | Zenodo | RNA | Germinal center-mediated broadening of B cell responses to SARS-CoV-2 booster immunization | exclude | single cell data |
| zenodo.15147562 | Zenodo | RNA | Biological modifications of the immune response to COVID-19 vaccine in patients treated wi | exclude | vaccination study |
| zenodo.15385493 | Zenodo | RNA | Analysis of COVID-19 Patient PBMCs Using Single-cell Multiome Profiling | exclude | single cell data |
| zenodo.15428751 | Zenodo | RNA | Intragenomic conflict associated with extreme phenotypic plasticity in queen-worker caste  | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| zenodo.15711415 | Zenodo | RNA | Targeting ARPC1B+ Cancer Stem Cells to Sensitize Pancreatic Cancer to Gemcitabine Treatmen | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| zenodo.16959681 | Zenodo | RNA | XChrom: a cross-cell chromatin accessibility prediction model integrating genomic sequence | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| zenodo.18047227 | Zenodo | RNA | Supplementary Material for "stCAMBL: Biased Multi-view Contrastive Learning with Attentive | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| zenodo.18850871 | Zenodo | RNA | Human Hematopoietic Stem Cells Remember Inflammatory Stress  -  Analysis Notebooks & Datas | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| zenodo.18945379 | Zenodo | RNA | Dataset and Codes for "Analysis of Plasma Extracellular Vesicles in Normal-Weight and Over | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| zenodo.19791315 | Zenodo | RNA | COSMIC Single-Cell and Spatial Atlas: Lung, Blood, Nasal, and Imaging Mass Cytometry Data  | exclude | single cell data |
| zenodo.20026819 | Zenodo | RNA | Derived tables and supplementary data: Integrative human transcriptomics and composition-a | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| zenodo.20053750 | Zenodo | RNA | Multimodal analysis reveals polyclonality and inflammatory programs specific to TRBV11-2+  | exclude | single cell data |
| zenodo.20113753 | Zenodo | RNA | Suco - PBMC | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| zenodo.20337344 | Zenodo | RNA | Single-cell multiomics reveals microglia-neuron granulin signaling associated with brain a | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| zenodo.20631946 | Zenodo | RNA | Integrating Mutation, Methylation, and Expression Data in Acute Myeloid Leukemia | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| zenodo.21398550 | Zenodo | RNA | Predicting Long COVID Persistence Using Mitochondrial Gene Expression Profiles: An XGBoost | exclude | post-acute only (no acute samples by its record) |
| zenodo.23693856 | Zenodo | RNA | Multi-layer evidence matrix for the 123 shared genes of Ménière's disease and systemic lup | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| zenodo.24155547 | Zenodo | RNA | PROACT: Prior-informed inference of context-specific transcription factor cooperation and  | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| zenodo.3357167 | Zenodo | RNA | A comparison of automatic cell identification methods for single-cell RNA-sequencing data | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| zenodo.3842289 | Zenodo | RNA | JTI | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| zenodo.4312046 | Zenodo | RNA | Data from: Uncovering genetic mechanisms of hypertension through multi-omic analysis of th | exclude | other assay or analyte |
| zenodo.4408458 | Zenodo | RNA | Gestational Cd exposure in the CD-1 mouse induces sex-specific hepatic insulin insensitivi | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| zenodo.4938413 | Zenodo | RNA | Data from: "Blood transcriptome sequencing of Common Crossbills (Loxia curvirostra) experi | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| zenodo.4978699 | Zenodo | RNA | Data from: "Polar bear (Ursus maritimus) transcriptome assembly and SNP discovery" in Geno | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| zenodo.4998681 | Zenodo | RNA | Data from: Genetic basis for variation in salinity tolerance between stickleback ecotypes | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| zenodo.5148862 | Zenodo | RNA | Single-cell epitope and RNA-sequencing reveals the systemic immune associations with type  | exclude | single cell data |
| zenodo.5153528 | Zenodo | RNA | Longitudinal Single Cell Atlas Identifies Type I Interferon Response as an Early Transient | exclude | single cell data |
| zenodo.5172943 | Zenodo | RNA | memento: Hypothesis testing of single cell RNA-sequencing data with method of moments esti | exclude | single cell data |
| zenodo.5564625 | Zenodo | RNA | moFluMemB - Dataset : scRNA-seq from Lymph node, Spleen and Lung - 10x_190712_m_moFluMemB | exclude | single cell data |
| zenodo.5565864 | Zenodo | RNA | moFluMemB - Dataset : scRNA-seq from Lymph node, Spleen and Lung - 10x_191105_m_moFluMemB | exclude | single cell data |
| zenodo.5566675 | Zenodo | RNA | moFluMemB - Dataset : scRNA-seq from Lymph node, Spleen and Lung | exclude | single cell data |
| zenodo.5895181 | Zenodo | RNA | Germinal centre-driven maturation of B cell response to SARS-CoV-2 mRNA vaccination | exclude | single cell data |
| zenodo.7506637 | Zenodo | RNA | Down-regulation of MALAT1 is a hallmark of tissue and peripheral proliferative T cells in  | exclude | single cell data |
| zenodo.7543245 | Zenodo | RNA | Data for: Corticosteroid-binding globulin (SERPINA6) establishes postpubertal sex differen | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| zenodo.7577699 | Zenodo | RNA | Influence of RNA-Seq library construction, sampling methods, and tissue harvesting time on | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| zenodo.7719030 | Zenodo | RNA | SARS-CoV-2 Omicron Boosting Induces De Novo B Cell Response in Humans | exclude | other assay or analyte |
| zenodo.7851053 | Zenodo | RNA | Genetic control of the dynamic transcriptional response to immune stimuli and glucocortico | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| zenodo.7983426 | Zenodo | RNA | Source code repo for `Cell-type-specific co-expression inference from single cell RNA-sequ | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| zenodo.7983559 | Zenodo | RNA | Additional raw data in `Cell-type-specific co-expression inference from single cell RNA-se | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| zenodo.8005706 | Zenodo | RNA | Supplementary datasets: sciCSR infers B cell state transition and predicts class-switch re | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| zenodo.8060366 | Zenodo | RNA | Dataset related to article "Transcriptomic profile of TNFhigh MAIT cells is linked to B ce | exclude | single cell data |
| zenodo.8305391 | Zenodo | RNA | Patterns of multiomic chromatin accessibility and gene expression in peripheral blood leuk | exclude | single cell data |
| zenodo.8417887 | Zenodo | RNA | Spatiotemporal Profiling Unveiling the Cellular Organization Patterns and Local Protumoral | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| EGAD00001006550 | controlled access archive | RNA | Single-cell transcriptomics reveal a unique memory-like NK cell subset that accumulates wi | exclude | not open access (application needed) |
| EGAD00001007982 | controlled access archive | RNA | Epigenetic and transcriptomic reprogramming in monocytes of severe COVID-19 patients refle | exclude | not open access (application needed) |
| EGAD00001008015 | controlled access archive | RNA | The contribution of genetic determinants of blood gene expression and splicing to molecula | exclude | not open access (application needed) |
| EGAD50000001405 | controlled access archive | RNA | Transcriptome changes in circulating immune cells of critical COVID-19 patients predict a  | exclude | not open access (application needed) |
| EGAD50000001406 | controlled access archive | RNA | Transcriptome changes in circulating immune cells of critical COVID-19 patients predict a  | exclude | not open access (application needed) |
| EGAS0000001004571 | controlled access archive | RNA | Identification of COVID-19-Associated DNA Methylation Variations by Integrating Methylatio | exclude | not open access (application needed) |
| EGAS00001003215 | controlled access archive | RNA | Immune cell identifier and classifier (ImmunIC) for single cell transcriptomic readouts. | exclude | not open access (application needed) |
| EGAS00001004481 | controlled access archive | RNA | Bulk and Single-Cell RNA Sequencing Elucidate the Etiology of Severe COVID-19. | exclude | not open access (application needed) |
| syn21041850 | controlled access archive | RNA | A gene regulatory network-aware graph learning method for cell identity annotation in sing | exclude | not open access (application needed) |
| syn21261143 | controlled access archive | RNA | Cell-type-specific co-expression inference from single cell RNA-sequencing data. | exclude | not open access (application needed) |
| syn22079621 | controlled access archive | RNA | Cell-type-specific co-expression inference from single cell RNA-sequencing data. | exclude | not open access (application needed) |
| E-MTAB-11671 | other | RNA | Diagnosis of Multisystem Inflammatory Syndrome in Children by a Whole-Blood Transcriptiona | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| E-MTAB-3929 | other | RNA | Profiling the transcriptomic age of single-cells in humans. | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| E-MTAB-9005 | other | RNA | Immune cell identifier and classifier (ImmunIC) for single cell transcriptomic readouts. | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| E-MTAB-9467 | other | RNA | Intracellular microbial shifts during COVID-19 infection and longitudinal recovery reveale | exclude | single cell data |
| PRJNA1335791 | BioProject | metabolites | Multiple-omic analyses reveals limited vertical transmission but alter placenta immune mic | exclude | single cell data |
| PRJNA680599 | BioProject | metabolites | Metabolic reprogramming and epigenetic changes of vital organs in SARS-CoV-2 induced syste | exclude | duplicate: BioProject of GSE162113 |
| PRJNA756623 | BioProject | metabolites | organismal metagenomes Raw sequence reads | exclude | not human (organismal metagenomes) |
| PRJNA828431 | BioProject | metabolites | Personalized system-wide multiomics analysis of COVID19 diease | exclude | single cell data |
| PRJNA905565 | BioProject | metabolites | Alterations in the gut microbiome and metabolome profiles of septic mice treated with Shen | exclude | not human |
| GSE168797 | GEO | metabolites | Cannabidiol Inhibits SARS-CoV-2 Replication through Induction of the Host ER Stress and In | exclude | not human |
| GSE32418 | GEO | metabolites | Integrated analysis of the transcriptome-lipidome reveals the interactions of inositol-cho | exclude | not human (Saccharomyces cerevisiae; Schizosaccharomyces pombe) |
| GSE46924 | GEO | metabolites | 27-Hydroxycholesterol links cholesterol and breast cancer pathophysiology. | exclude | not human |
| GSE7002 | GEO | metabolites | Gene Expression Changes in the Rat Nasal Epithelium Following Formaldehyde Exposure | exclude | not human (Rattus norvegicus) |
| GSE7664 | GEO | metabolites | Analysis of human cells response to benzene metabolites | exclude | 16 samples: cannot reach 10 participants per class |
| 10.17632/sjnv8kjxpb | Mendeley Data | metabolites | Clinical and metabolomics characterization of hospitalized COVID-19 patients during the De | exclude | vaccination study |
| 10.17632/tzydswhhb5 | Mendeley Data | metabolites | Multi-omics resolves a sharp disease-state shift between mild and moderate COVID-19. Su et | exclude | duplicate: Su et al. 2020 INCOV cohort, the participants of Su2022 already in CHORUS |
| ST001853 | Metabolomics Workbench | metabolites | Longitudinal Metabolomics of Human Plasma Reveals Robust Prognostic Markers of COVID-19 Di | exclude | not human (Mesocricetus auratus) |
| ST001921 | Metabolomics Workbench | metabolites | An Airway Organoid-Based Screen Identifies a Role for the HIF1α-Glycolysis Axis in SARS-Co | exclude | in vitro or cell culture |
| ST002086 | Metabolomics Workbench | metabolites | Serum NMR profiling reveals differential alterations in the lipoproteome induced by Pfizer | exclude | vaccination study |
| ST002476 | Metabolomics Workbench | metabolites | High body temperature increases gut microbiota-dependent host resistance to influenza A vi | exclude | not human (Mus musculus) |
| ST002479 | Metabolomics Workbench | metabolites | High body temperature increases gut microbiota-dependent host resistance to influenza A vi | exclude | not human (Mesocricetus Auratus) |
| ST003039 | Metabolomics Workbench | metabolites | A Non-Targeted Metabolomics Comparative Study on Plasma of Pfizer and Sinopharm COVID- 19  | exclude | vaccination study |
| ST003103 | Metabolomics Workbench | metabolites | Reinforcing the Evidence of Mitochondrial Dysfunction in Long COVID Patients using a Multi | exclude | post-acute only (no acute samples by its record) |
| ST004607 | Metabolomics Workbench | metabolites | Longitudinal NMR-Based Metabolomic and Immune Profiling of Perfluorooctanoic Acid (PFOA) E | exclude | not human (Mustela putorius furo) |
| ST004667 | Metabolomics Workbench | metabolites | Metabolomics analysis of human primary astrocytes treated with astrocyte-derived exosomes  | exclude | other assay or analyte |
| zenodo.10784399 | Zenodo | metabolites | NMR based clinical metabolomics study revealed aberrant proline biosynthesis and mitochond | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| zenodo.14044561 | Zenodo | metabolites | High-throughput untargeted metabolomics reveals metabolites and metabolic pathways that di | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| zenodo.14046073 | Zenodo | metabolites | Merging metabolomics and genomics provides a catalog of genetic factors that infuence mole | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| zenodo.14046986 | Zenodo | metabolites | Description of metabolic differences between castrated males and intact gilts obtained fro | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| zenodo.14560461 | Zenodo | metabolites | Understanding Metabolic Alterations in Advanced Stage Chronic Kidney Disease Patients by N | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| zenodo.15194775 | Zenodo | metabolites | Effects of Short-Term Exercise on Plasma Metabolic and Lipidomic Profiles of Individuals w | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| zenodo.15461218 | Zenodo | metabolites | NMR based Metabolomics Analysis of fasting serum samples of Type-II Diabetes Mellitus | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| zenodo.16531492 | Zenodo | metabolites | Plasma NMR Metabolomics Dataset Used in the Study: "Obesity-driven metabolic dysregulation | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| zenodo.16903869 | Zenodo | metabolites | International Summer School on Non-Target Metabolomics - Group 2 - Skin Swabs | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| zenodo.17556858 | Zenodo | metabolites | NMR based clinical metabolomics for predicting onset of Gestational diabetes mellitus (GDM | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| zenodo.18921694 | Zenodo | metabolites | Raw data used in the manuscript titled 'NMR-based serum metabolomics in patients with low- | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| zenodo.21787316 | Zenodo | metabolites | Summary Statistics from a Whole-Genome Sequencing Metabolomic GWAS and Whole-Blood cis-eQT | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| zenodo.21840009 | Zenodo | metabolites | Supplementary Data for "Abdominal Massage Improves Gut Barrier Function and Is Associated  | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| zenodo.22228513 | Zenodo | metabolites | Transcriptomic profiles of brain tumors (glioblastoma, GBM) and corresponding orthotopic x | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| zenodo.23788488 | Zenodo | metabolites | Source data and code for Figure 2 of ZT Scan DIA spectral data mining for metabolomics usi | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| zenodo.4009418 | Zenodo | metabolites | Data from: Metabolomics - Emory Cardiovascular Biobank | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| zenodo.4898767 | Zenodo | metabolites | Supplemental tables 1, 2 and 3 to the paper: "Impact of combined consumption of fish oil a | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| zenodo.5363485 | Zenodo | metabolites | Neonatal exposure to BPA, BDE-99, and PCB produces persistent changes in hepatic transcrip | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| zenodo.5513937 | Zenodo | metabolites | GTTAtlas - Metabolomics atlas of oral 13C-glucose tolerance test in mice | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| zenodo.5893528 | Zenodo | metabolites | Moose omics dataframe | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| zenodo.6597902 | Zenodo | metabolites | Normalized NMR integration values from the metabolomic and lipidomic analysis of red blood | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| zenodo.7385508 | Zenodo | metabolites | Urinary volatilome database | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| E-MTAB-12777 | ArrayExpress | proteins | Discovery and systematic assessment of early biomarkers that predict progression to severe | exclude | antibody array; severity class 'critical or severe' merges severe (oxygen) with critical, so it cannot be mapped |
| E-MTAB-14495 | ArrayExpress | proteins | Investigation of IgG-reactivity Profiles and Epitope Landscape in Divergent COVID-19 Traje | exclude | not a plasma or serum proteome of patients (cells, tissue, other fluid or viral proteins) |
| PRJDB18846 | BioProject | proteins | Epigenetic regulation by nitrosative stress as a potential mechanism of long COVID | exclude | other assay or analyte |
| PRJDB18847 | BioProject | proteins | Epigenetic regulation by nitrosative stress as a potential mechanism of long COVID | exclude | other assay or analyte |
| PRJEB42534 | BioProject | proteins | Blood RNA-seq data from individuals during the UK Covid-19 lockdown, spanning sequential w | exclude | not a plasma or serum proteome of patients (cells, tissue, other fluid or viral proteins) |
| PRJEB43418 | BioProject | proteins | HLA Class I-associated expansion of TRBV11-2 T cells via a CDR3-independent mechanism in M | exclude | not adults |
| PRJEB45262 | BioProject | proteins | Impact of high platelet turnover on the platelet transcriptome: Results from platelet RNA- | exclude | not a plasma or serum proteome of patients (cells, tissue, other fluid or viral proteins) |
| PRJEB46564 | BioProject | proteins | We characterized the expression profile of post-mortem human spleen samples by Ion AmpliSe | exclude | not a plasma or serum proteome of patients (cells, tissue, other fluid or viral proteins) |
| PRJEB53053 | BioProject | proteins | Longitudinal scRNAseq profiling of memory B cells in severe COVID19 patients | exclude | single cell data |
| PRJEB71810 | BioProject | proteins | Spike-specific T cell response in vaccinated | exclude | in vitro or cell culture |
| PRJNA1006628 | BioProject | proteins | Human RP105 monoclonal antibody enables the production of antigen-specific antibody | exclude | duplicate: BioProject of GSE241135 |
| PRJNA1009393 | BioProject | proteins | NULISA: a novel proteomic liquid biopsy platform with attomolar sensitivity and high multi | exclude | duplicate: BioProject of GSE241717 |
| PRJNA1010223 | BioProject | proteins | SARS-CoV-2 spike protein induces the cytokine release syndrome by stimulating T cells to p | exclude | duplicate: BioProject of GSE241843 |
| PRJNA1011802 | BioProject | proteins | A Hub Gene Signature as a Therapeutic Target and Biomarker for Sepsis and Geriatric Sepsis | exclude | duplicate: BioProject of GSE242127 |
| PRJNA1026530 | BioProject | proteins | Longitudinal Immune Profiling and Corticosteroid Modulation of the Immediate Innate Immune | exclude | vaccination study |
| PRJNA1036593 | BioProject | proteins | Monocyte-mediated thrombosis linked to circulating tissue factor and immune-paralysis in C | exclude | duplicate: BioProject of GSE247186 |
| PRJNA1045548 | BioProject | proteins | Single-cell characteristics of longitudinal immune profiling in subjects with breakthrough | exclude | duplicate: BioProject of GSE248669 |
| PRJNA1049421 | BioProject | proteins | High-dimensional comparison of monocytes and T cells in post-COVID and idiopathic pulmonar | exclude | duplicate: BioProject of GSE249513 |
| PRJNA1063318 | BioProject | proteins | Recruitment of plasma cells from follicular and extrafollicular immune reactions to the bo | exclude | duplicate: BioProject of GSE252959 |
| PRJNA1063433 | BioProject | proteins | Hybrid immunity to SARS-CoV-2 arises from serological recall of IgG antibodies distinctly  | exclude | vaccination study |
| PRJNA1067578 | BioProject | proteins | Recruitment of plasma cells from follicular and extrafollicular immune reactions to the bo | exclude | duplicate: BioProject of GSE253857 |
| PRJNA1070999 | BioProject | proteins | Inhibition of CXCL10 and IFN-γ ameliorates myocarditis in preclinical models of SARS-CoV-  | exclude | duplicate: BioProject of GSE254537 |
| PRJNA1071588 | BioProject | proteins | Inhibition of PCSK9 with evolocumab modulates lipoproteins and monocyte activation in high | exclude | duplicate: BioProject of GSE254732 |
| PRJNA1084870 | BioProject | proteins | Maintenance and functional regulation of immune memory to COVID-19 vaccines in tissues | exclude | duplicate: BioProject of GSE261278 |
| PRJNA1089238 | BioProject | proteins | Prior vaccination prevents overactivation of innate immune responses during COVID-19 break | exclude | duplicate: BioProject of GSE261862 |
| PRJNA1101439 | BioProject | proteins | A 50-gene high-risk profile predictive of COVID-19 and Idiopathic Pulmonary Fibrosis morta | exclude | duplicate: BioProject of GSE264196 |
| PRJNA1117008 | BioProject | proteins | Two-dimensional high-throughput on-cell screening of immunoglobulins against broad antigen | exclude | not a plasma or serum proteome of patients (cells, tissue, other fluid or viral proteins) |
| PRJNA1134144 | BioProject | proteins | Mapping of human monoclonal antibody responses to XBB.1.5 COVID-19 monovalent vaccines: a  | exclude | vaccination study |
| PRJNA1136554 | BioProject | proteins | Transcriptional response to the SARS-CoV-2 Delta variant in patients with severe COVID-19 | exclude | duplicate: BioProject of GSE272392 |
| PRJNA1140207 | BioProject | proteins | Dasatinib Overcomes Acquired Resistance to BRAF/MEK Inhibitor Therapy | exclude | duplicate: BioProject of GSE273120 |
| PRJNA1149196 | BioProject | proteins | Cardiovascular symptoms of PASC are associated with trace-level cytokines that affect the  | exclude | in vitro or cell culture |
| PRJNA1155161 | BioProject | proteins | Unveiling Alternate Pathways for SARS-CoV-2 Infection via Extracellular Vesicle-Mediated T | exclude | duplicate: BioProject of GSE276140 |
| PRJNA1157815 | BioProject | proteins | MVA-based SARS-CoV-2 vaccine candidates encoding different spike protein conformations ind | exclude | duplicate: BioProject of GSE276544 |
| PRJNA1162484 | BioProject | proteins | SARS-CoV-2 nsp15 enhances viral virulence by subverting host antiviral defenses | exclude | in vitro or cell culture |
| PRJNA1182358 | BioProject | proteins | The identification of a SARs-CoV2 S1 protein derived peptide with super-antigen-like stimu | exclude | duplicate: BioProject of GSE281135 |
| PRJNA1190999 | BioProject | proteins | Targeting G9a Regulated m6A Modification Pathway to Hinder SARS-CoV-2 Replication | exclude | in vitro or cell culture |
| PRJNA1191027 | BioProject | proteins | The vaccine platform used for COVID-19 primary immunization shapes the quality of the huma | exclude | single cell data |
| PRJNA1194117 | BioProject | proteins | Plasma proteomics and organ dysfunction in COVID-19 ARDS | exclude | duplicate: BioProject of GSE283501 |
| PRJNA1205767 | BioProject | proteins | Early innate immune response and evolution of a SARS-CoV-2 furin cleavage site inactive va | exclude | duplicate: BioProject of GSE285756 |
| PRJNA1235184 | BioProject | proteins | Macrophages direct location-dependent recall of B cell memory to vaccination [scRNA-seq an | exclude | duplicate: BioProject of GSE291745 |
| PRJNA1251309 | BioProject | proteins | Homo sapiens | exclude | not a plasma or serum proteome of patients (cells, tissue, other fluid or viral proteins) |
| PRJNA1268552 | BioProject | proteins | Metabolic quiescence of naïve-like memory T cells precedes and maintains antigen-specific  | exclude | duplicate: BioProject of GSE298237 |
| PRJNA1279089 | BioProject | proteins | Whole blood gene expression data of thrombotic and non-thrombotic critically ill COVID-19  | exclude | duplicate: BioProject of GSE300129 |
| PRJNA1291127 | BioProject | proteins | Single-cell multi-omics reveals dyssynchrony of the innate and adaptive immune system in p | exclude | single cell data |
| PRJNA1298582 | BioProject | proteins | Dosing interval is a major factor determining the quality of T cells induced by SARS-CoV-2 | exclude | duplicate: BioProject of GSE303984 |
| PRJNA1299289 | BioProject | proteins | Glycosylphosphatidylinositol biosynthesis functions as a conserved host defense pathway ag | exclude | not a plasma or serum proteome of patients (cells, tissue, other fluid or viral proteins) |
| PRJNA1303629 | BioProject | proteins | Synthetic mRNA Vaccines and Transcriptomic Dysregulation: Evidence from New-Onset Adverse  | exclude | duplicate: BioProject of GSE304973 |
| PRJNA1331709 | BioProject | proteins | Pulmonary fibrosis after COVID-19 is characterized by airway abnormalities and elevated cl | exclude | duplicate: BioProject of GSE308611 |
| PRJNA1405735 | BioProject | proteins | Extracellular vesicle protein and miRNA signatures as biomarkers for post-infectious ME/CF | exclude | other assay or analyte |
| PRJNA1416848 | BioProject | proteins | Integrative Epigenomic and Transcriptomic Analysis of Peripheral Blood Across Different SA | exclude | duplicate: BioProject of GSE318176 |
| PRJNA1420870 | BioProject | proteins | A comparative multiomic atlas of human immunity to multiple COVID-19 vaccines | exclude | duplicate: BioProject of GSE318946 |
| PRJNA1428820 | BioProject | proteins | A comparative multiomic atlas of human immunity to multiple COVID-19 vaccines | exclude | duplicate: BioProject of GSE320554 |
| PRJNA1508924 | BioProject | proteins | Isoform-level analysis reveals conserved early changes in transcript usage during human va | exclude | duplicate: BioProject of GSE342810 |
| PRJNA1511296 | BioProject | proteins | Narsoplimab (anti-MASP2) modulates long COVID patient plasma-induced transcriptional respo | exclude | duplicate: BioProject of GSE343379 |
| PRJNA1551418 | BioProject | proteins | Early Fc-effector antibody signatures impact COVID-19 disease trajectory | exclude | not a plasma or serum proteome of patients (cells, tissue, other fluid or viral proteins) |
| PRJNA1551419 | BioProject | proteins | Early Fc-effector antibody signatures impact COVID-19 disease trajectory | exclude | not a plasma or serum proteome of patients (cells, tissue, other fluid or viral proteins) |
| PRJNA277352 | BioProject | proteins | Defective structural RNA processing in relapsing-remitting multiple sclerosis | exclude | duplicate: BioProject of GSE66573 |
| PRJNA521077 | BioProject | proteins | RNA-seq and ribosome footprint profiling of human and mouse platelets with sepsis | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| PRJNA574628 | BioProject | proteins | VHL shapes the pre-immune B cell-repertoire by controlling HIF-1α | exclude | not human |
| PRJNA574906 | BioProject | proteins | VHL shapes the pre-immune B cell-repertoire by controlling HIF-1α (run2) | exclude | not human |
| PRJNA574931 | BioProject | proteins | VHL shapes the pre-immune B cell-repertoire by controlling HIF-1α (run3) | exclude | not human |
| PRJNA579053 | BioProject | proteins | Identification of differentially expressed lncRNAs and mRNAs in luminal-B breast cancer by | exclude | in vitro or cell culture |
| PRJNA615032 | BioProject | proteins | Transcriptional response to SARS-CoV-2 infection | exclude | duplicate: BioProject of GSE147507 |
| PRJNA634489 | BioProject | proteins | RNA-seq of platelets from SARS-CoV-2 Covid-19 | exclude | not a plasma or serum proteome of patients (cells, tissue, other fluid or viral proteins) |
| PRJNA646224 | BioProject | proteins | Transcriptional and proteomic insights into the host response in fatal COVID-19 cases | exclude | not blood |
| PRJNA648120 | BioProject | proteins | Transcriptional Profiling of Leukocytes in Critically Ill COVID19 Patients | exclude | duplicate: BioProject of GSE154998 |
| PRJNA648957 | BioProject | proteins | Multi-omic data integration allows baseline immune signatures to predict hepatitis B vacci | exclude | duplicate: BioProject of GSE155198 |
| PRJNA648991 | BioProject | proteins | Single-Cell Omics Reveals Dyssynchrony of the Innate and Adaptive Immune System in Progres | exclude | duplicate: BioProject of GSE155223 |
| PRJNA648992 | BioProject | proteins | Single-Cell Omics Reveals Dyssynchrony of the Innate and Adaptive Immune System in Progres | exclude | single cell data |
| PRJNA660067 | BioProject | proteins | Large-scale Multi-omic Analysis of COVID-19 Severity | exclude | duplicate: BioProject of GSE157103 |
| PRJNA662785 | BioProject | proteins | Dexamethasone modulates immature neutrophils and interferon programming in severe COVID-19 | exclude | duplicate: BioProject of GSE157789 |
| PRJNA669691 | BioProject | proteins | Human Tubular Epithelial Cells Activate a Coordinated Stress Response after Serum Exposure | exclude | duplicate: BioProject of GSE159554 |
| PRJNA669838 | BioProject | proteins | Human Tubular Epithelial Cells Activate a Coordinated Stress Response after Serum Exposure | exclude | duplicate: BioProject of GSE159586 |
| PRJNA674755 | BioProject | proteins | Age-determined expression of priming protease TMPRSS2 and localization of SARS-CoV-2 in lu | exclude | not human (Mus musculus) |
| PRJNA675120 | BioProject | proteins | Multi-omic data integration allows baseline immune signatures to predict hepatitis B vacci | exclude | vaccination study |
| PRJNA680293 | BioProject | proteins | Time-resolved Systems Immunology Reveals a Late Juncture Linked to Fatal COVID-19 | exclude | duplicate: BioProject of GSE161918 |
| PRJNA694174 | BioProject | proteins | Transcriptional analysis of SARS-CoV-2 specific CD4 T cells | exclude | duplicate: BioProject of GSE165373 |
| PRJNA701047 | BioProject | proteins | Immune dysregulation and autoreactivity correlate with disease severity in SARS-CoV-2-asso | exclude | duplicate: BioProject of GSE166489 |
| PRJNA701725 | BioProject | proteins | Clonal expansion and activation of tissue-resident memory-like Th17 cells expressing GM-CS | exclude | duplicate: BioProject of GSE167118 |
| PRJNA707445 | BioProject | proteins | Longitudinal single-cell epitope and RNA-sequencing reveals the immunological impact of ty | exclude | duplicate: BioProject of GSE168453 |
| PRJNA717250 | BioProject | proteins | Prevalent, protective, and convergent IgG recognition of SARS-CoV-2 non-RBD spike epitopes | exclude | not a plasma or serum proteome of patients (cells, tissue, other fluid or viral proteins) |
| PRJNA717310 | BioProject | proteins | Divergent and self-reactive immune responses in the CNS of COVID-19 patients with neurolog | exclude | single cell data |
| PRJNA722046 | BioProject | proteins | Identification of driver genes for severe forms of COVID-19 in a deeply phenotyped young p | exclude | duplicate: BioProject of GSE172114 |
| PRJNA724033 | BioProject | proteins | Nopp140-chaperoned 2’-O-methylation of small nuclear RNAs in Cajal bodies ensures splicing | exclude | single cell data |
| PRJNA725183 | BioProject | proteins | Blood transcriptomes of anti-SARS-CoV2 antibody  positive healthy individuals with prior   | exclude | duplicate: BioProject of GSE173317 |
| PRJNA726234 | BioProject | proteins | Longitudinal trajectories of mild/moderate SARS-CoV2 infection reveals key precursors and  | exclude | duplicate: BioProject of GSE173590 |
| PRJNA728070 | BioProject | proteins | Increased Systemic HSP70B Levels in Spinal Muscular Atrophy Infants | exclude | not adults |
| PRJNA728100 | BioProject | proteins | Multi-omic profiling reveals widespread dysregulation of innate immunity and hematopoiesis | exclude | duplicate: BioProject of GSE174072 |
| PRJNA728117 | BioProject | proteins | Large-scale genomic study reveals robust activation of the immune system following advance | exclude | duplicate: BioProject of GSE174083 |
| PRJNA729513 | BioProject | proteins | Prevalent, protective, and convergent IgG recognition of SARS-CoV-2 non-RBD spike epitopes | exclude | not a plasma or serum proteome of patients (cells, tissue, other fluid or viral proteins) |
| PRJNA730810 | BioProject | proteins | BNT162b2 vaccination enhances interferon-JAK-STAT-regulated antiviral programs in COVID-19 | exclude | duplicate: BioProject of GSE174621 |
| PRJNA733168 | BioProject | proteins | Bidirectional genome-wide CRISPR screens reveal host factors regulating SARS-CoV-2,  MERS- | exclude | in vitro or cell culture |
| PRJNA737991 | BioProject | proteins | Whole transcriptome miRNA profiling of Coronavirus Disease 2019 (COVID-19) patient plasma. | exclude | other assay or analyte |
| PRJNA738502 | BioProject | proteins | mRNA profiling in pooled human umbilical vein endothelial cells stimulated with Coronaviru | exclude | duplicate: BioProject of GSE178331 |
| PRJNA738876 | BioProject | proteins | Th2/Th1 cytokine imbalance is associated with higher COVID-19 risk mortality | exclude | duplicate: BioProject of GSE178399 |
| PRJNA743882 | BioProject | proteins | Homo sapiens Transcriptome or Gene expression | exclude | single cell data |
| PRJNA747055 | BioProject | proteins | Pulmonary epithelial and endothelial cell damage, loss of repair, inhibition of fibrinolys | exclude | duplicate: BioProject of GSE180226 |
| PRJNA749310 | BioProject | proteins | New-Onset IgG Autoantibodies in Hospitalized Patients with COVID-19 [Unpaired Cytokine and | exclude | duplicate: BioProject of GSE180741 |
| PRJNA749312 | BioProject | proteins | New-Onset IgG Autoantibodies in Hospitalized Patients with COVID-19 [Paired Cytokine] | exclude | duplicate: BioProject of GSE180740 |
| PRJNA749313 | BioProject | proteins | New-Onset IgG Autoantibodies in Hospitalized Patients with COVID-19 [Autoantigen V2] | exclude | duplicate: BioProject of GSE180738 |
| PRJNA749314 | BioProject | proteins | New-Onset IgG Autoantibodies in Hospitalized Patients with COVID-19 [Unpaired Autoantigen  | exclude | duplicate: BioProject of GSE180739 |
| PRJNA749315 | BioProject | proteins | New-Onset IgG Autoantibodies in Hospitalized Patients with COVID-19 [Autoantigen V1] | exclude | duplicate: BioProject of GSE180737 |
| PRJNA754796 | BioProject | proteins | Clinical significance of interferon signal based on mRNA-microRNA integration and plasma p | exclude | duplicate: BioProject of GSE182152 |
| PRJNA758448 | BioProject | proteins | Longitudinal immune dynamics of mild COVID-19 define signatures of recovery and persistenc | exclude | single cell data |
| PRJNA758449 | BioProject | proteins | Longitudinal immune dynamics of mild COVID-19 define signatures of recovery and persistenc | exclude | single cell data |
| PRJNA765050 | BioProject | proteins | Enteric coronavirus infection and treatment modeled with an immunocompetent human intestin | exclude | duplicate: BioProject of GSE184569 |
| PRJNA771040 | BioProject | proteins | Transcriptomic profiles revealed down regulation of low-density lipoprotein particle recep | exclude | duplicate: BioProject of GSE185863 |
| PRJNA781246 | BioProject | proteins | Prior BNT162b vaccination enhances interferon-JAK-STAT regulated antiviral programs in hos | exclude | duplicate: BioProject of GSE189039 |
| PRJNA789649 | BioProject | proteins | Comprehensive comparison of adaptive immune response by inactivated SARS-CoV-2 vaccine bet | exclude | duplicate: BioProject of GSE191089 |
| PRJNA795179 | BioProject | proteins | SARS-CoV-2 diverges from other betacoronaviruses in only partially activating the IRE1a/XB | exclude | duplicate: BioProject of GSE193169 |
| PRJNA798057 | BioProject | proteins | Bidirectional genome-wide CRISPR screens reveal host factors regulating SARS-CoV-2, MERS-C | exclude | in vitro or cell culture |
| PRJNA800394 | BioProject | proteins | Influenza vaccination reveals sex dimorphic imprints of prior mild COVID-19 [RNA-Seq] | exclude | duplicate: BioProject of GSE194378 |
| PRJNA803375 | BioProject | proteins | SARS-CoV-2 disrupts respiratory vascular barriers by suppressing Claudin-5 expression | exclude | duplicate: BioProject of GSE196113 |
| PRJNA804953 | BioProject | proteins | Viral mimicry of Interleukin-17A by SARS-COV-2 ORF8 | exclude | duplicate: BioProject of GSE196455 |
| PRJNA807364 | BioProject | proteins | A public antibody class recognizes an S2 epitope exposed on open conformations of SARS-CoV | exclude | duplicate: BioProject of GSE196820 |
| PRJNA814386 | BioProject | proteins | Proliferative Receptor Binding Domain (RBD)-specific CD4+ T cells correlate with neutraliz | exclude | duplicate: BioProject of GSE198281 |
| PRJNA828379 | BioProject | proteins | Development of Peripheral Spondyloarthritis after COVID-19 in HLA-B27-Positive Monozygotic | exclude | single cell data |
| PRJNA838497 | BioProject | proteins | Single-Cell RNA Sequencing to Decipher the Immunogenicity of ChAdOx1 nCoV-19/AZD1222 and m | exclude | duplicate: BioProject of GSE203081 |
| PRJNA848805 | BioProject | proteins | Transcriptomes analysis of COVID-19 Inactivated Vaccine and Protein Subunit Vaccine as Thi | exclude | duplicate: BioProject of GSE206023 |
| PRJNA849284 | BioProject | proteins | Augmenting SARS-CoV-2 vaccine induced spike protein-specific CD4+ T cell responses through | exclude | duplicate: BioProject of GSE206134 |
| PRJNA849921 | BioProject | proteins | Differential chromatin accessibility in peripheral blood mononuclear cells underlies COVID | exclude | duplicate: BioProject of GSE206263 |
| PRJNA849922 | BioProject | proteins | Influenza vaccination reveals sex dimorphic imprints of prior mild COVID-19 [single-cell C | exclude | duplicate: BioProject of GSE206265 |
| PRJNA849934 | BioProject | proteins | Differential chromatin accessibility in peripheral blood mononuclear cells underlies COVID | exclude | single cell data |
| PRJNA849954 | BioProject | proteins | Differential chromatin accessibility in peripheral blood mononuclear cells underlies COVID | exclude | duplicate: BioProject of GSE206283 |
| PRJNA851014 | BioProject | proteins | Differential chromatin accessibility in peripheral blood mononuclear cells underlies COVID | exclude | single cell data |
| PRJNA851623 | BioProject | proteins | SARS-CoV-2 disrupt respiratory vascular barrier by suppressing Claudin-5 expression (addit | exclude | duplicate: BioProject of GSE206635 |
| PRJNA853326 | BioProject | proteins | COVID-19 drives a distinct plasma proteome in pregnant and non-pregnant individuals | exclude | duplicate: BioProject of GSE207015 |
| PRJNA864802 | BioProject | proteins | CD8+ T-cell memory induced by successive SARS-CoV-2 mRNA vaccinations is characterized by  | exclude | single cell data |
| PRJNA864803 | BioProject | proteins | CD8+ T-cell memory induced by successive SARS-CoV-2 mRNA vaccinations is characterized by  | exclude | single cell data |
| PRJNA864804 | BioProject | proteins | CD8+ T-cell memory induced by successive SARS-CoV-2 mRNA vaccinations is characterized by  | exclude | duplicate: BioProject of GSE210229 |
| PRJNA873691 | BioProject | proteins | Longitudinal characterization of circulating neutrophils uncovers phenotypes associated wi | exclude | duplicate: BioProject of GSE212041 |
| PRJNA877667 | BioProject | proteins | Dynamics of gene expression profiling by RNA-Seq and identification of high-risk patients  | exclude | duplicate: BioProject of GSE212861 |
| PRJNA877669 | BioProject | proteins | Dynamics of gene expression profiling by microarrays and identification of high-risk patie | exclude | duplicate: BioProject of GSE212865 |
| PRJNA890633 | BioProject | proteins | Adaptive immune responses to SARS-CoV-2 persist in the pharyngeal lymphoid tissue of child | exclude | duplicate: BioProject of GSE215802 |
| PRJNA908773 | BioProject | proteins | Associations of miRNAs with blood phenotypes and ICU admission in COVID-19 patients [mRNA- | exclude | duplicate: BioProject of GSE220076 |
| PRJNA908778 | BioProject | proteins | Associations of miRNAs with blood phenotypes and ICU admission in COVID-19 patients [miRNA | exclude | duplicate: BioProject of GSE220077 |
| PRJNA910934 | BioProject | proteins | Activation of Coagulation and Proinflammatory Pathways in Thrombosis with Thrombocytopenia | exclude | duplicate: BioProject of GSE220659 |
| PRJNA912320 | BioProject | proteins | COVID-19 severity microbiome | exclude | not a plasma or serum proteome of patients (cells, tissue, other fluid or viral proteins) |
| PRJNA923335 | BioProject | proteins | Autoantibodies are highly prevalent in non-SARS-CoV-2 respiratory infections and critical  | exclude | duplicate: BioProject of GSE222760 |
| PRJNA923336 | BioProject | proteins | Autoantibodies are highly prevalent in non-SARS-CoV-2 respiratory infections and critical  | exclude | duplicate: BioProject of GSE222764 |
| PRJNA923688 | BioProject | proteins | Primary ChAdOx1 vaccination does not reactivate pre-existing, cross-reactive immunity | exclude | duplicate: BioProject of GSE222872 |
| PRJNA925441 | BioProject | proteins | Cytotoxic T cells targeting spike glycoprotein are associated with hybrid immunity to SARS | exclude | duplicate: BioProject of GSE223236 |
| PRJNA931952 | BioProject | proteins | Long COVID manifests with T cell dysregulation, inflammation, and an uncoordinated adaptiv | exclude | duplicate: BioProject of GSE224615 |
| PRJNA935327 | BioProject | proteins | Proteomic Characterization of Acute Kidney Injury in Patients Hospitalized with SARS-CoV2  | exclude | duplicate: BioProject of GSE225349 |
| PRJNA938055 | BioProject | proteins | Neuro-PASC is characterized by enhanced CD4+ and diminished CD8+ T cell responses to SARS- | exclude | not a plasma or serum proteome of patients (cells, tissue, other fluid or viral proteins) |
| PRJNA939253 | BioProject | proteins | Long COVID involves activation of proinflammatory and immune exhaustion pathways | exclude | duplicate: BioProject of GSE226260 |
| PRJNA941336 | BioProject | proteins | TGF-β1 alter the SARS-CoV2 pathogenesis and modulates ACE2 expression by miRNA dependent m | exclude | in vitro or cell culture |
| PRJNA971080 | BioProject | proteins | "Secreted ORF8 is a pathogenic cause of severe COVID-19 and is potentially targetable with | exclude | in vitro or cell culture |
| PRJNA971962 | BioProject | proteins | CD8+ T-cell memory induced by successive SARS-CoV-2 mRNA vaccinations is characterized by  | exclude | single cell data |
| PRJNA976394 | BioProject | proteins | Microfluidic immuno-serolomic assay reveals systems level association with COVID-19 pathol | exclude | duplicate: BioProject of GSE233506 |
| PRJNA979185 | BioProject | proteins | Transcriptomes analysis of COVID-19 Inactivated Vaccine and Protein Subunit Vaccine as Thi | exclude | duplicate: BioProject of GSE234008 |
| PRJNA983740 | BioProject | proteins | The Transcriptional and Phenotypic Characteristics that Define Alveolar Macrophage Subsets | exclude | duplicate: BioProject of GSE234918 |
| PRJNA984257 | BioProject | proteins | scRNA-seq characterization of immune cell types among individuals with or without Long COV | exclude | duplicate: BioProject of GSE235050 |
| PRJNA985596 | BioProject | proteins | Transcriptome analysis of transplanted vascular organoid grafts under cranial window at th | exclude | duplicate: BioProject of GSE235232 |
| GSE100054 | GEO | proteins | Expression profiling of peripheral blood mononuclear cells in normal controls and patients | exclude | 19 samples: cannot reach 10 participants per class |
| GSE100866 | GEO | proteins | CITE-seq: Large scale simultaneous measuremnt of epitopes and transcriptomes in single cel | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE103412 | GEO | proteins | Mucositis in patients with tonsil squamous cell carcinoma: gene expression profile of oral | exclude | no characteristic splits 10 COVID samples from 10 others and no severity field (automatic; best : 0 against 0) |
| GSE104214 | GEO | proteins | BPIFA1 Regulates Interferon-Mediated Lung Neutrophil Recruitment During Acute Inflammation | exclude | not human (Mus musculus) |
| GSE110496 | GEO | proteins | Single cell transcriptional dynamics of flavivirus infection | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE114297 | GEO | proteins | Pseudotime Ordering of Single Human Beta-Cells Reveals States of Insulin Production and Un | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE124263 | GEO | proteins | Neonatal and adult human testis defined at the single-cell level | exclude | single cell data |
| GSE124334 | GEO | proteins | Single-Cell Transcriptome Analysis of CD34+ Stem Cell-Derived Myeloid Cells Infected With  | exclude | single cell data |
| GSE12654 | GEO | proteins | Gene expression from human prefrontal cortex (BA10) | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE128639 | GEO | proteins | Comprehensive integration of single-cell data | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE130117 | GEO | proteins | Simultaneous measurement of DNA repair activities and mRNA abundance in single cells | exclude | single cell data |
| GSE13044 | GEO | proteins | Gene expression profiling in the lung and liver of low and high dose Perfluorooctanoic Aci | exclude | not human (Mus musculus) |
| GSE134347 | GEO | proteins | Protein-coding and non-coding RNA expression in whole blood leukocytes of healthy subjects | exclude | no characteristic splits 10 COVID samples from 10 others and no severity field (automatic; best : 0 against 0) |
| GSE135710 | GEO | proteins | Single-cell RNA-seq of B cell polustion in breast cancer patients | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE138064 | GEO | proteins | Interferon-β corrects massive gene dysregulation in multiple sclerosis:  Short-term and lo | exclude | no characteristic splits 10 COVID samples from 10 others and no severity field (automatic; best : 0 against 0) |
| GSE139369 | GEO | proteins | Single-cell, multi-omic analysis identifies regulatory programs in mixed phenotype acute l | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE140231 | GEO | proteins | A human single cell atlas of the substantia nigra reveals novel cell specific pathways ass | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE141578 | GEO | proteins | CD4+ T cells contribute to neurodegeneration in Lewy body dementia [CSF_and_PBMCs_Healthy  | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE145443 | GEO | proteins | An atlas of cell types in the mammalian epididymis and vas deferens [single-cell RNA-Seq] | exclude | not human (Mus musculus) |
| GSE14630 | GEO | proteins | Effect of Mycophenolic Acid on renal transplant recipients | exclude | 10 samples: cannot reach 10 participants per class |
| GSE149973 | GEO | proteins | Decoding SARS-CoV-2 coding capacity | exclude | no characteristic splits 10 COVID samples from 10 others and no severity field (automatic; best : 0 against 0) |
| GSE150819 | GEO | proteins | Generation of human bronchial organoids for SARS-CoV-2 research | exclude | in vitro or cell culture |
| GSE15129 | GEO | proteins | Coenzyme Q10-dependent gene expression in SAMP1 mice tissues | exclude | not human (Mus musculus) |
| GSE151973 | GEO | proteins | A gateway for SARS-CoV-2 infection in the human olfactory neuroepithelium | exclude | single cell data |
| GSE153684 | GEO | proteins | SARS-CoV-2 infection and replication in human fetal and pediatric gastric organoids | exclude | in vitro or cell culture |
| GSE153712 | GEO | proteins | Meta-analysis of genome-wide DNA methylation identifies shared associations across neurode | exclude | single cell data |
| GSE15490 | GEO | proteins | Sequential gene expression profiling in CLL during treatment | exclude | no characteristic splits 10 COVID samples from 10 others and no severity field (automatic; best : 0 against 0) |
| GSE156473 | GEO | proteins | High throughput joint profiling of chromatin accessibility and protein levels in single ce | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE156701 | GEO | proteins | Vascular Disease and Thrombosis in SARS-CoV-2 Infected Humans and Rhesus Macaques | exclude | not human (Macaca mulatta) |
| GSE159787 | GEO | proteins | Temporal and Spatial Heterogeneity of Host Response to SARS-CoV-2 Pulmonary Infection [gen | exclude | not a plasma or serum proteome of patients (cells, tissue, other fluid or viral proteins) |
| GSE162208 | GEO | proteins | Innate and adaptive immune responses to SARS-CoV-2 in Syrian hamsters disclosed by single- | exclude | not human (Mesocricetus auratus) |
| GSE162629 | GEO | proteins | A novel cell culture system modeling the SARS-CoV-2 life cycle | exclude | in vitro or cell culture |
| GSE164593 | GEO | proteins | The role of XIST in X-inactivation maintenance in B cells (CRISPRi screen) | exclude | single cell data |
| GSE164594 | GEO | proteins | The role of XIST in X-inactivation maintenance in B cells (ChIP-seq, Cut&Run)) | exclude | single cell data |
| GSE166990 | GEO | proteins | Modeling SARS-CoV-2 infection and its individual differences with ACE2-expressing human iP | exclude | in vitro or cell culture |
| GSE17078 | GEO | proteins | Cell Adhesion Molecule 1 (CADM1): A Novel Risk Factor for Venous Thrombosis | exclude | no characteristic splits 10 COVID samples from 10 others and no severity field (automatic; best : 0 against 0) |
| GSE171370 | GEO | proteins | Genome-wide Analyses Revealed the Global Impacts of SARS-CoV-2 Viral Gene Orf9c on Human C | exclude | not blood |
| GSE173507 | GEO | proteins | Discovery and functional interrogation of the virus and host RNA interactome of SARS-CoV-2 | exclude | in vitro or cell culture |
| GSE17400 | GEO | proteins | Dynamic Innate Immune Responses of Human Bronchial Epithelial Cells against SARS-CoV and D | exclude | in vitro or cell culture |
| GSE178399 | GEO | proteins | Th2/Th1 cytokine imbalance is associated with higher COVID-19 risk mortality | exclude | COVID-19 patients only; outcome (died or survived) has no mild or moderate class |
| GSE180737 | GEO | proteins | New-Onset IgG Autoantibodies in Hospitalized Patients with COVID-19 [Autoantigen V1] | exclude | not a plasma or serum proteome of patients (cells, tissue, other fluid or viral proteins) |
| GSE180738 | GEO | proteins | New-Onset IgG Autoantibodies in Hospitalized Patients with COVID-19 [Autoantigen V2] | exclude | not a plasma or serum proteome of patients (cells, tissue, other fluid or viral proteins) |
| GSE180739 | GEO | proteins | New-Onset IgG Autoantibodies in Hospitalized Patients with COVID-19 [Unpaired Autoantigen  | exclude | not a plasma or serum proteome of patients (cells, tissue, other fluid or viral proteins) |
| GSE180740 | GEO | proteins | New-Onset IgG Autoantibodies in Hospitalized Patients with COVID-19 [Paired Cytokine] | exclude | not a plasma or serum proteome of patients (cells, tissue, other fluid or viral proteins) |
| GSE180741 | GEO | proteins | New-Onset IgG Autoantibodies in Hospitalized Patients with COVID-19 [Unpaired Cytokine and | exclude | not a plasma or serum proteome of patients (cells, tissue, other fluid or viral proteins) |
| GSE180742 | GEO | proteins | New-Onset IgG Autoantibodies in Hospitalized Patients with COVID-19 [viral] | exclude | not a plasma or serum proteome of patients (cells, tissue, other fluid or viral proteins) |
| GSE182152 | GEO | proteins | Clinical significance of interferon signal based on mRNA-microRNA integration and plasma p | exclude | other assay or analyte |
| GSE186650 | GEO | proteins | A virus-specific monocyte inflammatory phenotype is induced by SARS-CoV2 at the immune-epi | exclude | no characteristic splits 10 COVID samples from 10 others and no severity field (automatic; best : 0 against 0) |
| GSE187420 | GEO | proteins | Deciphering COVID-19 host transcriptomic complexity and variations for therapeutic discove | exclude | 9 samples: cannot reach 10 participants per class |
| GSE189706 | GEO | proteins | RNA sequencing of SARS-CoV and SARS-CoV-2 N proteins expressing cell lines | exclude | in vitro or cell culture |
| GSE190193 | GEO | proteins | Next Generation Sequencing of SARS-CoV-2 infected lung epithelial cells derived from iPSC | exclude | in vitro or cell culture |
| GSE190850 | GEO | proteins | Aerosol delivery, but not intramuscular injection, of adenovirus-vectored tuberculosis vac | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE19392 | GEO | proteins | Dynamic responses of primary human bronchial epithelial cells to influenza virus, viral RN | exclude | single cell data |
| GSE194122 | GEO | proteins | A sandbox for prediction and integration of DNA, RNA, and proteins in single cells | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE19587 | GEO | proteins | Imaging-guided microarray: Identifies molecular markers in the pathogenesis of Parkinson’s | exclude | not human |
| GSE207945 | GEO | proteins | Evaluating the utility of proteomics for the identification of circulating pharmacodynamic | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE20940 | GEO | proteins | Lactobacillus rhamnosus LGG and LC705 effects on human primary macrophages | exclude | 18 samples: cannot reach 10 participants per class |
| GSE20948 | GEO | proteins | The Effect of Hepatitis C Virus Infection on Host Gene Expression | exclude | no characteristic splits 10 COVID samples from 10 others and no severity field (automatic; best : 0 against 0) |
| GSE20966 | GEO | proteins | Gene expression profiles of beta-cell enriched tissue obtained by Laser Capture Microdisse | exclude | no characteristic splits 10 COVID samples from 10 others and no severity field (automatic; best : 0 against 0) |
| GSE21618 | GEO | proteins | Expression profiling of ligand-stimulated wild type and tamoxifen-resistant MCF-7 | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE220077 | GEO | proteins | Associations of miRNAs with blood phenotypes and ICU admission in COVID-19 patients [miRNA | exclude | other assay or analyte |
| GSE222760 | GEO | proteins | Autoantibodies are highly prevalent in non-SARS-CoV-2 respiratory infections and critical  | exclude | not a plasma or serum proteome of patients (cells, tissue, other fluid or viral proteins) |
| GSE222764 | GEO | proteins | Autoantibodies are highly prevalent in non-SARS-CoV-2 respiratory infections and critical  | exclude | not a plasma or serum proteome of patients (cells, tissue, other fluid or viral proteins) |
| GSE225349 | GEO | proteins | Proteomic Characterization of Acute Kidney Injury in Patients Hospitalized with SARS-CoV2  | exclude | COVID-negative participants are hospitalised patients (other disease) and no severity field: no contrast |
| GSE230566 | GEO | proteins | Chronic inflammatory activation of arterial endothelial cells in mice under normal lipopro | exclude | not human (Mus musculus) |
| GSE231910 | GEO | proteins | Respiratory SARS-CoV-2 Infection Induces Skeletal Muscle Atrophy and Long-Lasting Energy M | exclude | not human (Mesocricetus auratus) |
| GSE2431 | GEO | proteins | microarray profiling of potential calorie restriction mimetics | exclude | not human (Mus musculus) |
| GSE243633 | GEO | proteins | Runt-related Transcription Factor 1 (RUNX1) is a mediator of acute kidney injury | exclude | not human (Mus musculus) |
| GSE26484 | GEO | proteins | Sema4A, a novel serum marker of multiple sclerosis, implicates Th17 pathology and efficacy | exclude | in vitro or cell culture |
| GSE27131 | GEO | proteins | Peripheral blood cells expression data form 7 patients with severe pdm(H1N1) influensa and | exclude | no characteristic splits 10 COVID samples from 10 others and no severity field (automatic; best : 0 against 0) |
| GSE283501 | GEO | proteins | Plasma proteomics and organ dysfunction in COVID-19 ARDS | exclude | COVID-19 ARDS patients only: no contrast |
| GSE300488 | GEO | proteins | Characterizing circulating pharmacodynamic biomarkers of IFNβ-1a biologics [inter] | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE300489 | GEO | proteins | Characterizing circulating pharmacodynamic biomarkers of IFNβ-1a biologics [low] | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE30351 | GEO | proteins | Gene expression analysis of the antiviral effect of chloroquine and derivatives on in vitr | exclude | in vitro or cell culture |
| GSE30903 | GEO | proteins | Expression data from untreated and ATP treated AML blasts | exclude | not blood |
| GSE3248 | GEO | proteins | Contribution of Nuclear and Extranuclear PolyQ to Neurological Phenotypes in Mouse Models  | exclude | not human (Mus musculus) |
| GSE325178 | GEO | proteins | Time-Resolved Multi-Omics Identify Biomarkers of Immediate Reactions to mRNA Vaccination [ | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE339049 | GEO | proteins | Severe long-COVID is driven by monocyte reprogramming and broad immune dysregulation linke | exclude | single cell data |
| GSE35765 | GEO | proteins | Effect of fluoxetine treatment on translational profiles of S100a10 cortical pyramidal cel | exclude | not human (Mus musculus) |
| GSE3606 | GEO | proteins | Gene expression profiles in white blood cells in response to exercise | exclude | no characteristic splits 10 COVID samples from 10 others and no severity field (automatic; best : 0 against 0) |
| GSE3737 | GEO | proteins | Omega-6 fatty acids, arachidonic acid (AA) activates PI3K signaling and induces gene expre | exclude | 8 samples: cannot reach 10 participants per class |
| GSE3920 | GEO | proteins | EC_interferon | exclude | not human |
| GSE42026 | GEO | proteins | Transcriptomic profiling in childhood H1N1/09 influenza reveals reduced expression of prot | exclude | not adults |
| GSE48424 | GEO | proteins | Severe preeclampsia: new blood circulating biomarkers | exclude | not adults |
| GSE48575 | GEO | proteins | Pandemic H1N1 Influenza A Directly Induces a Robust and Acute Inflammatory Gene Signature  | exclude | in vitro or cell culture |
| GSE51808 | GEO | proteins | Systems biological analysis of immunity to dengue | exclude | not human |
| GSE59185 | GEO | proteins | Host responses contributing to the attenuation of severe acute respiratory syndrome corona | exclude | not human (Mus musculus) |
| GSE63042 | GEO | proteins | An integrated transcriptome and expressed variant analysis of sepsis survival and death | exclude | not a plasma or serum proteome of patients (cells, tissue, other fluid or viral proteins) |
| GSE66573 | GEO | proteins | Defective structural RNA processing in relapsing-remitting multiple sclerosis | exclude | single cell data |
| GSE68004 | GEO | proteins | Whole blood Transcriptional Profiles as a Prognostic and Diagnostic Tool in Complete and I | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE6930 | GEO | proteins | Expression data from A673 cell line treated with compounds at several doses and time point | exclude | in vitro or cell culture |
| GSE72377 | GEO | proteins | Response of peripheral blood mononuclear cells from  20 healthy donor and 15 patients with | exclude | in vitro or cell culture |
| GSE76925 | GEO | proteins | Functional interactors of genomewide association study genes are differentially expressed  | exclude | single cell data |
| GSE79209 | GEO | proteins | Airway epithelial cells from smokers with and without bronchial premalignant lesions | exclude | single cell data |
| GSE89042 | GEO | proteins | Leucocyte small nucleolar RNA and microRNA expression data from healthy control subjects o | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| GSE93777 | GEO | proteins | Multi-omics monitoring of drug response in rheumatoid arthritis. | exclude | not a plasma or serum proteome of patients (cells, tissue, other fluid or viral proteins) |
| MSV000086946 | MassIVE | proteins | Persistent immune and clotting dysfunction detected in saliva and blood plasma after COVID | exclude | not a plasma or serum proteome of patients (cells, tissue, other fluid or viral proteins) |
| MSV000087733 | MassIVE | proteins | Combined Transcriptome and Proteome Leukocyte's Profiling Reveals Up-Regulated Module of G | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| MSV000093012 | MassIVE | proteins | Integrated Immunopeptidomic and Proteomic Analysis of COVID-19 lung biopsies. | exclude | not a plasma or serum proteome of patients (cells, tissue, other fluid or viral proteins) |
| 10.17632/2mc6rrc5j3 | Mendeley Data | proteins | COVIDome Datasets | exclude | duplicate: COVIDome datasets, version 2 of 10.17632/73j3nzjxd3 (no contrast) |
| 10.17632/44p9ck395x | Mendeley Data | proteins | Supplementary information for "Plasma Proteomics identify Biomarkers and Pathogenesis of C | exclude | fewer than 10 per class: 8 healthy; severity 5 fatal against 10 mild (severe without intensive care not assigned) |
| 10.17632/6z5zkzb3hv | Mendeley Data | proteins | Omics-driven systems interrogation of metabolic dysregulationin COVID-19 pathogenesis | exclude | other assay or analyte |
| 10.17632/734559nz2k | Mendeley Data | proteins | Systems-level immunomonitoring from acute to recovery phase of severe COVID-19 | exclude | sample list only; no data deposited |
| 10.17632/73j3nzjxd3 | Mendeley Data | proteins | The COVIDome Explorer Researcher Portal. Sullivan et. al. | exclude | COVIDome: COVID-negative participants are symptomatic patients (other disease) and no severity field: no contrast |
| 10.17632/8zfdjsypd8 | Mendeley Data | proteins | The plasma metabolome of  long COVID-19 patients two years after infection | exclude | post-acute only (no acute samples by its record) |
| 10.17632/b34tpw33vh | Mendeley Data | proteins | Long COVID metabolomics | exclude | post-acute only (no acute samples by its record) |
| 10.17632/csx49c78c6 | Mendeley Data | proteins | Quantitative, multiplexed, targeted proteomics for ascertaining variant specific SARS-CoV- | exclude | vaccination study |
| 10.17632/ds6g796xyg | Mendeley Data | proteins | The Immunology of Multisystem Inflammatory Syndrome in Children with COVID-19. Consiglio e | exclude | not adults |
| 10.17632/dvf6yvrg4x | Mendeley Data | proteins | Persistent complement dysregulation with signs of thromboinflammation in active Long Covid | exclude | single cell data |
| 10.17632/dx69469n8x | Mendeley Data | proteins | Proteomic analysis of circulating immune cells identifies novel cellular phenotypes associ | exclude | not a plasma or serum proteome of patients (cells, tissue, other fluid or viral proteins) |
| 10.17632/f8kbg4798h | Mendeley Data | proteins | Cross-platform Clinical Proteomics using the Charité Open Standard for Plasma Proteomics ( | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| 10.17632/fwn3kbbh6y | Mendeley Data | proteins | Landscape of epitopes targeted by T cells in 852 convalescent COVID-19 patients: Meta-anal | exclude | not a plasma or serum proteome of patients (cells, tissue, other fluid or viral proteins) |
| 10.17632/h7zhh3s9jj | Mendeley Data | proteins | COVID project | exclude | in vitro or cell culture |
| 10.17632/hz4g8yddpk | Mendeley Data | proteins | Supplementary Data from: Proteomic analysis of circulating immune cells identifies cellula | exclude | not a plasma or serum proteome of patients (cells, tissue, other fluid or viral proteins) |
| 10.17632/jv9jvm3r98 | Mendeley Data | proteins | Flow cytometry data (2) from: Proteomic analysis of circulating immune cells identifies ce | exclude | not a plasma or serum proteome of patients (cells, tissue, other fluid or viral proteins) |
| 10.17632/mdnb359tp9 | Mendeley Data | proteins | The systemic inflammatory landscape of COVID19 in pregnancy: Extensive serum proteomic pro | exclude | not adults |
| 10.17632/nzd3f5c3j4 | Mendeley Data | proteins | TNFα-producing CD4+ T cells dominate the SARS-CoV-2-specific T cell response in COVID-19 o | exclude | not a plasma or serum proteome of patients (cells, tissue, other fluid or viral proteins) |
| 10.17632/s69k229z2h | Mendeley Data | proteins | Soluble Biomarkers during accute MIS-C, Kawasaky disease and Still's disease before treate | exclude | not adults |
| 10.17632/syvmvw2k2f | Mendeley Data | proteins | LONG COVID METABOLOME | exclude | post-acute only (no acute samples by its record) |
| 10.17632/th35tt4zwm | Mendeley Data | proteins | Immuno-proteomic profiling reveals aberrant immune cell regulation in the airways of indiv | exclude | post-acute only (no acute samples by its record) |
| 10.17632/v65cjwvsjz | Mendeley Data | proteins | Flow cytometry data (3) from: Proteomic analysis of circulating immune cells identifies ce | exclude | not a plasma or serum proteome of patients (cells, tissue, other fluid or viral proteins) |
| 10.17632/v9dkw64s6c | Mendeley Data | proteins | Immune dysregulation and autoreactivity correlate with disease severity in SARS-CoV-2-asso | exclude | single cell data |
| 10.17632/xd9g443btr | Mendeley Data | proteins | Predictive models of patient severity in intensive care units based on serum cytokine prof | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| 10.17632/xss9pdtc4f | Mendeley Data | proteins | Integrated plasma proteomic and single-cell immune signaling network signatures demarcate  | exclude | no per sample protein values deposited (summary tables only) |
| 10.17632/yj5rm3k3py | Mendeley Data | proteins | Flow cytometry data (1) from: Proteomic analysis of circulating immune cells identifies ce | exclude | not a plasma or serum proteome of patients (cells, tissue, other fluid or viral proteins) |
| PXD017710 | PRIDE | proteins | Proteome and Translatome of SARS-CoV-2 infected cells | exclude | in vitro or cell culture |
| PXD020216 | PRIDE | proteins | Proteomics of lung tissues from COVID19 patients | exclude | not blood |
| PXD020601 | PRIDE | proteins | Serum proteomics in COVID-19 patients: Altered coagulation and complement status as a func | exclude | not a plasma or serum proteome of patients (cells, tissue, other fluid or viral proteins) |
| PXD021144 | PRIDE | proteins | Mapping the Host Circulating Exosome Response to SARS-CoV-2 | exclude | other assay or analyte |
| PXD021197 | PRIDE | proteins | COVID-19 proteome profiling reveals resolution of inflammatory pathways and respiratory di | exclude | not a plasma or serum proteome of patients (cells, tissue, other fluid or viral proteins) |
| PXD022475 | PRIDE | proteins | Proteomics and Machine Learning Approaches Reveal a Set of Prognostic Markers for COVID-19 | exclude | no files listed in PRIDE (not public in PRIDE, or hosted elsewhere) |
| PXD022990 | PRIDE | proteins | Synthetic Siglec-9 agonists inhibit neutrophil activation associated with COVID-19 | exclude | not a plasma or serum proteome of patients (cells, tissue, other fluid or viral proteins) |
| PXD023419 | PRIDE | proteins | MS-based sequencing of the anti-FLAG-M2 antibody | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| PXD024707 | PRIDE | proteins | Is there a serum proteome signature to predict mortality in severe COVID-19 patients? | exclude | not adults |
| PXD024728 | PRIDE | proteins | SOD1 is a possible predictor of COVID-19 progression as revealed by plasma proteomics | exclude | MS raw files only (94 files); no abundance matrix deposited |
| PXD025752 | PRIDE | proteins | A time-resolved proteomic and prognostic map of COVID-19 | exclude | MS raw files only (100 files); no abundance matrix deposited |
| PXD025884 | PRIDE | proteins | De novo sequencing of antibody light chain proteoforms from patients with multiple myeloma | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| PXD029143 | PRIDE | proteins | Noninvasive Nasopharyngeal Proteomics Identify Coagulation Abnormalities in SARS-CoV-2 | exclude | vaccination study |
| PXD029376 | PRIDE | proteins | Large-scale multi-omics analysis identifies novel early markers of COVID-19 severity | exclude | glycopeptide identifications only; no per sample protein matrix |
| PXD030662 | PRIDE | proteins | Proteomic and metabolomic profiling of urine uncovers immune responses in patients with CO | exclude | not a plasma or serum proteome of patients (cells, tissue, other fluid or viral proteins) |
| PXD031941 | PRIDE | proteins | Template-based assembly of proteomic short reads for de novo antibody sequencing and reper | exclude | not a plasma or serum proteome of patients (cells, tissue, other fluid or viral proteins) |
| PXD032808 | PRIDE | proteins | Proteome of extracellular vesicles from plasma of COVID-19 patients | exclude | other assay or analyte |
| PXD033335 | PRIDE | proteins | Integrated Multiomics analysis to infer covid 19 biological insights | exclude | MS raw files only (100 files); no abundance matrix deposited |
| PXD034391 | PRIDE | proteins | Human sepsis with plasma samples | exclude | not human |
| PXD036491 | PRIDE | proteins | Proteomic analysis identifies a signature of disease severity in the plasma of COVID-19 pn | exclude | not a plasma or serum proteome of patients (cells, tissue, other fluid or viral proteins) |
| PXD036594 | PRIDE | proteins | prm-PASEF and g-dia-PASEF analysis of plasma samples from colorectal cancer and control pa | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| PXD036969 | PRIDE | proteins | Multi-omics provide evidence for an anti-inflammatory immune signature and metabolic alter | exclude | vaccination study |
| PXD036970 | PRIDE | proteins | Multi-omics provide evidence for an anti-inflammatory immune signature and metabolic alter | exclude | vaccination study |
| PXD036972 | PRIDE | proteins | A multi-omics based anti-inflammatory immune signature characterizes long COVID-19 syndrom | exclude | post-acute only (no acute samples by its record) |
| PXD037676 | PRIDE | proteins | Plasma proteomics in mouse sepsis | exclude | not human (Mus musculus (mouse)) |
| PXD039618 | PRIDE | proteins | Integrative proteomic profiling of lung tissues and blood in acute respiratory distress sy | exclude | not human (Mus musculus) |
| PXD040914 | PRIDE | proteins | Polymorphic Histidine-rich Glycoprotein  LC-MS/MS | exclude | not a plasma or serum proteome of patients (cells, tissue, other fluid or viral proteins) |
| PXD041931 | PRIDE | proteins | Exploring Diverse Clinical Presentations in COVID-19 through Circulating Extracellular Ves | exclude | not human |
| PXD044040 | PRIDE | proteins | Protein modifications in plasma during SARS-CoV-2 infection | exclude | not a plasma or serum proteome of patients (cells, tissue, other fluid or viral proteins) |
| PXD044733 | PRIDE | proteins | Characterization of MPro substrate specificity by HTPS | exclude | in vitro or cell culture |
| PXD046242 | PRIDE | proteins | Identification of direct cleavages on FVII and FXII by MPro | exclude | in vitro or cell culture |
| PXD048139 | PRIDE | proteins | Multi-omic characteristics of longitudinal immune profiling after breakthrough infections  | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| PXD054443 | PRIDE | proteins | Regulatory Role of the VCL/ICAM-1 Pathway in Inflammatory Exudative Damage to the Gas-Bloo | exclude | not human |
| PXD055973 | PRIDE | proteins | Human plasma proteome immune response to COVID-19 | exclude | hospitalised patients only; acquisition setup file, no abundance matrix |
| PXD072061 | PRIDE | proteins | Proteomic profiling of circulating extracellular vesicles from COVID-19 patients and their | exclude | in vitro or cell culture |
| PXD018094 | ProteomeXchange | proteins | Proteomics Characterization of COVID-19 lung tissues | exclude | not a plasma or serum proteome of patients (cells, tissue, other fluid or viral proteins) |
| PXD018117 | ProteomeXchange | proteins | A SARS-CoV-2-Human Protein-Protein Interaction Map Reveals Drug Targets and Potential Drug | exclude | vaccination study |
| PXD018307 | ProteomeXchange | proteins | The Analysis of Plasma Metabolome and Lipidome of Patients with COVID-19 | exclude | no files listed in PRIDE (not public in PRIDE, or hosted elsewhere) |
| PXD018357 | ProteomeXchange | proteins | Inhibition of growth factor signaling prevents SARS-CoV-2 replication | exclude | in vitro or cell culture |
| PXD018594 | ProteomeXchange | proteins | Shotgun proteomics of Vero E6 cells infected by Italy-INMI1 SARS-CoV-2 virus | exclude | not human (Betacoronavirus sp.; Chlorocebus aethiops aethiops) |
| PXD018760 | ProteomeXchange | proteins | TARGETED PROTEOMICS FOR THE DETECTION OF SARS-COV-2 PROTEINS | exclude | not human (Cercopithecus aethiops (Green monkey) (Grivet)) |
| PXD018804 | ProteomeXchange | proteins | Extensive proteomic dataset of Vero E6 cells infected by Italy-INMI1 SARS-CoV-2 virus at D | exclude | not human (Betacoronavirus sp.; Chlorocebus aethiops aethiops) |
| PXD018970 | ProteomeXchange | proteins | Urine proteomics profiling of the COVID-19 | exclude | not blood |
| PXD019106 | ProteomeXchange | proteins | Plasma proteomics dataset of COVID-19 patients | exclude | no files listed in PRIDE (not public in PRIDE, or hosted elsewhere) |
| PXD019113 | ProteomeXchange | proteins | The Global Phosphorylation Landscape of SARS-CoV-2 Infection | exclude | not human (Cercopithecus aethiops (Green monkey) (Grivet)) |
| PXD019221 | ProteomeXchange | proteins | [HPP-2.0] Proteomic and Metabolomic Characterization of COVID-19 Patient Sera | exclude | not a plasma or serum proteome of patients (cells, tissue, other fluid or viral proteins) |
| PXD019232 | ProteomeXchange | proteins | [HPP-2.0] Proteomic and Metabolomic Characterization of COVID-19 Patient Sera MRMHR | exclude | no files listed in PRIDE (not public in PRIDE, or hosted elsewhere) |
| PXD019300 | ProteomeXchange | proteins | SARS-CoV-2_PRM_SRM_targeted_acquisitions | exclude | not human (Severe acute respiratory syndrome coronavirus 2) |
| PXD019645 | ProteomeXchange | proteins | Data, reagents, assays and merits of proteomics for SARS-CoV-2 research and testing | exclude | in vitro or cell culture |
| PXD019648 | ProteomeXchange | proteins | LC-MSMSMass Spectrometric detection of SARS-CoV-2 virus in scrapings of the epithelium of  | exclude | not human (Severe acute respiratory syndrome coronavirus 2) |
| PXD019686 | ProteomeXchange | proteins | Proteotyping of SARS-CoV-2 viral peptides in nasal swab simili | exclude | not human (Betacoronavirus sp.) |
| PXD019968 | ProteomeXchange | proteins | Integrated Proteomics Characterization of SARS-CoV-2 Associated Liver Injury | exclude | not a plasma or serum proteome of patients (cells, tissue, other fluid or viral proteins) |
| PXD019970 | ProteomeXchange | proteins | Tissue proteomics dataset of COVID-19 patients | exclude | not a plasma or serum proteome of patients (cells, tissue, other fluid or viral proteins) |
| PXD020067 | ProteomeXchange | proteins | interactome and proteomic survey of PMBCs from COVID-19 patients | exclude | not a plasma or serum proteome of patients (cells, tissue, other fluid or viral proteins) |
| PXD020183 | ProteomeXchange | proteins | Phosphoproteomic Survey of Human iPSC-derived Lung Alveolar Type 2 Cells Reveals Cytopatho | exclude | in vitro or cell culture |
| PXD020339 | ProteomeXchange | proteins | Integrative Vectors for Regulated Expression of SARS-CoV-2 Proteins Implicated in RNA Meta | exclude | in vitro or cell culture |
| PXD020354 | ProteomeXchange | proteins | In-depth plasma proteome profiling of COVID-19 patients | exclude | not a plasma or serum proteome of patients (cells, tissue, other fluid or viral proteins) |
| PXD020394 | ProteomeXchange | proteins | Quantitative proteomic dataset from oro- and naso-pharyngeal swabs used for COVID-19 diagn | exclude | not blood |
| PXD020444 | ProteomeXchange | proteins | Urinary peptides enable assessment of disease severity in COVID-19 | exclude | not blood |
| PXD020580 | ProteomeXchange | proteins | Proteomic Investigation Reveals Dominant Alterations of Neutrophil Degranulation and mRNA  | exclude | not blood |
| PXD020994 | ProteomeXchange | proteins | Discovery of drugs to treat cytokine storm-induced cardiac dysfunction | exclude | not human |
| PXD021073 | ProteomeXchange | proteins | Metabolic Remodeling in Recovered COVID-19 Patients | exclude | no files listed in PRIDE (not public in PRIDE, or hosted elsewhere) |
| PXD021103 | ProteomeXchange | proteins | Omics study reveals abnormal alterations of breastmilk proteins and  metabolites in puerpe | exclude | not adults |
| PXD021152 | ProteomeXchange | proteins | C4PR_LIV: N-terminomic analysis of SARS-CoV-2-infected ACE2-A549 cells (Enriched) | exclude | in vitro or cell culture |
| PXD021153 | ProteomeXchange | proteins | C4PR_LIV: N-terminomic analysis of SARS-CoV-2-infected Vero E6 cells (Enriched) | exclude | not human (Cercopithecus aethiops (Green monkey) (Grivet)) |
| PXD021154 | ProteomeXchange | proteins | C4PR_LIV: N-terminomic analysis of SARS-CoV-2-infected Vero E6 cells (Unenriched) | exclude | not human (Cercopithecus aethiops (Green monkey) (Grivet)) |
| PXD021193 | ProteomeXchange | proteins | Breastmilk proteomics dataset of puerperant women with COVID-19 | exclude | not a plasma or serum proteome of patients (cells, tissue, other fluid or viral proteins) |
| PXD021388 | ProteomeXchange | proteins | Detection of SARS-CoV-2 in nasal swabs using MALDI-MS | exclude | not blood |
| PXD021402 | ProteomeXchange | proteins | C4PR_LIV: N-terminomic analysis of SARS-CoV-2-infected ACE2-A549 cells (Enriched - Variabl | exclude | in vitro or cell culture |
| PXD021403 | ProteomeXchange | proteins | C4PR_LIV: N-terminomic analysis of SARS-CoV-2-infected Vero E6 cells (Enriched - Variable  | exclude | not human (Cercopithecus aethiops (Green monkey) (Grivet)) |
| PXD021702 | ProteomeXchange | proteins | Identification of serum prognostic biomarkers of severe COVID-19 using a quantitative prot | exclude | no files listed in PRIDE (not public in PRIDE, or hosted elsewhere) |
| PXD021954 | ProteomeXchange | proteins | Serum proteomic features of active and recovered COVID-19 patients: potential diagnostic a | exclude | no files listed in PRIDE (not public in PRIDE, or hosted elsewhere) |
| PXD022017 | ProteomeXchange | proteins | Comparative multiplexed interactomics of SARS-CoV-2 and homologous coronavirus non-structu | exclude | not a plasma or serum proteome of patients (cells, tissue, other fluid or viral proteins) |
| PXD022085 | ProteomeXchange | proteins | Proteomic analysis of bronchoalveolar lavage fluid in COVID-19 patients | exclude | not blood |
| PXD022086 | ProteomeXchange | proteins | A proximity labeling map of SARS-CoV-2 and human | exclude | not a plasma or serum proteome of patients (cells, tissue, other fluid or viral proteins) |
| PXD022273 | ProteomeXchange | proteins | Aberrant metabolism is associated with fast antibody fading in COVID-19 convalescent patie | exclude | vaccination study |
| PXD022296 | ProteomeXchange | proteins | Mass-spectrometry based plasma proteomics identifies biomarkers for the detection of COVID | exclude | MS raw files only (76 files); no abundance matrix deposited |
| PXD022751 | ProteomeXchange | proteins | Proteomic Analysis of Inactivated and Sterilized Sera might be Applied to Monitoring the D | exclude | no files listed in PRIDE (not public in PRIDE, or hosted elsewhere) |
| PXD022817 | ProteomeXchange | proteins | COVIDomics_plasma_2020_DDA-PASEF | exclude | COVIDome (Colorado) MS plasma: negatives are symptomatic patients (other disease) and no severity field |
| PXD022889 | ProteomeXchange | proteins | Proteomic signature of host response to SARS-CoV-2 infection in the nasopharynx | exclude | not blood |
| PXD023016 | ProteomeXchange | proteins | Proteomic Investigation Reveals Dominant Alterations of Neutrophil Degranulation and mRNA  | exclude | not blood |
| PXD023099 | ProteomeXchange | proteins | Molecular Phenomics of the Gut Ecosystem in COVID-19 Patients | exclude | not a plasma or serum proteome of patients (cells, tissue, other fluid or viral proteins) |
| PXD023154 | ProteomeXchange | proteins | Lipidomics data of cured COVID-19 patients | exclude | no files listed in PRIDE (not public in PRIDE, or hosted elsewhere) |
| PXD023175 | ProteomeXchange | proteins | A blood atlas of COVID-19 defines hallmarks of disease severity and specificity | exclude | not a plasma or serum proteome of patients (cells, tissue, other fluid or viral proteins) |
| PXD023305 | ProteomeXchange | proteins | Proteomic analysis identifies prolonged disturbances in pathways related to cholesterol me | exclude | no files listed in PRIDE (not public in PRIDE, or hosted elsewhere) |
| PXD023422 | ProteomeXchange | proteins | Fluvastatin mitigates SARS-CoV-2 infection in human lung cells | exclude | in vitro or cell culture |
| PXD023521 | ProteomeXchange | proteins | A longitudinal study of COVID-19 patients using proteomics approach | exclude | SDRF only; no abundance matrix deposited |
| PXD023538 | ProteomeXchange | proteins | C4PR_LIV: N-terminomic analysis of SARS-CoV-2-infected Vero E6 cells treated with protease | exclude | not human (Cercopithecus aethiops (Green monkey) (Grivet)) |
| PXD023539 | ProteomeXchange | proteins | C4PR_LIV: N-terminomic analysis of SARS-CoV-2-infected Vero E6 cells treated with protease | exclude | not human (Cercopithecus aethiops (Green monkey) (Grivet)) |
| PXD023834 | ProteomeXchange | proteins | A type I IFN, prothrombotic hyperinflammatory neutrophil signature is distinct for COVID-1 | exclude | not a plasma or serum proteome of patients (cells, tissue, other fluid or viral proteins) |
| PXD024026 | ProteomeXchange | proteins | SARS-CoV-2 RNAemia and proteomic biomarker trajectory inform prognostication in COVID-19 p | exclude | not a plasma or serum proteome of patients (cells, tissue, other fluid or viral proteins) |
| PXD024087 | ProteomeXchange | proteins | Kinome and phosphoproteome reprogramming underlie the immunological abnormalities in sever | exclude | not a plasma or serum proteome of patients (cells, tissue, other fluid or viral proteins) |
| PXD024089 | ProteomeXchange | proteins | SARS-CoV-2 RNAemia and proteomic biomarker trajectory inform prognostication in COVID-19 p | exclude | other assay or analyte |
| PXD024410 | ProteomeXchange | proteins | Distinct Molecular Signature of Asymptomatic COVID-19 Revealed by Lipidomic Characterizati | exclude | no files listed in PRIDE (not public in PRIDE, or hosted elsewhere) |
| PXD024549 | ProteomeXchange | proteins | Serum quantitative proteomics identifies prognostic biomarkers and physiological processes | exclude | vaccination study |
| PXD024566 | ProteomeXchange | proteins | Comparative host interactomes of the SARS-CoV-2 nonstructural protein 3 and human coronavi | exclude | not a plasma or serum proteome of patients (cells, tissue, other fluid or viral proteins) |
| PXD024674 | ProteomeXchange | proteins | The Trans-omics Landscape of COVID-19 | exclude | MS raw files only (100 files); no abundance matrix deposited |
| PXD024776 | ProteomeXchange | proteins | Proteomic characterization of COVID-19 patient Breast milk | exclude | not a plasma or serum proteome of patients (cells, tissue, other fluid or viral proteins) |
| PXD025015 | ProteomeXchange | proteins | Divergent and self-reactive immune responses in the CNS of COVID-19 patients with neurolog | exclude | single cell data |
| PXD025016 | ProteomeXchange | proteins | Divergent and self-reactive immune responses in the CNS of COVID-19 patients with neurolog | exclude | single cell data |
| PXD025035 | ProteomeXchange | proteins | Divergent and self-reactive immune responses in the CNS of COVID-19 patients with neurolog | exclude | single cell data |
| PXD025130 | ProteomeXchange | proteins | Proteomics characterization of purified SARS-CoV-2 viral particles for search of mutations | exclude | not human (Severe acute respiratory syndrome coronavirus 2) |
| PXD025131 | ProteomeXchange | proteins | Proteomics characterization of purified SARS-CoV-2 viral particles. | exclude | not human (Severe acute respiratory syndrome coronavirus 2) |
| PXD025265 | ProteomeXchange | proteins | Drivers genes in severe forms of COVID-19 | exclude | not a plasma or serum proteome of patients (cells, tissue, other fluid or viral proteins) |
| PXD025294 | ProteomeXchange | proteins | Targeted proteomics as a tool to detect SARS-CoV-2 proteins in clinical specimens | exclude | in vitro or cell culture |
| PXD025453 | ProteomeXchange | proteins | Clinical and immunological benefits of convalescent plasma therapy in severe COVID-19: ins | exclude | not a plasma or serum proteome of patients (cells, tissue, other fluid or viral proteins) |
| PXD025554 | ProteomeXchange | proteins | Glycoproteomics profiling of recombinant SARS-CoV-2 stabilized Spike construct | exclude | not human (Severe acute respiratory syndrome coronavirus 2) |
| PXD026079 | ProteomeXchange | proteins | Differential interferon-Î± subtype immune signatures suppress SARS-CoV-2 infection | exclude | not blood |
| PXD026091 | ProteomeXchange | proteins | COVID-19 patient plasma exosome proteomics | exclude | other assay or analyte |
| PXD026412 | ProteomeXchange | proteins | FIRST DETECTION OF SARS-COV-2 PROTEINS IN WASTEWATER SAMPLES BY MASS SPECTROMETRY | exclude | not human (Coronavirus PREDICT CoV-2) |
| PXD026431 | ProteomeXchange | proteins | Protein N-glycosylation of SARS-CoV-2 infection | exclude | not human (Severe acute respiratory syndrome coronavirus 2) |
| PXD026442 | ProteomeXchange | proteins | Immune response pattern across the asymptomatic, symptomatic and convalescent periods of C | exclude | not a plasma or serum proteome of patients (cells, tissue, other fluid or viral proteins) |
| PXD026510 | ProteomeXchange | proteins | Plasma gradient of soluble urokinase-type plasminogen activator receptor is linked to path | exclude | MS raw files only (81 files); no abundance matrix deposited |
| PXD026568 | ProteomeXchange | proteins | DeHart_RBD_Manuscript_Resubmission_Dataset | exclude | not human (Sars-CoV-2) |
| PXD026703 | ProteomeXchange | proteins | Viral fallout: Semen proteome reveals COVID-19 dysregulates several biological processes l | exclude | not a plasma or serum proteome of patients (cells, tissue, other fluid or viral proteins) |
| PXD026848 | ProteomeXchange | proteins | Validation of proteins showing significant dysregulation from shotgun proteomics experimen | exclude | no files listed in PRIDE (not public in PRIDE, or hosted elsewhere) |
| PXD026849 | ProteomeXchange | proteins | Targeted validation of significantly altered semen proteins in COVID-19 patients | exclude | not a plasma or serum proteome of patients (cells, tissue, other fluid or viral proteins) |
| PXD026962 | ProteomeXchange | proteins | Multi-omic Profiling of Plasma Identify Biomarkers and Pathogenesis of COVID-19 in Childre | exclude | not adults |
| PXD027168 | ProteomeXchange | proteins | Affinity Purification Mass Spectrometry Analysis of Mutant and Wildtype SARS-CoV-2 Nucleoc | exclude | in vitro or cell culture |
| PXD027557 | ProteomeXchange | proteins | Pulmonary and Renal Long COVID at Two-year Revisit | exclude | post-acute only (no acute samples by its record) |
| PXD027582 | ProteomeXchange | proteins | Absence of COVID-19-associated Changes in Coagulation and Thrombosis in the Ferret Model | exclude | not human |
| PXD027605 | ProteomeXchange | proteins | Structural O-Glycoform Heterogeneity of the SARS-CoV-2 Spike Protein Receptor-Binding Doma | exclude | not a plasma or serum proteome of patients (cells, tissue, other fluid or viral proteins) |
| PXD027959 | ProteomeXchange | proteins | Targeted MS based multi-omic analysis of blood plasma from COVID-19 patients to predict su | exclude | no files listed in PRIDE (not public in PRIDE, or hosted elsewhere) |
| PXD027989 | ProteomeXchange | proteins | Proteomics analysis of  mass spectrometry data | exclude | no files listed in PRIDE (not public in PRIDE, or hosted elsewhere) |
| PXD028429 | ProteomeXchange | proteins | Proteomics of serum samples from COVID-19 ARDS and bacterial ARDS patients | exclude | single cell data |
| PXD028505 | ProteomeXchange | proteins | A SISCAPA-based approach for detection of SARS-CoV-2 viral antigens from clinical samples. | exclude | not blood |
| PXD029007 | ProteomeXchange | proteins | Plasma-derived Exosome-Like Vesicles (ELVs) of COVID-19 patients expose SARS-CoV-2-Spike-d | exclude | other assay or analyte |
| PXD029023 | ProteomeXchange | proteins | Proteomics dataset of COVID-19 and control patients | exclude | no files listed in PRIDE (not public in PRIDE, or hosted elsewhere) |
| PXD029291 | ProteomeXchange | proteins | N1-methylpseudouridine found within COVID-19 mRNA vaccines produces faithful protein produ | exclude | not human (Severe acute respiratory syndrome coronavirus 2) |
| PXD029300 | ProteomeXchange | proteins | Detection of mutated peptides of SARS-CoV-2 from COVID-19-infected patients | exclude | not a plasma or serum proteome of patients (cells, tissue, other fluid or viral proteins) |
| PXD029321 | ProteomeXchange | proteins | Anti-SARS-CoV-2 antibodies elicited by COVID-19mRNA vaccine exhibits a unique glycosylatio | exclude | vaccination study |
| PXD029426 | ProteomeXchange | proteins | SARS-CoV-2 Infection Impacts Glutamate/Glutamine Metabolism in Syrian Hamster Central Nerv | exclude | not human (Mesocricetus auratus (Golden hamster)) |
| PXD029437 | ProteomeXchange | proteins | Longitudinal Plasma Proteomics Analysis of Acute COVID-19 | exclude | acute COVID-19 patients only (no controls, no severity classes); peptide-level MRM table |
| PXD029547 | ProteomeXchange | proteins | Neutrophil and monocyte dysfunctional effector response towards bacterial challenge in cri | exclude | not a plasma or serum proteome of patients (cells, tissue, other fluid or viral proteins) |
| PXD029567 | ProteomeXchange | proteins | Increased soluble HLA in COVID-19 present a disease-related, diverse immunopeptidome assoc | exclude | not a plasma or serum proteome of patients (cells, tissue, other fluid or viral proteins) |
| PXD029740 | ProteomeXchange | proteins | Proteomics of Cerebrospinal fluid from COVID-19 Patients with Neurological Manifestations | exclude | not a plasma or serum proteome of patients (cells, tissue, other fluid or viral proteins) |
| PXD029756 | ProteomeXchange | proteins | Protein Posttranslational Signatures Identified in COVID-19 Patient Plasma | exclude | not a plasma or serum proteome of patients (cells, tissue, other fluid or viral proteins) |
| PXD029800 | ProteomeXchange | proteins | Identification of Endogenous Peptides in MALDI-TOF-MS based COVID-19 diagnostics | exclude | not blood |
| PXD030012 | ProteomeXchange | proteins | Detection SARS-CoV-2 infection in mice by ATR-FTIR | exclude | not human |
| PXD030753 | ProteomeXchange | proteins | The cGAS-STING pathway drives type I IFN immunopathology in COVID-19 | exclude | not a plasma or serum proteome of patients (cells, tissue, other fluid or viral proteins) |
| PXD030785 | ProteomeXchange | proteins | Protective immune trajectories in early viral containment of non-pneumonic SARS-CoV-2 infe | exclude | not a plasma or serum proteome of patients (cells, tissue, other fluid or viral proteins) |
| PXD031241 | ProteomeXchange | proteins | Post-covid surgical complications in adipose tissue | exclude | post-acute only (no acute samples by its record) |
| PXD031249 | ProteomeXchange | proteins | Spatial-region-resolved proteomic map reveals mechanism of COVID-19-associated heart injur | exclude | not blood |
| PXD031251 | ProteomeXchange | proteins | PLATELET PROTEOME REVEALS FEATURES OF CELL DEATH, ANTIVIRAL RESPONSE AND VIRAL REPLICATION | exclude | not a plasma or serum proteome of patients (cells, tissue, other fluid or viral proteins) |
| PXD031969 | ProteomeXchange | proteins | Serum proteomics covid-19 vaccine | exclude | vaccination study |
| PXD032127 | ProteomeXchange | proteins | Proteomic and metabolomic signatures associated with the immune response in healthy indivi | exclude | vaccination study |
| PXD032150 | ProteomeXchange | proteins | Proteome Map of Human nasopharyngeal/oropharyngeal swab samples in VTM media following SAR | exclude | not blood |
| PXD032677 | ProteomeXchange | proteins | MultiOMICs landscape of SARS-CoV2-induced host responses in human lung epithelial cells | exclude | in vitro or cell culture |
| PXD032937 | ProteomeXchange | proteins | Cortisol binding sites on SARS-CoV-2 S1 identified by limited proteolysis-coupled mass spe | exclude | not human (Severe acute respiratory syndrome coronavirus 2) |
| PXD033015 | ProteomeXchange | proteins | SARS-CoV-2 phosphoproteomics - ACE2-A549 and Vero | exclude | in vitro or cell culture |
| PXD033293 | ProteomeXchange | proteins | Proteomic analysis neutrophils from COVID-19 patients | exclude | not a plasma or serum proteome of patients (cells, tissue, other fluid or viral proteins) |
| PXD033558 | ProteomeXchange | proteins | Multi-omics Characterization of Neutrophil Extracellular Trap Formation in Severe and Mild | exclude | not a plasma or serum proteome of patients (cells, tissue, other fluid or viral proteins) |
| PXD033661 | ProteomeXchange | proteins | Persistent Dysfunction of Innate Immune and Clotting Factors Detected in Saliva and Blood  | exclude | not a plasma or serum proteome of patients (cells, tissue, other fluid or viral proteins) |
| PXD033790 | ProteomeXchange | proteins | SARS-CoV-2 PRM-MS Assay with Clinical Patient Samples | exclude | not blood |
| PXD033979 | ProteomeXchange | proteins | COVID-19 PATIENTS WITH NEUROLOGICAL SYMPTOMS PRESENT BROAD CENTRAL NERVOUS SYSTEM MOLECULA | exclude | not a plasma or serum proteome of patients (cells, tissue, other fluid or viral proteins) |
| PXD034045 | ProteomeXchange | proteins | Shotgun LC-MS analysis of human serum from COVID-19 patients | exclude | vaccination study |
| PXD034086 | ProteomeXchange | proteins | Structural- and Site-specific N-glycosylation Characterization of COVID-19 Virus Spike wit | exclude | not a plasma or serum proteome of patients (cells, tissue, other fluid or viral proteins) |
| PXD034510 | ProteomeXchange | proteins | SARS-CoV-2 infects adipose tissue in a fat depot- and viral lineage-dependent manner | exclude | not a plasma or serum proteome of patients (cells, tissue, other fluid or viral proteins) |
| PXD034648 | ProteomeXchange | proteins | SARS-CoV-2 spike sequence dictates host TP53 activity and COVID-19 pathogenicity | exclude | not human (Chlorocebus sabaeus) |
| PXD034949 | ProteomeXchange | proteins | urine proteome of QIV and COVID-19 | exclude | vaccination study |
| PXD035168 | ProteomeXchange | proteins | Serum Proteome of COVID-19 Patients with Different Severities | exclude | no files listed in PRIDE (not public in PRIDE, or hosted elsewhere) |
| PXD035247 | ProteomeXchange | proteins | Cerebrospinal fluid protein markers indicate neuro-damage in SARS-CoV-2-infected non-human | exclude | not human (Macaca mulatta (Rhesus macaque); Cercopithecus aethiops (Green monkey) (Grivet)) |
| PXD035451 | ProteomeXchange | proteins | SARS-CoV-2 hijacks p38Î²/MAPK11 to promote virus replication | exclude | in vitro or cell culture |
| PXD035471 | ProteomeXchange | proteins | Temporal proteomics characterization of COVID-19 patient sera | exclude | no files listed in PRIDE (not public in PRIDE, or hosted elsewhere) |
| PXD035474 | ProteomeXchange | proteins | Multi-organ proteomic and phosphoproteomic landscape of K18-hACE2 mice infected with SARS- | exclude | not human (Mus musculus) |
| PXD035510 | ProteomeXchange | proteins | urine metabolomics of COVID-19 patients | exclude | not blood |
| PXD035709 | ProteomeXchange | proteins | Spike S1 domain interactome in non-pulmonary systems: a role beyond the receptor recogniti | exclude | in vitro or cell culture |
| PXD035904 | ProteomeXchange | proteins | Distinct Core Glycan and O-Glycoform Utilization of SARS-CoV-2 Omicron Variant Spike Prote | exclude | vaccination study |
| PXD036159 | ProteomeXchange | proteins | Serum Peptidome Profiles Immune Response of COVID-19 Vaccine Administration | exclude | vaccination study |
| PXD036193 | ProteomeXchange | proteins | Temporal proteomic characterization of COVID-19 infected mouse lungs | exclude | not human (Mus musculus (Mouse)) |
| PXD036325 | ProteomeXchange | proteins | A yeast-based system to study SARS-CoV-2 Mpro structure and to identify nirmatrelvir resis | exclude | not human (Caenorhabditis elegans; Saccharomyces cerevisiae) |
| PXD036571 | ProteomeXchange | proteins | COVID-19 salivary protein profile: unravelling molecular aspects of SARS-COV-2 infection | exclude | not blood |
| PXD036590 | ProteomeXchange | proteins | Plasma proteomic signature predicts who will get persistent symptoms following SARS-CoV-2  | exclude | not a plasma or serum proteome of patients (cells, tissue, other fluid or viral proteins) |
| PXD036798 | ProteomeXchange | proteins | SARS-CoV-2 variants evolve convergent strategies to remodel the host response: Affinity Pu | exclude | in vitro or cell culture |
| PXD036821 | ProteomeXchange | proteins | The accessory protein ORF6 contributes to SARS-CoV-2 pathogenesis by regulating host and v | exclude | in vitro or cell culture |
| PXD036862 | ProteomeXchange | proteins | Quantitative proteomics analysis of COVID-19 infected hamster lung tissue | exclude | not human (Mesocricetus auratus (Golden hamster)) |
| PXD036968 | ProteomeXchange | proteins | SARS-CoV-2 variants evolve convergent strategies to remodel the host response: Infection p | exclude | not human (Severe acute respiratory syndrome coronavirus 2) |
| PXD037204 | ProteomeXchange | proteins | A simple protein-based model for predicting liver dysfunction in COVID-19 patients with no | exclude | not a plasma or serum proteome of patients (cells, tissue, other fluid or viral proteins) |
| PXD037451 | ProteomeXchange | proteins | Mass-spectrometry proteomics of Calu-3 cell lines infected with SARS-CoV-2, in comparison  | exclude | in vitro or cell culture |
| PXD037486 | ProteomeXchange | proteins | In-depth serum proteomics of hospitalized COVID-19 patients and PCR-negative healthy contr | exclude | 7 healthy controls (fewer than 10 per class) |
| PXD038220 | ProteomeXchange | proteins | In Depth Characterization of Protein N-Glycosylation for a COVID-19 Variant-Design Vaccine | exclude | not human (Cricetulus griseus) |
| PXD038629 | ProteomeXchange | proteins | Selennoprotein S interacts with the replication and transcription complex of the SARS-CoV- | exclude | not blood |
| PXD038693 | ProteomeXchange | proteins | The central nervous systemâ€™s proteogenomic and spatial imprint upon systemic viral infec | exclude | not a plasma or serum proteome of patients (cells, tissue, other fluid or viral proteins) |
| PXD038748 | ProteomeXchange | proteins | COVID-19-associated monocytic encephalitis (CAME): histological and proteomic evidence fro | exclude | not a plasma or serum proteome of patients (cells, tissue, other fluid or viral proteins) |
| PXD039058 | ProteomeXchange | proteins | Proteomic characteristics and diagnostic potential of exhaled breath particles in patients | exclude | MS raw files only (36 files); no abundance matrix deposited |
| PXD039095 | ProteomeXchange | proteins | Proteomics Characterization of human COVID-19 tissues | exclude | not a plasma or serum proteome of patients (cells, tissue, other fluid or viral proteins) |
| PXD039598 | ProteomeXchange | proteins | Metaproteomic analysis of nasopharyngeal swab samples to identify microbial peptides and p | exclude | not blood |
| PXD039722 | ProteomeXchange | proteins | Proteomics Landscape of Diverse Serological Patterns in COVID-19 | exclude | no files listed in PRIDE (not public in PRIDE, or hosted elsewhere) |
| PXD040175 | ProteomeXchange | proteins | Clinical and pulmonary function analysis in long-COVID revealed that long-term pulmonary d | exclude | not a plasma or serum proteome of patients (cells, tissue, other fluid or viral proteins) |
| PXD040245 | ProteomeXchange | proteins | TMT labeling and LC-MS/MS strategy to quantify differentially abundant proteins among COVI | exclude | not a plasma or serum proteome of patients (cells, tissue, other fluid or viral proteins) |
| PXD040703 | ProteomeXchange | proteins | Proteomic analysis of circulating immune cells identifies novel cellular phenotypes associ | exclude | not a plasma or serum proteome of patients (cells, tissue, other fluid or viral proteins) |
| PXD040982 | ProteomeXchange | proteins | Mass-spectrometry phosphoproteomics of Calu-3 cell lines infected with SARS-CoV-2, in comp | exclude | in vitro or cell culture |
| PXD041142 | ProteomeXchange | proteins | Urine metabolomics analysis of patients recovered of the Omicron variant of COVID-19 | exclude | not blood |
| PXD041438 | ProteomeXchange | proteins | Cross-linking mass spectrometry uncovers interactions between high-density lipoproteins an | exclude | not a plasma or serum proteome of patients (cells, tissue, other fluid or viral proteins) |
| PXD041681 | ProteomeXchange | proteins | Changes in the proteome of platelets from patients with critical progression of COVID-19 | exclude | not a plasma or serum proteome of patients (cells, tissue, other fluid or viral proteins) |
| PXD041762 | ProteomeXchange | proteins | PROTEOMICS EVOLUTION FROM ACUTE TO POST-COVID-19 CONDITION | exclude | not a plasma or serum proteome of patients (cells, tissue, other fluid or viral proteins) |
| PXD041831 | ProteomeXchange | proteins | Cross-linking mass spectrometry uncovers interactions between high-density lipoproteins an | exclude | not a plasma or serum proteome of patients (cells, tissue, other fluid or viral proteins) |
| PXD041854 | ProteomeXchange | proteins | Plitidepsin is a broad-spectrum antiviral that transiently inhibits viral 1 protein transl | exclude | not human (Chlorocebus) |
| PXD041904 | ProteomeXchange | proteins | Semen proteome analysis comparing the persistent effects of different SARS-COV-2 variants  | exclude | vaccination study |
| PXD041909 | ProteomeXchange | proteins | Serum proteomics hint at an early T-cell response and modulation of SARS-CoV-2-related pat | exclude | no files listed in PRIDE (not public in PRIDE, or hosted elsewhere) |
| PXD042106 | ProteomeXchange | proteins | Plasma proteomics-based biomarkers for predicting response to mesenchymal stem cell therap | exclude | no files listed in PRIDE (not public in PRIDE, or hosted elsewhere) |
| PXD042148 | ProteomeXchange | proteins | A viral pan-end RNA element and host complex define a SARS-CoV-2 regulon | exclude | not a plasma or serum proteome of patients (cells, tissue, other fluid or viral proteins) |
| PXD042938 | ProteomeXchange | proteins | Proteomic changes of CD19+ B Lymphocytes in COVID-19 Patients | exclude | not a plasma or serum proteome of patients (cells, tissue, other fluid or viral proteins) |
| PXD043058 | ProteomeXchange | proteins | Extracellular Vesicles in COVID-19 convalescence can regulate T cell metabolism and functi | exclude | other assay or analyte |
| PXD043740 | ProteomeXchange | proteins | Personalized serum proteome profiles of COVID-19-associated acute respiratory distress syn | exclude | not a plasma or serum proteome of patients (cells, tissue, other fluid or viral proteins) |
| PXD043807 | ProteomeXchange | proteins | Deep spatial proteomics reveals region-specific features of severe COVID-19-related pulmon | exclude | not a plasma or serum proteome of patients (cells, tissue, other fluid or viral proteins) |
| PXD045427 | ProteomeXchange | proteins | Proteomic study of lung tissues from COVID-19 | exclude | not blood |
| PXD046103 | ProteomeXchange | proteins | Protein and cytokine profiles associated with COVID-19 disease severity in Puerto Rico | exclude | 22 patients across severity groups and few negative controls (fewer than 10 per class) |
| PXD046411 | ProteomeXchange | proteins | Mapping the immunopeptidome of several SARS-CoV-2 antigens across common HLA haplotypes | exclude | in vitro or cell culture |
| PXD047427 | ProteomeXchange | proteins | SARS-CoV-2-Infected hiPSC-Derived Cardiomyocytes Reveal Dynamic Changes in the COVID-19 He | exclude | in vitro or cell culture |
| PXD047446 | ProteomeXchange | proteins | Proteomic analysis of salivary extracellular vesicles from COVID-19 patients reveals a spe | exclude | other assay or analyte |
| PXD047923 | ProteomeXchange | proteins | Deep-proteomic of 52 serum samples from 52 COVID-19 patients and 30 matched healthy contro | exclude | in vitro or cell culture |
| PXD048080 | ProteomeXchange | proteins | Mutation T9I in Envelope conveys autophagy resistance of SARS-CoV-2 Omicron | exclude | not a plasma or serum proteome of patients (cells, tissue, other fluid or viral proteins) |
| PXD048192 | ProteomeXchange | proteins | Plasma proteomic signatures of COVID-19 at the end of 2022 in China | exclude | no files listed in PRIDE (not public in PRIDE, or hosted elsewhere) |
| PXD050303 | ProteomeXchange | proteins | Clinical features and Multiomics Profiles Indicate Coagulation and Platelet Dysfunction as | exclude | not a plasma or serum proteome of patients (cells, tissue, other fluid or viral proteins) |
| PXD050384 | ProteomeXchange | proteins | SARS-CoV-2 Spike Protein in the Skull-Meninges-Brain Axis | exclude | not human |
| PXD050392 | ProteomeXchange | proteins | Glycosylation  analysis of SARS-CoV-2 spike protein | exclude | not human (Severe acute respiratory syndrome coronavirus 2) |
| PXD050512 | ProteomeXchange | proteins | COVID-19 infection on human semen quality | exclude | not a plasma or serum proteome of patients (cells, tissue, other fluid or viral proteins) |
| PXD051205 | ProteomeXchange | proteins | Integrated multiomics implicates dysregulation of ECM and cell adhesion pathways as driver | exclude | single cell data |
| PXD051702 | ProteomeXchange | proteins | Proteomic insights into the impact of COVID-19 on peripheral blood mononuclear cells in ki | exclude | not a plasma or serum proteome of patients (cells, tissue, other fluid or viral proteins) |
| PXD051948 | ProteomeXchange | proteins | Differential antibody capture pneumonia samples Munich cohort | exclude | not a plasma or serum proteome of patients (cells, tissue, other fluid or viral proteins) |
| PXD051982 | ProteomeXchange | proteins | Differential antibody capture pneumonia samples Chicago cohort | exclude | not a plasma or serum proteome of patients (cells, tissue, other fluid or viral proteins) |
| PXD051984 | ProteomeXchange | proteins | Functional Heterogeneity Analysis of Mesenchymal Stem Cell Lines and their Therapeutic Pot | exclude | not human |
| PXD052187 | ProteomeXchange | proteins | Integrating Machine Learning-Enhanced Immunopeptidomics and SARS-CoV-2 Population-Scale An | exclude | vaccination study |
| PXD052222 | ProteomeXchange | proteins | Proteomics and metabolomics profiles of unvaccinated nonagenarian patients with severe SAR | exclude | not a plasma or serum proteome of patients (cells, tissue, other fluid or viral proteins) |
| PXD053240 | ProteomeXchange | proteins | Time-based quantitative proteomic and phosphoproteomic analysis of A549-ACE2 cells during  | exclude | in vitro or cell culture |
| PXD053895 | ProteomeXchange | proteins | Serum Proteomics Identifies Biomarkers for Predicting Non-survivors in Elderly COVID-19 Pa | exclude | not a plasma or serum proteome of patients (cells, tissue, other fluid or viral proteins) |
| PXD054133 | ProteomeXchange | proteins | 4D-DIA Proteomics Uncovers New Insights into Host Salivary Response Following SARS-CoV-2 O | exclude | not blood |
| PXD054533 | ProteomeXchange | proteins | Proteomic analysis of 442 clinical plasma samples from healthy donors and COVID-19 covales | exclude | no files listed in PRIDE (not public in PRIDE, or hosted elsewhere) |
| PXD054718 | ProteomeXchange | proteins | Glycoprofile Comparison of SARS-CoV-2 Spike Proteins Expressed in CHO and HEK Cell Lines | exclude | in vitro or cell culture |
| PXD055810 | ProteomeXchange | proteins | Neutrophils display novel cytosolic proliferating cell nuclear antigen (PCNA) partners inv | exclude | in vitro or cell culture |
| PXD057019 | ProteomeXchange | proteins | Glycosylation signature of plasma IgA of critically ill COVID-19 patients | exclude | not a plasma or serum proteome of patients (cells, tissue, other fluid or viral proteins) |
| PXD057158 | ProteomeXchange | proteins | Metaproteomics reveals age-specific alterations of gut microbiome in hamsters with SARS-Co | exclude | not human (Mesocricetus auratus (Golden hamster)) |
| PXD057656 | ProteomeXchange | proteins | Caspase-8â€™s non-apoptotic role is critical for orchestrating exaggerated inflammation du | exclude | not human (Mus musculus (Mouse)) |
| PXD058234 | ProteomeXchange | proteins | Unveiling Multi-Omics and Immune Landscapes in COVID-19 Patients' Ovaries Microenvironment | exclude | not a plasma or serum proteome of patients (cells, tissue, other fluid or viral proteins) |
| PXD058743 | ProteomeXchange | proteins | Targeting G9a Regulated m6A Modification Pathway to Hinder SARS-CoV-2 Replication | exclude | in vitro or cell culture |
| PXD058750 | ProteomeXchange | proteins | Structural and Functional Insights into the Evolution of SARS-CoV-2 KP.3.1.1 Spike Protein | exclude | not a plasma or serum proteome of patients (cells, tissue, other fluid or viral proteins) |
| PXD060880 | ProteomeXchange | proteins | SARS-CoV-2 non-spike structural proteins hijack mucosa epithelial cell fate | exclude | not a plasma or serum proteome of patients (cells, tissue, other fluid or viral proteins) |
| PXD061500 | ProteomeXchange | proteins | Biomarker identification in hospitalised patients with COVID-19: a cohort study of acute i | exclude | MS raw files only (32 files); no abundance matrix deposited |
| PXD063172 | ProteomeXchange | proteins | integrating transcriptomic, proteomic, and interactome identify potential therapeutic targ | exclude | not a plasma or serum proteome of patients (cells, tissue, other fluid or viral proteins) |
| PXD066966 | ProteomeXchange | proteins | Proteomic characteristics of diabetic patients infected with SARS-CoV-2 Omicron variant | exclude | not a plasma or serum proteome of patients (cells, tissue, other fluid or viral proteins) |
| PXD067176 | ProteomeXchange | proteins | Acute Neurological Manifestations in COVIDâ€19 Are Associated With Complement Activation  | exclude | not a plasma or serum proteome of patients (cells, tissue, other fluid or viral proteins) |
| PXD067986 | ProteomeXchange | proteins | Comparative proteomics of human naso-oropharyngeal swabs from COVID19 patients for identif | exclude | not blood |
| PXD068257 | ProteomeXchange | proteins | Antigen-specific IgG Glycosylation Profiles in Hamsters and Macaques Following COVID-19 Va | exclude | not human (Mesocricetus auratus (Golden hamster); Macaca mulatta (Rhesus macaque)) |
| PXD069164 | ProteomeXchange | proteins | Platelet proteomics on less than a drop of previously frozen, non-citrate plasma | exclude | not a plasma or serum proteome of patients (cells, tissue, other fluid or viral proteins) |
| PXD070575 | ProteomeXchange | proteins | Proteomics of A549 and A549 cells transducted with ORF7a SARS-CoV-2 accessory protein | exclude | in vitro or cell culture |
| PXD071515 | ProteomeXchange | proteins | Dysregulation of adhesion and immune pathways in Omicron-associated neurological complicat | exclude | not blood |
| PXD073010 | ProteomeXchange | proteins | Plasma Proteomic Profiling of Comorbid and Non-Comorbid COVID-19 Patients in ICU | exclude | not a plasma or serum proteome of patients (cells, tissue, other fluid or viral proteins) |
| PXD077914 | ProteomeXchange | proteins | Efferocytosis of apoptotic bodies drives SARS-CoV-2 infection and macrophage inflammation | exclude | in vitro or cell culture |
| PXD080937 | ProteomeXchange | proteins | Metabolic reprogramming is associated with symptomatic COVID-19: a serum proteomics | exclude | no files listed in PRIDE (not public in PRIDE, or hosted elsewhere) |
| PXD083254 | ProteomeXchange | proteins | Study of symptoms in patients with COVID-19 | exclude | no files listed in PRIDE (not public in PRIDE, or hosted elsewhere) |
| PXD083958 | ProteomeXchange | proteins | Cross-platform clinical proteomics using the Charité open standard for plasma proteomics ( | exclude | no files listed in PRIDE (not public in PRIDE, or hosted elsewhere) |
| RPXD034683 | ProteomeXchange | proteins | PXD019113-1 The Global Phosphorylation Landscape of SARS-CoV-2 Infection [Reanalysis: PXD0 | exclude | not human (cell culture) |
| RPXD034684 | ProteomeXchange | proteins | PXD019113-2 The Global Phosphorylation Landscape of SARS-CoV-2 Infection [Reanalysis: PXD0 | exclude | not human (cell culture) |
| RPXD034685 | ProteomeXchange | proteins | PXD019645-1 Data, reagents, assays and merits of proteomics for SARS-CoV-2 research and te | exclude | not human (cell culture) |
| RPXD034686 | ProteomeXchange | proteins | PXD019645-2 Data, reagents, assays and merits of proteomics for SARS-CoV-2 research and te | exclude | not human (cell culture) |
| RPXD034687 | ProteomeXchange | proteins | PXD019645-3 Data, reagents, assays and merits of proteomics for SARS-CoV-2 research and te | exclude | not human (cell culture) |
| RPXD034688 | ProteomeXchange | proteins | PXD019645-4 Data, reagents, assays and merits of proteomics for SARS-CoV-2 research and te | exclude | not human (cell culture) |
| RPXD034689 | ProteomeXchange | proteins | PXD019645-5 Data, reagents, assays and merits of proteomics for SARS-CoV-2 research and te | exclude | not human (cell culture) |
| RPXD034690 | ProteomeXchange | proteins | PXD019645-6 Data, reagents, assays and merits of proteomics for SARS-CoV-2 research and te | exclude | not human (cell culture) |
| RPXD034691 | ProteomeXchange | proteins | PXD019645-7 Data, reagents, assays and merits of proteomics for SARS-CoV-2 research and te | exclude | not human (cell culture) |
| RPXD034692 | ProteomeXchange | proteins | PXD019423-1 MS analysis of SARS-CoV2 proteins from patient samples [Reanalysis: PXD019423] | exclude | not human (Body Fluid or Substance) |
| RPXD034693 | ProteomeXchange | proteins | PXD018804-1 Extensive proteomic dataset of Vero E6 cells infected by Italy-INMI1 SARS-CoV- | exclude | not human (cell culture) |
| RPXD034694 | ProteomeXchange | proteins | PXD018594-1 Shotgun proteomics of Vero E6 cells infected by Italy-INMI1 SARS-CoV-2 virus [ | exclude | not human (cell culture) |
| RPXD034695 | ProteomeXchange | proteins | PXD018357-1 Inhibition of growth factor signaling prevents SARS-CoV-2 replication [Reanaly | exclude | not human (cell culture) |
| RPXD034696 | ProteomeXchange | proteins | PXD018357-2 Inhibition of growth factor signaling prevents SARS-CoV-2 replication [Reanaly | exclude | not human (cell culture) |
| RPXD034697 | ProteomeXchange | proteins | PXD018117-1 A SARS-CoV-2-Human Protein-Protein Interaction Map Reveals Drug Targets and Po | exclude | not human (cell culture) |
| zenodo.10065067 | Zenodo | proteins | Effects of high dietary zinc supplementation timing on the biological responses of gestati | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| zenodo.10076491 | Zenodo | proteins | Mass spectrometry data for: A small protein coded within the mitochondrial canonical gene  | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| zenodo.10161903 | Zenodo | proteins | A cohort-based study of host gene expression: tumor suppressor and innate immune/inflammat | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| zenodo.10223608 | Zenodo | proteins | Binding and sequestration of poison frog alkaloids by a plasma globulin | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| zenodo.10456585 | Zenodo | proteins | Reactivating PTEN to impair glioma stem cells by inhibiting cytosolic iron-sulfur assembly | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| zenodo.10499011 | Zenodo | proteins | Synaptopodin regulates denervation-induced plasticity at hippocampal mossy fiber synapses | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| zenodo.10570139 | Zenodo | proteins | E2P4 COVID19 Trial Cytokine profile data | exclude | not a plasma or serum proteome of patients (cells, tissue, other fluid or viral proteins) |
| zenodo.10609655 | Zenodo | proteins | Dataset related to the article: HCG18, LEF1AS1 and lncCEACAM21 as biomarkers of disease se | exclude | single cell data |
| zenodo.10653443 | Zenodo | proteins | Supplementary Data to "Acute and short-term fluctuations in gravity are associated with ch | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| zenodo.10784727 | Zenodo | proteins | Genome-wide analysis of the effects of siRNA-mediated removal of components of ESCRT-I or  | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| zenodo.10882661 | Zenodo | proteins | The bone marrow endothelial progenitor cell response to septic infection: Dataset | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| zenodo.10885474 | Zenodo | proteins | PI3Kg inhibition circumvents inflammation and mortality in SARS-CoV-2 and other infections | exclude | not a plasma or serum proteome of patients (cells, tissue, other fluid or viral proteins) |
| zenodo.11167377 | Zenodo | proteins | Evidence of molecular mimicry in multisystem inflammatory syndrome in children (MIS-C) | exclude | not adults |
| zenodo.12110542 | Zenodo | proteins | Cross-species analysis identifies mitochondrial dysregulation as a functional consequence  | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| zenodo.14241784 | Zenodo | proteins | Replication Data for: Precision Oncology, Signaling Pathways Reprogramming and Targeted Th | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| zenodo.14680520 | Zenodo | proteins | APNet, an explainable sparse deep learning model to discover differentially active drivers | exclude | single cell data |
| zenodo.14719593 | Zenodo | proteins | Comprehensive Multimodal Immune Response Dataset for LAIV Vaccination in Pediatric Cohorts | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| zenodo.14835054 | Zenodo | proteins | Serum proteomics reveals high-affinity and convergent antibodies by tracking SARS-COV-2 hy | exclude | single cell data |
| zenodo.15050130 | Zenodo | proteins | Blood metabolome profiling for patient stratification and assessment of disease severity a | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| zenodo.15199328 | Zenodo | proteins | SARS-CoV-2 plasma antigen, plasma proteomics, and CyTOF | exclude | post-acute only (no acute samples by its record) |
| zenodo.15237087 | Zenodo | proteins | Bansal_Olechnowicz dataset | exclude | post-acute only (no acute samples by its record) |
| zenodo.15701829 | Zenodo | proteins | MAMSI Supplementary Dataset | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| zenodo.15883238 | Zenodo | proteins | Clinical outcomes and proteomics data of IVIG-treated patients | exclude | IVIG-treated COVID-19 patients only: no contrast |
| zenodo.17202421 | Zenodo | proteins | Plasma Global Metabolomics for dengue and malaria Part 2 | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| zenodo.17341280 | Zenodo | proteins | Lake Sturgeon Genome and Annotation | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| zenodo.17551035 | Zenodo | proteins | Host proteins associated with strong neutralizing SARS-CoV-2 antibody responses in a South | exclude | not a plasma or serum proteome of patients (cells, tissue, other fluid or viral proteins) |
| zenodo.17654510 | Zenodo | proteins | Systems Immunology Dataset: Innate Predictors of TSCM Response to Influenza Vaccination (O | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| zenodo.18409437 | Zenodo | proteins | Data Repository for Beat AML Drug Combo | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| zenodo.18805277 | Zenodo | proteins | Raw values for plasma SOMAScan data of CALERIE participants at baseline and 2 years after  | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| zenodo.19124653 | Zenodo | proteins | A NOVEL RNA INTERFERENCE THERAPY TO REDUCE CARDIOVASCULAR RISK BY LOWERING LIPOPROTEIN(A)  | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| zenodo.21283519 | Zenodo | proteins | Small extracellular vesicles proteome reveals persistent inflammatory and coagulopathic dy | exclude | other assay or analyte |
| zenodo.21387173 | Zenodo | proteins | Dataset: Does pathologically seeded TDP-43 exit the retina via anterograde axonal transpor | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| zenodo.21534848 | Zenodo | proteins | Dataset: Parkinson's Disease versus Alzheimer's Disease: Shared and Exclusive Pathways - P | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| zenodo.21785855 | Zenodo | proteins | Enhancing metabolite coverage using dedicated mobile phases for individual polarity modes  | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| zenodo.21789566 | Zenodo | proteins | Suppression of Ribosome Biogenesis Underlies Aging-Like Symptoms Following COVID-19: Multi | exclude | not a plasma or serum proteome of patients (cells, tissue, other fluid or viral proteins) |
| zenodo.22271531 | Zenodo | proteins | Research data on Newcastle disease virus surveillance in backyard poultry in São Paulo Sta | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| zenodo.24010209 | Zenodo | proteins | PROGNOSTIC SIGNIFICANCE OF THE ALANINE AMINOTRANSFERASE TO ASPARTATE AMINOTRANSFERASE RATI | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| zenodo.3743880 | Zenodo | proteins | Proteomics of S and N proteins | exclude | not a plasma or serum proteome of patients (cells, tissue, other fluid or viral proteins) |
| zenodo.3895886 | Zenodo | proteins | hamelpatel/COVID19_proteomic: Proteomic blood profiling in mild, severe and critical COVID | exclude | analysis scripts only; no data deposited |
| zenodo.3964503 | Zenodo | proteins | Protein Subcellular localization prediction data used in the article entitled "MSclassifie | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| zenodo.4110560 | Zenodo | proteins | The spatial landscape of lung pathology during COVID-19 progression - raw IMC data | exclude | single cell data |
| zenodo.4304890 | Zenodo | proteins | Sphaerospora molnari (Myxozoa) SMBS transcriptome data | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| zenodo.437469 | Zenodo | proteins | The organophosphate pesticide methamidophos opens the blood-testis barrier and covalently  | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| zenodo.4452813 | Zenodo | proteins | Data from: Acute necrotizing encephalopathy with SARS-CoV-2 RNA confirmed in Cerebrospinal | exclude | not a plasma or serum proteome of patients (cells, tissue, other fluid or viral proteins) |
| zenodo.4533702 | Zenodo | proteins | Data set from Vianello E, Dozio E, Bandera F, Froldi M, Micaglio E, Lamont J, Tacchini L,  | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| zenodo.4716063 | Zenodo | proteins | Two lipidomics data sources for the publication "Multi-omics profiling of living human pan | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| zenodo.4951201 | Zenodo | proteins | Data from: Metatranscriptomics yields new genomic resources and sensitive detection of inf | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| zenodo.4969242 | Zenodo | proteins | Data from: Untargeted metabolomics profiles delineate metabolic alterations in mouse plasm | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| zenodo.4977081 | Zenodo | proteins | Data from: Large scale plasma metabolome analysis reveals alterations in HDL metabolism in | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| zenodo.5012436 | Zenodo | proteins | Data from: The role of plasma and urine metabolomics in identifying new biomarkers in seve | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| zenodo.5012796 | Zenodo | proteins | The utility of reptile blood transcriptomes in molecular ecology | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| zenodo.5031173 | Zenodo | proteins | Data from: Molecular physiology of chemical defenses in a poison frog | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| zenodo.5476120 | Zenodo | proteins | Identification and evaluation of serum protein biomarkers which differentiate psoriatic fr | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| zenodo.5642579 | Zenodo | proteins | Proteomic and metabolomic profiling of urine uncovers immune responses in COVID-19 | exclude | vaccination study |
| zenodo.5771937 | Zenodo | proteins | Complement activation induces excessive T cell cytotoxicity in severe COVID-19: Analysis o | exclude | single cell data |
| zenodo.5844957 | Zenodo | proteins | ArtemisPapadaki/BAL_proteomics: | exclude | post-acute only (no acute samples by its record) |
| zenodo.5873701 | Zenodo | proteins | Dataset related to article "A 'Multiomic' Approach of Saliva Metabolomics, Microbiota, and | exclude | not a plasma or serum proteome of patients (cells, tissue, other fluid or viral proteins) |
| zenodo.6497251 | Zenodo | proteins | Multi-omics identify LRRC15 as a COVID-19 severity predictor and persistent pro-thrombotic | exclude | not blood |
| zenodo.6911339 | Zenodo | proteins | Mass spectrometry analysis of DENV-2 Infected HeLa Supernatant | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| zenodo.6924738 | Zenodo | proteins | Untargeted saliva metabolomics by liquid chromatography - mass spectrometry reveals marker | exclude | not blood |
| zenodo.7045902 | Zenodo | proteins | Integrated plasma proteomic and single-cell immune signaling network signatures demarcate  | exclude | single cell data |
| zenodo.7688661 | Zenodo | proteins | Data for manuscript,  "An optimized workflow for MS-based quantitative proteomics of chall | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| zenodo.7785543 | Zenodo | proteins | CIHR-CORONABIO DATASET | exclude | not a plasma or serum proteome of patients (cells, tissue, other fluid or viral proteins) |
| zenodo.7872791 | Zenodo | proteins | Persistent serum protein signatures define an inflammatory subcategory of long COVID | exclude | not a plasma or serum proteome of patients (cells, tissue, other fluid or viral proteins) |
| zenodo.8252959 | Zenodo | proteins | Multimodal immune profiling of SARS-CoV-2 in Uganda – soluble immune mediators | exclude | 48-plex cytokine panel, not a proteomics platform; no healthy controls; 'severe' without intensive care not assigned |
| zenodo.8279780 | Zenodo | proteins | Global effects of deletions of the sitABCD, mntH, cbiMNQO and corA genes, encoding transpo | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| E-MTAB-5273 | other | proteins | Combined Transcriptome and Proteome Leukocyte's Profiling Reveals Up-Regulated Module of G | exclude | not a COVID-19 dataset (no mention of COVID-19 or SARS-CoV-2 in its record) |
| E-MTAB-10026 | ArrayExpress | unknown | Deciphering the molecular immune response to COVID-19 using single cell multi-omics | exclude | not an RNA, protein or metabolite dataset (by its record) |
| E-MTAB-10129 | ArrayExpress | unknown | Multiple Early Factors Anticipate Post-Acute COVID-19 Sequelae | exclude | not an RNA, protein or metabolite dataset (by its record) |
| E-MTAB-10740 | ArrayExpress | unknown | Linear epitope-specific IgG responses to SARS-CoV-2 and related human coronaviruses in hum | exclude | not an RNA, protein or metabolite dataset (by its record) |
| E-MTAB-10970 | ArrayExpress | unknown | SnoRNAs and miRNAs expression profiles underlying COVID-19 disease severity | exclude | not an RNA, protein or metabolite dataset (by its record) |
| E-MTAB-11011 | ArrayExpress | unknown | Single Cell Analysis of B cells in COVID-19 comprising active and recovered disease | exclude | not an RNA, protein or metabolite dataset (by its record) |
| E-MTAB-11084 | ArrayExpress | unknown | Pre-existing cross-reactive antibody responses do not significantly impact inactivated COV | exclude | not an RNA, protein or metabolite dataset (by its record) |
| E-MTAB-11117 | ArrayExpress | unknown | Pre-existing cross-reactive antibody responses do not significantly impact inactivated COV | exclude | not an RNA, protein or metabolite dataset (by its record) |
| E-MTAB-11229 | ArrayExpress | unknown | TCR repertoire profiling of CMV specific and SARS-CoV2 Spike specific CD8+ T cells using s | exclude | not an RNA, protein or metabolite dataset (by its record) |
| E-MTAB-11240 | ArrayExpress | unknown | Trancriptome analysis of peripheral leukocytes from SARS-CoV-2 infected patients and healt | exclude | not an RNA, protein or metabolite dataset (by its record) |
| E-MTAB-11523 | ArrayExpress | unknown | Coronavirus antibody profiling in convalescent sera of RT-PCR SARS-CoV-2 confirmed hospita | exclude | not an RNA, protein or metabolite dataset (by its record) |
| E-MTAB-12161 | ArrayExpress | unknown | Longitudinal clonal IgG cross-reactivity profiles of coronavirus-specific B-cells in COVID | exclude | not an RNA, protein or metabolite dataset (by its record) |
| E-MTAB-12392 | ArrayExpress | unknown | Single cell RNA, BCR and CITE-Seq of longitudinal COVID-19 patients | exclude | not an RNA, protein or metabolite dataset (by its record) |
| E-MTAB-12505 | ArrayExpress | unknown | CRISPR screen of human embryonic stem cell line h-pES10 for the identification of host fac | exclude | not an RNA, protein or metabolite dataset (by its record) |
| E-MTAB-12779 | ArrayExpress | unknown | Discovery and systematic assessment of early biomarkers that predict progression to severe | exclude | not an RNA, protein or metabolite dataset (by its record) |
| E-MTAB-13040 | ArrayExpress | unknown | B cell epitope mapping of SARS-CoV-2 spike vaccinated rabbits | exclude | not an RNA, protein or metabolite dataset (by its record) |
| E-MTAB-13321 | ArrayExpress | unknown | single cell RNAseq of peripheral CD4+ T cells in the context SARS-CoV-2 infection | exclude | not an RNA, protein or metabolite dataset (by its record) |
| E-MTAB-13745 | ArrayExpress | unknown | A distinct effector B cell population drives autoantibody production in SARS-CoV-2 infecti | exclude | not an RNA, protein or metabolite dataset (by its record) |
| E-MTAB-13809 | ArrayExpress | unknown | Comparison between healthy individuals and Alzheimer's disease patients with covid | exclude | not an RNA, protein or metabolite dataset (by its record) |
| E-MTAB-13916 | ArrayExpress | unknown | Targeting ATP2B1 impairs PI3K/Akt/Forkhead box (Fox) transcription factor signaling and re | exclude | not an RNA, protein or metabolite dataset (by its record) |
| E-MTAB-14002 | ArrayExpress | unknown | APMAT analysis reveals associations between antigen-TCR sequences and CD8 T cell phenotype | exclude | not an RNA, protein or metabolite dataset (by its record) |
| E-MTAB-14933 | ArrayExpress | unknown | SmartSeq single-cell RNAseq of SARS-CoV-2 spike-specific CD4+ T cells | exclude | not an RNA, protein or metabolite dataset (by its record) |
| E-MTAB-15201 | ArrayExpress | unknown | Linear epitope-specific IgG pan-coronavirus responses in ferrets pre- and post ancestral S | exclude | not an RNA, protein or metabolite dataset (by its record) |
| E-MTAB-15470 | ArrayExpress | unknown | Polymorphisms associated with post-COVID-19 hyperglycemia using the Axiom Human Genotyping | exclude | not an RNA, protein or metabolite dataset (by its record) |
| E-MTAB-15820 | ArrayExpress | unknown | Clonal hematopoiesis activates pro-calcific pathways in macrophages and promotes aortic va | exclude | not an RNA, protein or metabolite dataset (by its record) |
| E-MTAB-16106 | ArrayExpress | unknown | HTGTS-Rep-Seq for Humanized Antibodies that Broadly Neutralize Omicron Sub-variants | exclude | not an RNA, protein or metabolite dataset (by its record) |
| E-MTAB-16472 | ArrayExpress | unknown | Atherosclerosis licenses for an exceeding immune response in COVID-19 disease by interfero | exclude | not an RNA, protein or metabolite dataset (by its record) |
| E-MTAB-16855 | ArrayExpress | unknown | Oncolytic Reovirus mediates innate-driven SARS-CoV-2 elimination in the absence of cell to | exclude | not an RNA, protein or metabolite dataset (by its record) |
| E-MTAB-9357 | ArrayExpress | unknown | Multi-omic analysis resolves coordinated immune responses to mild, moderate and severe COV | exclude | not an RNA, protein or metabolite dataset (by its record) |
| E-MTAB-9638 | ArrayExpress | unknown | Functional genomic screens in human cells to identify host factors for SARS-CoV-2 and comm | exclude | not an RNA, protein or metabolite dataset (by its record) |
| E-MTAB-9995 | ArrayExpress | unknown | Longitudinal scRNAseq profiling of memory B cells in severe COVID19 patients | exclude | not an RNA, protein or metabolite dataset (by its record) |
| PRJEB38339 | BioProject | unknown | Next Generation Sequencing of T and B cell receptor repertoires from COVID-19 patients sho | exclude | not an RNA, protein or metabolite dataset (by its record) |
| PRJEB70840 | BioProject | unknown | TCR Rep-Seq of COVID-19 affected individuals | exclude | not an RNA, protein or metabolite dataset (by its record) |
| PRJEB82466 | BioProject | unknown | Peripheral Blood Mononuclear Cell Electrophysiology as a Marker for Immune Response to Exp | exclude | not an RNA, protein or metabolite dataset (by its record) |
| PRJNA1005648 | BioProject | unknown | Circulating miRNA Profiles in COVID-19 Patients and Meta-Analysis: Implications for Diseas | exclude | not an RNA, protein or metabolite dataset (by its record) |
| PRJNA1014436 | BioProject | unknown | Stereotypic persistent B cell receptor clonotypes in Alzheimer’s Disease | exclude | not an RNA, protein or metabolite dataset (by its record) |
| PRJNA1015194 | BioProject | unknown | SNP (A>G - rs13057211) but not GT(n) polymorphism in HMOX-1 promoter gene is associated wi | exclude | not an RNA, protein or metabolite dataset (by its record) |
| PRJNA1021596 | BioProject | unknown | Single-cell sequencing of 12 PBMC samples to study the different immune response after SAR | exclude | duplicate: BioProject of GSE244173 |
| PRJNA1032283 | BioProject | unknown | Impact of SARS-CoV-2 vaccination on FcgammaRIIIA/CD16 dynamcis in Natural Killer cells: re | exclude | not an RNA, protein or metabolite dataset (by its record) |
| PRJNA1040842 | BioProject | unknown | An epigenome-wide DNA methylation study conducted on 96 patients six months after COVID-19 | exclude | not an RNA, protein or metabolite dataset (by its record) |
| PRJNA1046698 | BioProject | unknown | The impact of primary SARS-CoV-2 infection on human anelloviruses | exclude | not an RNA, protein or metabolite dataset (by its record) |
| PRJNA1047084 | BioProject | unknown | Altered expression of PBMC miRNA profile following the third dose of inactivated COVID-19  | exclude | not an RNA, protein or metabolite dataset (by its record) |
| PRJNA1051817 | BioProject | unknown | Analysis of KIR and HLA polymorphisms in the Chinese individuals with COVID-19 | exclude | not an RNA, protein or metabolite dataset (by its record) |
| PRJNA1059130 | BioProject | unknown | Identification of host genetic factors that are determinant for the development of severe  | exclude | not an RNA, protein or metabolite dataset (by its record) |
| PRJNA1059131 | BioProject | unknown | Identification of host genetic factors that are determinant for the development of severe  | exclude | not an RNA, protein or metabolite dataset (by its record) |
| PRJNA1071186 | BioProject | unknown | RNA expression differences in innate immune response to COVID19 with advancing age. | exclude | duplicate: BioProject of GSE254647 |
| PRJNA1081927 | BioProject | unknown | Tracking protective inflammation signatures in lungs after SARS-CoV-2 omicron BA.1 infecti | exclude | not an RNA, protein or metabolite dataset (by its record) |
| PRJNA1085011 | BioProject | unknown | REDS-IV-P Epidemiology, Surveillance and Preparedness of the Novel SARS-CoV-2 Epidemic (RE | exclude | not an RNA, protein or metabolite dataset (by its record) |
| PRJNA1085012 | BioProject | unknown | REDS-IV-P Epidemiology, Surveillance and Preparedness of the Novel SARS-CoV-2 Epidemic (RE | exclude | not an RNA, protein or metabolite dataset (by its record) |
| PRJNA1099166 | BioProject | unknown | Homo sapiens Raw sequence reads | exclude | not an RNA, protein or metabolite dataset (by its record) |
| PRJNA1130853 | BioProject | unknown | Rewired type I IFN signaling is linked to age-dependent differences in COVID-19 | exclude | duplicate: BioProject of GSE271284 |
| PRJNA1138814 | BioProject | unknown | Paired analysis of host and pathogen genomes identifies determinants of human tuberculosis | exclude | not an RNA, protein or metabolite dataset (by its record) |
| PRJNA1138815 | BioProject | unknown | Paired analysis of host and pathogen genomes identifies determinants of human tuberculosis | exclude | not an RNA, protein or metabolite dataset (by its record) |
| PRJNA1144424 | BioProject | unknown | Antigen-specific T helper cells predict intensity and longevity of immune responses to SAR | exclude | duplicate: BioProject of GSE273994 |
| PRJNA1193318 | BioProject | unknown | Adaptation of PCR-based library preparation for DNBSEQ platform for cancer mutation testin | exclude | not an RNA, protein or metabolite dataset (by its record) |
| PRJNA1193547 | BioProject | unknown | A diverse, naive-biased B cell compartment is associated with vaccine responsiveness in SO | exclude | not an RNA, protein or metabolite dataset (by its record) |
| PRJNA1199552 | BioProject | unknown | Human ACE2 Transgenic Mouse Line | exclude | not an RNA, protein or metabolite dataset (by its record) |
| PRJNA1199949 | BioProject | unknown | CRISPR-LRS mapping of BAC-hACE2 | exclude | not an RNA, protein or metabolite dataset (by its record) |
| PRJNA1214709 | BioProject | unknown | Investigating the Risk of Arrhythmogenesis Associated with Fentanyl Abuse Using Mice and H | exclude | duplicate: BioProject of GSE287796 |
| PRJNA1242341 | BioProject | unknown | Associations of epigenetic aging and COVID-19: a 3-year longitudinal study | exclude | not an RNA, protein or metabolite dataset (by its record) |
| PRJNA1257454 | BioProject | unknown | Comparison of single-cell profiles in patients with Long Covid (LC) and those who recovere | exclude | not an RNA, protein or metabolite dataset (by its record) |
| PRJNA1269494 | BioProject | unknown | DNA methylation profiles of hospitalized COVID-19 patients | exclude | not an RNA, protein or metabolite dataset (by its record) |
| PRJNA1288046 | BioProject | unknown | Prevention and Early Treatment of Acute Lung Injury (PETAL) Network – Functional, imaging, | exclude | not an RNA, protein or metabolite dataset (by its record) |
| PRJNA1288047 | BioProject | unknown | Prevention and Early Treatment of Acute Lung Injury (PETAL) Network – Functional, imaging, | exclude | not an RNA, protein or metabolite dataset (by its record) |
| PRJNA1291097 | BioProject | unknown | Epigenetic liquid biopsies reveal elevated vascular endothelial cell turnover and erythrop | exclude | not an RNA, protein or metabolite dataset (by its record) |
| PRJNA1356723 | BioProject | unknown | Homo sapiens Epigenomics | exclude | not an RNA, protein or metabolite dataset (by its record) |
| PRJNA1368929 | BioProject | unknown | Persistent immune dysregulation and metabolic alterations following SARS-CoV-2 infection | exclude | not an RNA, protein or metabolite dataset (by its record) |
| PRJNA1368930 | BioProject | unknown | Persistent immune dysregulation and metabolic alterations following SARS-CoV-2 infection | exclude | not an RNA, protein or metabolite dataset (by its record) |
| PRJNA1368937 | BioProject | unknown | T cell immune responses to SARS-CoV-2 | exclude | not an RNA, protein or metabolite dataset (by its record) |
| PRJNA1368938 | BioProject | unknown | T cell immune responses to SARS-CoV-2 | exclude | not an RNA, protein or metabolite dataset (by its record) |
| PRJNA1381519 | BioProject | unknown | Blood PAXgene Hospitalized and Lethal COVID-19 relative to controls (vaccinated and unvacc | exclude | not an RNA, protein or metabolite dataset (by its record) |
| PRJNA1417731 | BioProject | unknown | Analysis of the Relationship Between Naive T Cells and HLA Composition in Peripheral Blood | exclude | not an RNA, protein or metabolite dataset (by its record) |
| PRJNA1465773 | BioProject | unknown | The CRAC channel inhibitor Auxora interrupts inflammatory circuits between alveolar macrop | exclude | duplicate: BioProject of GSE330843 |
| PRJNA1499456 | BioProject | unknown | Single-cell whole-genome DNA methylation profiling of whole blood from PASC patients and h | exclude | not an RNA, protein or metabolite dataset (by its record) |
| PRJNA234047 | BioProject | unknown | Clinical validation of next-generation sequencing approaches for pathogen diagnosis | exclude | not an RNA, protein or metabolite dataset (by its record) |
| PRJNA252396 | BioProject | unknown | Eucalyptus grandis Brasuz1 X M43D1 genome sequencing | exclude | not human (Eucalyptus grandis) |
| PRJNA306908 | BioProject | unknown | Comprehensive Omics Analysis of Pediatric Tumors and Establishment of a Repository for Rel | exclude | not an RNA, protein or metabolite dataset (by its record) |
| PRJNA341405 | BioProject | unknown | Transcription profiling of PfSPZ Malaria Vaccine trial subjects in a malaria endemic regio | exclude | not an RNA, protein or metabolite dataset (by its record) |
| PRJNA352062 | BioProject | unknown | A blood RNA signature for predicting the treatment outcome in the Tuberculosis Treatment R | exclude | duplicate: BioProject of GSE89403 |
| PRJNA369684 | BioProject | unknown | A blood RNA signature for tuberculosis disease risk in household contact study - GC6 cohor | exclude | duplicate: BioProject of GSE94438 |
| PRJNA397222 | BioProject | unknown | Gene responses to RTS,S malaria vaccination in controlled human malaria study | exclude | not an RNA, protein or metabolite dataset (by its record) |
| PRJNA459725 | BioProject | unknown | Targeted sequencing based maternal whole blood expression changes with gestational age and | exclude | not an RNA, protein or metabolite dataset (by its record) |
| PRJNA494572 | BioProject | unknown | Sequencing the Peripheral Blood B and T cell Repertoire - Quantifying robustness and limit | exclude | not an RNA, protein or metabolite dataset (by its record) |
| PRJNA526839 | BioProject | unknown | Whole blood RNAseq following live attenuated influenza vaccine in Gambian children | exclude | not an RNA, protein or metabolite dataset (by its record) |
| PRJNA552286 | BioProject | unknown | Next generation sequencing for temporal analysis at high altitude in Indian and  Kyrgyz he | exclude | not an RNA, protein or metabolite dataset (by its record) |
| PRJNA560990 | BioProject | unknown | SINGLE CELL ANALYSIS OF THE HUMAN PLACENTA | exclude | not an RNA, protein or metabolite dataset (by its record) |
| PRJNA605983 | BioProject | unknown | Severe acute respiratory syndrome coronavirus 2 Raw sequence reads | exclude | not human (Severe acute respiratory syndrome coronavirus 2) |
| PRJNA628125 | BioProject | unknown | Human B cell clonal and convergent responses to acute SARS-CoV-2 infection | exclude | not an RNA, protein or metabolite dataset (by its record) |
| PRJNA629682 | BioProject | unknown | Attenuation of SARS-CoV-2 Infection by Losartan in Human Kidney Organoids | exclude | duplicate: BioProject of GSE149687 |
| PRJNA630455 | BioProject | unknown | Comprehensive mapping of immune perturbations associated with severe COVID-19 | exclude | not an RNA, protein or metabolite dataset (by its record) |
| PRJNA636824 | BioProject | unknown | Alterations of the Gut Microbiota in Patients with COVID-19 or H1N1 Influenza | exclude | not an RNA, protein or metabolite dataset (by its record) |
| PRJNA637034 | BioProject | unknown | Gut mycobiota alterations in patients with COVID-19 and H1N1 and associations with immune  | exclude | not an RNA, protein or metabolite dataset (by its record) |
| PRJNA638224 | BioProject | unknown | BCR repertoire sequencing from COVID-19 patients | exclude | not an RNA, protein or metabolite dataset (by its record) |
| PRJNA649279 | BioProject | unknown | Capturing and Recreating Diverse Antibody Repertoires as Multivalent Recombinant Polyclona | exclude | not an RNA, protein or metabolite dataset (by its record) |
| PRJNA655740 | BioProject | unknown | Systems biological assessment of immunity to severe and mild COVID-19 infections | exclude | duplicate: BioProject of GSE155673 |
| PRJNA655952 | BioProject | unknown | Immunological Biomarkers of Fatal COVID-19: A Study of 868 Patients | exclude | duplicate: BioProject of GSE155897 |
| PRJNA656660 | BioProject | unknown | Pharyngeal microbiota sequencing | exclude | not an RNA, protein or metabolite dataset (by its record) |
| PRJNA667083 | BioProject | unknown | A customizable purification assay designed for small-RNA biomarker identification and func | exclude | duplicate: BioProject of GSE158948 |
| PRJNA670712 | BioProject | unknown | High resolution chromosome conformation capture from gene promoters at COVID-19, T1D, AS a | exclude | not an RNA, protein or metabolite dataset (by its record) |
| PRJNA671427 | BioProject | unknown | Genetic Determinants of Susceptibility to Severe COVID-19 Infection | exclude | not an RNA, protein or metabolite dataset (by its record) |
| PRJNA671428 | BioProject | unknown | Genetic Determinants of Susceptibility to Severe COVID-19 Infection | exclude | not an RNA, protein or metabolite dataset (by its record) |
| PRJNA674312 | BioProject | unknown | Primary human trophoblasts exposed to normoxia or hypoxia | exclude | not an RNA, protein or metabolite dataset (by its record) |
| PRJNA674329 | BioProject | unknown | Primary human trophoblasts exposed to normoxia or hypoxia three replicates | exclude | not an RNA, protein or metabolite dataset (by its record) |
| PRJNA674366 | BioProject | unknown | Primary human trophoblasts exposed to DMSO in normoxia | exclude | not an RNA, protein or metabolite dataset (by its record) |
| PRJNA679336 | BioProject | unknown | Longitudinal multi-omics identifies responses of megakaryocytes, erythroid cells and plasm | exclude | duplicate: BioProject of GSE161678 |
| PRJNA679920 | BioProject | unknown | Dynamics of B-cell repertoires in response to SARS-CoV-2 infection in COVID-19 patients wi | exclude | not an RNA, protein or metabolite dataset (by its record) |
| PRJNA680276 | BioProject | unknown | Epigenetic clocks are not accelerated in COVID-19 patients | exclude | not an RNA, protein or metabolite dataset (by its record) |
| PRJNA687910 | BioProject | unknown | cell-free DNA sequencing of plasma | exclude | not an RNA, protein or metabolite dataset (by its record) |
| PRJNA688684 | BioProject | unknown | Global Absence and Targeting of Protective Immune States in Severe COVID-19 | exclude | duplicate: BioProject of GSE163668 |
| PRJNA689961 | BioProject | unknown | Gut microbiome sequencing | exclude | not an RNA, protein or metabolite dataset (by its record) |
| PRJNA690513 | BioProject | unknown | Outcome of SARS-CoV-2 infection is linked to MAIT cell activation and cytotoxicity | exclude | duplicate: BioProject of GSE164386 |
| PRJNA693881 | BioProject | unknown | Viral Bronchiolitis | exclude | not an RNA, protein or metabolite dataset (by its record) |
| PRJNA699486 | BioProject | unknown | Identification of distinct miRNAs associated with SARS-CoV-2 infection | exclude | not an RNA, protein or metabolite dataset (by its record) |
| PRJNA703929 | BioProject | unknown | Host methylation predicts SARS-CoV-2 infection and clinical outcome | exclude | not an RNA, protein or metabolite dataset (by its record) |
| PRJNA714242 | BioProject | unknown | Microbiome regulation of mucosal immune responses to SARS-CoV-2 | exclude | not an RNA, protein or metabolite dataset (by its record) |
| PRJNA718009 | BioProject | unknown | Single-cell epigenomic landscape of peripheral immune cells reveals establishment of train | exclude | not an RNA, protein or metabolite dataset (by its record) |
| PRJNA718557 | BioProject | unknown | Non-severe SARS-CoV-2 infection is characterised by very early T cell proliferation indepe | exclude | not an RNA, protein or metabolite dataset (by its record) |
| PRJNA723367 | BioProject | unknown | Predicting COVID-19 Severity with a Specific Nucleocapsid Antibody plus Disease Risk Facto | exclude | duplicate: BioProject of GSE172471 |
| PRJNA732198 | BioProject | unknown | Genome-wide DNA methylation analysis of COVID-19 severity | exclude | duplicate: BioProject of GSE174818 |
| PRJNA733666 | BioProject | unknown | Identification of LZTFL1 as a candidate effector gene at a COVID-19 risk locus (ChIP-seq) | exclude | not an RNA, protein or metabolite dataset (by its record) |
| PRJNA733669 | BioProject | unknown | Identification of LZTFL1 as a candidate effector gene at a COVID-19 risk locus (ATAC-seq) | exclude | not an RNA, protein or metabolite dataset (by its record) |
| PRJNA742421 | BioProject | unknown | Up-regulation of ACE2, the SARS-CoV-2 Receptor, in Asthmatics on Maintenance Inhaled Corti | exclude | duplicate: BioProject of GSE179156 |
| PRJNA743287 | BioProject | unknown | COVID-19 DNA methylation analysis | exclude | not an RNA, protein or metabolite dataset (by its record) |
| PRJNA744567 | BioProject | unknown | Efficient discovery of potently neutralizing SARS-CoV-2 antibodies using LIBRA-seq with li | exclude | not an RNA, protein or metabolite dataset (by its record) |
| PRJNA746586 | BioProject | unknown | Characterization of B-cell receptor repertoire associated with asymptomatic and symptomati | exclude | not an RNA, protein or metabolite dataset (by its record) |
| PRJNA748731 | BioProject | unknown | People critically ill with COVID-19 exhibit peripheral immune profiles predictive of morta | exclude | duplicate: BioProject of GSE180578 |
| PRJNA748746 | BioProject | unknown | Human TLR7 and plasmacytoid dendritic cells are essential for type I IFN pulmonary immunit | exclude | duplicate: BioProject of GSE181787 |
| PRJNA750531 | BioProject | unknown | Rapid Acceleration of Diagnostics - Radical (RADx-rad): Bioinformatics Framework for Waste | exclude | not an RNA, protein or metabolite dataset (by its record) |
| PRJNA750532 | BioProject | unknown | Rapid Acceleration of Diagnostics - Radical (RADx-rad): Bioinformatics Framework for Waste | exclude | not an RNA, protein or metabolite dataset (by its record) |
| PRJNA750902 | BioProject | unknown | The salivary microbiome of COVID-19 inpatients and outpatients | exclude | not an RNA, protein or metabolite dataset (by its record) |
| PRJNA753626 | BioProject | unknown | The differential analysis about N6-methyladenosine of COVID-19 and the health | exclude | not an RNA, protein or metabolite dataset (by its record) |
| PRJNA755276 | BioProject | unknown | Analysis of circulating miRNAs altered in the acute respiratory distress syndrome related  | exclude | not an RNA, protein or metabolite dataset (by its record) |
| PRJNA755325 | BioProject | unknown | DNA methylation profiling of immune cells from tuberculosis-exposed individuals overlaps w | exclude | not an RNA, protein or metabolite dataset (by its record) |
| PRJNA762387 | BioProject | unknown | Rapid Acceleration of Diagnostics - Radical (RADx-rad): Identifying Biomarker Signatures o | exclude | not an RNA, protein or metabolite dataset (by its record) |
| PRJNA762389 | BioProject | unknown | Rapid Acceleration of Diagnostics - Radical (RADx-rad): Identifying Biomarker Signatures o | exclude | not an RNA, protein or metabolite dataset (by its record) |
| PRJNA770903 | BioProject | unknown | Monoclonal antibody-mediated neutralization of SARS-CoV-2 in an IRF9-deficient child | exclude | not an RNA, protein or metabolite dataset (by its record) |
| PRJNA780093 | BioProject | unknown | Next-Generation Sequencing of Cerebrospinal fluid from COVID-19 Patients with Neurological | exclude | not an RNA, protein or metabolite dataset (by its record) |
| PRJNA780877 | BioProject | unknown | A new tractable method for generating Human Alveolar Macrophage Like cells in vitro to stu | exclude | duplicate: BioProject of GSE188945 |
| PRJNA782883 | BioProject | unknown | SARS-CoV-2 reactive and neutralizing antibodies discovered by single cell sequencing of pl | exclude | not an RNA, protein or metabolite dataset (by its record) |
| PRJNA783883 | BioProject | unknown | Nanostring plasma microRNA profiling in hospitalised patients with COVID-19 | exclude | not an RNA, protein or metabolite dataset (by its record) |
| PRJNA785568 | BioProject | unknown | Durable Alterations of Hematopoiesis and Chromatin post-COVID-19 [Bulk ATAC-seq] | exclude | not an RNA, protein or metabolite dataset (by its record) |
| PRJNA791978 | BioProject | unknown | Clinical Trial of COVID-19 Convalescent Plasma in Outpatients (C3PO) | exclude | not an RNA, protein or metabolite dataset (by its record) |
| PRJNA791979 | BioProject | unknown | Clinical Trial of COVID-19 Convalescent Plasma in Outpatients (C3PO) | exclude | not an RNA, protein or metabolite dataset (by its record) |
| PRJNA792793 | BioProject | unknown | DNA methylation profiles in pneumonia patients reflect changes in cell types and pneumonia | exclude | not an RNA, protein or metabolite dataset (by its record) |
| PRJNA798304 | BioProject | unknown | Epigenetic Profiling Linked to Multisystem Inflammatory Syndrome in Children (MIS-C): A Mu | exclude | not an RNA, protein or metabolite dataset (by its record) |
| PRJNA799434 | BioProject | unknown | Rapid Acceleration of Diagnostics - Radical (RADx-rad): Microfluidic Isolation and Charact | exclude | not an RNA, protein or metabolite dataset (by its record) |
| PRJNA799435 | BioProject | unknown | Rapid Acceleration of Diagnostics - Radical (RADx-rad): Microfluidic Isolation and Charact | exclude | not an RNA, protein or metabolite dataset (by its record) |
| PRJNA802573 | BioProject | unknown | microRNA profiling in the plasma of hospitalized COVID-19 patients | exclude | not an RNA, protein or metabolite dataset (by its record) |
| PRJNA808247 | BioProject | unknown | Durable Alterations of Hematopoiesis and Chromatin post-COVID-19 [snATAC-seq] | exclude | not an RNA, protein or metabolite dataset (by its record) |
| PRJNA808249 | BioProject | unknown | Durable Alterations of Hematopoiesis and Chromatin post-COVID-19 [MultiomeATAC-seq] | exclude | not an RNA, protein or metabolite dataset (by its record) |
| PRJNA809228 | BioProject | unknown | Genome-wide DNA methylation analysis one year post-COVID-19 infection | exclude | not an RNA, protein or metabolite dataset (by its record) |
| PRJNA835962 | BioProject | unknown | Collaborative Cohort of Cohorts for COVID-19 Research (C4R): Genetic Epidemiology of COPD  | exclude | not an RNA, protein or metabolite dataset (by its record) |
| PRJNA835963 | BioProject | unknown | Collaborative Cohort of Cohorts for COVID-19 Research (C4R): Genetic Epidemiology of COPD  | exclude | not an RNA, protein or metabolite dataset (by its record) |
| PRJNA836688 | BioProject | unknown | Collaborative Cohort of Cohorts for COVID-19 Research (C4R): Severe Asthma Research Progra | exclude | not an RNA, protein or metabolite dataset (by its record) |
| PRJNA836689 | BioProject | unknown | Collaborative Cohort of Cohorts for COVID-19 Research (C4R): Severe Asthma Research Progra | exclude | not an RNA, protein or metabolite dataset (by its record) |
| PRJNA839182 | BioProject | unknown | mouse gut metagenome Genome sequencing and assembly | exclude | not human (mouse gut metagenome) |
| PRJNA851448 | BioProject | unknown | IMMUNE DYSREGULATION FOLLOWING INFECTION BY THE NEW CORONAVIRUS SARS-CoV-2 (COVID-19); A N | exclude | not an RNA, protein or metabolite dataset (by its record) |
| PRJNA854732 | BioProject | unknown | Collaborative Cohort of Cohorts for COVID-19 Research (C4R): Mediators of Atherosclerosis  | exclude | not an RNA, protein or metabolite dataset (by its record) |
| PRJNA854733 | BioProject | unknown | Collaborative Cohort of Cohorts for COVID-19 Research (C4R): Mediators of Atherosclerosis  | exclude | not an RNA, protein or metabolite dataset (by its record) |
| PRJNA855896 | BioProject | unknown | Serological mechanisms driving cardiotoxicity in COVID-19 | exclude | duplicate: BioProject of GSE207473 |
| PRJNA867002 | BioProject | unknown | TFEB and TFE3 are required for efficient cellular response against beta-coronavirus [ChIP- | exclude | not an RNA, protein or metabolite dataset (by its record) |
| PRJNA869780 | BioProject | unknown | COVID-19 Low Responder Study | exclude | not an RNA, protein or metabolite dataset (by its record) |
| PRJNA870070 | BioProject | unknown | Whole Blood profiling of COVID convalescent and Healthy donors with nCounter | exclude | duplicate: BioProject of GSE211378 |
| PRJNA870094 | BioProject | unknown | TCR profiling of COVID convalescent and healthy donor whole blood with nCounter | exclude | duplicate: BioProject of GSE211394 |
| PRJNA870911 | BioProject | unknown | SARS-CoV-2 infection of Lung Organoids Reveals Conserved Use of Tetraspanin-8 by Ancestral | exclude | not an RNA, protein or metabolite dataset (by its record) |
| PRJNA877691 | BioProject | unknown | TCR beta chain sequencing of SARS-CoV-2 peptide stimulated CD4 T cells | exclude | not an RNA, protein or metabolite dataset (by its record) |
| PRJNA883706 | BioProject | unknown | MIS-C: A COVID-19-associated condition between hypoimmunity and hyperimmunity | exclude | not an RNA, protein or metabolite dataset (by its record) |
| PRJNA886870 | BioProject | unknown | Comparative evolution of SARS-CoV-2 under immune suppression in the presence and absence o | exclude | not an RNA, protein or metabolite dataset (by its record) |
| PRJNA887739 | BioProject | unknown | Bisulfite Sequencing PCR of leukocytes from whole blood of COVID-19 patients | exclude | not an RNA, protein or metabolite dataset (by its record) |
| PRJNA904484 | BioProject | unknown | Rapid Acceleration of Diagnostics - Radical (RADx-rad):  Exosome-based Non-traditional Tec | exclude | not an RNA, protein or metabolite dataset (by its record) |
| PRJNA904485 | BioProject | unknown | Rapid Acceleration of Diagnostics - Radical (RADx-rad):  Exosome-based Non-traditional Tec | exclude | not an RNA, protein or metabolite dataset (by its record) |
| PRJNA906138 | BioProject | unknown | single-nucleus ATAC-seq in PBMCs with COVID-19 | exclude | not an RNA, protein or metabolite dataset (by its record) |
| PRJNA906621 | BioProject | unknown | Methylation data to study immune response to SARS-CoV-2 virus infection in COVID-19 Health | exclude | not an RNA, protein or metabolite dataset (by its record) |
| PRJNA908301 | BioProject | unknown | Gastrointestinal symptoms of long COVID relate to the ectopic colonization of specific bac | exclude | not an RNA, protein or metabolite dataset (by its record) |
| PRJNA922682 | BioProject | unknown | Single molecule resolution of near full-length human immunoglobulin heavy chain repertoire | exclude | not an RNA, protein or metabolite dataset (by its record) |
| PRJNA925228 | BioProject | unknown | The enhanced immune response of healthy individuals after vaccination with SARS CoV-2 inac | exclude | not an RNA, protein or metabolite dataset (by its record) |
| PRJNA934454 | BioProject | unknown | Nucleic acid biomarkers of immune response and cell and tissue damage in children with COV | exclude | duplicate: BioProject of GSE225220 |
| PRJNA934458 | BioProject | unknown | Nucleic acid biomarkers of immune response and cell and tissue damage in children with COV | exclude | not an RNA, protein or metabolite dataset (by its record) |
| PRJNA934459 | BioProject | unknown | Nucleic acid biomarkers of immune response and cell and tissue damage in children with COV | exclude | duplicate: BioProject of GSE225221 |
| PRJNA934665 | BioProject | unknown | Single Cell DNA Methylation Analysis for Forensic Epigenetics | exclude | not an RNA, protein or metabolite dataset (by its record) |
| PRJNA934666 | BioProject | unknown | Single Cell DNA Methylation Analysis for Forensic Epigenetics | exclude | not an RNA, protein or metabolite dataset (by its record) |
| PRJNA939179 | BioProject | unknown | Longitudinal study of DNA methylation changes in patients with severe COVID-19 | exclude | not an RNA, protein or metabolite dataset (by its record) |
| PRJNA946749 | BioProject | unknown | Characterization of all small RNAs in and comparisons across cultured megakaryocytes and p | exclude | not an RNA, protein or metabolite dataset (by its record) |
| PRJNA957447 | BioProject | unknown | Altered plasma microRNA signature in hospitalized COVID-19 patients requiring oxygen suppo | exclude | not an RNA, protein or metabolite dataset (by its record) |
| PRJNA967329 | BioProject | unknown | Functional CVIDs phenotype clusters identified by the integration of immune parameters aft | exclude | not an RNA, protein or metabolite dataset (by its record) |
| PRJNA970759 | BioProject | unknown | IMPACT OF MATERNAL SARS-COV-2 INFECTION ON UCB | exclude | not an RNA, protein or metabolite dataset (by its record) |
| PRJNA970789 | BioProject | unknown | IMPACT OF MATERNAL SARS-COV-2 INFECTION ON VILLOUS | exclude | not an RNA, protein or metabolite dataset (by its record) |
| PRJNA986577 | BioProject | unknown | Association of inflammatory mediators with mitochondrial DNA mutations in geriatric COVID- | exclude | not an RNA, protein or metabolite dataset (by its record) |
| PRJNA988119 | BioProject | unknown | Salivary miRNA profiles in COVID-19 patients with different disease severity [Plasma] | exclude | not an RNA, protein or metabolite dataset (by its record) |
| PRJNA988125 | BioProject | unknown | Salivary miRNA profiles in COVID-19 patients with different disease severity [Saliva] | exclude | not an RNA, protein or metabolite dataset (by its record) |
| PRJNA999322 | BioProject | unknown | SARS-CoV-2 infection elucidates unique features of pregnancy-specific immunity [ADT] | exclude | not an RNA, protein or metabolite dataset (by its record) |
| PRJNA999323 | BioProject | unknown | SARS-CoV-2 infection elucidates unique features of pregnancy-specific immunity [GEX] | exclude | duplicate: BioProject of GSE239450 |
| PRJNA999324 | BioProject | unknown | SARS-CoV-2 infection elucidates unique features of pregnancy-specific immunity [VDJ] | exclude | not an RNA, protein or metabolite dataset (by its record) |
| PRJNA999592 | BioProject | unknown | PostCov | exclude | not an RNA, protein or metabolite dataset (by its record) |
| GSE100467 | GEO | unknown | A comprehensive, cell specific microRNA catalogue of human peripheral blood | exclude | not an RNA, protein or metabolite dataset (by its record) |
| GSE101702 | GEO | unknown | Host response to IAV infections in human patients | exclude | not an RNA, protein or metabolite dataset (by its record) |
| GSE102592 | GEO | unknown | Hedgehog gradient in the stroma maintains niche diversity and organ function | exclude | not an RNA, protein or metabolite dataset (by its record) |
| GSE103842 | GEO | unknown | Respiratory Syncytial Virus Genotypes and Disease Severity in Young Children Hospitalized  | exclude | not an RNA, protein or metabolite dataset (by its record) |
| GSE104600 | GEO | unknown | Early Changes In B Cell Subsets Predict Risk Of Autoimmunity Following Combination Immune  | exclude | not an RNA, protein or metabolite dataset (by its record) |
| GSE11056 | GEO | unknown | Expression data from mouse lung | exclude | not human (Mus musculus) |
| GSE112087 | GEO | unknown | Systemic Lupus Erythematosus patient blood with controls | exclude | not an RNA, protein or metabolite dataset (by its record) |
| GSE115198 | GEO | unknown | Distinct immune responses to type I and III interferons are regulated by interferon regula | exclude | not an RNA, protein or metabolite dataset (by its record) |
| GSE11903 | GEO | unknown | Effective treatment of psoriasis with etanercept is linked to suppression of IL17 signalin | exclude | not an RNA, protein or metabolite dataset (by its record) |
| GSE122692 | GEO | unknown | Immune parameters and outcome during Ebola virus disease | exclude | not an RNA, protein or metabolite dataset (by its record) |
| GSE124220 | GEO | unknown | Human Bacille Calmette-Guérin vaccination elicits trained immunity via the hematopoietic p | exclude | not an RNA, protein or metabolite dataset (by its record) |
| GSE12484 | GEO | unknown | Periodontitis is associated with a type-1 interferon signature in peripheral blood neutrop | exclude | not an RNA, protein or metabolite dataset (by its record) |
| GSE124898 | GEO | unknown | Single-cell profiling of cutaneous T-cell lymphoma after treatment with vorinostat and pho | exclude | not an RNA, protein or metabolite dataset (by its record) |
| GSE127472 | GEO | unknown | High-throughput tissue dissection and cell purification with digital cytometry | exclude | not an RNA, protein or metabolite dataset (by its record) |
| GSE128889 | GEO | unknown | Identification of a mesenchymal progenitor cell hierarchy in adipose tissue | exclude | not an RNA, protein or metabolite dataset (by its record) |
| GSE134080 | GEO | unknown | RNASeq whole blood of Dutch 500FG cohort | exclude | not an RNA, protein or metabolite dataset (by its record) |
| GSE134425 | GEO | unknown | Inflammatory Cytokines shape a Changing DNA Methylome in Monocytes Mirroring Disease Activ | exclude | not an RNA, protein or metabolite dataset (by its record) |
| GSE134809 | GEO | unknown | Single-cell analysis of Crohn’s disease lesions identifies a pathogenic cellular module as | exclude | not an RNA, protein or metabolite dataset (by its record) |
| GSE13485 | GEO | unknown | Profiling 20,077 genes in 25 vaccine recipients: U133 Plus 2.0, Custom CDF Version 9 | exclude | not an RNA, protein or metabolite dataset (by its record) |
| GSE135929 | GEO | unknown | scRNASeq of splenic CD3+ cells from immunodeficient mice injected with progenitor T-cells  | exclude | not an RNA, protein or metabolite dataset (by its record) |
| GSE136831 | GEO | unknown | IPF Cell Atlas | exclude | not an RNA, protein or metabolite dataset (by its record) |
| GSE138074 | GEO | unknown | Inflammatory cytokines and organ dysfunction associate with the aberrant DNA methylome of  | exclude | not an RNA, protein or metabolite dataset (by its record) |
| GSE139522 | GEO | unknown | Single-cell analysis of olfactory neurogenesis and differentiation in adult humans | exclude | not an RNA, protein or metabolite dataset (by its record) |
| GSE1420 | GEO | unknown | Barrett's esophagus, Barrett's-associated adenocarcinomas and normal esophageal epithelium | exclude | not an RNA, protein or metabolite dataset (by its record) |
| GSE145918 | GEO | unknown | TNF-signaling modulates neutrophil-mediated immunity at the feto-maternal interface  durin | exclude | not human (Macaca mulatta) |
| GSE147143 | GEO | unknown | Single-cell sequencing of BALF reveals immune response of COV2019 infection | exclude | not an RNA, protein or metabolite dataset (by its record) |
| GSE1488290 | GEO | unknown | Altered ACE2 and interferon landscape in the COVID-19 microenvironment correlate with the  | exclude | not an RNA, protein or metabolite dataset (by its record) |
| GSE149273 | GEO | unknown | RV infections in asthmatics increase ACE2 expression and stimulate cytokine pathways impli | exclude | not an RNA, protein or metabolite dataset (by its record) |
| GSE14949 | GEO | unknown | SARS-CoV-2 infection impairs the insulin/IGF signaling pathway in the lung, liver, adipose | exclude | not an RNA, protein or metabolite dataset (by its record) |
| GSE149497 | GEO | unknown | IFN-b treatment of primary murine hepatocytes | exclude | not human (Mus musculus) |
| GSE152254 | GEO | unknown | A feedback loop between the androgen receptor and 6-phosphogluoconate dehydrogenase (6PGD) | exclude | not an RNA, protein or metabolite dataset (by its record) |
| GSE154104 | GEO | unknown | SARS-CoV-2 infection in the lungs of human ACE2 transgenic mice causes severe inflammation | exclude | not human (Mus musculus) |
| GSE1552641 | GEO | unknown | Genomic Biomarker Heterogeneities between SARS-CoV-2 and COVID-19. | exclude | not an RNA, protein or metabolite dataset (by its record) |
| GSE155286 | GEO | unknown | SARS-CoV-2 infection is effectively treated and prevented by EIDD-2801 | exclude | not an RNA, protein or metabolite dataset (by its record) |
| GSE155518 | GEO | unknown | SARS-CoV-2 infection of primary human lung epithelium for COVID-19 modeling and drug disco | exclude | not an RNA, protein or metabolite dataset (by its record) |
| GSE155817 | GEO | unknown | Single keratinocyte analysis uncovers Foxm1 as a Yap-dependent regulator of human epiderma | exclude | not an RNA, protein or metabolite dataset (by its record) |
| GSE158712 | GEO | unknown | Rapid and efficacious humoral responses are required for recovery from hemorrhagic fever w | exclude | not an RNA, protein or metabolite dataset (by its record) |
| GSE159252 | GEO | unknown | Single-cell analyses identify circulating anti-tumor CD8 T cells and markers for their enr | exclude | not an RNA, protein or metabolite dataset (by its record) |
| GSE160435 | GEO | unknown | bulk RNAseq of SARS-COV2 infected human lung organoid | exclude | not an RNA, protein or metabolite dataset (by its record) |
| GSE16160 | GEO | unknown | Expression data of HL-60 leukemia cells after treatment with Curcumin and TCDD | exclude | not an RNA, protein or metabolite dataset (by its record) |
| GSE161678 | GEO | unknown | Longitudinal multi-omics identifies responses of megakaryocytes, erythroid cells and plasm | exclude | not an RNA, protein or metabolite dataset (by its record) |
| GSE163431 | GEO | unknown | Developmental programs in neuroblastoma | exclude | not an RNA, protein or metabolite dataset (by its record) |
| GSE163530 | GEO | unknown | A single-cell and spatial atlas of autopsy tissues reveals pathology and cellular targets  | exclude | not an RNA, protein or metabolite dataset (by its record) |
| GSE16524 | GEO | unknown | Expression data from skin fibroblasts derived from Setleis Syndrome patients and normal co | exclude | not an RNA, protein or metabolite dataset (by its record) |
| GSE165890 | GEO | unknown | Deleterious effects of SARS-CoV-2 infection on human pancreatic cells | exclude | not an RNA, protein or metabolite dataset (by its record) |
| GSE1676 | GEO | unknown | DNA damage response (NCS treatment) in several siRNA-manipulated cellular systems | exclude | not an RNA, protein or metabolite dataset (by its record) |
| GSE168739 | GEO | unknown | Epigenome-Wide Association Study of COVID-19 Severity with Respiratory Failure | exclude | not an RNA, protein or metabolite dataset (by its record) |
| GSE1719 | GEO | unknown | Macular Degeneration | exclude | not an RNA, protein or metabolite dataset (by its record) |
| GSE172471 | GEO | unknown | Predicting COVID-19 Severity with a Specific Nucleocapsid Antibody plus Disease Risk Facto | exclude | not an RNA, protein or metabolite dataset (by its record) |
| GSE173707 | GEO | unknown | ORAI1 establishes resistance to SARS-CoV-2 infection by regulating tonic type I interferon | exclude | not an RNA, protein or metabolite dataset (by its record) |
| GSE1739 | GEO | unknown | SARS | exclude | not an RNA, protein or metabolite dataset (by its record) |
| GSE174818 | GEO | unknown | Genome-wide DNA methylation analysis of COVID-19 severity | exclude | not an RNA, protein or metabolite dataset (by its record) |
| GSE1751 | GEO | unknown | Human blood expression for Huntington's disease versus control | exclude | not an RNA, protein or metabolite dataset (by its record) |
| GSE175943 | GEO | unknown | High-fat high-sugar diet-induced changes in the lipid metabolism are associated with mildl | exclude | not human (Mesocricetus auratus) |
| GSE178333 | GEO | unknown | Human organoid systems reveal in vitro correlates of fitness for SARS-CoV-2 B.1.1.7 | exclude | not an RNA, protein or metabolite dataset (by its record) |
| GSE180743 | GEO | unknown | New-Onset IgG Autoantibodies in Hospitalized Patients with COVID-19 | exclude | not an RNA, protein or metabolite dataset (by its record) |
| GSE182299 | GEO | unknown | Molecular Profiling of COVID-19 Autopsies Uncovers Novel Disease Mechanisms | exclude | not an RNA, protein or metabolite dataset (by its record) |
| GSE183904 | GEO | unknown | Single-cell atlas of lineage states, tumor microenvironment and subtype specific expressio | exclude | not an RNA, protein or metabolite dataset (by its record) |
| GSE185555 | GEO | unknown | Single-cell Architecture and Functional Requirement of Alternative Splicing during Hematop | exclude | not human (Mus musculus) |
| GSE1922 | GEO | unknown | Leukemia study / K562 cell line | exclude | not an RNA, protein or metabolite dataset (by its record) |
| GSE196830 | GEO | unknown | Single-cell eQTL mapping identifies cell type specific genetic control of autoimmune disea | exclude | not an RNA, protein or metabolite dataset (by its record) |
| GSE1987 | GEO | unknown | Non Small Cell Lung Cancer | exclude | not an RNA, protein or metabolite dataset (by its record) |
| GSE199123 | GEO | unknown | Single cell RNAseq reveal covid19 vaccine response difference | exclude | not an RNA, protein or metabolite dataset (by its record) |
| GSE2052 | GEO | unknown | IPF versus Control | exclude | not an RNA, protein or metabolite dataset (by its record) |
| GSE211394 | GEO | unknown | TCR profiling of COVID convalescent and healthy donor whole blood with nCounter | exclude | not an RNA, protein or metabolite dataset (by its record) |
| GSE213168 | GEO | unknown | RNAseq analysis of blood from influenza infected patients and healthy controls, project SI | exclude | not an RNA, protein or metabolite dataset (by its record) |
| GSE2136 | GEO | unknown | Effects of mutation in a soluble guanylyl cyclase gene | exclude | not human (Drosophila melanogaster) |
| GSE2189 | GEO | unknown | A549 teatement with MGd | exclude | not an RNA, protein or metabolite dataset (by its record) |
| GSE219015 | GEO | unknown | Deciphering cell states and genealogies of human hematopoiesis with single-cell multi-omic | exclude | not an RNA, protein or metabolite dataset (by its record) |
| GSE22147 | GEO | unknown | microRNAs-449 control vertebrate multi-ciliogenesis by repressing Notch signalling | exclude | not an RNA, protein or metabolite dataset (by its record) |
| GSE221856 | GEO | unknown | Microglia-mediated T cell Infiltration Drives Neurodegeneration in Tauopathy | exclude | not human (Mus musculus) |
| GSE222765 | GEO | unknown | Autoantibodies are highly prevalent in non-SARS-CoV-2 respiratory infections and critical  | exclude | not an RNA, protein or metabolite dataset (by its record) |
| GSE235646 | GEO | unknown | Single cell multiome profiles of bone marrow hematopoietic stem and progenitor cells from  | exclude | not an RNA, protein or metabolite dataset (by its record) |
| GSE23832 | GEO | unknown | Expression profiles of PBMC from multiple sclerosis patients | exclude | not an RNA, protein or metabolite dataset (by its record) |
| GSE249479 | GEO | unknown | Identification of a human hematopoietic stem cell subset that retains memory of inflammato | exclude | not an RNA, protein or metabolite dataset (by its record) |
| GSE253635 | GEO | unknown | Interferon-lambda uniquely promotes CD8 T cell immunity against SARS-CoV-2 relative to typ | exclude | not human (Mus musculus) |
| GSE2624 | GEO | unknown | TNF time course study | exclude | not an RNA, protein or metabolite dataset (by its record) |
| GSE274086 | GEO | unknown | Projecting development of TB-IRIS using RNA velocity of whole blood | exclude | not an RNA, protein or metabolite dataset (by its record) |
| GSE2742 | GEO | unknown | Genomic Strategies Identify the Antitumor Agent Apratoxin A as a Potent Antagonist of FGF  | exclude | not an RNA, protein or metabolite dataset (by its record) |
| GSE278673 | GEO | unknown | LSK immunophenotype hematopoietic stem cells after irradiation at day 2 | exclude | not human (Mus musculus) |
| GSE280145 | GEO | unknown | LK immunophenotype cells after sleep deprivation. | exclude | not human (Mus musculus) |
| GSE280270 | GEO | unknown | UCBCs cultured with TPO | exclude | not an RNA, protein or metabolite dataset (by its record) |
| GSE29614 | GEO | unknown | Time Course of Young Adults Vaccinated with Influenza TIV Vaccine during 2007/08 Flu Seaso | exclude | not an RNA, protein or metabolite dataset (by its record) |
| GSE29615 | GEO | unknown | Time Course of Young Adults Vaccinated with Influenza LAIV Vaccine during 2008/09 Flu Seas | exclude | not an RNA, protein or metabolite dataset (by its record) |
| GSE29617 | GEO | unknown | Time Course of Young Adults Vaccinated with Influenza TIV Vaccine during 2008/09 Flu Seaso | exclude | not an RNA, protein or metabolite dataset (by its record) |
| GSE30719 | GEO | unknown | Microarray Analysis of West Nile Virus infected Human Retinal Pigment Epithelium | exclude | not an RNA, protein or metabolite dataset (by its record) |
| GSE32138 | GEO | unknown | The response of human primary airway epithelial cells to Influenza or RSV infection [hAECs | exclude | not an RNA, protein or metabolite dataset (by its record) |
| GSE3311 | GEO | unknown | Effect of chronic ethanol consumption on rat pancreas | exclude | not human (Rattus norvegicus) |
| GSE35182 | GEO | unknown | Sex differences during  acute myocarditis against coxsackievirus B3 (CVB3)-induced myocard | exclude | not human (Mus musculus) |
| GSE37571 | GEO | unknown | Host Regulatory Network Response in Calu-3 cell Infection with A/CA/04/2009 Influenza viru | exclude | not an RNA, protein or metabolite dataset (by its record) |
| GSE3807 | GEO | unknown | CD3+ T cells derived from bone marrow or peripheral blood of aplastic anemia patients | exclude | not an RNA, protein or metabolite dataset (by its record) |
| GSE38206 | GEO | unknown | Transcription profiling of major depression peripheral blood mononuclear cells (PBMCs) at  | exclude | not an RNA, protein or metabolite dataset (by its record) |
| GSE40885 | GEO | unknown | Data expression in alveolar macrophages induced by lipopolysaccharide in humans | exclude | not an RNA, protein or metabolite dataset (by its record) |
| GSE42832 | GEO | unknown | Microarray study using whole blood and purified leucocyte populations from healthy control | exclude | not an RNA, protein or metabolite dataset (by its record) |
| GSE4302 | GEO | unknown | Genome-Wide Profiling of Airway Epithelial Cells in Asthmatics, Smokers and Healthy Contro | exclude | not an RNA, protein or metabolite dataset (by its record) |
| GSE43696 | GEO | unknown | Asthma | exclude | not an RNA, protein or metabolite dataset (by its record) |
| GSE44555 | GEO | unknown | Multiple tissue expression data from inbreds and F1 of CAST, PWK, and WSB | exclude | not human (Mus musculus) |
| GSE4748 | GEO | unknown | A cyanobacterial LPS antagonist prevents endotoxin shock and blocks sustained TLR4 stimula | exclude | not an RNA, protein or metabolite dataset (by its record) |
| GSE47960 | GEO | unknown | SHAE002: SARS-CoV, SARS-dORF6 and SARS-BatSRBD infection of HAE cultures. | exclude | not an RNA, protein or metabolite dataset (by its record) |
| GSE47961 | GEO | unknown | SHAE003: SARS-CoV, SARS-dORF6 and SARS-BatSRBD infection of HAE cultures. | exclude | not an RNA, protein or metabolite dataset (by its record) |
| GSE47962 | GEO | unknown | SHAE004: SARS-CoV, SARS-dORF6 and SARS-BatSRBD infection of HAE cultures. | exclude | not an RNA, protein or metabolite dataset (by its record) |
| GSE47963 | GEO | unknown | SARS-CoV, SARS-dORF6 and SARS-BatSRBD infection of HAE cultures. | exclude | not an RNA, protein or metabolite dataset (by its record) |
| GSE48027 | GEO | unknown | Host directed activity of Pyrazinamide in Mycobacterium tuberculosis infection | exclude | not human (Mus musculus) |
| GSE48280 | GEO | unknown | Expression data from inflammatory myopathies | exclude | not an RNA, protein or metabolite dataset (by its record) |
| GSE5509 | GEO | unknown | Expression data from Rat liver 48 hours after treated with different toxic compounds. | exclude | not human (Rattus norvegicus) |
| GSE55100 | GEO | unknown | Peripheral blood mononuclear cell in patients with type 1 diabetes compared with normal co | exclude | not an RNA, protein or metabolite dataset (by its record) |
| GSE5542 | GEO | unknown | Differential Gene Induction by Type I And Type II Interferons and Their Combination | exclude | not an RNA, protein or metabolite dataset (by its record) |
| GSE57253 | GEO | unknown | An activating NLRC4 inflammasome mutation causes autoinflammation with recurrent macrophag | exclude | not an RNA, protein or metabolite dataset (by its record) |
| GSE66099 | GEO | unknown | Unique Patients from the Genomics of Pediatric SIRS and Septic Shock Investigators (GPSSSI | exclude | not an RNA, protein or metabolite dataset (by its record) |
| GSE67198 | GEO | unknown | The effect of Sendai virus infection using 2 concentrations on U937 cells | exclude | not an RNA, protein or metabolite dataset (by its record) |
| GSE6878 | GEO | unknown | Expression data from human liver cells | exclude | not an RNA, protein or metabolite dataset (by its record) |
| GSE7111 | GEO | unknown | Resveratrol treatment of 3T3-L1 cells | exclude | not human (Mus musculus) |
| GSE72829 | GEO | unknown | Diagnosis of childhood bacterial and viral infection using host RNA expression | exclude | not an RNA, protein or metabolite dataset (by its record) |
| GSE73464 | GEO | unknown | Diagnosis of Kawasaki Disease in children using host RNA expression | exclude | not an RNA, protein or metabolite dataset (by its record) |
| GSE74235 | GEO | unknown | Hydroxychloroquine inhibits responses to group A streptococcus in peripheral blood mononuc | exclude | not an RNA, protein or metabolite dataset (by its record) |
| GSE82050 | GEO | unknown | IFI27 is a gene-expression biomarker that identifies influenza-specific immune response in | exclude | not an RNA, protein or metabolite dataset (by its record) |
| GSE84465 | GEO | unknown | Single-Cell RNAseq analysis of diffuse neoplastic infiltrating cells at the migrating fron | exclude | not an RNA, protein or metabolite dataset (by its record) |
| GSE87197 | GEO | unknown | BLUEPRINT hematopoietic progenitors | exclude | not an RNA, protein or metabolite dataset (by its record) |
| GSE89403 | GEO | unknown | A blood RNA signature for predicting the treatment outcome in the Tuberculosis Treatment R | exclude | not an RNA, protein or metabolite dataset (by its record) |
| GSE8961 | GEO | unknown | Identification of human metapneumovirus-induced gene networks in airway epithelial cells b | exclude | not an RNA, protein or metabolite dataset (by its record) |
| GSE90732 | GEO | unknown | Correlates of self-reported symptom scores in a human H1N1 influenza challenge trial | exclude | not an RNA, protein or metabolite dataset (by its record) |
| GSE94438 | GEO | unknown | A blood RNA signature for tuberculosis disease risk in household contact study - GC6 cohor | exclude | not an RNA, protein or metabolite dataset (by its record) |
| GSE9476 | GEO | unknown | Abnormal Expression Changes in AML | exclude | not an RNA, protein or metabolite dataset (by its record) |
| GSE9692 | GEO | unknown | Validation of Genome-wide Expression patterns in Pediatric Septic Shock | exclude | not an RNA, protein or metabolite dataset (by its record) |
| GSE98600 | GEO | unknown | BCG reprogramming of hematopoietic stem cells generates protective innate immunity against | exclude | not human (Mus musculus) |
| 10.17632/5xfzcb4kn8 | Mendeley Data | unknown | Whole-cortex in situ sequencing reveals peripheral input-dependent cellular and area ident | exclude | not an RNA, protein or metabolite dataset (by its record) |
| 10.17632/7686ww5z33 | Mendeley Data | unknown | Longitudinal multi-omics analyses identify responses of megakaryocytes, erythroid cells an | exclude | not an RNA, protein or metabolite dataset (by its record) |
| 10.17632/7xmbpz6p2f | Mendeley Data | unknown | Supplemental tables for "Pre-infection antiviral innate immunity contributes to sex differ | exclude | not an RNA, protein or metabolite dataset (by its record) |
| 10.17632/8kd3xjfrh4 | Mendeley Data | unknown | Medini 2021 COVID19 Supplements | exclude | not an RNA, protein or metabolite dataset (by its record) |
| 10.17632/9y93mjff4b | Mendeley Data | unknown | Characteristics of mice infected with SARS-CoV-2 prototype and Omicron BA.1 strain | exclude | not an RNA, protein or metabolite dataset (by its record) |
| 10.17632/c2vfpvbrbz | Mendeley Data | unknown | Evidence of premature lymphocyte aging in people with low anti-spike antibody levels after | exclude | not an RNA, protein or metabolite dataset (by its record) |
| 10.17632/fm7v43jrbk | Mendeley Data | unknown | Longitudinal profiling of respiratory and systemic immune response reveals myeloid cell-dr | exclude | not an RNA, protein or metabolite dataset (by its record) |
| 10.17632/r83csstphc | Mendeley Data | unknown | Critically ill COVID-19 patients exhibit peripheral immune profiles predictive of mortalit | exclude | not an RNA, protein or metabolite dataset (by its record) |
| 10.17632/sbv5rnmp8t | Mendeley Data | unknown | Clinically distinct COVID-19 cases share notably similar immune response progression: A fo | exclude | not an RNA, protein or metabolite dataset (by its record) |
| 10.17632/wkhswx3nf8 | Mendeley Data | unknown | Poly(I:C) enhances mesenchymal stem cell control of myeloid cells from COVID-19 patients. | exclude | not an RNA, protein or metabolite dataset (by its record) |
| 10.17632/zyzt62gbrw | Mendeley Data | unknown | Multi-omics and Machine Learning Analysis of Human Plasma to Identify Biomarkers in Patien | exclude | not an RNA, protein or metabolite dataset (by its record) |
| zenodo.10021757 | Zenodo | unknown | JRosenheim/SARS-CoV-2_challenge: Version 1.0 | exclude | not an RNA, protein or metabolite dataset (by its record) |
| zenodo.10038121 | Zenodo | unknown | Dynamic evolution of locomotor performance independent of changes in extended phenotype us | exclude | not an RNA, protein or metabolite dataset (by its record) |
| zenodo.10064877 | Zenodo | unknown | Divergent molecular networks program functionally distinct CD8+ skin-resident memory T cel | exclude | not an RNA, protein or metabolite dataset (by its record) |
| zenodo.10215356 | Zenodo | unknown | Interspecies chimerism with human embryonic stem cells generates functional human dopamine | exclude | not an RNA, protein or metabolite dataset (by its record) |
| zenodo.10432060 | Zenodo | unknown | An evolution-based framework for describing human gut bacteria | exclude | not an RNA, protein or metabolite dataset (by its record) |
| zenodo.10671954 | Zenodo | unknown | Low-dose biliatresone treatment of pregnant mice causes subclinical biliary disease in the | exclude | not an RNA, protein or metabolite dataset (by its record) |
| zenodo.10797098 | Zenodo | unknown | Chelysheva/COVID_multiomics_codes: v1.0 | exclude | not an RNA, protein or metabolite dataset (by its record) |
| zenodo.10817195 | Zenodo | unknown | Data Repository: LRRK2 G2019S mutation affects human iPSC-derived endothelial cells (line  | exclude | not an RNA, protein or metabolite dataset (by its record) |
| zenodo.10837057 | Zenodo | unknown | Data Repository: LRRK2 G2019S mutation affects human iPSC-derived endothelial cells (line  | exclude | not an RNA, protein or metabolite dataset (by its record) |
| zenodo.10914402 | Zenodo | unknown | Biogeography and evolution of social parasitism in Australian Myrmecia bulldog ants reveal | exclude | not an RNA, protein or metabolite dataset (by its record) |
| zenodo.10928039 | Zenodo | unknown | Widespread gene-environment interactions shape the immune response to SARS-CoV-2 infection | exclude | not an RNA, protein or metabolite dataset (by its record) |
| zenodo.10959358 | Zenodo | unknown | Transcobalamin receptor antibodies in autoimmune vitamin B12 central deficiency | exclude | not an RNA, protein or metabolite dataset (by its record) |
| zenodo.11154013 | Zenodo | unknown | From hidden hearing loss to supranormal auditory processing by neurotrophin 3-mediated mod | exclude | not an RNA, protein or metabolite dataset (by its record) |
| zenodo.11228952 | Zenodo | unknown | Data from: Puerarin alleviates atherosclerosis via the inhibition of Prevotella copri and  | exclude | not an RNA, protein or metabolite dataset (by its record) |
| zenodo.11446124 | Zenodo | unknown | Genomic insights into the critically endangered King Island scrubtit | exclude | not an RNA, protein or metabolite dataset (by its record) |
| zenodo.12205202 | Zenodo | unknown | Identification of age-related CpG sites from longitudinal avian methylomes | exclude | not an RNA, protein or metabolite dataset (by its record) |
| zenodo.12528938 | Zenodo | unknown | Software of "Blood DNA virome associates with autoimmune diseases and COVID-19" | exclude | not an RNA, protein or metabolite dataset (by its record) |
| zenodo.12647132 | Zenodo | unknown | Data from: Museomics help resolving the phylogeny of snowfinches (Aves, Passeridae, Montif | exclude | not an RNA, protein or metabolite dataset (by its record) |
| zenodo.1284334 | Zenodo | unknown | Joint Longitudinal And Time-To-Event Models Via Stan | exclude | not an RNA, protein or metabolite dataset (by its record) |
| zenodo.13896251 | Zenodo | unknown | jishnu-lab/ImmModCoV: v1 | exclude | not an RNA, protein or metabolite dataset (by its record) |
| zenodo.13953637 | Zenodo | unknown | Scripts for "Blood DNA virome associates with autoimmune diseases and COVID-19" | exclude | not an RNA, protein or metabolite dataset (by its record) |
| zenodo.14015194 | Zenodo | unknown | INTERVAL-RNAseq/manuscript-scripts: v1.0.0 | exclude | not an RNA, protein or metabolite dataset (by its record) |
| zenodo.14046594 | Zenodo | unknown | BALLESTARLAB/CVID-COVID: COVIDCVID | exclude | not an RNA, protein or metabolite dataset (by its record) |
| zenodo.14174632 | Zenodo | unknown | Scopeplus open architecture source code | exclude | not an RNA, protein or metabolite dataset (by its record) |
| zenodo.14447221 | Zenodo | unknown | Reproducible codes for Whole genome transcript intensive profiling from self-collected fin | exclude | not an RNA, protein or metabolite dataset (by its record) |
| zenodo.14846551 | Zenodo | unknown | Inflammatory citokines in tears and plasma of patients with type 2 diabetes and healthy co | exclude | not an RNA, protein or metabolite dataset (by its record) |
| zenodo.14991132 | Zenodo | unknown | tatsuhikonaito/OASIS_HLATCR: v1.0.0 | exclude | not an RNA, protein or metabolite dataset (by its record) |
| zenodo.15025566 | Zenodo | unknown | A single-cell atlas characterizes dysregulation of the bone marrow immune microenvironment | exclude | not an RNA, protein or metabolite dataset (by its record) |
| zenodo.15349126 | Zenodo | unknown | Data from: "Comparison of enzymatic- and bisulfite conversion of circulating cell-free tum | exclude | not an RNA, protein or metabolite dataset (by its record) |
| zenodo.15358026 | Zenodo | unknown | scCOVID-19/COVID_Imperial: Cell Genomics Publication | exclude | not an RNA, protein or metabolite dataset (by its record) |
| zenodo.15595324 | Zenodo | unknown | ZadorLaboratory/colormycells: Initial release | exclude | not an RNA, protein or metabolite dataset (by its record) |
| zenodo.15877644 | Zenodo | unknown | REdahiro/OASIS_project: v1.1 | exclude | not an RNA, protein or metabolite dataset (by its record) |
| zenodo.16986111 | Zenodo | unknown | Integrating Developmental Traits, Morphometrics, and Machine Learning for Sex Classificati | exclude | not an RNA, protein or metabolite dataset (by its record) |
| zenodo.17162225 | Zenodo | unknown | scPortrait Golgi Morphology Benchmark – Raw Files | exclude | not an RNA, protein or metabolite dataset (by its record) |
| zenodo.17527930 | Zenodo | unknown | Molecular signature of COVID-19 prior to its exacerbation by multi-omics survey | exclude | not an RNA, protein or metabolite dataset (by its record) |
| zenodo.18265276 | Zenodo | unknown | scTELL TE Accessibility Matrices: Single-cell ATAC-seq Transposable Element Analysis Data | exclude | not an RNA, protein or metabolite dataset (by its record) |
| zenodo.19188548 | Zenodo | unknown | MICROCOSM in vivo study data | exclude | not an RNA, protein or metabolite dataset (by its record) |
| zenodo.20514089 | Zenodo | unknown | SARS-CoV-2 infection induced alterations in ADAR editing patterns differ between patients  | exclude | not an RNA, protein or metabolite dataset (by its record) |
| zenodo.20734147 | Zenodo | unknown | Efficient globin production during terminal erythropoiesis depends on the cooperative acti | exclude | not an RNA, protein or metabolite dataset (by its record) |
| zenodo.21061789 | Zenodo | unknown | From viral pressure to private oligoclonality: a TCRβ repertoire sequencing study of post‑ | exclude | not an RNA, protein or metabolite dataset (by its record) |
| zenodo.21283425 | Zenodo | unknown | Catalytic P-O Bond Hydrolysis: The Mechanism Underlying Most Heavy Metal Poisoning? | exclude | not an RNA, protein or metabolite dataset (by its record) |
| zenodo.2590196 | Zenodo | unknown | Hoohm/CITE-seq-Count: 1.4.2 | exclude | not an RNA, protein or metabolite dataset (by its record) |
| zenodo.3254438 | Zenodo | unknown | Transferrin gene ChIP-seq data | exclude | not an RNA, protein or metabolite dataset (by its record) |
| zenodo.3629446 | Zenodo | unknown | mwaskom/seaborn: v0.10.0 (January 2020) | exclude | not an RNA, protein or metabolite dataset (by its record) |
| zenodo.3893991 | Zenodo | unknown | immunomind/immunarch: 0.6.5: Basic single-cell support | exclude | not an RNA, protein or metabolite dataset (by its record) |
| zenodo.4030220 | Zenodo | unknown | pietrobarbiero/digital-patient: Absolutno | exclude | not an RNA, protein or metabolite dataset (by its record) |
| zenodo.4030228 | Zenodo | unknown | pietrobarbiero/computational-patient: Absolutno | exclude | not an RNA, protein or metabolite dataset (by its record) |
| zenodo.4484341 | Zenodo | unknown | anfoss/HTPS_workflow v1.0 | exclude | not an RNA, protein or metabolite dataset (by its record) |
| zenodo.4624526 | Zenodo | unknown | DongshengChen-TY/COVID19: | exclude | not an RNA, protein or metabolite dataset (by its record) |
| zenodo.4932273 | Zenodo | unknown | Data from: Perinatal exposure to glyphosate and a glyphosate-based herbicide affect sperma | exclude | not an RNA, protein or metabolite dataset (by its record) |
| zenodo.4938153 | Zenodo | unknown | Data from: Preservation of RNA and DNA from mammal samples under field conditions | exclude | not an RNA, protein or metabolite dataset (by its record) |
| zenodo.4940432 | Zenodo | unknown | Data from: Epithelial magnesium transport by TRPM6 is essential for prenatal development a | exclude | not an RNA, protein or metabolite dataset (by its record) |
| zenodo.4959981 | Zenodo | unknown | Performance evaluation of a laboratory developed PCR test for quantitation of HIV-2 viral  | exclude | not an RNA, protein or metabolite dataset (by its record) |
| zenodo.5018162 | Zenodo | unknown | Data from: "Identification of SNP markers for the endangered Ugandan red colobus (Procolob | exclude | not an RNA, protein or metabolite dataset (by its record) |
| zenodo.5036924 | Zenodo | unknown | Cutaneous and systemic hyperinflammation drives maculopapular drug exanthema in severely i | exclude | not an RNA, protein or metabolite dataset (by its record) |
| zenodo.5081091 | Zenodo | unknown | cusom/CUSOM.COVIDome.Shiny-Apps: Version 2.3.0 | exclude | not an RNA, protein or metabolite dataset (by its record) |
| zenodo.5081093 | Zenodo | unknown | cusom/CUSOM.ShinyHelpers: Version 1.4.3 | exclude | not an RNA, protein or metabolite dataset (by its record) |
| zenodo.5127899 | Zenodo | unknown | FelixKrueger/TrimGalore: v0.6.7 - DOI via Zenodo | exclude | not an RNA, protein or metabolite dataset (by its record) |
| zenodo.5182825 | Zenodo | unknown | Inflammatory Responses in the Placenta upon SARS-CoV-2 Infection Late in Pregnancy - IHC d | exclude | not an RNA, protein or metabolite dataset (by its record) |
| zenodo.56794 | Zenodo | unknown | Pi_1.0.0 | exclude | not an RNA, protein or metabolite dataset (by its record) |
| zenodo.5748636 | Zenodo | unknown | jamrute/2021_COVID_amrute_steed: | exclude | not an RNA, protein or metabolite dataset (by its record) |
| zenodo.6568981 | Zenodo | unknown | Datasets from "Circulating miRNA and Lung Cancer: - a More Comprehensive Analysis of Avail | exclude | not an RNA, protein or metabolite dataset (by its record) |
| zenodo.6706499 | Zenodo | unknown | neogilab/COVID_GSMM: COVID_V1 | exclude | not an RNA, protein or metabolite dataset (by its record) |
| zenodo.7030528 | Zenodo | unknown | lasalletj/COVID_Neutrophils: COVID-19 Neutrophils | exclude | not an RNA, protein or metabolite dataset (by its record) |
| zenodo.7034435 | Zenodo | unknown | SUMMIT | exclude | not an RNA, protein or metabolite dataset (by its record) |
| zenodo.7076472 | Zenodo | unknown | Longitudinal characterization of circulating neutrophils uncovers distinct phenotypes asso | exclude | not an RNA, protein or metabolite dataset (by its record) |
| zenodo.7114840 | Zenodo | unknown | kuanrongchan/COVID19-severity: The codes to make the COVID-19 database website is now publ | exclude | not an RNA, protein or metabolite dataset (by its record) |
| zenodo.7333789 | Zenodo | unknown | jackgisby/covid-longitudinal-multi-omics: Published version of the code | exclude | not an RNA, protein or metabolite dataset (by its record) |
| zenodo.7384115 | Zenodo | unknown | nf-core/atacseq: nf-core/atacseq v2.0 - Iron Iguana | exclude | not an RNA, protein or metabolite dataset (by its record) |
| zenodo.7455070 | Zenodo | unknown | zj-zhang/CHARM-AlternativeSplicing: Version freeze for publication | exclude | not an RNA, protein or metabolite dataset (by its record) |
| zenodo.7455172 | Zenodo | unknown | zhanglab-aim/JEMM: Version freeze for publication | exclude | not an RNA, protein or metabolite dataset (by its record) |
| zenodo.7555525 | Zenodo | unknown | cxli233/customized_upset_plots: customized_upset_plots_v3 | exclude | not an RNA, protein or metabolite dataset (by its record) |
| zenodo.8394171 | Zenodo | unknown | Alamar-Biosciences/NULISAseqR: v1.0.0 | exclude | not an RNA, protein or metabolite dataset (by its record) |
| EGAD00001007931 | controlled access archive | unknown | A blood atlas of COVID-19 defines hallmarks of disease severity and specificity. | exclude | not open access (application needed) |
| EGAD00001007932 | controlled access archive | unknown | A blood atlas of COVID-19 defines hallmarks of disease severity and specificity. | exclude | not open access (application needed) |
| EGAD00001007957 | controlled access archive | unknown | A blood atlas of COVID-19 defines hallmarks of disease severity and specificity. | exclude | not open access (application needed) |
| EGAD00001007959 | controlled access archive | unknown | A blood atlas of COVID-19 defines hallmarks of disease severity and specificity. | exclude | not open access (application needed) |
| EGAD00001007960 | controlled access archive | unknown | A blood atlas of COVID-19 defines hallmarks of disease severity and specificity. | exclude | not open access (application needed) |
| EGAD00001007961 | controlled access archive | unknown | A blood atlas of COVID-19 defines hallmarks of disease severity and specificity. | exclude | not open access (application needed) |
| EGAD00001007962 | controlled access archive | unknown | A blood atlas of COVID-19 defines hallmarks of disease severity and specificity. | exclude | not open access (application needed) |
| EGAD00001007963 | controlled access archive | unknown | A blood atlas of COVID-19 defines hallmarks of disease severity and specificity. | exclude | not open access (application needed) |
| EGAD00001007964 | controlled access archive | unknown | A blood atlas of COVID-19 defines hallmarks of disease severity and specificity. | exclude | not open access (application needed) |
| EGAD00001007965 | controlled access archive | unknown | A blood atlas of COVID-19 defines hallmarks of disease severity and specificity. | exclude | not open access (application needed) |
| EGAD00001007995 | controlled access archive | unknown | Longitudinal single cell atlas identifies complex temporal relationship between type I int | exclude | not open access (application needed) |
| EGAD00001008007 | controlled access archive | unknown | A blood atlas of COVID-19 defines hallmarks of disease severity and specificity. | exclude | not open access (application needed) |
| EGAD00001009331 | controlled access archive | unknown | Single-cell immune aging clocks reveal inter-individual heterogeneity during infection and | exclude | not open access (application needed) |
| EGAD50000000543 | controlled access archive | unknown | COVID-19 progression and convalescence in common variable immunodeficiency patients show d | exclude | not open access (application needed) |
| EGAD50000000684 | controlled access archive | unknown | Predicting trajectories of illness using RNA velocity of whole blood. | exclude | not open access (application needed) |
| EGAD50000000942 | controlled access archive | unknown | Predicting trajectories of illness using RNA velocity of whole blood. | exclude | not open access (application needed) |
| EGAD50000000956 | controlled access archive | unknown | Predicting trajectories of illness using RNA velocity of whole blood. | exclude | not open access (application needed) |
| EGAS00001004085 | controlled access archive | unknown | Single-Cell RNAseq Profiling of Human γδ T Lymphocytes in Virus-Related Cancers and COVID- | exclude | not open access (application needed) |
| EGAS00001004419 | controlled access archive | unknown | ELF5 is a potential respiratory epithelial cell-specific risk gene for severe COVID-19. | exclude | not open access (application needed) |
| EGAS00001004503 | controlled access archive | unknown | Hypoxia and Activation of Neutrophil Degranulation-Related Genes in the Peripheral Blood o | exclude | not open access (application needed) |
| EGAS00001004571 | controlled access archive | unknown | Longitudinal characterization of circulating neutrophils uncovers phenotypes associated wi | exclude | not open access (application needed) |
| EGAS00001004689 | controlled access archive | unknown | ELF5 is a potential respiratory epithelial cell-specific risk gene for severe COVID-19. | exclude | not open access (application needed) |
| EGAS00001004717 | controlled access archive | unknown | Severe COVID-19 patients display hyper-activated NK cells and NK cell-platelet aggregates. | exclude | not open access (application needed) |
| EGAS00001004928 | controlled access archive | unknown | Longitudinal characterization of circulating neutrophils uncovers phenotypes associated wi | exclude | not open access (application needed) |
| EGAS00001005332 | controlled access archive | unknown | Coagulation factor V is a T-cell inhibitor expressed by leukocytes in COVID-19. | exclude | not open access (application needed) |
| EGAS00001005446 | controlled access archive | unknown | Single-cell immune aging clocks reveal inter-individual heterogeneity during infection and | exclude | not open access (application needed) |
| EGAS00001005461 | controlled access archive | unknown | ELF5 is a potential respiratory epithelial cell-specific risk gene for severe COVID-19. | exclude | not open access (application needed) |
| EGAS00001005493 | controlled access archive | unknown | COMBATdb: a database for the COVID-19 Multi-Omics Blood ATlas. | exclude | not open access (application needed) |
| EGAS00001005529 | controlled access archive | unknown | Single-cell immune aging clocks reveal inter-individual heterogeneity during infection and | exclude | not open access (application needed) |
| EGAS00001005545 | controlled access archive | unknown | Longitudinal single cell atlas identifies complex temporal relationship between type I int | exclude | not open access (application needed) |
| EGAS00001005594 | controlled access archive | unknown | Single-cell immune aging clocks reveal inter-individual heterogeneity during infection and | exclude | not open access (application needed) |
| EGAS00001005634 | controlled access archive | unknown | Longitudinal characterization of circulating neutrophils uncovers phenotypes associated wi | exclude | not open access (application needed) |
| EGAS00001006778 | controlled access archive | unknown | Multi-omics identify falling LRRC15 as a COVID-19 severity marker and persistent pro-throm | exclude | not open access (application needed) |
| EGAS00001006787 | controlled access archive | unknown | Single-cell immune aging clocks reveal inter-individual heterogeneity during infection and | exclude | not open access (application needed) |
| EGAS00001006990 | controlled access archive | unknown | Single-cell immune aging clocks reveal inter-individual heterogeneity during infection and | exclude | not open access (application needed) |
| EGAS00001007358 | controlled access archive | unknown | Human haematopoietic stem cells remember inflammatory stress. | exclude | not open access (application needed) |
| EGAS00001007886 | controlled access archive | unknown | Blood DNA virome associates with autoimmune diseases and COVID-19. | exclude | not open access (application needed) |
| EGAS00001007955 | controlled access archive | unknown | Blood DNA virome associates with autoimmune diseases and COVID-19. | exclude | not open access (application needed) |
| EGAS00001008016 | controlled access archive | unknown | Deciphering state-dependent immune features from multi-layer omics data at single-cell res | exclude | not open access (application needed) |
| EGAS50000000142 | controlled access archive | unknown | Single-cell immune aging clocks reveal inter-individual heterogeneity during infection and | exclude | not open access (application needed) |
| EGAS50000001623 | controlled access archive | unknown | Human haematopoietic stem cells remember inflammatory stress. | exclude | not open access (application needed) |
| EGAS50000001624 | controlled access archive | unknown | Human haematopoietic stem cells remember inflammatory stress. | exclude | not open access (application needed) |
| syn35874390 | controlled access archive | unknown | Mast cell activation may contribute to adverse health transitions in COVID-19 patients wit | exclude | not open access (application needed) |
| syn51228232 | controlled access archive | unknown | ACE2-independent SARS-CoV-2 virus entry through cell surface GRP78 on monocytes - evidence | exclude | not open access (application needed) |
| syn51730532 | controlled access archive | unknown | A single-cell compendium of human cerebrospinal fluid identifies disease-associated immune | exclude | not open access (application needed) |
| MGH_Olink_Filbin2021 | memory | unknown | candidate from memory (to verify) | exclude | not an RNA, protein or metabolite dataset (by its record) |
| E-MTAB-14728 | other | unknown | Predicting trajectories of illness using RNA velocity of whole blood. | exclude | not an RNA, protein or metabolite dataset (by its record) |
| E-MTAB-2836 | other | unknown | The Landscape of Expressed Chimeric Transcripts in the Blood of Severe COVID-19 Infected P | exclude | not an RNA, protein or metabolite dataset (by its record) |
| E-MTAB-4421 | other | unknown | A conserved monocyte activation program links brain injury to systemic immune adaptation a | exclude | not an RNA, protein or metabolite dataset (by its record) |
| E-MTAB-5195 | other | unknown | Multi-cohort analysis of host immune response identifies conserved protective and detrimen | exclude | not an RNA, protein or metabolite dataset (by its record) |
| E-MTAB-8077 | other | unknown | Inflammatory risk contributes to post-COVID endothelial dysfunction through anti-ACKR1 aut | exclude | not an RNA, protein or metabolite dataset (by its record) |
| PRJNA63448 | other | unknown | Platelets of COVID-19 patients display mitochondrial dysfunction, oxidative stress, and en | exclude | not an RNA, protein or metabolite dataset (by its record) |
| PRJNA76339 | other | unknown | Role of miR-2392 in driving SARS-CoV-2 infection. | exclude | not an RNA, protein or metabolite dataset (by its record) |
