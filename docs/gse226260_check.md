# GSE226260: how the PBMC and whole blood sample sets were defined

This file compares the GEO series record with the paper. It does not interpret any benchmark result.

**Sources (checked 13 November 2026).**
- GEO series GSE226260 and the GEO sample records, including the per-sample protocols of GSM7069899 (PBMC),
  GSM8989185 (whole blood) and GSM8260984 (`Blood`).
- Paper: Aid M, et al. "Long COVID involves activation of proinflammatory and immune exhaustion pathways." Nat Immunol
  27, 61–71 (2026; online 12 December 2025). PMID 41388153, PMC12764429, doi 10.1038/s41590-025-02353-x.
  - Its author correction (PMID 41578033, PMC12956568) only fixes an author's surname.
- The paper's code repository, github.com/Barouch-and-BCGE-Labs/PASC-Project. Its scripts read `Cohort1.csv` and
  `Cohort2.csv` sample sheets (SampleID, SubjectID, Time, PASCSTATUS), which are not in the repository.

The config comments elsewhere in this repository call the paper "Aid 2025"; it appeared online in 2025 and in the 2026
volume.

## The three GEO subsets

GSE226260 holds 331 samples in three subsets, which were deposited at different times.

| GEO source_name | samples | title pattern | first submitted | used in CHORUS as |
|---|---|---|---|---|
| PBMC | 86 | `p22179-s001_116521-02-Apr-2021` (project, subject, collection date) | 27 February 2023 | GSE226260_PBMC |
| Blood | 142 | `p20169.s001_1_PT023` | 9 May 2024 | not used (no LC / recovered label) |
| Whole blood | 103 | `BI1419.140.Recovered` (subject, days, status) | 16 May 2025 | GSE226260_WB |

The GEO series summary is the paper abstract. The GEO overall design describes only PBMCs: "Blood samples were collected
from 178 patients ... PBMCs were harvested and RNA was extracted, and used for total RNA-seq". It does not mention the
whole blood subset.

## They are two cohorts, enrolled in different periods

The paper has two cohorts.

| | 2020–2021 cohort | 2023–2024 cohort |
|---|---|---|
| paper role | discovery | "validation", "independent" |
| enrolment | MassCPR observational cohort and BIDMC, Boston, April 2020 to October 2021 | BIDMC clinical site of NIH RECOVER and the BIDMC biorepository, October 2022 to December 2024 |
| groups (paper) | uninfected 35, acute 54, convalescent controls (CC) 24, LC 28 | CC 20, LC 18 |
| symptom ascertainment | self-reported. GEO adds: "Severity of long COVID was determined by self-reporting, and validated by a clinician's assessment" | RECOVER standardised symptom questionnaire (the paper cites Thaweethai et al. JAMA 2023) |
| sample type in the paper | PBMC | "peripheral blood" |

The GEO subsets map onto the two cohorts as follows.
- **PBMC subset = 2020–2021 cohort.**
  - The collection dates run from 13 April 2020 to 18 November 2021.
  - The extraction protocol matches the paper's PBMC methods word for word.
  - There are 28 PASC subjects, the paper's LC n = 28.
- **Whole blood subset = 2023–2024 cohort.**
  - There are 18 PASC and 20 Recovered subjects, exactly the paper's LC n = 18 and CC n = 20 (Table 2).
- **Blood subset (p20169).** Its protocol text also describes PBMCs. It carries only a `pasc severity` field: critical 61,
  severe 23, moderate 7, mild 7, asymptomatic 6, empty 38. It has 99 distinct `PT` ids. GEO does not say which paper group
  these samples belong to.

**Are they different participants?**
- The paper describes two separately enrolled cohorts from different periods and calls the second one independent.
- The two subsets use different id systems: numeric ids such as `116521` for PBMC and `BIxxxx` for whole blood.
- Neither GEO nor the paper links ids across the cohorts, and the paper does not say whether anyone was in both.
- So the overlap cannot be checked from public data. `deduplicate_subjects` cannot link people across the two subsets
  (`data/meta/curation_notes.md`).

## What COVID_recovered means in each subset

CHORUS maps the GEO labels to COVID_recovered as follows (`scripts/build_cohort_meta.py`):
- PBMC: GEO `pasc status` NOPASC.
- Whole blood: GEO `current long_covid_status` Recovered.

**PBMC (2020–2021).**
- The paper calls this group "convalescent controls" and "individuals who recovered". In the acute analysis it says
  "acutely infected patients who fully recovered".
- Symptoms were self-reported. The paper gives no symptom duration or questionnaire criterion for the LC versus CC split.
- The status is constant within each subject in GEO.

**Whole blood (2023–2024).**
- The paper calls this group CC. Symptoms come from the RECOVER standardised questionnaire.
- The paper does not say whether the RECOVER PASC index threshold or another rule set the LC versus CC split.
- The GEO field is named *current* long COVID status. It is constant within each subject across all of that subject's
  draws.

**GEO and paper counts do not match for the PBMC controls.**

| | GEO PBMC subset | paper, 2020–2021 cohort |
|---|---|---|
| LC / PASC subjects | 28 | 28 |
| CC / NOPASC subjects | 8 | 24 |
| LC / PASC PBMC samples | 72 | 53 (acute 6, 90–180 days 26, >180 days 21) |
| CC / NOPASC PBMC samples | 14 | 31 (acute 5, 90–180 days 21, >180 days 5) |

The paper's bulk RNA-seq comparison used "LC (n = 26) and CCs (n = 21) at day 90–180". GEO has 14 NOPASC PBMC samples in
total. GEO does not explain the difference, and the paper's per-sample sheet is not public.

## When the PBMC controls were collected relative to the cases

GEO gives each PBMC sample a **collection date** but **no infection date**, so days since infection cannot be derived. The
paper groups PBMC draws into acute (<30 days), 90–180 days and >180 days after infection, but without per-sample
assignments.

Collection dates by status (all draws):

| quarter | NOPASC | PASC |
|---|---|---|
| 2020 Q2 | 2 | 5 |
| 2020 Q3 | 5 | 23 |
| 2020 Q4 | 3 | 21 |
| 2021 Q1 | 1 | 14 |
| 2021 Q2 | 2 | 5 |
| 2021 Q3 | 1 | 2 |
| 2021 Q4 | 0 | 2 |

| | NOPASC (8 subjects) | PASC (28 subjects) |
|---|---|---|
| draws per subject | 14 samples | 72 samples |
| first draw, median date | 2 September 2020 | 25 July 2020 |
| first draw, range | 13 April 2020 to 18 May 2021 | 28 April 2020 to 12 May 2021 |
| span from first to last draw, median days | 81 | 122 |
| span from first to last draw, max days | 185 | 359 |
| age at first draw, median (range) | 39 (27–49), 7 with age | 49 (27–70), 27 with age |
| female | 6 of 8 | 24 of 28 |

- **Calendar time.** Controls and cases were collected over the same period, April 2020 to mid 2021.
- **Time since infection.** Unknown per sample. According to the paper, 5 of the 31 CC PBMC samples, against 21 of the
  53 LC samples, were drawn more than 180 days after infection.
- **Acute draws among the kept samples.** The paper reports acute (<30 days) PBMC draws in both groups: LC 6, CC 5.
  - CHORUS keeps each PBMC subject's earliest collection date (`one_sample_per_subject`).
  - GEO has no infection date, so it cannot tell whether that earliest draw is acute.
  - Some of the PBMC samples used could therefore be acute draws.
- **One inconsistent date.** GSM7069947 (`p22179-s049_601096-14-Dec-2020`) has `timepoint` 7/1/2020, which differs from the
  date in its title. It is the sample that carried the malformed tail of the combined matrix (config comment).

**Whole blood, for comparison.** The second title field agrees with the GEO draw window in every sample, so it is read as
days since infection (`data/meta/curation_notes.md`).

| GEO draw window | PASC | Recovered |
|---|---|---|
| Acute | 4 | 0 |
| 3_6M | 15 | 16 |
| 6_12M | 11 | 20 |
| >12M | 22 | 15 |

- Days range from 15 to 1678 for PASC and from 100 to 954 for Recovered. The paper states "days 15–700 after infection".
- CHORUS keeps each subject's earliest post-acute draw:

  | | PASC (17 subjects) | Recovered (20 subjects) |
  |---|---|---|
  | median days | 143 | 126 |
  | range (days) | 79–833 | 100–368 |
  | 3_6M / 6_12M / >12M | 10 / 4 / 3 | 13 / 6 / 1 |

## Library preparation differs between the two subsets

| | PBMC (2020–2021) | whole blood (2023–2024) |
|---|---|---|
| extraction | TRIzol, miRNeasy Mini with on-column DNase | RNeasy Plus Mini |
| library | 10 ng total RNA, SMART-Seq v4 Ultra Low Input, Nextera XT | poly(A) selection (oligo d(T)25), NEBNext Ultra II Directional |
| sequencing | NovaSeq 6000, PE100, about 30 M reads | NovaSeq X Plus 10B, PE150, about 40 M reads |
| processing (GEO) | STAR 2.7.3a, htseq-count, GRCh38.p14 plus SARS-CoV-2 | Trimmomatic, STAR 2.7.1a, StringTie, featureCounts or HTSeq, GRCh38.p14 plus SARS-CoV-2 |

So in GSE226260, sample matrix (PBMC versus whole blood) goes together with:
- the cohort and enrolment period;
- the symptom ascertainment method;
- the library chemistry (full-length SMART-Seq total RNA versus poly(A) stranded);
- the sequencer and read length;
- the processing pipeline.
