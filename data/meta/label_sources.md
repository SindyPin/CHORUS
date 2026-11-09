# Label sources: Long COVID versus COVID recovered

Compiled on 2026-11-09. For each transcriptomic cohort, this file says where the per-sample `lc_status` label comes
from, what has been extracted, and what still has to be requested. The labels extracted so far are in
`data/meta/labels_geo_and_taco.csv` (one row per GSM: study_id, gsm, title, subset, lc_status). The GEO sample
tables are `data/meta/<GSE>_geo_samples.csv`.

## Summary

| cohort | label source | per-sample labels in hand | LC / recovered (other) | status |
|---|---|---|---|---|
| GSE226260 whole blood | GEO characteristic `current long_covid_status` | yes, 103 GSM | 52 / 51 samples; 18 / 20 subjects | done |
| GSE226260 PBMC | GEO characteristic `pasc status` | yes, 86 GSM | 72 / 14 samples; 28 / 8 subjects | done |
| GSE224615 | GEO characteristic `lc status` | yes, 36 GSM | 23 / 13 | done |
| GSE251849 | GEO characteristic `group` | yes, 23 GSM | 11 / 5 (7 never infected) | done |
| GSE215865 | TACO clinical table (`Long_COVID`) joined on GEO `subject id` | yes, 1,370 of 1,392 GSM | subjects 303 / 186 (incl. 64 COVID negative) | done, with caveats |
| GSE270045 | GEO title + paper (controls are recovered) | yes, 36 GSM | 19 / 17 | done, one discrepancy |
| GSE169687 | paper: referral to a long COVID clinic | no | 21 / 48 subjects (14 never infected) | **request from authors** |
| GSE267625 | paper supplement / P4O2 data committee | no | 77 / 18 (config) | **request** |
| GSE222253 | paper / BQC19 biobank | no | 10 / 14 patients | **request** |
| PRJNA1184005 | authors (no label in BioSample or paper) | no; 3 of 7 recovered inferable by age | 7 / 7 | **request** |

## Per cohort

### GSE226260 (Aid 2025)
GEO characteristics. Whole blood (`tissue = Whole blood`): `current long_covid_status` PASC → LC, Recovered →
COVID_recovered. There are repeated samples per subject (38 subjects: 18 PASC, 20 recovered) at 3–6 M, 6–12 M,
>12 M and four acute draws. Use subject-level de-duplication or a timepoint choice, and drop the 4 `Acute` draws.
PBMC (`cell type = PBMC`): `pasc status` PASC / NOPASC (36 subjects: 28 / 8). The 142 `Blood` samples carry only
acute severity and no LC label.

### GSE224615 (Yin 2024)
GEO `lc status` LC / Non-LC, 23 / 13 at 8 months. The paper cohort is 27 / 16. The remaining participants have
no RNA-seq in GEO.

### GSE251849 (Greene 2024)
GEO `group`: Long COVID (6) + Long COVID brain fog (5) → LC; COVID convalescent (5) → COVID_recovered; Healthy (7)
→ never_infected. The count matrix uses its own ids (`LongCOVID001`…), and the brain fog split within LC cannot be
mapped (see `docs/count_matrices.md`). That split is not needed for LC vs recovered.

### GSE215865 (Thompson 2023, Mount Sinai) via TACO
Source: `github.com/SindyPin/TACO`, commit cac2afb (2 Apr 2026), `Data/411_DCEG_Paths_average.csv`, columns
`Subject_ID` and `Long_COVID`. 489 subjects (303 with `Long_COVID` = 1, 186 with 0). `Subject_ID` uses the same
`Subj_xxxxxxxx` form as the GEO characteristic `subject id`, so it joins directly.
- **1,370 of 1,392 GSM** (489 of 508 subjects) map to a TACO label. 22 GSM from 19 subjects have no label (16
  COVID positive, 6 negative). Saved per GSM in `data/meta/GSE215865_labels_from_TACO.csv`.
- `Long_COVID = 0` is not the same as COVID recovered. Among subjects, 0 contains 122 COVID-positive and 64
  COVID-negative (by GEO `covid-19 positive`) subjects. In `labels_geo_and_taco.csv`, 0 + positive is mapped to
  COVID_recovered and 0 + negative to `never_infected_or_negative`. One subject is labelled LC but COVID negative
  in GEO and should be checked.
- Samples are acute hospital draws (T0–T13), so this cohort is the acute transfer set, not a post-acute cohort.
  The original clinical table is Synapse syn35874390 (registration and data use terms), and the TACO labels
  should be checked against it before publication.

### GSE270045 (Shahbaz 2025, Cell Rep Med, PMID 41205594, PMC12711683)
The paper settles the control question: "We randomly selected 18 patients (15 LCF and 3 LCM) and 17 control
individuals (12 RF and 5 RM) for bulk RNA-seq", where R are "COVID-19 patients who had recovered from the disease
without any obvious symptoms and complications". **The 17 GEO 'Healthy Control' samples are COVID recovered**, not
never infected. Discrepancy: GEO has 19 'Long Covid' samples while the paper says 18 patients. Ask the lead
contact which sample is extra. The LC group is LC with ME/CFS (a stricter phenotype).
Contact: Shokrollah Elahi, elahi@ualberta.ca.

### GSE169687 (Ryan 2022, BMC Med, PMID 35027067, PMC8758383)
- Controls: 14 seronegative healthy controls "with no history of prior SARS-CoV-2 infection" (never_infected).
- LC label: the paper defines LC as referral to the SA Health long COVID clinic, which applies to 21 of 69
  convalescents (30.4%). A retrospective questionnaire was completed by 57. The paper's DE contrast is referred
  vs not referred at 24 wpi. **This cohort therefore does have an LC vs COVID recovered contrast** (21 vs 48
  subjects), unlike the config comment.
- Downloaded open access supplements: Table S1 (MOESM4: per-subject sample availability, severity, timeline),
  Table S4 (MOESM7: per-library sample statistics, SampleID ↔ subject CV id, severity, timepoint), Table S5
  (MOESM8: referred vs not DE results). **None has the per-subject referral flag.** The GEO `subject id`
  (e.g. 15-063) also differs from the paper's `CV` numbers. Table S4 maps library ids to CV ids.
- Request: the per-subject referral status (and questionnaire answers) keyed by GEO subject id or CV number.
  Contact: David Lynn, david.lynn@sahmri.com (also simon.barry@adelaide.edu.au, branka.grubor@adelaide.edu.au).

### GSE267625 (Blankestijn 2024, J Allergy Clin Immunol, doi 10.1016/j.jaci.2024.04.032, PMID 38830512)
- Not open access (no PMC copy; publisher and Elsevier API return 400/403), so the supplement could not be
  downloaded. GEO carries no LC label: `treatment` and `genotype` are placeholders, only `time` (3 or 12 months) is
  informative.
- Request: per-sample LC / recovered status (and the clinical variables) keyed by the P4O2 subject id in the GEO
  titles (e.g. `AA002`, `VV001`) and timepoint. GEO says: "Additional clinical variables can be requested
  through the P4O2 data committee (p4o2@amsterdamumc.nl)". Corresponding authors are in the P4O2 consortium
  (Blankestijn JM, Maitland-van der Zee AH, Amsterdam UMC). If library access works, check the JACI supplement
  first.

### GSE222253 (PMID 37680625, Front Immunol 2023, doi 10.3389/fimmu.2023.1243689, PMC10482103)
- Open access. Supplement DataSheet_1.docx was downloaded from Frontiers. It contains questionnaire definitions
  and group-level tables only: 14 patients without and 10 with post-COVID symptoms, and endotype tables
  (Resolved 12 / Suppressive 7 / Unresolved 5). **There are no per-sample labels.** GEO gives only the endotype
  cluster, hospital vs discharged location and batch.
- All 24 patients were hospitalised (BQC19). Endotype is not a proxy for post-COVID symptoms.
- Request: per-patient post-COVID symptom status keyed by GEO sample title (`sample125` …). Contact: Robert
  (Bob) Hancock, bob@hancocklab.com, or the BQC19 biobank (https://en.quebeccovidbiobank.ca/) via its access process.
- Related series GSE221234 (124 GSM, COVID-19 severity endotypes, same group) is acute only and not needed here.

### PRJNA1184005 (Missailidis 2024, Front Immunol 15:1450853, PMID 39691709, PMC11649547)
- BioSample attributes (`data/meta/PRJNA1184005_biosamples.csv`): age, sex, tissue PBMC, collection date,
  isolate. **No LC field.** The paper has no supplementary files and gives group summaries only (Table 1:
  LC aged 23–53, recovered aged 23–63).
- By age alone, three samples must be COVID recovered: vacSG34 (63), vacSG102 (62), vacSG112 (62). The other 11
  cannot be assigned.
- Request: the LC / recovered key for the 14 isolates (`SG22`, `vacSG21` … `vacSG117`). Contact: Sarah Annesley,
  S.Annesley@latrobe.edu.au.
