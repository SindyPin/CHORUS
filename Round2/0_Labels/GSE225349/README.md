# GSE225349: SomaScan proteomics of the Mount Sinai COVID-19 Biobank

**Data.** SomaScan v4 serum proteomics (Paranjpe et al., Commun Med 2023, PMID 37308534; GEO GSE225349, public).
- **Same people as GSE215865** (RNA-seq) and the Synapse clinical table. Samples are linked by `Blood_Sample_ID`.
- **Built as an acute proteomics cohort. No model has been run on it.**

## Build

**Matrix.** The "covid only" matrix: 1,105 samples, 571 subjects, every time point, log2 RFU.
- It agrees with log2 of the raw "covid AKI" file: r = 0.97 on a shared sample.
- The covid AKI file holds 437 samples, one per patient.

**Protein mapping.**
- The matrix writes each SeqId as a number (`10000.28`), which drops trailing zeros (`10063-10` becomes `10063.1`).
  So aptamers are matched to the annotation on the numeric form. That form is unique in the annotation, and the two
  files are also in the same row order.
- All 4,930 aptamers are matched. The Round1 `somascan` adapter gives 4,446 UniProt features; 37 aptamers have no
  UniProt id.

**Metadata.**
- Subject ids are `GSE215865:Subj_…`, the same as GSE215865.
- `lc_status` (Amendment 7): the post COVID checklist label per subject (`../GSE215865_checklist`), with COVID-19
  status from this series: 170 LC, 32 recovered, 82 never infected, 287 unknown (never answered). Symptoms only
  (sensitivity): 167 / 35 / 82 / 287. The TACO labels are withdrawn.

**Files.** Written to `data/` (matrix, features, metadata; not committed). `output/build_summary.csv` has the
aggregate counts.

## Overlap on the first acute sample (`output/overlap.csv`)

The first acute sample is the earliest RNA-seq sample of each Round2 early detection participant. Labels of
Amendment 7.

| | main: LC | main: recovered | symptoms only: LC | symptoms only: recovered |
|---|---|---|---|---|
| participants with an RNA-seq first acute sample | 159 | 29 | 156 | 32 |
| with Synapse clinical data | 159 | 29 | 156 | 32 |
| with SomaScan on any sample | 157 | 28 | 154 | 31 |
| **RNA-seq, SomaScan and clinical data on the same first acute sample** | **135** | **21** | **132** | **24** |
| that sample is also their first SomaScan sample | 125 | 20 | 122 | 23 |
| GSE225349 subjects with any acute SomaScan sample | 170 | 32 | 167 | 35 |

## Script

`build.py`, run from `Code/`.
