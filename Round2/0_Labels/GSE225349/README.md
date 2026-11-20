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
- `lc_status` is the current GSE215865 (TACO) label per subject: 290 LC, 116 recovered, 63 never infected, 102
  unknown.

**Files.** Written to `data/` (matrix, features, metadata; not committed). `output/build_summary.csv` has the
aggregate counts.

## Overlap on the first acute sample (`output/overlap.csv`)

The first acute sample is the earliest RNA-seq sample of each Round2 early detection participant.

| | LC | recovered |
|---|---|---|
| participants with an RNA-seq first acute sample | 302 | 122 |
| with Synapse clinical data | 302 | 122 |
| with SomaScan on any sample | 290 | 116 |
| **RNA-seq, SomaScan and clinical data on the same first acute sample** | **236** | **93** |
| that sample is also their first SomaScan sample | 223 | 88 |
| of the 236 / 93: checklist PASC | 110 | 21 |
| of the 236 / 93: checklist no PASC | 12 | 8 |
| of the 236 / 93: no checklist answer | 114 | 64 |

The labels are TACO's; the checklist rows are for information (see `../GSE215865_checklist`).

## Script

`build.py`, run from `Code/`.
