# Author count matrices: verdicts

Downloaded on 2026-11-09 from the GEO supplementary folders into `data/counts/<GSE>/` (git-ignored). The logCPM and
voom weight matrices of GSE215865 were not downloaded because they are derived from the raw count file. Per-GSM sample
tables (title, source, characteristics, BioSample, SRX) were written to `data/meta/<GSE>_geo_samples.csv` and used
for the id checks. No sequencing reads were downloaded.

## One line per file

| file | verdict |
|---|---|
| GSE215865_rnaseq_raw_count_matrix.csv.gz | integer counts; genes x samples; Ensembl versioned (ENSG…\.5); 58,929 genes x 1,392 samples; ids are not GSM titles but match `blood_sample_id` + `_` + `library_prep_plate` 1,392/1,392 (e.g. `Subj_16ae238fT0_Plate_2`), unique, so the mapping to GSM is exact. |
| GSE226260_AdditionalSamples.rawCounts.csv.gz | integer counts; genes x samples; Ensembl unversioned; 78,932 genes x 103 samples; ids equal GSM titles 103/103, all of them the whole blood subset. |
| GSE226260_AdditionalSamples.logNormalizedCounts.csv.gz | normalised (log scale, negative values); genes x samples; Ensembl; 78,932 x 103; ids equal GSM titles 103/103 (whole blood). Not needed: use the raw counts. |
| GSE226260_combined.ExprData.updated.csv.gz | integer counts; genes x samples; Ensembl unversioned; 24,203 genes x 228 samples (86 PBMC + 142 'Blood' acute); ids equal GSM titles 228/228. **Malformed tail**: 17,345 extra rows with an empty gene id, each holding a single value in column `p22179-s049_601096-14-Dec-2020`; dropped. Treat sample s049 (GSM of that title) with suspicion until checked. |
| GSE267625_P4O2_COVID_blood_transcriptome.csv.gz | not integer (24.7% of values fractional, e.g. 4.008; non-negative, column sums ~31 M): estimated counts (Salmon/RSEM style), not normalised; genes x samples; Ensembl unversioned; 34,515 genes x 111 samples; ids are `<subject>_<index>` (e.g. `AA002_1`) and match no GSM title, but the subject sequence of the 111 columns is identical to the subject sequence of the GSM titles in GSM order, so a positional mapping is plausible. 25 subjects have two samples (3 and 12 months), so the timepoint can only be recovered by position. Prefer re-quantification (reads are public). |
| GSE251849_Counts.txt.gz | not integer (10.2% fractional, non-negative, column sums ~15 M): estimated counts; genes x samples; Ensembl unversioned; 62,710 genes x 23 samples; ids `Control001-007`, `Convalescent001-005`, `LongCOVID001-011` match no GSM title. Groups map by count to Healthy 7 / COVID convalescent 5 / Long COVID 6 + brain fog 5, but which LongCOVID columns are brain fog is unknown. LC vs recovered labels are unaffected by that. |
| GSE270045_LC_counts.tsv.gz | not integer (48.5% fractional, column sums ~3 M): normalised or estimated, not raw counts; genes x samples; **gene symbols** (incl. `5S_rRNA`); 28,889 genes x 36 samples; ids match the `sample id` characteristic 36/36 (not the GSM titles). Re-quantify from SRA. |
| GSE169687_RAW.tar (152 featureCounts files) | integer counts; one file per GSM, merged to genes x samples `GSE169687_merged_counts.csv.gz`; Ensembl unversioned (GRCh38.101, HISAT2 + featureCounts, stranded); 60,671 genes x 152 samples; file names carry the GSM accession, 152/152 matched. |
| GSE222253_RAW.tar (35 count files) | integer counts (HTSeq style, five `__` summary rows dropped); one file per GSM, merged to `GSE222253_merged_counts.csv.gz`; Ensembl unversioned; 60,664 genes x 35 samples; file names carry the GSM accession, 35/35 matched. |

## GSE226260 split

GEO puts three cohorts into this series, distinguished by `source_name` and the `tissue` and `cell type`
characteristics. They were written as separate matrices (genes x GSM, integer counts, Ensembl):

| subset | rule | file | genes x samples | labels in GEO |
|---|---|---|---|---|
| whole blood | tissue = Whole blood | `GSE226260_wholeblood_counts.csv.gz` | 78,932 x 103 | current long_covid_status: 52 PASC / 51 Recovered; time of blood draw 3–6 M, 6–12 M, >12 M, 4 Acute |
| PBMC | cell type = PBMC | `GSE226260_PBMC_counts.csv.gz` | 24,203 x 86 | pasc status: 72 PASC / 14 NOPASC |
| acute blood | source = Blood (no cell type) | `GSE226260_blood_acute_counts.csv.gz` | 24,203 x 142 | pasc severity only (acute), no LC field |

The whole blood and PBMC subsets share 23,864 genes. The whole blood subset comes from a different count file
with a different annotation, so the subsets must be harmonised (or re-quantified together) before they are compared.

## PRJNA1184005 (Missailidis 2024, Front Immunol 15:1450853, PMID 39691709, PMC11649547)

- 14 BioSamples (SAMN44656608–SAMN44656621) and 14 SRA runs, all PBMC, NovaSeq 6000, 13 paired-end plus 1
  single-end (vacSG112). Attributes are age, sex, tissue, collection date, geo_loc_name and isolate (`SG22`,
  `vacSG21` … `vacSG117`); written to `data/meta/PRJNA1184005_biosamples.csv`. **There is no LC or recovered attribute.**
- No processed matrix: the BioProject is raw reads only and has no GEO series. The paper's data statement points
  only to PRJNA1184005, and the article has no supplementary files. Table 1 gives only group summaries (7 LC aged
  23–53, 7 recovered aged 23–63). The paper's tables are a group summary and a 70-gene DE list.
- Consequence: re-quantify from SRA (`requantify`). Labels must come from the authors (see `data/meta/label_sources.md`).
