# Provenance of re-quantified matrices and the run environment

## GSE224615 (Yin 2024, LIINC PBMC): salmon on a laptop

**Source reads.**
- ENA project PRJNA931952: 36 paired-end runs, about 230 GB of FASTQ.
- Downloaded over HTTPS from the ENA `fastq_ftp` links. Each file was checked against the ENA `fastq_md5` before use and
  deleted after quantification.
- ENA limits a single connection to about 1 MB/s and answers HTTP 403 when loaded. `scripts/salmon_quant_ena.py`
  therefore downloads each file in 12 parallel byte ranges, retrying and resuming each range with back-off.
- The ENA file report used is saved as `data/quant/GSE224615/ena_filereport.tsv`.

**Library preparation (GEO sample protocol, for example GSM7027478).**
- Molecule: `polyA RNA`. PBMCs were incubated in RPMI 1640 with 10% FBS, then total RNA was extracted with the Qiagen
  AllPrep kit.
- Library: "following oligo dT enrichment, fragmentation and random priming, cDNA synthesis", then end repair,
  dA-tailing, adaptor ligation and PCR. This is a **poly(A)-selected, unstranded** library, not rRNA depletion.
- Sequencing: Illumina HiSeq 4000, paired-end 150 bp.
- Salmon's automatic library type detection agrees: `IU` (inward, unstranded) in every sample checked so far.
- The authors' own processing was Trimmomatic, STAR on GRCh37/hg19 and featureCounts. Only a DEG table with normalised
  counts was deposited, which is why the reads are re-quantified here.

**Reference.**
- GENCODE release 47 (GRCh38), from https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/:
  - `gencode.v47.transcripts.fa.gz`, md5 9c3ce49c6b6f13f3e4f4deb53ca80bf9
  - `gencode.v47.annotation.gtf.gz`, md5 2459f58694fd330740627a8c371fed6a
- Both md5 values match the release MD5SUMS file.

**Index.**
- `salmon index -t gencode.v47.transcripts.fa.gz -i salmon_index_v47_k31 --gencode -k 31 -p 12`.
- Built in 116 s; 1.1 GB on disk at `data/ref/gencode_v47/`.
- **Transcriptome only, without decoys.** The genome is not added as a decoy, to fit a laptop. Reads from unannotated or
  intronic sequence can therefore be assigned to similar transcripts more often than with a decoy-aware index (Srivastava
  et al. 2020). This is a known, small difference from nf-core/rnaseq's default salmon setup.

**Quantification.**
- Command: `salmon quant -i salmon_index_v47_k31 -l A -1 R1 -2 R2 --validateMappings --gcBias -p 10`.
- Salmon 1.10.3 from bioconda, in the conda env `chorus-salmon` (osx-64 under Rosetta).
- Samples are named by the ENA `sample_alias`, which is the GSM accession, and these match `sample_id` in
  `data/meta/GSE224615_meta.csv`.

**Gene level.**
- `NumReads` is summed per gene with a transcript-to-gene table built from the GENCODE v47 GTF
  (`tx2gene_gencode_v47.tsv`: transcript_id, gene_id, gene_name, gene_type).
- Gene ids have the version removed.
- Output: `data/quant/GSE224615/salmon.merged.gene_counts.tsv`, estimated counts (genes x samples). The loader rounds
  them to integers and keeps provenance `requantify`.

**QC.**
- `data/quant/GSE224615/salmon_qc.tsv` gives, per sample:
  - processed and mapped fragments and the mapping rate;
  - the mitochondrial fraction (share of NumReads on chrM transcripts);
  - the detected library type;
  - the cohort median and MAD of the mapping rate.
- A sample is **flagged** when its mapping rate is below 40%, or more than 2.5 MAD (raw, unscaled) below the cohort median.
- Flags are not exclusions. Every sample enters the benchmark, and flagged samples are listed in `results/README.md`
  so that excluding them can be a later sensitivity run.
- With a poly(A) library and a transcriptome-only index, low mapping rates mostly reflect intronic, intergenic and
  mitochondrial reads. In the first sample (GSM7027478, 52.5% mapped), 58% of the quantified reads were on chrM
  transcripts, against 4.5% in the second.

**QC summary (all 36 runs, 12 November 2026).**
- 36 of 36 runs quantified; no run failed. Library type `IU` in all 36.
- Gene matrix: 78,277 genes x 36 samples.
- Mapping rate: cohort median 81.57%, MAD 1.77 (raw), so the MAD lower bound is 77.13%. No sample is below 40%.
- Mitochondrial fraction: median 0.045 (range 0.027 to 0.580).
- Flagged (more than 2.5 MAD below the median), all kept in the benchmark:

  | sample | mapping rate (%) | mitochondrial fraction | processed fragments |
  |---|---|---|---|
  | GSM7027478 | 52.53 | 0.580 | 28,651,080 |
  | GSM7027490 | 59.07 | 0.457 | 16,969,427 |
  | GSM7027501 | 56.84 | 0.052 | 76,655,739 |
- Two downloads needed a retry (SRR23348807: an incomplete byte range; SRR23348806: an md5 mismatch); both passed the md5
  check on the next attempt.

## Environments

| purpose | interpreter | key versions |
|---|---|---|
| curation, smoke tests | `.venv` (miniconda python 3.12.9, x86_64 under Rosetta) | numpy 2.5.4, pandas 3.1.1, scikit-learn 1.9.1 |
| full benchmark runs | `~/.venvs/chorus-arm64` (Homebrew python 3.12.16, arm64), `pip install -e ".[all,dev]"` | numpy 2.5.4, pandas 3.1.1, scikit-learn 1.9.1, torch 2.15.0, tabpfn 9.2.0, inmoose 0.9.1, xgboost 3.4.1 |
| salmon | conda env `chorus-salmon` | salmon 1.10.3 |

**Benchmarks run in the arm64 environment.** Every benchmark in `results/` was run with
`~/.venvs/chorus-arm64/bin/python`, called through `arch -arm64`: the three real-data runs, both simulation runs, and
the smoke tests from 2026-11-11 on. The exact package versions are pinned in `docs/requirements-chorus-arm64.txt` (`pip
freeze`, editable install of this package excluded). To recreate the environment:
`python3.12 -m venv ~/.venvs/chorus-arm64 && ~/.venvs/chorus-arm64/bin/pip install -r docs/requirements-chorus-arm64.txt && ~/.venvs/chorus-arm64/bin/pip install -e .`
on an arm64 Python 3.12. Earlier smoke tests and all curation scripts used `.venv`.

**Why a second environment.** The project `.venv` is an x86_64 Python. PyTorch publishes no x86_64 macOS wheels after
2.2.2, which is too old for TabPFN (it fails with `module 'torch' has no attribute 'OutOfMemoryError'`). That torch
build is also compiled against NumPy 1.x. The arm64 environment has the same numpy, pandas and scikit-learn versions as
`.venv`, and the 16 tests pass in both.

## Proteomics cohorts (downloaded 2026-11-14)

`scripts/proteomics/download.sh` fetches every file below into `data/proteomics/<study>/raw/` (local, not committed)
and writes `data/proteomics/downloads.tsv`: URL, size, the checksum the repository publishes, the checksum observed, and
the SHA-256 of each file. Every published checksum matched. The repositories published:
- Zenodo: MD5;
- Mendeley Data: SHA-256;
- PRIDE: SHA-1 for PXD036969, PXD036970 and PXD036590, none for PXD047623, so only size and SHA-256 are recorded there;
- GitHub: commit hash.

| study | source | accession | licence | files taken | checked |
|---|---|---|---|---|---|
| Buggert2025 (Gao 2025 Nat Immunol) | Zenodo | 10.5281/zenodo.14772494 (published 2025-01-30) | CC-BY-4.0 | `plate1_olink.csv`, `plate2_olink.csv`, `plate3_olink.csv`, `meta_UKsamples.csv`, `meta_SEsamples.csv` (all 5) | MD5 |
| Zoodsma2022 (Front Immunol) | GitHub | CiiM-Bioinformatics-group/Postcovid-targeted-proteomics, commit d45e36077d735c242788e72b871ce4a3d1e728a9 | MIT | whole repository (`git clone`, 110 MB) | commit |
| Cervia2024 (Cervia-Hasler 2024 Science) | Mendeley Data | 10.17632/dvf6yvrg4x.1 (published 2024-01-22) | CC BY 4.0 | the 5 proteomics files: `..._Proteomics_Clinical_Data_6M_timepoint.xlsx`, `..._Acute_plus_healthy-81EAkK.xlsx`, `..._Labels-TAklSO.xlsx`, `..._SomaScan_Aptamer_Info-KnqtVJ.xlsx`, `..._Mass_Spectrometry_6M_timepoint-VOZ1Hu.xlsx`. Not taken: scRNA-seq (1.4 GB rds), flow cytometry, VirScan, images | SHA-256 |
| Su2022 (Cell, INCOV) | Mendeley Data | 10.17632/96v329bg7g.1 (published 2022-01-17) | CC BY 4.0 | `Table S1.xlsx` (clinical), `Table S2.xlsx` (proteomics, metabolomics). Not taken: Tables S3 to S7 | SHA-256 |
| Kovarik2023 (iScience) | PRIDE | PXD036969 (plasma), PXD036970 (macrophage supernatants), published 2023-03-10 | CC0 | `MQ_Search.zip` of each (MaxQuant search output). Not taken: RAW, `MQ_Peak.zip` | SHA-1 |
| Captur2022 (eBioMedicine) | PRIDE | PXD036590 (published 2022-10-11) | CC0 | 8 `Plate_*_TransitionResults_QuantQuali_YE_IS.csv`, `Supplementary_Data_1.xlsx`. Not taken: 32 RAW, two docx | SHA-1 |
| Chowdhury2024 (Front Physiol) | PRIDE | PXD047623 (published 2024-10-17) | CC0 | `diann.7z` (DIA-NN output, 1,302,826,886 bytes, equal to the PRIDE size), `metadata.sdrf.tsv`, `PXD047623_community_annotated.sdrf.tsv`. Not taken: 52 raw `.d` archives (349 GB) | size, SHA-256 |

- **Chowdhury2024 extraction:** from `diann.7z` only the DIA-NN summary tables were extracted (`report.pg_matrix.tsv`,
  `report.pr_matrix.tsv`, `report.gg_matrix.tsv`, `report.unique_genes_matrix.tsv`, `report.stats.tsv`,
  `report.log.txt`), with py7zr 1.1.4 in a separate virtual environment. 7-Zip is not installed.
- **Vijayakumar2022 (Dryad 10.5061/dryad.2ngf1vhq3, CC0) was not downloaded.**
  - The Dryad file API answers 401 without a token, and the landing-page download link answers 403 to scripted
    requests.
  - The dataset record lists 4 files with SHA-256: `plasma_long_personal_clinical_flow.csv` 031ffed7…aa82,
    `bal_long_personal_clinical_flow.csv` 5ddc5bf6…ba5, `data_description.xlsx` f2b26abe…2d98, `README.md` c7c49eb3…c01e.
    The full hashes are in the record.
  - These files can be checked once they are downloaded through a browser.
- **Paper text** for the curation was read from Europe PMC full text where available (for example PMC12043503 Gao 2025,
  PMC9670186 Zoodsma 2022, PMC9719844 Kovarik 2023). Cervia-Hasler 2024 is not in Europe PMC, so its methods were read
  on the Science site. Kovarik 2023 Table S1 (mmc2.xlsx) was read from the Europe
  PMC supplementary files to map MaxQuant sample codes to groups; it is not stored in the repository.
- **Kovarik 2023 supplement (downloaded 2026-11-15)** for Test 2: `mmc2.xlsx` (Table S1, proteomics), `mmc3.xlsx`
  (Table S2, eicosanoids, nAUC), `mmc4.xlsx` (Table S3, Biocrates MxP Quant 500, log2 µM), `mmc5.xlsx` (Table S4,
  triglyceride heatmap values) and the record `PMC9719844.2.json`, from the PMC Open Access dataset on AWS
  (`https://pmc-oa-opendata.s3.amazonaws.com/PMC9719844.2/`, article version 2, CC BY), into
  `data/proteomics/Kovarik2023_plasma/raw/supplement/` (local). The MD5 of each file equals the one in the record:
  mmc2 967c9d2b168cc807c9feaa6276a75790, mmc3 42686e41f05ac4edf3737ee19a51217a, mmc4 d5468dba9c7df74d94e9c70ffd7e079e,
  mmc5 d71a151d7e8f2895b74fc083be07e216. The Europe PMC supplementary files service first answered 503 and then broke
  off the transfer, so it was not used. Why the tables were not used is in `docs/paired_cohorts.md`.
- **Su2022 metabolomics:** sheet "S2.2 Metabolomics" of the `Table S2.xlsx` above (same file, same SHA-256), built by
  `scripts/proteomics/build_Su2022_metabolomics.py` into `data/metabolomics/Su2022_metabolomics/` (local).
- **R for DIABLO:** R 4.5.3 (`/usr/local/bin/Rscript`) with mixOmics 6.32.0 from Bioconductor, installed on
  2026-11-15 into the separate library `~/Library/R/chorus-mixomics` (with BiocManager and its dependencies), so the
  system R library is unchanged. `chorus.eval.paired` reads that path from `CHORUS_R_LIB` (default as above).
- **Reference:** `data/ref/hgnc_complete_set.txt` is the HGNC complete set already in use. It maps symbols to UniProt for
  Su2022 and Captur2022.
- **MSigDB gene sets (downloaded 2026-11-15)** for Test 3 (`docs/mosaic.md`): release 2026.2 human, from
  `https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2026.2.Hs/` into `data/ref/msigdb/` (local), under the
  Creative Commons Attribution 4.0 International licence (MSigDB v2022.1 and later; the extra terms of the licence page
  concern KEGG, BioCarta and AAAS/STKE sets only). MD5: `h.all.v2026.2.Hs.entrez.gmt` f945c99da9f9c996387a1e0f9a88a00e,
  `h.all.v2026.2.Hs.symbols.gmt` c945afd59d949398471c8768fe4d41b0, `c2.cp.reactome.v2026.2.Hs.entrez.gmt`
  6db9a240e0edd0171bd7c349da51a895, `c2.cp.reactome.v2026.2.Hs.symbols.gmt` 74705d0ea35325b3d56a3412801e5592. The
  Entrez files are used, mapped through the HGNC complete set above.

## Metabolomics cohorts (downloaded 2026-11-15 and 2026-11-16)

For the metabolomics arm (`docs/metabolomics_cohorts.md`), into `data/metabolomics/<study>/raw/` (local). The raw MS
files were not downloaded.
- **Martinez2024:** Metabolomics Workbench ST003103 (project PR001926), licence CC BY-NC-ND 4.0 as stated by the
  Workbench study record. The mwTab text of the five analyses was taken from the REST service
  (`/rest/study/analysis_id/<AN>/mwtab/txt`), plus the sample factors (`/rest/study/study_id/ST003103/factors`). The
  Workbench publishes no checksum; SHA-256 of the files as downloaded: AN005077
  964546da90dc5d12a789048ba8ccb6b11694d70c5cf4d346254361cd067da705, AN005078
  cca14df184de90c7a769e64f0b974a5af803367f1eb3610b49460eeb7e1bd510, AN005079
  7165c6161393adca57549c3f8aa9d6047ddfe5e1769ed409b6f48412be9759a5, AN005080
  feea1bb94b67ce60e4e6f3f06406a4a9b6682ecdea609c81e294068de53ebe80, AN005081
  2fd70f4128b6e312bd5c7f90db95c2f849c24a0b8425d09a2131a0af9cc633e5, factors
  151d34e3a25850acce9b055575b09d89c1cabb2a684353581a72d40df58144ce.
- **Lopez2023:** Mendeley Data 10.17632/8zfdjsypd8.1, CC BY 4.0, `DATA AVAILABILITY.xls` (the only file). Its
  SHA-256 equals the one in the Mendeley file record: 5554301ca0e8a893215e26ed340f3480f5ac431c0bdbe7b86c26d1811b6e558c.
  Read with xlrd 2.0.2, installed into `~/.venvs/chorus-arm64` for this file (`docs/requirements-chorus-arm64.txt`).
- **Holmes2021:** Zenodo 10.5281/zenodo.4713391 (concept 10.5281/zenodo.4713390, the accession 4713390 of the
  catalogue), CC BY 4.0. The README, `sampleDescriptionList.txt` and the four `targetedMS_*.txt` tables were downloaded;
  the MD5 of each equals the Zenodo record (`record_4713391.json`). The 32 mzML files were not downloaded.
- **Saito2024:** supplementary material of Saito S et al. 2024 Front Immunol (PMC10830702), CC BY 4.0, from the Europe
  PMC supplementary files service. `DataSheet_1.pdf` (Supplementary Table 1) SHA-256
  2378fc7d20742eb444b827afc99b304f62d9a6d182bba79bf38e8628caf34019, `DataSheet_2.xlsx` (Supplementary Table 2)
  361588d72fef930dfdf0cc2128d7ecd8ee429301058b1502eac13cf237631e14. No publisher checksum is available.
- **Mikuteit2023:** supplementary material of Mikuteit M et al. 2023 Amino Acids (PMC10564820) from Europe PMC:
  `726_2023_3305_MOESM1_ESM.docx` (SHA-256 df7e04a85fce35d42ab72df1297745e6319aa16ac5777de36c8fbc31e938ec48) holds QC
  results and group summaries only, so no cohort was built and the file is not stored.
- **RefMet:** the full RefMet table from `https://www.metabolomicsworkbench.org/databases/refmet/refmet_download.php`,
  saved as `data/ref/refmet/refmet.csv` (208,173 entries; SHA-256
  e2c630f9e9c5a5378eb8687f0bdfe02e4e93f5020c6c6226aa2ee72ce1ed2952). Name matching uses the RefMet match REST service;
  its answers are cached in `data/ref/refmet/refmet_match_cache.csv` with the query date.

## Cohorts added 16 November 2026

Downloaded 2026-11-16 into `data/` (local, not committed). Definitions are in `docs/proteomics_cohorts.md` and
`docs/count_matrices.md` ("Cohorts added 16 November 2026").

| study | source | accession | licence | files taken | checked |
|---|---|---|---|---|---|
| Chowdhury2026 (Sci Rep, PMID 41998026) | PRIDE | PXD066724 (released 2026-03-24) | CC0 (PRIDE) | `diann.7z` (1,498,231,135 bytes, equal to the PRIDE size), `metadata.sdrf.tsv`, `PXD066724_community_annotated.sdrf.tsv`. Not taken: the raw `.d.7z` archives | size, SHA-256 (PRIDE publishes no checksum for this project) |
| Fineschi2026 (Front Immunol, PMID 41808838) | SciLifeLab figshare | 10.17044/scilifelab.31316779 (2026-02-13) | CC BY 4.0 | all 4 files: `XL-4180_raw_data.csv`, `XL-4180_Analysis_Report_2025-06-04.pdf`, `README.txt`, `MANIFEST.txt` | MD5 (figshare) |
| Fineschi2025 PBMC labels (Front Immunol, PMID 40519915) | SciLifeLab figshare | 10.17044/scilifelab.28832492 (2025-04-29) | CC BY 4.0 | `SamplesPC.txt`, `README.txt`, `MANIFEST.txt`, into `data/counts/Fineschi2025_check/`, for the overlap check only. Not taken: `salmon.merged.gene_counts.tsv` | MD5 (figshare) |
| Huapaya2024 (Front Immunol, PMID 39324144) | NHLBI figshare | 10.25444/nhlbi.26383798 (2025-07-22) | public domain | `V2 SOMA COVID_project_data.xlsx` (no per sample PASC label; not built) | MD5 (figshare) |
| GSE297499 (J Infect Dis, PMID 40920586) | GEO | GSE297499 (public 2026-01-21) | GEO | `GSE297499_final_processed_tpm_values_pbmc.txt.gz`, `GSE297499_series_matrix.txt.gz`. Not taken: the CSF TPM file and the reads (SRA PRJNA1265093) | SHA-256 recorded (GEO publishes none) |
| PRJNA1149196 (Nat Microbiol, PMID 39478108) | ENA / SRA | PRJNA1149196 | INSDC | leading bytes of the FASTQ files of all 101 runs and the whole FASTQ files of the 9 QLD runs (below); BioSample records (`efetch db=biosample`) as `data/meta/raw_PRJNA1149196/biosamples.xml` and `.csv` | see below |

- **Chowdhury2026 extraction.** From `diann.7z` only the DIA-NN summary tables (`report.pg_matrix.tsv`,
  `report.pr_matrix.tsv`, `report.gg_matrix.tsv`, `report.unique_genes_matrix.tsv`, `report.stats.tsv`,
  `report.log.txt`) were extracted, with py7zr 1.1.4 in a separate virtual environment, as for Chowdhury2024.
- **SHA-256 of the files as downloaded:**
  - `diann.7z` 8529344794ec95eff511cf8e49c358a2441397268dc5c1eb020af0fcadb30603; `metadata.sdrf.tsv`
    9c562822260554cfcacc904b6018062aa1b5e5db28ca62eecbcc288563ab307f; `PXD066724_community_annotated.sdrf.tsv`
    353f448e4b1f1c9ed80020ab802a3472eda9a21eb5f69a9a4fc9f97ba3aea7ad.
  - `XL-4180_raw_data.csv` 44c04cff0df7661c5107e4cefcb025ba50a5ecf705ccbdddd924f4900d065e5b (MD5
    bf3536286e0d31ea3b0df82183ad1bad, equal to figshare); `SamplesPC.txt`
    7d22f7581ef2ea8f9196b62f3ebfdeac2e3007fdbb8404b3ff95634b8a18b4f9 (MD5 d9e2c54105686e68d6d0ac59c51534ef, equal);
    `V2 SOMA COVID_project_data.xlsx` 76465868c305cd65088763dc426486c0bdcf5f056a54da8dd7f01d002841211f (MD5
    721ff0a0fa5818098a6dd7f320fca432, equal). The README, MANIFEST and PDF MD5s also equal the figshare values.
  - `GSE297499_final_processed_tpm_values_pbmc.txt.gz` 37b6c414a64cac36a070bc46ed7494f431f376e8ec1e4e6b3cfce3390010e290;
    `GSE297499_series_matrix.txt.gz` 746b14e48a6302156f9bcc8f25c7219fa563db52588abade668600878ce7d982.
  - `biosamples.xml` 05e8530b6a3028c186a429f5f7e2b9621c1df55ad44b97f57c66831da12d8bc0.
- **GSE297499 matrix.** The PBMC TPM file is keyed by gene symbol (56,609 rows). `scripts/build_rna_cohorts_2026_11.py`
  maps symbols to Ensembl with the HGNC table (approved symbol, then unique previous symbol; 36,599 rows mapped, the
  others are mostly StringTie or unannotated ids) and sums rows of one gene, giving
  `data/counts/GSE297499/GSE297499_pbmc_tpm_ensembl.tsv.gz` (36,499 genes x 26 samples).

## PRJNA1149196 (Sinclair 2024): SA on the first 5,000,000 read pairs of each run, QLD on all reads

**First pass, all 101 runs on a read prefix** (`scripts/salmon_quant_ena.py PRJNA1149196
data/ref/gencode_v47/salmon_index_v47_k31 data/quant/PRJNA1149196 --max-pairs 5000000 --threads 8 --segments 16`,
salmon 1.10.3 from the conda env `chorus-salmon`, 16 November 2026).
- Runs: the ENA file report (`data/quant/PRJNA1149196/ena_filereport.tsv`), 101 paired-end runs, 16.2 M to 90.7 M read
  pairs per run (median 49.2 M); 1,113.6 GB of gzipped FASTQ in total.
- For each mate file only the leading bytes expected to hold 5,000,000 reads (from the ENA `fastq_bytes` and
  `read_count`, plus 15%) are downloaded in 16 byte-range segments: 129.3 GB in total.
- The prefixes are decompressed, cut to 5,000,000 complete records and checked: both mates keep the same number of
  records and the read names of the last record agree. A partial file cannot be checked against `fastq_md5`.
- Index and quantification options as for GSE224615 (GENCODE v47 transcripts, k 31; `salmon quant -l A --validateMappings --gcBias`).
- Reads used per run: `data/quant/PRJNA1149196/prefix_reads.tsv` (5,000,000 pairs in all 101 runs).
- **Failed runs and retries.** The first pass (log `quant.log`) quantified 98 runs. Three failed:
  - SRR30515205 (23-00316) and SRR30515189 (23-00330): a download segment was still incomplete after 10 attempts.
  - SRR30565305 (23-00343): the last read names differed between mates (`@SRR30565305.5000000` and `+`), so the run
    was not quantified.
  The restart (`quant_retry.log`) skips runs that already have `quant.sf`; on the first retry all three were downloaded
  again and quantified (mapping 39.7%, 41.1% and 35.0%). 101 of 101 runs quantified.
- Aggregation (`--aggregate --gtf data/ref/gencode_v47/gencode.v47.annotation.gtf.gz`): 78,277 genes x 101 samples in
  `salmon.merged.gene_counts.tsv`, QC in `salmon_qc.tsv`.

**The QLD prefix is not a sample of the reads.** In the prefix counts of the 9 QLD runs, chr1 genes hold 75.8% to 95.7%
of the counts and chr10 most of the rest (median 11.3%); in SA chr1 holds 8.7% to 10.8%. The QLD FASTQ files are
ordered by alignment position, so their first reads come from the first chromosomes in lexical order. Only 5,914 genes
passed the expression filter (CPM at least 1 in at least half of the samples) in the QLD prefix counts, against 18,004
in SA. The SA prefix counts are used: their chr1 share (8.7% to 10.8%) is in the range of the QLD counts on all reads
(7.9% to 10.6%, below).

**QLD on all reads** (`scripts/salmon_quant_ena.py PRJNA1149196 data/ref/gencode_v47/salmon_index_v47_k31
data/quant/PRJNA1149196_QLD_full --runs SRR30280882,...,SRR30280890 --threads 8 --segments 16`, 16 and 17 November
2026, log `quant.log`).
- The 9 QLD runs downloaded whole (27.4 GB), each file checked against the ENA `fastq_md5`; 9 of 9 quantified on the
  first attempt (44.5 M to 53.4 M fragments).
- Aggregation as above: 78,277 genes x 9 samples, `data/quant/PRJNA1149196_QLD_full/salmon.merged.gene_counts.tsv`
  and `salmon_qc.tsv`. chr1 holds 7.9% to 10.6% of the counts; 16,053 genes pass the expression filter.
- `scripts/build_rna_cohorts_2026_11.py` writes `PRJNA1149196_SA_gene_counts.tsv.gz` (92 samples) from the prefix
  matrix and `PRJNA1149196_QLD_gene_counts.tsv.gz` (9 samples) from the full-depth matrix, both into
  `data/quant/PRJNA1149196/`.
- The version 4 transcriptomic runs made with the QLD prefix matrix (16 November) were not kept in `results/`; all
  version 4 transcriptomic results use the full-depth QLD matrix.

**QC.**

| site | runs | reads quantified | library type | mapping rate median (range) | mitochondrial fraction median (range) | flagged |
|---|---|---|---|---|---|---|
| SA | 92 | first 5 M pairs | ISF (all) | 39.56% (34.72 to 47.05) | 0.0174 (0.0071 to 0.0369) | 56 |
| QLD | 9 | all (44.5 M to 53.4 M) | ISR (all) | 59.37% (50.81 to 67.33) | 0.0072 (0.0046 to 0.0104) | 0 |

- The QC rule is the one used for GSE224615, applied within each aggregation. SA: in the 101-run prefix aggregation,
  median 39.71%, MAD 1.5 (raw), lower bound 35.96%. QLD: in the full-depth aggregation, median 59.4%, MAD 4.2, lower
  bound 48.7%.
- 56 SA samples are flagged: 51 below 40% only, and 5 below 40% and more than 2.5 MAD below the median (23-00304
  35.65%, 23-00343 35.04%, 23-00408 34.72%, 23-00409 35.07%, 23-00414 35.95%). The flagged ids are listed by
  `aggregate.log` and in `salmon_qc.tsv` (`flag`, `flag_reason`). No QLD sample is flagged.
- As for GSE224615, flags are not exclusions: every sample enters the benchmark.

## Cohorts added 17 November 2026

Downloaded 2026-11-17 into `data/` (local, not committed). Definitions and the counts reported before building are in
`docs/proteomics_cohorts.md` and `docs/count_matrices.md` ("Cohorts added 17 November 2026" and "GSE228320"). The two
Zenodo files are fetched by `scripts/proteomics/download.sh` (run as `ONLY="Woodruff2023 Bansal2025"`; `ONLY` limits a
run to the named studies), which logs them in `data/proteomics/downloads.tsv`.

| study | source | accession | licence | files taken | checked |
|---|---|---|---|---|---|
| Woodruff2023 (Nat Commun, PMID 37452024) | Zenodo | 10.5281/zenodo.8092298 (2023-06-28) | CC BY 4.0 | `Covid-v0.4.0.tar.gz` (26,071,732 bytes); `data/complete.csv` extracted from it (34,243,700 bytes) | MD5 (Zenodo) |
| Bansal2025 (Commun Med, PMID 41975235) | Zenodo | 10.5281/zenodo.15237087, files of version 15237088 (2025-04-17) | CC BY 4.0 | `baseline57.xlsx` (1,234,901 bytes). Not taken: `Longitudinal_paired34.xlsx`, `baseline_182proteins_univariate.xlsx`, `rcode_longcovid_2024.R` | MD5 (Zenodo) |
| GSE228320 (Biomed Pharmacother, PMID 36058144) | GEO | GSE228320 | GEO | `GSE228320_Raw_counts.txt.gz` (1,763,745 bytes), `GSE228320_series_matrix.txt.gz` (5,899 bytes), into `data/counts/GSE228320/`. Not taken: `GSE228320_Counts_Normalized_postCOVID.txt.gz`, the reads | SHA-256 recorded (GEO publishes none) |

- **Checksums of the files as downloaded:**

  | file | MD5 | SHA-256 |
  |---|---|---|
  | `Covid-v0.4.0.tar.gz` | e2ad5be6185f7f14d9bf6d22c85a3730, equal to Zenodo | 723d12248a60c8e7c73b3d5d3477dea2d1ff633f2618a50eb016b58f67bcc1c6 |
  | `complete.csv` (extracted) | not recorded | 1d7537e9a90372e31b6d09fc52797ecff450bc62f5e32151a44ebc0cc80256a8 |
  | `baseline57.xlsx` | d46a796a46c6ed883e5eac8adac249f7, equal to Zenodo | d19681bf7c6297b0978101bdb314e9cebe4b7af4b6b307e83b17cbec0481246a |
  | `GSE228320_Raw_counts.txt.gz` | not recorded | 407226d486ddbbf81cb1cd5bba886b0fca1cbcfaa8748799c12533b64d929642 |
  | `GSE228320_series_matrix.txt.gz` | not recorded | 0a8d3d81780a65e506f3fadf342d9d725b1a85e0bf45fcb820735e8c3d066d2f |

- **Inspected, not downloaded into `data/`:**
  - PRIDE PXD051740 and PXD041762: project records and file lists from the PRIDE API (v3). For PXD041762, the first
    3 MB of `REC-2349_Russell_v2.sky` were read with an HTTP range request.
  - The Berezhnoy 2023 supplement (Europe PMC, PMC10203989). It and the partial Skyline file were kept in the session
    scratch space only.
