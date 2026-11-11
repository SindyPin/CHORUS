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
