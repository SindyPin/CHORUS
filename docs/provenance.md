# Provenance of re-quantified matrices and the run environment

## GSE224615 (Yin 2024, LIINC PBMC): salmon on a laptop

**Source reads.**
- ENA project PRJNA931952: 36 paired-end runs, about 230 GB of FASTQ.
- Downloaded over HTTPS from the ENA `fastq_ftp` links. Each file was checked against the ENA `fastq_md5` before use and
  deleted after quantification.
- ENA limits a single connection to about 1 MB/s and answers HTTP 403 when loaded. `scripts/salmon_quant_ena.py`
  therefore downloads each file in 12 parallel byte ranges, retrying and resuming each range with back-off.
- The ENA file report used is saved as `data/quant/GSE224615/ena_filereport.tsv`.

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
- `data/quant/GSE224615/salmon_qc.tsv` gives the mapping rate per sample and flags samples below 60%.
- The summary is below.

QC summary: _filled in when the run completes._

## Environments

| purpose | interpreter | key versions |
|---|---|---|
| curation, smoke tests | `.venv` (miniconda python 3.12.9, x86_64 under Rosetta) | numpy 2.5.4, pandas 3.1.1, scikit-learn 1.9.1 |
| full benchmark runs | `~/.venvs/chorus-arm64` (Homebrew python 3.12.16, arm64), `pip install -e ".[all,dev]"` | numpy 2.5.4, pandas 3.1.1, scikit-learn 1.9.1, torch 2.15.0, tabpfn 9.2.0, inmoose 0.9.1, xgboost 3.4.1 |
| salmon | conda env `chorus-salmon` | salmon 1.10.3 |

**Why a second environment.** The project `.venv` is an x86_64 Python. PyTorch publishes no x86_64 macOS wheels after
2.2.2, which is too old for TabPFN (it fails with `module 'torch' has no attribute 'OutOfMemoryError'`). That torch
build is also compiled against NumPy 1.x. The arm64 environment has the same numpy, pandas and scikit-learn versions as
`.venv`, and the 16 tests pass in both.
