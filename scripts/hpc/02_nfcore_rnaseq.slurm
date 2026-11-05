#!/bin/bash
#SBATCH --job-name=chorus_rnaseq
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=logs/nfcore_%j.out
# Uniform re-quantification of one cohort with nf-core/rnaseq (Salmon pseudo-alignment, one GENCODE release).
# The Nextflow head job is light; executors are submitted by Nextflow through the Slurm profile.
# Usage: sbatch 02_nfcore_rnaseq.slurm GSE275334 /scratch/$USER/chorus
set -euo pipefail
GSE=${1:?GSE}
ROOT=${2:-/scratch/$USER/chorus}
module load nextflow singularity 2>/dev/null || module load Nextflow Apptainer
export NXF_SINGULARITY_CACHEDIR=$ROOT/singularity
mkdir -p $ROOT/quant/$GSE logs

# Same reference for every cohort: GRCh38 primary assembly + GENCODE v47 (pin once, record in provenance)
REF=$ROOT/ref
mkdir -p $REF
[ -f $REF/GRCh38.primary_assembly.genome.fa ] || (cd $REF && wget -q https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/GRCh38.primary_assembly.genome.fa.gz && gunzip GRCh38.primary_assembly.genome.fa.gz)
[ -f $REF/gencode.v47.primary_assembly.annotation.gtf ] || (cd $REF && wget -q https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.primary_assembly.annotation.gtf.gz && gunzip gencode.v47.primary_assembly.annotation.gtf.gz)

nextflow run nf-core/rnaseq -r 3.19.0 \
  -profile singularity,slurm \
  --input $ROOT/raw/$GSE/samplesheet.csv \
  --outdir $ROOT/quant/$GSE \
  --fasta $REF/GRCh38.primary_assembly.genome.fa \
  --gtf $REF/gencode.v47.primary_assembly.annotation.gtf \
  --gencode \
  --pseudo_aligner salmon --skip_alignment \
  --save_reference \
  --max_cpus 16 --max_memory 64.GB \
  -resume

# Outputs used by CHORUS: $ROOT/quant/$GSE/salmon/salmon.merged.gene_counts.tsv (gene x sample, integer counts)
# and salmon.merged.gene_lengths.tsv; MultiQC report for QC thresholds (mapping rate, rRNA, globin).
