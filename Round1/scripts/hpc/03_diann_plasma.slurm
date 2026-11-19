#!/bin/bash
#SBATCH --job-name=chorus_diann
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --output=logs/diann_%j.out
# Uniform re-search of one PRIDE DIA plasma proteomics dataset with DIA-NN (library free, one FASTA, 1% FDR).
# For DDA datasets use FragPipe (MSFragger + IonQuant) with the same FASTA; or run both through nf-core/quantms.
# Usage: sbatch 03_diann_plasma.slurm PXD026442 /scratch/$USER/chorus
set -euo pipefail
PXD=${1:?PXD}
ROOT=${2:-/scratch/$USER/chorus}
mkdir -p $ROOT/prot/$PXD/raw $ROOT/prot/$PXD/out logs
cd $ROOT/prot/$PXD

# 1. fetch raw files listed by the PRIDE API (requires network on the login/data node; run once)
if [ ! -f files.json ]; then
  curl -s "https://www.ebi.ac.uk/pride/ws/archive/v2/projects/$PXD/files?pageSize=1000" > files.json
  python3 -c "import json;[print(f['publicFileLocations'][0]['value']) for f in json.load(open('files.json')) if f['fileCategory']['value']=='RAW']" > raw_urls.txt
  (cd raw && aria2c -c -x 8 -j 4 -i ../raw_urls.txt)
fi

# 2. one FASTA for every dataset (UniProt human reviewed + isoforms, pinned release; record in provenance)
FASTA=$ROOT/ref/uniprot_human_sp_2025_04.fasta
[ -f $FASTA ] || (mkdir -p $ROOT/ref && curl -s "https://rest.uniprot.org/uniprotkb/stream?query=organism_id:9606+AND+reviewed:true&format=fasta&includeIsoform=true" > $FASTA)

# 3. convert vendor raw to mzML if needed (ThermoRawFileParser container), then DIA-NN
module load singularity 2>/dev/null || true
DIANN=${DIANN:-$ROOT/bin/diann-2.1.0/diann-linux}
$DIANN --dir raw --out out/report.tsv --out-lib out/lib.tsv \
  --fasta $FASTA --fasta-search --predictor --gen-spec-lib \
  --threads $SLURM_CPUS_PER_TASK --qvalue 0.01 --matrices --reanalyse \
  --min-pep-len 7 --max-pep-len 30 --missed-cleavages 1 --unimod4 --var-mods 1 --var-mod UniMod:35,15.994915,M \
  --mass-acc 15 --mass-acc-ms1 15 --window 0 --smart-profiling --pg-level 2

# Outputs used by CHORUS: out/report.pg_matrix.tsv (protein groups x runs, MaxLFQ intensities, NaN = not quantified).
# Record sample to run mapping from the PRIDE sdrf file (files.json: fileCategory EXPERIMENTAL_DESIGN) into metadata.
