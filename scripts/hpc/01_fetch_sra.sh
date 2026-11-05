#!/usr/bin/env bash
# Fetch FASTQ for a GEO series through ENA (faster and more reliable than SRA toolkit on HPC).
# Usage: bash 01_fetch_sra.sh GSE275334 /scratch/$USER/chorus/raw
# Produces <out>/<GSE>/samplesheet.csv ready for nf-core/rnaseq and a metadata skeleton.
set -euo pipefail
GSE=${1:?GSE accession}
OUT=${2:-raw}
mkdir -p "$OUT/$GSE"
cd "$OUT/$GSE"

# 1. GEO series -> BioProject -> ENA run table (includes fastq_ftp links and sample titles)
PRJ=$(curl -s "https://ftp.ncbi.nlm.nih.gov/geo/series/${GSE:0:${#GSE}-3}nnn/${GSE}/soft/${GSE}_family.soft.gz" | zcat | grep -m1 "BioProject" | grep -o "PRJNA[0-9]*" || true)
if [ -z "$PRJ" ]; then echo "no BioProject for $GSE: raw reads may be controlled access (dbGaP). Use the author count matrix and set is_counts: true with provenance author_matrix."; exit 2; fi
echo "$GSE -> $PRJ"
curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${PRJ}&result=read_run&fields=run_accession,sample_accession,sample_title,library_layout,library_strategy,read_count,fastq_ftp,fastq_md5&format=tsv&download=true" > ena_runs.tsv
N=$(($(wc -l < ena_runs.tsv) - 1)); echo "$N runs"

# 2. download with aria2c (fallback wget), verify md5
mkdir -p fastq
tail -n +2 ena_runs.tsv | awk -F'\t' '{print $7}' | tr ';' '\n' | sed 's#^#ftp://#' > urls.txt
if command -v aria2c >/dev/null; then aria2c -c -x 8 -j 4 -d fastq -i urls.txt; else wget -c -P fastq -i urls.txt; fi

# 3. nf-core/rnaseq samplesheet: sample,fastq_1,fastq_2,strandedness
python3 - <<'EOF'
import csv, os
rows = list(csv.DictReader(open("ena_runs.tsv"), delimiter="\t"))
with open("samplesheet.csv", "w", newline="") as fh:
    w = csv.writer(fh); w.writerow(["sample", "fastq_1", "fastq_2", "strandedness"])
    for r in rows:
        files = [os.path.join("fastq", os.path.basename(u)) for u in r["fastq_ftp"].split(";") if u]
        f1 = files[0]; f2 = files[1] if len(files) > 1 else ""
        w.writerow([r["run_accession"], os.path.abspath(f1), os.path.abspath(f2) if f2 else "", "auto"])
with open("metadata_skeleton.csv", "w", newline="") as fh:
    w = csv.writer(fh); w.writerow(["sample_id", "run_accession", "sample_title", "subject_id", "study_id", "lc_status", "case_definition", "sample_matrix", "omics_layer", "days_since_infection", "time_anchor", "acute_severity", "hospitalised", "age", "sex"])
    for r in rows:
        w.writerow([r["run_accession"], r["run_accession"], r["sample_title"], "", os.path.basename(os.getcwd()), "", "", "", "bulk_rna", "", "", "", "", "", ""])
print("wrote samplesheet.csv and metadata_skeleton.csv (fill lc_status etc. from GEO characteristics; see 00_curate_geo_metadata.py)")
EOF
