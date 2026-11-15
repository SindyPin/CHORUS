"""Quantify an ENA project with salmon on a laptop, one run at a time, then aggregate to genes.

For each run in the ENA filereport the FASTQ files are downloaded over HTTPS from the fastq_ftp
links (in parallel byte ranges), checked against fastq_md5, quantified with ``salmon quant --validateMappings --gcBias -l A``
and deleted, so disk use stays at about one run (the next run downloads while the current one is
quantified). Finished runs (quant.sf present) are skipped, so the script can be restarted.
Samples are named by the ENA sample_alias (the GSM for GEO submissions).

``--max-pairs N`` quantifies only the first N read pairs (spots) of each run, for projects too large to download
at full depth: only the leading bytes of each FASTQ file that hold about N reads (from the ENA fastq_bytes and
read_count, with a 15 percent margin) are fetched, decompressed, and cut to the same number of complete records in
both mates. A partial file cannot be checked against fastq_md5, so the check is that the prefix decompresses, that
both mates keep the same number of records and that the read names of the last record agree. The reads kept per
run are written to ``prefix_reads.tsv``. A run whose files are smaller than the cap is downloaded whole and md5
checked as usual.

``--aggregate`` sums NumReads per gene with a tx2gene table built from the GENCODE GTF and writes
genes x samples counts (Ensembl gene ids without version) plus a QC table with the mapping rate
from each run's aux_info/meta_info.json and the mitochondrial (chrM) share of NumReads. A sample is
flagged when its mapping rate is below 40 percent or more than 2.5 MAD (raw, unscaled) below the cohort
median; flags are for sensitivity analyses, not automatic exclusion.

  python scripts/salmon_quant_ena.py PRJNA931952 data/ref/gencode_v47/salmon_index_v47_k31 data/quant/GSE224615
  python scripts/salmon_quant_ena.py PRJNA931952 ... data/quant/GSE224615 --aggregate --gtf data/ref/gencode_v47/gencode.v47.annotation.gtf.gz
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import queue
import re
import subprocess
import sys
import threading
import time
import urllib.request
from pathlib import Path

import pandas as pd

FIELDS = "run_accession,sample_alias,library_layout,read_count,fastq_bytes,fastq_md5,fastq_ftp"
MIN_MAPPING_RATE = 40.0  # flag below this, or more than MAD_FLAG median absolute deviations below the cohort median
MAD_FLAG = 2.5


def filereport(project: str) -> pd.DataFrame:
    url = f"https://www.ebi.ac.uk/ena/portal/api/filereport?accession={project}&result=read_run&fields={FIELDS}&format=tsv&limit=0"
    with urllib.request.urlopen(url, timeout=120) as r:
        from io import StringIO

        return pd.read_csv(StringIO(r.read().decode()), sep="\t").sort_values("sample_alias")


def md5sum(path: Path) -> str:
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 22), b""):
            h.update(chunk)
    return h.hexdigest()


def _segment(url: str, part: Path, start: int, end: int, attempts: int = 10) -> None:
    """Fetch bytes start..end into part, resuming after errors; ENA answers 403 under load, so back off."""
    want = end - start + 1
    for i in range(attempts):
        have = part.stat().st_size if part.exists() else 0
        if have == want:
            return
        if have > want:
            part.unlink()
            have = 0
        with open(part, "ab") as fh:  # a failed transfer leaves a valid prefix, resumed on the next attempt
            res = subprocess.run(["curl", "-sSfL", "-r", f"{start + have}-{end}", url], stdout=fh, stderr=subprocess.DEVNULL)
        if res.returncode != 0:
            time.sleep(min(300, 10 * 2 ** i))
    if not part.exists() or part.stat().st_size != want:
        raise RuntimeError(f"segment {part.name} incomplete after {attempts} attempts")


def download(url: str, dest: Path, size: int, segments: int) -> None:
    """Download in parallel byte ranges (ENA throttles single connections to about 1 MB/s)."""
    step = -(-size // segments)
    ranges = [(i * step, min(size, (i + 1) * step) - 1) for i in range(segments) if i * step < size]
    parts = [dest.with_name(f"{dest.name}.part{i}") for i in range(len(ranges))]
    errors = []
    threads = [threading.Thread(target=lambda p=p, r=r: _try(errors, _segment, url, p, *r)) for p, r in zip(parts, ranges)]
    for t in threads:
        t.start()
    for t in threads:
        t.join()
    if errors:
        for p in parts:
            p.unlink(missing_ok=True)
        raise RuntimeError(errors[0])
    with open(dest, "wb") as out:
        for p in parts:
            with open(p, "rb") as fh:
                while chunk := fh.read(1 << 24):
                    out.write(chunk)
            p.unlink()


def _try(errors: list, fn, *args) -> None:
    try:
        fn(*args)
    except Exception as e:
        errors.append(str(e))


def _read_name(line: str) -> str:
    """Read name without the mate suffix (/1, /2) or the comment after the first space."""
    return re.sub(r"/[12]$", "", line.split()[0])


def fetch_prefix(row, tmp: Path, segments: int, max_pairs: int, log: Path) -> list:
    """Leading bytes of each FASTQ file holding about max_pairs reads, decompressed and cut to the same number of
    complete records in both mates (see --max-pairs)."""
    links, sizes = _mates(row.fastq_ftp.split(";"), [int(x) for x in str(row.fastq_bytes).split(";")])
    caps = [min(size, int(1.15 * max_pairs * size / row.read_count) + 1) for size in sizes]
    if all(c == size for c, size in zip(caps, sizes)):
        return fetch(row, tmp, segments)
    gz, fq = [], []
    for link, cap in zip(links, caps):
        dest = tmp / Path(link).name
        download("https://" + link, dest, cap, segments)
        gz.append(dest)
    n_lines = []
    for g in gz:
        out = g.with_suffix("")  # .fastq.gz -> .fastq
        with open(out, "w") as fh:  # a truncated gzip stream decompresses up to the cut; the error at the end is expected
            zc = subprocess.Popen(["gzip", "-dc", str(g)], stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
            subprocess.run(["head", "-n", str(4 * max_pairs)], stdin=zc.stdout, stdout=fh)
            zc.stdout.close()
            zc.wait()
        g.unlink()
        n_lines.append(sum(1 for _ in open(out)))
        fq.append(out)
    n = min(x // 4 for x in n_lines)
    if n == 0:
        raise RuntimeError(f"{row.run_accession}: no complete record in the downloaded prefix")
    names = []
    for f in fq:
        if sum(1 for _ in open(f)) != 4 * n:
            tmpf = f.with_suffix(".cut")
            subprocess.run(["head", "-n", str(4 * n)], stdin=open(f), stdout=open(tmpf, "w"))
            tmpf.replace(f)
        with open(f) as fh:
            last = None
            for i, line in enumerate(fh):
                if i % 4 == 0:
                    last = line
        names.append(_read_name(last))
    if len(set(names)) != 1:
        raise RuntimeError(f"{row.run_accession}: last read names differ between mates: {names}")
    new = not log.exists()
    with open(log, "a") as fh:
        if new:
            fh.write("run_accession\tsample_alias\tread_pairs_in_run\tread_pairs_used\tbytes_downloaded\tbytes_in_run\n")
        fh.write(f"{row.run_accession}\t{row.sample_alias}\t{row.read_count}\t{n}\t{sum(caps)}\t{sum(sizes)}\n")
    return fq


def _mates(links: list, *other: list):
    """For a paired run that ENA lists with an extra unpaired file (<run>.fastq.gz next to <run>_1 and <run>_2), keep
    the two mate files; the same positions of the other lists are kept."""
    keep = [i for i, l in enumerate(links) if re.search(r"_[12]\.fastq\.gz$", l)]
    if len(keep) != 2:
        keep = list(range(len(links)))
    return ([links[i] for i in keep], *([o[i] for i in keep] for o in other))


def fetch(row, tmp: Path, segments: int) -> list:
    """Download and md5 check the FASTQ files of one run; retry each file up to three times."""
    files = []
    for link, md5, size in zip(*_mates(row.fastq_ftp.split(";"), row.fastq_md5.split(";"), str(row.fastq_bytes).split(";"))):
        dest = tmp / Path(link).name
        for attempt in range(3):
            try:
                download("https://" + link, dest, int(size), segments)
            except RuntimeError as e:
                print(f"[retry] {row.run_accession} {dest.name}: {e}", flush=True)
                continue
            if md5sum(dest) == md5:
                break
            print(f"[retry] {row.run_accession} {dest.name}: md5 mismatch", flush=True)
            dest.unlink(missing_ok=True)
        else:
            raise RuntimeError(f"{row.run_accession}: md5 mismatch for {dest.name} after 3 attempts")
        files.append(dest)
    return files


def quantify(project: str, index: Path, out: Path, threads: int, segments: int, max_pairs: int = 0) -> None:
    runs = filereport(project)
    (out / "salmon").mkdir(parents=True, exist_ok=True)
    runs.to_csv(out / "ena_filereport.tsv", sep="\t", index=False)
    todo = [r for r in runs.itertuples() if not (out / "salmon" / r.sample_alias / "quant.sf").exists()]
    print(f"{len(runs)} runs, {len(runs) - len(todo)} already quantified, {len(todo)} to do", flush=True)
    tmp = out / "fastq_tmp"
    tmp.mkdir(exist_ok=True)
    ready: queue.Queue = queue.Queue(maxsize=1)

    def producer():
        for r in todo:
            try:
                files = fetch_prefix(r, tmp, segments, max_pairs, out / "prefix_reads.tsv") if max_pairs else fetch(r, tmp, segments)
                ready.put((r, files, None))
            except Exception as e:  # pass the failure to the consumer
                ready.put((r, [], e))
        ready.put(None)

    threading.Thread(target=producer, daemon=True).start()
    while (item := ready.get()) is not None:
        r, files, err = item
        if err:
            print(f"[fail] {r.run_accession} {r.sample_alias}: {err}", flush=True)
            continue
        t0 = time.time()
        reads = ["-1", str(files[0]), "-2", str(files[1])] if len(files) == 2 else ["-r", str(files[0])]
        dest = out / "salmon" / r.sample_alias
        res = subprocess.run(["salmon", "quant", "-i", str(index), "-l", "A", *reads, "--validateMappings", "--gcBias",
                              "-p", str(threads), "-o", str(dest)], capture_output=True, text=True)
        for f in files:
            f.unlink(missing_ok=True)
        if res.returncode != 0:
            print(f"[fail] {r.run_accession} {r.sample_alias}: salmon exit {res.returncode}\n{res.stderr[-1500:]}", flush=True)
            continue
        rate = json.loads((dest / "aux_info" / "meta_info.json").read_text())["percent_mapped"]
        print(f"[done] {r.run_accession} {r.sample_alias} mapping {rate:.1f}% in {time.time() - t0:.0f} s", flush=True)


def tx2gene(gtf: Path) -> pd.DataFrame:
    rows = []
    with gzip.open(gtf, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t", 8)
            if f[2] != "transcript":
                continue
            a = dict(re.findall(r'(\S+) "([^"]*)"', f[8]))
            rows.append((a["transcript_id"], a["gene_id"], a.get("gene_name", ""), a.get("gene_type", ""), f[0]))
    return pd.DataFrame(rows, columns=["transcript_id", "gene_id", "gene_name", "gene_type", "chrom"])


def aggregate(out: Path, gtf: Path) -> None:
    t2g = tx2gene(gtf)
    t2g.to_csv(out / "tx2gene_gencode_v47.tsv", sep="\t", index=False)
    m = t2g.set_index("transcript_id")["gene_id"].str.split(".").str[0]
    mito = set(t2g.loc[t2g["chrom"] == "chrM", "transcript_id"])
    counts, qc = {}, []
    for d in sorted((out / "salmon").iterdir()):
        if not (d / "quant.sf").exists():
            continue
        q = pd.read_csv(d / "quant.sf", sep="\t", usecols=["Name", "NumReads"])
        genes = q["Name"].map(m)
        if genes.isna().any():
            raise ValueError(f"{d.name}: {genes.isna().sum()} transcripts missing from the tx2gene table")
        counts[d.name] = q.groupby(genes.values)["NumReads"].sum()
        info = json.loads((d / "aux_info" / "meta_info.json").read_text())
        lib = json.loads((d / "lib_format_counts.json").read_text())
        qc.append({"sample": d.name, "num_processed": info["num_processed"], "num_mapped": info["num_mapped"],
                   "mapping_rate": round(info["percent_mapped"], 2),
                   "mito_fraction": round(q.loc[q["Name"].isin(mito), "NumReads"].sum() / q["NumReads"].sum(), 4),
                   "library_type": lib.get("expected_format"), "salmon_version": info.get("salmon_version")})
    X = pd.DataFrame(counts).fillna(0.0)
    X.index.name = "gene_id"
    X.to_csv(out / "salmon.merged.gene_counts.tsv", sep="\t", float_format="%.3f")
    qc = pd.DataFrame(qc)
    med = qc["mapping_rate"].median()
    mad = (qc["mapping_rate"] - med).abs().median()  # raw MAD, not scaled to a normal SD
    low = med - MAD_FLAG * mad
    qc["cohort_median_mapping"], qc["cohort_mad_mapping"], qc["mad_lower_bound"] = round(med, 2), round(mad, 2), round(low, 2)
    reason = []
    for r in qc["mapping_rate"]:
        why = ([f"below {MIN_MAPPING_RATE:g}%"] if r < MIN_MAPPING_RATE else []) + ([f"more than {MAD_FLAG} MAD below median"] if r < low else [])
        reason.append("; ".join(why))
    qc["flag"] = [bool(x) for x in reason]
    qc["flag_reason"] = reason
    qc.to_csv(out / "salmon_qc.tsv", sep="\t", index=False)
    print(f"{X.shape[0]} genes x {X.shape[1]} samples; mapping rate median {med:.1f}%, MAD {mad:.1f}, lower bound {low:.1f}%; "
          f"flagged: {qc.loc[qc['flag'], 'sample'].tolist() or 'none'}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("project")
    ap.add_argument("index")
    ap.add_argument("out")
    ap.add_argument("--threads", type=int, default=10)
    ap.add_argument("--segments", type=int, default=16, help="parallel byte ranges per FASTQ download")
    ap.add_argument("--max-pairs", type=int, default=0, help="quantify only the first N read pairs of each run (0: all)")
    ap.add_argument("--aggregate", action="store_true")
    ap.add_argument("--gtf")
    a = ap.parse_args()
    out = Path(a.out)
    if a.aggregate:
        aggregate(out, Path(a.gtf))
    else:
        quantify(a.project, Path(a.index), out, a.threads, a.segments, a.max_pairs)


if __name__ == "__main__":
    sys.exit(main())
