#!/usr/bin/env python3
"""Fetch the long reads each release was built from, at the placements only.

  fetch_long_reads.py DATA [--only RUN] [--limit N]

HPRC publishes, for each release, its HiFi reads aligned to that release, indexed. The address is
the BAM variable of the driver script of each run (`scripts/run_na19240_heldout.sh`,
`scripts/run_hg01943.sh`); nothing is fetched from anywhere else. For every strong missing-allele
placement of a run, and for every other row of the cluster table that belongs to the same event
(`twin_group`), the records that overlap the placement's interval and 100 bases on either side
are read through the index and kept in

  DATA/work/missing_alleles_long_reads/<run>/long_reads.bam   (with .bai)

without base qualities and without the `cs` tag, which this record does not use and which make up
two thirds of the size. Beside it: `regions.bed` (what was asked for), `fetch.json` (address, times,
counts, sizes) and the index of the remote file as samtools cached it. A fetch that is complete
(`fetch.json` exists and names the same regions) is not made again.
"""

import argparse
import hashlib
import json
import subprocess
import sys
import time
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import placements as P  # noqa: E402

PAD = P.PAD  # bases on either side of a placement's interval: every read that covers the window of
# placements.long_reads_at overlaps the stretch asked for here
STRIP = r"""BEGIN { FS = OFS = "\t" } /^@/ { print; next } { $11 = "*"; print }"""

ap = argparse.ArgumentParser()
ap.add_argument("data")
ap.add_argument("--only", help="one run key")
ap.add_argument("--limit", type=int, help="the first N regions only, into a folder of its own (a trial)")
args = ap.parse_args()
DATA = Path(args.data)


def regions_of(table: pd.DataFrame) -> pd.DataFrame:
    """The intervals to fetch: every row of the events that have a strong missing-allele placement."""
    events = set(P.strong_alleles(table)["twin_group"])
    rows = table[table["twin_group"].isin(events)]
    r = pd.DataFrame({"contig": rows["contig"], "start": (rows["start"] - 1 - PAD).clip(lower=0), "end": rows["end"] + PAD})
    r = r.sort_values(["contig", "start", "end"])
    merged = []
    for contig, s, e in r.itertuples(index=False):
        if merged and merged[-1][0] == contig and s <= merged[-1][2]:
            merged[-1][2] = max(merged[-1][2], e)
        else:
            merged.append([contig, int(s), int(e)])
    return pd.DataFrame(merged, columns=["contig", "start", "end"]), len(rows)


for key, _, rel, script in P.RUNS:
    if script is None or (args.only and key != args.only):
        continue
    address = P.bam_address(script)
    out = DATA / "work" / P.WORK / (key if not args.limit else f"{key}.trial{args.limit}")
    out.mkdir(parents=True, exist_ok=True)
    regions, n_rows = regions_of(P.load_table(DATA / "work" / rel))
    if args.limit:
        regions = regions.head(args.limit)
    bed = "".join(f"{c}\t{s}\t{e}\n" for c, s, e in regions.itertuples(index=False))
    stamp = hashlib.md5((address + bed).encode()).hexdigest()
    done = out / "fetch.json"
    if done.is_file() and json.loads(done.read_text()).get("stamp") == stamp and (out / "long_reads.bam.bai").is_file():
        print(f"{key}: already fetched ({json.loads(done.read_text())['records']} records)")
        continue
    (out / "regions.bed").write_text(bed)
    t0 = time.time()
    started = time.strftime("%F %T %Z")
    # -M: every record once, however many regions it overlaps; unmapped records are left out
    cmd = (f"samtools view -h -M -L regions.bed -F 0x4 -x cs '{address}' | awk '{STRIP}' "
           f"| samtools view -b -o long_reads.tmp.bam -")
    with open(out / "fetch.log", "w") as log:
        subprocess.run(["bash", "-o", "pipefail", "-c", cmd], check=True, cwd=out, stderr=log)
    subprocess.run(["samtools", "quickcheck", "long_reads.tmp.bam"], check=True, cwd=out)
    (out / "long_reads.tmp.bam").rename(out / "long_reads.bam")
    subprocess.run(["samtools", "index", "long_reads.bam"], check=True, cwd=out)
    stats = subprocess.run(["samtools", "view", "long_reads.bam"], capture_output=True, text=True, check=True, cwd=out).stdout
    n = bases = primary = 0
    for line in stats.splitlines():
        t = line.split("\t", 10)
        n += 1
        if not int(t[1]) & 0x900:
            primary += 1
            bases += len(t[9])
    info = {
        "run": key, "address": address, "stamp": stamp, "started": started, "seconds": round(time.time() - t0, 1),
        "samtools": subprocess.run(["samtools", "--version"], capture_output=True, text=True).stdout.splitlines()[0],
        "cluster_table_rows": n_rows, "regions": len(regions), "bases_asked_for": int((regions["end"] - regions["start"]).sum()),
        "records": n, "records_not_secondary_or_supplementary": primary, "read_bases": bases,
        "bam_bytes": (out / "long_reads.bam").stat().st_size,
    }
    done.write_text(json.dumps(info, indent=1) + "\n")
    print(json.dumps(info, indent=1))
