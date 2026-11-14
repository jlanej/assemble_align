#!/usr/bin/env python3
"""Run copyproof on one assembly with existing read k-mer databases (and optionally the reads,
to localise missing k-mers). A thin wrapper over copyproof.pipeline for the experiments here.

  run_sample.py --out DIR --asm LABEL=FASTA ... [--group-map TSV] --db sample=X.cpkdb
                [--db mat=M.cpkdb --db pat=P.cpkdb] [--truth-db T.cpkdb] [--reads FASTQ ...]
"""

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from copyproof.pipeline import Config, run  # noqa: E402
from copyproof.reads import ReadSource  # noqa: E402


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--asm", action="append", required=True)
    ap.add_argument("--group-map")
    ap.add_argument("--group-regex", action="append", default=[])
    ap.add_argument("--db", action="append", required=True, help="LABEL=database; the first is the sample")
    ap.add_argument("--reads", nargs="*")
    ap.add_argument("--reference")
    ap.add_argument("--threads", type=int, default=10)
    ap.add_argument("--no-missing", action="store_true")
    ap.add_argument("--long-read-bam")
    ap.add_argument("--cards", type=int, default=0)
    a = ap.parse_args()
    dbs = [d.split("=", 1) for d in a.db]
    cfg = Config(
        out=a.out, asm=[tuple(x.split("=", 1)) for x in a.asm], group_map=a.group_map,
        group_regex=[tuple(x.split("=", 1)) for x in a.group_regex if "=" in x],
        reads=ReadSource(a.reads, a.reference) if a.reads else None, reads_db=dbs[0][1],
        parents={k: v for k, v in dbs[1:]}, threads=a.threads, missing=not a.no_missing, long_read_bam=a.long_read_bam,
        cards=a.cards,
    )
    s = run(cfg)
    print(json.dumps({k: s[k] for k in ("depth_per_copy", "assessable_fraction", "issues", "cn_segments", "missing", "timings")}, indent=1))


# Worker processes re-import this file; without the guard each of them would start the run again.
if __name__ == "__main__":
    main()
