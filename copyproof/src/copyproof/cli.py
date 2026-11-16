"""Command line of copyproof."""

from __future__ import annotations

import argparse
import json
import sys

from . import __version__
from .pipeline import Config, run
from .reads import ReadSource


def _pairs(values, what):
    out = []
    for v in values or []:
        if "=" not in v:
            sys.exit(f"copyproof: {what} expects LABEL=VALUE, got '{v}'")
        a, b = v.split("=", 1)
        out.append((a, b))
    return out


def _cards(a) -> int:
    """Draw cards after the fact: a run draws only its best-scoring issues."""
    from pathlib import Path

    import pandas as pd

    from .cards import write_cards
    from .cn import CnParams, assessable
    from .model import Model
    from .track import Scan
    from .windows import Windows

    run_dir = Path(a.run)
    issues = pd.read_csv(run_dir / "issues.tsv", sep="\t")
    pick = issues["id"].isin(a.id)
    for reg in a.region:
        seq, span = reg.rsplit(":", 1)
        start, end = (int(x.replace(",", "")) for x in span.split("-"))
        pick |= (issues["contig"] == seq) & (issues["start"] < end) & (issues["end"] > start)
    if not pick.any():
        sys.exit("copyproof cards: no issue matches --id or --region")
    sc = Scan(str(run_dir / "scan"))
    wfile = sorted(run_dir.glob("windows.w*.npz"))[0]
    w = Windows.load(str(wfile))
    model = Model.from_json((run_dir / "model.json").read_text())
    ok = assessable(w, model, CnParams())
    parents = (sc.reads[1], sc.reads[2]) if w.hap is not None and len(sc.reads) >= 3 else None
    n = write_cards(issues, w, model, ok, sc.gaps, a.out or run_dir / "cards", parents, only=issues[pick])
    print(f"wrote {n} cards to {a.out or run_dir / 'cards'}")
    return 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        prog="copyproof",
        description="Proofread a genome assembly with the sample's own short reads: copies that are missing, doubled, "
        "unsupported or on the wrong haplotype, and read sequence the assembly lacks.",
    )
    ap.add_argument("--version", action="version", version=f"copyproof {__version__}")
    sub = ap.add_subparsers(dest="cmd", required=True)

    r = sub.add_parser("run", help="run every stage for one sample")
    r.add_argument("--out", required=True, help="output directory")
    r.add_argument("--asm", action="append", required=True, metavar="LABEL=FASTA",
                   help="one FASTA of the release with its group label (hap1, hap2, other); repeat for every file, "
                        "including unplaced, EBV and mitochondrial sequence")
    r.add_argument("--group-regex", action="append", metavar="LABEL=REGEX",
                   help="assign sequences to a group by name instead of by file (first match wins)")
    r.add_argument("--group-map", help="TSV of sequence name and group label")
    r.add_argument("--reads", nargs="+", metavar="FILE", help="the sample's reads: FASTQ/FASTA (.gz) or BAM/CRAM")
    r.add_argument("--reads-db", help="an existing k-mer database of the sample's reads (skips counting; "
                                      "missing k-mers are then counted but not localised)")
    r.add_argument("--parent", action="append", metavar="NAME=FILE[,FILE...]",
                   help="reads (or a .cpkdb database) of one parent; give two for the hapmer audit")
    r.add_argument("--reference", help="FASTA that decodes CRAM inputs")
    r.add_argument("-k", type=int, default=31)
    r.add_argument("--sample", type=int, default=4, help="track 1 in SAMPLE k-mers (a power of two; default 4)")
    r.add_argument("--threads", type=int, default=8)
    r.add_argument("--window", type=int, default=1000)
    r.add_argument("--cap", type=int, default=8, help="highest k-mer multiplicity in the release that is used")
    r.add_argument("--genome-size", default="3.1g", help="haploid genome size, for sizing the counting table")
    r.add_argument("--no-missing", action="store_true", help="skip the localisation of missing k-mers (S2)")
    r.add_argument("--cards", type=int, default=200, help="evidence cards to draw, best-scoring first (0: none)")
    r.add_argument("--long-read-bam", help="long reads aligned to the release (indexed BAM, path or URL): audit every "
                                           "run of N for reads that link its two sides, and ask whether long reads "
                                           "hold the assembly k-mers that no short read holds")

    c = sub.add_parser("cards", help="draw evidence cards for chosen issues of a finished run")
    c.add_argument("--run", required=True, help="output directory of `copyproof run`")
    c.add_argument("--id", nargs="*", default=[], help="issue ids (cp000123 ...)")
    c.add_argument("--region", nargs="*", default=[], metavar="SEQ:START-END", help="every issue that overlaps a region")
    c.add_argument("--out", help="where to write (default: RUN/cards)")

    a = ap.parse_args(argv)
    if a.cmd == "cards":
        return _cards(a)
    if a.cmd == "run":
        if not a.reads and not a.reads_db:
            sys.exit("copyproof run: give --reads or --reads-db")
        parents = {}
        for name, val in _pairs(a.parent, "--parent"):
            paths = val.split(",")
            parents[name] = paths[0] if paths[0].endswith(".cpkdb") else ReadSource(paths, a.reference)
        if len(parents) not in (0, 2):
            sys.exit("copyproof run: give no parent or exactly two")
        cfg = Config(
            out=a.out, asm=_pairs(a.asm, "--asm"), reads=ReadSource(a.reads, a.reference) if a.reads else None,
            reads_db=a.reads_db, parents=parents, group_regex=_pairs(a.group_regex, "--group-regex"), group_map=a.group_map,
            k=a.k, sample=a.sample, threads=a.threads, window=a.window, cap=a.cap, genome_size=a.genome_size,
            missing=not a.no_missing, cards=a.cards, long_read_bam=a.long_read_bam,
        )
        summary = run(cfg)
        json.dump({k: summary[k] for k in ("depth_per_copy", "assessable_fraction", "issues", "issues_by_class")}, sys.stdout, indent=1)
        print()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
