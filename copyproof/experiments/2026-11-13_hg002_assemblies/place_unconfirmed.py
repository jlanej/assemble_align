#!/usr/bin/env python3
"""Put every unconfirmed high-confidence copy-number issue into the truth genome's coordinates.

  place_unconfirmed.py DATA assemblies.tsv unconfirmed_high.tsv [--queue-behind PATTERN] [--threads N]

DATA is the data root: it holds the assemblies (hg002/gqc/NAME/), the truth genome
(hg002/asm/hg002v1.1.fasta.gz) and one copyproof run per assembly (work/hg002_assemblies/NAME).
unconfirmed_high.tsv is the output of unconfirmed.py; its rows are numbered u01, u02, ... in order.

This is the part that loads genomes. It writes only under DATA/work/hg002_assemblies/unconfirmed/
and skips every step whose output is there already:

  NAME.regions.bed, NAME.calls.fa   the sequence of each call, cut from its assembly (cpk extract);
                                    every call is cut 30 bases longer, so that the k-mers starting
                                    in its last 30 bases are whole
  NAME.track.npz                    the tracked k-mers of each call from that assembly's run: position,
                                    copies in the assembly, read count, copies in the truth genome
  ebv.fa                            the Epstein-Barr virus sequence that the draft hg002v1.0.1 carries
                                    (the truth genome v1.1 holds no such sequence)
  calls.fa                          all calls
  calls.truth.hits.tsv              cpk locate of the calls in the truth genome: tracked k-mers
                                    (1 in 4) that occur at most 64 times there
  calls.ebv.hits.tsv                cpk locate of the calls in ebv.fa

sort_unconfirmed.py turns these into the committed tables.

--queue-behind PATTERN waits, before each step that loads a genome, while a process whose command
line matches PATTERN (a regular expression, as for pgrep -f) is running: for a disk that serves
one heavy job at a time. Write the pattern so that it cannot match this command line, which holds
it: '/cp[k] |run_sampl[e]\\.py' and not '/cpk |run_sample\\.py'. Otherwise another job that waits
in the same way sees this one as the thing to wait for.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
from copyproof.engine import engine_path  # noqa: E402
from copyproof.track import Scan  # noqa: E402

K = 31
EBV = ("hg002v1.0.1", "chrEBV")  # the assembly and sequence the virus genome is taken from
SAMPLE, MAX_OCC = 4, 64  # the k-mers the runs track; those found more than 64 times in the truth genome anchor nothing


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", file=sys.stderr, flush=True)


def wait_for(pattern: str | None, what: str) -> None:
    """Sleep while another process matching `pattern` runs; it must be gone twice, 10 s apart.
    This script and the shell that started it carry the pattern in their own command lines, so
    command lines that name this script are not counted."""
    if not pattern:
        return
    clear = 0
    told = False
    me = Path(__file__).name
    while clear < 2:
        found = subprocess.run(["pgrep", "-fl", pattern], capture_output=True, text=True).stdout.splitlines()
        busy = [x for x in found if me not in x]
        if busy:
            clear = 0
            if not told:
                log(f"waiting before {what}: {len(busy)} process(es) match '{pattern}'")
                told = True
            time.sleep(10)
        else:
            clear += 1
            if clear < 2:
                time.sleep(10)
    if told:
        log(f"the other job is gone; starting {what}")


def read_calls(path: str) -> pd.DataFrame:
    u = pd.read_csv(path, sep="\t")
    u.insert(0, "call", [f"u{i + 1:02d}" for i in range(len(u))])
    return u


def assemblies(path: str) -> dict:
    out = {}
    for line in Path(path).read_text().splitlines():
        if line and not line.startswith("#"):
            t = line.split("\t")
            out[t[0]] = t[1]
    return out


def cpk(args: list, logfile: Path) -> None:
    with open(logfile, "a") as lf:
        subprocess.run([engine_path()] + [str(a) for a in args], check=True, stderr=lf)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("data")
    ap.add_argument("assemblies")
    ap.add_argument("unconfirmed")
    ap.add_argument("--queue-behind")
    ap.add_argument("--threads", type=int, default=8)
    a = ap.parse_args()
    D = Path(a.data)
    runs = D / "work/hg002_assemblies"
    W = runs / "unconfirmed"
    W.mkdir(parents=True, exist_ok=True)
    elog = W / "engine.log"
    calls = read_calls(a.unconfirmed)
    fasta = assemblies(a.assemblies)

    for name, g in calls.groupby("assembly", sort=False):
        run = runs / name
        contigs = pd.read_csv(run / "scan.contigs.tsv", sep="\t", dtype={"name": str}).set_index("name")["length"]
        bed = "".join(f"{r.contig}\t{int(r.start)}\t{min(int(r.end) + K - 1, int(contigs[r.contig]))}\t{r.call}\n"
                      for r in g.itertuples(index=False))
        extra = f"{EBV[1]}\t0\t{int(contigs[EBV[1]])}\t{EBV[1]}\n" if name == EBV[0] else ""
        bedf, fa, npz = W / f"{name}.regions.bed", W / f"{name}.calls.fa", W / f"{name}.track.npz"
        if not (bedf.is_file() and fa.is_file() and bedf.read_text() == bed + extra):
            wait_for(a.queue_behind, f"cpk extract on {name}")
            fa.unlink(missing_ok=True)
            bedf.write_text(bed + extra)
            log(f"{name}: cutting {len(g)} calls out of {fasta[name]}")
            tmp = W / f"{name}.calls.tmp"
            cpk(["extract", "--regions", bedf, "--asm", f"asm={D / 'hg002/gqc' / name / fasta[name]}", "--out", tmp], elog)
            tmp.rename(fa)
        if not npz.is_file():
            # reads the records of the calls' sequences only: a few hundred MB, so no waiting
            log(f"{name}: tracked k-mers of {len(g)} calls from the run's track")
            sc = Scan(str(run / "scan"))
            ti = sc.reads.index("truth")
            lam = json.loads((run / "summary.json").read_text())["depth_per_copy"]
            arrays = {"lam": np.float64(lam)}
            for contig, gc in g.groupby("contig", sort=False):
                ci = sc.contig_index(contig)
                gi = sc.group_index(ci)
                pos, m, c = sc.records(ci)
                for r in gc.itertuples(index=False):
                    i, j = np.searchsorted(pos, [int(r.start), int(r.end)])
                    arrays[f"{r.call}.pos"] = pos[i:j].astype(np.int64)
                    arrays[f"{r.call}.m"] = m[i:j].sum(axis=1).astype(np.int32)  # copies in the whole assembly
                    arrays[f"{r.call}.m_own"] = m[i:j, gi].astype(np.int32)  # copies in the sequence's own group
                    arrays[f"{r.call}.reads"] = c[i:j, 0].astype(np.int32)
                    arrays[f"{r.call}.truth"] = c[i:j, ti].astype(np.int32)
                del pos, m, c
            np.savez(str(npz).removesuffix(".npz") + ".tmp.npz", **arrays)
            Path(str(npz).removesuffix(".npz") + ".tmp.npz").rename(npz)

    # all calls in one file; the virus genome apart
    all_fa, ebv_fa = W / "calls.fa", W / "ebv.fa"
    text, ebv = [], []
    for name in calls["assembly"].unique():
        keep = None
        for line in (W / f"{name}.calls.fa").read_text().splitlines():
            if line.startswith(">"):
                keep = ebv if line[1:] == EBV[1] else text
            keep.append(line)
    new = "\n".join(text) + "\n"
    if not (all_fa.is_file() and all_fa.read_text() == new):
        for f in ("calls.truth.hits.tsv", "calls.ebv.hits.tsv"):
            (W / f).unlink(missing_ok=True)
        all_fa.write_text(new)
    ebv_fa.write_text("\n".join(ebv) + "\n")

    truth = D / "hg002/asm/hg002v1.1.fasta.gz"
    for out, ref in ((W / "calls.truth.hits.tsv", f"truth={truth}"), (W / "calls.ebv.hits.tsv", f"ebv={ebv_fa}")):
        if out.is_file():
            continue
        if ref.startswith("truth="):  # the virus genome is 173 kb: nothing to wait for
            wait_for(a.queue_behind, "cpk locate into the truth genome")
        log(f"locating the calls: {out.name}")
        tmp = Path(str(out) + ".tmp")
        cpk(["locate", "--query", all_fa, "--asm", ref, "--out", tmp, "--k", K, "--sample", SAMPLE, "--max-occ", MAX_OCC,
             "--threads", a.threads], elog)
        tmp.rename(out)
    log(f"done: {W}")


if __name__ == "__main__":
    main()
