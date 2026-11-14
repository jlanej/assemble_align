#!/usr/bin/env python3
"""Compare cpk's counts of tracked k-mers with jellyfish's counts of all k-mers on the same reads.

  compare_jellyfish.py READS.fastq.gz WORKDIR [N_READS] [THREADS]

Takes the first N_READS reads (default 4,000,000), counts canonical 31-mers with jellyfish and with
`cpk count --adapter none` (jellyfish does not trim adapters), and compares every k-mer whose hash
makes it a tracked k-mer (1 in 4). Needs jellyfish 2 on the PATH. Uses about 6 GB of scratch.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parents[2]
CPK = HERE / "engine/target/release/cpk"
K, SAMPLE = 31, 4
M1, M2 = np.uint64(0xFF51AFD7ED558CCD), np.uint64(0xC4CEB9FE1A85EC53)
CODE = np.full(256, 255, dtype=np.uint8)
for i, ch in enumerate(b"ACGT"):
    CODE[ch] = i


def fmix64(x: np.ndarray) -> np.ndarray:
    x = x.copy()
    with np.errstate(over="ignore"):
        x ^= x >> np.uint64(33)
        x *= M1
        x ^= x >> np.uint64(33)
        x *= M2
        x ^= x >> np.uint64(33)
    return x


def encode(kmers: pd.Series) -> np.ndarray:
    """2-bit codes of canonical k-mers given as text (jellyfish -C and cpk dump both print the
    canonical form: the smaller of a k-mer and its reverse complement)."""
    b = np.frombuffer("".join(kmers).encode(), dtype=np.uint8).reshape(-1, K)
    c = CODE[b].astype(np.uint64)
    x = np.zeros(len(c), dtype=np.uint64)
    for j in range(K):
        x = (x << np.uint64(2)) | c[:, j]
    return x


def main() -> None:
    reads, work = sys.argv[1], Path(sys.argv[2])
    n_reads = int(sys.argv[3]) if len(sys.argv) > 3 else 4_000_000
    threads = sys.argv[4] if len(sys.argv) > 4 else "8"
    work.mkdir(parents=True, exist_ok=True)
    fq = work / "sub.fq"
    if not fq.is_file():
        subprocess.run(f"gzip -dc '{reads}' | head -n {4 * n_reads} > '{fq}'", shell=True, check=False)
    if not (work / "jf.txt").is_file():
        subprocess.run(["jellyfish", "count", "-m", str(K), "-C", "-s", "1G", "-t", threads, "-o", str(work / "sub.jf"), str(fq)], check=True)
        with open(work / "jf.txt", "w") as f:
            subprocess.run(["jellyfish", "dump", "-c", "-L", "2", str(work / "sub.jf")], check=True, stdout=f)
    if not (work / "sub.cpkdb").is_file():
        subprocess.run([str(CPK), "count", "--out", str(work / "sub"), "--k", str(K), "--sample", str(SAMPLE), "--threads", threads,
                        "--adapter", "none", str(fq)], check=True)
    with open(work / "cpk.txt", "w") as f:
        subprocess.run([str(CPK), "dump", "--db", str(work / "sub.cpkdb")], check=True, stdout=f)

    cpk = pd.read_csv(work / "cpk.txt", sep="\t", names=["kmer", "count"])
    cpk_count = pd.Series(cpk["count"].values, index=encode(cpk["kmer"]))
    assert cpk_count.index.is_unique
    n_jf = n_tracked = same = plus1 = saturated = lower = other = missing = 0
    seen = []
    for chunk in pd.read_csv(work / "jf.txt", sep=" ", names=["kmer", "count"], chunksize=4_000_000):
        code = encode(chunk["kmer"])
        tracked = (fmix64(code) % np.uint64(SAMPLE)) == 0
        n_jf += len(chunk)
        n_tracked += int(tracked.sum())
        code, jf = code[tracked], chunk["count"].values[tracked]
        got = cpk_count.reindex(code).values
        present = ~np.isnan(got)
        missing += int((~present).sum())
        d = got[present] - jf[present]
        sat = (got[present] == 65535) & (jf[present] > 65535)
        same += int((d == 0).sum())
        plus1 += int((d == 1).sum())
        saturated += int(sat.sum())
        lower += int(((d < 0) & ~sat).sum())
        other += int((d > 1).sum())
        seen.append(code[present])
    seen = np.concatenate(seen)
    only = cpk_count.drop(index=seen)
    rows = [
        ("jellyfish k-mers with count >= 2", n_jf),
        ("of those, tracked (hash multiple of 4)", n_tracked),
        ("present in the cpk database", n_tracked - missing),
        ("identical count", same),
        ("cpk higher by 1 (Bloom false positive at first sighting)", plus1),
        ("cpk higher by more than 1", other),
        ("saturated at 65,535 in cpk", saturated),
        ("cpk lower", lower),
        ("missing from cpk", missing),
        ("in cpk only", len(only)),
        ("... of those with count 2 (singletons that hit a Bloom false positive)", int((only == 2).sum())),
    ]
    out = pd.DataFrame(rows, columns=["what", "tracked_kmers"])
    out.to_csv(work / "compare_jellyfish.tsv", sep="\t", index=False)
    print(out.to_string(index=False))


if __name__ == "__main__":
    main()
