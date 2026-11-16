#!/usr/bin/env python3
"""At a locus of the truth genome: how many copies of its k-mers do the short reads, the parents' short reads and the long reads hold?

  locus_depth.py DATA REGIONS.tsv OUT.tsv [--bam URL_OR_PATH] [--tables-exist]

REGIONS.tsv lists stretches of the truth genome (locus, role, contig, start, end). The stretches of
one locus are all the copies of one sequence, so that a k-mer whose every copy lies inside them can
be counted in reads that were fetched by alignment to them. Role `control` marks ordinary sequence
near the locus, from which the depth of one copy is taken for each read set.

For every stretch the tracked k-mers come from the truth genome's run with both parents
(DATA/work/hg002_v1.1/run35_trio): copies in the maternal haplotype, the paternal haplotype and the
mitochondrial genome, and the counts in the child's, the mother's and the father's short reads.
With --bam, the long reads aligned to the stretches of a locus are fetched (primary alignments, each
read once) and the same k-mers are counted in them.

Depth of one copy, per read set: half the median count of the control k-mers that the truth genome
holds once in each haplotype (long reads: the median count of control k-mers found once in the
genome, because reads are fetched from one haplotype's stretch). OUT.tsv has one row per stretch
and k-mer class: copies in the truth genome, whether the mitochondrial genome holds the k-mer,
k-mers, and the median copies in each read set (count / depth of one copy). Long-read copies use
only k-mers whose every copy lies inside the fetched stretches (`kmers_closed`).

OUT.stretches.tsv has one row per stretch: the long reads whose primary alignment lies on it, how
many of them run through the whole stretch with 500 bases to spare on either side, their mean
mapping quality, and the depth of their alignments over the stretch as a multiple of the same
depth over the control stretches. This second table rests on alignment and says which copy the
aligner gave each read to; the k-mer counts do not depend on that.

OUT.loci.tsv is the summary the record quotes, one row per locus: at its first stretch, over the
classes with 50 k-mers or more (all copies inside the fetched stretches, none in the mitochondrial
genome), the lowest and the highest ratio of each read set's copies to the truth genome's; the
copies of the k-mers that the mitochondrial genome shares; and over the locus's stretches the
fewest long reads running through one and the median alignment depth. --tables-exist writes it
again from OUT.tsv and OUT.stretches.tsv without reading any reads.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
from copyproof.track import Scan  # noqa: E402

K = 31
_CODE = np.full(256, 4, np.uint8)
for _i, _ch in enumerate(b"ACGT"):
    _CODE[_ch] = _i
    _CODE[_ch + 32] = _i


def kmer_codes(seq: bytes, k: int = K):
    """Canonical 2-bit code of every k-mer of `seq` (as copyproof.kmers.canonical), and which are free of N."""
    c = _CODE[np.frombuffer(seq, np.uint8)]
    n = len(c) - k + 1
    if n <= 0:
        return np.zeros(0, np.uint64), np.zeros(0, bool)
    bad = c > 3
    cc = np.where(bad, 0, c).astype(np.uint64)
    fw = np.zeros(n, np.uint64)
    rv = np.zeros(n, np.uint64)
    for j in range(k):
        fw = (fw << np.uint64(2)) | cc[j : j + n]
        rv |= (np.uint64(3) - cc[j : j + n]) << np.uint64(2 * j)
    nb = np.concatenate([[0], np.cumsum(bad)])
    return np.minimum(fw, rv), (nb[k : k + n] - nb[:n]) == 0


def records_upto(sc: Scan, contig: str, end: int, step: int = 1 << 20):
    """Tracked k-mers of a sequence from its start up to position `end`: as Scan.records, but the
    track is read only as far as needed (the stretches here lie near the starts of chromosomes)."""
    row = sc.contigs.iloc[sc.contig_index(contig)]
    a, n = int(row.rec_start), int(row.n_rec)
    got, last, parts = 0, 0, []
    while got < n and last < end:
        rec = sc.track[a + got : a + min(got + step, n)]
        pos = last + np.cumsum(rec["gap"], dtype=np.int64)
        parts.append((pos, np.asarray(rec["m"]), np.asarray(rec["c"])))
        got += len(rec)
        last = int(pos[-1])
    pos, m, c = (np.concatenate([p[i] for p in parts]) for i in range(3))
    real = m.any(axis=1)
    return pos[real], m[real], c[real]


def ref_span(cigar: bytes) -> int:
    """Bases of the reference that an alignment covers."""
    n, num = 0, 0
    for ch in cigar:
        if 48 <= ch <= 57:
            num = num * 10 + ch - 48
        else:
            if ch in b"MDN=X":
                n += num
            num = 0
    return n


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", file=sys.stderr, flush=True)


def summarise(t: pd.DataFrame, st: pd.DataFrame | None, regions: pd.DataFrame, min_kmers: int = 50) -> pd.DataFrame:
    """One row per locus: each read set's copies as a ratio to the truth genome's, at the locus's first stretch."""
    rows = []
    for locus, g in regions[regions["role"] != "control"].groupby("locus", sort=False):
        first = f"{g.iloc[0].contig}:{int(g.iloc[0].start)}-{int(g.iloc[0].end)}"
        here = t[(t["locus"] == locus) & (t["region"] == first)]
        x = here[~here["in_mito"] & (here["kmers_closed"] >= min_kmers)]
        row = {"locus": locus, "stretch": first, "stretches": len(g), "truth_copies": ",".join(map(str, x["truth_copies"])),
               "kmers": int(x["kmers_closed"].sum())}
        for name in ("child", "mother", "father", "long_reads"):
            ratio = x[name] / x["truth_copies"]
            row[f"{name}_low"], row[f"{name}_high"] = round(float(ratio.min()), 2), round(float(ratio.max()), 2)
        mito = here[here["in_mito"]]
        row["kmers_also_in_chrM"] = int(mito["kmers"].sum())
        for name in ("child", "mother", "father"):
            row[f"{name}_copies_of_those"] = round(float((mito[name] * mito["kmers"]).sum() / mito["kmers"].sum())) if len(mito) else np.nan
        if st is not None:
            c = st[(st["locus"] == locus) & (st["role"] == "copy")]
            row["run_through_fewest"] = int(c["run_through"].min())
            row["alignment_depth_vs_control_median"] = round(float(c["depth_vs_control"].median()), 2)
        rows.append(row)
    return pd.DataFrame(rows)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("data")
    ap.add_argument("regions")
    ap.add_argument("out")
    ap.add_argument("--bam")
    ap.add_argument("--tables-exist", action="store_true")
    a = ap.parse_args()
    pd.set_option("display.width", 250)
    pd.set_option("display.max_rows", 500)
    stem = str(a.out).removesuffix(".tsv")
    if a.tables_exist:
        st_path = Path(stem + ".stretches.tsv")
        loci = summarise(pd.read_csv(a.out, sep="\t"), pd.read_csv(st_path, sep="\t") if st_path.is_file() else None,
                         pd.read_csv(a.regions, sep="\t", comment="#"))
        loci.to_csv(stem + ".loci.tsv", sep="\t", index=False)
        print(loci.to_string(index=False))
        return
    D = Path(a.data)
    fasta = D / "hg002/asm/hg002v1.1.fasta.gz"
    cache = D / "work/hg002_v1.1/allele"  # samtools keeps the index of a remote BAM in its working directory
    sc = Scan(str(D / "work/hg002_v1.1/run35_trio/scan"))
    sets = {"child": sc.reads.index("sample"), "mother": sc.reads.index("mother"), "father": sc.reads.index("father")}
    gm, gp, go = sc.groups.index("mat"), sc.groups.index("pat"), sc.groups.index("other")
    regions = pd.read_csv(a.regions, sep="\t", comment="#")

    # tracked k-mers of every stretch
    parts = []
    for contig, g in regions.groupby("contig", sort=False):
        pos, m, c = records_upto(sc, contig, int(g["end"].max()))
        for r in g.itertuples(index=False):
            i, j = np.searchsorted(pos, [int(r.start), int(r.end)])
            fa = subprocess.run(["samtools", "faidx", str(fasta), f"{contig}:{int(r.start) + 1}-{int(r.end) + K - 1}"],
                                capture_output=True, text=True, check=True).stdout
            seq = "".join(fa.split("\n")[1:]).encode()
            codes, valid = kmer_codes(seq)
            rel = pos[i:j] - int(r.start)
            assert valid[rel].all(), "a tracked position holds no k-mer"
            t = pd.DataFrame({"locus": r.locus, "role": r.role, "region": f"{contig}:{int(r.start)}-{int(r.end)}", "pos": pos[i:j],
                              "code": codes[rel], "m_mat": m[i:j, gm], "m_pat": m[i:j, gp], "m_mito": m[i:j, go]})
            for name, ri in sets.items():
                t[name] = c[i:j, ri].astype(np.int64)
            parts.append(t)
        del pos, m, c
    km = pd.concat(parts, ignore_index=True)
    km["m"] = km["m_mat"].astype(int) + km["m_pat"] + km["m_mito"]
    km["inside"] = km.groupby(["locus", "code"])["code"].transform("size")
    km["closed"] = km["inside"] == km["m"]
    log(f"{len(km)} tracked k-mers in {len(regions)} stretches")

    # long reads: each read once per locus, the stretches' k-mers counted in it
    km["long"] = np.nan
    stretches = []
    if a.bam:
        for locus, g in regions.groupby("locus", sort=False):
            sel = km["locus"] == locus
            uniq = np.unique(km.loc[sel, "code"].values)
            n = np.zeros(len(uniq), dtype=np.int64)
            seen, bases = set(), 0
            for r in g.itertuples(index=False):
                # a remote file fails now and then: the stretch is counted only when its fetch is whole
                for attempt in range(5):
                    got = subprocess.run(["samtools", "view", "-F", "0x904", a.bam, f"{r.contig}:{int(r.start) + 1}-{int(r.end)}"],
                                         capture_output=True, cwd=cache)
                    if got.returncode == 0:
                        break
                    log(f"fetch of {r.contig}:{r.start}-{r.end} failed ({got.stderr.decode()[-80:].strip()}); trying again")
                    time.sleep(20)
                else:
                    raise RuntimeError(f"samtools view failed on {r.contig}:{r.start}-{r.end}")
                here, through, mapq, covered = 0, 0, 0, 0
                for line in got.stdout.splitlines():
                    f = line.split(b"\t", 10)
                    p0 = int(f[3]) - 1
                    p1 = p0 + ref_span(f[5])
                    here += 1
                    mapq += int(f[4])
                    through += p0 <= int(r.start) - 500 and p1 >= int(r.end) + 500
                    covered += max(0, min(p1, int(r.end)) - max(p0, int(r.start)))
                    if f[0] in seen or len(f[9]) < K:
                        continue
                    seen.add(f[0])
                    bases += len(f[9])
                    codes, valid = kmer_codes(f[9])
                    codes = codes[valid]
                    at = np.searchsorted(uniq, codes)
                    at[at == len(uniq)] = 0
                    hit = uniq[at] == codes
                    np.add.at(n, at[hit], 1)
                del got
                log(f"{locus}: {r.contig}:{r.start}-{r.end} fetched; {len(seen)} reads, {bases / 1e6:.1f} Mb so far")
                stretches.append({"locus": locus, "role": r.role, "region": f"{r.contig}:{int(r.start)}-{int(r.end)}", "long_reads": here,
                                  "run_through": int(through), "mean_mapq": round(mapq / max(here, 1), 1),
                                  "alignment_depth": covered / max(int(r.end) - int(r.start), 1)})
            km.loc[sel, "long"] = n[np.searchsorted(uniq, km.loc[sel, "code"].values)]

    # depth of one copy, per read set, from the control stretches
    ctl = km[km["role"] == "control"]
    two = ctl[(ctl["m_mat"] == 1) & (ctl["m_pat"] == 1) & (ctl["m_mito"] == 0)]
    one = ctl[ctl["m"] == 1]
    lam = {name: float(np.median(two[name])) / 2 for name in sets}
    lam["child_from_single"] = float(np.median(one["child"]))
    lam["long"] = float(np.median(one["long"])) if a.bam else np.nan
    print(f"control k-mers: {len(two)} found once in each haplotype, {len(one)} found once in the genome")
    print("depth of one copy:", {k: round(v, 2) for k, v in lam.items()})

    rows = []
    for (locus, role, region, m, mito), g in km[km["m"] <= 8].groupby(["locus", "role", "region", "m", km["m_mito"] > 0], sort=False):
        row = {"locus": locus, "role": role, "region": region, "truth_copies": int(m), "in_mito": bool(mito), "kmers": len(g)}
        for name in sets:
            row[name] = round(float(np.median(g[name])) / lam[name], 2)
        gc = g[g["closed"]]
        row["kmers_closed"] = len(gc)
        row["long_reads"] = round(float(np.median(gc["long"])) / lam["long"], 2) if a.bam and len(gc) else np.nan
        rows.append(row)
    t = pd.DataFrame(rows).sort_values(["locus", "role", "region", "in_mito", "truth_copies"])
    t = t[t["kmers"] >= 20]
    t.to_csv(a.out, sep="\t", index=False)
    print(t.to_string(index=False))
    st = None
    if stretches:
        st = pd.DataFrame(stretches)
        ctl_depth = float(st.loc[st["role"] == "control", "alignment_depth"].mean())
        st["depth_vs_control"] = (st["alignment_depth"] / ctl_depth).round(2)
        st["alignment_depth"] = st["alignment_depth"].round(1)
        st.to_csv(stem + ".stretches.tsv", sep="\t", index=False)
        print(f"\nlong-read alignments per stretch (depth over the control stretches: {ctl_depth:.1f}):")
        print(st.to_string(index=False))
        arms = st[(st["role"] == "copy") & st["region"].str.match(r"chr(13|14|15|21|22)_")]
        print(f"stretches on chromosomes 13, 14, 15, 21 and 22: {len(arms)}; alignment depth against the control, median "
              f"{arms['depth_vs_control'].median():.2f}, quartiles {arms['depth_vs_control'].quantile(0.25):.2f} and "
              f"{arms['depth_vs_control'].quantile(0.75):.2f}; fewest reads running through one: {int(arms['run_through'].min())}")
    loci = summarise(t, st, regions)
    loci.to_csv(stem + ".loci.tsv", sep="\t", index=False)
    print("\nper locus, copies as a ratio to the truth genome's count:")
    print(loci.to_string(index=False))


if __name__ == "__main__":
    main()
