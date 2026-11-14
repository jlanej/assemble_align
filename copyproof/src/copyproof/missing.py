"""S2: read k-mers that are missing from the release, localised.

Solid read k-mers that occur nowhere in the release are mostly expected (recurrent sequencing
errors, contaminants). The rest mark sequence the sample has and the assembly lacks. Reads that
carry such k-mers are recruited, their k-mers that do exist in the release are located, and each
cluster of missing k-mers is classified by where the flanks of its reads fall:

MISSING_ALLELE    both flanks adjacent on one sequence: an allele or small variant is absent
MISSING_JUNCTION  the flanks lie far apart or on different sequences: a breakpoint is absent
MISSING_SEQUENCE  reads run from an anchored flank into sequence found nowhere in the release
UNPLACED          reads with no flank in the release at all (missing sequence, or a contaminant)
"""

from __future__ import annotations

import hashlib
import subprocess
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from .engine import engine_path
from .reads import ReadSource
from .track import Scan


@dataclass
class MissingParams:
    min_count_frac: float = 0.5  # a missing k-mer is solid at this fraction of the depth per copy
    min_count: int = 4
    max_count_frac: float = 3.5  # missing k-mers above this many copies are high-copy variants, not localised
    min_hits: int = 2  # missing k-mers a read needs to be recruited
    max_reads: int = 1_000_000
    max_repeat_frac: float = 0.5  # reads with a larger share of high-copy k-mers are not recruited
    max_occ: int = 4  # flank k-mers with more occurrences in the release are not used as anchors
    anchor_sample: int = 2  # 1 in this many flank k-mers is located (halves the work; flanks stay well anchored)
    min_anchor: int = 3  # co-linear anchor k-mers a flank needs
    min_reads: int = 3  # reads a cluster needs to be reported
    min_kmers: int = 3  # tracked missing k-mers a cluster needs
    local_max: int = 10_000  # flanks closer than this on one sequence make a local allele
    diag_tol: int = 30  # anchors within this many bases of one diagonal are co-linear
    locus_tol: int = 40  # placements of one locus are at most this far from the next
    min_strand: int = 1  # reads a locus needs on each strand ...
    min_strand_frac: float = 0.15  # ... and as a share of its reads


# (period, shortest span in bases) of the tandem repeats that make a k-mer "low complexity"
TANDEM_RULES = ((1, 7), (2, 10), (3, 12), (4, 14), (5, 15), (6, 18))


def tandem_mask(kmers: np.ndarray, k: int, period: int, min_span: int) -> np.ndarray:
    """True for 2-bit encoded k-mers that hold a tandem repeat of `period` spanning `min_span` bases."""
    x = kmers.astype(np.uint64)
    d = x ^ (x >> np.uint64(2 * period))
    eq = ~(d | (d >> np.uint64(1))) & np.uint64(0x5555555555555555)
    eq &= np.uint64((1 << (2 * (k - period))) - 1)  # only pairs of bases inside the k-mer
    r = eq.copy()
    for i in range(1, min_span - period):
        r &= eq >> np.uint64(2 * i)
    return r != 0


def low_complexity(kmers: np.ndarray, k: int) -> np.ndarray:
    """k-mers that span a short tandem repeat or homopolymer. Polymerase slippage makes reads with
    one repeat unit more or fewer at a high rate, and every k-mer that spans such a run looks
    like sequence the assembly lacks."""
    out = np.zeros(len(kmers), dtype=bool)
    for period, span in TANDEM_RULES:
        out |= tandem_mask(kmers, k, period, span)
    return out


def solid_missing(sc: Scan, lam: float, p: MissingParams):
    """Solid read k-mers absent from the release.

    Returns (k-mers, counts in every read set, threshold, candidate mask). Candidates are the
    k-mers worth localising: not low-complexity, and not so frequent that they must come from a
    high-copy repeat (variants of rDNA, satellites, a virus)."""
    kmers, counts = sc.novel()
    thr = max(p.min_count, int(round(p.min_count_frac * lam)))
    keep = counts[:, 0] >= thr
    kmers, counts = kmers[keep], counts[keep]
    cand = (counts[:, 0] <= p.max_count_frac * lam) & ~low_complexity(kmers, sc.k)
    return kmers, counts, thr, cand


def summarise(sc: Scan, lam: float, kmers: np.ndarray, counts: np.ndarray, thr: int, cand: np.ndarray,
              parents: tuple[int, int] | None, p: MissingParams) -> dict:
    """Genome-wide tally of solid missing k-mers. Each tracked k-mer stands for `sample` k-mers."""
    c = counts[:, 0].astype(float)
    lc = low_complexity(kmers, sc.k)
    high = c > p.max_count_frac * lam
    out = {
        "threshold": int(thr),
        "tracked_kmers": int(len(kmers)),
        "estimated_kmers": int(len(kmers)) * sc.sample,
        "low_complexity": int(lc.sum()) * sc.sample,
        "high_copy": int((high & ~lc).sum()) * sc.sample,
        "candidates": int(cand.sum()) * sc.sample,
        "candidates_about_one_copy": int((cand & (c < 1.5 * lam)).sum()) * sc.sample,
        "candidates_about_two_copies": int((cand & (c >= 1.5 * lam)).sum()) * sc.sample,
    }
    if parents is not None and counts.shape[1] > max(parents):
        a, b = counts[:, parents[0]] > 0, counts[:, parents[1]] > 0
        out["candidates_in_both_parents"] = int((cand & a & b).sum()) * sc.sample
        out["candidates_in_parent1_only"] = int((cand & a & ~b).sum()) * sc.sample
        out["candidates_in_parent2_only"] = int((cand & ~a & b).sum()) * sc.sample
        out["candidates_in_neither_parent"] = int((cand & ~a & ~b).sum()) * sc.sample
    return out


class _UnionFind:
    def __init__(self, n: int):
        self.p = list(range(n))

    def find(self, x: int) -> int:
        p = self.p
        while p[x] != x:
            p[x] = p[p[x]]
            x = p[x]
        return x

    def union(self, a: int, b: int) -> None:
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.p[rb] = ra


def read_recruits(path) -> tuple[list[str], list[str], list[list[tuple[int, int]]]]:
    """Recruited reads: names, sequences and, per read, the (k-mer index, read position) hits."""
    names, seqs, hits = [], [], []
    with open(path) as f:
        for header in f:
            seq = f.readline().strip()
            name, rest = header[1:].split(" ", 1)
            field = rest.strip().split("k=")[1]
            hs = [tuple(map(int, h.split(":"))) for h in field.split(",")]
            names.append(name)
            seqs.append(seq)
            hits.append(hs)
    return names, seqs, hits


def build_chains(anchors: pd.DataFrame, names: list[str], p: "MissingParams") -> dict:
    """Group every read's anchors into co-linear chains.

    An anchor is a read k-mer found in the release; anchors of one read on one sequence and
    strand whose diagonals (assembly position minus read position, or plus on the reverse strand)
    lie within `diag_tol` form a chain. Returns {read index: [chain, ...]} with chains as dicts.
    """
    a = anchors[anchors["occ"] <= p.max_occ]
    if not len(a):
        return {}
    read = pd.Categorical(a["query"], categories=names).codes.astype(np.int64)
    contig_code, contig_names = pd.factorize(a["contig"])
    strand = (a["strand"].values == "-").astype(np.int64)
    q = a["qpos"].values.astype(np.int64)
    pos = a["pos"].values.astype(np.int64)
    diag = np.where(strand == 0, pos - q, pos + q)
    order = np.lexsort((diag, strand, contig_code, read))
    r, c, s, d, qq = read[order], contig_code[order], strand[order], diag[order], q[order]
    brk = np.r_[True, (r[1:] != r[:-1]) | (c[1:] != c[:-1]) | (s[1:] != s[:-1]) | (d[1:] - d[:-1] > p.diag_tol)]
    starts = np.flatnonzero(brk)
    ends = np.r_[starts[1:], len(r)]
    n = ends - starts
    qmin = np.minimum.reduceat(qq, starts)
    qmax = np.maximum.reduceat(qq, starts)
    dmed = d[(starts + ends - 1) // 2]  # anchors are sorted by diagonal within a chain
    keep = (n >= p.min_anchor) & (r[starts] >= 0)
    out: dict = {}
    for ri, ci, si, di, ni, lo, hi in zip(r[starts][keep].tolist(), c[starts][keep].tolist(), s[starts][keep].tolist(),
                                          dmed[keep].tolist(), n[keep].tolist(), qmin[keep].tolist(), qmax[keep].tolist()):
        out.setdefault(ri, []).append({"contig": contig_names[ci], "strand": "-" if si else "+", "diag": di, "n": ni,
                                       "qmin": lo, "qmax": hi})
    return out


def _left_cut(chain: dict, q0: int, k: int) -> int:
    """Assembly coordinate of the cut where a read's left flank ends. q0 is the read position of
    the first missing k-mer, so the flank is read bases [0, q0 + k - 1). On the forward strand read
    base r lies at diag + r; on the reverse strand at diag + k - 1 - r."""
    return chain["diag"] + q0 + k - 1 if chain["strand"] == "+" else chain["diag"] - q0 + 1


def _right_cut(chain: dict, q1: int, k: int) -> int:
    """Assembly coordinate of the cut where a read's right flank starts: read bases [q1 + 1, end),
    where q1 is the read position of the last missing k-mer."""
    return chain["diag"] + q1 + 1 if chain["strand"] == "+" else chain["diag"] + k - 1 - q1


def classify_reads(names, hits, anchors: pd.DataFrame, k: int, p: MissingParams) -> pd.DataFrame:
    """Where each recruited read's flanks anchor, and what that implies.

    A read gets one row per placement. Sequence present on both haplotypes gives a read a placement
    on each, so an allele missing from both is reported at both. Columns: ri (read index), read,
    kind, contig, pos, end, size, contig2, pos2, strand (of the read on the assembly).
    """
    chains = build_chains(anchors, names, p)
    rows = []
    for ri, (name, hs) in enumerate(zip(names, hits)):
        q0 = min(h[1] for h in hs)
        q1 = max(h[1] for h in hs)
        ch = chains.get(ri)
        if not ch:
            rows.append((ri, name, "UNPLACED", None, -1, -1, 0, None, -1, ""))
            continue
        left = [c for c in ch if c["qmin"] < q0]
        right = [c for c in ch if c["qmax"] > q1]
        # a local allele: one sequence and strand carries both flanks close together
        local = {}
        for lc in left:
            for rc in right:
                if lc is rc:
                    cand = (lc["n"], lc, rc, 0)  # one chain spans the missing k-mers
                elif lc["contig"] == rc["contig"] and lc["strand"] == rc["strand"] and abs(rc["diag"] - lc["diag"]) <= p.local_max:
                    shift = rc["diag"] - lc["diag"] if lc["strand"] == "+" else lc["diag"] - rc["diag"]
                    cand = (lc["n"] + rc["n"], lc, rc, int(shift))
                else:
                    continue
                if lc["contig"] not in local or cand[0] > local[lc["contig"]][0]:
                    local[lc["contig"]] = cand
        if local:
            for contig, (_, lc, rc, shift) in local.items():
                a, b = _left_cut(lc, q0, k), _right_cut(rc, q1, k)
                lo, hi = min(a, b), max(a, b)
                rows.append((ri, name, "MISSING_ALLELE", contig, int(lo), int(max(hi, lo + 1)), int(shift), None, -1, lc["strand"]))
        elif left and right:
            lc = max(left, key=lambda c: c["n"])
            rc = max(right, key=lambda c: c["n"])
            # the two ends are reported in a fixed order, so reads of either strand agree
            e1 = (lc["contig"], _left_cut(lc, q0, k), lc["strand"])
            e2 = (rc["contig"], _right_cut(rc, q1, k), rc["strand"])
            if e2[:2] < e1[:2]:
                e1, e2 = e2, e1  # each end keeps the strand of the read where it anchors there
            rows.append((ri, name, "MISSING_JUNCTION", e1[0], int(e1[1]), int(e1[1]) + 1, 0, e2[0], int(e2[1]), e1[2]))
        else:
            side = left or right
            best = {}
            for c in side:
                if c["contig"] not in best or c["n"] > best[c["contig"]]["n"]:
                    best[c["contig"]] = c
            for contig, c in best.items():
                pa = _left_cut(c, q0, k) if left else _right_cut(c, q1, k)
                rows.append((ri, name, "MISSING_SEQUENCE", contig, int(pa), int(pa) + 1, 0, None, -1, c["strand"]))
    return pd.DataFrame(rows, columns=["ri", "read", "kind", "contig", "pos", "end", "size", "contig2", "pos2", "strand"])


def cluster_reads(names, hits, reads: pd.DataFrame, kmers: np.ndarray, counts: np.ndarray, sample: int, p: MissingParams,
                  parents: tuple[int, int] | None = None, labels: list[str] | None = None) -> pd.DataFrame:
    """Group read placements into loci and summarise each.

    Placements on one sequence less than `locus_tol` apart form a cluster. Clustering by position,
    not by shared k-mers, matters at genome scale: k-mers that are errors of a repeat recur in
    reads from many places and would chain unrelated loci together. A cluster is kept when enough
    different reads agree on it.
    """
    cols = ["cluster", "kind", "contig", "start", "end", "size", "contig2", "pos2", "n_reads", "n_reads_kind", "n_reads_locus",
            "reads_fwd", "reads_rev", "tracked_kmers", "estimated_kmers", "median_count", "frac_in_parent1", "frac_in_parent2",
            "other_loci"]
    placed = reads[reads["kind"] != "UNPLACED"]
    rows = []
    members: list = []
    skipped = Counter()
    if len(placed):
        code, cnames = pd.factorize(placed["contig"])
        pos = placed["pos"].values.astype(np.int64)
        order = np.lexsort((pos, code))
        code_s, pos_s = code[order], pos[order]
        brk = np.r_[True, (code_s[1:] != code_s[:-1]) | (pos_s[1:] - pos_s[:-1] > p.locus_tol)]
        starts = np.flatnonzero(brk)
        ends = np.r_[starts[1:], len(order)]
        ri_s = placed["ri"].values[order]
        kind_s = placed["kind"].values[order]
        end_s = placed["end"].values[order]
        size_s = placed["size"].values[order]
        c2_s = placed["contig2"].values[order]
        p2_s = placed["pos2"].values[order]
        fwd_s = placed["strand"].values[order] == "+"
        # where else each read is placed (the homologous haplotype, usually)
        n_place = np.bincount(placed["ri"].values, minlength=len(names))
        by_read = None
        cid = 0
        for a, b in zip(starts, ends):
            if b - a < p.min_reads:
                continue
            ris = np.unique(ri_s[a:b])
            if len(ris) < p.min_reads:
                continue
            # A real allele is read from both strands. Sequencing errors that recur at one place
            # (a base miscalled in a fixed sequence context) come from one strand only.
            first = np.unique(ri_s[a:b], return_index=True)[1] + a
            n_fwd = int(fwd_s[first].sum())
            n_rev = len(ris) - n_fwd
            if min(n_fwd, n_rev) < max(p.min_strand, p.min_strand_frac * len(ris)):
                skipped["one_strand"] += 1
                continue
            ks = np.unique(np.fromiter((h[0] for r in ris for h in hits[r]), dtype=np.int64))
            if len(ks) < p.min_kmers:
                continue
            kinds = Counter(kind_s[a:b].tolist())
            n = b - a
            kind = "MISSING_SEQUENCE"
            for cand in ("MISSING_JUNCTION", "MISSING_ALLELE"):
                if kinds[cand] >= max(2, n / 3):
                    kind = cand
                    break
            sel = np.flatnonzero(kind_s[a:b] == kind) + a
            row = {"cluster": cid, "kind": kind, "contig": cnames[code_s[a]], "start": int(np.median(pos_s[sel])),
                   "reads_fwd": n_fwd, "reads_rev": n_rev,
                   "end": int(np.median(end_s[sel])), "size": int(np.median(size_s[sel])) if kind == "MISSING_ALLELE" else 0,
                   "contig2": None, "pos2": -1, "n_reads": len(ris), "n_reads_kind": int(len(np.unique(ri_s[sel]))),
                   "n_reads_locus": len(ris), "tracked_kmers": len(ks), "estimated_kmers": len(ks) * sample,
                   "median_count": float(np.median(counts[ks, 0])), "other_loci": ""}
            if kind == "MISSING_JUNCTION":
                mates = Counter(zip(c2_s[sel].tolist(), (p2_s[sel] // max(p.locus_tol, 1)).tolist()))
                (c2, pbin), nm = mates.most_common(1)[0]
                agree = sel[(c2_s[sel] == c2) & (p2_s[sel] // max(p.locus_tol, 1) == pbin)]
                if nm < p.min_reads:
                    kind = row["kind"] = "MISSING_SEQUENCE"  # the reads do not agree on the other side
                else:
                    row.update(contig2=c2, pos2=int(np.median(p2_s[agree])), n_reads_kind=int(len(np.unique(ri_s[agree]))))
            if parents is not None and counts.shape[1] > max(parents):
                row["frac_in_parent1"] = round(float((counts[ks, parents[0]] > 0).mean()), 3)
                row["frac_in_parent2"] = round(float((counts[ks, parents[1]] > 0).mean()), 3)
            for i, label in enumerate(labels or []):
                if 0 < i < counts.shape[1]:
                    row[f"frac_in_{label}"] = round(float((counts[ks, i] > 0).mean()), 3)
            if (n_place[ris] > 1).any():
                if by_read is None:
                    by_read = placed.groupby("ri")["contig"].agg(list)
                other = Counter(c for r in ris for c in by_read[r] if c != row["contig"])
                row["other_loci"] = ";".join(f"{c}:{k}" for c, k in other.most_common(3))
            rows.append(row)
            members.append(ris)
            cid += 1
    df = pd.DataFrame(rows)
    for c in cols:
        if c not in df.columns:
            df[c] = np.nan
    if len(df):
        df = _absorb_junction_ends(df, p)
        df["twin_group"] = _twin_groups(df, members)
    else:
        df["twin_group"] = pd.Series(dtype=int)
    cols = cols + ["twin_group"] + [c for c in df.columns if c.startswith("frac_in_") and c not in cols]
    out = df[cols].sort_values(["contig", "start"]).reset_index(drop=True)
    out.attrs["one_strand"] = int(skipped["one_strand"])
    return out


def _absorb_junction_ends(df: pd.DataFrame, p: MissingParams) -> pd.DataFrame:
    """Reads that reach a junction from its far side with a short flank anchor on one side only.
    Their cluster sits at the junction's other end; it is that junction, not a separate locus."""
    j = df[df["kind"] == "MISSING_JUNCTION"]
    if not len(j):
        return df
    drop = np.zeros(len(df), dtype=bool)
    seq = df["kind"].values == "MISSING_SEQUENCE"
    for c2, p2 in zip(j["contig2"], j["pos2"]):
        drop |= seq & (df["contig"].values == c2) & (np.abs(df["start"].values - p2) <= p.locus_tol)
    return df[~drop]


def _twin_groups(df: pd.DataFrame, members: list) -> np.ndarray:
    """Clusters that hold mostly the same reads (the same event anchored on both haplotypes) share
    a group number."""
    uf = _UnionFind(len(members))
    seen: dict = {}
    for ci in df["cluster"].tolist():
        for r in members[ci].tolist():
            seen.setdefault(r, []).append(ci)
    pair = Counter()
    for cs in seen.values():
        for i in range(len(cs)):
            for j in range(i + 1, len(cs)):
                pair[(cs[i], cs[j])] += 1
    for (a, b), n in pair.items():
        if n >= 0.5 * min(len(members[a]), len(members[b])):
            uf.union(a, b)
    return np.array([uf.find(ci) for ci in df["cluster"].tolist()], dtype=np.int64)


def find_missing(sc: Scan, lam: float, source: ReadSource, asm_args: list[str], workdir, p: MissingParams | None = None,
                 parents: tuple[int, int] | None = None, threads: int = 8, log=None):
    """Run the whole channel. Returns (summary, clusters, reads)."""
    p = p or MissingParams()
    wd = Path(workdir)
    wd.mkdir(parents=True, exist_ok=True)
    kmers, counts, thr, cand = solid_missing(sc, lam, p)
    summary = summarise(sc, lam, kmers, counts, thr, cand, parents, p)
    kmers, counts = kmers[cand], counts[cand]
    if len(kmers) == 0:
        return summary, pd.DataFrame(), pd.DataFrame()
    kfile, rfile, afile = wd / "missing.kmers.u64", wd / "missing.reads.fa", wd / "missing.anchors.tsv"
    done = wd / "missing.done"
    blob = kmers.astype("<u8").tobytes()
    stamp = f"{hashlib.md5(blob).hexdigest()} min_hits={p.min_hits} max_reads={p.max_reads} repeat={p.max_repeat_frac} " \
            f"anchor_sample={p.anchor_sample} max_occ={p.max_occ}"
    reuse = done.is_file() and done.read_text() == stamp and rfile.is_file() and afile.is_file()
    if not reuse:
        kfile.write_bytes(blob)
        args = ["recruit", "--kmers", kfile, "--out", rfile, "--k", sc.k, "--min-hits", p.min_hits, "--max-reads", p.max_reads,
                "--threads", threads, "--set-sample", sc.sample]
        repeats = Path(f"{sc.prefix}.repeat.u64")
        if repeats.is_file():
            args += ["--repeat-hashes", repeats, "--sample", sc.sample, "--max-repeat-frac", p.max_repeat_frac]
        source.run(args, log=log)
    names, seqs, hits = read_recruits(rfile)
    summary["recruited_reads"] = len(names)
    summary["recruit_cap_reached"] = len(names) >= p.max_reads
    if not names:
        return summary, pd.DataFrame(), pd.DataFrame()
    if not reuse:
        subprocess.run([engine_path(), "locate", "--query", str(rfile), "--out", str(afile), "--k", str(sc.k),
                        "--sample", str(p.anchor_sample), "--max-occ", str(p.max_occ), "--threads", str(threads)]
                       + [str(a) for a in asm_args], check=True, stderr=log)
        done.write_text(stamp)
    anchors = pd.read_csv(afile, sep="\t", dtype={"query": str, "contig": str, "strand": str})
    reads = classify_reads(names, hits, anchors, sc.k, p)
    clusters = cluster_reads(names, hits, reads, kmers, counts, sc.sample, p, parents, sc.reads)
    unplaced = reads.loc[reads["kind"] == "UNPLACED", "ri"].values
    summary["unplaced_reads"] = int(len(unplaced))
    summary["unplaced_tracked_kmers"] = int(len({h[0] for r in unplaced for h in hits[r]}))
    summary["loci_on_one_strand_only"] = int(clusters.attrs.get("one_strand", 0))
    summary["clusters"] = int(len(clusters))
    for kind, n in Counter(clusters["kind"]).items():
        summary[f"clusters_{kind}"] = int(n)
    return summary, clusters, reads
