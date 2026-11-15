"""Shared by the scripts of this record: the strong missing-allele placements of a run, with the
reads and the k-mers behind each, recovered from the run's `missing/` folder.

A run's `missing_clusters.tsv` gives one row per placement and does not store which reads and which
k-mers made it. `recover` gets them again: it reads the recruited reads, their anchors in the release
and the k-mer list of the run's `missing/` folder, places every read with the package's own
`classify_reads`, and then groups the placements exactly as `cluster_reads` in
`src/copyproof/missing.py` does, keeping what that function drops: the reads and the k-mers of each
cluster. `compare` holds the result against every row of the table.

Terms as in the record of 2026-11-14 (`../2026-11-14_missing_alleles_by_parent`): placement, event,
strong, in a parent's reads, read count.
"""

from __future__ import annotations

import hashlib
import json
import re
import subprocess
import sys
from collections import Counter
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "src"))
from copyproof.kmers import fmix64, revcomp  # noqa: E402
from copyproof.missing import MissingParams, classify_reads, read_recruits, solid_missing  # noqa: E402
from copyproof.track import Scan  # noqa: E402

CPK = REPO / "engine" / "target" / "release" / "cpk"
WORK = "missing_alleles_long_reads"  # this record's folder under DATA/work
ANCHORS = ("kmerdb", "work", "hprc_r2", "hg002", "hgsvc3", "nygc")  # first-level folders of the data root

# key, what is audited, run directory under DATA/work, the script whose BAM variable names the
# release's own long reads aligned to the release (None: no long reads are looked at)
RUNS = [
    ("NA19240_release", "NA19240: HPRC release 2", "rung6/NA19240_short", "scripts/run_na19240_heldout.sh"),
    ("HG01943_release", "HG01943: HPRC release 2", "hprc_r2/HG01943", "scripts/run_hg01943.sh"),
    ("HG002_hprc_year2", "HG002: the assembly hprc_year2", "rung6/hg002_calibration/hprc_year2", None),
    ("HG002_truth_genome", "HG002: the truth genome itself", "hg002_v1.1/run35_trio", None),
]


def under_data(data: Path, path: str) -> Path:
    """A path written by a run on the machine it ran on, found again under the data root."""
    parts = Path(path).parts
    at = [i for i, p in enumerate(parts) if p in ANCHORS]
    return data.joinpath(*parts[min(at):]) if at else Path(path)


def bam_address(script: str) -> str:
    """The address of the release's long reads: the BAM variable of a driver script, and nothing else."""
    lines = [ln for ln in (REPO / script).read_text().splitlines() if ln.startswith("BAM=")]
    assert len(lines) == 1, script
    return lines[0].split("=", 1)[1].strip().strip("'\"")


def load_table(run: Path) -> pd.DataFrame:
    """The run's cluster table with the columns of the record of 2026-11-14."""
    c = pd.read_csv(run / "missing_clusters.tsv", sep="\t")
    assert c[["frac_in_mother", "frac_in_father"]].notna().all().all(), run
    c["strength"] = np.where((c["n_reads_locus"] >= 6) & (c["tracked_kmers"] >= 5), "strong", "weak")
    c["parent_max"] = c[["frac_in_mother", "frac_in_father"]].max(axis=1)
    c["in_a_parents_reads"] = np.where(c["parent_max"] >= 0.5, "yes", "no")
    return c


def strong_alleles(c: pd.DataFrame) -> pd.DataFrame:
    return c[(c["kind"] == "MISSING_ALLELE") & (c["strength"] == "strong")]


def read_count_group(v: float) -> str:
    return "under 9" if v < 9 else ("9 to 12.5" if v < 13 else "13 or more")


def decode_kmers(codes: np.ndarray, k: int) -> list[str]:
    """2-bit codes (A=0, C=1, G=2, T=3, first base most significant) as sequence."""
    codes = np.asarray(codes, dtype=np.uint64)
    out = np.empty((len(codes), k), dtype="S1")
    letters = np.array([b"A", b"C", b"G", b"T"])
    for i in range(k):
        out[:, k - 1 - i] = letters[((codes >> np.uint64(2 * i)) & np.uint64(3)).astype(np.int64)]
    return [row.tobytes().decode() for row in out]


_CODE = {"A": 0, "C": 1, "G": 2, "T": 3}


def tracked_starts(seq: str, k: int, sample: int) -> list[int]:
    """Starts of the k-mers of `seq` that the engine tracks (the hash of the canonical k-mer is a
    multiple of `sample`)."""
    mask, shift = (1 << (2 * k)) - 1, 2 * (k - 1)
    fw = rv = valid = 0
    out = []
    for i, ch in enumerate(seq):
        c = _CODE.get(ch)
        if c is None:
            valid = 0
            continue
        fw = ((fw << 2) | c) & mask
        rv = (rv >> 2) | ((3 - c) << shift)
        valid += 1
        if valid >= k and fmix64(min(fw, rv)) % sample == 0:
            out.append(i + 1 - k)
    return out


@dataclass
class Recovered:
    k: int
    sample: int
    lam: float
    kmers: np.ndarray  # the run's candidate k-mers (codes), in the order of missing.kmers.u64
    counts: np.ndarray  # their counts in every read set of the scan
    seqs: list  # recruited reads
    hits: list  # per read, (index into kmers, position in the read)
    reads: pd.DataFrame  # one row per read and placement, from classify_reads
    members: dict  # cluster number -> (read indices, k-mer indices, rows of `reads`, kind, start)
    kmer_list_is_the_runs: bool


def recover(run: Path) -> Recovered:
    """The reads and k-mers of every cluster of a run, collected the way cluster_reads collects them."""
    j = json.loads((run / "run.json").read_text())
    p = MissingParams(**j["parameters"]["missing"])
    lam = json.loads((run / "model.json").read_text())["lam"]
    sc = Scan(str(run / "scan"))
    kmers, counts, _, cand = solid_missing(sc, lam, p)
    kmers, counts = kmers[cand], counts[cand]
    m = run / "missing"
    # the stamp that find_missing leaves: the k-mer list the reads were recruited with
    same = (m / "missing.done").read_text().split()[0] == hashlib.md5(kmers.astype("<u8").tobytes()).hexdigest()
    names, seqs, hits = read_recruits(m / "missing.reads.fa")
    anchors = pd.read_csv(m / "missing.anchors.tsv", sep="\t", dtype={"query": str, "contig": str, "strand": str})
    reads = classify_reads(names, hits, anchors, sc.k, p)
    del anchors
    placed = reads[reads["kind"] != "UNPLACED"]
    # from here on: the loop of cluster_reads, line for line, without the summary it writes
    code, _ = pd.factorize(placed["contig"])
    pos = placed["pos"].values.astype(np.int64)
    order = np.lexsort((pos, code))
    code_s, pos_s = code[order], pos[order]
    brk = np.r_[True, (code_s[1:] != code_s[:-1]) | (pos_s[1:] - pos_s[:-1] > p.locus_tol)]
    starts = np.flatnonzero(brk)
    ends = np.r_[starts[1:], len(order)]
    ri_s = placed["ri"].values[order]
    kind_s = placed["kind"].values[order]
    fwd_s = placed["strand"].values[order] == "+"
    row_s = placed.index.values[order]
    members, cid = {}, 0
    for a, b in zip(starts, ends):
        if b - a < p.min_reads:
            continue
        ris = np.unique(ri_s[a:b])
        if len(ris) < p.min_reads:
            continue
        first = np.unique(ri_s[a:b], return_index=True)[1] + a
        n_fwd = int(fwd_s[first].sum())
        if min(n_fwd, len(ris) - n_fwd) < max(p.min_strand, p.min_strand_frac * len(ris)):
            continue
        ks = np.unique(np.fromiter((h[0] for r in ris for h in hits[r]), dtype=np.int64))
        if len(ks) < p.min_kmers:
            continue
        kinds = Counter(kind_s[a:b].tolist())
        kind = "MISSING_SEQUENCE"
        for cand_kind in ("MISSING_JUNCTION", "MISSING_ALLELE"):
            if kinds[cand_kind] >= max(2, (b - a) / 3):
                kind = cand_kind
                break
        members[cid] = (ris, ks, row_s[a:b], kind, int(np.median(pos_s[a:b][kind_s[a:b] == kind])))
        cid += 1
    return Recovered(sc.k, sc.sample, lam, kmers, counts, seqs, hits, reads, members, same)


def compare(table: pd.DataFrame, rec: Recovered) -> pd.DataFrame:
    """Every row of the cluster table beside what was recovered for its cluster number."""
    rows = []
    for r in table.itertuples(index=False):
        ris, ks, idx, kind, start = rec.members[int(r.cluster)]
        sub = rec.reads.loc[idx]
        first = sub.drop_duplicates("ri")
        fwd = int((first["strand"] == "+").sum())
        rows.append({
            "cluster": int(r.cluster), "kind": r.kind, "strength": r.strength,
            "sequence": bool((sub["contig"] == r.contig).all()),
            # a junction whose reads do not agree on the other side is written as one-sided sequence
            "kind_recovered": kind == r.kind or (kind == "MISSING_JUNCTION" and r.kind == "MISSING_SEQUENCE"),
            "reads": len(ris) == int(r.n_reads_locus),
            "strands": fwd == int(r.reads_fwd) and len(ris) - fwd == int(r.reads_rev),
            "tracked_kmers": len(ks) == int(r.tracked_kmers),
            "read_count": float(np.median(rec.counts[ks, 0])) == float(r.median_count),
            "start": start == int(r.start),
            "in_the_mother": round(float((rec.counts[ks, 1] > 0).mean()), 3) == float(r.frac_in_mother),
            "in_the_father": round(float((rec.counts[ks, 2] > 0).mean()), 3) == float(r.frac_in_father),
        })
    return pd.DataFrame(rows)


def release_sequences(data: Path, run: Path, regions: pd.DataFrame, prefix: Path) -> dict[str, str]:
    """Stretches of the release (columns name, contig, start, end), with `cpk extract --regions`.
    The FASTA is kept beside `prefix` and read again when the same regions are asked for."""
    j = json.loads((run / "run.json").read_text())
    bed, fa = Path(f"{prefix}.regions.bed"), Path(f"{prefix}.regions.fa")
    text = "".join(f"{r.contig}\t{int(r.start)}\t{int(r.end)}\t{r.name}\n" for r in regions.itertuples(index=False))
    if not (fa.is_file() and bed.is_file() and bed.read_text() == text):
        prefix.parent.mkdir(parents=True, exist_ok=True)
        fa.unlink(missing_ok=True)
        bed.write_text(text)
        args = [str(CPK), "extract", "--regions", str(bed), "--out", str(fa) + ".tmp"]
        for a in j["asm"]:
            args += ["--asm", f"{a['label']}={under_data(data, a['path'])}"]
        with open(f"{prefix}.extract.log", "w") as log:
            subprocess.run(args, check=True, stderr=log)
        Path(str(fa) + ".tmp").rename(fa)
    out, name = {}, None
    for line in fa.read_text().splitlines():
        if line.startswith(">"):
            name = line[1:]
        elif name is not None:
            out[name] = line.upper()
    return out


# ---- the reads' allele against the release ------------------------------------------------------------

FLANK = 300  # bases of the release taken on either side of a placement


def _consensus(rec: Recovered, ris, idx):
    """The placement's reads laid over one another by the k-mers they share.

    Every read is turned to the release's strand. The first read fixes a coordinate for each of
    its missing k-mers; every further read is shifted so that its k-mers fall on those coordinates.
    Returns (consensus bases by coordinate, reads used, reads that fit no shift, coordinate of
    every k-mer)."""
    sub = rec.reads.loc[idx].drop_duplicates("ri").set_index("ri")
    k = rec.k
    oriented = {}
    for ri in ris:
        seq = rec.seqs[ri]
        if sub.at[ri, "strand"] == "-":
            oriented[ri] = (revcomp(seq), {h: len(seq) - k - q for h, q in rec.hits[ri]})
        else:
            oriented[ri] = (seq, {h: q for h, q in rec.hits[ri]})
    todo = sorted(oriented, key=lambda r: -len(oriented[r][1]))
    u, shift, odd = dict(oriented[todo[0]][1]), {todo[0]: 0}, set()
    progress = True
    while progress:
        progress = False
        for ri in todo:
            if ri in shift or ri in odd:
                continue
            known = Counter(u[h] - q for h, q in oriented[ri][1].items() if h in u)
            if not known:
                continue
            (s, n), = known.most_common(1)
            if n < sum(known.values()):  # the read's k-mers do not all sit at one shift
                odd.add(ri)
                continue
            shift[ri] = s
            for h, q in oriented[ri][1].items():
                u.setdefault(h, q + s)
            progress = True
    cols: dict = {}
    for ri, s in shift.items():
        for i, ch in enumerate(oriented[ri][0]):
            cols.setdefault(i + s, []).append(ch)
    cons = {}
    for x, bases in cols.items():
        (b, n), = Counter(bases).most_common(1)
        # a column counts when two reads or more agree on it and they are three quarters of its reads
        if n >= 2 and n >= 0.75 * len(bases) and b in "ACGT":
            cons[x] = b
    return cons, len(shift), len(ris) - len(shift), u


MAX_MISMATCH = 3  # a k-mer laid on the release without gaps may differ from it at this many bases
SLACK = 10  # bases beyond the placement's interval within which a k-mer may still be laid


def _lay(window: str, kmers: list[str], s0: int, e0: int, k: int) -> list[list]:
    """Each k-mer laid on the release without gaps where it fits best.

    A k-mer of the placement holds the changed base, so it overlaps the interval [s0, e0] of the
    window. Both strands are tried at every start from s0 - k + 1 - SLACK to e0 + SLACK. Returns, per
    k-mer, the places with the fewest mismatches when those are MAX_MISMATCH or fewer: tuples of
    (window position, base of the k-mer there), one tuple per place."""
    w = np.frombuffer(window.encode(), dtype=np.uint8)
    lo, hi = max(s0 - k + 1 - SLACK, 0), min(e0 + SLACK, len(w) - k)
    view = np.lib.stride_tricks.sliding_window_view(w[lo:hi + k], k)
    out = []
    for s in kmers:
        found, best = [], MAX_MISMATCH
        for t in (s, revcomp(s)):
            a = np.frombuffer(t.encode(), dtype=np.uint8)
            differ = view != a
            h = differ.sum(axis=1)
            m = int(h.min())
            if m > best:
                continue
            if m < best:
                found, best = [], m
            for x in np.flatnonzero(h == m):
                found.append(tuple((lo + int(x) + int(i), t[int(i)]) for i in np.flatnonzero(differ[x])))
        out.append(found)
    return out


def changed_base(rec: Recovered, cluster: int, start: int, end: int, window: str, window_start: int) -> dict:
    """What the placement's reads hold where the release differs.

    `window` is the release's sequence around the placement. Each tracked k-mer of the placement is
    laid on it without gaps where it fits best. `resolved` is "one changed base" when every k-mer
    differs from the release at one and the same position by one and the same base. When the k-mers
    differ at more positions, those are counted. When a k-mer does not fit without gaps, the reads
    are laid over one another and compared with the release between two anchors (`between_anchors`),
    which gives the difference in length."""
    ris, ks, idx = rec.members[cluster][:3]
    k = rec.k
    out = {"resolved": "", "changes": 0, "release_position": -1, "release_base": "", "reads_base": "", "before": "", "after": "",
           "release_has": "", "reads_have": "", "kmers": len(ks), "kmers_laid_without_gaps": 0}
    laid = _lay(window, decode_kmers(rec.kmers[ks], k), start - 1 - window_start, end - window_start, k)
    out["kmers_laid_without_gaps"] = sum(1 for f in laid if f)
    if all(laid):
        sure = set().union(*[set(f[0]) for f in laid if len(f) == 1])
        chosen = []
        for f in laid:
            # a k-mer that fits at several places takes the one that the other k-mers point to
            fit = f if len(f) == 1 else [x for x in f if set(x) <= sure]
            if len(fit) != 1:
                out["resolved"] = "no: a k-mer fits the release at more than one place"
                return out
            chosen.append(fit[0])
        changes = sorted(set().union(*[set(c) for c in chosen]))
        places = sorted({x for x, _ in changes})
        out["changes"] = len(places)
        if len(changes) == 1:
            at, base = changes[0]
            out.update(resolved="one changed base", release_position=window_start + at, release_base=window[at], reads_base=base,
                       before=window[at - 1], after=window[at + 1], release_has=window[at], reads_have=base)
            return out
        out.update(release_position=window_start + places[0], release_has=window[places[0]:places[-1] + 1],
                   reads_have="".join(dict(changes).get(x, window[x]) for x in range(places[0], places[-1] + 1)))
        if len(changes) > len(places):
            out["resolved"] = "no: the k-mers hold two different bases at one position"
        else:
            out["resolved"] = f"no: {len(places)} changed bases within {places[-1] - places[0] + 1} bases"
        return out
    ref_mid, alt_mid, at = between_anchors(rec, ris, idx, start, window, window_start)
    if at is None:
        out["resolved"] = "no: the k-mers do not fit without gaps, and the reads give no sequence anchored in the release on both sides"
        return out
    out.update(release_position=window_start + at, release_has=ref_mid, reads_have=alt_mid)
    d = len(alt_mid) - len(ref_mid)
    if d:
        out["resolved"] = f"no: the reads have {abs(d)} base{'s' if abs(d) > 1 else ''} {'more' if d > 0 else 'fewer'} than the release"
    else:
        diff = sum(1 for x, y in zip(ref_mid, alt_mid) if x != y)
        out["resolved"] = f"no: {diff} changed bases within {len(ref_mid)} bases, more than {MAX_MISMATCH} in one k-mer"
    return out


def between_anchors(rec: Recovered, ris, idx, start: int, window: str, window_start: int):
    """The reads' consensus between two anchors in the release: the 31 bases before its first
    missing k-mer and the 31 bases after its last one (moved outwards base by base, up to 40, when
    such a stretch is not in the release). Returns (the release's bases, the reads' bases, window
    position) of the stretch that differs, with what is common at either end taken off."""
    k = rec.k
    cons, _, _, u = _consensus(rec, ris, idx)
    lo, hi = min(u.values()), max(u.values()) + k  # the span of the missing k-mers, in the reads' coordinates
    left = right = None
    for g in range(0, 40):
        if left is None and all(x in cons for x in range(lo - g - k, lo - g)):
            at = _find_all(window, "".join(cons[x] for x in range(lo - g - k, lo - g)))
            if at:
                left = (g, at)
        if right is None and all(x in cons for x in range(hi + g, hi + g + k)):
            at = _find_all(window, "".join(cons[x] for x in range(hi + g, hi + g + k)))
            if at:
                right = (g, at)
    if left is None or right is None or not all(x in cons for x in range(lo - left[0], hi + right[0])):
        return "", "", None
    alt = "".join(cons[x] for x in range(lo - left[0], hi + right[0]))
    best = None  # of the places where the anchors are found, the pair nearest to the placement
    for xl in left[1]:
        for xr in right[1]:
            if xr >= xl + k:
                d = abs(xl + k - (start - window_start)) + abs((xr - xl - k) - len(alt))
                if best is None or d < best[0]:
                    best = (d, xl, xr)
    if best is None:
        return "", "", None
    _, xl, xr = best
    ref = window[xl + k:xr]
    n = 0
    while n < min(len(ref), len(alt)) and ref[n] == alt[n]:
        n += 1
    m = 0
    while m < min(len(ref), len(alt)) - n and ref[len(ref) - 1 - m] == alt[len(alt) - 1 - m]:
        m += 1
    return ref[n:len(ref) - m], alt[n:len(alt) - m], xl + k + n


def _find_all(text: str, what: str) -> list[int]:
    out, i = [], text.find(what)
    while i >= 0:
        out.append(i)
        i = text.find(what, i + 1)
    return out


_PAIR = {"A": "T", "C": "G", "G": "C", "T": "A"}


def substitution(release_base: str, reads_base: str, before: str, after: str) -> tuple[str, str, str]:
    """A change read from the release to the reads, written from the pyrimidine: (kind, the base
    before it, the base after it), as in the usual 96 classes."""
    if release_base in "CT":
        return f"{release_base}>{reads_base}", before, after
    return f"{_PAIR[release_base]}>{_PAIR[reads_base]}", _PAIR.get(after, "N"), _PAIR.get(before, "N")


# ---- the release's own long reads at a placement ----------------------------------------------------------
# Everything from here to `event_group` follows the block "Sorting" of NOTEBOOK entry u, which was
# written before any long read was searched.

PAD = 100  # bases on either side of a placement's interval that a long read must cover
ENOUGH = 10  # reads that hold the allele or the release's sequence, below which nothing is said
MINORITY, NEARLY_ALL = 0.25, 0.75  # share under the first: a minority; above the second: nearly all
GROUPS = ["no long read holds the allele", "a minority do", "about half do", "nearly all do", "too few reads to say"]
SHARE_STEPS = ["none", "under 0.05", "0.05 to under 0.1", "0.1 to under 0.25", "0.25 to under 0.4", "0.4 to under 0.6",
               "0.6 to 0.75", "over 0.75 to under 0.9", "0.9 to under 1", "all"]

_CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")


def both_strands(kmers) -> set[str]:
    out = set(kmers)
    out.update(revcomp(s) for s in kmers)
    return out


def release_kmers_there(rec: Recovered, cluster: int, start: int, end: int, window: str, window_start: int) -> list[str]:
    """The release's sequence at a placement, as k-mers: the tracked k-mers of the release that
    overlap the interval from the base before `start` to `end` and that none of the placement's
    short reads holds."""
    k = rec.k
    ris = rec.members[cluster][0]
    in_reads = set()
    for ri in ris:
        s = rec.seqs[ri]
        in_reads.update(s[i:i + k] for i in range(len(s) - k + 1))
    in_reads = both_strands(in_reads)
    a, b = max(start - 1 - window_start - (k - 1), 0), min(end - window_start + k, len(window))
    piece = window[a:b]
    out = []
    for i in tracked_starts(piece, k, rec.sample):
        s = piece[i:i + k]
        if s not in in_reads:
            out.append(s)
    return out


def _query_at(cigar: str, pos: int, targets: list[int]) -> list:
    """Where in the read each position of the release lies, for an alignment that starts at `pos`:
    the index of the read's base aligned there, or of the next base of the read when the position
    is deleted in the read (then with False)."""
    out = [None] * len(targets)
    r, q = pos, 0
    for n, op in _CIGAR.findall(cigar):
        n = int(n)
        if op in "M=X":
            for i, t in enumerate(targets):
                if out[i] is None and r <= t < r + n:
                    out[i] = (q + t - r, True)
            r += n
            q += n
        elif op in "DN":
            for i, t in enumerate(targets):
                if out[i] is None and r <= t < r + n:
                    out[i] = (q, False)
            r += n
        elif op in "IS":
            q += n
    return out


def long_reads_at(bam: Path, contig: str, start: int, end: int, allele: list[str], release: list[str], k: int,
                  changed_at: int = -1) -> dict:
    """The long reads at one placement, counted as the block "Sorting" says.

    `allele` are the placement's tracked k-mers and `release` the release's k-mers there. With
    `changed_at` (a position of the release), the base each covering read has there by its
    alignment is counted too."""
    lo, hi = start - 1 - PAD, end + 1 + PAD  # the window [lo, hi) that a read must cover: bases start - 1 to end, and PAD beyond
    out = {"records": 0, "left_out_secondary_or_supplementary": 0, "not_covering": 0, "covering": 0, "hold_the_allele": 0,
           "hold_the_release": 0, "hold_both": 0, "hold_neither": 0, "allele_only_outside_the_window": 0,
           "hold_every_allele_kmer": 0, "hold_under_half_of_the_allele_kmers": 0, "covering_with_mapping_quality_under_20": 0,
           "allele_holders_with_mapping_quality_under_20": 0, "bases": ""}
    if lo < 0:
        return out
    alt, ref = both_strands(allele), both_strands(release)
    proc = subprocess.run(["samtools", "view", str(bam), f"{contig}:{lo + 1}-{hi}"], capture_output=True, text=True, check=True)
    bases = Counter()
    for line in proc.stdout.splitlines():
        t = line.split("\t", 10)
        out["records"] += 1
        if int(t[1]) & 0x900 or t[9] == "*":
            out["left_out_secondary_or_supplementary"] += 1
            continue
        pos = int(t[3]) - 1
        ref_len = sum(int(n) for n, op in _CIGAR.findall(t[5]) if op in "MDN=X")
        if pos > lo or pos + ref_len < hi:
            out["not_covering"] += 1
            continue
        seq = t[9].upper()
        (q0, _), (q1, _), *rest = _query_at(t[5], pos, [lo, hi - 1] + ([changed_at] if changed_at >= 0 else []))
        stretch = seq[q0:q1 + 1]
        found = {stretch[i:i + k] for i in range(len(stretch) - k + 1)}
        n_alt = sum(1 for s in allele if s in found or revcomp(s) in found)
        has_alt, has_ref = n_alt > 0, not found.isdisjoint(ref)
        low = int(t[4]) < 20
        out["covering"] += 1
        out["covering_with_mapping_quality_under_20"] += low
        out["hold_the_allele"] += has_alt
        out["hold_the_release"] += has_ref
        out["hold_both"] += has_alt and has_ref
        out["hold_neither"] += not has_alt and not has_ref
        out["hold_every_allele_kmer"] += n_alt == len(allele)
        out["hold_under_half_of_the_allele_kmers"] += 0 < n_alt < 0.5 * len(allele)
        out["allele_holders_with_mapping_quality_under_20"] += has_alt and low
        if not has_alt and any(s in seq for s in alt):
            out["allele_only_outside_the_window"] += 1
        if rest:
            q, aligned = rest[0]
            bases[(seq[q] if aligned else "-", "allele" if has_alt else ("release" if has_ref else "neither"))] += 1
    out["bases"] = ";".join(f"{b}:{what}={n}" for (b, what), n in sorted(bases.items()))
    return out


def share_of(hold_the_allele: int, hold_the_release: int, hold_both: int) -> tuple[int, float]:
    """(reads that hold the allele or the release's sequence, share of them that hold the allele)."""
    either = hold_the_allele + hold_the_release - hold_both
    return either, (hold_the_allele / either if either else float("nan"))


def group_of(hold_the_allele: int, either: int) -> str:
    if either < ENOUGH:
        return GROUPS[4]
    if hold_the_allele == 0:
        return GROUPS[0]
    share = hold_the_allele / either
    return GROUPS[1] if share < MINORITY else (GROUPS[2] if share <= NEARLY_ALL else GROUPS[3])


def event_group(groups: list[str]) -> str:
    """An event from the groups of all its placements: the highest of nearly all, about half, a
    minority; "no long read" only when every placement has enough reads and none holds the allele."""
    for g in (GROUPS[3], GROUPS[2], GROUPS[1]):
        if g in groups:
            return g
    return GROUPS[0] if all(g == GROUPS[0] for g in groups) else GROUPS[4]


def step_of(share: float) -> str:
    """The share in finer steps; the cuts of the groups (0.25 and 0.75) are among the bounds."""
    bounds = [0.05, 0.1, 0.25, 0.4, 0.6]
    if share == 0:
        return SHARE_STEPS[0]
    for i, b in enumerate(bounds):
        if share < b:
            return SHARE_STEPS[i + 1]
    if share <= 0.75:
        return SHARE_STEPS[6]
    return SHARE_STEPS[7] if share < 0.9 else (SHARE_STEPS[8] if share < 1 else SHARE_STEPS[9])
