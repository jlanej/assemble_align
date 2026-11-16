"""S4: join audit with long reads aligned to the assembly.

For every run of N in the release, does any long read link the two flanks? A read links a join
when one alignment runs through it, or when two alignments of the read sit on either side with
the same orientation, each anchored by enough bases. Reads that end at the join with a long clip
are counted too: many of them, and no link, is what an unsupported scaffold join looks like.

Only the neighbourhood of each gap is read, through the BAM index, so the alignment file can stay
remote (HPRC publishes read-to-assembly BAMs for every release-2 sample).
"""

from __future__ import annotations

import re
import subprocess
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

_CIG = re.compile(r"(\d+)([MIDNSHP=X])")


@dataclass
class JoinParams:
    min_gap: int = 10  # runs of N at least this long are joins
    flank: int = 10_000  # bases fetched on each side
    anchor: int = 5_000  # bases a read must align on each side to link a join
    min_mapq: int = 20
    clip_min: int = 500  # a clip at least this long, ending near the join, counts as a clip
    clip_window: int = 300
    min_depth: int = 3  # reads within `anchor` of the join, on each side, needed to trust the absence of a link


def _parse(cigar: str):
    """Reference length and leading/trailing clip of an alignment."""
    ops = _CIG.findall(cigar)
    ref = sum(int(n) for n, o in ops if o in "MDN=X")
    lead = int(ops[0][0]) if ops and ops[0][1] in "SH" else 0
    trail = int(ops[-1][0]) if len(ops) > 1 and ops[-1][1] in "SH" else 0
    return ref, lead, trail


def bam_sequences(bam: str, workdir) -> set[str]:
    out = subprocess.run(["samtools", "view", "-H", bam], capture_output=True, text=True, check=True, cwd=workdir).stdout
    return {f.split(":", 1)[1] for line in out.splitlines() if line.startswith("@SQ") for f in line.split("\t") if f.startswith("SN:")}


def fetch(bam: str, contig: str, start: int, end: int, workdir):
    """Alignments overlapping a region: (name, flag, start, end, mapq, leading clip, trailing clip)."""
    region = f"{contig}:{max(start, 0) + 1}-{end}"
    # secondary alignments are kept: where two copies are identical, half the reads are secondary here
    proc = subprocess.run(["samtools", "view", "-F", "0x4", bam, region], capture_output=True, text=True, check=True, cwd=workdir)
    rows = []
    for line in proc.stdout.splitlines():
        t = line.split("\t", 6)
        ref, lead, trail = _parse(t[5])
        pos = int(t[3]) - 1
        rows.append((t[0], int(t[1]), pos, pos + ref, int(t[4]), lead, trail))
    return rows


def audit_joins(gaps: pd.DataFrame, bam: str, workdir, p: JoinParams | None = None, lengths: dict | None = None) -> pd.DataFrame:
    """One row per join: reads that link it, depth and clips on either side, and a verdict."""
    p = p or JoinParams()
    wd = Path(workdir)
    wd.mkdir(parents=True, exist_ok=True)
    have = bam_sequences(bam, wd)
    rows = []
    for g in gaps.itertuples(index=False):
        if g.end - g.start < p.min_gap:
            continue
        row = {"contig": g.contig, "gap_start": int(g.start), "gap_end": int(g.end), "linking_reads": 0, "linking_reads_mapq": 0,
               "depth_left": 0, "depth_right": 0, "primary_left": 0, "primary_right": 0, "clips_left": 0, "clips_right": 0,
               "verdict": "not_assessable"}
        L = (lengths or {}).get(g.contig)
        at_end = g.start < p.anchor or (L is not None and L - g.end < p.anchor)
        if g.contig not in have or at_end:
            row["note"] = "sequence absent from the alignment file" if g.contig not in have else "gap at a sequence end"
            rows.append(row)
            continue
        alns = fetch(bam, g.contig, g.start - p.flank, g.end + p.flank, wd)
        left = defaultdict(list)  # read -> alignments anchored left of the gap
        right = defaultdict(list)
        linked, linked_q = set(), set()
        for name, flag, s, e, mq, lead, trail in alns:
            strand = "-" if flag & 16 else "+"
            # reads that come within `anchor` of the join on either side: were the join real,
            # about half of them would run through it
            near_l = s < g.start and e > g.start - p.anchor
            near_r = e > g.end and s < g.end + p.anchor
            row["depth_left"] += near_l
            row["depth_right"] += near_r
            row["primary_left"] += near_l and not flag & 0x100
            row["primary_right"] += near_r and not flag & 0x100
            if s <= g.start - p.anchor and e >= g.end + p.anchor:
                linked.add(name)
                if mq >= p.min_mapq:
                    linked_q.add(name)
                continue
            if e <= g.start + 50 and e - max(s, g.start - p.flank) >= p.anchor:
                left[name].append((strand, mq))
            if s >= g.end - 50 and min(e, g.end + p.flank) - s >= p.anchor:
                right[name].append((strand, mq))
            if abs(e - g.start) <= p.clip_window and trail >= p.clip_min:
                row["clips_left"] += 1
            if abs(s - g.end) <= p.clip_window and lead >= p.clip_min:
                row["clips_right"] += 1
        for name in set(left) & set(right):
            for sl, ql in left[name]:
                for sr, qr in right[name]:
                    if sl == sr:
                        linked.add(name)
                        if min(ql, qr) >= p.min_mapq:
                            linked_q.add(name)
        row["linking_reads"], row["linking_reads_mapq"] = len(linked), len(linked_q)
        if linked:
            row["verdict"] = "linked"
        elif min(row["depth_left"], row["depth_right"]) >= p.min_depth:
            row["verdict"] = "unsupported"
        rows.append(row)
    return pd.DataFrame(rows)


def issues_from_joins(joins: pd.DataFrame, pad: int = 500) -> pd.DataFrame:
    from .calls import ISSUE_COLUMNS

    rows = []
    if not len(joins):  # a release without runs of N
        return pd.DataFrame(columns=ISSUE_COLUMNS[1:])
    for r in joins[joins["verdict"] == "unsupported"].itertuples(index=False):
        rows.append(
            {
                "contig": r.contig, "start": max(int(r.gap_start) - pad, 0), "end": int(r.gap_end) + pad, "class": "UNSUPPORTED_JOIN",
                "role": "locus", "confidence": "high" if r.clips_left + r.clips_right >= 5 else "low",
                "score": float(r.clips_left + r.clips_right), "channel": "S4",
                "summary": f"no long read links the two sides of this {int(r.gap_end - r.gap_start)}-N join "
                           f"({r.depth_left} and {r.depth_right} reads on the flanks, {r.clips_left + r.clips_right} clipped at it)",
                "evidence": f"linking_reads=0;depth_left={r.depth_left};depth_right={r.depth_right};clips_left={r.clips_left};clips_right={r.clips_right}",
            }
        )
    return pd.DataFrame(rows, columns=ISSUE_COLUMNS[1:])


NEXT_TO_A_JOIN = ("MISPLACED_PIECE", "PHASE_SWITCH", "COLLAPSE", "EXPANSION", "UNSUPPORTED", "FALSE_JOIN", "MISSING_JUNCTION",
                  "MISSING_HAPLOTYPE")


def corroborate_joins(join_issues: pd.DataFrame, others: list[pd.DataFrame], pad: int = 20_000) -> pd.DataFrame:
    """A scaffold join that no long read links is ordinary: that is what a gap is (34 of the 42
    joins of one HPRC release-2 sample). It becomes a firm issue when something else is wrong next
    to it: a piece on the wrong haplotype, a copy too many or too few, a missing junction. Those
    joins get high confidence and name their neighbour; the rest stay low."""
    if not len(join_issues):
        return join_issues
    parts = [t[t["class"].isin(NEXT_TO_A_JOIN)] for t in others if t is not None and len(t)]
    near = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame(columns=join_issues.columns)
    by = {c: g for c, g in near.groupby("contig")}
    rows = join_issues.to_dict("records")
    for row in rows:
        g = by.get(row["contig"])
        hit = g[(g["start"] < row["end"] + pad) & (g["end"] > row["start"] - pad)] if g is not None else []
        if len(hit):
            classes = sorted(set(hit["class"]))
            row["confidence"] = "high"
            row["summary"] += "; next to " + ", ".join(classes)
            row["evidence"] += ";next_to=" + ",".join(classes)
        else:
            row["confidence"] = "low"
    return pd.DataFrame(rows, columns=join_issues.columns)


_COMP = str.maketrans("ACGTacgt", "TGCAtgca")


def extract_regions(regions: pd.DataFrame, asm_args: list, prefix, log=None) -> dict[str, str]:
    """Sequences of stretches of the release (columns name, contig, start, end), through the engine."""
    from .engine import engine_path

    bed, fa = Path(f"{prefix}.regions.bed"), Path(f"{prefix}.regions.fa")
    bed.write_text("".join(f"{r.contig}\t{int(r.start)}\t{int(r.end)}\t{r.name}\n" for r in regions.itertuples(index=False)))
    subprocess.run([engine_path(), "extract", "--regions", str(bed), "--out", str(fa)] + [str(x) for x in asm_args],
                   check=True, stderr=log)
    out, name = {}, None
    for line in fa.read_text().splitlines():
        if line.startswith(">"):
            name = line[1:]
        elif name is not None:
            out[name] = line
    return out


def unread_kmers(sc, regions: pd.DataFrame, seqs: dict[str, str], cap: int = 8, limit: int = 200) -> dict[str, list[str]]:
    """The tracked assembly k-mers of each region that no short read holds, as sequence (at most
    `limit` per region, evenly spread)."""
    import numpy as np

    out = {}
    for contig, g in regions.groupby("contig", sort=False):
        pos, m, c = sc.records(contig)  # one read of the sequence's track for all its regions
        unread = (c[:, 0] == 0) & (m.sum(axis=1) <= cap)
        for r in g.itertuples(index=False):
            seq = seqs.get(str(r.name), "")
            a, b = np.searchsorted(pos, [int(r.start), int(r.end)])
            sel = pos[a:b][unread[a:b]]
            if len(sel) > limit:
                sel = sel[:: -(-len(sel) // limit)]
            ks = [seq[p - int(r.start) : p - int(r.start) + sc.k] for p in sel]
            out[str(r.name)] = [k for k in ks if len(k) == sc.k and "N" not in k]
    return out


def audit_points(points: pd.DataFrame, bam: str, workdir, margin: int = 500, clip_min: int = 200, window: int = 100,
                 reach: int = 2000, threads: int = 8, kmers_of: dict | None = None, min_reads: int = 5) -> pd.DataFrame:
    """What do long reads say at each short interval (name, contig, start, end)?

    Short reads leave some correct sequence unread, and a run of assembly k-mers without short-read
    support is one of three things: such a stretch, a wrong base, or a false join. Long reads
    separate them, in two ways.

    By alignment: `spanning` counts alignments that cover the interval and `margin` bases on both
    sides; `clipped` counts alignments that stop within `window` of it with a clip of at least
    `clip_min`. `verdict` is `spanned` (at least 3 spanning), `broken` (none spanning, at least 3
    clipped) or `unclear`.

    By content, when `kmers_of` gives the unread assembly k-mers of each interval: `kmers_carried`
    counts those found in at least two long reads. `long_read_kmers` is `carried` (at least 80%:
    the sequence is real and short reads fail on it), `not_carried` (at most 20%, with at least
    `min_reads` long reads there: the release's sequence is in no read of either kind) or
    `unclear`. An alignment can run through a stretch the read does not hold; its k-mers cannot.
    """
    from concurrent.futures import ThreadPoolExecutor

    wd = Path(workdir)
    wd.mkdir(parents=True, exist_ok=True)
    have = bam_sequences(bam, wd)  # also fetches the index once, before the threads start

    def one(r):
        row = {"name": str(getattr(r, "name", "")), "contig": r.contig, "start": int(r.start), "end": int(r.end), "spanning": 0,
               "spanning_mapq": 0, "clipped": 0, "reads": 0, "verdict": "unclear", "kmers": 0, "kmers_carried": 0,
               "long_read_kmers": "unclear"}
        if r.contig not in have:
            row["verdict"] = row["long_read_kmers"] = "not_in_alignment_file"
            return row
        region = f"{r.contig}:{max(int(r.start) - reach, 0) + 1}-{int(r.end) + reach}"
        try:
            proc = subprocess.run(["samtools", "view", "-F", "0x4", bam, region], capture_output=True, text=True, check=True, cwd=wd)
        except subprocess.CalledProcessError:
            row["verdict"] = row["long_read_kmers"] = "fetch_failed"
            return row
        seqs = []
        for line in proc.stdout.splitlines():
            t = line.split("\t", 10)
            ref, lead, trail = _parse(t[5])
            s = int(t[3]) - 1
            e = s + ref
            row["reads"] += 1
            if s <= r.start - margin and e >= r.end + margin:
                row["spanning"] += 1
                row["spanning_mapq"] += int(t[4]) >= 20
            if (abs(e - r.start) <= window or abs(e - r.end) <= window) and trail >= clip_min:
                row["clipped"] += 1
            if (abs(s - r.start) <= window or abs(s - r.end) <= window) and lead >= clip_min:
                row["clipped"] += 1
            # the read's own sequence, once (secondary alignments carry none; supplementary ones a part)
            if not int(t[1]) & 0x900 and s < r.end and e > r.start and len(t) > 9 and t[9] != "*":
                seqs.append(t[9].upper())
        if row["spanning"] >= 3:
            row["verdict"] = "spanned"
        elif row["spanning"] == 0 and row["clipped"] >= 3:
            row["verdict"] = "broken"
        ks = (kmers_of or {}).get(row["name"], [])
        row["kmers"] = len(ks)
        row["long_reads_over"] = len(seqs)
        if ks:
            carried = 0
            for k in ks:
                rc = k.translate(_COMP)[::-1]
                n = 0
                for q in seqs:
                    if k in q or rc in q:
                        n += 1
                        if n >= 2:
                            break
                carried += n >= 2
            row["kmers_carried"] = carried
            if len(ks) >= 4:
                if carried >= 0.8 * len(ks):
                    row["long_read_kmers"] = "carried"
                elif carried <= 0.2 * len(ks) and len(seqs) >= min_reads:
                    row["long_read_kmers"] = "not_carried"
        return row

    with ThreadPoolExecutor(max_workers=threads) as ex:
        rows = list(ex.map(one, list(points.itertuples(index=False))))
    return pd.DataFrame(rows)


def issues_from_points(points: pd.DataFrame) -> pd.DataFrame:
    """Runs of assembly k-mers that neither short nor long reads hold.

    By content: long reads lie over the run and do not carry its k-mers. With long reads running
    through, that is a wrong base or small consensus error; with none, a false join. By
    alignment, where content cannot say (at a wrong join the reads' primary alignments go to the
    other haplotype, and secondary alignments carry no sequence): no long read runs through and
    several are clipped at the run. A run whose k-mers long reads do carry is never an issue."""
    from .calls import ISSUE_COLUMNS

    rows = []
    for r in points.itertuples(index=False):
        content = getattr(r, "long_read_kmers", "unclear")
        if content == "not_carried":
            ev = (f"unread_kmers={r.kmers};in_long_reads={r.kmers_carried};long_reads_over={r.long_reads_over};"
                  f"spanning_long_reads={r.spanning};clipped_long_reads={r.clipped}")
            if r.spanning >= 1:
                cls, conf = "BASE_ERROR", "high" if r.spanning >= 3 and r.kmers >= 6 else "low"
                what = (f"{r.kmers} assembly k-mers here are in no short read and in no long read, though {r.spanning} long "
                        "reads run through: a wrong base or small consensus error")
            else:
                cls, conf = "FALSE_JOIN", "high" if r.clipped >= 3 or r.long_reads_over >= 10 else "low"
                what = (f"{r.kmers} assembly k-mers here are in no short read and in no long read, and no long read runs "
                        f"through ({r.clipped} are clipped here): the two sides are not adjacent in the reads")
            score = float(r.kmers)
        elif r.verdict == "broken" and content != "carried":
            cls, conf, score = "FALSE_JOIN", "high" if r.clipped >= 5 else "low", float(r.clipped)
            what = (f"no short read supports the assembly across this point, no long read runs through it, and "
                    f"{r.clipped} long reads are clipped at it")
            ev = f"spanning_long_reads=0;clipped_long_reads={r.clipped};long_reads_nearby={r.reads}"
        else:
            continue
        rows.append({"contig": r.contig, "start": int(r.start), "end": int(r.end), "class": cls, "role": "locus",
                     "confidence": conf, "score": score, "channel": "S1+S4", "summary": what, "evidence": ev})
    return pd.DataFrame(rows, columns=ISSUE_COLUMNS[1:])


def adjudicate_cn(cn_issues: pd.DataFrame, points: pd.DataFrame) -> pd.DataFrame:
    """Let long reads settle copy-number issues whose own k-mers no short read holds. Where long
    reads do not carry those k-mers either, the issue stands with high confidence. Where they do,
    the sequence is real: a short issue is dropped (one allele that short reads fail on), a long
    one is kept at low confidence (sequence in the long reads' DNA and not in the short reads')."""
    if not len(cn_issues) or not len(points):
        return cn_issues
    verdict = {(r.contig, int(r.start), int(r.end)): r for r in points.itertuples(index=False)}
    keep, rows = [], cn_issues.to_dict("records")
    for row in rows:
        r = verdict.get((row["contig"], int(row["start"]), int(row["end"])))
        if r is None or r.long_read_kmers not in ("carried", "not_carried"):
            keep.append(row)
            continue
        numbers = f";unread_kmers={r.kmers};in_long_reads={r.kmers_carried};long_reads_over={r.long_reads_over}"
        if r.long_read_kmers == "carried":
            if "low_because=short" in row["evidence"]:
                continue  # a short stretch one allele of which short reads fail on
            row["confidence"] = "low"
            row["channel"] = "S1+S4"
            row["summary"] += ("; long reads carry these k-mers, so the sequence is real in their DNA: short reads fail on "
                               "it, or the two DNA samples differ")
            row["evidence"] += ";low_because=long_reads_carry" + numbers
            keep.append(row)
            continue
        row["confidence"] = "high"
        row["channel"] = "S1+S4"
        row["summary"] = row["summary"].split("; only ")[0] + f"; long reads do not carry these k-mers either ({r.kmers_carried} of {r.kmers})"
        row["evidence"] = row["evidence"].replace(";low_because=short", "") + numbers
        keep.append(row)
    return pd.DataFrame(keep, columns=cn_issues.columns)
