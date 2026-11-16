#!/usr/bin/env python3
"""Sort the unconfirmed high-confidence copy-number issues by where they lie in the truth genome.

  sort_unconfirmed.py DATA assemblies.tsv unconfirmed_high.tsv OUT_PREFIX

Reads what place_unconfirmed.py left under DATA/work/hg002_assemblies/unconfirmed/, the Q100
annotation (DATA/hg002/annot), the benchmark's own alignments of each assembly to the truth genome
(DATA/hg002/gqc/NAME/) and the truth genome's own copyproof run (DATA/work/hg002_v1.1/run35).
Writes OUT_PREFIX_placements.tsv (one row per placement of a call), OUT_PREFIX_sorted.tsv (one row
per call), OUT_PREFIX_loci.tsv (one row per locus), OUT_PREFIX_by_sort.tsv, OUT_PREFIX_by_assembly.tsv
(beside the counts of scoreboard.tsv, which lies next to assemblies.tsv) and OUT_PREFIX_summary.tsv
(every count the record quotes, with its denominator).

Nothing is rescored: every call stays unconfirmed. The tables say where each call lies, what is
known about that place, and which evidence speaks for which explanation.

Placement. The tracked k-mers of a call (1 in 4) that occur at most 64 times in the truth genome
are located there; hits on one diagonal are one placement (copyproof.copies.copy_loci), kept when
it holds 20 k-mers and 5% of the located ones. The placement with the most k-mers is the `best`
one; placements with at least half as many are `principal` (the other haplotype's copy, paralogous
copies). The virus genome that the truth genome leaves out (chrEBV of the draft hg002v1.0.1) is
one more sequence to place on.

Locus. Two calls belong to one locus when the best placement of each overlaps a principal
placement of the other over half of the shorter of the two intervals.

Sort (the first rule that applies):
  virus                  the best placement is chrEBV: DNA of the cell line that the truth genome
                         leaves out
  mitochondrion          the best placement is chrM: the truth genome holds it once, a cell
                         hundreds of times or more
  mito_kmers             a call for more copies on nuclear sequence whose surplus sits on the
                         k-mers it shares with chrM (a fifth of its k-mers or more): those are
                         read at ten times the depth of its other k-mers or more, and the others at
                         no more than 1.5 times the assembly's count
  truth_differs          the truth genome holds more copies (collapse) or fewer (expansion,
                         unsupported) than the assembly of at least half of the call's k-mers,
                         counted over all of them or over those the scoring rule uses (k-mers the
                         assembly holds at most 8 times). The rule also asks that such k-mers be
                         the majority in 1-kb windows covering half of the call, which is where
                         these calls fell short
  truth_differs_in_part  the same for a fifth to a half of the k-mers
  reads_only_recurrent   the assembly holds the truth genome's copy numbers and only the reads
                         depart; and the truth genome's own run reports an issue of the same class
                         at a principal placement, or such calls come from two independent
                         assemblies or more at this locus
  reads_only_alone       the same, in one assembly (or one assembly and its other versions)
Known issues and v1.2 changes are reported beside the sort, not as a step of it.
"""

from __future__ import annotations

import gzip
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
from copyproof.copies import copy_loci  # noqa: E402
from copyproof.intervals import IntervalSet  # noqa: E402

K = 31
CAP = 8  # k-mers the assembly holds more often are not scored (as in eval_truth.py)
MIN_KMERS, MIN_FRAC = 20, 0.05
PRINCIPAL = 0.5
EBV = "chrEBV"
SORTS = ["virus", "mitochondrion", "mito_kmers", "truth_differs", "truth_differs_in_part", "reads_only_recurrent", "reads_only_alone"]
# Assemblies that are versions of one another are one witness, not several. The three drafts are
# stages of the truth genome. hprc_curated carries hprc_year1's sequence under other names: calls
# u31 and u38 are the same 7,030 bases, and every 31-mer of u41 is in u29.
LINEAGE = {"hg002v0.7": "Q100 drafts", "hg002v0.9": "Q100 drafts", "hg002v1.0.1": "Q100 drafts",
           "hprc_year1": "HPRC year 1", "hprc_curated": "HPRC year 1"}


def read_calls(path: str) -> pd.DataFrame:
    u = pd.read_csv(path, sep="\t")
    u.insert(0, "call", [f"u{i + 1:02d}" for i in range(len(u))])
    return u


def typed_bed(path: Path) -> dict:
    """contig -> (starts, ends, types) of a BED file with a type in column 4."""
    acc = defaultdict(list)
    for line in path.read_text().splitlines():
        t = line.split("\t")
        if len(t) >= 4 and not line.startswith(("#", "track", "browser")):
            acc[t[0]].append((int(t[1]), int(t[2]), t[3]))
    return {c: (np.array([x[0] for x in v]), np.array([x[1] for x in v]), np.array([x[2] for x in v], dtype=object)) for c, v in acc.items()}


def typed_overlap(tab: dict, c: str, s: int, e: int) -> dict:
    """Bases of [s, e) under each type, largest first."""
    if c not in tab:
        return {}
    a, b, t = tab[c]
    sel = (a < e) & (b > s)
    out = defaultdict(int)
    for x, y, z in zip(a[sel], b[sel], t[sel]):
        out[z] += int(min(y, e) - max(x, s))
    return dict(sorted(out.items(), key=lambda kv: -kv[1]))


def overlap(a: tuple, b: tuple) -> int:
    return max(0, min(a[2], b[2]) - max(a[1], b[1])) if a[0] == b[0] else 0


def half_of_shorter(a: tuple, b: tuple) -> bool:
    return overlap(a, b) >= 0.5 * min(a[2] - a[1], b[2] - b[1])


def items(series) -> list:
    """The comma-separated entries of a column of strings, empty cells skipped."""
    return [x for s in series if isinstance(s, str) for x in s.split(",") if x]


class Annotation:
    """What is known about a stretch of the truth genome."""

    def __init__(self, annot: Path, truth_run: Path):
        b = annot / "v4_issue_breakdown"
        self.known = IntervalSet.from_bed(str(annot / "hg002v1.1_issues.v4.bed"))
        self.known_typed = typed_bed(annot / "hg002v1.1_issues.v4.bed")
        self.submitted = typed_bed(annot / "hg002v1.1_submitted_issues.bed")
        self.sd = IntervalSet.from_bed(str(b / "hg002v1.1.SDs.013025.merged.bed"))
        self.cen = IntervalSet.from_bed(str(b / "hg002v1.1.centromeres.bed"))
        self.sat = typed_bed(b / "hg002v1.1.hsats.bed")
        rows = []
        with gzip.open(annot / "changes/hg002v1.1_to_hg002v1.2.vcf.gz", "rt") as f:
            for line in f:
                if not line.startswith("#"):
                    t = line.split("\t", 5)
                    rows.append((t[0], int(t[1]) - 1, int(t[1]) - 1 + len(t[3]), len(t[4]) - len(t[3])))
        v = pd.DataFrame(rows, columns=["contig", "start", "end", "net"])
        self.v12 = {c: g[["start", "end", "net"]].values for c, g in v.groupby("contig")}
        self.n_v12, self.n_known, self.known_bp = len(v), len(self.known), self.known.total()
        self.mask = IntervalSet.from_bed(str(truth_run / "cannot_assess.bed"))
        seg = pd.read_csv(truth_run / "cn_segments.tsv", sep="\t")
        self.seg = {c: g for c, g in seg.groupby("contig")}
        iss = pd.read_csv(truth_run / "issues.tsv", sep="\t")
        iss = iss[(iss["channel"] == "S1") & iss["class"].isin(["COLLAPSE", "EXPANSION", "UNSUPPORTED"])]
        self.iss = {c: g for c, g in iss.groupby("contig")}

    def of(self, c: str, s: int, e: int) -> dict:
        L = max(e - s, 1)
        kt = typed_overlap(self.known_typed, c, s, e)
        sat = typed_overlap(self.sat, c, s, e)
        v = self.v12.get(c)
        v = v[(v[:, 0] < e) & (v[:, 1] > s)] if v is not None else np.zeros((0, 3), dtype=int)
        g = self.seg.get(c)
        segs = [f"{r.direction}:{r.kind}:{min(r.end, e) - max(r.start, s)}"
                for r in g[(g["start"] < e) & (g["end"] > s)].itertuples(index=False)] if g is not None else []
        i = self.iss.get(c)
        issues = [f"{r._4}:{r.confidence}:{r.contig}:{r.start}-{r.end}"
                  for r in i[(i["start"] < e) & (i["end"] > s)].itertuples(index=False)] if i is not None else []
        return {
            "known_issue_bp": self.known.overlap(c, s, e), "known_issue_share": round(self.known.overlap(c, s, e) / L, 3),
            "known_issue_types": ",".join(f"{k}:{n}" for k, n in kt.items()),
            "submitted_issue": "; ".join(typed_overlap(self.submitted, c, s, e)),
            "segdup_share": round(self.sd.overlap(c, s, e) / L, 3), "centromere_share": round(self.cen.overlap(c, s, e) / L, 3),
            "satellite_share": round(sum(sat.values()) / L, 3),
            "satellites": ",".join(f"{k.replace(',', '/')}:{n}" for k, n in sat.items()),  # names such as cenSat(CER,HSAT5) hold commas
            "v12_changes": ",".join(f"{c}:{int(a) + 1}:{int(n):+d}" for a, _, n in v),
            "truth_run_masked": round(self.mask.overlap(c, s, e) / L, 3),
            "truth_run_segments": ",".join(segs), "truth_run_issues": ",".join(issues),
        }


def gqc_tables(folder: Path, name: str):
    """The benchmark's own alignments of one assembly: blocks (assembly interval -> truth intervals)
    and the structural differences it lists."""
    blocks = []
    for hap in ("mat", "pat"):
        for line in (folder / f"test{hap}covered.{name}.merged.bed").read_text().splitlines():
            t = line.split("\t")
            blocks.append((t[0], int(t[1]), int(t[2]), t[3]))
    sv = []
    for line in (folder / f"{name}.svs.bed").read_text().splitlines():
        t = line.split("\t")
        if t[5].isdigit() and t[6].isdigit():
            sv.append((t[4], int(t[5]), int(t[6]), t[3], t[9]))
    return blocks, sv


def gqc_of(blocks: list, sv: list, contig: str, s: int, e: int) -> tuple[list, str, int]:
    """Where the benchmark's alignments put [s, e) of an assembly sequence (truth intervals, each
    projected along its block when the block maps to one stretch), which of its structural
    differences lie inside, and how many bases of the call the alignments cover."""
    out, covered = [], []
    for c, a, b, names in blocks:
        if c != contig or a >= e or b <= s:
            continue
        covered.append((max(a, s), min(b, e)))
        parts = names.split("|")
        if len(parts) == 1:
            p = parts[0].rsplit(".", 3)
            tc, ts, te, strand = p[0], int(p[1]), int(p[2]), p[3]
            x0, x1 = max(a, s) - a, min(b, e) - a
            out.append((tc, ts + x0, ts + x1) if strand == "F" else (tc, te - x1, te - x0))
        else:
            for part in parts:
                p = part.rsplit(".", 3)
                out.append((p[0], int(p[1]), int(p[2])))
    cov = int(IntervalSet.from_rows(("x", a, b) for a, b in covered).total()) if covered else 0
    svs = ",".join(f"{t}:{size}" for c, a, b, t, size in sv if c == contig and a < e and b > s)
    return out, svs, cov


def main() -> None:
    data, unconfirmed, out = Path(sys.argv[1]), sys.argv[3], sys.argv[4]
    runs = data / "work/hg002_assemblies"
    W = runs / "unconfirmed"
    pd.set_option("display.width", 250)
    pd.set_option("display.max_rows", 400)
    pd.set_option("display.max_colwidth", 90)
    ann = Annotation(data / "hg002/annot", data / "work/hg002_v1.1/run35")

    calls = read_calls(unconfirmed)
    calls["direction"] = np.where(calls["class"] == "COLLAPSE", "gain", "loss")
    calls["lineage"] = [LINEAGE.get(a, a) for a in calls["assembly"]]

    # ---- placements in the truth genome and on the virus genome -----------------------------
    dt = {"query": str, "contig": str, "strand": str}
    hits = pd.concat([pd.read_csv(W / "calls.truth.hits.tsv", sep="\t", dtype=dt), pd.read_csv(W / "calls.ebv.hits.tsv", sep="\t", dtype=dt)])
    by = {q: g for q, g in hits.groupby("query", sort=False)}

    gap = {}  # share of each call's bases that are N: a scaffold gap has no k-mers to score
    name = None
    for line in (W / "calls.fa").read_text().splitlines():
        if line.startswith(">"):
            name = line[1:]
        else:
            gap[name] = round(line.count("N") / max(len(line), 1), 3)

    place, per_call = [], []
    for name, g in calls.groupby("assembly", sort=False):
        z = np.load(W / f"{name}.track.npz")
        lam = float(z["lam"])
        blocks, sv = gqc_tables(data / "hg002/gqc" / name, name)
        scored = pd.read_csv(runs / f"{name}.eval.cn_issues_scored.tsv", sep="\t").set_index(["contig", "start", "end"])
        for r in g.itertuples(index=False):
            e = dict(x.split("=", 1) for x in str(scored.loc[(r.contig, r.start, r.end), "evidence"]).split(";") if "=" in x)
            pos, m, reads, truth = z[f"{r.call}.pos"], z[f"{r.call}.m"], z[f"{r.call}.reads"], z[f"{r.call}.truth"]
            n = max(len(m), 1)
            el = m <= CAP
            same = truth > m if r.direction == "gain" else truth < m
            opp = truth < m if r.direction == "gain" else truth > m
            per_copy = reads / (lam * np.maximum(m, 1))
            h = by.get(r.call)
            located = int(h["qpos"].nunique()) if h is not None else 0
            on = {}  # k-mers of the call that chrM, or the virus genome, also holds
            for tag, seq in (("chrM", "chrM"), ("virus", EBV)):
                q = h.loc[h["contig"] == seq, "qpos"].unique() + int(r.start) if h is not None else np.zeros(0, dtype=np.int64)
                on[tag] = np.isin(pos, q)
            mito = on["chrM"]
            loci = copy_loci(h, "", 0, 0, K, min_kmers=MIN_KMERS, min_frac=MIN_FRAC) if h is not None else []
            best = loci[0]["kmers"] if loci else 0
            gq, gsv, gcov = gqc_of(blocks, sv, r.contig, int(r.start), int(r.end))
            agrees = False
            for i, d in enumerate(loci):
                iv = (d["contig"], d["start"], d["end"])
                principal = d["kmers"] >= PRINCIPAL * best
                agrees |= principal and any(overlap(iv, x) > 0 for x in gq)
                row = {"call": r.call, "assembly": r.assembly, "contig": r.contig, "start": r.start, "end": r.end, "class": r._6,
                       "rank": i + 1, "principal": principal, "truth_contig": d["contig"], "truth_start": d["start"],
                       "truth_end": d["end"], "strand": d["strand"], "kmers": d["kmers"], "share_of_call": round(d["kmers"] / n, 3),
                       "call_from": d["q_start"], "call_to": d["q_end"]}
                row.update(ann.of(*iv) if d["contig"] != EBV else {})
                place.append(row)
            per_call.append({
                "call": r.call, "delta": int(e["delta"]), "windows_assessed": int(e.get("windows", 0)), "gap_share": gap[r.call],
                "copies_by_multiplicity": e.get("copies_by_multiplicity", ""),
                "kmers": int(len(m)), "kmers_scored": int(el.sum()), "asm_copies_median": float(np.median(m)),
                "truth_copies_median": float(np.median(truth)), "reads_copies_median": round(float(np.median(reads / lam)), 1),
                "reads_per_asm_copy_median": round(float(np.median(per_copy)), 2),
                "truth_same_direction": round(float(same.mean()), 3),
                "truth_same_direction_scored": round(float(same[el].mean()), 3) if el.any() else np.nan,
                "truth_opposite": round(float(opp.mean()), 3), "truth_absent": round(float((truth == 0).mean()), 3),
                "located": round(located / n, 3), "on_virus": round(float(on["virus"].mean()), 3), "on_chrM": round(float(mito.mean()), 3),
                "reads_per_copy_on_chrM": round(float(np.median(per_copy[mito])), 1) if mito.sum() >= 20 else np.nan,
                "reads_per_copy_off_chrM": round(float(np.median(per_copy[~mito])), 2) if (~mito).sum() >= 20 else np.nan,
                "placements": len(loci), "principal_placements": sum(d["kmers"] >= PRINCIPAL * best for d in loci),
                "best_placement": f"{loci[0]['contig']}:{loci[0]['start']}-{loci[0]['end']}" if loci else "",
                "best_share_of_call": round(best / n, 3),
                "benchmark_alignment": ";".join(f"{c}:{a}-{b}" for c, a, b in gq[:4]) + (f";+{len(gq) - 4}" if len(gq) > 4 else ""),
                "benchmark_aligned_share": round(gcov / max(int(r.end) - int(r.start), 1), 3),
                "benchmark_alignment_agrees": agrees if gq else np.nan, "benchmark_svs": gsv,
            })
    place = pd.DataFrame(place)
    calls = calls.merge(pd.DataFrame(per_call), on="call")

    # ---- what is known at the principal placements of each call ------------------------------
    pp = place[place["principal"]]
    agg = []
    for r in calls.itertuples(index=False):
        g = pp[pp["call"] == r.call]
        b = g.iloc[0]
        want = "COLLAPSE" if r.direction == "gain" else ("EXPANSION", "UNSUPPORTED")
        v12 = sorted(set(items(g["v12_changes"])))
        agg.append({
            "call": r.call, "known_issue_share": b.get("known_issue_share", np.nan),
            "known_issue_share_max": g["known_issue_share"].max(), "known_issue_types": ",".join(sorted(set(items(g["known_issue_types"])))),
            "submitted_issue": "; ".join(sorted({x for x in g["submitted_issue"] if isinstance(x, str) and x})),
            "v12_changes": ",".join(v12), "segdup_share": b.get("segdup_share", np.nan), "centromere_share": b.get("centromere_share", np.nan),
            "satellite_share": b.get("satellite_share", np.nan), "satellites": b.get("satellites", ""),
            "truth_run_masked": b.get("truth_run_masked", np.nan),
            "truth_run_issue": ",".join(sorted({x for x in items(g["truth_run_issues"]) if x.startswith(want)})),
            "truth_run_segment_same_direction": ",".join(sorted({x.split(":")[1] for x in items(g["truth_run_segments"]) if x.startswith(r.direction)})),
        })
    calls = calls.merge(pd.DataFrame(agg), on="call")

    # ---- loci ----------------------------------------------------------------------------------
    best = {r.call: (r.truth_contig, r.truth_start, r.truth_end) for r in place[place["rank"] == 1].itertuples(index=False)}
    princ = defaultdict(list)
    for r in pp.itertuples(index=False):
        princ[r.call].append((r.truth_contig, r.truth_start, r.truth_end))
    parent = {c: c for c in calls["call"]}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    ids = list(best)
    for i, a in enumerate(ids):
        for b in ids[i + 1:]:
            if any(half_of_shorter(best[a], p) for p in princ[b]) and any(half_of_shorter(best[b], p) for p in princ[a]):
                parent[find(a)] = find(b)
    order = {}
    for c in calls["call"]:
        order.setdefault(find(c), f"L{len(order) + 1:02d}")
    calls["locus"] = [order[find(c)] for c in calls["call"]]
    calls["locus_assemblies"] = calls["locus"].map(calls.groupby("locus")["assembly"].nunique())
    calls["locus_lineages"] = calls["locus"].map(calls.groupby("locus")["lineage"].nunique())

    # ---- the sort ------------------------------------------------------------------------------
    calls["truth_differs_share"] = calls[["truth_same_direction", "truth_same_direction_scored"]].max(axis=1)

    def first_sort(r) -> str:
        if r.best_placement.startswith(EBV + ":"):
            return "virus"
        if r.best_placement.startswith("chrM:"):
            return "mitochondrion"
        if (r.direction == "gain" and r.on_chrM >= 0.2 and pd.notna(r.reads_per_copy_on_chrM) and pd.notna(r.reads_per_copy_off_chrM)
                and r.reads_per_copy_on_chrM >= 10 * r.reads_per_copy_off_chrM and r.reads_per_copy_off_chrM <= 1.5):
            return "mito_kmers"
        if r.truth_differs_share >= 0.5:
            return "truth_differs"
        if r.truth_differs_share >= 0.2:
            return "truth_differs_in_part"
        return "reads_only"

    calls["sort"] = [first_sort(r) for r in calls.itertuples(index=False)]
    # a locus recurs when its reads-only calls come from two independent assemblies or more
    ro = calls[calls["sort"] == "reads_only"]
    calls["locus_reads_only_lineages"] = calls["locus"].map(ro.groupby("locus")["lineage"].nunique()).fillna(0).astype(int)
    recurs = (calls["truth_run_issue"] != "") | (calls["locus_reads_only_lineages"] >= 2)
    calls.loc[calls["sort"] == "reads_only", "sort"] = np.where(recurs[calls["sort"] == "reads_only"], "reads_only_recurrent", "reads_only_alone")

    # ---- tables --------------------------------------------------------------------------------
    loci = []
    for locus, g in calls.groupby("locus", sort=False):
        bp = place[(place["call"].isin(g["call"])) & (place["rank"] == 1)]
        hull = bp.groupby("truth_contig").agg(s=("truth_start", "min"), e=("truth_end", "max"), n=("call", "size")).sort_values("n", ascending=False)
        loci.append({
            "locus": locus, "calls": len(g), "assemblies": g["assembly"].nunique(), "lineages": g["lineage"].nunique(),
            "in_assemblies": ",".join(sorted(g["assembly"].unique())), "classes": ",".join(f"{k}:{v}" for k, v in g["class"].value_counts().items()),
            "kb_min": int(g["kb"].min()), "kb_max": int(g["kb"].max()),
            "best_placements": "; ".join(f"{c}:{int(r.s)}-{int(r.e)} ({int(r.n)})" for c, r in hull.iterrows()),
            "sort": ",".join(f"{k}:{v}" for k, v in g["sort"].value_counts().items()),
            "truth_run_issue": ",".join(sorted(set(items(g["truth_run_issue"])))),
            "truth_run_segment_same_direction": ",".join(sorted(set(items(g["truth_run_segment_same_direction"])))),
            "known_issue_share_max": g["known_issue_share_max"].max(), "v12_changes": ",".join(sorted(set(items(g["v12_changes"])))),
            "satellites": ",".join(sorted({x.split(":")[0] for x in items(g["satellites"])})),
            "truth_differs_share_median": round(float(g["truth_differs_share"].median()), 3), "calls_list": ",".join(g["call"]),
        })
    loci = pd.DataFrame(loci).sort_values(["lineages", "calls"], ascending=False)

    first = ["call", "assembly", "contig", "start", "end", "kb", "class", "truth_overlap", "sort", "locus", "locus_assemblies", "locus_lineages",
             "locus_reads_only_lineages", "best_placement", "best_share_of_call", "placements", "principal_placements"]
    calls = calls[first + [c for c in calls.columns if c not in first and c not in ("summary", "lineage")]]
    place.to_csv(f"{out}_placements.tsv", sep="\t", index=False)
    calls.to_csv(f"{out}_sorted.tsv", sep="\t", index=False)
    loci.to_csv(f"{out}_loci.tsv", sep="\t", index=False)

    # ---- what the record quotes: every count with its denominator, printed and written ---------
    N = len(calls)
    facts = []

    def fact(what: str, value, of="") -> None:
        facts.append({"quantity": what, "value": value, "of": of})

    is_sort = {s: calls["sort"] == s for s in SORTS}
    fact("unconfirmed high-confidence copy-number issues", N)
    fact("... with a placement", int((calls["placements"] > 0).sum()), N)
    fact("placements", len(place))
    fact("... principal", int(place["principal"].sum()), len(place))
    al = calls[calls["benchmark_alignment_agrees"].notna()]
    fact("calls that the benchmark's own alignment of the assembly covers at all", len(al), N)
    fact("... and it overlaps a principal placement", int(al["benchmark_alignment_agrees"].astype(bool).sum()), len(al))
    for srt in SORTS:
        fact(f"sort: {srt}", int(is_sort[srt].sum()), N)

    g = calls[is_sort["virus"]]
    fact("virus: the truth genome holds none of the call's k-mers", int((g["truth_absent"] >= 0.99).sum()), len(g))
    fact("virus: lowest share of a call's k-mers that lie on chrEBV", g["on_virus"].min())
    fact("virus: copies in the reads (median over a call's k-mers), lowest", g["reads_copies_median"].min())
    fact("virus: copies in the reads, highest", g["reads_copies_median"].max())
    fact("virus: copies the assembly holds (median over a call's k-mers), values", ",".join(str(int(x)) for x in sorted(g["asm_copies_median"].unique())))
    g = calls[is_sort["mitochondrion"]]
    fact("mitochondrion: copies in the reads, lowest", g["reads_copies_median"].min())
    fact("mitochondrion: copies in the reads, highest", g["reads_copies_median"].max())
    g = calls[is_sort["mito_kmers"]]
    fact("mito_kmers: the truth genome's own run reports the same class at a principal placement", int((g["truth_run_issue"] != "").sum()), len(g))
    fact("mito_kmers: share of a call's k-mers that chrM also holds, lowest", g["on_chrM"].min())
    fact("mito_kmers: share of a call's k-mers that chrM also holds, highest", g["on_chrM"].max())
    fact("mito_kmers: reads per assembly copy on those k-mers, lowest", g["reads_per_copy_on_chrM"].min())
    fact("mito_kmers: reads per assembly copy on those k-mers, highest", g["reads_per_copy_on_chrM"].max())
    fact("mito_kmers: reads per assembly copy on the other k-mers, lowest", g["reads_per_copy_off_chrM"].min())
    fact("mito_kmers: reads per assembly copy on the other k-mers, highest", g["reads_per_copy_off_chrM"].max())

    full = calls[is_sort["truth_differs"]]
    td = calls[is_sort["truth_differs"] | is_sort["truth_differs_in_part"]]
    fact("truth_differs: by all the call's k-mers", int((full["truth_same_direction"] >= 0.5).sum()), len(full))
    fact("truth_differs: only by the k-mers the scoring rule uses", int((full["truth_same_direction"] < 0.5).sum()), len(full))
    fact("truth_differs: lowest share of the call that a truth discordance covers", full["truth_overlap"].min())
    fact("truth_differs: highest share of the call that a truth discordance covers (0.5 is under 0.5 before rounding)", full["truth_overlap"].max())
    high_m, gappy = full["asm_copies_median"] > CAP, full["gap_share"] >= 0.5
    fact("truth_differs: the median multiplicity of the call's k-mers is above 8", int(high_m.sum()), len(full))
    fact("truth_differs: at least half of the call is a gap (N)", int((gappy & ~high_m).sum()), len(full))
    fact("truth_differs: ... smallest share of gap among those", full.loc[gappy, "gap_share"].min())
    fact("truth_differs: ... largest share of gap among those", full.loc[gappy, "gap_share"].max())
    rest = full[~high_m & ~gappy]
    fact("truth_differs: neither", len(rest), len(full))
    fact("truth_differs: ... lowest share of the call that a truth discordance covers, among those", rest["truth_overlap"].min())
    fact("truth_differs: ... highest, among those", rest["truth_overlap"].max())
    fact("truth_differs, in whole or in part: the benchmark's alignment lists a structural difference inside the call",
         int((td["benchmark_svs"] != "").sum()), len(td))

    nuc = calls[~(is_sort["virus"] | is_sort["mitochondrion"])]
    touch = nuc[nuc["known_issue_share_max"] > 0]
    fact("calls whose best placement is nuclear sequence of the truth genome", len(nuc), N)
    fact("known-issue intervals after merging", ann.n_known)
    fact("known issues: a principal placement touches an interval", len(touch), len(nuc))
    fact("known issues: they cover half of a principal placement", int((nuc["known_issue_share_max"] >= 0.5).sum()), len(nuc))
    fact("known issues: they cover a tenth of a principal placement", int((nuc["known_issue_share_max"] >= 0.1).sum()), len(nuc))
    fact("known issues: largest share of a principal placement", nuc["known_issue_share_max"].max())
    fact("known issues: touched over less than 1% of the placement", int((touch["known_issue_share_max"] < 0.01).sum()), len(touch))
    fact("known issues: a curators' submitted issue at a principal placement", int((nuc["submitted_issue"] != "").sum()), len(nuc))
    changes = sorted(set(items(nuc["v12_changes"])))
    fact("v1.2: changes that lie inside a principal placement", len(changes), ann.n_v12)
    fact("v1.2: those changes (position, bases added or removed)", "; ".join(changes))
    fact("v1.2: calls with such a change", int((nuc["v12_changes"] != "").sum()), len(nuc))
    fact("truth genome's own run: an issue of the same class at a principal placement", int((nuc["truth_run_issue"] != "").sum()), len(nuc))

    ro = calls[is_sort["reads_only_recurrent"] | is_sort["reads_only_alone"]]
    fact("reads_only (recurrent or alone)", len(ro), N)
    fact("reads_only: expansions", int((ro["class"] == "EXPANSION").sum()), len(ro))
    fact("reads_only: largest share of a call's k-mers of which the truth genome holds another number, in the call's direction", ro["truth_differs_share"].max())
    fact("reads_only: truth genome's own run has an issue of the same class at a principal placement", int((ro["truth_run_issue"] != "").sum()), len(ro))
    fact("reads_only: ... no issue, but a segment of the same direction that it does not report",
         int(((ro["truth_run_issue"] == "") & (ro["truth_run_segment_same_direction"] != "")).sum()), len(ro))
    fact("reads_only: ... neither", int(((ro["truth_run_issue"] == "") & (ro["truth_run_segment_same_direction"] == "")).sum()), len(ro))
    fact("reads_only: best placement in the first 20 Mb of chromosome 13, 14, 15, 21 or 22",
         int(ro["best_placement"].str.match(r"chr(13|14|15|21|22)_[A-Z]+:(\d{1,7}|1\d{7})-").sum()), len(ro))
    rec = loci[loci["lineages"] >= 2]
    fact("loci", len(loci))
    fact("loci called in two or more independent assemblies", len(rec), len(loci))
    fact("calls at those loci", int(rec["calls"].sum()), N)
    summary = pd.DataFrame(facts)
    summary.to_csv(f"{out}_summary.tsv", sep="\t", index=False)
    print(summary.to_string(index=False))

    by_sort = calls.groupby("sort").agg(
        calls=("call", "size"), assemblies=("assembly", "nunique"), loci=("locus", "nunique"),
        collapse=("class", lambda x: int((x == "COLLAPSE").sum())), expansion=("class", lambda x: int((x == "EXPANSION").sum())),
        unsupported=("class", lambda x: int((x == "UNSUPPORTED").sum())), no_truth_discordance=("truth_overlap", lambda x: int((x == 0).sum())),
        kb_median=("kb", "median")).reindex(SORTS).fillna(0).astype(int)
    by_sort.to_csv(f"{out}_by_sort.tsv", sep="\t")
    print("\nby sort:")
    print(by_sort.to_string())

    # per assembly, beside the scoreboard's counts of high-confidence issues
    board = pd.read_csv(Path(sys.argv[2]).parent / "scoreboard.tsv", sep="\t").set_index("assembly")[["high", "high_true"]]
    board.columns = ["high_confidence", "confirmed"]
    per = board.join(calls.groupby("assembly")["sort"].value_counts().unstack(fill_value=0).reindex(columns=SORTS, fill_value=0)).fillna(0).astype(int)
    per.insert(2, "unconfirmed", per["high_confidence"] - per["confirmed"])
    per.loc["all"] = per.sum()
    per.loc["all but verkko_2023"] = per.loc["all"] - per.loc["verkko_2023"]
    per.to_csv(f"{out}_by_assembly.tsv", sep="\t")
    print("\nby assembly:")
    print(per.to_string())

    print(f"\nloci with two calls or more:")
    print(loci[loci["calls"] >= 2].drop(columns=["calls_list", "satellites", "v12_changes"]).to_string(index=False))
    print("\nper call:")
    cols = ["call", "assembly", "contig", "start", "kb", "class", "truth_overlap", "sort", "locus", "locus_lineages", "best_placement",
            "best_share_of_call", "truth_same_direction", "truth_same_direction_scored", "truth_opposite", "truth_absent",
            "reads_copies_median", "known_issue_share_max", "truth_run_issue", "truth_run_segment_same_direction", "benchmark_svs"]
    print(calls[cols].to_string(index=False))

if __name__ == "__main__":
    main()
