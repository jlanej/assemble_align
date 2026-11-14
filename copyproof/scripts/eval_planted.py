#!/usr/bin/env python3
"""Score a copyproof run on a planted genome (GUIDE section 8, rung 2).

  eval_planted.py --expected plant.expected.bed --run RUN_DIR [--baseline BASELINE_RUN_DIR] --out PREFIX

For every expected signal, was an issue of the right class reported at the right place? Issues away
from every planted edit are compared with the baseline run on the unedited genome: those also in
the baseline are not caused by planting.
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--expected", required=True)
ap.add_argument("--run", required=True)
ap.add_argument("--baseline")
ap.add_argument("--edits", help="the planting tool's edits table: maps edited coordinates back for the baseline comparison")
ap.add_argument("--out", required=True)
ap.add_argument("--min-cover", type=float, default=0.5)
ap.add_argument("--point-tol", type=int, default=200)
ap.add_argument("--pad", type=int, default=5000, help="issues within this distance of a planted edit are attributed to it")
a = ap.parse_args()

exp = pd.read_csv(a.expected, sep="\t", names=["contig", "start", "end", "signal", "tag"])
exp[["id", "kind", "size", "stratum"]] = exp["tag"].str.split("|", expand=True)
exp["size"] = exp["size"].astype(int)
iss = pd.read_csv(Path(a.run) / "issues.tsv", sep="\t")
iss["mate_contig"] = iss["evidence"].str.extract(r"mate=([^:;]+):")[0]
iss["mate_pos"] = pd.to_numeric(iss["evidence"].str.extract(r"mate=[^:;]+:(-?\d+)")[0], errors="coerce")
by = {c: g for c, g in iss.groupby("contig")}

# which issue classes answer each expected signal
ANSWER = {
    "EXPANSION": (["EXPANSION", "UNSUPPORTED"], ["locus", "unresolved", "partner"]),
    "EXPANSION_PARTNER": (["EXPANSION"], ["partner", "unresolved", "locus"]),
    "COLLAPSE": (["COLLAPSE"], ["locus", "unresolved", "partner"]),
    "MISSING_HAPLOTYPE": (["MISSING_HAPLOTYPE", "MISSING_ALLELE"], ["locus"]),
    "MISSING_JUNCTION": (["MISSING_JUNCTION"], ["locus"]),
    "FALSE_JOIN": (["FALSE_JOIN", "BASE_ERROR", "UNSUPPORTED"], ["locus"]),
    "PHASE_SWITCH": (["PHASE_SWITCH", "MISPLACED_PIECE"], ["locus"]),
}


def cover(g, s, e):
    if not len(g):
        return 0.0
    iv = sorted(zip(np.maximum(g["start"].values, s), np.minimum(g["end"].values, e)))
    tot, cur_s, cur_e = 0, None, None
    for x, y in iv:
        if y <= x:
            continue
        if cur_e is None or x > cur_e:
            tot += 0 if cur_e is None else cur_e - cur_s
            cur_s, cur_e = x, y
        else:
            cur_e = max(cur_e, y)
    tot += 0 if cur_e is None else cur_e - cur_s
    return tot / max(e - s, 1)


rows = []
for r in exp.itertuples(index=False):
    classes, roles = ANSWER[r.signal]
    g = by.get(r.contig)
    found, how, conf = False, "", ""
    if g is not None:
        g = g[g["class"].isin(classes)]
        if r.signal in ("MISSING_JUNCTION", "FALSE_JOIN"):
            mid = (r.start + r.end) // 2
            half = (r.end - r.start) // 2 + a.point_tol
            near = g[(g["start"] <= mid + half) & (g["end"] >= mid - half)]
            if r.signal == "MISSING_JUNCTION":
                # the junction is also found if it is listed from its other end
                allj = iss[(iss["class"] == "MISSING_JUNCTION") & (iss["mate_contig"] == r.contig) & ((iss["mate_pos"] - mid).abs() <= half)]
                found = len(near) > 0 or len(allj) > 0
                pick = near if len(near) else allj
            else:
                found = len(near) > 0
                pick = near
            if found:
                how, conf = pick.iloc[0]["class"], pick.iloc[0]["confidence"]
        else:
            ov = g[(g["start"] < r.end) & (g["end"] > r.start)]
            c = cover(ov, r.start, r.end)
            if r.signal == "MISSING_HAPLOTYPE":
                found = len(ov) > 0
            else:
                found = c >= a.min_cover
            if found:
                best = ov.sort_values("score", ascending=False).iloc[0]
                how, conf = f"{best['class']}/{best['role']}", best["confidence"]
            how = how or (f"cover={c:.2f}" if len(ov) else "")
    rows.append({**r._asdict(), "found": found, "as": how, "confidence": conf})
res = pd.DataFrame(rows)
res.to_csv(f"{a.out}.expected_scored.tsv", sep="\t", index=False)

tab = res.groupby(["kind", "signal", "stratum", "size"])["found"].agg(["sum", "count"]).reset_index()
tab["sensitivity"] = (tab["sum"] / tab["count"]).round(2)
tab = tab.rename(columns={"sum": "found", "count": "planted"})
tab.to_csv(f"{a.out}.sensitivity.tsv", sep="\t", index=False)
pd.set_option("display.width", 200)
pd.set_option("display.max_rows", 200)
print(tab.to_string(index=False))

# issues away from every planted edit
regions = {c: g[["start", "end"]].values for c, g in exp.groupby("contig")}


def near_edit(c, s, e):
    r = regions.get(c)
    return r is not None and bool(np.any((r[:, 0] - a.pad < e) & (r[:, 1] + a.pad > s)))


iss["planted"] = [near_edit(c, s, e) for c, s, e in zip(iss["contig"], iss["start"], iss["end"])]
# An event whose flanks occur at several places is listed at each of them as a twin, with the
# place it was assigned to; a twin belongs to a planted edit if that place does.
tw = iss["evidence"].str.extract(r"twin_of=([^:;]+):(-?\d+)")
iss["planted"] |= np.array([isinstance(c, str) and near_edit(c, int(x), int(x) + 1) for c, x in zip(tw[0], tw[1].fillna(0))], dtype=bool)
# The other copies of a copy-number issue (paralogues of an edit planted in a duplication).
oc = iss["evidence"].str.extract(r"other_copies=([^;]+)")[0]


def copy_at_edit(text):
    if not isinstance(text, str):
        return False
    for item in text.split(","):
        c, span = item.rsplit(":", 2)[:2]
        a, b = span.split("-")
        if near_edit(c, int(a), int(b)):
            return True
    return False


iss["planted"] |= oc.map(copy_at_edit).astype(bool)
away = iss[~iss["planted"]].copy()
print(f"\nissues: {len(iss)} total, {int(iss['planted'].sum())} at planted edits, {len(away)} elsewhere")
if a.baseline:
    base = pd.read_csv(Path(a.baseline) / "issues.tsv", sep="\t")
    bb = {(c, k): g[["start", "end"]].values for (c, k), g in base.groupby(["contig", "class"])}
    shifts = {}
    if a.edits:
        ed = pd.read_csv(a.edits, sep="\t")
        ed["delta"] = (ed["new_end"] - ed["new_start"]) - (ed["end"] - ed["start"])
        ed.loc[ed["kind"] == "extra", "delta"] = 0
        for c, g in ed.sort_values("new_start").groupby("contig"):
            shifts[c] = (g["new_end"].values, np.cumsum(g["delta"].values))

    def to_original(c, pos):
        """Coordinate in the unedited genome of a position outside every edit."""
        t = shifts.get(c)
        if t is None:
            return pos
        i = int(np.searchsorted(t[0], pos, side="right"))
        return pos - (int(t[1][i - 1]) if i else 0)

    def in_base(c, k, s, e):
        r = bb.get((c, k))
        s, e = to_original(c, s), to_original(c, e)
        return r is not None and bool(np.any((r[:, 0] - 2000 < e) & (r[:, 1] + 2000 > s)))

    away["in_baseline"] = [in_base(c, k, s, e) for c, k, s, e in zip(away["contig"], away["class"], away["start"], away["end"])]
    print(f"of those elsewhere, {int(away['in_baseline'].sum())} are also in the baseline run; {int((~away['in_baseline']).sum())} are new")
    new = away[~away["in_baseline"]].copy()
    if a.edits:
        # the stretch of the other haplotype that is homologous to an edit (original coordinates)
        hom = {c: g[["hom_start", "hom_end"]].values for c, g in ed.dropna(subset=["hom_contig"]).groupby("hom_contig")}

        def at_homologue(c, s, e):
            r = hom.get(c)
            s, e = to_original(c, s), to_original(c, e)
            return r is not None and bool(np.any((r[:, 0] - 4 * a.pad < e) & (r[:, 1] + 4 * a.pad > s)))

        new["at_homologue"] = [at_homologue(c, s, e) for c, s, e in zip(new["contig"], new["start"], new["end"])]
        print(f"of the new ones, {int(new['at_homologue'].sum())} lie within {4 * a.pad} bp of the stretch homologous to an edit")
        new = new[~new["at_homologue"]]
    print(new.groupby(["class", "role", "confidence"]).size().to_string())
    new.to_csv(f"{a.out}.new_elsewhere.tsv", sep="\t", index=False)
print(away.groupby(["class", "role", "confidence"]).size().to_string())
away.to_csv(f"{a.out}.elsewhere.tsv", sep="\t", index=False)
