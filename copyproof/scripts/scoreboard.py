#!/usr/bin/env python3
"""Rung 4: one table and one figure over all benchmark assemblies scored against the truth genome.

  scoreboard.py --work DIR --list assemblies.tsv --out PREFIX

DIR holds NAME.eval.truth_eval.json and NAME.eval.cn_by_size.tsv from scripts/eval_truth.py.
Writes PREFIX.tsv (one row per assembly), PREFIX.by_size.tsv (recall and precision by segment
size, summed over assemblies), PREFIX.md (both as Markdown tables, and every sum the experiment
record quotes), PREFIX.json (the figure's numbers) and PREFIX.png.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--work", required=True)
ap.add_argument("--list", required=True)
ap.add_argument("--out", required=True)
ap.add_argument("--k", type=int, default=31)
a = ap.parse_args()

names = [line.split("\t") for line in Path(a.list).read_text().splitlines() if line and not line.startswith("#")]
rows, sizes = [], []
for t in names:
    name, what = t[0], t[3] if len(t) > 3 else ""
    f = Path(a.work) / f"{name}.eval.truth_eval.json"
    if not f.is_file():
        continue
    r = json.loads(f.read_text())
    cn, mk, uk = r["copy_number"], r["missing_kmers"], r["unsupported_kmers"]
    high = cn["gain_high_issues"] + cn["loss_high_issues"]
    high_true = cn["gain_high_true"] + cn["loss_high_true"]
    low = cn["gain_low_issues"] + cn["loss_low_issues"]
    low_true = cn["gain_low_true"] + cn["loss_low_true"]
    t_ok = cn["gain_truth_assessable"] + cn["loss_truth_assessable"]
    found = cn["gain_found_assessable"] + cn["loss_found_assessable"]
    amb = cn["gain_CN_OR_DEPTH_not_reported"] + cn["loss_CN_OR_DEPTH_not_reported"]
    amb_true = cn["gain_CN_OR_DEPTH_not_reported_true"] + cn["loss_CN_OR_DEPTH_not_reported_true"]
    anom = sum(cn[f"{d}_{k}_not_reported"] for d in ("gain", "loss") for k in ("DEPTH", "IRREGULAR"))
    anom_true = sum(cn[f"{d}_{k}_not_reported_true"] for d in ("gain", "loss") for k in ("DEPTH", "IRREGULAR"))
    wrong = uk["wrong_kmer_rate"]
    rows.append({
        "assembly": name, "what": what, "Gb": round(r["assembly_bases"] / 1e9, 3), "sequences": r["sequences"],
        "assessable": round(r["assessable_fraction"], 4),
        "truth_segments": cn["truth_segments"], "truth_Mb": round(cn["truth_bp"] / 1e6, 2),
        "truth_assessable": t_ok, "truth_found": found, "recall": round(found / t_ok, 3) if t_ok else float("nan"),
        "truth_assessable_Mb": round((cn["gain_truth_bp_assessable"] + cn["loss_truth_bp_assessable"]) / 1e6, 2),
        "truth_found_Mb": round((cn["gain_found_bp_assessable"] + cn["loss_found_bp_assessable"]) / 1e6, 2),
        "high": high, "high_true": high_true, "precision_high": round(high_true / high, 3) if high else float("nan"),
        "low": low, "low_true": low_true, "precision_low": round(low_true / low, 3) if low else float("nan"),
        "ambiguous_not_reported": amb, "ambiguous_not_reported_true": amb_true,
        "depth_anomalies": anom, "depth_anomalies_true": anom_true,
        # consensus quality from the truth's k-mers, and what the reads alone would say
        "qv_truth": round(-10 * math.log10(1 - (1 - wrong) ** (1 / a.k)), 1) if 0 < wrong < 1 else float("nan"),
        "qv_reads": r.get("kmer_support", {}).get("qv", float("nan")),
        "wrong_kmers_found_by_reads": uk["recall"], "unsupported_kmers_wrong": uk["precision"],
        "truth_kmers_missing": mk["estimated_missing_truth_kmers"], "missing_seen_in_reads": mk["recall_all_solid"],
    })
    b = pd.read_csv(Path(a.work) / f"{name}.eval.cn_by_size.tsv", sep="\t")
    b["assembly"] = name
    sizes.append(b)

tab = pd.DataFrame(rows)
tab.to_csv(f"{a.out}.tsv", sep="\t", index=False)
pd.set_option("display.width", 250)
pd.set_option("display.max_columns", 40)
print(tab.drop(columns=["what"]).to_string(index=False))
by = pd.concat(sizes, ignore_index=True)
agg = by.groupby(["direction", "size"], sort=False)[["truth_assessable", "found", "found_high", "issues", "issues_true", "high", "high_true"]].sum().reset_index()
agg["recall"] = (agg["found"] / agg["truth_assessable"]).round(3)
agg["precision"] = (agg["issues_true"] / agg["issues"]).round(3)
agg["precision_high"] = (agg["high_true"] / agg["high"]).round(3)
agg.to_csv(f"{a.out}.by_size.tsv", sep="\t", index=False)
print()
print(agg.to_string(index=False))

# ---- the same as a Markdown table, for the experiment record
def pct(a, b):
    return f"{a:,} of {b:,} ({100 * a / b:.0f}%)" if b else "none"


md = ["| Assembly | What it is | Size (Gb), sequences | Discordant with the truth (Mb) | Assessable discordances found | "
      "High-confidence issues confirmed | Low-confidence issues confirmed | QV by truth k-mers / by read k-mers |",
      "|---|---|---|---|---|---|---|---|"]
for r in tab.sort_values("truth_Mb", ascending=False).itertuples(index=False):
    md.append(f"| `{r.assembly}` | {r.what} | {r.Gb:.2f}, {r.sequences:,} | {r.truth_Mb:.1f} | {pct(r.truth_found, r.truth_assessable)} | "
              f"{pct(r.high_true, r.high)} | {pct(r.low_true, r.low)} | {r.qv_truth} / {r.qv_reads} |")

SIZE = {"2-5 kb": "under 5 kb", "5-20 kb": "5 to 20 kb", "20-100 kb": "20 to 100 kb", ">100 kb": "over 100 kb"}
DIRECTION = {"gain": "The truth holds more (collapse)", "loss": "The truth holds fewer (surplus copy)"}
md += ["", f"Summed over the {len(tab)} assemblies, by size of the segment:", "",
       "| Direction | Size | Assessable discordances | Found | Issues | Confirmed | High-confidence issues | Confirmed |",
       "|---|---|---|---|---|---|---|---|"]
for d in ("gain", "loss"):
    for i, r in enumerate(agg[agg["direction"] == d].itertuples(index=False)):
        md.append(f"| {DIRECTION[d] if i == 0 else ''} | {SIZE.get(r.size, r.size)} | {r.truth_assessable:,} | "
                  f"{r.found:,} ({100 * r.found / max(r.truth_assessable, 1):.0f}%) | {r.issues:,} | "
                  f"{pct(r.issues_true, r.issues).split(' of ')[0]} ({100 * r.issues_true / max(r.issues, 1):.0f}%) | {r.high:,} | "
                  + (f"{r.high_true:,} ({100 * r.high_true / r.high:.0f}%) |" if r.high else " |"))


def span(v, digits=0):
    v = [x for x in v if x == x]
    return f"{min(v):.{digits}f} to {max(v):.{digits}f}" if v else "none"


def share_span(agg, d, sizes, num, den):
    x = agg[(agg["direction"] == d) & agg["size"].isin(sizes)]
    return span(100 * x[num] / x[den].clip(lower=1))


big = tab.sort_values("high", ascending=False).iloc[0]
rest = tab[tab["assembly"] != big.assembly]
many = tab[tab["high"] >= 20]
few = tab[tab["high"] < 20]
unconfirmed = tab["high"] - tab["high_true"]
ge5 = ["5-20 kb", "20-100 kb", ">100 kb"]
lt5 = agg[agg["size"] == "2-5 kb"].set_index("direction")
md += ["", "Sums the record quotes:", "",
       f"- High-confidence copy-number issues confirmed: {pct(int(tab.high_true.sum()), int(tab.high.sum()))}; "
       f"without `{big.assembly}`, {pct(int(rest.high_true.sum()), int(rest.high.sum()))}.",
       f"- Per assembly, among the {len(many)} with 20 such issues or more: {span(100 * many.high_true / many.high)}% confirmed.",
       f"- The {len(few)} with fewer than 20 ({', '.join('`' + n + '`' for n in few.assembly)}): "
       f"{pct(int(few.high_true.sum()), int(few.high.sum()))} confirmed.",
       f"- High-confidence issues not confirmed, per assembly: {span(unconfirmed)} "
       f"({int(unconfirmed.sum())} in all), while the confirmed ones range from {int(tab.high_true.min()):,} to {int(tab.high_true.max()):,}.",
       f"- Found, of assessable discordances of 5 kb and more: surplus copies {share_span(agg, 'loss', ge5, 'found', 'truth_assessable')}%, "
       f"collapses {share_span(agg, 'gain', ge5, 'found', 'truth_assessable')}% (the three size classes).",
       f"- Under 5 kb: found {100 * lt5.loc['gain', 'found'] / lt5.loc['gain', 'truth_assessable']:.0f}% (collapses) and "
       f"{100 * lt5.loc['loss', 'found'] / lt5.loc['loss', 'truth_assessable']:.0f}% (surplus copies); "
       f"{pct(int(lt5.loc['loss', 'issues_true']), int(lt5.loc['loss', 'issues']))} surplus-copy issues confirmed.",
       f"- Segments that fit a depth factor or nothing regular (not reported): "
       f"{pct(int(tab.depth_anomalies_true.sum()), int(tab.depth_anomalies.sum()))} confirmed; "
       f"without `{big.assembly}`, {pct(int(rest.depth_anomalies_true.sum()), int(rest.depth_anomalies.sum()))}.",
       f"- Ambiguous segments without support (not reported): "
       f"{pct(int(tab.ambiguous_not_reported_true.sum()), int(tab.ambiguous_not_reported.sum()))} confirmed.",
       f"- Sequence discordant with the truth: {span(tab.truth_Mb, 1)} Mb.",
       f"- QV by truth k-mers minus QV by read k-mers: {span(tab.qv_truth - tab.qv_reads, 1)}; highest QV by read k-mers {tab.qv_reads.max():.2f}.",
       f"- Unread assembly k-mers that are wrong by the truth: {span(100 * tab.unsupported_kmers_wrong)}%; "
       f"k-mers wrong by the truth that the reads find: {span(100 * tab.wrong_kmers_found_by_reads)}%."]
Path(f"{a.out}.md").write_text("\n".join(md) + "\n")

# ---- figure: for each assembly, the high-confidence copy-number issues split by what the truth says
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

BLUE, ORANGE, INK, MUTED, GRID = "#2a78d6", "#eb6834", "#1f2328", "#59636e", "#d9dee3"
# Counts differ a hundredfold between assemblies, so each bar is a share and the counts are its label.
order = tab.sort_values("truth_assessable", ascending=True)
fig, axes = plt.subplots(1, 2, figsize=(11.5, 0.44 * len(order) + 2.1), sharey=True, gridspec_kw={"wspace": 0.34})
y = list(range(len(order)))
panels = (
    ("High-confidence copy-number issues:\nshare confirmed by the truth genome", "high", "high_true", "confirmed", "not confirmed"),
    ("Truth discordances in assessable sequence:\nshare found (any confidence)", "truth_assessable", "truth_found", "found", "missed"),
)
for ax, (title, total, good, lab_good, lab_rest) in zip(axes, panels):
    g, n = order[good].values.astype(float), order[total].values.astype(float)
    share = 100 * g / np.maximum(n, 1)
    ax.barh(y, share, color=BLUE, height=0.6, label=lab_good, edgecolor="white", linewidth=1.0)
    ax.barh(y, np.where(n > 0, 100 - share, 0), left=share, color=ORANGE, height=0.6, label=lab_rest, edgecolor="white", linewidth=1.0)
    for yi, gi, ni in zip(y, g, n):
        ax.text(102, yi, f"{int(gi):,} of {int(ni):,}" if ni else "none", va="center", ha="left", fontsize=8.5, color=INK)
    ax.set_title(title, fontsize=10, color=INK, loc="left")
    ax.set_xlim(0, 100)
    ax.set_xticks([0, 25, 50, 75, 100])
    ax.set_xticklabels(["0", "25", "50", "75", "100%"])
    ax.xaxis.grid(True, color=GRID, linewidth=0.7)
    ax.set_axisbelow(True)
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.spines["bottom"].set_color(GRID)
    ax.tick_params(axis="x", colors=MUTED, labelsize=8.5, length=0)
    ax.tick_params(axis="y", colors=INK, labelsize=9, length=0)
    ax.legend(loc="upper left", bbox_to_anchor=(0, -0.09), ncol=2, frameon=False, fontsize=8.5, labelcolor=INK, handlelength=1.2)
axes[0].set_yticks(y)
axes[0].set_yticklabels(order["assembly"])
fig.text(0.01, 1.0, "copyproof on HG002 assemblies, scored against the k-mers of T2T-HG002 Q100 v1.1", fontsize=11.5, color=INK, ha="left", va="bottom")
fig.savefig(f"{a.out}.png", dpi=170, bbox_inches="tight", facecolor="white")
Path(f"{a.out}.json").write_text(json.dumps(order[["assembly", "high", "high_true", "truth_assessable", "truth_found"]].to_dict("records"), indent=1))
print(f"\nwrote {a.out}.tsv, .by_size.tsv, .json, .png")
