#!/usr/bin/env python3
"""Score a copyproof run against other assemblies of the same person (GUIDE section 8, rungs 5 and 6).

Without a truth genome, other assemblies of the same individual can label an assembly's k-mers
the way a truth genome does in scripts/eval_truth.py: scan the release with each other assembly's
exact k-mer counts as an extra read set (cpk count --exact true), and every tracked k-mer carries
the number of copies that assembly holds. These labels are a proxy, not truth: two assemblies can
share an error, and where they differ from the release either side can be the one that is wrong.

Label sets, for two reference assemblies A and B:

  A, B     the release against one reference: the rule of eval_truth.py, unchanged
  both     a k-mer is contradicted when A and B hold the same number of copies of it and that
           number differs from the release's
  either   (issues only) confirmed by A or by B
  truth    (with --truth-label) the truth genome's own labels, to measure the proxy where a
           truth genome exists

As in eval_truth.py: a discordance is a run of windows in which at least half of the k-mers are
contradicted in one direction (gain: the labels hold more copies than the release; loss: fewer);
an issue is confirmed when a discordance of its direction covers at least half of it; a
discordance is found when copy-number issues cover at least half of it; it is assessable when
under half of it lies in the run's cannot-assess mask.

  eval_assemblies.py --run RUN_DIR --refs A B --out PREFIX [--truth-label truth] [--scan PREFIX]

Writes PREFIX.cn_scores.tsv, PREFIX.cn_issues.tsv, PREFIX.segments.tsv, PREFIX.missing_clusters.tsv,
PREFIX.unread_kmers.tsv, PREFIX.json and, with --truth-label, PREFIX.calibration.tsv.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from copyproof.track import Scan  # noqa: E402
from copyproof.truth import TruthWindows, match_segments, truth_segments  # noqa: E402

CN_CLASSES = ["COLLAPSE", "EXPANSION", "UNSUPPORTED"]
SIZES = [("under 5 kb", 0, 5_000), ("5 to 20 kb", 5_000, 20_000), ("20 to 100 kb", 20_000, 100_000),
         ("over 100 kb", 100_000, 10**12), ("5 kb and more", 5_000, 10**12), ("all", 0, 10**12)]
KEYS = ("contig", "start", "n", "gain", "loss", "absent", "net")


def label_windows(sc: Scan, sets: dict, sample: int | None, W: int, cap: int):
    """One pass over the track. `sets` maps a label-set name to the read-set indices it uses
    (one index: that assembly alone; two: both must agree). Returns per-set window labels and
    per-set counts of assembly k-mers that no read of the sample holds."""
    acc = {name: {k: [] for k in KEYS} for name in sets}
    unread = {name: {"positions": 0, "unread": 0, "absent": 0, "unread_and_absent": 0} for name in sets}
    for ci in range(len(sc.contigs)):
        L = int(sc.contigs.at[ci, "length"])
        nw = -(-L // W)
        pos, m, c = sc.records(ci)
        mt = m.sum(axis=1, dtype=np.int64)
        elig = mt <= cap
        win = pos[elig] // W
        me = mt[elig]
        ce = c[elig]
        no_reads = (ce[:, sample] == 0) if sample is not None else None
        for name, idx in sets.items():
            t = ce[:, idx[0]].astype(np.int64)
            if len(idx) == 1:
                agree = np.ones(len(t), dtype=bool)
                absent = t == 0
            else:
                u = ce[:, idx[1]].astype(np.int64)
                agree = t == u
                absent = (t == 0) & (u == 0)
            d = np.where(agree, t - me, 0)
            o = acc[name]
            o["contig"].append(np.full(nw, ci, dtype=np.int32))
            o["start"].append(np.arange(nw, dtype=np.int64) * W)
            o["n"].append(np.bincount(win, minlength=nw).astype(np.int32))
            o["gain"].append(np.bincount(win[d > 0], minlength=nw).astype(np.int32))
            o["loss"].append(np.bincount(win[d < 0], minlength=nw).astype(np.int32))
            o["absent"].append(np.bincount(win[absent], minlength=nw).astype(np.int32))
            o["net"].append(np.bincount(win, weights=np.clip(d, -3, 3), minlength=nw).astype(np.float32))
            if no_reads is not None:
                s = unread[name]
                s["positions"] += int(len(t))
                s["unread"] += int(no_reads.sum())
                s["absent"] += int(absent.sum())
                s["unread_and_absent"] += int((no_reads & absent).sum())
    tws = {name: TruthWindows(W=W, cap=cap, **{k: np.concatenate(v) for k, v in o.items()}) for name, o in acc.items()}
    return tws, unread


def masked_fraction(mask: pd.DataFrame, seg: pd.DataFrame) -> list[float]:
    mby = {c: g[["start", "end"]].values for c, g in mask.groupby("contig")}
    out = []
    for c, s, e in zip(seg["contig"], seg["start"], seg["end"]):
        m = mby.get(c)
        if m is None:
            out.append(0.0)
            continue
        ov = np.minimum(m[:, 1], e) - np.maximum(m[:, 0], s)
        out.append(float(ov[ov > 0].sum()) / max(e - s, 1))
    return out


def match(calls: pd.DataFrame, truth: pd.DataFrame):
    """match_segments, with empty tables allowed."""
    if not len(calls) or not len(truth):
        return (calls.assign(truth_overlap=0.0, true=False), truth.assign(call_overlap=0.0, found=False))
    return match_segments(calls, truth)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run", required=True)
    ap.add_argument("--refs", nargs=2, required=True, metavar=("A", "B"), help="read-set labels of the two other assemblies")
    ap.add_argument("--out", required=True)
    ap.add_argument("--truth-label", help="read-set label of a truth genome's k-mers, if the scan holds one")
    ap.add_argument("--sample-label", default="sample")
    ap.add_argument("--scan", help="scan prefix that carries the labels (default: RUN_DIR/scan)")
    a = ap.parse_args()

    run = Path(a.run)
    sc = Scan(a.scan or str(run / "scan"))
    summary = json.loads((run / "summary.json").read_text())
    W, cap = int(summary["window"]), int(summary["multiplicity_cap"])
    A, B = a.refs
    sets = {A: (sc.reads.index(A),), B: (sc.reads.index(B),), "both": (sc.reads.index(A), sc.reads.index(B))}
    if a.truth_label:
        sets["truth"] = (sc.reads.index(a.truth_label),)
    sample = sc.reads.index(a.sample_label) if a.sample_label in sc.reads else None
    tws, unread = label_windows(sc, sets, sample, W, cap)

    mask_path = run / "cannot_assess.bed"
    mask = (pd.read_csv(mask_path, sep="\t", names=["contig", "start", "end"]) if mask_path.stat().st_size
            else pd.DataFrame(columns=["contig", "start", "end"]))
    issues = pd.read_csv(run / "issues.tsv", sep="\t")
    iss = issues[(issues["channel"] == "S1") & issues["class"].isin(CN_CLASSES)].copy()
    iss["direction"] = np.where(iss["class"] == "COLLAPSE", "gain", "loss")
    iss["len"] = iss["end"] - iss["start"]
    high = iss[iss["confidence"] == "high"]

    # ---- copy number: discordances per label set, issues scored against each ---------------
    segs = {}
    for name, tw in tws.items():
        t = truth_segments(tw, sc.contigs)
        t["masked"] = masked_fraction(mask, t)
        t["len"] = t["end"] - t["start"]
        scored, t = match(iss, t)
        iss[f"overlap_{name}"] = scored["truth_overlap"].values
        iss[f"confirmed_{name}"] = scored["true"].values
        _, th = match(high, t.drop(columns=["call_overlap", "found"]))
        t["found_high"] = th["found"].values
        segs[name] = t
    iss["confirmed_either"] = iss[f"confirmed_{A}"] | iss[f"confirmed_{B}"]

    rows = []
    for name in list(sets) + ["either"]:
        t = segs.get(name)
        for direction in ("gain", "loss", "both directions"):
            c = iss if direction == "both directions" else iss[iss["direction"] == direction]
            for size, lo, hi in SIZES:
                cc = c[(c["len"] >= lo) & (c["len"] < hi)]
                ch, cl = cc[cc["confidence"] == "high"], cc[cc["confidence"] == "low"]
                row = {"labels": name, "direction": direction, "size": size}
                if t is not None:
                    td = t if direction == "both directions" else t[t["direction"] == direction]
                    td = td[(td["len"] >= lo) & (td["len"] < hi)]
                    ok = td[td["masked"] < 0.5]
                    row.update(discordances=len(td), assessable=len(ok), found=int(ok["found"].sum()),
                               found_high=int(ok["found_high"].sum()))
                col = f"confirmed_{name}"
                row.update(issues=len(cc), issues_confirmed=int(cc[col].sum()), high=len(ch), high_confirmed=int(ch[col].sum()),
                           low=len(cl), low_confirmed=int(cl[col].sum()))
                rows.append(row)
    scores = pd.DataFrame(rows)
    for col in ("discordances", "assessable", "found", "found_high"):
        scores[col] = scores[col].astype("Int64")
    scores.to_csv(f"{a.out}.cn_scores.tsv", sep="\t", index=False)
    keep = ["id", "contig", "start", "end", "class", "role", "confidence", "score", "direction", "len"]
    iss[keep + [c for c in iss.columns if c.startswith(("overlap_", "confirmed_"))]].to_csv(f"{a.out}.cn_issues.tsv", sep="\t", index=False)
    pd.concat([t.assign(labels=name) for name, t in segs.items()], ignore_index=True).to_csv(f"{a.out}.segments.tsv", sep="\t", index=False)

    # ---- with a truth genome: how does the proxy compare with it? ---------------------------
    if a.truth_label:
        rows = []
        for name in [A, B, "both", "either"]:
            for conf in ("high", "low"):
                c = iss[iss["confidence"] == conf]
                p, t = c[f"confirmed_{name}"], c["confirmed_truth"]
                rows.append({"what": "issues", "labels": name, "group": conf, "n": len(c), "by_truth": int(t.sum()),
                             "by_proxy": int(p.sum()), "by_both": int((p & t).sum()), "proxy_only": int((p & ~t).sum()),
                             "truth_only": int((~p & t).sum()), "by_neither": int((~p & ~t).sum())})
        # the proxy's own discordances: is each one a discordance by the truth as well?
        tt = segs["truth"][["contig", "start", "end", "direction"]]
        for name in [A, B, "both"]:
            d = segs[name]
            d = d[d["masked"] < 0.5]
            real, _ = match(d[["contig", "start", "end", "direction", "len", "found"]], tt)
            for size, lo, hi in SIZES:
                r = real[(real["len"] >= lo) & (real["len"] < hi)]
                rows.append({"what": "assessable discordances", "labels": name, "group": size, "n": len(r),
                             "by_truth": int(r["true"].sum()), "by_proxy": len(r), "by_both": int(r["true"].sum()),
                             "proxy_only": int((~r["true"]).sum()), "found": int(r["found"].sum()),
                             "found_and_by_truth": int((r["found"] & r["true"]).sum())})
        pd.DataFrame(rows).to_csv(f"{a.out}.calibration.tsv", sep="\t", index=False)

    # ---- read sequence the release lacks: do the other assemblies hold it? ------------------
    mc = run / "missing_clusters.tsv"
    miss_rows = []
    if mc.is_file():
        c = pd.read_csv(mc, sep="\t")
        if len(c) and {f"frac_in_{A}", f"frac_in_{B}"} <= set(c.columns):
            c["strength"] = np.where((c["n_reads_locus"] >= 6) & (c["tracked_kmers"] >= 5), "strong", "weak")  # as in the issue list
            ha, hb = c[f"frac_in_{A}"] >= 0.5, c[f"frac_in_{B}"] >= 0.5
            c["held_both"], c["held_either"] = ha & hb, ha | hb
            ht = c[f"frac_in_{a.truth_label}"] >= 0.5 if a.truth_label else None
            for (kind, strength), g in c.groupby(["kind", "strength"]):
                i = g.index
                row = {"kind": kind, "strength": strength, "clusters": len(g), f"held_by_{A}": int(ha[i].sum()),
                       f"held_by_{B}": int(hb[i].sum()), "held_by_both": int(g["held_both"].sum()),
                       "held_by_either": int(g["held_either"].sum())}
                if ht is not None:
                    row.update(held_by_truth=int(ht[i].sum()), both_and_truth=int((g["held_both"] & ht[i]).sum()),
                               either_and_truth=int((g["held_either"] & ht[i]).sum()))
                miss_rows.append(row)
    pd.DataFrame(miss_rows).to_csv(f"{a.out}.missing_clusters.tsv", sep="\t", index=False)

    # ---- assembly k-mers that no read holds: do the other assemblies lack them too? ---------
    urows = [{"labels": name, **s, "unread_share_absent": round(s["unread_and_absent"] / max(s["unread"], 1), 4),
              "absent_share_unread": round(s["unread_and_absent"] / max(s["absent"], 1), 4)} for name, s in unread.items()]
    pd.DataFrame(urows).to_csv(f"{a.out}.unread_kmers.tsv", sep="\t", index=False)

    pick = scores[(scores["direction"] == "both directions") & (scores["size"] == "all")].set_index("labels")
    res = {"run": str(run), "refs": [A, B], "truth_label": a.truth_label, "window": W, "cap": cap,
           "high_confidence_cn_issues": int(len(high)),
           "high_confirmed": {n: int(pick.at[n, "high_confirmed"]) for n in pick.index},
           "assessable_discordances": {n: int(pick.at[n, "assessable"]) for n in pick.index if n != "either"},
           "found": {n: int(pick.at[n, "found"]) for n in pick.index if n != "either"}}
    Path(f"{a.out}.json").write_text(json.dumps(res, indent=1))
    pd.set_option("display.width", 250)
    pd.set_option("display.max_rows", 500)
    print(json.dumps(res, indent=1))
    print(scores[scores["size"].isin(["all", "5 kb and more", "under 5 kb"])].to_string(index=False))


if __name__ == "__main__":
    main()
