"""Development aid: re-classify the cached segments of a run with the classifier as it is in the
source tree, and summarise (see cache_segments.py). Used for the revision in NOTEBOOK entry j.

  try_classifier.py truth   RUN_DIR
  try_classifier.py planted RUN_DIR plant.expected.bed
  try_classifier.py asm     RUN_DIR NAME.eval.truth_cn_segments_raw.tsv
"""
import sys, pickle, json
from pathlib import Path
import numpy as np, pandas as pd
from copyproof.windows import Windows
from copyproof.model import fit_model
from copyproof.cn import CnParams, assessable, classify_segment
from copyproof import calls as C
from copyproof.truth import match_segments

def build(run):
    run = Path(run)
    w = Windows.load(str(run / "windows.w1000.npz")); model = fit_model(w); p = CnParams()
    ok = assessable(w, model, p)
    cache = pickle.load(open(run / "segcache.pkl", "rb"))
    rows = []
    for ci, a, b, d, llr, parts in cache:
        sl = w.slice_of(ci); name = w.contigs.at[ci, "name"]; L = int(w.contigs.at[ci, "length"])
        def add(a, b, d, llr):
            w0, w1 = sl.start + a, sl.start + b
            row = {"contig": name, "start": a * w.W, "end": min(b * w.W, L), "direction": "gain" if d > 0 else "loss", "llr": round(llr, 2),
                   "n_win": b - a, "n_assessed": int(ok[w0:w1].sum())}
            row.update(classify_segment(w, model, p, w0, w1, ok))
            return row
        row = add(a, b, d, llr)
        if row["kind"] == "IRREGULAR":
            ps = [add(x, y, d, l) for x, y, l in parts]
            if len(ps) > 1 and [r for r in ps if r["kind"] != "IRREGULAR"]:
                rows.extend(ps); continue
        rows.append(row)
    segs = pd.DataFrame(rows)
    if "own_classes" not in segs: segs["own_classes"] = ""
    return run, segs, p

def main():
    mode, run = sys.argv[1], sys.argv[2]
    run, segs, p = build(run)
    iss = C.issues_from_cn(segs, False, None, p.deep_ratio, p.deep_windows)
    iss["direction"] = np.where(iss["class"] == "COLLAPSE", "gain", "loss")
    out = [f"== {mode} {run.parent.name}/{run.name}: segments {len(segs)}", "  kinds: " + str({f"{k}/{d}": int(n) for (k, d), n in segs.groupby(["kind", "direction"]).size().items()}),
           "  S1 issues: " + str({f"{c}/{k}/{r}": int(n) for (c, k, r), n in iss.groupby(["class", "confidence", "role"]).size().items()})]
    if mode == "planted":
        exp = pd.read_csv(sys.argv[3], sep="\t", names=["contig", "start", "end", "signal", "tag"])
        exp[["id", "kind", "size", "stratum"]] = exp["tag"].str.split("|", expand=True)
        by = {c: g for c, g in iss.groupby("contig")}
        res = []
        for r in exp[exp.signal.isin(["EXPANSION", "EXPANSION_PARTNER", "COLLAPSE"])].itertuples(index=False):
            g = by.get(r.contig)
            want = "COLLAPSE" if r.signal == "COLLAPSE" else "EXPANSION"
            cov, hi = 0.0, False
            if g is not None:
                o = g[(g["class"].isin([want, "UNSUPPORTED"] if want == "EXPANSION" else [want])) & (g.start < r.end) & (g.end > r.start)]
                cov = float((np.minimum(o.end, r.end) - np.maximum(o.start, r.start)).clip(lower=0).sum()) / max(r.end - r.start, 1)
                hi = bool((o.confidence == "high").any())
            res.append((r.kind, r.signal, r.stratum, int(r.size), cov >= 0.5, cov >= 0.5 and hi))
        t = pd.DataFrame(res, columns=["kind", "signal", "stratum", "size", "found", "high"])
        tab = t.groupby(["kind", "signal", "stratum", "size"]).agg(found=("found", "sum"), high=("high", "sum"), n=("found", "size")).reset_index()
        out.append("  " + " | ".join(f"{k}/{s[:9]}/{st[:2]}/{sz//1000}k {f}({h})/{n}" for k, s, st, sz, f, h, n in tab.itertuples(index=False)))
        regs = {c: g[["start", "end"]].values for c, g in exp.groupby("contig")}
        near = [(c in regs) and bool(np.any((regs[c][:, 0] - 5000 < e) & (regs[c][:, 1] + 5000 > s))) for c, s, e in zip(iss.contig, iss.start, iss.end)]
        away = iss[~np.array(near, dtype=bool)]
        out.append(f"  S1 issues away from planted edits: {len(away)} (high {int((away.confidence=='high').sum())})")
    if mode == "asm":
        t = pd.read_csv(sys.argv[3], sep="\t")
        mask = pd.read_csv(run / "cannot_assess.bed", sep="\t", names=["contig", "start", "end"])
        mby = {c: g[["start", "end"]].values for c, g in mask.groupby("contig")}
        def mf(c, s, e):
            m = mby.get(c)
            if m is None: return 0.0
            ov = np.minimum(m[:, 1], e) - np.maximum(m[:, 0], s)
            return float(ov[ov > 0].sum()) / max(e - s, 1)
        t["masked"] = [mf(c, s, e) for c, s, e in zip(t.contig, t.start, t.end)]
        iss2, t2 = match_segments(iss, t)
        for d in ("gain", "loss"):
            c = iss2[iss2.direction == d]; tt = t2[(t2.direction == d) & (t2.masked < 0.5)]
            h, l = c[c.confidence == "high"], c[c.confidence == "low"]
            out.append(f"  {d}: high {int(h['true'].sum())}/{len(h)} true, low {int(l['true'].sum())}/{len(l)} true; truth assessable found {int(tt.found.sum())}/{len(tt)}")
        rest = segs[segs.kind.isin(["CN_OR_DEPTH", "DEPTH", "IRREGULAR"])]
        key = set(zip(iss.contig, iss.start, iss.end))
        rest = rest[[k not in key for k in zip(rest.contig, rest.start, rest.end)]]
        rest, _ = match_segments(rest, t)
        out.append("  not reported, true/all: " + str({f"{k}/{d}": f"{int(g['true'].sum())}/{len(g)}" for (k, d), g in rest.groupby(["kind", "direction"])}))
    segs.to_pickle(str(run / "segs_try.pkl"))
    print("\n".join(out))

if __name__ == "__main__":
    main()
