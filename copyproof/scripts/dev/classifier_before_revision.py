"""The segment classifier as it was before the revision of NOTEBOOK entry j, kept so that the
revision can be re-run (scripts/dev/sensitivity.py). Not used by copyproof.

It differs from copyproof.cn.classify_segment in three ways: every k-mer class was allowed the
full regional wander of depth on its own, under the concordant and the whole-copy hypothesis; a
hypothesis fitted only if no class at all was further than 3.5 sd from it; and depth factors were
tried in steps of 0.125.
"""

from __future__ import annotations

import numpy as np

from copyproof.cn import ADD_GRID, CnParams, _fmt_classes
from copyproof.model import Model
from copyproof.windows import Windows

MULT_GRID = (0.0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.125, 1.25, 1.375, 1.5, 1.75, 2.0, 2.5, 3.0)


def classify_segment_before(w: Windows, model: Model, p: CnParams, w0: int, w1: int, ok: np.ndarray,
                     base: np.ndarray | None = None, known: np.ndarray | None = None) -> dict:
    """Decide what kind of departure a segment shows, from the summed evidence of its windows."""
    use_w = ok[w0:w1]
    n = w.n[w0:w1][use_w].astype(float)
    S = w.S[w0:w1][use_w].astype(float)
    V = w.V[w0:w1][use_w].astype(float)
    n_sp = w.n_sp[w0:w1][use_w].astype(float)
    S_sp = w.S_sp[w0:w1][use_w].astype(float)
    cap = w.cap
    f = np.repeat(model.gc_factor(w.gc[w0:w1][use_w].astype(float))[:, None], cap + 1, axis=1)
    psi_class = np.full(cap + 1, model.psi)
    if base is not None:
        f[:, 3:] *= base[w0:w1][use_w][:, None]
        if known is not None and known[w0:w1][use_w].mean() >= 0.5:
            psi_class[3:] = min(model.psi, p.psi_local)  # the regional wander is taken out by the baseline
    j = np.arange(cap + 1, dtype=float)
    om = model.omega_vec(cap)
    N = n.sum(axis=0)
    F = (f * n).sum(axis=0)  # depth-weighted k-mer counts: E = lam * copies * F
    FV = (f * V).sum(axis=0)
    Sj = S.sum(axis=0)
    info = (N >= p.min_seg_class) & (j >= 1)
    res = {
        "kind": "IRREGULAR", "delta": 0, "factor": 1.0, "llr_cn": 0.0, "llr_depth": 0.0, "n_kmers": int(N[1:].sum()),
        "classes": "", "own": "na", "own_llr": 0.0, "dominant_m": 0, "copies_dominant": float("nan"),
        "fit_cn": float("nan"), "fit_depth": float("nan"), "share_cn": 0.0,
    }
    if not info.any():
        return res
    js = j[info]
    lam = model.lam

    Ninfo = N[info]
    copies = Sj[info] / (lam * F[info])
    di = int(np.argmax(Ninfo))

    def fit(expected, rel):
        """Log-likelihood and z-score of each informative class at the expected copies, with a
        relative uncertainty `rel` of the expectation."""
        expected = np.asarray(expected, dtype=float)
        mu = lam * (expected + model.eps)
        E = mu * F[info]
        psi = np.where(expected == 0.0, 1.0, rel)
        sd = np.sqrt(om[info] * mu * FV[info] + (psi * E) ** 2 + 1.0)
        z = (Sj[info] - E) / sd
        return -np.log(sd) - 0.5 * (p.nu + 1.0) * np.log1p(z * z / p.nu), z

    wide = psi_class[info]  # how far depth wanders by region
    l_null, z_null = fit(js, wide)
    null_ll = float(l_null.sum())
    # additive: each class either follows the change or stays concordant
    cands = []
    for d in ADD_GRID:
        l, z = fit(np.maximum(js + d, 0.0), wide)
        part = l > l_null
        if not part.any():
            continue
        tot = float(np.where(part, l, l_null).sum())
        worst = float(np.abs(np.where(part, z, z_null)).max())
        share = float(Ninfo[part].sum() / Ninfo.sum())
        # the commonest class must itself follow the change, and land near a whole copy number
        dom_ok = bool(part[di]) and abs(copies[di] - max(js[di] + d, 0.0)) <= p.dominant_tol
        cands.append((worst <= p.max_z and dom_ok, tot, d, worst, share, dom_ok))
    best_add, best_add_ll, add_worst, add_share, add_dom = 0, -np.inf, np.inf, 0.0, False
    if cands:
        _, best_add_ll, best_add, add_worst, add_share, add_dom = max(cands, key=lambda t: (t[0], t[1]))
    # multiplicative: one factor for every class. Whatever scales all the k-mers of a stretch
    # alike (uneven depth, a fraction of cells) scales them by the same factor, so the classes
    # must agree with it much more closely than depth agrees with its expectation.
    best_mult, best_mult_ll, mult_worst = 1.0, -np.inf, np.inf
    for r in MULT_GRID:
        _, z = fit(js * r, p.psi_local)
        l, _ = fit(js * r, wide)
        cand = (float(np.abs(z).max()) <= p.max_z, float(l.sum()))
        if cand > (mult_worst <= p.max_z, best_mult_ll):
            best_mult, best_mult_ll, mult_worst = r, float(l.sum()), float(np.abs(z).max())
    dom = int(js[di])
    res.update(
        delta=int(best_add), factor=float(best_mult), llr_cn=round(best_add_ll - null_ll, 2),
        llr_depth=round(best_mult_ll - null_ll, 2), classes=_fmt_classes(js, copies), dominant_m=dom,
        copies_dominant=round(float(copies[di]), 3), fit_cn=round(add_worst, 2),
        fit_depth=round(mult_worst, 2), share_cn=round(add_share, 3),
    )
    add_ok = add_dom and add_worst <= p.max_z and add_share >= p.min_share and best_add_ll - null_ll >= p.margin
    mult_ok = mult_worst <= p.max_z and best_mult_ll - null_ll >= p.margin
    # Short stretches that the reads hardly cover at all are common in a correct assembly (sequence
    # that the sequencing chemistry reads poorly), so only a long one counts as absent sequence.
    absent = mult_ok and best_mult <= p.absent_factor and int(use_w.sum()) >= p.absent_min_windows
    if add_ok and mult_ok:
        diff = best_add_ll - best_mult_ll
        res["kind"] = "CN" if diff >= p.margin else (("ABSENT" if absent else "DEPTH") if diff <= -p.margin else "CN_OR_DEPTH")
    elif add_ok:
        res["kind"] = "CN"
    elif absent:
        res["kind"] = "ABSENT"  # every class far below its depth: not bias, the reads lack this sequence
    elif mult_ok:
        res["kind"] = "DEPTH"
    else:
        res["kind"] = "IRREGULAR"
    # Does this sequence's own group carry the change? k-mers found only in this group follow an
    # additive change when the change is here, and stay concordant when it is on the homologue.
    if best_add != 0 and res["kind"] in ("CN", "CN_OR_DEPTH"):
        Nsp = n_sp.sum(axis=0)
        sel = (Nsp >= max(10, p.min_seg_class // 2)) & (j >= 1)
        if sel.any():
            Fsp = (f * n_sp).sum(axis=0)[sel]
            FVsp = (FV * np.divide(Nsp, np.maximum(N, 1.0)))[sel]
            Ssp = S_sp.sum(axis=0)[sel]

            def ll_sp(copies):
                copies = np.asarray(copies, dtype=float)
                mu = lam * (copies + model.eps)
                E = mu * Fsp
                psi = np.where(copies == 0.0, 1.0, psi_class[sel])
                sd = np.sqrt(om[sel] * mu * FVsp + (psi * E) ** 2 + 1.0)
                z = (Ssp - E) / sd
                return -np.log(sd) - 0.5 * (p.nu + 1.0) * np.log1p(z * z / p.nu)

            own_llr = float((ll_sp(np.maximum(j[sel] + best_add, 0.0)) - ll_sp(j[sel])).sum())
            res["own_llr"] = round(own_llr, 2)
            res["own"] = "own" if own_llr >= p.margin else ("other" if own_llr <= -p.margin else "unresolved")
            res["own_classes"] = _fmt_classes(j[sel], Ssp / (lam * Fsp))
    return res
