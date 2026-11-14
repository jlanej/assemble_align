"""S1: positional k-mer copy number.

For every window the read counts of the assembly's k-mers are compared with what the assembly
implies. A k-mer that occurs m times in the release should have a read count of about
depth x m. Two kinds of departure are modelled:

additive (d = +1, +2, -1)
    The genome holds d more copies of this window's sequence than the release does, so a k-mer
    of multiplicity m is expected at depth x (m + d). This is what a collapsed or a falsely
    duplicated copy looks like. A copy that is missing or extra need not be identical to the
    copies that remain, so each multiplicity class may either follow d or stay concordant.

multiplicative (factor 0.5, 0.75, 1.25, 1.5)
    Every class is scaled by one factor. This is what uneven sequencing depth looks like, and
    also sequence present in only a fraction of cells. It is not an assembly error by itself.

A hidden Markov model over windows finds the segments; each segment is then classified from the
summed evidence of all its windows. In the model an additive change involves all k-mer classes
or only those the release holds at least 2, 3 or 5 times (the k-mers shared between copies). k-mers unique to one haplotype are what separates the two
kinds: at a real one-copy change they sit at 0, 1 or 2 copies, under a depth factor f at f.

A segment is labelled CN when an integer change fits every class and a single factor does not,
DEPTH when a factor fits and an integer change does not, CN_OR_DEPTH when both fit (typically a
segment with a single class, where one copy fewer and half the depth are the same thing), ABSENT
when every class sits at a quarter of its depth or less, and IRREGULAR when nothing fits. CN and
ABSENT segments are candidate assembly errors. A CN_OR_DEPTH segment is one when something else
supports it: k-mers unique to the other copies of its sequence (see `copies`), k-mers of its own
haplotype, another channel, or its sheer length and depth (`deep_ratio`, `deep_windows`).

Depth wanders by region, but all k-mers of a stretch share that wander. So in the classification
every hypothesis has one regional factor for all classes (prior width: the fitted wander), and
each class may then differ from its expectation by 4% plus sampling noise.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
import pandas as pd

from .model import Model
from .windows import Windows

ADD_GRID = (-1, 1, -2, 2, -3, 3, -4, 4)  # tried in this order; ties keep the smaller change
# depth factors tried in the classification: 0, then 0.05 to 3 in steps of 3% (finer than `psi_local`)
MULT_GRID = np.concatenate([[0.0], np.exp(np.arange(np.log(0.05), np.log(3.0) + 1e-9, 0.03))])


@dataclass
class CnParams:
    add: tuple = (1, -1, 2)  # additive hypotheses of the segmentation
    mult: tuple = (0.0, 0.25, 0.5, 0.75, 1.25, 1.5, 2.0)  # multiplicative hypotheses of the segmentation
    add_from: tuple = (1, 2, 3, 5)  # an additive change may involve only classes of at least this multiplicity
    nu: float = 5.0  # degrees of freedom of the Student-t noise
    min_class: int = 4  # k-mers a class needs to be used in a window
    min_kmers: int = 30  # usable k-mers a window needs to be assessed
    max_high_copy: float = 0.5  # windows with a larger share of k-mers above the multiplicity cap are not assessed
    min_gc_factor: float = 0.5  # nor windows whose GC predicts less than this share of the depth
    min_acgt: float = 0.5  # nor windows with a smaller share of A/C/G/T bases
    min_island: int = 5  # assessable runs shorter than this, between unassessable windows, are dropped
    p_enter: float = 1e-6  # prior probability that an event starts at a window
    p_leave: float = 0.02  # prior probability that an event ends at a window
    min_llr: float = 20.0  # log-likelihood ratio (nats) a segment needs to be reported
    merge_gap: int = 2  # segments of the same direction this many windows apart are joined
    min_seg_class: int = 20  # k-mers a class needs to count in the classification of a segment
    margin: float = 3.0  # log-likelihood margin to prefer one kind of explanation
    max_z: float = 3.5  # a hypothesis fits a segment if its classes are no further than this from it ...
    fit_share: float = 0.9  # ... counting the classes that hold this share of the segment's k-mers
    min_share: float = 0.5  # share of a segment's k-mers that must follow an additive change
    dominant_tol: float = 0.35  # the commonest class must land this close to a whole copy number
    use_baseline: bool = False  # see `local_baseline`: off, because it invents gains where local depth is low
    base_span: int = 10  # windows on either side used for the local depth baseline
    base_min: int = 100  # one- and two-copy k-mers the baseline needs ...
    base_windows: int = 7  # ... in at least this many windows of the span ...
    base_win_min: int = 4  # ... each with at least this many of them
    base_range: tuple = (0.8, 1.25)  # a baseline outside this range is an event, not bias, and is not applied
    psi_local: float = 0.04  # relative uncertainty of one class against another in the same segment
    max_wander: float = 2.0  # the regional factor of a segment may be this many sd of the fitted wander from 1
    absent_factor: float = 0.25  # every class at or below this share of its depth: sequence absent from the reads
    absent_min_windows: int = 20  # ... over at least this many assessed windows
    # An ambiguous segment (a depth factor fits as well as a whole copy) is reported without other
    # support when it is long and deep: a truth genome has none in 5.3 Gb, against 142 that are
    # short or shallow (105 of them dips of one or two windows to half the depth).
    deep_ratio: float = 1.5  # the whole-copy change alters the commonest class's depth by at least this factor
    deep_windows: int = 8  # ... over at least this many assessed windows
    # A whole-copy call on fewer assessed windows than this is of low confidence. Stretches of 1 to
    # 3 kb where one haplotype's own k-mers have no reads are as common in a truth genome (13) as
    # in an ordinary assembly of the same person (18, of which the truth confirms 1): mostly
    # sequence that short reads fail on, not errors.
    high_windows: int = 5

    def states(self):
        """(kind, value, lowest multiplicity that takes part). An additive change may involve every
        k-mer class, or only those shared with other copies: k-mers the release holds at least
        2, 3 or 5 times. Which classes take part is a property of the whole event, so it is part
        of the state and not drawn anew in every window."""
        return ([("N", 0.0, 0)] + [("A", float(d), t) for d in self.add for t in self.add_from]
                + [("M", float(r), 0) for r in self.mult])


def assessable(w: Windows, model: Model, p: CnParams) -> np.ndarray:
    """Windows whose k-mers can support a copy-number statement (the complement is the
    CANNOT_ASSESS mask of this channel)."""
    n = w.n.astype(np.int64)
    use = n >= p.min_class
    use[:, 0] = False
    usable = np.where(use, n, 0).sum(axis=1)
    n_all = w.n_all.astype(np.int64)
    high = 1.0 - n.sum(axis=1) / np.maximum(n_all, 1)
    f = model.gc_factor(w.gc.astype(float))
    ok = (usable >= p.min_kmers) & (high <= p.max_high_copy) & (f >= p.min_gc_factor) & (w.acgt >= p.min_acgt * w.W)
    if p.min_island > 1:
        for ci in range(len(w.contigs)):
            sl = w.slice_of(ci)
            a = ok[sl]
            if not a.any():
                continue
            cut = np.flatnonzero(np.diff(a.astype(np.int8))) + 1
            starts = np.r_[0, cut]
            ends = np.r_[cut, len(a)]
            for s, e in zip(starts, ends):
                # short islands strictly inside the sequence; sequence ends count as unassessable
                if a[s] and e - s < p.min_island:
                    a[s:e] = False
    return ok


def local_baseline(w: Windows, model: Model, ok: np.ndarray, p: CnParams) -> np.ndarray:
    """Local depth relative to the genome-wide model, from the k-mers the release holds once or
    twice, over `base_span` windows on either side.

    Sequencing depth wanders by region beyond what GC predicts, and where a sequence has six or
    eight copies, one copy more or fewer is a smaller step than that wander. But the wander is
    shared by all k-mers of a region, so the low-copy k-mers nearby measure it. The baseline is
    applied to classes of multiplicity 3 and above only: for the low-copy classes themselves it
    would hide the very changes it is measured from. Where it leaves the plausible range of bias
    it is not applied at all (returned as 1), and where there are too few low-copy k-mers it is
    unknown.

    Off by default (`CnParams.use_baseline`). It recovers real collapses at multiplicity 6 to 8
    where a whole region is read at 0.9 of the expected depth, but it rests on the copies of a
    multi-copy k-mer all sharing the local depth. Often they do not: most copies lie elsewhere,
    the k-mer is read at the genome's average depth, and dividing by a locally low baseline makes
    a gain out of nothing. On the truth genome it raised confident calls from 17 to 40.
    """
    base = np.ones(len(w.contig))
    known = np.zeros(len(w.contig), dtype=bool)
    lam_f = model.lam * model.gc_factor(w.gc.astype(float))
    width = 2 * p.base_span + 1
    for ci in range(len(w.contigs)):
        sl = w.slice_of(ci)
        okc = ok[sl]
        n1, n2 = w.n[sl, 1].astype(float), w.n[sl, 2].astype(float)
        N = n1 + n2
        E = lam_f[sl] * (n1 + 2.0 * n2)
        S = w.S[sl, 1].astype(float) + w.S[sl, 2]
        with np.errstate(invalid="ignore", divide="ignore"):
            r = np.where(okc & (N >= p.base_win_min), S / E, np.nan)
        ser = pd.Series(r)
        # the median over the span ignores an event that fills less than half of it, so the
        # baseline next to an event is not pulled towards the event's own depth
        med = ser.rolling(width, center=True, min_periods=p.base_windows).median().values
        cnt = pd.Series(np.where(np.isnan(r), 0.0, N)).rolling(width, center=True, min_periods=1).sum().values
        have = ~np.isnan(med) & (cnt >= p.base_min)
        inside = have & (med >= p.base_range[0]) & (med <= p.base_range[1])
        base[sl] = np.where(inside, med, 1.0)
        known[sl] = inside
    return base, known


def _class_ll(S, n, V, mu, omega, psi, nu):
    """Student-t log-density of class sums S (constants dropped). mu is the expected count of one
    k-mer, V the pair weight that carries the overlap of k-mers, omega the extra-Poisson factor,
    psi the relative uncertainty of the expectation."""
    E = mu * n
    sd = np.sqrt(omega * mu * V + (psi * E) ** 2 + 1.0)
    z = (S - E) / sd
    return -np.log(sd) - 0.5 * (nu + 1.0) * np.log1p(z * z / nu)


def emissions(w: Windows, model: Model, p: CnParams, sl: slice, ok: np.ndarray, base: np.ndarray | None = None) -> np.ndarray:
    """Log-likelihood of each window of `sl` under each state; zero rows where not assessable."""
    n = w.n[sl].astype(float)
    S = w.S[sl].astype(float)
    V = w.V[sl].astype(float)
    lam_w = np.repeat((model.lam * model.gc_factor(w.gc[sl].astype(float)))[:, None], w.cap + 1, axis=1)
    if base is not None:
        lam_w[:, 3:] *= base[sl][:, None]  # classes of multiplicity 3 and above follow the local baseline
    j = np.arange(w.cap + 1, dtype=float)[None, :]
    om = model.omega_vec(w.cap)[None, :]
    use = n >= p.min_class
    use[:, 0] = False
    states = p.states()
    ll = np.zeros((n.shape[0], len(states)))
    l0 = _class_ll(S, n, V, lam_w * (j + model.eps), om, model.psi, p.nu)
    ll[:, 0] = np.where(use, l0, 0.0).sum(axis=1)
    cache = {}
    for si, (kind, v, t) in enumerate(states[1:], start=1):
        if kind == "A":
            if v not in cache:
                copies = np.maximum(j + v, 0.0)
                # what "no copies" leaves behind (stray error k-mers) is known only roughly
                psi = np.where(copies == 0.0, 1.0, model.psi)
                cache[v] = _class_ll(S, n, V, lam_w * (copies + model.eps), om, psi, p.nu)
            l = np.where(j >= t, cache[v], l0)  # classes below t stay concordant
        else:
            psi = 1.0 if v == 0.0 else model.psi
            l = _class_ll(S, n, V, lam_w * (j * v + model.eps), om, psi, p.nu)
        ll[:, si] = np.where(use, l, 0.0).sum(axis=1)
    ll[~ok[sl]] = 0.0
    return ll


def viterbi_star(ll: np.ndarray, p_enter: float, p_leave: float) -> np.ndarray:
    """Viterbi path for a model where state 0 connects to every event state and event states
    connect only to themselves and to state 0. Plain Python floats: fast enough for a genome."""
    n, k = ll.shape
    if n == 0:
        return np.zeros(0, dtype=np.int8)
    stay0 = math.log(1.0 - p_enter)
    enter = math.log(p_enter / (k - 1))
    stay = math.log(1.0 - p_leave)
    leave = math.log(p_leave)
    rows = ll.tolist()
    score = [rows[0][0]] + [rows[0][s] + enter for s in range(1, k)]
    back0 = np.zeros(n, dtype=np.int8)  # predecessor of state 0
    backs = np.zeros((n, k), dtype=np.bool_)  # for event state s: True if it came from itself
    rng = range(1, k)
    for t in range(1, n):
        row = rows[t]
        s0 = score[0]
        best, arg = s0 + stay0, 0
        for s in rng:
            v = score[s] + leave
            if v > best:
                best, arg = v, s
        new = [best + row[0]]
        back0[t] = arg
        from0 = s0 + enter
        bt = backs[t]
        for s in rng:
            a = score[s] + stay
            if a >= from0:
                new.append(a + row[s])
                bt[s] = True
            else:
                new.append(from0 + row[s])
        score = new
    path = np.zeros(n, dtype=np.int8)
    cur = max(range(k), key=lambda s: score[s])
    for t in range(n - 1, -1, -1):
        path[t] = cur
        if t == 0:
            break
        if cur == 0:
            cur = int(back0[t])
        elif not backs[t, cur]:
            cur = 0
    return path


def _direction(state) -> int:
    kind, v = state[0], state[1]
    if kind == "N":
        return 0
    if kind == "A":
        return 1 if v > 0 else -1
    return 1 if v > 1 else -1


def segments_of_path(path: np.ndarray, ll: np.ndarray, ok: np.ndarray, p: CnParams):
    """Event segments as (first window, end window, direction, llr). Runs of event states of one
    direction are joined across short stretches of concordant windows, and across unassessable
    ones.

    An event that involves only k-mers of multiplicity 3 and above says nothing about a window
    that holds none: there its likelihood equals the concordant one, and staying in the event
    costs the model nothing. Such windows are therefore not part of a segment (they count as
    concordant), or an event would run on until the next window that contradicts it."""
    states = p.states()
    dirs = np.array([_direction(s) for s in states], dtype=np.int8)[path]
    n = len(path)
    untested = ok & (np.abs(ll[np.arange(n), path] - ll[:, 0]) < 1e-9)
    dirs = np.where(untested, 0, dirs).astype(np.int8)
    out = []
    i = 0
    while i < n:
        d = int(dirs[i])
        if d == 0:
            i += 1
            continue
        j = i
        llr = 0.0
        end = i
        while j < n:
            if dirs[j] == d:
                llr += ll[j, path[j]] - ll[j, 0]
                end = j + 1
                j += 1
                continue
            # look ahead over a short gap
            k = j
            while k < n and k - j < p.merge_gap and dirs[k] != d and (dirs[k] == 0):
                k += 1
            if k < n and k > j and dirs[k] == d:
                j = k
                continue
            break
        a, b = i, end
        while a < b and not ok[a]:
            a += 1
        while b > a and not ok[b - 1]:
            b -= 1
        if b > a:
            out.append((a, b, d, float(llr)))
        i = max(end, i + 1)
    return out


def state_runs(path: np.ndarray, a: int, b: int, ok: np.ndarray, ll: np.ndarray):
    """Maximal runs of one event state inside windows [a, b), trimmed to assessable windows, as
    (start, end, llr)."""
    out = []
    i = a
    while i < b:
        s = int(path[i])
        j = i
        while j < b and path[j] == s:
            j += 1
        if s != 0:
            x, y = i, j
            while x < y and not ok[x]:
                x += 1
            while y > x and not ok[y - 1]:
                y -= 1
            if y > x:
                out.append((x, y, float((ll[x:y, s] - ll[x:y, 0]).sum())))
        i = j
    return out


def _fmt_classes(js, copies) -> str:
    return ",".join(f"{int(j)}:{c:.2f}" for j, c in zip(js, copies))


def classify_segment(w: Windows, model: Model, p: CnParams, w0: int, w1: int, ok: np.ndarray,
                     base: np.ndarray | None = None, known: np.ndarray | None = None) -> dict:
    """Decide what kind of departure a segment shows, from the summed evidence of its windows.

    Depth wanders by region beyond what GC predicts (by `model.psi`), but the wander is shared by
    all k-mers of a stretch. So every hypothesis has one regional factor g, common to all classes,
    with a prior of that width around 1, and each class may then differ from its expectation only
    by `psi_local` plus sampling noise:

      concordant   class m at g x m
      whole copy   class m at g x (m + d) if it follows the change, at g x m if it does not
      depth factor class m at r x m, for any r (no prior: the factor is whatever it is)
    """
    use_w = ok[w0:w1]
    n = w.n[w0:w1][use_w].astype(float)
    S = w.S[w0:w1][use_w].astype(float)
    V = w.V[w0:w1][use_w].astype(float)
    n_sp = w.n_sp[w0:w1][use_w].astype(float)
    S_sp = w.S_sp[w0:w1][use_w].astype(float)
    cap = w.cap
    f = np.repeat(model.gc_factor(w.gc[w0:w1][use_w].astype(float))[:, None], cap + 1, axis=1)
    if base is not None:
        f[:, 3:] *= base[w0:w1][use_w][:, None]
    psi_class = np.full(cap + 1, model.psi)
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

    def fit(expected):
        """Log-likelihood and z-score of each informative class at the expected copies; `expected`
        may carry leading axes (one row per regional factor or depth factor)."""
        expected = np.asarray(expected, dtype=float)
        mu = lam * (expected + model.eps)
        E = mu * F[info]
        psi = np.where(expected == 0.0, 1.0, p.psi_local)
        sd = np.sqrt(om[info] * mu * FV[info] + (psi * E) ** 2 + 1.0)
        z = (Sj[info] - E) / sd
        return -np.log(sd) - 0.5 * (p.nu + 1.0) * np.log1p(z * z / p.nu), z

    def misfit(z):
        """How far the classes are from a hypothesis: the |z| that classes holding `fit_share` of
        the k-mers stay within. Minor classes are mixtures (a k-mer that lies twice in a surplus
        copy loses two copies, its classmates one), so they may not veto a fit; they still count
        in the likelihood. Works along the last axis."""
        az = np.abs(z)
        order = np.argsort(az, axis=-1)
        cum = np.cumsum(np.take_along_axis(np.broadcast_to(Ninfo, az.shape), order, axis=-1), axis=-1)
        k = (cum < p.fit_share * Ninfo.sum() - 1e-9).sum(axis=-1, keepdims=True)
        return np.take_along_axis(np.take_along_axis(az, order, axis=-1), k, axis=-1)[..., 0]

    # regional factors within two standard deviations of the fitted wander: beyond that, a
    # departure is a depth factor in its own right and not the background of a whole-copy change
    G = np.exp(np.linspace(-p.max_wander * model.psi, p.max_wander * model.psi, 17))
    g = G[:, None]
    prior = -0.5 * (np.log(G) / model.psi) ** 2
    l_null, z_null = fit(g * js[None, :])
    tot_null = l_null.sum(axis=1) + prior
    null_ll = float(tot_null.max())
    # whole-copy change: each class either follows it or stays concordant, at the same regional factor
    cands = []
    for d in ADD_GRID:
        l, z = fit(g * np.maximum(js + d, 0.0)[None, :])
        part = l > l_null
        tot = np.where(part.any(axis=1), np.where(part, l, l_null).sum(axis=1) + prior, -np.inf)
        gi = int(np.argmax(tot))
        if not np.isfinite(tot[gi]):
            continue
        zz = np.where(part[gi], z[gi], z_null[gi])
        worst = float(misfit(zz))
        share = float(Ninfo[part[gi]].sum() / Ninfo.sum())
        # the commonest class must itself follow the change, and land near a whole copy number
        dom_ok = (bool(part[gi, di]) and abs(zz[di]) <= p.max_z
                  and abs(copies[di] / G[gi] - max(js[di] + d, 0.0)) <= p.dominant_tol)
        cands.append((worst <= p.max_z and dom_ok, float(tot[gi]), d, worst, share, dom_ok))
    best_add, best_add_ll, add_worst, add_share, add_dom = 0, -np.inf, np.inf, 0.0, False
    if cands:
        _, best_add_ll, best_add, add_worst, add_share, add_dom = max(cands, key=lambda t: (t[0], t[1]))
    # depth factor: one factor for every class
    lm, zm = fit(MULT_GRID[:, None] * js[None, :])
    tot_m, worst_m = lm.sum(axis=1), misfit(zm)
    fits = (worst_m <= p.max_z) & (np.abs(zm[:, di]) <= p.max_z)
    ri = int(np.argmax(np.where(fits, tot_m, -np.inf))) if fits.any() else int(np.argmax(tot_m))
    best_mult, best_mult_ll, mult_worst = float(MULT_GRID[ri]), float(tot_m[ri]), float(worst_m[ri])
    dom = int(js[di])
    res.update(
        delta=int(best_add), factor=round(best_mult, 3), llr_cn=round(best_add_ll - null_ll, 2),
        llr_depth=round(best_mult_ll - null_ll, 2), classes=_fmt_classes(js, copies), dominant_m=dom,
        copies_dominant=round(float(copies[di]), 3), fit_cn=round(add_worst, 2),
        fit_depth=round(mult_worst, 2), share_cn=round(add_share, 3),
    )
    add_ok = add_dom and add_worst <= p.max_z and add_share >= p.min_share and best_add_ll - null_ll >= p.margin
    mult_ok = bool(fits[ri])
    # Short stretches that the reads hardly cover at all are common in a correct assembly (sequence
    # that the sequencing chemistry reads poorly), so only a long one counts as absent sequence.
    absent = mult_ok and best_mult <= p.absent_factor and int(use_w.sum()) >= p.absent_min_windows
    if add_ok and mult_ok:
        diff = best_add_ll - best_mult_ll
        res["kind"] = "CN" if diff >= p.margin else (("ABSENT" if absent else "DEPTH") if diff <= -p.margin else "CN_OR_DEPTH")
    elif add_ok:
        # no factor fits; but if one is far more likely than the whole-copy change, neither is right
        res["kind"] = "CN" if best_add_ll - best_mult_ll > -p.margin else "IRREGULAR"
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


def call_segments(w: Windows, model: Model, p: CnParams | None = None, progress=None):
    """Segments where read copy number departs from the release, over all sequences.

    Returns (segments, assessable): a table with one row per segment and the boolean window mask.
    """
    p = p or CnParams()
    ok = assessable(w, model, p)
    base, known = local_baseline(w, model, ok, p) if p.use_baseline else (None, None)
    rows = []
    for ci in range(len(w.contigs)):
        sl = w.slice_of(ci)
        ll = emissions(w, model, p, sl, ok, base)
        path = viterbi_star(ll, p.p_enter, p.p_leave)
        name = w.contigs.at[ci, "name"]
        L = int(w.contigs.at[ci, "length"])

        def add(a, b, d, llr):
            w0, w1 = sl.start + a, sl.start + b
            row = {
                "contig": name, "start": a * w.W, "end": min(b * w.W, L), "direction": "gain" if d > 0 else "loss",
                "llr": round(llr, 2), "n_win": b - a, "n_assessed": int(ok[w0:w1].sum()),
            }
            row.update(classify_segment(w, model, p, w0, w1, ok, base, known))
            row.update(ci=ci, w0=w0, w1=w1)
            return row

        for a, b, d, llr in segments_of_path(path, ll, ok[sl], p):
            if llr < p.min_llr:
                continue
            row = add(a, b, d, llr)
            if row["kind"] == "IRREGULAR":
                # a merged segment can mix two different events; try its runs of one state each
                parts = [add(x, y, d, l) for x, y, l in state_runs(path, a, b, ok[sl], ll) if l >= p.min_llr]
                regular = [r for r in parts if r["kind"] != "IRREGULAR"]
                if len(parts) > 1 and regular:
                    rows.extend(parts)
                    continue
            rows.append(row)
        if progress:
            progress(ci, len(w.contigs))
    cols = [
        "contig", "start", "end", "direction", "kind", "delta", "factor", "llr", "llr_cn", "llr_depth", "own", "own_llr",
        "dominant_m", "copies_dominant", "fit_cn", "fit_depth", "share_cn", "classes", "own_classes", "n_kmers", "n_win",
        "n_assessed", "ci", "w0", "w1",
    ]
    df = pd.DataFrame(rows)
    for c in cols:
        if c not in df.columns:
            df[c] = pd.Series(dtype=object)
    return df[cols], ok


def window_ratio(w: Windows, model: Model, cls: int | None = None) -> np.ndarray:
    """Copies per window implied by the reads for one multiplicity class (default: the base class);
    NaN where the class has no k-mers."""
    cls = model.base_class if cls is None else cls
    n = w.n[:, cls].astype(float)
    lam_w = model.lam * model.gc_factor(w.gc.astype(float))
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(n > 0, w.S[:, cls] / (lam_w * n), np.nan)
