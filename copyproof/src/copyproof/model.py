"""Depth model: depth per copy, its dependence on GC, and how noisy a window's k-mer sum is.

Everything is fitted from the windows themselves, using the commonest multiplicity class (two
copies in a diploid release). The fit is robust (medians and median absolute deviations), so the
small fraction of windows that are truly discordant does not move it.
"""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field

import numpy as np

from .windows import Windows


@dataclass
class Model:
    lam: float  # read depth per copy at the reference GC content
    base_class: int  # multiplicity class the fit used
    gc_x: list = field(default_factory=list)  # GC fractions ...
    gc_f: list = field(default_factory=list)  # ... and relative depth there
    # Var(S) = omega[m] * depth-per-k-mer * V + (psi * E)^2 for the sum S (mean E) of a class with
    # multiplicity m. omega is extra-Poisson noise, one value per class; psi is the slow, relative
    # wander of depth along the genome that the GC curve does not explain.
    omega: list = field(default_factory=list)
    psi: float = 0.05
    eps: float = 0.02  # expected count per k-mer, in units of lam, for sequence absent from the genome
    n_fit: int = 0
    notes: dict = field(default_factory=dict)

    def gc_factor(self, gc: np.ndarray) -> np.ndarray:
        if not self.gc_x:
            return np.ones_like(gc, dtype=float)
        return np.interp(gc, self.gc_x, self.gc_f)

    def omega_vec(self, cap: int) -> np.ndarray:
        """omega for classes 0..cap."""
        om = np.ones(cap + 1)
        k = min(len(self.omega), cap + 1)
        om[:k] = self.omega[:k]
        if k and k < cap + 1:
            om[k:] = self.omega[k - 1]
        return om

    def to_json(self) -> str:
        return json.dumps(asdict(self), indent=1)

    @staticmethod
    def from_json(text: str) -> "Model":
        return Model(**json.loads(text))


def mad_sd(x: np.ndarray) -> float:
    """Standard deviation estimated from the median absolute deviation."""
    med = np.median(x)
    return float(1.4826 * np.median(np.abs(x - med)))


def ok_all(w: Windows, min_kmers: int) -> np.ndarray:
    """Windows with enough base-class k-mers and few non-ACGT bases."""
    base = int(np.argmax(w.n.sum(axis=0)))
    return (w.n[:, base] >= min_kmers) & (w.acgt >= 0.9 * w.W)


def fit_model(w: Windows, min_kmers: int = 100, gc_step: float = 0.01, min_bin: int = 200) -> Model:
    """Fit depth per copy, the GC curve and the variance of window sums."""
    base = int(np.argmax(w.n.sum(axis=0)))
    n = w.n[:, base].astype(float)
    S = w.S[:, base].astype(float)
    ok = (n >= min_kmers) & (w.acgt >= 0.9 * w.W)
    if ok.sum() < 1000:
        ok = n >= max(10, min_kmers // 5)
    r = S[ok] / (base * n[ok] * w.lam)
    gc = w.gc[ok].astype(float)
    # GC curve: median ratio per GC bin, relative to the overall median
    centre = float(np.median(r))
    bins = np.floor(gc / gc_step).astype(int)
    xs, fs = [], []
    for b in np.unique(bins):
        sel = bins == b
        if sel.sum() >= min_bin:
            xs.append((b + 0.5) * gc_step)
            fs.append(float(np.median(r[sel])) / centre)
    lam = w.lam * centre
    model = Model(lam=lam, base_class=base, gc_x=xs, gc_f=fs, n_fit=int(ok.sum()))
    # Variance of a class sum: omega[m] * mu * V + (psi * E)^2, where mu is the expected count of
    # one k-mer and V counts overlapping k-mer pairs (windows.py). omega is read from windows where
    # a class is sparse, so that the relative term is negligible; psi from the dense base class.
    f_all = model.gc_factor(w.gc.astype(float))
    good = ok_all(w, min_kmers)
    om = {}
    for j in range(1, w.cap + 1):
        if j == base:
            continue
        nj = w.n[:, j].astype(float)
        sel = (nj >= 8) & (nj <= 0.25 * np.maximum(n, 1)) & good
        if sel.sum() < 2000:
            continue
        mu = lam * f_all[sel] * j
        z = (w.S[sel, j] - mu * nj[sel]) / np.sqrt(mu * w.V[sel, j])
        om[j] = (int(sel.sum()), float(np.clip(mad_sd(z) ** 2, 1.0, 10.0)))
    # the base class takes the value of the best-sampled sparse class (heterozygous k-mers in a
    # diploid release); classes without enough windows take the largest fitted value
    ref = max(om.values(), key=lambda t: t[0])[1] if om else 1.0
    worst = max((v for _, v in om.values()), default=ref)
    omega_list = [1.0] + [om[j][1] if j in om else (ref if j == base else worst) for j in range(1, w.cap + 1)]
    omega = omega_list[base]
    mu = lam * model.gc_factor(gc) * base
    E = mu * n[ok]
    resid = S[ok] - E
    core = omega * mu * w.V[ok, base]
    order = np.argsort(E)
    nb = max(5, min(40, len(E) // 2000))
    psi2 = []
    for chunk in np.array_split(order, nb):
        v = mad_sd(resid[chunk]) ** 2
        psi2.append((v - float(np.median(core[chunk]))) / float(np.median(E[chunk])) ** 2)
    model.omega = [round(float(x), 4) for x in omega_list]
    model.psi = float(np.sqrt(max(np.median(psi2), 1e-6)))
    model.notes = {
        "ratio_sd_per_window": mad_sd(r / model.gc_factor(gc)) / centre,
        "median_kmers_per_window": float(np.median(n[ok])),
        "median_expected_sum": float(np.median(E)),
        "gc_bins": len(xs),
        "omega_windows": {int(j): int(v[0]) for j, v in om.items()},
        "median_pair_weight_per_kmer": float(np.median(w.V[ok, base] / n[ok])),
    }
    return model
