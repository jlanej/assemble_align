"""S3: parental hapmer audit.

A hapmer is a k-mer of the assembly that the child inherited from exactly one parent: present in
the child's reads and in that parent's reads, absent from the other parent's. One variant makes up
to k adjacent hapmers, so evidence is counted in sites (runs of adjacent hapmers), not k-mers.

Each sequence is decoded with a two-state hidden Markov model (parent 1 / parent 2). Stretches
that follow the other parent than the rest of their scaffold are phase switches; a whole
gap-delimited piece that does so is a misplaced piece. Windows without hapmers carry no evidence,
so runs of homozygosity are "not assessable", never "consistent".
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
import pandas as pd

from .windows import Windows


@dataclass
class HapParams:
    err: float = 0.02  # probability that a hapmer site points to the wrong parent
    p_switch: float = 1e-3  # prior probability of a phase switch between adjacent windows
    min_sites: int = 4  # sites of the other parent a switch segment needs
    min_ratio: float = 3.0  # and its ratio of other-parent to own-parent sites
    min_piece_sites: int = 20  # sites a piece needs to be assessed


def viterbi_two_state(a: np.ndarray, b: np.ndarray, err: float, p_switch: float) -> np.ndarray:
    """0 where the window follows parent 1, 1 where it follows parent 2."""
    n = len(a)
    if n == 0:
        return np.zeros(0, dtype=np.int8)
    lr = math.log((1.0 - err) / err)
    d = ((a - b) * lr).tolist()  # log-likelihood of state 0 minus state 1
    stay, sw = math.log(1.0 - p_switch), math.log(p_switch)
    s0, s1 = d[0], 0.0
    back = np.zeros((n, 2), dtype=np.bool_)  # True: came from the other state
    for t in range(1, n):
        a0, b0 = s0 + stay, s1 + sw
        a1, b1 = s1 + stay, s0 + sw
        if b0 > a0:
            n0 = b0
            back[t, 0] = True
        else:
            n0 = a0
        if b1 > a1:
            n1 = b1
            back[t, 1] = True
        else:
            n1 = a1
        s0, s1 = n0 + d[t], n1
    path = np.zeros(n, dtype=np.int8)
    cur = 0 if s0 >= s1 else 1
    for t in range(n - 1, -1, -1):
        path[t] = cur
        if back[t, cur]:
            cur = 1 - cur
    return path


def pieces_of(length: int, gaps: list[tuple[int, int]], min_gap: int = 10) -> list[tuple[int, int]]:
    """Split a sequence at runs of at least `min_gap` non-ACGT bases."""
    out, prev = [], 0
    for s, e in sorted(gaps):
        if e - s < min_gap:
            continue
        if s > prev:
            out.append((prev, s))
        prev = e
    if prev < length:
        out.append((prev, length))
    return out


def audit(w: Windows, gaps: pd.DataFrame, parent_names: tuple[str, str] = ("parent1", "parent2"), p: HapParams | None = None):
    """Phase audit of every sequence.

    Returns (contigs, pieces, switches): per-sequence totals, per-piece calls, and the stretches
    that follow the other parent than their sequence does overall.
    """
    if w.hap is None:
        raise ValueError("the windows carry no hapmer counts: scan with both parents' read sets")
    p = p or HapParams()
    lr = math.log((1.0 - p.err) / p.err)
    crow, prow, srow = [], [], []
    gap_by = {c: list(zip(g["start"], g["end"])) for c, g in gaps.groupby("contig")} if len(gaps) else {}
    for ci in range(len(w.contigs)):
        sl = w.slice_of(ci)
        name = w.contigs.at[ci, "name"]
        L = int(w.contigs.at[ci, "length"])
        a = w.hap[sl, 0].astype(np.int64)
        b = w.hap[sl, 1].astype(np.int64)
        ta, tb = int(a.sum()), int(b.sum())
        major = 0 if ta >= tb else 1
        path = viterbi_two_state(a, b, p.err, p.p_switch)
        own, oth = (a, b) if major == 0 else (b, a)
        assessed = ta + tb >= p.min_piece_sites
        minority = int(oth[path == major].sum()) + int(own[path != major].sum())  # sites that disagree with the decoded phase
        crow.append(
            {
                "contig": name, "group": w.contigs.at[ci, "group"], "length": L,
                f"sites_{parent_names[0]}": ta, f"sites_{parent_names[1]}": tb,
                "parent": parent_names[major] if assessed else "not_assessable",
                "purity": round(max(ta, tb) / (ta + tb), 4) if ta + tb else float("nan"),
                "switched_bp": int((path != major).sum()) * w.W if assessed else 0,
                "discordant_sites": minority,
                "sites_neither": int(w.hap[sl, 2].sum()) if w.hap.shape[1] > 2 else 0,
            }
        )
        if not assessed:
            continue
        # switch segments: runs decoded as the other parent
        cut = np.flatnonzero(np.diff(path)) + 1
        starts = np.r_[0, cut]
        ends = np.r_[cut, len(path)]
        for s, e in zip(starts, ends):
            if path[s] == major:
                continue
            n_oth, n_own = int(oth[s:e].sum()), int(own[s:e].sum())
            if n_oth < p.min_sites or n_oth < p.min_ratio * max(n_own, 1):
                continue
            hit = np.flatnonzero(oth[s:e] > 0)
            s2, e2 = s + int(hit[0]), s + int(hit[-1]) + 1  # trim to the first and last informative window
            srow.append(
                {
                    "contig": name, "start": s2 * w.W, "end": min(e2 * w.W, L), "scaffold_parent": parent_names[major],
                    "segment_parent": parent_names[1 - major], "sites_segment_parent": n_oth, "sites_scaffold_parent": n_own,
                    "llr": round((n_oth - n_own) * lr, 1),
                }
            )
        # gap-delimited pieces
        pcs = pieces_of(L, gap_by.get(name, []))
        if len(pcs) > 1:
            for i, (s, e) in enumerate(pcs):
                ws, we = s // w.W, -(-e // w.W)
                n_own, n_oth = int(own[ws:we].sum()), int(oth[ws:we].sum())
                if n_own + n_oth < p.min_piece_sites:
                    call = "not_assessable"
                elif n_oth >= p.min_ratio * max(n_own, 1):
                    call = "other_parent"
                elif n_own >= p.min_ratio * max(n_oth, 1):
                    call = "consistent"
                else:
                    call = "mixed"
                prow.append(
                    {
                        "contig": name, "piece": i, "start": s, "end": e, "scaffold_parent": parent_names[major],
                        "sites_scaffold_parent": n_own, "sites_other_parent": n_oth, "call": call,
                    }
                )
    return pd.DataFrame(crow), pd.DataFrame(prow), pd.DataFrame(srow)
