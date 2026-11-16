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
    # A hapmer says "this parent's reads hold the k-mer and the other's do not". The second half
    # is only evidence where the other parent's reads cover the place: where a parent's depth on
    # the window's two-copy k-mers is outside this range of its usual depth (a dropout, or a
    # deletion in that parent's cell line, as at the immunoglobulin loci), the window's sites are
    # not used.
    parent_depth: tuple = (0.7, 1.5)
    min_depth_kmers: int = 20  # two-copy k-mers a window needs for that check
    # A real switch carries the other parent's sites at about the usual density. A few stray
    # sites spread over a long stretch that has almost none of its own are noise.
    min_density: float = 0.1  # of the sequence's own density of sites


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


def readable_windows(w: Windows, p: HapParams) -> np.ndarray:
    """Windows where both parents' reads cover the release at about their usual depth. Windows
    with too few two-copy k-mers to tell (a hemizygous X, for instance) count as readable."""
    ok = np.ones(len(w.contig), dtype=bool)
    if w.par is None:
        return ok
    n2 = w.n[:, 2].astype(float)
    known = n2 >= p.min_depth_kmers
    for i in range(w.par.shape[1]):
        per = w.par[:, i] / np.maximum(n2, 1.0)
        usual = float(np.median(per[known])) if known.any() else 0.0
        if usual > 0:
            rel = per / usual
            ok &= ~(known & ((rel < p.parent_depth[0]) | (rel > p.parent_depth[1])))
    return ok


def audit(w: Windows, gaps: pd.DataFrame, parent_names: tuple[str, str] = ("parent1", "parent2"), p: HapParams | None = None):
    """Phase audit of every sequence.

    Returns (contigs, pieces, switches): per-sequence totals, per-piece calls, and the stretches
    that follow the other parent than their sequence does overall.
    """
    if w.hap is None:
        raise ValueError("the windows carry no hapmer counts: scan with both parents' read sets")
    p = p or HapParams()
    lr = math.log((1.0 - p.err) / p.err)
    readable = readable_windows(w, p)
    crow, prow, srow = [], [], []
    gap_by = {c: list(zip(g["start"], g["end"])) for c, g in gaps.groupby("contig")} if len(gaps) else {}
    for ci in range(len(w.contigs)):
        sl = w.slice_of(ci)
        name = w.contigs.at[ci, "name"]
        L = int(w.contigs.at[ci, "length"])
        a = np.where(readable[sl], w.hap[sl, 0], 0).astype(np.int64)
        b = np.where(readable[sl], w.hap[sl, 1], 0).astype(np.int64)
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
                "switched_bp": 0,  # filled below from the switch segments that are kept
                "discordant_sites": minority,
                "sites_neither": int(w.hap[sl, 2].sum()) if w.hap.shape[1] > 2 else 0,
                "windows_unreadable": int((~readable[sl]).sum()),
            }
        )
        if not assessed:
            continue
        own_density = int(own.sum()) / max(len(own), 1)  # own-parent sites per window
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
            if n_oth / (e2 - s2) < p.min_density * own_density:
                continue
            crow[-1]["switched_bp"] += (min(e2 * w.W, L) - s2 * w.W)
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
