"""A small hidden Markov model over windows: Viterbi path and posterior state probabilities.

Emissions are given as log-likelihoods (windows x states); a window with no evidence has equal
values for every state and simply inherits its neighbours' state.
"""

from __future__ import annotations

import numpy as np
from scipy.special import logsumexp


def transition_matrix(n_states: int, p_enter: float, p_leave: float) -> np.ndarray:
    """State 0 is "concordant". Leaving it has total probability p_enter, split evenly; an event
    state returns to 0 with probability p_leave and never jumps straight to another event state
    (two adjacent events are separated by at least one concordant window or reached via 0)."""
    T = np.zeros((n_states, n_states))
    T[0, 0] = 1.0 - p_enter
    T[0, 1:] = p_enter / (n_states - 1)
    for s in range(1, n_states):
        T[s, s] = 1.0 - p_leave
        T[s, 0] = p_leave
    return T


def viterbi(loglik: np.ndarray, T: np.ndarray, start: np.ndarray | None = None) -> np.ndarray:
    """Most probable state path."""
    n, k = loglik.shape
    if n == 0:
        return np.zeros(0, dtype=np.int8)
    with np.errstate(divide="ignore"):
        logT = np.log(T)
        logp = np.log(start) if start is not None else np.log(np.r_[1.0 - 1e-6, np.full(k - 1, 1e-6 / (k - 1))])
    back = np.zeros((n, k), dtype=np.int8)
    score = logp + loglik[0]
    for t in range(1, n):
        cand = score[:, None] + logT  # from state (rows) to state (cols)
        back[t] = np.argmax(cand, axis=0)
        score = cand[back[t], np.arange(k)] + loglik[t]
    path = np.zeros(n, dtype=np.int8)
    path[-1] = int(np.argmax(score))
    for t in range(n - 1, 0, -1):
        path[t - 1] = back[t, path[t]]
    return path


def posterior(loglik: np.ndarray, T: np.ndarray, start: np.ndarray | None = None) -> np.ndarray:
    """Posterior probability of each state in each window (forward-backward in log space)."""
    n, k = loglik.shape
    if n == 0:
        return np.zeros((0, k))
    with np.errstate(divide="ignore"):
        logT = np.log(T)
        logp = np.log(start) if start is not None else np.log(np.r_[1.0 - 1e-6, np.full(k - 1, 1e-6 / (k - 1))])
    fwd = np.zeros((n, k))
    fwd[0] = logp + loglik[0]
    for t in range(1, n):
        fwd[t] = logsumexp(fwd[t - 1][:, None] + logT, axis=0) + loglik[t]
    bwd = np.zeros((n, k))
    for t in range(n - 2, -1, -1):
        bwd[t] = logsumexp(logT + (loglik[t + 1] + bwd[t + 1])[None, :], axis=1)
    post = fwd + bwd
    post -= logsumexp(post, axis=1, keepdims=True)
    return np.exp(post)


def runs(path: np.ndarray):
    """Yield (start, end, state) for each maximal run of equal states; end is exclusive."""
    n = len(path)
    if n == 0:
        return
    cut = np.flatnonzero(np.diff(path)) + 1
    starts = np.r_[0, cut]
    ends = np.r_[cut, n]
    for a, b in zip(starts, ends):
        yield int(a), int(b), int(path[a])
