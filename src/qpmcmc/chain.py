"""Metropolis chains built from an arbitrary (possibly non-symmetric) proposal matrix."""

from __future__ import annotations

import numpy as np


def single_flip_proposal(n: int) -> np.ndarray:
    """Classical baseline: flip one uniformly chosen link."""
    dim = 1 << n
    q = np.zeros((dim, dim))
    s = np.arange(dim)
    for k in range(n):
        q[s, s ^ (1 << k)] += 1.0 / n
    return q


def uniform_proposal(n: int) -> np.ndarray:
    """Classical baseline: independent uniform proposal over all configurations."""
    dim = 1 << n
    return np.full((dim, dim), 1.0 / dim)


def metropolis_matrix(q: np.ndarray, pi: np.ndarray) -> np.ndarray:
    """Metropolis kernel that *assumes* a symmetric proposal: ``a = min(1, pi(s')/pi(s))``."""
    acc = np.minimum(1.0, pi[None, :] / pi[:, None])
    p = q * acc
    np.fill_diagonal(p, 0.0)
    np.fill_diagonal(p, 1.0 - p.sum(axis=1))
    return p


def stationary_and_gap(p: np.ndarray) -> tuple[np.ndarray, float, bool]:
    """Stationary distribution, spectral gap ``1 - |lambda_2|`` and a reducibility flag."""
    w, v = np.linalg.eig(p.T)
    order = np.argsort(-np.abs(w))
    w, v = w[order], v[:, order]
    reducible = bool(np.abs(w[1]) > 1 - 1e-9)
    if reducible:
        lazy = 0.5 * (p + np.eye(len(p)))
        for _ in range(40):
            lazy = lazy @ lazy
        pi_st = np.full(len(p), 1.0 / len(p)) @ lazy
    else:
        pi_st = np.abs(np.real(v[:, 0]))
    return pi_st / pi_st.sum(), float(1.0 - np.abs(w[1])), reducible


def proposal_defect(q: np.ndarray) -> float:
    """Relative asymmetry ``||Q - Q^T||_1 / ||Q||_1`` (zero iff the proposal is symmetric)."""
    return float(np.abs(q - q.T).sum() / np.abs(q).sum())


def chain_metrics(q: np.ndarray, pi: np.ndarray, energies: np.ndarray) -> dict[str, float]:
    """Exact bias and mixing of the Metropolis chain driven by proposal matrix ``q``.

    ``energy_bias`` is the stationary-mean energy error in units of the thermal standard deviation of
    the energy, so ``1/sqrt(N)`` is the statistical error of ``N`` independent samples.
    """
    p = metropolis_matrix(q, pi)
    pi_st, gap, reducible = stationary_and_gap(p)
    sigma = float(np.sqrt(pi @ energies**2 - (pi @ energies) ** 2))
    return {
        "tv_bias": float(0.5 * np.abs(pi_st - pi).sum()),
        "energy_bias": float((pi_st - pi) @ energies / sigma),
        "gap": gap,
        "defect": proposal_defect(q),
        "reducible": float(reducible),
    }


def simulate_chain(q: np.ndarray, pi: np.ndarray, steps: int, rng: np.random.Generator) -> np.ndarray:
    """Run the Metropolis chain with sampled proposals and return the empirical distribution."""
    p_cdf = np.cumsum(metropolis_matrix(q, pi), axis=1)
    counts = np.zeros(len(pi))
    s = int(rng.integers(len(pi)))
    u = rng.random(steps)
    for k in range(steps):
        s = min(int(np.searchsorted(p_cdf[s], u[k])), len(pi) - 1)
        counts[s] += 1
    return counts / steps
