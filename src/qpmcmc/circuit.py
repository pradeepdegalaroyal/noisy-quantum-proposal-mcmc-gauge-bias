"""Exact transition matrix of a Trotterised, noisy quantum-proposal circuit.

The proposal prepares basis state ``|s>``, evolves it under a palindromic (Strang) Trotter circuit of
``H = (1 - g) alpha E(s)/E_scale + g sum_i X_i`` and measures in the computational basis. Noise is
applied exactly to batched density matrices, so no sampling error enters the transition matrix.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import torch

from . import channels as ch
from .lattice import GaugeModel


@dataclass(frozen=True)
class Proposal:
    """Quantum proposal parameters: mixer weight ``g``, total time ``t``, Trotter steps, scale ``alpha``."""

    g: float = 0.5
    t: float = 4.0
    steps: int = 4
    alpha: float = 4.0


@dataclass(frozen=True)
class NoiseSpec:
    """Per-qubit noise per CNOT touch (``p_dep``, ``gamma``, ``p_z``) and readout flips ``e0``, ``e1``.

    ``e0 = P(read 1 | true 0)`` and ``e1 = P(read 0 | true 1)``. Single-qubit gates cost
    ``one_q_factor`` of one CNOT touch.
    """

    p_dep: float = 0.0
    gamma: float = 0.0
    p_z: float = 0.0
    e0: float = 0.0
    e1: float = 0.0
    one_q_factor: float = 0.1

    def with_(self, **kw) -> "NoiseSpec":
        return NoiseSpec(**{**self.__dict__, **kw})


def _trotter_angles(prop: Proposal) -> tuple[float, float]:
    dt = prop.t / prop.steps
    return prop.g * dt, dt * (1 - prop.g) * prop.alpha


def ideal_transition_matrix(model: GaugeModel, energies: np.ndarray, prop: Proposal) -> np.ndarray:
    """Noiseless transition matrix ``|<s'|U|s>|^2`` (dense, float64); symmetric because ``U^T = U``."""
    theta, phi = _trotter_angles(prop)
    e_norm = energies / np.abs(energies).max()
    rx_all = np.array([[1.0 + 0j]])
    for _ in range(model.n_links):
        rx_all = np.kron(rx_all, ch.rx(theta))
    diag = np.exp(-1j * phi * e_norm)
    step = rx_all @ (diag[:, None] * rx_all)
    u = np.linalg.matrix_power(step, prop.steps)
    return np.abs(u.T) ** 2        # row = input state


def _layer_superop(noise: NoiseSpec, scale: float, touches: float, frac: float, twirl: bool) -> np.ndarray:
    f = scale * touches * frac
    s = ch.noise_superop(noise.p_dep * f, noise.gamma * f, noise.p_z * f)
    return ch.pauli_twirl(s) if twirl else s


def _build_ops(model, energies, prop, noise, scales, twirl, placement, device):
    """One Trotter step as a list of ('s', qubit, 2x2x2x2 superop) / ('d', phase matrix) operations."""
    n = model.n_links
    theta, phi = _trotter_angles(prop)
    e_norm = energies / np.abs(energies).max()
    touches = model.cnot_touches()
    cx = lambda a: torch.tensor(a.reshape(2, 2, 2, 2), dtype=torch.complex64, device=device)
    ops: list = []

    def xlayer():
        u = ch.unitary_superop(ch.rx(theta))
        for q in range(n):
            if placement == "post":
                s = _layer_superop(noise, scales[q], noise.one_q_factor, 1.0, twirl) @ u
            else:
                h = _layer_superop(noise, scales[q], noise.one_q_factor, 0.5, twirl)
                s = h @ u @ h
            ops.append(("s", q, cx(s)))

    xlayer()
    if placement == "sandwich":
        for q in range(n):
            ops.append(("s", q, cx(_layer_superop(noise, scales[q], touches[q], 0.5, twirl))))
    ph = np.exp(-1j * phi * (e_norm[:, None] - e_norm[None, :]))
    ops.append(("d", torch.tensor(ph, dtype=torch.complex64, device=device)))
    frac = 1.0 if placement == "post" else 0.5
    for q in range(n):
        ops.append(("s", q, cx(_layer_superop(noise, scales[q], touches[q], frac, twirl))))
    xlayer()
    return ops


def _apply_super(rho: torch.Tensor, s: torch.Tensor, q: int, n: int) -> torch.Tensor:
    b, dim = rho.shape[0], rho.shape[1]
    r = rho.view(b, 1 << (n - 1 - q), 2, 1 << q, 1 << (n - 1 - q), 2, 1 << q)
    return torch.einsum("pqij,zaibcjd->zapbcqd", s, r).reshape(b, dim, dim)


def noisy_transition_matrix(
    model: GaugeModel,
    energies: np.ndarray,
    prop: Proposal,
    noise: NoiseSpec,
    scales: np.ndarray | None = None,
    *,
    twirl: bool = False,
    placement: str = "post",
    device: str = "cpu",
    budget_bytes: float = 1.2e9,
) -> np.ndarray:
    """Gate-noise transition matrix ``T[s, s']`` (readout error is applied by :func:`apply_readout`).

    ``placement='post'`` puts the noise of a layer after the layer (usual circuit-level model);
    ``'sandwich'`` splits it half before and half after (time-symmetric, used as a diagnostic).
    ``twirl=True`` replaces every noise channel by its Pauli twirl (idealised randomised compiling).
    """
    n, dim = model.n_links, model.dim
    scales = np.ones(n) if scales is None else scales
    ops = _build_ops(model, energies, prop, noise, scales, twirl, placement, device) * prop.steps
    chunk = int(max(1, min(dim, budget_bytes // (dim * dim * 8))))
    out = np.empty((dim, dim))
    for start in range(0, dim, chunk):
        b = min(chunk, dim - start)
        rho = torch.zeros(b, dim, dim, dtype=torch.complex64, device=device)
        idx = torch.arange(b, device=device)
        rho[idx, start + idx, start + idx] = 1.0
        for op in ops:
            rho = _apply_super(rho, op[2], op[1], n) if op[0] == "s" else rho * op[1]
        out[start:start + b] = torch.diagonal(rho, dim1=1, dim2=2).real.double().cpu().numpy()
        del rho
    out = np.clip(out, 0.0, None)
    return out / out.sum(axis=1, keepdims=True)


def readout_matrix(n: int, e0: float, e1: float, twirl: bool = False) -> np.ndarray:
    """Confusion matrix ``C[true, read]``; ``twirl`` symmetrises it via random pre-measurement X flips."""
    if twirl:
        e0 = e1 = 0.5 * (e0 + e1)
    c1 = np.array([[1 - e0, e0], [e1, 1 - e1]])
    c = np.array([[1.0]])
    for _ in range(n):
        c = np.kron(c1, c)          # qubit 0 is the least-significant bit
    return c


def apply_readout(t: np.ndarray, n: int, e0: float, e1: float, twirl: bool = False,
                  pad_input: bool = False) -> np.ndarray:
    """Right-multiply the transition matrix by the (optionally flip-twirled) readout confusion matrix.

    ``pad_input`` also left-multiplies by the same symmetric matrix, i.e. the prepared state is
    randomly bit-flipped with the readout flip probability before the circuit. ``C T C`` is symmetric
    whenever ``T`` is, which restores the symmetry that ``T C`` alone lacks (``C`` need not commute with ``T``).
    """
    if e0 == 0 and e1 == 0:
        return t
    c = readout_matrix(n, e0, e1, twirl or pad_input)
    return c @ t @ c if pad_input else t @ c
