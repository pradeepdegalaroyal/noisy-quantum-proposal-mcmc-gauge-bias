"""Single-qubit noise channels as 4x4 superoperators acting on row-major vec(rho)."""

from __future__ import annotations

import numpy as np

I2 = np.eye(2, dtype=complex)
PX = np.array([[0, 1], [1, 0]], dtype=complex)
PY = np.array([[0, -1j], [1j, 0]], dtype=complex)
PZ = np.array([[1, 0], [0, -1]], dtype=complex)
PAULIS = (I2, PX, PY, PZ)


def superop(kraus: list[np.ndarray]) -> np.ndarray:
    """Superoperator ``sum_k K_k (x) conj(K_k)`` so that ``vec(K rho K^dag) = S vec(rho)``."""
    return sum(np.kron(k, k.conj()) for k in kraus)


def unitary_superop(u: np.ndarray) -> np.ndarray:
    return np.kron(u, u.conj())


def rx(theta: float) -> np.ndarray:
    """Rotation exp(-i theta X / 2)."""
    return np.cos(theta / 2) * I2 - 1j * np.sin(theta / 2) * PX


def pauli_channel(probs: tuple[float, float, float, float]) -> np.ndarray:
    return superop([np.sqrt(max(p, 0.0)) * s for p, s in zip(probs, PAULIS)])


def depolarizing(p: float) -> np.ndarray:
    """rho -> (1 - p) rho + p I / 2."""
    return pauli_channel((1 - 0.75 * p, p / 4, p / 4, p / 4))


def dephasing(p: float) -> np.ndarray:
    return pauli_channel((1 - p, 0.0, 0.0, p))


def amplitude_damping(gamma: float) -> np.ndarray:
    """T1-type decay |1> -> |0>; the non-unital channel that biases the proposal."""
    k0 = np.array([[1, 0], [0, np.sqrt(1 - gamma)]], dtype=complex)
    k1 = np.array([[0, np.sqrt(gamma)], [0, 0]], dtype=complex)
    return superop([k0, k1])


def noise_superop(p_dep: float, gamma: float, p_z: float) -> np.ndarray:
    """Depolarising, then amplitude damping, then dephasing (rates clipped to valid ranges)."""
    p_dep, gamma, p_z = (float(np.clip(v, 0.0, 1.0)) for v in (p_dep, gamma, p_z))
    return dephasing(p_z) @ amplitude_damping(gamma) @ depolarizing(p_dep)


def pauli_probabilities(s: np.ndarray) -> np.ndarray:
    """Probabilities of the Pauli-twirled version of the channel with superoperator ``s``."""
    lam = np.empty(4)
    for a, sig in enumerate(PAULIS):
        out = (s @ sig.reshape(-1)).reshape(2, 2)
        lam[a] = np.real(np.trace(sig @ out)) / 2
    chi = np.array([[1 if (a == 0 or k == 0 or a == k) else -1 for k in range(4)] for a in range(4)])
    return np.clip(chi.T @ lam / 4, 0.0, None)


def pauli_twirl(s: np.ndarray) -> np.ndarray:
    """Keep only the diagonal of the Pauli transfer matrix: a stochastic (unital) Pauli channel."""
    probs = pauli_probabilities(s)
    return pauli_channel(tuple(probs / probs.sum()))
