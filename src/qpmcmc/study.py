"""Study helpers: random instances, the proposal variants compared, and calibration."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .chain import chain_metrics, single_flip_proposal, uniform_proposal
from .circuit import NoiseSpec, Proposal, apply_readout, ideal_transition_matrix, noisy_transition_matrix
from .lattice import GaugeModel

NOMINAL_NOISE = NoiseSpec(p_dep=0.002, gamma=0.002, p_z=0.002, e0=0.01, e1=0.03)
BETA = 1.0
KAPPA = 0.3


@dataclass
class Instance:
    """One random problem: disordered plaquette couplings, qubit-to-qubit noise scales, target."""

    model: GaugeModel
    couplings: np.ndarray
    scales: np.ndarray
    kappa: float
    beta: float
    energies: np.ndarray
    pi: np.ndarray


def make_instance(model: GaugeModel, seed: int, kappa: float = KAPPA, beta: float = BETA) -> Instance:
    """Couplings ``J_p = 1 + 0.25 xi`` and lognormal qubit noise scales (sigma 0.3, mean 1)."""
    rng = np.random.default_rng(seed)
    j = np.clip(1.0 + 0.25 * rng.standard_normal(len(model.plaquettes)), 0.3, None)
    scales = np.exp(0.3 * rng.standard_normal(model.n_links))
    scales /= scales.mean()
    e = model.energies(j, kappa)
    pi = np.exp(-beta * (e - e.min()))
    return Instance(model, j, scales, kappa, beta, e, pi / pi.sum())


def gate_matrix(inst: Instance, prop: Proposal, noise: NoiseSpec, device: str, *,
                twirl: bool = False, placement: str = "post") -> np.ndarray:
    return noisy_transition_matrix(inst.model, inst.energies, prop, noise, inst.scales,
                                   twirl=twirl, placement=placement, device=device)


def measure(inst: Instance, q: np.ndarray) -> dict[str, float]:
    return chain_metrics(q, inst.pi, inst.energies)


def classical_metrics(inst: Instance) -> dict[str, dict[str, float]]:
    n = inst.model.n_links
    return {"classical single-flip": measure(inst, single_flip_proposal(n)),
            "classical uniform": measure(inst, uniform_proposal(n))}


def calibrate(model: GaugeModel, seed: int = 99) -> tuple[Proposal, list[dict]]:
    """Pick (g, t, alpha) maximising the noiseless spectral gap on a held-out instance."""
    inst = make_instance(model, seed)
    rows, best = [], (-1.0, Proposal())
    for g in (0.3, 0.5, 0.7):
        for t in (2.0, 4.0, 8.0):
            for a in (2.0, 4.0, 8.0):
                prop = Proposal(g=g, t=t, steps=4, alpha=a)
                m = measure(inst, ideal_transition_matrix(model, inst.energies, prop))
                rows.append({"g": g, "t": t, "alpha": a, "gap": m["gap"]})
                if m["gap"] > best[0]:
                    best = (m["gap"], prop)
    return best[1], rows


def readout_variants(t: np.ndarray, n: int, noise: NoiseSpec, twirl: bool, pad_input: bool = False) -> np.ndarray:
    return apply_readout(t, n, noise.e0, noise.e1, twirl=twirl, pad_input=pad_input)
