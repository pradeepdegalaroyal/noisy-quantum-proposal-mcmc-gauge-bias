"""Exact analysis of noisy quantum-proposal Markov chain Monte Carlo on a Z2 lattice gauge target."""

from .lattice import GaugeModel, ladder, torus2
from .channels import noise_superop, pauli_twirl
from .circuit import NoiseSpec, Proposal, ideal_transition_matrix, noisy_transition_matrix, apply_readout
from .chain import metropolis_matrix, stationary_and_gap, chain_metrics, single_flip_proposal, uniform_proposal

__all__ = [
    "GaugeModel", "ladder", "torus2", "noise_superop", "pauli_twirl", "NoiseSpec", "Proposal",
    "ideal_transition_matrix", "noisy_transition_matrix", "apply_readout", "metropolis_matrix",
    "stationary_and_gap", "chain_metrics", "single_flip_proposal", "uniform_proposal",
]
__version__ = "0.1.0"
