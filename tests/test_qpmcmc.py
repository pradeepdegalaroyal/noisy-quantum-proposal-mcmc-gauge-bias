import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from qpmcmc import (NoiseSpec, Proposal, apply_readout, ideal_transition_matrix, ladder,  # noqa: E402
                    metropolis_matrix, noisy_transition_matrix, pauli_twirl, torus2)
from qpmcmc.channels import amplitude_damping, noise_superop, pauli_probabilities  # noqa: E402
from qpmcmc.chain import chain_metrics, single_flip_proposal, stationary_and_gap  # noqa: E402
from qpmcmc.study import make_instance  # noqa: E402

PROP = Proposal()
NOISE = NoiseSpec(p_dep=0.004, gamma=0.01, p_z=0.004, e0=0.01, e1=0.04)


def _trace_preserving(s: np.ndarray) -> bool:
    # vec(rho) -> trace: components 0 and 3 (row-major) must sum to 1 for any input
    return np.allclose(np.array([1, 0, 0, 1]) @ s, [1, 0, 0, 1])


def test_channels_are_trace_preserving_and_twirl_is_pauli():
    s = noise_superop(0.01, 0.02, 0.01)
    assert _trace_preserving(s) and _trace_preserving(amplitude_damping(0.3))
    probs = pauli_probabilities(s)
    assert np.isclose(probs.sum(), 1.0, atol=1e-9) and probs.min() >= 0
    assert _trace_preserving(pauli_twirl(s))


def test_lattice_sizes_and_cnot_touches():
    assert ladder(1).n_links == 4 and ladder(2).n_links == 7 and torus2().n_links == 8
    assert torus2().cnot_touches().sum() == 4 * 6 * 2  # six CNOTs per plaquette, two qubits per CNOT


def test_noiseless_density_matrix_matches_dense_unitary():
    m = ladder(1)
    inst = make_instance(m, 0)
    ideal = ideal_transition_matrix(m, inst.energies, PROP)
    dm = noisy_transition_matrix(m, inst.energies, PROP, NoiseSpec())
    assert np.allclose(ideal, ideal.T, atol=1e-12)
    assert np.abs(dm - ideal).max() < 1e-5


def test_ideal_proposal_gives_unbiased_chain():
    m = ladder(1)
    inst = make_instance(m, 1)
    q = ideal_transition_matrix(m, inst.energies, PROP)
    met = chain_metrics(q, inst.pi, inst.energies)
    assert met["tv_bias"] < 1e-9 and met["defect"] < 1e-9


def test_amplitude_damping_biases_but_twirled_time_symmetric_noise_does_not():
    m = ladder(1)
    inst = make_instance(m, 2)
    raw = noisy_transition_matrix(m, inst.energies, PROP, NOISE, inst.scales)
    sym = noisy_transition_matrix(m, inst.energies, PROP, NOISE, inst.scales, twirl=True, placement="sandwich")
    m_raw = chain_metrics(apply_readout(raw, 4, NOISE.e0, NOISE.e1), inst.pi, inst.energies)
    m_sym = chain_metrics(sym, inst.pi, inst.energies)
    assert m_raw["tv_bias"] > 1e-3
    assert m_sym["tv_bias"] < 1e-5 and m_sym["defect"] < 1e-5


def test_metropolis_rows_are_stochastic_and_flip_baseline_is_exact():
    m = ladder(1)
    inst = make_instance(m, 3)
    p = metropolis_matrix(single_flip_proposal(4), inst.pi)
    assert np.allclose(p.sum(axis=1), 1.0) and p.min() >= 0
    pi_st, gap, reducible = stationary_and_gap(p)
    assert np.allclose(pi_st, inst.pi, atol=1e-9) and gap > 0 and not reducible
