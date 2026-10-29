# Research gap: does hardware noise break detailed balance in quantum-proposal sampling?

## State of the art, from the supplied reading list

Quantum simulation for high-energy physics (HEP) is now part of the community planning process, and the
roadmap of Bauer et al. (2023, digital object identifier (DOI) 10.1103/PRXQuantum.4.027001) asks for algorithm,
theory and hardware co-design for lattice-gauge-type problems. A natural quantum primitive for sampling a
Boltzmann distribution is a Markov chain Monte Carlo (MCMC) step whose proposal is produced by a short
time-evolution circuit: prepare the current configuration, evolve under a Hamiltonian built from the target
energy plus a transverse field, measure, and accept or reject with the ordinary Metropolis rule. Because the
noiseless evolution operator is complex-symmetric, the proposal is exactly symmetric and the chain has the
right stationary distribution. **The reading list contains no paper on quantum MCMC itself, so this
construction is described here but not cited.**

The hardware papers in the list show what noise a proposal circuit will really meet: structured errors, not
white depolarising noise. Ma et al. (2023, DOI 10.1038/s41586-023-06438-1) report that most gate errors in
a ytterbium qubit are decays out of the qubit subspace, Sivak et al. (2023, DOI 10.1038/s41586-023-05782-6)
stabilise a logical qubit with engineered dissipation, Evered et al. (2023, DOI 10.1038/s41586-023-06481-y)
reach 99.5% two-qubit fidelity and characterise the remaining error sources, and Burkard et al. (2023,
DOI 10.1103/RevModPhys.95.025003) describe the impact of noise on spin qubits through the filter-function
formalism. Energy relaxation (decay of the excited
state, called T1 decay) is **non-unital**: it pushes every state towards the ground state.

## The gap

Papers on the list judge noisy hardware by gate fidelity, logical coherence or learnability (Jerbi et al., 2023,
DOI 10.1038/s41467-023-36159-y). None asks what the noise does to the **stationary distribution** of a
sampling algorithm. A noisy proposal matrix is generally not symmetric, the Metropolis rule then no longer
satisfies detailed balance, and the chain converges to a biased distribution whose bias no amount of
additional sampling removes. Whether that bias matters for a gauge-theory target, and whether a
hardware-compatible fix exists, is unmeasured.

**Hypotheses (falsifiable), tested on a Z2 gauge-Higgs model with 8 links (8 qubits), simulated exactly:**

1. **H1.** At a nominal noise level (rates per touch by a controlled-NOT (CNOT) gate of 0.2% for depolarising, amplitude damping and
   dephasing, plus asymmetric readout error of 1% and 3%), the stationary energy of the unmitigated quantum
   chain is biased by at least 0.02 thermal standard deviations, which is twice the statistical error of 10^4
   independent samples.
2. **H2.** Pauli-twirling the gate noise and flip-twirling the readout (together, "tailoring") reduces that
   bias at least tenfold, below the 10^4-sample statistical error, without reducing the spectral gap below that
   of the classical single-link-flip chain.
3. **H3 (negative control).** For the pure gauge theory (Higgs hopping zero) the energy is blind to gauge
   orbits, so the same noise leaves the energy essentially unbiased even though the full link-configuration
   distribution is biased.

These hypotheses were written after small exploratory probes on a few instances but before the full runs. In H2, "tailoring" originally meant the two twirls only; a third component (matched input flips) was added after the first result and is reported as a post-hoc extension.

**Success metrics.** Total variation (TV) distance between the exact stationary and Boltzmann distributions;
energy bias in thermal standard deviations; spectral gap of the chain; relative asymmetry of the proposal
matrix. All are computed exactly from the transition matrix, with no sampling noise.

## Related work

The reading list carries no role labels, so roles were assigned here: *recent* = dated 2026-2027, *landmark* =
experimental or methods papers with more than 200 citations, *review* = reviews and roadmaps, *adjacent* =
other platforms or protocols.

| Paper | Year | Role | What it does | What it does NOT cover |
|---|---|---|---|---|
| Bauer et al., DOI 10.1103/PRXQuantum.4.027001 | 2023 | review | Roadmap for quantum simulation in HEP | How noise changes the distribution a sampler returns |
| Jerbi et al., DOI 10.1038/s41467-023-36159-y | 2023 | landmark | Compares parametrised-circuit learning models under near-term resource limits | Monte Carlo sampling; noise-induced bias |
| Fisher et al., DOI 10.1146/annurev-conmatphys-031720-030658 | 2023 | review | Random and monitored circuits as tractable many-body dynamics | Stationary laws of noisy proposal kernels |
| Sivak et al., DOI 10.1038/s41586-023-05782-6 | 2023 | landmark | Error correction beyond break-even with engineered dissipation | Application-level effect of residual noise |
| Evered et al., DOI 10.1038/s41586-023-06481-y | 2023 | landmark | 99.5% parallel two-qubit gates on a neutral-atom array | Algorithmic consequences of the remaining error structure |
| Ma et al., DOI 10.1038/s41586-023-06438-1 | 2023 | landmark | Gates and erasure conversion in an atomic qubit; errors mostly leave the qubit space | A noise model for algorithms: unital or not |
| Burkard et al., DOI 10.1103/RevModPhys.95.025003 | 2023 | review | Spin qubits; noise through filter functions | Effect on sampling algorithms |
| Storz et al., DOI 10.1038/s41586-023-05885-0 | 2023 | adjacent | Loophole-free Bell test with fast, high-fidelity readout | Readout asymmetry in algorithms |
| Azuma et al., DOI 10.1103/RevModPhys.95.045006 | 2023 | adjacent | Repeaters; loss and noise in quantum channels | Computation on noisy processors |
| Hu et al., DOI 10.1038/s42254-023-00588-x | 2023 | review | Teleportation progress and its role in distributed computation | Sampling under noise |
| Das and Dash, DOI 10.1007/978-3-032-35516-4_32 | 2027 | recent | Quantum Hamming encoder and decoder with T-count reporting in Qiskit | Noise effect on output statistics |
| Sbert et al., DOI 10.1007/978-3-032-37976-4_5 | 2027 | recent | Entropy families and a data-processing inequality for quantum divergences | Distances between sampled distributions under noise |
| Mehra and Arnon, DOI 10.1109/JSTQE.2026.3737986 | 2027 | recent | Shows gate fidelity depends on the choice of computational basis under hardware limits | Sampling-algorithm correctness |

Where they meet: the roadmap asks for hardware-aware algorithms, the landmark hardware papers show that the
dominant errors are structured (decay, erasure-like leakage), and none of the reviews' open problems asks
whether a noisy *sampler* still samples the right law.

## Novelty in one sentence

We show, by exact transition-matrix analysis, that the *unitality* and time-symmetric placement of hardware
noise, rather than its strength, decides whether a quantum-proposal Metropolis chain keeps its gauge-Higgs
Boltzmann target, and we quantify how much of the bias a compiler-level Pauli-twirl plus readout flip-twirl removes.

## Scale and scope

Everything is classical simulation of up to 10 qubits with exact density matrices. No hardware data, no
quantum advantage and no claim about larger lattices is made.
