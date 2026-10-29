# Method

## Target

Z2 gauge-Higgs model in unitary gauge. One qubit per lattice link, basis state `s`, link spins `s_l = +1`
(bit 0) or `-1` (bit 1):

`E(s) = - sum_p J_p prod_{l in p} s_l - kappa sum_l s_l`,  `pi(s) ~ exp(-beta E(s))`,  `beta = 1`, `kappa = 0.3`.

The first sum runs over plaquettes (four-link loops) with couplings `J_p = 1 + 0.25 xi_p`, clipped below at 0.3.
The second term is the Higgs hopping; `kappa = 0` gives the pure gauge theory. Lattices: open strips of 1, 2 and
3 plaquettes (4, 7 and 10 links) and the 2 x 2 periodic lattice (8 links, the main case).

## Quantum proposal (the standard quantum baseline)

From the current state `s`: prepare `|s>`, apply four palindromic (Strang) Trotter steps of
`H = (1 - g) alpha E(s)/max|E| + g sum_i X_i`, measure in the computational basis, and take the outcome `s'`
as the proposed move. Accept with `min(1, pi(s')/pi(s))`. This acceptance rule is correct only if
`Q(s -> s') = Q(s' -> s)`. For the noiseless circuit `U` is a product of complex-symmetric matrices in a
palindromic order, so `U^T = U` and `Q(s -> s') = |<s'|U|s>|^2` is symmetric. The parameters `(g, t, alpha)` are
chosen by a grid search that maximises the noiseless spectral gap on a held-out instance (seed 99), then frozen.

## Noise model (hardware-realistic, exact)

Noise is applied to batched density matrices, so the transition matrix is exact.

* Each ZZZZ plaquette term is compiled as a controlled-NOT (CNOT) ladder with 6 CNOTs: the target qubit is touched 6 times, each
  control twice. Every touch adds, on that qubit, depolarising `p`, amplitude damping `gamma` (decay of
  `|1>` to `|0>`) and dephasing `p_z`, each 0.2% at the nominal point. A single-qubit gate costs 0.1 of a touch.
* Per-qubit noise scales are lognormal (spread 0.3, mean 1) to mimic calibration spread.
* Readout confusion `e0 = P(read 1 | 0) = 1%`, `e1 = P(read 0 | 1) = 3%`.
* State-preparation error is not modelled; readout is the only classical noise.

## Why noise breaks detailed balance, and the fix

Write the noisy circuit as a mixture over inserted Pauli operators `W_k`. Then
`T(s -> s') = sum_k p_k |<s'|W_k|s>|^2` and, since each unitary is complex-symmetric and `P^T = +-P`,
`T(s' -> s)` is the same sum over the *reversed* insertion pattern. Two conditions make the sum symmetric:

1. **Unital Pauli noise.** Amplitude damping is non-unital and not a Pauli mixture, so it breaks the argument.
   Pauli-twirling (randomised compiling) replaces each channel by the diagonal of its Pauli transfer matrix
   (PTM), a stochastic Pauli channel with `p_X = p_Y`. In the simulator the twirl is idealised: the exact Pauli
   twirl of each layer's channel.
2. **Time-symmetric placement.** The reversed pattern must have the same slot strengths; noise lumped *after*
   each layer does not, noise split half before, half after each layer does.

Readout asymmetry is fixed by random pre-measurement bit flips, which turn `(e0, e1)` into the symmetric
`(e0 + e1)/2`. A symmetric readout matrix `C` still does not commute with `T`, so `T C` is only approximately
symmetric; the residual is measured, not assumed away.

**Tailoring** = gate Pauli-twirl + readout flip-twirl + matched input flips (the prepared state is randomly
bit-flipped with the same symmetric probability as the readout, so the effective proposal is `C T C`, which is
symmetric whenever `T` is). The input flips were added after the first full run showed that the two twirls
alone left a residual just above the target (see the README); both are realisable on hardware.
**Time-symmetric placement** is a diagnostic of the theory (an idealised noise model), not a mitigation.

## Baselines

* Classical single-link-flip Metropolis (symmetric, exact).
* Classical uniform independence proposal (symmetric, exact).
* Quantum proposal, noiseless (exact).
* Quantum proposal with raw noise: the standard quantum baseline on hardware.

## Metrics (all exact)

* Energy bias: `(E_stationary - E_target)/sigma_E` where `sigma_E` is the thermal standard deviation of the
  energy. `1/sqrt(N)` is the statistical error of `N` independent samples, so 0.01 matches `N = 10^4`.
* Total variation (TV) distance between stationary and target distributions.
* Spectral gap `1 - |lambda_2|` of the Metropolis kernel (mixing per step; a quantum step also costs a circuit).
* Proposal asymmetry `||Q - Q^T||_1 / ||Q||_1`.

## Ablation, sweeps, seeds

* **Ablation (waterfall):** raw noise, then readout flip-twirl, gate Pauli-twirl, matched input flips and
  time-symmetric placement added in order; the first two are also reported alone. A row with only unital noise isolates the non-unital part.
* **Noise sweep:** amplitude-damping rate 0 to 1.6% per touch, raw versus tailored (line plot with band).
* **Two-parameter map:** amplitude damping versus Higgs hopping `kappa` (heatmap), which includes the pure gauge
  case.
* **Scaling:** 4, 7, 8 and 10 qubits.
* **Check:** a sampled 300,000-step chain compared with the exact stationary distribution.
* **Seeds:** 5 random instances (couplings and qubit noise scales) for the main table and the sweep, 3 for the
  map, 3 for scaling (1 at 10 qubits). Mean and sample standard deviation are reported.

## Numerics

Density matrices use single precision on the graphics processing unit (GPU); everything downstream is double
precision. The single-precision floor (about 1e-7 in the asymmetry) is far below every reported effect.
