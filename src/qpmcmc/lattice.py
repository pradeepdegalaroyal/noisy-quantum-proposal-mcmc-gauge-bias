"""Small Z2 lattice gauge models: one qubit per link, one four-body term per plaquette."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class GaugeModel:
    """Z2 gauge model defined by its plaquettes (tuples of link indices).

    Link ``l`` is qubit ``l``; basis state ``s`` has bit ``l`` of ``s`` equal to the link value
    (spin ``+1`` for bit 0, ``-1`` for bit 1). The classical energy is
    ``E(s) = -sum_p J_p prod_{l in p} spin_l - kappa sum_l spin_l``.
    """

    name: str
    n_links: int
    plaquettes: tuple[tuple[int, ...], ...]

    @property
    def dim(self) -> int:
        return 1 << self.n_links

    def spins(self) -> np.ndarray:
        """Return the (2**n, n) array of link spins in {+1, -1}."""
        states = np.arange(self.dim)[:, None]
        bits = (states >> np.arange(self.n_links)[None, :]) & 1
        return 1.0 - 2.0 * bits

    def plaquette_values(self) -> np.ndarray:
        """Return the (2**n, n_plaquettes) array of plaquette products."""
        sp = self.spins()
        return np.stack([np.prod(sp[:, list(p)], axis=1) for p in self.plaquettes], axis=1)

    def energies(self, couplings: np.ndarray, kappa: float = 0.0) -> np.ndarray:
        """Classical energy of every basis state.

        ``E = -sum_p J_p P_p - kappa sum_l s_l``: the second term is the Higgs hopping of the Z2
        gauge-Higgs model in unitary gauge (``kappa = 0`` is the pure gauge theory, whose energy does
        not distinguish gauge-equivalent link configurations).
        """
        e = -(self.plaquette_values() @ np.asarray(couplings, dtype=float))
        return e - kappa * self.spins().sum(axis=1)

    def cnot_touches(self) -> np.ndarray:
        """CNOTs touching each qubit when every four-body ZZZZ term is compiled as a CNOT ladder.

        A weight-4 term uses 6 CNOTs: the target qubit is touched by all 6, each control by 2.
        """
        touches = np.zeros(self.n_links)
        for p in self.plaquettes:
            touches[list(p)] += 2.0
            touches[p[-1]] += 4.0
        return touches


def ladder(n_plaquettes: int) -> GaugeModel:
    """Open ``n_plaquettes x 1`` strip of plaquettes: ``3 p + 1`` links."""
    p = n_plaquettes
    h = lambda x, y: y * p + x          # horizontal links: 2 rows of p
    v = lambda x: 2 * p + x             # vertical links: p + 1
    plaq = tuple((h(x, 0), v(x + 1), h(x, 1), v(x)) for x in range(p))
    return GaugeModel(f"ladder{p}", 3 * p + 1, plaq)


def torus2() -> GaugeModel:
    """2 x 2 periodic lattice: 8 links, 4 plaquettes (the product of all plaquettes is +1)."""
    L = 2
    link = lambda x, y, mu: 2 * ((x % L) + L * (y % L)) + mu
    plaq = tuple(
        (link(x, y, 0), link(x + 1, y, 1), link(x, y + 1, 0), link(x, y, 1))
        for y in range(L) for x in range(L)
    )
    return GaugeModel("torus2x2", 2 * L * L, plaq)
