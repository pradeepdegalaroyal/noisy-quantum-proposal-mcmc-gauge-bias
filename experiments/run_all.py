"""Reproduce every result and figure: ``python experiments/run_all.py`` (``--quick`` for a smoke run)."""

from __future__ import annotations

import argparse
import json
import platform
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import torch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from qpmcmc import ladder, torus2  # noqa: E402
from qpmcmc.chain import simulate_chain, metropolis_matrix, stationary_and_gap  # noqa: E402
from qpmcmc.circuit import Proposal, ideal_transition_matrix  # noqa: E402
from qpmcmc.study import (KAPPA, NOMINAL_NOISE, calibrate, classical_metrics, gate_matrix,  # noqa: E402
                          make_instance, measure, readout_variants)
from qpmcmc import viz  # noqa: E402

torch.set_num_threads(6)
MC_ERR = 0.01   # statistical error (in thermal sigmas) of 10**4 independent samples


def agg(df: pd.DataFrame, by: list[str], cols: list[str]) -> pd.DataFrame:
    g = df.groupby(by, sort=False)[cols]
    out = g.mean().add_suffix("_mean").join(g.std(ddof=1).add_suffix("_std"))
    return out.reset_index()


def main_table(model, seeds, prop, device) -> pd.DataFrame:
    """Compare classical, ideal-quantum and noisy-quantum proposals plus the mitigation ablation."""
    n, noise, rows = model.n_links, NOMINAL_NOISE, []
    unital = noise.with_(gamma=0.0, e0=0.5 * (noise.e0 + noise.e1), e1=0.5 * (noise.e0 + noise.e1))
    for seed in seeds:
        inst = make_instance(model, seed)
        res = dict(classical_metrics(inst))
        res["quantum ideal (no noise)"] = measure(inst, ideal_transition_matrix(model, inst.energies, prop))
        raw = gate_matrix(inst, prop, noise, device)
        tw = gate_matrix(inst, prop, noise, device, twirl=True)
        sym = gate_matrix(inst, prop, noise, device, twirl=True, placement="sandwich")
        res["quantum raw noise"] = measure(inst, readout_variants(raw, n, noise, False))
        res["quantum unital noise only"] = measure(
            inst, readout_variants(gate_matrix(inst, prop, unital, device), n, unital, False))
        res["+ readout flip-twirl"] = measure(inst, readout_variants(raw, n, noise, True))
        res["+ gate Pauli-twirl"] = measure(inst, readout_variants(tw, n, noise, False))
        res["+ both (twirled)"] = measure(inst, readout_variants(tw, n, noise, True))
        res["+ matched input flips (tailored)"] = measure(inst, readout_variants(tw, n, noise, True, True))
        res["tailored + time-symmetric placement"] = measure(inst, readout_variants(sym, n, noise, True, True))
        res["gate-twirl + time-symmetric, no readout"] = measure(inst, sym)
        for name, m in res.items():
            rows.append({"seed": seed, "method": name, **m})
        print(f"  main seed {seed} done", flush=True)
    return pd.DataFrame(rows)


def gamma_sweep(model, seeds, prop, device, gammas) -> pd.DataFrame:
    """Raw versus tailored bias as the amplitude-damping (T1-type) rate grows; also stores raw matrices' metrics."""
    n, rows = model.n_links, []
    for seed in seeds:
        inst = make_instance(model, seed)
        for g in gammas:
            noise = NOMINAL_NOISE.with_(gamma=g)
            raw = gate_matrix(inst, prop, noise, device)
            tw = gate_matrix(inst, prop, noise, device, twirl=True)
            for name, q in (("raw", readout_variants(raw, n, noise, False)),
                            ("tailored", readout_variants(tw, n, noise, True, True))):
                rows.append({"seed": seed, "gamma": g, "variant": name, **measure(inst, q)})
        print(f"  gamma sweep seed {seed} done", flush=True)
    return pd.DataFrame(rows)


def kappa_heatmap(model, seeds, prop, device, gammas, kappas) -> pd.DataFrame:
    """Raw-noise bias over (amplitude damping, Higgs hopping); kappa = 0 is the pure gauge theory."""
    n, rows = model.n_links, []
    for seed in seeds:
        for kappa in kappas:
            inst = make_instance(model, seed, kappa=kappa)
            for g in gammas:
                noise = NOMINAL_NOISE.with_(gamma=g)
                q = readout_variants(gate_matrix(inst, prop, noise, device), n, noise, False)
                rows.append({"seed": seed, "gamma": g, "kappa": kappa, **measure(inst, q)})
        print(f"  heatmap seed {seed} done", flush=True)
    return pd.DataFrame(rows)


def scaling(models, seeds_by_model, prop, device) -> pd.DataFrame:
    rows = []
    for model in models:
        n = model.n_links
        for seed in seeds_by_model[model.name]:
            t0 = time.time()
            inst = make_instance(model, seed)
            res = dict(classical_metrics(inst))
            res["quantum ideal (no noise)"] = measure(inst, ideal_transition_matrix(model, inst.energies, prop))
            raw = gate_matrix(inst, prop, NOMINAL_NOISE, device)
            tw = gate_matrix(inst, prop, NOMINAL_NOISE, device, twirl=True)
            res["quantum raw noise"] = measure(inst, readout_variants(raw, n, NOMINAL_NOISE, False))
            res["+ matched input flips (tailored)"] = measure(inst, readout_variants(tw, n, NOMINAL_NOISE, True, True))
            for name, m in res.items():
                rows.append({"model": model.name, "n": n, "seed": seed, "method": name, **m})
            print(f"  scaling {model.name} seed {seed}: {time.time() - t0:.1f}s", flush=True)
    return pd.DataFrame(rows)


def chain_check(model, prop, device, steps) -> dict:
    """Sampled Metropolis chain versus the exact stationary distribution of the same noisy proposal."""
    inst = make_instance(model, 0)
    q = readout_variants(gate_matrix(inst, prop, NOMINAL_NOISE, device), model.n_links, NOMINAL_NOISE, False)
    pi_st, _, _ = stationary_and_gap(metropolis_matrix(q, inst.pi))
    emp = simulate_chain(q, inst.pi, steps, np.random.default_rng(0))
    return {"steps": steps,
            "exact_tv_bias": float(0.5 * np.abs(pi_st - inst.pi).sum()),
            "empirical_tv_to_target": float(0.5 * np.abs(emp - inst.pi).sum()),
            "empirical_tv_to_exact_stationary": float(0.5 * np.abs(emp - pi_st).sum()),
            "exact_energy_bias_sigma": float((pi_st - inst.pi) @ inst.energies / np.sqrt(
                inst.pi @ inst.energies ** 2 - (inst.pi @ inst.energies) ** 2)),
            "empirical_energy_bias_sigma": float((emp - inst.pi) @ inst.energies / np.sqrt(
                inst.pi @ inst.energies ** 2 - (inst.pi @ inst.energies) ** 2))}


def make_figures(main, sweep, heat, scal, model_name, fig_dir: Path, gammas, kappas) -> None:
    ab = main.groupby("method")["energy_bias"].apply(lambda s: s.abs().mean())
    stages = ["quantum raw noise", "+ readout flip-twirl", "+ both (twirled)", "+ matched input flips (tailored)",
              "tailored + time-symmetric placement"]
    vals = [ab[s] for s in stages]
    labels = ["raw", "readout flip-twirl", "gate Pauli-twirl", "matched input flips", "time-symmetric placement", "remaining bias"]
    labels = [l.replace(" ", chr(10), 1) if len(l) > 10 else l for l in labels]
    viz.waterfall(labels,
                  vals[0], [vals[i + 1] - vals[i] for i in range(4)], fig_dir / "fig1_waterfall_bias_ablation.png",
                  "mean |energy bias| (thermal standard deviations)",
                  f"Where the noise-induced bias goes ({model_name}, nominal noise)")
    s = sweep.assign(a=sweep.energy_bias.abs()).groupby(["variant", "gamma"])["a"].agg(["mean", "std"]).reset_index()
    ser = {}
    for v, label in (("raw", "raw noisy proposal"), ("tailored", "tailored (twirls + matched input flips)")):
        d = s[s.variant == v].sort_values("gamma")
        ser[label] = (d["mean"].to_numpy(), d["std"].fillna(0).to_numpy())
    viz.line_band(np.array(gammas), ser, fig_dir / "fig2_bias_vs_amplitude_damping.png",
                  "amplitude damping probability per CNOT touch", "|energy bias| (thermal standard deviations)",
                  "Bias grows with non-unital (T1-type) noise; tailoring flattens it",
                  hline=(MC_ERR, "statistical error of 10^4 independent samples"))
    h = heat.assign(a=heat.energy_bias.abs()).groupby(["kappa", "gamma"])["a"].mean().unstack("gamma")
    viz.heatmap(h.to_numpy(), list(h.columns), list(h.index), fig_dir / "fig3_heatmap_bias_gamma_kappa.png",
                "amplitude damping per CNOT touch", "Higgs hopping kappa (0 = pure gauge)",
                "Raw-noise energy bias over damping and Higgs hopping", "|energy bias| (sigma)", "{:.3f}")
    sc = scal[scal.method.isin(["quantum raw noise", "+ matched input flips (tailored)"])].assign(a=lambda d: d.energy_bias.abs())
    sizes = sorted(sc.n.unique())
    ser = {}
    for m, label in (("quantum raw noise", "raw noisy proposal"), ("+ matched input flips (tailored)", "tailored proposal")):
        g = sc[sc.method == m].groupby("n")["a"].agg(["mean", "std"]).reindex(sizes)
        ser[label] = (g["mean"].to_numpy(), g["std"].fillna(0).to_numpy())
    viz.line_band(np.array(sizes), ser, fig_dir / "fig4_bias_vs_qubits.png", "number of qubits (links)",
                  "|energy bias| (thermal standard deviations)", "Bias versus lattice size at nominal noise",
                  hline=(MC_ERR, "statistical error of 10^4 independent samples"))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--device", default="auto")
    ap.add_argument("--figures-only", action="store_true", help="redraw figures from saved results/*_raw.csv")
    args = ap.parse_args()
    device = ("cuda" if torch.cuda.is_available() else "cpu") if args.device == "auto" else args.device
    res_dir, fig_dir = ROOT / ("results/quick" if args.quick else "results"), ROOT / ("figures/quick" if args.quick else "figures")
    res_dir.mkdir(parents=True, exist_ok=True)
    fig_dir.mkdir(parents=True, exist_ok=True)
    t_start = time.time()
    if args.figures_only:
        d = {k: pd.read_csv(res_dir / f"{k}_raw.csv") for k in ("main", "gamma_sweep", "heatmap", "scaling")}
        make_figures(d["main"], d["gamma_sweep"], d["heatmap"], d["scaling"], "torus2x2", fig_dir,
                     sorted(d["gamma_sweep"].gamma.unique()), None)
        return

    if args.quick:
        model, seeds, gammas, kappas = ladder(1), [0, 1], [0.0, 0.004, 0.016], [0.0, 0.3]
        scal_models, scal_seeds, chain_steps = [ladder(1), ladder(2)], {"ladder1": [0, 1], "ladder2": [0]}, 20000
        prop, cal_rows = Proposal(), []
    else:
        model, seeds, gammas, kappas = torus2(), [0, 1, 2, 3, 4], [0.0, 0.001, 0.002, 0.004, 0.008, 0.016], [0.0, 0.15, 0.3, 0.45, 0.6]
        scal_models = [ladder(1), ladder(2), torus2(), ladder(3)]
        scal_seeds = {"ladder1": [0, 1, 2], "ladder2": [0, 1, 2], "torus2x2": [0, 1, 2], "ladder3": [0]}
        chain_steps = 300000
        print("calibrating proposal on a held-out instance", flush=True)
        prop, cal_rows = calibrate(model)
        pd.DataFrame(cal_rows).to_csv(res_dir / "calibration.csv", index=False)
    print(f"device={device} model={model.name} n={model.n_links} proposal={prop}", flush=True)

    timings = {}
    for name, fn in (
        ("main", lambda: main_table(model, seeds, prop, device)),
        ("gamma_sweep", lambda: gamma_sweep(model, seeds, prop, device, gammas)),
        ("heatmap", lambda: kappa_heatmap(model, seeds[:3], prop, device, gammas, kappas)),
        ("scaling", lambda: scaling(scal_models, scal_seeds, prop, device)),
    ):
        t0 = time.time()
        print(f"[{name}]", flush=True)
        globals()[f"df_{name}"] = fn()
        timings[name] = time.time() - t0
        globals()[f"df_{name}"].to_csv(res_dir / f"{name}_raw.csv", index=False)
    t0 = time.time()
    chk = chain_check(ladder(1) if args.quick else ladder(2), prop, device, chain_steps)
    timings["chain_check"] = time.time() - t0
    (res_dir / "chain_check.json").write_text(json.dumps(chk, indent=2))

    metric_cols = ["tv_bias", "energy_bias", "gap", "defect"]
    main_df = df_main.assign(energy_bias=df_main.energy_bias)
    summary = agg(main_df.assign(abs_energy_bias=main_df.energy_bias.abs()), ["method"],
                  ["tv_bias", "abs_energy_bias", "gap", "defect"])
    summary.to_csv(res_dir / "main_summary.csv", index=False)
    agg(df_gamma_sweep.assign(abs_energy_bias=df_gamma_sweep.energy_bias.abs()), ["variant", "gamma"],
        ["tv_bias", "abs_energy_bias", "gap"]).to_csv(res_dir / "gamma_sweep_summary.csv", index=False)
    agg(df_heatmap.assign(abs_energy_bias=df_heatmap.energy_bias.abs()), ["kappa", "gamma"],
        ["tv_bias", "abs_energy_bias", "gap"]).to_csv(res_dir / "heatmap_summary.csv", index=False)
    agg(df_scaling.assign(abs_energy_bias=df_scaling.energy_bias.abs()), ["n", "method"],
        ["tv_bias", "abs_energy_bias", "gap"]).to_csv(res_dir / "scaling_summary.csv", index=False)
    make_figures(df_main, df_gamma_sweep, df_heatmap, df_scaling, model.name, fig_dir, gammas, kappas)

    meta = {"quick": args.quick, "device": device, "model": model.name, "n_qubits": model.n_links, "seeds": seeds,
            "proposal": prop.__dict__, "kappa": KAPPA, "nominal_noise": NOMINAL_NOISE.__dict__,
            "stage_seconds": timings, "total_seconds": time.time() - t_start,
            "python": platform.python_version(), "torch": torch.__version__,
            "gpu": torch.cuda.get_device_name(0) if device == "cuda" else None}
    (res_dir / "meta.json").write_text(json.dumps(meta, indent=2))
    print(summary.to_string(index=False))
    print(f"total {time.time() - t_start:.1f}s")


if __name__ == "__main__":
    main()
