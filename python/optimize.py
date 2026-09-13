"""Multi-objective optimisation of the elution phase of the HGMS hIgG process.

Five parameters (elution buffer pH, number of elution cycles, and three step durations that are
shared by all cycles) are searched with NSGA-II for the best trade-off between native hIgG yield
and elution buffer consumption, subject to a product-quality constraint. Loading and washing are
identical for every candidate, so they are simulated once and each candidate is warm-started from
that state. Run with ``python optimize.py``.
"""

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List

import matplotlib.pyplot as plt
import numpy as np
import optuna
import pandas as pd
from scipy.optimize import brentq

import fs3
from experiment import build_experiment, state
from reactions import DENATURATION
from recipe import FlowSheetConfig, RecipeStep, default_recipe
from solutions import build_solutions
from unit_operations import _FLOW_CURVE, _PUMP_PERCENTAGES

OUTPUT_DIR = Path(__file__).resolve().parent.parent / "data" / "optimization"
PREFIX_END = "Wash 3 pause"  # last step shared by all candidates
FILL_PUMP, LOOP_PUMP, REGEN_PUMP = 40.0, 30.0, 25.0
REFERENCE = dict(n_cycles=5, t_fill=52.0, t_resuspend_loop=300.0, t_recapture_loop=200.0)

DISCRETIZATION_FACTOR, TAU_REACTION = 0.2, 0.1
HOLD_TIME_S = 1800.0  # the pooled eluate sits at the elution pH before it is neutralised
MAX_DENATURED = 0.05  # of the pooled eluate
N_TRIALS, N_JOBS, SEED = 600, 12, 0


def flow_rate(pump_percentage: float) -> float:
    return float(np.interp(pump_percentage, _PUMP_PERCENTAGES, _FLOW_CURVE))


def elution_recipe(n_cycles: int, t_fill: float, t_resuspend_loop: float, t_recapture_loop: float) -> List[RecipeStep]:
    """Elution phase only, starting at t = 0. Cycle n is pushed out (and collected) by fill n+1.
    All cycles share the same durations, so the recipe has four degrees of freedom."""
    steps = []
    for n in range(1, n_cycles + 1):
        collect = f"Elution {n - 1}" if n > 1 else "Wash 3"
        steps += [
            RecipeStep(f"Elution {n} fill", t_fill, FILL_PUMP, 0, FlowSheetConfig.LINE, "Buffer 3", collect),
            RecipeStep(f"Elution {n} resuspend", 90, 0, 60, FlowSheetConfig.NO_FLOW, "Water (B5)", None),
            RecipeStep(f"Elution {n} resuspend loop", t_resuspend_loop, LOOP_PUMP, 40, FlowSheetConfig.LOOP, "Water (B5)", None),
            RecipeStep(f"Elution {n} recapture", 60, 0, 0, FlowSheetConfig.NO_FLOW, "Water (B5)", None),
            RecipeStep(f"Elution {n} recapture loop", t_recapture_loop, LOOP_PUMP, 0, FlowSheetConfig.LOOP, "Water (B5)", None),
            RecipeStep(f"Elution {n} pause", 8, 0, 0, FlowSheetConfig.NO_FLOW, "Water (B5)", None),
        ]
    steps.append(RecipeStep("Regeneration", 52, REGEN_PUMP, 0, FlowSheetConfig.LINE, "Buffer 1", f"Elution {n_cycles}"))
    return steps


@dataclass
class Setup:
    """Everything shared by all candidates: solutions, the state after washing, and the amount of
    hIgG on the particles at that point (the denominator of the elution yield)."""

    cs: fs3.ComponentSystem
    rs: fs3.ReactionSystem
    solutions: Dict[str, np.ndarray]
    prefix_state: np.ndarray
    hIgG_captured_g: float
    denaturation: dict

    def buffer_at_pH(self, pH: float) -> np.ndarray:
        """200 mM sodium acetate titrated with HCl to the requested (equilibrated) pH."""
        base = np.zeros(self.cs.n_components)
        base[self.cs.get_idx("H₂O")] = 1000.0 / self.cs["H₂O"].molar_mass
        base[self.cs.get_idx("Na⁺")] = base[self.cs.get_idx("Ac⁻")] = 0.2e3

        def equilibrated(c_hcl):
            b = base.copy()
            b[self.cs.get_idx("H⁺")] = b[self.cs.get_idx("Cl⁻")] = c_hcl
            b, _ = fs3.one_cell_reaction(self.rs, b, 100.0, fs3.SolverType.ERK)
            return b

        def error(c_hcl):
            return -np.log10(equilibrated(c_hcl)[self.cs.get_idx("H⁺")] * 1e-3) - pH

        return equilibrated(brentq(error, 150.0, 400.0, xtol=1e-3))

    @property
    def reference_pH(self) -> float:
        return -np.log10(self.solutions["Buffer 3"][self.cs.get_idx("H⁺")] * 1e-3)


def build_setup(denaturation: dict = DENATURATION) -> Setup:
    """Simulate the shared loading + washing phases once."""
    from components import build_component_system
    from reactions import build_reaction_system

    cs = build_component_system(denaturation=True)
    rs, _ = build_reaction_system(cs, TAU_REACTION, denaturation)
    solutions = build_solutions(cs, rs, fs3.SolverType.ERK)

    prefix = default_recipe()
    prefix = prefix[: 1 + next(i for i, s in enumerate(prefix) if s.name == PREFIX_END)]
    e = build_experiment(prefix, discretization_factor=DISCRETIZATION_FACTOR, tau_reaction=TAU_REACTION,
                         dt_obs=sum(s.t_duration for s in prefix) / 2, solutions=solutions,
                         denaturation=denaturation)
    e.solve()
    y = state(e).reshape(-1, cs.n_components)
    captured = 1000.0 * sum(y[:, cs.get_idx(n)].sum() for n in ("MNP-hIgG", "hIgG", "hIgG-den"))
    return Setup(cs, rs, solutions, state(e), captured, denaturation)


def evaluate(setup: Setup, pH: float, n_cycles: int, t_fill: float, t_resuspend_loop: float,
             t_recapture_loop: float, timeout_seconds: float = 60.0) -> Dict[str, float]:
    """Simulate one elution design. ``reached_end`` is False if the solver hit the wall-clock
    timeout, which happens for the few designs that make the ODE degenerate."""
    recipe = elution_recipe(n_cycles, t_fill, t_resuspend_loop, t_recapture_loop)
    duration = sum(s.t_duration for s in recipe)
    solutions = {**setup.solutions, "Buffer 3": setup.buffer_at_pH(pH)}
    e = build_experiment(recipe, discretization_factor=DISCRETIZATION_FACTOR, tau_reaction=TAU_REACTION,
                         dt_obs=duration / 2, solutions=solutions, denaturation=setup.denaturation,
                         initial_state=setup.prefix_state, timeout_seconds=timeout_seconds)
    e.solve()

    native, denatured = e.fraction_masses("hIgG"), e.fraction_masses("hIgG-den")
    pooled = [f"Elution {n}" for n in range(1, n_cycles + 1)]
    m_native = sum(native[f] for f in pooled)
    m_denatured = sum(denatured[f] for f in pooled)

    # FS³ runs no reactions in the collection vessels, so the standard low-pH hold of the pooled
    # eluate is applied afterwards, with the same kinetics.
    d = setup.denaturation
    surviving = np.exp(-d["k_ref"] * 10 ** ((d["pH_ref"] - pH) / d["pH_per_decade"]) * HOLD_TIME_S)
    m_native, m_denatured = m_native * surviving, m_denatured + m_native * (1 - surviving)

    return {
        "yield": m_native / setup.hIgG_captured_g,
        "buffer_L": 1e3 * n_cycles * t_fill * flow_rate(FILL_PUMP),
        "time_s": duration,
        "native_fraction": m_native / (m_native + m_denatured) if m_native + m_denatured > 0 else 1.0,
        "hIgG_native_g": m_native,
        "hIgG_denatured_g": m_denatured,
        "reached_end": e.t >= duration - 1e-6,
    }


# ============================== optimisation ==============================


def search_space(trial) -> dict:
    """The five decision variables."""
    return dict(
        pH=trial.suggest_float("pH", 2.0, 4.0),
        n_cycles=trial.suggest_int("n_cycles", 1, 5),  # one per available fraction vessel
        t_fill=trial.suggest_float("t_fill", 20.0, 120.0),
        t_resuspend_loop=trial.suggest_float("t_resuspend_loop", 30.0, 600.0),
        t_recapture_loop=trial.suggest_float("t_recapture_loop", 60.0, 300.0),
    )


def optimise(setup: Setup) -> pd.DataFrame:
    """NSGA-II on (yield, buffer), with the quality spec as a constraint. Returns every trial."""

    def objective(trial):
        m = evaluate(setup, **search_space(trial))
        if not m.pop("reached_end"):
            raise optuna.TrialPruned("solver did not reach the end of the recipe")
        trial.set_user_attr("metrics", m)
        trial.set_user_attr("constraint", (1.0 - m["native_fraction"] - MAX_DENATURED,))
        return m["yield"], m["buffer_L"]

    study = optuna.create_study(
        directions=["maximize", "minimize"],
        sampler=optuna.samplers.NSGAIISampler(seed=SEED, constraints_func=lambda t: t.user_attrs["constraint"]),
    )
    study.optimize(objective, n_trials=N_TRIALS, n_jobs=N_JOBS)

    df = pd.DataFrame([{**t.params, **t.user_attrs["metrics"], "feasible": t.user_attrs["constraint"][0] <= 0}
                       for t in study.trials if t.values is not None])
    return df.assign(pareto=_is_pareto(df))


def _is_pareto(df: pd.DataFrame) -> np.ndarray:
    """Feasible designs that no other feasible design beats in both yield and buffer."""
    best, on_front = -np.inf, np.zeros(len(df), dtype=bool)
    for i in df[df.feasible].sort_values("buffer_L").index:
        if df.loc[i, "yield"] > best:
            on_front[df.index.get_loc(i)], best = True, df.loc[i, "yield"]
    return on_front


def plot_pareto(df: pd.DataFrame, reference: Dict[str, float], path: Path):
    front = df[df.pareto].sort_values("buffer_L")
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.scatter(df[~df.feasible].buffer_L, df[~df.feasible]["yield"], s=8, c="lightgrey",
               label=f"> {MAX_DENATURED:.0%} denatured")
    pts = ax.scatter(df[df.feasible].buffer_L, df[df.feasible]["yield"], s=14,
                     c=df[df.feasible].pH, cmap="RdYlBu", label="feasible")
    ax.plot(front.buffer_L, front["yield"], "k-", lw=1.5, label="Pareto front")
    ax.plot(reference["buffer_L"], reference["yield"], "*", c="black", ms=16, label="reference recipe")
    fig.colorbar(pts, label="elution buffer pH")
    ax.set(xlabel="elution buffer consumption [L]", ylabel="native hIgG yield [-]",
           title="Elution optimisation: yield vs. buffer consumption")
    ax.legend(loc="lower right", fontsize=8)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(path, dpi=200)


def plot_pH_sensitivity(path: Path, factors=(0.0, 1.0, 3.0)):
    """How far the best pH moves if the assumed denaturation rate is wrong. Uses the reference
    recipe, so only the pH changes."""
    pH_values = np.arange(2.1, 3.81, 0.1)
    fig, (ax_yield, ax_quality) = plt.subplots(1, 2, figsize=(10, 4))
    for factor, style in zip(factors, (":", "-", "--")):
        setup = build_setup({**DENATURATION, "k_ref": factor * DENATURATION["k_ref"]})
        runs = [evaluate(setup, pH=pH, **REFERENCE) for pH in pH_values]
        label = "no denaturation" if factor == 0 else f"{factor:g}x assumed rate"
        ax_yield.plot(pH_values, [r["yield"] for r in runs], style, label=label)
        ax_quality.plot(pH_values, [r["native_fraction"] for r in runs], style, label=label)
        print(f"  {label:20s} best pH {pH_values[int(np.argmax([r['yield'] for r in runs]))]:.1f}, "
              f"yield {max(r['yield'] for r in runs):.3f}")
    ax_yield.set(xlabel="elution buffer pH", ylabel="native hIgG yield [-]", title="Yield")
    ax_quality.axhline(1 - MAX_DENATURED, color="grey", lw=1, ls="-.")
    ax_quality.set(xlabel="elution buffer pH", ylabel="native fraction of eluate [-]", title="Product quality")
    for ax in (ax_yield, ax_quality):
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=200)


if __name__ == "__main__":
    optuna.logging.set_verbosity(optuna.logging.WARNING)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    columns = ["pH", "n_cycles", "t_fill", "t_resuspend_loop", "t_recapture_loop",
               "yield", "buffer_L", "time_s", "native_fraction"]

    setup = build_setup()
    reference = evaluate(setup, pH=setup.reference_pH, **REFERENCE)
    print(f"hIgG on the particles after washing: {setup.hIgG_captured_g:.3f} g")
    print(f"reference recipe (pH {setup.reference_pH:.2f}): "
          + "  ".join(f"{k}={v:.4g}" for k, v in reference.items()))

    df = optimise(setup)
    df.to_csv(OUTPUT_DIR / "trials.csv", index=False)
    plot_pareto(df, reference, OUTPUT_DIR / "pareto_front.png")
    print(f"\n{len(df)} trials, {df.feasible.sum()} within the quality spec, {df.pareto.sum()} on the front")
    print(df[df.pareto].sort_values("buffer_L")[columns].to_string(index=False, float_format=lambda v: f"{v:.3f}"))

    print("\npH sensitivity to the assumed denaturation rate (reference recipe):")
    plot_pH_sensitivity(OUTPUT_DIR / "pH_sensitivity.png")
