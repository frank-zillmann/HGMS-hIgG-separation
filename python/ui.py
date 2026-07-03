import threading
from datetime import datetime
from pathlib import Path
from typing import Callable, List

import pandas as pd
import streamlit as st

import fs3
from components import build_component_system
from experiment import Experiment, build_experiment
from observers import save_observations
from plot_fractions import EXPERIMENT, build_fractions_figure
from plot_pH import build_ph_figure, load_experimental
from plot_time_step_sizes import build_time_step_figure
from reactions import build_reaction_system
from recipe import (
    FLOW_SHEET_LABELS,
    FRACTIONS,
    RecipeStep,
    default_recipe,
    phase_transitions,
    total_duration,
)
from solutions import build_solutions

REPO_ROOT = Path(__file__).resolve().parent.parent
OUTPUT_ROOT = REPO_ROOT / "data" / "ui_runs"
EXPERIMENTAL_PATH = REPO_ROOT / "data" / "experimental_data.ods"

FLOW_SHEET_LABEL_TO_ENUM = {label: enum for enum, label in FLOW_SHEET_LABELS.items()}
NO_FRACTION = "None"
FRACTION_OPTIONS = [NO_FRACTION] + FRACTIONS
SOLVERS = {"ERK": fs3.SolverType.ERK, "ARK": fs3.SolverType.ARK, "ADAMS": fs3.SolverType.ADAMS, "BDF": fs3.SolverType.BDF}


@st.cache_resource(show_spinner="Equilibrating solutions...")
def get_solutions(tau_reaction: float, solver_name: str):
    cs = build_component_system()
    rs, _ = build_reaction_system(cs, tau_reaction)
    return build_solutions(cs, rs, SOLVERS[solver_name])


@st.cache_data(show_spinner=False)
def get_experimental():
    return load_experimental(str(EXPERIMENTAL_PATH))


# Solution names are independent of tau/solver, so build once for the recipe dropdown.
SOLUTION_NAMES = list(get_solutions(0.1, "ERK"))


def solve_streaming(experiment: Experiment, on_update: Callable[[Experiment], None], interval: float = 1.0) -> None:
    """Solve in a background thread; call ``on_update(experiment)`` every ``interval`` s and once at the end."""
    done = threading.Event()
    error: list = []

    def worker():
        try:
            experiment.solve()
        except BaseException as exc:  # re-raised on the main thread below
            error.append(exc)
        finally:
            done.set()

    thread = threading.Thread(target=worker, daemon=True)
    thread.start()
    while not done.wait(interval):
        on_update(experiment)
    thread.join()
    on_update(experiment)
    if error:
        raise error[0]


# =============================================================
# ========================== RECIPE <-> DF ====================
# =============================================================
def recipe_to_dataframe(recipe_steps: List[RecipeStep]) -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "Name": step.name,
                "Duration_s": step.t_duration,
                "Pump_%": step.pump_percentage,
                "Mixing_%": step.mixing_percentage,
                "Flow_Config": FLOW_SHEET_LABELS[step.flow_sheet_configuration],
                "Inlet": step.inlet,
                "Fraction": step.fraction or NO_FRACTION,
            }
            for step in recipe_steps
        ]
    )


def dataframe_to_recipe(df: pd.DataFrame) -> List[RecipeStep]:
    return [
        RecipeStep(
            name=str(row["Name"]).strip(),
            t_duration=float(row["Duration_s"]),
            pump_percentage=float(row["Pump_%"]),
            mixing_percentage=float(row["Mixing_%"]),
            flow_sheet_configuration=FLOW_SHEET_LABEL_TO_ENUM[row["Flow_Config"]],
            inlet=row["Inlet"],
            fraction=None if row["Fraction"] == NO_FRACTION else row["Fraction"],
        )
        for row in df.to_dict(orient="records")
    ]


def validate_recipe(df: pd.DataFrame) -> List[str]:
    errors = []
    durations = pd.to_numeric(df["Duration_s"], errors="coerce")
    pump = pd.to_numeric(df["Pump_%"], errors="coerce")
    mixing = pd.to_numeric(df["Mixing_%"], errors="coerce")

    if durations.isna().any() or (durations <= 0).any():
        errors.append("All durations must be positive numbers.")
    if pump.isna().any() or (pump < 0).any() or (pump > 100).any():
        errors.append("Pump percentage must be within 0-100.")
    if mixing.isna().any() or (mixing < 0).any() or (mixing > 100).any():
        errors.append("Mixing percentage must be within 0-100.")
    if df["Name"].isna().any() or (df["Name"].astype(str).str.strip() == "").any():
        errors.append("Every step needs a name.")
    if not df["Inlet"].isin(SOLUTION_NAMES).all():
        errors.append(f"Inlet must be one of: {', '.join(SOLUTION_NAMES)}.")
    if not df["Fraction"].isin(FRACTION_OPTIONS).all():
        errors.append("Fraction contains an unknown value.")
    return errors


@st.dialog("Save run")
def save_run_dialog(experiment: Experiment):
    name = st.text_input("Run name", value=datetime.now().strftime("run_%Y-%m-%d_%H-%M-%S"))
    if st.button("Save", type="primary", disabled=not name.strip()):
        path = save_observations(experiment, OUTPUT_ROOT / name.strip() / "obs")
        st.session_state.saved_path = str(path)
        st.rerun()


# =============================================================
# ========================== LAYOUT ===========================
# =============================================================
st.set_page_config(page_title="HGMS hIgG Simulation Studio", page_icon="🧪", layout="wide")
st.title("HGMS hIgG Simulation Studio")
st.caption("Shape the recipe, launch a run, and inspect the core figures in one place.")

if "editor_version" not in st.session_state:
    st.session_state.editor_version = 0
if "experiment" not in st.session_state:
    st.session_state.experiment = None

if st.session_state.get("saved_path"):
    st.success(f"Saved run to {st.session_state.saved_path}")
    st.session_state.saved_path = None

# ----- Recipe editor -----------------------------------------------------------
st.subheader("Recipe editor")
edited_df = st.data_editor(
    recipe_to_dataframe(default_recipe()),
    key=f"recipe_editor_{st.session_state.editor_version}",
    num_rows="dynamic",
    width="stretch",
    column_config={
        "Duration_s": st.column_config.NumberColumn("Duration [s]", min_value=0.0, step=1.0),
        "Pump_%": st.column_config.NumberColumn("Pump [%]", min_value=0.0, max_value=100.0, step=1.0),
        "Mixing_%": st.column_config.NumberColumn("Mixing [%]", min_value=0.0, max_value=100.0, step=1.0),
        "Flow_Config": st.column_config.SelectboxColumn("Flow config", options=list(FLOW_SHEET_LABEL_TO_ENUM), required=True),
        "Inlet": st.column_config.SelectboxColumn("Inlet solution", options=SOLUTION_NAMES, required=True),
        "Fraction": st.column_config.SelectboxColumn("Fraction", options=FRACTION_OPTIONS, required=True),
    },
)

errors = validate_recipe(edited_df)
if errors:
    st.error(" ".join(errors))
    st.markdown("**Total duration:** –")
else:
    st.markdown(f"**Total duration:** {total_duration(dataframe_to_recipe(edited_df)):.1f} s")

if st.button("Reset to default recipe"):
    st.session_state.editor_version += 1
    st.rerun()

# ----- Run ---------------------------------------------------------------------
st.subheader("Run")
opt1, opt2, opt3 = st.columns(3)
with opt1:
    discretization_factor = st.number_input("Discretization factor", min_value=0.01, max_value=2.0, value=0.1, step=0.05)
with opt2:
    tau_reaction = st.number_input("τ reaction [s]", min_value=0.0, max_value=100.0, value=0.1, step=0.1, format="%.2f")
with opt3:
    solver_name = st.selectbox("Solver type", options=list(SOLVERS), index=0)

run_button = st.button("Run simulation", type="primary", disabled=bool(errors), width="stretch")

if run_button:
    recipe_steps = dataframe_to_recipe(edited_df)
    transitions = phase_transitions(recipe_steps)
    experiment = build_experiment(
        recipe_steps,
        discretization_factor=discretization_factor,
        tau_reaction=tau_reaction,
        solver_type=SOLVERS[solver_name],
        solutions=get_solutions(tau_reaction, solver_name),
    )
    span = (0.0, experiment.total_duration)

    st.markdown("##### Live")
    ph_slot = st.empty()
    step_slot = st.empty()
    progress = st.progress(0.0, text="Starting simulation…")
    frame = {"i": 0}

    def on_update(experiment: Experiment):
        frame["i"] += 1
        times, ph = experiment.live_pH()
        ph_slot.plotly_chart(
            build_ph_figure(times, ph, transitions=transitions, x_range=span, title="pH (live)"),
            width="stretch", key=f"live_ph_{frame['i']}",
        )
        step_slot.plotly_chart(
            build_time_step_figure(*experiment.step_sizes(), x_range=span, title="Solver time steps (live)"),
            width="stretch", key=f"live_ts_{frame['i']}",
        )
        progress.progress(experiment.progress, text=f"Simulating… {experiment.t:.0f} / {experiment.total_duration:.0f} s")

    solve_streaming(experiment, on_update, interval=1.0)
    progress.progress(1.0, text="Simulation finished.")
    ph_slot.empty()
    step_slot.empty()

    st.session_state.experiment = experiment
    st.session_state.transitions = transitions
    st.success("Run finished.")

# ----- Results -----------------------------------------------------------------
st.subheader("Results")
experiment = st.session_state.experiment
if experiment is None:
    st.info("Run a simulation to view results.")
else:
    span = (0.0, experiment.total_duration)
    st.markdown("##### pH profile")
    st.plotly_chart(
        build_ph_figure(
            experiment.times(), experiment.outlet_pH(True), experiment.outlet_pH(False),
            experimental=get_experimental(), transitions=st.session_state.transitions, x_range=span,
        ),
        width="stretch",
    )
    st.markdown("##### Fraction masses")
    st.plotly_chart(
        build_fractions_figure({"Simulation": experiment.fraction_masses(), "Experiment": EXPERIMENT}),
        width="stretch",
    )
    st.markdown("##### Solver time-step sizes")
    st.plotly_chart(build_time_step_figure(*experiment.step_sizes(), x_range=span), width="stretch")

    st.divider()
    if st.button("Save run…", type="primary"):
        save_run_dialog(experiment)
