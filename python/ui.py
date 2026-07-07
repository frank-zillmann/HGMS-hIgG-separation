import threading
from datetime import datetime
from pathlib import Path
from typing import List

import pandas as pd
import streamlit as st

import fs3
from components import build_component_system
from experiment import Experiment, build_experiment
from observers import save_observations
from plot_fractions import build_fractions_figure, default_datasets, load_reference_fractions
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

OUTPUT_ROOT = Path(__file__).resolve().parent.parent / "data" / "ui_runs"

FLOW_SHEET_LABEL_TO_ENUM = {label: enum for enum, label in FLOW_SHEET_LABELS.items()}
FRACTION_OPTIONS = ["None"] + FRACTIONS
SOLVERS = {"ERK": fs3.SolverType.ERK, "ARK": fs3.SolverType.ARK, "ADAMS": fs3.SolverType.ADAMS, "BDF": fs3.SolverType.BDF}


@st.cache_resource(show_spinner="Equilibrating solutions...")
def get_solutions(tau_reaction: float, solver_name: str):
    cs = build_component_system()
    rs, _ = build_reaction_system(cs, tau_reaction)
    return build_solutions(cs, rs, SOLVERS[solver_name])


get_experimental = st.cache_data(show_spinner=False)(load_experimental)
get_references = st.cache_data(show_spinner=False)(load_reference_fractions)


# Solution names are independent of tau/solver, so build once for the recipe dropdown.
SOLUTION_NAMES = list(get_solutions(0.1, "ERK"))


def _solve(experiment: Experiment, holder: dict):
    """Run the (blocking) solve; store any error so the UI thread can surface it afterwards."""
    try:
        experiment.solve()
    except BaseException as exc:
        holder["error"] = exc


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
                "Fraction": step.fraction or "None",
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
            fraction=None if row["Fraction"] == "None" else row["Fraction"],
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
    # Streamlit runs server-side, so it can't open a native OS "Save As" dialog; instead we let
    # the user pick a folder (server filesystem) and a name. Keys keep edits stable across reruns.
    folder = st.text_input("Folder", key="save_run_dir")
    name = st.text_input("Run name", key="save_run_name")
    target = Path(folder) / name.strip() / "obs" if name.strip() else None
    st.caption(f"Saves to `{target}`" if target else "Enter a run name.")
    if st.button("Save", type="primary", disabled=target is None):
        st.session_state.saved_path = str(save_observations(experiment, target))
        st.rerun()


# =============================================================
# ========================== LAYOUT ===========================
# =============================================================
st.set_page_config(page_title="HGMS hIgG Simulation Studio", page_icon="🧪", layout="wide")
st.title("HGMS hIgG Simulation Studio")
st.caption("Shape the recipe, launch a run, and inspect the core figures in one place.")

st.session_state.setdefault("editor_version", 0)
st.session_state.setdefault("experiment", None)  # last finished (or currently running) experiment
st.session_state.setdefault("running", False)

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
    discretization_factor = st.number_input("Discretization factor", min_value=0.0, value=0.5, format="%g", step=0.1)
with opt2:
    tau_reaction = st.number_input("τ reaction [s]", min_value=0.001, value=0.1, format="%g", step=0.01)
with opt3:
    solver_name = st.selectbox("Solver type", options=list(SOLVERS), index=0)

run_col, abort_col = st.columns(2)
run_clicked = run_col.button("Run simulation", type="primary", disabled=bool(errors) or st.session_state.running, width="stretch")
abort_clicked = abort_col.button("Abort simulation", disabled=not st.session_state.running, width="stretch")

if run_clicked:
    recipe_steps = dataframe_to_recipe(edited_df)
    experiment = build_experiment(
        recipe_steps,
        discretization_factor=discretization_factor,
        tau_reaction=tau_reaction,
        solver_type=SOLVERS[solver_name],
        solutions=get_solutions(tau_reaction, solver_name),
    )
    holder: dict = {}
    worker = threading.Thread(target=_solve, args=(experiment, holder), daemon=True)
    worker.start()
    st.session_state.update(
        experiment=experiment, worker=worker, holder=holder,
        transitions=phase_transitions(recipe_steps), running=True,
    )
    st.rerun()

if abort_clicked:
    # The C++ solve can't be cancelled, so the background thread finishes on its own; we
    # simply stop waiting on it and discard the (partial) experiment.
    st.session_state.update(running=False, experiment=None, worker=None)
    st.rerun()

# ----- Live (only while a run is in progress) ----------------------------------
if st.session_state.running:
    @st.fragment(run_every=1.0)
    def live_view():
        experiment = st.session_state.experiment
        span = (0.0, experiment.total_duration)
        times, ph_activity, ph_concentration = experiment.live_pH()
        st.plotly_chart(
            build_ph_figure(times, ph_activity, ph_concentration, experimental=get_experimental(), transitions=st.session_state.transitions, x_range=span, title="pH (live)"),
            width="stretch", key="live_ph",
        )
        st.plotly_chart(
            build_time_step_figure(*experiment.step_sizes(), x_range=span, title="Solver time steps (live)"),
            width="stretch", key="live_ts",
        )
        st.progress(experiment.progress, text=f"Simulating… {experiment.t:.0f} / {experiment.total_duration:.0f} s")
        if not st.session_state.worker.is_alive():
            st.session_state.running = False
            if st.session_state.holder.get("error"):
                st.session_state.update(experiment=None, run_error=str(st.session_state.holder["error"]))
            st.rerun(scope="app")

    live_view()

# ----- Results (hidden while a run is streaming above) -------------------------
if not st.session_state.running:
    st.subheader("Results")
    if st.session_state.get("run_error"):
        st.error(f"Simulation failed: {st.session_state.run_error}")
        st.session_state.run_error = None

    experiment = st.session_state.experiment
    experimental = get_experimental()
    st.markdown("##### pH profile")
    if experiment is None:
        # No run yet (or aborted): still show the experimental pH so it is always comparable.
        x_range = (0.0, float(experimental["time_s"].max())) if experimental is not None else None
        st.plotly_chart(build_ph_figure([], [], experimental=experimental, x_range=x_range), width="stretch")
        st.info("Run a simulation to overlay the model results.")
    else:
        span = (0.0, experiment.total_duration)
        st.plotly_chart(
            build_ph_figure(
                experiment.times(), experiment.outlet_pH(True), experiment.outlet_pH(False),
                experimental=experimental, transitions=st.session_state.transitions, x_range=span,
            ),
            width="stretch",
        )
        st.markdown("##### Fraction masses")
        st.plotly_chart(
            build_fractions_figure(default_datasets(experiment.fraction_masses(), get_references())),
            width="stretch",
        )
        st.markdown("##### Solver time-step sizes")
        st.plotly_chart(build_time_step_figure(*experiment.step_sizes(), x_range=span), width="stretch")

        st.divider()
        if st.button("Save run…", type="primary"):
            st.session_state.save_run_dir = str(OUTPUT_ROOT)
            st.session_state.save_run_name = datetime.now().strftime("run_%Y-%m-%d_%H-%M-%S")
            save_run_dialog(experiment)
