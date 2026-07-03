import threading
from datetime import datetime
from pathlib import Path
from typing import Callable, List

import pandas as pd
import streamlit as st

import fs3
from components import build_component_system
from experiment import Simulation, build_simulation
from figures import result_figures
from observers import save_observations
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

# Fixed solver settings for this prototype UI.
SOLVER_TYPE = fs3.SolverType.ERK
TAU_REACTION = 0.1
DISCRETIZATION_FACTOR = 0.1


@st.cache_resource(show_spinner="Building model and equilibrating solutions...")
def get_solutions():
    cs = build_component_system()
    rs, _ = build_reaction_system(cs, TAU_REACTION)
    return build_solutions(cs, rs, SOLVER_TYPE)


@st.cache_data(show_spinner=False)
def get_experimental():
    return load_experimental(str(EXPERIMENTAL_PATH))


SOLUTIONS = get_solutions()
SOLUTION_NAMES = list(SOLUTIONS)


def solve_streaming(sim: Simulation, on_update: Callable[[Simulation], None], interval: float = 1.0) -> None:
    """Solve in a background thread; call ``on_update(sim)`` every ``interval`` seconds and once at the end."""
    done = threading.Event()
    error: list = []

    def worker():
        try:
            sim.solve()
        except BaseException as exc:  # re-raised on the main thread below
            error.append(exc)
        finally:
            done.set()

    thread = threading.Thread(target=worker, daemon=True)
    thread.start()
    while not done.wait(interval):
        on_update(sim)
    thread.join()
    on_update(sim)
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


# =============================================================
# ========================== LAYOUT ===========================
# =============================================================
st.set_page_config(page_title="HGMS hIgG Simulation Studio", page_icon="🧪", layout="wide")
st.title("HGMS hIgG Simulation Studio")
st.caption("Shape the recipe, launch a run, and inspect the core figures in one place.")

if "editor_version" not in st.session_state:
    st.session_state.editor_version = 0
if "results" not in st.session_state:
    st.session_state.results = None

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
run_col_left, run_col_right = st.columns([3, 1])
with run_col_left:
    run_name = st.text_input("Run name", placeholder="e.g. baseline_recipe")
    save_to_disk = st.checkbox("Save observations to disk", value=False)
with run_col_right:
    st.markdown("<div style='height: 1.75rem'></div>", unsafe_allow_html=True)
    run_button = st.button("Run simulation", type="primary", disabled=bool(errors), width="stretch")

if run_button:
    recipe_steps = dataframe_to_recipe(edited_df)
    transitions = phase_transitions(recipe_steps)
    sim = build_simulation(recipe_steps, discretization_factor=DISCRETIZATION_FACTOR, solutions=SOLUTIONS)

    st.markdown("##### Live")
    ph_slot = st.empty()
    step_slot = st.empty()
    progress = st.progress(0.0, text="Starting simulation…")
    frame = {"i": 0}

    def on_update(sim: Simulation):
        frame["i"] += 1
        times, ph = sim.live_pH()
        ph_slot.plotly_chart(
            build_ph_figure(times, ph, transitions=transitions, x_range=(0.0, sim.total_duration), title="pH (live)"),
            width="stretch", key=f"live_ph_{frame['i']}",
        )
        step_times, step_sizes = sim.step_sizes()
        step_slot.plotly_chart(
            build_time_step_figure(step_times, step_sizes, title="Solver time steps (live)"),
            width="stretch", key=f"live_ts_{frame['i']}",
        )
        progress.progress(sim.progress, text=f"Simulating… {sim.t:.0f} / {sim.total_duration:.0f} s")

    solve_streaming(sim, on_update, interval=1.0)
    progress.progress(1.0, text="Simulation finished.")
    ph_slot.empty()
    step_slot.empty()

    if save_to_disk:
        tag = (run_name.strip() or "run") + datetime.now().strftime("_%Y-%m-%d_%H-%M-%S")
        save_observations(sim, OUTPUT_ROOT / tag / "obs")

    st.session_state.results = sim.results()
    st.success("Run finished.")

# ----- Results -----------------------------------------------------------------
st.subheader("Results")
results = st.session_state.results
if results is None:
    st.info("Run a simulation to view results.")
else:
    figures = result_figures(results, get_experimental())
    st.markdown("##### pH profile")
    st.plotly_chart(figures["pH"], width="stretch")
    st.markdown("##### Fraction masses")
    st.plotly_chart(figures["fractions"], width="stretch")
    st.markdown("##### Solver time-step sizes")
    st.plotly_chart(figures["time_steps"], width="stretch")
