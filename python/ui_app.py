from datetime import datetime
from pathlib import Path
from typing import List, Optional

import pandas as pd
import streamlit as st

import fs3
from HGMS_hIgG_separation import run_HGMS_hIgG_separation
from model import (
    FLOW_SHEET_LABELS,
    FRACTIONS,
    RecipeStep,
    build_component_system,
    build_reaction_system,
    build_solutions,
    default_recipe,
    total_duration,
)
from plot_fractions import plot_fractions
from plot_pH import build_ph_figure, plot_pH
from plot_time_step_sizes import plot_time_step_sizes

APP_ROOT = Path(__file__).resolve().parent
REPO_ROOT = APP_ROOT.parent
OUTPUT_ROOT = REPO_ROOT / "data" / "ui_runs"
EXPERIMENTAL_PATH = REPO_ROOT / "data" / "experimental_data.ods"

FLOW_SHEET_LABEL_TO_ENUM = {label: enum for enum, label in FLOW_SHEET_LABELS.items()}
NO_FRACTION = "None"
FRACTION_OPTIONS = [NO_FRACTION] + FRACTIONS

# Fixed solver settings for this prototype UI.
SOLVER_TYPE = fs3.SolverType.ERK
TAU_REACTION = 0.1
DISCRETIZATION_FACTOR = 0.1
KF_ION = 1.0e3
TIMEOUT_SECONDS = float("inf")


@st.cache_resource(show_spinner="Building model and equilibrating solutions...")
def get_solution_names() -> List[str]:
    cs = build_component_system()
    rs, _ = build_reaction_system(cs, tau_reaction=TAU_REACTION)
    return list(build_solutions(cs, rs, SOLVER_TYPE))


SOLUTION_NAMES = get_solution_names()


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


def phase_transitions(recipe_steps: List[RecipeStep]) -> List[float]:
    """Cumulative end time of every step -> vertical guide lines matching the recipe."""
    ends, total = [], 0.0
    for step in recipe_steps:
        total += step.t_duration
        ends.append(total)
    return ends


# =============================================================
# ========================== LAYOUT ===========================
# =============================================================
st.set_page_config(page_title="HGMS hIgG Simulation Studio", page_icon="🧪", layout="wide")
st.title("HGMS hIgG Simulation Studio")
st.caption("Shape the recipe, launch a run, and inspect the core figures in one place.")

if "editor_version" not in st.session_state:
    st.session_state.editor_version = 0
if "last_run" not in st.session_state:
    st.session_state.last_run = None

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
with run_col_right:
    st.markdown("<div style='height: 1.75rem'></div>", unsafe_allow_html=True)
    run_button = st.button("Run simulation", type="primary", disabled=bool(errors), width="stretch")

if run_button:
    recipe_steps = dataframe_to_recipe(edited_df)
    total = total_duration(recipe_steps)
    transitions = phase_transitions(recipe_steps)
    run_tag = run_name.strip() or None
    if run_tag is not None:
        run_tag = f"{run_tag}_{datetime.now().strftime('%Y-%m-%d_%H-%M-%S')}"

    st.markdown("##### Live pH")
    live_chart = st.empty()
    progress = st.progress(0.0, text="Starting simulation…")
    frame = {"i": 0}

    def on_progress(times, ph, t_now, total_time):
        frame["i"] += 1
        live_chart.plotly_chart(
            build_ph_figure(times, ph, transitions=transitions, x_range=(0.0, total_time), title="pH (live)"),
            width="stretch",
            key=f"live_ph_{frame['i']}",
        )
        progress.progress(min(t_now / total_time, 1.0), text=f"Simulating… {t_now:.0f} / {total_time:.0f} s")

    result, run_dir, obs_dir = run_HGMS_hIgG_separation(
        kf_ion=KF_ION,
        tau_reaction=TAU_REACTION,
        save_obs=True,
        timeout_seconds=TIMEOUT_SECONDS,
        solver_type=SOLVER_TYPE,
        discretization_factor=DISCRETIZATION_FACTOR,
        recipe_steps=recipe_steps,
        output_base_dir=OUTPUT_ROOT,
        run_tag=run_tag,
        return_paths=True,
        progress_callback=on_progress,
    )
    progress.progress(1.0, text="Simulation finished.")
    live_chart.empty()  # the Results section below renders the full figure
    st.session_state.last_run = {"run_dir": str(run_dir), "transitions": transitions}
    st.success("Run finished.")

# ----- Results -----------------------------------------------------------------
st.subheader("Results")
last_run = st.session_state.last_run
if last_run is None:
    st.info("Run a simulation to view results.")
else:
    obs_dir = Path(last_run["run_dir"]) / "obs"
    if not obs_dir.exists():
        st.warning("Observation data not found for this run.")
    else:
        st.markdown("##### pH profile")
        fig, _ = plot_pH(
            path_to_obs=str(obs_dir),
            experimental_path=str(EXPERIMENTAL_PATH),
            phase_transitions=last_run["transitions"],
            save_plots=False,
        )
        st.plotly_chart(fig, width="stretch")

        st.markdown("##### Fraction masses")
        fig, _ = plot_fractions(
            path_to_obs_dict={Path(last_run["run_dir"]).name: str(obs_dir)},
            save_plots=False,
        )
        st.plotly_chart(fig, width="stretch")

        st.markdown("##### Solver time-step sizes")
        fig, _ = plot_time_step_sizes(path_to_obs=str(obs_dir), save_plots=False)
        st.plotly_chart(fig, width="stretch")
