from datetime import datetime
from pathlib import Path
from typing import List, Optional

import altair as alt
import numpy as np
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

APP_ROOT = Path(__file__).resolve().parent
REPO_ROOT = APP_ROOT.parent
OUTPUT_ROOT = REPO_ROOT / "data" / "ui_runs"
EXPERIMENTAL_PATH = REPO_ROOT / "data" / "experimental_data.ods"

FLOW_SHEET_LABEL_TO_ENUM = {label: enum for enum, label in FLOW_SHEET_LABELS.items()}
NO_FRACTION = "None"
FRACTION_OPTIONS = [NO_FRACTION] + FRACTIONS

# The time-step scatter has tens of thousands of points; let Altair embed them all.
alt.data_transformers.disable_max_rows()


# =============================================================
# ========================== MODEL ============================
# =============================================================
@st.cache_resource(show_spinner="Building model and equilibrating solutions...")
def get_model():
    cs = build_component_system()
    rs, _ = build_reaction_system(cs, tau_reaction=0.1)
    solutions = build_solutions(cs, rs, fs3.SolverType.ERK)
    return cs, solutions


COMPONENT_SYSTEM, SOLUTIONS = get_model()
SOLUTION_NAMES = list(SOLUTIONS)


# Fixed solver settings for this prototype UI.
SOLVER_TYPE = fs3.SolverType.ERK
TAU_REACTION = 0.1
DISCRETIZATION_FACTOR = 0.1
KF_ION = 1.0e3
TIMEOUT_SECONDS = float("inf")


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
# ========================== CHARTS ===========================
# =============================================================
def ph_chart(obs_dir: Path, transitions: List[float]) -> alt.Chart:
    activity = np.load(obs_dir / "pipe_outlet_middle_cell_pH_activity.npy")
    concentration = np.load(obs_dir / "pipe_outlet_middle_cell_pH_concentration.npy")
    time = np.arange(0.0, len(activity) * 2.0, 2.0)

    frames = [
        pd.DataFrame({"time": time, "pH": activity, "series": "Simulation (activity)"}),
        pd.DataFrame({"time": time, "pH": concentration, "series": "Simulation (concentration)"}),
    ]

    experimental = load_experimental()
    if experimental is not None:
        frames.append(
            pd.DataFrame(
                {"time": experimental["time_s"], "pH": experimental["pH"] - 0.35, "series": "Experimental"}
            )
        )

    df = pd.concat(frames, ignore_index=True)
    lines = (
        alt.Chart(df)
        .mark_line()
        .encode(
            x=alt.X("time:Q", title="Time [s]"),
            y=alt.Y("pH:Q", title="pH [-]", scale=alt.Scale(zero=False)),
            color=alt.Color("series:N", title=None),
            tooltip=[
                alt.Tooltip("series:N", title="Series"),
                alt.Tooltip("time:Q", title="Time [s]", format=".0f"),
                alt.Tooltip("pH:Q", format=".2f"),
            ],
        )
    )
    rules = (
        alt.Chart(pd.DataFrame({"t": transitions}))
        .mark_rule(color="gray", opacity=0.35, strokeDash=[3, 3])
        .encode(x="t:Q")
    )
    return alt.layer(rules, lines).properties(height=420).interactive()


def fractions_chart(obs_dir: Path, component: str = "hIgG") -> Optional[alt.Chart]:
    npz_path = obs_dir / "unit_operations.npz"
    if not npz_path.exists():
        return None
    data = np.load(npz_path)
    component_idx = COMPONENT_SYSTEM.get_idx(component)

    rows = [
        {"scenario": "Simulation", "fraction": name, "mass": float(data[name][-1, 0, component_idx]) * 1000.0}
        for name in FRACTIONS
        if name in data
    ]
    experiment = {"Elution 1": 0.32, "Elution 2": 0.63, "Elution 3": 0.42, "Elution 4": 0.32, "Elution 5": 0.06}
    rows += [{"scenario": "Experiment", "fraction": name, "mass": mass} for name, mass in experiment.items()]

    df = pd.DataFrame(rows)
    return (
        alt.Chart(df)
        .mark_bar()
        .encode(
            x=alt.X("scenario:N", title=None),
            y=alt.Y("mass:Q", title=f"Mass {component} [g]"),
            color=alt.Color("fraction:N", title="Fraction", sort=FRACTIONS),
            order=alt.Order("fraction:N", sort="ascending"),
            tooltip=[
                alt.Tooltip("scenario:N", title="Scenario"),
                alt.Tooltip("fraction:N", title="Fraction"),
                alt.Tooltip("mass:Q", title="Mass [g]", format=".3f"),
            ],
        )
        .properties(height=420)
    )


def time_step_chart(obs_dir: Path) -> alt.Chart:
    timestamps = np.load(obs_dir / "internal_timestamps.npy")
    df = pd.DataFrame({"time": timestamps[:-1], "dt": np.diff(timestamps)})
    return (
        alt.Chart(df)
        .mark_circle(size=8, opacity=0.5)
        .encode(
            x=alt.X("time:Q", title="Time [s]"),
            y=alt.Y("dt:Q", title="Step size [s]"),
            tooltip=[alt.Tooltip("time:Q", title="Time [s]", format=".2f"), alt.Tooltip("dt:Q", title="Step size [s]", format=".4f")],
        )
        .properties(height=380)
        .interactive()
    )


@st.cache_data(show_spinner=False)
def load_experimental() -> Optional[pd.DataFrame]:
    if not EXPERIMENTAL_PATH.exists():
        return None
    try:
        df = pd.read_excel(EXPERIMENTAL_PATH, engine="odf")
    except Exception:
        return None
    return df[["time_s", "pH"]].sort_values("time_s")


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
        "Flow_Config": st.column_config.SelectboxColumn("Flow config", options=list(FLOW_SHEET_LABEL_TO_ENUM)),
        "Inlet": st.column_config.SelectboxColumn("Inlet solution", options=SOLUTION_NAMES),
        "Fraction": st.column_config.SelectboxColumn("Fraction", options=FRACTION_OPTIONS),
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
    run_tag = run_name.strip() or None
    if run_tag is not None:
        run_tag = f"{run_tag}_{datetime.now().strftime('%Y-%m-%d_%H-%M-%S')}"
    with st.spinner("Running simulation..."):
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
        )
    st.session_state.last_run = {
        "run_dir": run_dir,
        "transitions": phase_transitions(recipe_steps),
    }
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
        st.altair_chart(ph_chart(obs_dir, last_run["transitions"]), width="stretch")

        st.markdown("##### Fraction masses")
        frac = fractions_chart(obs_dir)
        if frac is None:
            st.warning("Fraction data not found for this run.")
        else:
            st.altair_chart(frac, width="stretch")

        st.markdown("##### Solver time-step sizes")
        st.altair_chart(time_step_chart(obs_dir), width="stretch")
