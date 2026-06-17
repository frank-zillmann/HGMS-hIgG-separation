import math
from datetime import datetime
from pathlib import Path
from typing import List, Optional

import pandas as pd
import streamlit as st

import fs3
from HGMS_hIgG_separation import run_HGMS_hIgG_separation
from model import build_component_system, build_reaction_system, build_solutions
from plot_fractions import plot_fractions
from plot_pH import plot_pH
from plot_time_step_sizes import plot_time_step_sizes
from shared_data import (
    FLOW_SHEET_LABELS,
    FRACTIONS,
    RecipeStep,
    default_recipe,
    total_duration,
)

APP_ROOT = Path(__file__).resolve().parent
REPO_ROOT = APP_ROOT.parent
OUTPUT_ROOT = REPO_ROOT / "data" / "ui_runs"

FLOW_SHEET_LABEL_TO_ENUM = {label: enum for enum, label in FLOW_SHEET_LABELS.items()}
NO_FRACTION = "None"


@st.cache_resource
def get_solutions():
    cs = build_component_system()
    rs, _ = build_reaction_system(cs, tau_reaction=0.1)
    return build_solutions(cs, rs, fs3.SolverType.ERK)


SOLUTIONS = get_solutions()


st.set_page_config(
    page_title="HGMS hIgG Simulation Studio",
    page_icon="HG",
    layout="wide",
)

st.markdown(
    """
<style>
@import url('https://fonts.googleapis.com/css2?family=Space+Grotesk:wght@400;600;700&family=IBM+Plex+Mono:wght@400;600&display=swap');

:root {
  --ink: #0b1a1a;
  --muted: #4b5b5b;
  --accent: #0f8a6a;
  --accent-2: #f0b429;
  --panel: rgba(255, 255, 255, 0.85);
  --panel-border: rgba(13, 33, 33, 0.12);
  --bg: radial-gradient(circle at 15% 20%, #e6f7f2 0%, #f7f5ea 45%, #f2efe7 100%);
}

html, body, [class*="stApp"] {
  background: var(--bg);
  color: var(--ink);
  font-family: 'Space Grotesk', sans-serif;
}

section[data-testid="stSidebar"], div[data-testid="stSidebar"] {
    display: none;
}

h1, h2, h3 {
  letter-spacing: -0.02em;
}

.small-note {
  font-family: 'IBM Plex Mono', monospace;
  color: var(--muted);
  font-size: 0.85rem;
}

.stButton>button {
  background: var(--accent);
  color: white;
  border-radius: 999px;
  padding: 0.6rem 1.4rem;
  border: none;
  box-shadow: 0 10px 20px rgba(15, 138, 106, 0.2);
  transition: transform 0.15s ease, box-shadow 0.15s ease;
}

.stButton>button:hover {
  transform: translateY(-1px);
  box-shadow: 0 12px 26px rgba(15, 138, 106, 0.25);
}

.panel {
  background: var(--panel);
  border: 1px solid var(--panel-border);
  border-radius: 18px;
    margin: 1rem 0;
  padding: 1.2rem 1.4rem;
  box-shadow: 0 18px 40px rgba(15, 24, 24, 0.08);
}

@keyframes fadeInUp {
  from { opacity: 0; transform: translateY(8px); }
  to { opacity: 1; transform: translateY(0); }
}

.fade-in {
  animation: fadeInUp 0.6s ease both;
}
</style>
""",
    unsafe_allow_html=True,
)

st.markdown(
        """
<div class="panel fade-in">
    <h1>HGMS hIgG Simulation Studio</h1>
    <p class="small-note">Shape the recipe, launch a run, and inspect the core figures in one place.</p>
</div>
""",
        unsafe_allow_html=True,
)


def recipe_to_dataframe(recipe_steps: List[RecipeStep]) -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "Name": step.name,
                "Duration_s": step.t_duration,
                "Pump_%": step.pump_percentage,
                "Mixing_%": step.mixing_percentage,
                "Flow_Config": FLOW_SHEET_LABELS[step.flow_sheet_configuration],
                "Inlet": step.inlet.name,
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
            inlet=SOLUTIONS[row["Inlet"]],
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
    return errors


def recipe_phase_transitions(recipe_steps: List[RecipeStep]) -> List[float]:
    transitions = [0.0]
    total = 0.0
    for step in recipe_steps:
        total += step.t_duration
        transitions.append(total)
    return transitions


solver_type = fs3.SolverType.ERK
tau_reaction = 0.1
discretization_factor = 0.1
kf_ion = 1.0e3
timeout_seconds = float("inf")
save_obs = True
save_plots = True


if "recipe_df" not in st.session_state:
    st.session_state.recipe_df = recipe_to_dataframe(default_recipe(SOLUTIONS))

if "last_run" not in st.session_state:
    st.session_state.last_run = None

st.markdown("\n")

st.markdown("<div class=\"panel\"><h2>Recipe Editor</h2></div>", unsafe_allow_html=True)
edited_df = st.data_editor(
    st.session_state.recipe_df,
    num_rows="fixed",
    use_container_width=True,
    column_config={
        "Flow_Config": st.column_config.SelectboxColumn(
            "Flow Config",
            options=list(FLOW_SHEET_LABEL_TO_ENUM.keys()),
        ),
        "Inlet": st.column_config.SelectboxColumn(
            "Inlet Solution",
            options=list(SOLUTIONS),
        ),
        "Fraction": st.column_config.SelectboxColumn(
            "Fraction",
            options=[NO_FRACTION] + FRACTIONS,
        ),
    },
)
st.session_state.recipe_df = edited_df

errors = validate_recipe(edited_df)
if errors:
    st.error(" ".join(errors))

total_time = total_duration(dataframe_to_recipe(edited_df)) if not errors else math.nan
st.markdown(f"**Total duration:** {total_time:.1f} s")

if st.button("Reset to default recipe"):
    st.session_state.recipe_df = recipe_to_dataframe(default_recipe(SOLUTIONS))
    st.experimental_rerun()

st.markdown("\n")

st.markdown("<div class=\"panel\"><h2>Run</h2></div>", unsafe_allow_html=True)
run_col_left, run_col_right = st.columns([3, 1])
with run_col_left:
    run_name = st.text_input("Run name", placeholder="e.g. baseline_recipe")
with run_col_right:
    run_button = st.button("Run simulation")

if run_button:
    if errors:
        st.error("Fix the recipe issues before running.")
    else:
        recipe_steps = dataframe_to_recipe(edited_df)
        run_tag = run_name.strip() or None
        if run_tag is not None:
            run_tag = f"{run_tag}_{datetime.now().strftime('%Y-%m-%d_%H-%M-%S')}"
        with st.spinner("Running simulation..."):
            result, run_dir, obs_dir = run_HGMS_hIgG_separation(
                kf_ion=kf_ion,
                tau_reaction=tau_reaction,
                save_obs=save_obs,
                timeout_seconds=timeout_seconds,
                solver_type=solver_type,
                discretization_factor=discretization_factor,
                recipe_steps=recipe_steps,
                output_base_dir=OUTPUT_ROOT,
                run_tag=run_tag,
                return_paths=True,
            )
        st.session_state.last_run = {
            "result": result,
            "run_dir": run_dir,
            "obs_dir": obs_dir,
        }
        st.success("Run finished.")


st.markdown("\n")

st.markdown("\n")

st.markdown("<div class=\"panel\"><h2>Results</h2></div>", unsafe_allow_html=True)

selected_run_dir: Optional[Path] = None
if st.session_state.last_run is not None:
    selected_run_dir = st.session_state.last_run["run_dir"]

if selected_run_dir is None:
    st.info("Run a simulation to view results.")
else:
    obs_dir = selected_run_dir / "obs"
    plots_dir = selected_run_dir / "plots"

    if obs_dir.exists():
        phase_transitions = None
        if not errors:
            phase_transitions = recipe_phase_transitions(dataframe_to_recipe(edited_df))

        st.markdown("<div class=\"panel\"><h3>pH Profile</h3></div>", unsafe_allow_html=True)
        fig, _ = plot_pH(
            path_to_obs=str(obs_dir),
            path_to_save=str(plots_dir),
            include_experimental=True,
            phase_transitions=phase_transitions,
            save_plots=save_plots,
            close_fig=False,
        )
        st.pyplot(fig, use_container_width=True)

        st.markdown("<div class=\"panel\"><h3>Fraction Masses</h3></div>", unsafe_allow_html=True)
        fig, _ = plot_fractions(
            path_to_obs_dict={selected_run_dir.name: str(obs_dir)},
            path_to_save=str(plots_dir),
            save_plots=save_plots,
            close_fig=False,
        )
        st.pyplot(fig, use_container_width=True)

        st.markdown("<div class=\"panel\"><h3>Time Step Sizes</h3></div>", unsafe_allow_html=True)
        fig, _ = plot_time_step_sizes(
            path_to_obs=str(obs_dir),
            path_to_save=str(plots_dir),
            save_plots=save_plots,
            close_fig=False,
        )
        st.pyplot(fig, use_container_width=True)
    else:
        st.warning("Observation data not found for this run.")
