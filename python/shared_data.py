import numpy as np

# used by scripts that can only take one observation e.g. animate_unit_operations.py
path_to_obs = "data/run_fitted_no_MNP_hydroxyl_reactions_ERK_for_time_steps/obs"

# used by scripts that can take multiple observations e.g. plot_fractions.py

path_to_obs_dict = {#"FS³, fitted Langmuir, main simulation" : "data/run_fitted_main_simulation/obs",
                    #"FS³, fitted Langmuir, no Gouy-Chapman Model" : "data/run_fitted_no_gouy_chapman/obs",
                    "FS³, fitted Langmuir, no MNP hydroxyl reactions" : "data/run_fitted_no_MNP_hydroxyl_reactions/obs",
                    "FS³, fitted Langmuir, no MNP hydroxyl reactions, tau=1, discretization factor = 0.5" : "data/run_fitted_no_MNP_hydroxyl_reactions_tau_1_discretization_factor_0,5_ERK/obs",
                    "FS³, fitted Langmuir, no MNP hydroxyl reactions, tau=5, discretization factor = 0.2" : "data/run_fitted_no_MNP_hydroxyl_reactions_tau_5_discretization_factor_0,2_ERK/obs"
                    #"FS³, fitted Langmuir, liquid-slurry dispersion 1e-8" : "data/run_fitted_liquid_slurry_dispersion_1e-8/obs",
                    #"FS³, isotherm Langmuir, main simulation" : "data/run_isotherm_main_simulation/obs",
                    #"FS³, isotherm Langmuir, no MNP hydroxyl reactions" : "data/run_isotherm_no_MNP_hydroxyl_reactions/obs",
                    #"FS³, isotherm Langmuir, liquid-slurry dispersion 1e-8" : "data/run_isotherm_liquid_slurry_dispersion_1e-8/obs"
                    }

path_to_save = "data/run_fitted_no_MNP_hydroxyl_reactions_ERK_for_time_steps/plots"

show_plots = False

time_stamps = np.arange(0, 6600 + 1e-5, 2)

idx = {
    "H₂O": 0,
    "H⁺": 1,
    "OH⁻": 2,
    "Na⁺": 3,
    "Cl⁻": 4,
    "TrisH⁺": 5,
    "Tris": 6,
    "AcH": 7,
    "Ac⁻": 8,
    "GlyH₂⁺": 9,
    "GlyH": 10,
    "Gly⁻": 11,
    "hIgG": 12,
    "MNP-OH₂⁺": 13,
    "MNP-OH": 14,
    "MNP-O⁻": 15,
    "MNP-hIgG": 16,
    "MNP1": 17,
    "MNP2": 18,
    "MNP3": 19,
    "MNP4": 20,
    "MNP5": 21,
    "MNP6": 22,
    "MNP7": 23,
    "MNP8": 24,
    "MNP9": 25,
    "MNP10": 26
}