#include <sundials/sundials_types.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <memory>
#include <numeric>
#include <ostream>
#include <sstream>
#include <stdexcept>
#include <string>

#include "Components/Component.hpp"
#include "EigenDataTypes.hpp"
#include "Interpolators.tpp"
#include "Logger.hpp"
#include "Observers/NumpyIO.hpp"
#include "Observers/PHConverter.hpp"
#include "Observers/SnapshotObserver.hpp"
#include "Observers/TimeSeriesObserver.hpp"
#include "PhysicalConstants.hpp"
#include "Process.hpp"
#include "Reactions/ActivityModels.hpp"
#include "Reactions/MassActionLaw.hpp"
#include "Reactions/OneCellReaction.hpp"
#include "Reactions/Reaction.hpp"
#include "Reactions/ReactionSystem.hpp"
#include "Solver.hpp"
#include "UnitOperations/Inlet.hpp"
#include "UnitOperations/Pipe.hpp"
#include "UnitOperations/RS_MagneticCaptureProcessChamber.hpp"
#include "UnitOperations/Volume.hpp"

std::tuple<double, double, size_t>
run_HGMS_hIgG_separation(realtype kf_ion, realtype tau_reaction, bool save_obs,
                         realtype timeout_seconds, SolverType solverType,
                         realtype discretization_factor) {
#ifdef NDEBUG
  std::cout << "Release mode" << std::endl;
#else
  std::cout << "Debug mode" << std::endl;
#endif

  // Generally, everything happens in SI units, i.e., mol/m³ for concentrations,
  // m³ for volumes, m² for areas, etc. The only exceptions are the
  // Truesdell-Jones activity model parameters, which are in nm and L/mol and
  // that the values in the y vector which are in kg are converted to g/L and
  // mol/L for the reactions

  // =============================================================
  // ========================== COMPONENTS =======================
  // =============================================================
  ComponentSystem cs{
      {// Format:
       // Component(name).setCharge(...).setMolarMass(...).setTruesdellJonesAlpha(...).setTruesdellJonesBeta(...)
       Component("H₂O").setMolarMass(18.01528e-3),
       Component("H⁺")
           .setCharge(1)
           .setMolarMass(1.008e-3)
           .setTruesdellJonesParameters(4.78e-10, 0.24e-3),
       //  Component("Mg²⁺").setCharge(2).setTruesdellJonesParameters(5.5e-10,
       //  0.2e-3),  // Only for testing of Truessdell-Jones model -> plays no
       //  role in process
       Component("OH⁻")
           .setCharge(-1)
           .setMolarMass(17.007e-3)
           .setTruesdellJonesParameters(10.65e-10, 0.21e-3),
       Component("Na⁺")
           .setCharge(1)
           .setMolarMass(22.990e-3)
           .setTruesdellJonesParameters(4.32e-10, 0.06e-3),
       Component("Cl⁻")
           .setCharge(-1)
           .setMolarMass(35.453e-3)
           .setTruesdellJonesParameters(3.71e-10, 0.01e-3),
       Component("TrisH⁺")
           .setCharge(1)
           .setMolarMass(121.14e-3)
           .setTruesdellJonesParameters(4.0e-10, 0.0e-3),
       Component("Tris").setMolarMass(121.14e-3),
       Component("AcH").setMolarMass(60.052e-3),
       Component("Ac⁻")
           .setCharge(-1)
           .setMolarMass(59.044e-3)
           .setTruesdellJonesParameters(4.5e-10, 0.0e-3),
       Component("GlyH₂⁺").setCharge(1).setMolarMass(
           75.067e-3), // TODO: No truesdell_jones paramters available -> use
                       // Debye-Hückel model and see if it makes a difference
       Component("GlyH").setMolarMass(75.067e-3),
       Component("Gly⁻")
           .setCharge(-1)
           .setMolarMass(74.059e-3)
           .setTruesdellJonesParameters(4.0e-10, 0.0e-3),
       Component("hIgG").setMolarMass(
           150), // Human Immunoglobulin G (IgG) antibody, roughly 150 kDa
       Component("MNP-OH₂⁺")
           .setType(Type::MagneticNanoParticleGroup)
           .setMolarMass(
               18.015e-3), // MNP hydroxyl groups are not noted as charged
                           // because they should not be considered in the ionic
                           // strength of the activity model and in gouy chapman
                           // equation
       Component("MNP-OH")
           .setType(Type::MagneticNanoParticleGroup)
           .setMolarMass(17.007e-3),
       Component("MNP-O⁻")
           .setType(Type::MagneticNanoParticleGroup)
           .setMolarMass(
               15.999e-3), // MNP hydroxyl groups are not noted as charged
                           // because they should not be considered in the ionic
                           // strength of the activity model and in gouy chapman
                           // equation
       Component("MNP-hIgG")
           .setType(Type::MagneticNanoParticleGroup)
           .setMolarMass(150000e-3),
       Component("MNP1")
           .setType(Type::MagneticNanoParticle)
           .setRadius(1039e-9)
           .setDensity(5170)
           .setMagneticSaturation(3.5e5),
       Component("MNP2")
           .setType(Type::MagneticNanoParticle)
           .setRadius(1209e-9)
           .setDensity(5170)
           .setMagneticSaturation(3.5e5),
       Component("MNP3")
           .setType(Type::MagneticNanoParticle)
           .setRadius(1406e-9)
           .setDensity(5170)
           .setMagneticSaturation(3.5e5),
       Component("MNP4")
           .setType(Type::MagneticNanoParticle)
           .setRadius(1635e-9)
           .setDensity(5170)
           .setMagneticSaturation(3.5e5),
       Component("MNP5")
           .setType(Type::MagneticNanoParticle)
           .setRadius(1901e-9)
           .setDensity(5170)
           .setMagneticSaturation(3.5e5),
       Component("MNP6")
           .setType(Type::MagneticNanoParticle)
           .setRadius(2211e-9)
           .setDensity(5170)
           .setMagneticSaturation(3.5e5),
       Component("MNP7")
           .setType(Type::MagneticNanoParticle)
           .setRadius(2571e-9)
           .setDensity(5170)
           .setMagneticSaturation(3.5e5),
       Component("MNP8")
           .setType(Type::MagneticNanoParticle)
           .setRadius(2990e-9)
           .setDensity(5170)
           .setMagneticSaturation(3.5e5),
       Component("MNP9")
           .setType(Type::MagneticNanoParticle)
           .setRadius(3477e-9)
           .setDensity(5170)
           .setMagneticSaturation(3.5e5),
       Component("MNP10")
           .setType(Type::MagneticNanoParticle)
           .setRadius(4043e-9)
           .setDensity(5170)
           .setMagneticSaturation(3.5e5)}};

  // Other physical parameters
  const realtype MNP_Ns =
      1.08 * 1e-5; // [mol/m²] - number of hydroxyl groups per area
  // TODO: Where does this come from? Is radius correct or only nm scale? Then
  // it would fit better with MNP_specific_surface
  const realtype MNP_specific_surface = 1e5; // m²/kg

  // =============================================================
  // ========================== REACTIONS ========================
  // =============================================================

  // NoActivityModel acitivityModel(cs);
  TruesdellJonesActivityModel acitivityModel(cs);
  ReactionSystem reactionSystem(cs, acitivityModel);

  realtype K_W = std::pow(10, -14) * 1e6; // [mol²/(m³)²] - ion product of water
  realtype c_W =
      1000 / cs("H₂O").molarMass; // [mol/m³] - Concentration of water

  realtype Ks_W =
      K_W /
      c_W; // [mol/m³] - Equilibrium constant for water dissociation in mol/m³
  realtype Ks_TrisH_plus = std::pow(10, -8.3) * 1e3;  // [mol/m³]
  realtype Ks_AcH = std::pow(10, -4.75) * 1e3;        // [mol/m³]
  realtype Ks_GlyH2_plus = std::pow(10, -2.35) * 1e3; // [mol/m³]
  realtype Ks_GlyH = std::pow(10, -9.78) * 1e3;       // [mol/m³]
  realtype Ks_MNP_OH2_plus = std::pow(10, -(7.9 - 0.5 * 2.5)) * 1e3; // [mol/m³]
  realtype Ks_MNP_OH = std::pow(10, -(7.9 + 0.5 * 2.5)) * 1e3; // [mol/m³]

  realtype kb_W = (kf_ion / Ks_W);                       // [m³/(mol*s)]
  realtype kb_TrisH_plus = (kf_ion / Ks_TrisH_plus);     // [m³/(mol*s)]
  realtype kb_AcH = (kf_ion / Ks_AcH);                   // [m³/(mol*s)]
  realtype kb_GlyH2_plus = (kf_ion / Ks_GlyH2_plus);     // [m³/(mol*s)]
  realtype kb_GlyH = (kf_ion / Ks_GlyH);                 // [m³/(mol*s)]
  realtype kb_MNP_OH2_plus = (kf_ion / Ks_MNP_OH2_plus); // [m³/(mol*s)]
  realtype kb_MNP_OH = (kf_ion / Ks_MNP_OH);             // [m³/(mol*s)]

  // reactionSystem.add(massActionLaw(cs, "H₂O <=> H⁺ + OH⁻", kf_ion, kb_W,
  // cs.getIdx("H⁺"))); reactionSystem.add(massActionLaw(cs, "TrisH⁺ <=> Tris +
  // H⁺", kf_ion, kb_Tris, cs.getIdx("H⁺")));
  // reactionSystem.add(massActionLaw(cs, "AcH <=> Ac⁻ + H⁺", kf_ion, kb_AcH,
  // cs.getIdx("H⁺"))); reactionSystem.add(massActionLaw(cs, "GlyH₂⁺ <=> GlyH +
  // H⁺", kf_ion, kb_GlyH2, cs.getIdx("H⁺")));
  // reactionSystem.add(massActionLaw(cs, "GlyH <=> Gly⁻ + H⁺", kf_ion, kb_GlyH,
  // cs.getIdx("H⁺")));

  // reactionSystem.add(massActionLaw_damped(cs, "H₂O <=> H⁺ + OH⁻", kf_ion,
  // kb_W, tau_reaction, cs.getIdx("H⁺")));
  // reactionSystem.add(massActionLaw_damped(cs, "TrisH⁺ <=> Tris + H⁺", kf_ion,
  // kb_Tris, tau_reaction, cs.getIdx("H⁺")));
  // reactionSystem.add(massActionLaw_damped(cs, "AcH <=> Ac⁻ + H⁺", kf_ion,
  // kb_AcH, tau_reaction, cs.getIdx("H⁺")));
  // reactionSystem.add(massActionLaw_damped(cs, "GlyH₂⁺ <=> GlyH + H⁺", kf_ion,
  // kb_GlyH2, tau_reaction, cs.getIdx("H⁺")));
  // reactionSystem.add(massActionLaw_damped(cs, "GlyH <=> Gly⁻ + H⁺", kf_ion,
  // kb_GlyH, tau_reaction, cs.getIdx("H⁺")));

  reactionSystem.add(massActionLaw_inverseRatePrediction(
      cs, "H₂O <=> H⁺ + OH⁻", Ks_W, tau_reaction, cs.getIdx("H⁺")));
  reactionSystem.add(massActionLaw_inverseRatePrediction(
      cs, "TrisH⁺ <=> Tris + H⁺", Ks_TrisH_plus, tau_reaction,
      cs.getIdx("H⁺")));
  reactionSystem.add(massActionLaw_inverseRatePrediction(
      cs, "AcH <=> Ac⁻ + H⁺", Ks_AcH, tau_reaction, cs.getIdx("H⁺")));
  reactionSystem.add(massActionLaw_inverseRatePrediction(
      cs, "GlyH₂⁺ <=> GlyH + H⁺", Ks_GlyH2_plus, tau_reaction,
      cs.getIdx("H⁺")));
  reactionSystem.add(massActionLaw_inverseRatePrediction(
      cs, "GlyH <=> Gly⁻ + H⁺", Ks_GlyH, tau_reaction, cs.getIdx("H⁺")));

  // TODO: maybe try smaller time constant for gouy chapman reaction
  // auto reaction_OH2_plus = massActionLaw(cs, "MNP-OH₂⁺ <=> MNP-OH + H⁺",
  // kf_ion, kb_MNP_OH2_plus, cs.getIdx("H⁺")); auto reaction_OH =
  // massActionLaw(cs, "MNP-OH <=> MNP-O⁻ + H⁺", kf_ion, kb_MNP_OH,
  // cs.getIdx("H⁺"));

  // auto reaction_OH2_plus = massActionLaw_damped(cs, "MNP-OH₂⁺ <=> MNP-OH +
  // H⁺", kf_ion, kb_MNP_OH2_plus, tau_reaction,
  //                                               cs.getIdx("H⁺"));
  // auto reaction_OH = massActionLaw_damped(cs, "MNP-OH <=> MNP-O⁻ + H⁺",
  // kf_ion, kb_MNP_OH, tau_reaction, cs.getIdx("H⁺"));

  //     auto reaction_OH2_plus = massActionLaw_inverseRatePrediction(cs,
  //     "MNP-OH₂⁺ <=> MNP-OH + H⁺", Ks_MNP_OH2_plus,
  //                                                                  tau_reaction,
  //                                                                  cs.getIdx("H⁺"));
  //     auto reaction_OH = massActionLaw_inverseRatePrediction(cs, "MNP-OH <=>
  //     MNP-O⁻ + H⁺", Ks_MNP_OH, tau_reaction,
  //                                                            cs.getIdx("H⁺"));

  //     auto reaction_H_plus_gouy_chapman =
  //         [  // Precompute indices directly in capture list
  //             idx_MNP_OH2_plus = cs.getIdx("MNP-OH₂⁺"), idx_MNP_OH =
  //             cs.getIdx("MNP-OH"), idx_MNP_O_minus = cs.getIdx("MNP-O⁻"),
  //             idx_H_plus = cs.getIdx("H⁺"), idx_MNP1 = cs.getIdx("MNP1"),
  //             idx_MNP10 = cs.getIdx("MNP10"),
  //             // Precompute charge mask - only works for ions with z = 1 or z
  //             = -1 charges_nonzero_mask = (cs.getCharges() !=
  //             0).cast<realtype>(),
  //             // Precompute physical constants
  //             surface_charge_factor = constants::elementary_charge *
  //             constants::avogadro * MNP_Ns, normalization_factor = 2 *
  //             cs.relative_permittivity * constants::vacuum_permittivity *
  //                                    constants::boltzmann * cs.temperature,
  //             MNP_specific_surface, &reaction_OH2_plus,
  //             &reaction_OH](realtype t, const ConstArrayMap& concentrations,
  //             ArrayMap& activities) { const auto c_MNP_OH2_plus =
  //             concentrations.col(idx_MNP_OH2_plus); const auto c_MNP_OH =
  //             concentrations.col(idx_MNP_OH); const auto c_MNP_O_minus =
  //             concentrations.col(idx_MNP_O_minus); const auto c_H_plus =
  //             concentrations.col(
  //                 idx_H_plus);  // Using H⁺ concentration here because
  //                 original Gouy Chapman was derived like this -> TODO:
  //                 discuss newer approaches of Martin Bazant et al.

  //             const auto surface_charge_raw = surface_charge_factor *
  //             ((c_MNP_OH2_plus - c_MNP_O_minus) /
  //                                                                      (c_MNP_OH2_plus
  //                                                                      +
  //                                                                      c_MNP_OH
  //                                                                      +
  //                                                                      c_MNP_O_minus));

  //             // Take absolute value here has the same effect as choosing the
  //             right sign for the sqrt later const auto surface_charge =
  //             surface_charge_raw.isFinite().select(surface_charge_raw,
  //                                                                              0.0);  // Replace non-finite with 0.0

  //             const auto c_MNP = activities.middleCols(idx_MNP1, idx_MNP10 -
  //             idx_MNP1 + 1).rowwise().sum(); const auto
  //             surface_charge_MNP_alternative_raw =
  //             (constants::elementary_charge * constants::avogadro /
  //                                                              MNP_specific_surface)
  //                                                              *
  //                                                             ((c_MNP_OH2_plus
  //                                                             -
  //                                                             c_MNP_O_minus)
  //                                                             / c_MNP);
  //             const auto surface_charge_MNP_alternative =
  //             surface_charge_MNP_alternative_raw.isFinite()
  //                                                             .select(surface_charge_MNP_alternative_raw,
  //                                                                     0.0);
  //                                                                     //
  //                                                                     Replace
  //                                                                     non-finite
  //                                                                     with
  //                                                                     0.0

  //             // TODO: Discuss if n_0 -> 0, f -> 0 is that correct?

  //             const auto n_0 = constants::avogadro *
  //             (concentrations.rowwise() *
  //             charges_nonzero_mask).rowwise().sum();

  //             const auto normalized_surface_charge = 0.5 * (surface_charge /
  //             ((normalization_factor * n_0).sqrt()));

  //             const auto f_raw = (-normalized_surface_charge +
  //                                 (normalized_surface_charge *
  //                                 normalized_surface_charge + 1).sqrt())
  //                                    .square();

  //             const auto f =
  //             normalized_surface_charge.isFinite().select(f_raw, 0.0);

  //             const auto c_H_plus_surface = f * c_H_plus;

  //             activities.col(idx_H_plus) = c_H_plus_surface;

  // #if LOG_ENABLED
  //             static int call_count = 0;
  //             if (call_count < LOG_FIRST_N_CALLS || call_count %
  //             LOG_EVERY_N_CALLS == 0) {
  //                 std::ostringstream oss;
  //                 oss << "Call count: " << call_count << " at t= " << t <<
  //                 "\n"; oss << "H⁺ concentration: " << c_H_plus.transpose()
  //                 << "\n"; oss << "MNP-OH₂⁺ concentration: " <<
  //                 c_MNP_OH2_plus.transpose() << "\n"; oss << "MNP-OH
  //                 concentration: " << c_MNP_OH.transpose() << "\n"; oss <<
  //                 "MNP-O⁻ concentration: " << c_MNP_O_minus.transpose() <<
  //                 "\n"; oss << "Raw Surface charge: " <<
  //                 surface_charge_raw.transpose() << "\n"; oss << "Suface
  //                 charge: " << surface_charge.transpose() << "\n"; oss <<
  //                 "Surface Charge MNP Alternative: " <<
  //                 surface_charge_MNP_alternative.transpose() << "\n"; oss <<
  //                 "N_0: " << n_0.transpose() << "\n"; oss << "Normalized
  //                 surface charge: " << normalized_surface_charge.transpose()
  //                 << "\n"; oss << "f_raw: " << f_raw.transpose() << "\n"; oss
  //                 << "f: " << f.transpose() << "\n\n";

  //                 LOG("MNP_H+_gouy_chapman.log", oss.str());
  //             }
  //             call_count++;
  // #endif
  //         };

  //     // amphoteric MNP hydroxl group reactions
  //     reactionSystem.add(Reaction(
  //         [=](realtype t, const ConstArrayMap& concentrations, const
  //         ConstArrayMap& activities, ArrayMap& dc_dt) {
  //             Array modified_activities{activities};
  //             ArrayMap modified_activities_map(modified_activities.data(),
  //             activities.rows(), activities.cols(),
  //                                              PhaseStride(activities.cols(),
  //                                              1));
  //             ConstArrayMap
  //             modified_activities_const_map(modified_activities.data(),
  //             activities.rows(),
  //                                                         activities.cols(),
  //                                                         PhaseStride(activities.cols(),
  //                                                         1));

  //             reaction_H_plus_gouy_chapman(t, concentrations,
  //             modified_activities_map); reaction_OH2_plus.rhs(t,
  //             concentrations, modified_activities_const_map, dc_dt);
  //         },

  //         [=](realtype t, const ConstArrayMap& concentrations, const
  //         ConstArrayMap& activities) -> realtype {
  //             Array modified_activities{activities};
  //             ArrayMap modified_activities_map(modified_activities.data(),
  //             activities.rows(), activities.cols(),
  //                                              PhaseStride(activities.cols(),
  //                                              1));
  //             ConstArrayMap
  //             modified_activities_const_map(modified_activities.data(),
  //             activities.rows(),
  //                                                         activities.cols(),
  //                                                         PhaseStride(activities.cols(),
  //                                                         1));

  //             reaction_H_plus_gouy_chapman(t, concentrations,
  //             modified_activities_map); auto error_OH2_plus =
  //             reaction_OH2_plus.errorFunction(t, concentrations,
  //             modified_activities_const_map); return error_OH2_plus;
  //         }));

  //     reactionSystem.add(Reaction(
  //         [=](realtype t, const ConstArrayMap& concentrations, const
  //         ConstArrayMap& activities, ArrayMap& dc_dt) {
  //             Array modified_activities{activities};
  //             ArrayMap modified_activities_map(modified_activities.data(),
  //             activities.rows(), activities.cols(),
  //                                              PhaseStride(activities.cols(),
  //                                              1));
  //             ConstArrayMap
  //             modified_activities_const_map(modified_activities.data(),
  //             activities.rows(),
  //                                                         activities.cols(),
  //                                                         PhaseStride(activities.cols(),
  //                                                         1));

  //             reaction_H_plus_gouy_chapman(t, concentrations,
  //             modified_activities_map); reaction_OH.rhs(t, concentrations,
  //             modified_activities_const_map, dc_dt);
  //         },

  //         [=](realtype t, const ConstArrayMap& concentrations, const
  //         ConstArrayMap& activities) -> realtype {
  //             Array modified_activities{activities};
  //             ArrayMap modified_activities_map(modified_activities.data(),
  //             activities.rows(), activities.cols(),
  //                                              PhaseStride(activities.cols(),
  //                                              1));
  //             ConstArrayMap
  //             modified_activities_const_map(modified_activities.data(),
  //             activities.rows(),
  //                                                         activities.cols(),
  //                                                         PhaseStride(activities.cols(),
  //                                                         1));

  //             reaction_H_plus_gouy_chapman(t, concentrations,
  //             modified_activities_map); auto error_OH =
  //             reaction_OH.errorFunction(t, concentrations,
  //             modified_activities_const_map);

  //             return error_OH;
  //         }));

  // Inverse fitted Langmuir parameters
  ColVector langmuir_pH(7);
  langmuir_pH << 2.0, 2.5, 2.6, 3.5, 3.7, 4.5, 7.0;
  ColVector langmuir_K_B(7);
  langmuir_K_B << 0.03, 2.14, 2.28, 3.22, 3.84, 35.15, 36.52; // [m³/mol]
  langmuir_K_B = langmuir_K_B;
  ColVector langmuir_q_max(7);
  langmuir_q_max << 0.00, 0.01, 0.02, 0.05, 0.13, 0.28, 0.31; // [g/g]

  // // Small scale isotherms
  // ColVector langmuir_pH(5);
  // langmuir_pH << 2.9, 3.5, 4.5, 7.0, 7.5;
  // ColVector langmuir_K_B(5);
  // langmuir_K_B << 53.28, 53.52, 93.5, 95.98, 42.63;  // [m³/mol]
  // ColVector langmuir_q_max(5);
  // langmuir_q_max << 0.04, 0.03, 0.18, 0.13, 0.15;  // [g/g]
  LinearInterpolator<ColVector> langmuir_K_B_interpolator(langmuir_pH,
                                                          langmuir_K_B);
  LinearInterpolator<ColVector> langmuir_q_max_interpolator(langmuir_pH,
                                                            langmuir_q_max);

  // Langmuir adsorption reaction for MNP-HIgG
  reactionSystem.add(Reaction(
      [
          // Precompute indices directly in capture list
          idx_H_plus = cs.getIdx("H⁺"), idx_MNP1 = cs.getIdx("MNP1"),
          idx_MNP10 = cs.getIdx("MNP10"), idx_MNP_hIgG = cs.getIdx("MNP-hIgG"),
          idx_hIgG = cs.getIdx("hIgG"), M_hIgG = cs("hIgG").molarMass,
          langmuir_K_B_interpolator, langmuir_q_max_interpolator,
          tau_reaction](realtype t, const ConstArrayMap &concentrations,
                        const ConstArrayMap &activities, ArrayMap &dc_dt) {
        // TODO: Maybe change to true activity based pH calculation later
        const auto pH =
            -(concentrations.col(idx_H_plus)).log10() +
            3; // +3 because pH is defined w.r.t. mol/L instead of mol/m³
        const auto pH_activity = -(activities.col(idx_H_plus)).log10() + 3;
        LOG("pH_conc_vs_pH_activity.log",
            "Time: " + std::to_string(t) +
                ", pH_conc: " + std::to_string(pH(0)) +
                ", pH_activity: " + std::to_string(pH_activity(0)) + "\n");
        const auto c_MNP =
            activities.middleCols(idx_MNP1, idx_MNP10 - idx_MNP1 + 1)
                .rowwise()
                .sum();
        const auto c_hIgG = activities.col(idx_hIgG);
        const auto c_MNP_hIgG = activities.col(idx_MNP_hIgG);

        auto dc_dt_MNP_hIgG = dc_dt.col(idx_MNP_hIgG);
        auto dc_dt_hIgG = dc_dt.col(idx_hIgG);

        // k_ads/k_des = K_B must be in m³/mol (not in L/g)
        const auto K_B_raw =
            pH.unaryExpr(langmuir_K_B_interpolator); // in L/g = m³/kg
        const auto K_B = K_B_raw * M_hIgG;           // in m³/mol

        const auto q_max = pH.unaryExpr(langmuir_q_max_interpolator);

        // Devide by molar mass of Protein to account for q being m_MNP_hIgG /
        // m_MNP -> (q_max - c_MNP_hIgG / c_MNP) * c_MNP / M_hIgG gives molar
        // concentration of binding sites rate = k_ads * c_hIgG * (q_max - q) *
        // (c_MNP / M_hIgG) - k_des * c_MNP_hIgG = 0 at equilibrium k_ads +
        // (c_hIgG - r) * (q_max - (C_MNP_hIgG + r) / c_MNP) * (c_MNP / M_hIgG)
        // - k_des * (c_MNP_hIgG + r) = 0

        const auto p = c_MNP_hIgG - c_hIgG - c_MNP * q_max - M_hIgG / K_B;
        const auto q_tilde =
            c_hIgG * (c_MNP * q_max - c_MNP_hIgG) - (M_hIgG * c_MNP_hIgG / K_B);

        const auto radicant = (p.square() - 4.0 * q_tilde).sqrt();

        const auto delta_xi_2 = (-p - radicant) / 2.0;

        // Assuming delta_xi_2 is the physically correct solution
        const auto net_rate = delta_xi_2 / tau_reaction;
        dc_dt_MNP_hIgG += net_rate;
        dc_dt_hIgG -= net_rate;

#if LOG_ENABLED || !defined(NDEBUG)
        // const auto residual_current = K_B * c_hIgG * (q_max - c_MNP_hIgG /
        // c_MNP) * (c_MNP / M_hIgG) - c_MNP_hIgG;
        const auto Q_current =
            c_MNP_hIgG /
            (c_hIgG * (q_max - c_MNP_hIgG / c_MNP) * (c_MNP / M_hIgG));

        const auto q_raw = c_MNP_hIgG / c_MNP;
        const auto q = (q_raw.isFinite() && (q_raw >= 0.0) && (q_raw <= q_max))
                           .select(q_raw, q_max);

        const auto delta_xi_1 = (-p + radicant) / 2.0;

        const auto c_hIgG_eq_1 = c_hIgG - delta_xi_1;
        const auto c_MNP_hIgG_eq_1 = c_MNP_hIgG + delta_xi_1;
        // const auto residual_eq_1 = K_B * c_hIgG_eq_1 * (q_max -
        // c_MNP_hIgG_eq_1 / c_MNP) * (c_MNP / M_hIgG) -
        //                            c_MNP_hIgG_eq_1;
        const auto Q_eq_1 =
            c_MNP_hIgG_eq_1 / (c_hIgG_eq_1 * (q_max - c_MNP_hIgG_eq_1 / c_MNP) *
                               (c_MNP / M_hIgG));
        // Set NaN/Inf values to 0 to avoid assertion failures when dividing by
        // zero
        const auto Q_relative_error_1_raw = (Q_eq_1 - K_B).cwiseAbs() / K_B;
        const auto Q_relative_error_1 =
            Q_relative_error_1_raw.isFinite().select(Q_relative_error_1_raw,
                                                     0.0);

        const auto c_hIgG_eq_2 = c_hIgG - delta_xi_2;
        const auto c_MNP_hIgG_eq_2 = c_MNP_hIgG + delta_xi_2;
        // const auto residual_eq_2 = K_B * c_hIgG_eq_2 * (q_max -
        // c_MNP_hIgG_eq_2 / c_MNP) * (c_MNP / M_hIgG) -
        //                            c_MNP_hIgG_eq_2;
        const auto Q_eq_2 =
            c_MNP_hIgG_eq_2 / (c_hIgG_eq_2 * (q_max - c_MNP_hIgG_eq_2 / c_MNP) *
                               (c_MNP / M_hIgG));
        // Set NaN/Inf values to 0 to avoid assertion failures when dividing by
        // zero
        const auto Q_relative_error_2_raw = (Q_eq_2 - K_B).cwiseAbs() / K_B;
        const auto Q_relative_error_2 =
            Q_relative_error_2_raw.isFinite().select(Q_relative_error_2_raw,
                                                     0.0);
#endif
#if LOG_ENABLED

        static int call_count = 0;
        if (call_count < LOG_FIRST_N_CALLS ||
            call_count % LOG_EVERY_N_CALLS == 0) {
          std::ostringstream oss;
          oss.precision(6);
          oss << std::scientific;
          oss << "MNP-hIgG adsorption reaction call " << call_count
              << " at t = " << t << " seconds\n";
          oss << "pH:\n" << pH.transpose() << "\n";
          oss << "c_MNP:\n" << c_MNP.transpose() << "\n";
          oss << "c_hIgG:\n" << c_hIgG.transpose() << "\n";
          oss << "c_MNP_hIgG:\n" << c_MNP_hIgG.transpose() << "\n";
          oss << "q_raw:\n" << q_raw.transpose() << "\n";
          oss << "q:\n" << q.transpose() << "\n";
          oss << "K_B:\n" << K_B.transpose() << "\n";
          oss << "q_max:\n" << q_max.transpose() << "\n";
          // oss << "residual_current:\n" << residual_current.transpose() <<
          // "\n";
          oss << "Q_current:\n" << Q_current.transpose() << "\n";
          oss << "p:\n" << p.transpose() << "\n";
          oss << "q_tilde:\n" << q_tilde.transpose() << "\n";
          oss << "Radicant:\n" << radicant.transpose() << "\n";
          oss << "delta_xi_1:\n" << delta_xi_1.transpose() << "\n";
          oss << "c_hIgG_eq_1:\n" << c_hIgG_eq_1.transpose() << "\n";
          oss << "c_MNP_hIgG_eq_1:\n" << c_MNP_hIgG_eq_1.transpose() << "\n";
          // oss << "residual_eq_1:\n" << residual_eq_1.transpose() << "\n";
          oss << "Q_eq_1:\n" << Q_eq_1.transpose() << "\n";
          oss << "Relative error Q_eq_1 vs K_B:\n"
              << Q_relative_error_1.transpose() << "\n";
          oss << "delta_xi_2:\n" << delta_xi_2.transpose() << "\n";
          oss << "c_hIgG_eq_2:\n" << c_hIgG_eq_2.transpose() << "\n";
          oss << "c_MNP_hIgG_eq_2:\n" << c_MNP_hIgG_eq_2.transpose() << "\n";
          // oss << "residual_eq_2:\n" << residual_eq_2.transpose() << "\n";
          oss << "Q_eq_2:\n" << Q_eq_2.transpose() << "\n";
          oss << "Relative error Q_eq_2 vs K_B:\n"
              << Q_relative_error_2.transpose() << "\n";
          oss << "dr/dt:\n" << net_rate.transpose() << "\n\n";

          LOG("MNP_hIgG_adsorption_rhs.log", oss.str());
        }
        call_count++;
#endif
        assert(((c_hIgG_eq_2 + 1e-12) >= 0).all() &&
               ((c_MNP_hIgG_eq_2 + 1e-12) >= 0).all() &&
               "Equilibrium concentrations (solution 2) must be non-negative");

        assert(((c_hIgG_eq_1 - 1e-12) < 0).all() ||
               ((c_MNP_hIgG_eq_1 - 1e-12) < 0).all() &&
                   "Equilibrium concentrations (solution 1) should contain "
                   "negative values");

        // Skip these assertions because they are fragile e.g. for extremely
        // dilute solutions, etc.

        // assert((Q_relative_error_2 < 1e-1).all() &&
        //        "Equilibrium concentrations (solution 2) must satisfy
        //        equilibrium condition");

        // assert((residual_eq_1.abs() < 1e-6).all() &&
        //        "Equilibrium concentrations (solution 1) must satisfy
        //        equilibrium condition");

        // assert((residual_eq_2.abs() < residual_current.abs()).all() &&
        //        "Equilibrium concentrations (solution 2) must be closer to
        //        equilibrium than current concentrations");
      }));

  // =============================================================
  // ========================== SOLUTIONS ========================
  // =============================================================

  const realtype buffer_eps =
      0.0; // To avoid numerical issues with zero concentrations
  const realtype t_reaction_duration =
      100.0; // Time for oneCellReaction to reach equilibrium

  // FEED
  realtype c_MNP_feed = 16.5;      // kg/m³
  realtype V_MNPs_feed = 0.000588; // m³

  realtype c_hIgG_feed = 0.45;   // kg/m³ // TODO: vs 0.47 in Matlab?
  realtype V_hIgG_feed = 0.0046; // m³ // TODO: vs 0.0044 in Matlab?

  realtype c_HCl_feed = 0.01835e3;                               // mol/m³
  realtype c_Tris_feed = 0.02e3;                                 // mol/m³
  realtype c_NaCl_feed = 0.15e3;                                 // mol/m³
  realtype V_buffer_feed = 0.000135 + 0.006 - 0.6e-3 - 0.873e-3; // m³

  realtype V_feed = V_MNPs_feed + V_hIgG_feed + V_buffer_feed;

  // Update concentrations due to dilution
  c_MNP_feed *= V_MNPs_feed / V_feed;
  c_hIgG_feed *= V_hIgG_feed / V_feed;
  c_HCl_feed *= V_buffer_feed / V_feed;
  c_Tris_feed *= V_buffer_feed / V_feed;
  c_NaCl_feed *= V_buffer_feed / V_feed;

  RowVector feed = RowVector(cs.n_components);
  feed.setConstant(buffer_eps);
  feed(cs.getIdx("H₂O")) = c_W; // ~55.5 mol/L
  feed(cs.getIdx("H⁺")) = c_HCl_feed;
  feed(cs.getIdx("Na⁺")) = c_NaCl_feed;
  feed(cs.getIdx("Cl⁻")) = c_HCl_feed + c_NaCl_feed;
  feed(cs.getIdx("Tris")) = c_Tris_feed;
  feed(cs.getIdx("hIgG")) = c_hIgG_feed;

  realtype mnp_distribution[10] = {0.04852, 0.1853,  0.25496, 0.14803, 0.03495,
                                   0.03743, 0.09137, 0.10777, 0.068,   0.02367};

  if (std::abs(std::accumulate(mnp_distribution, mnp_distribution + 10, 0.0) -
               1.0) > 1e-3) {
    throw std::runtime_error(
        "Error: MNP distribution does not sum to 1 but to " +
        std::to_string(
            std::accumulate(mnp_distribution, mnp_distribution + 10, 0.0)) +
        "!");
  }

  for (int i = 0; i < 10; ++i) {
    feed(cs.getIdx("MNP" + std::to_string(i + 1))) =
        c_MNP_feed * mnp_distribution[i];
  }

  feed(cs.getIdx("MNP-OH")) =
      MNP_Ns * MNP_specific_surface *
      c_MNP_feed; // MNP_Ns * MNP_specific_surface is in mol/kg, multiply with
                  // c_MNP_feed in kg/m³ to get mol/m³

  std::cout << "Feed before reaction: " << feed << std::endl;
  realtype feed_error;
  std::tie(feed, feed_error) =
      oneCellReaction(reactionSystem, feed, t_reaction_duration * 10,
                      solverType, timeout_seconds);
  std::cout << "Feed after reaction: " << feed << std::endl;
  const auto pH_feed = -std::log10(feed(cs.getIdx("H⁺")) * 1e-3);
  std::cout << "pH feed: " << pH_feed << " (should be ~7.3)" << std::endl;
  const auto q_feed =
      feed(cs.getIdx("MNP-hIgG")) /
      (feed.middleCols(cs.getIdx("MNP1"),
                       cs.getIdx("MNP10") - cs.getIdx("MNP1") + 1)
           .rowwise()
           .sum());
  const auto q_max = langmuir_q_max_interpolator(pH_feed);
  std::cout << "q feed: " << q_feed << " q_max at this pH: " << q_max
            << std::endl;
  std::cout << "Feed reaction error: " << feed_error << std::endl;

  // BUFFER 1 (TBS pH 7)
  realtype c_HCl_buffer1 = 0.019156e3; // mol/m³
  realtype c_Tris_buffer1 = 0.02e3;    // mol/m³
  realtype c_NaCl_buffer1 = 0.15e3;    // mol/m³

  RowVector buffer1 = RowVector(cs.n_components);
  buffer1.setConstant(buffer_eps);
  buffer1(cs.getIdx("H₂O")) = c_W;
  buffer1(cs.getIdx("H⁺")) = c_HCl_buffer1;
  buffer1(cs.getIdx("Na⁺")) = c_NaCl_buffer1;
  buffer1(cs.getIdx("Cl⁻")) = c_HCl_buffer1 + c_NaCl_buffer1;
  buffer1(cs.getIdx("Tris")) = c_Tris_buffer1;

  std::cout << "B1 before reaction: " << buffer1 << std::endl;
  realtype buffer1_error;
  std::tie(buffer1, buffer1_error) =
      oneCellReaction(reactionSystem, buffer1, t_reaction_duration, solverType,
                      timeout_seconds);
  std::cout << "B1 after reaction: " << buffer1 << std::endl;
  auto pH_B1 = -std::log10(buffer1(cs.getIdx("H⁺")) * 1e-3);
  std::cout << "pH B1: " << pH_B1 << " (should be ~7)" << std::endl;
  std::cout << "B1 reaction error: " << buffer1_error << std::endl;

  // BUFFER 2 (NaAc pH 4.75 (Goal/Paper) / 4.77 (Matlab))
  realtype c_HCl_buffer2 = 0.072107e3; // mol/m³
  realtype c_NaAc_buffer2 = 0.2e3;     // mol/m³

  RowVector buffer2 = RowVector(cs.n_components);
  buffer2.setConstant(buffer_eps);
  buffer2(cs.getIdx("H₂O")) = c_W;
  buffer2(cs.getIdx("H⁺")) = c_HCl_buffer2;
  buffer2(cs.getIdx("Na⁺")) = c_NaAc_buffer2;
  buffer2(cs.getIdx("Cl⁻")) = c_HCl_buffer2;
  buffer2(cs.getIdx("Ac⁻")) = c_NaAc_buffer2;

  std::cout << "B2 before reaction: " << buffer2 << std::endl;
  realtype buffer2_error;
  std::tie(buffer2, buffer2_error) =
      oneCellReaction(reactionSystem, buffer2, t_reaction_duration, solverType,
                      timeout_seconds);
  std::cout << "B2 after reaction: " << buffer2 << std::endl;
  auto pH_B2 = -std::log10(buffer2(cs.getIdx("H⁺")) * 1e-3);
  std::cout << "pH B2: " << pH_B2
            << " (should be ~4.75 (Goal/Paper) / 4.77 (Matlab))" << std::endl;
  std::cout << "B2 reaction error: " << buffer2_error << std::endl;

  // BUFFER 3 (NaAc+Gly pH 2.5 (Goal/Paper) / 2.41 (Matlab))
  realtype c_HCl_buffer3 = 0.20221e3; // mol/m³
  realtype c_NaAc_buffer3 = 0.2e3;    // mol/m³

  RowVector buffer3 = RowVector(cs.n_components);
  buffer3.setConstant(buffer_eps);
  buffer3(cs.getIdx("H₂O")) = c_W;
  buffer3(cs.getIdx("H⁺")) = c_HCl_buffer3;
  buffer3(cs.getIdx("Na⁺")) = c_NaAc_buffer3;
  buffer3(cs.getIdx("Cl⁻")) = c_HCl_buffer3;
  buffer3(cs.getIdx("Ac⁻")) = c_NaAc_buffer3;

  std::cout << "B3 before reaction: " << buffer3 << std::endl;
  realtype buffer3_error;
  std::tie(buffer3, buffer3_error) =
      oneCellReaction(reactionSystem, buffer3, t_reaction_duration, solverType,
                      timeout_seconds);
  std::cout << "B3 after reaction: " << buffer3 << std::endl;
  auto pH_B3 = -std::log10(buffer3(cs.getIdx("H⁺")) * 1e-3);
  std::cout << "pH B3: " << pH_B3
            << " (should be ~2.5 (Goal/Paper) / 2.41 (Matlab))" << std::endl;
  std::cout << "B3 reaction error: " << buffer3_error << std::endl;

  // BUFFER WATER (Empty/Pure Water)
  realtype pH_buffer5 = 6.5; // TODO: Maybe adjust for pH value?

  RowVector buffer5_water = RowVector(cs.n_components);
  buffer5_water.setConstant(buffer_eps);
  buffer5_water(cs.getIdx("H₂O")) = c_W;

  std::cout << "B5 before reaction: " << buffer5_water << std::endl;
  realtype buffer5_water_error;
  std::tie(buffer5_water, buffer5_water_error) =
      oneCellReaction(reactionSystem, buffer5_water, t_reaction_duration,
                      solverType, timeout_seconds);
  std::cout << "B5 after reaction: " << buffer5_water << std::endl;
  auto pH_B5 = -std::log10(buffer5_water(cs.getIdx("H⁺")) * 1e-3);
  std::cout << "pH B5: " << pH_B5 << " (should be ~7.0 / 6.5)" << std::endl;
  std::cout << "B5 reaction error: " << buffer5_water_error << std::endl;

  reactionSystem.max_error =
      0.0; // Set back to default for following main simulation

  // =============================================================
  // ========================== RECIPE ===========================
  // =============================================================

  auto frac_feed_1 = std::make_shared<Volume>(reactionSystem);
  auto frac_feed_2 = std::make_shared<Volume>(reactionSystem);
  auto frac_wash_1 = std::make_shared<Volume>(reactionSystem);
  auto frac_wash_2 = std::make_shared<Volume>(reactionSystem);
  auto frac_wash_3 = std::make_shared<Volume>(reactionSystem);
  auto frac_elution_1 = std::make_shared<Volume>(reactionSystem);
  auto frac_elution_2 = std::make_shared<Volume>(reactionSystem);
  auto frac_elution_3 = std::make_shared<Volume>(reactionSystem);
  auto frac_elution_4 = std::make_shared<Volume>(reactionSystem);
  auto frac_elution_5 = std::make_shared<Volume>(reactionSystem);

  enum FlowSheetConfiguration { NoFlow, Line, Loop };

  struct TimeSection {
    std::string name;
    realtype t_duration;
    realtype pump_percentage;
    realtype mixing_percentage;
    FlowSheetConfiguration flow_sheet_configuration;
    RowVector inlet_solution;
    std::shared_ptr<Volume> fraction;
  };
  std::vector<TimeSection> recipe;

  recipe.push_back({"Load 1", 1199, 20, 0, Line, feed, frac_feed_1});
  recipe.push_back({"Load 1 pause", 34, 0, 0, NoFlow, buffer5_water, nullptr});
  recipe.push_back({"Load 2", 96, 20, 0, Line, feed, frac_feed_1});
  recipe.push_back({"Load 2 pause", 33, 0, 0, NoFlow, buffer5_water, nullptr});

  recipe.push_back({"Wash 1 fill", 50, 40, 0, Line, buffer1, frac_feed_2});
  recipe.push_back(
      {"Wash 1 resuspend", 65, 0, 60, NoFlow, buffer5_water, nullptr});
  recipe.push_back(
      {"Wash 1 resuspend loop", 51, 40, 50, Loop, buffer5_water, nullptr});
  recipe.push_back(
      {"Wash 1 recapture", 69, 0, 0, NoFlow, buffer5_water, nullptr});
  recipe.push_back(
      {"Wash 1 recapture loop", 196, 30, 0, Loop, buffer5_water, nullptr});
  recipe.push_back({"Wash 1 pause", 6, 0, 0, NoFlow, buffer5_water, nullptr});

  recipe.push_back({"Wash 2 fill", 50, 40, 0, Line, buffer1, frac_wash_1});
  recipe.push_back(
      {"Wash 2 resuspend", 65, 0, 60, NoFlow, buffer5_water, nullptr});
  recipe.push_back(
      {"Wash 2 resuspend loop", 51, 40, 50, Loop, buffer5_water, nullptr});
  recipe.push_back(
      {"Wash 2 recapture", 69, 0, 0, NoFlow, buffer5_water, nullptr});
  recipe.push_back(
      {"Wash 2 recapture loop", 196, 30, 0, Loop, buffer5_water, nullptr});
  recipe.push_back({"Wash 2 pause", 14, 0, 0, NoFlow, buffer5_water, nullptr});

  recipe.push_back({"Wash 3 fill", 50, 40, 0, Line, buffer2, frac_wash_2});
  recipe.push_back(
      {"Wash 3 resuspend", 65, 0, 60, NoFlow, buffer5_water, nullptr});
  recipe.push_back(
      {"Wash 3 resuspend loop", 51, 40, 50, Loop, buffer5_water, nullptr});
  recipe.push_back(
      {"Wash 3 recapture", 69, 0, 0, NoFlow, buffer5_water, nullptr});
  recipe.push_back(
      {"Wash 3 recapture loop", 196, 30, 0, Loop, buffer5_water, nullptr});
  recipe.push_back({"Wash 3 pause", 13, 0, 0, NoFlow, buffer5_water, nullptr});

  // 6s with Buffer2 in Matlab code, but not in the Appendix / Thesis

  recipe.push_back({"Elution 1 fill", 52, 40, 0, Line, buffer3, frac_wash_3});
  recipe.push_back(
      {"Elution 1 resuspend", 90, 0, 60, NoFlow, buffer5_water, nullptr});
  recipe.push_back(
      {"Elution 1 resuspend loop", 300, 30, 40, Loop, buffer5_water, nullptr});
  recipe.push_back(
      {"Elution 1 recapture", 60, 0, 0, NoFlow, buffer5_water, nullptr});
  recipe.push_back(
      {"Elution 1 recapture loop", 200, 30, 0, Loop, buffer5_water, nullptr});
  recipe.push_back(
      {"Elution 1 pause", 7, 0, 0, NoFlow, buffer5_water, nullptr});

  recipe.push_back(
      {"Elution 2 fill", 52, 40, 0, Line, buffer3, frac_elution_1});
  recipe.push_back(
      {"Elution 2 resuspend", 90, 0, 60, NoFlow, buffer5_water, nullptr});
  recipe.push_back(
      {"Elution 2 resuspend loop", 300, 30, 40, Loop, buffer5_water, nullptr});
  recipe.push_back(
      {"Elution 2 recapture", 60, 0, 0, NoFlow, buffer5_water, nullptr});
  recipe.push_back(
      {"Elution 2 recapture loop", 200, 30, 0, Loop, buffer5_water, nullptr});
  recipe.push_back(
      {"Elution 2 pause", 9, 0, 0, NoFlow, buffer5_water, nullptr});

  recipe.push_back(
      {"Elution 3 fill", 52, 40, 0, Line, buffer3, frac_elution_2});
  recipe.push_back(
      {"Elution 3 resuspend", 90, 0, 60, NoFlow, buffer5_water, nullptr});
  recipe.push_back(
      {"Elution 3 resuspend loop", 300, 30, 40, Loop, buffer5_water, nullptr});
  recipe.push_back(
      {"Elution 3 recapture", 60, 0, 0, NoFlow, buffer5_water, nullptr});
  recipe.push_back(
      {"Elution 3 recapture loop", 200, 30, 0, Loop, buffer5_water, nullptr});
  recipe.push_back(
      {"Elution 3 pause", 10, 0, 0, NoFlow, buffer5_water, nullptr});

  recipe.push_back(
      {"Elution 4 fill", 52, 40, 0, Line, buffer3, frac_elution_3});
  recipe.push_back(
      {"Elution 4 resuspend", 90, 0, 60, NoFlow, buffer5_water, nullptr});
  recipe.push_back(
      {"Elution 4 resuspend loop", 300, 30, 40, Loop, buffer5_water, nullptr});
  recipe.push_back(
      {"Elution 4 recapture", 60, 0, 0, NoFlow, buffer5_water, nullptr});
  recipe.push_back(
      {"Elution 4 recapture loop", 200, 30, 0, Loop, buffer5_water, nullptr});
  recipe.push_back(
      {"Elution 4 pause", 154, 0, 0, NoFlow, buffer5_water, nullptr});

  recipe.push_back(
      {"Elution 5 fill", 52, 40, 0, Line, buffer3, frac_elution_4});
  recipe.push_back(
      {"Elution 5 resuspend", 90, 0, 60, NoFlow, buffer5_water, nullptr});
  recipe.push_back(
      {"Elution 5 resuspend loop", 300, 30, 40, Loop, buffer5_water, nullptr});
  recipe.push_back(
      {"Elution 5 recapture", 60, 0, 0, NoFlow, buffer5_water, nullptr});
  recipe.push_back(
      {"Elution 5 recapture loop", 200, 30, 0, Loop, buffer5_water, nullptr});
  recipe.push_back(
      {"Elution 5 pause", 170, 0, 0, NoFlow, buffer5_water, nullptr});

  recipe.push_back({"Regeneration", 52, 25, 0, Line, buffer1, frac_elution_5});

  realtype total_time = 0.0;
  for (const auto &section : recipe) {
    std::ostringstream oss;
    oss << "From t=" << total_time
        << " to t=" << total_time + section.t_duration << " : " << section.name
        << " with " << section.pump_percentage << "% pump and "
        << section.mixing_percentage << "% mixing in configuration "
        << (section.flow_sheet_configuration == NoFlow
                ? "NoFlow"
                : (section.flow_sheet_configuration == Line ? "Line" : "Loop"))
        << "\n";
    LOG("recipe.log", oss.str());
    total_time += section.t_duration;
  }

  // =============================================================
  // ========================== DERIVED RECIPE FUNCTIONS =========
  // =============================================================

  ColVector pump_percentages(19);
  pump_percentages << 0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70,
      75, 80, 85, 90;
  ColVector upper_flowRate_curve(19);
  upper_flowRate_curve << 0, 91, 220, 330, 465, 621, 783, 916, 1075, 1236, 1420,
      1495, 1750, 1830, 1858, 1865, 1860, 1865, 1895;
  ColVector lower_flowRate_curve(19);
  lower_flowRate_curve << 0, 117, 260, 380, 505, 654, 783, 920, 1080, 1238,
      1424, 1501, 1758, 1838, 1861, 1875, 1868, 1875, 1910;
  ColVector flowRate_curve =
      (upper_flowRate_curve + lower_flowRate_curve) / 2.0 / 6e7;
  // Create the linear interpolator for pump percentages and flow rates
  LinearInterpolator<ColVector> pump_flowRate_Interpolator(pump_percentages,
                                                           flowRate_curve);

  auto flowRateFunction = [&](realtype t) -> realtype {
    realtype t_cache = 0.0;
    for (const auto &section : recipe) {
      t_cache += section.t_duration;
      if (t_cache > t) {
        return pump_flowRate_Interpolator(section.pump_percentage);
      }
    }
    return pump_flowRate_Interpolator(0.0);
    // throw std::runtime_error("Time " + std::to_string(t) + " exceeds total
    // duration of the recipe.");
  };

  realtype D_ax_pipes =
      0.0; // TODO: 0 vs 6*e-6? Both values are in Matlab and inconsistent
  ColVector mixing_percentages(5);
  mixing_percentages << 0, 20, 40, 60, 80;
  ColVector D_ax_curve(5);
  D_ax_curve << D_ax_pipes, 1.3e-5, 5e-5, 6.9e-5,
      8.4e-5; // TODO: inconsistent with 0.5e-6 for 60% in paper
  LinearInterpolator<ColVector> mixing_D_ax_Interpolator(mixing_percentages,
                                                         D_ax_curve);

  auto pc_D_ax_function = [&](realtype t) -> realtype {
    realtype t_cache = 0.0;
    for (const auto &section : recipe) {
      t_cache += section.t_duration;
      if (t_cache > t) {
        return mixing_D_ax_Interpolator(section.mixing_percentage);
      }
    }
    return mixing_D_ax_Interpolator(0.0);
    // throw std::runtime_error("Time " + std::to_string(t) + " exceeds total
    // duration of the recipe.");
  };

  // Capture function: 1.0 for magnetic capture, negative values for uncapture
  // i.e. dm_dt = capture_function(t) * m
  auto capture_function = [&](realtype t) -> realtype {
    realtype t_cache = 0.0;
    for (const auto &section : recipe) {
      t_cache += section.t_duration;
      if (t_cache > t) {
        if (section.mixing_percentage == 0.0) {
          return 1.0;
        } else {
          auto t_start = t_cache - section.t_duration;
          auto t_in_section = t - t_start;
          auto T_halfLife_uncapture = 5.0; // Half life of uncapture in seconds
          realtype uncapture_rate = std::log(2) / T_halfLife_uncapture;
          // To avoid sharp discontinuities, use a quadratic ramp over the first
          // 10s of the section
          if (t_in_section < 10.0) {
            return -(t_in_section / 10.0) * (t_in_section / 10.0) *
                   uncapture_rate;
          } else {
            return -uncapture_rate;
          }
        }
      }
    }
    return 1.0;
    // throw std::runtime_error("Time " + std::to_string(t) + " exceeds total
    // duration of the recipe.");
  };

  auto inlet_function = [&](realtype t) -> RowVector {
    realtype t_cache = 0.0;
    for (const auto &section : recipe) {
      t_cache += section.t_duration;
      if (t_cache > t) {
        return section.inlet_solution;
      }
    }
    return recipe.back().inlet_solution;
    // throw std::runtime_error("Time " + std::to_string(t) + " exceeds total
    // duration of the recipe.");
  };

  auto flowRate_function = [&](realtype t) -> realtype {
    realtype t_cache = 0.0;
    for (const auto &section : recipe) {
      t_cache += section.t_duration;
      if (t_cache > t) {
        if (section.flow_sheet_configuration == Loop ||
            section.flow_sheet_configuration == Line) {
          return pump_flowRate_Interpolator(section.pump_percentage);
        } else {
          assert(section.flow_sheet_configuration == NoFlow &&
                 section.pump_percentage == 0.0 &&
                 "If NoFlow, pump_percentage must be 0!");
          return 0.0;
        }
      }
    }
    return pump_flowRate_Interpolator(0.0);
    // throw std::runtime_error("Time " + std::to_string(t) + " exceeds total
    // duration of the recipe.");
  };

  auto flowRate_function_line = [&](realtype t) -> realtype {
    realtype t_cache = 0.0;
    for (const auto &section : recipe) {
      t_cache += section.t_duration;
      if (t_cache > t) {
        if (section.flow_sheet_configuration == Line) {
          return pump_flowRate_Interpolator(section.pump_percentage);
        } else {
          return 0.0;
        }
      }
    }
    return pump_flowRate_Interpolator(0.0);
    // throw std::runtime_error("Time " + std::to_string(t) + " exceeds total
    // duration of the recipe.");
  };

  auto flowRate_function_loop = [&](realtype t) -> realtype {
    realtype t_cache = 0.0;
    for (const auto &section : recipe) {
      t_cache += section.t_duration;
      if (t_cache > t) {
        if (section.flow_sheet_configuration == Loop) {
          return pump_flowRate_Interpolator(section.pump_percentage);
        } else {
          return 0.0;
        }
      }
    }
    return pump_flowRate_Interpolator(0.0);
    // throw std::runtime_error("Time " + std::to_string(t) + " exceeds total
    // duration of the recipe.");
  };

  auto flowRate_function_frac =
      [&](realtype t, const std::shared_ptr<Volume> &fraction) -> realtype {
    realtype t_cache = 0.0;
    for (const auto &section : recipe) {
      t_cache += section.t_duration;
      if (t_cache > t) {
        if (section.fraction == fraction) {
          return pump_flowRate_Interpolator(section.pump_percentage);
        } else {
          return 0.0;
        }
      }
    }
    return pump_flowRate_Interpolator(0.0);
    // throw std::runtime_error("Time " + std::to_string(t) + " exceeds total
    // duration of the recipe.");
  };

  // =============================================================
  // ========================== UNIT OPERATIONS ==================
  // =============================================================

  auto inlet = std::make_shared<Inlet>(reactionSystem, inlet_function);

  realtype d_pipes = 0.0096;
  realtype A_cross_pipes = M_PI / 4 * d_pipes * d_pipes;

  realtype l_inlet_pipe = 0.67 + (0.00005 / A_cross_pipes);

  auto scale_n = [&](int base) {
    return std::max(1,
                    static_cast<int>(std::round(base * discretization_factor)));
  };
  auto pipe_inlet =
      std::make_shared<Pipe>(reactionSystem, scale_n(10), A_cross_pipes,
                             l_inlet_pipe, flowRateFunction, D_ax_pipes);
  pipe_inlet->setConstInitialConcentration(buffer5_water);

  auto a_eff_function = [](realtype um_u0_ratio) {
    auto a1 = 2.035;
    auto a2 = 107.1;
    auto a3 = -0.00808;
    auto p = 0.07477;
    auto q = 1.083;
    return 1 / (1 + a1 * std::exp(-p * um_u0_ratio) +
                a2 * std::exp(-q * um_u0_ratio) + a3);
  };

  realtype l_pc = 9.755e-02;
  realtype d_pc = 9.031e-2;
  realtype d_shaft_pc = 6.691e-2;
  realtype A_cross_pc = (M_PI / 4) * (d_pc * d_pc - d_shaft_pc * d_shaft_pc);
  realtype V_int_pc = 2.5558e-4;
  realtype porosity_pc = V_int_pc / (A_cross_pc * l_pc);
  realtype magnetic_field_strength = 2.227 * 1e5;
  realtype MNP_capacity = 40e-3;

  // RS_MagneticCaptureProcessChamber with all parameters
  auto pc = std::make_shared<RS_MagneticCaptureProcessChamber>(
      reactionSystem,
      scale_n(30), // n_cells (scaled)
      A_cross_pc,  // crossSectionArea - Cross-sectional area of the chamber
      l_pc,        // length - Height of the discs in the chamber
      porosity_pc, // empty_porosity - Porosity of the chamber when empty
      25,          // n_disks - Number of disks in the chamber
      8.0e-4,      // disk_height - Height of the discs in the chamber
      19, // alpha_fluid_particle_volume_ratio - Slurry factor ratio of fluid
          // volume to MNP volume in the slurry
      2, // deposition_rate - Empirical parameter in G to account for deposition
         // rate
      MNP_capacity, // MNP_capacity - Maximum mass of magnetic particles with
                    // that can be captured
      magnetic_field_strength, // T_halfLife_uncapture - In the uncapture case,
                               // half of the MNPs in slurry phase decay back to
                               // liquid phase every T_halfLife seconds
      a_eff_function,          // effective capture area function
      capture_function,        // magnetic_field_strength_function
      flowRateFunction,        // flowRateFunction
      pc_D_ax_function);       // dispersion_coefficient_function
  pc->setConstInitialConcentration(buffer5_water);

  realtype V_inletRegion_pc = 4.113e-5;
  realtype V_outletRegion_pc = 8.381e-5;
  realtype V_dead = V_inletRegion_pc + V_outletRegion_pc;

  auto dead_volume = std::make_shared<Volume>(reactionSystem, V_dead);
  dead_volume->setConstInitialConcentration(buffer5_water);

  realtype l_pipe_outlet = 1.02;
  auto pipe_outlet =
      std::make_shared<Pipe>(reactionSystem, scale_n(10), A_cross_pipes,
                             l_pipe_outlet, flowRateFunction, D_ax_pipes);
  pipe_outlet->setConstInitialConcentration(buffer5_water);

  realtype l_pipe_loop = 1.73;
  auto pipe_loop =
      std::make_shared<Pipe>(reactionSystem, scale_n(10), A_cross_pipes,
                             l_pipe_loop, flowRate_function_loop, D_ax_pipes);
  pipe_loop->setConstInitialConcentration(buffer5_water);

  // Calculate total duration
  realtype total_duration = 0.0;
  for (const auto &section : recipe) {
    total_duration += section.t_duration;
  }
  std::cout << "Total duration of the process: " << total_duration
            << " seconds." << std::endl;

  Process process(cs, {inlet, pipe_inlet, pc, dead_volume, pipe_outlet,
                       pipe_loop, frac_feed_1, frac_feed_2, frac_wash_1,
                       frac_wash_2, frac_wash_3, frac_elution_1, frac_elution_2,
                       frac_elution_3, frac_elution_4,
                       frac_elution_5}); // TODO: t_end remove, make functions
                                         // robust to t>total_duration

  // =============================================================
  // ========================== CONNECTIONS ======================
  // =============================================================
  // WARNING: You need to check consistency of connections and flow rates for
  // all times yourself!

  // Inlet / Line connections
  process.addConnection(inlet->out(), pipe_inlet->in(), flowRate_function_line);

  // Always connected
  process.addConnection(pipe_inlet->out(), pc->in(), flowRateFunction);
  process.addConnection(pc->out(), dead_volume->in(), flowRateFunction);
  process.addConnection(dead_volume->out(), pipe_outlet->in(),
                        flowRateFunction);
  // process.addConnection(pipe_outlet->out(), outlet->in(),
  // flowRate_function_line);

  // Connections to fractions
  process.addConnection(pipe_outlet->out(), frac_feed_1->in(), [&](realtype t) {
    return flowRate_function_frac(t, frac_feed_1);
  });
  process.addConnection(pipe_outlet->out(), frac_feed_2->in(), [&](realtype t) {
    return flowRate_function_frac(t, frac_feed_2);
  });
  process.addConnection(pipe_outlet->out(), frac_wash_1->in(), [&](realtype t) {
    return flowRate_function_frac(t, frac_wash_1);
  });
  process.addConnection(pipe_outlet->out(), frac_wash_2->in(), [&](realtype t) {
    return flowRate_function_frac(t, frac_wash_2);
  });
  process.addConnection(pipe_outlet->out(), frac_wash_3->in(), [&](realtype t) {
    return flowRate_function_frac(t, frac_wash_3);
  });
  process.addConnection(
      pipe_outlet->out(), frac_elution_1->in(),
      [&](realtype t) { return flowRate_function_frac(t, frac_elution_1); });
  process.addConnection(
      pipe_outlet->out(), frac_elution_2->in(),
      [&](realtype t) { return flowRate_function_frac(t, frac_elution_2); });
  process.addConnection(
      pipe_outlet->out(), frac_elution_3->in(),
      [&](realtype t) { return flowRate_function_frac(t, frac_elution_3); });
  process.addConnection(
      pipe_outlet->out(), frac_elution_4->in(),
      [&](realtype t) { return flowRate_function_frac(t, frac_elution_4); });
  process.addConnection(
      pipe_outlet->out(), frac_elution_5->in(),
      [&](realtype t) { return flowRate_function_frac(t, frac_elution_5); });

  // Loop Connections
  process.addConnection(pipe_outlet->out(), pipe_loop->in(),
                        flowRate_function_loop);
  process.addConnection(pipe_loop->out(), pipe_inlet->in(),
                        flowRate_function_loop);

  // =============================================================
  // ========================== OBSERVER =========================
  // =============================================================

  // Compute maximum flow rate and minimum cell_volume
  // realtype max_flow_rate = pump_flowRate_Interpolator(40);  // Max at 40%
  // pump percentage realtype min_cell_volume =
  // std::min({pipe_inlet->cell_volume, pc->cell_volume_l_and_sl,
  // dead_volume->cell_volume,
  //                                      pipe_outlet->cell_volume,
  //                                      pipe_loop->cell_volume});
  // std::cout << "Max flow rate: " << max_flow_rate << " m³/s" << std::endl;
  // std::cout << "Min cell volume: " << min_cell_volume << " m³" << std::endl;
  // std::cout << "Cell volume l and sl pc: " << pc->cell_volume_l_and_sl << "
  // m³" << std::endl;

  Solver solver(process, solverType);

  std::size_t dt_obs = 2.0; // Observe every dt_obs seconds
  std::size_t n_obs_time_steps =
      static_cast<std::size_t>(total_duration / dt_obs) + 1;
  std::cout << "Number of observation time steps: " << n_obs_time_steps
            << std::endl;

  bool compute_errors = true;
  TimeSeriesObserver obs_pipe_inlet(0, total_duration, n_obs_time_steps,
                                    pipe_inlet->all(), solver, save_obs,
                                    compute_errors);
  TimeSeriesObserver obs_pc_liquid(0, total_duration, n_obs_time_steps,
                                   pc->liquid(), solver, save_obs,
                                   compute_errors);
  TimeSeriesObserver obs_pc_slurry(0, total_duration, n_obs_time_steps,
                                   pc->slurry(), solver, save_obs,
                                   compute_errors);
  TimeSeriesObserver obs_pipe_oulet(0, total_duration, n_obs_time_steps,
                                    pipe_outlet->all(), solver, save_obs,
                                    compute_errors);
  TimeSeriesObserver obs_pipe_loop(0, total_duration, n_obs_time_steps,
                                   pipe_loop->all(), solver, save_obs,
                                   compute_errors);
  TimeSeriesObserver obs_frac_feed_1(0, total_duration, n_obs_time_steps,
                                     frac_feed_1->all(), solver, save_obs,
                                     compute_errors);
  TimeSeriesObserver obs_frac_feed_2(0, total_duration, n_obs_time_steps,
                                     frac_feed_2->all(), solver, save_obs,
                                     compute_errors);
  TimeSeriesObserver obs_frac_wash_1(0, total_duration, n_obs_time_steps,
                                     frac_wash_1->all(), solver, save_obs,
                                     compute_errors);
  TimeSeriesObserver obs_frac_wash_2(0, total_duration, n_obs_time_steps,
                                     frac_wash_2->all(), solver, save_obs,
                                     compute_errors);
  TimeSeriesObserver obs_frac_wash_3(0, total_duration, n_obs_time_steps,
                                     frac_wash_3->all(), solver, save_obs,
                                     compute_errors);
  TimeSeriesObserver obs_frac_elution_1(0, total_duration, n_obs_time_steps,
                                        frac_elution_1->all(), solver, save_obs,
                                        compute_errors);
  TimeSeriesObserver obs_frac_elution_2(0, total_duration, n_obs_time_steps,
                                        frac_elution_2->all(), solver, save_obs,
                                        compute_errors);
  TimeSeriesObserver obs_frac_elution_3(0, total_duration, n_obs_time_steps,
                                        frac_elution_3->all(), solver, save_obs,
                                        compute_errors);
  TimeSeriesObserver obs_frac_elution_4(0, total_duration, n_obs_time_steps,
                                        frac_elution_4->all(), solver, save_obs,
                                        compute_errors);
  TimeSeriesObserver obs_frac_elution_5(0, total_duration, n_obs_time_steps,
                                        frac_elution_5->all(), solver, save_obs,
                                        compute_errors);

  TimeSeriesObserver obs_complete_y_vector(
      0, total_duration, n_obs_time_steps,
      ArrayMapper(inlet.get(), solver.getYSize() / cs.n_components,
                  cs.n_components),
      solver, save_obs);
  TimeSeriesObserver obs_pc_outlet(0, total_duration, n_obs_time_steps,
                                   pc->out(), solver, save_obs);

  // To observe:
  // Elution Protein
  // ph/Salts

  // To optimize:
  // Mass MP
  // Elution steps

  // =============================================================
  // ========================== SOLVING ==========================
  // =============================================================

  solver.solve(total_duration, timeout_seconds);
  realtype t_solve_duration = solver.getSolveTime();

  // =============================================================
  // ========================== RESULTS ==========================
  // =============================================================

  std::cout << "Max pH error: " << reactionSystem.max_error << std::endl;
  const auto &internal_time_stamps = solver.getInternalTimeStamps();
  // Save internal time stamps to .npy file
  auto obs_dir = logger::obs_dir();
  FS3::npy_save(obs_dir + "/internal_timestamps.npy", &internal_time_stamps[0],
                {internal_time_stamps.size()}, "w");

  if (save_obs) {
    auto obs_unitOperations_filename = obs_dir + "/unit_operations.npz";
    obs_pipe_inlet.save_to_npz(obs_unitOperations_filename, "pipe_inlet", "a");
    obs_pc_liquid.save_to_npz(obs_unitOperations_filename, "pc_liquid", "a");
    obs_pc_slurry.save_to_npz(obs_unitOperations_filename, "pc_slurry", "a");
    obs_pipe_oulet.save_to_npz(obs_unitOperations_filename, "pipe_outlet", "a");
    obs_pipe_loop.save_to_npz(obs_unitOperations_filename, "pipe_loop", "a");
    obs_frac_feed_1.save_to_npz(obs_unitOperations_filename, "frac_feed_1",
                                "a");
    obs_frac_feed_2.save_to_npz(obs_unitOperations_filename, "frac_feed_2",
                                "a");
    obs_frac_wash_1.save_to_npz(obs_unitOperations_filename, "frac_wash_1",
                                "a");
    obs_frac_wash_2.save_to_npz(obs_unitOperations_filename, "frac_wash_2",
                                "a");
    obs_frac_wash_3.save_to_npz(obs_unitOperations_filename, "frac_wash_3",
                                "a");
    obs_frac_elution_1.save_to_npz(obs_unitOperations_filename,
                                   "frac_elution_1", "a");
    obs_frac_elution_2.save_to_npz(obs_unitOperations_filename,
                                   "frac_elution_2", "a");
    obs_frac_elution_3.save_to_npz(obs_unitOperations_filename,
                                   "frac_elution_3", "a");
    obs_frac_elution_4.save_to_npz(obs_unitOperations_filename,
                                   "frac_elution_4", "a");
    obs_frac_elution_5.save_to_npz(obs_unitOperations_filename,
                                   "frac_elution_5", "a");

    auto obs_complete_y_vector_filename = obs_dir + "/complete_y_vector.npy";
    obs_complete_y_vector.save_to_npy(obs_complete_y_vector_filename);
    auto obs_pc_outlet_filename = obs_dir + "/pc_outlet.npy";
    obs_pc_outlet.save_to_npy(obs_pc_outlet_filename);

    // Convert pipe outlet data to pH values for each cell and save
    std::size_t cell_index = pipe_outlet->n_cells / 2; // middle cell
    ColVector pH_activity_pipe_outlet =
        convert_to_pH(obs_pipe_oulet, cell_index, pipe_outlet->cell_volume, cs,
                      acitivityModel);
    FS3::npy_save(obs_dir + "/pipe_outlet_middle_cell_pH_activity.npy",
                  pH_activity_pipe_outlet.data(),
                  {static_cast<size_t>(pH_activity_pipe_outlet.size())}, "w");
    ColVector pH_concentration_pipe_outlet =
        convert_to_pH(obs_pipe_oulet, cell_index, pipe_outlet->cell_volume, cs,
                      NoActivityModel(cs));
    FS3::npy_save(obs_dir + "/pipe_outlet_middle_cell_pH_concentration.npy",
                  pH_concentration_pipe_outlet.data(),
                  {static_cast<size_t>(pH_concentration_pipe_outlet.size())},
                  "w");
  }

  std::cout << "Solved with " << internal_time_stamps.size()
            << " internal time stamps in " << t_solve_duration << " seconds."
            << std::endl;

  if (t_solve_duration > timeout_seconds) {
    std::cout << "WARNING: Solver timed out after " << timeout_seconds
              << " seconds at t=" << solver.getT() << " / " << total_duration
              << std::endl;
  }

  return std::make_tuple(reactionSystem.max_error, t_solve_duration,
                         internal_time_stamps.size());
}
