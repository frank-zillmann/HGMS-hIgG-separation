#ifndef HGMS_HIGG_SEPARATION_SHARED_HPP
#define HGMS_HIGG_SEPARATION_SHARED_HPP

#include <limits>
#include <tuple>

#include "Solver.hpp"
#include "sundials/sundials_types.h"

// Extended signature: add discretization_factor (scales n_cells of unit operations). Default keeps behavior unchanged.
std::tuple<double, double, size_t> run_HGMS_hIgG_separation(realtype kf_ion = 1e10,
                                                            realtype tau_reaction = 1,
                                                            bool save_obs = true,
                                                            realtype timeout_seconds = std::numeric_limits<realtype>::infinity(),
                                                            SolverType solverType = SolverType::ERK,
                                                            realtype discretization_factor = 1.0);

#endif  // HGMS_HIGG_SEPARATION_SHARED_HPP