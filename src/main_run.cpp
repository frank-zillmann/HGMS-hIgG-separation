#include <sundials/sundials_types.h>

#include <Solver.hpp>
#include <limits>

#include "HGMS_hIgG_separation_shared.hpp"

int main() {
  run_HGMS_hIgG_separation(1e3, 0.1, true,
                           std::numeric_limits<realtype>::infinity(),
                           SolverType::ERK, 0.1);
  return 0;
}