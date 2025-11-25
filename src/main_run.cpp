#include <sundials/sundials_types.h>

#include <Solver.hpp>
#include <limits>

#include "HGMS_hIgG_separation_shared.hpp"

int main() {
    run_HGMS_hIgG_separation(1e3, 0.1, true, std::numeric_limits<realtype>::infinity(), SolverType::ERK, 1.0);
    return 0;
}