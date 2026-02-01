#include <sundials/sundials_types.h>

#include <chrono>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

#include "HGMS_hIgG_separation_shared.hpp"
#include "Observers/NumpyIO.hpp"

int main() {
  // Define the hyperparameter grid

  // for jacobian approx method
  // std::vector<realtype> kf_ion_values = {1e1, 1e10};
  // std::vector<realtype> tau_reaction_values = {1e-1, 1e0, 1e1, 1e2, 1e3,
  // 1e4};

  // For this search we vary the discretization factor (scales n_cells) instead
  // of k_f, and keep the tau values as before.
  std::vector<realtype> discretization_factors = {0.5, 1.0, 2.0};
  std::vector<realtype> tau_reaction_values = {0.01, 0.1, 1.0};

  size_t n_disc = discretization_factors.size();
  size_t n_tau = tau_reaction_values.size();

  // Timeout in seconds for each simulation
  const realtype timeout_seconds = 600.0;

  // Prepare arrays to store results
  std::vector<realtype> max_errors;
  std::vector<realtype> solve_times;
  std::vector<realtype> num_steps;

  max_errors.reserve(n_disc * n_tau);
  solve_times.reserve(n_disc * n_tau);
  num_steps.reserve(n_disc * n_tau);

  std::cout << "Starting hyperparameter grid search..." << std::endl;
  std::cout << "Grid size: " << n_disc << " x " << n_tau << " = "
            << n_disc * n_tau << " runs" << std::endl;
  std::cout << "Timeout per run: " << timeout_seconds << " seconds ("
            << timeout_seconds / 60.0 << " minutes)" << std::endl;

  int run_counter = 0;
  auto total_start = std::chrono::high_resolution_clock::now();

  // Grid search: iterate through all parameter combinations
  for (size_t i = 0; i < n_disc; ++i) {
    for (size_t j = 0; j < n_tau; ++j) {
      run_counter++;
      realtype discretization_factor = discretization_factors[i];
      realtype tau_reaction = tau_reaction_values[j];

      std::cout << "\n========================================" << std::endl;
      std::cout << "Run " << run_counter << "/" << n_disc * n_tau << std::endl;
      std::cout << "discretization_factor = " << discretization_factor
                << ", tau_reaction = " << tau_reaction << std::endl;
      std::cout << "========================================" << std::endl;

      try {
        // Run simulation with timeout built into the solver
        const realtype kf_ion_fixed =
            1e10; // keep k_f fixed; inverse method does not use it
        auto [max_error, t_solve, internal_steps] = run_HGMS_hIgG_separation(
            kf_ion_fixed, tau_reaction, false, timeout_seconds,
            SolverType::ADAMS, discretization_factor);

        std::cout << "Results: max_error = " << max_error
                  << ", t_solve = " << t_solve
                  << " s, internal_steps = " << internal_steps << std::endl;

        max_errors.push_back(max_error);
        solve_times.push_back(t_solve);
        num_steps.push_back(static_cast<realtype>(internal_steps));

      } catch (const std::exception &e) {
        // Handle simulation failures
        std::cerr << "Simulation failed with exception: " << e.what()
                  << std::endl;
        max_errors.push_back(std::numeric_limits<realtype>::quiet_NaN());
        solve_times.push_back(std::numeric_limits<realtype>::quiet_NaN());
        num_steps.push_back(std::numeric_limits<realtype>::quiet_NaN());
      }
    }
  }

  auto total_end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<realtype> total_time = total_end - total_start;

  std::cout << "\n========================================" << std::endl;
  std::cout << "Grid search completed!" << std::endl;
  std::cout << "Total time: " << total_time.count() << " seconds" << std::endl;
  std::cout << "========================================" << std::endl;

  // Generate timestamp for unique filename
  auto now = std::chrono::system_clock::now();
  auto now_c = std::chrono::system_clock::to_time_t(now);
  std::tm now_tm = *std::localtime(&now_c);

  std::ostringstream filename_stream;
  filename_stream << "hyperparameter_search_results_"
                  << std::put_time(&now_tm, "%d-%m-%Y_%H-%M") << ".npz";
  std::string filename = filename_stream.str();

  // Save results to npz file
  std::cout << "Saving results to " << filename << "..." << std::endl;

  // Save arrays - results are already in row-major order (iterate
  // discretization_factors outer, tau_reaction inner)
  FS3::npz_save(filename, "max_errors", max_errors.data(), {n_disc, n_tau},
                "w");
  FS3::npz_save(filename, "solve_times", solve_times.data(), {n_disc, n_tau},
                "a");
  FS3::npz_save(filename, "num_steps", num_steps.data(), {n_disc, n_tau}, "a");
  FS3::npz_save(filename, "discretization_factors",
                discretization_factors.data(), {n_disc}, "a");
  FS3::npz_save(filename, "tau_reaction_values", tau_reaction_values.data(),
                {n_tau}, "a");

  std::cout << "Results saved successfully!" << std::endl;
  std::cout << "\nSummary:" << std::endl;
  std::cout << "- max_errors: shape (" << n_disc << ", " << n_tau
            << ") - 2D array" << std::endl;
  std::cout << "- solve_times: shape (" << n_disc << ", " << n_tau
            << ") - 2D array" << std::endl;
  std::cout << "- num_steps: shape (" << n_disc << ", " << n_tau
            << ") - 2D array" << std::endl;
  std::cout << "- discretization_factors: shape (" << n_disc
            << ") - values for axis 0" << std::endl;
  std::cout << "- tau_reaction_values: shape (" << n_tau
            << ") - values for axis 1" << std::endl;

  return 0;
}
