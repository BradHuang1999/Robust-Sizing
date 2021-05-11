//
// Created by Brad Huang on 8/16/20.
//

#ifndef ROBUST_SIZING_CHEBY_MULTIROOF_H
#define ROBUST_SIZING_CHEBY_MULTIROOF_H

#include "params_multiroof.h"
#include <stdexcept>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

// VARIABLES

extern double lambda2;

// CONSTANTS

/**
 * lambda_factor: this is how much excess lambda^2 we tolerate
 */
constexpr double static CHEBYSHEV_BETA = 0.1;

/**
 * cheby_num_steps: number of steps in each dimension to search for
 */
constexpr size_t static cheby_num_steps = 1000;

// FUNCTIONS

/**
 * update_chebyshev_params: update before chebyshev simulation
 */
void update_chebyshev_params(const valarray<bool> &is_zeros);

/**
 * convert_simulation_result_to_matrix: convert vector of SimulationMultiRoofResult's
 *   to Eigen matrix. can be refactored to other modules
 */
MatrixXd convert_simulation_result_to_matrix(
        const vector<SimulationMultiRoofResult> &adagrad_sims,
        const valarray<bool> &is_zeros, bool normalize_battery = true);

/**
 * get_chebyshev_bound: calculates chebyshev bound based on the above variables
 * @param adagrad_sims adagra simulation results (min_cost points)
 * @return non-dominated, pareto efficient, chebyshev bound
 */
vector<SimulationMultiRoofResult> get_chebyshev_bound(
        const vector<SimulationMultiRoofResult> &adagrad_sims,
        const valarray<bool> &is_zeros);

#endif //ROBUST_SIZING_CHEBY_MULTIROOF_H
