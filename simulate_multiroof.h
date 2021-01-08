//
// Created by Brad Huang on 8/16/20.
//

#ifndef ROBUST_SIZING_SIMULATE_MULTIROOF_H
#define ROBUST_SIZING_SIMULATE_MULTIROOF_H

#include <vector>
#include "params_multiroof.h"

using namespace std;

double static num_cells = 200.0; // just a default value that will be updated every time we check a new battery size
double static nominal_voltage_c = 3.8793;
double static nominal_voltage_d = 3.5967;
double static a1_slope = 0.1920;
double static a2_slope = -0.4865;
double static a1_intercept = 0.0*num_cells;
double static a2_intercept = kWh_in_one_cell*num_cells;
double static eta_d = 1/0.9; // taking reciprocal so that we don't divide by eta_d when updating the battery energy content
double static eta_c = 0.9942;
double static alpha_d = a2_intercept*1.0; // the 1 indicates the maximum discharging C-rate
double static alpha_c = a2_intercept*1.0; // the 1 indicates the maximum charging C-rate

double static default_adagrad_step_size = 1;
double static cost_threshold = 50;
size_t static mean_window_size = 200;
size_t static adagrad_max_it = 20000;

extern size_t total_sim_called;

double sim(vector <double> &load_trace, vector<vector<double>> &solar_traces,
           size_t start_index, size_t end_index, double cells, valarray<double> pvs, double b_0);

vector<SimulationMultiRoofResult> simulate(
        vector<double> &load_trace, vector<vector<double>> &solar_traces,
        size_t start_index, size_t end_index, double b_0);

SimulationMultiRoofResult binary_search_result(
        vector<double> &load_trace, vector<vector<double>> &solar_traces,
        size_t start_index, size_t end_index, valarray<double> &pvs, double b_0 = 0,
        double cells_U = cells_max, double cells_L = cells_min);

vector<SimulationMultiRoofResult> simulate_deterministic_adagrad(
        vector<double> &load_trace, vector<vector<double>> &solar_traces, size_t start_index, size_t end_index,
        valarray<double> &init_pv, double b_0 = 0, double step_size = default_adagrad_step_size, double fudge_factor = 1e-6);

double random_simulate_cheroff(
        vector <double> &load_trace, vector<vector<double>> &solar_traces,
        double cells, valarray<double> &pvs, double b_0 = 0);

double deterministic_simulate_cheroff(
        vector <double> &load_trace, vector<vector<double>> &solar_traces,
        double cells, valarray<double> &pvs, double b_0 = 0);

double get_cheroff(double p_tilde);

double chernoff_result(
        vector <double> &load_trace, vector<vector<double>> &solar_traces,
        double cells, valarray<double> &pvs, double b_0 = 0);

vector<SimulationMultiRoofResult>
tabu_cheroff(vector<double> &load_trace, vector<vector<double>> &solar_traces, double target_p, double b_0 = 0);

#endif //ROBUST_SIZING_SIMULATE_MULTIROOF_H
