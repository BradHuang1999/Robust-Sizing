//
// Created by Brad Huang on 8/16/20.
//

#include <cmath>
#include <random>
#include <deque>
#include <unordered_map>
#include <valarray>

#include "simulate_multiroof.h"

using namespace std;

size_t total_sim_called = 0;

// parameters specified for an NMC cell with operating range of 1 C charging and discharging

void update_parameters(double n) {
    num_cells = n;

    a1_intercept = 0.0*num_cells;
    a2_intercept = kWh_in_one_cell*num_cells;

    alpha_d = a2_intercept*1.0;
    alpha_c = a2_intercept*1.0;
}

// decrease the applied (charging) power by increments of (1/30) until the power is
// low enough to avoid violating the upper energy limit constraint.
double calc_max_charging(double power, double b_prev) {

    double step = power/30.0;

    for (double c = power; c >= 0; c -= step) {
        double upper_lim = a2_slope*(c/nominal_voltage_c) + a2_intercept;
        double b = b_prev + c*eta_c*T_u;
        if (b <= upper_lim) {
            return c;
        }
    }
    return 0;
}


// decrease the applied (discharging) power by increments of (1/30) until the power is
// low enough to avoid violating the lower energy limit constraint.
double calc_max_discharging(double power, double b_prev) {

    double step = power/30.0;

    for (double d = power; d >= 0; d -= step) {
        double lower_lim = a1_slope*(d/nominal_voltage_d) + a1_intercept;
        double b = b_prev - d*eta_d*T_u;
        if (b >= lower_lim) {
            return d;
        }
    }
    return 0;
}


// Note: sim_year calls procedures calc_max_charging and calc_max_discharging.
// You could potentially speed up the computation by expanding these functions into sim_year
// to avoid procedure calls in this inner loop.
double sim(vector <double> &load_trace, vector<vector<double>> &solar_traces,
           size_t start_index, size_t end_index, double cells, valarray<double> pvs, double b_0)
{
    ++total_sim_called;

    update_parameters(cells);

    // set the battery
    double b = b_0*cells*kWh_in_one_cell; //0.5*a2_intercept

    int loss_events = 0;

    double load_deficit = 0;
    double load_sum = 0;

    size_t trace_length_load = load_trace.size();

    for (size_t t = start_index; t < end_index; t++)
    {
        // wrap around to the start of the trace if we hit the end.
        const double total_load = load_trace[t % trace_length_load];

        double total_solar = 0;
        for (size_t i = 0; i < n_solars; ++i) {
            vector<double>& solar_trace = solar_traces[i];
            size_t t_index = t % solar_trace.size();
            total_solar += solar_trace[t_index] * pvs[i];
        }

        load_sum += total_load;

        // first, calculate how much power is available for charging, and how much is needed to discharge
        if (total_solar > total_load) {
            const double c = total_solar - total_load;
            const double max_c = fmin(calc_max_charging(c,b), alpha_c);
            b = b + max_c*eta_c*T_u;
        } else {
            const double d = total_load - total_solar;
            const double max_d = fmin(calc_max_discharging(d,b), alpha_d);
            b = b - max_d*eta_d*T_u;

            // if we didnt get to discharge as much as we wanted, there is a loss
            if (max_d < d) {
                loss_events += 1;
                load_deficit += (d - max_d);
            }
        }
    }

    if (metric == 0) {
        // lolp
        return loss_events / (double)(end_index - start_index);
    } else {
        // metric == 1, eue
        return load_deficit / (double)load_sum;
    }
}

// Run simulation for provides solar and load trace to find cheapest combination of
// load and solar that can meet the epsilon target
vector<SimulationMultiRoofResult> simulate(vector<double> &load_trace, vector<vector<double>> &solar_traces,
                                           size_t start_index, size_t end_index, double b_0)
{
    vector<SimulationMultiRoofResult> pareto_set;

    valarray<double> pv_values(pv_maxs);
    vector<size_t> pv_nonzero_indexes;
    for (size_t i = 0; i < n_solars; ++i) {
        if (pv_values[i] > 0) {
            pv_nonzero_indexes.push_back(i);
        }
    }

    // first, find the lowest value of cells that will get us epsilon loss when the PV is maximized
    // use binary search
    double cells_U = cells_max;
    double cells_L = cells_min;

    while (cells_U - cells_L > cells_step)
    {
        double mid_cells = (cells_L + cells_U) / 2.0;
        double loss = sim(load_trace, solar_traces, start_index, end_index, mid_cells, pv_values, b_0);

        //cout << "sim result with " << a2_intercept << " kWh and " << pv_max << " pv: " << loss << endl;
        if (loss > epsilon)
        {
            cells_L = mid_cells;
        }
        else
        {
            // (loss <= epsilon)
            cells_U = mid_cells;
        }
    }

    default_random_engine ran_gen;
    uniform_int_distribution<size_t> rdm(0, pv_nonzero_indexes.size() - 1);

    double b_value;
    double b_init = cells_U;

    while (!pv_nonzero_indexes.empty())
    {
        bool found = false;
        for (b_value = b_init; b_value < cells_max && !found; b_value += cells_step)
        {
            double loss = sim(load_trace, solar_traces, start_index, end_index, b_value, pv_values, b_0);
            if (loss < epsilon) {
                found = true;
                pareto_set.emplace_back(b_value, pv_values, start_index);
                b_init = b_value;
            }
        }

        if (b_value >= cells_max) {
            return pareto_set;
        } else {
            size_t rindex_a_nonzero = rdm(ran_gen);
            size_t reduce_index = pv_nonzero_indexes[rindex_a_nonzero];

            pv_values[reduce_index] -= pv_steps[reduce_index];
            if (pv_values[reduce_index] <= 0) {
                // erase reduce_index, at rindex_a_nonzero, from pv_nonzero_indexes
                pv_nonzero_indexes.erase(pv_nonzero_indexes.begin() + rindex_a_nonzero);
                // replace rdm
                rdm = uniform_int_distribution<size_t>(0, pv_nonzero_indexes.size() - 1);
            }
        }
    }

    return pareto_set;
}

SimulationMultiRoofResult
binary_search_result(vector<double> &load_trace, vector<vector<double>> &solar_traces, size_t start_index,
                     size_t end_index, valarray<double> &pvs, double b_0, double cells_U, double cells_L)
{
    double loss_U = INFTY;

    while (cells_U - cells_L > cells_step)
    {
        double mid_cells = (cells_L + cells_U) / 2.0;
        double loss = sim(load_trace, solar_traces, start_index, end_index, mid_cells, pvs, 0);

        if (loss > epsilon)
        {
            cells_L = mid_cells;
        }
        else
        {
            // (loss <= epsilon)
            loss_U = loss;
            cells_U = mid_cells;
        }
    }

    if (loss_U <= epsilon) {
        return SimulationMultiRoofResult(cells_U, pvs, start_index);
    } else {
        return SimulationMultiRoofResult(start_index);
    }
}

vector<SimulationMultiRoofResult> simulate_deterministic_adagrad(
        vector<double> &load_trace, vector<vector<double>> &solar_traces, size_t start_index, size_t end_index,
        valarray<double> &init_pv, double b_0, double step_size, double fudge_factor)
{
    size_t pv_n = init_pv.size();
    valarray<bool> is_zero = init_pv == 0;

    default_random_engine generator;

    vector<normal_distribution<double>> distribution;
    for (size_t t = 0; t < pv_n; ++t) {
        distribution.emplace_back(0, pv_steps[t]);
    }

    valarray<double> pv_values(init_pv);
    valarray<double> g_ti(0.0, pv_n);

    deque<double> mean_window;
    double window_sum;

    vector<SimulationMultiRoofResult> ret;

    SimulationMultiRoofResult last_run_result = binary_search_result(
            load_trace, solar_traces, start_index, end_index, pv_values, b_0);
    if (!last_run_result.feasible) {
        return ret;
    }

    for (size_t it_count = 0; it_count < adagrad_max_it; ++it_count) {
        ret.push_back(last_run_result);

        bool has_feasible = false;
        valarray<double> grad(0.0, pv_n);

        for (size_t t = 0; t < pv_n; ++t) {
            if (!is_zero[t]) {
                pv_values[t] += pv_steps[t];
                SimulationMultiRoofResult next_step_result = binary_search_result(
                        load_trace, solar_traces, start_index, end_index, pv_values, b_0);

                if (next_step_result.feasible) {
                    has_feasible = true;
                    grad[t] = (next_step_result.cost - last_run_result.cost) / pv_steps[t];
                }
                pv_values[t] -= pv_steps[t];
            }
        }

        if (!has_feasible) {
            return ret;
        }

        g_ti += grad * grad;
        pv_values -= step_size * (grad / (fudge_factor + sqrt(g_ti)));

        for (size_t t = 0; t < pv_n; ++t) {
            if (!is_zero[t]) {
                pv_values[t] += distribution[t](generator);
            }
        }

        pv_values[(pv_values != 0) && (pv_values < pv_steps)] = pv_steps[(pv_values != 0) && (pv_values < pv_steps)];

        SimulationMultiRoofResult curr_run_result = binary_search_result(
                load_trace, solar_traces, start_index, end_index, pv_values, b_0);
        if (!curr_run_result.feasible) {
            return ret;
        }
        if (it_count >= mean_window_size) {
            if (window_sum / mean_window_size - curr_run_result.cost < cost_threshold) {
                ret.push_back(curr_run_result);
                return ret;
            }
            window_sum -= mean_window.front();
            mean_window.pop_front();
        }
        window_sum += curr_run_result.cost;
        mean_window.push_back(curr_run_result.cost);

        last_run_result = curr_run_result;
    }

    ret.push_back(last_run_result);
    return ret;
}


default_random_engine generator;

double random_simulate_cheroff(vector <double> &load_trace, vector<vector<double>> &solar_traces,
        double cells, valarray<double> &pvs, double b_0) {

    uniform_int_distribution<size_t> start_index_dist(0, chunk_total);


    size_t n_target_satisfied = 0;

    for (size_t i = 0; i < number_of_chunks; ++i) {
        size_t start_index = start_index_dist(generator);
        double loss = sim(load_trace, solar_traces,
                start_index, start_index + chunk_size, cells, pvs, b_0);
        if (loss < epsilon) {
            ++n_target_satisfied;
        }
//        cout << start_index << "," << loss << "," << (loss < epsilon) << endl;
    }

    return (double) n_target_satisfied / (double) number_of_chunks;

}

double deterministic_simulate_cheroff(vector <double> &load_trace, vector<vector<double>> &solar_traces,
        double cells, valarray<double> &pvs, double b_0) {

    size_t n_target_satisfied = 0;

//    cout << "deterministic_simulate_cheroff"
//         << ", cells=" << cells
//         << ", pv1=" << pvs[0]
//         << ", pv2=" << pvs[1]
//         << endl
//         << "chunk_start,loss,satisfied"
//         << endl;

    for (size_t chunk_num = 0; chunk_num < number_of_chunks; chunk_num += 1) {
        size_t chunk_start = chunk_num * chunk_step;
        double loss = sim(load_trace, solar_traces,
                          chunk_start, chunk_start + chunk_size, cells, pvs, b_0);
        if (loss < epsilon) {
            ++n_target_satisfied;
        }
//        cout << chunk_start << "," << loss << "," << (loss < epsilon) << endl;
    }

    return (double) n_target_satisfied / (double) number_of_chunks;
}

double get_cheroff(double p_tilde) {
    return p_tilde - sqrt( -3 * p_tilde * log(1 - confidence) / number_of_chunks);
}

double chernoff_result(vector <double> &load_trace, vector<vector<double>> &solar_traces,
                       double cells, valarray<double> &pvs, double b_0) {
    double p_tilde = deterministic_simulate_cheroff(load_trace, solar_traces, cells, pvs, b_0);
    double p_delta = get_cheroff(p_tilde);
    return p_delta;
}

vector<SimulationMultiRoofResult>
tabu_cheroff(vector<double> &load_trace, vector<vector<double>> &solar_traces, double target_p, double b_0) {
    // search for an upper, pareto-efficient surface that
    // satisfies the cheroff bound

    // first perform binary search to determine one point
    bool viable = false;
    int l = 0, u = num_steps;
    while (l <= u) {
        int m = (l + u) / 2;
        valarray<double> pv_M = pv_steps * m;
        double cells_M = cells_step * m;

        double p_delta = chernoff_result(load, solar, cells_M, pv_M, b_0);
        SimulationMultiRoofResult result_M(cells_M, pv_M);
//        cout << result_M << "," << p_delta << endl;

        if (p_delta >= target_p) {
            viable = true;
            u = m - 1;
        }
        else
        {
            l = m + 1;
        }
    }

    if (!viable) {
        return vector<SimulationMultiRoofResult>();
    }

    // if viable, (pv_U, cells_U) is the first point on surface.
    // perform tabu search originating from this point

    vector<SimulationMultiRoofResult> boundary;
    unordered_map<string, bool> seen;
    deque<pair<SimulationMultiRoofResult, bool>> search_q;

    double first_cells = cells_step * l;
    valarray<double> first_PV = pv_steps * l;
    SimulationMultiRoofResult first_result(first_cells, first_PV);

    seen[first_result.cells_pv_serialize()] = true;
    search_q.emplace_back(first_result, true);

    while (!search_q.empty()) {
        auto& curr_pair = search_q.front();
        search_q.pop_front();

        SimulationMultiRoofResult curr_result = curr_pair.first;
        bool curr_b = curr_pair.second;

        if (curr_b) {
            // explore all neighbors with less value
            bool is_boundary = true;

            for (size_t i = 0; i < n_solars; ++i) {
                valarray<double> neighbor_pvs = curr_result.PVs;
                neighbor_pvs[i] -= pv_steps[i];
                if (neighbor_pvs[i] < pv_mins[i] - numeric_limits<double>::epsilon()) {
                    continue;
                }

                SimulationMultiRoofResult neighbor_result(curr_result.B, neighbor_pvs);
                string neighbor_serial = neighbor_result.cells_pv_serialize();
                if (seen.count(neighbor_serial)) {
                    if (seen[neighbor_serial]) {
                        is_boundary = false;
                    }
                    continue;
                }

                double p_delta = chernoff_result(load, solar, curr_result.B, neighbor_pvs, b_0);
//                cout << neighbor_serial << "," << p_delta << endl;

                if (p_delta >= target_p) {
                    seen[neighbor_serial] = true;
                    search_q.emplace_back(neighbor_result, true);
                    is_boundary = false;
                } else {
                    seen[neighbor_serial] = false;
                    search_q.emplace_back(neighbor_result, false);
                }
            }

            // also check for one unit less battery
            double neighbor_cells = curr_result.B - cells_step;
            if (neighbor_cells >= cells_min) {
                valarray<double> neighbor_pvs = curr_result.PVs;
                SimulationMultiRoofResult neighbor_result(neighbor_cells, neighbor_pvs);
                string neighbor_serial = neighbor_result.cells_pv_serialize();
                if (seen.count(neighbor_serial)) {
                    if (seen[neighbor_serial]) {
                        is_boundary = false;
                    }
                } else {
                    double p_delta = chernoff_result(load, solar, neighbor_cells, neighbor_pvs, b_0);
//                    cout << neighbor_serial << "," << p_delta << endl;

                    if (p_delta >= target_p) {
                        seen[neighbor_serial] = true;
                        search_q.emplace_back(neighbor_result, true);
                        is_boundary = false;
                    } else {
                        seen[neighbor_serial] = false;
                        search_q.emplace_back(neighbor_result, false);
                    }
                }
            }

            if (is_boundary) {
                boundary.push_back(move(curr_result));
            }
        } else {
            for (size_t i = 0; i < n_solars; ++i) {
                valarray<double> neighbor_pvs = curr_result.PVs;
                neighbor_pvs[i] += pv_steps[i];
                if (neighbor_pvs[i] > pv_maxs[i]) {
                    continue;
                }

                SimulationMultiRoofResult neighbor_result(curr_result.B, neighbor_pvs);
                string neighbor_serial = neighbor_result.cells_pv_serialize();
                if (seen.count(neighbor_serial)) {
                    continue;
                }

                double p_delta = chernoff_result(load, solar, curr_result.B, neighbor_pvs, b_0);
//                cout << neighbor_serial << "," << p_delta << endl;

                bool b = p_delta >= target_p;
                seen[neighbor_serial] = b;
                search_q.emplace_back(neighbor_result, b);
            }

            // also check for one unit more battery
            double neighbor_cells = curr_result.B + cells_step;
            if (neighbor_cells <= cells_max) {
                valarray<double> neighbor_pvs = curr_result.PVs;
                SimulationMultiRoofResult neighbor_result(neighbor_cells, neighbor_pvs);
                string neighbor_serial = neighbor_result.cells_pv_serialize();
                if (seen.count(neighbor_serial)) {
                    continue;
                }

                double p_delta = chernoff_result(load, solar, neighbor_cells, neighbor_pvs, b_0);
//                cout << neighbor_serial << "," << p_delta << endl;

                bool b = p_delta >= target_p;
                seen[neighbor_serial] = b;
                search_q.emplace_back(neighbor_result, b);
            }
        }
    }

    return boundary;
}

