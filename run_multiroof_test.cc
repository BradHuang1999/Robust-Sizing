//
// Created by Brad Huang on 9/8/20.
//

#include "simulate_multiroof.h"
#include "params_multiroof.h"
#include <memory>
#include <valarray>
#include <ctime>
#include <iomanip>

struct curr_time {
    const char *fmt;
    curr_time(const char *fmt) : fmt(fmt) {}
};

tm *get_local_time() {
    auto t = time(nullptr);
    return localtime(&t);
}

ostream& operator<<( ostream& os, const curr_time ct) {
    tm *local_time = get_local_time();
    os << put_time(local_time, ct.fmt);
    return os;
}

void adagrad_search(size_t start_index, valarray<double> init_pv, double step_size) {
    vector<SimulationMultiRoofResult> adagrad_results = simulate_deterministic_adagrad(
            load, solar, start_index, start_index + chunk_size, init_pv, 0, step_size);

    stringstream filename_ss;
    filename_ss << "results/adagrad_search/" << curr_time("%m_%d_%H_%M")
                << "_start=" << start_index
                << "_stepsize=" << step_size
                << "_initpv=";
    for (double p: init_pv) {
        filename_ss << ((p > 0) ? 1 : 0);
    }
    filename_ss << ".csv";

    ofstream os(filename_ss.str());
    os << "battery,pv1,pv2,cost,chunk_start" << endl;
    for (SimulationMultiRoofResult& r: adagrad_results) {
        os << r << endl;
    }
}

void grid_search(size_t start_index) {
    vector<SimulationMultiRoofResult> results = vector<SimulationMultiRoofResult>();

    stringstream filename_ss;
    filename_ss << "results/grid_search/" << curr_time("%m_%d_%H_%M")
                << "_start=" << start_index
                << ".csv";

    for (size_t pv1_idx = 0; pv1_idx < 350; ++pv1_idx) {
        double pv1 = pv_mins[0] + pv1_idx * pv_steps[0];
        valarray<double> pvs(2);
        pvs[0] = pv1;
        for (double pv2 = pv_mins[1]; pv2 <= pv_maxs[1]; pv2 += pv_steps[1]) {
            pvs[1] = pv2;
            SimulationMultiRoofResult result = binary_search_result(
                    load, solar, start_index, start_index + chunk_size, pvs);
            if (result.feasible) {
                results.push_back(result);
            }
        }
    }

    ofstream os(filename_ss.str());
    os << "battery,pv1,pv2,cost,chunk_start" << endl;
    for (SimulationMultiRoofResult& r: results) {
        os << r << endl;
    }
}

SimulationMultiRoofResult min(const vector<SimulationMultiRoofResult> &results) {
    SimulationMultiRoofResult min_result = results[0];
    double min_val = results[0].cost;
    for (const auto& s: results) {
        if (s.cost < min_val) {
            min_val = s.cost;
            min_result = s;
        }
    }
    return min_result;
}

void full_search(const string& foldername="") {
    valarray<double> pv1_maxs = {pv_maxs[0], 0};
    valarray<double> pv2_maxs = {0, pv_maxs[1]};

    vector<SimulationMultiRoofResult> adagrad_mins, adagrad_pv1_mins, adagrad_pv2_mins;

    stringstream filename_ss;
    filename_ss << "results/full_search/" << foldername
                << curr_time("%m_%d_%H_%M")
                << "_daysinchunk=" << days_in_chunk
                << "_conf=" << confidence
                << "_epsilon=" << epsilon
                << ".csv";

    ofstream os(filename_ss.str());
    os << "type,battery,pv1,pv2,cost,chunk_start" << endl;

    for (size_t chunk_num = 0; chunk_num < number_of_chunks; chunk_num += 1)
    {
        size_t chunk_start = chunk_num * chunk_step;

        vector<SimulationMultiRoofResult> adagrad_results = simulate_deterministic_adagrad(
                load, solar, chunk_start, chunk_start + chunk_size, pv_maxs, 0, 3);
        adagrad_mins.push_back(min(adagrad_results));

        vector<SimulationMultiRoofResult> adagrad_pv1_results = simulate_deterministic_adagrad(
                load, solar, chunk_start, chunk_start + chunk_size, pv1_maxs, 0, 3);
        adagrad_pv1_mins.push_back(min(adagrad_pv1_results));

        vector<SimulationMultiRoofResult> adagrad_pv2_results = simulate_deterministic_adagrad(
                load, solar, chunk_start, chunk_start + chunk_size, pv2_maxs, 0, 3);
        adagrad_pv2_mins.push_back(min(adagrad_pv2_results));
    }

    for (SimulationMultiRoofResult& s: adagrad_mins) {
        os << "0," << s << endl;
    }
    for (SimulationMultiRoofResult& s: adagrad_pv1_mins) {
        os << "1," << s << endl;
    }
    for (SimulationMultiRoofResult& s: adagrad_pv2_mins) {
        os << "2," << s << endl;
    }

    os.close();
}

void chernoff_grid_map() {
    valarray<double> pvs(2);

    stringstream filename_ss;
    filename_ss << "results/chernoff_grid_map/" << curr_time("%m_%d_%H_%M")
                << "_conf=" << confidence
                << ".csv";

    ofstream os(filename_ss.str());
    os << "battery,pv1,pv2,p_tilde,p_delta" << endl;

    for (double pv1 = pv_mins[0]; pv1 < pv_maxs[0]; pv1 += pv_steps[0]) {
        pvs[0] = pv1;
        for (double pv2 = pv_mins[1]; pv2 <= pv_maxs[1]; pv2 += pv_steps[1]) {
            pvs[1] = pv2;
            for (double cell = cells_min; cell <= cells_max; cell += cells_step) {
                double p_tilde = random_simulate_cheroff(load, solar, cell, pvs);
                double p_delta = get_cheroff(p_tilde);
                os << cell << ","
                   << pv1 << ","
                   << pv2 << ","
                   << p_tilde << ","
                   << p_delta << endl;
            }
        }
    }

    os.close();
}

void chernoff_search(double target_p) {
    valarray<double> pvs(2);

    stringstream filename_ss;
    filename_ss << "results/chernoff_search/" << curr_time("%m_%d_%H_%M")
                << "_daysinchunk=" << days_in_chunk
                << "_conf=" << confidence
                << "_epsilon=" << epsilon
                << "_targetp=" << target_p
                << "_n=" << number_of_chunks
                << ".csv";

    ofstream os(filename_ss.str());
    os << "battery,pv1,pv2,p_delta" << endl;

    for (double pv1 = pv_mins[0]; pv1 < pv_maxs[0]; pv1 += pv_steps[0]) {
        pvs[0] = pv1;
        for (double pv2 = pv_mins[1]; pv2 <= pv_maxs[1]; pv2 += pv_steps[1]) {
            pvs[1] = pv2;

            double cells_U = cells_max, cells_L = cells_min;
            double p_delta_U = 0;

            while (cells_U - cells_L > cells_step)
            {
                double mid_cells = (cells_L + cells_U) / 2.0;
                double p_delta = chernoff_result(load, solar, mid_cells, pvs);

                if (p_delta >= target_p)
                {
                    cells_U = mid_cells;
                    p_delta_U = p_delta;
                }
                else
                {
                    // (loss < epsilon)
                    cells_L = mid_cells;
                }
            }

            if (p_delta_U >= target_p) {
                os << cells_L << ","
                   << pv1 << ","
                   << pv2 << ","
                   << p_delta_U << endl;
            }
        }
    }

    os.close();
}

void chernoff_tabu_search(double target_p, const string& foldername="") {
    stringstream filename_ss;
    filename_ss << "results/chernoff_tabu_search/" << foldername << curr_time("%m_%d_%H_%M")
                << "_daysinchunk=" << days_in_chunk
                << "_conf=" << confidence
                << "_epsilon=" << epsilon
                << "_targetp=" << target_p
                << ".csv";

    vector<SimulationMultiRoofResult> ret = tabu_cheroff(load, solar, target_p);

    ofstream os(filename_ss.str());
    os << "battery,pv1,pv2,cost" << endl;

    for (const SimulationMultiRoofResult& s: ret) {
        os << s << endl;
    }

    if (ret.empty()) {
        cout << "(no viable result)";
    } else {
        SimulationMultiRoofResult min_search = min(ret);
        cout << min_search;
    }
}

void experiment(const string& foldername="") {
    vector<double> epsilon_vals {0.1, 0.3};
//    vector<double> epsilon_vals {0.1};
//    vector<double> epsilon_vals {0.3};

//    vector<double> p_vals {0.8, 0.9};
    vector<double> p_vals {0.5};
//    vector<double> p_vals {0.8};
//    vector<double> p_vals {0.9};

    vector<double> conf_vals {0.85, 0.95};
//    vector<double> conf_vals {0.85};
//    vector<double> conf_vals {0.95};

    vector<size_t> days_in_chunk_vals {100, 200, 365};
//    vector<size_t> days_in_chunk_vals {365};

    cout << "\nsim_type"
         << ",days_in_chunk"
         << ",confidence"
         << ",epsilon"
         << ",p"
         << ",effective_epsilon"
         << ",battery"
         << ",pv1"
         << ",pv2"
         << ",cost"
         << ",duration"
         << ",total_sim_called"
         << endl;

    for (double epsilon_val: epsilon_vals) {
        for (double p_val: p_vals) {
            for (double conf_val: conf_vals) {
                for (size_t days_in_chunk_val: days_in_chunk_vals) {

                    confidence = conf_val;
                    days_in_chunk = days_in_chunk_val;

                    double effective_epsilon = 1 - (p_val * (1 - epsilon_val));

                    // RUN CHEROFF EXPERIMENT
                    cout << "chernoff"
                         << "," << days_in_chunk_val
                         << "," << conf_val
                         << "," << epsilon_val
                         << "," << p_val
                         << "," << effective_epsilon
                         << ",";

                    epsilon = epsilon_val;

                    time_t t_before_cheroff = time(nullptr);
                    update_number_of_chunks(1000);
                    total_sim_called = 0;
                    chernoff_tabu_search(p_val, foldername);
                    time_t t_after_cheroff = time(nullptr);
                    time_t t_duration_cheroff = t_after_cheroff - t_before_cheroff;

                    cout << "," << t_duration_cheroff
                         << "," << total_sim_called
                         << endl;

//                    // RUN CHEBYSHEV EXPERIMENT
//                    epsilon = effective_epsilon;
//
//                    cout << "chebyshev"
//                         << "," << days_in_chunk_val
//                         << "," << conf_val
//                         << "," << epsilon_val
//                         << "," << p_val
//                         << "," << effective_epsilon
//                         << ",";
//
//                    time_t t_before_cheby = time(nullptr);
//                    update_number_of_chunks(calc_chebyshev_number_of_chunks());
//                    total_sim_called = 0;
//                    full_search(foldername);
//                    time_t t_after_cheby = time(nullptr);
//                    time_t t_duration_cheby = t_after_cheby - t_before_cheby;
//
//                    cout << "-1,-1,-1,-1," << t_duration_cheby
//                         << "," << total_sim_called
//                         << endl;
                }
            }
        }
    }
}

int main(int argc, char **argv) {

    int input_process_status = process_input(argc, argv);

    cout << "Operation started at " << curr_time("%m/%d %H:%M:%S") << endl;

//    calc_chebyshev_number_of_chunks();
//    update_number_of_chunks(32);

    if (input_process_status) {
        cerr << "Illegal input" << endl;
        return 1;
    }

//    grid_search(32240);
//    adagrad_search(15000, pv_maxs, 3);
//    adagrad_search(32240, {pv_maxs[0], 0}, 3);

//    chernoff_grid_map(1000);
//    chernoff_search(0.9, 1000);

//    12_23_02_18_daysinchunk=100_conf=0.85_epsilon=0.1_targetp=0.8.csv
//    days_in_chunk = 100;
//    confidence = 0.85;
//    epsilon = 0.1;
//    double p = 0.8;
//    update_number_of_chunks(1000);

//    chernoff_search(p);
//    chernoff_tabu_search(p, "load=3482_pv=7989 6423/");

    string folder_name = "load=3482_pv=7989 6423_numsteps=20/";

    cout << "num_steps = " << num_steps << endl;
    cout << "folder_name = " << folder_name << endl;

    experiment(folder_name);

//    confidence = 0.95;
//    epsilon = 0.65;
//    days_in_chunk = 365;
//    update_number_of_chunks(1000);
//
//    valarray<double> pvs = {0, 5.9914981617647065};
//    double cell = 7.952769535294118 / kWh_in_one_cell;
//
//    double p = random_simulate_cheroff(load, solar, cell, pvs);
//    double pdelta = get_cheroff(p);
//    cout << p << ", " << pdelta << endl;

//    confidence = 0.95;
//    epsilon = 1 - (0.9 * (1 - 0.1));
//    days_in_chunk = 365;
//    update_number_of_chunks(calc_chebyshev_number_of_chunks());
//
//    full_search("load=6423_pv=7989 6423/");

    cout << "\nOperation ended at " << curr_time("%m/%d %H:%M:%S") << endl;

    return 0;
}