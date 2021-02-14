//
// Created by Brad Huang on 9/8/20.
//

#include "params_multiroof.h"
#include "cheby_multiroof.h"
#include "simulate_multiroof.h"

#include <memory>
#include <ctime>
#include <iomanip>
#include <filesystem>

struct curr_time {
    const char *fmt;
    const time_t t;

    curr_time(const char *fmt) : fmt(fmt), t(time(nullptr)) {}

    _Put_time<char> put() const {
        tm *loc_t = localtime(&t);
        return put_time(loc_t, fmt);
    }
};

ostream& operator<<( ostream& os, const curr_time& ct) {
    return os << ct.put();
}

template<typename T>
std::ostream& operator<<( std::ostream& os, const std::valarray<T>& v ) {
    os << "{ ";
    for (const T& d: v) {
        os << d << " ";
    }
    os << "}";
    return os;
}

void adagrad_search(size_t start_index) {
    stringstream filename_ss;
    filename_ss << "results/adagrad_search/" << output_folder_path << "/";
    filesystem::create_directory(filename_ss.str());
    filename_ss << curr_time("%m_%d_%H_%M")
                << "_start=" << start_index
                << "_daysinchunk=" << days_in_chunk
                << "_conf=" << confidence
                << "_epsilon=" << epsilon
                << "_num_steps=" << num_steps
                << "_adagrad_stepsize=" << adagrad_stepsize
                << ".csv";

    ofstream os(filename_ss.str());
    if (os.fail()) {
        throw runtime_error("os failed");
    }
    cout << "Starting adagrad search, file stored at " << filename_ss.str() << endl;
    os << "type,battery,";
    for (size_t i = 0; i < n_solars; ++i) {
        os << "pv" << i << ",";
    }
    os << "cost" << endl;

    // combinatorial explosion
    for (size_t type = 1; type < (1u << n_solars); ++type) {
        if (type_mode != 0 && type_mode != type) {
            continue;
        }

        cout << "\n--------------------- start of type " << type << " ---------------------" << endl;
        cout << endl << "Operation started at " << curr_time("%m/%d %H:%M:%S") << endl;

        valarray<double> pv_start = pv_maxs;
        valarray<bool> is_zeros = (pv_maxs == (double) 0);

        bool skip = false;
        for (size_t i = 0; i < n_solars; ++i) {
            if ((type & (1u << i)) == 0) {
                if (is_zeros[i]) {
                    skip = true;
                    break;
                }

                pv_start[i] = 0;
                is_zeros[i] = true;
            }
        }
        if (skip) {
            continue;
        }

        cout << "training len: " << load.size() << endl;
        cout << "pv_start: " << pv_start << endl;
        cout << "is_zeros: " << is_zeros << endl;

        update_chebyshev_params(is_zeros);
        update_number_of_chunks();

        vector<SimulationMultiRoofResult> adagrad_results = simulate_deterministic_adagrad(
                load, solar, start_index, start_index + chunk_size, pv_start, is_zeros);

        for (SimulationMultiRoofResult& r: adagrad_results) {
            os << type << "," << r << endl;
        }
    }
}

void grid_search(size_t start_index) {
    valarray<bool> is_zeros(false, n_solars);
    update_chebyshev_params(is_zeros);
    update_number_of_chunks();

    vector<SimulationMultiRoofResult> results = vector<SimulationMultiRoofResult>();

    stringstream filename_ss;
    filename_ss << "results/grid_search/" << output_folder_path << "/";
    filesystem::create_directory(filename_ss.str());
    filename_ss << curr_time("%m_%d_%H_%M")
                << "_start=" << start_index
                << "_daysinchunk=" << days_in_chunk
                << "_conf=" << confidence
                << "_epsilon=" << epsilon
                << "_num_steps=" << num_steps
                << ".csv";

    ofstream os(filename_ss.str());
    if (os.fail()) {
        throw runtime_error("os failed");
    }
    os << "battery,pv1,pv2,cost" << endl;
    cout << "Starting grid search, file stored at " << filename_ss.str() << endl;

    for (double pv1 = pv_mins[0]; pv1 <= pv_maxs[0]; pv1 += pv_steps[0]) {
        valarray<double> pvs(2);
        pvs[0] = pv1;
        for (double pv2 = pv_mins[1]; pv2 <= pv_maxs[1]; pv2 += pv_steps[1]) {
            pvs[1] = pv2;
            SimulationMultiRoofResult result = binary_search_result(
                    load, solar, start_index, start_index + chunk_size, pvs);
            if (result.feasible) {
                os << result << endl;
            }
        }
    }
}

SimulationMultiRoofResult min_cost(const vector<SimulationMultiRoofResult> &results) {
    SimulationMultiRoofResult min_result;
    double min_val = INFTY;
    for (const auto& s: results) {
        if (s.cost < min_val) {
            min_val = s.cost;
            min_result = s;
        }
    }
    return min_result;
}

SimulationMultiRoofResult full_search(const string& foldername, const curr_time& ct) {
    filesystem::create_directory(foldername);

    stringstream cheby_ss;
    cheby_ss << foldername << "cheby/";
    filesystem::create_directory(cheby_ss.str());

    cout << "Starting full search, files stored at " << foldername << endl;

    vector<SimulationMultiRoofResult> min_chebys;

    // combinatorial explosion
    for (size_t type = 1; type < (1u << n_solars); ++type) {
        if (type_mode != 0 && type_mode != type) {
            continue;
        }

        cout << "\n--------------------- start of type " << type << " ---------------------" << endl;
        cout << endl << "Operation started at " << curr_time("%m/%d %H:%M:%S") << endl;

        stringstream adagrad_results_ss;
        adagrad_results_ss << foldername
                           << ct
                           << "_type=" << type
                           << "_daysinchunk=" << days_in_chunk
                           << "_conf=" << confidence
                           << "_epsilon=" << epsilon
                           << ".csv";
        ofstream adagrad_results_os(adagrad_results_ss.str());
        if (adagrad_results_os.fail()) {
            throw runtime_error("adagrad_results_os failed");
        }
        adagrad_results_os << "type,battery,";
        for (size_t i = 0; i < n_solars; ++i) {
            adagrad_results_os << "pv" << i << ",";
        }
        adagrad_results_os << "cost" << endl;

        stringstream cheby_results_ss;
        cheby_results_ss << cheby_ss.str()
                         << ct
                         << "_type=" << type
                         << "_daysinchunk=" << days_in_chunk
                         << "_conf=" << confidence
                         << "_epsilon=" << epsilon
                         << ".csv";
        ofstream cheby_results_os(cheby_results_ss.str());
        if (cheby_results_os.fail()) {
            throw runtime_error("adagrad_results_os failed");
        }
        cheby_results_os << "type,battery,";
        for (size_t i = 0; i < n_solars; ++i) {
            cheby_results_os << "pv" << i << ",";
        }
        cheby_results_os << "cost" << endl;

        valarray<double> pv_start = pv_maxs;
        valarray<bool> is_zeros = (pv_maxs == (double)0);

        bool skip = false;
        for (size_t i = 0; i < n_solars; ++i) {
            if ((type & (1u << i)) == 0) {
                if (is_zeros[i]) {
                    skip = true;
                    break;
                }

                pv_start[i] = 0;
                is_zeros[i] = true;
            }
        }
        if (skip) {
            continue;
        }

        cout << "training len: " << load.size() << endl;
        cout << "pv_start: " << pv_start << endl;
        cout << "is_zeros: " << is_zeros << endl;

        update_chebyshev_params(is_zeros);
        update_number_of_chunks();

        cout << "days_in_chunk = " << days_in_chunk << endl
             << "conf = " << confidence << endl
             << "epsilon = " << epsilon << endl
             << "lambda2 = " << lambda2 << endl
             << "number_of_chunks = " << number_of_chunks << endl
             << endl;

        time_t t_before = time(nullptr);

        vector<SimulationMultiRoofResult> min_adagrads;

        for (size_t chunk_num = 0; chunk_num < number_of_chunks; chunk_num += 1) {
            size_t chunk_start = chunk_num * chunk_step;
            vector<SimulationMultiRoofResult> adagrad_results = simulate_deterministic_adagrad(
                    load, solar, chunk_start, chunk_start + chunk_size, pv_start, is_zeros);

            SimulationMultiRoofResult min_adagrad_result = min_cost(adagrad_results);
            min_adagrads.push_back(min_adagrad_result);
        }

        // print out min_adagrad points
        for (const SimulationMultiRoofResult& r: min_adagrads) {
            adagrad_results_os << type << "," << r << endl;
        }
        adagrad_results_os.close();

        cout << "Adagrad ended at " << curr_time("%m/%d %H:%M:%S") << ", starting cheby bound search" << endl << endl;

        // get cheby bound
        vector<SimulationMultiRoofResult> cheby_results = get_chebyshev_bound(min_adagrads, is_zeros);

        // print cheby bound
        for (const SimulationMultiRoofResult& r: cheby_results) {
            cheby_results_os << type << "," << r << endl;
        }
        cheby_results_os.close();

        SimulationMultiRoofResult min_cheby_result = min_cost(cheby_results);
        min_chebys.push_back(min_cheby_result);

        cout << "min_cheby_result: " << min_cheby_result << endl;
        cout << "\nOperation ended at " << curr_time("%m/%d %H:%M:%S") << endl;

        time_t t_after = time(nullptr);

        stringstream type_summary_ss;
        type_summary_ss << foldername << "type_summary" << n_solars << ".csv";
        ofstream type_summary_os;
        type_summary_os.open(type_summary_ss.str(), ios_base::app);
        type_summary_os << ct << ","
                        << output_folder_path << ","
                        << type << ","
                        << days_in_chunk << ","
                        << confidence << ","
                        << epsilon << ","
                        << min_cheby_result.feasible << ","
                        << min_cheby_result << ","
                        << (t_after - t_before) << ","
                        << total_sim_called
                        << endl;
        type_summary_os.close();
    }

    SimulationMultiRoofResult ret_result = min_cost(min_chebys);
    cout << "\n--------------------- final result ---------------------" << endl;
    cout << ret_result << endl;

    return ret_result;
}

void full_search_with_cross_validation(const string& foldername) {
    curr_time ct("%m_%d_%H_%M_%S");

    for (size_t s = 0; s < n_solars; ++s) {
        if (solar[s].size() != load.size()) {
            throw range_error("solar size is not equal to load size. cannot perform cross validation");
        }
    }
    size_t tot_yrs = load.size() / T_yr;

    vector<double> load_bk((tot_yrs - 1) * T_yr);
    swap(load, load_bk);

    vector<vector<double>> solar_bk(n_solars, vector<double>((tot_yrs - 1) * T_yr));
    swap(solar, solar_bk);

    vector<double> load_validation(T_yr);
    vector<vector<double>> solar_validation(n_solars, vector<double>(T_yr));

    for (size_t validation_yr = 0; validation_yr < tot_yrs; ++validation_yr) {
        cout << "\n\n------------------------------- validation year " << validation_yr << " -------------------------------\n\n" << endl;

        size_t validation_begin = validation_yr * T_yr;
        size_t validation_end = (validation_yr + 1) * T_yr;

        // copy from ..v_begin to load (training set)
        copy(load_bk.begin(),
             load_bk.begin() + validation_begin,
             load.begin());
        // copy from v_begin..v_end to load_validation
        copy(load_bk.begin() + validation_begin,
             load_bk.begin() + validation_end,
             load_validation.begin());
        // copy from v_end.. to load (training set)
        copy(load_bk.begin() + validation_end,
             load_bk.end(),
             load.begin() + validation_begin);

        for (size_t s = 0; s < n_solars; ++s) {
            // copy from ..v_begin to load (training set)
            copy(solar_bk[s].begin(),
                 solar_bk[s].begin() + validation_begin,
                 solar[s].begin());
            // copy from v_begin..v_end to load_validation
            copy(solar_bk[s].begin() + validation_begin,
                 solar_bk[s].begin() + validation_end,
                 solar_validation[s].begin());
            // copy from v_end.. to load (training set)
            copy(solar_bk[s].begin() + validation_end,
                 solar_bk[s].end(),
                 solar[s].begin() + validation_begin);
        }

        stringstream full_search_foldername_ss;
        full_search_foldername_ss << "results/full_search/" << foldername;
        filesystem::create_directory(full_search_foldername_ss.str());

        full_search_foldername_ss << "/yr" << validation_yr << "/";
        filesystem::create_directory(full_search_foldername_ss.str());

        total_sim_called = 0;
        time_t t_before = time(nullptr);
        SimulationMultiRoofResult train_result = full_search(full_search_foldername_ss.str(), ct);
        time_t t_after = time(nullptr);

        double validation_loss;
        if (train_result.feasible) {
            cout << "sizing feasible" << endl;
            validation_loss = sim(load_validation, solar_validation,
                        0, T_yr, train_result.B, train_result.PVs);
            cout << "feasible, loss = " << validation_loss << endl;
//            double sum_validation_loss = 0;
//            vector<double> load_validation(T_yr);
//            for (size_t yr = 0; yr < tot_yrs; ++yr) {
//                size_t load_validation_begin = yr * T_yr;
//                size_t load_validation_end = (yr + 1) * T_yr;
//                copy(load.begin() + load_validation_begin,
//                     load.begin() + load_validation_end,
//                     load_validation.begin());
//                double loss = sim(load_validation, solar_validation,
//                                  0, T_yr, train_result.B, train_result.PVs);
//                cout << "yr " << yr << " load validation loss = " << loss << endl;
//                sum_validation_loss += loss;
//            }
//            validation_loss = sum_validation_loss / tot_yrs;
//            cout << "avg validation loss = " << validation_loss << endl;
        } else {
            validation_loss = INFTY;
            cout << "infeasible sizing" << endl;
        }

        // print out a the returned cost line on summary2.csv
        stringstream cheby_summary_ss;
        cheby_summary_ss << "results/full_search/cheby_summary/summary" << n_solars << ".csv";
        ofstream cheby_summary_os;
        cheby_summary_os.open(cheby_summary_ss.str(), ios_base::app);
        cheby_summary_os << ct << ","
                         << foldername << ","
                         << validation_yr << ","
                         << days_in_chunk << ","
                         << confidence << ","
                         << epsilon << ","
                         << train_result.feasible << ","
                         << train_result << ","
                         << validation_loss << ","
                         << (t_after - t_before) << ","
                         << total_sim_called << endl;

        cout << "\n\n------------------------------- validation year " << validation_yr
             << " finished -------------------------------\n\n" << endl;
    }
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
    if (os.fail()) {
        throw runtime_error("os failed");
    }
    os << "battery,pv1,pv2,cost" << endl;

    for (const SimulationMultiRoofResult& s: ret) {
        os << s << endl;
    }

    if (ret.empty()) {
        cout << "(no viable result)";
    } else {
        SimulationMultiRoofResult min_search = min_cost(ret);
        cout << min_search;
    }
}

void experiment(const string& foldername="") {
//    vector<double> epsilon_vals {0.1, 0.3};
//    vector<double> epsilon_vals {0.1};
//    vector<double> epsilon_vals {0.3};
    vector<double> epsilon_vals {0.5};

    vector<double> p_vals {0.9, 0.8, 0.5};
//    vector<double> p_vals {0.5};
//    vector<double> p_vals {0.8};
//    vector<double> p_vals {0.9};

    vector<double> conf_vals {0.95, 0.85};
//    vector<double> conf_vals {0.85};
//    vector<double> conf_vals {0.95};

//    vector<size_t> days_in_chunk_vals {100, 200, 365};
    vector<size_t> days_in_chunk_vals {365};

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

                    // RUN CHEBYSHEV EXPERIMENT
                    epsilon = effective_epsilon;

                    cout << "chebyshev"
                         << "," << days_in_chunk_val
                         << "," << conf_val
                         << "," << epsilon_val
                         << "," << p_val
                         << "," << effective_epsilon
                         << ",";

                    time_t t_before_cheby = time(nullptr);
                    total_sim_called = 0;
                    full_search(foldername, curr_time("%m_%d_%H_%M_%S"));
                    time_t t_after_cheby = time(nullptr);
                    time_t t_duration_cheby = t_after_cheby - t_before_cheby;

                    cout << "-1,-1,-1,-1," << t_duration_cheby
                         << "," << total_sim_called
                         << endl;
                }
            }
        }
    }
}

int main(int argc, char **argv) {

    int input_process_status = process_input(argc, argv);

//    // manually change parameters. not recommended.
//    days_in_chunk = 100;
//    confidence = 0.85;
//    epsilon = 0.1;
//    double p = 0.8;
//    update_number_of_chunks(1000);

    cout << "num_steps = " << num_steps << endl;
    cout << "output_folder_path = " << output_folder_path << endl << endl;

    cout << "Operation started at " << curr_time("%m/%d %H:%M:%S") << endl;

    if (input_process_status) {
        cerr << "Illegal input" << endl;
        return 1;
    }

    stringstream foldername_ss;

    switch (search_mode) {
        case 'g':
            grid_search(0);
            break;
        case 'a':
            adagrad_search(0);
            break;
        case 'r':
            chernoff_grid_map();
            break;
        case 'c':
            chernoff_search(0.7);
            break;
        case 't':
            chernoff_tabu_search(0.7, output_folder_path);
            break;
        case 'f':
            foldername_ss << "results/full_search/" << output_folder_path << "/";
            full_search(foldername_ss.str(), curr_time("%m_%d_%H_%M_%S"));
            break;
        case 'v':
            full_search_with_cross_validation(output_folder_path);
            break;
        default:
            throw runtime_error("search mode not found");
    }

    cout << "\nOperation ended at " << curr_time("%m/%d %H:%M:%S") << endl;

    return 0;
}
