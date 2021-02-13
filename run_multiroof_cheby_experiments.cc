//
// Created by Brad Huang on 1/10/21.
//

#include "simulate_multiroof.h"
#include "params_multiroof.h"
#include "cheby_multiroof.h"

#include <iostream>
#include <ctime>
#include <iomanip>
#include <filesystem>

using namespace std;

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

vector<vector<SimulationMultiRoofResult>> parse_csv(const string& file_name,
        size_t type_col = 0, size_t number_of_types = 3, bool normalize_battery = true,
        size_t battery_col = 1, size_t pv_col_start = 2, size_t pv_col_end = 3,
        bool ignore_header = true) {

    if (pv_col_start > pv_col_end) {
        throw range_error("required: pv_col_start <= pv_col_end");
    }
    if (battery_col == type_col) {
        throw range_error("battery_col overlaps with type_col");
    }
    if (pv_col_start <= battery_col && battery_col <= pv_col_end) {
        throw range_error("battery_col overlapping between pv_col_start and pv_col_end");
    }
    if (pv_col_start <= type_col && type_col <= pv_col_end) {
        throw range_error("battery_col overlapping between pv_col_start and pv_col_end");
    }

    size_t pv_col_size = pv_col_end - pv_col_start + 1;
    vector<vector<SimulationMultiRoofResult>> ret(number_of_types);

    ifstream fs(file_name);
    if (fs.fail()) {
        throw runtime_error("file failed at path");
    } else if (!fs) {
        throw runtime_error("fs is empty");
    }

    string curr_line, curr_token;

    if (ignore_header) {
        // ignore the first line as header
        getline(fs, curr_line);
    }

    for (size_t l = 1; getline(fs, curr_line); ++l) {
        istringstream line_iss( curr_line );

        size_t type;
        double battery;
        valarray<double> pvs(pv_col_size);

        size_t i = 0;
        while (getline(line_iss, curr_token, ',')) {
            if (i == type_col) {
                type = stoi(curr_token);
            }

            if (i == battery_col) {
                if (normalize_battery) {
                    battery = stod(curr_token) / kWh_in_one_cell;
                } else {
                    battery = stod(curr_token);
                }
            }

            if (pv_col_start <= i && i <= pv_col_end) {
                pvs[i - pv_col_start] = stod(curr_token);
            }

            ++i;
        }

        if (i <= type_col || i <= battery_col || i <= pv_col_end) {
            stringstream err_msg_ss;
            err_msg_ss << "line " << l << " too short with only " << i << " values";
            throw range_error(err_msg_ss.str());
        }

        ret[type].emplace_back(battery, pvs);
    }

    return ret;
}

SimulationMultiRoofResult min_cost(const vector<SimulationMultiRoofResult> &results) {
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

int main(int argc, char **argv)
{
    int input_process_status = process_input(argc, argv);
    if (input_process_status) {
        cerr << "Illegal input" << endl;
        return 1;
    }

    stringstream filename_ss;

    filename_ss << "results/full_search/" << output_folder_path << "/";
    const string adagrad_sims_path = filename_ss.str();

    filename_ss << "cheby/";
    const string cheby_path = filename_ss.str();

    ofstream summary_os(cheby_path + "summary2.csv");
    summary_os << "days_in_chunk,confidence,epsilon,battery,pv1,pv2,cost" << endl;

    for (const auto& adagrad_sims_file: filesystem::directory_iterator(adagrad_sims_path)) {
        const string adagrad_sims_filename = adagrad_sims_file.path().filename();
        cout << adagrad_sims_filename << endl;

        if (adagrad_sims_filename == "cheby" || adagrad_sims_filename != "02_09_01_57_09_daysinchunk=365_conf=0.75_epsilon=0.45.csv") {
            cout << "skipping" << endl << endl;
            continue;
        }

        const string days_in_chunk_string = adagrad_sims_filename.substr(27, 3);
        cout << days_in_chunk_string << ", ";
        days_in_chunk = stoi(days_in_chunk_string);

        const string conf_string = adagrad_sims_filename.substr(36, 4);
        cout << conf_string << ", ";
        confidence = stod(conf_string);

        const string epsilon_string = adagrad_sims_filename.substr(49, 4);
        cout << epsilon_string << ", " << endl;
        epsilon = stod(epsilon_string);

        cout << "days_in_chunk = " << days_in_chunk << endl
             << "conf = " << confidence << endl
             << "epsilon = " << epsilon << endl
             << endl;

        ofstream os(cheby_path + adagrad_sims_filename);
        os << "type,battery,pv1,pv2,cost" << endl;

        cout << endl << "Operation started at " << curr_time("%m/%d %H:%M:%S") << endl;

        vector<vector<SimulationMultiRoofResult>> adagrad_sims =
                parse_csv(adagrad_sims_file.path().relative_path(),
                        0, (1u << n_solars), true,
                        1, 2, 1 + n_solars, true);

        vector<SimulationMultiRoofResult> min_results;

        for (size_t type = 1; type < (1u << n_solars); ++type) {
            cout << "\n--------------------- start of type " << type << " ---------------------" << endl;
            valarray<bool> is_zeros(n_solars);
            for (size_t i = 0; i < n_solars; ++i) {
                if ((type & (1u << i)) == 0) {
                    is_zeros[i] = true;
                }
                cout << boolalpha << is_zeros[i] << ", ";
            }
            cout << endl;

            vector<SimulationMultiRoofResult> cheby_bound = get_chebyshev_bound(adagrad_sims[type], is_zeros);
            for (const SimulationMultiRoofResult& cb: cheby_bound) {
                os << type + 1 << "," << cb << endl;
            }

            if (cheby_bound.empty()) {
                cout << "CHEBY BOUND EMPTY" << endl;
            } else {
                SimulationMultiRoofResult min_result = min_cost(cheby_bound);
                cout << "\nmin_cost" << endl << min_result << endl << endl;
                min_results.push_back(move(min_result));
            }
            cout << "--------------------- end of type " << type << " ---------------------\n" << endl;
        }

        os.close();

        if (min_results.empty()) {
            cout << "BAD: MIN RESULTS EMPTY" << endl;
        } else {
            SimulationMultiRoofResult min_min_result = min_cost(min_results);
            cout << "min_min_cost" << endl << min_min_result << endl << endl;
            summary_os << days_in_chunk << ","
                       << confidence << ","
                       << epsilon << ","
                       << min_min_result << endl;
        }

        cout << "\nOperation ended at " << curr_time("%m/%d %H:%M:%S") << endl;
    }

    summary_os.close();

    return 0;
}
