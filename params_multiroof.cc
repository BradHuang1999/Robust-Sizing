//
// Created by Brad Huang on 8/17/20.
//

#include "params_multiroof.h"

using namespace std;

size_t n_solars;

valarray<double> pv_mins;
valarray<double> pv_maxs;
valarray<double> pv_steps; // search in steps of x kW

valarray<double> PV_fix_costs;
valarray<double> PV_invs;

vector<vector<double>> solar;

size_t calc_chebyshev_number_of_chunks() {
    double asymptote = (double)(n_solars + 1) / (1 - confidence);
    double lambda2 = asymptote * lambda_factor;
    double beta = lambda_factor - 1;
    double eta = (lambda2 + sqrt(lambda2 * lambda2 - 4 * beta)) / (2 * beta);
    return (size_t)(eta + 1);
}

void update_number_of_chunks(size_t nchunks) {
    number_of_chunks = nchunks;
    chunk_size = days_in_chunk * 24 / T_u;

    chunk_total = load.size();
    for (auto& solar_vec: solar) {
        chunk_total *= solar_vec.size() / T_yr;
    }
    chunk_step = chunk_total / number_of_chunks;

#ifdef DEBUG
    cout << "chunk_size = " << chunk_size
         << ", total_yrs = " << total_yrs
         << ", chunk_step = " << chunk_step
         << endl;
#endif
}

int process_input(int argc, char **argv) {

    int i = 0;

    // n_solars
    string n_solars_string = argv[++i];
    n_solars = stoi(n_solars_string);

#ifdef DEBUG
    cout << "n_solars_string = " << n_solars_string
         << ", n_solars = " << n_solars << endl;
#endif

    // metric
    string metric_string = argv[++i];
    metric = stoi(metric_string);

#ifdef DEBUG
    cout << "metric_string = " << metric_string
        << ", metric = " << metric << endl;
#endif

    // epsilon
    string epsilon_string = argv[++i];
    epsilon = stod(epsilon_string);

#ifdef DEBUG
    cout << "epsilon_string = " << epsilon_string
         << ", epsilon = " << epsilon << endl;
#endif

    // confidence
    string confidence_string = argv[++i];
    confidence = stod(confidence_string);

#ifdef DEBUG
    cout << "confidence_string = " << confidence_string
         << ", confidence = " << confidence << endl;
#endif

    // days_in_chunk
    string days_in_chunk_string = argv[++i];
    days_in_chunk = stoi(days_in_chunk_string);

#ifdef DEBUG
    cout << "days_in_chunk_string = " << days_in_chunk_string
         << ", days_in_chunk = " << days_in_chunk << endl;
#endif

    // B_inv
    string inv_B_string = argv[++i];
    B_inv = stod(inv_B_string) * kWh_in_one_cell; // convert from per-kWh to per-cell cost

#ifdef DEBUG
    cout << "inv_B_string = " << inv_B_string
         << ", B_inv = " << B_inv << endl;
#endif

    // cells_max
    string cells_max_string = argv[++i];
    cells_max = stod(cells_max_string) / kWh_in_one_cell;

    // set default cells_min and cells_step
    cells_min = 0;
    cells_step = (cells_max - cells_min) / num_steps;

    // load
    string loadfile = argv[++i];

#ifdef DEBUG
    cout << "loadfile = " << loadfile << endl;
#endif

    if (loadfile == string("--")) {
        // read from cin
        int limit = stoi(argv[++i]);

#ifdef DEBUG
        cout << "reading load data from stdin. limit = " << limit << endl;
#endif

        load = read_data_from_file(cin, limit);
    } else {

#ifdef DEBUG
        cout << "reading load file" << endl;
#endif

        // read in data into vector
        ifstream loadstream(loadfile.c_str());
        load = read_data_from_file(loadstream);
    }

#ifdef DEBUG
    cout << "checking for errors in load file..." << endl;
#endif

    if (load[0] < 0) {
        cerr << "error reading load file " << loadfile << endl;
        return 1;
    } else if (load.size() % T_yr > 0) {
        cerr << "load file length needs to be multiple of " << T_yr << endl;
        return 1;
    }

    // pv params
    PV_fix_costs = valarray<double>(n_solars);
    PV_invs = valarray<double>(n_solars);
    pv_mins = valarray<double>(n_solars);
    pv_maxs = valarray<double>(n_solars);
    pv_steps = valarray<double>(n_solars);
    solar = vector<vector<double>>(n_solars);

    for (size_t s = 0; s < n_solars; ++s) {
#ifdef DEBUG
        cout << "solar file " << s << endl;
#endif

        string fix_PV_string = argv[++i];
        PV_fix_costs[s] = stod(fix_PV_string);

#ifdef DEBUG
        cout << "fix_PV_string = " << fix_PV_string
             << ", PV_fix_cost = " << PV_fix_costs[s] << endl;
#endif

        string inv_PV_string = argv[++i];
        PV_invs[s] = stod(inv_PV_string);

#ifdef DEBUG
        cout << "inv_PV_string = " << inv_PV_string
         << ", PV_inv = " << PV_invs[s] << endl;
#endif

        string pv_max_string = argv[++i];
        pv_maxs[s] = stod(pv_max_string);

        // set default pv_min and pv_step
        pv_mins[s] = 0;
        pv_steps[s] = (pv_maxs[s] - pv_mins[s]) / num_steps;

#ifdef DEBUG
        cout << "pv_max_string = " << pv_max_string
         << ", pv_max = " << pv_maxs[s]
         << ", pv_min = " << pv_mins[s]
         << ", pv_step = " << pv_steps[s]
         << endl;
#endif

        string solarfile = argv[++i];

#ifdef DEBUG
        cout << "solarfile = " << solarfile << endl;
#endif

        if (solarfile == string("--")) {

#ifdef DEBUG
            cout << "reading solar file" << endl;
#endif

            // read from cin
            int limit = stoi(argv[++i]);

#ifdef DEBUG
            cout << "reading solar data from stdin. limit = " << limit << endl;
#endif

            solar[s] = read_data_from_file(cin, limit);
        } else {
            // read in data into vector
            ifstream solarstream(solarfile.c_str());
            solar[s] = read_data_from_file(solarstream);
        }

#ifdef DEBUG
        cout << "checking for errors in solar file..." << endl;
#endif

        if (solar[s][0] < 0) {
            cerr << "error reading solar file " << solarfile << endl;
            return 1;
        } else if (solar[s].size() % T_yr > 0) {
            cerr << "solar file length needs to be multiple of " << T_yr << endl;
            return 1;
        }
    }

    update_number_of_chunks(calc_chebyshev_number_of_chunks());

    return 0;
}
