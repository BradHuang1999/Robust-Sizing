#include <fstream>
#include <sstream>
#include <cstring>
#include <cstdlib>
#include <vector>
#include <iostream>
#include <iomanip>
#include "snc_eue_pertrace.h"
#include "params.h"

using namespace std;

// run_snc_eue
// load_filename: filename, each line in file contains electricity consumption value
// solar_filename: filename, each line in file contains solar generation value
// id: request id
// metric: 0 for LOLP, 1 for unmet load
// epsilon: number in range [0,1] representing LOLP or unmet load fraction.
// chunk_size: length of time (in days)
SimulationResult run_snc_eue(vector<double> &load, vector<double> &solar, double epsilon, double confidence) {
    vector<int> chunk_starts;
    vector<int> chunk_ends;

    for (size_t chunk_num = 0, chunk_start = 0;
         chunk_num < number_of_chunks; chunk_num += 1, chunk_start += chunk_step) {
        size_t chunk_end = chunk_start + chunk_size;
        chunk_starts.push_back(chunk_start);
        chunk_ends.push_back(chunk_end);

        double sum_solar = 0;
        double sum_load = 0;
        for (size_t index = chunk_start; index < chunk_end; index++) {
            sum_solar += solar[index % solar.size()];
            sum_load += load[index % solar.size()];
        }
    }

    return snc_eue(load, solar, chunk_starts, chunk_ends, epsilon, confidence, chunk_size);
}

int main(int argc, char **argv) {

    int input_process_status = process_input(argc, argv, false);

    if (input_process_status) {
        cerr << "Illegal input" << endl;
        return 1;
    }

    SimulationResult sr = run_snc_eue(load, solar, epsilon, confidence);
    cout << sr.B << "\t" << sr.C << "\t" << sr.cost << endl;

    return 0;
}
