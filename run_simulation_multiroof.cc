//
// Created by Brad Huang on 8/16/20.
//

#include "simulate_multiroof.h"
#include "params_multiroof.h"
#include "cheby_multiroof.h"

using namespace std;

// run_simulations
// load_filename: filename, each line in file contains electricity consumption value
// solar_filename: filename, each line in file contains solar generation value
// id: request id
// metric: 0 for LOLP, 1 for unmet load
// epsilon: number in range [0,1] representing LOLP or unmet load fraction.
// chunk_size: length of time (in days)
vector<vector<SimulationMultiRoofResult>> run_simulations()
{
    vector<vector<SimulationMultiRoofResult>> results;

    // get random start times and run simulation on this chunk of data
    for (size_t chunk_num = 0, chunk_start = 0; chunk_num < number_of_chunks; chunk_num += 1, chunk_start += chunk_step)
    {
        for (size_t run = 0; run < runs_per_chunk; ++run) {
            vector<SimulationMultiRoofResult> pareto_set = simulate(
                    load, solar, chunk_start, chunk_start + chunk_size, 0);
            results.push_back(move(pareto_set));
        }
    }

#ifdef TRACE
    for (auto& pareto_set: results) {
        for (SimulationMultiRoofResult s: pareto_set) {
            cout << s.B << ',';
            for (auto pv: s.PVs) {
                cout << pv << ',';
            }
            cout << s.cost << endl;
        }
    }
#endif

    return results;

    // calculate the chebyshev curves, find the cheapest system along their upper envelope, and return it
    // calculate_sample_bound(results, epsilon, confidence);
}

int main(int argc, char ** argv) {

    int input_process_status = process_input(argc, argv);

    if (input_process_status) {
        cerr << "Illegal input" << endl;
        return 1;
    }

    run_simulations();
//    cout << sr.B << "\t" << sr.C << "\t" << sr.cost << endl;

    return 0;
}
