#include "simulate_system.h"
#include "params.h"
#include "cheby.h"

using namespace std;

// run_simulations
// load_filename: filename, each line in file contains electricity consumption value
// solar_filename: filename, each line in file contains solar generation value
// id: request id
// metric: 0 for LOLP, 1 for unmet load
// epsilon: number in range [0,1] representing LOLP or unmet load fraction.
// chunk_size: length of time (in days)
SimulationResult run_simulations()
{
	vector<vector<SimulationResult>> results;

	// get random start times and run simulation on this chunk of data
	for (size_t chunk_num = 0, chunk_start = 0; chunk_num < number_of_chunks; chunk_num += 1, chunk_start += chunk_step)
	{
		results.push_back(simulate(load, solar, chunk_start, chunk_start + chunk_size, 0));
	}

#ifdef DEBUG
	// print all of the curves
	int chunk_index = 1;
	cout << "DEBUG: sizing_curves" << endl;
	for (auto it = results.begin(); it != results.end(); ++it, ++chunk_index) {
		cout << "chunk_" << chunk_index << endl;
		for (auto & it2 : *it) {
			cout << it2.B << "\t" << it2.C << "\t" << it2.cost << endl;
		}
	}
	cout << "DEBUG: sizing_curves_end" << endl;
#endif

	// calculate the chebyshev curves, find the cheapest system along their upper envelope, and return it
	return calculate_sample_bound(results, epsilon, confidence);
}

int main(int argc, char ** argv) {
	int input_process_status = process_input(argc, argv, true);

	if (input_process_status) {
		cerr << "Illegal input" << endl;
		return 1;
	}
	
	SimulationResult sr = run_simulations();
	cout << sr.B << "\t" << sr.C << "\t" << sr.cost << endl;

	return 0;
}
