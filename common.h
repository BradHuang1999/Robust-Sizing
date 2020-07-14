// common.h
#ifndef COMMON_H
#define COMMON_H

#include <vector>
#include <limits>

using namespace std;

// INPUTS

extern double B_inv; // cost per cell
extern double PV_inv; // cost per unit (kW) of PV
extern double epsilon;
extern double confidence;
extern int metric;

extern size_t days_in_chunk;
extern size_t chunk_size;
extern size_t chunk_step;

extern vector<double> load;
extern vector<double> solar;

// define the upper and lower values to test for battery cells and pv,
// as well as the step size of the search
extern double cells_min;
extern double cells_max;
extern double cells_step; // search in step of x cells

extern double pv_min;
extern double pv_max;
extern double pv_step; // search in steps of x kW

// CONSTANTS

// defines the number of samples, set via command line input
size_t static number_of_chunks = 100;

/**
 * T_u: this is the time unit, representing the number of hours in
 *      each time slot of the load and solar traces
 */
size_t static T_u = 1;

/**
 * T_yr: this is year unit, representing the number of traces that constitutes a year.
 *       Inputs must have multiples of this size.
 */
size_t static T_yr = 365 * 24 / T_u;

double static kWh_in_one_cell = 0.011284;
double static num_cells_steps = 400; // search in total of n steps for cells
double static num_pv_steps = 350; // search in total of n steps for pv

double static INFTY = numeric_limits<double>::infinity();

struct SimulationResult {

	double B;
	double C;
	double cost;

	SimulationResult(double B_val, double C_val, double cost_val) : 
					B(B_val), C(C_val), cost(cost_val) {}

};

vector<double> read_data_from_file(string);

int process_input(int argc, char **argv, bool process_metric_input);

#endif
