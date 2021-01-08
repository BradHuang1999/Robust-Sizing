// params.h
#ifndef PARAMS_H
#define PARAMS_H

#include "params_common.h"

extern double PV_inv; // cost per unit (kW) of PV
extern vector<double> solar;

// define the upper and lower values to test for pv,
// as well as the step size of the search
extern double pv_min;
extern double pv_max;
extern double pv_step; // search in steps of x kW

struct SimulationResult {

	double B;
	double C;
	double cost;

	SimulationResult(double B_val, double C_val, double cost_val) : 
					B(B_val), C(C_val), cost(cost_val) {}

};

int process_input(int argc, char **argv, bool process_metric_input);

#endif
