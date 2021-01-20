//
// Created by Brad Huang on 8/17/20.
//

#ifndef ROBUST_SIZING_PARAMS_MULTIROOF_H
#define ROBUST_SIZING_PARAMS_MULTIROOF_H

#include "params_common.h"
#include <utility>
#include <valarray>

const size_t static runs_per_chunk = 5;

extern string output_folder_path; // an optional path to output process file

extern size_t n_solars;

extern valarray<double> PV_fix_costs; // initial cost for installing the PV, regardless of unit
extern valarray<double> PV_invs; // cost per unit (kW) of PV
extern vector<vector<double>> solar;

// define the upper and lower values to test for pv,
// as well as the step size of the search
extern valarray<double> pv_mins;
extern valarray<double> pv_maxs;
extern valarray<double> pv_steps; // search in steps of x kW

struct SimulationMultiRoofResult {

    bool feasible;
    double B;
    valarray<double> PVs;
    double cost;

    SimulationMultiRoofResult() :
            feasible(false), B(INFTY), PVs(n_solars), cost(INFTY) {
        PVs = INFTY;
    }

    SimulationMultiRoofResult(double B_val, valarray<double> PVs_val):
            feasible(true), B(B_val), PVs(move(PVs_val))
    {
        cost = B * B_inv;
        for (size_t i = 0; i < n_solars; ++i) {
            double pv_value = PVs[i];
            if (pv_value >= numeric_limits<double>::epsilon()) {
                cost += pv_value * PV_invs[i] + PV_fix_costs[i];
            }
        }
    }

public:

    string cells_pv_serialize() const;

    friend ostream& operator<< (ostream& os, const SimulationMultiRoofResult& result);
};

void update_number_of_chunks(size_t nchunks=number_of_chunks);

int process_input(int argc, char **argv);

#endif //ROBUST_SIZING_PARAMS_MULTIROOF_H
