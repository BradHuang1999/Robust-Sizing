//
// Created by Brad Huang on 8/17/20.
//

#ifndef ROBUST_SIZING_PARAMS_MULTIROOF_H
#define ROBUST_SIZING_PARAMS_MULTIROOF_H

#include "params_common.h"
#include <valarray>

const size_t static runs_per_chunk = 5;

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
    size_t chunk_start;

    SimulationMultiRoofResult(size_t chunk_start = -1): feasible(false), chunk_start(chunk_start) {}

    SimulationMultiRoofResult(double B_val, valarray<double>& PVs_val, size_t chunk_start = -1):
            feasible(true), B(B_val), PVs(PVs_val), chunk_start(chunk_start)
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

    string cells_pv_serialize() const {
        stringstream ss;
        ss.precision(8);
        ss << (B * kWh_in_one_cell);
        for (const double &pv: PVs) {
            ss << "," << pv;
        }
        return ss.str();
    }

    friend ostream& operator<<( ostream& os, const SimulationMultiRoofResult& result) {
        if (result.feasible) {
            os << result.cells_pv_serialize() << "," << result.cost;
        } else {
            os << "infeasible";
        }

        if (result.chunk_start != -1) {
            os << "," << result.chunk_start;
        }

        return os;
    }

};

size_t calc_chebyshev_number_of_chunks();

void update_number_of_chunks(size_t nchunks=number_of_chunks);

int process_input(int argc, char **argv);

#endif //ROBUST_SIZING_PARAMS_MULTIROOF_H
