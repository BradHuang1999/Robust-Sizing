//
// Created by Brad Huang on 8/16/20.
//

#include "cheby_multiroof.h"
#include <unordered_map>
#include <deque>

double lambda2;

void update_chebyshev_params(const valarray<bool>& is_zeros) {
    size_t non_zero_cols = 0;
    for (size_t i = 0; i < n_solars; ++i) {
        if (!is_zeros[i]) {
            ++non_zero_cols;
        }
    }

    double asymptote = (double)(non_zero_cols + 1) / (1 - confidence);
    lambda2 = asymptote * (CHEBYSHEV_BETA + 1);

    double eta = (lambda2 + sqrt(lambda2 * lambda2 - 4 * CHEBYSHEV_BETA)) / (2 * CHEBYSHEV_BETA);
    number_of_chunks = (size_t)(eta + 1);
}

RowVectorXd convert_to_rowvec(
        const SimulationMultiRoofResult& result,
        const valarray<bool>& is_zeros) {

    valarray<double> result_pvs = result.PVs[!is_zeros];
    RowVectorXd ret(result_pvs.size() + 1);
    for (size_t i = 0; i < result_pvs.size(); ++i) {
        ret(i) = result_pvs[i];
    }
    ret(result_pvs.size()) = result.B * kWh_in_one_cell;
    return ret;
}

SimulationMultiRoofResult convert_to_result(
        const RowVectorXd& vec,
        const valarray<bool>& is_zeros) {

    size_t arr_size = vec.size() - 1;
    valarray<double> arr(arr_size);
    for (size_t i = 0; i < arr_size; ++i) {
        arr[i] = vec(i);
    }
    valarray<double> result_arr(is_zeros.size());
    result_arr[!is_zeros] = arr;
    result_arr[is_zeros] = 0;

    return SimulationMultiRoofResult(vec(arr_size) / kWh_in_one_cell, result_arr);
}

MatrixXd convert_simulation_result_to_matrix(
        const vector<SimulationMultiRoofResult>& adagrad_sims,
        const valarray<bool>& is_zeros, bool normalize_battery) {

    if (adagrad_sims.empty()) {
        throw range_error("adagrad_sims is empty");
    }

    size_t adagrad_size_x = adagrad_sims.size();
    size_t adagrad_size_y = adagrad_sims[0].PVs[!is_zeros].size();
    MatrixXd ret(adagrad_size_x, adagrad_size_y + 1);

    for (size_t l = 0; l < adagrad_size_x; ++l) {
        ret.row(l) = convert_to_rowvec(adagrad_sims[l], is_zeros);
    }

    return ret;
}

inline double get_l(
        const RowVectorXd& xi,
        const RowVectorXd& mu_eta,
        const MatrixXd& sigma_eta_inv) {
    RowVectorXd xi_mueta_diff = xi - mu_eta;
    return xi_mueta_diff * sigma_eta_inv * xi_mueta_diff.transpose();
}

RowVectorXd get_cheby_steps(const valarray<bool>& is_zeros) {
    valarray<double> pv_diffs = (pv_maxs - pv_mins) / cheby_num_steps;
    valarray<double> non_zero_pv_diffs = pv_diffs[!is_zeros];
    size_t non_zero_pv_diffs_size = non_zero_pv_diffs.size();
    RowVectorXd ret(non_zero_pv_diffs_size + 1);
    for (size_t l = 0; l < non_zero_pv_diffs_size; ++l) {
        ret(l) = non_zero_pv_diffs[l];
    }
    ret(non_zero_pv_diffs_size) = (cells_max - cells_min) * kWh_in_one_cell / cheby_num_steps;
    return ret;
}

inline string to_string(const RowVectorXd& rv) {
    stringstream ss;
    ss << rv;
    return ss.str();
}

inline bool lt(const RowVectorXd& a, const RowVectorXd& b) {
    RowVectorXd diff = a - b;
    for (size_t l = 0; l < diff.size(); ++l) {
        if (diff(l) >= numeric_limits<double>::epsilon()) {
            return false;
        }
    }
    return true;
}

vector<SimulationMultiRoofResult> get_chebyshev_bound(
        const vector<SimulationMultiRoofResult>& adagrad_sims,
        const valarray<bool>& is_zeros) {

    size_t non_zero_cols = 0;
    for (size_t i = 0; i < n_solars; ++i) {
        if (!is_zeros[i]) {
            ++non_zero_cols;
        }
    }

    MatrixXd sigma = convert_simulation_result_to_matrix(adagrad_sims, is_zeros);

    size_t eta = sigma.rows();
    cout << "eta=" << eta << endl;

    double l2 = (double)(eta * eta - 1) / ((1 - confidence) * eta * eta / (double)(non_zero_cols + 1) - eta);
    cout << "l2=" << l2 << endl;

    RowVectorXd mu_eta = sigma.colwise().mean();
    cout << "mu_eta=" << mu_eta << endl;

    MatrixXd sigma_diff = sigma - mu_eta.replicate(eta, 1);
    MatrixXd sigma_eta = sigma_diff.transpose() * sigma_diff / (eta - 1);
    MatrixXd sigma_eta_inv = sigma_eta.inverse();

    cout << "sigma_eta=" << endl << sigma_eta << endl << endl;

    RowVectorXd cheby_mins = convert_to_rowvec(SimulationMultiRoofResult(cells_min, pv_mins), is_zeros);
    RowVectorXd cheby_maxs = convert_to_rowvec(SimulationMultiRoofResult(cells_max, pv_maxs), is_zeros);
    RowVectorXd cheby_steps = get_cheby_steps(is_zeros);

    // Tabu search for the non-dominated, pareto-efficient boundary
    vector<SimulationMultiRoofResult> ret;
    deque<pair<SimulationMultiRoofResult, bool>> search_q;
    unordered_map<string, bool> seen;

    auto is_outside = [&](const RowVectorXd& lxi) -> pair<bool, bool> {
        string lxi_str = to_string(lxi);
        if (seen.count(lxi_str)) {
//            cout << lxi_str << " -> (cached) " << boolalpha << seen[lxi_str] << endl;
            return {seen[lxi_str], true};
        }
        double lval = get_l(lxi, mu_eta, sigma_eta_inv);
        bool loutside = lval >= l2;
//        cout << lxi_str << " -> " << lval << ", " << boolalpha << loutside << endl;
        seen[lxi_str] = loutside;
        return {loutside, false};
    };

    // First binary search at mu_eta points, increasing one dimension at once
    for (size_t i = 0; i <= non_zero_cols; ++i) {
        double col_L = mu_eta[i], col_U = cheby_maxs(i);
        bool is_viable = false;
        RowVectorXd bxi = mu_eta;

        while (col_U - col_L > cheby_steps(i)) {
            double col_M = (col_L + col_U) / 2;
            bxi(i) = col_M;

            if (is_outside(bxi).first) {
                is_viable = true;
                col_U = col_M;
            } else {
                col_L = col_M;
            }
        }

        if (is_viable) {
            cout << "col " << i << " is viable at " << col_U << endl
                 << "- viable point is " << bxi << endl;
            bxi(i) = col_U;
            search_q.emplace_back(convert_to_result(bxi, is_zeros), true);
        } else {
            cout << "col " << i << " not viable" << endl;
        }
    }

    // by now, xi is the first point to be in the search_q
    // perform tabu search. the criteria is that a point is:
    //   - outside the ellipse
    //   - all of its reduced-by-one neighbours are inside the ellipse
    while (!search_q.empty()) {
        const auto& curr_pair = search_q.front();
        search_q.pop_front();

        SimulationMultiRoofResult curr_result = curr_pair.first;
        RowVectorXd curr_xi = convert_to_rowvec(curr_result, is_zeros);
        bool curr_outside = curr_pair.second;

        bool all_lower_neighbor_inside = true;
        bool all_lower_neighbor_outside = true;
        bool all_higher_neighbor_inside = true;
        bool all_higher_neighbor_outside = true;

        vector<pair<SimulationMultiRoofResult, bool>> neighbors;

        // test for all its neighbors
        for (size_t i = 0; i <= non_zero_cols; ++i) {
            // test for lower neighbor
            RowVectorXd lower_xi = curr_xi;
            lower_xi(i) -= cheby_steps(i);
            if (lt(mu_eta, lower_xi)) {
                auto lower_outside = is_outside(lower_xi);
                if (lower_outside.first) { // outside
                    all_lower_neighbor_inside = false;
                } else { // inside
                    all_lower_neighbor_outside = false;
                }
                if (!lower_outside.second) {
//                    cout << "(^pushed)" << endl;
                    neighbors.emplace_back(convert_to_result(lower_xi, is_zeros), lower_outside.first);
                }
            } else {
                all_lower_neighbor_inside = false;
            }

            // test for higher neighbor
            RowVectorXd higher_xi = curr_xi;
            higher_xi(i) += cheby_steps(i);
            auto higher_outside = is_outside(higher_xi);
            if (higher_outside.first) { // outside
                all_higher_neighbor_inside = false;
            } else { // inside
                all_higher_neighbor_outside = false;
            }
            if (!higher_outside.second) {
//                cout << "(^pushed)" << endl;
                neighbors.emplace_back(convert_to_result(higher_xi, is_zeros), higher_outside.first);
            }
        }

        if (curr_outside && lt(curr_xi, cheby_maxs) && all_lower_neighbor_inside) {
            ret.push_back(convert_to_result(curr_xi, is_zeros));
        }

        if (!(curr_outside && all_lower_neighbor_outside && all_higher_neighbor_outside) &&
            !(!curr_outside && all_lower_neighbor_inside && all_higher_neighbor_inside)) {
            for (auto& n: neighbors) {
                search_q.push_back(n);
            }
        }
    }

    return ret;
}
