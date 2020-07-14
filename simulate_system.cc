#include <cmath>
#include <iostream>
#include <vector>
#include "simulate_system.h"

using namespace std;

// parameters specified for an NMC cell with operating range of 1 C charging and discharging

void update_parameters(double n) {

	num_cells = n;

	a1_intercept = 0.0*num_cells;
	
	a2_intercept = kWh_in_one_cell*num_cells;
	
	alpha_d = a2_intercept*1.0;
	alpha_c = a2_intercept*1.0;
	return;
}

// decrease the applied (charging) power by increments of (1/30) until the power is 
// low enough to avoid violating the upper energy limit constraint.
double calc_max_charging(double power, double b_prev) {

	double step = power/30.0;

	for (double c = power; c >= 0; c -= step) {
		double upper_lim = a2_slope*(c/nominal_voltage_c) + a2_intercept;
		double b = b_prev + c*eta_c*T_u;
		if (b <= upper_lim) {
			return c;
		}
	}
	return 0;
}


// decrease the applied (discharging) power by increments of (1/30) until the power is
// low enough to avoid violating the lower energy limit constraint.
double calc_max_discharging(double power, double b_prev) {

	double step = power/30.0;

	for (double d = power; d >= 0; d -= step) {
		double lower_lim = a1_slope*(d/nominal_voltage_d) + a1_intercept;
		double b = b_prev - d*eta_d*T_u;
		if (b >= lower_lim) {
			return d;
		}
	}
	return 0;
}


// Note: sim_year calls procedures calc_max_charging and calc_max_discharging.
// You could potentially speed up the computation by expanding these functions into sim_year
// to avoid procedure calls in this inner loop.
double sim(vector <double> &load_trace, vector <double> &solar_trace,
           size_t start_index, size_t end_index, double cells, double pv, double b_0)
{
	update_parameters(cells);

	// set the battery
	double b = b_0*cells*kWh_in_one_cell; //0.5*a2_intercept

	int loss_events = 0;

	double load_deficit = 0;
	double load_sum = 0;

    size_t trace_length_solar = solar_trace.size();
    size_t trace_length_load = load_trace.size();

	double c = 0.0;
	double d = 0.0;
	double max_c = 0.0;
	double max_d = 0.0;
    size_t index_t_solar;
    size_t index_t_load;

	for (size_t t = start_index; t < end_index; t++)
	{
		// wrap around to the start of the trace if we hit the end.
		index_t_solar = t % trace_length_solar;
		index_t_load = t % trace_length_load;

		load_sum += load_trace[index_t_load];

		// first, calculate how much power is available for charging, and how much is needed to discharge
		c = fmax(solar_trace[index_t_solar]*pv - load_trace[index_t_load],0);
		d = fmax(load_trace[index_t_load] - solar_trace[index_t_solar]*pv, 0);

		// constrain the power
		max_c = fmin(calc_max_charging(c,b), alpha_c);
		max_d = fmin(calc_max_discharging(d,b), alpha_d);

		b = b + max_c*eta_c*T_u - max_d*eta_d*T_u;

		// if we didnt get to discharge as much as we wanted, there is a loss
		if (max_d < d) {
			loss_events += 1;
			load_deficit += (d - max_d);
		}
	}

	if (metric == 0) {
		// lolp
		return loss_events / (double)(end_index - start_index);
	} else {
		// metric == 1, eue
		return load_deficit / (double)load_sum;
	}
}


// Run simulation for provides solar and load trace to find cheapest combination of
// load and solar that can meet the epsilon target
vector <SimulationResult> simulate(vector <double> &load_trace, vector <double> &solar_trace,
                                   size_t start_index, size_t end_index, double b_0) {

	// first, find the lowest value of cells that will get us epsilon loss when the PV is maximized
	// use binary search
	double cells_U = cells_max;
	double cells_L = cells_min;

	while (cells_U - cells_L > cells_step)
	{
        double mid_cells = (cells_L + cells_U) / 2.0;
        double loss = sim(load_trace, solar_trace, start_index, end_index, mid_cells, pv_max, b_0);

		//cout << "sim result with " << a2_intercept << " kWh and " << pv_max << " pv: " << loss << endl;
		if (loss > epsilon)
		{
			cells_L = mid_cells;
		}
		else
		{
		 	// (loss <= epsilon)
			cells_U = mid_cells;
		}
	}

	// set the starting number of battery cells to be the upper limit that was converged on
	double starting_cells = cells_U;
	double starting_cost = B_inv*starting_cells + PV_inv * pv_max;
	double lowest_feasible_pv = pv_max;

	vector <SimulationResult> curve;
	curve.emplace_back(starting_cells * kWh_in_one_cell, lowest_feasible_pv, starting_cost);

	for (double cells = starting_cells; cells <= cells_max; cells += cells_step)
	{
		// for each value of cells, find the lowest pv that meets the epsilon loss constraint
        bool binary_search = true;

		if (curve.size() >= 2)
		{
            double lastC1 = curve.end()[-1].C;
            double lastC2 = curve.end()[-2].C;

            if (lastC1 - lastC2 < 10 * pv_step)
            {
                // use linear search if last two pv values are close
                binary_search = false;
            }
        }

		if (binary_search)
        {
            double pv_U = lowest_feasible_pv;
            double pv_L = 0;

            while (pv_U - pv_L > pv_step)
            {
                double mid_pv = (pv_L + pv_U) / 2.0;
                double loss = sim(load_trace, solar_trace, start_index, end_index, cells, mid_pv, b_0);

                if (loss > epsilon)
                {
                    pv_L = mid_pv;
                }
                else
                {
                    // (loss <= epsilon)
                    pv_U = mid_pv;
                }
            }

            lowest_feasible_pv = pv_U;
        }
		else
        {
		    while (true)
		    {
		        double loss = sim(load_trace, solar_trace, start_index, end_index, cells, lowest_feasible_pv - pv_step, b_0);

		        if (loss < epsilon)
		        {
		            lowest_feasible_pv -= pv_step;
		        }
		        else
		        {
		            break;
		        }

		        // this only happens if the trace is very short, since the battery starts half full
		        // and can prevent loss without pv for a short time
		        if (lowest_feasible_pv <= 0)
		        {
		            lowest_feasible_pv = 0;
		            break;
		        }
		    }
        }

		double cost = B_inv * cells + PV_inv * lowest_feasible_pv;

		curve.emplace_back(cells * kWh_in_one_cell, lowest_feasible_pv, cost);
	} 

	return curve;
}
