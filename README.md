# Robust Sizing

Statistically robust, simulation & optimization based algorithms for **sizing solar panels and storage** for single and multiple roofs.

See publications for [single-roof algorithms](https://www.cs.stanford.edu/~fiodar/pubs/robust-sizing-journal.pdf) and [multiple-roof algorithms](http://pv.bradhuang.me).

The repo is developed by the [ISS4E lab](https://iss4e.ca/) at the University of Waterloo.
Feel free to reach out or contribute.

### Requirements

Multi-roof algorithms require C++ 17 and g++ 9.

Please also download the [Eigen](https://eigen.tuxfamily.org/) library downloaded on path `~/eigen`.
This can be tweaked in [Makefile](Makefile).

Single-roof algorithms only require C++ 14.

### To compile all programs
```
make
```

### To run programs
[NEW] Multi-roof Simulation method. Currently supports up to 5 roof segments.
```
./bin/multiroof_test f 0 <folder_name> <n_solar> <metric> <epsilon> <conf> <n_days> \
                     <battery_cost> <battery_max> <load_file> \
                     <pv_fixed_cost> <pv_marginal_cost> <pv_max> <pv_file> \
                     [...repeat for more roofs]
```

Simulation method for sizing:
```
./bin/sim <pv_cost> <battery_cost> <pv_max> <battery_max> <metric> \
          <epsilon> <conf> <n_days> <load_file> <pv_file>
```

Stochastic network calculus method for sizing with LOLP target:
```
./bin/snc_lolp <pv_cost> <battery_cost> <pv_max> <battery_max> \
               <epsilon> <conf> <n_days> <load_file> <pv_file>
```

Stochastic network calculus method for sizing with EUE target:
```
./bin/snc_eue <pv_cost> <battery_cost> <pv_max> <battery_max>\
              <epsilon> <conf> <n_days> <load_file> <pv_file>
```

Where:

`n_solar`: (multi-roof only) number of roof segments available for solar installation

`folder_name`: (multi-roof only) folder name, under results/full_search, for intermediate outputs

`pv_cost`: (single-roof only) cost per unit of PV panel

`pv_fixed_cost`: (multi-roof only) fixed cost for installing panels (regardless how many) on the specified roof

`pv_marginal_cost`: (multi-roof only) cost per unit of panel on the specified roof

`battery_cost`: cost per kWh of battery

`pv_max`: max units of PV panels available, in terms of kW

`battery_max`: max kWh of battery available, in terms of kWh

`metric`: 0 for LOLP metric, 1 for EUE metric

`epsilon`: epsilon (for LOLP target) or theta (for EUE target) value

`conf`: confidence interval

`n_days`: number of days in each sample (sample period)

`load_file`: name of load file, values in kW (for example, see [load.txt](example_inputs/load.txt))

`pv_file`: name of PV generation file, values in kW (for example, see [pv.txt](example_inputs/pv.txt))

For example, suppose we wish to run the simulation method with:

- LOLP target of **20%** (or EUE target of **45%**)
- Confidence of **85%** over any **365 day** period.
- The price of PV panels is **$2000/unit**, the fixed cost is **$10000/roof**, the price of battery storage is **$500/kWh**
- There are at most **70 units** of PV panels and **225 kWh** of batteries available
- The names of load and pv files are "example_inputs/load.txt" and "example_inputs/pv.txt" respectively.

Then after downloading the code and going into the directory into which it was downloaded, we would write:
```bash
make

# Multi-roof Simulation EUE
./bin/multiroof_test \
    f \
    0 \
    "rothera_pv=nbh admirals hangar giants" \ # folder_name
    2 \ # n_solar
    0 \ # metric - LOLP
    0.2 \ # epsilon
    0.85 \ # conf
    365 \ # n_days
    500 \ # battery_cost
    225 \ # battery_max
    example_inputs/rothera/load_data_cleaned.txt \ # load_file
    10000 \ # pv_fixed_cost for the FIRST roof
    2000 \ # pv_marginal_cost
    70 \ # pv_max
    example_inputs/rothera/NBH_effective.txt \ # pv_file
    10000 \ # pv_fixed_cost for the SECOND roof
    2000 \ # pv_marginal_cost
    70 \ # pv_max
    example_inputs/rothera/Admirals_effective.txt  # pv_file

# Simulation Method LOLP
./bin/sim 2000 500 70 225 0 0.2 0.85 365 example_inputs/load.txt example_inputs/pv.txt

# SNC Method LOLP
./bin/snc_lolp 2000 500 70 225 0.02 0.85 365 example_inputs/load.txt example_inputs/pv.txt

# Multi-roof Simulation EUE
./bin/multiroof_test \
    f \
    0 \
    "rothera_pv=nbh admirals hangar giants" \ # folder_name
    2 \ # n_solar
    1 \ # metric
    0.45 \ # epsilon
    0.85 \ # conf
    365 \ # n_days
    500 \ # battery_cost
    225 \ # battery_max
    example_inputs/rothera/load_data_cleaned.txt \ # load_file
    10000 \ # pv_fixed_cost for the FIRST roof
    2000 \ # pv_marginal_cost
    70 \ # pv_max
    example_inputs/rothera/NBH_effective.txt \ # pv_file
    10000 \ # pv_fixed_cost for the SECOND roof
    2000 \ # pv_marginal_cost
    70 \ # pv_max
    example_inputs/rothera/Admirals_effective.txt  # pv_file

# Simulation Method EUE
./bin/sim 2000 500 70 225 1 0.45 0.85 365 example_inputs/load.txt example_inputs/pv.txt

# SNC Method EUE
./bin/snc_eue 2000 500 70 225 0.45 0.85 365 example_inputs/load.txt example_inputs/pv.txt
```

We can also take electricity load and pv generation metrics through stdin. Instead of entering the file name in the command, enter `--` followed by number of lines to be taken in stdin. For example,
```bash
./bin/sim 2000 500 70 225 0 0.2 0.85 100 -- 8760 -- 8760
```
would first take 8760 lines of load metrics, then 8760 lines of pv metrics.

The output values will be in the standard ouput. It contains a line with three values in the following order, separated by tabs: # of kWh of battery, # of kW of PV, and total cost of the system.

To output the results into a file, add `> output.file` at the end of the command:
```bash
./bin/sim 2000 500 70 225 0 0.2 0.85 100 example_inputs/load.txt \
          example_inputs/pv.txt > results/output.txt
```

Additional configuration parameters, such as the size of the search space, or the number of samples taken from the input files, can be found in the corresponding .h file of each method, or `params_common.h`.
