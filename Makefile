all: multiroof_test sim snc_lolp snc_eue

multiroof_test:
	g++-9 -std=c++17 -I ~/eigen/ cheby_multiroof.cc params_common.cc params_multiroof.cc run_multiroof_test.cc simulate_multiroof.cc -o bin/multiroof_test

sim: 
	g++ -std=c++14 -O3 run_simulation.cc cheby.cc simulate_system.cc params.cc params_common.cc -o bin/sim

sim_multiroof:
	g++ -std=c++14 -O3 run_simulation_multiroof.cc cheby_multiroof.cc simulate_multiroof.cc params_multiroof.cc params_common.cc -o bin/sim_multiroof

snc_lolp: 
	g++ -std=c++14 -O3 run_snc_lolp.cc snc_lolp_pertrace.cc params.cc params_common.cc -o bin/snc_lolp

snc_eue: 
	g++ -std=c++14 -O3 run_snc_eue.cc snc_eue_pertrace.cc params.cc params_common.cc -o bin/snc_eue

debug: debug_sim debug_sim_multiroof debug_snc_lolp debug_snc_eue

debug_sim: 
	g++ -std=c++14 -O0 -ggdb -D DEBUG run_simulation.cc cheby.cc simulate_system.cc params.cc params_common.cc -o bin/debug/sim

debug_sim_multiroof:
	g++ -std=c++14 -O0 -ggdb -D DEBUG run_simulation_multiroof.cc cheby_multiroof.cc simulate_multiroof.cc params_multiroof.cc params_common.cc -o bin/debug/sim_multiroof

debug_snc_lolp: 
	g++ -std=c++14 -O0 -ggdb -D DEBUG run_snc_lolp.cc snc_lolp_pertrace.cc params.cc params_common.cc -o bin/debug/snc_lolp

debug_snc_eue: 
	g++ -std=c++14 -O0 -ggdb -D DEBUG run_snc_eue.cc snc_eue_pertrace.cc params.cc params_common.cc -o bin/debug/snc_eue

clean: 
	rm bin/sim bin/sim_multiroof bin/snc_lolp bin/snc_eue bin/debug/sim bin/debug/sim_multiroof bin/debug/snc_lolp bin/debug/snc_eue
