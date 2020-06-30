cmake_minimum_required(VERSION 3.16)
project(Robust_Sizing)

set(CMAKE_CXX_STANDARD 14)

include_directories(.)

add_executable(sim_eue
        cheby.cc
        cheby.h
        common.cc
        common.h
        run_simulation.cc
        simulate_system.cc
        simulate_system.h)

add_executable(sim_lolp
        cheby.cc
        cheby.h
        common.cc
        common.h
        run_simulation.cc
        simulate_system.cc
        simulate_system.h)

add_executable(snc_eue
        cheby.cc
        cheby.h
        common.cc
        common.h
        run_snc_eue.cc
        snc_eue_pertrace.cc
        snc_eue_pertrace.h)

add_executable(snc_lolp
        cheby.cc
        cheby.h
        common.cc
        common.h
        run_snc_lolp.cc
        snc_lolp_pertrace.cc
        snc_lolp_pertrace.h)
