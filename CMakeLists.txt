cmake_minimum_required(VERSION 3.16)
project(Robust_Sizing)

find_package(Eigen3 3.3 REQUIRED NO_MODULE)

set(CMAKE_C_COMPILER /usr/local/bin/gcc-10)
set(CMAKE_CXX_COMPILER /usr/local/bin/g++-10)
set(CMAKE_CXX_STANDARD 17)

add_compile_definitions(TRACE)

include_directories(.)

add_executable(run_multiroof_cheby_experiments
        cheby_multiroof.cc
        cheby_multiroof.h
        params_common.cc
        params_common.h
        params_multiroof.cc
        params_multiroof.h
        run_multiroof_cheby_experiments.cc
        simulate_multiroof.cc
        simulate_system.h
        )
target_link_libraries(run_multiroof_cheby_experiments Eigen3::Eigen)

add_executable(run_multiroof_test
        cheby_multiroof.cc
        cheby_multiroof.h
        params_common.cc
        params_common.h
        params_multiroof.cc
        params_multiroof.h
        run_multiroof_test.cc
        simulate_multiroof.cc
        simulate_system.h
        )
target_link_libraries(run_multiroof_test Eigen3::Eigen)

add_executable(sim_eue
        cheby.cc
        cheby.h
        params_common.cc
        params_common.h
        params.cc
        params.h
        run_simulation.cc
        simulate_system.cc
        simulate_system.h
        )

add_executable(sim_lolp
        cheby.cc
        cheby.h
        params_common.cc
        params_common.h
        params.cc
        params.h
        run_simulation.cc
        simulate_system.cc
        simulate_system.h
        )

add_executable(sim_eue_multiroof
        cheby_multiroof.cc
        cheby_multiroof.h
        params_common.cc
        params_common.h
        params_multiroof.cc
        params_multiroof.h
        run_simulation_multiroof.cc
        simulate_multiroof.cc
        simulate_system.h
        )
target_link_libraries(sim_eue_multiroof Eigen3::Eigen)

add_executable(sim_lolp_multiroof
        cheby_multiroof.cc
        cheby_multiroof.h
        params_common.cc
        params_common.h
        params_multiroof.cc
        params_multiroof.h
        run_simulation_multiroof.cc
        simulate_multiroof.cc
        simulate_system.h
        )
target_link_libraries(sim_lolp_multiroof Eigen3::Eigen)

add_executable(snc_eue
        cheby.cc
        cheby.h
        params_common.cc
        params_common.h
        params.cc
        params.h
        run_snc_eue.cc
        snc_eue_pertrace.cc
        simulate_system.h
        )

add_executable(snc_lolp
        cheby.cc
        cheby.h
        params_common.cc
        params_common.h
        params.cc
        params.h
        run_snc_lolp.cc
        snc_lolp_pertrace.cc
        simulate_system.h
        )
