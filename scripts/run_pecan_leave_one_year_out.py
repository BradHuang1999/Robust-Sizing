#!/usr/bin/env python3.7

from glob import glob
import os

import re
import itertools

import subprocess
import multiprocessing as mp

from sys import argv

def run_command(cmd):
    print('running:', cmd)
    result = subprocess.check_output(cmd, shell=True)
    return result.decode()

if __name__ == '__main__':
    # parse arguments
    _, load_paths_start, load_paths_end, confs, epsilons, num_processes = argv
    load_paths_start = int(load_paths_start)
    load_paths_end = int(load_paths_end)
    confs = confs.split(',')
    epsilons = epsilons.split(',')
    num_processes = int(num_processes)

    # get load/pv paths
    os.chdir('..')
    load_paths = sorted(glob('example_inputs/pecan/normed/load_[0-9]*'))[load_paths_start:load_paths_end]
    pv_paths = ['example_inputs/pecan/normed/PV_7989.txt',
                'example_inputs/pecan/normed/PV_6423.txt']

    # fixed params
    param_search_mode = 'v'
    param_type_mode = 0
    param_binary = 'bin/multiroof_test'
    param_n_solars = len(pv_paths)
    param_metric = 1 # 1=EUE, 0=LOLP
    param_days_in_chunk = 365
    param_battery_varcost = 500
    param_battery_max = 120
    param_pv_fixcost = 2000
    param_pv_varcost = 500
    param_pv_max = 60

    # get ids
    load_ids = [re.search(r"example_inputs/pecan/normed/load_([0-9]+).txt", path).group(1) for path in load_paths]
    pv_ids = [re.search(r"example_inputs/pecan/normed/PV_([0-9]+).txt", path).group(1) for path in pv_paths]
    loads = list(zip(load_ids, load_paths))

    # generate commands
    load_pv_ids = f'load=%s_pv={" ".join(pv_ids)}'
    pv_command = " ".join(f'{param_pv_fixcost} {param_pv_varcost} {param_pv_max} {pv_path}' for pv_path in pv_paths)

    commands = [
        f'{param_binary} {param_search_mode} {param_type_mode} "{load_pv_ids % load_id}" {param_n_solars} {param_metric} {epsilon} {conf} {param_days_in_chunk} {param_battery_varcost} {param_battery_max} {load_path} {pv_command}' \
        for (load_id, load_path), conf, epsilon in itertools.product(loads, confs, epsilons)
    ]

    # print params
    print('Load ID start:', load_paths_start, ', end (non-inclusive):', load_paths_end)
    print()
    print('Load IDs being run:')
    print('\t'.join(load_ids))
    print()
    print('Confidence Intervals:')
    print('\t'.join(confs))
    print()
    print('Epsilons:')
    print('\t'.join(epsilons))
    print()
    print('Sample command:')
    print(commands[0])
    print()

    pool = mp.Pool(processes=num_processes)
    results = pool.map(run_command, commands)
    pool.close()
    pool.join()

    for command, result in zip(commands, results):
        print(command, end='\n\n')
        print(result, end='\n\n')
