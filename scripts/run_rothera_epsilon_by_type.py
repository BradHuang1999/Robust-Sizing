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
    _, confs, epsilons, num_processes  = argv
    confs = confs.split(',')
    epsilons = epsilons.split(',')
    num_processes = int(num_processes)

    param_binary = 'bin/multiroof_test'
    n_solar = 5
    n_types = 1 << n_solar
    
    os.chdir('..')

    commands = [
        f'{param_binary} f {type_num} "rothera_pv=nbh admirals bonner hangar giants" 5 1 {epsilon} {conf} ' + \
        '365 670 4000 example_inputs/rothera/load_data_cleaned.txt ' + \
        '10000 188.16 450.15 example_inputs/rothera/NBH_effective.txt ' + \
        '10000 188.16 1043.86 example_inputs/rothera/Admirals_effective.txt ' + \
        '10000 188.16 69.25 example_inputs/rothera/Bonner_effective.txt ' + \
        '10000 188.16 1048.84 example_inputs/rothera/Hangar_effective.txt ' + \
        '10000 188.16 317.32 example_inputs/rothera/Giants_effective.txt'
        for type_num, conf, epsilon in itertools.product(range(n_types), confs, epsilons)
    ]

    # print params
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
