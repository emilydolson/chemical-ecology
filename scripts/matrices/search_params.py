import subprocess
import pandas as pd
import random
import sys
import csv
import numpy as np

from matrices import write_matrix
from klemm_lfr_graph import create_matrix as klemm_create_matrix


def test_matrix_size(seed): #5
    params = []
    for n_types in [9, 15, 21, 27, 33]:
        params.append([n_types, 4, 0.9, 0.5, 1, 0.75, 0.25, seed])
    return params


def create_klemm_params(seed): #9,600
    params = []
    for clique_size in range(2, 10, 2): #4
        for clique_linkage in np.arange(0, 1.2, 0.2): #5
            for muw in np.arange(0, 1.2, 0.2): #5
                for beta in np.arange(0, 3.5, 0.5): #6
                    for pct_pos_in in np.arange(0, 1.25, 0.25): #4
                        for pct_pos_out in np.arange(0, 1.25, 0.25): #4
                            params.append([9, clique_size, clique_linkage, muw, beta, pct_pos_in, pct_pos_out, seed])
    return params


def search_params(matrix_scheme, seed):
    abiotic_range = [0.001, 0.333, 0.666, 0.999]
    interaction_matrix_file = 'interaction_matrix.dat'

    if matrix_scheme == 'klemm':
        params = create_klemm_params(int(seed))
    elif matrix_scheme == 'test_sizing':
        params = test_matrix_size(int(seed))
    else:
        print('invalid matrix scheme')
        exit()

    for param_list in params:
        write_matrix(klemm_create_matrix(*param_list), interaction_matrix_file)
        fitnesses = []
        abiotics = []
        for diffusion in abiotic_range:
            for seeding in abiotic_range:
                for clear in abiotic_range:
                    chem_eco = subprocess.Popen(
                        [(f'./chemical-ecology '
                        f'-DIFFUSION {diffusion} '
                        f'-SEEDING_PROB {seeding} '
                        f'-PROB_CLEAR {clear} ' 
                        f'-INTERACTION_SOURCE {interaction_matrix_file} '
                        f'-REPRO_THRESHOLD {100000000} '
                        f'-MAX_POP {10000} '
                        f'-WORLD_X {10} '
                        f'-WORLD_Y {10} '
                        f'-UPDATES {1000} '
                        f'-N_TYPES {9}')],
                        shell=True, 
                        stdout=subprocess.DEVNULL)
                    return_code = chem_eco.wait()
                    if return_code != 0:
                        print("Error in a-eco, return code:", return_code)
                        sys.stdout.flush()
                    df = pd.read_csv('scores.csv')
                    b_score = df['Biomass_Score'].values
                    g_score = df['Growth_Rate_Score'].values
                    h_score = df['Heredity_Score'].values
                    i_score = df['Invasion_Ability_Score'].values
                    r_score = df['Resiliance_Score'].values
                    a_score = df['Assembly_Score'].values
                    fitnesses.append({
                        'Biomass_Score': b_score[0] - a_score[0],
                        'Growth_Rate_Score': g_score[0] - a_score[0],
                        'Heredity_Score': h_score[0] - a_score[0],
                        'Invasion_Ability_Score': i_score[0] - a_score[0],
                        'Resiliance_Score': r_score[0] - a_score[0]
                    })
                    abiotics.append([diffusion, seeding, clear])
        with open('results.txt', 'a') as f:
            for i in range(len(fitnesses)):
                f.write(f'{fitnesses[i]} | {abiotics[i]} | {param_list}\n')


def main():
    search_params(sys.argv[1], sys.argv[2])


main()