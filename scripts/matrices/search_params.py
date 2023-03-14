import subprocess
import random
import sys
import csv
import os
import pandas as pd
import numpy as np
from scipy.stats import qmc
from itertools import product

from matrices import write_matrix
from klemm_lfr_graph import create_matrix as klemm_create_matrix
from klemm_eguiliuz import create_matrix as klemm_random_create_matrix
from motifs import create_matrix as motif_create_matrix


#upper_bounds is not inclusive for integers
def sample_params(num_samples, lower_bounds, upper_bounds, ints, constant_ntypes, seed, sampling_scheme):
    if sampling_scheme == 'lhs':
        lower_bounds = [0, 0, 0] + lower_bounds
        upper_bounds = [1, 0.1, 1] + upper_bounds

        sampler = qmc.LatinHypercube(d=len(lower_bounds), seed=seed)
        unscaled_sample = sampler.random(n=num_samples)
        sample = qmc.scale(unscaled_sample, lower_bounds, upper_bounds).tolist()

        matrix_sample = [s[3:] for s in sample]
        abiotic_params = [[round(x, 3) for x in s[0:3]] for s in sample]
        if constant_ntypes:
            matrix_params = [[9]+[int(s[i]) if ints[i] else round(s[i], 3) for i in range(len(s))]+[seed] for s in matrix_sample]
        else:
            matrix_params = [[int(s[i]) if ints[i] else round(s[i], 3) for i in range(len(s))]+[seed] for s in matrix_sample]
    
    return abiotic_params, matrix_params


def sample_matrix_params(num_samples, lower_bounds, upper_bounds, ints, constant_ntypes, seed, sampling_scheme):
    if sampling_scheme == 'lhs':
        sampler = qmc.LatinHypercube(d=len(lower_bounds), seed=seed)
        unscaled_sample = sampler.random(n=num_samples)
        sample = qmc.scale(unscaled_sample, lower_bounds, upper_bounds).tolist()
        if constant_ntypes:
            params = [[9]+[int(s[i]) if ints[i] else round(s[i], 3) for i in range(len(s))]+[seed] for s in sample]
        else:
            params = [[int(s[i]) if ints[i] else round(s[i], 3) for i in range(len(s))]+[seed] for s in sample]
    return params


def get_mangal_params():
    params = []
    for f in os.listdir('mangal_matrices'):
        species_count = sum(1 for line in open(f'mangal_matrices/{f}'))
        if species_count <= 20:
            params.append([species_count, f])
    return params


def read_mangal(num_types, name):
    matrix = np.loadtxt(f'mangal_matrices/{name}', delimiter=',').tolist()
    return matrix


def random_create_matrix(num_types, prob_edge, matrix_seed, seed):
    random.seed(matrix_seed)
    matrix = [[0 for _ in range(num_types)] for _ in range(num_types)]
    for row in range(len(matrix)):
        for col in range(len(matrix)):
            if random.random() < prob_edge:
                matrix[row][col] = random.uniform(-1, 1)
    return matrix


def search_params(matrix_scheme, seed):
    matrix_scheme = matrix_scheme[:-1]
    interaction_matrix_file = 'interaction_matrix.dat'
    sampling_scheme = 'lhs'
    num_samples = 50000

    if matrix_scheme == 'klemm':
        abiotic_params, matrix_params = sample_params(num_samples, [2, 0, 0, 0, 0, 0], [9, 1, 3, 1, 1, 1], [True, False, False, False, False, False], True, int(seed), sampling_scheme)
        create_matrix_func = klemm_create_matrix
    elif matrix_scheme == 'klemm_random':
        abiotic_params, matrix_params = sample_params(num_samples, [2, 0], [9, 1], [True, False], True, int(seed), sampling_scheme)
        create_matrix_func = klemm_random_create_matrix
    elif matrix_scheme == 'mangal':
        matrix_params = get_mangal_params()
        abiotic_params = list(product([0.001, 0.333, 0.666, 0.999], repeat=3))*len(matrix_params)
        matrix_params = matrix_params*64
        num_samples = len(matrix_params)
        create_matrix_func = read_mangal
    elif matrix_scheme == 'motif':
        abiotic_params, matrix_params = sample_params(num_samples, [0]*18, ([14]*9)+([3]*9), [True]*18, True, int(seed), sampling_scheme)
        matrix_params = [[x[0], 3, x[1:10], x[10:-1], x[-1]] for x in matrix_params]
        create_matrix_func = motif_create_matrix
    elif matrix_scheme == 'random':
        abiotic_params, matrix_params = sample_params(num_samples, [0, 0], [1, 100000000], [False, True], True, int(seed), sampling_scheme)
        create_matrix_func = random_create_matrix
    else:
        print('invalid matrix scheme')
        exit()

    fitnesses = []
    abiotics = []
    matrices = []
    for i in range(num_samples):
        diffusion = abiotic_params[i][0]
        seeding = abiotic_params[i][1]
        clear = abiotic_params[i][2]
        param_list = matrix_params[i]

        if seeding > 0.1:
            continue
        
        write_matrix(create_matrix_func(*param_list), interaction_matrix_file)
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
            f'-N_TYPES {param_list[0]}')],
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
        num_communities = df['Num_Final_Communities'].values
        df_aeco = pd.read_csv('a-eco_data.csv')
        biomasses = df_aeco['mean_Biomass'].values
        growth_rates = df_aeco['mean_Growth_Rate'].values
        heredities = df_aeco['mean_Heredity'].values
        fitnesses.append({
            'Biomass': biomasses[-1] - biomasses[0],
            'Growth_Rate': growth_rates[-1] - growth_rates[0],
            'Heredity': heredities[-1],
            'Biomass_Score': b_score[0] / a_score[0] if a_score[0] != 0 else 0,
            'Growth_Rate_Score': g_score[0] / a_score[0] if a_score[0] != 0 else 0,
            'Heredity_Score': h_score[0] / a_score[0] if a_score[0] != 0 else 0,
            'Invasion_Ability_Score': i_score[0] / a_score[0] if a_score[0] != 0 else 0,
            'Resilience_Score': r_score[0] / a_score[0] if a_score[0] != 0 else 0,
            'Final_Communities_Present': df['Assembly_Final_Communities'].values[0] / num_communities[0] if num_communities[0] != 0 else 0
        })
        abiotics.append([diffusion, seeding, clear])
        matrices.append(param_list)

        if (i+1) % 100 == 0 or (i+1) == num_samples:
            with open('results.txt', 'a') as f:
                for i in range(len(fitnesses)):
                    f.write(f'{fitnesses[i]} | {abiotics[i]} | {matrices[i]}\n')
            fitnesses = []
            abiotics = []
            matrices = []


def main():
    search_params(sys.argv[1], sys.argv[2])


if __name__ == '__main__':
    main()