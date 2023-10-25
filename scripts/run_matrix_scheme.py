import pandas as pd
import numpy as np
import subprocess
import inspect
import random
import csv
import sys
import os
from scipy.stats import qmc


def sample_params(num_samples, lower_bounds, upper_bounds, ints, seed):
    ''' Latin Hypercube Sampling of params
    Arguments:
        num_samples (int): how many samples to take
        lower_bounds (list): the lower bound of each param
        upper_bounds (list): the upper bound of each param, not inclusive
        ints (list of bools): whether each param should be an int or a float
        seed: seed for the sample
    Returns:
        a list of lists, outer list is num_samples long and inner list is length of params
    '''

    if (len(lower_bounds) != len(upper_bounds) or len(lower_bounds) != len(ints)):
        print('Error in sample_params input')
        print(f'\tlower_bounds: {lower_bounds}\n\tupper_bounds: {lower_bounds}\n\tints: {lower_bounds}')
        exit()

    sampler = qmc.LatinHypercube(d=len(lower_bounds), seed=seed)
    unscaled_sample = sampler.random(n=num_samples)
    sample = qmc.scale(unscaled_sample, lower_bounds, upper_bounds).tolist()

    params = [[int(s[i]) if ints[i] else round(s[i], 3) for i in range(len(s))] for s in sample]
    return params


def get_scheme_bounds(scheme_name):
    ''' Record and access the bounds of parameters/arguments for each matrix generation scheme
    Arguments:
        scheme_name (string): the name of the matrix scheme function
    Returns:
        a list consisting of a list of lower_bounds, a list of upper_bounds, and a list of bools
        list of bools marks if the param should be sampled as an int or float
    '''

    bounds = {
        'random_matrix': [[0], [1], [False]]
    }
    return bounds[scheme_name]


def random_matrix(ntypes, prob_edge, seed=0):
    random.seed(seed)
    return [[round(random.uniform(-1, 1), 3) if random.random() < prob_edge else 0 for _ in range(ntypes)] for _ in range(ntypes)]


def write_matrix(interactions, output_name):
    with open(output_name, 'w') as f:
        wr = csv.writer(f)
        wr.writerows(interactions)


def search_params(scheme, rep):
    matrix_file_name = 'interaction_matrix.dat'
    num_samples = 10

    scheme_name = scheme.__name__
    scheme_args = inspect.getfullargspec(scheme)[0]
    header = ['scheme', 'replicate', 'score', 'diffusion', 'seeding', 'clear'] + scheme_args[:-1] #remove seed
    with open('results.txt', 'w') as f:
        f.write(f"{','.join(header)}\n")

    scheme_param_bounds = get_scheme_bounds(scheme_name)
    lower_bounds = [0, 0, 0] + [9] + scheme_param_bounds[0] #abiotic + ntypes + matrix
    upper_bounds = [1, 1, 1] + [10] + scheme_param_bounds[1]
    ints = [False, False, False] + [True] + scheme_param_bounds[2]

    samples = sample_params(num_samples, lower_bounds, upper_bounds, ints, rep)
    results = []
    for i in range(len(samples)):
        sample = samples[i]

        diffusion = sample[0]
        seeding = sample[1]
        clear = sample[2]
        ntypes = sample[3]
        matrix_params = sample[3:] + [rep]
        
        matrix = scheme(*matrix_params)
        write_matrix(matrix, matrix_file_name)

        chem_eco = subprocess.Popen(
            [(f'./chemical-ecology '
            f'-DIFFUSION {diffusion} '
            f'-SEEDING_PROB {seeding} '
            f'-PROB_CLEAR {clear} '
            f'-INTERACTION_SOURCE {matrix_file_name} '
            f'-SEED {rep} '
            f'-MAX_POP {10000} '
            f'-WORLD_WIDTH {10} '
            f'-WORLD_HEIGHT {10} '
            f'-UPDATES {1000} '
            f'-N_TYPES {ntypes}')],
            shell=True, 
            stdout=subprocess.DEVNULL)
        return_code = chem_eco.wait()
        if return_code != 0:
            print("Error in a-eco, return code:", return_code)
            sys.stdout.flush()
        
        df = pd.read_csv('output/world_summary_pwip.csv')
        df = df.loc[df['update'] == 1000]
        df = df.replace('ERR', 0)
        score = np.prod(df['proportion']*pd.to_numeric(df['adaptive_assembly_ratio']))

        result = [scheme_name, rep, score] + sample
        results.append(result)

        if (i+1) % 100 == 0 or (i+1) == num_samples:
            with open('results.txt', 'a') as f:
                [f.write(f"{','.join(map(str, x))}\n") for x in results]
            results = []


def main():
    schemes = [random_matrix]

    if len(sys.argv) == 3:
        try:
            scheme = schemes[int(sys.argv[1])]
            rep = int(sys.argv[2])
        except:
            print('Please give the valid arguments (scheme id, replicate).')
            print('Scheme id options:')
            for i in range(len(schemes)):
                print(f'\t{i}: {schemes[i]}')
            exit()
        search_params(scheme, rep)
    else:
        print('Please give the valid arguments (scheme id, replicate).')
        exit()


if __name__ == '__main__':
    main()