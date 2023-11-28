from common import get_matrix_function, get_common_columns, get_scheme_bounds
from scipy.stats import qmc
import networkx as nx
import pandas as pd
import numpy as np
import subprocess
import inspect
import csv
import sys


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


def is_valid(matrix):
    graph = nx.DiGraph.reverse(nx.DiGraph(np.array(matrix)))
    return nx.is_weakly_connected(graph) and not nx.is_empty(graph)

    
def write_matrix(interactions, output_name):
    with open(output_name, 'w') as f:
        wr = csv.writer(f)
        wr.writerows(interactions)


def search_params(scheme, rep):
    matrix_file_name = 'interaction_matrix.dat'
    num_samples = 100

    scheme_name = scheme.__name__
    scheme_args = inspect.getfullargspec(scheme)[0]
    header = get_common_columns() + scheme_args[1:-1] #remove ntypes and seed
    with open('results.txt', 'w') as f:
        f.write(f"{','.join(header)}\n")

    scheme_param_bounds = get_scheme_bounds(scheme_name)
    lower_bounds = [0, 0, 0] + scheme_param_bounds[0] #abiotic + ntypes + matrix
    upper_bounds = [1, 1, 0.5] + scheme_param_bounds[1]
    ints = [False, False, False] + scheme_param_bounds[2]

    samples = sample_params(num_samples, lower_bounds, upper_bounds, ints, rep)
    results = []
    for ntypes in [5, 10, 15, 20, 25]:
        for world_size in [10]:
            for i in range(len(samples)):
                sample = samples[i]

                diffusion = sample[0]
                seeding = sample[1]
                clear = sample[2]
                matrix_params = [ntypes] + sample[3:] + [rep]

                matrix = scheme(*matrix_params)
                if is_valid(matrix):
                    write_matrix(matrix, matrix_file_name)

                    chem_eco = subprocess.Popen(
                        [(f'./chemical-ecology '
                        f'-DIFFUSION {diffusion} '
                        f'-SEEDING_PROB {seeding} '
                        f'-PROB_CLEAR {clear} '
                        f'-INTERACTION_SOURCE {matrix_file_name} '
                        f'-SEED {rep} '
                        f'-MAX_POP {10000} '
                        f'-WORLD_WIDTH {world_size} '
                        f'-WORLD_HEIGHT {world_size} '
                        f'-UPDATES {1000} '
                        f'-N_TYPES {ntypes}')],
                        shell=True, 
                        stdout=subprocess.DEVNULL)
                    return_code = chem_eco.wait()
                    if return_code != 0:
                        print("Error in a-eco, return code:", return_code)
                        sys.stdout.flush()
                        continue
                    
                    df = pd.read_csv('output/recorded_communities_scores_pwip.csv')
                    score = df['logscore'][0]
                    num_communities = len(df)

                    results.append([scheme_name, rep, score, num_communities, world_size, ntypes] + sample)
                else:
                    results.append([scheme_name, rep, 'NA', 0, world_size, ntypes] + sample)

                if (i+1) % 100 == 0 or (i+1) == num_samples:
                    with open('results.txt', 'a') as f:
                        [f.write(f"{','.join(map(str, x))}\n") for x in results]
                    results = []


def main():
    if len(sys.argv) == 3:
        try:
            scheme = get_matrix_function(sys.argv[1])
            rep = int(sys.argv[2])
        except:
            print('Please give the valid arguments (scheme, replicate).')
            exit()
        search_params(scheme, rep)
    else:
        print('Please give the valid arguments (scheme, replicate).')
        exit()


if __name__ == '__main__':
    main()