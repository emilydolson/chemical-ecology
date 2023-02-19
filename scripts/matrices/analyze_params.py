import json
import ast
import csv
import sys
import os
import numpy as np
import matplotlib.pyplot as plt


def fitness_correlation(params, fitnesses, score, matrix_scheme):
    figure, axis = plt.subplots(1, 1)
    for i in range(len(params)):
        scores = [y[score] for y in fitnesses[i]]
        correlations = [np.corrcoef([x[j] for x in params[i]], scores)[0,1] for j in range(len(params[0][0]))]
        axis.scatter(correlations, range(len(params[0][0])))
    figure.suptitle(f'Correlation between parameters and {score}')
    figure.supxlabel(f'Correlation with {score}')
    plt.savefig(f'plots/{matrix_scheme}/{score}_correlation.png')


def is_adaptive(fitnesses):
    adaptiveness = {x:[] for x in fitnesses[0].keys() if 'Score' in x}
    for fitness in fitnesses:
        for score in adaptiveness.keys():
            if fitness[score] > 0 and fitness['Biomass'] > 5000:
                adaptiveness[score].append(True)
            else:
                adaptiveness[score].append(False)
    return adaptiveness


def get_results(file_location):
    results = open(f'{file_location}/results.txt').read()
    param_list = []
    fitness_list = []
    for line in results.split('\n'):
        if len(line) > 0:
            if line[0] == '{':
                fitnesses, abiotic, matrix = line.split(' | ')
                fitnesses = json.loads(fitnesses.replace('\'', '\"'))
                abiotic = ast.literal_eval(abiotic)
                matrix = ast.literal_eval(matrix)
                param_list.append(abiotic + matrix)
                fitness_list.append(fitnesses)
    return param_list, fitness_list


def main(file_name):
    file_path = f'/mnt/gs21/scratch/{os.getlogin()}/chemical-ecology/data/{file_name}'
    params = []
    fitnesses = []
    adaptiveness = []
    for f in os.listdir(file_path):
        full_path = os.path.join(file_path, f)
        if os.path.exists(os.path.join(full_path, 'results.txt')):
            param_list, fitness_list = get_results(full_path)
            params.append(param_list)
            fitnesses.append(fitness_list)
            adaptiveness.append(is_adaptive(fitness_list))
        else:
            print(f'results not found for seed {f}')
    fitness_correlation(params, fitnesses, 'Biomass', file_name)
    

if __name__ == '__main__':
    main(sys.argv[1])
