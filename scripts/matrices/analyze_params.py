import json
import ast
import csv
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sb


def get_network_properties(adaptive, score, matrix_scheme):
    if matrix_scheme == 'klemm':
        create_matrix_func = klemm_create_matrix
    elif matrix_scheme == 'klemm_random':
        create_matrix_func = klemm_random_create_matrix
    elif matrix_scheme == 'mangal':
        create_matrix_func = read_mangal
    elif matrix_scheme == 'motif':
        create_matrix_func = motif_create_matrix
    elif matrix_scheme == 'random':
        create_matrix_func = random_create_matrix
    else:
        print('invalid matrix scheme')
        return

    properties = []
    for genome in adaptive[score]:
        matrix = create_matrix_func(*genome[3:])
        connectance = sum([sum([1 for y in x if y != 0]) for x in matrix]) / len(matrix)**2


def correlation_heatmap(params, score, matrix_scheme):
    df = pd.DataFrame(params)
    df.drop(df.columns[[-1]], axis=1, inplace=True)
    corr = df.corr()
    fig = sb.heatmap(corr, center=0, vmin=-1, vmax=1, cmap='coolwarm').get_figure()
    fig.savefig(f'plots/{matrix_scheme}/{score}_heatmap.png')
    plt.close()


def param_histograms(params, score, matrix_scheme):
    figure, axis = plt.subplots(5, 2, figsize=(10,10))

    if len(params[0]) > 12:
        print('too many parameters to visualize in histogram figure')
        return

    row = 0
    col = 0
    for i in range(len(params[0])-1):
        param_i = []
        for j in range(len(params)):
            param_i.append(params[j][i])
        bins = set(param_i)
        if len(bins) < 10:
            bins = sorted(list(bins))
            h = np.histogram(param_i, bins=[bins[0]-0.001]+[x+0.001 for x in bins])
            axis[row][col].bar(range(len(bins)), h[0], label=i, alpha=0.66, width=1)
            axis[row][col].set_xticks(np.arange(0, len(bins), 1), bins)
        else:
            if round(max(param_i)) == 1:
                axis[row][col].hist(param_i, bins=np.arange(0, 1.1, 0.1), label=i, alpha=0.66)
            else:
                axis[row][col].hist(param_i, bins=np.arange(0, round(max(param_i))+0.5, 0.5), label=i, alpha=0.66)
        axis[row][col].set_ylim(0, len(params))
        axis[row][col].legend()
        if i > 1:
            row += 1
        if row == 5:
            col += 1
            row = 0
    figure.suptitle(f'{matrix_scheme} {score}')
    figure.supxlabel('Value')
    figure.supylabel('Count')

    plt.savefig(f'plots/{matrix_scheme}/{score}_histogram.png')
    plt.close()


def fitness_correlation(params, fitnesses, score, matrix_scheme):
    figure, axis = plt.subplots(1, 1)
    for i in range(len(params)):
        scores = [y[score] for y in fitnesses[i]]
        correlations = [np.corrcoef([x[j] for x in params[i]], scores)[0,1] for j in range(len(params[0][0]))]
        axis.scatter(correlations, range(len(params[0][0])))
    figure.suptitle(f'Correlation between parameters and {score}')
    figure.supxlabel(f'Correlation with {score}')
    plt.savefig(f'plots/{matrix_scheme}/{score}_correlation.png')
    plt.close()


def adaptive_params(fitnesses, params):
    adaptive = {x:[] for x in fitnesses[0][0].keys() if 'Score' in x}
    for i in range(len(params)):
        for j in range(len(params[i])):
            for score in adaptive.keys():
                if fitnesses[i][j][score] > 0 and fitnesses[i][j]['Biomass'] > 5000:
                    adaptive[score].append(params[i][j])
    return adaptive


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
    adaptive = []
    for f in os.listdir(file_path):
        full_path = os.path.join(file_path, f)
        if os.path.exists(os.path.join(full_path, 'results.txt')):
            param_list, fitness_list = get_results(full_path)
            params.append(param_list)
            fitnesses.append(fitness_list)
        else:
            print(f'results not found for seed {f}')

    print(f'{sum([len(x) for x in params])} combinations')
    adaptive = adaptive_params(fitnesses, params)
    for score in [x for x in fitnesses[0][0].keys() if 'Score' in x]:
        print(f'Adaptive {score}: {len(adaptive[score])}')
        fitness_correlation(params, fitnesses, score, file_name)
        param_histograms(adaptive[score], score, file_name)
        correlation_heatmap(adaptive[score], score, file_name)

    all_params = [val for sublist in params for val in sublist]
    param_histograms(all_params, 'All', file_name)
    correlation_heatmap(all_params, 'All', file_name)


if __name__ == '__main__':
    main(sys.argv[1])
