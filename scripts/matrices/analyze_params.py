import json
import ast
import csv
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sb


def correlation_heatmap(adaptive, score, matrix_scheme):
    df = pd.DataFrame(adaptive[score])
    corr = df.corr()
    fig = sb.heatmap(corr, center=0, vmin=-1, vmax=1, cmap='coolwarm').get_figure()
    fig.savefig(f'plots/{matrix_scheme}/{score}_heatmap.png')
    plt.close()


def param_histograms(adaptive, score, matrix_scheme):
    figure, axis = plt.subplots(5, 2, figsize=(10,10))

    if len(adaptive[score][0]) > 12:
        print('too many parameters to visualize in histogram figure')
        return

    row = 0
    col = 0
    for i in range(len(adaptive[score][0])):
        param_i = []
        for j in range(len(adaptive[score])):
            param_i.append(adaptive[score][j][i])
        bins = sorted(list(set(param_i)))
        h = np.histogram(param_i, bins=[bins[0]-0.001]+[x+0.001 for x in bins])
        axis[row][col].bar(range(len(bins)), h[0], label=i, alpha=0.66, width=1)
        axis[row][col].set_xticks(np.arange(0, len(bins), 1), bins)
        axis[row][col].set_ylim(0, len(adaptive[score]))
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
                    adaptive[score].append(params[i][j][:-1])
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
    print(f'{len(params[0])*len(params)} combinations')
    adaptive = adaptive_params(fitnesses, params)
    for score in [x for x in fitnesses[0][0].keys() if 'Score' in x]:
        print(f'Adaptive {score}: {len(adaptive[score])}')
        fitness_correlation(params, fitnesses, score, file_name)
        param_histograms(adaptive, score, file_name)
        correlation_heatmap(adaptive, score, file_name)


if __name__ == '__main__':
    main(sys.argv[1])
