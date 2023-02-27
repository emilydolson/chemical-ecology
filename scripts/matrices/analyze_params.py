import json
import ast
import csv
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sb

from klemm_lfr_graph import create_matrix as klemm_create_matrix
from klemm_eguiliuz import create_matrix as klemm_random_create_matrix
from motifs import create_matrix as motif_create_matrix
from search_params import read_mangal, random_create_matrix


def visualize_decision_tree():
    return


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
    if len(params[0]) > 12:
        print('too many parameters to visualize in histogram figure')
        return
    
    figure, axis = plt.subplots(5, 2, figsize=(10,10))
    
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


def individual_correlation(params, fitnesses, score, param_id, matrix_scheme, adaptive, replicate=0):
    figure, axis = plt.subplots(1, 1)

    if adaptive:
        scores = [y[score] for y in fitnesses[replicate] if y['Biomass'] > 5000 and y[score] > 0]
        params_id = [params[replicate][j][param_id] for j in range(len(params[replicate])) if fitnesses[replicate][j]['Biomass'] > 5000 and fitnesses[replicate][j][score] > 0]
        colors = [y['Biomass'] for y in fitnesses[replicate] if y['Biomass'] > 5000 and y[score] > 0]
        cmap = 'Oranges'
    else:
        scores = [y[score] for y in fitnesses[replicate]]
        params_id = [params[replicate][j][param_id] for j in range(len(params[replicate]))]
        colors_temp = [y['Biomass'] for y in fitnesses[replicate]]
        colors = []
        min_c = min(colors_temp)
        max_c = max(colors_temp)
        for c in colors_temp:
            if c < 0:
                colors.append(-1)
            elif c > 0 and c < 5000:
                colors.append(-0.5)
            else:
                colors.append((c - min_c) / (max_c - min_c))
        cmap = 'PuOr_r'

    sc = axis.scatter(params_id, scores, c=colors, cmap=cmap)
    cb = figure.colorbar(sc, ax=axis)
    cb.set_label('Biomass')
    figure.suptitle(f'{matrix_scheme} param {param_id} and {score} (adaptive = {adaptive})')
    figure.supxlabel(f'Parameter {param_id}')
    figure.supylabel(score)

    if adaptive:
        plt.savefig(f'plots/{matrix_scheme}/{score}_{param_id}_adpt_correlation.png')
    else:
        plt.savefig(f'plots/{matrix_scheme}/{score}_{param_id}_all_correlation.png')
    plt.close()


def fitness_correlation_custom(params, fitnesses, score, matrix_scheme):
    figure, axis = plt.subplots(1, 1, figsize=(8,4))
    for i in range(len(params)):
        scores = [y[score] for y in fitnesses[i] if y['Biomass'] > 5000]
        params_i = [params[i][j] for j in range(len(params[i])) if fitnesses[i][j]['Biomass'] > 5000]
        correlations = [np.corrcoef([x[j] for x in params_i], scores)[0,1] for j in range(len(params[0][0]))]
        axis.scatter(correlations, range(len(params[0][0])))
    figure.suptitle(f'Correlation between parameters and {score}')
    figure.supxlabel(f'Correlation with {score}')
    plt.savefig(f'plots/{matrix_scheme}/temp_{score}_correlation.png')
    plt.close()


def fitness_correlation(params, fitnesses, score, matrix_scheme):
    figure, axis = plt.subplots(1, 1, figsize=(8,4))
    for i in range(len(params)):
        scores = [y[score] for y in fitnesses[i]] 
        correlations = [np.corrcoef([x[j] for x in params[i]], scores)[0,1] for j in range(len(params[0][0]))]
        axis.scatter(correlations, range(len(params[0][0])))
    figure.suptitle(f'Correlation between parameters and {score}')
    figure.supxlabel(f'Correlation with {score}')
    #axis.yaxis.set_ticks(list(range(10)), ['diffusion', 'seeding', 'clear', 'num_nodes', 'clique_size', 'clique_linkage', 'muw', 'beta', 'pct_pos_in', 'pct_pos_out'])
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


def get_results(file_location, matrix_scheme):
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
                if matrix_scheme == 'motif':
                    new_matrix = []
                    for i in range(len(matrix)):
                        if i == 0 or i == 4:
                            new_matrix.append(matrix[i])
                        elif i == 2:
                            counts = [0]*14
                            for m in matrix[i]:
                                counts[m] += 1
                            new_matrix += counts
                        elif i == 3:
                            counts = [0]*3
                            for m in matrix[i]:
                                counts[m] += 1
                            new_matrix += counts
                    param_list.append(abiotic + new_matrix)
                else:
                    param_list.append(abiotic + matrix)
                fitness_list.append(fitnesses)
    return param_list, fitness_list


def main(file_name):
    file_path = f'/mnt/gs21/scratch/{os.getlogin()}/chemical-ecology/data/{file_name}'

    params = []
    fitnesses = []
    for f in os.listdir(file_path):
        full_path = os.path.join(file_path, f)
        if os.path.exists(os.path.join(full_path, 'results.txt')):
            param_list, fitness_list = get_results(full_path, file_name)
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
