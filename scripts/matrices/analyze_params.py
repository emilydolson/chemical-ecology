import json
import ast
import csv
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sb
import networkx as nx
from scipy.stats import kstest, chisquare

from klemm_lfr_graph import create_matrix as klemm_create_matrix
from klemm_eguiliuz import create_matrix as klemm_random_create_matrix
from motifs import create_matrix as motif_create_matrix
from search_params import read_mangal, random_create_matrix


#https://towardsdatascience.com/how-to-compare-two-or-more-distributions-9b06ee4d30bf
def cdf(df, score, param_names, matrix_scheme):
    if len(param_names) <= 10:
        figure, axis = plt.subplots(5, 2, figsize=(10,10))
    elif len(param_names) > 10 and len(param_names) <= 15:
        figure, axis = plt.subplots(5, 3, figsize=(12,12))
    elif len(param_names) > 15 and len(param_names) <= 20:
        figure, axis = plt.subplots(5, 4, figsize=(15,15))
    else:
        print('too many parameters to visualize in histogram figure')
        return

    row = 0
    col = 0
    for parameter in param_names:
        adaptive_params = df.loc[df[score[:-5]+'Adaptive']==True, parameter].values
        successional_params = df.loc[df[score[:-5]+'Adaptive']==False, parameter].values

        df_ks = pd.DataFrame()
        df_ks[parameter] = np.sort(df[parameter].unique())
        df_ks['F_adaptive'] = df_ks[parameter].apply(lambda x: np.mean(adaptive_params<=x))
        df_ks['F_successional'] = df_ks[parameter].apply(lambda x: np.mean(successional_params<=x))
        k = np.argmax(np.abs(df_ks['F_successional'] - df_ks['F_adaptive']))
        ks_stat = np.abs(df_ks['F_adaptive'][k] - df_ks['F_successional'][k])
        stat, p_value = kstest(adaptive_params, successional_params)

        y = (df_ks['F_adaptive'][k] + df_ks['F_successional'][k])/2
        axis[row][col].plot(parameter, 'F_successional', data=df_ks, label=parameter)
        axis[row][col].plot(parameter, 'F_adaptive', data=df_ks, label='adaptive')
        axis[row][col].errorbar(x=df_ks[parameter][k], y=y, yerr=ks_stat/2, color='k',
                    capsize=5, mew=3, label=f'Test staistic: {stat:.3f} \n p-value: {p_value:.3f}')
        axis[row][col].legend(loc='lower right', fontsize='x-small')

        row += 1
        if row % 5 == 0:
            col += 1
            row = 0

    figure.suptitle('Kolmogorov-Smirnov Test')
    figure.supylabel('Cumulative Probability')

    plt.savefig(f'plots/{matrix_scheme}/{score}_cdf.png')
    plt.close()


def visualize_decision_tree():
    return


def correlation_heatmap(df, score, param_names, matrix_scheme):
    df_params = df[param_names]
    corr = df_params.corr()
    fig = sb.heatmap(corr, center=0, vmin=-1, vmax=1, cmap='coolwarm').get_figure()
    fig.savefig(f'plots/{matrix_scheme}/{score}_heatmap.png')
    plt.close()


def param_histograms(df, score, param_names, matrix_scheme):
    if len(param_names) <= 10:
        figure, axis = plt.subplots(5, 2, figsize=(10,10))
    elif len(param_names) > 10 and len(param_names) <= 15:
        figure, axis = plt.subplots(5, 3, figsize=(12,12))
    elif len(param_names) > 15 and len(param_names) <= 20:
        figure, axis = plt.subplots(5, 4, figsize=(15,15))
    else:
        print('too many parameters to visualize in histogram figure')
        return

    row = 0
    col = 0
    for param in param_names:
        params = df[param]
        bins = set(params)
        if len(bins) < 10:
            bins = sorted(list(bins))
            h = np.histogram(params, bins=[bins[0]-0.001]+[x+0.001 for x in bins])
            axis[row][col].bar(range(len(bins)), h[0], label=param, alpha=0.66, width=1)
            axis[row][col].set_xticks(np.arange(0, len(bins), 1), bins)
        else:
            if round(max(params)) == 1:
                axis[row][col].hist(params, bins=np.arange(0, 1.1, 0.1), label=param, alpha=0.66)
            else:
                axis[row][col].hist(params, label=param, alpha=0.66)
        axis[row][col].set_ylim(0, len(df)//2)
        axis[row][col].legend()
        row += 1
        if row % 5 == 0:
            col += 1
            row = 0
    figure.suptitle(f'{matrix_scheme} {score}')
    figure.supxlabel('Value')
    figure.supylabel('Count')

    plt.savefig(f'plots/{matrix_scheme}/{score}_histogram.png')
    plt.close()


def individual_scatter(df, score, param, color, matrix_scheme):
    figure, axis = plt.subplots(1, 1)

    sc = axis.scatter(df[param], df[score], c=df[color], cmap='Oranges')
    cb = figure.colorbar(sc, ax=axis)
    cb.set_label(color)
    axis.set_ylim(df[score].min(), 5)
    figure.suptitle(f'{matrix_scheme} {param} and {score}')
    figure.supxlabel(param)
    figure.supylabel(score)

    plt.savefig(f'plots/{matrix_scheme}/{score}_{param}_scatter.png')
    plt.close()


def fitness_correlation(df, score, param_names, matrix_scheme):
    figure, axis = plt.subplots(1, 1, figsize=(8,4))
    for i in df['replicate'].unique():
        df_replicate = df.loc[df['replicate'] == i]
        scores = df_replicate[score].values
        correlations = [np.corrcoef(df_replicate[param].values, scores)[0,1] for param in param_names]
        axis.scatter(correlations, range(len(param_names)))
    axis.yaxis.set_ticks(list(range(len(param_names))), param_names)
    figure.suptitle(f'Correlation between parameters and {score}')
    figure.supxlabel(f'Correlation with {score}')
    plt.savefig(f'plots/{matrix_scheme}/{score}_correlation.png')
    plt.close()


def params_to_dataframe(fitnesses, params, param_names):
    fitness_names = list(fitnesses[0][0].keys())
    combined_lists = []
    for i in range(len(fitnesses)):
        for j in range(len(fitnesses[i])):
            fitnesses_flattened = []
            adaptive_or_not = []
            for fit in fitness_names:
                fitnesses_flattened.append(fitnesses[i][j][fit])
                if 'Score' in fit:
                    adaptive_or_not.append(True if fitnesses[i][j][fit] > 1 and fitnesses[i][j]['Biomass'] > 5000 else False)
            combined_lists.append(fitnesses_flattened+params[i][j]+adaptive_or_not)

    adaptive_names = [x[:-5]+'Adaptive' for x in fitness_names if 'Score' in x]
    return pd.DataFrame(combined_lists, columns=fitness_names+param_names+['replicate']+adaptive_names)


def get_network_properties(params, matrix_scheme):
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
        print('undefined matrix scheme')
        return

    properties = []
    for genome in params:
        if matrix_scheme == 'motif':
            matrix = create_matrix_func(*genome[3:8])
        else:
            matrix = create_matrix_func(*genome[3:-1])
        graph = nx.DiGraph(np.array(matrix))

        num_nodes = len(matrix)
        num_interactions = sum([sum([1 for y in x if y != 0]) for x in matrix])
        if num_interactions == 0:
            properties.append([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
            continue

        #connectedness
        if matrix_scheme == 'klemm' or matrix_scheme == 'klemm_random': #no intrinstic interactions
            connectance = num_interactions / (num_nodes*(num_nodes-1))
            hetero_range = range(1,num_nodes)
        else:
            connectance = num_interactions / num_nodes**2
            hetero_range = range(1,num_nodes+1)
        
        #heterogeneity
        in_degree_k, in_degree_count = np.unique(list(d for n, d in graph.in_degree()), return_counts=True)
        out_degree_k, out_degree_count = np.unique(list(d for n, d in graph.out_degree()), return_counts=True)
        in_degree_P = {in_degree_k[i]: in_degree_count[i]/num_nodes for i in range(len(in_degree_k))}
        out_degree_P = {out_degree_k[i]: out_degree_count[i]/num_nodes for i in range(len(out_degree_k))}
        in_heterogeneity = np.sqrt(sum([(1-in_degree_P[k])**2 for k in set(in_degree_k)])/num_nodes) \
                         / np.sqrt(sum([(1-in_degree_P[k])**2 if k in in_degree_P else 1 for k in hetero_range])/num_nodes)
        out_heterogeneity = np.sqrt(sum([(1-out_degree_P[k])**2 for k in set(out_degree_k)])/num_nodes) \
                          / np.sqrt(sum([(1-out_degree_P[k])**2 if k in out_degree_P else 1 for k in hetero_range])/num_nodes)

        #network diameter
        shortest_path = dict(nx.shortest_path_length(graph))
        diameter = max([max(shortest_path[i].values()) for i in range(len(shortest_path))])

        #motifs
        double_pos_motif = 0
        double_neg_motif = 0
        for i in range(num_nodes):
            for j in range(num_nodes):
                if i != j:
                    if matrix[i][j] < 0 and matrix[j][i] < 0:
                        double_neg_motif += 1
                    if matrix[i][j] > 0 and matrix[j][i] > 0:
                        double_pos_motif += 1
        double_pos_motif = double_pos_motif / num_interactions
        double_neg_motif = double_neg_motif / num_interactions

        #cluster coefficient
        cluster_coeff = sum(dict(nx.clustering(graph)).values()) / num_nodes

        #weak interactions
        weak = sum([sum([1 for y in x if y < 0.25 and y > -0.25 and y != 0]) for x in matrix]) / num_interactions

        #mutualistic interactions
        mutual = sum([sum([1 for y in x if y > 0]) for x in matrix]) / num_interactions

        properties.append(genome[0:3] +
                          [connectance, in_heterogeneity, out_heterogeneity, diameter,
                          double_pos_motif, double_neg_motif, cluster_coeff, weak, mutual] +
                          [genome[-1]])
    return properties


def get_results(file_location, matrix_scheme, replicate):
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
                        if i == 2:
                            counts = [0]*14
                            for m in matrix[i]:
                                counts[m] += 1
                            new_matrix += counts
                        if i == 3:
                            counts = [0]*3
                            for m in matrix[i]:
                                counts[m] += 1
                            new_matrix += counts
                    param_list.append(abiotic + matrix + new_matrix + [replicate])
                else:
                    param_list.append(abiotic + matrix + [replicate])
                fitness_list.append(fitnesses)
    return param_list, fitness_list


def run_analysis(df, param_names, file_name):
    print(f'{len(df)} combinations')
    for score in [x for x in df.columns if 'Score' in x]:
        df_adaptive = df.loc[df[score[:-5]+'Adaptive'] == True]
        print(f'Adaptive {score}: {len(df_adaptive)}')

        fitness_correlation(df, score, param_names, file_name)
        param_histograms(df_adaptive, score, param_names, file_name)
        correlation_heatmap(df_adaptive, score, param_names, file_name)
        cdf(df, score, param_names, file_name)

    param_histograms(df, 'All', param_names, file_name)
    correlation_heatmap(df, 'All', param_names, file_name)


def scheme_analysis(file_name):
    file_path = f'/mnt/gs21/scratch/{os.getlogin()}/chemical-ecology/data/{file_name}'

    abiotic_params = ['diffusion', 'seeding', 'clear']
    if file_name == 'klemm':
        param_names = ['num_nodes', 'clique_size', 'clique_linkage', 'muw', 'beta', 'pct_pos_in', 'pct_pos_out', 'seed']
        analysis_param_names = param_names[1:-1]
    elif file_name == 'klemm_random':
        param_names = ['num_nodes', 'clique_size', 'clique_linkage', 'seed']
        analysis_param_names = param_names[1:-1]
    elif file_name == 'mangal':
        param_names = ['num_nodes', 'name']
        analysis_param_names = []
    elif file_name == 'motif':
        param_names = ['num_nodes', 'chunk_size', 'motifs', 'strengths', 'seed'] + [f'motif{i}' for i in range(14)] + [f'strength{i}' for i in range(3)]
        analysis_param_names = param_names[5:]
    elif file_name == 'random':
        param_names = ['num_nodes', 'prob_edge', 'matrix_seed', 'seed']
        analysis_param_names = [param_names[1]]
    else:
        print('undefined matrix scheme')
        return
    network_param_names = ['connectance', 'in_heterogeneity', 'out_heterogeneity', 'diameter', \
                           'double_pos_motifs', 'double_neg_motifs', 'cluster_coeff', 'weak', 'mutual']

    params = []
    fitnesses = []
    network_params = []
    for f in os.listdir(file_path):
        full_path = os.path.join(file_path, f)
        if os.path.exists(os.path.join(full_path, 'results.txt')):
            param_list, fitness_list = get_results(full_path, file_name, int(f))
            params.append(param_list)
            fitnesses.append(fitness_list)
            network_params.append(get_network_properties(param_list, file_name))
        else:
            print(f'results not found for replicate {f}')
        print(f'Finished reading replicate {f}')
    df_generation = params_to_dataframe(fitnesses, params, abiotic_params+param_names)
    df_network = params_to_dataframe(fitnesses, network_params, abiotic_params+network_param_names)

    run_analysis(df_generation, abiotic_params+analysis_param_names, file_name)
    run_analysis(df_network, abiotic_params+network_param_names, file_name+'/network')


if __name__ == '__main__':
    argument = sys.argv[1]
    if argument != None:
        scheme_analysis(argument)
    else:
        print('none')
