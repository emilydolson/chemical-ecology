import json
import ast
import csv
import sys
import os
import numpy as np
import pandas as pd
import networkx as nx

from klemm_lfr_graph import create_matrix as klemm_create_matrix
from klemm_eguiliuz import create_matrix as klemm_random_create_matrix
from motifs import create_matrix as motif_create_matrix
from search_params import read_mangal, random_create_matrix
from matrices import visualize_network


def get_matrix_func(matrix_scheme):
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
    return create_matrix_func


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
    create_matrix_func = get_matrix_func(matrix_scheme)

    properties = []
    for genome in params:
        if matrix_scheme == 'motif':
            matrix = create_matrix_func(*genome[3:8])
        elif matrix_scheme == 'mangal':
            matrix = create_matrix_func(*genome[3:5])
        else:
            matrix = create_matrix_func(*genome[3:-1])
        graph = nx.DiGraph(np.array(matrix))

        num_nodes = len(matrix)
        num_interactions = sum([sum([1 for y in x if y != 0]) for x in matrix])
        if num_interactions == 0:
            properties.append([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
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
        num_neg_stronger = 0
        num_pairs = 0
        for i in range(num_nodes):
            for j in range(num_nodes):
                if i != j:
                    if matrix[i][j] < 0 and matrix[j][i] < 0:
                        double_neg_motif += 1
                    if matrix[i][j] > 0 and matrix[j][i] > 0:
                        double_pos_motif += 1
                    if matrix[i][j] < 0 and matrix[j][i] > 0:
                        if abs(matrix[i][j]) > abs(matrix[j][i]):
                            num_neg_stronger += 1
                    if matrix[j][i] < 0 and matrix[i][j] > 0:
                        if abs(matrix[i][j]) < abs(matrix[j][i]):
                            num_neg_stronger += 1
                    if matrix[i][j] != 0 and matrix[j][i] != 0:
                        num_pairs += 1
        if num_pairs > 0:
            double_pos_motif = double_pos_motif / num_pairs
            double_neg_motif = double_neg_motif / num_pairs
            num_neg_stronger = num_neg_stronger / num_pairs
        else:
            double_pos_motif = 0
            double_neg_motif = 0
            num_neg_stronger = 0

        #cluster coefficient
        cluster_coeff = sum(dict(nx.clustering(graph)).values()) / num_nodes

        #weak interactions
        weak = sum([sum([1 for y in x if y < 0.25 and y > -0.25 and y != 0]) for x in matrix]) / num_interactions

        #mutualistic interactions
        mutual = sum([sum([1 for y in x if y > 0]) for x in matrix]) / num_interactions

        properties.append(genome[0:3] +
                          [connectance, in_heterogeneity, out_heterogeneity, diameter,
                          double_pos_motif, double_neg_motif, cluster_coeff, weak, mutual, num_neg_stronger] +
                          [genome[-1]])
    return properties


def read_results(file_location, matrix_scheme, replicate):
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
                elif matrix_scheme == 'mangal':
                    matrix_num = int(matrix[-1].split('_')[0])
                    new_matrix = [matrix_num, 0, 0, 0, 0]
                    if matrix_num == 89:
                        new_matrix[1] = (int(matrix[-1].split('_')[-1][0]))
                    elif matrix_num == 1464:
                        new_matrix[2] = (int(matrix[-1].split('_')[-1][0]))
                    elif matrix_num == 1473:
                        new_matrix[3] = (int(matrix[-1].split('_')[-1][0]))
                    elif matrix_num == 1512:
                        new_matrix[4] = (int(matrix[-1].split('_')[-1][0]))
                    param_list.append(abiotic + matrix + new_matrix + [replicate])
                else:
                    param_list.append(abiotic + matrix + [replicate])
                fitness_list.append(fitnesses)
    return param_list, fitness_list


def get_data(file_name, include_network=True):
    file_path = f'/mnt/gs21/scratch/{os.getlogin()}/chemical-ecology/data/{file_name}'

    abiotic_params = ['diffusion', 'seeding', 'clear']
    if file_name == 'klemm':
        param_names = ['num_nodes', 'clique_size', 'clique_linkage', 'muw', 'beta', 'pct_pos_in', 'pct_pos_out', 'seed']
        analysis_param_names = param_names[1:-1]
    elif file_name == 'klemm_random':
        param_names = ['num_nodes', 'clique_size', 'clique_linkage', 'seed']
        analysis_param_names = param_names[1:-1]
    elif file_name == 'mangal':
        param_names = ['num_nodes', 'name'] + ['matrix_num', 'matrix89', 'matrix1464', 'matrix1473', 'matrix1512']
        analysis_param_names = param_names[2:]
    elif file_name == 'motif':
        param_names = ['num_nodes', 'chunk_size', 'motifs', 'strengths', 'seed'] + [f'motif{i}' for i in range(14)] + [f'strength{i}' for i in range(3)]
        analysis_param_names = param_names[5:]
    elif file_name == 'random':
        param_names = ['num_nodes', 'prob_edge', 'matrix_seed', 'seed']
        analysis_param_names = [param_names[1]]
    else:
        print('undefined matrix scheme')
        exit()
    network_param_names = ['connectance', 'in_heterogeneity', 'out_heterogeneity', 'diameter', \
                           'double_pos_motifs', 'double_neg_motifs', 'cluster_coeff', 'weak', 'mutual', 'num_neg_stronger']

    params = []
    fitnesses = []
    network_params = []
    for f in os.listdir(file_path):
        full_path = os.path.join(file_path, f)
        if os.path.exists(os.path.join(full_path, 'results.txt')):
            param_list, fitness_list = read_results(full_path, file_name, int(f))
            params.append(param_list)
            fitnesses.append(fitness_list)
            if include_network:
                network_params.append(get_network_properties(param_list, file_name))
        else:
            print(f'results not found for replicate {f}')
        print(f'Finished reading replicate {f}')
        #break
    df_generation = params_to_dataframe(fitnesses, params, abiotic_params+param_names)
    print(f'Read {len(df_generation)} samples')
    
    if include_network:
        df_network = params_to_dataframe(fitnesses, network_params, abiotic_params+network_param_names)
        return df_generation, df_network, abiotic_params, param_names, analysis_param_names, network_param_names
    else:
        return df_generation, abiotic_params, param_names, analysis_param_names