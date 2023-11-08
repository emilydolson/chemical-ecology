from common import get_matrix_function, get_common_columns, get_raw_data_path, get_processed_data_path
import networkx as nx
import pandas as pd
import numpy as np
import os


def save_data(scheme, df):
    df.to_pickle(f'{get_processed_data_path()}{scheme}.pkl')

    
def get_topological_properties(row, matrix_scheme):
    create_matrix_func = get_matrix_function(matrix_scheme)
    scheme_specific_params = row[~row.index.isin(get_common_columns())].tolist()
    matrix_params = [row['ntypes']] + scheme_specific_params + [row['replicate']]
    matrix = create_matrix_func(*matrix_params)

    num_nodes = len(matrix)
    num_interactions = sum([sum([1 for y in x if y != 0]) for x in matrix])
    
    graph = nx.DiGraph.reverse(nx.DiGraph(np.array(matrix)))
    
    #structural properties
    row['connectance'] = num_interactions / num_nodes**2
    row['cluster_coeff'] = sum(dict(nx.clustering(graph)).values()) / num_nodes

    #weight properties
    num_positive = sum([sum([1 for y in x if y > 0]) for x in matrix])
    num_negative = sum([sum([1 for y in x if y < 0]) for x in matrix])
    sum_positive = sum([sum([y for y in x if y > 0]) for x in matrix])
    sum_negative = sum([sum([y for y in x if y < 0]) for x in matrix])
    row['weak_proportion'] = sum([sum([1 for y in x if y < 0.25 and y > -0.25 and y != 0]) for x in matrix]) / num_interactions
    row['strong_proportion'] = sum([sum([1 for y in x if y > 0.75 or y < -0.75]) for x in matrix]) / num_interactions
    row['pos_proportion'] = num_positive / num_interactions
    row['neg_proportion'] = num_negative / num_interactions
    row['avg_pos_strength'] = sum_positive / num_positive if num_positive > 0 else 0
    row['avg_neg_strength'] = sum_negative / num_negative if num_negative > 0 else 0
    row['avg_strength_diff'] = row['avg_pos_strength'] - row['avg_neg_strength']
    row['mutualistic_pairs'] = sum([sum([1 for j in range(i+1, num_nodes) if matrix[i][j] > 0 and matrix[j][i] > 0]) for i in range(num_nodes)]) / num_interactions
    row['pos_self_strength'] = sum([1 for i in range(num_nodes) if matrix[i][i] > 0]) / num_nodes

    return row


def read_data(file_path, matrix_scheme):
    scheme_dfs = []
    for replicate in os.listdir(file_path+matrix_scheme):
        if os.path.exists(file_path+matrix_scheme+'/'+replicate+'/results.txt'):
            df_rep = pd.read_csv(file_path+matrix_scheme+'/'+replicate+'/results.txt')
            scheme_dfs.append(df_rep)
        else:
            print(f'Results not found for {matrix_scheme} replicate {replicate}.')
    scheme_df = pd.concat(scheme_dfs)
    return scheme_df


def main():
    file_path = get_raw_data_path()
    for matrix_scheme in os.listdir(file_path):
        scheme_df = read_data(file_path, matrix_scheme)
        print(f"{matrix_scheme}: removing {scheme_df['score'].isna().sum()} NA scores")
        scheme_df = scheme_df.dropna()
        scheme_df = scheme_df.apply(get_topological_properties, args=(matrix_scheme,), axis=1)
        save_data(matrix_scheme, scheme_df)


if __name__ == '__main__':
    main()