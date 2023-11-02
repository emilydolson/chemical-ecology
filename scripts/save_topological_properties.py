from common import get_matrix_function, get_common_columns, get_file_path
import networkx as nx
import pandas as pd
import numpy as np
import os


def save_data(df):
    return

    
def get_topological_properties(row, matrix_scheme):
    create_matrix_func = get_matrix_function(matrix_scheme)
    scheme_specific_params = row[~row.index.isin(get_common_columns())].tolist()
    matrix_params = [row['ntypes']] + scheme_specific_params + [row['replicate']]
    matrix = create_matrix_func(*matrix_params)

    num_nodes = len(matrix)
    num_interactions = sum([sum([1 for y in x if y != 0]) for x in matrix])
    
    graph = nx.DiGraph.reverse(nx.DiGraph(np.array(matrix)))
    
    row['cluster_coeff'] = sum(dict(nx.clustering(graph)).values()) / num_nodes

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
    file_path = get_file_path()
    for matrix_scheme in os.listdir(file_path):
        scheme_df = read_data(file_path, matrix_scheme)
        scheme_df = scheme_df.apply(get_topological_properties, args=(matrix_scheme,), axis=1)
        print(scheme_df)


if __name__ == '__main__':
    main()