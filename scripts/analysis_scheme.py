from common import get_processed_data_path, get_plots_path, get_common_columns, get_common_param_columns, get_schemes
from analysis import histograms, fitness_correlation, param_correlation_heatmap, summary_statistics
import pandas as pd
import numpy as np
import sys
import os


def create_folder_structure(scheme):
    path = f'{get_plots_path()}{scheme}'
    if not os.path.exists(path):
        os.makedirs(path)


def read_scheme_data(scheme):
    file_path = f'{get_processed_data_path()}{scheme}.pkl'
    if os.path.exists(file_path):
        df = pd.read_pickle(file_path)
    else:
        print(f'{scheme} pickle not found')
        exit()
    param_names = list(set(df.columns) - set(get_common_columns())) + get_common_param_columns()
    df = df.reset_index()
    df['adaptive'] = np.where(df['score'] > 0.5, True, False)
    return df, param_names


def scheme_analysis(scheme):
    df, param_names = read_scheme_data(scheme)
    create_folder_structure(scheme)
    
    fitness_correlation(df, param_names, scheme)
    histograms(df, param_names, scheme)
    param_correlation_heatmap(df, param_names, scheme)
    summary_statistics(df, scheme)


def main():
    if len(sys.argv) == 2:
        scheme = sys.argv[1]
        scheme_analysis(scheme)
    else:
        for scheme in get_schemes():
            print(scheme)
            scheme_analysis(scheme)


if __name__ == '__main__':
    main()