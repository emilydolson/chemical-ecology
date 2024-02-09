import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from common import get_non_property_column_names, get_plots_path, get_processed_data_path

'''
Plots for either individual schemes or all combined
'''
def histograms(df, param_names, exp_name):
    figure, axis = plt.subplots(5, 5, figsize=(15,15))
    row = 0
    col = 0
    for param in param_names:
        sns.histplot(data=df, x=param, hue='adaptive', multiple='stack', stat='proportion', ax=axis[row][col])
        axis[row][col].set_title(param)
        row += 1
        if row % 5 == 0:
            col += 1
            row = 0
    figure.tight_layout(rect=[0, 0.03, 1, 0.95])
    figure.suptitle(f'{exp_name} parameter histograms')
    plt.savefig(f'{get_plots_path()}{exp_name}/histogram_param.png')
    plt.close()


def boxplots(df, param_names, exp_name):
    figure, axis = plt.subplots(5, 4, figsize=(15,15))
    row = 0
    col = 0
    for param in param_names:
        sns.boxplot(data=df, x='adaptive', y=param, ax=axis[row][col])
        axis[row][col].set_title(param)
        row += 1
        if row % 5 == 0:
            col += 1
            row = 0
    figure.tight_layout(rect=[0, 0.03, 1, 0.95])
    figure.suptitle(f'{exp_name} parameter boxplots')
    plt.savefig(f'{get_plots_path()}{exp_name}/boxplot_param.png')
    plt.close()


def fitness_correlation(df, param_names, exp_name, grouping):
    figure, axis = plt.subplots(1, 1, figsize=(8,4))
    for i in df[grouping].unique():
        df_group = df.loc[df[grouping] == i]
        scores = df_group["score"].values
        correlations = [np.corrcoef(df_group[param].values, scores)[0,1] for param in param_names]
        axis.scatter(correlations, range(len(param_names)), label=i)
    axis.yaxis.set_ticks(list(range(len(param_names))), param_names)
    figure.tight_layout(rect=[0, 0.03, 1, 0.95])
    figure.suptitle(f"{exp_name} correlation between parameters and score")
    figure.supxlabel("Correlation with score")
    figure.legend()
    plt.grid(axis="y")
    plt.savefig(f'{get_plots_path()}{exp_name}/correlation_score.png')
    plt.close()


'''
Main functions
'''
def read_data(exp_name):
    file_path = get_processed_data_path()
    if exp_name == "combined":
        dfs = []
        for exp_df_file in os.listdir(file_path):
            exp_df = pd.read_pickle(file_path+exp_df_file)
            dfs.append(exp_df)
        df = pd.concat(dfs)
    else:
        df = pd.read_pickle(f"{file_path}{exp_name}.pkl")
    param_names = list(set(df.columns) - set(get_non_property_column_names()))
    df = df.reset_index()
    df["adaptive"] = np.where(df["score"] > 2, True, False)
    return df, param_names


def main(exp_name):
    df, param_names = read_data(exp_name)

    path = f"{get_plots_path()}{exp_name}"
    if not os.path.exists(path):
        os.makedirs(path)

    if exp_name == "combined":
        fitness_correlation(df, param_names, exp_name, "experiment")
        histograms(df, param_names, exp_name)
        boxplots(df, param_names, exp_name)
    else:
        df_exp = df.loc[df["experiment"] == exp_name]
        fitness_correlation(df_exp, param_names, exp_name, "config_num")
        histograms(df_exp, param_names, exp_name)
        boxplots(df_exp, param_names, exp_name)

        for config_num in df_exp["config_num"].unique():
            if not os.path.exists(f"{path}/{config_num}"):
                os.makedirs(f"{path}/{config_num}")
            cond_name = f"{exp_name}/{config_num}"
            df_config = df_exp.loc[df_exp["config_num"] == config_num]
            fitness_correlation(df_config, param_names, cond_name, "replicate")
            histograms(df_config, param_names, cond_name)
            boxplots(df_config, param_names, cond_name)


if __name__ == "__main__":
    if len(sys.argv) == 2:
        exp_name = sys.argv[1]
        main(exp_name)
    else:
        print("Please provide an experiment name or \"combined\".")