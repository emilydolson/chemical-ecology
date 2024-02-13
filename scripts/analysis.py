import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from common import get_non_property_column_names, get_plots_path, get_processed_data_path

sns.set_palette(sns.color_palette(["#ef7c8e", "#4c956c", "#d68c45", "#809bce", "#66786a", "#f2d91a", "#ab46d2"]))

"""
Data Exploration Plots
"""
def histograms(df, param_names, exp_name):
    figure, axis = plt.subplots(5, 5, figsize=(15,15))
    row = 0
    col = 0
    for param in param_names:
        sns.histplot(data=df, x=param, hue="adaptive", multiple="stack", stat="proportion", ax=axis[row][col])
        axis[row][col].set_title(param)
        row += 1
        if row % 5 == 0:
            col += 1
            row = 0
    figure.tight_layout(rect=[0, 0.03, 1, 0.95])
    figure.suptitle(f"{exp_name} parameter histograms")
    plt.savefig(f"{get_plots_path()}{exp_name}/histogram_param.png")
    plt.close()


def parameter_boxplots(df, param_names, exp_name):
    figure, axis = plt.subplots(5, 4, figsize=(15,15))
    row = 0
    col = 0
    for param in param_names:
        sns.boxplot(data=df, x="adaptive", y=param, ax=axis[row][col])
        axis[row][col].set_title(param)
        row += 1
        if row % 5 == 0:
            col += 1
            row = 0
    figure.tight_layout(rect=[0, 0.03, 1, 0.95])
    figure.suptitle(f"{exp_name} parameter boxplots")
    plt.savefig(f"{get_plots_path()}{exp_name}/boxplot_param.png")
    plt.close()


def adaptive_parameter_boxplots(df, param_names, exp_name):
    figure, axis = plt.subplots(5, 4, figsize=(15,15))
    row = 0
    col = 0
    for param in param_names:
        sns.boxplot(data=df, x="config_num", y=param, hue="adaptive", ax=axis[row][col])
        axis[row][col].set_title(param)
        row += 1
        if row % 5 == 0:
            col += 1
            row = 0
    figure.tight_layout(rect=[0, 0.03, 1, 0.95])
    figure.suptitle(f"{exp_name} adaptive parameter boxplots")
    plt.savefig(f"{get_plots_path()}{exp_name}/adaptive_boxplot_param.png")
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
    plt.savefig(f"{get_plots_path()}{exp_name}/correlation_score.png")
    plt.close()


"""
Other Functions
"""
def param_correlation_heatmap(df, param_names, exp_name):
    figure = plt.figure()
    df_params = df[param_names]
    corr = df_params.corr()
    sns.heatmap(corr, center=0, vmin=-1, vmax=1, cmap="coolwarm")
    figure.tight_layout()
    plt.savefig(f"{get_plots_path()}{exp_name}/correlation_param.png")
    plt.close()


def plot_score_correlated_properties(df, param_names, exp_name):
    configs = sorted(df["config_num"].unique())
    num_params = len(param_names)
    highly_correlated = set()
    correlations = dict()
    for i in configs:
        df_group = df.loc[df["config_num"] == i]
        scores = df_group["score"].values
        correlations[i] = {param:np.corrcoef(df_group[param].values, scores)[0,1] for param in param_names}
        for j in range(num_params):
            param = param_names[j]
            correlation = correlations[i][param]
            if abs(correlation) > 0.2:
                highly_correlated.add(param)

    for p in highly_correlated:
        max_p = max(df[p])
        min_p = min(df[p])
        max_score = max(df["score"])
        min_score = min(df["score"])
        num_configs = len(configs)
        fig_col_cnt = 2 if num_configs <= 4 else 4
        fig_row_cnt = int(np.ceil(num_configs/fig_col_cnt))
        figure, axis = plt.subplots(fig_row_cnt, fig_col_cnt, figsize=(5*fig_row_cnt, 3*fig_col_cnt), squeeze=False)
        fig_row = 0
        fig_col = 0
        for i in configs:
            df_group = df.loc[df["config_num"] == i]
            axis[fig_row][fig_col].scatter(df_group[p], df_group["score"])
            axis[fig_row][fig_col].set_title(f"Config {i}, r = {round(correlations[i][p], 2)}")
            axis[fig_row][fig_col].set_xlim(min_p, max_p)
            axis[fig_row][fig_col].set_ylim(min_score, max_score)
            fig_col += 1
            if fig_col % fig_row_cnt == 0:
                fig_row += 1
                fig_col = 0
        figure.tight_layout(rect=[0, 0.03, 1, 0.95])
        figure.suptitle(f"{exp_name} scatter plots")
        figure.supylabel("Score")
        figure.supxlabel(p)
        plt.savefig(f"{get_plots_path()}{exp_name}/score_{p}_scatter.png")
        plt.close()


"""
Main functions
"""
def read_data(exp_name):
    file_path = get_processed_data_path()
    df = pd.read_pickle(f"{file_path}{exp_name}.pkl")
    param_names = list(set(df.columns) - set(get_non_property_column_names()))
    df = df.reset_index()
    df["adaptive"] = np.where(df["score"] > 2, True, False)
    return df, param_names


def main(exp_name, config_level_analysis=False):
    df, param_names = read_data(exp_name)

    path = f"{get_plots_path()}{exp_name}"
    if not os.path.exists(path):
        os.makedirs(path)

    fitness_correlation(df, param_names, exp_name, "config_num")
    histograms(df, param_names, exp_name)
    adaptive_parameter_boxplots(df, param_names, exp_name)
    plot_score_correlated_properties(df, param_names, exp_name)

    if exp_name == "uniform":
        param_correlation_heatmap(df, param_names, exp_name)
    
    if config_level_analysis:
        config_nums = df["config_num"].unique()
        for config_num in config_nums:
            if not os.path.exists(f"{path}/{config_num}"):
                os.makedirs(f"{path}/{config_num}")
            cond_name = f"{exp_name}/{config_num}"
            df_config = df.loc[df["config_num"] == config_num]
            fitness_correlation(df_config, param_names, cond_name, "replicate")
            histograms(df_config, param_names, cond_name)
            parameter_boxplots(df_config, param_names, cond_name)


if __name__ == "__main__":
    if len(sys.argv) == 2:
        exp_name = sys.argv[1]
        main(exp_name)
    else:
        print("Please provide an experiment name.")