from common import get_processed_data_path, get_plots_path, get_common_columns, get_common_param_columns
from itertools import combinations
import matplotlib.pyplot as plt
from scipy import stats
import seaborn as sns
import pandas as pd
import numpy as np
import os


'''
Plots for either individual schemes or all combined
'''
def histograms(df, param_names, scheme):
    figure, axis = plt.subplots(5, 4, figsize=(15,15))
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
    figure.suptitle(f'{scheme} parameter histograms')
    plt.savefig(f'{get_plots_path()}{scheme}/histogram_param.png')
    plt.close()


def boxplots(df, param_names, scheme):
    figure, axis = plt.subplots(5, 4, figsize=(15,15))
    row = 0
    col = 0
    for param in param_names:
        sns.violinplot(data=df, x='adaptive', y=param, cut=0, ax=axis[row][col])
        axis[row][col].set_title(param)
        row += 1
        if row % 5 == 0:
            col += 1
            row = 0
    figure.tight_layout(rect=[0, 0.03, 1, 0.95])
    figure.suptitle(f'{scheme} parameter boxplots')
    plt.savefig(f'{get_plots_path()}{scheme}/boxplot_param.png')
    plt.close()


def fitness_correlation(df, param_names, scheme):
    figure, axis = plt.subplots(1, 1, figsize=(8,4))
    for i in df['replicate'].unique():
        df_replicate = df.loc[df['replicate'] == i]
        scores = df_replicate['score'].values
        correlations = [np.corrcoef(df_replicate[param].values, scores)[0,1] for param in param_names]
        axis.scatter(correlations, range(len(param_names)))
    axis.yaxis.set_ticks(list(range(len(param_names))), param_names)
    figure.tight_layout(rect=[0, 0.03, 1, 0.95])
    figure.suptitle(f'{scheme} correlation between parameters and score')
    figure.supxlabel('Correlation with score')
    plt.grid(axis='y')
    plt.savefig(f'{get_plots_path()}{scheme}/correlation_score.png')
    plt.close()


def param_correlation_heatmap(df, param_names, scheme):
    figure = plt.figure()
    df_params = df[param_names]
    corr = df_params.corr()
    sns.heatmap(corr, center=0, vmin=-1, vmax=1, cmap='coolwarm')
    figure.tight_layout()
    plt.savefig(f'{get_plots_path()}{scheme}/correlation_param.png')
    plt.close()


'''
Plots for all schemes combined
'''
def score_boxplot(df):
    sns.violinplot(data=df, x='scheme', y='score', cut=0)
    plt.title(f'score box plots')
    plt.savefig(f'{get_plots_path()}boxplot_score.png')
    plt.close()


'''
Other analysis
'''
def significance(df):
    scheme_combos = list(combinations(df['scheme'].unique(), 2))
    for rep in df['replicate'].unique():
        print(f'Replicate {rep}')
        df_rep = df.loc[df['replicate'] == rep]

        print('\tKruskal-Wallace')
        kruskal_stat, kruskal_p = stats.kruskal(*[df_rep.loc[df_rep['scheme'] == x]['score'] for x in df_rep['scheme'].unique()])
        print(f'\t\tstatistic: {kruskal_stat}, p-value: {kruskal_p}')
        print()

        print('\tWilcoxon Rank Sums')
        for scheme1, scheme2 in scheme_combos:
            print('\t', scheme1, scheme2)
            df_scheme1 = df_rep.loc[df_rep['scheme'] == scheme1]
            df_scheme2 = df_rep.loc[df_rep['scheme'] == scheme2]
            rank_stat, rank_p = stats.ranksums(df_scheme1['score'], df_scheme2['score'])
            print(f'\t\tstatistic: {rank_stat}, p-value: {rank_p}')
    print(f'Bonferroni Correction: {0.05/len(scheme_combos)}')


'''
Main functions
'''
def read_data():
    file_path = get_processed_data_path()
    dfs = []
    for scheme_df_file in os.listdir(file_path):
        scheme_df = pd.read_pickle(file_path+scheme_df_file)
        dfs.append(scheme_df)
    if len(dfs) > 1:
        overlapping_cols = list(set.intersection(*[set(x.columns) for x in dfs]))
        dfs = [x[overlapping_cols] for x in dfs]
    df = pd.concat(dfs)
    param_names = list(set(df.columns) - set(get_common_columns())) + get_common_param_columns()
    df = df.reset_index()
    df['adaptive'] = np.where(df['score'] > 0.5, True, False)
    return df, param_names


def main():
    df, param_names = read_data()

    fitness_correlation(df, param_names, 'combined')
    histograms(df, param_names, 'combined')
    param_correlation_heatmap(df, param_names, 'combined')
    score_boxplot(df)
    significance(df)


if __name__ == '__main__':
    main()