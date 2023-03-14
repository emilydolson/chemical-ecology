import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.cbook import boxplot_stats
from scipy import stats
from itertools import combinations

from get_results import get_data


def significance(df):
    for scheme in pd.unique(df['scheme']): #ignore mangal and random results for kruskal
        print(scheme)
        df_scheme = df.loc[df['scheme'] == scheme]
        for score in [x for x in df.columns if 'Score' in x]:
            kruskal_stat, kruskal_p = stats.kruskal(*[df_scheme.loc[df_scheme['replicate'] == x][score] for x in pd.unique(df_scheme['replicate'])])
            print(f'\t{score} statistic: {kruskal_stat}, p-value: {kruskal_p}')
        print()
    for scheme1, scheme2 in combinations(pd.unique(df['scheme']), 2):
        print(scheme1, scheme2)
        df_scheme1 = df.loc[df['scheme'] == scheme1]
        df_scheme2 = df.loc[df['scheme'] == scheme2]
        for score in [x for x in df.columns if 'Score' in x]:
            rank_stat, rank_p = stats.ranksums(df_scheme1[score], df_scheme2[score])
            print(f'\t{score} statistic: {rank_stat}, p-value: {rank_p}')
    print(f"Bonferroni Correction: {0.05/len(list(combinations(pd.unique(df['scheme']), 2)))}")


def nonboxplots(df, group):
    for score in [x for x in df.columns if 'Score' in x]:
        df.loc[df[score] > 1e10, score] = 1e10
        df.loc[df[score] < 1e-10, score] = 1e-10
        df_a = df.loc[df[score[:-5]+'Adaptive'] == False]

        plt.figure()
        axis = sns.boxplot(data=df_a, x=group, y=score)
        axis.set(yscale='log', ylim=(1e-10, 1e10))
        plt.title(f'Non-Adaptive {score} Box Plots')
        plt.savefig(f'plots/combined/{score}_nonboxplot.png')
        plt.close()


def boxplots(df, group):
    for score in [x for x in df.columns if 'Score' in x]:
        df.loc[df[score] > 1e3, score] = 1e3
        df_a = df.loc[df[score[:-5]+'Adaptive'] == True]

        plt.figure()
        axis = sns.boxplot(data=df_a, x=group, y=score)
        axis.set(yscale='log', ylim=(1, 1e3))
        plt.title(f'Adaptive {score} Box Plots')
        plt.savefig(f'plots/combined0/{score}_boxplot.png')
        plt.close()


def scatter(df, param): #TODO this should compare network+abiotic params against each other (x=param, y=score, colored by scheme)
    figure, axis = plt.subplots(1, 5)

    sc = axis.scatter(df[param], df[score], c=df[color], cmap=cmap)
    cb = figure.colorbar(sc, ax=axis)
    cb.set_label(color)
    axis.set_ylim(ylim_min, ylim_max)
    figure.suptitle(f'{matrix_scheme} {param} and {score}')
    figure.supxlabel(param)
    figure.supylabel(score)

    plt.savefig(f'plots/{matrix_scheme}/{score}_{param}_scatter.png')
    plt.close()


def count_adaptive(df, group):
    df_count = df[[group]]
    df_count.insert(1, 'total', df.groupby(group)[group].transform('count'))
    for score in [x for x in df.columns if 'Adaptive' in x]:
        df_count = df_count.join(df.groupby(group)[score].sum(), on=[group])
    df_count = df_count.drop_duplicates()
    df_count.to_csv(f'num_adaptive_{group}.csv')


def main():
    schemes = ['klemm', 'klemm_random', 'motif', 'mangal', 'random']
    df = pd.DataFrame()
    for scheme in schemes:
        #df_generation, df_network, abiotic_params, param_names, analysis_param_names, network_param_names = get_data(scheme)
        df_network, abiotic_params, param_names, analysis_param_names = get_data(scheme, False)
        df_network.drop(columns=param_names, inplace=True)
        df_network['scheme'] = scheme
        df_network['scheme_rep'] = df_network['scheme'] + df_network['replicate'].map(str)
        df = pd.concat([df, df_network], ignore_index=True)

    significance(df)
    count_adaptive(df, 'scheme')
    boxplots(df, 'scheme')


if __name__ == '__main__':
    main()