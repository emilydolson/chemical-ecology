import sys
import random
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sb
import networkx as nx
from scipy.stats import kstest, chisquare
from sklearn import tree
from scipy import stats

from get_results import get_data, get_matrix_func
from matrices import visualize_network


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


def visualize_decision_tree(df, score, param_names, matrix_scheme):
    X = df.loc[:, param_names]
    y = df.loc[:, score]
    clf = tree.DecisionTreeClassifier(max_depth=5, min_samples_split=1000).fit(X, y)

    plt.figure(figsize=(35,10))
    tree.plot_tree(clf, filled=True, fontsize=7, feature_names=param_names)
    plt.savefig(f'plots/{matrix_scheme}/decisiontree.png', dpi=75)
    plt.close()


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
            axis[row][col].set_xticks(np.arange(0, len(bins), 1), [round(x, 2) for x in bins])
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


def individual_scatter(df, score, param, color, matrix_scheme, ylim_min=0, ylim_max=5, cmap='Oranges'):
    figure, axis = plt.subplots(1, 1)

    sc = axis.scatter(df[param], df[score], s=2, c=df[color], cmap=cmap)
    cb = figure.colorbar(sc, ax=axis)
    cb.set_label(color)
    #axis.set_ylim(ylim_min, ylim_max)
    axis.set_yscale('log')
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


def mangal(df_generation, df_network):
    df = pd.merge(df_generation, df_network, on=['Biomass', 'Growth_Rate', 'Heredity', 'Biomass_Score',
       'Growth_Rate_Score', 'Heredity_Score', 'Invasion_Ability_Score',
       'Resilience_Score', 'Final_Communities_Present', 'diffusion', 'seeding', 'clear', 'Biomass_Adaptive',
       'Growth_Rate_Adaptive', 'Heredity_Adaptive', 'Invasion_Ability_Adaptive', 'Resilience_Adaptive'])
    for score in [x for x in df.columns if 'Score' in x]:
        for matrix_num in sorted(pd.unique(df['name'])):
            neg = df.loc[df['name'] == matrix_num]['num_neg_stronger'].values[0]
            print(f"Adaptive {score} for {matrix_num}: {len(df.loc[(df[score[:-5]+'Adaptive'] == True) & (df['name'] == matrix_num)])}, {round(neg, 3)}")
        print()
        for matrix_num in sorted(pd.unique(df['matrix_num'])):
            print(f"Adaptive {score} for {matrix_num}: {len(df.loc[(df[score[:-5]+'Adaptive'] == True) & (df['matrix_num'] == matrix_num)])}")
        print()

    df['replicate2'] = df['name'].str[-5].astype(int)
    for score in [x for x in df.columns if 'Score' in x]:
        kruskal_stat, kruskal_p = stats.kruskal(*[df.loc[df['replicate2'] == x][score] for x in pd.unique(df['replicate2'])])
        print(f'\t{score} statistic: {kruskal_stat}, p-value: {kruskal_p}')

    group = 'name'
    df_count = df[[group]]
    df_count.insert(1, 'total', df.groupby(group)[group].transform('count'))
    for score in [x for x in df.columns if 'Adaptive' in x]:
        df_count = df_count.join(df.groupby(group)[score].sum(), on=[group])
    df_count = df_count.drop_duplicates()
    df_count['replicate'] = df_count['name'].str[-5].astype(int)
    df_count['Matrix'] = df_count['name'].str.split('_').str[0]
    plt.figure()
    axis = sb.barplot(data=df_count, x='Matrix', y='Biomass_Adaptive', hue='replicate')
    plt.ylabel('Count of Adaptive Biomass Results')
    plt.title('Count of Adaptive Biomass Results by Mangal Matrix and Replicate')
    plt.savefig(f'plots/mangal/barplot.png')
    plt.close()


def random_stats(df):
    for score in [x for x in df.columns if 'Score' in x]:
        kruskal_stat, kruskal_p = stats.kruskal(*[df.loc[df['matrix_seed'] == x][score] for x in pd.unique(df['matrix_seed'])])
        print(f'\t{score} statistic: {kruskal_stat}, p-value: {kruskal_p}')


def sample_visualizations(df, param_names, matrix_scheme):
    create_matrix_func = get_matrix_func(matrix_scheme)

    params = df[param_names+[x for x in df.columns if 'Score' in x]].values.tolist()
    if matrix_scheme == 'mangal':
        for matrix_name in sorted(os.listdir('mangal_matrices')):
            matrix = np.loadtxt(f'mangal_matrices/{matrix_name}', delimiter=',').tolist()
            title = matrix_name[0:2] if '89' in matrix_name else matrix_name[0:4]
            visualize_network(matrix, f'samples/{matrix_name[:-4]}.png', f'Matrix {title}')
    else:
        for i,genome in enumerate(random.sample(params, 10)):
            if matrix_scheme == 'motif':
                matrix = create_matrix_func(*genome[0:len(param_names)][3:8])
                genome = [round(x,3) if type(x) == float else x for x in genome]
            else:
                genome = [int(x) if x.is_integer() else round(x,3) for x in genome]
                matrix = create_matrix_func(*genome[0:len(param_names)][3:])
            visualize_network(matrix, f'samples/{matrix_scheme}_{i}.png', f'{genome[0:len(param_names)]}\n{genome[len(param_names):]}')


def scheme_plots(df, param_names, file_name):
    param_histograms(df, 'All', param_names, file_name)
    correlation_heatmap(df, 'All', param_names, file_name)
    for score in [x for x in df.columns if 'Score' in x]:
        df_adaptive = df.loc[df[score[:-5]+'Adaptive'] == True]
        fitness_correlation(df, score, param_names, file_name)
        param_histograms(df_adaptive, score, param_names, file_name)
        correlation_heatmap(df_adaptive, score, param_names, file_name)
        cdf(df, score, param_names, file_name)


def main(file_name):
    df_generation, df_network, abiotic_params, param_names, analysis_param_names, network_param_names = get_data(file_name)

    scheme_plots(df_generation, abiotic_params+analysis_param_names, file_name)
    scheme_plots(df_network, abiotic_params+network_param_names, file_name+'/network')

    sample_visualizations(df_generation.loc[df_generation['Biomass_Adaptive'] == True], abiotic_params+param_names, file_name)
    visualize_decision_tree(df_generation, abiotic_params+analysis_param_names, file_name, regression=False)
    individual_scatter(df.loc[df['Biomass_Adaptive'] == True], 'Biomass_Score', 'clear', 'seeding', file_name, ylim_min=1, ylim_max=None, cmap='plasma_r')


if __name__ == '__main__':
    main(sys.argv[1])
