from common import get_plots_path, get_matrix_function, get_common_columns
from analysis_scheme import read_scheme_data
from matplotlib.colors import BoundaryNorm
import matplotlib.pyplot as plt
from analysis import read_data
import networkx as nx
import seaborn as sns
import numpy as np
import inspect
import csv
import sys


def visualize_graph(file_name):
    matrix = np.loadtxt(file_name, delimiter=',')
    G = nx.DiGraph(matrix)

    cmap = plt.get_cmap('PuOr')
    cmaplist = [cmap(i) for i in range(cmap.N)]
    cmap = cmap.from_list('Custom cmap', cmaplist, cmap.N)
    bounds = [-1, -0.75, -0.5, -0.25, -.0001, .0001, 0.25, 0.5, 0.75, 1]
    norm = BoundaryNorm(bounds, cmap.N)

    weights = nx.get_edge_attributes(G, 'weight')
    color_map = [cmap(norm(weights[edge])) for edge in G.edges()]

    plt.figure(figsize=(8,8))
    pos = nx.circular_layout(G)
    nx.draw(G, pos=pos, with_labels=True, node_size=1000, edge_color=color_map, node_color='lightgrey')

    plt.savefig(f'{file_name[:-4]}_graph.png')
    plt.close()


def visualize_matrix(file_name):
    matrix = np.loadtxt(file_name, delimiter=',')

    cmap = plt.get_cmap('PuOr')
    cmaplist = [cmap(i) for i in range(cmap.N)]
    cmap = cmap.from_list('Custom cmap', cmaplist, cmap.N)
    bounds = [-1, -0.75, -0.5, -0.25, -.0001, .0001, 0.25, 0.5, 0.75, 1]
    norm = BoundaryNorm(bounds, cmap.N)

    print(cmap(0.5))
    print(norm(cmap(0.2)))

    plt.imshow(matrix, cmap=cmap, interpolation='none', norm=norm)
    plt.colorbar(ticks=[-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1])
    plt.title(file_name, fontsize=15)
    plt.savefig(f'{file_name[:-4]}.png')
    plt.close()


def write_matrix(interactions, output_name):
    with open(output_name, 'w') as f:
        wr = csv.writer(f)
        wr.writerows(interactions)


def get_highest_scoring_matrices(df, n):
    df_best = df.nlargest(n, 'score')
    for index, row in df_best.iterrows():
        create_matrix_func = get_matrix_function(row['scheme'])
        scheme_args = inspect.getfullargspec(create_matrix_func)[0]
        scheme_args[-1] = 'replicate'
        matrix_params = [row[x] for x in scheme_args]
        matrix = create_matrix_func(*matrix_params)
        write_matrix(matrix, f'matrix_{index}.dat')
        print(f'Matrix {index}')
        for col in get_common_columns():
            try:
                print(f"\t{col} {row[col]}")
            except:
                continue
        print(f"./chemical-ecology " +
              f"-DIFFUSION {row['diffusion']} -SEEDING_PROB {row['seeding']} -PROB_CLEAR {row['clear']} " +
              f"-INTERACTION_SOURCE matrix_{index}.dat -SEED {row['replicate']} -N_TYPES {row['ntypes']} " +
              f"-WORLD_WIDTH {row['world_size']} -WORLD_HEIGHT {row['world_size']}")


def count_ntypes(df):
    for ntypes in sorted(df['ntypes'].unique()):
        print(f"{ntypes} - {len(df.loc[df['ntypes'] == ntypes])}")


def individual_scatter(df, x, y, hue, scheme):
    plt.figure()
    sns.scatterplot(x=x, y=y, data=df, hue=hue)
    plt.xlabel(x)
    plt.ylabel(y)
    plt.savefig(f'{get_plots_path()}{scheme}/scatter_{x}_{y}.png')
    plt.close()


def main():
    if len(sys.argv) == 2:
        scheme = sys.argv[1]
        df, param_names = read_scheme_data(scheme)
        get_highest_scoring_matrices(df.loc[df['ntypes'] == 10], 1)
    else:
        scheme = 'combined'
        df, param_names = read_data()
    #individual_scatter(df.loc[df['ntypes'] == 10], 'score', 'num_communities', 'adaptive', scheme)
    #print(len(df.loc[(df['ntypes'] == 10) & (df['num_communities'] >= 80)]))


if __name__ == '__main__':
    main()