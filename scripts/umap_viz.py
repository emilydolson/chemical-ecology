import matplotlib.pyplot as plt
from umap import UMAP
import pandas as pd
import os


code_location = '/mnt/home/leithers/grant/chemical-ecology/'


def get_plots_path():
    return code_location+'scripts/output/plots/'


def get_processed_data_path():
    return code_location+'scripts/output/data/'


def get_common_param_columns():
    return ['diffusion', 'seeding', 'clear']


def get_common_columns():
    return ['scheme', 'replicate', 'score', 'ntypes'] + get_common_param_columns()


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
    param_names = list(set(df.columns) - set(get_common_columns()))
    df = df.reset_index()
    df['adaptive'] = False
    df.loc[df['score'] > 0.5, 'adaptive'] = True
    return df, param_names


def umap_plot(df, param_names):
    colors = ['pink', 'palegreen', 'wheat']
    darkcolors = ['crimson', 'green', 'darkorange']
    params = df[param_names].values
    mapped = UMAP(n_neighbors=15, min_dist=0.1, n_components=2, random_state=42).fit_transform(params)
    figure, axis = plt.subplots(1, 1, figsize=(8,4))
    for i,scheme in enumerate(df['scheme'].unique()):
        for adaptive in [False, True]:
            (X,Y) = list(zip(*[d for d,s,a in zip(mapped, df['scheme'].values, df['adaptive'].values) if s == scheme and a == adaptive]))
            if adaptive:
                axis.scatter(X, Y, c=darkcolors[i])
            else:
                axis.scatter(X, Y, label=scheme, c=colors[i])
    axis.legend()
    figure.suptitle('UMAP Visualization')
    plt.savefig(f'{get_plots_path()}umap.png')
    plt.close()


def main():
    df, param_names = read_data()
    umap_plot(df, param_names)


if __name__ == '__main__':
    main()