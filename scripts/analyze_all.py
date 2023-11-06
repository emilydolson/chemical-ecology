#analysis across all schemes (looking for scheme-independent patterns)

from common import get_processed_data_path, get_plots_path
from itertools import combinations
from scipy import stats
import pandas as pd
import os


def significance(df):
    for rep in df['replicate'].unique():
        print(f'Replicate {rep}')
        df_rep = df.loc[df['replicate'] == rep]
        kruskal_stat, kruskal_p = stats.kruskal(*[df_rep.loc[df_rep['scheme'] == x]['score'] for x in df_rep['scheme'].unique()])
        print(f'\tstatistic: {kruskal_stat}, p-value: {kruskal_p}')
        print()
    scheme_combos = combinations(df['scheme'].unique(), 2)
    for scheme1, scheme2 in scheme_combos:
        print(scheme1, scheme2)
        df_scheme1 = df.loc[df['scheme'] == scheme1]
        df_scheme2 = df.loc[df['scheme'] == scheme2]
        rank_stat, rank_p = stats.ranksums(df_scheme1['score'], df_scheme2['score'])
        print(f'\tstatistic: {rank_stat}, p-value: {rank_p}')
    print(f'Bonferroni Correction: {0.05/len(list(scheme_combos))}')


def read_data(file_path):
    dfs = []
    for scheme_df_file in os.listdir(file_path):
        scheme_df = pd.read_pickle(file_path+scheme_df_file)
        dfs.append(scheme_df)
    if len(dfs) > 1: #not tested
        common_columns = list(set(dfs[0].columns).intersection(set(dfs[1].columns)))
        dfs = [x[[common_columns]] for x in dfs]
    df = pd.concat(dfs)
    return df


def main():
    file_path = get_processed_data_path()
    df = read_data(file_path)


if __name__ == '__main__':
    main()