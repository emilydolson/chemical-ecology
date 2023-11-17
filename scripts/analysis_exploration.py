from common import get_plots_path, get_matrix_function, get_common_columns
from analysis_scheme import read_scheme_data
import matplotlib.pyplot as plt
from analysis import read_data
import seaborn as sns
import inspect
import csv
import sys


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
              f"-WORLD_WIDTH {10} -WORLD_HEIGHT {10}")


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
        #get_highest_scoring_matrices(df, 1)
    else:
        scheme = 'combined'
        df, param_names = read_data()
    individual_scatter(df.loc[df['score'] > -5], 'ntypes', 'score', 'adaptive', scheme)


if __name__ == '__main__':
    main()