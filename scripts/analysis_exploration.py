from analysis_scheme import read_scheme_data
from common import get_plots_path
import matplotlib.pyplot as plt
from analysis import read_data
import seaborn as sns
import sys


def individual_scatter(df, x, y, scheme):
    plt.figure()
    sns.scatterplot(x=x, y=y, data=df, hue='adaptive')
    plt.xlabel(x)
    plt.ylabel(y)
    plt.savefig(f'{get_plots_path()}{scheme}/scatter_{x}_{y}.png')
    plt.close()


def main():
    if len(sys.argv) == 2:
        scheme = sys.argv[1]
        df, param_names = read_scheme_data(scheme)
    else:
        scheme = 'combined'
        df, param_names = read_data()
    individual_scatter(df, 'connectance', 'clear', scheme)


if __name__ == '__main__':
    main()