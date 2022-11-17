import pandas as pd
import matplotlib.pyplot as plt
import sys


def main():
    file_location = sys.argv[1]

    df = pd.read_csv(f'{file_location}/all_data.csv')
    std_err = df.groupby('Time')['mean_Fitness'].sem()
    df = df.groupby('Time').mean()

    plt.errorbar(df.index, df['mean_Fitness']/max(df['mean_Fitness']), yerr=std_err/max(df['mean_Fitness']), color='blue', ecolor='lightblue')
    plt.xlabel('Time')
    plt.ylabel('Mean Growth Rate (normalized)')
    plt.savefig(f'{file_location}/graph.png')


if __name__ == '__main__':
    main()