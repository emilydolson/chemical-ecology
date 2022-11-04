import pandas as pd
import matplotlib.pyplot as plt
import sys


def main():
    file_location = sys.argv[1]

    df = pd.read_csv(f'{file_location}/all_data.csv')
    df = df.groupby('Time').mean()
    plt.plot(df.index, df['mean_Fitness']/max(df['mean_Fitness']))
    plt.xlabel('Time')
    plt.ylabel('Mean Growth Rate (normalized)')
    plt.savefig(f'{file_location}/graph.png')


if __name__ == '__main__':
    main()