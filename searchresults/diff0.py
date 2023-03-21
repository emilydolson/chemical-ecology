import numpy as np 
import csv 
import sys
import matplotlib 
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

#Plot for diffusion values of zero
def plot(dir, class_num):
    df = pd.read_csv(dir + '/Biomass.csv')
    df = df.round({"Diffusion": 1, "Seeding": 2, "Clear": 1})
    max_score = df['Score'].max()
    min_score = df['Score'].min() 

    df = df[df['Diffusion'] == 0]

    pivot_table = pd.pivot_table(df, values='Score', index='Seeding', columns='Clear')

    # Create the heatmap
    sns.heatmap(pivot_table, cmap="YlGnBu", vmin=-.2, vmax=.4)

    plt.suptitle(f"Diffusion values of 0 for class {class_num} matrix", fontsize=20)

    # norm = plt.Normalize(max_score, min_score)
    # sm = plt.cm.ScalarMappable(cmap="YlGnBu", norm=norm)

    plt.savefig(dir + "/diff0.png", bbox_inches="tight", pad_inches=0.1)


if __name__ == '__main__':
    plot(sys.argv[1], sys.argv[2])