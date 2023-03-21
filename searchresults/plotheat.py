import numpy as np 
import csv 
import sys
import matplotlib 
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter


def plot(dir, class_num):
    measures = ["Biomass", "Growth_Rate", "Heredity", "Invasion_Ability", "Resiliance"]
    # Need this to make color bar consistent accross graphs 
    best_score = -10
    worst_score = 10
    for measure in measures:
        data = pd.read_csv(dir + "/" + measure + ".csv")
        max_score = data['Score'].max()
        min_score = data['Score'].min() 
        if max_score > best_score:
            best_score = max_score
        if min_score < worst_score:
            worst_score = min_score

    for measure in measures:
        data = pd.read_csv(dir + "/" + measure + ".csv")
        data = data.round({"Diffusion": 1, "Seeding": 2, "Clear": 1}) 

        # Get the unique values of Clear
        clear_values = np.unique(data["Clear"])

        fig, axs = plt.subplots(nrows=2, ncols=6, figsize=(30, 10))
        fig.subplots_adjust(hspace=.02, wspace=.2)

        for i, clear in enumerate(clear_values):
            subset = data[data["Clear"] == clear]
            
            # Pivot the data to create a matrix for the heatmap
            heatmap_data = subset.pivot(index="Seeding", columns="Diffusion", values="Score")
            heatmap_data.index = heatmap_data.index.astype(str).str[1:]

            sns.heatmap(heatmap_data, cmap="YlGnBu", ax=axs.flat[i], vmin=-.3, vmax=.5, cbar=False)
            #Make heatmaps square
            axs.flat[i].set_aspect('equal')

            axs.flat[i].set_title(f"Clear = {clear}")

        fig.suptitle(f"Adaptive scores for a class {class_num} matrix", fontsize=20)

        #delete the empty figure that comes with the subplots
        for ax in axs.flat:
            if not bool(ax.has_data()):
                fig.delaxes(ax)

        # Create a colorbar
        norm = plt.Normalize(worst_score, best_score)
        sm = plt.cm.ScalarMappable(cmap="YlGnBu", norm=norm)
        fig.colorbar(sm, ax=axs, location="right", pad=.01)

        plt.savefig(dir + "/" + f"{measure}_hmap.png", bbox_inches="tight", pad_inches=0.1)


if __name__ == '__main__':
    plot(sys.argv[1], sys.argv[2])