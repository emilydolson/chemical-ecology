import numpy as np 
import csv 
import sys
import matplotlib 
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


def plot(dir):
    measures = ["Biomass", "Growth_Rate", "Heredity", "Invasion_Ability", "Resiliance"]
    for measure in measures:
        data = pd.read_csv(dir + "/" + measure + ".csv")
        data = data.round({"Diffusion": 1, "Seeding": 1, "Clear": 1})

        # Get the unique values of Clear
        clear_values = np.unique(data["Clear"])

        fig, axs = plt.subplots(nrows=2, ncols=6, figsize=(30, 10))
        fig.subplots_adjust(hspace=.02, wspace=.2)

        for i, clear in enumerate(clear_values):
            subset = data[data["Clear"] == clear]
            
            # Pivot the data to create a matrix for the heatmap
            heatmap_data = subset.pivot(index="Seeding", columns="Diffusion", values="Score")

            sns.heatmap(heatmap_data, cmap="YlGnBu", ax=axs.flat[i], vmin=-.3, vmax=.5, cbar=False)
            #Make heatmaps square
            axs.flat[i].set_aspect('equal')

            axs.flat[i].set_title(f"Clear = {clear}")

        fig.suptitle(f"Heatmaps for different Clear values of {measure}", fontsize=20)

        #delete the empty figure that comes with the subplots
        for ax in axs.flat:
            if not bool(ax.has_data()):
                fig.delaxes(ax)

        # Create a colorbar
        norm = plt.Normalize(-.3, .5)
        sm = plt.cm.ScalarMappable(cmap="YlGnBu", norm=norm)
        fig.colorbar(sm, ax=axs, location="right", pad=.01)

        plt.savefig(dir + "/" + f"{measure}_hmap.png", bbox_inches="tight", pad_inches=0.1)


if __name__ == '__main__':
    plot(sys.argv[1])