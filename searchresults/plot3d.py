#%matplotlib inline

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pylab import *
  
def main():
    dataset = pd.read_csv("Biomass.csv")
    x = dataset["Diffusion"].tolist()
    y = dataset["Seeding"].tolist()
    z = dataset["Clear"].tolist()
    
    colo = dataset["Score"].tolist()
    
    # creating 3d figures
    fig = plt.figure(figsize=(8, 5))
    ax = fig.add_subplot(111, projection='3d')
    
    # configuring colorbar
    color_map = cm.ScalarMappable(cmap=cm.gray)
    color_map.set_array(colo)
    
    # creating the heatmap
    img = ax.scatter(x, y, z, marker='s', 
                    s=99, color='gray')
    plt.colorbar(color_map)
    
    # adding title and labels
    ax.set_title("3D Heatmap")
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('')
    
    # displaying plot
    plt.savefig("test.png")

main()