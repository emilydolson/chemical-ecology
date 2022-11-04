import sys
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
import numpy as np


def main():
    file_name = sys.argv[1]
    experiment_name = sys.argv[2]
    matrix = np.loadtxt(file_name, delimiter=',')

    #https://stackoverflow.com/questions/23994020/colorplot-that-distinguishes-between-positive-and-negative-values
    cmap = plt.get_cmap('PuOr')
    cmaplist = [cmap(i) for i in range(cmap.N)]
    cmap = cmap.from_list('Custom cmap', cmaplist, cmap.N)
    bounds = [-1, -0.75, -0.5, -0.25, -.0001, .0001, 0.25, 0.5, 0.75, 1]
    norm = BoundaryNorm(bounds, cmap.N)

    plt.imshow(matrix, cmap=cmap, interpolation='none', norm=norm)
    plt.colorbar()
    plt.savefig(f'../experiments/{experiment_name}/interaction_matrix.png')


if __name__ == '__main__':
    main()