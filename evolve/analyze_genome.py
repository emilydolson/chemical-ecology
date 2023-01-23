import ast
import csv
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from lfr_graph import create_matrix

def visualize_network(input_name, output_name):
    matrix = np.loadtxt(input_name, delimiter=',')

    #https://stackoverflow.com/questions/23994020/colorplot-that-distinguishes-between-positive-and-negative-values
    cmap = plt.get_cmap('PuOr')
    cmaplist = [cmap(i) for i in range(cmap.N)]
    cmap = cmap.from_list('Custom cmap', cmaplist, cmap.N)
    bounds = [-1, -0.75, -0.5, -0.25, -.0001, .0001, 0.25, 0.5, 0.75, 1]
    norm = BoundaryNorm(bounds, cmap.N)

    plt.imshow(matrix, cmap=cmap, interpolation='none', norm=norm)
    plt.colorbar()
    plt.savefig(output_name)


def create_matrix_unformatted(genome):
    diffusion = genome[0]
    seeding = genome[1]
    clear = genome[2]
    average_k = genome[12]
    max_degree = genome[13]
    mut = genome[3]
    muw = genome[4]
    beta = genome[5]
    com_size_min = genome[10]
    com_size_max = genome[11]
    tau = genome[6]
    tau2 = genome[7]
    overlapping_nodes = genome[14]
    overlap_membership = genome[15]
    pct_pos_in = genome[8]
    pct_pos_out = genome[9]
    matrix = create_matrix(num_nodes=9, average_k=average_k, max_degree=max_degree, mut=mut, muw=muw, beta=beta, com_size_min=com_size_min, com_size_max=com_size_max, tau=tau, tau2=tau2, overlapping_nodes=overlapping_nodes, overlap_membership=overlap_membership, pct_pos_in=pct_pos_in, pct_pos_out=pct_pos_out)
    return matrix


def main(genome):
    interactions = create_matrix_unformatted(genome)
    with open("interaction_matrix.dat", "w") as f:
        wr = csv.writer(f)
        wr.writerows(interactions)
    visualize_network('interaction_matrix.dat', 'interaction_matrix.png')


if __name__ == '__main__':
    main()