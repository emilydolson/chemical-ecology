#https://github.com/sydleither/graphgen

import random
from graphgen.lfr_generators import klemm_as_adj

def create_matrix(num_nodes, clique_size, clique_linkage, muw, beta, pct_pos_in, pct_pos_out):
    random.seed(0)

    matrix, communities = klemm_as_adj(num_nodes=num_nodes, clique_size=clique_size, clique_linkage=clique_linkage, muw=muw, beta=beta, seed=0)

    min_all = min([min(x) for x in matrix])
    max_all = max([max(x) for x in matrix])

    for row in range(len(matrix)):
        for col in range(len(matrix)):
            if matrix[row][col] != 0:
                matrix[row][col] = (matrix[row][col] - min_all) / (max_all - min_all)
                if sum([x in communities[row] for x in communities[col]]) == 0:
                    if random.random() > pct_pos_in:
                        matrix[row][col] = -matrix[row][col]
                else:
                    if random.random() > pct_pos_out:
                        matrix[row][col] = -matrix[row][col]

    return matrix


def create_matrix_unformatted(genome):
    diffusion = genome[0]
    seeding = genome[1]
    clear = genome[2]
    clique_linkage = genome[3]
    pct_pos_in = genome[4]
    pct_pos_out = genome[5]
    muw = genome[6]
    beta = genome[7]
    clique_size = genome[8]
    interaction_matrix = create_matrix(num_nodes=9, clique_size=clique_size, clique_linkage=clique_linkage, muw=muw, beta=beta, pct_pos_in=pct_pos_in, pct_pos_out=pct_pos_out)
    return matrix