import numpy as np
import random


#from https://www.pnas.org/doi/full/10.1073/pnas.2005759118#data-availability
MOTIF_DICT = {0:[[0, 0, 0], [0, 0, 0], [0, 0, 0]],
              1:[[0, -1, 0], [0, 0, -1], [-1, 0, 0]], #intransitive competition
              2:[[0, 1, 1], [0, 0, -1], [0, -1, 0]], #facilitation-driven competition
              3:[[0, 0, 0], [-1, 0, 1], [-1, 1, 0]], #competition-driven facilitation
              4:[[0, -1, -1], [0, 0, 1], [0, 1, 0]], #competition-driven facilitation
              5:[[0, 1, 0], [0, 0, 1], [1, 0, 0]], #intransitive facilitation
              6:[[0, 1, -1], [0, 0, -1], [0, 0, 0]], #module 6
              7:[[0, -1, -1], [0, 0, 1], [0, 0, 0]], #module 7
              8:[[0, 1, 1], [0, 0, 1], [0, 0, 0]], #facilitation cascade
              9:[[0, -1, -1], [0, 0, -1], [0, -1, 0]], #dominance
              10:[[0, -1, -1], [0, 0, 0], [-1, -1, 0]], #codominance
              11:[[0, -1, -1], [0, 0, -1], [-1, 0, 0]], #quasiintransitive competition
              12:[[0, -1, -1], [-1, 0, -1], [-1, -1, 0]], #complete competition
              13:[[0, -1, -1], [0, 0, -1], [-1, -1, 0]]} #quasicomplete competition
LOW_RANGE_DICT = {0: 0, 1: 0.33, 2:0.66}
HIGH_RANGE_DICT = {0: 0.33, 1:0.66, 2:1}


def modify_strength(motif, strength, seed):
    random.seed(seed)

    new_motif = [[0 for _ in range(len(motif))] for _ in range(len(motif))]
    for row in range(len(motif)):
        for col in range(len(motif)):
            new_motif[row][col] = motif[row][col] * random.uniform(LOW_RANGE_DICT[strength], HIGH_RANGE_DICT[strength])
    return new_motif


#motifs: 0 to 13 in list of length 9
#strengths: 0 to 2 in list of length 9
def create_matrix(num_nodes, chunk_size, motifs, strengths, seed=0):
    random.seed(seed)
    interactions_temp = [[] for _ in range(chunk_size)]
    for i in range(chunk_size):
        interactions_temp[i] = modify_strength(MOTIF_DICT[motifs[i]], strengths[i], random.randint(0,100))
    for i in range(chunk_size, num_nodes):
        interactions_temp[i%3] = np.vstack([interactions_temp[i%3], modify_strength(MOTIF_DICT[motifs[i]], strengths[i], random.randint(0,100))])
    interactions = np.hstack([x for x in interactions_temp])
    return interactions