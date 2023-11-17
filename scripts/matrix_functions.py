import random
import networkx as nx


def scale_free(ntypes, split1, split2, delta_in, delta_out, seed):
    lower_split = split1 if split1 < split2 else split2
    higher_split = split1 if split1 > split2 else split2
    alpha = lower_split
    beta = higher_split - lower_split
    gamma = 1 - higher_split
    structure_nx = nx.scale_free_graph(ntypes, alpha, beta, gamma, delta_in, delta_out, seed)
    matrix = nx.to_numpy_array(structure_nx).tolist()
    for row in range(len(matrix)):
        for col in range(len(matrix)):
            if matrix[row][col] != 0:
                matrix[row][col] = round(random.uniform(-1, 1), 3)
    return matrix


def erdos_renyi_normal(ntypes, prob_edge, mean, var, seed):
    structure_nx = nx.erdos_renyi_graph(n=ntypes, p=prob_edge, seed=seed, directed=True)
    matrix = nx.to_numpy_array(structure_nx).tolist()
    for row in range(len(matrix)):
        for col in range(len(matrix)):
            if matrix[row][col] != 0:
                val = round(random.gauss(mean, var), 3)
                val = 1 if val > 1 else val
                val = -1 if val < -1 else val
                matrix[row][col] = val
    return matrix


def erdos_renyi_uniform(ntypes, prob_edge, seed):
    structure_nx = nx.erdos_renyi_graph(n=ntypes, p=prob_edge, seed=seed, directed=True)
    matrix = nx.to_numpy_array(structure_nx).tolist()
    for row in range(len(matrix)):
        for col in range(len(matrix)):
            if matrix[row][col] != 0:
                matrix[row][col] = round(random.uniform(-1, 1), 3)
    return matrix


def random_matrix(ntypes, prob_edge, seed):
    random.seed(seed)
    return [[round(random.uniform(-1, 1), 3) if random.random() < prob_edge else 0 for _ in range(ntypes)] for _ in range(ntypes)]