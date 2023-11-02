import random


def random_matrix(ntypes, prob_edge, seed):
    random.seed(seed)
    return [[round(random.uniform(-1, 1), 3) if random.random() < prob_edge else 0 for _ in range(ntypes)] for _ in range(ntypes)]