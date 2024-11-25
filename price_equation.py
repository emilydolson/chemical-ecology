from statistics import covariance

# N_SPECIES is the number of species
# DIFFUSION LIMIT is the amount of species required to diffused

N_SPECIES = 3
DIFFUSION_LIMIT = 10
INTS = [[0 for i in range(N_SPECIES)] for i in range(N_SPECIES)]

# Calculate group level fitness (assumes that all non-zero species are part of the community and must reach
# DIFFUSION_LIMIT for community to propogate)
def calc_W(community):
    community = community[:]
    W = 0
    time_limit = 1000
    while any([i > 0 and i < DIFFUSION_LIMIT for i in community]) and time_limit > 0:
        W += 1
        for i in range(N_SPECIES):
            for j in range(N_SPECIES):
                community[j] += INTS[i][j] * community[i]
        time_limit -= 1
    if W == 1000:
        return 0
    return 1/W


# Calculate individual level fitness (assumes that only the individual must reach DIFFUSION_LIMIT for it to propogate)
def calc_ind_W(community, ind):
    community = community[:]
    W = 0
    time_limit = 1000
    while community[ind] > 0 and community[ind] < DIFFUSION_LIMIT and time_limit > 0:
        W += 1
        for i in range(N_SPECIES):
            community[i] += INTS[ind][i] * community[ind]
        time_limit -= 1
    if W == 1000:
        return 0
    return 1/W

# Calculate the price equation for a population
def price_equation(population):
    W_vec = [calc_W(community) for community in population for ind in range(N_SPECIES) for _ in range(community[ind])]
    delta_w_vec = [calc_ind_W(community, ind) - calc_W(community) for community in population for ind in range(N_SPECIES) for _ in range(community[ind])]
    g_vecs = []

    print(W_vec)
    print(delta_w_vec)

    for trait in range(N_SPECIES):
        g_vecs.append([INTS[ind][trait] for community in population for ind in range(N_SPECIES) for _ in range(community[ind])])
        print(g_vecs[-1])

    results = {i:(0,0) for i in range(N_SPECIES)}
    for trait in range(N_SPECIES):

        between_groups = covariance(W_vec, g_vecs[trait])
        within_groups = covariance(delta_w_vec, g_vecs[trait])
        results[trait] = (between_groups, within_groups)

    return results

# Modify this to set up interaction nmatrix
INTS[0][0] = 1
INTS[0][1] = -1
INTS[0][2] = -1
INTS[1][0] = 0
INTS[1][1] = 0
INTS[1][2] = 1
INTS[2][0] = 0
INTS[2][1] = 1
INTS[2][2] = 0

# Modify this to set up population
# population = [[1, 0, 0], [0, 1, 1], [1, 0, 1], [1, 1, 1], [0, 1, 0], [0, 0, 1], [1, 1, 0]]
population = [[1, 0, 0], [0, 1, 1]]
print(price_equation(population))