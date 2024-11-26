from statistics import covariance
import itertools

# N_SPECIES is the number of species
# DIFFUSION LIMIT is the amount of species required to diffused

N_SPECIES = 5
DIFFUSION_LIMIT = 10

class Species:
    def __init__(self, id, interactions, self_interaction):
        self.id = id
        self.interactions = interactions
        self.self_interaction = self_interaction

    def __repr__(self):
        return f"Species {self.id}: {self.interactions}, {self.self_interaction}"

INTS = [Species(i, {j: 0 for j in range(N_SPECIES) if j != i}, 0) for i in range(N_SPECIES)]

# Calculate group level fitness (assumes that all non-zero species are part of the community and must reach
# DIFFUSION_LIMIT for community to propogate)
def calc_W(community):
    community = community[:]
    W = 0
    time_limit = 1000
    members = [i for i in range(N_SPECIES) if community[i] > 0]
    if len(members) == 0:
        return 0
    while any([community[i] < DIFFUSION_LIMIT for i in members]) and time_limit > 0:
        W += 1
        next_community = community[:]
        for i in range(N_SPECIES):
            modifier = 0
            for j in range(N_SPECIES):
                if i == j:
                    modifier += INTS[i].self_interaction * community[i]
                else:
                    modifier += INTS[j].interactions[i] * community[j]

            # print(i, modifier)
            next_community[i] += next_community[i] * modifier
            next_community[i] = min(DIFFUSION_LIMIT, next_community[i])
            next_community[i] = max(0, next_community[i])
        community = next_community
        # print(community)
        time_limit -= 1
    if W == 1000:
        return 0
    return 100/W


# Calculate individual level fitness (assumes that only the individual must reach DIFFUSION_LIMIT for it to propogate)
def calc_ind_W(community, ind):
    community = community[:]
    W = 0
    time_limit = 1000
    if community[ind] == 0:
        return 0
    while community[ind] < DIFFUSION_LIMIT and time_limit > 0:
        W += 1
        next_community = community[:]
        for i in range(N_SPECIES):
            modifier = 0
            for j in range(N_SPECIES):
                if i == j:
                    modifier += INTS[i].self_interaction * community[i]
                else:
                    modifier += INTS[j].interactions[i] * community[j]
            next_community[i] += next_community[i] * modifier
            next_community[i] = min(DIFFUSION_LIMIT, next_community[i])
            next_community[i] = max(0, next_community[i])            
        community = next_community
        time_limit -= 1
    if W == 1000:
        return 0
    return 100/W

# Calculate the price equation for a population
def price_equation(population):
    W_vec = [calc_W(community) for community in population for ind in range(N_SPECIES) for _ in range(community[ind])]
    delta_w_vec = [calc_ind_W(community, ind) - calc_W(community) for community in population for ind in range(N_SPECIES) for _ in range(community[ind])]
    g_vecs = []

    # print(W_vec)
    # print(delta_w_vec)
    for i in range(len(population)):
        print(population[i], calc_W(population[i]), [calc_ind_W(population[i], ind) - calc_W(population[i]) for ind in range(N_SPECIES)])

    for trait in range(N_SPECIES):
        g_vecs.append([INTS[ind].interactions[trait] if trait != ind else 0 for community in population for ind in range(N_SPECIES) for _ in range(community[ind])])
        # print(g_vecs[-1])
    g_vecs.append([INTS[ind].self_interaction for community in population for ind in range(N_SPECIES) for _ in range(community[ind])])

    results = {}
    for trait in range(N_SPECIES+1):
        between_groups = covariance(W_vec, g_vecs[trait])
        within_groups = covariance(delta_w_vec, g_vecs[trait])
        results[trait] = (between_groups, within_groups)

    return results



# Modify this to set up interaction nmatrix
# INTS[0][0] = 0
# INTS[1][0] = -1
# INTS[2][0] = -1

# INTS[0][1] = -1
# INTS[1][1] = 0
# INTS[2][1] = -1

# INTS[0][2] = -1
# INTS[1][2] = -1
# INTS[2][2] = 0

# INTS[3][4] = 1
# INTS[4][3] = 1
# INTS[0][3] = -.1
# # INTS[3][0] = -1
# INTS[0][N_SPECIES] = .5
# INTS[1][N_SPECIES] = .5
# INTS[2][N_SPECIES] = .5
# INTS[3][N_SPECIES] = 0
# INTS[4][N_SPECIES] = 0

INTS[0].interactions[1] = -1
INTS[0].interactions[2] = -1
INTS[0].interactions[3] = -1
INTS[0].interactions[4] = -1

INTS[1].interactions[0] = -1
INTS[1].interactions[2] = -1
INTS[1].interactions[3] = -1
INTS[1].interactions[4] = -1

# INTS[2].interactions[0] = 1
# INTS[2].interactions[1] = 1
# INTS[2].interactions[3] = 1
# INTS[2].interactions[4] = 1

# INTS[3].interactions[0] = -1
# INTS[3].interactions[1] = -1
INTS[3].interactions[4] = 1
INTS[4].interactions[3] = 1


INTS[0].self_interaction = .5
INTS[1].self_interaction = .5
INTS[2].self_interaction = .5

print(INTS)

# Modify this to set up population
# population = [[1, 0, 0], [0, 1, 1], [1, 0, 1], [1, 1, 1], [0, 1, 0], [0, 0, 1], [1, 1, 0]]
# population = [[1, 0, 0, 0, 0, 0, 0, 0, 0], [0, 1, 1, 0, 0, 0, 0, 0, 0], [1, 1, 1, 1, 0, 0, 0, 0, 0]]
population = [list(i) for i in itertools.product([0, 1], repeat=N_SPECIES)]
# print(population)
print(price_equation(population))

# print(calc_W([1, 1, 1, 0, 0]))