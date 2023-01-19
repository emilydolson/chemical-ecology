import json
import ast
import csv
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from lfr_graph import create_matrix


def histograms(population, file_name):
    figure, axis = plt.subplots(1, 1)

    param_names = ['diffusion', 'seeding', 'clear']
    for i in range(3):
        plt.hist([x[i] for x in population], bins=np.arange(0, 1 + 0.05, 0.05), stacked=True, label=param_names[i], alpha=0.66)
    axis.set_xlim(0, 1)
    figure.suptitle(file_name)
    figure.legend()

    plt.savefig(f'{file_name}_histograms.png')


def histograms_scores(fitnesses, file_name):
    figure, axis = plt.subplots(1, 1)
    rng = 0.1
    fitness_names = [x for x in fitnesses[0].keys() if 'Score' in x]
    for i in range(5):
        plt.hist([x[fitness_names[i]] for x in fitnesses], bins=np.arange(-rng, rng + 0.01, 0.01), stacked=True, label=fitness_names[i], alpha=0.66)
    axis.set_xlim(-rng, rng)
    figure.suptitle(file_name)
    figure.legend()

    plt.savefig(f'{file_name}_histograms_scores.png')


def pareto_front(fitnesses):
    biomass = [x['Biomass'] for x in fitnesses]
    growth_rate = [x['Growth_Rate'] for x in fitnesses]
    heredity = [x['Heredity'] for x in fitnesses]

    #fig = plt.figure()
    #ax = plt.axes(projection="3d")
    #ax.scatter(biomass, growth_rate, heredity)

    figure, axis = plt.subplots(3, 1, figsize=(9, 18))
    axis[0].scatter(biomass, growth_rate)
    axis[1].scatter(growth_rate, heredity)
    axis[2].scatter(heredity, biomass)
    axis[0].set_xlabel('biomass')
    axis[0].set_ylabel('growth rate')
    axis[1].set_xlabel('growth rate')
    axis[1].set_ylabel('heredity')
    axis[2].set_xlabel('heredity')
    axis[2].set_ylabel('biomass')

    plt.savefig('pareto.png')


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


def check_matrices(population):
    zeros = 0
    unique = []
    for genome in population:
        matrix = create_matrix_unformatted(genome)
        if sum([sum(x) for x in matrix]) == 0:
            zeros += 1
        if genome[3:] not in unique:
            unique.append(genome[3:])
    print(f'error matrices: {zeros}')
    print(f'unique matrices: {len(unique)}')


def get_final_pop(file_location):
    final_pop = open(f'{file_location}/final_population').read()
    pop_list = []
    fitness_list = []
    for line in final_pop.split('\n'):
        if len(line) > 0:
            if line[0] == '{':
                fitnesses, genome = line.split('} ')
                fitnesses = json.loads(fitnesses.replace('\'', '\"') + '}')
                genome = ast.literal_eval(genome[1:])
                pop_list.append(genome)
                fitness_list.append(fitnesses)
    return pop_list, fitness_list


def main(file_name):
    file_path = f'/mnt/gs21/scratch/leithers/chemical-ecology/data/{file_name}/evolve'
    population, fitnesses = get_final_pop(file_path)
    histograms(population, file_name)
    histograms_scores(fitnesses, file_name)
    #check_matrices(population)

    #interactions = create_matrix_unformatted(population[0])
    #with open("interaction_matrix.dat", "w") as f:
        #wr = csv.writer(f)
        #wr.writerows(interactions)
    #visualize_network('interaction_matrix.dat', 'interaction_matrix.png')


if __name__ == '__main__':
    main('evolve_ceil_scores_only')