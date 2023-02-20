import json
import ast
import csv
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from lfr_graph import create_matrix

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

def histograms_fitnesses(fitnesses, file_name):
    figure, axis = plt.subplots(1, 1)

    fitness_names = ['Biomass', 'Growth_Rate', 'Heredity']
    for i in range(len(fitness_names)):
        data = [x[fitness_names[i]] for x in fitnesses]
        axis.hist([d/max(data) for d in data], bins=np.arange(0, 1 + 0.1, 0.1), stacked=True, label=f'{fitness_names[i]} {round(max(data),2)}', alpha=0.66)
    axis.set_xlim(0, 1)
    axis.set_ylim(0, len(fitnesses))
    figure.legend()
    figure.suptitle(file_name)
    figure.supxlabel('value (normalized)')
    figure.supylabel('count')

    plt.savefig(f'{file_name}_histogram_fitnesses.png')


def histogram_params(population, file_name):
    figure, axis = plt.subplots(1, 1)

    param_names = ['diffusion', 'seeding', 'clear']
    for i in range(3):
        axis.hist([x[i] for x in population], bins=np.arange(0, 1 + 0.05, 0.05), stacked=True, label=param_names[i], alpha=0.66)
    axis.set_xlim(0, 1)
    axis.set_ylim(0, len(population))
    figure.suptitle(file_name)
    figure.legend()

    plt.savefig(f'{file_name}_histograms_params.png')


def histogram_scores(fitnesses, file_name):
    figure, axis = plt.subplots(1, 1)
    rng = 0.1
    fitness_names = [x for x in fitnesses[0].keys() if 'Score' in x]
    for i in range(len(fitness_names)):
        h = np.histogram([x[fitness_names[i]] for x in fitnesses], bins=(-1, 0, 1e-100, 1))
        axis.bar(range(3), h[0], label=fitness_names[i], alpha=0.66, width=1)
        axis.set_xticks(np.arange(0, 3, 1), ['[-1, 0)', '0', '(0, 1]'])
    axis.set_ylim(0, len(fitnesses))
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


def get_avg_fitness(file_location):
    avg_fitnesses = open(f'{file_location}/avg_fitness.txt').read()
    avg_fitness_list = [{}]
    for line in avg_fitnesses.split('\n'):
        if len(line) < 3:
            avg_fitness_list.append({})
        else:
            name, value = line.split(': ')
            avg_fitness_list[-1][name] = float(value)
    return avg_fitness_list[:-2]


def get_max_fitness(file_location):
    avg_fitnesses = open(f'{file_location}/max_fitness.txt').read()
    avg_fitness_list = [{}]
    for line in avg_fitnesses.split('\n'):
        if len(line) < 3:
            avg_fitness_list.append({})
        else:
            name, value = line.split(': ')
            avg_fitness_list[-1][name] = float(value)
    return avg_fitness_list[:-2]
            

def avg_fitness_over_time(avg_fitnesses, file_name):
    fitness_names = [x for x in avg_fitnesses[0].keys() if 'Score' in x]
    avg_fitness_sep = [[] for _ in range(len(fitness_names))]
    for fitnesses in avg_fitnesses:
        for i in range(len(fitness_names)):
            avg_fitness_sep[i].append(fitnesses[fitness_names[i]])
    
    figure, axis = plt.subplots(1, 1)
    for i in range(len(fitness_names)):
        plt.plot(avg_fitness_sep[i], label=fitness_names[i])
    figure.suptitle(file_name)
    figure.legend(loc='lower right')

    plt.savefig(f'{file_name}_avg_fitness.png')


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
    histogram_params(population, file_name)
    histogram_scores(fitnesses, file_name)
    histograms_fitnesses(fitnesses, file_name)
    avg_fitness_over_time(get_avg_fitness(file_path), file_name)
    check_matrices(population)


if __name__ == '__main__':
    main('temp')