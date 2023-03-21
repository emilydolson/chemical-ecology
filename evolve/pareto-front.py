import numpy as np
import sys
import os
from analyze_final_pop import get_final_pop, histogram_params
import matplotlib.pyplot as plt
from copy import deepcopy
from klemm_lfr_graph import create_matrix_unformatted
import matplotlib
from matplotlib.colors import BoundaryNorm
from ga import calc_all_fitness


def diffusion_correlation(pops):
    matrices = []
    for genome in pops:
        matrices.append(sum([sum(abs(x)) for x in create_matrix_unformatted(genome)]))

    figure, axis = plt.subplots(1, 1)

    axis.scatter([x[0] for x in pops], matrices)
    #figure.suptitle(file_name)
    figure.supxlabel('Diffusion')
    figure.supylabel('Sum of Matrix Weights')

    print(np.corrcoef([x[0] for x in pops], matrices))

    plt.savefig('test.png')


def matrices(genome, i, j):
    matrix = create_matrix_unformatted(genome)

    figure, axis = plt.subplots(1, 1)
    #https://stackoverflow.com/questions/23994020/colorplot-that-distinguishes-between-positive-and-negative-values
    cmap = plt.get_cmap('PuOr')
    cmaplist = [cmap(i) for i in range(cmap.N)]
    cmap = cmap.from_list('Custom cmap', cmaplist, cmap.N)
    bounds = [-1, -0.75, -0.5, -0.25, -.0001, .0001, 0.25, 0.5, 0.75, 1]
    norm = BoundaryNorm(bounds, cmap.N)

    figure.suptitle(genome)
    b = axis.imshow(matrix, cmap=cmap, interpolation='none', norm=norm)
    plt.colorbar(b)
    plt.savefig(f'matrices/{i}_{j}.png')


def find_pareto(pops, fitnesses):
    fitness_names = [x for x in fitnesses[0].keys() if 'Score' in x]
    for index, fit in enumerate(fitnesses):
        for fit2 in fitnesses:
            if(fit['Biomass_Score'] < fit2['Biomass_Score'] and fit['Growth_Rate_Score'] < fit2['Growth_Rate_Score'] and fit['Heredity_Score'] < fit2['Heredity_Score'] and fit['Invasion_Ability_Score'] < fit2['Invasion_Ability_Score'] and fit['Resiliance_Score'] < fit2['Resiliance_Score']):
                pops[index] = 0
    final_pops = []
    final_fit = []
    for i in range(len(pops)):
        if pops[i] != 0 and fitnesses[i]['Biomass'] > 5000:
            pos_count = 0
            for f in fitness_names:
                if fitnesses[i][f] > 0:
                    pos_count += 1
            if pos_count > 1:
                final_pops.append(pops[i])
                final_fit.append(fitnesses[i])
    return final_pops, final_fit


def main():
    user = os.getlogin()
    root_path = f"/mnt/home/{user}/grant/chemical-ecology/results/"
    pops = []
    fitnesses = []
    for i in range(10):
        full_path = root_path + str(i)
        population, fitness = get_final_pop(full_path+'/evolve')
        pops.append(population)
        fitnesses.append(fitness)
    pareto_pops = []
    for i in range(10):
        ppops, pfits = find_pareto(deepcopy(pops[i]), deepcopy(fitnesses[i]))
        pareto_pops.append(ppops)
    for i, replicate in enumerate(pareto_pops):
        print(len(pareto_pops[i]), "genomes found on the pareto front for replicate", i)
    all_pops = [val for sublist in pareto_pops for val in sublist]

    _, node_cnt, matrix_sum = calc_all_fitness(all_pops)
    print(np.corrcoef(node_cnt, matrix_sum))

    #for i in range(len(pareto_pops)):
    #    for j in range(len(pareto_pops[i])):
    #        matrices(pareto_pops[i][j], i, j)
    
    #histogram_params(all_pops, 'results')
    #diffusion_correlation(all_pops)


if __name__ == '__main__':
    main()