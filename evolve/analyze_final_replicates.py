import json
import ast
import csv
import os
import numpy as np
import sys
import matplotlib
import matplotlib.pyplot as plt
from analyze_final_pop import get_final_pop, get_avg_fitness, create_matrix_unformatted, get_max_fitness

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

def histograms_fitnesses(fitnesses, file_name):
    figure, axis = plt.subplots(5, 2, figsize=(10,10))

    fitness_names = ['Biomass', 'Growth_Rate', 'Heredity']
    k = 0
    for j in range(len(fitnesses)):
        for i in range(len(fitness_names)):
            data = [x[fitness_names[i]] for x in fitnesses[j]]
            axis[j%5][k%2].hist([d/max(data) for d in data], bins=np.arange(0, 1 + 0.1, 0.1), stacked=True, label=f'{fitness_names[i]} {round(max(data),2)}', alpha=0.66)
        axis[j%5][k%2].set_xlim(0, 1)
        axis[j%5][k%2].set_ylim(0, len(fitnesses[j]))
        axis[j%5][k%2].legend(loc='upper center')
        k += 1
    figure.suptitle(file_name)
    # figure.supxlabel('value (normalized)')
    # figure.supylabel('count')

    plt.savefig(f'{file_name}_histograms_fitnesses.png')


def histograms_mparams(populations, file_name):
    figure, axis = plt.subplots(5, 2, figsize=(10,10))

    param_names = ['clique_linkage', 'pct_pos_in', 'pct_pos_out', 'muw']
    k = 0
    for j in range(len(populations)):
        for i in range(len(param_names)):
            label = None if j > 0 else param_names[i]
            axis[j%5][k%2].hist([x[i+3] for x in populations[j]], bins=np.arange(0, 1 + 0.05, 0.05), stacked=True, label=label, alpha=0.66)
        axis[j%5][k%2].set_xlim(0, 1)
        axis[j%5][k%2].set_ylim(0, len(populations[j]))
        k += 1
    figure.suptitle(file_name)
    # figure.supxlabel('value')
    # figure.supylabel('count')
    figure.legend()

    plt.savefig(f'{file_name}_histograms_mparams.png')


def histograms_params(populations, file_name):
    figure, axis = plt.subplots(5, 2, figsize=(10,10))

    param_names = ['diffusion', 'seeding', 'clear']
    k = 0
    for j in range(len(populations)):
        for i in range(len(param_names)):
            label = None if j > 0 else param_names[i]
            axis[j%5][k%2].hist([x[i] for x in populations[j]], bins=np.arange(0, 1 + 0.05, 0.05), stacked=True, label=label, alpha=0.66)
        axis[j%5][k%2].set_xlim(0, 1)
        axis[j%5][k%2].set_ylim(0, len(populations[j]))
        k += 1
    figure.suptitle(file_name)
    # figure.supxlabel('value')
    # figure.supylabel('count')
    figure.legend()

    plt.savefig(f'{file_name}_histograms_params.png')


def histograms_scores(fitnesses, file_name):
    figure, axis = plt.subplots(5, 2, figsize=(10,10))
    fitness_names = [x for x in fitnesses[0][0].keys() if 'Score' in x]
    k = 0
    for j in range(len(fitnesses)):
        for i in range(len(fitness_names)):
            label = None if j > 0 else fitness_names[i]
            h = np.histogram([x[fitness_names[i]] for x in fitnesses[j]], bins=(-1, 0, 1e-100, 1))
            axis[j%5][k%2].bar(range(3), h[0], label=label, alpha=0.66, width=1)
            axis[j%5][k%2].set_xticks(np.arange(0, 3, 1), ['[-1, 0)', '0', '(0, 1]'])
        axis[j%5][k%2].set_ylim(0, len(fitnesses[j]))
        k += 1
    figure.suptitle(file_name)
    # figure.supxlabel('value')
    # figure.supylabel('count')
    figure.legend()

    plt.savefig(f'{file_name}_histograms_scores.png')


def x_fitnesses_over_time(avg_fitnesseses, file_name, max_or_avg, folder):
    figure, axis = plt.subplots(5, 2, figsize=(10,10))
    fitness_names = [x for x in avg_fitnesseses[0][0].keys() if 'Score' in x]
    k = 0
    for j,avg_fitnesses in enumerate(avg_fitnesseses):
        avg_fitness_sep = [[] for _ in range(5)]
        for fitnesses in avg_fitnesses:
            for i in range(len(fitness_names)):
                avg_fitness_sep[i].append(fitnesses[fitness_names[i]])
        
        for i in range(len(fitness_names)):
            label = None if j > 0 else fitness_names[i]
            axis[j%5][k%2].plot(avg_fitness_sep[i], label=label)
        axis[j%5][k%2].set_xlim(0, len(avg_fitnesses))
        axis[j%5][k%2].title.set_text(folder[j])
        k += 1
    figure.suptitle(file_name)
    figure.legend()

    plt.savefig(f'{file_name}_{max_or_avg}_fitnesses.png')


def get_adaptive_genomes(populations, fitnesses, print_genomes=True):
    fitness_names = [x for x in fitnesses[0][0].keys() if 'Score' in x]
    num_adaptive = 0
    for r in range(len(populations)):
        unique_genomes = []
        for g in range(len(populations[r])):
            num_adaptive_local = 0
            for f in fitness_names:
                if fitnesses[r][g][f] > 0:
                    num_adaptive_local += 1
            if num_adaptive_local >= 4:
                #matrix = create_matrix_unformatted(populations[r][g])
                #if sum([sum(x) for x in matrix]) != 0:
                if print_genomes:
                    print(populations[r][g], {x: fitnesses[r][g][x] for x in fitnesses[r][g] if 'Score' in x})
                num_adaptive += 1
            if populations[r][g] not in unique_genomes:
                unique_genomes.append(populations[r][g])
        print(f'Unique genomes: {len(unique_genomes)}')
    print(f'Adaptive genomes: {num_adaptive} out of {len(populations)*len(populations[0])}')


def main(user, file_name):
    file_path = f'/mnt/gs21/scratch/{user}/chemical-ecology/data/{file_name}'
    populations = []
    fitnesses = []
    avg_fitnesses = []
    max_fitnesses = []
    core_file_count = 0
    folders = []
    for f in os.listdir(file_path):
        full_path = os.path.join(os.path.join(file_path, f), 'evolve')
        if os.path.exists(os.path.join(full_path, 'final_population')):
            population, fitness = get_final_pop(full_path)
            populations.append(population)
            fitnesses.append(fitness)
            avg_fitnesses.append(get_avg_fitness(full_path))
            max_fitnesses.append(get_max_fitness(full_path))
            folders.append(f)
        for t in os.listdir(full_path):
            if 'core' in t and '.csv' not in t:
                core_file_count += 1
    print(f'Core files: {core_file_count}')
    if(len(populations) == 0):
        print("No final fitnesses found")
        exit()
    histograms_mparams(populations, file_name)
    histograms_scores(fitnesses, file_name)
    histograms_fitnesses(fitnesses, file_name)
    x_fitnesses_over_time(avg_fitnesses, file_name, 'avg', folders)
    x_fitnesses_over_time(max_fitnesses, file_name, 'max')
    get_adaptive_genomes(populations, fitnesses, False)


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])