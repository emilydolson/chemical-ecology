import json
import ast
import csv
import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from analyze_final_pop import get_final_pop, get_avg_fitness

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

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
    figure.supxlabel('value')
    figure.supylabel('count')
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
    figure.supxlabel('value')
    figure.supylabel('count')
    figure.legend()

    plt.savefig(f'{file_name}_histograms_scores.png')


def avg_fitnesses_over_time(avg_fitnesseses, file_name):
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
        axis[j%5][k%2].set_xlim(0, len(avg_fitnesseses))
        k += 1
    figure.suptitle(file_name)
    figure.legend()

    plt.savefig(f'{file_name}_avg_fitnesses.png')


def main(file_name):
    file_path = f'/mnt/gs21/scratch/leithers/chemical-ecology/data/{file_name}'
    populations = []
    fitnesses = []
    avg_fitnesses = []
    for f in os.listdir(file_path):
        full_path = os.path.join(file_path, f)
        population, fitness = get_final_pop(full_path+'/evolve')
        populations.append(population)
        fitnesses.append(fitness)
        avg_fitnesses.append(get_avg_fitness(full_path+'/evolve'))
    histograms_params(populations, file_name)
    histograms_scores(fitnesses, file_name)
    avg_fitnesses_over_time(avg_fitnesses, file_name)


if __name__ == '__main__':
    main('results_double_with_fitness')