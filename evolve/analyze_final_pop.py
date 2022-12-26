import json
import ast
import numpy as np
import matplotlib.pyplot as plt


def histograms(population):
    diffusion = [x[0] for x in population]
    seeding = [x[1] for x in population]
    clear = [x[2] for x in population]

    figure, axis = plt.subplots(3, 1, figsize=(9, 18))
    axis[0].hist(diffusion)
    axis[1].hist(seeding)
    axis[2].hist(clear)
    axis[0].set_xlabel('diffusion')
    axis[1].set_xlabel('seeding')
    axis[2].set_xlabel('clear')

    plt.savefig('histograms.png')


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


def get_final_pop():
    final_pop = open('final_population').read()
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


if __name__ == '__main__':
    population, fitnesses = get_final_pop()