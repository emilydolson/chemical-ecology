import numpy as np
import sys
import os
from analyze_final_pop import get_final_pop
import matplotlib.pyplot as plt

def histograms_params(populations):
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
    figure.suptitle("Parameters of the pareto front")
    figure.legend()

    plt.savefig('histograms_params.png')


# Based on https://stackoverflow.com/a/40239615
def find_pareto(pops, fitnesses):
    for index, fit in enumerate(fitnesses):
        for fit2 in fitnesses:
            if(fit['Biomass_Score'] < fit2['Biomass_Score'] and fit['Growth_Rate_Score'] < fit2['Growth_Rate_Score'] and fit['Heredity_Score'] < fit2['Heredity_Score'] and fit['Invasion_Ability_Score'] < fit2['Invasion_Ability_Score'] and fit['Resiliance_Score'] < fit2['Resiliance_Score']):
                pops[index] = 0
    pops = [i for i in pops if i != 0]
    return pops

def main():
    user = os.getlogin()
    root_path = f"/mnt/home/{user}/chemical-ecology/results/"
    pops = []
    fitnesses = []
    for i in range(10):
        full_path = root_path + str(i)
        population, fitness = get_final_pop(full_path+'/evolve')
        pops.append(population)
        fitnesses.append(fitness)
    pareto_pops = []
    for i in range(10):
        pareto_pops.append(find_pareto(pops[i], fitnesses[i]))
    for i, replicate in enumerate(pareto_pops):
        print(len(pareto_pops[i]), "genomes found on the pareto front for replicate", i)
    histograms_params(pareto_pops)
        

if __name__ == '__main__':
    main()