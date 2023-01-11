import numpy as np
import matplotlib.pyplot as plt
from ga import create_pop, calc_all_fitness
from lfr_graph import create_matrix

def main():
    fitness_types = ['Biomass', 'Growth_Rate', 'Heredity']#, 'Biomass_Score', 'Growth_Rate_Score', 'Heredity_Score', 'Invasion_Ability_Score', 'Resiliance_Score']
    population = create_pop(5)
    measures = [{'Mean':[], 'Standard Deviation':[], 'Range':[], 'Coefficient of Variation':[]} for _ in range(len(fitness_types))]
    figure, axis = plt.subplots(3, 1, figsize=(9, 18))

    for genome in population:
        print(genome)
        fitnesses = [[] for _ in range(len(fitness_types))]
        for trials in range(5):
            fitness = calc_all_fitness([genome])
            for k in range(len(fitness_types)):
                fitnesses[k].append(fitness[0][fitness_types[k]])
        for i in range(len(fitnesses)):
            measures[i]['Mean'].append(np.mean(fitnesses[i]))
            measures[i]['Standard Deviation'].append(np.std(fitnesses[i]))
            measures[i]['Range'].append(max(fitnesses[i]) - min(fitnesses[i]))
            measures[i]['Coefficient of Variation'].append(np.std(fitnesses[i])/np.mean(fitnesses[i]))
            axis[i].plot(range(len(fitnesses[i])), fitnesses[i])
    
    for i in range(len(fitness_types)):
        axis[i].set_title(fitness_types[i])
        print(fitness_types[i])
        for measure in measures[i].keys():
            print(f'Average {measure}: {np.mean(measures[i][measure])}')
        print()
        
    plt.savefig('graph.png')


if __name__ == '__main__':
    main()
    #print(create_matrix(9, 3, 5, 0.3, 0.1, 1.5, 4, 8, 2, 1, 8, 4, 0.1, 0.056))
    