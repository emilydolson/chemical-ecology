import subprocess
import pandas as pd
import random
import sys
import csv
import numpy as np
from copy import deepcopy
from klemm_lfr_graph import create_matrix

#REPRO_THRESHOLD 10000000000 -MAX_POP 10000 -WORLD_X 30 -WORLD_Y 30 -N_TYPES 9 -UPDATES 10000 <-- Used in poc
N_TYPES = 9
WORLD_X = 10
WORLD_Y = 10
UPDATES = 1000
MAX_POP = 10000
REPRO_THRESHOLD = 10000000000

def create_pop(size):
    pop = []
    for _ in range(size):
        diffusion = round(random.random(), 3)
        seeding = round(random.random(), 3)
        clear = round(random.random(), 3)
        clique_linkage = round(random.random(), 3)
        pct_pos_in = round(random.random(), 3)
        pct_pos_out = round(random.random(), 3)
        muw = round(random.random(), 3)
        beta = random.randint(1, 2)
        clique_size = random.randint(2, N_TYPES-1)
        genome = [diffusion, seeding, clear, clique_linkage, pct_pos_in, pct_pos_out, muw, beta, clique_size]
        pop.append(genome)
    return pop


def euclidean_distance(genome_x, genome_y):
    return np.linalg.norm(np.array(genome_x[0:6] + [b/10 for b in genome_x[6:]]) - np.array(genome_y[0:6] + [b/10 for b in genome_y[6:]]))


def calc_all_fitness(population, niched=False):
    fitness_lst = []
    for genome in population:
        diffusion = genome[0]
        seeding = genome[1]
        clear = genome[2]
        clique_linkage = genome[3]
        pct_pos_in = genome[4]
        pct_pos_out = genome[5]
        muw = genome[6]
        beta = genome[7]
        clique_size = genome[8]
        interaction_matrix = create_matrix(num_nodes=N_TYPES, clique_size=clique_size, clique_linkage=clique_linkage, muw=muw, beta=beta, pct_pos_in=pct_pos_in, pct_pos_out=pct_pos_out)
        interaction_matrix_file = 'interaction_matrix.dat'
        with open(interaction_matrix_file, 'w') as f:
            wr = csv.writer(f)
            wr.writerows(interaction_matrix)
        #run chem eco three times, store the fitness and assembly values for each run, then multiply the scores togetehr and average the fitnesses
        temp_fits = {
            "biomass": [],
            "growth_rate": [], 
            "heredity": [], 
            "biomass_score": [],
            "growth_rate_score": [], 
            "invasion_score": [], 
            "heredity_score": [], 
            "resiliance_score": [], 
            #assembly score is not a fitness value. It is used to calc fitness values. 
            "assembly_score": []
        }
        for _ in range(3):
            chem_eco = subprocess.Popen(
                [(f'../chemical-ecology '
                f'-DIFFUSION {diffusion} '
                f'-SEEDING_PROB {seeding} '
                f'-PROB_CLEAR {clear} ' 
                f'-INTERACTION_SOURCE {interaction_matrix_file} '
                f'-REPRO_THRESHOLD {REPRO_THRESHOLD} '
                f'-MAX_POP {MAX_POP} '
                f'-WORLD_X {WORLD_X} '
                f'-WORLD_Y {WORLD_Y} '
                f'-UPDATES {UPDATES} '
                f'-N_TYPES {N_TYPES}')],
                shell=True, 
                stdout=subprocess.DEVNULL)
            return_code = chem_eco.wait()
            if return_code != 0:
                print("Error in a-eco, return code:", return_code)
                sys.stdout.flush()
            df = pd.read_csv('a-eco_data.csv')
            biomasses = df['mean_Biomass'].values
            growth_rates = df['mean_Growth_Rate'].values
            heredities = df['mean_Heredity'].values
            df2 = pd.read_csv('scores.csv')
            b_score = df2['Biomass_Score'].values
            g_score = df2['Growth_Rate_Score'].values
            h_score = df2['Heredity_Score'].values
            i_score = df2['Invasion_Ability_Score'].values
            r_score = df2['Resiliance_Score'].values
            a_score =  df2['Assembly_Score'].values
            temp_fits["biomass"].append(biomasses[-1] - biomasses[0])
            temp_fits["growth_rate"].append(growth_rates[-1] - growth_rates[0])
            temp_fits["heredity"].append(heredities[-1])
            temp_fits["biomass_score"].append(b_score[0])
            temp_fits["growth_rate_score"].append(g_score[0])
            temp_fits["heredity_score"].append(h_score[0])
            temp_fits["invasion_score"].append(i_score[0])
            temp_fits["resiliance_score"].append(r_score[0])
            temp_fits["assembly_score"].append(a_score[0])
        new_fitness = {
            'Biomass': np.mean(temp_fits['biomass']),
            'Growth_Rate': np.mean(temp_fits['growth_rate']),
            'Heredity': np.mean(temp_fits['heredity']),
            'Biomass_Score': np.prod(temp_fits['biomass_score']) - np.prod(temp_fits['assembly_score']),
            'Growth_Rate_Score': np.prod(temp_fits['growth_rate_score']) - np.prod(temp_fits['assembly_score']),
            'Heredity_Score': np.prod(temp_fits['heredity_score']) - np.prod(temp_fits['assembly_score']),
            'Invasion_Ability_Score': np.prod(temp_fits['invasion_score']) - np.prod(temp_fits['assembly_score']),
            'Resiliance_Score': np.prod(temp_fits['resiliance_score']) - np.prod(temp_fits['assembly_score'])
        }
        '''if np.prod(temp_fits['assembly_score']) != 0:
            new_fitness = {
                'Biomass': np.mean(temp_fits['biomass']),
                'Growth_Rate': np.mean(temp_fits['growth_rate']),
                'Heredity': np.mean(temp_fits['heredity']),
                'Biomass_Score': np.prod(temp_fits['biomass_score']) / np.prod(temp_fits['assembly_score']),
                'Growth_Rate_Score': np.prod(temp_fits['growth_rate_score']) / np.prod(temp_fits['assembly_score']),
                'Heredity_Score': np.prod(temp_fits['heredity_score']) / np.prod(temp_fits['assembly_score']),
                'Invasion_Ability_Score': np.prod(temp_fits['invasion_score']) / np.prod(temp_fits['assembly_score']),
                'Resiliance_Score': np.prod(temp_fits['resiliance_score']) / np.prod(temp_fits['assembly_score'])
            }
        else:
            new_fitness = {
                'Biomass': np.mean(temp_fits['biomass']),
                'Growth_Rate': np.mean(temp_fits['growth_rate']),
                'Heredity': np.mean(temp_fits['heredity']),
                'Biomass_Score': 0,
                'Growth_Rate_Score': 0,
                'Heredity_Score': 0,
                'Invasion_Ability_Score': 0,
                'Resiliance_Score': 0
            }'''
        fitness_lst.append(new_fitness)
    #niching / fitness sharing
    #based on formula in https://www.sciencedirect.com/science/article/pii/S0304397518304882
    if niched:
        niched_fitness_lst = []
        sigma = len(population)/2
        alpha = 1
        for i in range(len(population)):
            genome_x = population[i]
            fitnesses = fitness_lst[i]
            sharing = sum([max(0, 1 - (euclidean_distance(genome_x, genome_y)/sigma)**(alpha)) for genome_y in population])
            new_fitness = {
                'Biomass': fitnesses['Biomass'] / sharing,
                'Growth_Rate': fitnesses['Growth_Rate'] / sharing,
                'Heredity': fitnesses['Heredity'] / sharing, 
                'Biomass_Score': fitnesses['Biomass_Score'] / sharing,
                'Growth_Rate_Score': fitnesses['Growth_Rate_Score'] / sharing,
                'Heredity_Score': fitnesses['Heredity_Score'] / sharing,
                'Invasion_Ability_Score': fitnesses['Invasion_Ability_Score'] / sharing,
                'Resiliance_Score': fitnesses['Resiliance_Score'] / sharing
            }
            niched_fitness_lst.append(new_fitness)
        return niched_fitness_lst
    return fitness_lst


def lexicase_select(pop, all_fitness, test_cases):
    random.shuffle(test_cases)
    while len(test_cases) > 0:
        # Select which test case to use 
        first = test_cases[0]
        # Initialize best to a negative value so it will lose the first comparison
        best = -9999999
        candidates = []
        for index, fitness in enumerate(all_fitness):
            if fitness[first] > best:
                best = fitness[first]
                # Save the index of the current best
                candidates = [index]
            # Also save the index of any ties
            elif fitness[first] == best:
                candidates.append(index)
        # If we have a winner
        if len(candidates) == 1:
            winner_index = candidates[0]
            winner = pop[winner_index]
            pop.pop(winner_index)
            all_fitness.pop(winner_index)
            # Return the parent, and a pop and fitness list without the winner
            return (winner, pop, all_fitness)
        # If we only have one test case left choose randomly from the ties
        if len(test_cases) == 1:
            winner_index = random.choice(candidates)
            winner = pop[winner_index]
            pop.pop(winner_index)
            all_fitness.pop(winner_index)
            # Return the parent, and a pop and fitness list without the winner
            return (winner, pop, all_fitness)
        # If there were ties and more than one test case, move to the next test case
        test_cases.pop(0)
    print("We should never get here. Find Max")


def crossover(genome1, genome2):
    # If we crossover, swap the world parameters and matrix parameters
    if (random.random() < .9):
        child1 = genome1[0:3] + genome2[3:]
        child2 = genome2[0:3] + genome1[3:]
        return (child1, child2)
    else:
        return(genome1, genome2)


# https://github.com/DEAP/deap/blob/master/deap/tools/mutation.py
def mutate(pop):
    for genome in pop:
        for i, param in enumerate(genome):
            if random.random() < 1/len(genome):
                if i == 1:
                    param += random.gauss(0, 0.1)
                    if param < 0.01:
                        param = 0.01
                    if param > 1:
                        param = 1
                if i < 4 and i != 1:
                    param += random.gauss(0, 0.1)
                    if param < 0:
                        param = 0
                    if param > 1:
                        param = 1
                if i == 4 or i == 5: #pct pos
                    param += random.gauss(0, 0.1)
                    if param < 0:
                        param = 0
                    if param > 0.75:
                        param = 0.75
                if i == 6: #muw
                    param += random.gauss(0, 0.1)
                    if param < 0:
                        param = 0
                    if param > 1:
                        param = 1
                elif i == 7: #beta
                    param += random.gauss(0, 0.1)
                    if param < 0:
                        param = 0
                    if param > 3:
                        param = 3
                elif i == 8: #clique size
                    param += random.choice([-1, 1])
                    if param < 2:
                        param = 2
                    if param > N_TYPES-1:
                        param = N_TYPES-1
            genome[i] = round(param, 3)
    return pop


def run():
    pop_size = 100
    generations = 250
    population = create_pop(pop_size)
    test_cases = ['Biomass', 'Growth_Rate', 'Heredity', 'Biomass_Score', 'Growth_Rate_Score', 'Heredity_Score', 'Invasion_Ability_Score', 'Resiliance_Score']
    #test_cases = ['Biomass_Score', 'Growth_Rate_Score', 'Heredity_Score', 'Invasion_Ability_Score', 'Resiliance_Score']
    for gen in range(generations):
        all_fitness = calc_all_fitness(population, False)
        track_avg_fitness(all_fitness, test_cases)
        track_max_fitness(all_fitness, test_cases)
        parent_tuple = (None, population, all_fitness)
        parents = []
        for _ in range(pop_size//2):
            #Returns (winner, pop, all_fitness)
            parent_tuple = lexicase_select(parent_tuple[1], parent_tuple[2], deepcopy(test_cases))
            parents.append(parent_tuple[0])
        #Crossover
        new_population = []
        for _ in range(pop_size//2):
            children = crossover(random.choice(parents), random.choice(parents))
            new_population.append(children[0])
            new_population.append(children[1])
        #Mutation
        new_population = mutate(new_population)
        population = new_population
        print("Finished generation", gen)
        sys.stdout.flush()
    final_fitness = calc_all_fitness(population, False)
    f = open("final_population", "w")
    f.write("types: " + str(N_TYPES) + "\nworld size: " + str(WORLD_X) + "\nupdates: " + str(UPDATES) + "\nmax pop: " + str(MAX_POP) + "\nrepro threshold: " + str(REPRO_THRESHOLD) + "\n")
    f.write("pop size: " + str(pop_size) + "\ngenerations: " + str(generations) + "\n")
    for i in range(len(final_fitness)):
        f.write(str(final_fitness[i]) + '  ' + str(population[i]) + '\n\n')
    f.close()


def track_avg_fitness(all_fitness, test_cases):
    f = open("avg_fitness.txt", "a")
    for test in test_cases:
        sum_fitness = 0
        for fitness in all_fitness:
            sum_fitness += fitness[test]
        avg = sum_fitness / len(all_fitness)
        f.write(test + ': ' + str(avg) + '\n')
    f.write('\n')
    f.close()


def track_max_fitness(all_fitness, test_cases):
    f = open("max_fitness.txt", "a")
    for test in test_cases:
        fitnesses = []
        for fitness in all_fitness:
            fitnesses.append(fitness[test])
        f.write(test + ': ' + str(max(fitnesses)) + '\n')
    f.write('\n')
    f.close()


if __name__ == '__main__':
    run()
