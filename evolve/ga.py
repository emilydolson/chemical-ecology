import subprocess
import pandas as pd
import random
import sys
import csv
import numpy as np
from klemm_eguiliuz import create_matrix

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
        mu = round(random.uniform(-0.5, 0.5), 3)
        sigma = round(random.uniform(0, 0.5), 3)
        clique_size = random.randint(1, N_TYPES-1)
        genome = [diffusion, seeding, clear, clique_linkage, mu, sigma, clique_size]
        pop.append(genome)
    return pop


def calc_all_fitness(population, niched=False):
    fitness_lst = []
    for genome in population:
        diffusion = genome[0]
        seeding = genome[1]
        clear = genome[2]
        interaction_matrix = create_matrix(num_nodes=N_TYPES, clique_size=genome[6], clique_linkage=genome[3], mu=genome[4], sigma=genome[5])
        interaction_matrix_file = 'interaction_matrix.dat'
        with open(interaction_matrix_file, 'w') as f:
            wr = csv.writer(f)
            wr.writerows(interaction_matrix)
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
        new_fitness = {
            # TODO use list slicing to subtract final pop from avg of first 5 pops
            'Biomass': biomasses[-1] - biomasses[0],
            'Growth_Rate': growth_rates[-1] - growth_rates[0],
            # Heredity always starts at 1 and goes down. Just take the final heredity value
            'Heredity': heredities[-1],
            # Scores
            'Biomass_Score': b_score[0],
            'Growth_Rate_Score': g_score[0],
            'Heredity_Score': h_score[0],
            'Invasion_Ability_Score': i_score[0],
            'Resiliance_Score': r_score[0]
        }
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
            sharing = sum([max(0, 1 - (np.linalg.norm(np.array(genome_x[0:6] + [b/10 for b in genome_x[6:]]) - np.array(genome_y[0:6] + [b/10 for b in genome_y[6:]]))/sigma)**(alpha)) for genome_y in population])
            new_fitness = {
                'Biomass': fitnesses['Biomass'] / sharing,
                'Growth_Rate': fitnesses['Growth_Rate'] / sharing,
                'Heredity': fitnesses['Heredity'] / sharing, 
                'Biomass_Score': b_score[0]/sharing,
                'Growth_Rate_Score': g_score[0]/sharing,
                'Heredity_Score': h_score[0]/sharing,
                'Invasion_Ability_Score': i_score[0]/sharing,
                'Resiliance_Score': r_score[0]/sharing
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
                if i < 4: #diffusion, seeding, clear, clique_linkage
                    param += random.gauss(0, 0.1)
                    if param < 0:
                        param = 0
                    if param > 1:
                        param = 1
                if i == 4: #mu
                    param += random.gauss(0, 0.1)
                    if param < -0.5:
                        param = -0.5
                    if param > 0.5:
                        param = 0.5
                if i == 5: #sigma
                    param += random.gauss(0, 0.1)
                    if param < 0:
                        param = 0
                    if param > 0.5:
                        param = 0.5
                if i >= 6: #clique_size, seed
                    param += random.choice([-1, 1])
                    if param < 1:
                        param = 1
                    if param > N_TYPES-1:
                        param = N_TYPES-1
            genome[i] = round(param, 3)
    return pop


def run():
    pop_size = 100
    generations = 100
    population = create_pop(pop_size)
    test_cases = ['Biomass', 'Growth_Rate', 'Heredity', 'Biomass_Score', 'Growth_Rate_Score', 'Heredity_Score', 'Invasion_Ability_Score', 'Resiliance_Score']
    for gen in range(generations):
        all_fitness = calc_all_fitness(population, True)
        parent_tuple = (None, population, all_fitness)
        parents = []
        for _ in range(pop_size//2):
            parent_tuple = lexicase_select(parent_tuple[1], parent_tuple[2], test_cases)
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


if __name__ == '__main__':
    run()
