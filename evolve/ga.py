import subprocess
import pandas as pd
import random
import sys
import csv
import numpy as np
from lfr_graph import create_matrix

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
        mut = 0.3 #round(random.random(), 3)
        muw = 0.1 #round(random.random(), 3)
        beta = 1.5 #round(random.random(), 3)
        tau = 2 #round(random.random(), 3)
        tau2 = 1 #round(random.random(), 3)
        pct_pos_in = round(random.random(), 3)
        pct_pos_out = round(random.random(), 3)
        com_size_min = random.randint(3, N_TYPES-3)
        com_size_max = random.randint(com_size_min+1, N_TYPES-1)
        average_k = random.randint(3, N_TYPES-3)
        max_degree = random.randint(average_k+1, N_TYPES-2)
        overlapping_nodes = random.randint(com_size_min, com_size_max)
        overlap_membership = com_size_min
        genome = [diffusion, seeding, clear, mut, muw, beta, tau, tau2, pct_pos_in, pct_pos_out, com_size_min, com_size_max, average_k, max_degree, overlapping_nodes, overlap_membership]
        pop.append(genome)
    return pop


def euclidean_distance(genome_x, genome_y):
    return np.linalg.norm(np.array(genome_x[0:5] + [b/10 for b in genome_x[5:]]) - np.array(genome_y[0:5] + [b/10 for b in genome_y[5:]]))


def calc_all_fitness(population, niched=False):
    fitness_lst = []
    for genome in population:
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
        interaction_matrix = create_matrix(num_nodes=N_TYPES, average_k=average_k, max_degree=max_degree, mut=mut, muw=muw, beta=beta, com_size_min=com_size_min, com_size_max=com_size_max, tau=tau, tau2=tau2, overlapping_nodes=overlapping_nodes, overlap_membership=overlap_membership, pct_pos_in=pct_pos_in, pct_pos_out=pct_pos_out)
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
            sharing = sum([max(0, 1 - (euclidean_distance(genome_x, genome_y)/sigma)**(alpha)) for genome_y in population])
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
                if i < 5:
                    param += random.gauss(0, 0.1)
                    if param < 0:
                        param = 0
                    if param > 1:
                        param = 1
                elif i >= 5 and i <= 7:
                    param += random.gauss(0, 0.1)
                    if param < 1:
                        param = 1
                    if param > 2:
                        param = 2
                elif i == 8 or i == 9:
                    param += random.gauss(0, 0.1)
                    if param < 0:
                        param = 0
                    if param > 0.75:
                        param = 0.75
                elif i == 10:
                    param += random.choice([-1, 1])
                    if param < 3:
                        param = 3
                    if param > N_TYPES-3:
                        param = N_TYPES-3
                elif i == 11:
                    param += random.choice([-1, 1])
                    if param < genome[10]+1:
                        param = genome[10]+1
                    if param > N_TYPES-1:
                        param = N_TYPES-1
                elif i == 12:
                    param += random.choice([-1, 1])
                    if param < 3:
                        param = 3
                    if param > N_TYPES-3:
                        param = N_TYPES-3
                elif i == 13:
                    param += random.choice([-1, 1])
                    if param < genome[12]+1:
                        param = genome[12]+1
                    if param > N_TYPES-2:
                        param = N_TYPES-2
                elif i == 14:
                    param += random.choice([-1, 1])
                    if param < genome[10]:
                        param = genome[10]
                    if param > genome[11]:
                        param = genome[11]
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
