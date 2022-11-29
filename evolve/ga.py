import subprocess
import pandas as pd
import random
import sys
import csv
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
        clique_size = random.randint(1,8)
        clique_linkage = round(random.random(), 3)
        genome = [diffusion, seeding, clear, clique_size, clique_linkage]
        pop.append(genome)
    return pop


#todo remove print statements from chemical-ecology, remove output of debug_file and interaction_matrix.dat
def calc_all_fitness(population):
    fitness_lst = []
    for genome in population:
        diffusion = genome[0]
        seeding = genome[1]
        clear = genome[2]
        interaction_matrix = create_matrix(N_TYPES, genome[3], genome[4])
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
        new_fitness = {
            # TODO this is redundant -- should we be averaging first 5(?) runs and then subtracting?
            'Biomass': biomasses[-1] - biomasses[0],
            'Growth_Rate': growth_rates[-1] - growth_rates[0],
            # Heredity always starts at 1 and goes down. Just take the final heredity value
            'Heredity': heredities[-1]
        }
        fitness_lst.append(new_fitness)
    return fitness_lst


def lexicase_select(pop, all_fitness, test_cases):
    random.shuffle(test_cases)
    while len(test_cases) > 0:
        # Select which test case to use 
        first = test_cases[0]
        # Initialize best to a negative value so it will lose the first comparison
        best = -999
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
    # If we mutate, swap the world parameters and matrix parameters
    if (random.random() < .9):
        child1 = genome1[0:3] + genome2[3:]
        child2 = genome2[0:3] + genome1[3:]
        return (child1, child2)
    else:
        return(genome1, genome2)


# https://github.com/DEAP/deap/blob/master/deap/tools/mutation.py
def mutate(pop):
    for genome in pop:
        params = genome[0:3]
        for param in params:
            if random.random() < 1/len(params):
                param += random.gauss(0, .1)
                # Temp error checking
                if param < 0:
                    param = 0
                if param > 1:
                    param = 1
        for i, param in enumerate(params):
            genome[i] = param
    return pop


def run():
    pop_size = 200
    generations = 100
    population = create_pop(pop_size)
    test_cases = ['Biomass', 'Growth_Rate', 'Heredity']
    for gen in range(generations):
        all_fitness = calc_all_fitness(population)
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
    final_fitness = calc_all_fitness(population)
    f = open("final_population", "w")
    for i in range(len(final_fitness)):
        f.write(str(final_fitness[i]) + '  ' + str(population[i]) + '\n')
    f.close()


if __name__ == '__main__':
    run()