import subprocess
import pandas as pd

N_TYPES = 10
WORLD_X = 10
WORLD_Y = 10
UPDATES = 1000
MAX_POP = 10000
SEEDING_PROB = 0.01
REPRO_THRESHOLD = 5
REPRO_DILUTION = 0.1
PROB_CLEAR = 0.1


#todo remove print statements from chemical-ecology, remove output of debug_file and interaction_matrix.dat
def fitness(diffusion, interaction_matrix):
    subprocess.run(['../chemical-ecology', f'-DIFFUSION {diffusion} -REPRO_THRESHOLD {REPRO_THRESHOLD} -MAX_POP {MAX_POP} -SEEDING_PROB {SEEDING_PROB} -PROB_CLEAR {PROB_CLEAR} -INTERACTION_SOURCE {interaction_matrix} -WORLD_X {WORLD_X} -WORLD_Y {WORLD_Y} -N_TYPES {N_TYPES}'], shell=True)
    df = pd.read_csv('a-eco_data.csv')
    #todo calculate fitness


if __name__ == '__main__':
    fitness(0.05, 'evolve/networks/proof_of_concept_interaction_matrix.dat')