import subprocess
import pandas as pd
import random
import sys
import csv
import numpy as np

def search_world_params(dat_file):
    interaction_matrix_file = f'../config/{dat_file}.dat'
    for diffusion in np.arange(0, 1, .1):
        for seeding in np.arange(0, 1, .1):
            for clear in np.arange(0, 1, .1):
                chem_eco = subprocess.Popen(
                    [(f'../chemical-ecology '
                    f'-DIFFUSION {diffusion} '
                    f'-SEEDING_PROB {seeding} '
                    f'-PROB_CLEAR {clear} ' 
                    f'-INTERACTION_SOURCE {interaction_matrix_file} '
                    f'-REPRO_THRESHOLD {100000000} '
                    f'-MAX_POP {10000} '
                    f'-WORLD_X {10} '
                    f'-WORLD_Y {10} '
                    f'-UPDATES {1000} '
                    f'-N_TYPES {9}')],
                    shell=True, 
                    stdout=subprocess.DEVNULL)
                return_code = chem_eco.wait()
                if return_code != 0:
                    print("Error in a-eco, return code:", return_code)
                    sys.stdout.flush()
                df = pd.read_csv('scores.csv')
                #df.next()
                scores = df.values
                with open("results.txt", "a") as f:
                    f.write(f"D{diffusion}, S{seeding}, C{clear}, {scores}\n")
        print(diffusion)


def search_world_params_strict():
    interaction_matrix_file = '../config/proof_of_concept_interaction_matrix.dat'
    for diffusion in np.arange(.3, .4, .01):
        for seeding in np.arange(.3, .4, .01):
            for clear in np.arange(.3, .4, .01):
                chem_eco = subprocess.Popen(
                    [(f'../chemical-ecology '
                    f'-DIFFUSION {diffusion} '
                    f'-SEEDING_PROB {seeding} '
                    f'-PROB_CLEAR {clear} ' 
                    f'-INTERACTION_SOURCE {interaction_matrix_file} '
                    f'-REPRO_THRESHOLD {100000000} '
                    f'-MAX_POP {10000} '
                    f'-WORLD_X {10} '
                    f'-WORLD_Y {10} '
                    f'-UPDATES {1000} '
                    f'-N_TYPES {9}')],
                    shell=True, 
                    stdout=subprocess.DEVNULL)
                return_code = chem_eco.wait()
                if return_code != 0:
                    print("Error in a-eco, return code:", return_code)
                    sys.stdout.flush()
                df = pd.read_csv('a-eco_data.csv')
                growth_rates = df['mean_Growth_Rate'].values
                gr = growth_rates[-1]
                if gr > 80000:
                    hits = replicate(diffusion, seeding, clear)
                    outfile = open("outfile.txt", "a")
                    outfile.write("D: " + str(diffusion) + "S: " + str(seeding) + "C: " + str(clear) + "hits: " + str(hits))
                    outfile.close()
        status_file = open("status.txt", "a")
        status_file.write(str(diffusion))
        status_file.close()

search_world_params(sys.argv[1])