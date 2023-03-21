import subprocess
import pandas as pd
import random
import sys
import csv
import numpy as np
import statistics

# Search a matrix over all parameters. Sbatch script is: sbatch search_params.sh class1test class1

def search_world_params(dat_file):
    interaction_matrix_file = f'../config/{dat_file}.dat'
    measures = ["Biomass", "Growth_Rate", "Heredity", "Invasion_Ability", "Resiliance"]
    fields = ["Diffusion", "Seeding", "Clear", "Variance"]
    # Create file and headers
    for m in measures:
        f_name = m + ".csv"
        with open(f_name, 'w') as csvf:
            csvwriter = csv.writer(csvf)
            csvwriter.writerow(fields)

    for diffusion in np.arange(0, 1.01, .1):
        for seeding in np.arange(0, .101, .01):
            for clear in np.arange(0, 1.01, .1):
                #Average over 5 runs
                dfs = []
                for _ in range(5):
                    chem_eco = subprocess.Popen(
                        [(f'../chemical-ecology '
                        f'-DIFFUSION {diffusion} '
                        f'-SEEDING_PROB {seeding} '
                        f'-PROB_CLEAR {clear} ' 
                        f'-INTERACTION_SOURCE {interaction_matrix_file} '
                        f'-REPRO_THRESHOLD {10000000} '
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
                    score = df.values[0][-1]
                    dfs.append(df)
                for i, m in enumerate(measures):
                    f_name = m + ".csv"
                    scores = []
                    for frame in dfs:
                        scores.append(df.values[0][i] - score)
                    variance = statistics.variance(scores)
                    with open(f_name, 'a') as csvf:
                        csvwriter = csv.writer(csvf)
                        row = [f"{diffusion}", f"{seeding}", f"{clear}", str(variance)]
                        csvwriter.writerow(row)
                # with open("results.txt", "a") as f:
                #     f.write(f"D{diffusion}, S{seeding}, C{clear}, {scores}\n")
        print(diffusion)


if __name__ == '__main__':
    search_world_params(sys.argv[1])