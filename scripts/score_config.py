import csv
import json
import pickle
import subprocess
import sys

import pandas as pd
from scipy.stats import qmc

from common import get_code_location
sys.path.insert(0, f"{get_code_location()}third-party/graph-evolution")
from ga import run
from eval_functions import Evaluation


def sample_abiotic_params(num_samples):
    sampler = qmc.LatinHypercube(d=3)
    unscaled_sample = sampler.random(n=num_samples)
    sample = qmc.scale(unscaled_sample, [0,0,0], [1,1,0.5]).tolist()
    abiotic_params = [[round(s[i], 3) for i in range(len(s))] for s in sample]
    return abiotic_params


def write_matrix(interactions, output_name):
    with open(output_name, "w") as f:
        wr = csv.writer(f)
        wr.writerows(interactions)


def is_properly_evolved(config, eval_obj, org):
    eval_funcs = config["eval_funcs"]
    for eval_func in eval_funcs.keys():
        target_val = eval_funcs[eval_func]["target"]
        if abs(org.getProperty(eval_func, eval_obj) - target_val) > 0.1:
            return False
    return True


def main(exp_name, config_name, rep):
    config = json.load(open(f"{get_code_location()}scripts/configs/{exp_name}/{config_name}.json"))
    eval_obj = Evaluation(config)
    config_num = config_name.split("_")[0]
    num_graphs = int(config["popsize"])
    ntypes = int(config["network_size"])

    abiotic_params = sample_abiotic_params(num_graphs)
    final_pop, _ = run(config)

    results = []
    for i in range(num_graphs):
        org = final_pop[i]
        diffusion = abiotic_params[i][0]
        seeding = abiotic_params[i][1]
        clear = abiotic_params[i][2]

        if not is_properly_evolved(config, eval_obj, org):
            continue
        matrix = org.adjacencyMatrix

        write_matrix(matrix, "interaction_matrix.dat")
        chem_eco = subprocess.Popen(
            [(f"./chemical-ecology "
            f"-DIFFUSION {diffusion} "
            f"-SEEDING_PROB {seeding} "
            f"-PROB_CLEAR {clear} "
            f"-INTERACTION_SOURCE interaction_matrix.dat "
            f"-SEED {rep} "
            f"-MAX_POP {10000} "
            f"-WORLD_WIDTH 10 "
            f"-WORLD_HEIGHT 10 "
            f"-UPDATES {1000} "
            f"-N_TYPES {ntypes}")],
            shell=True, 
            stdout=subprocess.DEVNULL)
        return_code = chem_eco.wait()
        if return_code != 0:
            print("Error in a-eco, return code:", return_code)
            sys.stdout.flush()
            continue

        df = pd.read_csv("output/ranked_threshold_communities_scores.csv")
        score = df["logged_mult_score"][0]
        num_communities = len(df)

        result = {
            "experiment": exp_name,
            "config": config_name,
            "config_num": config_num,
            "replicate": rep,
            "ntypes": ntypes,
            "seeding": seeding,
            "diffusion": diffusion,
            "clear": clear,
            "graph": org,
            "num_communities": num_communities,
            "score": score
        }
        results.append(result)

    with open("results.pkl", "wb") as f:
        pickle.dump(results, f)


if __name__ == "__main__":
    if len(sys.argv) == 4:
        exp_name = sys.argv[1]
        config_name = sys.argv[2]
        rep = sys.argv[3]
        main(exp_name, config_name, rep)
    else:
        print("Please give the valid arguments (experiment name, config name, replicate).")
        exit()