import json
import os
import sys

from common import get_code_location, get_configs_path, get_raw_data_path


def experiment_config(save_dir, exp_dir, exp_name, eval_funcs, network_size):
    config = {
        "data_dir": exp_dir,
        "name": exp_name,
        "reps": 1,
        "save_data": 1,
        "plot_data": 0,
        "scheme": "NSGAII",
        "popsize": 100,
        "mutation_rate": 0.005,
        "mutation_odds": [1,2,1],
        "crossover_odds": [1,2,2],
        "crossover_rate": 0.6,
        "weight_range": [-1,1],
        "network_size": network_size,
        "num_generations": 500,
        "eval_funcs": eval_funcs
    }

    config_path = f"{save_dir}{exp_dir}/{exp_name}.json"
    with open(config_path, "w") as f:
        json.dump(config, f, indent=4)


def topology_experiment(exp_dir):
    configs_path = get_configs_path()
    if not os.path.exists(configs_path+exp_dir):
        os.makedirs(configs_path+exp_dir)

    eval_funcs = []
    for connectance in [0.25, 0.5, 0.75, 1]:
        eval_funcs.append(
            {
            "weak_components": {"target": 1},
            "connectance": {"target": connectance}
            }
        )

    config_names = []
    for network_size in [10, 25, 50]:
        for i,eval_func in enumerate(eval_funcs):
            exp_name = f"{i}_{network_size}"
            experiment_config(configs_path, exp_dir, exp_name, eval_func, network_size)
            config_names.append(exp_name)
    
    return config_names


#https://stackoverflow.com/questions/312443/how-do-i-split-a-list-into-equally-sized-chunks
def chunks(lst, n):
    for i in range(0, len(lst), n):
        yield lst[i:i + n]


def generate_scripts(exp_dir, config_names):
    code_location = get_code_location()
    configs_path = get_configs_path()
    raw_data_path = get_raw_data_path()

    config_chunks = chunks(config_names, 99)
    for i,chunk in enumerate(config_chunks):
        with open(f"{configs_path}{exp_dir}/submit_experiments{i}", "w") as f:
            for config_name in chunk:
                f.write(f"sbatch {code_location}scripts/score_config.sb {exp_dir} {config_name}\n")


if __name__ == "__main__":
    experiment_name = sys.argv[1]
    if experiment_name == "topology":
        config_names = topology_experiment(experiment_name)
    else:
        print("Invalid experiment name.")
        exit()
    generate_scripts(experiment_name, config_names)