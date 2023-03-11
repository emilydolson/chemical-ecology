from get_results import get_data
import pandas as pd


def main():
    df = pd.DataFrame()
    for scheme in ['klemm', 'klemm_random', 'motif', 'mangal', 'random']:
        df_generation, df_network, abiotic_params, param_names, analysis_param_names, network_param_names = get_data(scheme)


if __name__ == '__main__':
    main()