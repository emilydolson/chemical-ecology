from matrix_functions import random_matrix


def get_plots_path():
    return '/mnt/home/leithers/grant/chemical-ecology/scripts/output/plots/'


def get_processed_data_path():
    return '/mnt/home/leithers/grant/chemical-ecology/scripts/output/data/'

   
def get_raw_data_path():
    return '/mnt/scratch/leithers/chemical-ecology/data/matrix_schemes/'


def get_common_columns():
    return ['scheme', 'replicate', 'score', 'ntypes', 'diffusion', 'seeding', 'clear']


def get_matrix_function(scheme_name):
    schemes = {'random_matrix':random_matrix}
    return schemes[scheme_name]