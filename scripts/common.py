from matrix_functions import random_matrix, erdos_renyi_uniform, scale_free, erdos_renyi_normal


schemes = {'random_matrix':random_matrix, 'erdos_renyi_uniform':erdos_renyi_uniform, 'scale_free':scale_free, 'erdos_renyi_normal':erdos_renyi_normal}
code_location = '/mnt/home/leithers/grant/chemical-ecology/'


def get_scheme_bounds(scheme_name):
    ''' Record and access the bounds of parameters/arguments for each matrix generation scheme
    Arguments:
        scheme_name (string): the name of the matrix scheme function
    Returns:
        a list consisting of a list of lower_bounds, a list of upper_bounds, and a list of bools
        list of bools marks if the param should be sampled as an int or float
    '''

    bounds = {
        'random_matrix': [[0], [1], [False]],
        'erdos_renyi_uniform': [[0], [1], [False]],
        'erdos_renyi_normal': [[0, -0.5, 0], [1, 0.5, 0.5], [False, False, False]],
        'scale_free': [[0.001, 0, 0, 0], [1, 1, 1, 1], [False, False, False, False]]
    }
    return bounds[scheme_name]


def get_plots_path():
    return code_location+'scripts/output/plots/'


def get_processed_data_path():
    return code_location+'scripts/output/data/'


def get_raw_data_path():
    return '/mnt/scratch/leithers/chemical-ecology/data/matrix_schemes/'


def get_common_param_columns():
    return ['num_communities', 'world_size', 'ntypes', 'diffusion', 'seeding', 'clear']


def get_common_columns():
    return ['scheme', 'replicate', 'score'] + get_common_param_columns()


def get_schemes():
    return list(schemes.keys())


def get_matrix_function(scheme_name):
    return schemes[scheme_name]