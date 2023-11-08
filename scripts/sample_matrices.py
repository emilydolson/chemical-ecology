from common import get_plots_path, get_schemes, get_matrix_function, get_scheme_bounds
from matplotlib.colors import BoundaryNorm
import matplotlib.pyplot as plt
from analysis_scheme import create_folder_structure
import random
import sys


def sample_plot(scheme):
    create_folder_structure(scheme)

    scheme_func = get_matrix_function(scheme)
    bounds = get_scheme_bounds(scheme)

    cmap = plt.get_cmap('PuOr')
    cmaplist = [cmap(i) for i in range(cmap.N)]
    cmap = cmap.from_list('Custom cmap', cmaplist, cmap.N)
    cm_bounds = [-1, -0.75, -0.5, -0.25, -.0001, .0001, 0.25, 0.5, 0.75, 1]
    norm = BoundaryNorm(cm_bounds, cmap.N)
    figure, axis = plt.subplots(5, 2, figsize=(12,12))

    row = 0
    col = 0
    for _ in range(10):
        matrix_params = [9]
        for i in range(len(bounds[0])):
            lower_bound = bounds[0][i]
            upper_bound = bounds[1][i]
            if bounds[2][i] == True:
                matrix_params.append(round(random.randint(lower_bound, upper_bound), 3))
            else:
                matrix_params.append(round(random.uniform(lower_bound, upper_bound), 3))
        matrix_params.append(0)
        matrix = scheme_func(*matrix_params)

        img = axis[row][col].imshow(matrix, cmap=cmap, interpolation='none', norm=norm)
        figure.colorbar(img, ticks=[-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1])
        axis[row][col].set_title(matrix_params)
  
        row += 1
        if row % 5 == 0:
            col += 1
            row = 0
    figure.tight_layout()
    plt.savefig(f'{get_plots_path()}{scheme}/samples.png')
    plt.close()


def main():
    if len(sys.argv) == 2:
        scheme = sys.argv[1]
        sample_plot(scheme)
    else:
        for scheme in get_schemes():
            sample_plot(scheme)


if __name__ == '__main__':
    main()