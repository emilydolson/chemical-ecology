import subprocess
import sys
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt

from scripts.matrices.matrices import write_matrix
from scripts.matrices.klemm_lfr_graph import create_matrix

def main():
    diffusion = 0.25
    seeding = 0.25
    clear = 0.1
    matrix_params = [3, 2, 0.5, 2.5, 0.75, 0.1, 0.2, 0]
    interaction_matrix_file = 'interaction_matrix.dat'

    write_matrix(create_matrix(*matrix_params), interaction_matrix_file)
    chem_eco = subprocess.Popen(
        [(f'./chemical-ecology '
        f'-DIFFUSION {diffusion} '
        f'-SEEDING_PROB {seeding} '
        f'-PROB_CLEAR {clear} '
        f'-INTERACTION_SOURCE {interaction_matrix_file} '
        f'-REPRO_THRESHOLD {100000000} '
        f'-MAX_POP {10000} '
        f'-WORLD_X {10} '
        f'-WORLD_Y {10} '
        f'-UPDATES {1000} '
        f'-N_TYPES {matrix_params[0]}')],
        shell=True)
    return_code = chem_eco.wait()
    if return_code != 0:
        print("Error in chemical-ecology, return code:", return_code)
        sys.stdout.flush()
        exit()

    assembly_graph_file = open('assembly_graph.txt').read()
    assembly_adj = []
    our_pageranks = []
    i = 0
    for line in assembly_graph_file.split('\n'):
        if len(line) == 0:
            i += 1
        if i == 0:
            row = [float(x) for x in line.split(',')[:-1]]
            assembly_adj.append(row)
        else:
            our_pageranks = [float(x) for x in line.split()]

    row_sums_to_one = [(sum(x) > 0.99 and sum(x) < 1.01) or sum(x) == 0 for x in assembly_adj]
    if row_sums_to_one.count(False) > 0:
        print(f'Invalid Assembly Nodes: {[i for i in range(len(assembly_adj)) if not row_sums_to_one[i]]}\n')

    G = nx.DiGraph(np.array(assembly_adj))
    pagerank = nx.pagerank(G, weight='weight', alpha=1)
    for i in range(len(our_pageranks)):
        print(f'{i}:\t{round(pagerank[i], 5)}\t\t{round(our_pageranks[i], 5)}')

    plt.figure(figsize=(10, 10))
    nx.draw(G, with_labels=True, edge_color='gray')
    plt.savefig('assembly_graph.png')
    plt.close()


if __name__ == '__main__':
    main()