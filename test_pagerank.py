import subprocess
import sys
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt

from scripts.matrices.matrices import write_matrix
from scripts.matrices.klemm_lfr_graph import create_matrix

def read_pagerank():
    weights = {}
    ranks = open('page_rank.txt').read().split('\n')
    curr_graph = ''
    for line in ranks:
        if len(line) > 0:
            if line[0] != '*':
                node, weight = line.split(' ')
                weight = round(float(weight), 5)
                weights[curr_graph][node] = weight
            else:
                graph_name = line.split('***')[1]
                weights[graph_name] = {}
                curr_graph = graph_name
    return weights


def read_graphs():
    graphs = {}
    communities = open('community_fitness.txt').read().split('\n')
    curr_graph = ''
    for line in communities:
        if len(line) > 0:
            if line[0] != '*':
                graphs[curr_graph].append(line.split(' '))
            else:
                graph_name = line.split('***')[1]
                graphs[graph_name] = []
                curr_graph = graph_name
    return graphs


def main():
    diffusion = 0.0
    seeding = 0.25
    clear = 0.1
    matrix_params = [9, 3, 0.5, 2.5, 0.75, 0.1, 0.2, 0]
    interaction_matrix_file = 'interaction_matrix.dat'

    write_matrix(create_matrix(*matrix_params), interaction_matrix_file)
    chem_eco = subprocess.Popen(
        [(f'./analyze_communities '
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
        print("Error in analyze_communities, return code:", return_code)
        sys.stdout.flush()
        exit()

    soup_world_file = open('soup_world.txt').read().split('\n')
    soup_world = {}
    for line in soup_world_file:
        if len(line) > 0:
            node, proportion = line.split(' ')
            soup_world[node] = float(proportion)

    assembly_adj = read_graphs()['Community Assembly']
    G = nx.DiGraph()
    G.add_weighted_edges_from(assembly_adj)
    our_pageranks = read_pagerank()['Community Assembly']
    pagerank = nx.pagerank(G, weight='weight', alpha=1)

    print('Node\tTrue Pagerank\tOur Pagerank\tSoup World Proportion')
    for i in our_pageranks:
        if i in soup_world:
            soup_world_proportion = soup_world[i]
        else: 
            soup_world_proportion = 0
        if i in pagerank:
            true_pagerank = round(pagerank[i], 5)
        else:
            true_pagerank = 0
        print(f'{i}:\t{true_pagerank}\t\t{round(our_pageranks[i], 5)}\t\t{soup_world_proportion}')

    #plt.figure(figsize=(10, 10))
    #nx.draw(G, with_labels=True, edge_color='gray')
    #plt.savefig('assembly_graph.png')
    #plt.close()


if __name__ == '__main__':
    main()