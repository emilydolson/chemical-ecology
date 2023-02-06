
#https://github.com/Nathaniel-Rodriguez/graphgen

import random
from multiprocessing import Process, Queue
import time
import networkx as nx
import matplotlib.pyplot as plt
from graphgen.lfr_generators import weighted_directed_lfr_as_adj

queue = Queue()

def matrix_process(num_nodes, average_k, max_degree, mut, muw, beta, com_size_min, com_size_max, tau, tau2, overlapping_nodes, overlap_membership, pct_pos_in, pct_pos_out):
    matrix, communities = weighted_directed_lfr_as_adj(
                                            num_nodes=num_nodes, 
                                            average_k=average_k, 
                                            max_degree=max_degree, 
                                            mut=mut, 
                                            muw=muw, 
                                            beta=beta,
                                            com_size_min=com_size_min, 
                                            com_size_max=com_size_max,
                                            tau=tau,
                                            tau2=tau2,
                                            overlapping_nodes=overlapping_nodes,
                                            overlap_membership=overlap_membership,
                                            seed=0
                                            )
    queue.put((matrix, communities))


def create_matrix(num_nodes, average_k, max_degree, mut, muw, beta, com_size_min, com_size_max, tau, tau2, overlapping_nodes, overlap_membership, pct_pos_in, pct_pos_out):
    '''
    :param num_nodes: Number of nodes in the network (starts id at 0)
    :param average_k: average degree of the nodes
    :param max_degree: largest degree of the nodes
    :param mut: mixing parameter, fraction of bridges
    :param muw: weight mixing parameter
    :param beta: minus exponent for weight distribution
    :param com_size_min: smallest community size
    :param com_size_max: largest community size
    :param tau: minus exponent for degree sequence
    :param tau2: minus exponent for community size distribution
    :param overlapping_nodes: number of overlapping nodes
    :param overlap_membership: number of memberships of overlapping nodes
    '''

    random.seed(0)

    p = Process(target=matrix_process, name='Create Matrix', args=(num_nodes, average_k, max_degree, mut, muw, beta, com_size_min, com_size_max, tau, tau2, overlapping_nodes, overlap_membership, pct_pos_in, pct_pos_out,))
    p.start()
    time.sleep(3)
    if p.is_alive():
        p.terminate()
        p.join()
        return [[0 for _ in range(num_nodes)] for _ in range(num_nodes)]
    matrix, communities = queue.get()

    min_all = min([min(x) for x in matrix])
    max_all = max([max(x) for x in matrix])

    if len(communities) > 0:
        for row in range(len(matrix)):
            for col in range(len(matrix)):
                if matrix[row][col] != 0:
                    matrix[row][col] = (matrix[row][col] - min_all) / (max_all - min_all)
                    if sum([x in communities[row] for x in communities[col]]) == 0:
                        if random.random() > 1:
                            matrix[row][col] = -matrix[row][col]
                    else:
                        if random.random() > 1:
                            matrix[row][col] = -matrix[row][col]

    return matrix


def create_matrix_unformatted(genome):
    diffusion = genome[0]
    seeding = genome[1]
    clear = genome[2]
    average_k = genome[12]
    max_degree = genome[13]
    mut = genome[3]
    muw = genome[4]
    beta = genome[5]
    com_size_min = genome[10]
    com_size_max = genome[11]
    tau = genome[6]
    tau2 = genome[7]
    overlapping_nodes = genome[14]
    overlap_membership = genome[15]
    pct_pos_in = genome[8]
    pct_pos_out = genome[9]
    matrix = create_matrix(num_nodes=9, average_k=average_k, max_degree=max_degree, mut=mut, muw=muw, beta=beta, com_size_min=com_size_min, com_size_max=com_size_max, tau=tau, tau2=tau2, overlapping_nodes=overlapping_nodes, overlap_membership=overlap_membership, pct_pos_in=pct_pos_in, pct_pos_out=pct_pos_out)
    return matrix


def save_graph(matrix):
    G = nx.DiGraph(matrix)

    widths = {(u, v): f'{d["weight"]*10:.2f}' for u, v, d in G.edges(data=True)}
    edge_labels = {(u, v): f'{d["weight"]:.2f}' for u, v, d in G.edges(data=True)}
    nodelist = G.nodes()

    weights = nx.get_edge_attributes(G, 'weight')
    color_map = []
    for edge in G.edges():
        if weights[edge] < 0:
            color_map.append('red')
        else: 
            color_map.append('blue')  

    plt.figure(figsize=(12,8))
    pos = nx.shell_layout(G)
    nx.draw_networkx_nodes(G,pos,
                        nodelist=nodelist,
                        node_size=1500,
                        node_color='black')
    nx.draw_networkx_edges(G,pos,
                        edgelist = widths.keys(),
                        width=list(widths.values()),
                        edge_color=color_map,
                        alpha=0.6)
    nx.draw_networkx_labels(G, pos=pos,
                            labels=dict(zip(nodelist,nodelist)),
                            font_color='white')
    nx.draw_networkx_edge_labels(G, pos, edge_labels)

    plt.savefig('lfr.png')