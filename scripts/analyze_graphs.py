import subprocess
import networkx as nx
import matplotlib.pyplot as plt

from matrices import write_matrix
from klemm_lfr_graph import create_matrix


def get_stats(graphs, weights, communities):
    print(f'{len(communities)} final communities')
    print()

    for graph_type, edges in graphs.items():
        graph = nx.DiGraph()
        graph.add_edges_from(edges)

        sink_nodes = [node for node in graph if graph.out_degree(node) == 0]
        sink_nodes_pageranks = [weights[graph_type][node] for node in sink_nodes]
        overlapping_sink_final = [node for node in sink_nodes if node in communities]

        print(f'{graph_type}') 
        print(f'\t{graph.number_of_nodes()} nodes, {len(edges)} edges')
        print(f'\t{len(sink_nodes)} sink nodes')
        print(f'\t{sum(weights[graph_type].values())} pagerank sum')
        print(f'\t{sum(sink_nodes_pageranks)} sink nodes pagerank sum')
        print(f'\t{len(overlapping_sink_final)} overlapping sink and final community nodes')
        print()


def visualize(edges, file_name, weights, communities):
    graph = nx.DiGraph()
    graph.add_edges_from(edges)

    color_map = []
    node_sizes = []
    labels = {}
    for node in graph:
        if graph.out_degree(node) == 0:
            node_sizes.append(8000)
        else:
            node_sizes.append(4500)
        if node in communities:
            color_map.append('hotpink')
            labels[node] = f'{node}\n{weights[node]}\n{communities[node]}'
        else:
            color_map.append('lightgray')
            labels[node] = f'{node}\n{weights[node]}'

    plt.figure(figsize=(25, 25))
    nx.draw(graph, with_labels=True, edge_color='gray', node_size=node_sizes, node_color=color_map, labels=labels)
    plt.savefig(f'{file_name}.png')
    plt.close()


def read_finalcommunities():
    communities = {}
    ranks = open('final_communities.txt').read().split('\n')
    for line in ranks:
        if len(line) > 0:
            node, proportion = line.split(' ')
            proportion = round(float(proportion), 3)
            communities[node] = proportion
    return communities


def read_pagerank():
    weights = {}
    ranks = open('page_rank.txt').read().split('\n')
    curr_graph = ''
    for line in ranks:
        if len(line) > 0:
            if line[0] != '*':
                node, weight = line.split(' ')
                weight = round(float(weight), 3)
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
    diffusion = 0.215
    seeding = 0.01
    clear = 0.5
    matrix_params = [9, 3, 0.92, 2.5, 0.94, 0.81, 0.82, 7]
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
        shell=True)#, 
        #stdout=subprocess.DEVNULL)
    return_code = chem_eco.wait()
    if return_code != 0:
        print("Error in analyze-communities, return code:", return_code)
        sys.stdout.flush()

    graphs = read_graphs()
    weights = read_pagerank()
    communities = read_finalcommunities()

    for graph in graphs:
        visualize(graphs[graph], graph, weights[graph], communities)

    get_stats(graphs, weights, communities)


if __name__ == '__main__':
    main()