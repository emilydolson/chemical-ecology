import subprocess
import sys
import networkx as nx
import matplotlib.pyplot as plt


def get_stats(edges, weights):
    graph = nx.DiGraph()
    graph.add_edges_from(edges)

    sink_nodes = [node for node in graph if graph.out_degree(node) == 0]
    sink_nodes_pageranks = [weights[node] for node in sink_nodes]

    print(f'\t{graph.number_of_nodes()} nodes, {len(edges)} edges')
    print(f'\t{len(sink_nodes)} sink nodes')
    print(f'\t{sum(weights.values())} pagerank sum')
    print(f'\t{sum(sink_nodes_pageranks)} sink nodes pagerank sum')
    print()


def visualize(edges, weights):
    graph = nx.DiGraph()
    graph.add_edges_from(edges)

    node_sizes = []
    color_map = []
    labels = {}
    for node in graph:
        if graph.out_degree(node) == 0:
            node_sizes.append(8000)
        else:
            node_sizes.append(4500)
        labels[node] = f'{node}\n{weights[node]}'
        color_map.append('lightgray')

    plt.figure(figsize=(15, 15))
    nx.draw(graph, with_labels=True, edge_color='gray', node_size=node_sizes, labels=labels, node_color=color_map)
    plt.savefig('graph.png')
    plt.close()


def read_pagerank():
    weights = {}
    ranks = open('custom_pagerank.txt').read().split('\n')
    for line in ranks:
        if len(line) > 0:
            node, weight = line.split(' ')
            weight = round(float(weight), 3)
            weights[node] = weight
    return weights


def read_graphs():
    graphs = []
    communities = open('custom_graph.txt').read().split('\n')
    for line in communities:
        if len(line) > 0:
            graphs.append(line.split(' '))
    return graphs


def main():
    chem_eco = subprocess.Popen([('./custom_graph')], shell=True)
    return_code = chem_eco.wait()
    if return_code != 0:
        print("Error in custom_graph, return code:", return_code)
        sys.stdout.flush()

    graph = read_graphs()
    weights = read_pagerank()

    visualize(graph, weights)
    get_stats(graph, weights)


if __name__ == '__main__':
    main()