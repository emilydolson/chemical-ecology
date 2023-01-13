import networkx as nx
import matplotlib.pyplot as plt


def with_pagerank(edges):
    weights = {}
    ranks = open('../pagerank.txt').read().split('\n')
    for line in ranks[22:]:
        if len(line) > 0:
            node, weight = line.split(' ')
            weight = round(float(weight), 3)
            weights[node] = weight

    graph = nx.DiGraph()
    graph.add_edges_from(edges)

    nx.set_node_attributes(graph, weights, name='weight')

    color_map = []
    for node in graph:
        if graph.out_degree(node) == 0:
            color_map.append('red')
        else:
            color_map.append('pink')

    labels = {n: n+'\n'+str(graph.nodes[n]['weight']) for n in graph.nodes}
    plt.figure(figsize=(25, 25))
    nx.draw(graph, with_labels=True, edge_color='gray', node_size=4200, node_color=color_map, labels=labels)
    plt.savefig('pagerank.png')


def visualize(edges, file_name):
    graph = nx.DiGraph()
    graph.add_edges_from(edges)

    print(file_name, ' nodes: ', graph.number_of_nodes())

    color_map = []
    for node in graph:
        if graph.out_degree(node) == 0:
            color_map.append('red')
        else:
            color_map.append('pink')

    plt.figure(figsize=(25, 25))
    nx.draw(graph, with_labels=True, edge_color='gray', node_size=4000, node_color=color_map)
    plt.savefig(f'{file_name}.png')


def main():
    graphs = {}
    communities = open('../community_fitness.txt').read().split('\n')
    curr_graph = ''
    for line in communities[22:]:
        if len(line) > 0:
            if line[0] != '*':
                graphs[curr_graph].append(line.split(' '))
            else:
                graph_name = line.split('***')[1]
                graphs[graph_name] = []
                curr_graph = graph_name

    #with_pagerank(graphs['Community Assembly'])

    for graph in graphs:
        print(graph, ' edges: ', len(graphs[graph]))
        visualize(graphs[graph], graph)


if __name__ == '__main__':
    main()