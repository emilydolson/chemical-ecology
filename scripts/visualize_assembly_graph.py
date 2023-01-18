import networkx as nx
import matplotlib.pyplot as plt

def visualize(edges, file_name, weights={}):
    graph = nx.DiGraph()
    graph.add_edges_from(edges)

    print('nodes: ', graph.number_of_nodes())

    color_map = []
    for node in graph:
        if graph.out_degree(node) == 0:
            color_map.append('red')
        else:
            color_map.append('pink')

    plt.figure(figsize=(25, 25))
    if len(weights) > 0:
        nx.set_node_attributes(graph, weights, name='weight')
        labels = {n: n+'\n'+str(graph.nodes[n]['weight']) for n in graph.nodes}
        nx.draw(graph, with_labels=True, edge_color='gray', node_size=4200, node_color=color_map, labels=labels)
    else:
        nx.draw(graph, with_labels=True, edge_color='gray', node_size=4000, node_color=color_map)
    plt.savefig(f'{file_name}.png')


def main(with_pagerank=False):
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

    weights = {}
    if with_pagerank:
        ranks = open('../pagerank.txt').read().split('\n')
        curr_graph = ''
        for line in ranks[22:]:
            if len(line) > 0:
                if line[0] != '*':
                    node, weight = line.split(' ')
                    weight = round(float(weight), 3)
                    weights[curr_graph][node] = weight
                else:
                    graph_name = line.split('***')[1]
                    weights[graph_name] = {}
                    curr_graph = graph_name

    for graph in graphs:
        print(graph)
        print('edges: ', len(graphs[graph]))
        if with_pagerank:
            visualize(graphs[graph], graph, weights[graph])
        else:
            visualize(graphs[graph], graph)
        print()


if __name__ == '__main__':
    main()