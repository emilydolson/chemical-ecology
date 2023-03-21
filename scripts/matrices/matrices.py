import csv
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
import networkx as nx


def write_matrix(interactions, output_name):
    with open(output_name, 'w') as f:
        wr = csv.writer(f)
        wr.writerows(interactions)


def visualize_network(matrix, output_name, title=''):
    #https://stackoverflow.com/questions/23994020/colorplot-that-distinguishes-between-positive-and-negative-values
    cmap = plt.get_cmap('PuOr')
    cmaplist = [cmap(i) for i in range(cmap.N)]
    cmap = cmap.from_list('Custom cmap', cmaplist, cmap.N)
    bounds = [-1, -0.75, -0.5, -0.25, -.0001, .0001, 0.25, 0.5, 0.75, 1]
    norm = BoundaryNorm(bounds, cmap.N)

    plt.imshow(matrix, cmap=cmap, interpolation='none', norm=norm)
    #plt.tick_params(which = 'both', size = 0, labelsize = 0)
    #plt.colorbar(ticks=[-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1], orientation='horizontal', location='bottom', fraction=0.052, pad=0.1)
    plt.colorbar(ticks=[-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1])
    plt.title(title, fontsize=15)
    plt.savefig(output_name)
    plt.close()


def visualize_graph(matrix, output_name):
    G = nx.DiGraph(matrix)

    widths = {(u, v): f'{d["weight"]*10:.2f}' for u, v, d in G.edges(data=True)}
    edge_labels = {(u, v): f'{d["weight"]:.2f}' for u, v, d in G.edges(data=True)}
    nodelist = G.nodes()

    weights = nx.get_edge_attributes(G, 'weight')
    color_map = []
    for edge in G.edges():
        if weights[edge] < 0:
            color_map.append('saddlebrown')
        elif weights[edge] > 0: 
            color_map.append('mediumblue')  

    plt.figure(figsize=(5,5))
    pos = nx.circular_layout(G)
    nx.draw(G, pos=pos, with_labels=True, node_size=1000, edge_color=color_map, node_color='lightgrey')

    plt.savefig(output_name)
    plt.close()