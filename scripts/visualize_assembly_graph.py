import networkx as nx
import matplotlib.pyplot as plt

graph = nx.DiGraph()
edges = []
communities = open('../community.txt').read().split('\n')
for line in communities[22:]:
    temp = line.split(', ')
    if len(temp) == 3:
        if temp[0] != temp[1]: #keep self-loops?
            edges.append(temp)
graph.add_weighted_edges_from(edges)

color_map = []
for node in graph:
    if graph.out_degree(node) == 0:
        color_map.append('red')
    else: 
        color_map.append('orange')  

plt.figure(figsize=(25, 25))
nx.draw(graph, with_labels=True, edge_color='gray', node_size=4000, node_color=color_map)
plt.savefig('community_graph.png')