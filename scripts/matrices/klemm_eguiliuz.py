#based on:
#https://www.frontiersin.org/articles/10.3389/fncom.2011.00011/full
#https://github.com/hallucigenia-sparsa/seqtime/blob/master/R/klemm.game.R

import random

def add_edges(interactions, node_degrees, i, j):
    interactions[i][j] = round(random.uniform(-1, 1), 3)
    interactions[j][i] = round(random.uniform(-1, 1), 3)
    node_degrees[i] += 2
    node_degrees[j] += 2
    return interactions, node_degrees


def create_matrix(num_nodes, clique_size, clique_linkage, seed=0):
    random.seed(seed)

    m = clique_size
    n = num_nodes
    u = clique_linkage #probability new edge will be connected to non-active nodes

    #create the fully connected initial network, and add its nodes to active nodes
    interactions = [[0 for _ in range(n)] for _ in range(n)]
    nodes = list(range(n))
    active_nodes = nodes[0:m]
    deactivated_nodes = []
    node_degrees = [0]*n

    for i in range(m):
        for j in range(i+1,m):
            interactions, node_degrees = add_edges(interactions, node_degrees, i, j)

    #iteratively connect remaining nodes to active nodes or a random node, then remove one active node
    for i in range(m, n):
        for j in active_nodes:
            if u > random.random() or len(deactivated_nodes) == 0:
                interactions, node_degrees = add_edges(interactions, node_degrees, i, j)
            else:
                connected = False
                no_infinite_loop = 0
                while (not connected) and no_infinite_loop < 100:
                    deactivated_degrees = sum(node_degrees[m:n])
                    node_j = random.choice(deactivated_nodes)
                    if deactivated_degrees != 0:
                        if node_degrees[node_j]/deactivated_degrees > random.random():
                            if interactions[i][node_j] == 0 and interactions[node_j][i] == 0:
                                interactions, node_degrees = add_edges(interactions, node_degrees, i, node_j)
                        connected = True
                    else:
                        interactions, node_degrees = add_edges(interactions, node_degrees, i, node_j)
                        connected = True
                    no_infinite_loop += 1
        #replace an active node with node i. Active nodes with lower degrees are more likely to be replaced
        active_nodes.append(i)
        chosen = False
        no_infinite_loop = 0
        while (not chosen) and no_infinite_loop < 100:
            node_j = random.choice(active_nodes)
            if node_degrees[node_j] != 0:
                p = (1/node_degrees[node_j]) / sum([1/k if k != 0 else 0 for k in node_degrees])
                if p > random.random():
                    active_nodes.remove(node_j)
                    deactivated_nodes.append(node_j)
                    chosen = True
            no_infinite_loop += 1

    return interactions


def adjust_connectance(interactions, connectance):
    num_nodes = len(interactions)
    curr_connectance = (sum([sum([1 for y in x if y != 0]) for x in interactions])-num_nodes) / (num_nodes**2 - num_nodes)
    #init_connectance = (sum([sum([1 for y in x if y != 0]) for x in interactions])) / (len(interactions)**2)
    
    if curr_connectance < connectance:
        edge_num_added = 0
        while(curr_connectance < connectance):
            xpos = random.randint(0, num_nodes-1)
            ypos = random.randint(0, num_nodes-1)
            if (xpos != ypos):
                if interactions[xpos][ypos] == 0:
                    edge_num_added += 1
                    interactions[xpos][ypos] = random.gauss(0, 0.25)
            curr_connectance = (sum([sum([1 for y in x if y != 0]) for x in interactions])-num_nodes) / (num_nodes**2 - num_nodes)
    if curr_connectance > connectance:
        edge_num_removed = 0
        while(curr_connectance > connectance):
            xpos = random.randint(0, num_nodes-1)
            ypos = random.randint(0, num_nodes-1)
            if (xpos != ypos):
                if interactions[xpos][ypos] != 0:
                    edge_num_removed += 1
                    interactions[xpos][ypos] = 0
            curr_connectance = (sum([sum([1 for y in x if y != 0]) for x in interactions])-num_nodes) / (num_nodes**2 - num_nodes)

    return interactions