# Sarah Baalbaki 
# Bioinformatics- Module 4 Code 

# 4.5- DPChange: 
import math  
def DPChange(money, Coins):
    minNumCoins = [float('inf')] * (money + 1)  
    minNumCoins[0] = 0  
    for m in range(1, money+1):  
        for coin_i in Coins:  
            if m >= coin_i:  
                if minNumCoins[m - coin_i] + 1 < minNumCoins[m]:
                    minNumCoins[m] = minNumCoins[m - coin_i] + 1  
    return minNumCoins[money]  

money = 7 
coins = [1,5]
result = DPChange(money, coins)
print(result) 

###############################
# 4.6- LongestPathLength: 
import sys
from typing import List, Dict, Iterable, Tuple
import numpy as np 

def LongestPathLength(n: int, m: int,
                      Down: List[List[int]], Right: List[List[int]]) -> int:
    s = np.zeros((n+1, m+1)) 
    for i in range(1, n+1):
        s[i][0] = s[i-1][0] + Down[i-1][0]
    for j in range(1, m+1):
        s[0][j] = s[0][j-1] + Right[0][j-1]
    for i in range(1, n+1):
        for j in range(1, m+1):
            term1 = s[i-1][j] + Down[i-1][j]  
            term2 = s[i][j-1] + Right[i][j-1]  
            s[i][j] = max(term1, term2)
    return int(s[n][m]) 

###############################
# 4.7- OutputLCS: 
def LCSBackTrack(v, w):
    s = np.zeros((len(v)+1, len(w)+1))
    backtrack = np.zeros((len(v)+1, len(w)+1))
    for i in range(len(v) + 1):
        s[i][0] = 0
    for j in range(len(w) + 1):
        s[0][j] = 0
    for i in range(1, len(v) + 1):
        for j in range(1, len(w) + 1):
            match = 0
            if v[i - 1] == w[j - 1]:
                match = 1
            s[i][j] = max(s[i-1][j], s[i][j-1], s[i-1][j-1] + match)
            if s[i][j] == s[i-1][j]:
                backtrack[i][j] = 1  # down
            elif s[i][j] == s[i][j-1]:
                backtrack[i][j] = 2  # right
            elif s[i][j] == s[i-1][j-1]:
                backtrack[i][j] = 3  # diagonal
    return backtrack

def OutputLCS(backtrack, v, i, j):
    if i == 0 or j == 0:
        return ""
    if backtrack[i][j] == 1:  # down
        return OutputLCS(backtrack, v, i - 1, j)
    elif backtrack[i][j] == 2:  # right
        return OutputLCS(backtrack, v, i, j - 1)
    else:
        return OutputLCS(backtrack, v, i - 1, j - 1) + v[i - 1]

def LCS(v, w):
    backtrack = LCSBackTrack(v, w)
    s = OutputLCS(backtrack, v, len(v), len(w))
    return s

v = "GACT"
w = "ATG"
print(LCS(v, w))

###############################
# 4.7- LongestPath: 
def LongestPath(s: int, t: int, E: Dict[int, List[Tuple[int, int]]]) -> Tuple[int, List[int]]:
    # make graph without weights 
    graph = make_graph(E)
    # get top order of nodes 
    top_order= topological_ordering(graph)
    #print(top_order)
    # make 2 dictionaries, one for distances (from source) and one for predecessors of the nodes 
    dist = {node: float('-inf') for node in top_order}
    pred = {node: None for node in top_order}
    # dist of source is set to 0 to start the path 
    dist[s] = 0

    # for all nodes in top order 
    for curr_node in top_order[:-1]:
        # print("node", node)
        # print(E.get(node, []))
        # go over all next poss nodes 
        for s_b, weight in E.get(curr_node, []):
            # if dist to next node is less than the current plus the weight, we found a longer path
            if dist[s_b] < dist[curr_node] + weight:
                # update distance 
                dist[s_b] = dist[curr_node] + weight
                # update parent of s_b
                pred[s_b] = curr_node

    # get longest path by going through parents of s_b
    longest_path = []
    # start at end node and go back 
    node = t
    # for all nodes until we start, get node and its parent
    while node is not None:
        longest_path.append(node)
        node = pred[node]
    # reverse after starting with t  
    longest_path.reverse()

    longest_dist = dist[t]

    return longest_dist, longest_path

def topological_ordering(graph):
    # array to store order 
    topological_order = []
    # store nodes w no incoming edges 
    no_incoming_edges= set()
    # get all nodes 
    nodes= get_nodes(graph)
    
    # for all nodes, 
    for n in nodes: 
        state= False
        # if node has no incoming edges, add it to no_incoming_edges 
        for key in graph: 
            entry= graph[key]
            if n in entry: 
                state = True 
        if state == False: 
            no_incoming_edges.add(n)
    
    # for all poss start nodes, go through graph to get top order 
    while no_incoming_edges:
        # remove node from candidates to start with
        node_1 = no_incoming_edges.pop() 
        # add node to top order 
        topological_order.append(node_1)

        # go over all neighbors of node, remove neighbor from  node 1 
        for neighbor in list(graph.get(node_1, [])):
            graph[node_1].remove(neighbor)
            # for remining nodes in graph, if neighbor doesnt have anymore incoming edges, add it to array
            for nodes in graph: 
                if neighbor not in graph[nodes]:
                    no_incoming_edges.add(neighbor)

    # return reverse of top ordering to get correct order 
    top_order_final= topological_order
    return topological_order

def get_nodes(E): 
    # func to get all nodes even if they are not keys in the graph (dest nodes, incoming edges, no outgoing)
    set_nodes= set()
    for key in E: 
        set_nodes.add(key)
        for entries in E[key]: 
            set_nodes.add(entries)
    return list(set_nodes)

def make_graph(graph):
    # make new graph with only nodes and no weights 
    new_graph = {}
    for node, neighbors in graph.items():
        new_graph[node] = [neighbor[0] for neighbor in neighbors]
    return new_graph