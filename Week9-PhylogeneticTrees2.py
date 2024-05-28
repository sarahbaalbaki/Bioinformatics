# Sarah Baalbaki 
# Bioinformatics- Module 9 

# 9.1- The Neighbor Joining Algorithm: 
import sys
from typing import List, Dict, Iterable, Tuple
import time
import numpy as np 

def NeighborJoining(distance_matrix: List[List[int]]) -> Dict[tuple, int]:
    # n= number of rows in D
    n = len(distance_matrix)
    clusters = [i for i in range(n+1)]

    T= dict()

    # counter for the number of nodes in tree
    counter = n 

    while True: 
        if n==2: 
            n_2= len(T)
            #print(len(T))
            #print("base")
            # tree consisting of a single edge of length D1,2
            T[n_2+1,n_2]= distance_matrix[0][1]
            T[n_2,n_2+1]= distance_matrix[0][1]
            #print("added", n_2, n_2+1)
            #print("T", T)
            T = FixDict(T)
            return T 
        
        # neighbor-joining matrix constructed from the distance matrix D
        D_star= Get_Neighbor_Joining_Matrix(distance_matrix, n)

        # find elements i and j such that D*i,j is a minimum non-diagonal element of D*
        i, j = GetClosest(D_star)
        #print("i, j:", i, j)

        # Δ ← (TotalDistanceD(i) - TotalDistanceD(j)) /(n - 2)
        delta = (Total_Distance_D(distance_matrix, i)- Total_Distance_D(distance_matrix, j))/(n-2)

        # limbLengthi ← (1/2)(Di,j + Δ)
        limbLengthi = (distance_matrix[i][j]+ delta)/2
        # limbLengthj ← (1/2)(Di,j - Δ)
        limbLengthj = (distance_matrix[i][j]- delta)/2

        # print(distance_matrix, type(distance_matrix))

        # add a new row/column m to D so that Dk,m = Dm,k = (1/2)(Dk,i + Dk,j - Di,j) for any k
        # m is the index of the new node, set it = counter 
        m= counter
        #print("M", m)
        new_row = []
        for k in range(len(distance_matrix)): 
            if k!= i and k!=j: 
                new_val = (distance_matrix[k][i]+distance_matrix[k][j]- distance_matrix[i][j])/2
                new_row.append(new_val)
                distance_matrix[k].append(new_val)
            else:
                new_row.append(0)
                distance_matrix[k].append(0)
        new_row.append(0)
        distance_matrix.append(new_row)
                
        # D ← D with rows i and j removed
        # D ← D with columns i and j removed
        # print(distance_matrix)
        distance_matrix= Remove_Rows_Columns(distance_matrix, [i,j])

        # add two new limbs (connecting node m with leaves i and j) to the tree T
        #print("added:", (clusters[i], m), limbLengthi)
        T[(m, clusters[i])] = limbLengthi
        #print("added:", (clusters[j], m), limbLengthj)
        T[(m, clusters[j])] = limbLengthj

        # increment counter of nodes since we added a new node 
        counter+=1

        #print(m, clusters[i], clusters[j], i, j)
        #print(clusters)
        # remove clusters in order of bigger one first so we do not remove the wrong cluster 
        if i < j:
            del clusters[j]
            del clusters[i]
        else:
            del clusters[i]
            del clusters[j]
        #print(clusters)
        if m not in clusters: 
            clusters.append(m)
        #print("CLUSTERS", clusters)
        
        # decrement length of distance matrix 
        n-=1
        
        #print(distance_matrix)
    #T = FixDict(T)
    return T 

# function to fix the dictionary to output desired nodes only once 
def FixDict(tree): 
    tree_new={}
    visited_pairs = set()
    for key, value in tree.items():
        #print(key, visited_pairs)
        if (key[0], key[1]) not in visited_pairs: 
            tree_new[key]= value
            #tree_new[(key[1], key[0])]=value
            visited_pairs.add((key[1], key[0]))
    tree_sorted = SortTree(tree_new)
    return tree_sorted

# function to sort the tree based on keys 
def SortTree(tree_new):
    sorted_keys = sorted(tree_new.keys(), key=lambda x: (x[0], x[1]))
    sorted_tree = {}
    for key in sorted_keys:
        sorted_tree[key] = tree_new[key]
    return sorted_tree

# function to remove rows and cols of nodes we join to new one 
def Remove_Rows_Columns(DistMatrix, indices):
    DistMatrix = np.array(DistMatrix, dtype= float)
    # reove rows
    DistMatrix = np.delete(DistMatrix, indices, axis=0)
    # remove columns
    DistMatrix = np.delete(DistMatrix, indices, axis=1)
    return DistMatrix.tolist()
    
# function to make the neighbor joining matrix 
def Get_Neighbor_Joining_Matrix(distance_matrix, n): 
    D_star = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i != j:
                D_star[i][j] = (n - 2) * distance_matrix[i][j] - Total_Distance_D(distance_matrix, i) - Total_Distance_D(distance_matrix, j)
            elif i==j: 
                D_star[i][j] = 0
    return D_star

# function to get the total distance 
def Total_Distance_D(distance_matrix, node_index): 
    total_distance = 0
    for dist in distance_matrix[node_index]:
        total_distance += dist
    return total_distance

# get closest clusters from distance matrix 
def GetClosest(distance_matrix): 
    min= np.inf 
    c1= np.inf 
    c2 = np.inf 
    for i in range(len(distance_matrix)): 
        for j in range(len(distance_matrix[0])): 
            if i !=j:  # do not take diag since = 0, not min 
                if distance_matrix[i][j]<min: 
                    min = distance_matrix[i][j]
                    c1 = i
                    c2 = j 
    return c1, c2 

##############################################################
# 9.3- Small Parsimony: 
import sys
from typing import List, Dict, Iterable, Tuple, Set, Union
import itertools
import re
from collections import deque
import numpy as np 

def SmallParsimony(leaves_labels: List[str],
                   edge_list: List[Tuple[int, int]]) -> Tuple[int, List[str]]:
    #pass
    # nodes_num = len(nodes)
    # tag = [0 for i in range(len(nodes_num))]
    # for v in T: 
    #     if IsLeaf(v)== True: 
    #         tag[v]=1

    # get num leaves 
    num_leaves = len(leaves_labels)
    # get total num nodes 
    total_nodes = 2*num_leaves-1

    # dictioanry to store tree, with corresponding node and children
    T= {}
    for edge in edge_list:
        node1 = edge[0]
        node2 = edge[1]
        if node1 in T:
            T[node1].append(node2)
        else:
            T[node1] = [node2]

    # store parents of each node, use to traverse and backtrack
    parent_tracker={}
    for edge in edge_list:
        parent = edge[0]
        child = edge[1]
        parent_tracker[child] = parent
    
    # get root candidates based on parents of node
    root_candidates = set()
    for value in parent_tracker.values():
        if value not in parent_tracker:
            root_candidates.add(value)
    # get the root 
    root_candidates_list = list(root_candidates)
    root = root_candidates_list[0]

    # length of labels 
    len_leaves_labels = len(leaves_labels[0])

    # define our alphabet
    alphabet= ('A', 'C', 'T', 'G')
    
    # define scoring matrix, num nodes x leaf label x alphabet 
    scores_matrix= np.empty((total_nodes, len_leaves_labels, len(alphabet)), dtype=int)
    for i in range(total_nodes):
        for j in range(len_leaves_labels ):
            for k in range(len(alphabet)):
                # scores_matrix[i, j, k] = (total_nodes - 1) * len_leaves_labels
                scores_matrix[i, j, k] = 100000000000000
    
    # labels we get after backtracking
    obtained_labels = []
    for i in range(total_nodes):
        inner_list = [' '] * len_leaves_labels
        obtained_labels.append(inner_list)

    # fill leaves' labels and scoring matrix, init
    for i in range(len(leaves_labels)):
        label = leaves_labels[i]
        obtained_labels[i] = list(label)
        for j in range(len(label)):
            character = label[j]
            scores_matrix[i, j, alphabet.index(character)] = 0

    # get scoring matrix 
    scores_matrix= GetScoringMatrix(root, num_leaves, T, scores_matrix, len_leaves_labels)
    
    #print(scores_matrix)
    
    # get parsimony score, min of scores 
    min_scores = scores_matrix[root].min(axis=1)
    parsimony_score = sum(min_scores)

    # get labels by back propagating scoring matrix 
    obtained_labels = BackPropagate(scores_matrix, obtained_labels, root, parent_tracker, T, alphabet, len_leaves_labels, num_leaves)

    # print(obtained_labels)

    # get node labels as strings
    node_labels = []
    for label in obtained_labels: 
        label_joined = ''.join(label)
        node_labels.append(label_joined)

    # return parsimony score and labeled node strings     
    return (parsimony_score,node_labels)

# function to get scoring matrix using dynamic programming 
def GetScoringMatrix(root, num_leaves, T, scores_matrix, len_leaves_labels):
    # go over all leaf labels 
    for i in range(len_leaves_labels):
        # init stack with root and flag for visited node or not 
        stack = [(root, False)]
        # until stack is empty
        while stack:
            # get node and visited status 
            # remove them from stack 
            node, visited = stack.pop()
            # if node not visited 
            if not visited:
                # if curr node not a leaf
                if node >= num_leaves:
                    # get left and right children
                    left = T[node][0]
                    right = T[node][1]
                    # set current mode to visited True 
                    stack.append((node, True))
                    # append them to stack to be visited 
                    stack.append((right, False))
                    stack.append((left, False))
            # if node visited 
            else:
                # get left and right children 
                left = T[node][0]
                right = T[node][1]
                # get scores for all poss labels for current node
                for j in range(4):
                    # check if letters match or not, calc ai_j 
                    ai_j = np.ones(4)
                    ai_j[j] = 0
                    # add ai_j to scores 
                    scores_matrix[node, i, j] = min(scores_matrix[left, i, :] + ai_j) + min(scores_matrix[right, i, :] + ai_j)
    return scores_matrix

# function to back propagate from scores to get chars 
def BackPropagate(scores_matrix, obtained_labels, root, parent_tracker, T, alphabet, len_leaves_labels, num_leaves):
# def BackPropagate(sk, s, root, parent, tree, d, alphabet, l, num_leaves):
    for i in range(len_leaves_labels):
        # stack for dfs 
        to_be_visited = deque()
        # start from root
        to_be_visited.append(root)

        # do until all nodes vsiited 
        while to_be_visited:
            # remove current node from stack
            node = to_be_visited.pop()

            # if curr node not a leaf
            if node >= num_leaves: 
                # get scores for node 
                node_scores = scores_matrix[node, i, :]

                # if root, get min score overall 
                if node == root:
                    index_min_char= np.argmin(node_scores)
                    obtained_labels[node][i] = alphabet[index_min_char]
                else:
                    parent_node = parent_tracker[node]
                    # get char and char index, use for ai_j
                    char = obtained_labels[parent_node][i]
                    char_index = alphabet.index(char)
                    # check if letters match or not, calc ai_j 
                    for ind in range(4):
                        ai_j = np.ones(4, dtype=int)  
                        ai_j[char_index] = 0
                    #print(ai_j)
                    # add ai_j to scores 
                    node_scores+= ai_j
                    # get min 
                    # min_score = min(node_scores)
                    index_min_char = np.argmin(node_scores)
                    #print(index_min_char)
                    #print("node scores", node_scores)
                    obtained_labels[node][i] = alphabet[index_min_char]

                # node children added to stack to visit next 
                to_be_visited.append(T[node][0]) # left
                to_be_visited.append(T[node][1]) # right 
    return obtained_labels

##############################################################
# 9.3- Small Parsimony for Unrooted Tree: 
import sys
from typing import List, Dict, Iterable, Tuple, Set, Union
from collections import defaultdict
import queue
import itertools
import re
from collections import deque
import numpy as np 

def SmallParsimonyUnrootedTree(leaves_labels: List[str],
                               edge_list: List[Tuple[int, int]]) -> Tuple[int, List[str]]:
    #pass
    min_parsimony_score = float('inf')
    node_labels = []

    # sort the edges from large -> small 
    edge_list_new= []
    for edge in edge_list: 
        new= (edge[1], edge[0])
        edge_list_new.append(new)

    # print(edge_list_new)

    # get root as the max node value+1     
    root = edge_list_new[-1][0]
    #print(root)
    root+=1

    # add root tpo the last edge only (any edge will result in same min parsimony score)
    for i, root_edge in enumerate(edge_list_new):
        if i != len(edge_list_new) - 1:
            continue
        # add the root and root edges 
        # root = len(leaves_labels)
        temp_edge_list = edge_list_new[:i] + [(root, root_edge[0])] + [(root, root_edge[1])] + edge_list_new[i+1:]
        # print(temp_edge_list)
        # run small parsimony on the tree with root (new edge list)
        score, labels = SmallParsimony(leaves_labels, temp_edge_list, root)
        # if for all edges: if score < min_parsimony_score: 
        # we are only doing this for one edge, no need to compare, assign min_parsimony_score and obtained_labels
        min_parsimony_score = score
        # remove the label of the last edge which connects to root node that we added, not needed anymore since we convert back to unrooted 
        node_labels = labels[:-1]

    return min_parsimony_score, node_labels

def SmallParsimony(leaves_labels: List[str],
                   edge_list: List[Tuple[int, int]], root) -> Tuple[int, List[str]]:
    # get num leaves 
    num_leaves = len(leaves_labels)
    # get total num nodes 
    total_nodes = 2*num_leaves-1

    # dictioanry to store tree, with corresponding node and children
    T= {}
    for edge in edge_list:
        node1 = edge[0]
        node2 = edge[1]
        if node1 in T:
            T[node1].append(node2)
        else:
            T[node1] = [node2]
    
    # print(T)

    # store parents of each node, use to traverse and backtrack
    parent_tracker={}
    for edge in edge_list:
        parent = edge[0]
        child = edge[1]
        parent_tracker[child] = parent
    
    #print(parent_tracker)
    
    # get root candidates based on parents of node
    root_candidates = set()
    for value in parent_tracker.values():
        if value not in parent_tracker:
            root_candidates.add(value)
    # get the root 
    root_candidates_list = list(root_candidates)
    # print("root_candidates_list", root_candidates_list)
    # root = root_candidates_list[0]
    root=root

    # length of labels 
    len_leaves_labels = len(leaves_labels[0])

    # define our alphabet
    alphabet= ('A', 'C', 'T', 'G')
    
    # define scoring matrix, num nodes x leaf label x alphabet 
    scores_matrix= np.empty((total_nodes, len_leaves_labels, len(alphabet)), dtype=int)
    for i in range(total_nodes):
        for j in range(len_leaves_labels ):
            for k in range(len(alphabet)):
                # scores_matrix[i, j, k] = (total_nodes - 1) * len_leaves_labels
                scores_matrix[i, j, k] = 100000000000000
    
    # labels we get after backtracking
    obtained_labels = []
    for i in range(total_nodes):
        inner_list = [' '] * len_leaves_labels
        obtained_labels.append(inner_list)

    # fill leaves' labels and scoring matrix, init
    for i in range(len(leaves_labels)):
        label = leaves_labels[i]
        obtained_labels[i] = list(label)
        for j in range(len(label)):
            character = label[j]
            scores_matrix[i, j, alphabet.index(character)] = 0

    # get scoring matrix 
    scores_matrix= GetScoringMatrix(root, num_leaves, T, scores_matrix, len_leaves_labels)
    
    #print(scores_matrix)
    
    # get parsimony score, min of scores 
    min_scores = scores_matrix[root].min(axis=1)
    #print(scores_matrix)
    parsimony_score = sum(min_scores)

    # get labels by back propagating scoring matrix 
    obtained_labels = BackPropagate(scores_matrix, obtained_labels, root, parent_tracker, T, alphabet, len_leaves_labels, num_leaves)

    # print(obtained_labels)

    # get node labels as strings
    node_labels = []
    for label in obtained_labels: 
        label_joined = ''.join(label)
        node_labels.append(label_joined)

    # return parsimony score and labeled node strings     
    return (parsimony_score,node_labels)

# function to get scoring matrix using dynamic programming 
def GetScoringMatrix(root, num_leaves, T, scores_matrix, len_leaves_labels):
    # go over all leaf labels 
    for i in range(len_leaves_labels):
        # init stack with root and flag for visited node or not 
        stack = [(root, False)]
        # until stack is empty
        while stack:
            # get node and visited status 
            # remove them from stack 
            node, visited = stack.pop()
            # if node not visited 
            if not visited:
                # if curr node not a leaf
                if node >= num_leaves:
                    # get left and right children
                    left = T[node][0]
                    right = T[node][1]
                    # set current mode to visited True 
                    stack.append((node, True))
                    # append them to stack to be visited 
                    stack.append((right, False))
                    stack.append((left, False))
            # if node visited 
            else:
                # get left and right children 
                left = T[node][0]
                right = T[node][1]
                # get scores for all poss labels for current node
                for j in range(4):
                    # check if letters match or not, calc ai_j 
                    ai_j = np.ones(4)
                    ai_j[j] = 0
                    # add ai_j to scores 
                    scores_matrix[node, i, j] = min(scores_matrix[left, i, :] + ai_j) + min(scores_matrix[right, i, :] + ai_j)
    return scores_matrix

# function to back propagate from scores to get chars 
def BackPropagate(scores_matrix, obtained_labels, root, parent_tracker, T, alphabet, len_leaves_labels, num_leaves):
# def BackPropagate(sk, s, root, parent, tree, d, alphabet, l, num_leaves):
    for i in range(len_leaves_labels):
        # stack for dfs 
        to_be_visited = deque()
        # start from root
        to_be_visited.append(root)

        # do until all nodes vsiited 
        while to_be_visited:
            # remove current node from stack
            node = to_be_visited.pop()

            # if curr node not a leaf
            if node >= num_leaves: 
                # get scores for node 
                node_scores = scores_matrix[node, i, :]

                # if root, get min score overall 
                if node == root:
                    index_min_char= np.argmin(node_scores)
                    obtained_labels[node][i] = alphabet[index_min_char]
                else:
                    parent_node = parent_tracker[node]
                    # get char and char index, use for ai_j
                    char = obtained_labels[parent_node][i]
                    char_index = alphabet.index(char)
                    # check if letters match or not, calc ai_j 
                    for ind in range(4):
                        ai_j = np.ones(4, dtype=int)  
                        ai_j[char_index] = 0
                    #print(ai_j)
                    # add ai_j to scores 
                    node_scores+= ai_j
                    # get min 
                    # min_score = min(node_scores)
                    index_min_char = np.argmin(node_scores)
                    #print(index_min_char)
                    #print("node scores", node_scores)
                    obtained_labels[node][i] = alphabet[index_min_char]

                # node children added to stack to visit next 
                to_be_visited.append(T[node][0]) # left
                to_be_visited.append(T[node][1]) # right 
    return obtained_labels

