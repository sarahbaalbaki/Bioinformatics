# Sarah Baalbaki 
# Module 7 - Bioinformatics 

# 7.4- GreedySorting: 
def GreedySorting(P: List[int]) -> List[List[int]]:
    sorted_steps=[]
    approxReversalDistance = []
    for k in range(0,len(P)):
        if abs(P[k]) != k + 1:
            # perform k-sorting 
            P = KSorted(P, k)
            # print(P)
            approxReversalDistance.append(list(P))
        if P[k] == -(k + 1):
            # fix sign 
            P[k] = abs(P[k]) 
            approxReversalDistance.append(list(P))  
    return approxReversalDistance

def KSorted(P, k):
    # segment up to k sorted 
    segment = P[:k]
    P_new = list(segment) 
    # find new index to know where to cut segment to be reversed at, pos of k+1 val 
    k1_index = k
    # loop over remaining part (unsorted)
    for i in range(k + 1, len(P)):
        # if element at index is equal to k+1, assign k+1 index to it 
        if abs(P[i]) == k + 1:
            k1_index = i
    # reverse seg from k to k1_index
    P_new.extend(reversed([-val for val in P[k:k1_index + 1]]))
    # add remaining portion of permutation 
    P_new.extend(P[k1_index + 1:])
    return P_new

l = add_commas_between_numbers("-3 +4 +1 +5 -2")
GreedySorting(l)
ans= GreedySorting(add_commas_between_numbers("+2 +6 -8 -17 +7 -14 +18 +3 -5 -16 -11 -19 -4 +10 +13 -20 -1 +9 -12 +15"))
print(len(ans))
for a in ans: 
    print(a, end= "\n")

############################################
# 7.5- Breakpoints
def BreakpointCount(P: List[int]) -> int:
    # pass
    breakpoints = 0
    n = len(P)
    # breakpoints for beg and end
    if P[0] != 1:
        breakpoints += 1
    if P[-1] != n:
        breakpoints += 1
    # loop in permutation and get breakpoints
    for i in range(1, n):
        if P[i] - P[i-1] != 1:
            breakpoints += 1
    return breakpoints

def add_commas_between_numbers(numbers_str):
    number_strings = numbers_str.split()
    integers_list = [int(num_str) for num_str in number_strings]
    return integers_list

l = add_commas_between_numbers("+3 +4 +5 -12 -8 -7 -6 +1 +2 +10 +9 -11 +13 +14")
BreakpointCount(l)

############################################
# 7.9- TwoBreak Distance: 
def TwoBreakDistance(P: List[List[int]], Q: List[List[int]]) -> int:
    # pass
    # get edges 
    p_red_edges = ColoredEdges(P)
    q_blue_edges = ColoredEdges(Q)
    colored_edges = p_red_edges + q_blue_edges
    
    # adj list to store graph 
    adj_list = defaultdict(list)
    
    # for each edge (u, v), 
    # add v to adj list of u 
    # add u to adj list of v 
    for u, v in colored_edges:
        adj_list[u].append(v)
        adj_list[v].append(u)

    # let p be 0, red and q be 1,blue
    cycles_found = []
    visited_nodes = set()
    # BFS from each node in graph 
    for node in adj_list:
        if node not in visited_nodes:
            # update visited nodes from BFS nodes 
            bfs_res= (BFS(node, 'p', adj_list))
            visited_nodes.update(bfs_res)
            # add visited BFS nodes as cycle 
            cycles_found.append(list(visited_nodes))
    
    # get 2 break distance 
    num_colored = len(colored_edges)
    two_break_dist = num_colored / 2 - len(cycles_found) 
    return int(two_break_dist)

# BFS function from given node, chnaging colors accordingly 
def BFS(node, seq, adj_list):
        # queue for first node and seq
        queue = [(node, seq)] 
        # keep track of visited nodes
        visited = set()  
        # first node visisted 
        visited.add(node)  
        while queue:
            # remove node 
            node, seq = queue.pop(0)  
            for neighbor in adj_list[node]:
                if neighbor not in visited:
                    # change seq 
                    if seq == 'p': #p
                        next_seq = 'q'  # change to q
                    else: # q (=1)
                        next_seq = 'p'  # change to p
                    queue.append((neighbor, next_seq))  
                    # update visited node 
                    visited.add(neighbor)  
        return visited
    
def ColoredEdges(P: List[List[int]]) -> List[Tuple[int, int]]:
    Edges=[]
    for Chromosome in P: 
        Nodes = ChromosomeToCycle(Chromosome)
        #print(Nodes)
        for j in range (len(Chromosome)):
            edge= (Nodes[(2*j)-1], Nodes[(2*j)])
            Edges.append(edge)
    return Edges 
    
def ChromosomeToCycle(Chromosome: List[int]) -> List[int]:
    #pass
    nodes=[]
    for j in range(len(Chromosome)): 
        i = Chromosome[j]
        if i>0: 
            #print(2*i-1)
            nodes.append((2*i-1))
            #print(2*i)
            nodes.append(2*i)
        else: 
            #print(-2*i)
            nodes.append(-(2*i))
            
            nodes.append(-(2*i+1))
    return nodes 

############################################
# 7.9- TwoBreakSorting: 
genome_cycles_t= List[List[int]]
genome_graph_t= List[Tuple[int, int]]

def TwoBreakSorting(P: genome_cycles_t, Q: genome_cycles_t) -> List[genome_cycles_t]:
    #pass
    # get colored edges of P and Q 
    p_red_edges = ColoredEdges(P)
    q_blue_edges = ColoredEdges(Q)
    
    # store the results 
    results = [] 
    results.append(P)

    new_p = P 
    # make sure 2 break dist > 0 to continue modifying and sorting 
    while TwoBreakDistance(new_p,Q) > 0:
        #print(path)
        #print("P:", P)
        #print("ColoredEdgesP", ColoredEdges(P))
        p_new_red_edges= ColoredEdges(new_p)
        # get the cycles from the edges 
        cycles_found = GetCycles(p_new_red_edges, q_blue_edges)
                        
        #print("cycles found", cycles_found)
        # go over all cycles 
        for c in cycles_found:
            if len(c) >= 4:
                # make sure we have more than 4 nodes in cycle to perform 2-break 
                i0, i1, i2, i3= c[0],c[1],c[2],c[3]
                new_p = TwoBreakOnGenome(new_p,i0, i1, i2, i3)
                results.append(new_p)
                break          
    return results

def GetCycles(p_new_red_edges, q_blue_edges): 
        cycles_found=[]
        total_edges = len(p_new_red_edges) + len(q_blue_edges)
        adjacency_matrix = [[0, 0] for _ in range(total_edges+1)] 
        nodes_visited = set() #on stepik, this timed out (for test case 6) so I used a boolean array 

        # p is in index 0 and q is in index 1 
        for edge in p_new_red_edges:
            node1 = edge[0]
            node2 = edge[1]
            adjacency_matrix[node1][0] = node2
            adjacency_matrix[node2][0] = node1

        for edge in q_blue_edges:
            node1 = edge[0]
            node2 = edge[1]
            adjacency_matrix[node1][1] = node2
            adjacency_matrix[node2][1] = node1
    
        #print(adjacency_matrix)
        # go over all nodes 
        for node in range(1, total_edges+1):
            # if we did not prev visit the node, we visit it 
            if node not in nodes_visited:
                nodes_visited.add(node)
                # start from this node 
                start = node
                # keep track of cycle c 
                c = [start]
                # start with the p sequence, p:0:red
                seq = 0
                while len(nodes_visited) < total_edges:
                    node = adjacency_matrix[node][seq]
                    if node == start: # cycle 
                        cycles_found.append(c)
                        break
                    # if no cycle, we add node to cycle so far and continue 
                    c.append(node)
                    nodes_visited.add(node)
                    # change to second sequence
                    if seq==0: # p 
                        seq=1 # q 
                    else: 
                        seq=0 # p 
        return cycles_found 

def TwoBreakOnGenomeGraph(GenomeGraph: genome_graph_t,
                          i1: int, i2: int, i3: int, i4: int) -> genome_graph_t:
    #pass
    # remove edges (i1, i2) & (i3, i4) from graph
    if (i1, i2) in GenomeGraph: 
        GenomeGraph.remove((i1, i2))
    else: 
        GenomeGraph.remove((i2, i1))
    
    if (i3, i4) in GenomeGraph: 
        GenomeGraph.remove((i3, i4))
    else: 
        GenomeGraph.remove((i4, i3))
    
    # add colored edges: (i1, i3) & (i2, i4) to graph
    GenomeGraph.append((i1, i3))
    GenomeGraph.append((i2, i4))
    
    return GenomeGraph

def TwoBreakDistance(P: List[List[int]], Q: List[List[int]]) -> int:
    # pass
    # get edges 
    p_red_edges = ColoredEdges(P)
    q_blue_edges = ColoredEdges(Q)
    colored_edges = p_red_edges + q_blue_edges
    
    # adj list to store graph 
    adj_list = defaultdict(list)
    
    # for each edge (u, v), 
    # add v to adj list of u 
    # add u to adj list of v 
    for u, v in colored_edges:
        adj_list[u].append(v)
        adj_list[v].append(u)

    # let p be 0, red and q be 1,blue
    cycles_found = []
    visited_nodes = set()
    # BFS from each node in graph 
    for node in adj_list:
        if node not in visited_nodes:
            # update visited nodes from BFS nodes 
            bfs_res= (BFS(node, 'p', adj_list))
            visited_nodes.update(bfs_res)
            # add visited BFS nodes as cycle 
            cycles_found.append(list(visited_nodes))
    
    # get 2 break distance 
    num_colored = len(colored_edges)
    two_break_dist = num_colored / 2 - len(cycles_found) 
    return int(two_break_dist)

def BFS(node, seq, adj_list):
        # queue for first node and seq
        queue = [(node, seq)] 
        # keep track of visited nodes
        visited = set()  
        # first node visisted 
        visited.add(node)  
        while queue:
            # remove node 
            node, seq = queue.pop(0)  
            for neighbor in adj_list[node]:
                if neighbor not in visited:
                    # change seq 
                    if seq == 'p': #p
                        next_seq = 'q'  # change to q
                    else: # q (=1)
                        next_seq = 'p'  # change to p
                    queue.append((neighbor, next_seq))  
                    # update visited node 
                    visited.add(neighbor)  
        return visited
    
def ColoredEdges(P: List[List[int]]) -> List[Tuple[int, int]]:
    Edges=[]
    for Chromosome in P: 
        Nodes = ChromosomeToCycle(Chromosome)
        #print(Nodes)
        for j in range (len(Chromosome)):
            edge= (Nodes[(2*j)-1], Nodes[(2*j)])
            Edges.append(edge)
    return Edges 
    
def ChromosomeToCycle(Chromosome: List[int]) -> List[int]:
    #pass
    nodes=[]
    for j in range(len(Chromosome)): 
        i = Chromosome[j]
        if i>0: 
            #print(2*i-1)
            nodes.append((2*i-1))
            #print(2*i)
            nodes.append(2*i)
        else: 
            #print(-2*i)
            nodes.append(-(2*i))
            
            nodes.append(-(2*i+1))
    return nodes 

def TwoBreakOnGenome(P: List[int],
                     i1: int, i2: int, i3: int, i4: int) -> List[List[int]]:
    #pass
    edges= ColoredEdges(P)
    graph = TwoBreakOnGenomeGraph(edges, i1, i2, i3, i4)
    genome = GraphToGenome(graph)
    return genome 

def TwoBreakOnGenomeGraph(GenomeGraph: genome_graph_t,
                          i1: int, i2: int, i3: int, i4: int) -> genome_graph_t:
    #pass
    # remove edges (i1, i2) & (i3, i4) from graph
    if (i1, i2) in GenomeGraph: 
        GenomeGraph.remove((i1, i2))
    else: 
        GenomeGraph.remove((i2, i1))
    
    if (i3, i4) in GenomeGraph: 
        GenomeGraph.remove((i3, i4))
    else: 
        GenomeGraph.remove((i4, i3))
    
    # add colored edges: (i1, i3) & (i2, i4) to graph
    GenomeGraph.append((i1, i3))
    GenomeGraph.append((i2, i4))
    #print("GENOMEGRAPH", GenomeGraph)
    return GenomeGraph

def GraphToGenome(GenomeGraph: List[Tuple[int, int]]) -> List[List[int]]:
    # pass
    P = []
    seen = []
    adjacency_matrix = np.zeros(len(GenomeGraph)*2)
    
    for edge in GenomeGraph:
        adjacency_matrix[int(edge[0]-1)] = edge[1]-1
        adjacency_matrix[int(edge[1]-1)] = edge[0]-1

    for edge in GenomeGraph:
        start = edge[0]
        if start not in seen:
            seen.append(start)
            # if odd 
            if start % 2 != 0:
                end = start + 1
            # if even 
            else:
                end = start - 1
            sequence = []
            while True:
                if start % 2 == 0:
                    sequence.append(int((start) // 2)) #-
                else:
                    sequence.append(int(-(start + 1) // 2)) # no - 
                go_to = adjacency_matrix[int(start - 1)] + 1
                seen.append(go_to)
                if go_to == end:
                    P.append(sequence)
                    break
                if go_to % 2 == 0:
                    start = go_to - 1
                else:
                    start = go_to + 1
                seen.append(start)
    return P

# Extra Problem: 
def CycleToChromosome(Nodes: List[int]) -> List[int]:
    n= len(Nodes)//2
    Chromosome = [0] * (n)
    for j in range(n):
        if Nodes[2*j] < Nodes[2*j + 1]:
            Chromosome[j] = Nodes[2*j + 1] // 2
        else:
            Chromosome[j] = -Nodes[2*j] // 2
    return Chromosome


