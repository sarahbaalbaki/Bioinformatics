# Sarah Baalbaki 
# Module 2 

# Exercise 1: (2.3.2)
def genome_path(path: List[str]) -> str:
    """Forms the genome path formed by a collection of patterns."""
    text= path[0]
    k= len(path[0])
    for kmer in path[1:]: 
        text+= kmer[-1]
    return text 

#########################
# Exercise 2: (2.5.8)
def de_bruijn_kmers(k_mers: List[str]) -> Dict[str, List[str]]:
    """Forms the de Bruijn graph of a collection of k-mers."""
    # dict to store result 
    mapped_dict = {}
    # kmers to remove 
    keys_to_remove = set()
    # go over all kmers 
    for kmer in k_mers: 
        prefix = kmer[:-1]
        # if kmer prefix doesnt exist in dict, add it 
        if prefix not in mapped_dict: 
            mapped_dict[prefix] = []
        # if suffix not in dict, add it also 
        suffix = kmer[1:]
        if suffix not in mapped_dict: 
            mapped_dict[suffix] = []
        # add corresponding suffix to each prefix 
        mapped_dict[prefix].append(suffix)
    # if we have no suffix for a prefix, remove from dict 
    for key in mapped_dict:
        if not mapped_dict[key]:
            keys_to_remove.add(key)
    for key in keys_to_remove:
        mapped_dict.pop(key)
    return mapped_dict

#########################
# Exercise 3: (2.8.2)
def eulerian_cycle(g: Dict[int, List[int]]) -> Iterable[int]:
    # cycle to be outputted 
    cycle = []
    # use stack for eff 
    stack = deque()
    # randomly choose a node 
    current_node = random.choice(list(g.keys()))
    # while we still have unexplored edges 
    while True:
        # if curr node still has undiscovered edges 
        if current_node in g and len(g[current_node])!=0:
            # add node to still be discovered to stack, get new node, assign current to new node 
            stack.append(current_node)
            next_node = g[current_node].pop()
            current_node = next_node
        else:
            # node has no more unexplored edges, add node to cycle 
            cycle.append(current_node)
            if len(stack)==0:
                break
            current_node = stack.pop()
    # fix node order, reverse cycle 
    cycle = cycle[::-1]
    return cycle

#########################
# Exercise 4: (2.8.6) 
def eulerian_path(graph: Dict[int, List[int]]) -> Iterable[int]:
    # get indegree and outdegree for each node
    flattened_dict = flatten_dict(graph)
    # in degree for all incoming edges of node 
    indegree = {node: 0 for node in flattened_dict} 
    for node in indegree: 
        for key in graph: 
            if node in graph[key]: 
                indegree[node] += 1
    # print(indegree)  
    # out degree for all outgoing edges of node 
    outdegree = {node: 0 for node in flattened_dict}
    for node, neighbors in graph.items(): 
        outdegree[node]= len(neighbors)
    #print(outdegree)
    
    # get path start and end points 
    start_point = []
    end_point = []
    for node in flattened_dict:
        # if indegree less than out, it is start point  
        if outdegree[node] == indegree[node] + 1:
            start_point.append(node)
        # if indegree more than out, it is end point  
        elif indegree[node] == outdegree[node] + 1:
            end_point.append(node)
    #print(start_points, end_points)
    
    # starck for path starting from start_point 
    stack = deque([start_point[0]])
    path = []
    
    # while we still have undiscovered edges 
    while stack:
        # get current node as last eleemnt from stack 
        current_node = stack[-1]
        # if out degree of node is 0, remove it and add node to path
        if outdegree.get(current_node, 0) == 0:
            path.append(stack.pop())
        # if we still have edges coming from node, get the next node, decrement the out degree, and add the next node to the stack 
        else:
            next_node = graph[current_node].pop()
            outdegree[current_node] -= 1
            stack.append(next_node)
            
    # flip the path 
    return path[::-1]

def flatten_dict(d: Dict[int, List[int]]) -> List[int]:
    # get flattened dict with unique vals to contain all node vals 
    flattened = []
    for key, value in d.items():
        flattened.append(key)
        flattened.extend(value) 
    return np.unique(flattened)

#########################
# Exercise 5: (2.8.7) 
def string_reconstruction(patterns: List[str], k: int) -> str:
    """Reconstructs a string from its k-mer composition."""
    dB = de_bruijn_kmers(patterns)
    path = eulerian_path(dB)
    text = genome_path(path) 
    return text 

def de_bruijn_kmers(k_mers: List[str]) -> Dict[str, List[str]]:
    """Forms the de Bruijn graph of a collection of k-mers."""
    mapped_dict = defaultdict(list)
    for kmer in k_mers:
        prefix, suffix = kmer[:-1], kmer[1:]
        mapped_dict[prefix].append(suffix)
    #print(mapped_dict)
    return {key: value for key, value in mapped_dict.items() if value}
    
def eulerian_path(graph: Dict[int, List[int]]) -> Iterable[int]:
    # get indegree and outdegree for each node
    flattened_dict = flatten_dict(graph)
    # in degree for all incoming edges of node 
    indegree = {node: 0 for node in flattened_dict} 
    for node in indegree: 
        for key in graph: 
            if node in graph[key]: 
                indegree[node] += 1
    #print(indegree)  
    # out degree for all outgoing edges of node 
    outdegree = {node: 0 for node in flattened_dict}
    for node, neighbors in graph.items(): 
        outdegree[node]= len(neighbors)
    #print(outdegree)
    
    # get path start and end points 
    start_point = []
    end_point = []
    for node in flattened_dict:
        # if indegree less than out, it is start point  
        if outdegree[node] == indegree[node] + 1:
            start_point.append(node)
        # if indegree more than out, it is end point  
        elif indegree[node] == outdegree[node] + 1:
            end_point.append(node)
     
    if len(start_point)==0: 
        start_point.append(random.choice(list(indegree.keys())))
    #print(start_points, end_points)
    
    #print(start_point)
    # starck for path starting from start_point 
    stack = deque([start_point[0]])
    path = []
    
    # while we still have undiscovered edges 
    while stack:
        # get current node as last eleemnt from stack 
        current_node = stack[-1]
        # if out degree of node is 0, remove it and add node to path
        if outdegree.get(current_node, 0) == 0:
            path.append(stack.pop())
        # if we still have edges coming from node, get the next node, decrement the out degree, and add the next node to the stack 
        else:
            next_node = graph[current_node].pop()
            outdegree[current_node] -= 1
            stack.append(next_node)
            
    # flip the path 
    return path[::-1]

def flatten_dict(d: Dict[int, List[int]]) -> List[int]:
    # get flattened dict with unique vals to contain all node vals 
    flattened = []
    for key, value in d.items():
        flattened.append(key)
        flattened.extend(value) 
    return np.unique(flattened)

def genome_path(path: List[str]) -> str:
    """Forms the genome path formed by a collection of patterns."""
    return ''.join([path[0]] + [kmer[-1] for kmer in path[1:]])

#########################
# Exercise 6: (2.8.11) 
def k_universal_string(k: int) -> str:
    """Generates a k-universal circular string."""
    # generate all possible kmers 
    kmers = generate_kmers(k)
    # get the de bruijn graph
    de_bruijn= de_bruijn_kmers(kmers)
    # get a eulerian cycle from the graph 
    k_universal_circular_str = eulerian_cycle(de_bruijn)
    # get the path 
    path = genome_path(k_universal_circular_str)
    return path[:-(k-1)]

# function from before, get the joint path given consecutive kmers 
def genome_path(path: List[str]) -> str:
    """Forms the genome path formed by a collection of patterns."""
    text = path[0]
    for kmer in path[1:]:
        text += kmer[-1]
    return text

# function from before, get eulerian cycle rep. circular string 
def eulerian_cycle(g: Dict[str, List[str]]) -> Iterable[str]:
    # cycle to be outputted 
    cycle = []
    # use stack for eff 
    stack = deque()
    # randomly choose a node 
    current_node = random.choice(list(g.keys()))
    # while we still have unexplored edges 
    while True:
        # if curr node still has undiscovered edges 
        if current_node in g and len(g[current_node])!=0:
            # add node to still be discovered to stack, get new node, assign current to new node 
            stack.append(current_node)
            next_node = g[current_node].pop()
            current_node = next_node
        else:
            # node has no more unexplored edges, add node to cycle 
            cycle.append(current_node)
            if len(stack)==0:
                break
            current_node = stack.pop()
    # fix node order, reverse cycle 
    cycle = cycle[::-1]
    return cycle

# from prev code, get the kmers in dict (as de bruijn graph)
def de_bruijn_kmers(k_mers: List[str]) -> Dict[str, List[str]]:
    """Forms the de Bruijn graph of a collection of k-mers."""
    mapped_dict = {}
    keys_to_remove = set()
    for kmer in k_mers: 
        prefix = kmer[:-1]
        if prefix not in mapped_dict: 
            mapped_dict[prefix] = []
        suffix = kmer[1:]
        if suffix not in mapped_dict: 
            mapped_dict[suffix] = []
        mapped_dict[prefix].append(suffix)
    for key in mapped_dict:
        if not mapped_dict[key]:
            keys_to_remove.add(key)
    for key in keys_to_remove:
        mapped_dict.pop(key)
    return mapped_dict

# generate all the possbile kmers 
def generate_kmers(k):
    kmers = []
    for seq in product('01', repeat=k):
        kmers.append(''.join(seq))
    return kmers

#########################
# Exercise 7: (2.9.5)
def ContigGeneration(Pattern: List[str]) -> List[str]:
    contigs= []
    graph = de_bruijn_kmers(Pattern)
    #print(graph)
    paths = maximal_nonbranching_paths(graph)
    #print(paths)
    for path in paths: 
        contigs.append(genome_path(path))
    #print(contigs)
    return contigs

# Insert your maximal_non_branching_paths function here, along with any subroutines you need
def maximal_nonbranching_paths(g: Dict[int, List[int]]) -> Iterable[Iterable[int]]:
    """Finds maximal non-branching paths in a graph."""
    paths = []

    flattened_dict = flatten_dict(g)

    #print(flattened_dict)

    indegree = {node: 0 for node in flattened_dict} 
    for node in indegree: 
        for key in g: 
            for n in g[key]:
                if n==node: 
                    indegree[node] += 1
    #print(indegree)  
    # out degree for all outgoing edges of node 
    outdegree = {node: 0 for node in flattened_dict}
    for node, neighbors in g.items(): 
        outdegree[node]= len(neighbors)
    #print(outdegree)

    for node in g.keys(): 
        if indegree[node] != 1 or outdegree[node] != 1: 
            #print("node", node, indegree(g, node), outdegree(g, node))
            if outdegree[node] > 0: 
                outdegree[node]-=1
                for vals in g[node]:
                    nonbranching_path = [node, vals]
                    #indegree[vals]-=1
                    outdegree[node]-=1
                    while indegree[vals] == 1 and outdegree[vals] == 1: 
                        "HERE"
                        outdegree[vals]-=1
                        #print("NONBRANCHING_PATH", nonbranching_path)
                        nonbranching_path.append(g[vals][0])
                        vals = g[vals][0]
                    paths.append(nonbranching_path)

    g_new = {}
    #print("indegree", indegree)
    #print("outdegree", outdegree)
    for node in g: 
        if indegree[node]!=0 and outdegree[node]!=0: 
            g_new[node]= g[node]
    #print("GNEW", g_new)
    cycles = find_cycles(g_new) 
    #print("cycles", cycles)
    for cycle in cycles: 
        if cycle not in paths: 
            paths.append(cycle)
    #print(paths) 
    return paths

def find_cycles(graph: Dict[int, List[int]]) -> Iterable[Iterable[int]]:
    """Finds cycles in a graph."""
    cycles = []
    visited = set()
    for node in graph:
        cycle = dfs(graph, node, visited)
        if cycle:
            cycles.append(cycle)
    return cycles

def flatten_dict(d: Dict[int, List[int]]) -> List[int]:
    # get flattened dict with unique vals to contain all node vals 
    flattened = []
    for key, value in d.items():
        flattened.append(key)
        flattened.extend(value) 
    return np.unique(flattened)


def dfs(graph, start_node, visited):
    stack = [(start_node, [start_node])]
    while stack:
        node, path = stack.pop()
        if node in visited:
            cycle_start = path.index(node)
            if cycle_start != len(path) - 1:
                return path[cycle_start:]
            else:
                continue
        visited.add(node)
        for neighbor in graph.get(node, []):
            stack.append((neighbor, path + [neighbor]))
    return []

def find_cycles(graph: Dict[int, List[int]]) -> Iterable[Iterable[int]]:
    """Finds cycles in a graph."""
    cycles = []
    visited = set()
    for node in graph:
        cycle = dfs(graph, node, visited)
        if cycle:
            cycles.append(cycle)
    return cycles
    
def genome_path(path: List[str]) -> str:
    """Forms the genome path formed by a collection of patterns."""
    text= path[0]
    k= len(path[0])
    for kmer in path[1:]: 
        text+= kmer[-1]
    return text
   
# from prev code, get the kmers in dict (as de bruijn graph)
def de_bruijn_kmers(k_mers: List[str]) -> Dict[str, List[str]]:
    """Forms the de Bruijn graph of a collection of k-mers."""
    mapped_dict = {}
    keys_to_remove = set()
    for kmer in k_mers: 
        prefix = kmer[:-1]
        if prefix not in mapped_dict: 
            mapped_dict[prefix] = []
        suffix = kmer[1:]
        if suffix not in mapped_dict: 
            mapped_dict[suffix] = []
        mapped_dict[prefix].append(suffix)
    for key in mapped_dict:
        if not mapped_dict[key]:
            keys_to_remove.add(key)
    for key in keys_to_remove:
        mapped_dict.pop(key)
    return mapped_dict