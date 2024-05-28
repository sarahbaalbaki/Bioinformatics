# Sarah Baalbaki 
# Bioinformatics- Module 10 

# 10.3- Trie Construction: 
def TrieConstruction(Patterns: List[str]) -> List[Tuple[int, int, str]]:
    # pass
    # trie will have tuple (node, parent, symbol)
    trie = [(0, 0, '')]  
    # trie = []  
    # index for next tree node
    node_next = 1
    # go through all patterns 
    for p in Patterns:
        # start rom root in each pattern 
        node_now = 0
        # go through all characters in each pattern 
        for character in p:
            # flag for if node is child node 
            found_child = False
            # go through already present children to see if we have a matrch 
            for child in range(len(trie)):
                # if we have node with same parent node and same symbol
                if trie[child][0] == node_now and trie[child][2] == character:
                    # print("here", node_now, child)
                    # update node to child node 
                    node_now = child
                    # set flag of found to true 
                    found_child = True
                    break
            # if no match found 
            if not found_child:
                # print("here: ", (node_now, node_next, character))
                # add the node to the tree 
                trie.append((node_now, node_next, character))
                # update new node to curr node 
                node_now = node_next
                # move next node by 1 
                node_next += 1
    # remove init element
    del trie[0]
    sorted_trie = SortTrie(trie)
    return sorted_trie

def SortTrie(trie): 
    sorted_trie = []
    # go over all trie entries 
    for entry in trie:
        # go over el in sorted list 
        for i in range(len(sorted_trie)):
            # if el shoould be before n element in sorted list 
            if entry[0] < sorted_trie[i][0]:
                # add it at correct index
                sorted_trie.insert(i, entry)
                break
        else:
            # add at end, bigger 
            sorted_trie.append(entry)
    return sorted_trie

###########################################
# 10.5- Suffix Tree: 
def suffix_tree(text: str) -> List[str]:
    """
    Construct a suffix tree from the given text.
    """
    #pass
    # get the trie 
    Trie = ModifiedSuffixTrieConstruction(text) 
    # get adj graph to pass into maximal nonbranching 
    graph= MakeGraph(Trie)
    # graph with symbols for edge reconstruction
    symbol_graph= MakeGraphWithLetters(Trie)
    # get the non-branching paths 
    non_branching= maximal_non_branching_paths(graph)
    # print(non_branching)

    # store edge labels 
    edge_labels=[]
    # for nodes on a non-branching path, compress to an edge 
    for edge in non_branching: 
        # label of the compressed edge 
        label=""
        # go through nodes on this path 
        for index in range(len(edge)-1):
            # add char to the label 
            for key, value in symbol_graph.items(): 
                for v in value: 
                    if key== edge[index] and v[0]==edge[index+1]: 
                        label+=v[1]
        edge_labels.append(label)
    edge_labels.append("$")
    # return edge labels 
    return edge_labels

# code from suffix trie construction, modified for Text instead of Patterns 
def ModifiedSuffixTrieConstruction(Text: str) -> List[Tuple[int, int, str]]:
    # trie will have tuple (node, parent, symbol)
    trie = [(0, 0, '$')]  
    # go through the text 
    for i in range(len(Text)):
        # start rom root in each position 
        node_now = 0
        # go through all characters 
        for j in range(i, len(Text)):
            # get the current character 
            symbol = Text[j]
            # flag if the edge exists or not in the tree already 
            edge_exists = False
            for edge_index in range(len(trie)):
                # if we have node with same parent node and same symbol
                if trie[edge_index][0] == node_now and trie[edge_index][2] == symbol:
                    # update node to child node 
                    node_now = trie[edge_index][1]
                    edge_exists = True
                    break
            # if no match found
            if not edge_exists:
                node_next = len(trie)
                # add the node to the tree 
                trie.append((node_now, node_next, symbol))
                # update curr node to new node 
                node_now = node_next
            # at end of text 
            if j == len(Text) - 1:
                for edge_index in range(len(trie)):
                    if trie[edge_index][0] == node_now and trie[edge_index][1] == len(Text):
                        trie[edge_index] = (node_now, len(Text), Text[i:])
    # print(SortTrie(trie[1:]))
    # return SortTrie(trie[1:])  
    # remove init element
    del trie[0]
    return trie

# graph as adj list for parent and child connections 
def MakeGraph(Tree): 
    g={}
    for entry in Tree: 
        if entry[0] not in g: 
            g[entry[0]]= [entry[1]]
        else: 
            g[entry[0]].append(entry[1])
    return(g)

# # graph as adj list for parent and child connections and symbols
def MakeGraphWithLetters(Tree): 
    g={}
    for entry in Tree: 
        ans=()
        if entry[0] not in g: 
            g[entry[0]]= [(entry[1], entry[2])]
        else: 
            g[entry[0]].append((entry[1], entry[2]))
    return(g)

# code from previous maximal non branching paths 
def maximal_non_branching_paths(g: Dict[int, List[int]]) -> Iterable[Iterable[int]]:
    """Finds maximal non-branching paths in a graph."""
    paths = []
    for node in g.keys(): 
        if indegree(g, node) != 1 or outdegree(g, node) != 1: 
            #print("node", node, indegree(g, node), outdegree(g, node))
            if outdegree(g, node) > 0: 
                for vals in g[node]:
                    nonbranching_path = [node, vals]
                    while indegree(g, vals) == 1 and outdegree(g, vals) == 1: 
                        nonbranching_path.append(g[vals][0])
                        vals = g[vals][0]
                    paths.append(nonbranching_path)
    return paths
        
def indegree(graph: Dict[int, List[int]], node: int) -> int:
    """Returns the indegree of a node in the graph."""
    indeg = 0
    for v in graph.values():
        if node in v:
            indeg += 1
    #print(node, indeg)
    return indeg

def outdegree(graph: Dict[int, List[int]], node: int) -> int:
    """Returns the outdegree of a node in the graph."""
    #print(node, len(graph.get(node, [])))
    return len(graph.get(node, []))

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

###########################################
# 10.5- Longest Repeat 
def longest_repeat(text: str) -> str:
    """
    Find the longest repeated substring in the given text.
    """
    #pass
    max_length=0
    longest_repeat=""
    path=[]
    text+='$'
    edges, paths, symbol_graph= suffix_tree(text)

    # loop over all edges 
    for i in range (len(edges)-1): 
        # if len(edges[i]) >= max_length and len(edges[i])>1 and edges[i][-1]!="$":
        # print("path[0]", paths[i][-1])
        # if the node has children (branching), get the number of children
        if paths[i][-1] in symbol_graph: 
            num_children= len(symbol_graph[paths[i][-1]])
        else: #leaf, no children 
            num_children=0
        # len> max length and has children and last char is not end sign '$'
        if len(GetFromRoot(paths[i], paths)) >= max_length and num_children>1 and edges[i][-1]!="$": 
            # update the max length path and length 
            path= paths[i]
            max_length= len(GetFromRoot(paths[i], paths))
    # if no path found, return 
    if path==[]: 
        return ""
    # get the longest repeat from the path 
    for index in range(len(path)-1):
        for key, value in symbol_graph.items(): 
            for v in value: 
                if key== path[index] and v[0]==path[index+1]: 
                    longest_repeat+=v[1]
    return longest_repeat

# get path from root to current node index
def GetFromRoot(path, paths): 
    # start at first node of one path, which will be the end of the next 
    index_last = path[0]
    # until we reach the root 
    while index_last!= 0: 
        # get the path from the root until the last node in path 
        for p in paths: 
            if p[-1]== index_last: 
                for el in p[::-1]: 
                    if el not in path:
                        path.insert(0, el)
                index_last= p[0]
    # return the path from root to curr node 
    return path

def suffix_tree(text: str) -> List[str]:
    """
    Construct a suffix tree from the given text.
    """
    #pass
    # get the trie 
    Trie = ModifiedSuffixTrieConstruction(text) 
    # get adj graph to pass into maximal nonbranching 
    graph= MakeGraph(Trie)
    # graph with symbols for edge reconstruction
    symbol_graph= MakeGraphWithLetters(Trie)
    # get the non-branching paths 
    non_branching= maximal_non_branching_paths(graph)
    # print(non_branching)

    # store edge labels 
    edge_labels=[]
    # for nodes on a non-branching path, compress to an edge 
    for edge in non_branching: 
        # label of the compressed edge 
        label=""
        # go through nodes on this path 
        for index in range(len(edge)-1):
            # add char to the label 
            for key, value in symbol_graph.items(): 
                for v in value: 
                    if key== edge[index] and v[0]==edge[index+1]: 
                        label+=v[1]
        edge_labels.append(label)
    edge_labels.append("$")
    return edge_labels, non_branching, symbol_graph

# code from suffix trie construction, modified for Text instead of Patterns 
def ModifiedSuffixTrieConstruction(Text: str) -> List[Tuple[int, int, str]]:
    # trie will have tuple (node, parent, symbol)
    trie = [(0, 0, '$')]  
    # go through the text 
    for i in range(len(Text)):
        # start rom root in each position 
        node_now = 0
        # go through all characters 
        for j in range(i, len(Text)):
            # get the current character 
            symbol = Text[j]
            # flag if the edge exists or not in the tree already 
            edge_exists = False
            for edge_index in range(len(trie)):
                # if we have node with same parent node and same symbol
                if trie[edge_index][0] == node_now and trie[edge_index][2] == symbol:
                    # update node to child node 
                    node_now = trie[edge_index][1]
                    edge_exists = True
                    break
            # if no match found
            if not edge_exists:
                node_next = len(trie)
                # add the node to the tree 
                trie.append((node_now, node_next, symbol))
                # update curr node to new node 
                node_now = node_next
            # at end of text 
            if j == len(Text) - 1:
                for edge_index in range(len(trie)):
                    if trie[edge_index][0] == node_now and trie[edge_index][1] == len(Text):
                        trie[edge_index] = (node_now, len(Text), Text[i:])
    # print(SortTrie(trie[1:]))
    # return SortTrie(trie[1:])  
    # remove init element
    del trie[0]
    return trie

# graph as adj list for parent and child connections 
def MakeGraph(Tree): 
    g={}
    for entry in Tree: 
        if entry[0] not in g: 
            g[entry[0]]= [entry[1]]
        else: 
            g[entry[0]].append(entry[1])
    return(g)

# # graph as adj list for parent and child connections and symbols
def MakeGraphWithLetters(Tree): 
    g={}
    for entry in Tree: 
        ans=()
        if entry[0] not in g: 
            g[entry[0]]= [(entry[1], entry[2])]
        else: 
            g[entry[0]].append((entry[1], entry[2]))
    return(g)

# code from previous maximal non branching paths 
def maximal_non_branching_paths(g: Dict[int, List[int]]) -> Iterable[Iterable[int]]:
    """Finds maximal non-branching paths in a graph."""
    paths = []
    for node in g.keys(): 
        if indegree(g, node) != 1 or outdegree(g, node) != 1: 
            #print("node", node, indegree(g, node), outdegree(g, node))
            if outdegree(g, node) > 0: 
                for vals in g[node]:
                    nonbranching_path = [node, vals]
                    while indegree(g, vals) == 1 and outdegree(g, vals) == 1: 
                        nonbranching_path.append(g[vals][0])
                        vals = g[vals][0]
                    paths.append(nonbranching_path)
    return paths
        
def indegree(graph: Dict[int, List[int]], node: int) -> int:
    """Returns the indegree of a node in the graph."""
    indeg = 0
    for v in graph.values():
        if node in v:
            indeg += 1
    #print(node, indeg)
    return indeg

def outdegree(graph: Dict[int, List[int]], node: int) -> int:
    """Returns the outdegree of a node in the graph."""
    #print(node, len(graph.get(node, [])))
    return len(graph.get(node, []))

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

###########################################
# 10.5- Longest Shared Substring: 
def longest_shared_substring(text1: str, text2: str) -> str:
    """
    Find the longest shared substring between two texts.
    """
    # pass

    # add end cahr to each string 
    text1+='#'
    text2+='$'
    text= text1+ text2

    max_length=0
    longest_repeat=""
    path=[]

    edges, paths, symbol_graph, edge_labels_colored = suffix_tree(text1, text2, text)
    # print("edges", edges)
    # print("paths", paths)
    # print("symbol_graph", symbol_graph)
    
    #edge_labels_colored=[]
    #for edge in edge_labels_c: 
        #if edge[1]==0: 
            #edge_labels_colored.append(edge)

    for i in range (len(edge_labels_colored)-1): 
        if paths[i][-1] in symbol_graph: 
            num_children= len(symbol_graph[paths[i][-1]])
        else: 
            num_children=0
        # get color, if =0 then in both strings 
        color = edge_labels_colored[i][1]
        if len(GetFromRoot(paths[i], paths)) > max_length and num_children>1 and edge_labels_colored[i][0][-1]!="$" and edge_labels_colored[i][0][-1]!="#" and color==0: 
            path= paths[i]
            # print("PATH", path)
            max_length= len(GetFromRoot(paths[i], paths))
    if path==[]: 
        return ""
    for index in range(len(path)-1):
        for key, value in symbol_graph.items(): 
            for v in value: 
                if key== path[index] and v[0]==path[index+1]: 
                    longest_repeat+=v[1]
    return longest_repeat

# get path from root to current node index
def GetFromRoot(path, paths): 
    index_last = path[0]
    while index_last!= 0: 
        for p in paths: 
            if p[-1]== index_last: 
                for el in p[::-1]: 
                    if el not in path:
                        path.insert(0, el)
                index_last= p[0]
    return path

def suffix_tree(text1, text2, text: str) -> List[str]:
    Trie = ModifiedSuffixTrieConstruction(text1, text2, text) 
    graph = MakeGraph(Trie)
    symbol_graph = MakeGraphWithLetters(Trie)
    non_branching = maximal_non_branching_paths(graph)

    # get symbol_graph 
    symbol_graph_lookup = {}
    for key, value in symbol_graph.items():
        for v in value:
            symbol_graph_lookup[(key, v[0])] = v[1]

    # get edge labels 
    edge_labels = []
    for edge in non_branching:
        label = ""
        for i in range(len(edge) - 1):
            label += symbol_graph_lookup.get((edge[i], edge[i + 1]), "")
        edge_labels.append(label)

    edge_labels.append("$")

    # get edge labels colored, for string diff 
    edge_labels_colored = []
    for edge in edge_labels: 
        color = -1
        # edge in both text 1 and text 2: color=0
        if edge in text1 and edge in text2: 
            color = 0
        # edge only in text 1: color=1
        elif edge in text1 and edge not in text2: 
            color = 1
        # edge only in text 2: color=2
        elif edge in text2 and edge not in text1: 
            color = 2
        edge_labels_colored.append((edge, color))

    return edge_labels, non_branching, symbol_graph, edge_labels_colored

# code for suffix trie construction
# modified prev implementation for time eff
def ModifiedSuffixTrieConstruction(text1: str, text2: str, text: str) -> List[Tuple[int, int, str]]:
    trie = {(0, '$'): 0}  
    end_node = len(text)  
    
    for i in range(len(text)):
        current_node = 0
        for j in range(i, len(text)):
            current_symbol = text[j]
            transition = (current_node, current_symbol)
            if transition in trie:
                current_node = trie[transition]
            else:
                new_node = len(trie)
                trie[transition] = new_node
                current_node = new_node
            if j == len(text) - 1:
                # update leaf node with suffix
                trie[(current_node, '$')] = end_node
    
    # change trie dict to a list of tuples
    trie_list = [(node, next_node, symbol) for (node, symbol), next_node in trie.items()]
    return trie_list


# graph as adj list for parent and child connections 
def MakeGraph(Tree): 
    g={}
    for entry in Tree: 
        if entry[0] not in g: 
            g[entry[0]]= [entry[1]]
        else: 
            g[entry[0]].append(entry[1])
    return(g)

# # graph as adj list for parent and child connections and symbols
def MakeGraphWithLetters(Tree): 
    g={}
    for entry in Tree: 
        ans=()
        if entry[0] not in g: 
            g[entry[0]]= [(entry[1], entry[2])]
        else: 
            g[entry[0]].append((entry[1], entry[2]))
    return(g)

# modified since prev code is time inefficient 
def maximal_non_branching_paths(g: Dict[int, List[int]]) -> Iterable[Iterable[int]]:
    paths = []
    indeg, outdeg = get_in_out_degrees(g)
    for node in g.keys():
        if indeg[node] != 1 or outdeg[node] != 1:
            if outdeg[node] > 0:
                for vals in g[node]:
                    nonbranching_path = [node, vals]
                    while indeg[vals] == 1 and outdeg[vals] == 1:
                        nonbranching_path.append(g[vals][0])
                        vals = g[vals][0]
                    paths.append(nonbranching_path)
    return paths

def get_in_out_degrees(graph):
    # get the in and out degrees of the graph
    indeg = defaultdict(int)
    outdeg = defaultdict(int)
    for source, targets in graph.items():
        outdeg[source] = len(targets)
        for target in targets:
            indeg[target] += 1
    return indeg, outdeg

###########################################
# 10.6- Suffix Array 
def suffix_array(text: str) -> List[int]:
    """
    Generate the suffix array for the given text.
    """
    # pass
    # text= text+'$'
    # dictionary to map each suffix to its index 
    arrays_with_indices={}
    for i in range(0, len(text)): 
        suffix= text[i:]
        arrays_with_indices[suffix]= i
    # print(arrays_with_indices)
    # sort the suffixes
    suffixes= list(arrays_with_indices.keys())
    suffixes.sort()
    # print(suffixes)
    # get the indices of the sorted suffixes 
    sorted_indices=[]
    for suffix in suffixes: 
        index = arrays_with_indices[suffix]
        sorted_indices.append(index)
    return sorted_indices