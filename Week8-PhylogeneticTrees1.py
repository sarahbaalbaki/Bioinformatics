# Sarah Baalbaki 
# Bioinformatics- Week 8 

# 8.2- Distances Between Leaves Problem: 
WeightedEdge = Tuple[int, int, int]
DistMatrix = List[List[int]]

def DistancesBetweenLeaves(n: int, edge_list: List[WeightedEdge]) -> DistMatrix:
    #pass
    # get max node index for num_nodes 
    max_node_index = 0
    for edge in edge_list:
        max_node_index = max(max_node_index, edge[0], edge[1])
    
    num_nodes= max_node_index+1
    
    # get distance matrix
    dist_matrix = [[float('inf')] * (num_nodes) for _ in range(num_nodes)]

    # set diag to 0
    for i in range(num_nodes):
        dist_matrix[i][i] = 0

    # fill mat with edge weights
    for node1, node2, weight in edge_list:
        dist_matrix[node1][node2] = weight
        dist_matrix[node2][node1] = weight

    # find shortest distances and set in matrix 
    for k in range(num_nodes):
        for i in range(num_nodes):
            for j in range(num_nodes):
                dist_matrix[i][j] = min(dist_matrix[i][j], dist_matrix[i][k] + dist_matrix[k][j])

    # get dist btw leaves
    leaf_distances = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            leaf_distances[i][j] = dist_matrix[i][j]
            leaf_distances[j][i] = dist_matrix[i][j]

    return leaf_distances

#########################################
# 8.3- Limb Length Problem: 
DistMatrix = List[List[float]]
# Insert your LimbLength function here, along with any subroutines you need
def LimbLength(j: int, D: DistMatrix) -> float:
    n= len(D)
    length = np.inf
    for i in range(n):
        if i != j:
            for k in range(n):
                if k != j:
                    length = min(length, (D[i][j] + D[j][k] - D[i][k]) // 2)
    return length

#########################################
# 8.6- UPGMA: 
DistMatrix = List[List[int]]
WeightedEdge = Tuple[int, int, float]

# Insert your UPGMA function here, along with any subroutines you need
def UPGMA(D: DistMatrix) -> List[WeightedEdge]:
        # get number of nodes initially 
        num_nodes = len(D)

        # n single-element clusters labeled 1, ... , n: [cluster index, num_elements]
        clusters = []
        for cluster in range(num_nodes):
            clusters.append([cluster, 1])

	# print(clusters)

        # fix dist matrix diagonal so it is not the min element (since 0) 
        DistMatrix= np.array(D)

        # graph T with n isolated nodes 
        T = [[i] for i in range(num_nodes)]

        # age dictionary for nodes
        age= {}
        for node in clusters: 
            age[node[0]]= 0 

        # while there is more than one cluster 
        while len(clusters)>0:
            # find the two closest clusters Ci and Cj
            c1, c2 = GetClosest(DistMatrix)
            # print(c1, c2)

            # merge Ci and Cj into a new cluster Cnew with |Ci| + |Cj| elements
            index_new_m = len(T)
            # new node Cnew (index, clusters)
            num_elements_m = clusters[c1][1] + clusters[c2][1]
            cluster_new = [index_new_m, num_elements_m]
            # print(cluster_new)

            # add a new node labeled by cluster Cnew to T
            C_new = [clusters[c1][0], clusters[c2][0]]
            T.append(C_new)
            # print(T)
            
            # remove the nodes from T since they are in new []             
            T[c1]= []
            T[c2]= []

            # Age(Cnew) ← DCi, Cj / 2
            age[index_new_m]= DistMatrix[c1, c2] / 2

            # if only 2 clusters, break as if we merge we get one last cluster  
            if 2 == len(DistMatrix):
                break
            
            # get new weights 
            row_c1= DistMatrix[c1,:]
            row_c2 = DistMatrix[c2,:]
            num_elements_c1= clusters[c1][1]
            num_elements_c2= clusters[c2][1]
            num_elements_total= num_elements_c1+ num_elements_c2
            new_row = (row_c1*num_elements_c1 + row_c2*num_elements_c2) / (num_elements_total)
            # print(d_new) 
            new_row = np.delete(new_row, [c1, c2], axis=0)
            # print(d_new)

            #remove the rows and columns of D corresponding to Ci and Cj
            DistMatrix= Remove_Rows_Columns(DistMatrix, [c1, c2])
            # add row to D for Cnew by computing D(Cnew, C) for each C in Clusters
            DistMatrix = np.vstack([DistMatrix, new_row])
            
            # add column to D for Cnew by computing D(Cnew, C) for each C in Clusters
            new_col = np.insert(new_row, len(new_row), np.inf, axis = 0)
            new_col = np.reshape(new_col, (-1, 1))
            DistMatrix = np.hstack([DistMatrix, new_col])

            # # remove Ci and Cj from Clusters
            removed_c1_c2 = []
            for idx, cluster in enumerate(clusters):
            # Check if the index is not equal to c1 or c2
                if idx not in [c1, c2]:
                    # If the condition is met, add the cluster to the new list
                    removed_c1_c2.append(cluster)
            clusters = removed_c1_c2

            # add Cnew to Clusters
            clusters.append(cluster_new)

        # length of (v, w) ← Age(v) - Age(w)
        T_fixed=[]
        for i in range(len(T)):
            new_arr = []
            for j in range(len(T[i])):
                v = T[i][j]
                age_new = abs(age[i] - age[v])
                if age_new > 0: 
                    new_arr.append((v, age_new))
            T_fixed.append(new_arr)
        
        #print("T FIXED", T_fixed)

        adj_list = FormatList(T_fixed)
        #print(adj_list)

        # format the tree printed as needed
        #Tree_Final = GetTree(T_fixed)
        # return Tree_Final
        
        return adj_list

# reove rows and cols of closest clusters c1 and c2
def Remove_Rows_Columns(DistMatrix, indices):
    # reove rows
    DistMatrix = np.delete(DistMatrix, indices, axis=0)
    # remove columns
    DistMatrix = np.delete(DistMatrix, indices, axis=1)
    return DistMatrix
     
# get closest clusters from distance matrix 
def GetClosest(DistMatrix): 
    min= np.inf 
    c1= np.inf 
    c2 = np.inf 
    for i in range(len(DistMatrix)): 
        for j in range(len(DistMatrix[0])): 
            if i !=j:  # do not take diag since = 0, not min 
                if DistMatrix[i][j]<min: 
                    min = DistMatrix[i][j]
                    c1 = i
                    c2 = j 
    return c1, c2 

# format based on desired output 
def FormatList(T): 
    adj_list=[]
    for i in range(len(T)):
        for j in range(len(T[i])):
            node, weight = T[i][j]
            adj_list.append([i, node, weight])
        #adj_list.append(new_arr)
    return adj_list

# print the tree out as output given, autograder does this 
#def GetTree(T):
    #for i in range(len(T)):
        #for j in range(len(T[i])):
            #node, weight = T[i][j]
            #print(str(i)+ '->' + str(node) + ':' + str(weight))