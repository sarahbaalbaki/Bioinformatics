# Sarah Baalbaki 
# Bioinformatics- Module 6

# 6.6- Decoding Problem: 
def Viterbi_Decoding(observed_seq, alphabet, states, initial_ptobabilities, transition_probabilities, emission_probabilities):
    V = [{}]
    for s in states:
      V[0][s] = {"log_prob": np.log(initial_ptobabilities[s]) + np.log(emission_probabilities[s][alphabet[0]]), "prev": None}

    # loop starting from second observation to find the state with the highest proability 
    # assign the probability and the state to the dictionary in V for the corresponding observation index 
    for t in range(1, len(observed_seq)):
        # add a dictionary to the viterbi matrix for the time step we are at 
        V.append({})
        # loop over the states at the step we are at 
        for s in states:
          # variables for max transition probabilities and previous states 
          max_transition_prob = -np.inf
          previous_state = None
          # go over the previous states 
          for prev_st in states:
            # get the transition log probabilities 
            transition_prob = V[t-1][prev_st]["log_prob"] + np.log(transition_probabilities[prev_st][s])
            # update the transition log probabilities 
            if transition_prob > max_transition_prob:
              max_transition_prob = transition_prob
              previous_state = prev_st
          # get the max probability that leads to current stae 
          max_prob = max_transition_prob + np.log(emission_probabilities[s][observed_seq[t]])
          # add results of max prob and previous state into the dictionary at the step we are at 
          V[t][s] = {"log_prob": max_prob, "prev": previous_state}

    # array to store obtained states 
    decoded_states = []
    # maximum log probability counter 
    max_log_prob = -np.inf
    # optimal state 
    optimal_state = None

    # get most probable state and backtrack
    # loop over states adn observations 
    for st, data in V[-1].items():
        if data["log_prob"] > max_log_prob:
            max_log_prob = data["log_prob"]
            optimal_state = st
    # append result to the decoded states 
    decoded_states.append(optimal_state)
    # reset previous state 
    previous = optimal_state

    # follow backtrack till first observation
    for val in range(len(V) - 2, -1, -1):
        decoded_states.insert(0, V[val + 1][previous]["prev"])
        previous = V[val + 1][previous]["prev"]

    print("States " + "".join(decoded_states))
    return "".join(decoded_states)
observed_seq= "xzxyxzyzzyzzyxzzzzxyzxxzzzyzzzzyyzxyzzyyzyyyxzzyxzzyyxyxxxzxyzyxxzzzyxxxyzyzyyzzzxxyxyyxzzyyxzzyxxyx"

alphabet= ("x", "y", "z")

states= ("A", "B", "C", "D")

initial_ptobabilities = {"A": 0.25, "B": 0.25, "C": 0.25, "D": 0.25}

transition_probabilities = {
    "A": {"A": 0.301, "B": 0.28, "C": 0.164, "D": 0.255},
    "B": {"A": 0.065, "B": 0.385, "C": 0.032, "D": 0.518},
	"C": {"A": 0.015, "B": 0.111, "C": 0.394, "D": 0.48},
	"D": {"A": 0.036, "B": 0.255, "C": 0.085, "D": 0.624},
}

emission_probabilities = {
    "A": {"x": 0.205, "y": 0.059, "z": 0.736},
    "B": {"x": 0.703, "y": 0.08, "z": 0.217},
	"C": {"x": 0.092, "y": 0.466, "z": 0.442},
	"D": {"x": 0.447, "y": 0.395, "z": 0.158},
}

Viterbi_Decoding(observed_seq, alphabet, states, initial_ptobabilities, transition_probabilities, emission_probabilities)

#################################
# 6.7- Profile HMM 
def ProfileHMM(threshold, alphabet, alignment):
    is_included= FilterColumns(threshold, alignment)
    state_dict = MakeStateDictionary(is_included)
    symbols= GetCharDict(alphabet)
    emission_matrix = GetEmissions(symbols, alignment, is_included, state_dict)
    transition_matrix= GetTransitions(symbols, alignment, is_included, state_dict)
    emission_matrix_normalized = Normalize(emission_matrix)
    transition_matrix_normalized = Normalize(transition_matrix) 
    return transition_matrix_normalized,emission_matrix_normalized

def Normalize(matrix): 
    # get row sums, divide is sum !=0 
    row_sums = matrix.sum(axis=1)
    for i in range(len(matrix)): 
        if row_sums[i]>0: 
            matrix[i,:]= matrix[i,:]/row_sums[i]
    return matrix

def GetTransitions(symbols, alignment, is_included, state_dict): 
    list_of_strings = [list(string) for string in alignment]
    list_of_strings= np.array(list_of_strings)
    n_strings= list_of_strings.shape[0] 
    total_n_states = list_of_strings.shape[1] 

    # init transition matrix 
    transition_matrix= np.zeros((len(state_dict),len(state_dict)))

    state_cols = np.zeros(total_n_states, dtype=int)
    count_gaps = 0
    for j in range(len(is_included)):
        if is_included[j]==0: 
            count_gaps += 1
        # add count to col_to_state
        state_cols[j] = count_gaps

    # loop over strings
    for i in range(n_strings):
        # start at state S
        start_state = 'S'
        # # loop over all chars in string 
        for j in range(total_n_states): 
            # state 
            s = state_cols[j]
            # if insertion at end 
            if is_included[j]==0 and list_of_strings[i, j] in symbols:
                state_now = 'M' + str(s)
            elif is_included[j]==1 and list_of_strings[i, j] in symbols:
                state_now = 'I' + str(s)
            elif is_included[j]==0 and list_of_strings[i, j] not in symbols:
                state_now = 'D' + str(s)
            elif is_included[j]==1 and list_of_strings[i][j] not in symbols: 
                if j == total_n_states-1:
                    state_now= 'E'
                
            s1=state_dict[start_state]
            s2= state_dict[state_now]
            # print(s1, s2)
            transition_matrix[s1, s2]+=1
            start_state= state_now 

            if j == total_n_states - 1 and not (is_included[j]==1 and list_of_strings[i][j] not in symbols):
                state_now='E'
                s1 = state_dict[start_state]
                transition_matrix[s1, state_dict[state_now]] += 1

    return transition_matrix

def GetEmissions(symbols, alignment, is_included, state_dict): 
    list_of_strings = [list(string) for string in alignment]
    list_of_strings= np.array(list_of_strings)
    n_strings= list_of_strings.shape[0] 
    total_n_states = list_of_strings.shape[1] 

    state_cols = np.zeros(total_n_states, dtype=int)
    count_gaps = 0
    for j in range(len(is_included)):
        if is_included[j]==0: 
            count_gaps += 1
        # add count to col_to_state
        state_cols[j] = count_gaps
    # print(state_cols)

    # init matrix 
    emission_matrix= np.zeros((len(state_dict),len(symbols)))

    # loop over strings
    for i in range(n_strings):
        for j in range(total_n_states):
            # get symbol
            symbol = list_of_strings[i][j]
            if symbol in symbols: 
                # index 
                # print("is_included", is_included)
                # print(is_included[:j])
                s = state_cols[j]
                c = symbols[symbol]
                # get if I or M 
                if is_included[j]==1: 
                    st= 'I'
                else: 
                    st= 'M'
                # entry to emission matrix, with value 
                # print("st", st, "s", s)
                entry = st+ str(s)
                # print("entry", entry)
                # print('c', c)
                emission_matrix[state_dict[entry], c]+=1
    return emission_matrix 

def FilterColumns(threshold, alignment): 
    list_of_strings = [list(string) for string in alignment]
    list_of_strings= np.array(list_of_strings)
    n_strings= list_of_strings.shape[0] 
    total_n_states = list_of_strings.shape[1] 

    # bool array to store if col included or not 
    is_included=[] 
    # count_included= 0 # store count of included cols 
    # go over all states 
    for j in range(total_n_states):
        # sum the gaps for each 
        sum_gaps = np.sum(list_of_strings[:, j] == '-')
        # if gaps less than threshold, col is included else it is not  
        if sum_gaps < n_strings * threshold:
            #count_gaps += 1
            is_included.append(0)
        else: 
            is_included.append(1)
    return is_included

def GetCharDict(alphabet): 
    chars = {}
    for i in range(len(alphabet)):
        s = alphabet[i]
        chars[s] = i
    return chars 

def MakeStateDictionary(is_included): 
    # get nuber of columns included after thresholding 
    # sum_cols= 0 
    # for col in is_included:
    #     sum_cols += col
    
    # sum_cols= np.cumsum(is_included)
    # sum_cols_max= max(sum_cols)
    state_cols = np.zeros(len(is_included), dtype=int)
    count_gaps = 0
    for j in range(len(is_included)):
        if is_included[j]==0: 
            count_gaps += 1
        # add count to col_to_state
        state_cols[j] = count_gaps
    sum_cols= max(state_cols)

    # init S, I0, Mi, Di, Ii, and E states 
    state_dict= {'S': 0, 'I0': 1}
    # add the M, D, and I states based on number of cols included 
    # state_dict.update({f'M{i}': 3 * i - 1 for i in range(1, sum_cols + 2)})
    # state_dict.update({f'D{i}': 3 * i for i in range(1, sum_cols + 2)})
    # state_dict.update({f'I{i}': 3 * i + 1 for i in range(1, sum_cols + 2)})   
    for i in range(1, sum_cols + 1):
        state_dict[f'M{i}'] = 3 * i - 1
        state_dict[f'D{i}'] = 3 * i
        state_dict[f'I{i}'] = 3 * i + 1
    state_dict['E']= 3 * sum_cols + 2

    return state_dict

#################################
# 7.8- HMM with Pseudocounts 
def ProfileHMM_Pseudocount(threshold, alphabet, alignment, pseudocount):
    is_included= FilterColumns(threshold, alignment)
    state_dict = MakeStateDictionary(is_included)
    symbols= GetCharDict(alphabet)
    emission_matrix = GetEmissions_Pseudo(symbols, alignment, is_included, state_dict, pseudocount)
    transition_matrix= GetTransitions_Pseudo(symbols, alignment, is_included, state_dict, pseudocount)
    emission_matrix_normalized = Normalize(emission_matrix)
    transition_matrix_normalized = Normalize(transition_matrix) 
    return transition_matrix_normalized,emission_matrix_normalized

def Normalize(matrix): 
    # get row sums, divide is sum !=0 
    row_sums = matrix.sum(axis=1)
    for i in range(len(matrix)): 
        if row_sums[i]>0: 
            matrix[i,:]= matrix[i,:]/row_sums[i]
    return matrix

def GetTransitions_Pseudo(symbols, alignment, is_included, state_dict, pseudocount): 
    list_of_strings = [list(string) for string in alignment]
    list_of_strings= np.array(list_of_strings)
    n_strings= list_of_strings.shape[0] 
    total_n_states = list_of_strings.shape[1] 

    # init transition matrix 
    transition_matrix= np.zeros((len(state_dict),len(state_dict)))

    state_cols = np.zeros(total_n_states, dtype=int)
    count_gaps = 0
    for j in range(len(is_included)):
        if is_included[j]==0: 
            count_gaps += 1
        # add count to col_to_state
        state_cols[j] = count_gaps

    # loop over strings
    for i in range(n_strings):
        # start at state S
        start_state = 'S'
        # # loop over all chars in string 
        for j in range(total_n_states): 
            # state 
            s = state_cols[j]
            # if insertion at end 
            if is_included[j]==1 and list_of_strings[i][j] in symbols:
                state_now = 'M' + str(s)
                print("1", state_now)
            elif is_included[j]==0 and list_of_strings[i][j] in symbols:
                state_now = 'I' + str(s)
                print("2", state_now)
            elif is_included[j]==1 and list_of_strings[i][j] not in symbols:
                state_now = 'D' + str(s)
                print("3", state_now)
            elif is_included[j]==0 and list_of_strings[i][j] not in symbols: 
                if j == total_n_states-1:
                    state_now= 'E'
                    print("4", state_now)
                
            s1=state_dict[start_state]
            print("state_now", state_now)
            s2= state_dict[state_now]
            print("s2", s2)
            print("here", s1, s2)
            transition_matrix[s1, s2]+=1
            start_state= state_now 

            if j == total_n_states - 1 and not (is_included[j]==0 and list_of_strings[i][j] not in symbols):
                state_now= 'E'
                print("5", state_now)
                s1 = state_dict[start_state]
                # s2 = state_now
                transition_matrix[s1, state_dict[state_now]] += 1
                print("done")
    
    transition_matrix= Normalize(transition_matrix)
    
   # add pseudocounts 
    for i in range(0,2):
        # for first 2, pseudocounts by 2*3 jump 
        a = 1
        b = a+2
        transition_matrix[i,a:b+1] += pseudocount
    # for rest of states, pseudocounts by 3*3 jumps 
    for i in range(2, len(transition_matrix)-1, 3):
        print("i", i)
        a = 4+(i-2)
        b = a+3
        print(a, b)
        transition_matrix[i:i+3, a:b] += pseudocount

    return transition_matrix

def GetEmissions_Pseudo(symbols, alignment, is_included, state_dict, pseudocount): 
    list_of_strings = [list(string) for string in alignment]
    list_of_strings= np.array(list_of_strings)
    n_strings= list_of_strings.shape[0] 
    total_n_states = list_of_strings.shape[1] 

    state_cols = np.zeros(total_n_states, dtype=int)
    count_gaps = 0
    for j in range(len(is_included)):
        if is_included[j]==0: 
            count_gaps += 1
        # add count to col_to_state
        state_cols[j] = count_gaps
    # print(state_cols)

    # init matrix 
    emission_matrix= np.zeros((len(state_dict),len(symbols)))

    # loop over strings
    for i in range(n_strings):
        for j in range(total_n_states):
            # get symbol
            symbol = list_of_strings[i][j]
            if symbol in symbols: 
                # index 
                # print("is_included", is_included)
                # print(is_included[:j])
                s = state_cols[j]
                c = symbols[symbol]
                # get if I or M 
                if is_included[j]==1: 
                    st= 'I'
                else: 
                    st= 'M'
                # entry to emission matrix, with value 
                # print("st", st, "s", s)
                entry = st+ str(s)
                # print("entry", entry)
                # print('c', c)
                emission_matrix[state_dict[entry], c]+=1
    
    print(emission_matrix)
    emission_matrix= Normalize(emission_matrix)
    print(emission_matrix)
    # loop over all states , except E
    for i in range(len(emission_matrix)-1):
        # if state is not % 3, it is not D (dont add pseudo to del) -> add pseudocounts 
        if i % 3 != 0:
            # add pseudocounts for all the emitted symbols 
            for j in range(len(emission_matrix[0])):
                emission_matrix[i, j] += pseudocount
    
    return emission_matrix 

def FilterColumns(threshold, alignment): 
    list_of_strings = [list(string) for string in alignment]
    list_of_strings= np.array(list_of_strings)
    n_strings= list_of_strings.shape[0] 
    total_n_states = list_of_strings.shape[1] 

    # bool array to store if col included or not 
    is_included=[] 
    # count_included= 0 # store count of included cols 
    # go over all states 
    for j in range(total_n_states):
        # sum the gaps for each 
        sum_gaps = np.sum(list_of_strings[:, j] == '-')
        # if gaps less than threshold, col is included else it is not  
        if sum_gaps < n_strings * threshold:
            #count_gaps += 1
            is_included.append(0)
        else: 
            is_included.append(1)
    return is_included

def GetCharDict(alphabet): 
    chars = {}
    for i in range(len(alphabet)):
        s = alphabet[i]
        chars[s] = i
    return chars 

def MakeStateDictionary(is_included): 
    # get nuber of columns included after thresholding 
    state_cols = np.zeros(len(is_included), dtype=int)
    count_gaps = 0
    for j in range(len(is_included)):
        if is_included[j]==0: 
            count_gaps += 1
        # add count to col_to_state
        state_cols[j] = count_gaps
    sum_cols= max(state_cols)

    # init S, I0, Mi, Di, Ii, and E states 
    state_dict= {'S': 0, 'I0': 1}
    # add the M, D, and I states based on number of cols included  
    for i in range(1, sum_cols + 1):
        state_dict[f'M{i}'] = 3 * i - 1
        state_dict[f'D{i}'] = 3 * i
        state_dict[f'I{i}'] = 3 * i + 1
    state_dict['E']= 3 * sum_cols + 2

    return state_dict

#################################
# 7.8- SeqAlignment with Profile HMM: 
import numpy as np 

def SeqAlignmentWithProfileHMM(string_x,threshold_theta,pseudocount,alphabet,alignment): 
    # get transition and meission matrix from profile HMM with pseudo, return the alphabet state dictionary too to obtain states 
    transition_matrix_normalized,emission_matrix_normalized,state_dict= ProfileHMM_Pseudocount(threshold_theta, alphabet, alignment, pseudocount)
    optimal_path = ViterbiGraph(string_x, alphabet, state_dict, transition_matrix_normalized, emission_matrix_normalized)
    return optimal_path

def ViterbiGraph(string_x, alphabet, state_dict, transition_matrix, emission_matrix): 
    # get state array for states in order from state dict 
    states = [None] * len(state_dict)
    for k,v in state_dict.items():
        states[v] = k
    
    # get len_str and num_states 
    len_str = len(string_x) 
    num_states= len(states)//3  
    
    # initialize scoring and backtrack matrices 
    scores_matrix= MakeScores(num_states, len_str)
    backtrack_matrix= MakeBacktrack(num_states, len_str)

    # follow topological order
    # initially column of dels 
    # then column of I0, M1, D1, I1, ... until end states 
    # we have columns = num_obs emitted 

    # set first column deletion nodes
    # if  HMM enters Deletion(1) from the initial state, it can move down through deletion states 
    # before transitioning to a match or insertion state in the first column
    scores_matrix[2][1,0] = 0
    backtrack_matrix[2][1,0] = None
    # set first column- scores and backtrack to deletion nodes 
    for l in range(2,num_states):
        prev_deletion = states.index('D'+str(l-1))
        current_d = states.index('D'+str(l))
        scores_matrix[2][l,0] = scores_matrix[2][l-1,0] + np.log(transition_matrix[prev_deletion,current_d])
        backtrack_matrix[2][l,0] = ('D',l-1,0)
        
    # update scores and backtrack with insertions intial
    scores_matrix, backtrack_matrix= SetInsertions(scores_matrix, backtrack_matrix, len_str, num_states, states, string_x, transition_matrix, emission_matrix)

    # update scores and backtrack with matches initial
    scores_matrix, backtrack_matrix= SetMatches(scores_matrix, backtrack_matrix, len_str, num_states, states, string_x, transition_matrix, emission_matrix)
    
    # update scores and backtrack with deletions initial 
    scores_matrix, backtrack_matrix= SetDeletions(scores_matrix, backtrack_matrix, len_str, num_states, states, string_x, transition_matrix, emission_matrix)

    # # update scores and backtrack for remaning cols observations 
    scores_matrix, backtrack_matrix= GetScoresAndBacktrackMatrices(scores_matrix, backtrack_matrix, len_str, num_states, states, string_x, transition_matrix, emission_matrix)

    # get path, backtrack through matrices
    path = Backtrack(scores_matrix, backtrack_matrix, len_str, num_states, transition_matrix, states)
    
    return path

def MakeScores(num_states, len_str): 
    scores_matrix = np.zeros(shape=(3, num_states, len_str+1), dtype=float)
    return scores_matrix

def MakeBacktrack(num_states, len_str): 
    backtrack_matrix = np.empty(shape = (3, num_states,len_str+1), dtype = tuple)
    return backtrack_matrix

# let 0 be insertion matrix and insertion scores
# let 1 be match matrix and match scores
# let 2 be deletion matrix and deletion scores

def SetInsertions(scores_matrix, backtrack_matrix, len_str, num_states, states, string_x, transition_matrix, emission_matrix): 
    i0 = states.index('I0')
    s = states.index('S')
    obs = alphabet.index(string_x[0])
    # initialize the scores and backtrack matrices for all cols for the top row for insertion- I0 in all cols 
    scores_matrix[0][0,1] = np.log(transition_matrix[s,i0]*emission_matrix[i0,obs])
    backtrack_matrix[0][0,1] = None  
    for j in range(2,len_str+1): 
        obs = alphabet.index(string_x[j-1])
        scores_matrix[0][0,j] = scores_matrix[0][0,j-1] + np.log(transition_matrix[i0,i0]*emission_matrix[i0,obs])
        backtrack_matrix[0][0,j] = ('I',0,j-1)
    
    # initialize second column - first insertions 
    for i in range(1,num_states):
        current_deletion = states.index('D'+str(i))
        current_insertion = states.index('I'+str(i))
        obs = alphabet.index(string_x[0])
        # update score and backtrack based on first observation  
        scores_matrix[0][i,1] = scores_matrix[2][i,0] + np.log(transition_matrix[current_deletion,current_insertion]*emission_matrix[current_insertion,obs])
        backtrack_matrix[0][i,1] = ('D',i,0)

    # set rest of insertion columns vals for score and backtrack 
    for j in range(2,len_str+1):
        i = 1
        current_deletion = states.index('D'+str(i))
        current_m = states.index('M'+str(i))
        current_insertion = states.index('I'+str(i))
        obs = alphabet.index(string_x[j-1])
        match = scores_matrix[1][i,j-1] + np.log(transition_matrix[current_m,current_insertion]*emission_matrix[current_m,obs])
        deletion = scores_matrix[2][i,j-1] + np.log(transition_matrix[current_deletion,current_insertion]*emission_matrix[current_m,obs])
        insertion = scores_matrix[0][i,j-1] + np.log(transition_matrix[current_insertion,current_insertion]*emission_matrix[current_m,obs])
        max_score = max(match,deletion,insertion)
        # update score and backtrack based on observation sequence 
        if max_score == match:
            scores_matrix[0][1,j]= match
            backtrack_matrix[0][1,j]= ('M',i,j-1)
        if max_score == deletion:
            scores_matrix[0][1,j]= deletion
            backtrack_matrix[0][1,j]= ('D',i,j-1)
        if max_score == insertion:
            scores_matrix[0][1,j]= insertion
            backtrack_matrix[0][1,j]= ('I',i,j-1)
    
    return scores_matrix, backtrack_matrix 

def SetMatches(scores_matrix, backtrack_matrix, len_str, num_states, states, string_x, transition_matrix, emission_matrix): 
    m1 = states.index('M1')
    s = states.index('S')
    obs = alphabet.index(string_x[0])
    scores_matrix[1][1,1] = np.log(transition_matrix[s,m1]*emission_matrix[m1,obs])
    backtrack_matrix[1][1,1] = None  
    i0= states.index('I0')
    # initialize the scores and backtrack matrices for all cols for the first M row- M1 in all cols 
    for j in range(2,len_str+1):
        obs = alphabet.index(string_x[j-1])
        scores_matrix[1][1,j] = scores_matrix[0][0,j-1] + np.log(transition_matrix[i0,m1]*emission_matrix[m1,obs])
        backtrack_matrix[1][1,j] = ('I',0,j-1)
    
    # initialize second column - first matches 
    for i in range(2,num_states):
        prev_deletion = states.index('D'+str(i-1))
        current_m = states.index("M" +str(i)) 
        obs = alphabet.index(string_x[0])
        # update score and backtrack based on first observation
        scores_matrix[1][i,1] = scores_matrix[2][i-1,0] + np.log(transition_matrix[prev_deletion,current_m]*emission_matrix[current_m,obs])
        backtrack_matrix[1][i,1] = ('D',i-1,0)
    return scores_matrix, backtrack_matrix

def SetDeletions(scores_matrix, backtrack_matrix, len_str, num_states, states, string_x, transition_matrix, emission_matrix): 
    # go over all states and set del for all rows, first col 
    d1 = states.index('D1')
    s = states.index('S')
    scores_matrix[2][1,0] = np.log(transition_matrix[s,d1])
    backtrack_matrix[2][1,0] = None   
    i0 = states.index('I0')
    # initialize the scores and backtrack matrices for all cols for the first D row- D1 in all cols 
    for j in range(1,len_str+1): 
        scores_matrix[2][1,j] = scores_matrix[0][0,j] + np.log(transition_matrix[i0,d1])
        backtrack_matrix[2][1,j] = ('I',0,j)
    
    for i in range(2,num_states):
        j=1 
        prev_deletion = states.index('D'+str(i-1))
        prev_m = states.index('M'+str(i-1))
        prev_insertion = states.index('I'+str(i-1))
        current_d = states.index('D'+str(i))
        match = scores_matrix[1][i-1,j]  + np.log(transition_matrix[prev_m,current_d])
        deletion = scores_matrix[2][i-1,j] + np.log(transition_matrix[prev_deletion,current_d])
        insertion = scores_matrix[0][i-1,j] + np.log(transition_matrix[prev_insertion,current_d])
        max_score = max(match,deletion,insertion)
        if max_score == match:
            scores_matrix[2][i,1]= match
            backtrack_matrix[2][i,1]= ('M',i-1,j)
        if max_score == deletion:
            scores_matrix[2][i,1]= deletion
            backtrack_matrix[2][i,1]= ('D',i-1,j)
        if max_score == insertion:
            scores_matrix[2][i,1]= insertion
            backtrack_matrix[2][i,1]= ('I',i-1,j)
    return scores_matrix, backtrack_matrix

def GetScoresAndBacktrackMatrices(scores_matrix, backtrack_matrix, len_str, num_states, states, string_x, transition_matrix, emission_matrix): 
    # go over all states and all observations 
    # fill in score and backtrack matrix 
    for i in range(2,num_states):
            for j in range(2,len_str+1):
                # get previous state indices for 3 cases 
                prev_deletion = states.index('D'+str(i-1))
                prev_m = states.index('M'+str(i-1))
                prev_insertion = states.index('I'+str(i-1))
                
                # get current state indices for 3 cases 
                current_deletion = states.index('D'+str(i))
                current_m = states.index("M" +str(i)) 
                current_insertion = states.index('I'+str(i))
                
                # get previous observation 
                obs = alphabet.index(string_x[j-1])

                # get scores for match 
                m_match = scores_matrix[1][i-1,j-1] + np.log(transition_matrix[prev_m,current_m]*emission_matrix[current_m,obs])
                m_deletion = scores_matrix[2][i-1,j-1] + np.log(transition_matrix[prev_deletion,current_m]*emission_matrix[current_m,obs])
                m_insertion = scores_matrix[0][i-1,j-1] + np.log(transition_matrix[prev_insertion,current_m]*emission_matrix[current_m,obs])
                max_score_matches = max(m_match,m_deletion,m_insertion)
                # update scores and backtrack (dim=1 for match) based on match
                if max_score_matches == m_match:
                    scores_matrix[1][i,j]=m_match
                    backtrack_matrix[1][i,j]= ('M',i-1,j-1)
                if max_score_matches == m_deletion:
                    scores_matrix[1][i,j]=m_deletion
                    backtrack_matrix[1][i,j]= ('D',i-1,j-1)
                if max_score_matches == m_insertion:
                    scores_matrix[1][i,j]=m_insertion
                    backtrack_matrix[1][i,j]= ('I',i-1,j-1)
                
                # get scores for deletion 
                d_match = scores_matrix[1][i-1,j]  + np.log(transition_matrix[prev_m,current_deletion])
                d_deletion = scores_matrix[2][i-1,j] + np.log(transition_matrix[prev_deletion,current_deletion])
                d_insertion = scores_matrix[0][i-1,j] + np.log(transition_matrix[prev_insertion,current_deletion])
                max_score_deletion = max(d_match,d_deletion,d_insertion)
                # update scores and backtrack (dim=2 for del) based on deletion
                if max_score_deletion == d_match:
                    scores_matrix[2][i,j]=d_match
                    backtrack_matrix[2][i,j]= ('M',i-1,j)
                if max_score_deletion == d_deletion:
                    scores_matrix[2][i,j]=d_deletion
                    backtrack_matrix[2][i,j]= ('D',i-1,j)
                if max_score_deletion == d_insertion:
                    scores_matrix[2][i,j]=d_insertion
                    backtrack_matrix[2][i,j]= ('I',i-1,j)
                
                # get scores for insertion 
                i_match = scores_matrix[1][i,j-1] + np.log(transition_matrix[current_m,current_insertion]*emission_matrix[current_m,obs])
                i_deletion = scores_matrix[2][i,j-1] + np.log(transition_matrix[current_deletion,current_insertion]*emission_matrix[current_m,obs])
                i_insertion = scores_matrix[0][i,j-1] + np.log(transition_matrix[current_insertion,current_insertion]*emission_matrix[current_m,obs])
                max_score_insertion = max(i_match,i_deletion,i_insertion)
                # update scores and backtrack (dim=1 for match) based on match
                if max_score_insertion == i_match:
                    scores_matrix[0][i,j]= i_match
                    backtrack_matrix[0][i, j]= ('M',i,j-1)
                if max_score_insertion == i_deletion:
                    scores_matrix[0][i,j]= i_deletion
                    backtrack_matrix[0][i, j]= ('D',i,j-1)
                if max_score_insertion == i_insertion:
                    scores_matrix[0][i,j]= i_insertion
                    backtrack_matrix[0][i, j]= ('I',i,j-1)
            
    return scores_matrix, backtrack_matrix

def Backtrack(scores_matrix, backtrack_matrix, len_str, num_states, transition_matrix, states): 
    # get state index before "E"
    before_last_deletion = states.index('D'+str(num_states-1))
    before_last_m = states.index('M'+str(num_states-1))
    before_last_insertion = states.index('I'+str(num_states-1))
    obs_end = states.index('E')
    
    # get scores for end match, del, insertion
    end_match = scores_matrix[1][num_states-1,len_str] + np.log(transition_matrix[before_last_m,obs_end])
    end_deletion = scores_matrix[2][num_states-1,len_str] + np.log(transition_matrix[before_last_deletion,obs_end])
    end_insertion = scores_matrix[0][num_states-1,len_str] + np.log(transition_matrix[before_last_insertion,obs_end])
    end_max_score = max(end_deletion,end_match,end_insertion)

    # array to store path 
    path = []
    if end_max_score == end_match:
        # ends with a match, add thos as end state and backtrack (set next step) based on prev state 
        path.append('M'+str(num_states-1))
        matrix_next = backtrack_matrix[1][num_states-1,len_str]
    elif end_max_score == end_deletion:
        # ends with a deletion, add thos as end state and backtrack (set next step) based on prev state 
        path.append('D'+str(num_states-1))
        matrix_next = backtrack_matrix[2][num_states-1,len_str]
    elif end_max_score == end_insertion:
        # ends with a insertion, add thos as end state and backtrack (set next step) based on prev state 
        path.append('I'+str(num_states-1))
        matrix_next = backtrack_matrix[0][num_states-1,len_str]
    
    state= matrix_next[0]
    i = matrix_next[1]
    j = matrix_next[2]
    
    # until we get to beginning
    while True:
        # get state and set next matrix to the matrix of our state 
        if state == 'I':
            emit = 'I'+str(i)
            matrix_next = backtrack_matrix[0]
        if state == 'M':
            emit = 'M'+str(i)
            matrix_next = backtrack_matrix[1]
        if state == 'D':
            emit = 'D'+str(i)
            matrix_next = backtrack_matrix[2]
        
        # append to the path the state that enitted the observation
        path.append(emit)
       
        if matrix_next[i,j] == None:
            break
        # print(matrix_next[i,j])
        
        # get next state, i, j based on "next_step"
        next_step = matrix_next[i,j]
        state= next_step[0]
        i = next_step[1]
        j = next_step[2]
    
    return path 

def ProfileHMM_Pseudocount(threshold, alphabet, alignment, pseudocount):
    is_included= FilterColumns(threshold, alignment)
    state_dict = MakeStateDictionary(is_included)
    symbols= GetCharDict(alphabet)
    emission_matrix = GetEmissions_Pseudo(symbols, alignment, is_included, state_dict, pseudocount)
    transition_matrix= GetTransitions_Pseudo(symbols, alignment, is_included, state_dict, pseudocount)
    emission_matrix_normalized = Normalize(emission_matrix)
    transition_matrix_normalized = Normalize(transition_matrix) 
    return transition_matrix_normalized,emission_matrix_normalized,state_dict

def Normalize(matrix): 
    # get row sums, divide is sum !=0 
    row_sums = matrix.sum(axis=1)
    for i in range(len(matrix)): 
        if row_sums[i]>0: 
            matrix[i,:]= matrix[i,:]/row_sums[i]
    return matrix

def GetTransitions_Pseudo(symbols, alignment, is_included, state_dict, pseudocount): 
    list_of_strings = [list(string) for string in alignment]
    list_of_strings= np.array(list_of_strings)
    n_strings= list_of_strings.shape[0] 
    total_n_states = list_of_strings.shape[1] 

    # init transition matrix 
    transition_matrix= np.zeros((len(state_dict),len(state_dict)))

    state_cols = np.zeros(total_n_states, dtype=int)
    count_gaps = 0
    for j in range(len(is_included)):
        if is_included[j]==0: 
            count_gaps += 1
        # add count to col_to_state
        state_cols[j] = count_gaps

    # loop over strings
    for i in range(n_strings):
        # start at state S
        start_state = 'S'
        # # loop over all chars in string 
        for j in range(total_n_states): 
        #     # # state 
            s = state_cols[j]
            # if insertion at end 
            if list_of_strings[i][j] not in symbols and is_included[j]==1: 
                if j == total_n_states-1:
                    s1= state_dict[start_state]
                    s2= state_dict['E']
                    transition_matrix[s1, s2]+=1
            else: 
                if is_included[j]==0 and list_of_strings[i, j] in symbols:
                    state_now = 'M' 
                elif is_included[j]==1 and list_of_strings[i, j] in symbols:
                    state_now = 'I' 
                elif is_included[j]==0 and list_of_strings[i, j] not in symbols:
                    state_now = 'D' 
                else:
                    break
                        # state 
                
                state_now = state_now + str(s)
                s1=state_dict[start_state]
                s2= state_dict[state_now]
                # print(s1, s2)
                transition_matrix[s1, s2]+=1
                start_state= state_now 

                if j == total_n_states - 1:
                    state_now=state_dict['E']
                    s1 = state_dict[start_state]
                    transition_matrix[s1, state_now] += 1
                        # state 
            # s = state_cols[j]
            # # if insertion at end 
            # if is_included[j]==1 and list_of_strings[i][j] in symbols:
            #     state_now = 'M' + str(s)
            #     print("1", state_now)
            # elif is_included[j]==0 and list_of_strings[i][j] in symbols:
            #     state_now = 'I' + str(s)
            #     print("2", state_now)
            # elif is_included[j]==1 and list_of_strings[i][j] not in symbols:
            #     state_now = 'D' + str(s)
            #     print("3", state_now)
            # elif is_included[j]==0 and list_of_strings[i][j] not in symbols: 
            #     if j == total_n_states-1:
            #         state_now= 'E'
            #         print("4", state_now)
                
            # s1=state_dict[start_state]
            # print("state_now", state_now)
            # s2= state_dict[state_now]
            # print("s2", s2)
            # print("here", s1, s2)
            # transition_matrix[s1, s2]+=1
            # start_state= state_now 

            # if j == total_n_states - 1 and not (is_included[j]==0 and list_of_strings[i][j] not in symbols):
            #     state_now= 'E'
            #     print("5", state_now)
            #     s1 = state_dict[start_state]
            #     # s2 = state_now
            #     transition_matrix[s1, state_dict[state_now]] += 1
            #     print("done")
    
    transition_matrix= Normalize(transition_matrix)
    
   # add pseudocounts 
    for i in range(0,2):
        # for first 2, pseudocounts by 2*3 jump 
        a = 1
        b = a+2
        transition_matrix[i,a:b+1] += pseudocount
    # for rest of states, pseudocounts by 3*3 jumps 
    for i in range(2, len(transition_matrix)-1, 3):
        # print("i", i)
        a = 4+(i-2)
        b = a+3
        # print(a, b)
        transition_matrix[i:i+3, a:b] += pseudocount
    # for i in name_state.values():
    #     csum = transition[i,:].sum()
    #     if (csum > 0):
    #         transition[i,:] /= csum

    return transition_matrix

def GetEmissions_Pseudo(symbols, alignment, is_included, state_dict, pseudocount): 
    list_of_strings = [list(string) for string in alignment]
    list_of_strings= np.array(list_of_strings)
    n_strings= list_of_strings.shape[0] 
    total_n_states = list_of_strings.shape[1] 

    state_cols = np.zeros(total_n_states, dtype=int)
    count_gaps = 0
    for j in range(len(is_included)):
        if is_included[j]==0: 
            count_gaps += 1
        # add count to col_to_state
        state_cols[j] = count_gaps
    # print(state_cols)

    # init matrix 
    emission_matrix= np.zeros((len(state_dict),len(symbols)))

    # loop over strings
    for i in range(n_strings):
        for j in range(total_n_states):
            # get symbol
            symbol = list_of_strings[i][j]
            if symbol in symbols: 
                # index 
                # print("is_included", is_included)
                # print(is_included[:j])
                s = state_cols[j]
                c = symbols[symbol]
                # get if I or M 
                if is_included[j]==1: 
                    st= 'I'
                else: 
                    st= 'M'
                # entry to emission matrix, with value 
                # print("st", st, "s", s)
                entry = st+ str(s)
                # print("entry", entry)
                # print('c', c)
                emission_matrix[state_dict[entry], c]+=1
    
    # print(emission_matrix)
    emission_matrix= Normalize(emission_matrix)
    # print(emission_matrix)
    # loop over all states , except E
    for i in range(len(emission_matrix)-1):
        # if state is not % 3, it is not D (dont add pseudo to del) -> add pseudocounts 
        if i % 3 != 0:
            # add pseudocounts for all the emitted symbols 
            for j in range(len(emission_matrix[0])):
                emission_matrix[i, j] += pseudocount
    
    return emission_matrix 

def FilterColumns(threshold, alignment): 
    list_of_strings = [list(string) for string in alignment]
    list_of_strings= np.array(list_of_strings)
    n_strings= list_of_strings.shape[0] 
    total_n_states = list_of_strings.shape[1] 

    # bool array to store if col included or not 
    is_included=[] 
    # count_included= 0 # store count of included cols 
    # go over all states 
    for j in range(total_n_states):
        # sum the gaps for each 
        sum_gaps = np.sum(list_of_strings[:, j] == '-')
        # if gaps less than threshold, col is included else it is not  
        if sum_gaps < n_strings * threshold:
            #count_gaps += 1
            is_included.append(0)
        else: 
            is_included.append(1)
    return is_included

def GetCharDict(alphabet): 
    chars = {}
    for i in range(len(alphabet)):
        s = alphabet[i]
        chars[s] = i
    return chars 

def MakeStateDictionary(is_included): 
    # get nuber of columns included after thresholding 
    state_cols = np.zeros(len(is_included), dtype=int)
    count_gaps = 0
    for j in range(len(is_included)):
        if is_included[j]==0: 
            count_gaps += 1
        # add count to col_to_state
        state_cols[j] = count_gaps
    sum_cols= max(state_cols)

    # init S, I0, and E states 
    state_dict= {'S': 0, 'I0': 1}
    # add the M, D, and I states based on number of cols included 
    for i in range(1, sum_cols + 1):
        state_dict[f'M{i}'] = 3 * i - 1
        state_dict[f'D{i}'] = 3 * i
        state_dict[f'I{i}'] = 3 * i + 1
    state_dict['E']= 3 * sum_cols + 2

    return state_dict

string_x = 'AEFDFDC'
theta = 0.4 
pseudocount = 0.01
alphabet = ['A','B','C','D','E','F']
alignment = ['ACDEFACADF',
             'AFDA---CCF',
             'A--EFD-FDC',
             'ACAEF--A-C',
             'ADDEFAAADF']
SeqAlignmentWithProfileHMM(string_x,theta,pseudocount,alphabet,alignment)

string_x = 'AEDEAEEEBADDEDEDECABDAEECEBCBCDAECBEBAECDADEACCAB'
theta = 0.356
pseudocount = 0.01
alphabet = ['A','B','C','D','E']
alignment = ["AEDAAEEEDAD--DED-ABBAAE-C--CCA-BEC-CB-ECEADEADDAB",
"AED-AECEB-DDEDE-BB-CCAEEA-BCB-DBE-ACAAE-CA-EAADAB",
"AE-CAEE-B-DDEDEDDB-BEA-ECABCBA---E-B-AEEE-DEA-DA-",
"AAD-AE-ECAD-EDEDBBABAA-E-EB--ADBECCCBA-EEDDEAA--E",
"-CBEAEEEBA---CEDDBABBC-ECEBCB--BEC-CDA-EEC-D-A-AB",
"-DDEAEC-EADDEAEDBBABA-EECEBCBAABECBCB-EEE-DEAADAE",
"CEDEAE-EBADD-DEDBBE-BD-BCCBC-ADBE-BCBBBE-AD-AEDAB",
"ACA--EEEB-ADDDEDBBBB-C-ECEE-EABB-C-CBACEEACAADDAB",
"-EDEAEEEBBDDEDEDB-ABEAEECEB-BADBEB-C-AEEAAD-CADAB"]
SeqAlignmentWithProfileHMM(string_x,theta,pseudocount,alphabet,alignment)
