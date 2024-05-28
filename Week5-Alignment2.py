# Sarah Baalbaki 
# Bioinformatics- Module 5 

# .2- Global Alignment: 
import sys
from typing import List, Dict, Iterable, Tuple

def GlobalAlignment(match_reward: int, mismatch_penalty: int, indel_penalty: int, s: str, t: str) -> Tuple[int, str, str]:
    scores, backtrack = LCSBackTrack(match_reward,mismatch_penalty, indel_penalty, s, t)
    modified_s, modified_t = OutputLCS(backtrack, s, t, len(s), len(t))
    score_alignemt= scores[len(s)][len(t)]
    return score_alignemt, modified_s, modified_t

def LCSBackTrack(match_reward: int, mismatch_penalty: int, indel_penalty: int, v, w):
    len_v = len(v)
    len_w = len(w)
    s = [[0] * (len_w + 1) for _ in range(len_v + 1)]
    backtrack = [[0] * (len_w + 1) for _ in range(len_v + 1)]
    
    # init row and col indels 
    for i in range(len(v)+1):
        s[i][0]= -i*indel_penalty
        backtrack[i][0]= 1 # down 
    for j in range(1,len(w)+1):
        s[0][j]= -j*indel_penalty
        backtrack[0][j]= 2 # right 

    print(s)
    print(backtrack)

    for i in range(1, len_v + 1):
        for j in range(1, len_w + 1):
            #match = 1 if v[i - 1] == w[j - 1] else 0
            # if s[i-1][j-1]== s[i][j]: 
            #      score = match_reward 
            if v[i - 1] == w[j - 1]:
                score = match_reward 
            else: 
                 score = - mismatch_penalty
            s[i][j] = max(s[i-1][j] - indel_penalty, s[i][j-1] - indel_penalty, s[i-1][j-1] + score)
            if s[i][j] == s[i-1][j] - indel_penalty:
                backtrack[i][j] = 1  # down
            elif s[i][j] == s[i][j-1] - indel_penalty:
                backtrack[i][j] = 2  # right
            else:
                backtrack[i][j] = 3  # diagonal

    return s, backtrack

def OutputLCS(backtrack, v, w, i, j):
    v_aligned = ""
    w_aligned = ""
    
    while i > 0 and j > 0:
        if backtrack[i][j] == 3:  # diagonal
            v_aligned = v[i - 1] + v_aligned
            w_aligned = w[j - 1] + w_aligned
            i -= 1
            j -= 1
        elif backtrack[i][j] == 1:  # down
            v_aligned = v[i - 1] + v_aligned
            w_aligned = "-" + w_aligned
            i -= 1
        else:  # right
            v_aligned = "-" + v_aligned
            w_aligned = w[j - 1] + w_aligned
            j -= 1
    
    # Handle remaining characters in v or w, if any
    while i > 0:
        v_aligned = v[i - 1] + v_aligned
        w_aligned = "-" + w_aligned
        i -= 1
    while j > 0:
        v_aligned = "-" + v_aligned
        w_aligned = w[j - 1] + w_aligned
        j -= 1
    
    return v_aligned, w_aligned

match_score = 1 
mismatch_score = 1 
indel_score = 2
s= "GAGA"
t = "GAT"
GlobalAlignment(match_score, mismatch_score, indel_score, s, t)

###################################
# 5.2- Local Alignment: 
import sys
from typing import List, Dict, Iterable, Tuple
import numpy as np 

def max_entry_indices(matrix):
    max_value = float('-inf')
    max_i = max_j = -1

    for i in range(len(matrix)):
        for j in range(len(matrix[0])):
            if matrix[i][j] > max_value:
                max_value = matrix[i][j]
                max_i = i
                max_j = j

    return max_i, max_j

def LocalAlignment(match_reward: int, mismatch_penalty: int, indel_penalty: int, s: str, t: str) -> Tuple[int, str, str]:
    scores, backtrack = LCSBackTrack_Local(match_reward,mismatch_penalty, indel_penalty, s, t)
    print(scores) 
    max_i, max_j = max_entry_indices(scores)
    modified_s, modified_t = OutputLCS_Local(scores, backtrack, s, t, max_i, max_j)
    # score_alignemt= scores[len(s)][len(t)]
    score_alignemt= scores[max_i][max_j]
    return score_alignemt, modified_s, modified_t

def LCSBackTrack_Local(match_reward: int, mismatch_penalty: int, indel_penalty: int, v, w):
    len_v = len(v)
    len_w = len(w)
    s = [[0] * (len_w + 1) for _ in range(len_v + 1)]
    backtrack = [[0] * (len_w + 1) for _ in range(len_v + 1)]
    
    # init row and col indels 
    for i in range(len(v)+1):
        # s[i][0]= -i*indel_penalty
        s[i][0]=0
        backtrack[i][0]= 0 # down 
    for j in range(1,len(w)+1):
        # s[0][j]= -j*indel_penalty
        s[0][j]=0
        backtrack[0][j]= 0 # right 

    print(s)
    print(backtrack)

    for i in range(1, len_v + 1):
        for j in range(1, len_w + 1):
            #match = 1 if v[i - 1] == w[j - 1] else 0
            # if s[i-1][j-1]== s[i][j]: 
            #      score = match_reward 
            if v[i - 1] == w[j - 1]:
                score = match_reward 
            else: 
                 score = - mismatch_penalty
            s[i][j] = max(0, s[i-1][j] - indel_penalty, s[i][j-1] - indel_penalty, s[i-1][j-1] + score)
            if s[i][j] == s[i-1][j] - indel_penalty:
                backtrack[i][j] = 1  # down
            elif s[i][j] == s[i][j-1] - indel_penalty:
                backtrack[i][j] = 2  # right
            elif s[i][j]==0: 
                backtrack[i][j] = 4 # free ride 
            else: #s[i][j] == s[i-1][j] - indel_penalty:
                backtrack[i][j] = 3  # diagonal
            # elif s[i][j]==0: 
            #     backtrack[i][j] = 4 # free ride 

            # if s[i][j] < 0: 
            #     s[i][j] = 0
            #     backtrack[i-1][j-1] = 4 # free ride 

    return s, backtrack

def OutputLCS_Local(s, backtrack, v, w, i, j):
    v_aligned = ""
    w_aligned = ""

    while i > 0 and j > 0 and s[i][j] != 0:
        if backtrack[i][j] == 3:  # diagonal
            v_aligned = v[i - 1] + v_aligned
            w_aligned = w[j - 1] + w_aligned
            i -= 1
            j -= 1
        elif backtrack[i][j] == 1:  # down
            v_aligned = v[i - 1] + v_aligned
            w_aligned = "-" + w_aligned
            i -= 1
        elif backtrack[i][j] == 2:  # right
            v_aligned = "-" + v_aligned
            w_aligned = w[j - 1] + w_aligned
            j -= 1
        elif backtrack[i][j] == 4:  # freeride
            i -= 1
            j -= 1

    return v_aligned, w_aligned

match_score = 1 
mismatch_score = 1 
indel_score = 2
s= "GAGA"
t = "GAT"
LocalAlignment(match_score, mismatch_score, indel_score, s, t)

###################################
# 5.3- Edit Distance 
def EditDistance(s: str, t: str) -> int:
#     pass
    len_s= len(s)
    len_t= len(t)

    # initialize dynamic prog matrix of len s * t 
    matrix= [[0]* (len_t+1) for _ in range(len_s+1)]

    for i in range(len_s+1):
        for j in range(len_t+1):
            # first str empty, add all second string char
            if i == 0: 
                matrix[i][j]=j # number of char of 2nd str so far 
            # second str empty, add all char of first string 
            elif j == 0: 
                matrix[i][j]=i # number of char of 1nd str so far 
            # if last char are diff, check for min from all poss 
            elif s[i-1]!=t[j-1]: 
                # 1 for edit dist and chcek remaining for mismatches and indel  
                matrix[i][j]= 1 + min(matrix[i - 1][j - 1], matrix[i-1][j], matrix[i][j-1])
            # if letters match, continue checking for remianing str 
            else: #s[i-1]==t[j-1]
                matrix[i][j]= matrix[i-1][j-1]
    
    # retrn last dp entry in matrix for edit dist between s and t
    return matrix[len_s][len_t]

EditDistance("GAGA", "GAT")

###################################
# 5.3- Fitting Alignment: 
def FittingAlignment(s: str, t: str,
                     BLOSUM: Dict[str, Dict[str, int]], indel_penalty) -> Tuple[int, str, str]:
    # pass
    matrix, backtrack= MakeScoringMatrix(s, t, BLOSUM, indel_penalty)
    answer= BacktrackFitting(matrix, backtrack, s, t)
    return answer 

def MakeScoringMatrix(s, t, BLOSUM, indel_penalty): 
    len_s = len(s)
    len_t = len(t)

    matrix = [[0]* (len_t+1) for _ in range(len_s+1)]
    backtrack = [[0]* (len_t+1) for _ in range(len_s+1)]

    for i in range(0, len_s+1): 
        for j in range(1, len_t+1): 
            # if equal char, maatch 
            # if s[i - 1] == t[j - 1]:
            #     score_1 = matrix[i - 1][j - 1] + BLOSUM[s[i-1]][t[j-1]]
            # elif s[i-1]!=t[j-1]: 
            #     # print(s[i-1])
            #     # print(t[j-1])
            #     # print(BLOSUM[s[i-1]][t[j-1]])
            #     score_1 = matrix[i - 1][j - 1] + BLOSUM[s[i-1]][t[j-1]]
            # get max from ins, del, match, mis 
            matrix[i][j]= max(matrix[i][j-1]-1, matrix[i-1][j] - 1, matrix[i - 1][j - 1] + BLOSUM[s[i-1]][t[j-1]])
            if matrix[i][j] == matrix[i-1][j] - 1:
                backtrack[i][j] = 1  # down
            elif matrix[i][j] == matrix[i][j-1] - 1:
                backtrack[i][j] = 2  # right
            else: 
                backtrack[i][j] = 3  # diagonal

    return matrix, backtrack 

def find_max_score_index(matrix, s, t):
    max_score = -float('inf')
    max_i = -1
    for i in range(len(s), len(t)+1):
        score = matrix[i][len(s)]  # Score in the end column
        if score > max_score:
            max_score = score
            max_i = i
    return max_i

def BacktrackFitting(matrix, backtrack, s, t): 

    # get indices of v and w for ideal fitting end 
    j = len(t)
    i = find_max_score_index(matrix, t, s)

    # cut v and w and backtrack from end optimal position 
    s_aligned = s[:i]
    t_aligned = t[:j]

    score = matrix[i][j]
    #print(matrix)

    while i>0 and j>0 and matrix[i][j]!=0: # i*j!=0
        if backtrack[i][j] == 3:  # diagonal
            # simply decrement the indices
            i -= 1
            j -= 1
        elif backtrack[i][j] == 1:  # down
            # substitute char in w with - for indel , decrement i 
            t_aligned = t_aligned[:j]+ "-" + t_aligned[j:]
            i -= 1
        elif backtrack[i][j] == 2:  # right
            # substitute char in v with - for indel, decrement w 
            s_aligned = s_aligned[:i]+ "-" + s_aligned[i:]
            j -= 1
    
    # make them same length again and to align 
    s_aligned = s_aligned[i:]

    return score, s_aligned, t_aligned

###################################
# 5.3- Overlap Alignment 
def OverlapAlignment(match_reward: int, mismatch_penalty: int, indel_penalty: int,
                    s: str, t: str) -> Tuple[int, str, str]:
    #pass
    max_score, indices, matrix, backtrack= MakeScoringMatrix(s, t, match_reward, mismatch_penalty, indel_penalty)
    answer= BacktrackFitting(matrix, backtrack, s, t, max_score, indices)
    return answer 

def MakeScoringMatrix(s, t, match_reward, mismatch_penalty, indel_penalty): 
    len_s = len(s)
    len_t = len(t)

    matrix = [[0]* (len_t+1) for _ in range(len_s+1)]
    for j in range(len_t+1):
        matrix[0][j] = -j*indel_penalty
    backtrack = [[0]* (len_t+1) for _ in range(len_s+1)]
    for j in range(1, len_t+1):
        backtrack[0][j] = 2

    best_score = float("-inf")

    for i in range(1, len_s+1): 
        for j in range(1, len_t+1): 
            if s[i-1] == t[j-1]: 
                score = match_reward 
            else: # s[i-1]!=t[j-1]:
                score = - mismatch_penalty
            matrix[i][j]= max(matrix[i][j-1]-indel_penalty, matrix[i-1][j] - indel_penalty, matrix[i - 1][j - 1] + score) # 0 here? 
            if matrix[i][j] == matrix[i-1][j] - indel_penalty:
                backtrack[i][j] = 1  # down
            elif matrix[i][j] == matrix[i][j-1] - indel_penalty:
                backtrack[i][j] = 2  # right
            else: 
                backtrack[i][j] = 3  # diagonal
            
            #print(matrix)
            

            if i == len_s: 
                if matrix[i][j] >= best_score: 
                    #print(best_score, matrix[i][j])
                    best_score = matrix[i][j]
                    indices = (i, j)
            # index_i = len_s
            # #j = np.array(matrix[len_s][:]).argmax()

            # index_j = max(range(len(matrix[-1])), key=lambda x: matrix[-1][x])
                    

    #print("i", i, "j", j)
    # indices=(i,j)
    # indices=(index_i,index_j)
    i = indices[0]
    j=indices[1]
    best_score = matrix[i][j]
    
    print(matrix)

    return best_score, indices, matrix, backtrack 

def BacktrackFitting(matrix, backtrack, s, t, max_score, indices): 

    # get indices of v and w for ideal fitting end 
    i= indices[0]
    j= indices[1]

    print(i,j)
    
    s_aligned = ""
    t_aligned= ""

    #print("BACKTRACK")
    print(backtrack)
    #print(i, j)
    # if i==0: 

    while (j >0):
        #print(i, j)
        #print(s_aligned, t_aligned)
#        assert backtrack[i,j] != '*'
        #print("backtrack[i][j] ", backtrack[i][j] )
        if backtrack[i][j] == 1:
            #print("here 1")
            t_aligned+='-'
            s_aligned+=s[i-1]
            i -= 1
        elif backtrack[i][j] == 2:
            #print("here 2")
            t_aligned+=(t[j-1])
            s_aligned+='-'
            j -= 1
        elif backtrack[i][j]==3:
            #print("here 3")
            t_aligned+=t[j-1]
            s_aligned+=s[i-1]
            j -= 1
            i -= 1
        else: 
            break 
        
        print(s_aligned, t_aligned)

    # reverse after backtrcak 
    s_aligned = s_aligned[::-1]
    t_aligned = t_aligned[::-1]

    return max_score, s_aligned, t_aligned

###################################
# 5.4- Affine Gap Penalties 
from typing import Tuple
import numpy as np

def AffineAlignment(match_reward: int, mismatch_penalty: int, gap_opening_penalty: int, gap_extension_penalty: int, s: str, t: str) -> Tuple[int, str, str]:
    backtrack, max_score = AffineGaps(match_reward, mismatch_penalty, gap_opening_penalty, gap_extension_penalty, s, t)
    result = GetAlignment(s, t, backtrack, max_score)
    return result 

def AffineGaps(match_reward, mismatch_penalty, gap_opening_penalty, gap_extension_penalty, s, t): 
    len_s = len(s)
    len_t = len(t)

    sigma = gap_opening_penalty
    epsilon= gap_extension_penalty
    
    # init matrices 
    middle = np.zeros(shape = (len_s+1,len_t+1))
    upper = np.zeros(shape = (len_s+1,len_t+1))
    upper.fill(-np.inf)
    lower = np.zeros(shape = (len_s+1,len_t+1))
    lower.fill(-np.inf)

    backtrack = np.full(shape=(len_s + 1, len_t + 1, 3), fill_value='', dtype='<U100') 

    # init backtrack 
    for i in range(1, len_s+1):
        for j in range(1, len_t+1): 
            backtrack[i][j][0]= "continue_vertical"
            backtrack[i][j][2]= "continue_horizontal"
            # print(backtrack)

    # init gaps 
    for i in range(1, len_s+1):
        middle[i][0] = -(sigma+ epsilon*(i-1))
        backtrack[i][0][1]= "close_vertical"
    for j in range(1, len_t+1):
        middle[0][j]= -(sigma+ epsilon*(j-1))
        backtrack[0][j][1]= "close_horizontal"
    
    # get scores and backtrack 
    for i in range(1, len_s+1): 
        for j in range(1, len_t+1): 
            if s[i-1] == t[j-1]:
                score = match_reward
            else: 
                score = -mismatch_penalty 
            lower[i][j] = max(lower[i-1][j]-epsilon, middle[i-1][j]-sigma)
            upper[i][j] = max(upper[i][j-1]-epsilon, middle[i][j-1]-sigma)
            middle[i][j] = max(lower[i][j], middle[i-1][j-1]+score, upper[i][j])
            
            # backtrack: 0-lower, 1-middle, 2-upper 
            if lower[i][j] == lower[i-1][j]-epsilon:  # continue vertical gap 
                backtrack[i][j][0] = "continue_vertical" #|
            elif lower[i][j] == middle[i-1][j]-sigma:  # open vertical gap 
                backtrack[i][j][0] = "open" # + 
            
            if upper[i][j] == upper[i][j-1]-epsilon:
                # continue horizontal gap
                backtrack[i][j][2] = "continue_horizontal" # '-'
            elif upper[i][j] == middle[i][j-1]-sigma:
                # open horizontal gap
                backtrack[i][j][2] = "open" #'+'
            
            if middle[i][j] == lower[i][j]:
                # close vertical gap
                backtrack[i][j][1] = "close_vertical" #'|'
            elif middle[i][j] == upper[i][j]:
                # close horizontal gap
                backtrack[i][j][1] = "close_horizontal" #'-'
            elif middle[i][j] == middle[i-1][j-1] + score:
                backtrack[i][j][1] = "match" if s[i-1] == t[j-1] else "mismatch" #/, *
            else: 
                break
    
    max_score = middle[len_s][len_t]
    # print(middle)
    # print(upper)
    # print(lower)
    return backtrack, max_score

def GetAlignment(s, t, backtrack, max_score): 
    s_aligned = []
    t_aligned = []

    i = len(s)
    j = len(t)
    level = 1

    # backtrack: 0-lower, 1-middle, 2-upper 
    while i > 0 and j > 0: 
        # print(i, j)
        # print(level)
        # print("print(backtrack[i][j][level])", backtrack[i][j][level])
        if level == 1: # middle 
            if backtrack[i][j][level] == "close_vertical":
                level = 0
            elif backtrack[i][j][level] == "close_horizontal":
                level = 2
            elif backtrack[i][j][level]== "match" or backtrack[i][j][level]== "mismatch":
                s_aligned.append(s[i-1])
                t_aligned.append(t[j-1])
                level=1
                i -= 1
                j -= 1
        elif level == 0:
            t_aligned.append('-')
            s_aligned.append(s[i-1])
            if backtrack[i][j][level]== "continue_vertical": 
                level = 0
            if backtrack[i][j][level] == 'open':
                level = 1
            i -= 1
        elif level == 2:
            t_aligned.append(t[j-1])
            s_aligned.append('-')
            if backtrack[i][j][level]== "continue_horizontal": 
                level = 2
            if backtrack[i][j][level] == 'open':
                level = 1
            j -= 1
    
    # print(i)
    # print(j)

    if i>0: 
        while i>0:
            s_aligned.append(s[i-1])
            t_aligned.append('-')
            i-=1

    if j>0: 
        while j>0:
            s_aligned.append('-')
            t_aligned.append(t[j-1])
            j-=1
    
    s_aligned.reverse()
    t_aligned.reverse()
    s_aligned = ''.join(s_aligned)
    t_aligned = ''.join(t_aligned)
    return max_score, s_aligned, t_aligned

###################################
# 5.5- Multiple Sequence Alignment 
import sys
from typing import List, Dict, Iterable, Tuple

# Insert your MultipleAlignment function here, along with any subroutines you need
def MultipleAlignment(s1: str, s2: str, s3: str) -> Tuple[int, str, str, str]:
    #pass
    # get score and backtrack matrices 
    score, backtrack = GetMatrices(s1, s2, s3)
    # get max, s1_aligned, s2_aligned, s3_aligned 
    max_score, s1_aligned, s2_aligned, s3_aligned= BacktrackMultiple(score, backtrack, s1, s2, s3)
    return max_score, s1_aligned, s2_aligned, s3_aligned

def GetMatrices(s1, s2, s3):
    # init score and backtrack mat to size of s1 x s2 x s3
    score = [[[0 for k in range(len(s3)+1)] for j in range(len(s2)+1)] for i in range(len(s1)+1)]
    backtrack = [[[0 for k in range(len(s3)+1)] for j in range(len(s2)+1)] for i in range(len(s1)+1)]

    # three nested for loops: 1 for each seq 
    for i in range(1, len(s1)+1):
        for j in range(1, len(s2)+1):
            for k in range(1, len(s3)+1):
                # if the three match, we have a score of 1 
                if s1[i-1] == s2[j-1] == s3[k-1]: 
                    s_three_match = 1 
                # else, we have a score of 0 
                else: 
                    s_three_match = 0 
                # get the max score from all possible scores 
                scores = max(score[i-1][j][k], 
                          score[i][j-1][k], 
                          score[i][j][k-1], 
                          score[i-1][j-1][k], 
                          score[i-1][j][k-1],
                          score[i][j-1][k-1],
                          score[i-1][j-1][k-1] + s_three_match)
                # set score at this index to the max 
                score[i][j][k]= scores
                
                # set backtrack value based on the max score 
                # char from s1 added, - in 2 and 3 
                if score[i][j][k]== score[i-1][j][k]: 
                    backtrack[i][j][k]= 0
                # char from s2 added, - in 1 and 3 
                elif score[i][j][k]== score[i][j-1][k]: 
                    backtrack[i][j][k]= 1
                # char from s3 added, - in 1 and 2 
                elif score[i][j][k]== score[i][j][k-1]: 
                    backtrack[i][j][k]= 2
                # char from s1 and s3 added, - in 2  
                elif score[i][j][k]== score[i-1][j][k-1]: 
                    backtrack[i][j][k]= 3
                # char from s2 and s3 added, - in 1  
                elif score[i][j][k]== score[i][j-1][k-1]: 
                    backtrack[i][j][k]= 4
                # char from s1 and s2 added, - in 3  
                elif score[i][j][k]== score[i][j-1][k-1]: 
                    backtrack[i][j][k]= 5
                else: 
                    backtrack[i][j][k]= 6
                
                # backtrack[i][j][k] = max_score_ind
                # score[i][j][k] = max_score

    for s in score: 
        print(s, '\n')       
    return score, backtrack 

def BacktrackMultiple(score, backtrack, s1, s2, s3): 
    s1_aligned= s1 
    s2_aligned= s2 
    s3_aligned = s3

    len_s1 = len(s1)
    len_s2= len(s2)
    len_s3 = len(s3)

    max_score = score[len_s1][len_s2][len_s3]

    i = len_s1
    j = len_s2
    k = len_s3 

    # while i*j*k != 0:
    while i>0 and j>0 and k>0:
        print(i,j,k)
        if backtrack[i][j][k] == 0: # char from s1 added, - in 1 and 2
            s2_aligned = s2_aligned[:j] + '-' + s2_aligned[j:]
            s3_aligned = s3_aligned[:k] + '-' + s3_aligned[k:]
            i -= 1
        elif backtrack[i][j][k] == 1: # char from s2 added, - in 1 and 3 
            s1_aligned = s1_aligned[:i] + '-' + s1_aligned[i:]
            s3_aligned = s3_aligned[:k] + '-' + s3_aligned[k:]
            j -= 1
        elif backtrack[i][j][k] == 2: # char from s3 added, - in 1 and 2 
            s1_aligned = s1_aligned[:i] + '-' + s1_aligned[i:]
            s2_aligned = s2_aligned[:j] + '-' + s2_aligned[j:]
            k -= 1
        elif backtrack[i][j][k] == 3: # char from s1 and s3 added, - in 2  
            s2_aligned = s2_aligned[:j] + '-' + s2_aligned[j:]
            i -= 1
            k -= 1
        elif backtrack[i][j][k] == 4: # char from s2 and s3 added, - in 1  
            s1_aligned = s1_aligned[:i] + '-' + s1_aligned[i:]
            j -= 1
            k -= 1
        elif backtrack[i][j][k] == 5: # char from s1 and s2 added, - in 3  
            s3_aligned = s3_aligned[:k] + '-' + s3_aligned[k:]
            i -= 1
            j -= 1
        else: # char from s1 and s2 and s3 added 
            i -= 1
            j -= 1
            k -= 1

    # get len of  max str and add gaps to pad rest of str 
    max_length = max(len(s1_aligned), len(s2_aligned), len(s3_aligned))
    alignments = [s1_aligned, s2_aligned, s3_aligned]

    for i in range(3):
        while len(alignments[i]) < max_length:
            alignments[i] = '-' + alignments[i]

    s1_aligned, s2_aligned, s3_aligned = alignments
    
    return max_score, s1_aligned, s2_aligned, s3_aligned