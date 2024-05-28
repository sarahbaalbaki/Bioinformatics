# Sarah Baalbaki 
# Bioinformatics- Week 1 Code 

####################
# Problem 1: 
def pattern_count(text: str, pattern: str) -> int:
    count = 0
    for i in range(len(text)-len(pattern)+1):
        if text[i:i+len(pattern)] == pattern:
            count += 1
    return count
####################

####################
# Problem 2: 
def profile_most_probable_kmer(text: str, k: int,
                               profile: list[dict[str, float]]) -> str:
    """Identifies the most probable k-mer according to a given profile matrix.

    The profile matrix is represented as a list of columns, where the i-th element is a map
    whose keys are strings ("A", "C", "G", and "T") and whose values represent the probability
    associated with this symbol in the i-th column of the profile matrix.
    """
    max_prob= -1
    profile_max_prob_kmer = ""
    for i in range(len(text)-k+1):
        kmer= text[i:i+k]
        prob= 1.0
        for j in range(k):
            val = kmer[j]
            prob = prob * profile[j][val]
        if prob> max_prob: 
            max_prob= prob
            profile_max_prob_kmer = kmer 
    return profile_max_prob_kmer
####################

####################
# Problem 3: 
import sys
import random 
import math 
from functools import reduce
      
def Motifs(profile, dna):
    motifs = []
    # for each text in dna, get kmers 
    for text in dna:
        # get the most prob kmer
        k = len(profile[list(profile.keys())[0]])
        most_prob_kmer = ProfileMostProbableKmer(text, k, profile)
        motifs.append(most_prob_kmer)
    # return motifs 
    return motifs

def randomized_motif_search(dna, k, t):
    obtained =[]
    # perform this loop 1000 times 
    for i in range(1000): 
        # get random motifs 
        random_motifs = RandomMotifs(dna, k, t)
        # run until optimum is found 
        obtained.append(RunUntilOptimum(random_motifs, dna))
    # return best scoring collection of motifs from set of all runs 
    return BestScoringMotif(obtained)

def BestScoringMotif(obtained_set): 
    # get best scoring set of motifs 
    best_score = math.inf 
    set_final = []
    # get score for each, return the one with the best score 
    for s in obtained_set: 
        score_of_set = Score(s)
        if score_of_set < best_score: 
            best_score = score_of_set
            set_final = s
    return set_final 

def RunUntilOptimum(current_motifs, dna):
    # get profile matrix 
    profile_matrix = ProfileMatrix(current_motifs)
    # get updated motifs 
    updated_motifs = Motifs(profile_matrix, dna)
    # if score is the same, return current 
    if Score(updated_motifs) == Score(current_motifs):
        return current_motifs
    # else, run until optimum
    return RunUntilOptimum(updated_motifs, dna)

def RandomMotifs(dna, k, t):
    # array for random motifs
    random_motifs = []
    # for each str text in dna set
    for text in dna: 
        # get random position
        kmer_start_pos = random.randint(0, len(text) - k)
        # append the kmer starting at that position
        random_motifs.append(text[kmer_start_pos: kmer_start_pos + k])
    return random_motifs

def ProfileMatrix(motifs):
    # get motifs count 
    motifs_count = Count(motifs)
    # prepare dict for profile with pseudocounts = 1 
    motifs_pseudocounts = {}
    for key, value in motifs_count.items():
        #add_pseudocount_to_array(value)
        motifs_pseudocounts[key] = [x + 1 for x in value] 
    # update the counts 
    divisor = 4 + len(motifs) 
    motifs_pseudocounts_updated = {}
    for key, value in motifs_pseudocounts.items():
        motifs_pseudocounts_updated[key] = [x / divisor for x in value]
    return motifs_pseudocounts_updated

def ProfileMostProbableKmer(text, k, profile):
    # get set of kmers for text
    kmers = [text[iterator:iterator + k] for iterator in range(len(text) - k + 1)]
    max_probability = 0
    most_probable_kmer = ""
    # for each kmer, get prob of occurence 
    for kmer in kmers:
        probability = 1
        for i in range(len(kmer)):
            probability *= profile[kmer[i]][i]       
        # if prob of kmer is more than prev, change the kmer max 
        if probability > max_probability:
            max_probability = probability
            most_probable_kmer = kmer   
    # return most probable one 
    return most_probable_kmer

def Count(motifs):
    # get results
    result = {}
    for char in "ACGT":
        result[char] = []
        for motif in motifs:
            result[char].append([1 if c == char else 0 for c in motif])
    # make dict for counts of each
    count_dict = {}
    for char in result:
        count_dict[char] = []
        for i in range(len(result[char][0])):
            count = 0
            for j in range(len(result[char])):
                count += result[char][j][i]
            count_dict[char].append(count)
    return count_dict

def Score(Motifs):
    # hamming dist 
    hamming_distance = lambda p, q: sum(int(x != y) for x, y in zip(p, q))
    # consensus
    count_matrix = {}
    for symbol in "ACGT":
        count_matrix[symbol] = [0] * len(Motifs[0])
    for i in range(len(Motifs[0])):
        for motif in Motifs:
            count_matrix[motif[i]][i] += 1
    symbols = "ACGT"
    string_length = len(Motifs[0])
    consensus_list = [''] * string_length
    for i in range(string_length):
        symbol_counts = {count_matrix[symbol][i]: symbol for symbol in symbols}
        max_count = max(symbol_counts.keys())
        consensus_list[i] = symbol_counts[max_count]
    consensus = ''.join(consensus_list)
    # Score calculation
    return sum(hamming_distance(motif, consensus) for motif in Motifs)
####################

####################
# Problem 4: 
def gibbs_sampler(dna, k, t, n):
    # get the ranom motifs and assign them to best so far to begin with 
    randomly_sel_motifs = RandomMotifs(dna, k, t)
    best_motifs = randomly_sel_motifs
    # go over all the gibbs sampling with 20 random starts, selecting the best one 
    for rounds in range(20):
        current_motifs = RandomMotifs(dna, k, t)
        for j in range(n):
            # random string to remove from profile
            index_to_remove = random.randint(0, t-1)
            profile = ProfileMatrix(current_motifs[:index_to_remove] + current_motifs[index_to_remove+1:]) # generate the profile with pseudocounts
            motif_at_i = GenerateKmer(dna[index_to_remove], profile, k) # get the motif of the removed string based on the profile 
            current_motifs = current_motifs[:index_to_remove] + [motif_at_i] + current_motifs[index_to_remove+1:] # return the motifs (old and new)
            # check if better than previously
            if Score(current_motifs) < Score(best_motifs):
                best_motifs = current_motifs
    return best_motifs

def RandomMotifs(dna, k, t):
    # array for random motifs
    random_motifs = []
    # for each str text in dna set
    for text in dna: 
        # get random position
        kmer_start_pos = random.randint(0, len(text) - k)
        # append the kmer starting at that position
        random_motifs.append(text[kmer_start_pos: kmer_start_pos + k])
    return random_motifs

def ProfileMatrix(motifs):
    # get motifs count 
    motifs_count = Count(motifs)
    # prepare dict for profile with pseudocounts = 1 
    motifs_pseudocounts = {}
    for key, value in motifs_count.items():
        #add_pseudocount_to_array(value)
        motifs_pseudocounts[key] = [x + 1 for x in value] 
    # update the counts 
    divisor = 4 + len(motifs) 
    motifs_pseudocounts_updated = {}
    for key, value in motifs_pseudocounts.items():
        motifs_pseudocounts_updated[key] = [x / divisor for x in value]
    return motifs_pseudocounts_updated

def Count(motifs):
    # get results
    result = {}
    for char in "ACGT":
        result[char] = []
        for motif in motifs:
            result[char].append([1 if c == char else 0 for c in motif])
    # make dict for counts of each
    count_dict = {}
    for char in result:
        count_dict[char] = []
        for i in range(len(result[char][0])):
            count = 0
            for j in range(len(result[char])):
                count += result[char][j][i]
            count_dict[char].append(count)
    return count_dict

def GenerateKmer(text, profile, k):
    # store probs of eack kmer 
    probabilities = {}
    for i in range(len(text) - k + 1):
        kmer = text[i:i+k]
        # get prob of generation of kmer
        prob = 1
        for i in range(len(kmer)):
            prob *= profile[kmer[i]][i]
        probabilities[kmer]= prob
    all_probs = sum(probabilities.values())
    weighted_probs = {}
    for key, value in probabilities.items():
        weighted_probs[key] = value / all_probs
    #weighted die
    key_to_range = KeyRanges(weighted_probs)
    random_fraction = random.uniform(0, 1)
    for key, value in key_to_range.items():
        if value['lower'] < random_fraction <= value['upper']:
            return key
    # return next(key for key, value in key_to_range.items() if value['lower'] < random_fraction <= value['upper'])

def KeyRanges(probabilities):
    keys = probabilities.keys()
    vals = map(lambda key: probabilities[key], keys)
    # upper bounds
    # upper = reduce(lambda result, element: result + [result[-1] + element], vals, [0])[1:]
    upper = [0]
    for val in vals:
        upper.append(upper[-1] + val)
    upper = upper[1:]
    # lower bounds
    lower = [0] + upper[:-1]
    lower_and_upper_bounds = list(map(lambda lower, upper: {'lower': lower, 'upper': upper}, lower, upper))
    # dictionary of keys with ranges 
    return dict(zip(keys, lower_and_upper_bounds))

def Score(Motifs):
    # hamming dist 
    hamming_distance = lambda p, q: sum(int(x != y) for x, y in zip(p, q))
    # consensus
    count_matrix = {}
    for symbol in "ACGT":
        count_matrix[symbol] = [0] * len(Motifs[0])
    for i in range(len(Motifs[0])):
        for motif in Motifs:
            count_matrix[motif[i]][i] += 1
    symbols = "ACGT"
    string_length = len(Motifs[0])
    consensus_list = [''] * string_length
    for i in range(string_length):
        count_to_symbols_tuple = {count_matrix[symbol][i]: symbol for symbol in symbols}
        max_count = max(count_to_symbols_tuple.keys())
        consensus_list[i] = count_to_symbols_tuple[max_count]
    consensus = ''.join(consensus_list)
    # Score calculation
    return sum(hamming_distance(motif, consensus) for motif in Motifs)
####################