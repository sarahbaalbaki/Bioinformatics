# Sarah Baalbaki 
# Week 3- Code 

AminoAcidMass = {
    'A': 71,
    'R': 156,
    'N': 114,
    'D': 115,
    'C': 103,
    'Q': 128,
    'E': 129,
    'G': 57,
    'H': 137,
    'I': 113,
    'L': 113,
    'K': 128,
    'M': 131,
    'F': 147,
    'P': 97,
    'S': 87,
    'T': 101,
    'W': 186,
    'Y': 163,
    'V': 99
}

###############################
# Exercise 3.4.2- Cyclospectrum: 
def Cyclospectrum(sequence):
    spectrum = [0]
    n = len(sequence) 
    # get prefix masses
    prefix_masses = [0]
    for char in sequence:
        prefix_masses.append(prefix_masses[-1] + AminoAcidMass[char])
    # get cyclic masses
    for i in range(n):
        for j in range(i+1, n+1):
            spectrum.append(prefix_masses[j] - prefix_masses[i])

    return sorted(spectrum)

Cyclospectrum("AG")

###############################
# Exercise 3.6- Cyclopeptide Sequencing 
from itertools import chain
import numpy as np

# takes in list of ints spectrum 
def CyclopeptideSequencing(Spectrum): 
    candidate_peptides = set()  
    candidate_peptides.add('')
    final_peptides = []
    iter = 1
    # i = 0
    
    while len(candidate_peptides) != 0: 
        print(iter)
        iter += 1
        candidate_peptides = Expand(candidate_peptides, Spectrum)
        # pep to remove 
        peptides_to_remove = set()  
        
        for peptide in candidate_peptides: 
            cyclospec = Cyclospectrum(peptide)
            mass_peptide = cyclospec[-1] 
            parent_mass = max(Spectrum)
            
            if mass_peptide == parent_mass: 
                if Cyclospectrum(peptide) == Spectrum and peptide not in final_peptides: 
                    final_peptides.append(peptide)
                peptides_to_remove.add(peptide)  
            # elif Cyclospectrum(peptide) != Spectrum: 
            elif mass_peptide not in Spectrum: 
                peptides_to_remove.add(peptide)  

        candidate_peptides.difference_update(peptides_to_remove)
        #print(len(peptides_to_remove))
        
        #print("candidate_peptides_updated", candidate_peptides)
        #print("final_peptides", final_peptides)
        
        # i += 1
        #print("HERE IS A NEW ONE")
    
    result = GetNumbers(final_peptides)
    return result 

def GetNumbers(peptides): 
    # get the output in the desired format - between each base mass 
    result= []
    for peptide in peptides: 
        number_seq=[]
        for char in peptide: 
            number_seq.append(AminoAcidMass[char])
        if number_seq not in result: 
            result.append(number_seq)
    
    result= sorted(result, reverse= True)
    result_str_format=[]
    for l in result: 
        result_str_format.append('-'.join(map(str, l)))

    return result_str_format

def Expand(candidate_peptides, Spectrum):
    expanded= []
    for pep in candidate_peptides: 
        for key in AminoAcidMass: 
            if AminoAcidMass[key] in Spectrum:
                expanded.append(pep+key)
    print("expanded set", expanded)
    return set(expanded)

def Cyclospectrum(sequence): 
    kmers= []
    for i in range(0, len(sequence)+1): 
        # print(i)
        resulting= build_kmers(sequence, i)
        # print(resulting)
        kmers.append(resulting)
    
    kmers= list(chain.from_iterable(kmers))

    spectrum= []
    for str in kmers: 
        if str==['']: 
            spectrum.append(0)
        else: 
            sum=0
            for char in str: 
                sum+= AminoAcidMass[char]
            spectrum.append(sum)
    return sorted(spectrum)

def build_kmers(sequence, ksize):
    kmers = []
    n_kmers = len(sequence) 
    if ksize==0: 
        return [""]
    if ksize== n_kmers: 
        return [sequence]
    for i in range(n_kmers):
        # gen kmers
        kmer = sequence[i:i + ksize] + sequence[:max(0, ksize - (n_kmers - i))]
        kmers.append(kmer)

    return kmers

# CyclopeptideSequencing([0, 113])
result= CyclopeptideSequencing([0, 113, 128, 186, 241, 299, 314, 427])
for r in result: 
    print(r, end= " ")

###############################
# Exercise 3.7.8- Leaderboard 
def Trim(Leaderboard, Spectrum, N): 
    linear_scores = []
    for j in range(0, len(Leaderboard)): 
        peptide = Leaderboard[j]
        #print(peptide)
        #print(Spectrum)
        theo_spec = theoretical_spectrum_linear(peptide)
        spec = Spectrum
        #print(get_score_linear(theo_spec, spec))
        linear_scores.append(get_score_linear(theo_spec, Spectrum))
        
    #print(linear_scores)
    combined = [(score, peptide) for score, peptide in zip(linear_scores, Leaderboard)]
    sorted_combined = sorted(combined, key=lambda x: x[0], reverse=True)
    print(sorted_combined)
    sorted_linear_scores = [entry[0] for entry in sorted_combined]
    #print(sorted_linear_scores)
    sorted_leaderboard = [entry[1] for entry in sorted_combined]
    #print(sorted_leaderboard)

    top_N = np.unique(sorted_linear_scores)
    top_N = sorted(top_N, reverse= True)[:N]
    print(top_N)
    top_peptides = []
    for s in top_N: 
        for pep in range(len(sorted_leaderboard)): 
            if sorted_linear_scores[pep]== s: 
                top_peptides.append(sorted_leaderboard[pep])
    return top_peptides

def theoretical_spectrum_linear(sequence):
    kmers = []
    for i in range(0, len(sequence) + 1): 
        resulting = build_kmers_linear(sequence, i)
        kmers.extend(resulting)
    
    spectrum = []
    for kmer in kmers: 
        if kmer == '': 
            spectrum.append(0)
        else: 
            mass_sum = 0
            for char in kmer: 
                mass_sum += AminoAcidMass[char]  # Assuming AminoAcidMass is defined elsewhere
            spectrum.append(mass_sum)
    return sorted(spectrum)

def build_kmers_linear(sequence, ksize):
    n_kmers = len(sequence)  

    if ksize == 0: 
        return [""]
    
    if ksize == n_kmers: 
        return [sequence]
    
    for i in range(n_kmers):
#         #kmer = sequence[i:i + ksize] + sequence[:max(0, ksize - (n_kmers - i))]
          kmers = [sequence[i:i+ksize] for i in range(len(sequence) - ksize + 1)]
    return kmers

def LeaderboardCyclopeptideSequencing(Spectrum, N):
    leaderboard = set()
    leaderboard.add('')
    leader_peptide = ""
    max_score = 0
    
    parent_mass = max(Spectrum)
    while len(leaderboard)>0:
        leaderboard = Expand_1(leaderboard, Spectrum)
        peptides_to_remove = []
        for peptide in leaderboard:
            cyclospec = theoretical_spectrum_linear(peptide)
            #cyclospec = Cyclospectrum(peptide)
            mass_peptide = cyclospec[-1]
            if mass_peptide == parent_mass:
                score = get_score_cyclo(peptide, Spectrum)
                #score = get_score_linear(peptide, Spectrum)
                if score > max_score:
                    leader_peptide = peptide
                    max_score = score
                # elif score == max_score: 
                #     peptides_to_remove.append(peptide)
            elif mass_peptide > parent_mass:
                peptides_to_remove.append(peptide)
        #leaderboard = [peptide for peptide in leaderboard if peptide not in peptides_to_remove]
        for p in peptides_to_remove:
            leaderboard.remove(p)
        leaderboard = Trim(leaderboard, Spectrum, N)
    return GetNumbers([leader_peptide])

def get_score_cyclo(peptide, Spectrum): 
    score = 0 
    cycloscore = Cyclospectrum(peptide)
    experimental_spectrum_copy = copy.deepcopy(Spectrum)
    for s in cycloscore: 
        if s in experimental_spectrum_copy: 
            score += 1 
            experimental_spectrum_copy.remove(s)  
    return score

def Cyclospectrum(sequence):
    spectrum = [0]
    n = len(sequence)
    
    # Calculate prefix masses
    prefix_masses = [0]
    for char in sequence:
        prefix_masses.append(prefix_masses[-1] + AminoAcidMass[char])
    
    # Calculate cyclic masses
    for i in range(n):
        for j in range(i+1, n+1):
            spectrum.append(prefix_masses[j] - prefix_masses[i])

    return sorted(spectrum)

def Expand_1(candidate_peptides, Spectrum):
    expanded= []
    for pep in candidate_peptides: 
        for key in AminoAcidMass: 
            if AminoAcidMass[key] in Spectrum:
                expanded.append(pep+key)
    print("expanded set", expanded)
    return set(expanded)

def get_score_linear(theoretical_spectrum, experimental_spectrum): 
    score = 0
    experimental_spectrum_copy = copy.deepcopy(experimental_spectrum)
    for s in theoretical_spectrum: 
        if s in experimental_spectrum_copy: 
            score += 1 
            experimental_spectrum_copy.remove(s)  # This might not be the desired behavior
    return score

def Trim(Leaderboard, Spectrum, N):
    theoretical_spectra = [theoretical_spectrum_linear(peptide) for peptide in Leaderboard]
    linear_scores = [(get_score_linear(theo_spec, Spectrum), peptide) for theo_spec, peptide in zip(theoretical_spectra, Leaderboard)]
    sorted_combined = sorted(linear_scores, key=lambda x: x[0], reverse=True)
    print(sorted_combined)
    top_N_scores = sorted(set(score for score, _ in sorted_combined), reverse=True)[:N]
    
    top_peptides = [peptide for score in top_N_scores for score_, peptide in sorted_combined if score_ == score]
    return top_peptides

def GetNumbers(peptides):
    result = []
    for peptide in peptides:
        number_seq = [AminoAcidMass[char] for char in peptide]
        if number_seq not in result:
            result.append(number_seq)

    result = sorted(result, reverse=True)

    result_str_format = ['-'.join(map(str, l)) for l in result]
    return result_str_format

def add_commas_between_numbers(numbers_str):
    number_strings = numbers_str.split()
    integers_list = [int(num_str) for num_str in number_strings]

    return integers_list

int_in= 10
input_in= add_commas_between_numbers("0 71 113 129 147 200 218 260 313 331 347 389 460")
#input_in= add_commas_between_numbers("0 97 99 113 114 115 128 128 147 147 163 186 227 241 242 244 244 256 260 261 262 283 291 309 330 333 340 347 385 388 389 390 390 405 435 447 485 487 503 504 518 544 552 575 577 584 599 608 631 632 650 651 653 672 690 691 717 738 745 770 779 804 818 819 827 835 837 875 892 892 917 932 932 933 934 965 982 989 1039 1060 1062 1078 1080 1081 1095 1136 1159 1175 1175 1194 1194 1208 1209 1223 1322")
print(LeaderboardCyclopeptideSequencing(input_in, int_in))

import copy 
AminoAcidMass = [57, 71, 87, 97, 99, 101, 103, 113, 113, 114, 115, 128, 128, 129, 131, 137, 147, 156, 163, 186]

def GetScore(theoretical_spectrum, experimental_spectrum):
    # theo spec set to remove el
    theoretical_spectrum_set = set(theoretical_spectrum)
    score = 0
    # scorre based on common masses 
    for mass in experimental_spectrum:
        if mass in theoretical_spectrum_set:
            theoretical_spectrum_set.remove(mass)
            score += 1
    return score

def TheoreticalSpectrum_Leaderboard(peptide):
    l = len(peptide)
    # loop fopr cyclic 
    cyclo = peptide + peptide
    specs = [0, sum(peptide)]
    for start in range(l):
        for length in range(1, l):
            mass = sum(cyclo[start:start+length])
            specs.append(mass)
    # print(subspectrums)
    return specs

def Trim_Leaderboard(Leaderboard, Spectrum, N):
    if len(Leaderboard) < N:
        return Leaderboard
    else:
        results = []
        for pep in Leaderboard:
            exp_spec = TheoreticalSpectrum_Leaderboard(pep)
            score= GetScore(Spectrum, exp_spec)
            results.append((pep, score))
        # sort res 
        sorted_results = sorted(results, key=lambda x: x[1], reverse=True)
        # get top N pep scores 
        top_N = sorted_results[N-1][1]
        # select top N scores with repeats 
        results = [res for res in sorted_results if res[1] >= top_N]
        # get the pep
        top_N_peptides = [peptide for peptide, _ in results]
        return top_N_peptides

def MakeDict(Spectrum):
    d= {} 
    for s in Spectrum: 
        if s in d: 
            d[s]+=1
        else: 
            d[s]=1
    return d 

def Expand_Leaderboard(leaderboard, Spectrum): 
    expanded_leaderboard = []
    for pep in leaderboard:
        for mass in AminoAcidMass:
            # if mass in Spectrum: 
            # new_pep= (pep, mass)
            # expanded_leaderboard.append(list(new_pep))
            dict_masses= MakeDict(Spectrum)
            # print(dict_masses)
            for char in pep: 
                # print(char)
                dict_masses[char]-=1
            # print(dict_masses)
            if mass in dict_masses and dict_masses[mass] >0: 
                pep_copy = copy.deepcopy(pep)
                pep_copy.append(mass)
                expanded_leaderboard.append(pep_copy)
                dict_masses[mass]-=1
    # print(expanded_leaderboard)
    return expanded_leaderboard

def UpdateLeaderboard(leaderboard, pep): 
    new_leaderboard = []
    flattened_pep = pep[0] + [pep[1]]
    for p in leaderboard:
        if p == pep:
        # if x is equal to pep, append a modified version of x to the new leaderboard
            new_leaderboard.append(flattened_pep)
        else:
            # if x is not equal to pep, append x unchanged to the new leaderboard
            new_leaderboard.append(p)
    return new_leaderboard, flattened_pep

def LeaderboardCyclopeptideSequencing_2(Spectrum,N):
    leaderboard = [[0]]
    leader_peptide = []
    while leaderboard != []:
        leaderboard = Expand_Leaderboard(leaderboard, Spectrum)
        # print(leaderboard)
        for pep in leaderboard:
            # print(pep)
            pep_spec = TheoreticalSpectrum_Leaderboard(pep)
            mass_pep= max(pep_spec)
            parentMass= max(Spectrum)
            if mass_pep == parentMass:
                leader_spec = TheoreticalSpectrum_Leaderboard(leader_peptide)
                if GetScore(Spectrum,pep_spec) > GetScore(Spectrum,leader_spec):
                    leader_peptide = pep
                leaderboard.remove(pep)
            elif mass_pep > parentMass:
                leaderboard.remove(pep)
        leaderboard = Trim_Leaderboard(leaderboard,Spectrum,N)
    return leader_peptide[1:]

def add_commas_between_numbers(numbers_str):
    number_strings = numbers_str.split()
    integers_list = [int(num_str) for num_str in number_strings]

    return integers_list

N=10
Spectrum = add_commas_between_numbers("0 71 113 129 147 200 218 260 313 331 347 389 460")
results = LeaderboardCyclopeptideSequencing_2(Spectrum,N)
print(results)

inp_num = 325 
val_in = add_commas_between_numbers("0 71 71 71 87 97 97 99 101 103 113 113 114 115 128 128 129 137 147 163 163 170 184 184 186 186 190 211 215 226 226 229 231 238 241 244 246 257 257 276 277 278 299 300 312 316 317 318 318 323 328 340 343 344 347 349 356 366 370 373 374 391 401 414 414 415 419 427 427 431 437 441 446 453 462 462 462 470 472 502 503 503 511 515 529 530 533 533 540 543 547 556 559 569 574 575 584 590 600 600 604 612 616 617 630 640 640 643 646 648 660 671 683 684 687 693 703 703 719 719 719 729 730 731 737 740 741 745 747 754 774 780 784 790 797 800 806 818 826 827 832 833 838 846 846 847 850 868 869 877 884 889 893 897 903 908 913 917 930 940 947 956 960 960 961 964 965 966 983 983 985 1002 1009 1010 1011 1021 1031 1031 1036 1053 1054 1058 1059 1062 1063 1074 1076 1084 1092 1103 1113 1122 1124 1130 1133 1134 1145 1146 1146 1149 1150 1155 1156 1171 1173 1174 1187 1191 1193 1200 1212 1221 1233 1240 1242 1246 1259 1260 1262 1277 1278 1283 1284 1287 1287 1288 1299 1300 1303 1309 1311 1320 1330 1341 1349 1357 1359 1370 1371 1374 1375 1379 1380 1397 1402 1402 1412 1422 1423 1424 1431 1448 1450 1450 1467 1468 1469 1472 1473 1473 1477 1486 1493 1503 1516 1520 1525 1530 1536 1540 1544 1549 1556 1564 1565 1583 1586 1587 1587 1595 1600 1601 1606 1607 1615 1627 1633 1636 1643 1649 1653 1659 1679 1686 1688 1692 1693 1696 1702 1703 1704 1714 1714 1714 1730 1730 1740 1746 1749 1750 1762 1773 1785 1787 1790 1793 1793 1803 1816 1817 1821 1829 1833 1833 1843 1849 1858 1859 1864 1877 1886 1890 1893 1900 1900 1903 1904 1918 1922 1930 1930 1931 1961 1963 1971 1971 1971 1980 1987 1992 1996 2002 2006 2006 2014 2018 2019 2019 2032 2042 2059 2060 2063 2067 2077 2084 2086 2089 2090 2093 2105 2110 2115 2115 2116 2117 2121 2133 2134 2155 2156 2157 2176 2176 2187 2189 2192 2195 2202 2204 2207 2207 2218 2222 2243 2247 2247 2249 2249 2263 2270 2270 2286 2296 2304 2305 2305 2318 2319 2320 2320 2330 2332 2334 2336 2336 2346 2362 2362 2362 2433")
ans= LeaderboardCyclopeptideSequencing_2(val_in, inp_num)
print(ans)

################################
# 3.9.1- Spectral Convolution: 
def spectral_convolution(int_list): 
    array1= copy.deepcopy(int_list)
    array2= copy.deepcopy(int_list)
    # array2[1]=0
    matrix = []
    for i in array1: 
        to_add= []
        for j in array2: 
            to_add.append(i-j)
        matrix.append(to_add)
    print(matrix)
    flattened= [item for sublist in matrix for item in sublist]
    flattened_abs= []
    for item in flattened: 
        if item>0: 
            flattened_abs.append(item)
    occurrences = {}
    for element in flattened_abs:
        if element in occurrences:
            occurrences[element] += 1
        else:
            occurrences[element] = 1
    print(flattened_abs)
    print(occurrences)
    # occurrences.pop(0)
    sorted_occurrences = sorted(occurrences.items(), key=lambda x: x[1], reverse= True)
    print(sorted_occurrences)
    results = []
    for key in sorted_occurrences: 
        for number in range(key[1]):
            results.append(key[0])
    return results

# 3.9.4- ConvolutionCyclopeptideSequencing: 
def ConvolutionCyclopeptideSequencing(M, N, Spectrum): 
    convolution_exp = spectral_convolution(Spectrum)
    print("convolution_exp", convolution_exp)
    filtered_convolution = [x for x in convolution_exp if 57 < x < 200]
    print("filtered_convolution", filtered_convolution)
    # Count frequencies
    frequency_counter = Counter(filtered_convolution)
    
    # Sort frequencies with ties
    sorted_frequencies = sorted(frequency_counter.items(), key=lambda x: (-x[1], x[0]))
    
    # Select top M elements with ties
    selected_elements = []
    max_frequency = sorted_frequencies[0][1]
    for element, frequency in sorted_frequencies:
        if M <= 0 or frequency < max_frequency:
            break
        selected_elements.append(element)
        M -= 1
    
    # Extract corresponding masses
    # aa_options=[]
    candidate_masses = sorted(selected_elements)
    print("candidate_masses", candidate_masses)
    # for c in candidate_masses: 
    #     for k in AminoAcidMass: 
    #         if AminoAcidMass[k]==c: 

    result = LeaderboardCyclopeptideSequencing(Spectrum, N, candidate_masses)
    print(result)
    return result

def LeaderboardCyclopeptideSequencing(Spectrum, N, AminoAcidOptions):
    leaderboard = set()
    leaderboard.add('')
    leader_peptide = ""
    max_score = 0
    
    #parent_mass = np.sum(AminoAcidOptions)
    parent_mass = max(Spectrum)
    #parent_mass = np.sum(AminoAcidOptions)
    while leaderboard:
        leaderboard = Expand_1(leaderboard, Spectrum)
        print("leaderboard", leaderboard)
        peptides_to_remove = set()
        for peptide in leaderboard.copy():
            cyclospec = theoretical_spectrum_linear(peptide)
            mass_peptide = cyclospec[-1]
            if mass_peptide == parent_mass:
                score = get_score_linear(cyclospec, AminoAcidOptions)
                if score > max_score:
                    leader_peptide = peptide
                    max_score = score
            elif mass_peptide > parent_mass:
                peptides_to_remove.add(peptide)
        #leaderboard = [peptide for peptide in leaderboard if peptide not in peptides_to_remove]
        leaderboard.difference_update(peptides_to_remove)
        leaderboard = Trim(leaderboard, Spectrum, N)
    return GetNumbers([leader_peptide])

def spectral_convolution(int_list): 
    array1= copy.deepcopy(int_list)
    array2= copy.deepcopy(int_list)
    # array2[1]=0
    matrix = []
    for i in array1: 
        to_add= []
        for j in array2: 
            to_add.append(i-j)
        matrix.append(to_add)
    print(matrix)
    flattened= [item for sublist in matrix for item in sublist]
    flattened_abs= []
    for item in flattened: 
        if item>0: 
            flattened_abs.append(item)
    occurrences = {}
    for element in flattened_abs:
        if element in occurrences:
            occurrences[element] += 1
        else:
            occurrences[element] = 1
    print(flattened_abs)
    print(occurrences)
    # occurrences.pop(0)
    sorted_occurrences = sorted(occurrences.items(), key=lambda x: x[1], reverse= True)
    print(sorted_occurrences)
    results = []
    for key in sorted_occurrences: 
        for number in range(key[1]):
            results.append(key[0])
    return results

M=20
N= 10
Spectrum= add_commas_between_numbers("57 57 71 99 129 137 170 186 194 208 228 265 285 299 307 323 356 364 394 422 493")
#Spectrum= add_commas_between_numbers("853 113 585 796 924 502 423 1210 342 186 761 391 593 1412 1152 1396 260 129 1381 229 242 356 990 1047 57 748 1176 730 990 1038 1119 294 339 114 696 1251 1267 617 567 357 471 163 1266 1281 0 536 1395 454 1104 1362 1039 892 1509 1086 129 649 1095 713 258 777 1394 753 299 599 648 876 414 1249 813 242 859 1305 552 1284 861 650 1249 261 520 470 519 957 1233 405 260 861 762 810 1248 891 916 1346 390 981 147 1323 390 732 618 1380 1038 756 989 225 633 910 204 1452 243 1119 860 1395 129 57 503 1267 1153 276 462 228 1215 114 1170 357 973 388 519 699 131 128 1120 648 1452 1055 632 333 1380 528 747 389 656 97 1167 779 1380 1280 942 115 1121 1152 1007 990 1006 1118 519 877 1378 471")
answer= ConvolutionCyclopeptideSequencing(M, N, Spectrum)
print("answer", answer)

# this took a long time to run for the test cases, re-wrote the code based on Leaderboard_2 for using integer representation 

from collections import Counter

def spectral_convolution(int_list): 
    array1= copy.deepcopy(int_list)
    array2= copy.deepcopy(int_list)
    # array2[1]=0
    matrix = []
    for i in array1: 
        to_add= []
        for j in array2: 
            to_add.append(i-j)
        matrix.append(to_add)
    # print(matrix)
    flattened= [item for sublist in matrix for item in sublist]
    flattened_abs= []
    for item in flattened: 
        if item>0: 
            flattened_abs.append(item)
    occurrences = {}
    for element in flattened_abs:
        if element in occurrences:
            occurrences[element] += 1
        else:
            occurrences[element] = 1
    # print(flattened_abs)
    # print(occurrences)
    # occurrences.pop(0)
    sorted_occurrences = sorted(occurrences.items(), key=lambda x: x[1], reverse= True)
    # print(sorted_occurrences)
    results = []
    for key in sorted_occurrences: 
        for number in range(key[1]):
            results.append(key[0])
    return results

def ConvolutionCyclopeptideSequencing(M, N, Spectrum): 
    convolution_exp = spectral_convolution(Spectrum)
    # print("convolution_exp", convolution_exp)
    filtered_convolution = [x for x in convolution_exp if 57 < x < 200]
    # print("filtered_convolution", filtered_convolution)
    # Count frequencies
    frequency_counter = Counter(filtered_convolution)
    # sort frequencies with ties
    sorted_frequencies = sorted(frequency_counter.items(), key=lambda x: (-x[1], x[0]))
    # get top M elements with ties
    selected_elements = []
    max_frequency = sorted_frequencies[0][1]
    for element, frequency in sorted_frequencies:
        if M <= 0:
            break
        selected_elements.append(element)
        M -= 1
    
    # get corresponding masses
    # aa_options=[]
    candidate_masses = sorted(selected_elements)
    print("candidate_masses", candidate_masses)
    # for c in candidate_masses: 
    #     for k in AminoAcidMass: 
    #         if AminoAcidMass[k]==c: 

    result = LeaderboardCyclopeptideSequencing_2_Convolution(Spectrum, N, candidate_masses)
    # print(result)
    return result

def LeaderboardCyclopeptideSequencing_2_Convolution(Spectrum, N, AminoAcidOptions):
    leaderboard = [[0]]
    leader_peptide = []
    while leaderboard != []:
        leaderboard = Expand_Leaderboard_Convolution(leaderboard, Spectrum, AminoAcidOptions)
        # print(leaderboard)
        for pep in leaderboard:
            # print(pep)
            pep_spec = TheoreticalSpectrum_Leaderboard(pep)
            mass_pep= max(pep_spec)
            parentMass= max(Spectrum)
            if mass_pep == parentMass:
                leader_spec = TheoreticalSpectrum_Leaderboard(leader_peptide)
                if GetScore(Spectrum,pep_spec) > GetScore(Spectrum,leader_spec):
                    leader_peptide = pep
                leaderboard.remove(pep)
            elif mass_pep > parentMass:
                leaderboard.remove(pep)
        leaderboard = Trim_Leaderboard(leaderboard,Spectrum,N)
    return leader_peptide[1:]

def Expand_Leaderboard_Convolution(leaderboard, Spectrum, AminoAcidOptions): 
    expanded_leaderboard = []
    Spectrum.append(0)
    print(Spectrum)
    for pep in leaderboard:
        for mass in AminoAcidOptions:
            # if mass in Spectrum: 
            # new_pep= (pep, mass)
            # expanded_leaderboard.append(list(new_pep))
            dict_masses= MakeDict(Spectrum)
            # print(dict_masses)
            for char in pep: 
                #print("pep", pep, "char", char)
                # print(char)
                dict_masses[char]-=1
            # print(dict_masses)
            if mass in dict_masses and dict_masses[mass] >0: 
                pep_copy = copy.deepcopy(pep)
                pep_copy.append(mass)
                expanded_leaderboard.append(pep_copy)
                dict_masses[mass]-=1
    # print(expanded_leaderboard)
    return expanded_leaderboard

print(spectral_convolution(add_commas_between_numbers("0 137 186 323")))
M=20
# N= 60
N=373
#Spectrum= add_commas_between_numbers("57 57 71 99 129 137 170 186 194 208 228 265 285 299 307 323 356 364 394 422 493")
Spectrum= add_commas_between_numbers("853 113 585 796 924 502 423 1210 342 186 761 391 593 1412 1152 1396 260 129 1381 229 242 356 990 1047 57 748 1176 730 990 1038 1119 294 339 114 696 1251 1267 617 567 357 471 163 1266 1281 0 536 1395 454 1104 1362 1039 892 1509 1086 129 649 1095 713 258 777 1394 753 299 599 648 876 414 1249 813 242 859 1305 552 1284 861 650 1249 261 520 470 519 957 1233 405 260 861 762 810 1248 891 916 1346 390 981 147 1323 390 732 618 1380 1038 756 989 225 633 910 204 1452 243 1119 860 1395 129 57 503 1267 1153 276 462 228 1215 114 1170 357 973 388 519 699 131 128 1120 648 1452 1055 632 333 1380 528 747 389 656 97 1167 779 1380 1280 942 115 1121 1152 1007 990 1006 1118 519 877 1378 471")
answer= ConvolutionCyclopeptideSequencing(M, N, Spectrum)
print("answer", answer)