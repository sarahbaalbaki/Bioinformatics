# Sarah Baalbaki 
# Bioinformatics- Module 11 

# 11.1- The Burrows Wheeler Transform: 
def burrows_wheeler_transform(text: str) -> str:
    """
    Generate the Burrows-Wheeler Transform of the given text.
    """
    # pass
    text+='$'
    text_rotations=[]
    # generate the rotations of text
    for i in range(len(text), 0, -1): 
        word= text[i:] + text[:i]
        text_rotations.append(word)
    # print(text_rotations)
    # sort them lexicographically 
    sorted_text_rotations= sorted(text_rotations)
    # print(sorted_text_rotations)
    # get end letters of sorted order of rotations 
    end_letters=[]
    for rotated_word in  sorted_text_rotations:
        end_letter= rotated_word[-1]
        end_letters.append(end_letter)
    return ''.join(end_letters)

###########################################
# 11.3- Inverse Burrows-Wheeler Transform: 
def inverse_burrows_wheeler_transform(transform: str) -> str:
    """
    Generate the inverse of the Burrows-Wheeler Transform.
    """
    transform= LabelRepeatedLetters(transform)
    transform_ordered= sorted(transform, key=SortCorrect)
    
    # index = transform_ordered.index('$')  
    # follow sorted, start from 0 which is $ character 
    index = transform_ordered.index('$')  
    # print(index)
    inverse = ''
    for _ in range(len(transform)):
        # add the char to the inversed string 
        inverse = inverse + transform_ordered[index][0]
        # get the current letter in the ordered transform 
        curr_index_letter = transform_ordered[index]
        # get the index of the current letter in the unsorted 
        # this corr to index of next char to add 
        next_index = transform.index(curr_index_letter) 
        index = next_index
    # move $ to end 
    fixed_dollar_sign= inverse[1:]+inverse[0]
    # fixed_dollar_sign= inverse[::-1]
    return fixed_dollar_sign

def LabelRepeatedLetters(string):
    # make dict of copunts for each letter 
    dict_counts={}
    for letter in string: 
        if letter not in dict_counts: 
            dict_counts[letter]=0
        else: 
            dict_counts[letter]+=1
    
    #print(dict_counts)
    # go over letters, increment by 1 if not 0 for correct indexing 
    for letter in dict_counts:
        if dict_counts[letter]!=0: 
            dict_counts[letter]+=1
    
    # reconstruct string 
    reconstructed_str=[]
    for letter in string[::-1]: 
        if dict_counts[letter]==0: 
            reconstructed_str.append(letter)
        else: 
            reconstructed_str.append(letter+ str(dict_counts[letter])) 
            dict_counts[letter]-=1
        # print(dict_counts)
    #print("reconst", reconstructed_str[::-1])
    return reconstructed_str[::-1]

# custom sort func to sort correctly
def SortCorrect(transform):
    # split to alpha and numeric to sort numeric correctly 
    letter = ''.join(filter(str.isalpha, transform))
    number = ''.join(filter(str.isdigit, transform))
    # get int 
    num_value = int(number) if number else 0
    # return a tuple with the alphabet part and the numeric value
    return (letter, num_value)

###########################################
# 11.4- BWMatching: 
def bw_matching(bwt: str, patterns: List[str]) -> List[int]:
    """
    Perform Burrows-Wheeler Matching for a set of patterns against the Burrows-Wheeler Transform of a text.
    """
    #pass
    last_column, LastToFirst = GetLastToFirst(bwt)

    # perform BWMatching for each pattern, add results into list 
    indices= []
    for pattern in patterns: 
        index= BWMatching(last_column, pattern, LastToFirst)
        indices.append(index)
    # return matched indices 
    return indices 

def BWMatching(LastColumn, Pattern, LastToFirst): 
    top= 0 
    bottom = len(LastColumn)

    # binary serach to get matching psoitions 
    while top <= bottom: 
        if len(Pattern)>0:
        # while Pattern:
            # get last symbol 
            symbol = Pattern[-1]
            # get pattern without last symbol 
            Pattern= Pattern[:-1]
            # print(symbol, LastColumn[top:bottom+1])
            # get symbols in the curr block from the last column 
            block_last = LastColumn[top:bottom+1]
            # print("block_last", block_last)
            symbols_last_col= [char[0] for char in block_last]
            # check if last symbol of  pattern is in current block
            if symbol in symbols_last_col: 
                # print(symbol)
                # get top and bottom indices and update them 
                topIndex, bottomIndex= GetTopAndBottomIndices(LastColumn, top, bottom, symbol)
                # print(topIndex, bottomIndex)
                top= LastToFirst[topIndex]
                bottom= LastToFirst[bottomIndex]
                # print(top, bottom, topIndex, bottomIndex)
            else: 
                return 0 
        else: 
            # print(topIndex, bottomIndex)
            # return count of occurrences btw top and bottom 
            return bottom-top+1

def GetTopAndBottomIndices(LastColumn, top, bottom, symbol):
    topIndex = None
    bottomIndex = None
    for i, char in enumerate(LastColumn):
        if i>=top and i<=bottom: 
        # print(i, char[0], symbol)
            if char[0] == symbol:
                # print("yes")
                if topIndex is None:
                    topIndex = i
                bottomIndex = i
    # print(topIndex, bottomIndex)
    return topIndex, bottomIndex

def GetLastToFirst(string):
    # make dict of copunts for each letter 
    dict_counts={}
    for letter in string: 
        if letter not in dict_counts: 
            dict_counts[letter]=0
        else: 
            dict_counts[letter]+=1
    
    # print(dict_counts)
    # go over letters, increment by 1 if not 0 for correct indexing 
    for letter in dict_counts:
        if dict_counts[letter]!=0: 
            dict_counts[letter]+=1
    
    # reconstruct string 
    reconstructed_str=[]
    for letter in string[::-1]: 
        if dict_counts[letter]==0: 
            reconstructed_str.append(letter)
        else: 
            reconstructed_str.append(letter+ str(dict_counts[letter])) 
            dict_counts[letter]-=1
        # print(dict_counts)
    # print("reconst", reconstructed_str[::-1])

    last = reconstructed_str[::-1]
    first =  sorted(last, key=SortCorrect)
    # print(last, first)

    LastToFirst=[]
    for i in range(len(last)):
        curr_index_letter = last[i]
        index = first.index(curr_index_letter) 
        LastToFirst.append(index)

    # print("lasttofirst", LastToFirst)
    return last, LastToFirst

# custom sort func to sort correctly
def SortCorrect(transform):
    # split to alpha and numeric to sort numeric correctly 
    letter = ''.join(filter(str.isalpha, transform))
    number = ''.join(filter(str.isdigit, transform))
    # get int 
    num_value = int(number) if number else 0
    # return a tuple with the alphabet part and the numeric value
    return (letter, num_value)

###########################################
# 11.7- Multiple Pattern Matching: 
def multiple_pattern_matching(text: str, patterns: List[str]) -> Dict[str, List[int]]:
    # get suff array 
    suff_array= suffix_array(text)
    found_indices = {}
    # for each pattern
    for pattern in patterns:
        # match it 
        first_index, last_index = PatternMatch(pattern, suff_array, text)
        if first_index!=-1 and last_index!=-1:
            # pattern found, get the pattern from orig str based on suff array ind
            found_indices[pattern] = suff_array[first_index:last_index+1]
        else: 
            # no pattern found
            found_indices[pattern]=[]
    return found_indices

def PatternMatch(pattern, suff_array, text):
    # init indices 
    minIndex = 0
    maxIndex = len(suff_array) - 1

    # init first and last occurences to return 
    first_index=-1
    last_index=-1

    # get the first occurrence of the pattern
    while minIndex <= maxIndex:
        midIndex = (minIndex + maxIndex) // 2
        # get pos of ind in str from suff array 
        position = suff_array[midIndex]
        # get suffix 
        suffix = text[position:position+len(pattern)]
        # compare aptter and suffix 
        if pattern == suffix:
            # first index= mid
            first_index = midIndex
            maxIndex = midIndex - 1  
        elif pattern > suffix:
            # take top half 
            minIndex = midIndex + 1
        else:
            # tale bottom half 
            maxIndex = midIndex - 1

    minIndex = 0
    maxIndex = len(suff_array) - 1

    # get the last occurrence of the pattern
    while minIndex <= maxIndex:
        midIndex = (minIndex + maxIndex) // 2
        position = suff_array[midIndex]
        suffix = text[position:position+len(pattern)]
        if pattern == suffix:
            last_index = midIndex
            minIndex = midIndex + 1  
        elif pattern > suffix:
            # take top half 
            minIndex = midIndex + 1
        else:
            # tale bottom half 
            maxIndex = midIndex - 1

    # return first and last index occurences 
    return first_index, last_index

    
def partial_suffix_array(text: str, k: int) -> List[Tuple[int, int]]:
    """
    Generate a partial suffix array for the given text and interval K.
    """
    # pass
    suff_array = suffix_array(text)
    retained_indices=[]
    for i in range(len(suff_array)): 
        if suff_array[i]%k==0: 
            retained_indices.append((i, suff_array[i]))
    return retained_indices
            
        
def suffix_array(text: str) -> List[int]:
    """
    Generate the suffix array for the given text.
    """
    # pass
    text= text+'$'
    arrays_with_indices={}
    for i in range(0, len(text)): 
        suffix= text[i:]
        arrays_with_indices[suffix]= i
    # return arrays_with_indices
    suffixes= list(arrays_with_indices.keys())
    suffixes.sort()
    # print(suffixes)
    sorted_indices=[]
    for suffix in suffixes: 
        index = arrays_with_indices[suffix]
        sorted_indices.append(index)
    return sorted_indices

###########################################
# 11.8- Multiple Approximate Pattern Matching: 
def multiple_approximate_pattern_matching(text: str, patterns: List[str], d: int) -> Dict[str, Iterable[int]]:
    """
    Find all starting positions in text where each string from patterns appears as a substring, allowing for up to d mismatches.
    """
    #pass

    # add end char to text, needed for suff array and bwt 
    text+='$'

    # get bwt 
    bwt_transform = burrows_wheeler_transform(text)
    # get suffix array 
    suff_array = suffix_array(text) 

    # our letters 
    letters= ('$', 'A', 'C', 'G', 'T')

    # get dict of count arrays for each letter 
    counts= GetCounts(bwt_transform, letters)
    #print("COUNTS", counts)
    # get dict of first occurences for each letter 
    first_occurence= GetFirstOcuurence(letters, counts)
    #print("first_occurence", first_occurence)

    # dict to store start positions for each pattern 
    start_positions= {}
    # loop over all patterns 
    for pattern in patterns: 
        starts=[]
        # get k to know the k-mer length for text division 
        k = int(len(pattern)/(d+1))

        # iterate by d for d mismatches 
        for i in range(d):
            # get the pattern sub segment of length k: seed 
            seed = pattern[i * k: (i + 1) * k]
            # get seed start index in original pattern 
            start_p = i * k
            # print(sub_p, start_p)
            # get top and bottom pointers and seed index 
            top = 0
            bottom = len(bwt_transform) 
            seed_index = len(seed) - 1
            # while top<bottom
            while top <= bottom:
                # if the seed index is > 0
                if seed_index >= 0:
                    # get char 
                    char = seed[seed_index]
                    # chek if we have an occurence of the character 
                    counts_top= counts[char][top]
                    counts_bottom= counts[char][bottom]
                    if counts_top != counts_bottom: 
                        # update top and bottom pointers 
                        top = first_occurence[char] + counts_top
                        bottom = first_occurence[char] + counts_bottom
                        # decrement current seed index to go to prev 
                        seed_index -= 1
                    else: # no occurence, exit loop 
                        top = bottom + 1
                # at beg of substring, check matches
                else: 
                    # go over range 
                    for i in range(top, bottom):
                        # get index in original text
                        to_add= suff_array[i] - start_p
                        # count mismatches 
                        mismatches=0 
                        substr= text[to_add:len(pattern)+to_add]
                        for i in range(len(substr)):
                            if pattern[i] != substr[i]:
                                mismatches += 1
                        # if mismatches <=d, return substr start pos 
                        if mismatches<=d:
                            for char in pattern: 
                                flag = True
                                if char not in bwt_transform:
                                    flag = False
                                if flag== True: 
                                    if to_add not in starts: 
                                        starts.append(to_add)
                            # if to_add not in starts: 
                                # starts.append(to_add)
                    break
        # account for end one 
        seed= pattern[d*k:len(pattern)]
        start_p= d*k
        top = 0
        bottom = len(bwt_transform) 
        seed_index = len(seed) - 1
        while top <= bottom:
            # if the seed index is > 0
            if seed_index >= 0:
                # get char 
                char = seed[seed_index]
                # chek if we have an occurence of the character 
                if counts[char][bottom] != counts[char][top]: 
                    # update top and bottom pointers 
                    top = first_occurence[char] + counts[char][top]
                    bottom = first_occurence[char] + counts[char][bottom] 
                    # decrement current seed index to go to prev 
                    seed_index -= 1
                else: # no occurence, end loop 
                    top = bottom + 1
            # at beg of substring 
            else:
                # go over range 
                for i in range(top, bottom):
                    # get index in original text
                    to_add= suff_array[i] - start_p
                    # count mismatches 
                    mismatches=0 
                    substr= text[to_add:len(pattern)+to_add]
                    for i in range(len(substr)):
                        if pattern[i] != substr[i]:
                            mismatches += 1
                    # if mismatches <=d, return substr start pos 
                    if mismatches<=d:
                        if to_add not in starts: 
                            flag = True
                            for char in pattern: 
                                if char not in bwt_transform:
                                    flag = False
                            if flag== True: 
                                starts.append(to_add)
                break
        # add indies to dictionary for this pattern 
        if starts==[]:
            mismatches=0 
            substr= text[0:len(pattern)]
            for i in range(len(substr)):
                if pattern[i] != substr[i]:
                    mismatches += 1
            if mismatches<=d:
                start_positions[pattern]= [0]
            else: 
                start_positions[pattern]= starts
        else: 
            start_positions[pattern]= starts
    # return dict 
    return start_positions

def GetCounts(bwt_transform, letters): 
    # init count dict 
    counts = {}
    # init arrays for wach letter 
    for char in letters:
        arr_counts=[0]
        for c in range(len(bwt_transform)):
            arr_counts.append(0)
        counts[char]=arr_counts
    # fill in arrays with counts, loop over bwt 
    for i in range(len(bwt_transform)):
        for char in counts:
            counts[char][i + 1] = counts[char][i]
        counts[bwt_transform[i]][i + 1] += 1
    return counts 

def GetFirstOcuurence(letters, counts): 
    first_occurence={}
    index = 0
    # go over chars 
    for char in letters:
        # get char first occ in dict 
        first_occurence[char] = index
        # offset the index by number of occurences of the char from count array 
        num_char= counts[char][-1]
        index += num_char
    return first_occurence

def burrows_wheeler_transform(text: str) -> str:
    """
    Generate the Burrows-Wheeler Transform of the given text.
    """
    # pass
    # text+='$'
    text_rotations=[]
    for i in range(len(text), 0, -1): 
        word= text[i:] + text[:i]
        text_rotations.append(word)
    # print(text_rotations)
    sorted_text_rotations= sorted(text_rotations)
    # print(sorted_text_rotations)
    end_letters=[]
    for rotated_word in  sorted_text_rotations:
        end_letter= rotated_word[-1]
        end_letters.append(end_letter)
    return ''.join(end_letters)

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

