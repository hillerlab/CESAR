'''
helpers functions 
@author: Anas Elghafari

'''

import params
import math

#def read_sequence(s):
#    s_stripped = s.strip().replace("-", "")
#    params.sequence = s_stripped
#    params.seq_as_codons = [params.sequence[i:i+3]
#                                 for i in range(0,len(params.sequence)-2, 3)]



def read_eth_matrix(filename):
    '''
    parses the substitution frequencies from the ETH matrix in the given
    file.
    '''
    
    f = file(filename)
    codons = params.order_of_codons_eth
    matrix = {}
    #each codon is a line in the file:
    for row_c in codons:
        l = f.readline()
        values = l.split()
        values = list(map(float, values[1:]))
        pairs = zip(codons, values)
        for p in pairs:
            matrix.setdefault(p[0], {})[row_c] = p[1]
    params.eth_codon_subs_probs = matrix
    params.eth_codon_insertion_probs = get_eth_insert_probs(matrix)
    #After having set the insertion probs based on the sub probs, we can modify
    #the sub probs to add stopcodons.

    #First we decrease the original values enough to allow us to add stop
    #codon emission probability while still having the emission probs sum to 1
    for c1 in params.eth_codon_subs_probs.keys():
        for c2 in params.eth_codon_subs_probs.keys():
            original_value = params.eth_codon_subs_probs[c1][c2]
            params.eth_codon_subs_probs[c1][c2] = _normalize(original_value)
    #second: we add the stop codon emission prob
    _add_stopcodon_probs("eth")




def get_eth_insert_probs(codon_sub_probs):
    '''
    Calculates the codon insertion probabilities based on the ETH matrix
    If this function is used then codon insertion probabilities will not be even
    instead, the insertion probabilities of codon XYZ will correspond to
    the probability of a random codon mutating into XYZ.
    Some tests have shown such insertion probs to be better (i.e. yield
    better alignments) than uniform insertion probs
    '''
    raw_insertion_scores = {}
    for codon in params.order_of_codons_eth:
        l = [codon_sub_probs[x][codon]
             for x in params.order_of_codons_eth]
        raw_insertion_scores[codon] = sum(l)
    normalized_probs = {}
    scores_sum = sum(raw_insertion_scores.values())
    for codon in params.order_of_codons_eth:
        normalized_probs[codon]= raw_insertion_scores[codon]/scores_sum

    #remove stop codons and get their propability to know by how much to increase the non-stop-codons props
    stop_codon_sum = normalized_probs.pop('TAA') + normalized_probs.pop('TAG') + normalized_probs.pop('TGA')
    #print("sum of stop codon probs from eth matrix", stop_codon_sum)
    fix_factor = stop_codon_sum/1.0
    for codon, p in normalized_probs.items():
        normalized_probs[codon] = p *  (1+fix_factor)
    #print(normalized_probs)
    #print("sum" + str(sum(normalized_probs.values())))
    return normalized_probs            


def read_BLOSUM_matrix(filename):
    '''
    Parses the aminoacid substitution frequencies from the BLOSUM matrix.
    Uses the aminoacid sub probs to calculate codon sub probs, with the
    assumption that, if an aminoacid A can be encoded by codons C1, C2, C3,
    then all those three encodings are equally likely.


    The function also adds the probability that a codon can be substituted
    by a stop codon. (see function _add_stopcodon_prob).
    '''
    
    f = open(filename, 'r')
    bmatrix_raw = []
    for line in f:
        bmatrix_raw.append(line.split())
    f.close()

    
    for i in range(len(bmatrix_raw)):
        params.BLOSUM_matrix.append([])
        for j in range(len(bmatrix_raw[i])):
            params.BLOSUM_matrix[i].append(float(bmatrix_raw[i][j]))
        for j in range(len(params.BLOSUM_matrix[i]), 20):
            params.BLOSUM_matrix[i].append(float(bmatrix_raw[j][i]))

    #print(BLOSUM_matrix)
    _calculate_codon_insertion_probs()
    _calculate_aminoacid_subs_probs()
    _calculate_codon_subs_probs()
    _add_stopcodon_probs('blosum')
    _set_uniform_distribution()
    _set_uniform_dist_stop_codons()



def read_profile(profile_file):
    '''
     function to parse a profile, such as acceptor profile

     function assumes that the first line in file sepcifies the header, i.e.
     the order of columns (e.g. A   C   G  T).

     The following lines specify the emission prob for each of A, C,G, T.
     Line x specifies the emission prob at position x (state x)

    '''
    dists = []
    lines = open(profile_file).readlines()
    header = lines[0].strip().split("\t")
    #print("header now is: " + str(header))
    a_index = header.index("A")
    t_index = header.index("T")
    c_index = header.index("C")
    g_index = header.index("G")
    for l in lines[1:]:
        lparts = l.strip().split("\t")
        #print("lparts now is: " + str(lparts))
        dists.append(_build_emissions_table(lparts[a_index], lparts[t_index],
                                           lparts[c_index], lparts[g_index]))

    return dists
        
        

def _build_emissions_table(a_prob, t_prob, c_prob, g_prob):
    table = dict()
    _add_trios(table, "A", a_prob)
    _add_trios(table, "T", t_prob)
    _add_trios(table, "C", c_prob)
    _add_trios(table, "G", g_prob)
    return table


def _add_trios(table, symbol, prob):
    p = float(prob)/16
    trios = [symbol + x + y for x in ["A", "T", "C", "G"]
             for y in ["A", "T", "C", "G"]]
    for trio in trios:
        table[trio] = p
    return

    
    
                              
def _calculate_codon_insertion_probs():

    '''
    This function calculates the insertion probability for each codon, using
    the BLOSUM matrix (the frequencies version, NOT THE SCORES version).

    IN the BLOSUM frequencies matrix, a cell ij stores a number that represents
    the probability of aminoacid_i being along with aminoacid_j. Those numbers sum up
    to  1 for the entire matrix, not for each single line.

    For the i_th line (row) in the matrix, the *sum* of the numbers in that line
    represent the probability of aminoacid_i being aligned with any other aminoacid. i.e.
    the sum represents the probability  that aminoacid occuring at all.

    Some aminoacids aremore common than others. e.g.

    So the sum of each line is the probability of the aminoacid being seen.
    We also assume it's the probability of the aminoacid being inserted.
    To get the codon insertion prob we divide the aminoacid insertion prob
    by the number of codons that encode that aminoacid.
    '''
    
    total_prob = 0
    for i in range(len(params.BLOSUM_matrix)):
        aminoacid = params.aminoacids[i]
        aminoacid_prob = sum(params.BLOSUM_matrix[i])
        total_prob += aminoacid_prob
        #print("current aminoacid prob: " + str(aminoacid_prob))
        #print("aminoacid is :" + str(aminoacid))
        num_codons = len(params.aminoacids_to_codons[aminoacid])
        codon_prob = aminoacid_prob/num_codons
        for codon in params.aminoacids_to_codons[aminoacid]:
            params.codon_insertion_probs[codon] = codon_prob
        #print("total probability: " + str(total_prob))
    '''
    If we need to change to uniform insertion probs:
    
    p = 1.0/61
    
    for c in params.codons:
        params.codon_insertion_probs[c] = p

    '''
    for sc in params.stop_codons:
        params.codon_insertion_probs[sc] = 0
    




def _calculate_aminoacid_subs_probs():
    #calculate the aminacid substitution probabilities from the
    #(aminoacid1, aminoacid1) co-occurence frequencies in the BLOSUM matrix

    #This is done thus:
    #if aminoacid A co-occurs with B at freq f1, co-occurs with C at freq f2...f20, then:
    #the the prob of A turning into B (A/B) = frequence (A/B)  / frequence(A/anything)
    #prob (A/B) = f1/ (f1+f2+.....f20)
    for i in range(len(params.BLOSUM_matrix)):
        aminoacid = params.aminoacids[i]
        aminoacid_prob = sum(params.BLOSUM_matrix[i])
        normalized_row = dict()
        for j, freq in enumerate(params.BLOSUM_matrix[i]):
            #j is the index of the "target" aminoacid.
            normalized_row[params.aminoacids[j]] = (freq/aminoacid_prob)
        params.aminoacid_subs_probs[aminoacid] =normalized_row
    #print("the aminoacid substitution

        #probabilities:")
    #for k,v in params.aminoacid_subs_probs.items():
    #    print(k, v)
    #    print("\n")



def _calculate_codon_subs_probs():
    #calculating the codon substitution probabilities from the aminoacid
    #substitution probabilities.
    for aa in params.aminoacids:
        for c in params.aminoacids_to_codons[aa]:
            c_probs = dict()
            for aa2 in params.aminoacids:
                aa_to_aa2_subs_prob = params.aminoacid_subs_probs[aa][aa2]
                num_codons_for_aa2 = len(params.aminoacids_to_codons[aa2])
                for c2 in params.aminoacids_to_codons[aa2]:
                    c_to_c2_prob = _normalize(aa_to_aa2_subs_prob/
                                              num_codons_for_aa2)
                    c_probs[c2] = c_to_c2_prob
            params.codon_subs_probs[c] = c_probs

            
    
            
def _normalize(p):
    #normalizes a codon substitution probability to account for
    #the probability of a stop codon being emitted.
    #Effectively, this function makes every codon sub prob smaller
    #to "make room" for the codon/stopcodon sub prob, such that
    #all the sub probabilities for a codon still sum to 1.
    return p * (1-params.stop_codon_emission_prob)





def _add_stopcodon_probs(matrix = {"blosum", "eth"}):
    '''
    Adds the probability of substitution into a stop codon to the
    specified matrix. the probability used in specified in
    params.stop_codon_emission_prob
    '''
    if matrix == "blosum":
        table = params.codon_subs_probs
    else: #eth matrix:
        table = params.eth_codon_subs_probs
    for c in params.codons:
        for sc in params.stop_codons:
            #since there are 3 stop codons
            table[c][sc] = params.stop_codon_emission_prob/3
            

    

def _set_uniform_distribution():
    p = 1.0/64
    for c in params.codons:
        params.uniform_distribution[c] = p
    for c in params.stop_codons:
        params.uniform_distribution[c] = p

    #params.uniform_distribution = {c:p for c in params.codons
    #                               or c in params.stop_codons}


def _set_uniform_dist_stop_codons():
    #distribution for the codon match states that correspong to stop codons
    #at those states only stop codons can be emitted.
    for stop_codon in params.stop_codons:
        params.uniform_dist_stop_codons[stop_codon] = float(1)/3

        

def modify_sequence(seq, insertions):
    '''
    Adds the alignment markers (---) to the given sequence at
    the places specified in the argument 'insertions'.
    '''
    #print("insertions: " + str(insertions))
    seq_modified = ""
    s1_index = 0
    
    for i, inser in enumerate(insertions):
        seq_modified += seq[s1_index: inser[0]]
        seq_modified += "-"*inser[1]
        s1_index = inser[0]

    seq_modified += seq[s1_index:]
    return seq_modified


def insertionsCount(path, i):
    '''
    function needed for printing out the alignment.

    locates where the insertions take place, and their length
    '''
    if path[i].startswith("ntI") or path[i].startswith("I"):
        return 1
    else:
        return 0

def deletionsCount(path, i):
    '''
    function needed for printing out the alignment

    locates where the deletions take place and their length
    '''
    if _isC(2, path[i]) and _isJumpsite(path[i-1]):
        return 1
    elif _isC(3, path[i]) and _isJumpsite(path[i-1]):
        return 2
    elif _isJumpsite(path[i]) and _isJumpsite(path[i-1]):
        return span_length(path[i-1], path[i])
    elif _isC(3, path[i]) and _isC(1, path[i-1]):
        return 1
    elif _isJumpsite(path[i]) and _isC(1, path[i-1]):
        return 2
    else:
        return 0


def matchCount(path, i):
    if path[i].startswith("C"):
        return 1
    else:
        return 0

def intronsCount(path, i):
    if (path[i].startswith("BEGIN") or path[i].startswith("acc") or
        path[i].startswith("do") or path[i].startswith("END")):
        return 1
    else:
        return 0


    
def leading_introns_count(path):
    return _count_introns(path, "leading")

def trailing_introns_count(path):
    return _count_introns(path, "trailing")


def _count_introns(path, side = {"leading", "trailing"}):
    if side=="leading":
        return path.count("BEGIN2") 
    if side=="trailing":
        return path.count("END1")

    

def acceptor_skipped(path):
    '''
    takes a path as an argument and returns true if and only if
    the "acc" states were skipped over
    '''
    return not "acc1" in path





def donor_skipped(path):
    '''
    takes a path as an argument and returns true if and only if
    the donor states were skipped over,
    i.e. if the path has the transition JUMPSITE_END1
    '''
    return not "do1" in path

    '''
    if "END1" in path:
        index_end = path.index("END1")
    else:
        #If the path does not have END1 state, we look at the model-end
        #state (end2 state), which is the last state.
        index_end = len(path)-1

    #if the state prior to the 2st END1 state is a jump state, then the
    #donor site has been skipped.
    return path[index_end-1].startswith("jump") 
    '''
    
    

def _isC(i, state):
    return state.startswith("C") and state[-1] == str(i)

def _isJumpsite(state):
    return state.startswith("jump")


def span_length(jumpsiteA, jumpsiteB):
    start = int(jumpsiteA[8:])
    end = int(jumpsiteB[8:])
    return (end-start) *3
    
'''
we are not using those two printout functions anymore
def print_sequences_wrapped_firstdraft(s1, s2):
    if len(s1) < 101:  #no wrapping needed
        print(s1)
        print(s2)
    else:
        lasti = 0
        for i in range(0, len(s1)-100, 100):
            print(s1[i: i+100])
            print(s2[i:i+100] + "\n")
            lasti =  i
        print(s1[lasti:])
        print(s2[lasti:] + "\n")

def print_sequences_wrapped(s1, s2):
    s1_segments = _segment(s1, 100)
    s2_segments = _segment(s2, 100)
    for seg1, seg2 in zip(s1_segments, s2_segments):
        print(seg1)
        print(seg2 + "\n")
    
'''

def _segment(s1, n):
    s1_segments = [s1[i-n:i] for i in range(n,len(s1), n)]
    s1_segments.append(s1[-(len(s1)%n):])
    return s1_segments
    



def is_frameshift( path, index):
    '''
    Function returns true if there is a frameshift at the specified
    location in the path, e.g. if there is a transition from jump-state to
    CM3 state.
    '''
    
    s1 = path[index]
    s2 = path[index+1]
    if ((s1[0:2] == "nt" and
        not _nt_insertions_multiple_of_3(path, index))
        or
        (s1.startswith("jump") and
        s2[0] == "C" and s2[-1] == "2")
        or
        (s1.startswith("jump") and
        s2[0] == "C" and s2[-1] == "3")
        or
        (s1[0] == "C" and s1[-1] == "1" and
        s2[0] == "C" and s2[-1] == "3")
        or
        (s1[0] == "C" and s1[-1] == "1" and
        s2[0:4] == "jump")):
        return True
    else:
        return False




def _nt_insertions_multiple_of_3(path, start_index):
    #this function does not consider 0 to be a multiple of 3.
    #If there are no nt_Insertion states at the specified index,
    #it return False

    #counts the nt states in both directions


    #case where there are 0 nts  starting from the current position
    if not path[start_index].startswith("nt"):
        return False

    
    nti_count = 0
    for s in path[start_index:]:
        if s.startswith("nt"):
            nti_count += 1
        else:
            break

    for s in reversed(path[0: start_index]):
        if s.startswith("nt"):
            nti_count += 1
        else:
            break
        
    return (nti_count%3) == 0




def map_state_to_query_seq_position(path, state_index):
    '''
    The function takes a path and an index in that path (state)
    and returns the location in the sequence that
    corresponds to the state.
    '''
    
    p =  path[0:state_index]
    silent_states = [x for x in p if (x.startswith("jump") or
                                      x.endswith("-start"))]
    return len(p)-len(silent_states)






def emits_stop_codon(queryseq, path, i):
    '''
    Return True iff the states at position i in the path emit a stop codon.
    '''
    if ((_isC(1, path[i]) and _isC(2, path[i+1]) and _isC(3, path[i+2]))
         or
        (path[i].startswith("ntI") and path[i+1].startswith("ntI") and
         path[i+2].startswith("ntI"))):
        
        seqindex = map_state_to_query_seq_position(path, i)
        emitted_segment = queryseq[seqindex: seqindex+3]
        if emitted_segment in params.stop_codons:
            return True
        else:
            return False
    else:
        return False





def lastCodonInTheExon(path, index):
    '''
    True if the index marks the start of the last codon in the sequence
    i.e. what comes after the codon is a donor site.
    '''
    return path[i+3].startswith("do1")


def parse_instance(filename):
    '''
    parses instance from a FASTA file.

    instance means: 1 refernce sequence (1st sequence in file) and 1 OR MORE
                    query sequences.

    The upstream and downstream intron phases for the reference sequence must
    be indicated by number of lowercase nucleotides at the begining/end of the
    ref sequence.


    
    function can handle sequences that are already aligned (i.e. contain ---)

    

    '''
    #print("\nCurrently parsing file: " + filename + "\n")
    #incoming_phase = int(filename[-3])
    #outgoing_phase = int(filename[-1])
    #print("\nincoming intron phase: "+ str(incoming_phase))
    #print("\noutgoing intron phase: "+ str(outgoing_phase))
    f = open(filename)
    lines = f.readlines()
    ref_original = lines[1].strip()
    queries_original = [lines[i].strip() for i in range(3, len(lines), 2)]
    queries_names = [lines[i].strip()[1:] for i in range(2, len(lines)-1, 2)]
    incoming_phase = (3- leading_intron_in_sequence(ref_original))%3
    outgoing_phase = trailing_intron_in_sequence(ref_original)
    #print("\nUpstream intron phase: "+ str(incoming_phase) +
    #      ". Downstream intron phase: "+ str(outgoing_phase))
    ref = ref_original.upper().replace("-", "").replace("\n", "")
    queries = [q.upper().replace("-", "").replace("\n", "")
               for q in queries_original]
    #print("\nIN PARSE_INSTANCE. REF: " + str(ref_original))
    #print("\nIN PARSE_INSTANCE. QUERIES: " + str(queries_original))
    return (ref_original, queries_original, ref,
            queries, incoming_phase, outgoing_phase,
            queries_names)




def mark_intronic_region(seq_r, incoming_phase, outgoing_phase):
    ''' a fix for instances where the intronic NTs occuring after the
        splice and before the donor site are not marked correctly
        ideally, those instannces should be fixed at the source'''

    leading_introns =  leading_intron_in_sequence(seq_r)
    trailing_introns = trailing_intron_in_sequence(seq_r)
    seq = seq_r[leading_introns: len(seq_r)-trailing_introns]
    seq_m = (seq[0:(3-incoming_phase)%3].lower() +
             seq[(3-incoming_phase)%3:len(seq)-outgoing_phase] +
             seq[len(seq)-outgoing_phase:].lower())
    seq_f = (seq_r[0:leading_introns] +
             seq_m +
             seq_r[len(seq_r)-trailing_introns:])
    return seq_f






def pair_to_path(model, pair):
    '''
    Function takes an aligned pair and deduces the path in the model that
    corresponds to that alignment.

    Purpose is to be able to score how good an alignment is by getting the
    likelihood of the corresponding path and comparing it to the Viterbi
    '''
    max_del_span = model.paramsconfig.maximum_cd
    model_name = model.name
    
    ref_raw = pair[0].strip()
    query_raw = pair[1].strip()
    #print("\nPAIR_TO_PATH, ref_raw: " + ref_raw)
    #print("\nPAIR_TO_PATH, query_raw: " + query_raw)
    path = [model_name + "-start"]
    incoming_intron_symbols = leading_intron_in_sequence(ref_raw)
    outgoing_intron_symbols = trailing_intron_in_sequence(ref_raw)
    #print("INCOMING INTRON SYMBOLS: " + str(incoming_intron_symbols))
    #print("OUTGOING INTRON SYMBOLS: " + str(outgoing_intron_symbols))
    intron_length = leading_intron_in_sequence(query_raw)


    
    acc_states_count = len(model.paramsconfig.acc_dists)
    path += ['BEGIN2'] * (intron_length - acc_states_count - incoming_intron_symbols)
    for i in range(acc_states_count):
        path.append('acc' + str(i+1))
        
    if incoming_intron_symbols == 2:
        path.append("incoming_intron1")
        path.append("incoming_intron2")
    if incoming_intron_symbols == 1:
        path.append("incoming_intron1")
        
    path_tail = [model_name + "-end"]
    trailing_intron_length  = trailing_intron_in_sequence(query_raw)
    do_states_count = len(model.paramsconfig.do_dists)
    path_tail += ['END1'] * (trailing_intron_length - do_states_count -outgoing_intron_symbols)

    #adding the do states in reverse, because this "tail" section will be reversed
    #before it gets appended to the path.
    for i in range(do_states_count, 0, -1): 
        path_tail.append("do" + str(i))
      
    if outgoing_intron_symbols == 1:
        path_tail.append("outgoing_intron1")
    if outgoing_intron_symbols == 2:
        path_tail.append("outgoing_intron2")
        path_tail.append("outgoing_intron1")


    ref_introns_removed = ref_raw[incoming_intron_symbols:
                                  len(ref_raw)-outgoing_intron_symbols]
    #print("\nPAIR_TO_PATH, ref introns removed: " + ref_introns_removed)
    ref = ref_introns_removed.strip().replace("A", "N").replace("G", "N").replace("C", "N").replace("T", "N")
    query = query_raw[intron_length: len(pair[1])-trailing_intron_length]
    #print("\nPAIR_TO_PATH, query introns removed: " + query)
    query = query.strip().replace("A", "N").replace("G", "N").replace("C", "N").replace("T", "N")
    
    #print("\n length of ref: "+ str(len(ref)) + "\nlength of query: "+ str(len(query)))
    
    #print("in pair_to_path function, ref sequence: "+ ref)
    #print("\nquery sequence: " + ref)
    #print("intron length: "+ str(intron_length))
    #print("trailing intron length "+ str(trailing_intron_length))
    #print("query now is: " + query)
    codonmatch_indx = 0
    jumpsite_indx = 0
    ri =0
    qi = 0

    while  ri <len(ref) and qi < len(query):
        segments = (ref[ri:ri+3], query[qi:qi+3])
        #print("SEGMENTS NOW: "+ str(segments))
        if segments == ("NNN", "NNN"):
            path.append("jumpsite" + str(jumpsite_indx))
            path.append("C" + str(codonmatch_indx) + "1")
            path.append("C" + str(codonmatch_indx) + "2")
            path.append("C" + str(codonmatch_indx) + "3")
            jumpsite_indx += 1
            codonmatch_indx +=1
            qi += 3
            ri += 3
        elif segments == ("---", "NNN"):
            path.append("I" + str(codonmatch_indx -1) + "1")
            path.append("I" + str(codonmatch_indx -1) + "2")
            path.append("I" + str(codonmatch_indx -1) + "3")
            qi += 3
            ri += 3
            subsequent_ci = _dashes_ahead(query, qi+3)/3
            for sc_i in range(subsequent_ci):
                path.append("I" + str(codonmatch_indx -1) + "1")
                path.append("I" + str(codonmatch_indx -1) + "2")
                path.append("I" + str(codonmatch_indx -1) + "3")
                qi +=3
                ri +=3
        elif segments == ("NNN", "---"):
            dashes_ahead = _dashes_ahead(query, qi)
            if dashes_ahead%3 == 0:
                num_deleted_codons = dashes_ahead/3
            elif dashes_ahead%3 == 1:
                num_deleted_codons = (dashes_ahead-1)/3
            else:
                num_deleted_codons = (dashes_ahead-2)/3
                       
            
            deletion_composition = get_deletion_composition(num_deleted_codons,
                                                            max_del_span)
            
            path.append("jumpsite" + str(jumpsite_indx))
            for i in deletion_composition[:-1]:
                path.append("jumpsite" + str(jumpsite_indx+i))
                jumpsite_indx += i
                codonmatch_indx += i
            jumpsite_indx += deletion_composition[-1]
            codonmatch_indx += deletion_composition[-1]
            qi += num_deleted_codons*3
            ri += num_deleted_codons*3
        elif segments == ("NNN", "-NN"):
            path.append("jumpsite" + str(jumpsite_indx))
            path.append("C" + str(codonmatch_indx) + "2")
            path.append("C" + str(codonmatch_indx) + "3")
            jumpsite_indx += 1
            codonmatch_indx +=1
            qi += 3
            ri += 3
        elif segments == ("NNN", "--N"):
            path.append("jumpsite" + str(jumpsite_indx))
            path.append("C" + str(codonmatch_indx) + "3")
            jumpsite_indx += 1
            codonmatch_indx +=1
            qi += 3
            ri += 3
        elif segments == ("NNN", "N--"):
            path.append("jumpsite" + str(jumpsite_indx))
            path.append("C" + str(codonmatch_indx) + "1")
            jumpsite_indx += 1
            codonmatch_indx +=1
            qi += 3
            ri += 3
        elif segments == ("NNN", "N-N"):
            path.append("jumpsite" + str(jumpsite_indx))
            path.append("C" + str(codonmatch_indx) + "1")
            path.append("C" + str(codonmatch_indx) + "3")
            jumpsite_indx += 1
            codonmatch_indx +=1
            qi += 3
            ri += 3
        elif (segments == ("-NN", "NNN") or
              segments == ("-N", "NN") or
              segments == ("-", "N")):
            path.append("ntI" + str(codonmatch_indx -1))
            ri += 1
            qi += 1
        elif (segments == ("--N", "NNN") or
              segments == ("--", "NN")):
            path.append("ntI" + str(codonmatch_indx -1))
            path.append("ntI" + str(codonmatch_indx -1))
            ri += 2
            qi += 2   
        else:
            print("Unkown sequence pattern: " + str(segments) +
                            "\nref index now is: " + str(ri) + "\nquery index now is: " + str(qi))
            break
        


    
    path_tail.append("jumpsite" + str(jumpsite_indx))
    path_tail.reverse()
    path += path_tail
    #print("\nDECODED PATH:  " + str(path))
    return path



def leading_intron_in_sequence(s):
    '''
    function counts the number of intronic nt occuring at the start
    of a sequence
    '''
    seq = s.strip()
    l = 0
    for i in seq:
        if i in ['a', 'g', 'c', 't']:
            l +=1
        else:
            break
    return l


def trailing_intron_in_sequence(s):
    '''
    function counts the number of intronic nt occuring at the end
    of a sequence
    '''
    return leading_intron_in_sequence(s[::-1])



def _dashes_ahead(seq, indx):
    count = 0
    while indx < len(seq) and seq[indx] == "-":
        count += 1
        indx += 1
    return count




def aggregate_ratios(ratios):
    results = [0, 0, 0, 0]
    for r in ratios:
        if r >= 0.99 and r <= 1.01:
            results[0] += 1
        elif r>1.01 and r <=10:
            results[1] += 1
        elif r>10 and r <= 100:
            results[2] += 1
        elif r> 100:
            results[3] += 1
        else:
            print("UNEXPECTED VALUE for ratio of viterbi path to original: " +
                            str(r))
      
    if len(ratios)!= 0:
        arith_avg =  float(sum(ratios)/len(ratios))
        #print("\nsum of ratios: "+ str(sum(ratios)))
        #print("\nnum of ratios: " + str(len(ratios)))
        #print("divided: " + str(arith_avg))
        if arith_avg != 0:
            log_arith_avg = math.log(arith_avg)
        else:
            log_arith_avg = 0  #or maybe set to non-defined
    else:
        arith_avg = 0
        log_arith_avg = -1
    results.append(arith_avg)
    results.append(log_arith_avg)
    results.append(harmonic_mean(ratios))
    return results



def get_ratio_class(r):
    
    '''
    Categorizing the rations into 5 classes.
    '''
    if r >= 0.99 and r <= 1.01:
        return 1
    elif r>1.01 and  r <=10:
        return 2
    elif r>10 and r <= 100:
        return 3
    elif r > 100:
        return 4
    else:
        return "5 NONDEF"

    


def calculate_offset(query_seq):
    '''
    Calculates the needed white space when printing out the alignment.
    '''
    offset = 0
    for i, v in enumerate(query_seq):
        if v in ['A', 'C', 'G', 'T']:
            offset = i
            break

    return offset


def harmonic_mean(l):
    '''
    Calculates harmonic mean, which is better suited for handling
    outlier values than the arithmatic mean.
    '''
    return -2

    '''
    inverse = [1.0/float(x) for x in l if x!=0]
    if sum(inverse) == 0:
	return 0
    else:
        return float(len(inverse))/sum(inverse)
    '''



def get_deletion_composition(num_codons_deleted, max_del_span):
    '''
    Fuction that calculates the deletion compsotion of a multi-codon deletion.

    e.g. Assuming the model allows a maximum of 10 codon deletions.

        for a deletion composition of 6 codons. The result is [6]

        for a deletion composition of 13 codons, the results is [10,3]

        for a deletion composition of 24 codons, the result is [10,10,4]

    function is used in deducing a path from a given alignment.
    '''
    if num_codons_deleted <= max_del_span:
        return [num_codons_deleted]
    else:
        composition = [max_del_span] * int(float(num_codons_deleted)/
                                           max_del_span)
        if num_codons_deleted% max_del_span !=0:
            composition.append(num_codons_deleted%max_del_span)
        return composition



def get_exon_start(seq):
    '''
    Finds the location in a sequence where the exon starts.
    '''
    if seq.count('A') != 0:
        a = seq.index('A')
    else:
        a = len(seq)

    if seq.count('T') != 0:
        t = seq.index('T')
    else:
        t = len(seq)

    if seq.count('C') != 0:
        c = seq.index('C')
    else:
        c = len(seq)

    if seq.count('G') != 0:
        g = seq.index('G')
    else:
        g = len(seq)
    return min([a, t, g, c])







def get_exon_end(seq):
    '''
    Finds the location in a sequence where the exon ends.
    '''
    
    seq = seq[::-1]
    return len(seq) - get_exon_start(seq) -1





def get_underlying_codons(c, exclude_stop_codons):
    codons = [c]
    if c[0] == 'N':
        codons = [x+ c[1]+c[2] for x in ['A', 'T', 'C', 'G']]
    for codon in codons:
        if codon[1] == 'N':
            codons.extend([codon[0] + x + codon[2] for x in ['A', 'T', 'C', 'G']])
    for codon in codons:
        if codon[2] == 'N':
            codons.extend([codon[0] + codon[1] + x for x in ['A','T','C','G']])

    if exclude_stop_codons:
        return [x for x in codons if x.count('N') == 0
                and x not in params.stop_codons]
    else:
        return [x for x in codons if x.count('N') == 0]






def calculate_average_emissions_table(codons, matrix = {"eth", "blosum"}):
    '''
    this function constructs emission table (containing substitution probs)
    which is the average of the emission tables of the argument 'codons'.

    it is meant to be used to calculate  the emissions table for a codon
    containing N.

    example: in the case of codon: ANG, this function will be called thus:
    (matrix, [AAG, ATG, ACG, AGG]). Suppose thesubsitition prob of AAG/CCC is
    p1,  ATG/CCC = p2, ACG/CCC = p3, AGG/CCC = p4
    then:
    in ther esulting table, the value associated with CCC will be the
    arith. average of p1, p2, p3 p4
    '''

    tables = []
    if matrix.lower() == 'blosum':
        subs_probs = params.codon_subs_probs
    else:
        subs_probs = params.eth_codon_subs_probs

    for c in codons:
        tables.append(subs_probs[c])



    result = dict()
    for aa in tables[0].keys():
        values_for_aa = [t[aa] for t in tables]
        avg_value = float(sum(values_for_aa))/len(values_for_aa)
        result[aa] = avg_value

    return result

        
def get_emissions_table_N_codon(codon, matrix):
    codons = get_underlying_codons(codon, True)
    table = calculate_average_emissions_table(codons, matrix)
    return table

    

def table_to_function(table):
    #print("KEYS: " + str(len(table.keys())))
    #if len(table.keys()) < 61:
    #    print(table.keys())
    def func(c):
        if c.count('N') == 0:
            if table.has_key(c) and table[c] > 0:
                return math.log(table[c])
            else:
                return float("-inf")
        else:
            underlying_codons = get_underlying_codons(c, False)
            probs = [table[x] for x in underlying_codons]
            avg_prob = float(sum(probs))/len(probs)
            return math.log(avg_prob)if avg_prob > 0 else float("-inf")
    return func
