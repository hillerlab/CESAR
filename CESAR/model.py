'''
@author: Anas Elghafari
building the model
'''
import math
import time
from decimal import Decimal
from yahmm import *
#importing yahmm this  way is probably not a good practice. but too time consuming to change
import params
import helpers
f = helpers.table_to_function  #we will use this funcion to wrap around the
                               #emissions table.


def _make_one_cluster(self, i, codon):
    '''
    Creates the states corresponding to 1 codon:
    3 codon match states, 3 codon insert states, 1 jumpsite,
    and the NT_i state.
    '''


    if codon in params.stop_codons:
        table = params.uniform_dist_stop_codons
        dist = LambdaDistribution(f(table))
    else:
        # if the codon has 'N', we disambiguate it into the underlying codons
        # and then average the emission tables of those underlying codons
        if self.paramsconfig.codon_sub_matrix == 'BLOSUM':
            if codon.count('N') > 0:
                table = helpers.get_emissions_table_N_codon(codon, 'BLOSUM')
                dist = LambdaDistribution(f(table))
            else:  #'normal' codon, not containing N.
                dist = LambdaDistribution(f(params.codon_subs_probs[codon]))
        else: #eth matrix
            if codon.count('N') > 0:
                table = helpers.get_emissions_table_N_codon(codon, 'ETH')
                dist = LambdaDistribution(f(table))
            else:
                dist = LambdaDistribution(f(params.eth_codon_subs_probs[codon]))
    
    udist = self.udist
    idist = self.idist
    cstate_name = "C" + str(i)
    istate_name = "I" + str(i)
    
    self.js_states.append(State(None, "jumpsite" + str(i)))
    self.cm1_states.append(State(dist, cstate_name + "1"))
    self.cm2_states.append(State(udist, cstate_name + "2"))
    #should we create a new Distribution object for cm_3?
    self.cm3_states.append(State(udist, cstate_name + "3"))
    self.ci1_states.append(State(idist, istate_name+"1"))
    self.ci2_states.append(State(udist, istate_name+"2"))
    self.ci3_states.append(State(udist, istate_name+"3"))
    self.nti_states.append(State(udist, "ntI" + str(i)))


def _make_initial_and_final_states(self):
    '''
    Creating the states and the start and end of the model:
        -the begin state,
        -the states for the acceptor
        -states for the incoming intron phase (none if incoming phase is 0)
        -first jumpsite
        -first codon insertion trio
    AT the end of the mode:
        -states for the outgoing intron phase (none if outgoing phase is 0)
        -states for the donor site
        -the end state
        
    '''
    udist = self.udist
    idist = self.idist
    first_i_trio = self.first_i_trio
    js_states = self.js_states
    upstream_cs = self.upstream_cs_states
    downstream_cs = self.downstream_cs_states
    acc_dists = self.paramsconfig.acc_dists
    do_dists = self.paramsconfig.do_dists
    
    self.begin2 = State(udist, "BEGIN2")
    for i, d in enumerate(acc_dists):
        self.acc_states.append(State(LambdaDistribution(f(d)), "acc" + str(i+1)))
    for i, d in enumerate(do_dists):
        self.do_states.append(State(LambdaDistribution(f(d)), "do" + str(i+1)))
    self.end1 = State(udist, "END1")

    if self.paramsconfig.incoming_phase == 1:
        upstream_cs.append(State(udist, "incoming_intron1"))
        upstream_cs.append(State(udist, "incoming_intron2"))
    elif self.paramsconfig.incoming_phase == 2:
        upstream_cs.append(State(udist, "incoming_intron1"))
    else: #phase 0
        pass

    if self.paramsconfig.outgoing_phase == 1:
        downstream_cs.append(State(udist, "outgoing_intron1"))
    elif self.paramsconfig.outgoing_phase == 2:
        downstream_cs.append(State(udist, "outgoing_intron1"))
        downstream_cs.append(State(udist, "outgoing_intron2"))
    else: #phase 0
        pass
        

    first_i_trio.append(State(idist, "I-11"))
    first_i_trio.append(State(udist, "I-12"))
    first_i_trio.append(State(udist, "I-13"))

    #adding the final js state:
    js_states.append(State(None, "jumpsite" + str(len(js_states))))


def _add_all_states(self):
    '''
    Adding all states to the model.
    '''
    m = self.m
    
    for s in self.cm1_states:
        m.add_state(s)
        
    for s in self.cm2_states:
        m.add_state(s)

    for s in self.cm3_states:
        m.add_state(s)

    for s in self.ci1_states:
        m.add_state(s)

    for s in self.ci2_states:
        m.add_state(s)

    for s in self.ci3_states:
        m.add_state(s)

    for s in self.nti_states:
        m.add_state(s)

    for s in self.js_states:
        m.add_state(s)

    for s in self.acc_states:
        m.add_state(s)
    for s in self.do_states:
        m.add_state(s)
    for s in self.upstream_cs_states:
        m.add_state(s)
    for s in self.downstream_cs_states:
        m.add_state(s)
    for s in self.first_i_trio:
        m.add_state(s)
    m.add_state(self.first_nti)
    m.add_state(self.begin2)
    m.add_state(self.end1)


    

def _add_all_transitions(self):
    '''
    Specifying the transitions between the states, according to the
    model diagram (see PDF manual)
    '''
    params = self.paramsconfig
    m = self.m
    cm1_states = self.cm1_states
    cm2_states = self.cm2_states
    cm3_states = self.cm3_states
    ci1_states = self.ci1_states
    ci2_states = self.ci2_states
    ci3_states = self.ci3_states
    nti_states = self.nti_states
    acc_states = self.acc_states
    do_states = self.do_states
    upstream_cs = self.upstream_cs_states
    downstream_cs =  self.downstream_cs_states
    first_i_trio = self.first_i_trio
    first_nti = self.first_nti
    begin2 = self.begin2
    end1 = self.end1
    js_states = self.js_states

    
    for i in range(len(params.seq_as_codons)):
        #transitions starting from js state
        m.add_transition(js_states[i], cm1_states[i], params.js_c1)
        m.add_transition(js_states[i], cm2_states[i], params.js_c2)
        m.add_transition(js_states[i], cm3_states[i], params.js_c3)
        m.add_transition(js_states[i], js_states[i+1], params.js_js)

        #transitions starting from cm1 states:
        m.add_transition(cm1_states[i], cm2_states[i], params.c1_c2)
        m.add_transition(cm1_states[i], cm3_states[i], params.c1_c3)
        m.add_transition(cm1_states[i], js_states[i+1], params.c1_js)

        #transitions from cm2 states:
        m.add_transition(cm2_states[i], cm3_states[i], params.c2_c3)

        #transitions from cm3 states:
        m.add_transition(cm3_states[i], js_states[i+1], params.c3_js)
        m.add_transition(cm3_states[i], ci1_states[i], params.c3_i1)
        m.add_transition(cm3_states[i], nti_states[i], params.c3_nti)

        #transitions from ci1 states:
        m.add_transition(ci1_states[i], ci2_states[i], params.i1_i2)

        #transitions from ci2 states: 
        m.add_transition(ci2_states[i], ci3_states[i], params.i2_i3)

        #transitions from ci3 states: 
        m.add_transition(ci3_states[i], ci1_states[i], params.i3_i1)                 
        m.add_transition(ci3_states[i], js_states[i+1], params.i3_js)

        #transitions from nti states:
        m.add_transition(nti_states[i], nti_states[i], params.nti_nti)
        m.add_transition(nti_states[i], js_states[i+1], params.nti_js)


    #multiple codon deletions:
    for span, multi_cd_prob in enumerate(params.multiple_cd_probs):
        for i in range(len(params.seq_as_codons)-span-1):
            m.add_transition(js_states[i],
                             js_states[i+span+2],
                             multi_cd_prob)

    #transitions in the edge states:
    m.add_transition(m.start, begin2, params.b1_b2)
    m.add_transition(m.start, acc_states[0], params.b1_acc1)
    m.add_transition(begin2, acc_states[0], params.b2_acc1)
    m.add_transition(begin2, begin2, params.b2_b2)
    m.add_transition(do_states[-1], end1, params.do2_e1)
    m.add_transition(end1, end1, params.e1_e1)
    if self.has_end: 
        m.add_transition(do_states[-1], m.end, params.do2_e2)
        m.add_transition(end1, m.end, params.e1_e2)
    for i in range(len(acc_states)-1):
        m.add_transition(acc_states[i], acc_states[i+1], params.acc_acc)
    for i in range(len(do_states)-1):
        m.add_transition(do_states[i], do_states[i+1], params.do_do)

    if params.incoming_phase == 0:
        m.add_transition(acc_states[-1], js_states[0], params.splice_js)
        m.add_transition(acc_states[-1], first_i_trio[0], params.splice_i1)
        m.add_transition(acc_states[-1], first_nti, params.splice_nti)
        #transition allowing the model to skip the acceptor site.
        m.add_transition(begin2, js_states[0], params.skip_acc)
    elif params.incoming_phase == 1:
        m.add_transition(begin2, upstream_cs[0], params.skip_acc)
        m.add_transition(acc_states[-1], upstream_cs[0], params.acc2_ii)
        m.add_transition(upstream_cs[0], upstream_cs[1], params.ii1_ii2)
        m.add_transition(upstream_cs[-1], js_states[0], params.splice_js)
        m.add_transition(upstream_cs[-1], first_i_trio[0], params.splice_i1)
        m.add_transition(upstream_cs[-1], first_nti, params.splice_nti)
    else: #phase 2
        m.add_transition(begin2, upstream_cs[0], params.skip_acc)
        m.add_transition(acc_states[-1], upstream_cs[0], params.acc2_ii)
        m.add_transition(upstream_cs[-1], js_states[0], params.splice_js)
        m.add_transition(upstream_cs[-1], first_i_trio[0], params.splice_i1)
        m.add_transition(upstream_cs[-1], first_nti, params.splice_nti)

    if  params.outgoing_phase == 0:
        m.add_transition(js_states[-1], do_states[0], params.js_donor)
        #transition allowing the model to skip the donor site:
        m.add_transition(js_states[-1], end1, params.skip_do)
    elif params.outgoing_phase == 1:
        m.add_transition(js_states[-1], downstream_cs[0], params.js_donor)
        m.add_transition(downstream_cs[-1], do_states[0], params.oi_do)
        m.add_transition(downstream_cs[-1], end1, params.skip_do)
    else: #phase 2
        m.add_transition(js_states[-1], downstream_cs[0], params.js_donor)
        m.add_transition(downstream_cs[0], downstream_cs[1], params.oi1_oi2)
        m.add_transition(downstream_cs[-1], do_states[0], params.oi_do)
        m.add_transition(downstream_cs[-1], end1, params.skip_do)
        
    

    #transitions concerning the first insert states:
    m.add_transition(first_i_trio[0], first_i_trio[1], params.i1_i2)
    m.add_transition(first_i_trio[1], first_i_trio[2], params.i2_i3)
    m.add_transition(first_i_trio[2], first_i_trio[0], params.i3_i1)
    m.add_transition(first_i_trio[2], js_states[0], params.i3_js)
    m.add_transition(first_nti, first_nti, params.nti_nti)
    m.add_transition(first_nti, js_states[0], params.nti_js)


def _viterbi(self, s, printout = True):
    '''
    s: query sequence. treated as the sequence of emissions from the HMM point of view
    This function finds the Viterbi path in the model. ie. the most likely path in
    the model that could have generated the sequence of observations.

    printout: if set to "true", the function will print out information about the
    viterbi path (score, calculation time) and will print the alignment specified
    by the path.
    '''
    params =  self.paramsconfig
    
    s1 = params.sequence  #original sequence (reference sequence)
    s2 = s.replace("-", "") #query sequence
    #Cutting the query sequence into triplets (trios) as this is
    #the "alphabet" of the model (see pdf documentation).
    s2_trios = [s2[i:i+3] for i in range(len(s2)-2)]
    #completing the last 2  trios.
    s2_trios.append(s2[-2:] + "T")
    s2_trios.append(s2[-1] + "TT")
    t1 = time.time()
    result = self.m.viterbi(s2_trios)
    viterbi_time = time.time() - t1
    logprob = result[0]
    path_probability = Decimal(math.e)**Decimal(logprob)
    path = [state.name for (n, state) in result[1]]
    if not printout:
        alignment = _align(self, s1, s2, path,
                          self.paramsconfig.incoming_phase,
                          self.paramsconfig.outgoing_phase,
                          False)
        return (path, alignment, path_probability,
                self.compile_time, viterbi_time)
    print("Viterbi path:")
    print(path)
    print("\n")
    s1_index = 0
    s2_index = 0
    #printing out alignment between each codon in the reference sequence
    #and each codon in the query sequence.
    for i, state in enumerate(path):
        if state.endswith("start") or state.endswith("end"):
            continue
        elif state.startswith("ntI"):
            print("-:" + s2[s2_index] + "  " +  state)
            s2_index += 1

        elif state.startswith("jump"):
            #what comes after JUMP might indicate  deletion.
            nxt_state_name = path[i+1]
            nxt_state_index = nxt_state_name[-1]
            if nxt_state_name.startswith("C")   and nxt_state_index == "1":
                continue
                #matching one full codon. handled at the next iteration
            elif nxt_state_name.startswith("C") and nxt_state_index == "2":
                print(s1[s1_index] + "  :-  " +
                      state + "->" + nxt_state_name)
                s1_index += 1
            elif nxt_state_name.startswith("C") and nxt_state_index == "3":
                #deleting 2 nts
                print(s1[s1_index:s1_index+2] + "  :--  " +
                      state + "->" + nxt_state_name)
                s1_index += 2
            elif nxt_state_name.startswith("jump"):
                
                #deleting one or more codons
                deletion_span_start = int(state[8:])
                deletion_span_end = int(nxt_state_name[8:])
                deletion_span = deletion_span_end - deletion_span_start
                
                print(s1[s1_index:s1_index+(deletion_span *3)] + ":" +
                      "---  " * deletion_span +
                     state + "->" + nxt_state_name)
                s1_index += deletion_span *3
        elif state.startswith("C")and state[-1] == "1":
            #C1 states might match the codon, or -if c2 is skipped-
            #then c1 macthes only the first nt.
            nxt_state = path[i+1]
            nxt_state_index = nxt_state[-1]
            if nxt_state.startswith("C") and nxt_state_index == "2":
                #matching full codon
                print(s1[s1_index:s1_index+3] + ":" +
                      s2[s2_index:s2_index+3] + "  " + state)
                s1_index += 3
                s2_index += 3
            if nxt_state.startswith("C") and nxt_state_index == "3":
                #matching the 1st nt, deleting the 2nd nt
                print(s1[s1_index] + ":" + s2[s2_index] + "  " + state)
                s1_index +=1
                s2_index +=1
                print(s1[s1_index] + ":-  " + state + "->" + nxt_state)
                s1_index += 1
            if nxt_state.startswith("jump"):
                #matching the 1st nt, deleting the 2nd & 3rd nts
                print(s1[s1_index] + ":" + s2[s2_index] + "  " + state)
                s1_index +=1
                s2_index +=1
                print(s1[s1_index:s1_index+2] + ":--  " + state + "->" + nxt_state)
                s1_index += 2
            
        elif state.startswith("C") and state[-1] == "2":
            prv_state = path[i-1]
            if prv_state.startswith("C") and prv_state[-1] == "1":
                continue
            elif prv_state.startswith("jump"):
                #C2 state is matching the 2nd  nt
                print(s1[s1_index] + ":" + s2[s2_index] + "  " + state)
                s1_index += 1
                s2_index += 1
            else:
                print("ERROR. THIS SHOULDN'T OCCUR HERE: " +state)

        elif state.startswith("C") and state[-1] == "3":
            prv_state = path[i-1]
            pre_prv_state = path[i-2]
            if prv_state.startswith("C") and prv_state[-1] == "2":
                if pre_prv_state.startswith("C") and pre_prv_state[-1] == "1":
                    continue
                else: #pre_prv_state was jump state, so c2 and c3 match codons
                    print(s1[s1_index] + ":" + s2[s2_index] + "  " + state)
                    s1_index += 1
                    s2_index += 1
            elif (prv_state.startswith("jump") or
            ((prv_state.startswith("C") and prv_state[-1] == "1"))):
                #C3 state is matching the 3rd  nt
                print(s1[s1_index] + ":" + s2[s2_index] + "  " + state)
                s1_index += 1
                s2_index += 1
            else:
                print("ERROR. THIS SHOULDN'T OCCUR HERE: " +state)
                      
                
            #C2 is be skipped,if the previous state was C1
        elif state.startswith("I") and state[-1] == "1":
            print("---" + ":" +
            s2[s2_index:s2_index+3] + "  " + state)
            s2_index += 3
        elif state.startswith("I") and (state[-1] == "2" or state[-1] == "3"):
            continue
        elif state.startswith("END") or state.startswith("BEGIN"):
            print(state + ":" + s2[s2_index])
            s2_index += 1
        elif (state.startswith("acc") or state.startswith("do")):
            print(state + ":" + s2[s2_index])
            s2_index += 1
        elif (state.startswith("incoming") or state.startswith("outgoing")):
            print(state + ":" + s2[s2_index])
            s1_index += 1
            s2_index += 1
        else:
            print("STATE NAME ERROR. UNKNOWN STATE NAME: " +state)

    print("\nTime required to compile the model: " + str(self.compile_time))
    print("Time required to calculate viterbi: " + str(viterbi_time))
    print("log prob of the Viterbi path: " + str(logprob))
    print("probability of the Viterbi path: " + str(Decimal(math.e)**Decimal(logprob)) + "\n")
    print("\nAlignment (based on the Viterbi path):")
    
    #printing out the traditional alignment between the ref and the query:
    alignment = _align(self, s1, s2, path, self.paramsconfig.incoming_phase,
          self.paramsconfig.outgoing_phase, True)
    return (path, alignment, path_probability,
                self.compile_time, viterbi_time)




                      
    
    

def _align(self, seq1, seq2, full_path,
          incoming_phase, outgoing_phase,
          printout = True):
    '''
    Prints out the alignment between seq1 and seq2, as specified by "full_path"
    '''
    insertions = []
    deletions = []
    path = full_path[1:-1] #excluding model-start and model-end.
    incoming_intron_symbols = (3-incoming_phase)%3
    leading_introns_count = helpers.leading_introns_count(path)
    trailing_introns_count = helpers.trailing_introns_count(path)
    acc_skipped = helpers.acceptor_skipped(path)
    do_skipped = helpers.donor_skipped(path)

    #if acceptor has not been skipped over, then we take acc_states_count
    #from the acc_profile. if acceptor has been skipped over, then the count is 0
    #similarily for donor states.
    acc_states_count = len(self.paramsconfig.acc_dists) if not acc_skipped else 0
    do_states_count = len(self.paramsconfig.do_dists) if not do_skipped else 0

    #print("acc/do states_count", acc_states_count, do_states_count)
    #setting the indices to where the exonic region starts in both sequences: 
    s1_index = incoming_intron_symbols
    s2_index = leading_introns_count + acc_states_count + incoming_intron_symbols
    
    for i, s in enumerate(path):
        icount = helpers.insertionsCount(path, i)
        dcount = helpers.deletionsCount(path, i)
        #print("state: " + s + ". icount: " + str(icount) + ". dcount: "+ str(dcount))
        if icount > 0:
            insertions.append((s1_index, icount))
            s2_index += icount
        if dcount >0:
            deletions.append((s2_index, dcount))
            s1_index += dcount
    
        mcount = helpers.matchCount(path, i)
        s1_index += mcount
        s2_index += mcount

    #print("insertion locations: " + str(insertions))
    #print("deletions locations: "  + str(deletions))
    seq1_modified = helpers.modify_sequence(seq1, insertions)
    seq2_modified = helpers.modify_sequence(seq2, deletions)
    leading_space = " " * (leading_introns_count + acc_states_count)
    exon_start = leading_introns_count + acc_states_count + incoming_intron_symbols
    exon_end = len(seq2_modified) - (trailing_introns_count
                                     + do_states_count
                                     + outgoing_phase)
    
    #print("exon start: " + str(exon_start))
    #print("exon end: " + str(exon_end))

    #first line of the alignment, i.e. the ref sequence. including the leading
    #space and marking the carried over phase bps as lower case.
    firstline = (
        leading_space +
        seq1_modified[:incoming_intron_symbols].lower() +
        seq1_modified[incoming_intron_symbols: len(seq1_modified)-outgoing_phase] +
        seq1_modified[len(seq1_modified)-outgoing_phase:].lower())

    #second line of the alignment, i.e. query sequence. 
    secondline = (
        seq2_modified[:exon_start].lower() +
        seq2_modified[exon_start:exon_end] +
        seq2_modified[exon_end:].lower())              
    if printout:
        print(firstline)
        print("\n")
        print(secondline)
    return (firstline, secondline)

    

def score_alignment_against_model(self, pair):

    '''
    Takes an alignment (as a pair of lines), translates that alignment into
    a path in the model (see "helpers.pair_to_path" function), and scores that
    path against the model.

    That  score is how likely the given alignment is, according to the model.
    '''
    
    path = helpers.pair_to_path(self, pair)
    #print("path of original alignment: " + str(path))
    path_objects = [self.names_to_states[x] for x in path]
    normalized_seq = pair[1].strip().upper().replace("-", "")
    normalized_seq_trios = [normalized_seq[i:i+3]
                            for i in range(len(normalized_seq)-2)]
    normalized_seq_trios.append(normalized_seq[-2:] + "A")
    normalized_seq_trios.append(normalized_seq[-1] + "AA")
    log_score = self.m.log_probability(normalized_seq_trios,
                                       path = path_objects)
    score = Decimal(math.e)** Decimal(log_score)
    return score




def _map_names_to_states(self):
    '''
    Creates a mapping between the states names and the state objects.
    Allows easy access to the states object for the purposes of constructing a
    path and scoring it
    '''     
    states = self.m.states
    for s in states:
        self.names_to_states[s.name]= s




class HMM:
    viterbi = _viterbi
    def __init__(self, name, paramsconfig, has_end = True):
        self.name = name
        self.m = Model(name)
        self.paramsconfig = paramsconfig
        self.has_end = True
        self.names_to_states = {}
        self.udist = LambdaDistribution(f(params.uniform_distribution))
        if self.paramsconfig.codon_sub_matrix == 'BLOSUM':
            self.idist = LambdaDistribution(f(params.codon_insertion_probs))
        else: #eth matrix
            self.idist = LambdaDistribution(f(params.eth_codon_insertion_probs))

        #states:
        self.cm1_states = []
        self.cm2_states = []
        self.cm3_states = []
        self.ci1_states =[]
        self.ci2_states = []
        self.ci3_states = []
        self.nti_states = []
        self.js_states = []
        self.first_i_trio = []  #first insert trio, let's index them as -1
        self.acc_states = []
        self.do_states = []
        self.upstream_cs_states = [] #upstream codon split
        self.downstream_cs_states = []
        self.begin2  = 0         #to be set later. 
        self.end1 = 0         #begin1 is model.start. end2 is model.end
        self.first_nti = State(self.udist, "ntI-1")
        #first nt insert state. let's number it -1

        
        
        for i, c in enumerate(self.paramsconfig.seq_as_codons):
            if (c in params.stop_codons and
                i != len(self.paramsconfig.seq_as_codons)-1):
                raise ValueError("Stop codon encountered before the end of " +
                                 "the reference sequence. Codon: " + str(i) +
                                  " is: " + c)
            else:
                _make_one_cluster(self, i, c)
        _make_initial_and_final_states(self)
        _add_all_states(self)
        _add_all_transitions(self)
        compilation_start_time = time.time()
        self.m.bake(verbose=False, merge = "None")
        self.compile_time = time.time() - compilation_start_time
        _map_names_to_states(self)



    
