#author: Anas Elghafari
#building a toy profile HMM in two libraries, GHMM and YAHMM,
#with the purpose of comparing the viterbi calculations in both libraries


from yahmm import *
import ghmm
random.seed(0)  #needed for yahmm ??



#ghmm model:
alphabet = ghmm.Alphabet(['A', 'C', 'G', 'T'])
initial_probs = ([1] + [0]*13)
start_em = end_em = [0.25, 0.25, 0.25, 0.25]
m1_em = [0.8, 0.1, 0.05, 0.05]
m2_em = [0.8, 0.1, 0.05, 0.05]
m3_em = [0.1, 0.8, 0.05, 0.05]
m4_em = [0.1, 0.8, 0.05, 0.05]
i1_em = i2_em = i3_em = i4_em = [0.25, 0.25, 0.25, 0.25]
d1_em = d2_em = d3_em = d4_em = [0, 0, 0, 0]
#transitions:
m1_trans = [0, 0, 0.2, 0.2, 0.6]+ ([0]*9)
m2_trans = [0]*5 + [0.2, 0.2, 0.6]+ [0]*6
m3_trans = [0]*8 + [0.2, 0.2, 0.6] + [0]*3
m4_trans = [0]*11 + [0.2, 0.2, 0.6]

i1_trans = [0, 0, 0.2, 0, 0.8] + [0]*9
i2_trans = [0]*5 + [0.2, 0, 0.8] + [0]*6
i3_trans = [0]*8 + [0.2, 0, 0.8] + [0]*3
i4_trans = [0]*11 + [0.2, 0, 0.8]

d1_trans = [0]*2 + [0.1, 0, 0.9] + [0]*9
d2_trans = [0]*5 + [0.1, 0, 0.9] + [0]*6
d3_trans = [0]*8 + [0.1, 0, 0.9] + [0]*3
d4_trans = [0]*11 + [0.1, 0, 1]

start_trans = [0.5, 0.5] + [0]*12
#end_trans = [0]*13 + [1]  doesn't match what we have for yahmm
end_trans = [0]*14

transitions =  [start_trans, m1_trans, i1_trans, d1_trans,
                m2_trans, i2_trans, d2_trans,
                m3_trans, i3_trans, d3_trans,
                m4_trans, i4_trans, d4_trans,
                end_trans]

emissions = [start_em, m1_em, i1_em, d1_em,
             m2_em, i2_em, d2_em,
             m3_em, i3_em, d3_em,
             m4_em,i4_em, d4_em,
             end_em]

ghmm_model = ghmm.HMMFromMatrices(alphabet, ghmm.DiscreteDistribution(alphabet),
                             transitions, emissions, initial_probs)


print(ghmm_model)
#print("A sample:\n\n")
#print(ghmm_model.sampleSingle(30))


#yahmm:
yahmm_model = yahmm.Model(name="ProfileHMM")
s_d = yahmm.DiscreteDistribution({'A':0.25, 'C':0.25, 'G':0.25,
                                'T':0.25})
m_dist1 = yahmm.DiscreteDistribution({'A':0.8, 'C':0.1, 'G':0.05,
                                'T':0.05})
m_dist2 = yahmm.DiscreteDistribution({'A':0.1, 'C':0.8, 'G':0.05,
                                'T':0.05})
start = State(s_d, "START")
end = State(s_d, "END")
m1 = State(m_dist1, name="m1")
m2 = State(m_dist1, name="m2")
m3 = State(m_dist2, name="m3")
m4 = State(m_dist2, name="m4")
i1 = State(s_d, name="i1")
i2 = State(s_d, name="i2")
i3 = State(s_d, name="i3")
i4 = State(s_d, name="i4")
d1 = State(None, name="d1")
d2 = State(None, name="d2")
d3 = State(None, name="d3")
d4 = State(None, name="d4")
#print(yahmm_model.states)
yahmm_model.add_state(start)
yahmm_model.add_state(m1)
yahmm_model.add_state(i1)
yahmm_model.add_state(d1)
yahmm_model.add_state(m2)
yahmm_model.add_state(i2)
yahmm_model.add_state(d2)
yahmm_model.add_state(m3)
yahmm_model.add_state(i3)
yahmm_model.add_state(d3)
#print(yahmm_model.states)
yahmm_model.add_state(m4)
yahmm_model.add_state(i4)
yahmm_model.add_state(d4)
#yahmm_model.add_state(State(None, "silent test"))
yahmm_model.add_state(end)
yahmm_model.add_transition(yahmm_model.start, start, 1.0)
yahmm_model.add_transition(start, start, 0.5)
yahmm_model.add_transition(start, m1, 0.5)
yahmm_model.add_transition(m1, i1, 0.2)
yahmm_model.add_transition(m1, d1, 0.2)
yahmm_model.add_transition(m1, m2, 0.6)
yahmm_model.add_transition(i1, i1, 0.2)
yahmm_model.add_transition(i1, m2, 0.8)
yahmm_model.add_transition(d1, i1, 0.1)
yahmm_model.add_transition(d1, m2, 0.9)

yahmm_model.add_transition(m2, i2, 0.2)
yahmm_model.add_transition(m2, d2, 0.2)
yahmm_model.add_transition(m2, m3, 0.6)
yahmm_model.add_transition(i2, i2, 0.2)
yahmm_model.add_transition(i2, m3, 0.8)
yahmm_model.add_transition(d2, i2, 0.1)
yahmm_model.add_transition(d2, m3, 0.9)

yahmm_model.add_transition(m3, i3, 0.2)
yahmm_model.add_transition(m3, d3, 0.2)
yahmm_model.add_transition(m3, m4, 0.6)
yahmm_model.add_transition(i3, i3, 0.2)
yahmm_model.add_transition(i3, m4, 0.8)
yahmm_model.add_transition(d3, i3, 0.1)
yahmm_model.add_transition(d3, m4, 0.9)


yahmm_model.add_transition(m4, i4, 0.2)
yahmm_model.add_transition(m4, d4, 0.2)
yahmm_model.add_transition(m4, end, 0.6)
yahmm_model.add_transition(i4, i4, 0.2)
yahmm_model.add_transition(i4, end, 0.8)
yahmm_model.add_transition(d4, i4, 0.1)
yahmm_model.add_transition(d4, end, 0.9)

#yahmm_model.add_transition(end, end, 0.5)
#this doesn't match what we have for  ghmm

#yahmm_model.add_transition(end, yahmm_model.end, 1)

yahmm_model.bake(merge='None', verbose=True)
print(yahmm_model.states)


#takes a string and returns EMissionSequence object,
#needed because ghmm.viterbi can't be called directly on a string or a list of symbols.
def ghmm_seq(s):
    gs = ghmm.EmissionSequence(alphabet, list(s))
    return gs


def printv(result):
    r = [obj.name for (state, obj) in result[1]]
    r.append(result[0])
    return r


#smaller model  
def y_test():
    y = Model("testing-yahmm")
    dist = DiscreteDistribution({'A':0.25, 'C':0.25, 'G':0.25,
                                'T':0.25})
    dist2 = DiscreteDistribution({'A':0.8, 'C':0.1, 'G':0.05,
                                'T':0.05})
    s1 = State(dist, "S1")
    s2 = State(dist2, "S2")
    s3 = State(None, "S3")
    y.add_state(s1)
    y.add_state(s2)
    y.add_state(s3)
    y.add_transition(y.start, s1, 1)
    y.add_transition(s3, s2, 1)
    y.add_transition(s1, s2, 0.2)
    y.add_transition(s1, s3, 0.3)
    y.add_transition(s1, y.end, 0.5)
    y.add_transition(s2, y.end, 1)
    y.bake(merge= 'None', verbose= True)
    return y


def y_test2():
    y = Model("testing-yahmm 2")
    dist = DiscreteDistribution({'A':0.24, 'C':0.25, 'G':0.25,
                                'T':0.25})
    s1 = State(dist, "S1")
    y.add_state(s1)
    y.add_transition(y.start, s1, 1)
    y.add_transition(s1, s1, 1)
    y.bake(merge= 'None', verbose= True)
    return y


def y_test3():
    y = Model("testing-yahmm")
    dist = DiscreteDistribution({'A':0.25, 'C':0.25, 'G':0.25,
                                'T':0.25})
    dist2 = DiscreteDistribution({'A':0.8, 'C':0.1, 'G':0.05,
                                'T':0.05})
    dist3 = DiscreteDistribution({'G': 1.0})
    s0 = State(dist3, "S0")
    s1 = State(dist, "S1")
    s2 = State(dist2, "S2")
    s3 = State(dist, "S3")
    y.add_state(s0)
    y.add_state(s1)
    y.add_state(s2)
    y.add_state(s3)
    y.add_transition(y.start, s0, 1)
    y.add_transition(s0, s1, 1)
    y.add_transition(s1, s2, 1)
    y.add_transition(s2, s3, 1)
    y.add_transition(s3, y.end, 1)
    y.bake(merge= 'None', verbose= True)
    return y


def BLOSUM_sum(filename):
    freq_ij= []
    freq_ii = []  #diagonal 
    f = open(filename, 'r')
    for line in f:
        freq_ii.append((line.split())[-1])
        freq_ij.append((line.split())[0:-1])
    return freq_ii, freq_ij


                       


    


