'''
@author: Anas Elghafari
module represnting the parameter configurations of the profile HMM.
Also has Parameters Configuration class to represent objects
each  object of this class represents a particular parameter configuration
'''

import helpers


#For each of the 8 files below, the variable 'MATRICES_PATH_PREFIX'
#specifies the path prefix up until the matix folder
#(e.g. C:\MY-PROGRAMS\PRFOFILE-HMM\MATRICES\). The varibale 'CLADE' specifies
#the clade used (e.g. human). The variables for the file names (BLOSUM_FILE,
#U12_ACC_PROFILE..etc)  specify the name of the file inside the clade folder.
#So the code will look for the file in: 'MATRICES_PATH_PREFIX/CLADE/file_name'
MATRICES_PATH_PREFIX = "CESAR/matrices/"
CLADE = "human"
BLOSUM_FILE = "BLOSUM_FREQ62"
ETH_FILE = "eth_codon_sub.txt"
ACC_PROFILE = "acc_profile.txt"
FIRST_CODON_PROFILE = "firstCodon_profile.txt"
DO_PROFILE = "do_profile.txt"
LAST_CODON_PROFILE = "lastCodon_profile.txt"
U12_ACC_PROFILE = "u12_acc_profile.txt"
U12_DONOR_PROFILE = "u12_donor_profile.txt"


BLOSUM_matrix = []  #will be read from a file. see helpers.py
aminoacid_subs_probs = dict()
codon_subs_probs = dict()
codon_insertion_probs = dict()
eth_codon_subs_probs = dict()
eth_codon_insertion_probs = dict()
uniform_distribution = dict()
uniform_dist_stop_codons = dict()
oi1_dist = dict()   #oi: outgoing intron
oi2_dist = dict()





#the aminoacids as they appear in the BLOSUM matrix: 
aminoacids =      ['Ala', 'Arg', 'Asn', 'Asp', 'Cys', 'Gln',
                   'Glu', 'Gly', 'His', 'Ile', 'Leu', 'Lys',
                   'Met', 'Phe', 'Pro', 'Ser', 'Thr', 'Trp',
                   'Tyr', 'Val']


codons = ['GCT', 'GCC', 'GCA', 'GCG', 'CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG',
          'AAT', 'AAC', 'GAT', 'GAC', 'TGT', 'TGC','CAA', 'CAG', 'GAA', 'GAG',
          'GGT', 'GGC', 'GGA', 'GGG', 'CAT', 'CAC', 'ATT', 'ATC', 'ATA',
          'TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG', 'AAA', 'AAG', 'ATG',
          'TTT', 'TTC', 'CCT', 'CCC', 'CCA', 'CCG',
          'TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC', 'ACT', 'ACC', 'ACA', 'ACG',
          'TGG', 'TAT', 'TAC', 'GTT', 'GTC', 'GTA', 'GTG']


stop_codons = {'TAA', 'TGA', 'TAG'}
stop_codon_emission_prob  = 0.00019625991262287552 #based on 1704/8,682,364


aminoacids_to_codons ={'Ala': ['GCT', 'GCC', 'GCA', 'GCG'],
                       'Arg': ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
                       'Asn': ['AAT', 'AAC'],
                       'Asp': ['GAT', 'GAC'],
                       'Cys': ['TGT', 'TGC'],
                       'Gln': ['CAA', 'CAG'],
                       'Glu': ['GAA', 'GAG'],
                       'Gly': ['GGT', 'GGC', 'GGA', 'GGG'],
                       'His': ['CAT', 'CAC'],
                       'Ile': ['ATT', 'ATC', 'ATA'],
                       'Leu': ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'],
                       'Lys': ['AAA', 'AAG'],
                       'Met': ['ATG'],
                       'Phe': ['TTT', 'TTC'],
                       'Pro': ['CCT', 'CCC', 'CCA', 'CCG'],
                       'Ser': ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'],
                       'Thr': ['ACT', 'ACC', 'ACA', 'ACG'],
                       'Trp': ['TGG'],
                       'Tyr': ['TAT', 'TAC'],
                       'Val': ['GTT', 'GTC', 'GTA', 'GTG']}


##ETH empirical codon substitution matrix. Order of codons as it appears in file
order_of_codons_eth = ['AAA', 'AAC', 'AAG', 'AAT', 'ACA', 'ACC',
                      'ACG', 'ACT', 'AGA', 'AGC', 'AGG', 'AGT', 
                      'ATA', 'ATC', 'ATG', 'ATT', 'CAA', 'CAC',
                      'CAG', 'CAT', 'CCA', 'CCC', 'CCG', 'CCT',
                      'CGA', 'CGC', 'CGG', 'CGT', 'CTA', 'CTC',
                      'CTG', 'CTT', 'GAA', 'GAC', 'GAG', 'GAT',
                      'GCA', 'GCC', 'GCG', 'GCT', 'GGA', 'GGC',
                      'GGG', 'GGT', 'GTA', 'GTC', 'GTG', 'GTT',
                      'TAA', 'TAC', 'TAG', 'TAT', 'TCA', 'TCC',
                      'TCG', 'TCT', 'TGA', 'TGC', 'TGG', 'TGT',
                      'TTA', 'TTC', 'TTG', 'TTT']




#the default transition probs:
#What the transition probabilities stand for:
#fs_prob: frameshifting probability
#ci_prob: codon insertion prob
#cd_prob: codon deletion  prob
#nti_nti_prob: prob of nti to nti transition
#ci2_prob: prob of inserting another codon.       
FS_PROB = 0.0001
CI_PROB = 0.01
CI2_PROB = 0.2
NTI_NTI_PROB = 0.25
TOTAL_CD_PROB = 0.025

#Prob of deleting 2 codons at once = 0.432 * 1CD_PROB
#prob of deleting 3 codons at once = 0.276 * 1CD_PROB...etc.
MULTIPLE_CD_FACTORS = [0.432, 0.276, 0.208, 0.164, 0.147,
                       0.133, 0.127, 0.123, 0.118]
#Prob of deleting 1 codon
CD_PROB = TOTAL_CD_PROB/(1+ sum(MULTIPLE_CD_FACTORS))
NO_LEADING_INTRONS_PROB = 0.01 #prob that there are no leading introns
NO_TRAILING_INTRONS_PROB = 0.01  #prob that there are no trailing introns
B1_ACC1 = NO_LEADING_INTRONS_PROB
B1_B2   = 1.0 - NO_LEADING_INTRONS_PROB 
B2_B2  = 0.90  #test more thoroghly. only matters with long introns
DO2_E1   = 1.0 - NO_TRAILING_INTRONS_PROB
DO2_E2   = NO_TRAILING_INTRONS_PROB
E1_E1  = 0.90 #should be more thoroghly tested. only matters with long introns
#transitions that are currently 1.0:
ACC_ACC = 1.0
ACC2_II = 1.0
DO_DO = 1.0
II1_II2 = 1.0
I1_I2 = 1.0
I2_I3 = 1.0
C2_C3 = 1.0
OI1_OI2 = 1.0





class Params_config:

    def  __init__(self):
        return
  
        

        

    def set_sequence(self, seq, incoming_phase, outgoing_phase,
                     first_exon = False, last_exon = False,
                     has_AC_acceptor = False, has_AT_donor = False):
        '''
        Takes as arguments:
            the reference sequence,
            the phase of the upstream intron,
            the phase of the dowstream intron
            optional arg: first_exon: must be set to true if the ref sequence
                is the 1st exon, then there will be no accpetor site, and
                the model will use the FIRST_CODON_PROFILE instead of ACC_PROFILE
            optional arg: last_exon: must be set to true if the ref sequence is
                the last exon. In that case the LAST_CODON_PROFILE will be
                used instead of the do_profile
            

        '''
        self.sequence = seq
        if incoming_phase not in {0,1,2} or outgoing_phase not in {0,1,2}:
            raise ValueError("Incoming and outgoing phases must be 0 or 1 or 2. "
                             +"Value as parsed is: " + str(incoming_phase))
        self.incoming_phase = incoming_phase
        self.outgoing_phase = outgoing_phase
        incoming_offset = (3-incoming_phase)%3
        self.seq_sans_introns = seq[incoming_offset:len(seq)-outgoing_phase]
        #print("\nreference sequence: " + self.sequence)
        #print("sequence sans introns: " + self.seq_sans_introns)
        self.seq_as_codons = [self.seq_sans_introns[i:i+3]
                              for i in range(0, len(self.seq_sans_introns)-2, 3)]
        
        #if self.seq_as_codons[-1] in {'TAA', 'TGA', 'TAG'}:
        #    self.seq_as_codons = self.seq_as_codons[0:len(self.seq_as_codons)-1]

        #print("\nreference sequence as codons: " + str(self.seq_as_codons))

           
        prf = MATRICES_PATH_PREFIX + CLADE + "/"
        helpers.read_BLOSUM_matrix(prf+BLOSUM_FILE)
        helpers.read_eth_matrix(prf+ETH_FILE)
        if first_exon:
            self.acc_dists = helpers.read_profile(prf+FIRST_CODON_PROFILE)
        elif has_AC_acceptor:
            self.acc_dists = helpers.read_profile(prf+U12_ACC_PROFILE)
        else: #usual case, usual acceptor profile:
            self.acc_dists = helpers.read_profile(prf+ACC_PROFILE)


        if last_exon:
            self.do_dists = helpers.read_profile(prf+LAST_CODON_PROFILE)
        elif has_AT_donor:
            self.do_dists = helpers.read_profile(prf+U12_DONOR_PROFILE)
        else: #usual case, usual donor profile
            self.do_dists = helpers.read_profile(prf+DO_PROFILE)

        

        
    def set_params(self, paramstuple):
        #THE CRUCIAL TRANSITIONS:
        #What the transition probabilities stand for:
        #fs_prob: frameshifting probability
        #ci_prob: codon insertion prob
        #cd_prob: codon deletion  prob
        #total_cd_prob: cd_prob + transitions that delete multiple codons
        #nti_nti_prob: prob of nti to nti transition
        #ci2_prob: prob of inserting another codon, right after inserting one.      
        self.fs_prob = paramstuple[0]
        self.ci_prob = paramstuple[1]
        self.ci2_prob = paramstuple[2]
        self.total_cd_prob = paramstuple[3]
        self.codon_sub_matrix = paramstuple[4]  
        self.maximum_cd = len(MULTIPLE_CD_FACTORS)
        self.cd_prob = self.total_cd_prob/(1+sum(MULTIPLE_CD_FACTORS))
        
        self.multiple_cd_probs = [self.cd_prob*factor for factor in MULTIPLE_CD_FACTORS]

        #the derived transition probs, calculated from the crucial ones.
        #Before changing any of this, make sure to understand the HMM schema
        #And make compensating changes.
        self.splice_nti = self.fs_prob
        self.splice_i1  = self.ci_prob
        self.splice_js  = 1- (self.fs_prob + self.ci_prob)
        self.js_c1   = 1- (2*self.fs_prob)- self.total_cd_prob    
        self.js_c2   = self.fs_prob
        self.js_c3   = self.fs_prob
        self.js_js   = self.cd_prob
        self.c1_c2   = 1- (2*self.fs_prob)
        self.c1_c3   = self.fs_prob
        self.c1_js   = self.fs_prob
        self.c3_js   = 1 - (self.fs_prob + self.ci_prob)
        self.c3_nti  = self.fs_prob
        self.c3_i1   = self.ci_prob
        self.i3_i1   = self.ci2_prob
        self.i3_js   = 1- self.ci2_prob
        self.nti_nti = NTI_NTI_PROB
        self.nti_js  = 1- self.nti_nti

        #Transitions at the edges of the HMM
        self.no_leading_introns_prob = NO_LEADING_INTRONS_PROB
        self.no_trailing_introns_prob = NO_TRAILING_INTRONS_PROB
        self.b1_acc1  = B1_ACC1
        self.b1_b2   = B1_B2
        self.b2_b2   = B2_B2
        self.skip_acc = self.fs_prob
        self.b2_acc1  = 1 - (self.b2_b2 +  self.skip_acc)

        
        self.skip_do = self.fs_prob
        self.js_donor  = 1.0 - self.skip_do
        self.do2_e1   =  DO2_E1
        self.do2_e2   = DO2_E2
        self.e1_e1   = E1_E1
        self.e1_e2   = 1.0 - self.e1_e1
        self.acc_acc = ACC_ACC
        self.acc2_ii = ACC2_II
        self.do_do = DO_DO
        self.oi_do = 1- self.skip_do
        self.ii1_ii2 = II1_II2
        self.i1_i2  = I1_I2
        self.i2_i3  = I2_I3
        self.c2_c3  = C2_C3
        self.oi1_oi2 = OI1_OI2

        #print("params now are: ", vars(self))
        



