#!/usr/bin/env python

#Author: Anas Elghafari, for MPI-CBG, Dresden, Germany.
#Contact:  anas.elghafari@gmail.com

"""
This module is designed to provide a command-line interface to the HMM aligner
It allows the user to align a pair of sequences from a FASTA file, and to specify the main params of the HMM.
It supports multiple query sequences in the FASTA file, where the first sequence is the reference exon and all subsequent sequences are genomic sequences of different query species.
NOTE: Each sequence must be given as a single line. 

For usage options, run this module with the argument '--help'
"""

import sys
import analytics
import helpers
import params
import model
import argparse



if __name__ == "__main__":

    #Defining the command line arguments:
    parser = argparse.ArgumentParser(description=
                                     'Align a pair of sequences using HMM')
    parser.add_argument('--fs_prob', default=params.FS_PROB, type=float,
                       help='Probability of a frameshift.')
    parser.add_argument('--ci_prob', default=params.CI_PROB,  type=float, 
                       help='Probability of a transition into codon insertion state')
    parser.add_argument('--ci2_prob', default=params.CI2_PROB,  type=float, 
                       help=('Probability of a subsequent codon insertion.'))
    parser.add_argument('--total_cd_prob', default=params.TOTAL_CD_PROB, type=float, 
                       help=('Probability of deleting one or more codons. '+
                             'I.e. the sum of the transition probability of deleting a single codon and all transition ' +
                             'probabilities of deleting multiple codons'))
    parser.add_argument('--matrix', default='ETH', choices=['BLOSUM', 'ETH'],
                        help=('Substitution matrix to use.' +
                        'Choices are BLOSUM matrix and ETH matrix'))
    
    parser.add_argument('--clade', help=('The clade of sequences to be aligned.' +
                        'To see the choices of clades available see the subfolders in the \'matrices\' folder.'))
    
    parser.add_argument('--is_first_exon', type=bool, default= False,
                        help = ("set this to true to use FIRST_CODON_PROFILE instead "
                                 + "of the usual acceptor profile. Default is False."))
    parser.add_argument('--is_last_exon', type=bool, default= False,
                        help = ("set this to true to use LAST_CODON_PROFILE instead "
                                 + "of the usual donor profile. Default is False."))
    parser.add_argument('--has_AC_acceptor', type=bool, default= False,
                        help = ("set this to true to use U12 acceptor profile instead "
                                 + "of the usual acceptor profile. Default is False."))
    parser.add_argument('--has_AT_donor', type=bool, default= False,
                        help = ("set this to true to use U12 donor profile instead "
                                 + "of the usual donor profile. Default is False."))
            
                              

    parser.add_argument('--verbosity', type=int, default=0,  choices=[0, 1, 2],
                        help=('Level of details to print out along with the alignment.' +
                              '\n0: only the alignment is printed out' +
                              '\n1: level 0 + the parameters used + time taken to compile the model' +
                              '\n2: level 1 + the viterbi path' +
                              '\nDefault value is: 1'))
    parser.add_argument('input_file',
                        help='FASTA file containing the reference sequence and the query sequence(s)')  

    args = parser.parse_args()

    if args.verbosity > 0:
        print("\nTHE HMM WILL BE CONSTRUCTED WITH THE FOLLOWING VALUES:")
        print("fs_prob = " + str(args.fs_prob))
        print("ci_prob = " + str(args.ci_prob))
        print("ci2_prob = " + str(args.ci2_prob))
        print("total_cd_prob = " + str(args.total_cd_prob))
        print("clade = " + args.clade)
        path_prefix = params.MATRICES_PATH_PREFIX + args.clade + "/"
        if args.matrix.upper() == 'BLOSUM':
            print("substitution matrix = " + str(args.matrix) +
                  ". FILE: " + path_prefix + params.BLOSUM_FILE)
        else:
            print("substitution matrix = " + str(args.matrix) +
                  ". FILE: " + path_prefix + params.ETH_FILE)
        if args.is_first_exon:
            print(path_prefix + params.FIRST_CODON_PROFILE + " used instead of the usual acceptor profile")
        if args.is_last_exon:
            print(path_prefix + params.LAST_CODON_PROFILE +  " used instead of donor profile")
        if args.has_AC_acceptor:
            print(path_prefix + params.U12_ACC_PROFILE + " used instead of the usual acceptor profile")
        if args.has_AT_donor:
            print(path_prefix + params.U12_DONOR_PROFILE + " used instead of the usual donor profile")
            

    #parsing the FASTA input file:    
    instance = helpers.parse_instance(args.input_file)
    params.CLADE = args.clade
    paramsconfig = params.Params_config()
    paramsconfig.set_params((args.fs_prob, args.ci_prob, args.ci2_prob, args.total_cd_prob, args.matrix.upper()))
    paramsconfig.set_sequence(instance[2], instance[4], instance[5],
                              args.is_first_exon, args.is_last_exon,
                              args.has_AC_acceptor, args.has_AT_donor)

    if args.verbosity> 0:
        print("\nReference sequene (raw): " + instance[0])
        print("\nCodons in the reference sequence: " + str(paramsconfig.seq_as_codons))
        print("\nUpstream intron phase: " + str(instance[4])) 
        print("\nDownstream intron phase: " + str(instance[5]))
        print("\nBuilding a profile HMM for the reference sequence...")
        
    m = model.HMM(args.input_file, paramsconfig)

    #iterating over the query sequences in the file
    for i in range(0, len(instance[3])):
        results = m.viterbi(instance[3][i], args.verbosity > 0)
        if args.verbosity> 0:
            print(results)
            if i== 0: print("\nHMM was compiled in: " + str(results[3]) + " seconds.")
            print("\nQuery sequence " + str(instance[6][i]) + ": ")
            if args.verbosity > 1:
                print("\nViterbi path: " + str(results[0]))
                print("\nProbability of viterbi path: " + str(results[2]))
                print("\nAlignment based on the Viterbi path:")
        print(">referenceExon\n"+ results[1][0])
        print(">" + str(instance[6][i]) + "\n" + results[1][1])
