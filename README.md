# CESAR
CESAR (Coding Exon Structure Aware Realigner) is a tool to realign coding exons with a Hidden-Markov-Model (HMM) that directly incorporates the reading frame and splice site annotation of each exon.

Its main purpose is to facilitate assessing exon conservation and thus transferring exon annotations from a reference to aligned (query) genomes. Given a whole genome alignment between a reference genome with annotated genes and a query genome, one can extract the query sequence that aligns to a coding exon in the reference genome. However, whole genome alignment is not aware of the reading frame and the splice sites. Thus, the aligned sequence may indicate a broken reading frame and/or non-functional splice sites, even if the exon is conserved in the query. To improve this, CESAR will align the query sequence again (“realign”), considering the reading frame and splice site position of the exon. The resulting alignment will preserve the reading frame and splice sites if the query sequence contains an intact exon.

#Installation
Except for Python 2.x, CESAR uses Jacob Schreiber's YAHMM (Yet Another Hidden Markov Model) library. 
Please follow the [these instructions](https://github.com/jmschrei/yahmm) to install YAHMM.

#Input & Output
The exon sequence of the reference genome and the query sequences are given in single FASTA file. The first fasta entry is the reference exon. 
Every codon that is entirely encoded by this exon should be in upper case letters. Bases belonging to codons that are split between this exon and the up/downstream exon should be in lower case. CESAR uses the number of lower case letters at the beginning and end of the exon to determine the reading frame. 

In this example, the first exon base is position 3 in a codon split by the upstream and this exon (upstream intron is in phase 2). Exon bases 2-4 encode the first full codon (CCT). The last two exon bases are position 1 and 2 in a codon split by this and the downstream exon (downstream intron is in phase 2). Thus, the last full codon encoded by this exon is ATG. 

  `gCCTGGGAACTTCACCTACCACATCCCTGTCAGTAGTGGCACCCCACTGCACCTCAGCCTGACTCTGCAGATGaa`

NOTE: Each sequence must be given in a single line (no 80 char line breaks). 

The output goes to stdout and contains the aligned sequences in FASTA format. 

#Usage
Given the input fasta file, just do 

`CESAR/CESAR.py fasta-file --clade human`

This will create an HMM for the reference exon and align it against every query sequence in the file. The mandatory parameter --clade determines from which species the donor and acceptor profile is taken (there must be a CESAR/matrices/$clade/ directory). Right now, CESAR provides profiles for human and mouse (both very similar). 

#Examples
The examples/ directory has two example*.fa files. 
Just do 

`CESAR/CESAR.py examples/example1.fa --clade human > example1.aligned.fa`

`CESAR/CESAR.py examples/example2.fa --clade human > example2.aligned.fa`

#Parameters
Default values for the transition probabilities are given in CESAR/params.py. The default values for codon substitution probabilities are taken from http://www.biomedcentral.com/1471-2105/6/134 (CESAR/matrices/eth_codon_sub.txt). All these values can be changed as command line parameters. See 

 `CESAR/CESAR.py --help`

By default, CESAR assumes you align an internal exon that has both a donor splice site upstream and an acceptor splice site downstream. Command line parameters allow to specify that the exon is the first coding exon (no donor site), the last coding exon (no acceptor site) or if the exon has non-canonical (U12) splice sites. 




