'''
@author: Anas Elghafari
some functions to compare various codon substitution matrices
'''

stop_codons= {'TAA', 'TGA', 'TAG'}

def order_of_codons_eth():
    codons = ['AAA', 'AAC', 'AAG', 'AAT', 'ACA', 'ACC',
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
    return codons


def parse_eth_matrix(filename):
    f = file(filename)
    codons = order_of_codons_eth()
    matrix = {}
    for row_c in codons:
        l = f.readline()
        values = l.split()
        values = values[1:]
        pairs = zip(codons, values)
        for p in pairs:
            matrix.setdefault(p[0], {})[row_c] = p[1]
    return matrix


def add_stop_codons(matrix):
    for c in stop_codons:
        matrix[c] = {}
        for c2 in order_of_codons_eth():
            matrix[c][c2] = 0

    for c in order_of_codons_eth():
        for c2 in stop_codons:
            matrix[c][c2] = 0

    return

def generate_spreadsheet(filename, matrices):
    f = file(filename, 'w+')
    codons = order_of_codons_eth()
    for c in codons:
        for c2 in codons:
            l = "(" + c + " -> " + c2 + "), "
            l += _get_codon_sub_probs(c, c2, matrices)
            l += "\n"
            f.write(l)
    f.close()


def _get_codon_sub_probs(c1, c2, matrices):
    fst_val = matrices[0][c1][c2]
    line = fst_val + ", "
    for m in matrices[1:]:
        val = m[c1][c2]
        ratio = ""
        if float(fst_val) != 0:
            ratio = float(val)/float(fst_val)
        else:
            ration = "NA"
        line += str(val) + ", " + str(ratio) + ", "
    return line[0:-2]



            
