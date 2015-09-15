# Author:  Rachel Ehrlich
# This program takes as input the folder with the output from roary and writes
# to standard out the input for cgs_supragenome.m which can be found at
# https://github.com/rehrlich/fsgm

import os
import sys
sys.path.append(os.getcwd())
import pairwise_table as pt
import numpy as np
from collections import Counter

# Input is a gene possession matrix as a numpy matrix
# Ouput is a dictionary, keys = number of strain, values = number of genes found
# in exactly that many strains
def get_gene_counts_dict(poss_mat):
    gene_counts = np.sum(poss_mat, axis=1)
    return Counter(gene_counts)

# Input is the gene_counts_dict and number of strains from roary
# Writes the number of genes ordered by the key numbers
def print_c_for_fsgm(gene_counts_dict, num_strains):
    ordered_gene_counts = []
    for i in range(1, num_strains + 1):
        ordered_gene_counts.append(gene_counts_dict[i])
    sys.stdout.write(' '.join(map(str, ordered_gene_counts)))
                   
def main():
    folder = sys.argv[1]
    poss_mat, cols = pt.get_pres_abs_mat(folder)
    gene_counts_dict = get_gene_counts_dict(poss_mat)
    print_c_for_fsgm(gene_counts_dict, len(cols))
    

if __name__ == "__main__":
    main()
