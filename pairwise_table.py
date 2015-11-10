#!/usr/bin/env python2
# Author:  Rachel Ehrlich
# This program takes as input the folder with the output from roary, an output
# folder and a nickname for the results.  It creates five files in the
# output directory with the prefix the supplied nickname with the pairwise
# tables and summary statistics.
# In the pairwise table files, each number in the matrix refers to the gene
# counts for the corresponding strains.
# similarity = present in both strains (includes core)
# difference = present in exactly one of the strains (xor)
# comparison = similarity - difference
# pair unique = present in only those two stains

import numpy as np
from collections import namedtuple
import sys

# Input is the folder with the roary output files
# This takes the gene_presence_absence.csv file and converts it to a gene
# presence and absence matrix
# Returns the matrix and the strain names
def get_pres_abs_mat(folder):
    with open(folder + "/gene_presence_absence.csv", 'rU') as f:
        lines_end_quotes = f.read().split('\n')
    lines = [x[1:-1] for x in lines_end_quotes]   
    
    # The file is a csv but the annotations contain commas...
    split_lines = [x.split('","') for x in lines if len(x) > 0]
    
    col_headings = split_lines[0][11:]
    rows = split_lines[1:]
    poss_data = np.zeros((len(rows), len(col_headings)), dtype=bool)

    for i in xrange(len(rows)):
        row = rows[i]
        gene = row[0]
        copy_num = row[11:]
        poss_data[i] = np.array([len(x) > 0 for x in copy_num], dtype=bool)

    return poss_data, col_headings

# Input: two numpy column vectors
# returns the number of positions where the value is true in both strains
def get_sim(vec1, vec2):
    return sum(np.logical_and(vec1, vec2))

# Input: two numpy column vectors
# returns the number of positions where the value is true in exactly one strain
def get_diff(vec1, vec2):
    return sum(np.logical_xor(vec1, vec2))

# Input:  the gene possession matrix and the column indices to consider
# Output: the number of rows where the value is true
# in only the two specified columns
def get_pair_unique(poss_mat, i, j):
    a = np.sum(poss_mat, axis=1) == 2
    b = poss_mat[:,i]
    c = poss_mat[:,j]
    return sum(np.logical_and(np.logical_and(a, b), c))

# Inputs: the gene possession matrix and number of strains it contains
# Output: the strain_pairs data structure which is a dict() whose keys are
# tuples of strain pair index numbers and whose values are a named tuple
# with counts for the similarity, difference, comparison and pair unique
# for the two strains
def compare_all_strain_pairs(poss_mat, num_strains):
    num_genes = namedtuple('num_genes', ['sim', 'diff', 'comp', 'pair_unique'])
    strain_pairs = dict()
    
    # Iterate over all unique pairs of strains
    for i in xrange(num_strains):
        for j in xrange(i + 1, num_strains):
            s = get_sim(poss_mat[:,i], poss_mat[:,j])
            d = get_diff(poss_mat[:,i], poss_mat[:,j])
            c = s - d
            p = get_pair_unique(poss_mat, i, j)
            strain_pairs[(i,j)] = num_genes(s, d, c, p)
    return strain_pairs

def tests1():
    a = np.array([True, True], dtype=bool)
    b = np.array([False, False], dtype=bool)
    assert(get_sim(a,a) == 2)
    assert(get_sim(a,b) == 0)
    assert(get_sim(b,a) == 0)
    assert(get_sim(b,b) == 0)
    
    assert(get_diff(a,a) == 0)
    assert(get_diff(a,b) == 2)
    assert(get_diff(b,a) == 2)
    assert(get_diff(b,b) == 0)
    
    m = np.zeros((5,5), dtype=bool)
    m[0,0] = True
    m[0,1] = True
    assert(get_pair_unique(m,0,1) == 1)
    assert(get_pair_unique(m,1,0) == 1)
    assert(get_pair_unique(m,0,2) == 0)
    assert(get_pair_unique(m,2,3) == 0)
    print 'tests pass'

# Makes a single table containing all four types of output
def make_output(strain_pairs, col_headings, num_strains):
    prefix = ''
    output = list()
    output.append("Strain\t" + '\t'.join(col_headings[1:]))
    row_labels = ['\tSimilarity', '\tDifference', '\tComparison', '\tPairUnique']

    for row in xrange(num_strains - 1):
        for i in xrange(4):
            curr = [prefix + col_headings[row]]
            for col in xrange(row + 1, num_strains):
                curr.append(strain_pairs[(row, col)][i])
            output.append('\t'.join(map(str, curr)) + row_labels[i])
        prefix = prefix + '\t'
    return '\n'.join(output)

# Makes three tables for the similarity, difference and pair unique output
def make_output_3(strain_pairs, col_headings, num_strains):
    row_labels = ['similarity', 'difference', 'comparison', 'pair_unique']
    rev_col_headings = list(reversed(col_headings[1:]))
    outputs = []
    
    for i in [0, 1, 3]:
        output = []
        output.append("Strain\t" + '\t'.join(rev_col_headings))
        for row in xrange(num_strains - 1):
            curr = [col_headings[row]]

            for col in list(reversed(xrange(row + 1, num_strains))):
                curr.append(strain_pairs[(row, col)][i])
            output.append('\t'.join(map(str, curr)))
        outputs.append(('\n'.join(output), row_labels[i]))

    return outputs

# Gets all data values for the output type specified by index
# example:  returns all pairwise similarity counts
def get_all_values(strain_pairs, index):
    return [val[index] for k, val in strain_pairs.items()]

# Computes the min, max, mean and standard deviation for the four output types
def calc_stats(strain_pairs):
    row_labels = ['Similarity', 'Difference', 'Comparison', 'PairUnique']
    col_labels = ['Min', 'Max', 'Mean', 'StdDev']
    output = list()
    output.append('\t' + '\t'.join(col_labels))
    for i in xrange(4):
        curr = [row_labels[i]]
        vals = get_all_values(strain_pairs, i)
        curr.append(min(vals))
        curr.append(max(vals))
        curr.append(np.mean(vals))
        curr.append(np.std(vals))
        output.append('\t'.join(map(str, curr)))       
    return '\n'.join(output)

def write_output(text, file_name):
    with open(file_name, 'w') as the_file:
        the_file.write(text)
        
def main():
    in_folder = sys.argv[1]
    out_folder = sys.argv[2]
    nickname = out_folder + '/' + sys.argv[3]
    
    poss_mat, col_headings = get_pres_abs_mat(in_folder)
    num_strains = len(col_headings)
    tests1()

    strain_pairs = compare_all_strain_pairs(poss_mat, num_strains)

    text = make_output(strain_pairs, col_headings, num_strains)
    write_output(text, nickname + "_pairwise_table.txt")
    single_tables = make_output_3(strain_pairs, col_headings, num_strains)
    for (data, table_type) in single_tables:
        write_output(data, nickname + "_pairwise_" + table_type + "_table.txt")
    text = calc_stats(strain_pairs) 
    write_output(text, nickname + "_pairwise_table_stats.txt")

if __name__ == "__main__":
    main()
