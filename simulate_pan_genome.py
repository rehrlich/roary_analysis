# Author:  Rachel Ehrlich
# This program takes as input a folder with roary output, an output directory
# for the results, a nickname for the analysis, a comma separated string of 
# cutoff values where 0 < value <= 1 and the number of simulations to run.
# An example cutoff argument is "0.15,0.95,0.99,1.0"
# This does the specified number of simulations to model how the gene 
# frequencies change when strains are sequenced in a different order.
# For each of the cutoffs, a tab separated .Rtab (for consistency 
# with roary) has one line for each simulation.  The nth number in a line
# is the number of genes whose frequency in the first n strains sampled
# is next smallest cutoff <= freqnecy < cutoff.  For the case of the 
# lowest cutoff, the lower number is 0.  For the case of 1.0, the upper 
# range is <= to include 100%.

from __future__ import division    
import numpy as np
import pairwise_table as pt
import sys


# Input is a boolean gene possession matrix
# Output is an array of gene frequencies
def get_gene_freq(poss_mat):
    gene_counts = np.sum(poss_mat, 1)
    num_strains = len(poss_mat[0,:])
    freq = gene_counts / num_strains
    return freq

# Input is a gene frequencies array amd a sorted list of cutoff frequencies
# Outputs a list of counts of genes whose frequencies are less than each cutoff
# but not the previous.
def get_counts_per_bin(freq, cutoffs):
    counts = []
    for val in cutoffs: 
        counts.append((freq < val).sum())
        
    bins = []
    for index, val in enumerate(counts):
        if index > 0:
            val -= counts[index - 1]
        bins.append(val)
    return bins

# Inputs are a gene possession matrix, an array of column headings,
# a list of cutoff frequencies and the number of simulations
# Returns a list of matrices where each matrix is the gene counts
# for each frequency bin for all simulations
def simulate_reordering(poss_mat, col_headings, cutoffs, num_iter):
    tot_strains = len(col_headings)
    
    # pre allocate array
    results = [np.empty((num_iter, tot_strains), dtype=int) for x in cutoffs]

    # Each iteration is one simulation
    for x in xrange(num_iter):
        
        # pre allocate array
        temp_results = [np.empty((tot_strains), dtype=int) for z in cutoffs]
        
        # reorder the strains
        order = np.random.permutation(range(tot_strains))
        poss_mat = poss_mat[:, order]
        
        # simulate adding one strain at a time
        for num_strains in xrange(1, tot_strains + 1):

            freq = get_gene_freq(poss_mat[:, :num_strains])
            nonzero_freq = freq[np.nonzero(freq)]
            counts_per_bin = get_counts_per_bin(nonzero_freq, cutoffs)
            
            assert(len(nonzero_freq) == sum(counts_per_bin))

            for i in xrange(len(cutoffs)):
                temp_results[i][num_strains - 1] = int(counts_per_bin[i])

        for i in xrange(len(cutoffs)):
            results[i][x,:] = temp_results[i]
            
    return results

# Input is a list of lists and an output file name
# connects each list with tabs and the list of list with newlines
# writes the string to the file
def make_rtab(data, file_name):
    txt = '\n'.join('\t'.join(map(str,x)) for x in data)
    pt.write_output(txt, file_name)


def main():   
    in_folder = sys.argv[1]
    out_dir = sys.argv[2]
    nickname = sys.argv[3]
    cutoffs = sys.argv[4]
    num_iter = int(sys.argv[5])
    
    cutoffs = [float(x) for x in cutoffs.split(',')]
    cutoffs = sorted([1.01 if x == 1 else x for x in cutoffs])

    poss_mat, col_headings = pt.get_pres_abs_mat(in_folder)
    results = simulate_reordering(poss_mat, col_headings, cutoffs, num_iter)
    
    cutoffs = [min(x, 1.0) for x in cutoffs]

    for result, cutoff in zip(results, cutoffs):
        file_name = out_dir + '/' + nickname + '_' + str(cutoff) + '.Rtab'
        make_rtab(result, file_name)


if __name__ == "__main__":
    main()
