# Author:  Rachel Ehrlich
# This program takes as input a folder with output from roary,
# an output directory that has the results from simulate_pan_genome.py
# and a nickname for the outputs.
# This makes two plots, one from the roary Rtab data and one from
# the simulated gene frequency data.

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import sys
import os

# Input is a list of lines from an rtab file output by roary
# Returns a list where each element is a list of integers from 
# the corresponding line
def rtab_to_list_of_lists(data):
    lol = [x.replace('\n','').split('\t') for x in data]
    lol_ints = [map(int, x) for x in lol]
    return lol_ints

# Input is a folder containing .Rtab files and a list of those files
# Returns a list where each element is a list of lists with rtab data
# from one file
def get_rtab_data(path, files):
    data = []
    for x in files:
        with open(path + '/' + x, 'rU') as f:
            data.append(rtab_to_list_of_lists(f.readlines()))
    return data

# Input is a list with the core, total, new_genes and unique rtab data and
# an output file
# Creates a pdf with a line plot of the data
def make_plots((core, total, new_genes, unique), out_file):
    num_genomes = range(1, len(core[0]) + 1)
    with PdfPages(out_file) as pdf:
        for i in range(len(new_genes)):
            plt.plot(num_genomes, new_genes[i], 'r', label = "new")
            plt.plot(num_genomes, total[i], 'b', label = 'total')
            plt.plot(num_genomes, core[i], 'g', label = 'core')
            plt.plot(num_genomes, unique[i], 'k', label = 'unique')
            
            if i == 0:
                plt.legend(loc=0)
                
        plt.xlabel("Number of genomes")
        plt.ylabel('Number of clusters')
        plt.title("Observed size of pan genome per strain sequenced")
        
        pdf.savefig()
        plt.close()

# Input is a list with the simulated rtab data, a list of labels,
# and an output file
# Creates a pdf with a line plot of the data
def make_plots2(data, labels, out_file):
    num_genomes = range(1, len(data[0][0]) + 1)
    num_simulations = len(data[0])
    num_bins = len(data)
    colors = ['m', 'g', 'b', 'r']

    with PdfPages(out_file) as pdf:
        for sim_num in range(num_simulations):
            all_bin_data = [x[sim_num] for x in data]
            
            for one_bin_data, label, color in zip(all_bin_data, labels, colors):
                plt.plot(num_genomes, one_bin_data, color, label = label)

            if sim_num == 0:
                plt.legend(loc=0)
             
        plt.xlabel("Number of genomes")
        plt.ylabel('Number of genes')
        plt.title("Observed gene frequency per strain sequenced")
        
        pdf.savefig()
        plt.close()

# Input is a directoy containing the output from simulate_pan_genome.py
# Returns a list of the .Rtab files and a list of their cutoffs sorted
# by cutoffs
def get_simulated_files(outdir):
    data = []
    for file1 in os.listdir(outdir):
        if file1.endswith(".Rtab"):
            percent = file1.rsplit('_', 1)[1][:-5]
            data.append((percent, file1))
            
    data.sort(key=lambda x: float(x[0]))
    
    files = [x[1] for x in data]
    cutoffs = [x[0] for x in data]
    return files, cutoffs
               
def main():
    roary_output = sys.argv[1]
    outdir = sys.argv[2]
    nickname = sys.argv[3]
    
    roary_files = ["/number_of_conserved_genes.Rtab",
                   "/number_of_genes_in_pan_genome.Rtab",
                   "/number_of_new_genes.Rtab", "/number_of_unique_genes.Rtab"]
    data = get_rtab_data(roary_output, roary_files)
    make_plots(data, outdir + '/' + nickname + '_observed_genome_size.pdf')
    
    sim_files, cutoffs = get_simulated_files(outdir)

    data = get_rtab_data(outdir, sim_files)
    plot_file =  outdir + '/' + nickname + '_observed_gene_frequencies.pdf'
    make_plots2(data, map(str, cutoffs), plot_file)

if __name__ == "__main__":
    main()
