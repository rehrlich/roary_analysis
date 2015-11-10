#!/usr/bin/env python2

# Author:  Rachel Ehrlich
# This program takes as input a directory of roary output directories and an
# output file name.  The roary directories are named with numbers corresponding
# to their blastP percentages.  The output file is a pdf which plots the
# number of core, soft core, shell and cloud genes for each blastP value.
# Example usage:
# python plot_blastp_comparison.py /home/rachel/Data/pg/roary_except_MB1843_no_delete cluster_freq_blastp
# The above command creates HFlu_blastP_plot.pdf in your current directory


import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import os
import sys

# Input is a line from roary's summary_statistics.txt file
# Output is the gene count for that line
def get_num(data):
    num = data.split(":")[1].strip()
    return int(num)

# Input is a folder and a file name
# Output is the number of lines (gene clusters) in the files
def get_num_lines(direc, file_name):
    num_lines = 0
    with open(direc + file_name, 'rU') as f:
        for line in f:
            num_lines += 1
    return num_lines

# Input is the directory that contains the folders of roary output
# Outputs lists of counts for core, soft, shell, cloud, unsplit, split groups
# and the folder names (blastp)
def get_all_summary_stats(roary_output_dir):
    core = []
    soft = []
    shell = []
    cloud = []
    unsplit = []
    split = []
    total = []
    blastp = []

    for folder in os.listdir(roary_output_dir):
        if not os.path.isdir(roary_output_dir + '/' + folder):
            continue
        direc = roary_output_dir + '/' + folder

        with open(direc + "/summary_statistics.txt", 'rU') as f:
            file_contents = f.read().split('\n')
            
        core.append(get_num(file_contents[0]))
        soft.append(get_num(file_contents[1]))
        shell.append(get_num(file_contents[2]))
        cloud.append(get_num(file_contents[3]))
        total.append(get_num(file_contents[4]))
        blastp.append(int(folder))

        unsplit.append(get_num_lines(direc, "/_inflated_unsplit_mcl_groups"))
        split.append(get_num_lines(direc, "/_inflated_mcl_groups"))
        assert(split == total)
        
    return (core, soft, shell, cloud, unsplit, split, blastp)
    

# Inputs are the gene counts and the output file name
# This plots the counts for each gene group vs the blastP percentage
def plot_counts((core, soft, shell, cloud, unsplit, split, blastp), plot_name):
    with PdfPages(plot_name + '.pdf') as pdf:
        plt.plot(blastp, core, 'ro', label='core (99-100% of strains)')
        plt.plot(blastp, soft, 'bs', label='soft core (95-99% of strains)')
        plt.plot(blastp, shell, 'gd', label='shell (15-95% of strains)')
        plt.plot(blastp, cloud, 'm<', label='cloud (0-15% of strains)')
        plt.plot(blastp, unsplit, 'kD', label='unsplit paralogs')
        plt.plot(blastp, split, 'y*', label='split paralogs')
        
        plt.xlabel('BlastP percent identity')
        plt.ylabel('Number of clusters')
        plt.title('Gene cluster frequency in sequenced strains')
        plt.legend(loc=0, numpoints=1)
        pdf.savefig()
        plt.close()  

def main():
    roary_output_dir = sys.argv[1]
    plot_name = sys.argv[2]

    summary_counts = get_all_summary_stats(roary_output_dir)
    plot_counts(summary_counts, plot_name)

if __name__ == "__main__":
    main()
