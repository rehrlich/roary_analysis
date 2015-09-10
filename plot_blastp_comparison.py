# Author:  Rachel Ehrlich
# This program takes as input a directory of roary output directories and an
# output file name.  The roary directories are named with numbers corresponding
# to their blastP percentages.  The output file is a pdf which plots the
# number of core, soft core, shell and cloud genes for each blastP value.
# Example usage:
# python plot_blastp_comparison.py /home/rachel/Data/HFluGenomes/roary_output HFlu_blastP_plot
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

# Input is the directory that contains the folders of roary output
# Outputs lists of counts for core, soft, shell, cloud
# and the folder names (blastp)
def get_all_summary_stats(roary_output_dir):
    core = []
    soft = []
    shell = []
    cloud = []
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
        blastp.append(int(folder))
        
    return (core, soft, shell, cloud, blastp)

# Inputs are the gene counts and the output file name
# This plots the counts for each gene group vs the blastP percentage
def plot_counts((core, soft, shell, cloud, blastp), plot_name):
    with PdfPages(plot_name + '.pdf') as pdf:
        plt.plot(blastp, core, 'ro', label='core (99-100)')
        plt.plot(blastp, soft, 'bs', label='soft core (95-99)')
        plt.plot(blastp, shell, 'g^', label='shell (15-95)')
        plt.plot(blastp, cloud, 'm<', label='cloud (0-15)')
        plt.xlabel('BlastP %')
        plt.ylabel('Number of genes')
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
