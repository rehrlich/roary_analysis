Scripts for analyzing the output of Roary which can be found at:  https://github.com/sanger-pathogens/Roary

plot_blastp_comparison.py

This program takes as input a directory of roary output directories and an
output file name.  The roary directories are named with numbers corresponding
to their blastP percentages.  The output file is a pdf which plots the
number of core, soft core, shell and cloud genes for each blastP value.
Example usage:
python plot_blastp_comparison.py /home/rachel/Data/HFluGenomes/roary_output HFlu_blastP_plot
The above command creates HFlu_blastP_plot.pdf in your current directory

pairwise_table.py

This program takes as input the folder with the output from roary, an output
folder and a nickname for the results.  It creates three files in the
output directory with the prefix the supplied nickname with the pairwise
tables and summary statistics.
In the two pairwise table files, each number in the matrix refers to the gene
counts for the corresponding strains.
similarity = present in both strains (includes core)
difference = present in exactly one of the strains (xor)
comparison = similarity - difference
pair unique = present in only those two stains
