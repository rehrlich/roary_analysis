## Scripts for analyzing the output of Roary which can be found at:  https://github.com/sanger-pathogens/Roary
Suggested use:  run plot_blastp_comparison.py, look at the graph to choose blastp values of interest and then run analyze_blastp_raory.sh


# plot_blastp_comparison.py

This program takes as input a directory of roary output directories and an
output file name.  The roary directories are named with numbers corresponding
to their blastP percentages.  The output file is a pdf which plots the
number of core, soft core, shell and cloud genes for each blastP value.

Example usage:

python plot_blastp_comparison.py /home/rachel/Data/HFluGenomes/roary_output HFlu_blastP_plot

The above command creates HFlu_blastP_plot.pdf in your current directory

# analyze_blastp_raory.sh

This program takes as input the folder with the fsgm .m files, the folder
containing all the roary output, a nickname for the analysis, an output
directory and the blastp roary folder of interest.
This creates outdir containing output from gene_counts_heat_map.r,
pairwise_table.py and the fsgm cgs_supragenome program
This assumes the roary output is structured:  roary_output/blastp to allow
for simple iteration over multiple folders.  If not, use the roary parent for
the roary_output and the roary folder for blastp.

# pairwise_table.py

This program takes as input the folder with the output from roary, an output
folder and a nickname for the results.  It creates five files in the
output directory with the prefix the supplied nickname with the pairwise
tables and summary statistics.
In the pairwise table files, each number in the matrix refers to the gene
counts for the corresponding strains.

* similarity = present in both strains (includes core)
* difference = present in exactly one of the strains (xor)
* comparison = similarity - difference
* pair unique = present in only those two stains

# get_fsgm_input.py
This program takes as input the folder with the output from roary and writes
to standard out the input for cgs_supragenome.m which can be found at https://github.com/rehrlich/fsgm

# gene_counts_heat_map.r
This program makes heat maps from Roary's output gene counts.  Command line inputs are the roary output folder, nickname for the run, and the output directory.

Example usage:

Rscript gene_counts_heat_map.r roary_dir nickname outdir

# plot_fsgm_results.py
This program takes as input an output directory (where it expects to find fsgm
output with the same nickname) and a nickname for the analysis.  It creates
two plots of the data, the likelihood of various N values and the number of
expected new genes per strain sequenced.

# simulate_pan_genome.py
This program takes as input a folder with roary output, an output directory
for the results, a nickname for the analysis, a comma separated string of 
cutoff values where 0 < value <= 1 and the number of simulations to run.
An example cutoff argument is "0.15,0.95,0.99,1.0"
This does the specified number of simulations to model how the gene 
frequencies change when strains are sequenced in a different order.
For each of the cutoffs, a tab separated .Rtab (for consistency 
with roary) has one line for each simulation.  The nth number in a line
is the number of genes whose frequency in the first n strains sampled
is next smallest cutoff <= freqnecy < cutoff.  For the case of the 
lowest cutoff, the lower number is 0.  For the case of 1.0, the upper 
range is <= to include 100%.
