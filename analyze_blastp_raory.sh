# Author:  Rachel Ehrlich
# This program takes as input the folder with the fsgm .m files, the folder
# containing all the roary output, a nickname for the analysis, an output
# directory and the blastp roary folder of interest.
# This creates outdir containing output from gene_counts_heat_map.r,
# pairwise_table.py and the fsgm cgs_supragenome program
# This assumes the roary output is structured:  roary_output/blastp to allow
# for simple iteration over multiple folders.  If not, use the roary parent for
# the roary_output and the roary folder for blastp.

#bash analyze_blastp_raory.sh /home/rachel/Documents/Finite_supragenome_model/fsgm/ /home/rachel/Data/pg/roary_except_MB1843 roary_except_MB1843 /home/rachel/Data/pg/roary_except_MB1843/results_matlab_test 85

#bash analyze_blastp_raory.sh /home/rachel/Documents/Finite_supragenome_model/fsgm/ /home/rachel/Documents/roary/roary_analysis roary_except_MB1843 /home/rachel/Documents/roary/roary_analysis

#bash analyze_blastp_raory.sh /home/rachel/Documents/Finite_supragenome_model/fsgm/ /home/rachel/mnt/shared/homes/rachel/borrelia/roary_75 Bbur /home/rachel/mnt/shared/homes/rachel/borrelia/roary_figures 75
fsgm=$1
roary_output=$2
nickname=$3
outdir=$4
blastp=$5

mkdir -p $outdir
roary_run=$roary_output
#roary_run=$roary_output"/"$blastp
name=$nickname"_"$blastp
name=$nickname

DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
cd $DIR


Rscript gene_counts_heat_map.r $roary_output $name $outdir

# So I don't lose data when unsplitting paralogs
cp $roary_run"/summary_statistics.txt"  $outdir"/"$name"_roary_summary_statistics.txt"
cp $roary_run"/gene_presence_absence.csv" $outdir"/"$name"_cluster_presence_absence.csv"

# Simulate adding strains to the pan genome in different orders
cutoffs="0.15,0.95,0.99,1.0" #graph needs more colors to increase this
num_simulations=5
python simulate_pan_genome.py $roary_run $outdir $name $cutoffs $num_simulations
python plot_rtab.py $roary_run $outdir $name


python pairwise_table.py $roary_run $outdir $name


tree=$roary_run'/core_gene_alignment.aln.newick'
Rscript pairwise_heat_map.r $outdir"/" $name $tree
Rscript make_trees.r $outdir $tree

exit 1
gene_counts=`python get_fsgm_input.py $roary_output/$blastp`

cd $fsgm
matlab -nodesktop -nojvm -nosplash -r "cgs_supragenome '$name' '$gene_counts' '$outdir';exit"

cd $DIR
python plot_fsgm_results.py $outdir $name
