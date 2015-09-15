# Author:  Rachel Ehrlich
# This program takes as input the folder with the fsgm .m files, the folder
# containing all the roary output, a nickname for the analysis, an output
# directory and the blastp roary folder of interest.
# This creates outdir containing output from gene_counts_heat_map.r,
# pairwise_table.py and the fsgm cgs_supragenome program
# This assumes the roary output is structured:  roary_output/blastp to allow
# for simple iteration over multiple folders.  If not, use the roary parent for
# the roary_output and the roary folder for blastp.


fsgm=$1
roary_output=$2
nickname=$3
outdir=$4
blastp=$5

mkdir -p $outdir

roary_run=$roary_output"/"$blastp
name=$nickname"_"$blastp

cd `dirname $0`

Rscript gene_counts_heat_map.r $roary_run $name $outdir

python pairwise_table.py $roary_run $outdir $name

gene_counts=`python get_fsgm_input.py $roary_output/$blastp`
cd $fsgm
matlab -nodesktop -nojvm -nosplash -r "cgs_supragenome '$name' '$gene_counts' '$outdir';exit"
