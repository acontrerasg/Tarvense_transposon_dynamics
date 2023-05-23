#$ -l h_vmem=4G  # amount of memory
#$ -pe parallel 120 # number of core needed
#$ -S /bin/bash
#$ -cwd
#$ -e /ebio/abt6_projects8/Tarvense_TE/data/Tarvense_ancestry/treemix/boostrap_runs/$JOB_NAME.$JOB_ID.ER # where to write stderr
#$ -o /ebio/abt6_projects8/Tarvense_TE/data/Tarvense_ancestry/treemix/boostrap_runs/$JOB_NAME.$JOB_ID.OU # where to write stdout
#$ -m beas # send status as email
#$ -N treeemix_bs_2.5k #job name

# Functions
run_treemix() {
                a=$RANDOM
                b=1$(date +%N)
                c=$(echo $a + $b | bc)
                treemix -i $2 -bootstrap -k $3 -se -root $4 -seed $c -o ${5}"treemix_boostrap_"$1
}
export -f run_treemix

# Inputs

treemix_file="/ebio/abt6_projects8/Tarvense_TE/data/Tarvense_ancestry/treemix/treemix.strat.frq.gz"
boostraps=2500
threads=120
k=100
outgroup="outgroup"
outdir="/ebio/abt6_projects8/Tarvense_TE/data/Tarvense_ancestry/treemix/boostrap_runs/"

seq 1 $boostraps | parallel -j $threads run_treemix {} $treemix_file  $k  $outgroup $outdir >> "logfile_treemix"$boostraps"_boostrap.log"
