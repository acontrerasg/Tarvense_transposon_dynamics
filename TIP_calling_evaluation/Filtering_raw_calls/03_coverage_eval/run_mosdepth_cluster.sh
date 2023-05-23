#$ -l h_vmem=8G  # amount of memory
#$ -pe parallel 64 # number of core needed
#$ -S /bin/bash
#$ -cwd
#$ -e /tmp/global2/acontreras/Tarvense_TE_pops/splitreader_bowtie_bams
## -o /tmp/global2/acontreras/Tarvense_TE_pops/splitreader_bowtie_bams
#$ -m beas # send status as email
#$ -N  mosdepth_run #job name

for pop in AM US CN DE FR NL SE ; do

parallel 'mosdepth -n {.} {}' ::: /tmp/global2/acontreras/Tarvense_TE_pops/splitreader_bowtie_bams/${pop}_bams/*bam
parallel 'mosdepth -m -b 10000 -n -t 8  {.}_by10kb {}' ::: /tmp/global2/acontreras/Tarvense_TE_pops/splitreader_bowtie_bams/${pop}_bams/*bam

done


#Unite the summary files:


#Stich together and format the stats:
for summary in  /tmp/global2/acontreras/Tarvense_TE_pops/splitreader_bowtie_bams/*_bams/*.mosdepth.summary.txt   ;
do
name=$(basename $summary .mosdepth.summary.txt) ;
tail -n +2  $summary | awk -v var=$name '{print $0"\t"var}'  ;
done | sed 's/^Chr/Scaffold/' | grep -v 'total' | grep -P "Scaffold_[1-7]\t" > /ebio/abt6_projects8/Tarvense_TE/doc/mosdepth_bam_summary.txt



#Stich together and format the stats:
for i in /tmp/global2/acontreras/Tarvense_TE_pops/splitreader_bowtie_bams/*_bams/*global.dist.txt ;
do
name=$(basename $i .mosdepth.global.dist.txt) ;
tail -n +2  $i | awk -v var=$name '{print $0"\t"var}'  ;
done |  grep "^total" | sed 's/^Chr/Scaffold/' > /ebio/abt6_projects8/Tarvense_TE/doc/mosdepth_bam_global.txt
