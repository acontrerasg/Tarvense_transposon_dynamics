##This script  first gets the .depthGCcorrected.regions  file from the raw reads and then runs the custom CNV pipeline in  the CN pops.
#hen it  appies a test statisct to compare them with the stats in the reference file for the TE annotation.

##Scripts_paths
get_cov="/ebio/abt6_projects8/Tarvense_TE/code/custom_cnv_pipeline/get_cov_sample.sh"
genome_file="/ebio/abt6_projects8/Tarvense_TE/data/Tarvense_TE_pops/Thlaspi_genome_2021/thlaspi_genome.fasta"
TE_CNV="/ebio/abt6_projects8/Tarvense_TE/code/custom_cnv_pipeline/TE_CNV.sh"
testR="/ebio/abt6_projects8/Tarvense_TE/code/custom_cnv_pipeline/Wilcoxon.test.R"
reference_file="/ebio/abt6_projects8/Tarvense_TE/data/Tarvense_TE_pops/CNV/RCcorrected_reference/8853_8874_avg_depthGCcorrected.per10_bp.bed.gz"
TE_anno="/ebio/abt6_projects8/Tarvense_TE/data/Tarvense_TE_pops/FINAL_ANALYSIS/Tarvense_acc_anno_files/sorted.final.TEs_denested_covF.gff3"


for pop in CN DE FR NL SE AM US ; do 
  ## First part, map to the gneome and get the GCcorrected.regions  files.
  cd /tmp/global2/acontreras/Tarvense_TE_pops/custom_CNV_output/${pop}_output/
  parallel "bash ${get_cov}  -g ${genome_file} -w 10 -f {} -t 8 " ::: /ebio/abt6_projects8/Tarvense_TE/data/Tarvense_TE_pops/Tarvense_TE_pops_raw_reads/Geng_etal/raw_reads/*1.fastq.gz
  ## Second part, get the statistics  comparing this  with the ref for the TE annnot file.
  #It uses one thread per sample but needs a big chunk of RAM.
  #Run the TE CNV cript in parallel.
  parallel "bash $TE_CNV -i {} -r ${reference_file} -a ${TE_anno} -t 300  -S ${testR}" ::: *.depthGCcorrected.per10_bp.bed.gz
done 
