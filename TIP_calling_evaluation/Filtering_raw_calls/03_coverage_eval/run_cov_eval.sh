#!/usr/bin/env bash

cd /ebio/abt6_projects8/Tarvense_TE/data/Tarvense_TE_pops/Splitreader/RESULTS/Results_filtering/03_coverage_eval/


#fixed paths
coverage_eval_insertions_script="./coverage_eval_insertions.sh"
final_NC_file="/ebio/abt6_projects8/Tarvense_TE/data/Tarvense_TE_pops/Splitreader/RESULTS/Results_filtering/02_acc_filtering/lifted_V4_V5/Tarvense-insertions.FINAL.NC.DP3_acc_lifted.bed"
path_mosdepth_results="/tmp/global2/acontreras/Tarvense_TE_pops/splitreader_bowtie_bams/mosdepth_results_all"
file_path="/ebio/abt6_projects8/Tarvense_TE/data/Tarvense_TE_pops/Splitreader/RESULTS/Results_filtering/02_acc_filtering/results_individual/lifted_files"

#Descrive function
cov_eval_ins() {
 coverage_eval_insertions=$1
 final_NC_file=$2
 individual_file=$3
 mosdepth_10kb_median_file=$4
 bash ${coverage_eval_insertions}  -f ${final_NC_file} -i ${individual_file} -r ${mosdepth_10kb_median_file}
}

export -f cov_eval_ins

#run in parallel:

find ${file_path} -name "*.bed" |
parallel --jobs 10  cov_eval_ins ${coverage_eval_insertions_script} ${final_NC_file} {} ${path_mosdepth_results}/{/.}_by10kb.regions.bed.gz
