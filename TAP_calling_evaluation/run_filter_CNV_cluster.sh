#!/usr/bin/env bin

# This function   filters a given popTAP output file  based on 0.05 pvalue and 1/10 of the calculated median of the samnple. 

filter_output() {
  input_file=$1
  #MEAN
  echo calculating overall median
  temp_median=$(zcat ${input_file} | tail -n +2 | awk ' { a[i++]=$5; } END { print a[int(i/2)]; }')
  echo DONE
  #GET NAME OUTPUT
  dirname=$(dirname ${input_file})
  name=$(basename $input_file .depthGCcorrected.per10_bp.bed.gz)
  #RUN THE FILTERING:
  zcat ${dirname}/${name}.depthGCcorrected.per10_bp.bed.output.tsv.gz |
  awk -v median=${temp_median} '{if(NR==1) print $0 ; else if($6 <= 0.05 &&  $3 < median/10 ) print $0 }' | bgzip > ${name}_final_TAPS.tsv.gz
 }
export -f filter_output

for pop in CN DE FR NL SE AM US ; do 
  cd /tmp/global2/acontreras/Tarvense_TE_pops/custom_CNV_output/${pop}_output/
  parallel --jobs 24  filter_output {} ::: *.depthGCcorrected.per10_bp.bed.gz
done
