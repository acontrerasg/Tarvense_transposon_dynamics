#!/usr/bin/env bash

usage() { echo "$0 Usage:
This script interacts over every bed file of the final call folder of the SPL raw filtering.
It  does two things. First it removes the DHH type TE calls. Second, it removes any TE call outside of the
accesibility regions file. For this, it swtiches back from the Chr_ notation splitreader uses to  the one
the reference genome has (Scaffold_). 
It also fixes the problem of files ending up with an extra tab at the end and also a extra space at the header." && grep " .)\ #" $0; exit 0; }
[ $# -eq 0 ] && usage
while getopts ":ha:o:i:" arg; do
  case $arg in
    a) #  accesibility bed file
      acc_file=${OPTARG}
      ;;
    i) # # path/to/final_call_folder/
      Input_dir=${OPTARG}
      ;;
    o) # # path/to/output
      output_dir=${OPTARG}
     ;;
    h | *) # Display help.
      usage
      exit 0
      ;;
  esac
done


output_dir="/ebio/abt6_projects8/Tarvense_TE/data/Tarvense_TE_pops/Splitreader/RESULTS/Results_filtering/acc_filtering"
input_dir="/ebio/abt6_projects8/Tarvense_TE/data/Tarvense_TE_pops/Splitreader/RESULTS/Results_filtering/baduel_filtering/final_call"
acc_file="/ebio/abt6_projects8/Tarvense_TE/data/Tarvense_TE_pops/Splitreader/RESULTS/Results_filtering/Regions_filtered_3Medians_merged.bed"


for bed_file in ${input_dir}/*bed ; 
  do
   name=$(basename $bed_file .bed) ;
   head -n 1  $bed_file | sed 's/\t$//;s/ //g' > header.txt ;
   sed 's/Chr_/Scaffold_/' ${bed_file} | tail -n +2 |
   sed 's/\t$//' |
   awk 'BEGIN{OFS=FS="\t"} {if($4 !~ /DHH/) print $0}' |
   bedtools intersect -wa -f 1 -a - -b ${acc_file} | 
   cat header.txt - > ${output_dir}/${name}_acc.bed ;
   rm header.txt ; 
  done 
