#!/bin/usr/env bash

# This script makes a general matrix of TIPS from the insertion coverage corrected file.
# 0 No insertion ; 1 insertion.

cd /ebio/abt6_projects8/Tarvense_TE/data/Tarvense_TE_pops/FINAL_ANALYSIS/TIPS/MATRIX/

TIPS="/ebio/abt6_projects8/Tarvense_TE/data/Tarvense_TE_pops/FINAL_ANALYSIS/TIPS/ALL_TIPS_coverage_corrected.tsv"
genome_index="/ebio/abt6_projects8/Tarvense_TE/data/Tarvense_TE_pops/Thlaspi_genome_2021/Thlaspi_paper_final_genome/modified.fasta.fai"
METADATA="/ebio/abt6_projects8/Tarvense_TE/doc/METADATA_modified_NONL_NOSE.tsv"

tmp_body=$(mktemp body_temp.XXXXXX)
tmp_positions=$(mktemp positions_temp.XXXXXX)
tmp_sample_list=$(mktemp sample_list_temp.XXXXXX)
tmp_tmp=$(mktemp tmp_tmp.XXXXXX)
tmp_region_list=$(mktemp region_list_temp.XXXXXX)

# Get all the possible positions in the file:

tail -n +2 ${TIPS}  |
cut -f1,2,3  |
sort -V -k1,1 -k2,2  |
bedtools  merge -i - -d -1  >  ${tmp_positions}



# MAIN MATRIX ALL SAMPLES:

# Get the list of all samples:
tail -n +2 ${TIPS} | cut -f11 | sort | uniq > ${tmp_sample_list}



for sample in $(cat  ${tmp_sample_list}) ; do
  awk  -v sample=$sample 'BEGIN{OFS=FS="\t"} {if( $11 == sample ) print $1,$2,$3,"Ins="$4}' ${TIPS} |
  sort -V -k1,1 -k2,2  |
  bedtools  merge -i - -d -1  -c 4 -o distinct |
  bedtools intersect -wao -F 0.5 -f 0.5 -a ${tmp_positions} -b - |
  awk '{if($7 == ".") print "0" ; else print "1"}' |
  sed '1 i\'${sample}' '  | paste ${tmp_body} - > ${tmp_tmp} ;
  mv ${tmp_tmp} ${tmp_body} ;
done



sed -i '1 i\#Chrom\tStart\tEnd' ${tmp_positions}

paste ${tmp_positions} ${tmp_body} > ALL_TIPS.matrix


# MATRIXES PER REGION

cat ${METADATA} | tail -n +2 | cut -f1,3 > ${tmp_region_list}

for region in $(cut -f2 ${tmp_region_list} | sort | uniq)  ; do
  # Initate a clean body
  touch  ${region}_tmp_body.bed
  for sample in $(cat  ${tmp_region_list} | grep ${region} | cut -f1) ; do
    awk  -v sample=$sample 'BEGIN{OFS=FS="\t"} {if( $11 == sample ) print $1,$2,$3,"Ins="$4}' ${TIPS} |
    sort -V -k1,1 -k2,2  |
    bedtools  merge -i - -d -1  -c 4 -o distinct |
    bedtools intersect -wao -F 0.5 -f 0.5 -a ${tmp_positions} -b - |
    awk '{if($7 == ".") print "0" ; else print "1"}' |
    sed '1 i\'${sample}' '  | paste ${region}_tmp_body.bed  - > ${tmp_tmp} ;
    mv ${tmp_tmp} ${region}_tmp_body.bed  ;
  done
  paste ${tmp_positions}  ${region}_tmp_body.bed  > ${region}_TIPS.matrix
  rm ${region}_tmp_body.bed ;
done

#Cleanup

rm ${tmp_positions}
rm ${tmp_body}
rm ${tmp_sample_list}
rm ${tmp_region_list}
