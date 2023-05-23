#!/bin/bash

# Get metadata
metadata="/ebio/abt6_projects8/Tarvense_TE/doc/METADATA.tsv"
# Path to individual insertions

ind_insertions_path="/ebio/abt6_projects8/Tarvense_TE/data/Tarvense_TE_pops/FINAL_ANALYSIS/TIPS/individual_TIPS"
# Set tmp files
tmp_region_info=$(mktemp region_info.XXXXXX)
tmp_inds=$(mktemp inds_temp.XXXXXX)

# MAIN
tail -n +2 ${metadata} | cut -f3 | sort | uniq > ${tmp_region_info}

while read region ;
        do
        grep "${region}" $metadata  | cut -f1 | shuf >  ${tmp_inds}
        counter=0
                for ind in $(cat  $tmp_inds) ;
                        do
                        counter=$((counter + 1))
                        cat ${ind_insertions_path}/${ind}*  |
                        awk  '{if($6==1) print $1"\t"$2"\t"$3}'  >> region_tmpfile.bed
                        bedtools sort -i region_tmpfile.bed | bedtools merge -i - | wc -l | echo -e ${counter}"\t"${region}"\t"$(</dev/stdin) ;
                        done
                rm region_tmpfile.bed
        done  <  ${tmp_region_info} > ./saturation_curve.tsv

#Cleanup

rm ${tmp_region_info}
rm ${tmp_inds}
