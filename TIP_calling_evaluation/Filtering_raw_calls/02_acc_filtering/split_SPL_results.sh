#!/usr/bin/env bash

# First argument is the boolean file obtained after splitreader filtering. 

##### Input arg:
results_file=$1
####

tmp_info=$(mktemp info_temp.XXXXXX)
tmp_ind=$(mktemp ind_temp.XXXXXX)

## Get number of lines in file

col_num=$(awk '{print NF}' ${results_file} | uniq)

## Get the info in the first five colmuns:

cut -f1-5 ${results_file} > ${tmp_info}

## Create folder (if it doesnt exists)

mkdir -p results_individual/


## Loop over every individual and add  the info fields

for i in $(seq 6 $col_num) ; do
        name=$(cut -f ${i}  ${results_file} | head -n 1) ;
        cut -f ${i} ${results_file} > ${tmp_ind} ;
        #by individuals
        paste  ${tmp_info} ${tmp_ind} | 
        awk '{if($6==1) print $0}' > results_individual/${name}.DP3.bed ;

done

# Cleanup
rm ${tmp_info}
rm ${tmp_ind}
