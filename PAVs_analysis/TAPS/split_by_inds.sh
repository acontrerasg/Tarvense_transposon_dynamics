#!/bin/usr/env bash

mkdir -p /ebio/abt6_projects8/Tarvense_TE/data/Tarvense_TE_pops/FINAL_ANALYSIS/TAPS/individual_TAPs/
cd /ebio/abt6_projects8/Tarvense_TE/data/Tarvense_TE_pops/FINAL_ANALYSIS/TAPS/

for sample in $(tail -n +2 ALL_TAPS_flank_filtered.gff | cut -f 9 | cut -f 3 -d ';' | sort | uniq ) ; do 
        name=$(echo $sample | sed 's/Sample=//' ) ; 
        echo $name
        grep $name  ALL_TAPS_flank_filtered.gff  > /ebio/abt6_projects8/Tarvense_TE/data/Tarvense_TE_pops/FINAL_ANALYSIS/TAPS/individual_TAPs/${name}_ind.gff3 ;
        wc -l /ebio/abt6_projects8/Tarvense_TE/data/Tarvense_TE_pops/FINAL_ANALYSIS/TAPS/individual_TAPs/${name}_ind.gff3 ;
done

echo "DONE!"
