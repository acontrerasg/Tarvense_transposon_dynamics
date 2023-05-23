#!/bin/bash
# This script realtes the output IDs od LTRpred with the TE IDs of the original annotation. 
LTRpred_gff="/ebio/abt6_projects8/Tarvense_TE/data/Tarvense_LTRpred/modified_LTRpred.gff"
TEanno="/ebio/abt6_projects8/Tarvense_TE/data/Tarvense_TE_pops/Thlaspi_genome_2021/Thlaspi_paper_final_genome/sorted.final.TEs_lifted.gff"

# Retrieve only Retrotransposns bgger than 4kb (as minimun size in LTRpred is ~4kb)

awk '{if($3 ~ /RL./ && $5-$4 > 4000) print $0 }' ${TEanno} > TEanno_LTRs_only.gff
 
bedtools intersect -f 0.8 -wb -a TEanno_LTRs_only.gff -b ${LTRpred_gff} | 
cut -f9,18 |
uniq |
cut -f1-1 -d ';' |
sed 's/;/\t/' |
sed 's/ID=//g;s/Target=//' > relation_TEanno_ltrpred.tsv
