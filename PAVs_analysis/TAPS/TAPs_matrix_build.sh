#!/bin/usr/env bash
lifted_anno_disjoined="/ebio/abt6_projects8/Tarvense_TE/data/Tarvense_TE_pops/Thlaspi_genome_2021/Thlaspi_paper_final_genome/sorted.final.TEs_lifted.gff"


mkdir -p /ebio/abt6_projects8/Tarvense_TE/data/Tarvense_TE_pops/FINAL_ANALYSIS/TAPS/MATRIX/
cd /ebio/abt6_projects8/Tarvense_TE/data/Tarvense_TE_pops/FINAL_ANALYSIS/TAPS/MATRIX/

#Get the > 300 bp TE file (we skip TEs smaller than 300 becasue we dont evaluate its absences)
awk 'BEGIN{OFS=FS="\t"}{ if($5-$4 > 299) print $0}' ${lifted_anno_disjoined} |
sed "1 i\Chrom\tSource\tFeature\tStart\tEnd\tScore\tStrand\tPhase\tAttributes" > disjoin_anno_300.gff

#  Iteare over the input file and add the  annotation file and add to column 10 a -1 if there is an absence and a 0 if its  present.

touch body_temp_XXXXXXXXXX
for TAP_file in /ebio/abt6_projects8/Tarvense_TE/data/Tarvense_TE_pops/FINAL_ANALYSIS/TAPS/individual_TAPs/*_ind.gff3 ; do
    #wc -l ${TAP_file} ;
    name=$(basename ${TAP_file} _ind.gff3) ;
    sed 's/;Sample=.*$//;' ${TAP_file} > temp_tap ;
    awk 'NR==FNR {h[$9] = $0; next} {if(h[$9] == $0) print h[$9]"\t-1" ; else print $0"\t0" }' temp_tap  disjoin_anno_300.gff |
    sed "1 i\Chrom\tSource\tFeature\tStart\tEnd\tScore\tStrand\tPhase\tAttributes\t${name}"  > ${name}_binary.gff ;
    cut -f10 ${name}_binary.gff | paste - body_temp_XXXXXXXXXX  > tmp ;
    mv tmp body_temp_XXXXXXXXXX ;
    awk '{print NF}' body_temp_XXXXXXXXXX  | sort | uniq -c  ;
    rm  temp_tap ;
done

paste disjoin_anno_300.gff body_temp_XXXXXXXXXX  |
sed '$ d' | sed 's/\t$//' > final_TAP_matrix.tsv

rm body_temp_XXXXXXXXXX
echo "DONE"
