#!/bin/bash

#Relate the LTR_predicted age LTR parent positions to the original TE names  present in the annotation:


#Create a LTR common annoation from the disjoined TE anno file:

cat disjoined/RLG.gff3 disjoined/RLA.gff3 disjoined/RLC.gff3 | sort -k1,1V -k2,2n -k3,3n  > LTR_disjoined.gff3
cat denested/RLG.gff3 denested/RLA.gff3 denested/RLC.gff3  | sort -k1,1V -k2,2n -k3,3n > LTR_denested.gff3

#Get the final file  from the LTR_age.R  ->  LTRs_tarvense_ages.txt

#Recorver coordinates from the  Parents ID name:

cat  LTRs_tarvense_ages.txt  | sed 's/Parent_LTRpos#//;s/_/\t/2;s/-/\t/' |
awk 'BEGIN{OFS=FS="\t"}{print $2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$1}'  |
sort -k1,1V -k2,2n -k3,3n >  LTRs_tarvense_ages.bed

#we deciend to go with denested for simplicity:
bedtools  intersect  -wao -a LTRs_tarvense_ages.bed  -b LTR_denested.gff3  |
awk 'BEGIN{OFS=FS="\t"}{if ($20 != ".") print $11,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$20,$21/($3-$2)}' > LTRs_tarvense_ages_denested_overlap.tsv

#Now  one can use the row key in the  1st fiel of the final file "LTRs_tarvense_ages_denested_overlap.tsv"  to asign a unique stimation value to the TE faimily with
#the most amount of overlap in the 13th field.
#This can be done in R with:
#group %>% group_by(row_index) %>% top_n(1, overlap)
