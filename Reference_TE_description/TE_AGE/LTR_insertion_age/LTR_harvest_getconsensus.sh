#!/bin/bash
#Retrieve the 3p and 5p sequence frm the output of LTR harvest with bedtools getfasta, then using mafft make the aligment between both ends and 
#make consensus fasta from it. 

genome="/ebio/abt6_projects8/Tarvense_TE/data/Tarvense_TE_pops/Thlaspi_genome_2021/thlaspi_genome.fasta"

ltr_harvest_output="/ebio/abt6_projects8/Tarvense_TE/data/Tarvense_TE_pops/Thlaspi_genome_2021/EDTA_thlaspi_2021/thlaspi_genome.fasta.mod.EDTA.raw/LTR/thlaspi_genome.fasta.mod.harvest.combine.scn"


#Scn format:
#scn is always in + sense
#s(ret) e(ret) l(ret) s(lLTR) e(lLTR) l(lLTR) s(rLTR) e(rLTR) l(rLTR) sim(LTRs) seq-nr chr
# where:
# s = starting position
# e = ending position
# l = length
# ret = LTR-retrotransposon
# lLTR = left LTR
# rLTR = right LTR
# sim = similarity #I am not going to use the similarity produced by ltr harvest as it is based as some empirical value for rice ltrs.
# seq-nr = sequence order
# and each line is a ltr 

#Create dirs to store results
mkdir -p LTR_ends_fasta/ aligned_ltrs/  consensus_ltrs/

#MAIN
grep -v '^#' $ltr_harvest_output  | while read -r line ;
do 
LTR_ID=$( echo ${line} | awk 'BEGIN{OFS="\t"} {print  "Parent_LTRpos#"$12"_"$1"-"$2}')
echo ${line} | awk -v name=$LTR_ID 'BEGIN{OFS="\t"} {print  $12,$4,$5,name";5pLTR",".","+"}' |
bedtools getfasta -fi ${genome} -bed -  -name >  LTR_ends_fasta/${LTR_ID}.fa ;
echo ${line} | awk -v name=$LTR_ID  'BEGIN{OFS="\t"} {print  $12,$7,$8,name";3pLTR",".","+"}' |
bedtools getfasta -fi ${genome} -bed -  -name >> LTR_ends_fasta/${LTR_ID}.fa ;
mafft LTR_ends_fasta/${LTR_ID}.fa > aligned_ltrs/${LTR_ID}.aln.fa 2> /dev/null ;
consambig -name ${LTR_ID} -sequence aligned_ltrs/${LTR_ID}.aln.fa -outseq consensus_ltrs/${LTR_ID}.cons.fa 2> /dev/null ;
done

#Create a tabel with LTR_harvest age prediciton values and  the new TE_IDS:

grep -v '^#' $ltr_harvest_output  | while read -r line ;
do 
LTR_ID=$( echo ${line} | awk 'BEGIN{OFS="\t"} {print  "Parent_LTRpos#"$12"_"$1"-"$2}')
echo ${line} | awk -v name=$LTR_ID 'BEGIN{OFS="\t"} {print  name,$10}' ;
done  | sed '1 i\LTR_ID\tAge' > LTR_harvest_age_predict.tsv


