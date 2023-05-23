#!/bin/bash

# This script aims to annotate
mkdir -p /ebio/abt6_projects8/Tarvense_TE/data/Tarvense_TE_pops/FINAL_ANALYSIS/TIPS/TIPS_closest/
cd /ebio/abt6_projects8/Tarvense_TE/data/Tarvense_TE_pops/FINAL_ANALYSIS/TIPS/TIPS_closest/
ln -s /ebio/abt6_projects8/Tarvense_TE/data/Tarvense_TE_pops/Thlaspi_genome_2021/Thlaspi_paper_final_genome/sorted.TE.genes_lifted.gff . 
ln -s /ebio/abt6_projects8/Tarvense_TE/data/Tarvense_TE_pops/FINAL_ANALYSIS/TIPS/ALL_TIPS_curated.gff3 . 
ln -s /ebio/abt6_projects8/Tarvense_TE/data/Tarvense_TE_pops/Thlaspi_genome_2021/Thlaspi_paper_final_genome/modified.fasta.fai .


# Sort both files:
bedtools  sort -i ALL_TIPS_curated.gff3  -g modified.fasta.fai > sorted_tips.gff
bedtools  sort -i sorted.TE.genes_lifted.gff -g modified.fasta.fai > sorted_anno.gff

bedtools  closest  -d -g modified.fasta.fai -a sorted_tips.gff  -b  sorted_anno.gff -t first |  
cut -f1-9,12,18,19  >  TIPS_close_features.tsv 

sed -i 's/TIP_FAM=//;s/;Accession=/\t/;s/;Region=/\t/;s/ID=//;s/dispersed_repeat.*;Target=//' TIPS_close_features.tsv 

# Merge the closest call with the allele frequencies of TIPS:

allele_fr_bg="/ebio/abt6_projects8/Tarvense_TE/data/Tarvense_TE_pops/FINAL_ANALYSIS/TIPS/allele_freqs_bg/ALL_TIPS_allele.freq.bedgraph"

tail -n +2 ${allele_fr_bg} |
bedtools intersect -wb -a - -b TIPS_close_features.tsv | 
sed 's/TA_/\tTA_/' | 
cut -f4,5,8,9,13- >  TIPS_close_features_freq.tsv
