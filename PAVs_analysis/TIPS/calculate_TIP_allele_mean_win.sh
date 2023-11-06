#!/usr/bin/env bash

# Reviewer recoomendaiton: 
# Outliers like the few almost fixed TE insertions near chromosome ends can create the impression of an overall trend. 
# It may be more informative to display the average allele frequency by genomic windows (e.g. 100/500Kb windows) instead. 

# Create the windowed allele frequency file from the already calculated "ALL_TIPS_allele.freq.bedgraph"

# Files:
TIP_allele_freq="/ebio/abt6_projects8/Tarvense_TE/data/Tarvense_TE_pops/FINAL_ANALYSIS/TIPS/allele_freqs_bg/ALL_TIPS_allele.freq.bedgraph"
genome_fai="/ebio/abt6_projects8/Tarvense_TE/data/Tarvense_TE_pops/Thlaspi_genome_2021/Thlaspi_paper_final_genome/modified.fasta.fai"
# Path 
cd /ebio/abt6_projects8/Tarvense_TE/data/Tarvense_TE_pops/FINAL_ANALYSIS/TIPS/allele_freqs_bg/

### Calculate the GENOMIC STEPS OF 100 JKB:
# Divide the genome in windows of 100 kb.
bedtools makewindows -g ${genome_fai} -w 100000 > genome_windowed100kb.bed

bedtools intersect -wo -a genome_windowed100kb.bed -b ${TIP_allele_freq} |
bedtools groupby -g 1,2,3 -c 7 -o count,mean >  ALL_TIPS_allele.freq.windowed.tsv 
sed -i '1 i\#Chrom\tStart\tEnd\tWindowID\tallele_count\tMean_allele_freq' ALL_TIPS_allele.freq.windowed.tsv 

# Prepare the count bedgraph and the mean frequency bedgraph: 
## counts
cut -f1,2,3,4 ALL_TIPS_allele.freq.windowed.tsv |
tail -n +2 > ALL_TIPS_allele.freq.windowed.counts.bedgraph
sed -i '1 i\track type=bedGraph name="Counts TIP allele per 100kb window " ' ALL_TIPS_allele.freq.windowed.counts.bedgraph
## means
cut -f1,2,3,5 ALL_TIPS_allele.freq.windowed.tsv |
tail -n +2 > ALL_TIPS_allele.freq.windowed.mean_freqs.bedgraph
sed -i '1 i\track type=bedGraph name="Mean freqs TIP allele per 100kb window" ' ALL_TIPS_allele.freq.windowed.mean_freqs.bedgraph
