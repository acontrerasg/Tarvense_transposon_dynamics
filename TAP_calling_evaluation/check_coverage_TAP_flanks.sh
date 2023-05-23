#!/bin/bash

# This script checks the flanking region of the TAPS of a given sample

# Main fun to run in parallel
retrieve_cov_flanks() {
                # Input
                sample_name=$1
                sample_gff=$2
                mosdepth_filepath=$3
                genome_fai=$4
                #MAIN
                # Expand flanking region 1kb both ends
                five_prime_temp=$(mktemp tempfileXXXX --suffix=.gff )
                bedtools flank -i ${sample_gff} -g $genome_fai  -l 1000 -r 0 > ${five_prime_temp}
                three_prime_temp=$(mktemp tempfileXXXX --suffix=.gff )
                bedtools flank -i ${sample_gff} -g $genome_fai  -l 0 -r 1000 > ${three_prime_temp}
                # Check coverage
                five_prime_out=$(mktemp tempfileXXXX --suffix=.tsv )
                three_prime_out=$(mktemp tempfileXXXX --suffix=.tsv )
                #Get coverage file
                cov_file=$(echo ${mosdepth_filepath}${sample_name}_per-base.bed.gz)
                # 5 prime
                zcat ${cov_file} | bedtools intersect -wb -a ${five_prime_temp} -b - |
                bedtools groupby -i - -g 9 -c 13 -o median |
                awk '{print $1"\t;5prime_cov="$2}' > ${five_prime_out}
                # 3 prime
                zcat ${cov_file} | bedtools intersect -wb -a ${three_prime_temp} -b - |
                bedtools groupby -i - -g 9 -c 13 -o median |
                awk '{print $1"\t;3prime_cov="$2}' > ${three_prime_out}
                #Merge outputs
                awk -F '\t' 'NR==FNR{a[$1]=$2;next} ($1 in a){b=$1;print $0,a[b] }'  OFS="\t" ${three_prime_out} ${five_prime_out} |
                sed 's/\t/;/g' >  ${sample_name}_1kb_flanking_median_cov.tsv
                # Cleanup
                rm ${five_prime_temp}
                rm ${three_prime_temp}
                rm ${five_prime_out}
                rm ${three_prime_out}

}

export -f retrieve_cov_flanks


# Variables

genome_fai="/ebio/abt6_projects8/Tarvense_TE/data/Tarvense_TE_pops/Thlaspi_genome_2021/Thlaspi_paper_final_genome/modified.fasta.fai"
mosdepth_filepath="/tmp/global2/acontreras/Tarvense_TE_pops/bwa_mem_aligments/mosdepth_per_base/"
sample_list="/ebio/abt6_projects8/Tarvense_TE/data/Tarvense_TE_pops/FINAL_ANALYSIS/TAPS/TAPS_flanking_cov/sample_names.tsv"
gff_path="/ebio/abt6_projects8/Tarvense_TE/data/Tarvense_TE_pops/FINAL_ANALYSIS/TAPS/individual_TAPs/"
threads=10

# Prepare files and move to folder
cd /ebio/abt6_projects8/Tarvense_TE/data/Tarvense_TE_pops/FINAL_ANALYSIS/TAPS/TAPS_flanking_cov
ln -s /ebio/abt6_projects8/Tarvense_TE/data/Tarvense_TE_pops/FINAL_ANALYSIS/TAPS/ALL_TAPS_sorted.gff  .
cat ALL_TAPS_sorted.gff  | cut -f9 | cut -f3 -d ';'  | sort | uniq  > sample_names.tsv 

# Run in parallel 
parallel --jobs $threads -a $sample_list "retrieve_cov_flanks {} ${gff_path}{}_ind.gff3 ${mosdepth_filepath} ${genome_fai}"
