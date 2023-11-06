#!/usr/bin/env bash

# This script maps illumina short reads geenrated for polishing: 88473_CGATGT_CDMKUANXX and 88473_CGATGT_CDMKUANXX

genome="/ebio/abt6_projects8/Tarvense_TE/data/Tarvense_TE_pops/Thlaspi_genome_2021/Thlaspi_paper_final_genome/modified.fasta"


for read in /ebio/abt6_projects8/Tarvense_TE/data/Tarvense_TE_pops/Splitreader/reference_reads_MN106/fastq/*1.fq.gz ; do
read_name=`echo $read | sed 's/_1.fq.gz//'`
out_name=`basename $read_name`
bwa mem -c 1 -t 8 $genome ${read_name}_1.fq.gz ${read_name}_2.fq.gz  | samtools sort -@ 8 -o ${out_name}.bam - ;
done
