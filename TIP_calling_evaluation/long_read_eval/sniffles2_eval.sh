#!/bin/bash

# Map armenian reads to  Nunn genome and call structural variants. 
# Conda env and location 

cd /ebio/abt6_projects8/Tarvense_TE/data/Armenian_thlaspi/SV_sniffles2/
conda activate /ebio/abt6_projects8/Tarvense_TE/code/conda_envs/sniffles2/

# File with raw  PACBIO reads
ln -s /ebio/seq/lra/runs/runs/64079/r64079_20220629_055159/1_A01/m64079_220629_071805.subreads.bam  .
# File with HiFi reads in fastq format 
ln -s /tmp/global2/amovilli/HiFi_Thlaspi/1_Q20-kinetics/Armenian_Thlaspi.q20.5mC.fastq.gz . 
# Genome
ln -s /ebio/abt6_projects8/Tarvense_TE/data/Tarvense_TE_pops/Thlaspi_genome_2021/Thlaspi_paper_final_genome/modified.fasta .
# Telb 
ln -s /ebio/abt6_projects8/Tarvense_TE/data/Tarvense_TE_pops/Thlaspi_genome_2021/Thlaspi_paper_final_genome/thlaspi.TE.fa .
#Prepare TE lib for blast:
makeblastdb -in thlaspi.TE.fa   -parse_seqids -dbtype nucl -out thlaspi.TE

# For either mapping long reads or computing whole-genome alignments, 
# Winnowmap requires pre-computing high frequency k-mers (e.g., top 0.02% most frequent) in a reference. 
# Winnowmap uses meryl k-mer counting tool for this purpose.

meryl count k=15 output merylDB modified.fasta
meryl print greater-than distinct=0.9998 merylDB > repetitive_k15.txt

winnowmap -W repetitive_k15.txt -t 16 -ax map-pb modified.fasta Armenian_Thlaspi.q20.5mC.fastq.gz |
samtools sort -@ 16 -o Armenian_Thlaspi.q20.5mC.bam - ;
samtools index Armenian_Thlaspi.q20.5mC.bam

# Sniffles2 variant calling
sniffles --reference modified.fasta \
         --output-rnames \
         -i Armenian_Thlaspi.q20.5mC.bam \
         -v Armenian_Thlaspi.q20.5mC.2.vcf

# Filter out 0/0  instances:
vawk --header '{ if ( S$SAMPLE$GT != "0/0") print }'  Armenian_Thlaspi.q20.5mC.vcf > Armenian_Thlaspi.q20.5mC.filt.vcf

# Extract IDs that overlap with TIPS within 2 Kbpp from each other:
TIPS="/ebio/abt6_projects8/Tarvense_TE/data/Tarvense_TE_pops/FINAL_ANALYSIS/TIPS/TA_AM_01_01_F3_CC0_M1_1_tips.gff3"

# Add an ID to the TIPs
cat $TIPS | awk '{print $0";ID_"NR}' > TIPS_IDs.gff3 

# Extract the list of SVs that are closer than 200 bp to a TIP loci.
vawk '{ print $0}' Armenian_Thlaspi.q20.5mC.filt.vcf  | 
grep -P "^Scaffold_[1-7]\t" | 
awk '{print $1"\t"$2"\t"$2+1"\t"$3}' |  
bedtools merge -i - -d -1 -c 4 -o collapse  | # Collapse SVs that span the same loci 
bedtools  closest -d  -a - -b TIPS_IDs.gff3  |
awk 'BEGIN{OFS="\t"}{if($14 <= 200) print $5,$6,$7,$8,$9,$10,$11,$12,$13";SVsupport="$4}' > TIPS_SV_support.gff3
