#!/bin/bash

# Install Last withg ananconda:
## mamba install -c bioconda/label/cf201901 last  at fresh env /ebio/abt6_projects8/Tarvense_TE/code/conda_envs/LAST_tool

#lastdb -P0 -uNEAR -R01 ../thlaspi_final-3-NEAR ../thlaspi.fa
#lastdb -P0 -uMAM4 -R01 ../thlaspi_final-3-MAM8 ../thlaspi.fa

#Path and env

conda activate /ebio/abt6_projects8/Tarvense_TE/code/conda_envs/LAST_tool
cd /ebio/abt6_projects8/Tarvense_TE/data/Tarvense_ancestry
ln -s /ebio/abt6_projects8/Tarvense_TE/data/Tarvense_TE_pops/Thlaspi_genome_2021/Thlaspi_paper_final_genome/modified.fasta  ./thlaspi_final.fasta
#Create a fasta index
samtools faidx ./thlaspi_final.fasta
#Vars
sparvula_genome="/ebio/abt6_projects8/Tarvense_TE/data/Tarvense_ancestry/closest_genomes/Schrenkiella_parvula/Sparvula_574_v2.0.fa"
esalsu_genome="/ebio/abt6_projects8/Tarvense_TE/data/Tarvense_ancestry/closest_genomes/eutrema_salsugineum/GCF_000478725.1_Eutsalg1_0_genomic.fna"
pop_vcf="/ebio/abt6_projects8/Tarvense_TE/data/Tarvense_TE_pops/FINAL_ANALYSIS/FINAL_PAVs_normal/genetic_analysis/Ta_v5_vrts_10NA_biall_1MAF_imp_withUSspls.vcf.gz"
# Make last databases:
mkdir -p lastdb/
lastdb -P0 -uNEAR -R01 lastdb/thlaspi_final-3-NEAR ./thlaspi_final.fasta
lastdb -P0 -uMAM4 -R01 lastdb/thlaspi_final-3-MAM4 ./thlaspi_final.fasta
lastdb -P0 -uMAM8 -R01  lastdb/thlaspi_final-3-MAM8 ./thlaspi_final.fasta


#Explanation:
# The last-train command finds the rates of deletion, insertion, and each kind of substitution between these sequences, and writes them to
# a file. Then lastal's -p option uses this file to get more-accurate alignments.

# NEAR : This DNA seeding scheme is good for finding short-and-strong (near-identical) similarities.
# It is also good for similarities with many gaps (insertions and deletions), because it can find the short matches between the gaps.
# (Long-and-weak seeding schemes allow for mismatches but not gaps.) (FASTEST)

nice last-train -P0 --revsym --matsym --gapsym -E0.05 -C2 lastdb/thlaspi_final-3-NEAR ${sparvula_genome} > thlaspi_sparvula_near.mat &&
lastal -m50 -E0.05 -C2 -p thlaspi_sparvula_near.mat lastdb/thlaspi_final-3-NEAR ${sparvula_genome} |
last-split -m1 > thlaspi_sparvula_near.maf

nice last-train -P0 --revsym --matsym --gapsym -E0.05 -C2 lastdb/thlaspi_final-3-NEAR ${esalsu_genome} > thlaspi_esalsu_near.mat &&
lastal -m50 -E0.05 -C2 -p thlaspi_esalsu_near.mat lastdb/thlaspi_final-3-NEAR ${esalsu_genome} |
last-split -m1 > thlaspi_esalsu_near.maf


for i in *_near.maf ; do
  maf-swap ${i} |
  last-split -m1 > ${i%.maf}_2.maf ;
done

for file in *_near_2.maf ; do
  maf-swap ${file} > ${file%_2.maf}_3.maf
done

for file in *_near_3.maf ; do
  maf-convert sam -r "ID:1  PL:ILLUMINA SM:${file%_near_3.maf}" ${file} |
  samtools view -bt thlaspi_final.fasta.fai  |
  samtools sort >  ${file%_3.maf}.sorted.bam ;
done

conda activate base

for file in *_near.sorted.bam ; do
  samtools coverage ${file} > ${file}.coverage ;
  samtools mpileup -ugf ./thlaspi_final.fasta -l ${pop_vcf} ${file} |
  bcftools call -m |
  bcftools view -Oz -o ${file%.sorted.bam}.allChrsHaplotypeData.vcf.gz ;
  tabix -p vcf ${file%.sorted.bam}.allChrsHaplotypeData.vcf.gz ;
done

# NEXT: treemix for phylogeniy
