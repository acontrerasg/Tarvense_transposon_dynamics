#!/usr/bin/env Rscript

library(rtracklayer)
library(vroom)
library(karyoploteR)

chromlengths <- read.table('~/Tarvense_TE/data/Tarvense_TE_pops/Thlaspi_genome_2021/Thlaspi_paper_final_genome/modified.fasta.fai')
chromlengths <- chromlengths[,c(1,2)]
tarvense_genome <- getGenome(chromlengths)
genes_granges <- import("~/Tarvense_TE/data/Tarvense_TE_pops/Thlaspi_genome_2021/Thlaspi_paper_final_genome/sorted.final.gene_lifted.gff", format="GFF")
TE_granges <-  import("~/Tarvense_TE/data/Tarvense_TE_pops/Thlaspi_genome_2021/Thlaspi_paper_final_genome/sorted.final.TEs_lifted.gff", format="GFF")
TIP_allele_freqs_granges <-  import("~/Tarvense_TE/data/Tarvense_TE_pops/FINAL_ANALYSIS/TIPS/allele_freqs_bg/ALL_TIPS_allele.freq.bedgraph", 
                                    format = "bedGraph")
# TIP allele means
TIP_allele_mean_freqs_100kbwin <- vroom("~/Tarvense_TE/data/Tarvense_TE_pops/FINAL_ANALYSIS/TIPS/allele_freqs_bg/ALL_TIPS_allele.freq.windowed.mean_freqs.bedgraph", 
                                        skip = 1, col_names = FALSE)
kp <- plotKaryotype(genome = tarvense_genome, 
                    chromosomes=c("Scaffold_1", "Scaffold_2", "Scaffold_3","Scaffold_4", "Scaffold_5", "Scaffold_6" , "Scaffold_7"), 
                    plot.type=4)
kp <- kpPlotDensity(kp,TE_granges, col="#636363", window.size = 1000000, r0=0, r1=0.20)
