#!/usr/bin/env Rscript 

# load tidyverse package
library(tidyverse)
library(vroom)
library(ggplot2)
#Then we will use a combination of readr and the standard scan function to read in the data.

#Load metatada with location and Region names and colors:
METADATA <-  vroom("/ebio/abt6_projects8/Tarvense_TE/doc/METADATA.tsv", col_names = TRUE)
sample_names_regions <- METADATA[,c(1,3,9)]

Region_colors_underscore <- c(
  "South_France" = "#008fc3",
  #"North-East_France" = "#008fc3", #tbdeleted
  "Middle_Sweden" = "#cbbb42",
  "South_Sweden" = "#0a75ad", 
  "North_Germany" = "#2d2f2d",
  "South_Germany" = "#730000",
  "The_Netherlands" = "#ea561b",
  "Armenia" = "#2e6525",
  "Hefei" = "#35a79c" ,
  "Mangkang" = "#009688" ,
  "Xi’an" = "#65c3ba", 
  "Zuogong" = "#54b2a9",
  "America" = "lightpink1"
)




#### ALL #####
# read in data
path <-  "/ebio/abt6_projects8/Tarvense_TE/data/Tarvense_TE_pops/FINAL_ANALYSIS/genetic_analysis/SNPS_PCA/10NA_biall_1MAF_imp_PCA/"
pca <- read_table(paste0(path,"Tarvense_SNPs_PCA.eigenvec"), 
                  col_names = FALSE)

eigenval <- scan(paste0(path, "Tarvense_SNPs_PCA.eigenval"))

# sort out the pca data
# remove nuisance column
pca <- pca[,-1]
# set names
names(pca)[1] <- "Accession"
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))
# join Accessions and locations:
pca <- dplyr::inner_join(sample_names_regions, pca ,by = "Accession")  

############## Plot the data ##############

#First we will plot the eigenvalues. It is quite straightforward to translate these into percentage variance explained 
#(although note, you could just plot these raw if you wished).

# first convert to percentage variance explained
pve <- data.frame(PC = 1:20, pve = eigenval/sum(eigenval)*100)
#then plot with GG plot
a <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity")
a + ylab("Percentage variance explained") + theme_light()

# calculate the cumulative sum of the percentage variance explained
cumsum(pve$pve)

# plot pca ALL
psize = 2
b <- ggplot(pca, aes(PC1, PC2, col = Region)) + 
  geom_point(size =  psize) + 
  scale_colour_manual(values = Region_colors_underscore) + 
  coord_equal() + 
  theme_bw(base_rect_size = 1, base_size = 12) +
  geom_vline(xintercept = 0, linetype = "dotted") + 
  geom_hline(yintercept = 0,  linetype = "dotted") +
  xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + 
  ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))

c <- ggplot(pca, aes(PC3, PC4, col = Region)) + geom_point(size = psize) + 
  scale_colour_manual(values = Region_colors_underscore) + 
  coord_equal() +
  theme_bw(base_rect_size = 1, base_size = 12) +
  geom_vline(xintercept = 0, linetype = "dotted") + 
  geom_hline(yintercept = 0,  linetype = "dotted") +
  xlab(paste0("PC3 (", signif(pve$pve[3], 3), "%)")) + ylab(paste0("PC4 (", signif(pve$pve[4], 3), "%)"))

d <- ggplot(pca, aes(PC5, PC6, col = Region)) + geom_point(size =  psize) + 
  scale_colour_manual(values = Region_colors_underscore) + 
  coord_equal() +
  theme_bw(base_rect_size = 1, base_size = 12) +
  geom_vline(xintercept = 0, linetype = "dotted") + 
  geom_hline(yintercept = 0,  linetype = "dotted") +
  xlab(paste0("PC5 (", signif(pve$pve[5], 3), "%)")) + ylab(paste0("PC6 (", signif(pve$pve[6], 3), "%)"))

e <- ggplot(pca, aes(PC7, PC8, col = Region)) + geom_point(size =  psize) + 
  scale_colour_manual(values = Region_colors_underscore) + 
  coord_equal() +
  theme_bw(base_rect_size = 1, base_size = 12) +
  geom_vline(xintercept = 0, linetype = "dotted") + 
  geom_hline(yintercept = 0,  linetype = "dotted") +
  xlab(paste0("PC7 (", signif(pve$pve[7], 3), "%)")) + ylab(paste0("PC8 (", signif(pve$pve[8], 3), "%)"))

plot_pcaALL <-  b + c + d + e 
