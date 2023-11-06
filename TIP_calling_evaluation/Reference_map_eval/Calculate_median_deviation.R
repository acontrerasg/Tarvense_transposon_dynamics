#!/usr/bin/env Rscript

# This script takes the  mean of both read depth corrrected files calculated with mosdepth and bwa mem and averages both of them,
# then it calculates median coverage. LAstly, assigns regions of the genome with 3 times lower than the median 
# coverage as LOW; 3 times higher than median as HIGH and the rest as PASS,

GCcorected_depth_means <- read.delim("./GCcorected_depth_means.bed", header=FALSE)

upper_limit <- median(GCcorected_depth_means$V4)*3
lower_limit <- median(GCcorected_depth_means$V4)/3

filtered_regions_bed <-  GCcorected_depth_means  %>% 
  dplyr::filter( !(V4 > upper_limit ))  %>% 
  dplyr::filter( !(V4 < lower_limit ))
write.table(filtered_bed[,1:3],  
            file="./Regions_filtered_3Medians.bed", 
            quote= FALSE, 
            sep='\t', 
            col.names = FALSE, 
            row.names = FALSE )






