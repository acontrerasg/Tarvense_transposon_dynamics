#!/usr/bin/Rscript

library(ComplexHeatmap)
library(dplyr)
library(vroom)
library(magick) # for better raster
library(tibble)


Region_colors_underscore <- c(
  "South_France" = "#008fc3",
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


# Matrixes load and metadata load. 
METADATA <-  vroom("/ebio/abt6_projects8/Tarvense_TE/doc/METADATA.tsv", col_names = TRUE)
# T A P S 
MatrixTAPS <-  vroom("/ebio/abt6_projects8/Tarvense_TE/data/Tarvense_TE_pops/FINAL_ANALYSIS/TAPS/MATRIX/final_TAP_matrix.tsv", col_names = TRUE)
MatrixTAPS <-  MatrixTAPS[ , !names(MatrixTAPS) %in% c("TA_SE_14_12_F1_HC0_M1_1")]
MatrixTAPS <- MatrixTAPS  %>%  mutate(name = paste(MatrixTAPS$Chrom,MatrixTAPS$Start,MatrixTAPS$End, sep="-"))

MatrixTAPS <- tibble::column_to_rownames(MatrixTAPS, var = "name")
MatrixTAPS  <- MatrixTAPS[,-c(1:9)]
row_sub = apply(MatrixTAPS, 1, function(row) any(row !=0))
# Thin the MATRIX_TAPS to only present the rows with at least one absence:
MatrixTAPS2 <- MatrixTAPS[row_sub,]
MatrixTIP <-  vroom("/ebio/abt6_projects8/Tarvense_TE/data/Tarvense_TE_pops/FINAL_ANALYSIS/TIPS/MATRIX/ALL_TIPS.matrix",  col_names = TRUE)
# T I P S 
MatrixTIP <-  MatrixTIP[ , !names(MatrixTIP) %in% c("TA_SE_14_12_F1_HC0_M1_1")]
MatrixTIP <- MatrixTIP  %>%  mutate(name = paste(MatrixTIP$`#Chrom`,MatrixTIP$Start,MatrixTIP$End, sep="-"))
MatrixTIP <- tibble::column_to_rownames(MatrixTIP, var = "name")
MatrixTIP <- MatrixTIP[,-c(1:3)]

MatrixTIP <- MatrixTIP[ , order(names(MatrixTIP))]
MatrixTAPS2 <- MatrixTAPS2[ , order(names(MatrixTAPS2))]

PAV_MATRIX <-  dplyr::bind_rows(MatrixTAPS2, MatrixTIP)

#Remove all but Scaffolds 1 to 7  and sort
PAV_MATRIX <-  PAV_MATRIX  %>%  filter(grepl("Scaffold_[1-7]-", rownames(PAV_MATRIX)))
PAV_MATRIX <-  PAV_MATRIX[stringi::stri_sort(rownames(PAV_MATRIX), numeric=TRUE),] 
#Replace Na to 0 
#PAV_MATRIX[is.na(PAV_MATRIX)] <- 0

mat <- as.matrix(PAV_MATRIX)

print("Final Matrix:")
print("Number of rows: ")
print(nrow(mat))
print("Number of cols: ")
print(ncol(mat))



# to set the metadata in the same order as the  mat columns


metadata_mat <- METADATA[match(colnames(mat), METADATA$Accession), ] 

column_ha <- HeatmapAnnotation(Region=metadata_mat$Region, col=list(Region=Region_colors_underscore))

heatmap_legend_param = list(
  title = "PAVs", at = c(-1, 0, 1), 
  labels = c("Absdence", "Fix", "Insertion")
)




PAV_subset <- dplyr::sample_n(PAV_MATRIX, 10000, replace = FALSE)
mat_subset <-  as.matrix(PAV_subset)

chr_subset <-  as.vector(stringr::str_remove(rownames(PAV_subset), pattern = "-.*$")) 
Heatmap(mat_subset, 
        col=c("#67a9cf", "#f7f7f7", "#E92424"), 
        show_row_names = FALSE, 
        show_column_names = FALSE, 
        cluster_columns = TRUE,
        row_split = chr_subset, cluster_rows = FALSE, show_column_dend = FALSE,
        use_raster = TRUE,
        top_annotation = column_ha,
        heatmap_legend_param = list(
          title = "PAVs", at = c(-1, 0, 1), 
          labels = c("Absent", "Fix", "Insertion")
        ))


print("Plotting the heatmap and saving it to an object...")
chr <-  as.vector(stringr::str_remove(rownames(PAV_MATRIX), pattern = "-.*$")) 
heatmap_plot <- Heatmap(mat, 
                        col=c("#67a9cf","#f7f7f7","#000000"), 
                        show_row_names = FALSE, 
                        row_split = chr, cluster_rows = FALSE,
                        show_column_names = FALSE, 
                        cluster_columns = TRUE,
                        use_raster = TRUE,
                        top_annotation = column_ha,
                        heatmap_legend_param = list(
                          title = "PAVs", at = c(-1, 0, 1), 
                          labels = c("Absence", "Fix", "Insertion")
                        ))

filetow <- "/ebio/abt6_projects8/Tarvense_TE/doc/heatmap_TEpolymorphisms.pdf"
print(paste0("start writting to file: ", filetow))
pdf(file=filetow, width = 8, height = 8)
draw(heatmap_plot)
dev.off()
##############



MatrixTIP <-  MatrixTIP  %>%  filter(grepl("Scaffold_[1-7]-", rownames(MatrixTIP)))
MatrixTIP <-  MatrixTIP[stringi::stri_sort(rownames(MatrixTIP), numeric=TRUE),]
chr_tip <-  as.vector(stringr::str_remove(rownames(MatrixTIP), pattern = "-.*$"))

MatrixTAPS2 <-  MatrixTAPS2  %>%  filter(grepl("Scaffold_[1-7]-", rownames(MatrixTAPS2)))
 MatrixTAPS2 <-  MatrixTAPS2[stringi::stri_sort(rownames(MatrixTAPS2), numeric=TRUE),]
chr_taps2 <-  as.vector(stringr::str_remove(rownames(MatrixTAPS2), pattern = "-.*$"))

mat_tip <- as.matrix(MatrixTIP)
mat_tap <- as.matrix(MatrixTAPS2)

heatmap_plot_tip <- Heatmap(mat_tip,
                            col=c("#f7f7f7", "#000000"), 
                            show_row_names = FALSE,
                            row_split = chr_tip, cluster_rows = FALSE,
                            show_column_names = FALSE,
                            cluster_columns = TRUE,
                            use_raster = TRUE,
                            top_annotation = column_ha,
                            heatmap_legend_param = list(
                              title = "PAVs", at = c(0, 1),
                              labels = c("Fix", "Insertion")
                            ))

heatmap_plot_tap <- Heatmap(mat_tap,
                            col=c("#67a9cf","#f7f7f7"), 
                            show_row_names = FALSE,
                            row_split = chr_taps2, cluster_rows = FALSE,
                            show_column_names = FALSE,
                            cluster_columns = TRUE,
                            use_raster = TRUE,
                            top_annotation = column_ha,
                            heatmap_legend_param = list(
                              title = "PAVs", at = c(-1, 0),
                              labels = c("Absence", "Fix")
                            ))

filetow <- "/ebio/abt6_projects8/Tarvense_TE/doc/heatmap_TEpolymorphisms_splitchrom_tip.pdf"
print(paste0("start writting to file: ", filetow))
pdf(file=filetow, width = 8, height = 8)
draw(heatmap_plot_tip)
dev.off()

filetow <- "/ebio/abt6_projects8/Tarvense_TE/doc/heatmap_TEpolymorphisms_splitchrom_tap.pdf"
print(paste0("start writting to file: ", filetow))
pdf(file=filetow, width = 8, height = 8)
draw(heatmap_plot_tap)
dev.off()

###### Subset Chrom 4 and 5 of TIPs and TAPS


MatrixTIP_45 <- MatrixTIP %>%  filter(grepl( "Scaffold_[4-5].*", rownames(MatrixTIP))) 
chr_tip <-  as.vector(stringr::str_remove(rownames(MatrixTIP_45), pattern = "-.*$"))
MatrixTAPS2_45 <- MatrixTAPS2 %>%  filter(grepl( "Scaffold_[4-5].*", rownames(MatrixTAPS2))) 
chr_taps2 <-  as.vector(stringr::str_remove(rownames(MatrixTAPS2_45), pattern = "-.*$"))

mat_tip_45 <- as.matrix(MatrixTIP_45)
mat_tap_45 <- as.matrix(MatrixTAPS2_45)


heatmap_plot_tip_45 <- Heatmap(mat_tip_45,
                            col=c("#f7f7f7", "#000000"), 
                            show_row_names = FALSE,
                            row_split = chr_tip, cluster_rows = FALSE,
                            show_column_names = FALSE,
                            cluster_columns = TRUE,
                            use_raster = TRUE,
                            top_annotation = column_ha,
                            heatmap_legend_param = list(
                              title = "", at = c(0, 1),
                              labels = c("Fix", "Insertion")
                            ))

heatmap_plot_tap_45 <- Heatmap(mat_tap_45,
                            col=c("#67a9cf","#f7f7f7"), 
                            show_row_names = FALSE,
                            row_split = chr_taps2, cluster_rows = FALSE,
                            show_column_names = FALSE,
                            cluster_columns = TRUE,
                            use_raster = TRUE,
                            top_annotation = column_ha,
                            heatmap_legend_param = list(
                              title = "", at = c(-1, 0),
                              labels = c("Absence", "Fix")
                            ))



filetow <- "/ebio/abt6_projects8/Tarvense_TE/doc/heatmap_TEpolymorphisms_splitchrom_tip_45.pdf"
print(paste0("start writting to file: ", filetow))
pdf(file=filetow, width = 8, height = 8)
draw(heatmap_plot_tip_45)
dev.off()

filetow <- "/ebio/abt6_projects8/Tarvense_TE/doc/heatmap_TEpolymorphisms_splitchrom_tap_45.pdf"
print(paste0("start writting to file: ", filetow))
pdf(file=filetow, width = 8, height = 8)
draw(heatmap_plot_tap_45)
dev.off()


