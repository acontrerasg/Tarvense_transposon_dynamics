## Author: Dario Galanti
## Aim: Plot methylation in the flanking regions of insertions for samples with vs without the insertion
## Run: Rscript --vanilla plotTIPflank_met.R Flks_50binned_met.txt fout.pdf

### SETUP
library(dplyr)
library(plyr)       #for round_any
library(data.table)
library(stringr)
library(ggplot2)

### DEFINE INPUT FILE
args <- commandArgs(trailingOnly = TRUE)

fin <- args[1]
#fin <- "RNA_TIPs/TIP_R_flks_met_50bin_2kbself_1kbTEs.csv"
fout <- args[2]
#fout <- "TIP_R_flks_met_50bin_2kbself_smooth.pdf"

### READING FIN
df <- fread(fin, data.table=FALSE)

### ADD MOVING MEANS
bins <- nrow(df)/2
## With insertion
flk1 <- c(mean(df$Samples_with_insertion[1:2]),df$Samples_with_insertion[1:bins],mean(df$Samples_with_insertion[(bins-1):bins]))
flk2 <- c(mean(df$Samples_with_insertion[(bins+1):(bins+2)]),df$Samples_with_insertion[(bins+1):(bins*2)],mean(df$Samples_with_insertion[(bins*2-1):(bins*2)]))
df$Insertion <- c(apply(embed(flk1, 3), 1, mean),apply(embed(flk2, 3), 1, mean))
# without insertion
flk1 <- c(mean(df$Samples_without_insertion[1:2]),df$Samples_without_insertion[1:bins],mean(df$Samples_without_insertion[(bins-1):bins]))
flk2 <- c(mean(df$Samples_without_insertion[(bins+1):(bins+2)]),df$Samples_without_insertion[(bins+1):(bins*2)],mean(df$Samples_without_insertion[(bins*2-1):(bins*2)]))
df$No.Insertion <- c(apply(embed(flk1, 3), 1, mean),apply(embed(flk2, 3), 1, mean))


### Stack df
#stack_df <- data.frame(df[1:4],stack(df[5:6]))
stack_df <- data.frame(df[1:4],stack(df[7:8]))  # If smoothed values are also present
colnames(stack_df)[5:6] <- c("Methylation", "Group")

top <- max(stack_df$Methylation)
top <- round_any(top, 20, f = ceiling)
#top <- max(top,40)

### Cleanup names
stack_df$Group <- gsub("No.Insertion","No Insertion", stack_df$Group)
#stack_df$Group <- gsub("Samples_with_insertion_smooth","Insertion",stack_df$Group)
#stack_df$Group <- gsub("Samples_without_insertion_smooth","No Insertion",stack_df$Group)
#stack_df$Group <- gsub("Samples_with_insertion","Insertion",stack_df$Group)
#stack_df$Group <- gsub("Samples_without_insertion","No Insertion",stack_df$Group)

### PLOT 1
#pdf("TIPflanks_met.pdf", width = 12, height = 8)
#ggplot(stack_df, aes(x = bp, y = Methylation, colour = Group) ) +
#  geom_line() +
#  theme_bw() +
#  facet_wrap(~Flk, scales = "free_x")
#dev.off()

### PLOT 2
pdf(fout, width = 9, height = 6)
#pdf("RNA_TIPs/TIP_R_flks_met_50bin_2kbself_1kbTEs.pdf", width = 10, height = 6)
ggplot(stack_df, aes(x = bp, y = Methylation, colour = Group) ) +
  geom_line(size=1) +
  theme_bw() +
  facet_wrap(~Flk, scales = "free_x") +
  ylim(0,top) +
  ylab("Methylation %") + 
  theme(legend.position = "top", legend.text = element_text(size=19),
        axis.text = element_text(size=19),
        axis.title = element_text(size=21),
        panel.grid.minor = element_line("white"),
        strip.background = element_blank(), strip.text.x = element_blank(),   # Remove facet_wrap lables
        #panel.background = element_blank(), legend.title=element_blank())
        legend.title=element_blank(),
        plot.margin = margin(8,20,8,8, "points"))  #default 5.5 everywhere. Order: top,right,bottom,left
dev.off()


