## Author: Dario Galanti Sep 2022
## Aim: Plot TIP intersection with genomic features from TIPs_features_intersect.sh


### SETUP
library("dplyr")
library("data.table")
library(stringr)
library(ggplot2)
library(scales)

Features_colors <- c(
  "global" = "grey40",
  #"genes" = ,
  "CDS" = "#D06465",
  "introns" = "#BF6B92",
  "TEs" = "#698838",
  "prom" = "#CB6E34",
  "interg" = "#7BA5D6")


### DEFINE INPUT FILE
fin <- "Features_intersect_SupFam_TIPs.txt"
fin <- "Features_intersect_Fam_TIPs.txt"

### READING FIN
df <- fread(fin, data.table=FALSE)
df$Feature <- gsub("intergenic","interg",df$Feature)
df$Feature <- gsub("promoters","prom",df$Feature)

## SELECT FEATURES TO KEEP
#df$Feature <- factor(df$Feature, levels = c('global','genes','CDS','introns','TEs', 'prom', 'interg'), ordered = TRUE) #all
#df$Feature <- factor(df$Feature, levels = c('genes','TEs', 'prom', 'interg'), ordered = TRUE) # v1
df$Feature <- factor(df$Feature, levels = c('global','CDS','introns','TEs', 'prom', 'interg'), ordered = TRUE) #v2 (width=7)
#df$Feature <- factor(df$Feature, levels = c('CDS','TEs', 'prom', 'interg'), ordered = TRUE) #v3
df <- df[complete.cases(df), ]
rownames(df) <- NULL


### Cleanup
df <- df[!(df$Group == "RIC"),]
df <- df[!(df$Group == "RIC" | df$Group == "DTT"),]
# Data are on very different scales

## Cleanup families with less than x global insertions
cleanup <- df[(df$Feature == "global" & df$Num_TIPs < 100),]
remove_fam <- cleanup$Group
for (fam in remove_fam){df <- df[!df$Group == fam,]}

### PLOT INSERTIONS PER MB INTO DIFFERENT FEATURES FOR EACH SUPERFAMILY
pdf("Insertions_x_Mb_byFam_100TIPs.pdf", width = 12, height = 14)
#pdf("Insertions_x_Mb_bySupFam.pdf", width = 10, height = 8)
ggplot(df, aes(fill=Feature, y=insertions_x_Mb, x=Feature)) + 
  geom_bar(stat="identity") +            # Normal stacked barplot
  guides(fill=guide_legend(title="Feature")) +
  facet_wrap(~ Group, ncol = 5) +
  scale_fill_manual(values=Features_colors, aesthetics = "fill") +
  xlab("Feature") + ylab("N° insertions / Mb") + 
  theme(legend.position = "right", legend.text= element_text(size=12), legend.title = element_blank(),
        #axis.text.x = element_text(size=12, face = "bold", hjust = 0.5, vjust = 0.5),
        axis.text.x = element_text(size=12, face="bold", angle=45, hjust=0.8, vjust=0.8), # Turned lables
        axis.text.y = element_text(size=12),
        axis.title.x = element_blank(), axis.title.y = element_text(size=14),
        panel.grid.major = element_line("white"), panel.grid.minor = element_line("white"),
        strip.background = element_rect(fill="white", color="white"),
        panel.background = element_rect(fill="white"), strip.text = element_text(color="grey20", size = 12))
dev.off()


### PLOT RATE OF INSERTION PER MB INTO DIFFERENT FEATURES FOR EACH SUPERFAMILY
pdf("Insertions_x_Mb_byFam_scaled1.pdf", width = 10, height = 20)
#pdf("Insertions_x_Mb_bySupFam_scaled1.pdf", width = 10, height = 8)
ggplot(df, aes(fill=Feature, y=Insertion_x_MB_global_scaled, x=Feature)) + 
  geom_bar(stat="identity") +            # Normal stacked barplot
  guides(fill=guide_legend(title="Feature")) +
  facet_wrap(~ Group, ncol = 4) +
  scale_fill_manual(values=Features_colors, aesthetics = "fill") +
  xlab("Feature") + ylab("N° insertions / Mb") + 
  theme(legend.position = "right", legend.text= element_text(size=12), legend.title = element_blank(),
        #axis.text.x = element_text(size=12, face = "bold", hjust = 0.5, vjust = 0.5),
        axis.text.x = element_text(size=12, face="bold", angle=45, hjust=0.8, vjust=0.8), # Turned lables
        axis.text.y = element_text(size=12),
        axis.title.x = element_blank(), axis.title.y = element_text(size=14),
        panel.grid.major = element_line("white"), panel.grid.minor = element_line("white"),
        strip.background = element_rect(fill="white", color="white"),
        panel.background = element_rect(fill="white"), strip.text = element_text(color="grey20", size = 12))
dev.off()



