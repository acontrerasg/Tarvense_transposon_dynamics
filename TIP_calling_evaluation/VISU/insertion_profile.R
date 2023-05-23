#!/usr/bin/env Rscript

library(vroom)
library(ggplot2)
library(dplyr)

superfam_colors <- c(
  "RLG" = "#9754d1",
  "RLC" = "#6fb445",
  "RLA" = "#cf46ac",
  "RIX" ="#4ea98a",
  "RIR" ="#da415f",
  "RIL" ="#5b98cd",
  "RII" ="#d65b2b",
  "RIC" ="#6474d8",
  "DTT" ="#bf9740",
  "DTM"= "#725b9a",
  "DTH"= "#657c38",
  "DTC"= "#d28bc6",
  "DTA"= "#b86d54",
  "DHH"= "#a94a70",
  "CHI" = "black"
)

#Load insertin files:
system_files  <- list.files(path = "~/Tarvense_TE/data/Tarvense_TE_pops/Splitreader/RESULTS/Results_analysis/insertion_analysis/", 
                               pattern = "*bed")

insertion_files <- list()
for (i in seq_along(system_files)) {
 insertion_files[[i]] <- assign(system_files[i], 
         vroom(paste0( "~/Tarvense_TE/data/Tarvense_TE_pops/Splitreader/RESULTS/Results_analysis/insertion_analysis/",
                       system_files[i]), col_names = c("chrom", "start", "end", "superfam", "fam", "indvidiual","pop")))
}


for(df in  insertion_files){
  name <-  unique(df$pop)
  plot <- df %>% filter(grepl("Scaffold_[1-7]$", chrom))%>%
  group_by(chrom, start, end, superfam, pop) %>%  tally() %>% 
  ggplot() + 
  geom_segment(aes(x = start, xend = end, y=0,yend=n, color=superfam), size = 1) + 
  scale_color_manual(values=superfam_colors)+
    facet_grid(superfam~chrom) + 
  ggtitle(name) + 
  theme_light() +
  ylab("#TIP") +
  theme(text = element_text(size=8), strip.background = element_rect(fill="black")) +
    theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.x = element_blank())
  ggsave(plot, filename = paste0("~/Tarvense_TE/doc/MOBILOME/insertion_profile_", paste(name, collapse = "_"),".pdf"), dpi=600 )
}
 
