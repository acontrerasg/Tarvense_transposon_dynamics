library(dplyr)
library(ggplot2)
library(vroom)
library(patchwork)
library(ggpubr)
library(regionR)

'%!in%' <- function(x,y)!('%in%'(x,y))

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
  #"DHH"= "#a94a70",
  "CHI" = "black"
)

Pops_colors_mn106 <- c(
  "FR"= "#008fc3",
  "SE" ="#cbbb42" ,
  "DE" = "#2d2f2d",
  "NL" = "#730000",
  "AM" = "#2e6525",
  "CH" = "#c51b8a",
  "MN106A" = "#9B59B6"
)


#Create the dataset in bash:
## sed 's/;.*$//' sorted.final_anno_CovF.gff |  bedtools  closest -d -a ALL_insertions_pop_unique_pos.bed  -b - |  
# cut -f1-7,9,10,16,17 >   TIPS_distancefeat.tsv


TIPs <- vroom::vroom("~/Tarvense_TE/data/Tarvense_TE_pops/Splitreader/RESULTS/Results_analysis/ALL_insertions_pop_unique_pos.bed", 
                                  col_names = FALSE) 

TIP_distance_feat <- vroom::vroom("~/Tarvense_TE/data/Tarvense_TE_pops/Splitreader/RESULTS/Results_analysis/genomic_relation/TIPS_distancefeat.tsv", 
                                  col_names = FALSE) 


total_TIPS_perPOP <- TIPs %>%  group_by(X4,X7) %>%  tally() 
colnames(total_TIPS_perPOP ) <- c("X4", "X7", "total")

TIP_distance_feat %>%  filter(X9 == "gene") %>%  filter(X11 < 2000)  %>%  
  group_by(X4,X7) %>% tally() %>% 
  dplyr::left_join(., total_TIPS_perPOP, by=c("X4","X7")) %>%
  mutate(percent =  n/total*100) %>% 
  ggplot(aes(x=X4, y=percent, fill=X4)) + 
  geom_col(alpha=0.8) + 
  facet_wrap(~X7) + 
  ylab("%TIPS within 2kb of Genes") + 
  xlab("") + 
  scale_y_continuous(expand = c(0,0), limits = c(0,110)) +  
  scale_fill_manual(values=superfam_colors) +
  theme_light(base_size = 15) + 
  theme(axis.text.x = element_blank(), axis.ticks = element_blank(), legend.title = element_blank(), 
        strip.background = element_rect(fill="black")) 


total_TIPS_perind <- TIPs %>%  group_by(X4,X6,X7) %>%  tally() 
colnames(total_TIPS_perind) <- c("X4", "X6", "X7", "total")


TIP_distance_feat %>%  filter(X9 == "gene") %>%  filter(X11 < 2000)  %>%  
  group_by(X4,X6,X7) %>% tally() %>% 
  dplyr::left_join(., total_TIPS_perind, by=c("X4", "X6", "X7")) %>%
  mutate(percent =  n/total*100) %>% 
  ggplot(aes(x=X4, y=percent, color=X4)) + 
  geom_boxplot(fill="white") +
  geom_jitter(alpha =0.6) + 
  facet_wrap(~X7) + 
  ylab("%TIPS  per individual within 2kb from gene") + 
  xlab("") + 
  scale_y_continuous(expand = c(0,0), limits = c(0,110)) +
  scale_color_manual(values=superfam_colors) +
  theme_light(base_size = 15) + 
  theme(axis.text.x = element_blank(), axis.ticks = element_blank(), legend.title = element_blank(), 
          strip.background = element_rect(fill="black")) 


TIP_distance_feat %>%  filter(X8 == "RepeatMasker") %>%  filter(X11 < 2000)  %>%  
  group_by(X4,X6,X7) %>% tally() %>% 
  dplyr::left_join(., total_TIPS_perind, by=c("X4", "X6", "X7")) %>%
  mutate(percent =  n/total*100) %>% 
  ggplot(aes(x=X4, y=percent, color=X4)) + 
  geom_boxplot(fill="white") +
  geom_jitter(alpha =0.6) + 
  facet_wrap(~X7) + 
  ylab("%TIPS  per individual within 2kb from TE") + 
  xlab("") + 
  scale_y_continuous(expand = c(0,0), limits = c(0,110)) +
  scale_color_manual(values=superfam_colors) +
  theme_light(base_size = 15) + 
  theme(axis.text.x = element_blank(), axis.ticks = element_blank(), legend.title = element_blank(), 
        strip.background = element_rect(fill="black")) 


TIP_distance_feat %>%  filter(X8 != ".") %>%  
  filter(X9 %!in% c("intron", "exon", "five_UTR", "three_UTR"))  %>%  
  group_by(X8)  %>%  
  ggplot(aes(x=X8, y=X11/1000)) + 
  geom_boxplot() + 
  xlab('') +
  ylab("Distance to gene (kb)") + 
  facet_grid(~X7) + 
  stat_compare_means(method = "t.test" , label.y = 2000) +
  theme_light(base_size = 15) + 
  theme(axis.text.x = element_text(angle=25, size=10, hjust =1), legend.title = element_blank(), 
        strip.background = element_rect(fill="black")) 


  
