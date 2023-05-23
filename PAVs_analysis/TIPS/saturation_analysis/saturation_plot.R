
saturation_curve <- read.delim("/ebio/abt6_projects8/Tarvense_TE/data/Tarvense_TE_pops/FINAL_ANALYSIS/TIPS/saturation_analysis/saturation_curve.tsv", header=FALSE)

Region_colors <- c(
  "South_France" = "#008fc3",
  #"North-East France" = "#008fc3", #tbdeleted
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
  "America" = "#ffb6c1"
)

plot <- saturation_curve  %>% 
  ggplot(aes(x=V1, y=log2(V3), color=V2)) + 
  geom_line() +
  scale_color_manual(values = Region_colors,name="Region") + 
  theme_light() + 
  facet_wrap(~V2) + 
  xlab("Individuals") +
  ylab("log2 Unique insertions") +
  theme(strip.background = element_rect(fill="white" , color="black"), 
  strip.text =  element_text(color="black"))
  
ggsave(plot, filename = "/ebio/abt6_projects8/Tarvense_TE/doc/MOBILOME/figures/Saturation_curves.pdf")
