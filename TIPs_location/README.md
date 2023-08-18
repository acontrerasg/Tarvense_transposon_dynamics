# Count TIPs and TIPs density in different genomic features
This can be useful for example to observe the insertion behaviour of different TE families. Whether they are more likely to insert in already present TEs or in genes.

[TIPs_features_intersect](https://github.com/acontrerasg/Tarvense_transposon_dynamics/blob/main/TIPs_location/TIPs_features_intersect.sh) <br/>
Simple script to count TIPs and TIPs density of different groups (families, superfamilies or even TIPs present in different geographic regions) in different genomic features (genes, annotated TEs, promoters, intergenic regions ...). Density is calculated as number of TIPs for each MB of genome covered by each genomic feature.

[plotTIPs_location](https://github.com/acontrerasg/Tarvense_transposon_dynamics/blob/main/TIPs_location/plotTIPs_location.R) <br/>
R script to plot the results.
