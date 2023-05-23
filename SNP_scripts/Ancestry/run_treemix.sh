#!/bin/usr/env bash

#TREEMIX run (after ancestry pops have run).
cd /ebio/abt6_projects8/Tarvense_TE/data/Tarvense_ancestry/
mkdir -p treemix/

esalsu_vcf="/ebio/abt6_projects8/Tarvense_TE/data/Tarvense_ancestry/thlaspi_esalsu_near.allChrsHaplotypeData.vcf.gz"
sparvula_vcf="/ebio/abt6_projects8/Tarvense_TE/data/Tarvense_ancestry/thlaspi_sparvula_near.allChrsHaplotypeData.vcf.gz"
METADATA="/ebio/abt6_projects8/Tarvense_TE/doc/METADATA.tsv" 
#remeber to taske care of the sample names in the pop vcf files. 
pop_vcf="/ebio/abt6_projects8/Tarvense_TE/data/Tarvense_ancestry/Ta_v5_vrts_10NA_biall_1MAF_imp_withUSspls.vcf.gz"
plink2treemix="/ebio/abt6_projects8/Tarvense_TE/code/genetic_analysis/plink2treemix_py2.py"
plot_treemixR="/ebio/abt6_projects8/Tarvense_TE/code/genetic_analysis/treemix_plotting_output.R"
make_dataclust="/ebio/abt6_projects8/Tarvense_TE/code/genetic_analysis/make_data.clust_treemix.sh"
#Merge pop_vcf with the  2 vcf files from last ancestry
bcftools merge --merge all  ${pop_vcf}  ${esalsu_vcf}  ${sparvula_vcf} | bgzip > treemix/pop_woutgroups.vcf.gz

cd treemix/
# Create an index for the newly created vcf 
tabix -p vcf pop_woutgroups.vcf.gz
#Annotate the  new VCF and index it 
bcftools annotate -x ID -I +'%CHROM:%POS:%REF:%ALT' pop_woutgroups.vcf.gz  |
bgzip > pop_woutgroups_annotated.vcf.gz

tabix -p vcf pop_woutgroups_annotated.vcf.gz 

#convert the annotated vcf file to plink bed format
plink --vcf  pop_woutgroups_annotated.vcf.gz --recode   \
  --double-id \
  --allow-no-sex \
  --allow-extra-chr \
  --make-bed \
  --out pop_woutgroups

#reformat the ped file
awk '{print $1, $1"_"$2 }'  pop_woutgroups.ped  > pop_ind.tmp
cut -d " " -f 3-   pop_woutgroups.ped  > geno.tmp
cp pop_woutgroups.ped  pop_woutgroups.ped.back
paste pop_ind.tmp geno.tmp  > pop_woutgroups.ped
rm *.tmp


# Add this two new samples to the data.clust file:

# Create data.clust file using the pop_woutgroups.fam.
# The data.clust format consist of 3 columns the first two are the sample names as present in the vcf file and the 3rd column is the group you want
# them to cluster (Like Armenia,Netherlands etc...) as:
# "Sample_name\tSample_name\tGroup"
bash ${make_dataclust} $METADATA pop_woutgroups.vcf.gz
# Add specfic outgrops not present in metadata 
echo -e "thlaspi_esalsu \tthlaspi_esalsu\toutgroup\nthlaspi_sparvula\tthlaspi_sparvula\toutgroup" | 
cat data.clust - > tmp
mv tmp data.clust 

# plink to freq tables
plink --allow-extra-chr  --bfile pop_woutgroups  --freq --missing --within data.clust
gzip plink.frq.strat

python2 ${plink2treemix} plink.frq.strat.gz  treemix.strat.frq.gz


## RUN TREEMIX  no migrations

treemix  -i treemix.strat.frq.gz -root outgroup -bootstrap -k 100 -o  treemix_results_nomig_k100 | tee ./treemix.log


#create unique pops  file for the polot residuals R function:
awk '{if ( $3 == "Middle_Sweden")    color="#cbbb42" ;
else if ( $3  == "South_France")     color="#008fc3" ;
else if ( $3  == "South_Sweden")     color="#0a75ad" ;
else if ( $3  == "North_Germany")    color="#2d2f2d" ;
else if ( $3  == "South_Germany")    color="#730000" ;
else if ( $3  == "The_Netherlands")  color="#ea561b" ;
else if ( $3  == "Armenia")          color="#2e6525" ;
else if ( $3  == "Hefei")            color="#35a79c" ;
else if ( $3  == "Mangkang")         color="#009688" ;
else if ( $3  == "Xi’an" )           color="#65c3ba" ;
else if ( $3  == "Zuogong")          color="#54b2a9" ;
else if ( $3  == "outgroup")         color="black" ;
else if ( $3  == "America")         color="lightpink1" ;
print $3"\t"color}' data.clust  | uniq > pop_order.txt


#Rerun treemix with 25 boostraps:
mkdir -p /ebio/abt6_projects8/Tarvense_TE/data/Tarvense_ancestry/treemix/boostrap_runs/
cd /ebio/abt6_projects8/Tarvense_TE/data/Tarvense_ancestry/treemix/boostrap_runs/

ln -s /ebio/abt6_projects8/Tarvense_TE/data/Tarvense_ancestry/treemix/treemix.strat.frq.gz .
ln -s /ebio/abt6_projects8/Tarvense_TE/data/Tarvense_ancestry/treemix/pop_order.txt .

# Run 2.5k boostraps with treemix:

qsub ./treemix_cluster.sh


# Merge trees
### Create a file with all the bootstrapped trees
for num in $(seq 1 2500) ; do 
  gunzip -c .//boostrap_runs/treemix_boostrap_${num}.treeout.gz  |
  head -n 1 ;
done  >>  boostraped_tree_merged.tre


# Activate phylip envirment to merge trees:
conda activate /ebio/abt6_projects8/Tarvense_TE/code/conda_envs/phylip

# Find the position of Outgroup population
posOutgroup=$(head -1  boostraped_tree_merged.tre | tr "," "\n" | grep outgroup -n | cut -d":" -f1)
# echo $posOutgroup
echo  boostraped_tree_merged.tre > boostraped_tree_merged_PhylipInputFile
echo "O" >> boostraped_tree_merged_PhylipInputFile
echo $posOutgroup >>boostraped_tree_merged_PhylipInputFile
echo "Y" >> boostraped_tree_merged_PhylipInputFile



# Run Phylip consense
consense < boostraped_tree_merged_PhylipInputFile > consense.log
# Output is written as outtree. Modify it to netwick
cat outtree | tr -d "\n"  > outtree.netwick
# Rerun treemix  with the new  tree:
treemix -i treemix.strat.frq.gz  -k 100 -root outgroup -se -tf outtree.netwick -o final_tree_merged_boostraps_treemix > logfile_final_merged_tree.log




