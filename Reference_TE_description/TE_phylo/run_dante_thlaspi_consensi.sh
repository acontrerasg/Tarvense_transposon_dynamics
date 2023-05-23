#!/usr/bin/env bash

conda activate /ebio/abt6_projects8/Tarvense_TE/code/conda_envs/LAST_tool
# mamba install -c conda-forge biopython
# mamba install -c anaconda numpy
# cd /ebio/abt6_projects8/Tarvense_TE/code/
# git clone https://github.com/kavonrtep/dante.git

mkdir -p /ebio/abt6_projects8/Tarvense_TE/data/Tarvense_TE_phylo/dante/
cd /ebio/abt6_projects8/Tarvense_TE/data/Tarvense_TE_phylo/dante/

#Python scripts
dante="/ebio/abt6_projects8/Tarvense_TE/code/dante/dante.py"
dante_filt="/ebio/abt6_projects8/Tarvense_TE/code/dante/dante_gff_output_filtering.py"

#Query
thlaspi_TE_consensi="/ebio/abt6_projects8/Tarvense_TE/data/Tarvense_TE_pops/Thlaspi_genome_2021/thlaspi.TE.fa"

#Dante DB
viriplantae="/ebio/abt6_projects8/Tarvense_TE/code/dante/tool-data/protein_domains/Viridiplantae_v3.0_pdb"
class_file="/ebio/abt6_projects8/Tarvense_TE/code/dante/tool-data/protein_domains/Viridiplantae_v3.0_class"


mkdir -p dante_tmp
python3 ${dante} -q  ${thlaspi_TE_consensi} -pdb ${viriplantae} --domain_gff thlaspi_consensi_class_dante.gff3 -dir dante_tmp  -cs ${class_file}


# Fitler the raw input with default parameters:

    # alignment length threshold 0.8 ;
    # proportion of alignment identity threshold  0.35 ;
    # threshold for alignment proportional similarity: 0.45.

# Except for interruptions (frameshifts + stop codons) tolerance. Threshold per 100 AA (default: 3) ;
# I set it to 10  because I am not working with full length copies but its conensus.

python3 ${dante_filt} --dom_gff ./dante_tmp/thlaspi_consensi_class_dante.gff3  \
                      -ouf ./thlaspi_consensi_class_dante_filt_ir10.gff3 \
                      -dir ./ \
                      -ir 10

# Retrieving TE family and its final classification.

for fam in $(cat thlaspi_consensi_class_dante_filt_ir10.gff3 | grep -v "^##" | cut -f1 | sort | uniq) ; do
  awk -v fam=$fam 'BEGIN{OFS="\t"}{if($1==fam) print $0}' thlaspi_consensi_class_dante_filt_ir10.gff3 |
  cut -f1,9 |
  tr ';' '\t' |
  # highest identity  and similarity up
  sort -rnk10,11 |
  head -n 1 |
  cut -f1,3 |
  sed 's/Final_Classification=//' |
  sed 's/|/\t/'
done   > tmp 

# Modify this tmp file to the final format: "PHYLO_CLASS PHYLO_SUBCLASS PYLO_ORDER PHYLO_LINEAGE:
cat  tmp |
sed 's/|/\t/' | 
sed 's/|/\t/' | 
awk '{if($3 == "pararetrovirus") print $0"\tpararetrovirus\tpararetrovirus" ; 
         else if($3 == "LINE") print $0"\tLINE\tLINE" ; 
         else print $0}' |
 awk '{if($5 == "") print $1,$2,$3,$4,"NA" ; else print $0}' | 
sed  '1 i\TEFAM\tPHYLO_CLASS\tPHYLO_SUBCLASS\tPYLO_ORDER\tPHYLO_LINEAGE'  >  thlaspi_TEfam_consensusClass.tsv

#cleanup
rm tmp 
