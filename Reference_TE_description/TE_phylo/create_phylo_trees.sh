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


#GET RT  from LTRs

for  LTR in Ty3 Ty1 LINE ; do 
  
  python3 ${dante_filt} --dom_gff ./dante_tmp/thlaspi_consensi_class_dante.gff3  \
                        -dps ./thlaspi_consensi_RT_${LTR}_ir10.fa \
                        -dir ./ \
                        -sd RT \
                        -el ${LTR} \
                        -ir 10 ;
done
                      

# sd options {All,GAG,INT,PROT,RH,RT,aRH,CHDCR,CHDII,TPase,YR,HEL1,HEL2,ENDO}
#el ( arbitrary substring of the element classification ('Final_Classification' attribute in GFF)
 #  I  will use Helitron  ; CACTA ; Mutator ; Harbinger ; Mariner ;  LINE
 

#  Get TPase from TIR elements 

for TIR in CACTA Mutator Harbinger Mariner ; do 
  python3 ${dante_filt} --dom_gff ./dante_tmp/thlaspi_consensi_class_dante.gff3  \
                        -dps ./thlaspi_consensi_TPase_${TIR}_ir10.fa \
                        -dir ./ \
                        -sd TPase \
                        -el ${TIR} \
                        -ir 10  ;
 done 

# GET HEL1 (more abundant (82) than HEL2 (54) ; from Helitrons

  python3 ${dante_filt} --dom_gff ./dante_tmp/thlaspi_consensi_class_dante.gff3  \
                        -dps ./thlaspi_consensi_HEl1_Helitron_ir10.fa \
                        -dir ./ \
                        -sd HEL1 \
                        -el Helitron \
                        -ir 10  ;

                      
 # Make  one MSA  for each superfamily using mafft 
 
parallel "mafft  --maxiterate 1000 --localpair --thread 4 {} > {.}.aln" ::: thlaspi_consensi_RT_*_ir10.fa


#Generate raxml.reduced.phy  trees by running raxml-ng with aln  sequences (throws and error and fixed input (spaces in names):
for msa in thlaspi_consensi*aln ; do 
  echo  " working with .... ${msa} " ; 
  raxml-ng -msa $msa  -model JTT+G ;
done 
# Once RAXML has  fixed the seq names for us (how nice) after  panic error messages (not so nice)....
 
 #Run RAXML with the JTT+Gamma model for proteins and 100 boostrap trees:
for reduced_phy in thlaspi_consensi*.aln.raxml.reduced.phy ; do 
  echo  " working with .... ${reduced_phy} " ; 
  raxml-ng -all -msa  ${reduced_phy} -model JTT+G  --bs-trees 100 ;
done 

          
