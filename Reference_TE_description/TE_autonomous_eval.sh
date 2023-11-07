#!/usr/bin/env bash

# Here we use the output of DANTE to evaluate wether a given TE family is
# either autonomous (ie: it codifies for all the machinery needed to its transpositon)  or 
# non-autonomous (ie : it doesnt)

## Path and env

conda activate /ebio/abt6_projects8/Tarvense_TE/code/conda_envs/LAST_tool
cd /ebio/abt6_projects8/Tarvense_TE/data/Tarvense_TE_pops/Thlaspi_genome_2021/Thlaspi_paper_final_genome/TE_proteins/
mkdir -p new_eval
cd new_eval
#Dante scripts
dante="/ebio/abt6_projects8/Tarvense_TE/code/dante/dante.py"
dante_filt="/ebio/abt6_projects8/Tarvense_TE/code/dante/dante_gff_output_filtering.py"
#Dante DB
viriplantae="/ebio/abt6_projects8/Tarvense_TE/code/dante/tool-data/protein_domains/Viridiplantae_v3.0_pdb"
class_file="/ebio/abt6_projects8/Tarvense_TE/code/dante/tool-data/protein_domains/Viridiplantae_v3.0_class"
#Thlaspi consensi and annotation
thlaspi_TE_consensi="/ebio/abt6_projects8/Tarvense_TE/data/Tarvense_TE_pops/Thlaspi_genome_2021/thlaspi.TE.fa"
annotation="/ebio/abt6_projects8/Tarvense_TE/data/Tarvense_TE_pops/Thlaspi_genome_2021/Thlaspi_paper_final_genome/sorted.final.TEs_lifted.gff"

# Dante run
mkdir -p dante_tmp
python3 ${dante} -q  ${thlaspi_TE_consensi} -pdb ${viriplantae} --domain_gff thlaspi_consensi_class_dante.gff3 -dir dante_tmp  -cs ${class_file}
python3 ${dante_filt} --dom_gff ./dante_tmp/thlaspi_consensi_class_dante.gff3  \
                      -ouf ./thlaspi_consensi_class_dante_filt_ir10.gff3 \
                      -dir ./ \
                      -ir 10


# Autonomy evaluation: 
# We further divided families as autonomous or non autonomous based on these results.
# DNA families with  their consensus sequence containing a TPase (or HEL1/HEL2 in case of helitrons) were copnsidered autonomous.
# NON-LTR families (LINES) just need a reverse transcriptase  (RT) )
# and  LTR families with GAG, protease (PROT), reverse transcriptase (RT), ribonuclease H (RH) and integrase (INT). were assigned as autonomous and the rest as non autonomous.

# Separate dante output in different file per superfam:
for superfam in $(cut -f3 $annotation | sort | uniq) ;  do 
  cat thlaspi_consensi_class_dante_filt_ir10.gff3 | 
  grep "^#" | 
  grep "$superfam" | 
  sed 's/##//'  >  dante_prot_${superfam}_out.txt ; 
done 



# Loop over the superfamilies and classify the families in autonomouns and non-autionomouns base on the presence of proteins:

for superfam in $(cut -f3 $annotation | sort | uniq) ;  do
   #Helitrons
      if  [[ $superfam == "DHH" ]] ; then
          echo "${superfam} is DHH"
              for fam in $(grep ${superfam} ${annotation}   | cut -f9 | cut -f2 -d ';' | sed 's/Target=//' | sort | uniq) ; do           
                Hel_found=$(grep -w $fam dante_prot_${superfam}_out.txt | grep -ic "hel")
     
                if [ ${Hel_found} -ne 0 ] ;  then
                  echo -e "${fam}\tautonomous\t${superfam}"
                else
                  echo -e "${fam}\tnon-autonomous\t${superfam}"
                fi
          done > ${superfam}_mobility_summary.txt    
  #Retrotransposons
      elif  [[ $superfam =~ (RL) ]] ; then
        echo "${superfam} is a LTR retrotransposon"
        for fam in  $(grep ${superfam} ${annotation} | cut -f9 | cut -f2 -d ';' | sed 's/Target=//' | sort | uniq) ; do
              #Search presence of candidate protein profiles:
              RT=$(grep -w $fam dante_prot_${superfam}_out.txt | grep -ic "RT")
              GAG=$(grep -w $fam dante_prot_${superfam}_out.txt  | grep -ic "GAG")
              RH=$(grep -w $fam dante_prot_${superfam}_out.txt  | grep -ic "RH")
              #Ty3 caspid is an pesky one as its similar to ASP and AP both the protease I am searching for
              PROT=$(grep -w $fam dante_prot_${superfam}_out.txt  | grep -ic "PROT")
              INT=$(grep -w $fam dante_prot_${superfam}_out.txt  | grep -ic "INT")
              if [ ${RT} -ne 0 ] && [ ${GAG} -ne 0 ] && [ ${RH} -ne 0 ] && [ ${PROT} -ne 0 ] && [ ${INT} -ne 0 ]  ;  then
                  echo -e "${fam}\tautonomous\t${superfam}"
              else
                  echo -e "${fam}\tnon-autonomous\t${superfam}"
              fi
       done > ${superfam}_mobility_summary.txt
  #LINES
     elif  [[ $superfam  =~ (RI) ]] ; then
      echo "${superfam} is a LINE"
      for fam in $(grep ${superfam} ${annotation} | cut -f9 | cut -f2 -d ';' | sed 's/Target=//' | sort | uniq) ; do
            #find the reverse transcriptase for Lines 
            RT=$(grep -w $fam dante_prot_${superfam}_out.txt | grep -ic "RT" )
            if [ ${RT} -ne 0 ] ;  then
              echo -e "${fam}\tautonomous\t${superfam}"
          else
              echo -e "${fam}\tnon-autonomous\t${superfam}"
          fi
      done  > ${superfam}_mobility_summary.txt
 #DNA transposons
    elif  [[ $superfam  =~ (DT) ]] ; then
      echo "${superfam} is DNA transposon"
      for fam in $(grep ${superfam} ${annotation} | cut -f9 | cut -f2 -d ';' | sed 's/Target=//' | sort | uniq) ; do
          #Search for transposase proteins
          transposase_found=$(grep -w $fam  dante_prot_${superfam}_out.txt | grep -ic "TPase")
          if [ ${transposase_found} -ne 0 ] ;  then
            echo -e "${fam}\tautonomous\t${superfam}"
          else
            echo -e "${fam}\tnon-autonomous\t${superfam}"
          fi
      done > ${superfam}_mobility_summary.txt
    fi
done


#Final dataset:
cat *mobility_summary.txt > ../TE_proteins_mobility_eval.tsv
