#!/usr/bin/env bash
set -e
set -u
set -o pipefail
set -v


# This script  is a hardcoded script that feeds IGV_snap with different bams and coordinates to evaulate TIPs from Thlaspi arvense.
# It also means I need to use gshuf gawk and stuff like that.
# It runs locally in my mac, therefore I need to mount first the project folders and the global_tmp folders where I have stored the data.
# Therefore run locally:
# $ mount_temp
# $ mountproject2

# Main script:
IGV_snap="/Users/acontreras/Tarvense_TE/code/IGV_snap.sh"

# Load files
genome='/Users/acontreras/Tarvense_TE/data/Tarvense_TE_pops/Thlaspi_genome_2021/Thlaspi_paper_final_genome/modified.fasta'
TE_FILE='/Users/acontreras/Tarvense_TE/data/Tarvense_TE_pops/Thlaspi_genome_2021/Thlaspi_paper_final_genome/sorted.final.TEs_lifted.gff'
GENE_FILE='/Users/acontreras/Tarvense_TE/data/Tarvense_TE_pops/Thlaspi_genome_2021/Thlaspi_paper_final_genome/sorted.final.gene_lifted.gff'
TIP_FILE='/Users/acontreras/Tarvense_TE/data/Tarvense_TE_pops/FINAL_ANALYSIS/TIPS/All_TIPS_cov_corrected.bed'
output_dir='/Users/acontreras/Tarvense_TE/doc/VISUAL_TIP_IGV_VALIDATION'

#Initiate tmp files
temp_TIPshut=$(mktemp temp_TIPshut.XXXXXXXX)

for pop  in FR CN AM DE SE NL US ; do
  # Make a list of sample names:
  mkdir -p  ${output_dir}/${pop} ;
  mkdir -p  ${output_dir}/${pop}/True_pos/ ;
  mkdir -p  ${output_dir}/${pop}/False_pos/  ;
  for bam in  /Users/acontreras/cluster_tmp/Tarvense_TE_pops/bwa_mem_aligments/${pop}_output/TA*bam ; do
    name=$(basename $bam .bam) ;
    mkdir -p ${output_dir}/${pop}/${name} ;
    echo "sample is ${name} " ;
    # Get the correspoing TIP insertions for the sample, shuffle them and store the first 10 lines in a temporal bed file :
    gawk -v name=$name 'BEGIN{OFS=FS="\t"} {if($7==name) print $0}' $TIP_FILE |
    gshuf | head | tee ${temp_TIPshut} ;
    # Do the screenshots
    bash ${IGV_snap} -p ${temp_TIPshut} \
        -z 100 \
        -g ${genome} \
        -b ${bam} \
        -a ${temp_TIPshut} \
        -o ${output_dir}/${pop}/${name}

    rm ${temp_TIPshut} ;
  done ;
done

#cleanup
rm ${temp_TIPshut}
