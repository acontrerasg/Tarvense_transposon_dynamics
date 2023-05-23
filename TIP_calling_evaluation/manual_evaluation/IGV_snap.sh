#!/bin/env bash
set -e
set -u
set -o pipefail
set -v

usage() { echo "$0 usage:
This script  makes snapshots on IGV  at each position given in the  bed file at a specific zoom level (given by nucleotides)
at the genome indicated. It also loads other kind of annotation for the visualization.
" && grep " .)\ #" $0; exit 0; }

[ $# -eq 0 ] && usage

while getopts ":hp:z:g:b:a:o:" arg; do
  case $arg in
    p) # path to  bed file where the screenshots have to be taken
      bed_file=${OPTARG}
      ;;
    z) # ZOOM level (in bp) it should be an integer. Recomended 100. 
      ZOOM=${OPTARG}
      ;;
    g)  # path to  genome file  to load
      GENOME=${OPTARG}
      ;;
    b) # path  the bam file with the reads.
      BAMFILE=${OPTARG}
      ;;
    a) # Annotation files to load (ideally given as an array of files separated by a commas)
      ANNO_FILES=${OPTARG}
      ;;
    o) ## Path  to output directory to save the snapshoots.
      LOCALDIR=${OPTARG}
      ;;
    h | *) # Display help.
      usage
      exit 0
      ;;
  esac
done

## Temporal files

#for zooming
BED_FILE_ZOOM=$(mktemp temp_zoom_bed.XXXXXXXX)
#for  batch file
temp_batch_file=$(mktemp temp_bash_file.XXXXXXXX)
temp_header_file=$(mktemp temp_header_file.XXXXXXXX)
temp_anno_loader_file=$(mktemp temp_anno_file.XXXXXXXX)
temp_body_file=$(mktemp temp_body_file.XXXXXXXX)
temp_end_file=$(mktemp temp_end_file.XXXXXXXX)


#APPLY ZOOM:
bedtools slop -i ${bed_file} -g ${GENOME}.fai -b ${ZOOM} > ${BED_FILE_ZOOM}

# Name of the genome
GENOME_NAME=$(basename ${GENOME} )


# WRITE THE BATCH  FILE FOR IGV
printf "new\ngenome ${GENOME_NAME}\n
colorBy INSERT_SIZE\n
load ${BAMFILE}\n
snapshotDirectory ${LOCALDIR}\n
" > ${temp_header_file}
#I want to change this as an array so you can include several annotations but for now...
printf "load ${bed_file}\n" > ${temp_anno_loader_file}
awk 'BEGIN{FS="\t"}{print "goto "$1":"$2"-"$3"\nsort base\nsnapshot snapshot_"$1":"$2"-"$3"_"NR".png"}' ${BED_FILE_ZOOM} > ${temp_body_file}
printf "exit" > ${temp_end_file}

cat ${temp_header_file} ${temp_anno_loader_file} ${temp_body_file}  ${temp_end_file} > ${temp_batch_file}



#RUN THE SCRIPT
bash /Applications/IGV_2.12.3/igv.sh -b ${temp_batch_file}


## CLEANOUT
rm ${BED_FILE_ZOOM}
rm ${temp_batch_file}
rm ${temp_header_file}
rm ${temp_anno_loader_file}
rm ${temp_body_file}
rm ${temp_end_file}
