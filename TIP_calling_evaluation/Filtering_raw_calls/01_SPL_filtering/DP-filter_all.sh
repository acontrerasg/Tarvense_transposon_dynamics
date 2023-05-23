#!/usr/bin/env bash

set -u

######### Filtering Non-referenceTE Presence Variants by Positive Coverage (DP-filter).  #######
##Adapted from:
# Efficient Detection of Transposable Element Insertion Polymorphisms Between Genomes Using Short-ReadSequencing Data.
# Authors Pierre Baduel, Leandro Quadrana, and Vincent Colot.

usage() { echo "$0 Usage:
This scripts assumes that " && grep " .)\ #" $0; exit 0; }
[ $# -eq 0 ] && usage

while getopts ":hB:T:P:D:I:N:" arg; do
  case $arg in
    B) # Path to   folder with the  raw bed output of splitreader. (${SAMPLE}-insertion-sites.bed)
      bed_path=${OPTARG}
      ;;
    T) # Path to  TE families and TSD information.
      TSD_info=${OPTARG}
      ;;
    P) # Path to   the Filter_insertions_splitreader.pl script
      Filter_pl=${OPTARG}
      ;;
   D) ## Minimun number of spolit or discordant reads that support an insertion.
      depth=${OPTARG}
      ;;
    I) ##The insertion size matches the expected TSD size + $TSDthresh bp.
      TSDthresh=${OPTARG}
      ;;
    N) # #The  insertion  size  is  not  longer  than $noTSDthresh bp for TE families that do not generateTSD upon insertion.
       noTSDthresh=${OPTARG}
      ;;
    h | *) # Display help.
      usage
      exit 0
      ;;
  esac
done

### Fix  the output of ther bed files and sort them:
for  bed_results in  ${bed_path}*-insertion-sites.bed ;
do
awk 'BEGIN{OFS=FS="\t"} {if(NF == 9) print $1 , $2 , $3, $4 , $5 , $6 , $7, $8, $9, "NA", "NA" ;  else print $0}' $bed_results |
bedtools  sort -i - >  $(basename  $bed_results .bed)_sorted.bed
done

#Prepare files by TE family to feed into the perl script one family at the time.
for TE_fam in $( cut -f1 ${TSD_info} ) ;
do
  for bed_results_sorted in   *-insertion-sites_sorted.bed ;
        do
        name=$(basename $bed_results_sorted  _ALL-insertion-sites_sorted.bed) ;
        grep -w ${TE_fam} ${bed_results_sorted} >  ${TE_fam}-${name}-insertion-sites.bed
        grep -w ${TE_fam} ${bed_results_sorted}  >>  ${TE_fam}_intersect.sort.bed ;
  done
  
  
bedtools sort -i ${TE_fam}_intersect.sort.bed | bedtools merge -i - >  ${TE_fam}_intersect.mrg.bed ;
multiIntersectBed -header -i ${TE_fam}-*-insertion-sites.bed >   ${TE_fam}_intersect.sort.bed
sed -i 's/_sorted-TE_ALL-insertion-sites_sorted.bed//g'  ${TE_fam}_intersect.sort.bed
sed -i "s/${TE_fam}-//g"  ${TE_fam}_intersect.sort.bed

mkdir -p ./out_filtered_DP${depth}/
mkdir -p ./out_filtered_DP${depth}_info/

perl ${Filter_pl}  ${depth} ${TSD_info} ${TE_fam} ${TSDthresh} ${noTSDthresh}  out_filtered_DP${depth}/  ${TE_fam}-*-insertion-sites.bed  >  ./out_filtered_DP${depth}_info/${TE_fam}_info.tsv;

##Cleanup
rm  ${TE_fam}-*insertion-sites.bed
rm  ${TE_fam}_intersect.mrg.bed
rm  ${TE_fam}_intersect.sort.bed
done

#remove intial transformed beds.
rm  TA*TE_ALL-insertion-sites_sorted.bed

#Work with the ouput data:

for i in out_filtered_DP${depth}/*_insertions.filt.DP3.bed;
do
head -n1  $i > header.tsv
break
done

for i in out_filtered_DP${depth}/*_insertions.filt.DP3.bed;
do
grep -v "^#chrom" $i >> tmp.tsv
done

cat header.tsv tmp.tsv >  Tarvense_filtered_DP${depth}.bed
#cleanup
rm header.tsv
rm tmp.tsv


#Preparing  the input data for negative coverage filter:
#Insertions sites
cat  Tarvense_filtered_DP${depth}.bed | cut -f 1-5 > Tarvense-insertion-sites.0.bed
#Insertion sites 100bp upstream (-100)
cat  Tarvense_filtered_DP${depth}.bed | cut -f 1-5 |
bedtools shift -i - -g /ebio/abt6_projects8/Tarvense_TE/data/Tarvense_TE_pops/Thlaspi_genome_2021/thlaspi_chrom_sizes.txt  -s -100 > Tarvense-insertion-sites.100up.bed
#Insertion sites 100bp downstream (-100)
cat  Tarvense_filtered_DP${depth}.bed | cut -f 1-5 |
bedtools shift -i - -g /ebio/abt6_projects8/Tarvense_TE/data/Tarvense_TE_pops/Thlaspi_genome_2021/thlaspi_chrom_sizes.txt  -s 100 > Tarvense-insertion-sites.100down.bed


#Work with the info:
cat ./out_filtered_DP${depth}_info/*_info.tsv | bgzip > ALL_info.tsv.gz
rm ./out_filtered_DP${depth}_info/*_info.tsv
rmdir ./out_filtered_DP${depth}_info/

echo "DONE"
