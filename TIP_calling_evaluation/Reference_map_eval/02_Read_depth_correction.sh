#!/usr/bin/env bash

#Calculate the Read depth with correction of GC content as stated in Yoon et al 2009.
#Uses bioawk, bedtools and mosdepth.


usage() { echo "$0 usage:" && grep " .)\ #" $0; exit 0; }

[ $# -eq 0 ] && usage

while getopts ":hg:w:b:t:" arg; do
  case $arg in
    g) # path to genome size
      genome=${OPTARG}
      ;;
    w) # Window size (in integers)
      window_size=${OPTARG}
      ;;
    b) # path to bam file
      bam=${OPTARG}
      ;;
    t) # threads (in integers)
      threads=${OPTARG}
      ;;
    h | *) # Display help.
      usage
      exit 0
      ;;
  esac
done

#Get names:
gen_name=`echo "${genome##*/}"`
gen_name=`echo "${gen_name%%.*}"`
bam_name=`echo "${bam##*/}"`
bam_name=`echo "${bam_name%%.*}"`
win_size_str=`echo $window_size"_bp"`

#Temp file names:
temp_chrom_sizes=$(mktemp --suffix=.txt temp_chromsizes.XXXXXXXX)
temp_win_gen=$(mktemp --suffix=.bed temp_windowed_genome.XXXXXXXX)
temp_GC_calc=$(mktemp --suffix=.bed temp_windowed_genome.XXXXXXXX)
temp_GC_perwin=$(mktemp --suffix=.bed temp_windowed_GCarray.XXXXXXXX)
temp_median=$(mktemp --suffix=.txt temp_median.XXXXXXXX)
temp_mCG=$(mktemp --suffix=.txt temp_mCG.XXXXXXXX)



#####Main    ######
#calcualte length of each chromosome
bioawk -c fastx '{ print $name, length($seq) }' ${genome} >  ${temp_chrom_sizes}
#Slice the chromosome in non overalpping windows of $window_size
bedtools makewindows -g ${temp_chrom_sizes} -w ${window_size} > ${temp_win_gen}
#Calculate GC content of genome 
bedtools  nuc -fi ${genome} -bed ${temp_win_gen} | tail -n +2 | awk '{printf "%.2f\n", $5}' >  ${temp_GC_calc}
paste ${temp_win_gen} ${temp_GC_calc} > ${gen_name}_${win_size_str}_GC.bed
#Calculate median read coverage per window:
mosdepth -t ${threads} -m -n -b ${gen_name}_${win_size_str}_GC.bed ${bam_name}.depthGCcorr $bam
#GetGCperCentwinds
cat ${gen_name}_${win_size_str}_GC.bed | cut -f4 | sort | uniq > ${temp_GC_perwin}
#calculate overall median
zcat ${bam_name}.depthGCcorr.regions.bed.gz  | cut -f5 | awk ' { a[i++]=$1; } END { print a[int(i/2)]; }' > ${temp_median}
#calculate median reads per GC fraction value:
for value in `cat ${temp_GC_perwin}` ;
do
zcat  ${bam_name}.depthGCcorr.regions.bed.gz |
awk -v value=$value  '{if($4 == value) print ($5+1) }' |
sort -n |
awk  -v GCvalue=$value '{ a[i++]=$1; }END { x=int((i+1)/2);
if (x < (i+1)/2) print GCvalue"\t"(a[x-1]+a[x])/2;
else print GCvalue"\t"a[x-1]; }' ;
done > ${temp_mCG}
#Get GC corrected  coverage:
zcat  ${bam_name}.depthGCcorr.regions.bed.gz |
while read line;
do
total_median=`cat ${temp_median}` ;
search_mGC=`echo $line | awk '{print $4}'` ;
medianGC=`grep -P "^${search_mGC}\t" ${temp_mCG} | cut -f2`
echo $line $medianGC |
awk -v total_median=$total_median -v medianGC=$medianGC 'BEGIN{OFS="\t" ; FS=" "}
{RCGCcorr=(total_median/medianGC) ;  print $1,$2,$3,$4,$5,($5*RCGCcorr)}'
done  |  sed '1 i #Chrom\tStart\tEnd\tGC\tRC\tRCcorr' > ${bam_name}.depthGCcorrected.regions.bed

##remove used temp files:
rm ${temp_chrom_sizes}
rm ${temp_win_gen}
rm ${temp_GC_calc}
rm ${temp_GC_perwin}
rm ${temp_median}
rm ${temp_mCG}
rm ${bam_name}.depthGCcorr.mosdepth.summary.txt
rm ${bam_name}.depthGCcorr.mosdepth.global.dist.txt
rm ${bam_name}.depthGCcorr.regions.bed.gz
