#!/usr/bin/env bash

usage() { echo "$0
# This scripts retrieves the coverage  at the insertion site, the 100pb upstream and the 100pb dowstream  of it  from the global
# FINAL.NC.DP3.bed file generated by the SPlitreader filtering step and the final, filtered insertion file of a given sample (.acc.filt.DP3.bed).
# Also adds the mean coverage by 10kb windows generated by mosdepth "mosdepth -m -b 10000 -n -t 8  $(basename $bamfie .bam)_by10kb $bamfile"
# to each insertion position, to generate a more local coverage value and it calcualtes the total sample mean to give a file to STDOUT like:
#
# chrom start end ins_TE_fam TIP_cov Ref_cov Upstream_100cov downstream_100cov local_10kb_cov sample_mean_cov
#
# Usage:" && grep " .)\ #" $0; exit 0; }

[ $# -eq 0 ] && usage
while getopts ":hf:i:r:" arg; do
  case $arg in
    f) # path to  the splitreader filtering FINAL.NC.DP3.bed file
      spl_filt_file=${OPTARG}
      ;;
    i) # Path of the  FILTERED insertion file of a sample ( ${sample_name} .lifted.DP3.bed).
      individual_insertions_file=${OPTARG}
      ;;
    r) # Path of the median coverage by 10kb wins regions file calculated by mosdepth (median_by10kb.regions.bed.gz)
      region_file10kb=${OPTARG}
      ;;
    h | *) # Display help.
      usage
      exit 0
      ;;
  esac
done


# Get date and names of the individual, set temp  files
ind_name=$(basename ${individual_insertions_file} .bed)
local_cov_file=$(mktemp --suffix=.tsv tmp_local_cov_file.XXXXXXXXXX)
final_cov_file=$(mktemp --suffix=.tsv tmp_final_cov_file.XXXXXXXXXX)
cov_file_merged=$(mktemp --suffix=.tsv tmp_merged_cov_file.XXXXXXXXXX)

# Sample queued:
echo " Working with sample ${ind_name}"

# Calcualte the median for the sample file
sample_median=$(zcat ${region_file10kb} |
cut -f4 |
sort -n |
awk ' { a[i++]=$1; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }' )


# Loop over the file with awk to find the matching column with the name ID and save chrom start end TEfam and cov info:

awk -v ind_name=$ind_name 'NR==1 { for (i=1;i<=NF;++i) if ($i==ind_name)  { n=i; break }} { print $1"\t"$2"\t"$3"\t"$4"\t"$n }' ${spl_filt_file} > ${local_cov_file}

# Local cov file is a slash-separated 4 field  column with  integers values as following:
# Coverage supporting the non-reference insertion /Coverage supporting the reference absence / Coverage 100bp up / Coverage100bp down.


# Retrieve the filtered allele coverage info and retain the first two fields in  a temp file
tail -n+2 ${local_cov_file} |
bedtools intersect -wb -f 1 -a ${individual_insertions_file} -b - |
awk  'BEGIN{OFS=FS="\t"} {if($4==$10 && $2==$8 && $3==$9 && $1==$7 && $6==1) print $1,$2,$3,$4,$5,$11}' |
sed 's/\//\t/g'  > ${final_cov_file}


# Add the local  coverage and the mean
bedtools intersect -wb -f 1 -a  ${final_cov_file} -b ${region_file10kb} |
awk -v median=$sample_median 'BEGIN{OFS=FS="\t"}{print $1,$2,$3,$4,$6,$7,$8,$9,$13,median}'  |
cat > ${cov_file_merged}

# Add name and  column names
awk -v ind_name=$ind_name 'BEGIN{OFS=FS="\t"} {print $0,ind_name}' ${cov_file_merged} |
sed '1 i\chrom\tstart\tend\tTIP_TE_fam\tTIP_cov\tref_cov\tupstream_100cov\tdownstream_100cov\tlocal_10kb_median_cov\tsample_median_cov\tAccession' |
cat >  ${ind_name}_cov_TIPs_eval.tsv

# Cleanup
rm ${local_cov_file}
rm ${cov_file_merged}
rm ${final_cov_file}