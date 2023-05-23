#!/bin/usr/env bash

#Arguments:

METADATA=$1 # Assumes first line is header.
vcf_file=$2 # VCF file has to be bgzipped and tabix indexed

vcf_samples=$(mktemp --suffix=.txt vcf-samples.XXXXXXXX)
meta_samples_region=$(mktemp --suffix=.txt meta_samples_region.XXXXXXXX)

#get the names of the samples in the same order as columns in the vcf (not sure if the order is needed  but just in case).
bcftools query -l ${vcf_file} > ${vcf_samples}

tail -n +2  ${METADATA} |  cut -f1,3  > ${meta_samples_region}

#Order by hash in awk
awk 'BEGIN{OFS=FS="\t"}
NR==FNR { a[$1]=$0; next }
$1 in a { print a[$1], $2, $3 }'  ${meta_samples_region}  ${vcf_samples} | sed 's/\t$//;s/\t$//' | awk '{print $1"\t"$1"\t"$2}' >  data.clust


rm $vcf_samples
rm $meta_samples_region

