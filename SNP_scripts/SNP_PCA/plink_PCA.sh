#!/bin/usr/env bash
set -e
set -o pipefail
set -v

VCF=$VCF 

plink   --vcf $VCF \
        --double-id \
        --allow-extra-chr  \
        --list-duplicate-vars suppress-first \
        --indep-pairwise 50 10 0.1 \
        --out Tarvense_SNPs_pruned_plink

plink   --vcf $VCF \
        --double-id \
        --allow-extra-chr  \
        --list-duplicate-vars suppress-first \
        --extract Tarvense_SNPs_pruned_plink.prune.in \
        --make-bed --pca --out Tarvense_SNPs_PCA
   
# Run R script to visu the PCA
Rscript ./visu_SNP_PCA.R
