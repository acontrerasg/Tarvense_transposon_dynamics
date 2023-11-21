SNP CALLING

Scripts for calling short variants in GATK4.<br/>
Workflow for read mapping (bwa), short variant calling (GATK4) and filtering from short-read Whole Genome Sequencing data, for large datasets. <br/>
System setup: linux-based cluster with PBS queueing system.
The original repository is deposited at [BinAC_varcalling]([https://github.com/Dario-Galanti/multipheno_GWAS/tree/main/gemmaGWAS](https://github.com/Dario-Galanti/BinAC_varcalling/tree/main))

The workflow is meant for the analysis of paired-end short reads. Based on the GATK4 best practices [for germline short variant discovery](https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels-), is specifically customized for fairly large WGS datasets of non-model species. It can handle a large sample numbers and large and fragmented genomes. <br/>
<br/> 

WORKFLOW DESCRIPTION:

[0_multi_trim_BinAC.sh](https://github.com/Dario-Galanti/BinAC_varcalling/blob/main/0_multi_trim_BinAC.sh) <br/>
Sample-parallelized adaptors trimming with [cutadapt](https://cutadapt.readthedocs.io/en/stable/).

[1_BWA_multi_align_BinAC.sh](https://github.com/Dario-Galanti/BinAC_varcalling/blob/main/1_BWA_multi_align_BinAC.sh) <br/>
Sample-parallelized read alignment with [bwa-mem](http://bio-bwa.sourceforge.net/bwa.shtml) and detection of duplicates with [MarkDuplicatesSpark](https://gatk.broadinstitute.org/hc/en-us/articles/360046221811-MarkDuplicatesSpark).

[2_HaploCaller_multi_BinAC.sh](https://github.com/Dario-Galanti/BinAC_varcalling/blob/main/2_HaploCaller_multi_BinAC.sh) <br/>
Sample-parallelized variant calling with local reassembly ([Haplotypecaller](https://gatk.broadinstitute.org/hc/en-us/articles/360036715891-HaplotypeCaller)) to obtain single-sample GVCF files.

[3_4_GenDB_GenoGVCFs_BinAC.sh](https://github.com/Dario-Galanti/BinAC_varcalling/blob/main/3_4_GenDB_GenoGVCFs_BinAC.sh) <br/>
Combine single-sample GVCF files in a multisample vcf file (actually files, one per scaffold). To avoid large datasets from causing RAM saturation, this step is parallelized by scaffold and is composed of two steps: 1) Importing all samples in a GenomicsDB ([GenomicsDBImport](https://gatk.broadinstitute.org/hc/en-us/articles/360036732771-GenomicsDBImport)) and 2) joint genotyping of all samples ([GenotypeGVCFs](https://gatk.broadinstitute.org/hc/en-us/articles/360036348452-GenotypeGVCFs)).
    
[5_GatherVcfs_BinAC.sh](https://github.com/Dario-Galanti/BinAC_varcalling/blob/main/5_GatherVcfs_BinAC.sh) <br/>
Combine all single-scaffold vcf files in a unique multisample vcf file ([GatherVcfs](https://gatk.broadinstitute.org/hc/en-us/articles/360036711811-GatherVcfs-Picard-)).

[6_filt_spls_BinAC.sh](https://github.com/Dario-Galanti/BinAC_varcalling/blob/main/6_filt_spls_BinAC.sh) <br/>
Samples filtering ([SelectVariants](https://gatk.broadinstitute.org/hc/en-us/articles/360040508071-SelectVariants)). This allows to run downsteam analysis with different sets of samples without repeating the variant calling.

[7_filt_variants_BinAC.sh](https://github.com/Dario-Galanti/BinAC_varcalling/blob/main/7_filt_variants_BinAC.sh) <br/>
Variants filtering: 1) remove low quality variants ([VariantFiltration](https://gatk.broadinstitute.org/hc/en-us/articles/360036350452-VariantFiltration)) using different paramenters for SNPs and other variants, 2) remove variants with missing genotypes in more than a user-defined proportion of samples ([VCFtools](https://vcftools.github.io/man_latest.html)), 3) remove multiallelic sites and filter for Minor Allele Frequency and 4) remove scaffolds harbouring less than 3 variants ([lonely_vcf_pos.sh](https://github.com/Dario-Galanti/BinAC_varcalling/blob/main/lonely_vcf_pos.sh)) as these cause problems during phasing and imputation.

[vcf_impute_BinAC.sh](https://github.com/Dario-Galanti/BinAC_varcalling/blob/main/vcf_impute_BinAC.sh) <br/>
Phasing and imputation of missing genotype calls with BEAGLE. In addition this script adds the reference genotype to the vcf file.
