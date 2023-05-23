## PERL SCRIPTS for SPL filtering. 

Here are my modified version of the perl scripts  from: 
https://github.com/baduelp/public/tree/master/SPLITREADER .

I also place here the bash wrappers that call them. 
Those modifications are mostly in  the  data flow managment part.

Except in the script Filter_negative_calls_splitreader.pl  where I  change the  minimun number of individuals to be covered.
From 10% of the  total cohort as preset, to a minum of 4 individuals to support an insertion. 
My reasoning was that in the more divergent population (Armenia)  I have only 7 individuals,
therefore it will never reach a 10% of the overall total. 


#######################################################

### Step by step  Protocol and rationale.

From: Efficient Detection of Transposable Element Insertion Polymorphisms Between Genomes Using Short-Read Sequencing Data. <br>

DOI: 10.1007/978-1-0716-1134-0_15.

####### <br>
**IMPORTANT**
Perl scripts presented in this filtering do not allow "-" symbols in the sample names. This will cause dowstream problems.
####### <br>

Putative non-reference TE presence variant sites called across each genome are then intersected and merged by TE family to define genomic regions where the presence variant calls overlap. These regions are further refined and eventually split into several non-reference TE presence variants based on splitreads, when available, as these give precise information (at a 1bp resolution) on the boundaries of non-reference TE presence variants.

 **RUN DP filter:** <br>

 > ./DP-filter_all.sh -B /all_raw_SPL_calls/*-insertion-sites.bed \
		-T  TSD_info.txt \
		-P Filter_insertions_splitreader.pl \
		-D 3 \
		-I 10 \
		-N 50


After the DP-filter, non-reference TE presence variants are then combined across TE-families and superfamilies to define the complete set of genomic sites for the second filter based this time on negative coverage NC-filter. For a diploid species with high homozygosity, the negative coverage in the alignment of a carrier genome onto the reference genome should drop by at least half where non-reference TE presence variants are present. <br>

 **Retrieve coverage at insertion sites, and surrounding regions:** <br>
 > ./readcount_wrapper.sh -i $sample_SPL1_bowtie2.bam \
 >  -g ${reference_genome} \
 >  -P Process_BAMrc_splitreader.pl \
 >  -d ${depth} \
 >   -F DP-filter_output_folder \
 >   -O RC_folder/  <br> 

**Run NC-filter  over the intial insertion sites.bed files for each sample: ** <br>

> ind_thresh=4 <br>
> depth=3 <br>
> perl Filter_negative_calls_splitreader.pl ${ind_thresh} ${depth} ${RC_folder} ${output_dir} DP_filter_output.bed *-insertion-sites.bed

Thus, the minimum negative coverage over the putative insertion site is compared to the average of the minimum negative coverage found within a similar-sized window 100 bp upstream and 100 bp downstream of the insertion site. <br>

If the negative coverage (NC) does not drop by at least half or if the total number of reads over a presence variant (both supporting a presence, i.e., positive coverage, and supporting the absence by aligning on the reference genome, i.e., negative coverage) is over 100 reads, the non-reference TE presence variant is considered as a false positive for the genome in question. <br>
If an individual genome is not a carrier but the negative coverage over the presence variant is under 5 reads, it is considered as missing information or NA. Non-reference TE presence variants are kept when they fulfill all criteria listed above and have a rate of missing individuals below 90%.

Final files explanation:

 1.  Tarvense-insertions.FINAL.NC.DP3.bed: <br>
        Contains all the NC information but is not filtered for the ratio of negative coverage (ratioNCfilt): 
        for each putative insertion sites it shows for each sample the coverage supporting the non-reference insertion / the coverage supporting the reference  absence / the coverage 100bp up / the coverage 100bp down.
        
 2.  Tarvense-insertions.ratioNC.FINAL.bool.DP3.bed: <br>
        Filters the putative insertions based on the ratioNCfilt with a 1 for samples where insertions have passed the filter, 
        a 0 where they have not and the reference coverage is sufficient (>$DPrefmin reads) to be confident that the insertion 
        is indeed absent, and - for the samples where the reference coverage is not sufficient to confirm the absence.
   
 3.  Tarvense-insertions.ratioNC.FINAL.DP3.bed: <br>
        Gives the coverage supporting the ratioNCfilt insertions.
        
 4.  Tarvense-insertions.ratioNC.FINAL.NConly.DP3.bed: <br>
        Gives the reference coverage over ratioNCfilt insertion sites.
        


If the negative coverage (NC) does not drop by at least half or if the total number of reads over a presence variant (both supporting a presence, i.e., positive coverage, and supporting the absence by aligning on the reference genome, i.e., negative coverage) is over 100 reads, the non-reference TE presence variant is considered as a false positive for the genome in question. If an individual genome is not a carrier but the negative coverage over the presence variant is under 5 reads, it is considered as missing information or NA. Non-reference TE presence variants are kept when they fulfill all criteria listed above and have a rate of missing individuals below 90%.


### LAST TWO costum filtering steps

On top of this, I also run two custom filtering steps not present in the Baduel protocol.

The first one, **acc_filtering step**, removes all those insertions which fall outside the accessibilty short reads regions evaluation as thosae are regions with abnormal short read mapping and therefore are a source of false positives. This script interacts over every bed file of the final call folder of the SPL raw filtering. It does two things:
1. It removes the DHH type TE calls. 
2. It removes any TE call outside of the accesibility regions file.  <br>
It also swtiches back from the Chr_ notation splitreader uses to the one the reference genome has (Scaffold_) and fixes the problem of files ending up with an extra tab at the end and also with a extra space at the header row.

The second one, **coverage eval**,  brings together coverage information from the different files produced by SPL and from a default aligment done with bwa mem using the raw reads. The idea is to filter out those TIPS whose local read coverage is higher than 100 reads.

