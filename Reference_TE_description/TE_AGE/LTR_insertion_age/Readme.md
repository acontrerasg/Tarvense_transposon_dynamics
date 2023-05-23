## Get the LTR insertion age estimation. 

The long terminal repeats at both the 5’ and 3 are identical upon insertion,  
therefore we  can estimate  their insertion time using the number of substitutions 
that occur between the two repeats. 

In order to do that, we need to assess the structre of each LTR copy. 
We first used LTR-harvest in their de novo mode against the genome of thlaspi arvense with the following parameters:

 > LTR_HARVEST args= -minlenltr 100 -maxlenltr 7000 -mintsd 4 -maxtsd 6 -motif TGCA -motifmis 1 -similar 85 -vic 10 -seed 20 -seqids yes

Then  we retrive the 5p and 3p LTR  fragments for each copy from the cns output file produced by ltr-harvest and We alignt both fragments  using  Mafft with the following script:

[LTR_harvest_getconsensus.sh](https://github.com/acontrerasg/Tarvense-TEpop/blob/main/TE_AGE/LTR_insertion_age/LTR_harvest_getconsensus.sh)

Then we calculate nucleotide divergence with several models using dna.dist in the ape package of R:
K2P model 
Raw model
tn93 model 
We also retrieve the LTR age prediction that comes by defaul with LTR harvest (but it is fine tuned for rice):
[LTR_insertion_age.R](https://github.com/acontrerasg/Tarvense-TEpop/blob/main/TE_AGE/LTR_insertion_age/LTR_insertion_age.R)

Finally, with the file that the above Rscripts outpus,  
we try to realate the  LTRs predicted by LTR-harvest with the LTR present in the original, denested annotation:

[LTR_age2anno.sh](https://github.com/acontrerasg/Tarvense-TEpop/blob/main/TE_AGE/LTR_insertion_age/LTR_age2anno.sh)

This last scripts outputs the amount of overlap between the denested annotaiton and the output of LTR harvest. 
So one can use the row key in the  1st fiel of the final file "LTRs_tarvense_ages_denested_overlap.tsv"  
to asign a unique stimation value to the TE faimily with  the most amount of overlap in the 13th field.

For all age measures, we relate nucleotide divergence to absolute time using a mutation rate from arabidopsis thaliana **not done yet**
(WE DONT HAVE IT YET FOR TARVENSE WE NEED TO GET IT).
