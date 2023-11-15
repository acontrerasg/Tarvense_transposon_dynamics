Age estimation of TIPs

Contains two scripts to estimate the age of TIPs in large sequencing datasets, starting from a matrix of presence absence of TIPs

[AgeDeterm_TIPs.sh](https://github.com/acontrerasg/Tarvense_transposon_dynamics/blob/main/TIPs_age/AgeDeterm_TIPs.sh) <br/>
Estimates age of TIPs with Minor Allele Count >= 2 (non singletones), based on the maximum pairwise divergence (number of SNPs) between each combination of two carriers, in the 70 kb region around the insertion (as in [Baduel et al. 2021](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02348-5))


[AgeDeterm_singleTIPs.sh](https://github.com/acontrerasg/Tarvense_transposon_dynamics/blob/main/TIPs_age/AgeDeterm_singleTIPs.sh) <br/>
Estimates age of singletone TIPs, based on the number of private SNPs in the 70 kb region around the insertion. Most singletones are recent insertions, but this can give an idea whether there might be some old ones.

