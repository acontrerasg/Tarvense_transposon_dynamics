#!/usr/bin/env python
#
# Thom Nelson's adaptation of Tomaz Berisa's script to convert plink clustered allele frequencies to treemix format
#   (https://bitbucket.org/nygcresearch/treemix/downloads/)
#
import sys
import pandas as pd

if len(sys.argv) != 2:
    sys.stderr.write("\n")
    sys.stderr.write("Usage: plink2treemix.py plink.frq.strat\n")
    sys.stderr.write(" (script will write uncompressed output to stdout)\n")
    sys.stderr.write("\n")
    sys.exit()

plink_freq_file = sys.argv[1]

#freq_file = "/home/thom_nelson/LewCard/seq/popgenetics/mapped2CE10chromonome/plink/diversity_GATK_CE10chromonome.8.variable.SNPs.frq.strat_toy"

freq_df = pd.read_csv(plink_freq_file, delim_whitespace=True)

snp_IDs = freq_df["SNP"].unique()
pop_IDs = freq_df["CLST"].unique()
nsnps = len(snp_IDs)
npops = len(pop_IDs)

sys.stderr.write("SNPs in file: %s\n" % nsnps)
sys.stderr.write("Populations in input (n = %s):\n" % npops)
[sys.stderr.write("  %s\n" % pop) for pop in pop_IDs]

pop_IDs_out = " ".join(pop_IDs)
sys.stdout.write("%s\n" % pop_IDs_out)

nsnps_out = 0
sys.stderr.write("\n")
for snp in snp_IDs:
    snp_df = freq_df[freq_df["SNP"] == snp]
    unique_freqs = snp_df["MAC"].unique().tolist()
    # enforce presence of site in all pops
    if snp_df.shape[0] != npops:
        continue
    nsnps_out += 1
    sys.stderr.write("SNPs written: %s\r" % nsnps_out)
    allele_counts_per_pop = []
    for pop in pop_IDs:
        pop_df = snp_df[snp_df["CLST"] == pop]
        nalleles = pop_df["NCHROBS"].tolist()[0]
        minor_count = pop_df["MAC"].tolist()[0]
        major_count = nalleles - minor_count
        minor_count = str(minor_count)
        major_count = str(major_count)
        allele_counts = ",".join([major_count,minor_count])
        allele_counts_per_pop.append(allele_counts)
    allele_counts_str = " ".join(allele_counts_per_pop)
    sys.stdout.write("%s\n" % allele_counts_str)
sys.stderr.write("\n")
