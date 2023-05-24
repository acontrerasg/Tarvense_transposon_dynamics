#!/bin/bash

# This script runs meme in ALESIA.FAM.7  full copies TEs.

# First we aligned all full copy alesia longer than 5kb:

full_copies="/ebio/abt6_projects8/Tarvense_TE/data/Tarvense_ALESIA_TE/RLC.FAM.7_fullcopies.fa"

mafft --adjustdirection --maxiterate 1000 --globalpair --thread 16  ${full_copies} > full_copes.aln

# WE manually removed sequences that presented small deletions and we end up with a final list of 39 copies:
curated_full_copies.name.list.txt
seqkit  grep -nf curated_full_copies.name.list.txt ../RLC.FAM.7_fullcopies.fa  > curated_full_copies.fa
# Then we ran meme with the  resulting final 39 copies:
conda activate /ebio/abt6_projects8/Tarvense_TE/code/conda_envs/meme_suite/
mkdir -p MEME_SUITE
cd MEME_SUITE
# Run meme  for the full ALESIA copies:
mkdir -p meme_full_copies 
meme curated_full_copies.fa -dna -oc meme_full_copies -mod zoops -nmotifs 4 -minw 6 -maxw 50 -objfun classic -revcomp -markov_order 01

# We also run it for all of them 
meme RLC.FAM.7_fullcopies.fa -dna -oc meme_full_copies -nmotifs 4 -minw 6 -maxw 50

