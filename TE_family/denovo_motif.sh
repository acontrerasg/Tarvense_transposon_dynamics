#!/bin/bash
# Here we perform a de novo search of putative HRE motifs using all Alesia.FAM.7 copies. We retrieve a candidate motif with a clsuter of four nGAAn, as reported in the literature. =
conda activate /ebio/abt6_projects8/Tarvense_TE/code/conda_envs/meme_suite/
mkdir -p ./MEME_SUITE/
cd MEME_SUITE/

meme ../RLC.FAM.7_allcopies.fa -dna -oc ./allcopies_memeout -nostatus -time 14400 -mod zoops -nmotifs 3 -minw 6 -maxw 50 -objfun classic -revcomp -markov_order 0 
