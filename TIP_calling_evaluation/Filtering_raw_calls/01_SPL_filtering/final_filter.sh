#!/bin/bash

# # Filter presence variant calls based on negative coverage ratios (NC-filter)
# # This is a wrapper around the perl script: Filter_negative_calls_splitreader.pl .

# Variables

Fncspl="/ebio/abt6_projects8/Tarvense_TE/code/SPLITREADER/SPLITREADER_runs/Filter_negative_calls_splitreader.pl"
depth=3 #Based  on the previous step at DP_filtering_wrapper_run.sh
filtname="/ebio/abt6_projects8/Tarvense_TE/data/Tarvense_TE_pops/Splitreader/RESULTS/Results_filtering/baduel_filtering/DP-filter/Tarvense_filtered_DP3.bed"
# #  Filtname contains the path to the first pass filter (positive coverage)  (Tarvense_filtered_DP3.bed in my case. )
path_to_RC="/ebio/abt6_projects8/Tarvense_TE/data/Tarvense_TE_pops/Splitreader/RESULTS/Results_filtering/baduel_filtering/NC-filter/all_RC/"
ind_thresh=4  # Minimun number of individualss to be covered. In my case I set it up to 4 as I  have 7 armenian individuals.
out_dir="/ebio/abt6_projects8/Tarvense_TE/data/Tarvense_TE_pops/Splitreader/RESULTS/Results_filtering/baduel_filtering/final_call/"

perl ${Fncspl} ${ind_thresh} ${depth} ${path_to_RC} ${out_dir} ${filtname} /ebio/abt6_projects8/Tarvense_TE/data/Tarvense_TE_pops/Splitreader/RESULTS/Results_filtering/baduel_filtering/all_raw_SPL_calls/*-insertion-sites.bed > final_filter.log
