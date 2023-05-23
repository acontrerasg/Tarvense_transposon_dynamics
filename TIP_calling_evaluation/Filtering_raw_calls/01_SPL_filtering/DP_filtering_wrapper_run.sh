#!/usr/bin/env bash

#This script runs the DP filtering evaulation of Baduel et all in the raw splitreader calls.


main_script="/ebio/abt6_projects8/Tarvense_TE/code/SPLITREADER/SPLITREADER_runs/DP-filter_all.sh"
Filter_insertions_splitreader="/ebio/abt6_projects8/Tarvense_TE/code/SPLITREADER/SPLITREADER_runs/Filter_insertions_splitreader.pl"
TSD_info="/ebio/abt6_projects8/Tarvense_TE/data/Tarvense_TE_pops/Splitreader/thlaspi_TE_TSD_info.txt"
bash $main_script -B /ebio/abt6_projects8/Tarvense_TE/data/Tarvense_TE_pops/Splitreader/RESULTS/Results_filtering/baduel_filtering/all_raw_SPL_calls/ \
		-T  $TSD_info \
		-P $Filter_insertions_splitreader \
		-D 3 \
		-I 10 \
		-N 50
