#!/bin/bash

#############################################
#                                           #
#               SPLITREADER                 #
#  calculate negative coverage by genome    #
#           Baduel et al. 2020              #
#                                           #
#############################################

# # Modified by Adrian Contreras Garirdo.

# # This script calculate the negative coverage (RC) for all sites within -insertion-sites.[region].bed
# # in the whole-genome alignment bam file (${in}) and stores the minimum over each interval processed by
# # Process_BAMrc_splitreader.pl for each $popname in output file: ${popname}-insertion.[region].rc .
# # If the read coverage was already calculated over a subset of the sites in -insertion-sites.[region].bed it does not re-analyze them to reduce
# # computational time.
# # The TOTAL_insertion-sites.[region].bed input files can be generated directly from the combined list of putative insertion sites
# # (all the $fam.$subsetname-insertions.filt.DP$depth.bed from the Filter_insertions_splitreader.pl).
# # and reformated directly as a bed-file (-insertion-sites.0.bed) or shifted 100bp upstream (-100 on both start and stop columns
# # in -insertion-sites.100up.bed) or 100bp downstream (-insertion-sites.100down.bed).

##Questions or comments to pbaduel(ar)bio.ens.psl.eu

usage() { echo "$0 Usage:
Needs bam-readcount (https://github.com/genome/bam-readcount) installed.
Installed at conda env: splitreader." && grep " .)\ #" $0; exit 0; }
[ $# -eq 0 ] && usage
while getopts ":hi:g:P:d:O:F:" arg; do
  case $arg in
    i) # Bam file
      in=${OPTARG}
      ;;
    g) # Path to reference genome file.
      genomefile=${OPTARG}
      ;;
    P) # Path to  Process_BAMrc_splitreader.pl (/users/a2e/pbaduel/myScripts/Process_BAMrc_splitreader.pl)
      Process_BAMrc_splitreader=${OPTARG}
      ;;
    d) # number of reads required to call an insertion on 1st pass
      depth=${OPTARG}
      ;;
    O) # # path/to/output
      OutputDir=${OPTARG}
      ;;
    F) # # path/to/Filter_insertions_splitreader/output/0_UP_and_DOWN_insertion-sites.bed
      InputDir=${OPTARG}
     ;;
    h | *) # Display help.
      usage
      exit 0
      ;;
  esac
done

##Define other variables:
popname=$(basename ${in} .bam)

####### BEGIN ########

echo "["$(date +"%y-%m-%d %T")"] Calculating negative coverage over insertions " | tee -a $OutputDir/${popname}_log.txt
echo ${in}

## calculate coverage over insertions
bam-readcount -w 1 -f $genomefile ${in} -l ${InputDir}/Tarvense-insertion-sites.0.bed > $OutputDir/${popname}-insertion.tmp.0.rc
# # process rc counts to keep only min DP over each insertion
perl ${Process_BAMrc_splitreader} Tarvense $depth $popname $InputDir $OutputDir  0
rm $OutputDir/${popname}-insertion.tmp.0.rc

echo "["$(date +"%y-%m-%d %T")"] Finished local BAMrc "  | tee -a $OutputDir/${popname}_${popname}_log.txt

## calculates coverage 100bp up and down

## UP
bam-readcount -w 1 -f $genomefile ${in} -l ${InputDir}/Tarvense-insertion-sites.100up.bed > $OutputDir/${popname}-insertion.tmp.100up.rc
# # process rc counts to keep only min DP over each insertion
perl ${Process_BAMrc_splitreader} Tarvense $depth $popname $InputDir $OutputDir 100up
rm $OutputDir/${popname}-insertion.tmp.100up.rc
echo "["$(date +"%y-%m-%d %T")"] Finished 100up BAMrc "  | tee -a $OutputDir/${popname}_log.txt

## DOWN
bam-readcount -w 1 -f $genomefile ${in} -l ${InputDir}/Tarvense-insertion-sites.100down.bed > $OutputDir/${popname}-insertion.tmp.100down.rc
# # process rc counts to keep only min DP over each insertion
perl ${Process_BAMrc_splitreader} Tarvense $depth $popname $InputDir $OutputDir 100down
rm $OutputDir/${popname}-insertion.tmp.100down.rc
echo "["$(date +"%y-%m-%d %T")"] Finished 100down BAMrc "  | tee -a $OutputDir/${popname}_log.txt

echo "["$(date +"%y-%m-%d %T")"] Finished BAMrc "  | tee -a $OutputDir/${popname}_log.txt
