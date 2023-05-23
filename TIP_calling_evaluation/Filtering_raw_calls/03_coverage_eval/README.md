#### Explanation 
This last round of filtering merges together different coverage metrics  the tool uses:
"TIP_cov ref_cov upstream_100cov downstream_100cov".
Adds  to this the local coverage at the 10kb window the insertion was detected plus the 
overall median coverage fo the sample: "local_10kb_median_cov  sample_median_cov".

To calculate the local and mean coverage we used raw reads of the sample and map them using
bwa mem with defaults over the genome. Then we used mosdepth to make the coverage calculations.  

This way  we cna make finer adjustments of wether an insertion is realiable or not. 

Care of use lifoff_v4_v5.sh over all bed files.


**Workflow**

`bash run_cov_eval.sh`
`mkdir -p files_individual`
`mv *_cov_TIPs_eval.tsv  files_individual`
###### Filter higher than 100 coverage at upstream aor downstream TIPS and filter out single headers

`cat files_individual/*tsv | head -n 1  > header`

`for files in  files_individual/*tsv ; do tail -n +2 $files ; done | 
 awk  'BEGIN{OFS=FS="\t"}{if ( $7 < 100 && $8 < 100) print $0 }' | 
cat  header  -  > All_TIPS_coverage_corrected.tsv`



Your final file is : <br>
**All_TIPS_coverage_corrected.tsv**
