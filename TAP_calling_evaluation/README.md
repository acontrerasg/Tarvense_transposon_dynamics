This directory contains the scripts used to run the custom pipeline popTAP to call Transposon abascences polymorphisms in a collection of short read libraries. 
The main pipeline can be found here:
https://github.com/acontrerasg/popTAP

After running the pipeline for all our samples, using the run_CNV_cluster wrapper scripts, we further filtered these putative TAPs based on mean coverage and pvalue with the filter_output.sh script. These are our real absences. 
However, we wanted to make sure called TAPs were not the product of biggger rearragments, therefore, we performed a second filtering round of these TAPs  based on flanking coverage. 


**WORKFLOW**

1. Run `run_CNV_cluster.sh` This calls the popTAP tool with our samples.
2. Run `run_filter_CNV_cluster.sh` This filters putative calls called by popTAP based on p-value of 0.05 and coverage lower than 1/10 of the median coverage of the sample.  Merge all individual TAPs into a single file.
3. Run `check_coverage_TAP_flanks.sh` to evaluate flakning coverage of the individual TAPs. This step needs a file generated by msodepth (you cna use the one generated by popTAP) in order to retireve the sample coverage. Ex: `mosdepth -n $sample_name $sample.bam`
4. We futhered filter these TAPs based on the flaking coverage,  Any TAP with  less than 5X or  more than 50X mean coverage in either flank (3_prime or 5_prime) were further remove for the analysis. This is the final TAP set. 
