# Plot DNA methylation around insertions

This directory contains scripts to plot DNA methylation in the flanking regions of TE insertions in samples with vs without the insertion.
It takes methylation information from unionbed files that can be generated with code available at [this](https://github.com/Dario-Galanti/WGBS_downstream/tree/main) repository. <br/>

[met_TIP_flanks.sh](https://github.com/acontrerasg/Tarvense_transposon_dynamics/blob/main/TIP_methylation/met_TIP_flanks.sh) <br/>
This bash script uses a bed file of TIPs and a methylation unionbed file (samples as columns and positions as rows, see [here](https://github.com/Dario-Galanti/WGBS_downstream/tree/main) for more explanation) to calculate methylation in the flanking regions of TE insertions, in samples with and without the insertion. More information on the input files are in the beginning of the script.
It produces a bed file of 50 bp bins covering 2kb up and downstream insertions, with average methylation in samples with and without the insertion. <br/>
Example output file: <br/>
st      end     bp      Flk     Samples_with_insertion  Samples_without_insertion <br/>
-1950   -1900   -1925   flk1    15.69   16.33 <br/>
-1900   -1850   -1875   flk1    12.99   13.29

[TIPflks_met.py](https://github.com/acontrerasg/Tarvense_transposon_dynamics/blob/main/TIP_methylation/TIPflks_met.py) <br/>
This python script is required to run the previous one. It calculates methylation in samples with and without the insertion.

[plot_TIPflank_met.R](https://github.com/acontrerasg/Tarvense_transposon_dynamics/blob/main/TIP_methylation/plot_TIPflank_met.R) <br/>
This R script plots methylation of bins covering the flanks of TE insertions, producing output similar to the one below:

![image](https://github.com/acontrerasg/Tarvense_transposon_dynamics/assets/58292612/fa4b5a06-e675-4891-b6f2-d4cf5e4c7639)

