
### SECOND FILTERING STEP

After the default SPL filtering steps, 
I also run a custom filtering step not present in the Baduel protocol to remove all those insertions which
fall outside the accessibilty  short reads regions  evaluation as thosae are regions with abnormal short read mapping 
and therefore are a source of false positives. 

`bash acc_spl_filter.sh`

This script interacts over every bed file of the final call folder of the SPL raw filtering. <br>
It  does two things:

1. It removes the DHH type TE calls. 
2. It removes any TE call outside of the accesibility regions file.  <br>

For the second, it swtiches back from the Chr_ notation splitreader uses to  the one
the reference genome has (Scaffold_).  <br>

It also fixes the problem of files ending up with an extra tab at the end and also with a extra space at the header row.

** Output files in folder specified in acc_spl_filter.sh output folder** 

After this I split files by individuals using the split_SPL_results.sh  script. 
