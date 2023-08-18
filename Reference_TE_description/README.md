Here we deposited the scripts used to analize the TEs present in the reference. The TE library used to retrieve the models and the TE annotation file can be found in: <br>
"https://onlinelibrary.wiley.com/doi/full/10.1111/pbi.13775" Chromosome-level Thlaspi arvense genome provides new tools for translational research and for a newly domesticated cash cover crop of the cooler climates. Nunn et al 2022. <br>


There are several folders and scripts:
- `TE_autonomous_eval.sh`  Uses the dante pipeline: https://github.com/kavonrtep/dante over the TE models of the TE library to classify each famly as autonomous or non-autonomous based on the protein domains present in it. 
- `TE_models_NCBI_blast.sh` is  the script used to check the level of conservation of long TEmodels in *T. arvense* (>4kb) in the NCBI database as per June of 2022. 
-  The folder TE_phylo has two scripts. `run_dante_thlaspi_consensi.sh` runs the dante pipeline over the TE models of the TE library  to classify each TE family as different clades. <br> The second script, retrieves the RT protein domain of these TE models to construct two phylogenetic trees, one for Ty1 and other for Ty3 elements.
-  The folder TE_intacts has the scripts necessary to install and run LTRpred :https://github.com/HajkD/LTRpred over the reference genome to retrieve young, intact TE elements. Then uwe used a custom script `intersect_LTRpred_anno.sh` to associate these LTRs to known TEs of the annotation.
-   We used the LTR age stimates calculated by LTRpred.
