#!/bin/bash                                                                                                        blast_ALESIA.sh                                                                                                                     
blast_db="/ebio/abt6_projects9/abt6_databases/db/blast_nt_JUN2022/nt"
consensi="/ebio/abt6_projects8/Tarvense_TE/data/Tarvense_ALESIA_TE/RLC.FAM.7_consensi.fa"


# Version: blast 2.13.0
blastn -task blastn \
       -db ${blast_db} \
       -query ${consensi} \
       -out ALESIA_nr.blast.out \
       -num_threads  16 \
       -outfmt '7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs stitle salltitles'

# Retrieve species  that got a good match in their genome:

grep -v "Thlaspi arvense" ALESIA_nr.blast.out  | \
grep -v "#" | \
awk '{if($13 >= 80 && $3 >=80) print $0}' | \
cut -f14 | cut -f1 -d ',' | cut -f1,2 -d ' ' | \
sort | uniq > Alesia_species_matches.txt

# Retrieve matches that likely represent full copy TEs 
grep -v "#"  ALESIA_nr.blast.out |
awk '{if($13 >= 80 && $3 >= 80 && $4 > 3000) print $0}' > ALESIA_nr.blast.filtered.out

# Retrieve sequene matches:
cat  ALESIA_nr.blast.filtered.out | awk '{if($10 < $9) print $2"\t"$10"-"$9 ; else print $2"\t"$9"-"$10}' > batch_list
blastdbcmd -db ${blast_db} -entry_batch  batch_list -out  ALESIA_nr.blast.filtered.out.fa -outfmt %f -long_seqids


# Modify names so they are easier to parse (and  do not collide with raxml tree syntax):
## NCBI GI, underscore start, underscore end, underscore species.
cat ALESIA_nr.blast.filtered.out.fa | 
sed 's/ genome.*//;s/ HDEM//;s/ haplo.*//;s/|emb.*|//;s/ /_/g;s/:/_/'  > ALESIA_nr.blast.filtered.out.fixnames.fa
# remove duplicated sequences ( there is only one ). 
seqkit rmdup -s ALESIA_nr.blast.filtered.out.fixnames.fa > ALESIA_nr.blast.filtered.out.fixnames.nodups.fa 
# Make aligment:
mafft --adjustdirection --maxiterate 1000 --globalpair --thread 16 ALESIA_nr.blast.filtered.out.fixnames.nodups.fa > ALESIA_nr.blast.filtered.out.aln
# Make tree:
raxml-ng -all -msa ALESIA_nr.blast.filtered.out.aln  -model JTT+G -model JTT+G  --bs-trees 100 --threads 16
# Visu and color tree with figtree. 
