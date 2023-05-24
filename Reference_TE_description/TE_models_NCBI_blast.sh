#!/bin/bash 
# Here we first run BLAST vs the NCBI nr to explore how conserved are our TE families:

blast_db="/ebio/abt6_projects9/abt6_databases/db/blast_nt_JUN2022/nt"
TElib="/ebio/abt6_projects8/Tarvense_TE/data/Tarvense_TE_pops/Thlaspi_genome_2021/thlaspi.TE.fa"
Taxonomy_folder="/ebio/abt6_projects9/abt6_databases/db/NCBI_taxonomy/24Jan2023/"

# Main blast
blastn -task blastn -db ${blast_db} \
	   -query ${TElib}  \
	   -out ./ALL_TElib_blastout.tsv \
	   -num_threads 12 \
	   -outfmt '7 qseqid qlen sgi sseqid slen pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs stitle salltitles staxids'   
# Compresss
gzip ./ALL_TElib_blastout.tsv
# FILTER OUTPUT, only retain those hits that are NOT thlaspi arvense and are 80% similarity over 80% of the TEmodel length, and 4Kb long minimun. 
zcat ./ALL_TElib_blastout.tsv.gz |
grep -v "^#" | 
awk 'BEGIN{FS=OFS="\t"}{perc_cov = ($7-$9)/$2 ; if( $6 >= 80  &&  perc_cov  > 0.8 && $2 > 4000) print $0,perc_cov}'  > ALL_TElib_blastout.filtered.tsv
taxonkit  lineage -d ';' -n -i 19 ALL_TElib_blastout.filtered.tsv > ALL_TElib_blastout.filtered.taxonomy.tsv
