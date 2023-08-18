## Author: Dario Galanti 2022
## Aim: Calculate number of TIPs overlapping with different genomic features and then correct it for genome space covered by each genomic feature.
## Input 1: List of TIPs of different Superfamilies or different geographic locations (1 bed file per group)
## Input 2: Bed files of genomic features

## Location
work=/scr/episan/RP07/TIPs
indir=${work}/insertions_withUS/individuals_TIPs
ref_index=/scr/episan/RP07/Thlaspi_genome/Ta_v5_genome/final.fasta.fai
#country_dir=${work}/insertions_withUS/TIPs_byCountry
SupFam_dir=${work}/insertions_withUS/TIPs_bySupFam_merged
Fam_dir=${work}/insertions_withUS/TIPs_byFam_merged
fout=${work}/insertions_withUS/Features_intersect_SupFam_TIPs.txt

## Define feature bed files
featureDir=/scr/episan/RP07/region_meth/Features_meth_v5/feature_beds
genes=${featureDir}/Ta_v5_genes.bed					# NB: Feature filenames must contain feature name just before file extention
CDS=${featureDir}/Ta_v5_CDS.bed					# NB: Feature filenames must contain feature name just before file extention
introns=${featureDir}/Ta_v5_introns.bed				# NB: Feature filenames must contain feature name just before file extention
TEs=${featureDir}/Ta_v5_TEs.bed						# NB: Feature filenames must contain feature name just before file extention
promoters=${featureDir}/Ta_v5_promoters.bed			# NB: Feature filenames must contain feature name just before file extention
intergenic=${featureDir}/Ta_v5_intergenic.bed		# NB: Feature filenames must contain feature name just before file extention

### TIPs BY SuperFamily
#All_TIPs=All_TIPS_cov_corrected.tsv
#SupFamArr=($(tail -n+2 $All_TIPs | cut -f14 | sort | uniq))
#FamArr=($(tail -n+2 $All_TIPs | cut -f6 | sort | uniq))
#ClassArr=($(tail -n+2 $All_TIPs | cut -f14 | cut -c1-2 | sort | uniq))
#TEgroupArr=("${SupFamArr[@]}" "${ClassArr[@]}" "R")
#
#for group in ${TEgroupArr[@]}
#do
# group_fout1=${SupFam_dir}/TA_${group}_insertions.bed
# group_fout2=${SupFam_dir}/Genes_2kbFlk_${group}_TIPs.txt
# tail -n+2 $All_TIPs | grep -P '\t'$group | awk 'OFS="\t"{print $3,$4,$5,$6,$1,$2}' | sort -k1,1V -k2,2n > $group_fout1
# bedtools merge -i $group_fout1 | bedtools slop -i stdin -g $ref_index -b 2000 | bedtools intersect -a $Genes -b stdin -u | cut -f5 > $group_fout2
#done

# MAKE FILE ARRAYS
reg_arr=( global $genes $CDS $introns $TEs $promoters $intergenic )
#Group_arr=($(ls ${SupFam_dir}/TA_*.bed))
Group_arr=($(ls ${SupFam_dir}/TA_*.bed))

# INTERSECT ALL CONTEXT AND REGION FILES
echo -e Feature"\t"Group"\t"Num_TIPs"\t"Feature_bp"\t"insertions_x_Mb > $fout
# First calculate for whole genome (no intersection)
for reg in ${reg_arr[@]};
do
 if [ $reg == "global" ];
 then
  feature=global
  bp=$(awk '{tot+=$2} END {print tot}' $ref_index)
 else
  feature=$(basename $reg .bed | cut -d"_" -f3)
  bp=$(awk '{len+=($3-$2)} END {print len}' $reg)
 fi
 for SF in ${Group_arr[@]};
 do
  group=$(basename $SF | cut -d"_" -f2)
  if [ $reg == "global" ];
  then
   num=$(bedtools merge -i $SF -d -1 | wc -l)
  else
   num=$(bedtools merge -i $SF -d -1 | bedtools intersect -a stdin -b $reg -u | wc -l)
  fi
  num_xBP=$(echo $num | awk -v bp=$bp '{printf "%.2f", ($0/(bp/1000000))}')
  echo -e $feature"\t"$group"\t"$num"\t"$bp"\t"$num_xBP | tr " " "\t" >> $fout
 done
done

## Add Insertion rate and Insertion density scaled to 1 insertion/MB globally
## NB: We need "global" to be the first line of each family, so we use a trick and temporarily change it to "aglobal"
tmp_fout=${work}/insertions_withUS/TMP_Features_intersect_SupFam_TIPs.txt
echo -e $(head -1 $fout)"\t"Insertion_rate"\t"Insertion_x_MB_global_scaled | tr " " "\t" > $tmp_fout
tail -n+2 $fout | sed 's/global/aglobal/g' | sort -k2,2V -k1,1 | sed 's/aglobal/global/g' | awk 'OFS="\t"{if($1=="global"){TIPs=$3;TIPdens=$5;print $0,1,1}else{rate=($3/TIPs);density_rate=($5/TIPdens);print $0,rate,density_rate}}' >> $tmp_fout
mv $tmp_fout $fout


