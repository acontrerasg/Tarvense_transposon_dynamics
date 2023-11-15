### Aim: Age Determination of singletone TIPs (similar Baduel et al. 2021)
### Author: Dario Galanti Oct 2022
### Input: TIPs matrix. Chr Start End Sample1 Sample2 ...
### Run: sbatch --partition test --cpus-per-task 4 --mem 10G --time 20:00:00 --wrap "bash AgeDeterm_singleTIPs_Epi.sh"
### Dependencies: plink2

### Process:
# NB: In Baduel et al. 2021, they simply assumed all singletones were young and used them for GWAS.
# NB: Also their method for TIPs age determination does not work with singletones, so we use a sligltly different one
# 1) Iterate through singletone TIPs
# 2) Count private singletone SNPs in the flanks.
# 5) Calculate TIP age: Divide n° mutations by the mutation rate expected on 70kb with 1 generations per year.

## Mutation rate = 7*10^(-9) Mutations/(bp*generation) (Ossowski et al. 2010)
## Thlaspi mutation rate = 7*10^(-9) Mutations/(bp*generation) * 70000 bp * 1 generations/year   =   0.00049 Mutations/year

## Define input
work=/scr/episan/RP07
wDir=${work}/TIPs/insertions_withUS
TIPs_mx=${wDir}/ALL_TIPS_distinct.matrix
vcf=${work}/Ta_SNPs/Ta_v4v5/withUSspls_sw_filt/Ta_v5_SNPs_10NA_bial_2MAC_withUSsplsSW.vcf.gz	# We use all SNPs with MAC>=2
fin_base=${work}/Ta_SNPs/Ta_v4v5/withUSspls_sw_filt/Ta_v5_SNPs_10NA_bial_2MAC_withUSsplsSW

## Define output
outdir=${wDir}/Age_Determin
mkdir -p ${outdir}
TIPs_single=${wDir}/TIPS_distinct_single.matrix
tmp_pwdist=${outdir}/PwDist_TIPs_single_tmp
TIPs_age=${outdir}/TIPs_distinct_single_age_freq.bed

## Define parameters
cis_reg=35000

## Define tools
plink=~/conda/plink2/bin/plink2


### 0) Recode vcf in plink
#zcat $fin | awk 'OFS="\t"{if(!/\#/){$3=$1":"$2;print} else {gsub("_","-");print}}' > ${outdir}/temp_ready.vcf
#zcat $fin | awk 'OFS="\t"{if(!/\#/){$3=$1":"$2;print} else {print}}' > ${outdir}/temp_ready.vcf
#$plink --vcf ${outdir}/temp_ready.vcf --allow-extra-chr --make-pgen --out $fin_base

### 1) Extract singletone TIPs
cat $TIPs_mx | awk '{OFS="\t";if(NR==1){print $0} else {for(i=4;i<=NF;i++){s+=$i};if(s==1){print $0};s=0}}' > $TIPs_single

### 2) Iterate through singletones TIPs matrix
spls_str=$(head -1 $TIPs_single | tr "\t" "," | tr "_" "-")
samples=($(head -1 $TIPs_single | cut -f4- | tr "\t" " " | tr "_" "-"))
spl_num=${#samples[@]}
echo -e chr"\t"st"\t"end"\t"n_carriers"\t"Areas"\t"max_SNPcount"\t"years > $TIPs_age

tail -n+2 $TIPs_single | while read line;
do
 chr=$(echo $line | cut -d" " -f1)
 st=$(echo $line | cut -d" " -f2)
 end=$(echo $line | cut -d" " -f3)
 cis_st=$(echo $line | awk -v cis_r=$cis_reg '{mid=int(($2+$3)/2);st=(mid-cis_r);cis_st=(st > 0 ? st : 0);print cis_st}')
 cis_end=$(echo $line | awk -v cis_r=$cis_reg '{mid=int(($2+$3)/2);cis_end=(mid+cis_r);print cis_end}')
 TIP=${chr}:${st}-${end}
 carriers=$(echo $line | awk -v spls=$spls_str 'BEGIN{split(spls,s,",")} {for(i=4;i<=NF;i++){if($i==1){printf "%s ", s[i]}}}' | tr -d " ")
 carriers_arr=($carriers)
 n_carriers=$(echo ${#carriers_arr[@]})
 Areas_arr=($(echo $carriers | awk '{OFS="\n";for(i=1;i<=NF;i++){gsub(/FR|DE|NL|SE/,"EU",$i);$i=substr($i,4,2)};print}' | sort -u))
 Areas_str=$(echo ${Areas_arr[@]//} | tr " " ",")
 spls_ordered=$(echo $carriers ${samples[@]/$carriers})
 $plink --bfile $fin_base --allow-extra-chr --chr ${chr} --from-bp ${cis_st} --to-bp ${cis_end} --sample-diff base=${spls_ordered} --out $tmp_pwdist
 report=${tmp_pwdist}.$(echo 0_${carriers}).sdiff
 ## Plink only reports variable positions for each pw comp, so we simply check whether there are >= 279 comparisons for a specific sample and position (with 1/1 in the sample and 0/0 in all others)
 ## Stringent option
 #max_SNPcount=$(tail -n+2 ${report} | awk 'OFS="\t"{print $1,($2-1),$2,$10,$11}' | bedtools merge -i stdin -d -1 -c 4,5,4 -o distinct,distinct,count | awk 'BEGIN{c=0} {if($4=="1/1" && $5=="0/0" && $6>=279){c++}} END {print c}')
 ## Reasonable option (allows for a few samples with missing values, but still not for other samples with the variant)
 max_SNPcount=$(tail -n+2 ${report} | awk 'OFS="\t"{print $1,($2-1),$2,$10,$11}' | bedtools merge -i stdin -d -1 -c 4,5,4 -o distinct,distinct,count | awk 'BEGIN{c=0} {if($4=="1/1" && $5=="0/0" && $6>275){c++}} END {print c}')
 age=$(echo ${max_SNPcount} | awk '{years=int($0/0.00049);print years}')
 echo -e ${chr}"\t"${st}"\t"${end}"\t"${n_carriers}"\t"${Areas_str}"\t"${max_SNPcount}"\t"${age} >> $TIPs_age
 rm ${tmp_pwdist}*
done





