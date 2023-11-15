### Aim: Age Determination of TIPs (see Baduel et al. 2021)
### Author: Dario Galanti Oct 2022
### Input: TIPs matrix. Chr Start End Sample1 Sample2 ...
### Run: sbatch --partition test --cpus-per-task 4 --mem 10G --time 20:00:00 --wrap "bash AgeDeterm_TIPs_Epi.sh"
### Dependencies: plink2

### Process:
# 1) Iterate through TIPs
# 2) For each TIP subset the plink input to the region of interest and samples with the insertions
# 3) Calculate the pairwise distance count (number of SNPs) between all samples with the insertion
# 4) Extract the highest pairwise distance between any two TIP carriers and
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
TIPs_MAC2=${wDir}/TIPS_distinct_2MAC.matrix
tmp_pwdist=${outdir}/PwDist_TIPs_2MAC_tmp
TIPs_age=${outdir}/TIPs_distinct_2MAC_age_freq.bed

## Define parameters
cis_reg=35000

## Define tools
plink=~/conda/plink2/bin/plink2

### 0) OPTIONAL: Recode vcf in plink
#zcat $fin | awk 'OFS="\t"{if(!/\#/){$3=$1":"$2;print} else {gsub("_","-");print}}' > ${outdir}/temp_ready.vcf
#zcat $fin | awk 'OFS="\t"{if(!/\#/){$3=$1":"$2;print} else {print}}' > ${outdir}/temp_ready.vcf
#$plink --vcf ${outdir}/temp_ready.vcf --allow-extra-chr --make-pgen --out $fin_base

### 1) Extract TIPs with MAC>=2
cat $TIPs_mx | awk '{OFS="\t";if(NR==1){print $0} else {for(i=4;i<=NF;i++){s+=$i};if(s>1){print $0};s=0}}' > $TIPs_MAC2

### 2) Iterate through TIPs matrix without singletones
spls_str=$(head -1 $TIPs_MAC2 | tr "\t" "," | tr "_" "-")
echo -e chr"\t"st"\t"end"\t"n_carriers"\t"Areas"\t"max_SNPcount"\t"years > $TIPs_age

tail -n+2 $TIPs_MAC2 | while read line;
do
 chr=$(echo $line | cut -d" " -f1)
 st=$(echo $line | cut -d" " -f2)
 end=$(echo $line | cut -d" " -f3)
 cis_st=$(echo $line | awk -v cis_r=$cis_reg '{mid=int(($2+$3)/2);st=(mid-cis_r);cis_st=(st > 0 ? st : 0);print cis_st}')
 cis_end=$(echo $line | awk -v cis_r=$cis_reg '{mid=int(($2+$3)/2);cis_end=(mid+cis_r);print cis_end}')
 TIP=${chr}:${st}-${end}
 carriers=$(echo $line | awk -v spls=$spls_str 'BEGIN{split(spls,s,",")} {for(i=4;i<=NF;i++){if($i==1){printf "%s ", s[i]}}}')
 carriers_arr=($carriers)
 n_carriers=$(echo ${#carriers_arr[@]})
 Areas_arr=($(echo $carriers | awk '{OFS="\n";for(i=1;i<=NF;i++){gsub(/FR|DE|NL|SE/,"EU",$i);$i=substr($i,4,2)};print}' | sort -u))
 Areas_str=$(echo ${Areas_arr[@]//} | tr " " ",")
 $plink --bfile $fin_base --allow-extra-chr --chr ${chr} --from-bp ${cis_st} --to-bp ${cis_end} --sample-diff ids=${carriers} --out $tmp_pwdist
 max_SNPcount=$(tail -n+2 ${tmp_pwdist}.sdiff.summary | awk 'BEGIN{max=0} {if($6>max){max=$6}} END{print max}')
 age=$(echo ${max_SNPcount} | awk '{years=int($0/0.00049);print years}')
 echo -e ${chr}"\t"${st}"\t"${end}"\t"${n_carriers}"\t"${Areas_str}"\t"${max_SNPcount}"\t"${age} >> $TIPs_age
done

