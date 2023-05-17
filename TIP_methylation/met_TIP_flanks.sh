### Aim: Plot methylation in the flanking regions of insertions
### Author: Dario Galanti Sep 2022
### Input 1): Bed file of TIPs, properly sorted (sort -k1,1V -k2,2n). Example line:
###		Scaffold_1   3041404 3041419 RLG.FAM.15   North_Germany,South_Sweden   8   TA_DE_16_01_F1_HC0_M1_1,TA_DE_16_09_F1_HC0_M1_1,TA_SE_04_05_F1_HC0_M1_1
### Input 2): unionbed file of methylation (headers: chrom start end sample1 sample2 ...)
### Input 3): List of samples with methylation information, no header (to further subset the unionbed file if necessary)
### Input 4): TEs bed file
### Input 5): genome index file
### Run: bash met_TIP_flanks.sh fin.bed unionbed.bed
### Run: bash met_TIP_flanks.sh ../TIPs_bySupFam/TA_R_insertions.bed unionbed_v5_3cov_25NAs/CpG_unionbed_v5_25NAs.bed
### Run Epi: sbatch --partition test --cpus-per-task 2 --mem 10G --time 04:00:00 --wrap "bash met_TIP_flanks.sh ../TIPs_bySupFam/TA_R_insertions.bed unionbed_v5_3cov_25NAs/CpG_unionbed_v5_25NAs.bed"
### Process:
### 	1) Select TIPs present in lines with methylation data (not singletones).
###		2) Remove TIPs closer than 2kb to each other
###		3) Remove TIPs closer than 1kb to annotated TEs
###		5) Extract flanks of each TIP
###		6) Intersect them with unionbed methylation file
###		7) Average methylation of samples with and without insertion/deletion. Calculate position relative to flank (flank_pos) (python script)
###		8) Sort according to flank_pos, average methylation of flanks in bins
###		9) Combine flank1 and 2 files
###		10) Plot: NB the R script is calculating moving mean for every 3 bins to smooth the curve.

## Define input
TEs=/scr/episan/RP07/region_meth/Features_meth_v5/feature_beds/Ta_v5_TEs.bed
ref_index=/scr/episan/RP03/thlaspi/REVISION/FINAL/final.fasta.fai
METspls=METspls_NCmore1.2.txt
Rscript=~/conda/R4/bin/Rscript
fin=$1					# /scr/episan/RP07/TIPs/insertions_withUS/TIPs_bySupFam/TA_R_insertions.bed
unionbed=$2				# /scr/episan/RP07/TIPs/insertions_withUS/met_TIP_flanks/unionbed_v5_3cov_25NAs/CpG_unionbed_v5_25NAs.bed
group=$(basename $fin | cut -d"_" -f2)
cont=$(basename $unionbed | grep -o 'C[pH][HG]')	# Extract context

## Define output
outdir=${group}
mkdir -p $outdir
flks_met=${outdir}/TIP_${group}_m${cont}_flks_50bin_2kbself_1kbfromTEs.txt
plot=${outdir}/TIP_${group}_m${cont}_flks_50bin_2kbself_1kbfromTEs_smooth.pdf

## Define intermediate files
file1=${outdir}/TIPs_${cont}_inMETspls.bed
file2=${outdir}/TIPs_${cont}_lonely.bed
file3=${outdir}/TIPs_${cont}_lonely_2kbfromTEs.bed
file4=${outdir}/TIPs_${cont}_lonely_2kbfromTEs_ID.bed
tmp2=${outdir}/tmp2_${group}_${cont}.txt
tmp3=${outdir}/tmp3_${group}_${cont}.txt
flanked_TIPs=${outdir}/TA_${group}_TIPs_lonely_noTEs_ID_flk.bed
union_flanks=${outdir}/${group}_${cont}_unionTIPflanks.bed
tmp_flanks_met=${outdir}/${group}_${cont}_TIPflanks_met.txt
flk1_met=${outdir}/Flk1_${group}_${cont}_50binned_met.txt
flk2_met=${outdir}/Flk2_${group}_${cont}_50binned_met.txt

echo Input TIPs: $fin
echo Input unionbed: $unionbed


## 1) Select TIPs present in samples with methylation and remove singletones or insertions present in all samples
fgrep -f $METspls $fin | awk 'OFS="\t"{if($6>1 && $6<200){print}}' > $file1
echo $(cat $file1 | wc -l) TIPs of group ${group} present in samples with methylation info

## 2) Remove TIPs closer than 2kb to each other
self_dist=1000		# TIPs closer than 2*$self_dist to each other, or closer than $self_dist to scaffolds start are excluded.
awk -v sf=$self_dist 'OFS="\t"{st=($2-sf);end=($3+sf);if(st>=0){print $1,st,end,$4,$5,$6,$7}}' $file1 | bedtools merge -i stdin -c 4,5,6,7 -o collapse,collapse,collapse,collapse | awk -v sf=$self_dist 'OFS="\t"{split($6,fam,",");if(length(fam)==1){st=($2+sf);end=($3-sf);print $1,st,end,$4,$5,$6,$7}}' > $file2
echo $(cat $file2 | wc -l) TIPs of group ${group} after removing the ones closer than 2kb to each other and closer than 1kb to scaffold starts

## 3) Remove TIPs in TEs
#bedtools intersect -a $file2 -b $TEs -v > $file3													# remove TIPs in TEs
bedtools slop -i $TEs -g $ref_index -b 1000 | bedtools intersect -a $file2 -b stdin -v > $file3		# remove TIPs closer than 1kb to TEs
echo $(cat $file3 | wc -l) TIPs of group ${group} after removing the ones closer than 1kb to TEs.

## 4) Add TIP ID
sort -k4,4V $file3 | awk -v fam="start" 'OFS="\t"{if($4==fam){n+=1}else{n=1};ind=$4"."n;print $1,$2,$3,ind,$4,$5,$6,$7;fam=$4}' | sort -k1,1V -k2,2n > $file4

## 5) Flank the TIPs (and add flk1/2)
bedtools flank -i $file4 -g $ref_index -b 2000 | awk 'OFS="\t"{if(NR %2 == 0){flk="flk2"}else{flk="flk1"};print $1,$2,$3,$4,$5,flk,$7,$8}' > $flanked_TIPs

## 6) Intersect with unionbed
## NB: Do not use -u option in bedtools intersect, because some positions might overlap with both the flk1 of one TIP and flk2 of another
tail -n+2 $unionbed | cut -f1-3 | bedtools intersect -a stdin -b $flanked_TIPs -wb > $tmp2
tail -n+2 $unionbed | bedtools intersect -a stdin -b $flanked_TIPs | cut -f4- > $tmp3
echo -e chr"\t"st"\t"end"\t"flk_chr"\t"flk_st"\t"flk_end"\t"ID"\t"fam"\t"flk"\t"n_spls_withDel"\t"spls_withDel"\t"$(head -1 $unionbed | cut -f4-) | tr " " "\t" > $union_flanks
paste $tmp2 $tmp3 >> $union_flanks

## 7) Average met of samples with and without insertion/deletion and calculate position relative to flank
python3 TIPflks_met.py $union_flanks $tmp_flanks_met
echo $(tail -n+2 $tmp_flanks_met | cut -f3 | uniq | wc -l) TIPs of group ${group} present in more than 1 sample. These will be used for plotting methylation of their flanks

## 8) Bin flank1 and flank2 separately and calculate average methylation of bins
## Possibly buggy, if there are bins without any positions. But very unlikely using many TEs.
bin=50
bin_end=$(expr -2000 + $bin)
tail -n+2 $tmp_flanks_met | grep flk1 | sort -k1,1n | awk -v bin=$bin -v bin_end=$bin_end '{if($1>bin_end){met_w=(met_with/c);met_n=(met_no/c); printf   "%s\t%s\t%s\t%s\t%.2f\t%.2f\n", (bin_end-bin),bin_end,(bin_end-int((bin/2)+0.999)),"flk1",met_w,met_n; bin_end=(bin_end+bin);met_with=$5;met_no=$6;c=1} else {met_with+=$5;met_no+=$6;c++}}' > $flk1_met

bin_end=$(expr 0 + $bin)
tail -n+2 $tmp_flanks_met | grep flk2 | sort -k1,1n | awk -v bin=$bin -v bin_end=$bin_end '{if($1>bin_end){met_w=(met_with/c);met_n=(met_no/c); printf   "%s\t%s\t%s\t%s\t%.2f\t%.2f\n", (bin_end-bin),bin_end,(bin_end-int((bin/2)+0.999)),"flk2",met_w,met_n; bin_end=(bin_end+bin);met_with=$5;met_no=$6;c=1} else {met_with+=$5;met_no+=$6;c++}}' > $flk2_met

## 9) Combine Flk1 and Flk2
echo -e st"\t"end"\t"bp"\t"Flk"\t"Samples_with_insertion"\t"Samples_without_insertion > $flks_met
awk FNR!=1 $flk1_met $flk2_met >> $flks_met

## 10) Plot with R
$Rscript --vanilla plot_TIPflank_met.R $flks_met $plot

## 11) Cleanup
rm $tmp2 $tmp3 $file2 $file3

