library(ape)
library(reshape2)
library(stringr)
library(plyr)

#This scripts calculates the LTR insertion age based on the structural information form LTR harvest with  different models:
#K2p
#raw
#tn93
#Then it  also retrives the LTR age stimation produced by ltr harvest and  compares it with the 3 newly created metrics.
#Latsly it saves this dataframe as text file in order to relete it to the original annotation LTRs with the "LTR_age2anno.sh" script. 

#path to folder where aligned_ltrs are:
path="/ebio/abt6_projects8/Tarvense_TE/data/Tarvense_TE_pops/Thlaspi_genome_2021/annotation/TEanno_fams/get_fastas_fams/LTR_structural/"
filenames=list.files(paste0(path,'aligned_ltrs'), pattern='.fa')
ages=data.frame(tename=stringr::str_split_fixed(filenames, '\\.', 3)[,1])
ages$k2p=NA
ages$raw=NA
ages$tn93=NA

for (i in 1:length(filenames)){
  ltr=read.FASTA(paste0(path,'aligned_ltrs/',filenames[i]))
  if( !is.null(ltr) ){
    d=dist.dna(ltr, model='K80', gamma=F)
    ages$k2p[i]=d
    raw=dist.dna(ltr, model='raw', gamma=F)
    ages$raw[i]=raw
    tn93=dist.dna(ltr, model='tn93', gamma=F)
    ages$tn93[i]=tn93
  }
}

#Get LTR age estimation from LTR harvest.
te <- read.delim(paste0(path,"LTR_harvest_age_predict.tsv"))
colnames(te) <-  c("LTR_ID","ltr_similarity")
ages$lk2p=100-(ages$k2p*100)
ages$lraw=100-(ages$raw*100)
ages$ltn93=100-(ages$tn93*100)
ages$gtsim=as.numeric(as.character(mapvalues(ages$tename, from=te$LTR_ID, to=te$ltr_similarity, warn_missing=F)))


#plot the differences in value  between ltr harvest  and the other age estimators. 
png(paste0(path,'ltrharvest_age_difference.%03d.png'))
plot(ages$lraw~ages$gtsim, xlab='LTR harvest age', ylab='RAW mafft age', pch=19, type='p')
plot(ages$lk2p~ages$gtsim, xlab='LTR harvest age', ylab='K2P mafft age', pch=19, type='p')
plot(ages$ltn93~ages$gtsim, xlab='LTR harvest age', ylab='TN93 mafft age', pch=19, type='p')
dev.off()

#Write final file: 
write.table(ages, paste0(path,'LTRs_tarvense_ages.txt'), quote=F, col.names=F, row.names=T, sep='\t')
