###########
# bash scripts to set up parallelization - 10 permutations, run on 10 seperate chunks of the data
###########

#!/bin/bash
module load R/3.5.2
Rscript perm.PERMNUMBER1.PERMNUMBER2.R

###

for f in {1..10} ; do cat perm.R | sed -e s/PERMNUMBER1/$f/g > perm.$f.R; 
	for k in {1..10}; do cat perm.$f.R| sed -e s/PERMNUMBER2/$k/g > perm.$f.$k.R;
done; done

for f in {1..10} ; do cat perm.sh | sed -e s/PERMNUMBER1/$f/g > perm.$f.sh; 
	for k in {1..10}; do cat perm.$f.sh| sed -e s/PERMNUMBER2/$k/g > perm.$f.$k.sh;
done; done

rm commands1.sh; touch commands1.sh
for f in {1..10} ; do
	for k in {1..10}; do echo 'sh perm.'${f}'.'${k}'.sh' >> commands1.sh;
done; done

###########
# Rscript
###########

library(data.table)

alt=fread('/Genomics/ayroleslab2/shared/longevity/final_files/29Jun20_merged_alt_counts_allCHR.txt')
both=fread('/Genomics/ayroleslab2/shared/longevity/final_files/29Jun20_merged_both_counts_allCHR.txt')
info=read.delim('/Genomics/ayroleslab2/shared/longevity/final_files/29Jun20_merged_sample_info.txt')

# permute all time points
info$condition2<-NA
info$condition2<-as.character(sample(info$condition))

ITER1<-(13577*PERMNUMBER1)-13577+1
ITER2<-ITER1+13577

library('aod',lib='~/R')

pvals<-matrix(nrow=dim(alt)[1],ncol=2)

for (i in ITER1:ITER2){
	info$y=(t(both[i,-1]))
	info$x=(t(alt[i,-1]))
	
	# permute all time points
	tmp2<-subset(info,(y>0) & sex!='unknown' )
	mod1<-betabin(cbind(x, y - x) ~ condition2 + meta_cage + factor(sequencing_batch)+sex, ~1,data=tmp2)
	mod1_res<-attributes(summary(mod1))$Coef
	pvals[i,1]<-mod1_res[2,4]
	pvals[i,2]<-mod1_res[3,4] }

write.table(pvals,'/Genomics/ayroleslab2/shared/longevity/final_files/20Oct20_permute_PERMNUMBER1_PERMNUMBER2.txt',row.names=F,sep='\t')

