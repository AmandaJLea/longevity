#!/bin/bash Rscript

# the full list of sites were broken up into parts and run in parallel
# for f in `seq 1 5000 280000`; do cat bb.sh | sed -e s/STARTNUMBER/$f/g > bb.$f.sh; done
# for f in `seq 1 5000 280000`; do cat bb2.sh | sed -e s/STARTNUMBER/$f/g > bb2.$f.sh; done
# rm commands1.sh; touch commands1.sh; for f in `seq 1 5000 280000`; do echo "sh bb2.$f.sh" >> commands1.sh; done

library(data.table)
library(aod)
ITER1=STARTNUMBER
ITER2=ITER1+5000

# count data to use for beta binomial modeling
alt=fread('/Genomics/ayroleslab2/shared/longevity/final_files/29Jun20_merged_alt_counts_allCHR.txt')
both=fread('/Genomics/ayroleslab2/shared/longevity/final_files/29Jun20_merged_both_counts_allCHR.txt')
info=read.delim('/Genomics/ayroleslab2/shared/longevity/final_files/29Jun20_merged_sample_info.txt')

# check that sample/individual IDs match up properly
identical(names(alt[,-1]),names(both[,-1]))
identical(as.character(info$sample_ID),names(both[,-1]))

############
# beta binomial models - main effects
############

pvals1<-matrix(nrow=dim(alt)[1],ncol=6)
coefs1<-matrix(nrow=dim(alt)[1],ncol=6)
se1<-matrix(nrow=dim(alt)[1],ncol=6)
lik1<-matrix(nrow=dim(alt)[1],ncol=2)

# this is run as a loop b/c every once in a while a site has convergence issues; it is very slow 
for (i in ITER1:ITER2){
	info$tmp_x=as.vector(t(alt[i,-1]))
	info$tmp_y=as.vector(t(both[i,-1]))
	# exclude samples with no reads at a given site, and samples that could not be sexed 
	tmp_info<-subset(info,tmp_y>0 & sex!='unknown')
	mod1<-tryCatch(betabin(cbind(tmp_x, tmp_y - tmp_x) ~ condition  + sequencing_batch + meta_cage + sex, ~1,data=tmp_info), error=function(x){})
	mod2<-tryCatch(betabin(cbind(tmp_x, tmp_y - tmp_x) ~ sequencing_batch + meta_cage + sex, ~1,data=tmp_info), error=function(x){})
	tryCatch(lik1[i,1]<-logLik(mod1), error=function(x){})
	tryCatch(lik1[i,2]<-logLik(mod2), error=function(x){})
	tryCatch(pvals1[i,]<-attributes(summary(mod1))$Coef[,4] , error=function(x){})
	tryCatch(se1[i,]<-attributes(summary(mod1))$Coef[,2] , error=function(x){})
	tryCatch(coefs1[i,]<-attributes(summary(mod1))$Coef[,1] , error=function(x){})
print(i) }

write.table(pvals1,'29Jun20_main_effects_pvals_STARTNUMBER.txt',row.names=F,sep='\t',col.names=F)
write.table(coefs1,'29Jun20_main_effects_coefs_STARTNUMBER.txt',row.names=F,sep='\t',col.names=F)
write.table(lik1,'29Jun20_main_effects_loglik_STARTNUMBER.txt',row.names=F,sep='\t',col.names=F)
write.table(se1,'29Jun20_main_effects_se_STARTNUMBER.txt',row.names=F,sep='\t',col.names=F)
