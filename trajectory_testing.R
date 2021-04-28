
###########
# read in count data - beginning and end time points, int time points
###########

library(data.table)

alt=fread('/Genomics/ayroleslab2/shared/longevity/final_files/29Jun20_merged_alt_counts_allCHR.txt')
both=fread('/Genomics/ayroleslab2/shared/longevity/final_files/29Jun20_merged_both_counts_allCHR.txt')
info=read.delim('/Genomics/ayroleslab2/shared/longevity/final_files/29Jun20_merged_sample_info.txt')

# check that sample/individual IDs match up properly
identical(names(alt[,-1]),names(both[,-1]))
identical(as.character(info$sample_ID),names(both[,-1]))

both2=fread('/Genomics/ayroleslab2/shared/longevity/final_files/27Jul20_tot_counts_int_time_pts.txt')
alt2=fread('/Genomics/ayroleslab2/shared/longevity/final_files/27Jul20_alt_counts_int_time_pts.txt')

info2=read.delim('/Genomics/ayroleslab2/shared/longevity/merged_counts/20Mar20_samples_int_time_points.txt')
info2=info2[order(info2$col_id),]
identical(as.character(info2$sample_id),as.character(names(both2)[-1]))
info2$condition<-'C'
info2$condition[which(info2$sub_cage %in% 'HS1')]<-'HS'
info2$sub_cage2=paste(info2$meta_cage,info2$sub_cage,sep='_')

sig=fread('6Jul20_main_effects_ALL_covfilt.txt')
sig$FDR_CTRL<-p.adjust(sig$pval_CTRL,method='BH')
sig$FDR_HS<-p.adjust(sig$pval_HS,method='BH')

info1_tmp<-info[,c('sample_ID','condition','timepoint','meta_cage','sequencing_batch','sex')]
info1_tmp$timepoint2<-0
info1_tmp$timepoint2[which(info1_tmp$timepoint=='TN')]<-10
info1_tmp$meta_cage<-paste('Cage',info1_tmp$meta_cage,sep='')
names(info1_tmp)[1]<-'sample_id'

info2$sequencing_batch<-3
all_info=rbind(info1_tmp[,c('sample_id','condition','timepoint2','meta_cage','sequencing_batch','sex')],info2[,c('sample_id','condition','timepoint2','meta_cage','sequencing_batch','sex')])


###########
# SNPs significant in both conditions
###########

library(aod)

sig$cat<- 'NS'
sig$cat[which( (sig$FDR_CTRL<0.01 & sig$pval_HS<0.05) | (sig$FDR_HS<0.01 & sig$pval_CTRL<0.05)  )] <- 'both'
sig$cat[which(sig$pval_CTRL>0.05 & sig$FDR_HS<0.01)] <- 'hs'
sig$cat[which(sig$pval_HS>0.05 & sig$FDR_CTRL<0.01)] <- 'ctrl'

both_sig<-subset(sig, cat=='both')
both_sig_alt<-alt2[which(alt2$site %in% both_sig$site),]
both_sig_tot<-both2[which(both2$site %in% both_sig$site),]

both_sig_alt1<-alt[which(alt$site %in% both_sig$site),]
both_sig_tot1<-both[which(both$site %in% both_sig$site),]

all_alt=(merge(both_sig_alt1,both_sig_alt,by='site'))
all_both=(merge(both_sig_tot1,both_sig_tot,by='site'))

identical(names(all_alt),names(all_both))
identical(names(all_alt)[-1],as.character(all_info$sample_id))

pvals<-matrix(nrow=dim(all_both)[1],ncol=4)
coefs_lik<-matrix(nrow=dim(all_both)[1],ncol=8)

for (i in 1:dim(all_both)[1]){
	all_info$y=(t(all_both[i,-1]))
	all_info$x=(t(all_alt[i,-1]))
	tmp2<-subset(all_info,(x>0 | y>0) & sex!='unknown')
	# model1 - linear
	mod1<-betabin(cbind(x, y - x) ~ timepoint2 + meta_cage + factor(sequencing_batch) + sex, ~1,data=tmp2)
	mod1_res<-attributes(summary(mod1))$Coef
	# model2 - quadratic
	mod2<-betabin(cbind(x, y - x) ~ timepoint2 + I(timepoint2^2) + meta_cage + factor(sequencing_batch) +sex , ~1,data=tmp2)
	mod2_res<-attributes(summary(mod2))$Coef
	# model3 - breakpoint
	# breaks <- c(1,2,4,6,7,8,9)
	breaks <- c(4,6,7,8,9)
	mse <- numeric(length(breaks))
		
		for (k in 1:length(breaks)){
		tmp2$timepoint_tmp<-tmp2$timepoint2
		tmp2$timepoint_tmp[which(tmp2$timepoint2<breaks[k])]<-breaks[k]
		piecewise1 <- betabin(cbind(x, y - x) ~ timepoint_tmp + meta_cage + factor(sequencing_batch) +sex, ~1,data=tmp2)
	 	mse[k] <-as.numeric(attributes(AIC(piecewise1))$istats[2]) }
	mse <- as.numeric(mse)

	pick<-breaks[which(mse==min(mse))]
	tmp2$timepoint_tmp[which(tmp2$timepoint2<pick)]<-pick
	mod3 <- betabin(cbind(x, y - x) ~ timepoint_tmp + meta_cage + factor(sequencing_batch) +sex, ~1,data=tmp2)
	mod3_res<-attributes(summary(mod3))$Coef

	pvals[i,1]<-mod1_res[2,4]
	pvals[i,2]<-mod2_res[2,4]
	pvals[i,3]<-mod2_res[3,4]
	pvals[i,4]<-mod3_res[2,4]

	coefs_lik[i,1]<-mod1_res[2,1]
	coefs_lik[i,2]<-mod2_res[2,1]
	coefs_lik[i,3]<-mod2_res[3,1]
	coefs_lik[i,4]<-mod3_res[2,1]
	coefs_lik[i,5]<-pick

	coefs_lik[i,6]<-as.numeric(attributes(AIC(mod1))$istats[2])
	coefs_lik[i,7]<-as.numeric(attributes(AIC(mod2))$istats[2])
	coefs_lik[i,8]<-as.numeric(attributes(AIC(mod3))$istats[2])

print(i) }

pvals<-as.data.frame(pvals)
coefs_lik<-as.data.frame(coefs_lik)
coefs_lik$site<-all_both$site
pvals$site<-all_both$site

write.table(pvals,'/Genomics/ayroleslab2/shared/longevity/final_files/27Jun20_pvals_both_SNPs.txt',row.names=F,sep='\t')
write.table(coefs_lik,'/Genomics/ayroleslab2/shared/longevity/final_files/27Jun20_coefs_both_SNPs.txt',row.names=F,sep='\t')

###########
# SNPs significant in HS
###########

both_sig<-subset(sig, cat=='hs')
both_sig_alt<-alt2[which(alt2$site %in% both_sig$site),]
both_sig_tot<-both2[which(both2$site %in% both_sig$site),]

both_sig_alt1<-alt[which(alt$site %in% both_sig$site),]
both_sig_tot1<-both[which(both$site %in% both_sig$site),]

all_alt=(merge(both_sig_alt1,both_sig_alt,by='site'))
all_both=(merge(both_sig_tot1,both_sig_tot,by='site'))

identical(names(all_alt),names(all_both))
identical(names(all_alt)[-1],as.character(all_info$sample_id))

pvals<-matrix(nrow=dim(all_both)[1],ncol=4)
coefs_lik<-matrix(nrow=dim(all_both)[1],ncol=8)

for (i in 1:dim(all_both)[1]){
	all_info$y=(t(all_both[i,-1]))
	all_info$x=(t(all_alt[i,-1]))
	tmp2<-subset(all_info,(x>0 | y>0) & sex!='unknown' & condition!='C')
	# model1 - linear
	mod1<-betabin(cbind(x, y - x) ~ timepoint2 + meta_cage + factor(sequencing_batch)+sex, ~1,data=tmp2)
	mod1_res<-attributes(summary(mod1))$Coef
	# model2 - quadratic
	mod2<-betabin(cbind(x, y - x) ~ timepoint2 + I(timepoint2^2) + meta_cage + factor(sequencing_batch)+sex , ~1,data=tmp2)
	mod2_res<-attributes(summary(mod2))$Coef
	# model3 - breakpoint
	# breaks <- c(1,2,4,6,7)
	breaks <- c(4,6,7)
	mse <- numeric(length(breaks))
		
		for (k in 1:length(breaks)){
		tmp2$timepoint_tmp<-tmp2$timepoint2
		tmp2$timepoint_tmp[which(tmp2$timepoint2<breaks[k])]<-breaks[k]
		piecewise1 <- betabin(cbind(x, y - x) ~ timepoint_tmp + meta_cage + factor(sequencing_batch)+sex, ~1,data=tmp2)
	 	mse[k] <-as.numeric(attributes(AIC(piecewise1))$istats[2]) }
	mse <- as.numeric(mse)

	pick<-breaks[which(mse==min(mse))]
	tmp2$timepoint_tmp[which(tmp2$timepoint2<pick)]<-pick
	mod3 <- betabin(cbind(x, y - x) ~ timepoint_tmp + meta_cage + factor(sequencing_batch) +sex, ~1,data=tmp2)
	mod3_res<-attributes(summary(mod3))$Coef

	pvals[i,1]<-mod1_res[2,4]
	pvals[i,2]<-mod2_res[2,4]
	pvals[i,3]<-mod2_res[3,4]
	pvals[i,4]<-mod3_res[2,4]

	coefs_lik[i,1]<-mod1_res[2,1]
	coefs_lik[i,2]<-mod2_res[2,1]
	coefs_lik[i,3]<-mod2_res[3,1]
	coefs_lik[i,4]<-mod3_res[2,1]
	coefs_lik[i,5]<-pick

	coefs_lik[i,6]<-as.numeric(attributes(AIC(mod1))$istats[2])
	coefs_lik[i,7]<-as.numeric(attributes(AIC(mod2))$istats[2])
	coefs_lik[i,8]<-as.numeric(attributes(AIC(mod3))$istats[2])
print(i) }

pvals<-as.data.frame(pvals)
coefs_lik<-as.data.frame(coefs_lik)
coefs_lik$site<-all_both$site
pvals$site<-all_both$site

write.table(pvals,'/Genomics/ayroleslab2/shared/longevity/final_files/27Jun20_pvals_HS_SNPs.txt',row.names=F,sep='\t')
write.table(coefs_lik,'/Genomics/ayroleslab2/shared/longevity/final_files/27Jun20_coefs_HS_SNPs.txt',row.names=F,sep='\t')

###########
# SNPs significant in CTRL
###########

both_sig<-subset(sig, cat=='ctrl')
both_sig_alt<-alt2[which(alt2$site %in% both_sig$site),]
both_sig_tot<-both2[which(both2$site %in% both_sig$site),]

both_sig_alt1<-alt[which(alt$site %in% both_sig$site),]
both_sig_tot1<-both[which(both$site %in% both_sig$site),]

all_alt=(merge(both_sig_alt1,both_sig_alt,by='site'))
all_both=(merge(both_sig_tot1,both_sig_tot,by='site'))

identical(names(all_alt),names(all_both))
identical(names(all_alt)[-1],as.character(all_info$sample_id))

pvals<-matrix(nrow=dim(all_both)[1],ncol=4)
coefs_lik<-matrix(nrow=dim(all_both)[1],ncol=8)

for (i in 1:dim(all_both)[1]){
	all_info$y=(t(all_both[i,-1]))
	all_info$x=(t(all_alt[i,-1]))
	tmp2<-subset(all_info,(x>0 | y>0)  & condition!='HS')
	# model1 - linear
	mod1<-betabin(cbind(x, y - x) ~ timepoint2 + meta_cage + factor(sequencing_batch), ~1,data=tmp2)
	mod1_res<-attributes(summary(mod1))$Coef
	# model2 - quadratic
	mod2<-betabin(cbind(x, y - x) ~ timepoint2 + I(timepoint2^2) + meta_cage + factor(sequencing_batch) , ~1,data=tmp2)
	mod2_res<-attributes(summary(mod2))$Coef
	# model3 - breakpoint
	# breaks <- c(1,2,4,6,7,8,9)
	breaks <- c(4,6,7,8,9)
	mse <- numeric(length(breaks))
		
		for (k in 1:length(breaks)){
		tmp2$timepoint_tmp<-tmp2$timepoint2
		tmp2$timepoint_tmp[which(tmp2$timepoint2<breaks[k])]<-breaks[k]
		piecewise1 <- betabin(cbind(x, y - x) ~ timepoint_tmp + meta_cage + factor(sequencing_batch), ~1,data=tmp2)
	 	mse[k] <-as.numeric(attributes(AIC(piecewise1))$istats[2]) }
	mse <- as.numeric(mse)

	pick<-breaks[which(mse==min(mse))]
	tmp2$timepoint_tmp[which(tmp2$timepoint2<pick)]<-pick
	mod3 <- betabin(cbind(x, y - x) ~ timepoint_tmp + meta_cage + factor(sequencing_batch) , ~1,data=tmp2)
	mod3_res<-attributes(summary(mod3))$Coef

	pvals[i,1]<-mod1_res[2,4]
	pvals[i,2]<-mod2_res[2,4]
	pvals[i,3]<-mod2_res[3,4]
	pvals[i,4]<-mod3_res[2,4]

	coefs_lik[i,1]<-mod1_res[2,1]
	coefs_lik[i,2]<-mod2_res[2,1]
	coefs_lik[i,3]<-mod2_res[3,1]
	coefs_lik[i,4]<-mod3_res[2,1]
	coefs_lik[i,5]<-pick

	coefs_lik[i,6]<-as.numeric(attributes(AIC(mod1))$istats[2])
	coefs_lik[i,7]<-as.numeric(attributes(AIC(mod2))$istats[2])
	coefs_lik[i,8]<-as.numeric(attributes(AIC(mod3))$istats[2])
print(i) }

pvals<-as.data.frame(pvals)
coefs_lik<-as.data.frame(coefs_lik)
coefs_lik$site<-all_both$site
pvals$site<-all_both$site

write.table(pvals,'/Genomics/ayroleslab2/shared/longevity/final_files/27Jun20_pvals_C_SNPs.txt',row.names=F,sep='\t')
write.table(coefs_lik,'/Genomics/ayroleslab2/shared/longevity/final_files/27Jun20_coefs_C_SNPs.txt',row.names=F,sep='\t')
