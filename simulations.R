library(MASS)
library(Matrix)
library(data.table)

setwd('~/Dropbox/Amanda_files/ayroles_lab/longevity_GXE/manuscript/final_files_and_code')
params_real=read.delim('~/Dropbox/Amanda_files/ayroles_lab/longevity_GXE/simulations/23Mar20_neg_bin_longevity.txt',header=T)
maf_real=fread('~/Dropbox/Amanda_files/ayroles_lab/longevity_GXE/manuscript/final_files_and_code/data_files/6Jul20_main_effects_ALL_covfilt_w_t0_freq.txt',header=T)

#########
# simulation function
#########

longevity_sim <- function(sites_n,n,af_diff_HS,af_diff_C,params2,maf) {
	
		maf_t0<-sample(maf,sites_n)
		
		# get genotypes per fly and locus
		geno_t0 = t(apply( as.matrix(maf_t0), 1, function(y) {rbinom(n,2,y) } ))
		geno_C = t(apply( as.matrix(maf_t0 - maf_t0*af_diff_C), 1, function(y) {rbinom(n,2,y) } ))
		geno_HS = t(apply( as.matrix(maf_t0 - maf_t0*af_diff_HS), 1, function(y) {rbinom(n,2,y) } ))
		 # get total coverage per fly and locus
		tot_counts_HS <- t(Reduce(cbind,lapply(1:sites_n,function(x) {
   		 params_tmp<-sample(1:dim(params2)[1],1)
   		 tot_counts <- rnegbin(n,params2[params_tmp,2],params2[params_tmp,1])*((x+1)-x)
  		 return(tot_counts) })))
   
  		 tot_counts_C <- t(Reduce(cbind,lapply(1:sites_n,function(x) {
   		 tot_counts <- rnegbin(n,params2[params_tmp,2],params2[params_tmp,1])*((x+1)-x)
  		 return(tot_counts) })))
   
   		 tot_counts_t0 <- t(Reduce(cbind,lapply(1:sites_n,function(x) {
   		 tot_counts <- rnegbin(n,params2[params_tmp,2],params2[params_tmp,1])*((x+1)-x)
  		 return(tot_counts) })))
   				
		# get reads mapped to alternate allele per fly and locus
		alt_counts_HS<-t(Reduce(cbind,lapply(1:sites_n, function(k) {
		out<-c()
		for (j in 1:n){out<-c(out,rbinom(1, tot_counts_HS[k,j], geno_HS[k,j]/2))}
		return(out) })))
		
		alt_counts_C<-t(Reduce(cbind,lapply(1:sites_n, function(k) {
		out<-c()
		for (j in 1:n){out<-c(out,rbinom(1, tot_counts_C[k,j], geno_C[k,j]/2))}
		return(out) })))
		
		alt_counts_t0<-t(Reduce(cbind,lapply(1:sites_n, function(k) {
		out<-c()
		for (j in 1:n){out<-c(out,rbinom(1, tot_counts_t0[k,j], geno_t0[k,j]/2))}
		return(out) })))
		return(list(alt_counts_HS=alt_counts_HS,alt_counts_C=alt_counts_C,alt_counts_t0=alt_counts_t0,tot_counts_HS=tot_counts_HS,tot_counts_C=tot_counts_C,tot_counts_t0=tot_counts_t0)) }
	
#########
# sanity checks
#########

sim1<-longevity_sim(500,100,0,0)

mean_cov_HS<-apply(sim1$tot_counts_HS,1,mean)
mean_cov_C<-apply(sim1$tot_counts_C,1,mean)
mean_cov_t0<-apply(sim1$tot_counts_t0,1,mean)

missing_HS<-apply(sim1$tot_counts_HS,1,function(x) length(which(x==0))/dim(sim1$tot_counts_HS)[2])
missing_C<-apply(sim1$tot_counts_C,1,function(x) length(which(x==0))/dim(sim1$tot_counts_C)[2])
missing_t0<-apply(sim1$tot_counts_t0,1,function(x) length(which(x==0))/dim(sim1$tot_counts_t0)[2])

af_HS<-apply(sim1$alt_counts_HS/sim1$tot_counts_HS,1,function(x) mean(x,na.rm=T))
af_C<-apply(sim1$alt_counts_C/sim1$tot_counts_C,1,function(x) mean(x,na.rm=T))
af_t0<-apply(sim1$alt_counts_t0/sim1$tot_counts_t0,1,function(x) mean(x,na.rm=T))

par(mfrow=c(3,3))
hist(mean_cov_HS,breaks=50,main='');hist(mean_cov_C,breaks=50,main='');hist(mean_cov_t0,breaks=50,main='')
hist(missing_HS,breaks=50,main='');hist(missing_C,breaks=50,main='');hist(missing_t0,breaks=50,main='')
hist(af_C-af_t0,breaks=50,main='');hist(af_HS-af_t0,breaks=50,main='')

#########
# run beta binomial
#########

diff<-c(seq(0,0.2,by=0.02))
tested_sites<-1000
tested_flies<-1000

c=c(rep(0,tested_flies),rep(1,tested_flies),rep(2,tested_flies))

out1<-matrix(nrow=tested_sites,ncol=length(diff))
out2<-matrix(nrow=tested_sites,ncol=length(diff))

for (q in 1:length(diff)) {
sim1<-longevity_sim(tested_sites,tested_flies,diff[q],diff[q],params_real,maf_real$freq_maf)
print(q)
library(aod)

# for CMH
#alt_tmp1<-cbind(sim1$alt_counts_C,sim1$alt_counts_t0[,1:500])
#alt_tmp2<-cbind(sim1$alt_counts_HS,sim1$alt_counts_t0[,501:1000])

#tot_tmp1<-cbind(sim1$tot_counts_C-sim1$alt_counts_C,sim1$tot_counts_t0[,1:500]-sim1$alt_counts_t0[,1:500])
	#tot_tmp2<-cbind(sim1$tot_counts_HS-sim1$alt_counts_HS,sim1$tot_counts_t0[,501:1000]-sim1$alt_counts_t0[,501:1000])

#write.table(alt_tmp1,paste('data_files/3Sep20_sim_alt1_forCMH_',diff[q],'.txt',sep=''),row.names=F,sep='\t')
#write.table(alt_tmp2,paste('data_files/3Sep20_sim_alt2_forCMH_',diff[q],'.txt',sep=''),row.names=F,sep='\t')

#write.table(tot_tmp1,paste('data_files/3Sep20_sim_ref1_forCMH_',diff[q],'.txt',sep=''),row.names=F,sep='\t')
#write.table(tot_tmp2,paste('data_files/3Sep20_sim_ref2_forCMH_',diff[q],'.txt',sep=''),row.names=F,sep='\t')

for (p in 1:tested_sites) {
	x=as.vector(c(sim1$alt_counts_t0[p,],sim1$alt_counts_HS[p,],sim1$alt_counts_C[p,]))
	y=as.vector(c(sim1$tot_counts_t0[p,],sim1$tot_counts_HS[p,],sim1$tot_counts_C[p,]))
	tmp<-as.data.frame(cbind(x,y,c))
	tmp<-subset(tmp,y>0)
	mod1<-betabin(cbind(x, y - x) ~ as.factor(c), ~1,data=tmp)
	out1[p,q]<-attributes(summary(mod1))$Coef[2,4]
	out2[p,q]<-attributes(summary(mod1))$Coef[3,4]

} }
	
	
plot(diff,apply(out1,2,function(x) length(which(x<10^-6))/tested_sites ),xlab='difference in AF (% decrease tN-t0)',ylab='prop p-values < 10^-6')
write.table(out1,'~/Dropbox/Amanda_files/ayroles_lab/longevity_GXE/manuscript/final_files_and_code/CMH_bb_compare/3Sep20_bb_simulations1.txt',row.names=F,sep='\t')
write.table(out2,'~/Dropbox/Amanda_files/ayroles_lab/longevity_GXE/manuscript/final_files_and_code/CMH_bb_compare/3Sep20_bb_simulations2.txt',row.names=F,sep='\t')


