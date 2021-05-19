
###########
# how often is the derived allele decreasing or increasing in a random cage?
###########

setwd("~/Dropbox/Amanda_files/ayroles_lab/longevity_GXE/manuscript/final_files_and_code")

data_new=read.delim('data_files/SummaryTable_allsites_12Nov20.txt',stringsAsFactors=FALSE)
data_new$maf_quant<-cut(data_new$MAF_T0, breaks = seq(0,0.5,0.1), 
      include.lowest = TRUE, labels = 1:5)
data_new$ref_is_anc<-0
data_new$ref_is_anc[which(data_new$REF==data_new$ANC)]<-1
data_new$ref_is_anc[which(is.na(data_new$ANC))]<-NA

data_new_ns<-subset(data_new,sig_cat=='NS' & q_BH_perm_CTRL>0.2 & q_BH_perm_HS>0.2)
data_new_sig<-subset(data_new,sig_cat!='NS')

af=read.delim('data_files/20Aug20_AF_by_cage.txt')
af=merge(af,data_new[,c('site','MAF_T0','sig_cat','maf_quant','ref_is_anc')],by='site')
tmp<-subset(data_new,ANC!='NA')$site
af2=subset(af,af$site %in% tmp)

af1=unique(subset(af2,site %in% data_new_sig$site ))
af2=unique(subset(af2,site %in% data_new_ns$site ))

# in HS cages
af2$alt_dec_HS1A<-0
af2$alt_dec_HS2A<-0
af2$alt_dec_HS3D<-0

af2$alt_dec_HS1A[which(af2$af_t0_A - af2$af_tN_HS1_A > 0)]<-1
af2$alt_dec_HS2A[which(af2$af_t0_A - af2$af_tN_HS2_A > 0)]<-1
af2$alt_dec_HS3D[which(af2$af_t0_D - af2$af_tN_HS1_D > 0)]<-1

# in C cages
af2$alt_dec_C1A<-0
af2$alt_dec_C2A<-0
af2$alt_dec_C3D<-0

af2$alt_dec_C1A[which(af2$af_t0_A - af2$af_tN_C1_A > 0)]<-1
af2$alt_dec_C2A[which(af2$af_t0_A - af2$af_tN_C2_A > 0)]<-1
af2$alt_dec_C3D[which(af2$af_t0_D - af2$af_tN_C1_D > 0)]<-1

# subsample and MAF match

# maf in sig
freqs=as.data.frame(table(af1$maf_quant))

bias_HS<-c()
bias_C<-c()

for (i in 1:1000){

af2$keep<-0

for (k in 1:5){
freq_tmp<-subset(af2,maf_quant==k )
freq_keep<-freq_tmp$site[sample(1:dim(freq_tmp)[1],freqs$Freq[k])]
af2$keep[which(af2$site %in% freq_keep)]<-1
}

tmp3<-subset(af2,keep==1)

tmp1perm<-as.data.frame(table(tmp3$alt_dec_C1A,tmp3$ref_is_anc))
tmp2perm<-as.data.frame(table(tmp3$alt_dec_C2A,tmp3$ref_is_anc))

m1<-sum(subset(tmp1perm, Var1==Var2)$Freq)/sum(tmp1perm$Freq)
m2<-sum(subset(tmp2perm, Var1==Var2)$Freq)/sum(tmp2perm$Freq)

tmp1perm<-as.data.frame(table(tmp3$alt_dec_HS1A,tmp3$ref_is_anc))
tmp2perm<-as.data.frame(table(tmp3$alt_dec_HS2A,tmp3$ref_is_anc))

n1<-sum(subset(tmp1perm, Var1==Var2)$Freq)/sum(tmp1perm$Freq)
n2<-sum(subset(tmp2perm, Var1==Var2)$Freq)/sum(tmp2perm$Freq)

bias_HS<-c(bias_HS,n1,n2)
bias_C<-c(bias_C,m1,m2)

}

###########
# plot realtive to significant sites/observed data
###########

par(mfrow=c(1,1))

plot(density(c(bias_C,bias_HS)),xlim=c(0.45,0.7),bty='n')

x1=0.59; arrows(x1, 5, x1, 0,col='plum',lwd=2)
x1=0.65; arrows(x1, 5, x1, 0,col='darkorange',lwd=2)

