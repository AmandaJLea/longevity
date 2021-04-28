setwd('~/Dropbox/Amanda_files/ayroles_lab/longevity_GXE/manuscript/final_files_and_code/')

data=read.delim('data_files/6Jul20_main_effects_coef_ALL.txt')
sig=read.delim('data_files/SummaryTable_allsites_12Nov20.txt')
sig=subset(sig,sig_cat!='NS')
af=read.delim('data_files/14Aug20_AF_summary.txt')

data$tmp_C_t0<-data$intercept+data$CTRL*0+data$batch*1+data$cage*1+data$sex*1
data$tmp_C_tN<-data$intercept+data$CTRL*1+data$batch*1+data$cage*1+data$sex*1

data$tmp_HS_t0<-data$intercept+data$HS*0+data$batch*1+data$cage*1+data$sex*1
data$tmp_HS_tN<-data$intercept+data$HS*1+data$batch*1+data$cage*1+data$sex*1

data$tmp2_C_t0<-1/(1+exp(-data$tmp_C_t0) )
data$tmp2_C_tN<-1/(1+exp(-data$tmp_C_tN) )

data$tmp2_HS_t0<-1/(1+exp(-data$tmp_HS_t0) )
data$tmp2_HS_tN<-1/(1+exp(-data$tmp_HS_tN) )

data=merge(data,af,by='site')

mean(abs( data$af_tN_C - data$af_t0))
mean(abs( data$af_tN_HS - data$af_t0))
sd(abs( data$af_tN_HS - data$af_t0))

null=fread('~/Desktop/30Dec20_subsample_AF_summary.txt')

#########
# plot
#########

par(mfrow=c(2,2))
data2=subset(data,site %in% subset(sig,sig_cat=='shared')$site)
plot(density(abs( data2$af_tN_C - data2$af_t0)),col='navyblue',lwd=2,ylim=c(0,50),main='',xlab='Allele frequency change',xlim=c(0,0.25),bty='n')
lines(density(abs( data2$af_tN_HS - data2$af_t0)),col='navyblue',lwd=2,lty=2)
data2=subset(data,site %in% subset(sig,sig_cat=='HS')$site)
lines(density(abs( data2$af_tN_HS - data2$af_t0)),col='darkorange',lwd=2)
lines(density(abs( data$af_tN_HS - data$af_t0)),col='grey',lwd=2,lty=2)
lines(density(abs( data$af_tN_C - data$af_t0)),col='grey',lwd=2)

legend('topright',c('NS: CTRL','NS: HS','Shared: CTRL','Shared: HS','GxE: HS'),col=c('grey','grey','navyblue','navyblue','darkorange'),lty=c(1,2,1,2,1,1),lwd=c(2,2,2,2,2,2),bty='n',cex=0.75)

data2=subset(data,site %in% subset(sig,sig_cat=='shared')$site)
plot(density(abs( data2$tmp2_C_tN - data2$tmp2_C_t0)),col='navyblue',lwd=2,ylim=c(0,50),main='',xlab='Allele frequency change',xlim=c(0,0.25),bty='n')
lines(density(abs( data2$tmp2_HS_tN - data2$tmp2_HS_t0)),col='navyblue',lwd=2,lty=2)
data2=subset(data,site %in% subset(sig,sig_cat=='HS')$site)
lines(density(abs( data2$tmp2_HS_tN - data2$tmp2_HS_t0)),col='darkorange',lwd=2)
lines(density(abs( data$tmp2_HS_tN - data$tmp2_HS_t0)),col='grey',lwd=2,lty=2)
lines(density(abs( data$tmp2_C_tN - data$tmp2_C_t0)),col='grey',lwd=2)

legend('topright',c('NS: CTRL','NS: HS','Shared: CTRL','Shared: HS','GxE: HS'),col=c('grey','grey','navyblue','navyblue','darkorange'),lty=c(1,2,1,2,1,1),lwd=c(2,2,2,2,2,2),bty='n',cex=0.75)

data2=subset(data,site %in% subset(sig,sig_cat=='shared')$site)
plot(density(abs( data2$af_tN_C - data2$af_t0)),col='navyblue',lwd=2,ylim=c(0,50),main='',xlab='Allele frequency change',xlim=c(0,0.25),bty='n')
lines(density(abs( data2$af_tN_HS - data2$af_t0)),col='navyblue',lwd=2,lty=2)
data2=subset(data,site %in% subset(sig,sig_cat=='HS')$site)
lines(density(abs( data2$af_tN_HS - data2$af_t0)),col='darkorange',lwd=2)

for (i in 1:50){
null_tmp<-as.matrix(abs( null[,i,with=F]))
lines(density( null_tmp),col='lightgrey',lwd=2)
}
legend('topright',c('Subsampled','Shared: CTRL','Shared: HS','GxE: HS'),col=c('grey','navyblue','navyblue','darkorange'),lty=c(1,1,2,1,1),lwd=c(2,2,2,2,2),bty='n',cex=0.75)


