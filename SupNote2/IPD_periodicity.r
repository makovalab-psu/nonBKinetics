require(IWTomics)

load('IPD_forward.RData')

regionsFeatures_forward_quantiles=lapply(regionsFeatures_forward@features$IPD_Forward,
                                         function(feature) apply(feature,1,quantile,probs=c(0.05,0.25,0.5,0.75,0.95),na.rm=TRUE))
region_names=nameRegions(regionsFeatures_forward)

save(regionsFeatures_forward_quantiles,region_names,file='IPD_forward_quantiles.RData')




load('IPD_forward_quantiles.RData')

# autocorrelation plots
require(TSA)
acf_median=lapply(regionsFeatures_forward_quantiles,function(feature) acf(feature[3,],plot=FALSE))
pdf('autocorrelation_forward_median.pdf',width=10,height=5)
mapply(function(acf_results,name){
         plot(acf_results,xlim=c(1,20),ylim=c(-1,1),main='',xaxt='n',cex.lab=1.5,cex.axis=1.5)
         title(name,cex.main=3)
         axis(1,at=1:20,cex.axis=1.5)
       },acf_median,region_names)
dev.off()

# fast fourier transform plots
fft_median=lapply(regionsFeatures_forward_quantiles,function(feature) fft(feature[3,]))
fft_amplitude_median=lapply(fft_median,function(feature) Mod(feature)/length(feature))
pdf('fouriertransform_forward_median.pdf',width=10,height=5)
ymax=max(Reduce(cbind,fft_amplitude_median)[-1,])
mapply(function(fft_results,name){
         plot((1:99)/100,fft_results[-1],type='h',ylim=c(0,ymax),xaxt='n',
              cex.lab=1.5,cex.axis=1.5,cex.main=3,main=name,ylab='Amplitude',xlab='Frequency')
         axis(1,at=seq(0,1,0.1),cex.axis=1.5)
         return()
       },fft_amplitude_median,region_names)
dev.off()

