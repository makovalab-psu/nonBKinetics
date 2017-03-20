require(IWTomics)
source('IWTomics_modified_functions.r')



#######################
##### PLUS STRAND #####
#######################

### Complete datasets ###

### Mean test statistics ###

# Test for difference using IWT
# 10 random subsample of maximum 10.000 regions for each region dataset
# see load_IPD_data.r
for(i in 1:10){
  load(paste0('IPD_forward_max10000_',i,'.RData'))
  
  # Test all features
  result_mean=IWTomicsTest(regionsFeatures_forward,
                           setdiff(idRegions(regionsFeatures_forward),'Control'),rep('Control',nRegions(regionsFeatures_forward)-1),
                           statistics='mean',B=10000)
  save(result_mean,file=paste0('IPD_forward_max10000_',i,'_results_mean.RData'))
  
  # Plot test results
  pdf(paste0('IPD_forward_max10000_',i,'_results_mean.pdf'),width=8.5,height=10)
  plotTest_same_sample_size(result_mean,ask=FALSE,col=c(rep('blue',nRegions(result_mean)-1),'red'),xlab='Coordinates',
                            average=FALSE,size=TRUE,size_perc=TRUE,position=TRUE,cex.axis=1.5,cex.lab=1.5,cex.main=2,probs=c(0.05,0.25,0.5,0.75,0.95))
  dev.off()
  plotSummary(result_mean,groupby='feature',filenames=paste0('IPD_forward_max10000_',i,'_results_mean_summary.pdf'),
              only_significant=FALSE,gaps_tests=c(7,11,17,37),align_lab='Feature center',ask=FALSE,xlab='Coordinates',cellwidth=10,cellheight=10)
  
  # Test all features
  # number of control regions equal to feature regions
  result_mean=IWTomicsTest_same_sample_size(regionsFeatures_forward,
                                            setdiff(idRegions(regionsFeatures_forward),'Control'),rep('Control',nRegions(regionsFeatures_forward)-1),
                                            statistics='mean',B=10000)
  save(result_mean,file=paste0('IPD_forward_max10000_',i,'_results_same_sample_size_mean.RData'))
  
  # Plot test results
  pdf(paste0('IPD_forward_max10000_',i,'_results_same_sample_size_mean.pdf'),width=8.5,height=10)
  plotTest_same_sample_size(result_mean,ask=FALSE,col=c(rep('blue',nRegions(result_mean)-1),'red'),xlab='Coordinates',
                            average=FALSE,size=TRUE,size_perc=TRUE,position=TRUE,cex.axis=1.5,cex.lab=1.5,cex.main=2,probs=c(0.05,0.25,0.5,0.75,0.95))
  dev.off()
  plotSummary(result_mean,groupby='feature',filenames=paste0('IPD_forward_max10000_',i,'_results_same_sample_size_mean_summary.pdf'),
              only_significant=FALSE,gaps_tests=c(7,11,17,37),align_lab='Feature center',ask=FALSE,xlab='Coordinates',cellwidth=10,cellheight=10)
}


# Plot only reproducible results
alpha=0.05
min_significant=1
load('IPD_forward_max10000_1_results_mean.RData')
reproducible_pval_matrix=lapply(seq_len(nTests(result_mean)),function(name) matrix(TRUE,100,100))
for(i in 1:10){
  load(paste0('IPD_forward_max10000_',i,'_results_mean.RData'))
  reproducible_pval_matrix=mapply(function(reproducible_pval,result){
    reproducible_pval=(reproducible_pval & (result$IPD_Forward$adjusted_pval_matrix<=alpha))
  },reproducible_pval_matrix,result_mean@test$result,SIMPLIFY=FALSE)
  load(paste0('IPD_forward_max10000_',i,'_results_same_sample_size_mean.RData'))
  reproducible_pval_matrix=mapply(function(reproducible_pval,result){
    reproducible_pval=(reproducible_pval & (result$IPD_Forward$adjusted_pval_matrix<=alpha))
  },reproducible_pval_matrix,result_mean@test$result,SIMPLIFY=FALSE)
}
reproducible_pval_matrix=mapply(function(reproducible_pval){
  reproducible_pval[rowSums(reproducible_pval)<min_significant,]=matrix(FALSE,nrow=sum(rowSums(reproducible_pval)<min_significant),ncol=100)
  return(reproducible_pval)
},reproducible_pval_matrix,SIMPLIFY=FALSE)
save(reproducible_pval_matrix,file='IPD_forward_max10000_reproducible_results_mean.RData')

load('IPD_forward_max10000_1_results_mean.RData')
plotSummary_reproducible(result_mean,groupby='feature',filenames=paste0('IPD_forward_max10000_reproducible_results_mean_summary.pdf'),
                         reproducible_pval_matrix=reproducible_pval_matrix,only_significant=TRUE,
                         gaps_tests=c(7,11,17,37),align_lab='Feature center',xlab='',ylab='',ask=FALSE,cellwidth=5,cellheight=12,fontsize=12)
load('IPD_forward.RData')
pdf('boxplot_vs_control_forward_reproducible_results_mean.pdf',width=9,height=6)
for(i in seq_len(nrow(regionsFeatures_forward)-1)){
  plot_no_size_control(regionsFeatures_forward,reproducible_pval=reproducible_pval_matrix[[i]],
                       type='boxplot',probs=c(0.05,0.25,0.5,0.75,0.95),size_perc=TRUE,position=TRUE,average=FALSE,cex=1.4,cex.axis=1.4,cex.lab=1.4,cex.main=2.5,
                       id_regions_subset=c('Control',idRegions(regionsFeatures_forward)[i]),col=c('red','blue'),ask=FALSE,xlab='Coordinates')
}
dev.off()





### Median test statistics ###

# Test for difference using IWT
# 10 random subsample of maximum 10.000 regions for each region dataset
# see load_IPD_data.r
for(i in 1:10){
  load(paste0('IPD_forward_max10000_',i,'.RData'))
  
  # Test all features
  result_median=IWTomicsTest(regionsFeatures_forward,
                             setdiff(idRegions(regionsFeatures_forward),'Control'),rep('Control',nRegions(regionsFeatures_forward)-1),
                             statistics='median',B=10000)
  save(result_median,file=paste0('IPD_forward_max10000_',i,'_results_median.RData'))
  
  # Plot test results
  pdf(paste0('IPD_forward_max10000_',i,'_results_median.pdf'),width=8.5,height=10)
  plotTest_same_sample_size(result_median,ask=FALSE,col=c(rep('blue',nRegions(result_median)-1),'red'),xlab='Coordinates',
                            average=FALSE,size=TRUE,size_perc=TRUE,position=TRUE,cex.axis=1.5,cex.lab=1.5,cex.main=2,probs=c(0.05,0.25,0.5,0.75,0.95))
  dev.off()
  plotSummary(result_median,groupby='feature',filenames=paste0('IPD_forward_max10000_',i,'_results_median_summary.pdf'),
              only_significant=FALSE,gaps_tests=c(7,11,17,37),align_lab='Feature center',ask=FALSE,xlab='Coordinates',cellwidth=10,cellheight=10)
  
  # Test all features
  # number of control regions equal to feature regions
  result_median=IWTomicsTest_same_sample_size(regionsFeatures_forward,
                                              setdiff(idRegions(regionsFeatures_forward),'Control'),rep('Control',nRegions(regionsFeatures_forward)-1),
                                              statistics='median',B=10000)
  save(result_median,file=paste0('IPD_forward_max10000_',i,'_results_same_sample_size_median.RData'))
  
  # Plot test results
  pdf(paste0('IPD_forward_max10000_',i,'_results_same_sample_size_median.pdf'),width=8.5,height=10)
  plotTest_same_sample_size(result_median,ask=FALSE,col=c(rep('blue',nRegions(result_median)-1),'red'),xlab='Coordinates',
                            average=FALSE,size=TRUE,size_perc=TRUE,position=TRUE,cex.axis=1.5,cex.lab=1.5,cex.main=2,probs=c(0.05,0.25,0.5,0.75,0.95))
  dev.off()
  plotSummary(result_median,groupby='feature',filenames=paste0('IPD_forward_max10000_',i,'_results_same_sample_size_median_summary.pdf'),
              only_significant=FALSE,gaps_tests=c(7,11,17,37),align_lab='Feature center',ask=FALSE,xlab='Coordinates',cellwidth=10,cellheight=10)
}


# Plot only reproducible results
alpha=0.05
min_significant=1
load('IPD_forward_max10000_1_results_median.RData')
reproducible_pval_matrix=lapply(seq_len(nTests(result_median)),function(name) matrix(TRUE,100,100))
for(i in 1:10){
  load(paste0('IPD_forward_max10000_',i,'_results_median.RData'))
  reproducible_pval_matrix=mapply(function(reproducible_pval,result){
    reproducible_pval=(reproducible_pval & (result$IPD_Forward$adjusted_pval_matrix<=alpha))
  },reproducible_pval_matrix,result_median@test$result,SIMPLIFY=FALSE)
  load(paste0('IPD_forward_max10000_',i,'_results_same_sample_size_median.RData'))
  reproducible_pval_matrix=mapply(function(reproducible_pval,result){
    reproducible_pval=(reproducible_pval & (result$IPD_Forward$adjusted_pval_matrix<=alpha))
  },reproducible_pval_matrix,result_median@test$result,SIMPLIFY=FALSE)
}
reproducible_pval_matrix=mapply(function(reproducible_pval){
  reproducible_pval[rowSums(reproducible_pval)<min_significant,]=matrix(FALSE,nrow=sum(rowSums(reproducible_pval)<min_significant),ncol=100)
  return(reproducible_pval)
},reproducible_pval_matrix,SIMPLIFY=FALSE)
save(reproducible_pval_matrix,file='IPD_forward_max10000_reproducible_results_median.RData')

load('IPD_forward_max10000_1_results_median.RData')
plotSummary_reproducible(result_median,groupby='feature',filenames=paste0('IPD_forward_max10000_reproducible_results_median_summary.pdf'),
                         reproducible_pval_matrix=reproducible_pval_matrix,only_significant=TRUE,
                         gaps_tests=c(7,11,17,37),align_lab='Feature center',xlab='',ylab='',ask=FALSE,cellwidth=5,cellheight=12,fontsize=12)
load('IPD_forward.RData')
pdf('boxplot_vs_control_forward_reproducible_results_median.pdf',width=9,height=6)
for(i in seq_len(nrow(regionsFeatures_forward)-1)){
  plot_no_size_control(regionsFeatures_forward,reproducible_pval=reproducible_pval_matrix[[i]],
                       type='boxplot',probs=c(0.05,0.25,0.5,0.75,0.95),size_perc=TRUE,position=TRUE,average=FALSE,cex=1.4,cex.axis=1.4,cex.lab=1.4,cex.main=2.5,
                       id_regions_subset=c('Control',idRegions(regionsFeatures_forward)[i]),col=c('red','blue'),ask=FALSE,xlab='Coordinates')
}
dev.off()





### Quantile test statistics (quantiles 0.05, 0.25, 0.50, 0.75 and 0.95) ###

# Test for difference using IWT
# 10 random subsample of maximum 10.000 regions for each region dataset
# see load_IPD_data.r
for(i in 1:10){
  load(paste0('IPD_forward_max10000_',i,'.RData'))
  
  # Test all features
  result_multi_quantile=IWTomicsTest(regionsFeatures_forward,
                                     setdiff(idRegions(regionsFeatures_forward),'Control'),rep('Control',nRegions(regionsFeatures_forward)-1),
                                     statistics='quantile',probs=c(0.05,0.25,0.5,0.75,0.95),B=10000)
  save(result_multi_quantile,file=paste0('IPD_forward_max10000_',i,'_results_quantiles_5_25_50_75_95.RData'))
  
  # Plot test results
  pdf(paste0('IPD_forward_max10000_',i,'_results_quantiles_5_25_50_75_95.pdf'),width=8.5,height=10)
  plotTest_same_sample_size(result_multi_quantile,ask=FALSE,col=c(rep('blue',nRegions(result_multi_quantile)-1),'red'),xlab='Coordinates',
                            average=FALSE,size=TRUE,size_perc=TRUE,position=TRUE,cex.axis=1.5,cex.lab=1.5,cex.main=2,probs=c(0.05,0.25,0.5,0.75,0.95))
  dev.off()
  plotSummary(result_multi_quantile,groupby='feature',filenames=paste0('IPD_forward_max10000_',i,'_results_quantiles_5_25_50_75_95_summary.pdf'),
              only_significant=FALSE,gaps_tests=c(7,11,17,37),align_lab='Feature center',ask=FALSE,xlab='Coordinates',cellwidth=10,cellheight=10)
  
  # Test all features
  # number of control regions equal to feature regions
  result_multi_quantile=IWTomicsTest_same_sample_size(regionsFeatures_forward,
                                                      setdiff(idRegions(regionsFeatures_forward),'Control'),rep('Control',nRegions(regionsFeatures_forward)-1),
                                                      statistics='quantile',probs=c(0.05,0.25,0.5,0.75,0.95),B=10000)
  save(result_multi_quantile,file=paste0('IPD_forward_max10000_',i,'_results_same_sample_size_quantiles_5_25_50_75_95.RData'))
  
  # Plot test results
  pdf(paste0('IPD_forward_max10000_',i,'_results_same_sample_size_quantiles_5_25_50_75_95.pdf'),width=8.5,height=10)
  plotTest_same_sample_size(result_multi_quantile,ask=FALSE,col=c(rep('blue',nRegions(result_multi_quantile)-1),'red'),xlab='Coordinates',
                            average=FALSE,size=TRUE,size_perc=TRUE,position=TRUE,cex.axis=1.5,cex.lab=1.5,cex.main=2,probs=c(0.05,0.25,0.5,0.75,0.95))
  dev.off()
  plotSummary(result_multi_quantile,groupby='feature',filenames=paste0('IPD_forward_max10000_',i,'_results_same_sample_size_quantiles_5_25_50_75_95_summary.pdf'),
              only_significant=FALSE,gaps_tests=c(7,11,17,37),align_lab='Feature center',ask=FALSE,xlab='Coordinates',cellwidth=10,cellheight=10)
}


# Plot only reproducible results
alpha=0.05
min_significant=1
load('IPD_forward_max10000_1_results_quantiles_5_25_50_75_95.RData')
reproducible_pval_matrix=lapply(seq_len(nTests(result_multi_quantile)),function(name) matrix(TRUE,100,100))
for(i in 1:10){
  load(paste0('IPD_forward_max10000_',i,'_results_quantiles_5_25_50_75_95.RData'))
  reproducible_pval_matrix=mapply(function(reproducible_pval,result){
    reproducible_pval=(reproducible_pval & (result$IPD_Forward$adjusted_pval_matrix<=alpha))
  },reproducible_pval_matrix,result_multi_quantile@test$result,SIMPLIFY=FALSE)
  load(paste0('IPD_forward_max10000_',i,'_results_same_sample_size_quantiles_5_25_50_75_95.RData'))
  reproducible_pval_matrix=mapply(function(reproducible_pval,result){
    reproducible_pval=(reproducible_pval & (result$IPD_Forward$adjusted_pval_matrix<=alpha))
  },reproducible_pval_matrix,result_multi_quantile@test$result,SIMPLIFY=FALSE)
}
reproducible_pval_matrix=mapply(function(reproducible_pval){
  reproducible_pval[rowSums(reproducible_pval)<min_significant,]=matrix(FALSE,nrow=sum(rowSums(reproducible_pval)<min_significant),ncol=100)
  return(reproducible_pval)
},reproducible_pval_matrix,SIMPLIFY=FALSE)
save(reproducible_pval_matrix,file='IPD_forward_max10000_reproducible_results_quantiles_5_25_50_75_95.RData')

load('IPD_forward_max10000_1_results_quantiles_5_25_50_75_95.RData')
plotSummary_reproducible(result_multi_quantile,groupby='feature',filenames=paste0('IPD_forward_max10000_reproducible_results_quantiles_5_25_50_75_95_summary.pdf'),
                         reproducible_pval_matrix=reproducible_pval_matrix,only_significant=TRUE,
                         gaps_tests=c(7,11,17,37),align_lab='Feature center',xlab='',ylab='',ask=FALSE,cellwidth=5,cellheight=12,fontsize=12)
load('IPD_forward.RData')
pdf('boxplot_vs_control_forward_reproducible_results_quantiles_5_25_50_75_95.pdf',width=9,height=6)
for(i in seq_len(nrow(regionsFeatures_forward)-1)){
  plot_no_size_control(regionsFeatures_forward,reproducible_pval=reproducible_pval_matrix[[i]],
                       type='boxplot',probs=c(0.05,0.25,0.5,0.75,0.95),size_perc=TRUE,position=TRUE,average=FALSE,cex=1.4,cex.axis=1.4,cex.lab=1.4,cex.main=2.5,
                       id_regions_subset=c('Control',idRegions(regionsFeatures_forward)[i]),col=c('red','blue'),ask=FALSE,xlab='Coordinates')
}
dev.off()




### G quadruplex datasets ###

### Quantile test statistics (quantiles 0.05, 0.25, 0.50, 0.75 and 0.95) ###

# Test for difference using IWT
# 10 random subsample of maximum 10.000 regions for each region dataset
# see load_IPD_data.r
for(i in 1:10){
  load(paste0('IPD_forward_Gquad_max10000_',i,'.RData'))
  
  # Test all features
  result_multi_quantile=IWTomicsTest(regionsFeatures_forward_GQuadruplexMotifs_subset,
                                     setdiff(idRegions(regionsFeatures_forward_GQuadruplexMotifs_subset),'Control'),
                                     rep('Control',nRegions(regionsFeatures_forward_GQuadruplexMotifs_subset)-1),
                                     statistics='quantile',probs=c(0.05,0.25,0.5,0.75,0.95),B=10000)
  save(result_multi_quantile,file=paste0('IPD_forward_Gquad_max10000_',i,'_results_quantiles_5_25_50_75_95.RData'))
  
  # Plot test results
  pdf(paste0('IPD_forward_Gquad_max10000_',i,'_results_quantiles_5_25_50_75_95.pdf'),width=8.5,height=10)
  plotTest_same_sample_size(result_multi_quantile,ask=FALSE,col=c(rep('blue',nRegions(result_multi_quantile)-1),'red'),xlab='Coordinates',
                            average=FALSE,size=TRUE,size_perc=TRUE,position=TRUE,cex.axis=1.5,cex.lab=1.5,cex.main=2,probs=c(0.05,0.25,0.5,0.75,0.95))
  dev.off()
  
  # Test all features
  # number of control regions equal to feature regions
  result_multi_quantile=IWTomicsTest_same_sample_size(regionsFeatures_forward_GQuadruplexMotifs_subset,
                                                      setdiff(idRegions(regionsFeatures_forward_GQuadruplexMotifs_subset),'Control'),
                                                      rep('Control',nRegions(regionsFeatures_forward_GQuadruplexMotifs_subset)-1),
                                                      statistics='quantile',probs=c(0.05,0.25,0.5,0.75,0.95),B=10000)
  save(result_multi_quantile,file=paste0('IPD_forward_Gquad_max10000_',i,'_results_same_sample_size_quantiles_5_25_50_75_95.RData'))
  
  # Plot test results
  pdf(paste0('IPD_forward_Gquad_max10000_',i,'_results_same_sample_size_quantiles_5_25_50_75_95.pdf'),width=8.5,height=10)
  plotTest_same_sample_size(result_multi_quantile,ask=FALSE,col=c(rep('blue',nRegions(result_multi_quantile)-1),'red'),xlab='Coordinates',
                            average=FALSE,size=TRUE,size_perc=TRUE,position=TRUE,cex.axis=1.5,cex.lab=1.5,cex.main=2,probs=c(0.05,0.25,0.5,0.75,0.95))
  dev.off()
}


# Plot only reproducible results
alpha=0.05
min_significant=1
load('IPD_forward_Gquad_max10000_1_results_quantiles_5_25_50_75_95.RData')
reproducible_pval_matrix=lapply(seq_len(nTests(result_multi_quantile)),function(name) matrix(TRUE,100,100))
for(i in 1:10){
  load(paste0('IPD_forward_Gquad_max10000_',i,'_results_quantiles_5_25_50_75_95.RData'))
  reproducible_pval_matrix=mapply(function(reproducible_pval,result){
    reproducible_pval=(reproducible_pval & (result$IPD_Forward$adjusted_pval_matrix<=alpha))
  },reproducible_pval_matrix,result_multi_quantile@test$result,SIMPLIFY=FALSE)
  load(paste0('IPD_forward_Gquad_max10000_',i,'_results_same_sample_size_quantiles_5_25_50_75_95.RData'))
  reproducible_pval_matrix=mapply(function(reproducible_pval,result){
    reproducible_pval=(reproducible_pval & (result$IPD_Forward$adjusted_pval_matrix<=alpha))
  },reproducible_pval_matrix,result_multi_quantile@test$result,SIMPLIFY=FALSE)
}
reproducible_pval_matrix=mapply(function(reproducible_pval){
  reproducible_pval[rowSums(reproducible_pval)<min_significant,]=matrix(FALSE,nrow=sum(rowSums(reproducible_pval)<min_significant),ncol=100)
  return(reproducible_pval)
},reproducible_pval_matrix,SIMPLIFY=FALSE)
save(reproducible_pval_matrix,file='IPD_forward_Gquad_max10000_reproducible_results_quantiles_5_25_50_75_95.RData')

load('IPD_forward_Gquad.RData')
pdf('boxplot_vs_control_forward_Gquad_reproducible_results_quantiles_5_25_50_75_95.pdf',width=9,height=6)
for(i in seq_len(nrow(regionsFeatures_forward_GQuadruplexMotifs_subset)-1)){
  plot_no_size_control(regionsFeatures_forward_GQuadruplexMotifs_subset,reproducible_pval=reproducible_pval_matrix[[i]],
                       type='boxplot',probs=c(0.05,0.25,0.5,0.75,0.95),size_perc=TRUE,position=TRUE,average=FALSE,cex=1.4,cex.axis=1.4,cex.lab=1.4,cex.main=2.5,
                       id_regions_subset=c('Control',idRegions(regionsFeatures_forward_GQuadruplexMotifs_subset)[i]),col=c('red','blue'),ask=FALSE,xlab='Coordinates')
}
dev.off()




### A-phased repeats datasets ###

### Quantile test statistics (quantiles 0.05, 0.25, 0.50, 0.75 and 0.95) ###

# Test for difference using IWT
# 10 random subsample of maximum 10.000 regions for each region dataset
# see load_IPD_data.r
for(i in 1:10){
  load(paste0('IPD_forward_Aphased_max10000_',i,'.RData'))
  
  # Test all features
  result_multi_quantile=IWTomicsTest(regionsFeatures_forward_APhasedRepeats_subset,
                                     setdiff(idRegions(regionsFeatures_forward_APhasedRepeats_subset),'Control'),
                                     rep('Control',nRegions(regionsFeatures_forward_APhasedRepeats_subset)-1),
                                     statistics='quantile',probs=c(0.05,0.25,0.5,0.75,0.95),B=10000)
  save(result_multi_quantile,file=paste0('IPD_forward_Aphased_max10000_',i,'_results_quantiles_5_25_50_75_95.RData'))
  
  # Plot test results
  pdf(paste0('IPD_forward_Aphased_max10000_',i,'_results_quantiles_5_25_50_75_95.pdf'),width=8.5,height=10)
  plotTest_same_sample_size(result_multi_quantile,ask=FALSE,col=c(rep('blue',nRegions(result_multi_quantile)-1),'red'),xlab='Coordinates',
                            average=FALSE,size=TRUE,size_perc=TRUE,position=TRUE,cex.axis=1.5,cex.lab=1.5,cex.main=2,probs=c(0.05,0.25,0.5,0.75,0.95))
  dev.off()
  
  # Test all features
  # number of control regions equal to feature regions
  result_multi_quantile=IWTomicsTest_same_sample_size(regionsFeatures_forward_APhasedRepeats_subset,
                                                      setdiff(idRegions(regionsFeatures_forward_APhasedRepeats_subset),'Control'),
                                                      rep('Control',nRegions(regionsFeatures_forward_APhasedRepeats_subset)-1),
                                                      statistics='quantile',probs=c(0.05,0.25,0.5,0.75,0.95),B=10000)
  save(result_multi_quantile,file=paste0('IPD_forward_Aphased_max10000_',i,'_results_same_sample_size_quantiles_5_25_50_75_95.RData'))
  
  # Plot test results
  pdf(paste0('IPD_forward_Aphased_max10000_',i,'_results_same_sample_size_quantiles_5_25_50_75_95.pdf'),width=8.5,height=10)
  plotTest_same_sample_size(result_multi_quantile,ask=FALSE,col=c(rep('blue',nRegions(result_multi_quantile)-1),'red'),xlab='Coordinates',
                            average=FALSE,size=TRUE,size_perc=TRUE,position=TRUE,cex.axis=1.5,cex.lab=1.5,cex.main=2,probs=c(0.05,0.25,0.5,0.75,0.95))
  dev.off()
}


# Plot only reproducible results
alpha=0.05
min_significant=1
load('IPD_forward_Aphased_max10000_1_results_quantiles_5_25_50_75_95.RData')
reproducible_pval_matrix=lapply(seq_len(nTests(result_multi_quantile)),function(name) matrix(TRUE,100,100))
for(i in 1:10){
  load(paste0('IPD_forward_Aphased_max10000_',i,'_results_quantiles_5_25_50_75_95.RData'))
  reproducible_pval_matrix=mapply(function(reproducible_pval,result){
    reproducible_pval=(reproducible_pval & (result$IPD_Forward$adjusted_pval_matrix<=alpha))
  },reproducible_pval_matrix,result_multi_quantile@test$result,SIMPLIFY=FALSE)
  load(paste0('IPD_forward_Aphased_max10000_',i,'_results_same_sample_size_quantiles_5_25_50_75_95.RData'))
  reproducible_pval_matrix=mapply(function(reproducible_pval,result){
    reproducible_pval=(reproducible_pval & (result$IPD_Forward$adjusted_pval_matrix<=alpha))
  },reproducible_pval_matrix,result_multi_quantile@test$result,SIMPLIFY=FALSE)
}
reproducible_pval_matrix=mapply(function(reproducible_pval){
  reproducible_pval[rowSums(reproducible_pval)<min_significant,]=matrix(FALSE,nrow=sum(rowSums(reproducible_pval)<min_significant),ncol=100)
  return(reproducible_pval)
},reproducible_pval_matrix,SIMPLIFY=FALSE)
save(reproducible_pval_matrix,file='IPD_forward_Aphased_max10000_reproducible_results_quantiles_5_25_50_75_95.RData')

load('IPD_forward_Aphased_max10000_reproducible_results_quantiles_5_25_50_75_95.RData')
load('IPD_forward_Aphased.RData')
pdf('boxplot_vs_control_forward_Aphased_reproducible_results_quantiles_5_25_50_75_95.pdf',width=9,height=6)
for(i in seq_len(nrow(regionsFeatures_forward_APhasedRepeats_subset)-1)){
  plot_no_size_control(regionsFeatures_forward_APhasedRepeats_subset,reproducible_pval=reproducible_pval_matrix[[i]],
                       type='boxplot',probs=c(0.05,0.25,0.5,0.75,0.95),size_perc=TRUE,position=TRUE,average=FALSE,cex=1.4,cex.axis=1.4,cex.lab=1.4,cex.main=2.5,
                       id_regions_subset=c('Control',idRegions(regionsFeatures_forward_APhasedRepeats_subset)[i]),col=c('red','blue'),ask=FALSE,xlab='Coordinates')
}
dev.off()















#######################
##### MINUS STRAND ####
#######################

### Complete datasets ###

### Mean test statistics ###

# Test for difference using IWT
# 10 random subsample of maximum 10.000 regions for each region dataset
# see load_IPD_data.r
for(i in 1:10){
  load(paste0('IPD_reverse_max10000_',i,'.RData'))
  
  # Test all features
  result_mean=IWTomicsTest(regionsFeatures_reverse,
                           setdiff(idRegions(regionsFeatures_reverse),'Control'),rep('Control',nRegions(regionsFeatures_reverse)-1),
                           statistics='mean',B=10000)
  save(result_mean,file=paste0('IPD_reverse_max10000_',i,'_results_mean.RData'))
  
  # Plot test results
  pdf(paste0('IPD_reverse_max10000_',i,'_results_mean.pdf'),width=8.5,height=10)
  plotTest_same_sample_size(result_mean,ask=FALSE,col=c(rep('blue',nRegions(result_mean)-1),'red'),xlab='Coordinates',
                            average=FALSE,size=TRUE,size_perc=TRUE,position=TRUE,cex.axis=1.5,cex.lab=1.5,cex.main=2,probs=c(0.05,0.25,0.5,0.75,0.95))
  dev.off()
  plotSummary(result_mean,groupby='feature',filenames=paste0('IPD_reverse_max10000_',i,'_results_mean_summary.pdf'),
              only_significant=FALSE,gaps_tests=c(7,11,17,37),align_lab='Feature center',ask=FALSE,xlab='Coordinates',cellwidth=10,cellheight=10)
  
  # Test all features
  # number of control regions equal to feature regions
  result_mean=IWTomicsTest_same_sample_size(regionsFeatures_reverse,
                                            setdiff(idRegions(regionsFeatures_reverse),'Control'),rep('Control',nRegions(regionsFeatures_reverse)-1),
                                            statistics='mean',B=10000)
  save(result_mean,file=paste0('IPD_reverse_max10000_',i,'_results_same_sample_size_mean.RData'))
  
  # Plot test results
  pdf(paste0('IPD_reverse_max10000_',i,'_results_same_sample_size_mean.pdf'),width=8.5,height=10)
  plotTest_same_sample_size(result_mean,ask=FALSE,col=c(rep('blue',nRegions(result_mean)-1),'red'),xlab='Coordinates',
                            average=FALSE,size=TRUE,size_perc=TRUE,position=TRUE,cex.axis=1.5,cex.lab=1.5,cex.main=2,probs=c(0.05,0.25,0.5,0.75,0.95))
  dev.off()
  plotSummary(result_mean,groupby='feature',filenames=paste0('IPD_reverse_max10000_',i,'_results_same_sample_size_mean_summary.pdf'),
              only_significant=FALSE,gaps_tests=c(7,11,17,37),align_lab='Feature center',ask=FALSE,xlab='Coordinates',cellwidth=10,cellheight=10)
}


# Plot only reproducible results
alpha=0.05
min_significant=1
load('IPD_reverse_max10000_1_results_mean.RData')
reproducible_pval_matrix=lapply(seq_len(nTests(result_mean)),function(name) matrix(TRUE,100,100))
for(i in 1:10){
  load(paste0('IPD_reverse_max10000_',i,'_results_mean.RData'))
  reproducible_pval_matrix=mapply(function(reproducible_pval,result){
    reproducible_pval=(reproducible_pval & (result$IPD_Reverse$adjusted_pval_matrix<=alpha))
  },reproducible_pval_matrix,result_mean@test$result,SIMPLIFY=FALSE)
  load(paste0('IPD_reverse_max10000_',i,'_results_same_sample_size_mean.RData'))
  reproducible_pval_matrix=mapply(function(reproducible_pval,result){
    reproducible_pval=(reproducible_pval & (result$IPD_Reverse$adjusted_pval_matrix<=alpha))
  },reproducible_pval_matrix,result_mean@test$result,SIMPLIFY=FALSE)
}
reproducible_pval_matrix=mapply(function(reproducible_pval){
  reproducible_pval[rowSums(reproducible_pval)<min_significant,]=matrix(FALSE,nrow=sum(rowSums(reproducible_pval)<min_significant),ncol=100)
  return(reproducible_pval)
},reproducible_pval_matrix,SIMPLIFY=FALSE)
save(reproducible_pval_matrix,file='IPD_reverse_max10000_reproducible_results_mean.RData')

load('IPD_reverse_max10000_1_results_mean.RData')
plotSummary_reproducible(result_mean,groupby='feature',filenames=paste0('IPD_reverse_max10000_reproducible_results_mean_summary.pdf'),
                         reproducible_pval_matrix=reproducible_pval_matrix,only_significant=TRUE,
                         gaps_tests=c(7,11,17,37),align_lab='Feature center',xlab='',ylab='',ask=FALSE,cellwidth=5,cellheight=12,fontsize=12)
load('IPD_reverse.RData')
pdf('boxplot_vs_control_reverse_reproducible_results_mean.pdf',width=9,height=6)
for(i in seq_len(nrow(regionsFeatures_reverse)-1)){
  plot_no_size_control(regionsFeatures_reverse,reproducible_pval=reproducible_pval_matrix[[i]],
                       type='boxplot',probs=c(0.05,0.25,0.5,0.75,0.95),size_perc=TRUE,position=TRUE,average=FALSE,cex=1.4,cex.axis=1.4,cex.lab=1.4,cex.main=2.5,
                       id_regions_subset=c('Control',idRegions(regionsFeatures_reverse)[i]),col=c('red','blue'),ask=FALSE,xlab='Coordinates')
}
dev.off()





### Median test statistics ###

# Test for difference using IWT
# 10 random subsample of maximum 10.000 regions for each region dataset
# see load_IPD_data.r
for(i in 1:10){
  load(paste0('IPD_reverse_max10000_',i,'.RData'))
  
  # Test all features
  result_median=IWTomicsTest(regionsFeatures_reverse,
                             setdiff(idRegions(regionsFeatures_reverse),'Control'),rep('Control',nRegions(regionsFeatures_reverse)-1),
                             statistics='median',B=10000)
  save(result_median,file=paste0('IPD_reverse_max10000_',i,'_results_median.RData'))
  
  # Plot test results
  pdf(paste0('IPD_reverse_max10000_',i,'_results_median.pdf'),width=8.5,height=10)
  plotTest_same_sample_size(result_median,ask=FALSE,col=c(rep('blue',nRegions(result_median)-1),'red'),xlab='Coordinates',
                            average=FALSE,size=TRUE,size_perc=TRUE,position=TRUE,cex.axis=1.5,cex.lab=1.5,cex.main=2,probs=c(0.05,0.25,0.5,0.75,0.95))
  dev.off()
  plotSummary(result_median,groupby='feature',filenames=paste0('IPD_reverse_max10000_',i,'_results_median_summary.pdf'),
              only_significant=FALSE,gaps_tests=c(7,11,17,37),align_lab='Feature center',ask=FALSE,xlab='Coordinates',cellwidth=10,cellheight=10)
  
  # Test all features
  # number of control regions equal to feature regions
  result_median=IWTomicsTest_same_sample_size(regionsFeatures_reverse,
                                              setdiff(idRegions(regionsFeatures_reverse),'Control'),rep('Control',nRegions(regionsFeatures_reverse)-1),
                                              statistics='median',B=10000)
  save(result_median,file=paste0('IPD_reverse_max10000_',i,'_results_same_sample_size_median.RData'))
  
  # Plot test results
  pdf(paste0('IPD_reverse_max10000_',i,'_results_same_sample_size_median.pdf'),width=8.5,height=10)
  plotTest_same_sample_size(result_median,ask=FALSE,col=c(rep('blue',nRegions(result_median)-1),'red'),xlab='Coordinates',
                            average=FALSE,size=TRUE,size_perc=TRUE,position=TRUE,cex.axis=1.5,cex.lab=1.5,cex.main=2,probs=c(0.05,0.25,0.5,0.75,0.95))
  dev.off()
  plotSummary(result_median,groupby='feature',filenames=paste0('IPD_reverse_max10000_',i,'_results_same_sample_size_median_summary.pdf'),
              only_significant=FALSE,gaps_tests=c(7,11,17,37),align_lab='Feature center',ask=FALSE,xlab='Coordinates',cellwidth=10,cellheight=10)
}


# Plot only reproducible results
alpha=0.05
min_significant=1
load('IPD_reverse_max10000_1_results_median.RData')
reproducible_pval_matrix=lapply(seq_len(nTests(result_median)),function(name) matrix(TRUE,100,100))
for(i in 1:10){
  load(paste0('IPD_reverse_max10000_',i,'_results_median.RData'))
  reproducible_pval_matrix=mapply(function(reproducible_pval,result){
    reproducible_pval=(reproducible_pval & (result$IPD_Reverse$adjusted_pval_matrix<=alpha))
  },reproducible_pval_matrix,result_median@test$result,SIMPLIFY=FALSE)
  load(paste0('IPD_reverse_max10000_',i,'_results_same_sample_size_median.RData'))
  reproducible_pval_matrix=mapply(function(reproducible_pval,result){
    reproducible_pval=(reproducible_pval & (result$IPD_Reverse$adjusted_pval_matrix<=alpha))
  },reproducible_pval_matrix,result_median@test$result,SIMPLIFY=FALSE)
}
reproducible_pval_matrix=mapply(function(reproducible_pval){
  reproducible_pval[rowSums(reproducible_pval)<min_significant,]=matrix(FALSE,nrow=sum(rowSums(reproducible_pval)<min_significant),ncol=100)
  return(reproducible_pval)
},reproducible_pval_matrix,SIMPLIFY=FALSE)
save(reproducible_pval_matrix,file='IPD_reverse_max10000_reproducible_results_median.RData')

load('IPD_reverse_max10000_1_results_median.RData')
plotSummary_reproducible(result_median,groupby='feature',filenames=paste0('IPD_reverse_max10000_reproducible_results_median_summary.pdf'),
                         reproducible_pval_matrix=reproducible_pval_matrix,only_significant=TRUE,
                         gaps_tests=c(7,11,17,37),align_lab='Feature center',xlab='',ylab='',ask=FALSE,cellwidth=5,cellheight=12,fontsize=12)
load('IPD_reverse.RData')
pdf('boxplot_vs_control_reverse_reproducible_results_median.pdf',width=9,height=6)
for(i in seq_len(nrow(regionsFeatures_reverse)-1)){
  plot_no_size_control(regionsFeatures_reverse,reproducible_pval=reproducible_pval_matrix[[i]],
                       type='boxplot',probs=c(0.05,0.25,0.5,0.75,0.95),size_perc=TRUE,position=TRUE,average=FALSE,cex=1.4,cex.axis=1.4,cex.lab=1.4,cex.main=2.5,
                       id_regions_subset=c('Control',idRegions(regionsFeatures_reverse)[i]),col=c('red','blue'),ask=FALSE,xlab='Coordinates')
}
dev.off()





### Quantile test statistics (quantiles 0.05, 0.25, 0.50, 0.75 and 0.95) ###

# Test for difference using IWT
# 10 random subsample of maximum 10.000 regions for each region dataset
# see load_IPD_data.r
for(i in 1:10){
  load(paste0('IPD_reverse_max10000_',i,'.RData'))
  
  # Test all features
  result_multi_quantile=IWTomicsTest(regionsFeatures_reverse,
                                     setdiff(idRegions(regionsFeatures_reverse),'Control'),rep('Control',nRegions(regionsFeatures_reverse)-1),
                                     statistics='quantile',probs=c(0.05,0.25,0.5,0.75,0.95),B=10000)
  save(result_multi_quantile,file=paste0('IPD_reverse_max10000_',i,'_results_quantiles_5_25_50_75_95.RData'))
  
  # Plot test results
  pdf(paste0('IPD_reverse_max10000_',i,'_results_quantiles_5_25_50_75_95.pdf'),width=8.5,height=10)
  plotTest_same_sample_size(result_multi_quantile,ask=FALSE,col=c(rep('blue',nRegions(result_multi_quantile)-1),'red'),xlab='Coordinates',
                            average=FALSE,size=TRUE,size_perc=TRUE,position=TRUE,cex.axis=1.5,cex.lab=1.5,cex.main=2,probs=c(0.05,0.25,0.5,0.75,0.95))
  dev.off()
  plotSummary(result_multi_quantile,groupby='feature',filenames=paste0('IPD_reverse_max10000_',i,'_results_quantiles_5_25_50_75_95_summary.pdf'),
              only_significant=FALSE,gaps_tests=c(7,11,17,37),align_lab='Feature center',ask=FALSE,xlab='Coordinates',cellwidth=10,cellheight=10)
  
  # Test all features
  # number of control regions equal to feature regions
  result_multi_quantile=IWTomicsTest_same_sample_size(regionsFeatures_reverse,
                                                      setdiff(idRegions(regionsFeatures_reverse),'Control'),rep('Control',nRegions(regionsFeatures_reverse)-1),
                                                      statistics='quantile',probs=c(0.05,0.25,0.5,0.75,0.95),B=10000)
  save(result_multi_quantile,file=paste0('IPD_reverse_max10000_',i,'_results_same_sample_size_quantiles_5_25_50_75_95.RData'))
  
  # Plot test results
  pdf(paste0('IPD_reverse_max10000_',i,'_results_same_sample_size_quantiles_5_25_50_75_95.pdf'),width=8.5,height=10)
  plotTest_same_sample_size(result_multi_quantile,ask=FALSE,col=c(rep('blue',nRegions(result_multi_quantile)-1),'red'),xlab='Coordinates',
                            average=FALSE,size=TRUE,size_perc=TRUE,position=TRUE,cex.axis=1.5,cex.lab=1.5,cex.main=2,probs=c(0.05,0.25,0.5,0.75,0.95))
  dev.off()
  plotSummary(result_multi_quantile,groupby='feature',filenames=paste0('IPD_reverse_max10000_',i,'_results_same_sample_size_quantiles_5_25_50_75_95_summary.pdf'),
              only_significant=FALSE,gaps_tests=c(7,11,17,37),align_lab='Feature center',ask=FALSE,xlab='Coordinates',cellwidth=10,cellheight=10)
}


# Plot only reproducible results
alpha=0.05
min_significant=1
load('IPD_reverse_max10000_1_results_quantiles_5_25_50_75_95.RData')
reproducible_pval_matrix=lapply(seq_len(nTests(result_multi_quantile)),function(name) matrix(TRUE,100,100))
for(i in 1:10){
  load(paste0('IPD_reverse_max10000_',i,'_results_quantiles_5_25_50_75_95.RData'))
  reproducible_pval_matrix=mapply(function(reproducible_pval,result){
    reproducible_pval=(reproducible_pval & (result$IPD_Reverse$adjusted_pval_matrix<=alpha))
  },reproducible_pval_matrix,result_multi_quantile@test$result,SIMPLIFY=FALSE)
  load(paste0('IPD_reverse_max10000_',i,'_results_same_sample_size_quantiles_5_25_50_75_95.RData'))
  reproducible_pval_matrix=mapply(function(reproducible_pval,result){
    reproducible_pval=(reproducible_pval & (result$IPD_Reverse$adjusted_pval_matrix<=alpha))
  },reproducible_pval_matrix,result_multi_quantile@test$result,SIMPLIFY=FALSE)
}
reproducible_pval_matrix=mapply(function(reproducible_pval){
  reproducible_pval[rowSums(reproducible_pval)<min_significant,]=matrix(FALSE,nrow=sum(rowSums(reproducible_pval)<min_significant),ncol=100)
  return(reproducible_pval)
},reproducible_pval_matrix,SIMPLIFY=FALSE)
save(reproducible_pval_matrix,file='IPD_reverse_max10000_reproducible_results_quantiles_5_25_50_75_95.RData')

load('IPD_reverse_max10000_1_results_quantiles_5_25_50_75_95.RData')
plotSummary_reproducible(result_multi_quantile,groupby='feature',filenames=paste0('IPD_reverse_max10000_reproducible_results_quantiles_5_25_50_75_95_summary.pdf'),
                         reproducible_pval_matrix=reproducible_pval_matrix,only_significant=TRUE,
                         gaps_tests=c(7,11,17,37),align_lab='Feature center',xlab='',ylab='',ask=FALSE,cellwidth=5,cellheight=12,fontsize=12)
load('IPD_reverse.RData')
pdf('boxplot_vs_control_reverse_reproducible_results_quantiles_5_25_50_75_95.pdf',width=9,height=6)
for(i in seq_len(nrow(regionsFeatures_reverse)-1)){
  plot_no_size_control(regionsFeatures_reverse,reproducible_pval=reproducible_pval_matrix[[i]],
                       type='boxplot',probs=c(0.05,0.25,0.5,0.75,0.95),size_perc=TRUE,position=TRUE,average=FALSE,cex=1.4,cex.axis=1.4,cex.lab=1.4,cex.main=2.5,
                       id_regions_subset=c('Control',idRegions(regionsFeatures_reverse)[i]),col=c('red','blue'),ask=FALSE,xlab='Coordinates')
}
dev.off()




