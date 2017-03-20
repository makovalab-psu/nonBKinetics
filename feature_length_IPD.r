require(IWTomics)
source('IWTomics_modified_functions.r')



#######################
##### PLUS STRAND #####
#######################

load('IPD_forward.RData')

# Feature length histograms
length=lapply(regionsFeatures_forward@regions[-84],function(region) region$length)
length100=lapply(length,pmin,100) # Set maximum feature length to 100 (window length)
index=lapply(length100,order)
length_100_range=length100
pdf("length_distribution_forward.pdf",width=7,height=5)
for(i in seq_len(nrow(regionsFeatures_forward)-1)){
  counts=hist(length100[[i]],breaks=0:100+0.5,xlim=c(0.5,100.5),col="blue",xlab="Length",main=nameRegions(regionsFeatures_forward)[i])$counts
  length_100_range[[i]]=which(counts>=10)
}
dev.off()

# Boxplot of log mean IPD in windows, grouped by feature lengths
meanIPD=lapply(regionsFeatures_forward@features$IPD_Forward[-84],function(feature) colMeans(feature,na.rm=TRUE))
pdf("length_vs_meanIPD_log_forward_boxplot.pdf",width=7,height=7)
for(i in seq_len(nrow(regionsFeatures_forward)-1)){
  box=boxplot(log(meanIPD[[i]][index[[i]]]+0.01)~factor(length100[[i]][index[[i]]]),plot=FALSE)
  box$stats=as.matrix(Reduce(cbind,tapply(log(meanIPD[[i]][index[[i]]]+0.01),factor(length100[[i]][index[[i]]]),quantile,probs=c(0.05,0.25,0.5,0.75,0.95),na.rm=TRUE,SIMPLIFY=FALSE)))
  perc=round(box$n/max(box$n)*100000,0)
  col_perc=rev(heat.colors(100000))[perc]
  bxp(box,at=sort(unique(length100[[i]])),xlim=c(0,100),outline=FALSE,pars=list(whisklty=3,staplewex=0.8,boxfill=col_perc),
      xlab="Length",ylab="Log Mean IPD",main=nameRegions(regionsFeatures_forward)[i])
}
dev.off()

# Boxplot of log mean IPD in features, grouped by feature lengths
meanIPD_feature=mapply(function(feature,region){
  pos=apply(as.matrix(mcols(region)),1,
            function(info){
              info=pmax(pmin(info[1:2],100),0)
              info[2]=min((100-info[1]),info[2])
              return(c(rep(NA,info[2]),rep(1,info[1]),rep(NA,100-sum(info))))})
  feature=feature*pos
  meanIPD=colMeans(feature,na.rm=TRUE)
  meanIPD[is.nan(meanIPD)]=NA
  return(meanIPD)
},
regionsFeatures_forward@features$IPD_Forward[-84],regionsFeatures_forward@regions[-84],SIMPLIFY=FALSE)
pdf("length_vs_meanIPD_log_feature_forward_boxplot.pdf",width=7,height=7)
for(i in seq_len(nrow(regionsFeatures_forward)-1)){
  box=boxplot(log(meanIPD_feature[[i]][index[[i]]]+0.01)~factor(length100[[i]][index[[i]]]),plot=FALSE)
  box$stats=as.matrix(Reduce(cbind,tapply(log(meanIPD_feature[[i]][index[[i]]]+0.01),factor(length100[[i]][index[[i]]]),quantile,probs=c(0.05,0.25,0.5,0.75,0.95),na.rm=TRUE,SIMPLIFY=FALSE)))
  perc=round(box$n/max(box$n)*100000,0)
  col_perc=rev(heat.colors(100000))[perc]
  bxp(box,at=sort(unique(length100[[i]])),xlim=c(0,100),outline=FALSE,pars=list(whisklty=3,staplewex=0.8,boxfill=col_perc),
      xlab="Length",ylab="Log Mean IPD",main=nameRegions(regionsFeatures_forward)[i])
}
dev.off()

# Pointwise boxplots, grouped by feature lengths
for(region in idRegions(regionsFeatures_forward)){
  if(length(length_100_range[[region]])>0){
    pdf(paste0('boxplot_vs_control_forward_',region,'.pdf'),width=9,height=6)
    plot_no_size_control(regionsFeatures_forward,lengths=length_100_range[[region]],
                         type='boxplot',probs=c(0.05,0.25,0.5,0.75,0.95),size_perc=TRUE,position=TRUE,average=FALSE,cex=1.4,cex.axis=1.4,cex.lab=1.4,cex.main=2.5,
                         id_regions_subset=c('Control',region),col=c('red','blue'),ask=FALSE,xlab='Coordinates')
    dev.off()
  }
}




### An datasets, grouped by feature lengths ###

# Create dataset
region='An'
lengths=length_100_range[[region]]
x=regionsFeatures_forward[c('Control',region),]
new_ids=paste0(region,lengths)
reg=x@regions[[region]]
feat=x@features[["IPD_Forward"]][[region]]
len=x@length_features[["IPD_Forward"]][[region]]
index=lapply(lengths,function(len) which(reg$length==len))
if(100 %in% lengths)
  index[[which(lengths==100)]]=which(reg$length>=100)
names(index)=new_ids
size_new=unlist(lapply(index,length))
feat_new=lapply(index,function(index) as.matrix(feat[,index]))
reg_new=lapply(index,function(index) reg[index])
len_new=lapply(index,function(index) len[index])
reg_datasets_new=as.data.frame(matrix(x@metadata$region_datasets[region,],nrow=length(lengths),ncol=3,byrow=TRUE),row.names=new_ids)
names(reg_datasets_new)=names(x@metadata$region_datasets)
reg_datasets_new$name=paste0(reg_datasets_new$name,' - Length ',lengths)
reg_datasets_new$size=size_new
x@metadata$region_datasets=rbind(x@metadata$region_datasets,reg_datasets_new)
x@regions=c(x@regions,GRangesList(reg_new))
mcols(x@regions[new_ids])=mcols(x@regions[region])
for(new_id in new_ids){
  NAbefore=pmax(unique(mcols(x@regions[[new_id]])$NAbefore),0)
  len=unique(pmin(unique(mcols(x@regions[[new_id]])$length),100))
  if(length(len)>1)
    stop('Something wrong with the length.')
  mcols(x@regions[new_id])[seq_len(min(NAbefore))]=NA
  mcols(x@regions[new_id])[100-seq_len(100-max(NAbefore)-len)+1]=NA
}
x@metadata$feature_datasets=x@metadata$feature_datasets[,c('name',paste0('file_',c(c('Control',region),rep(region,length(lengths)))),'resolution')]
names(x@metadata$feature_datasets)=c('name',paste0('file_',c('Control',region,new_ids)),'resolution')
x@features[["IPD_Forward"]]=c(x@features[["IPD_Forward"]],feat_new)
x@length_features[["IPD_Forward"]]=c(x@length_features[["IPD_Forward"]],len_new)
regionsFeatures_forward=x[setdiff(idRegions(x),region),]
save(regionsFeatures_forward,file='IPD_forward_An_different_lengths.RData')

# Subsample regions, maximum 10.000 for each region dataset
# 10 random subsamples
# to test for differences using IWT
maxsize=10000
regions_maxsize=rownames(regionsFeatures_forward@metadata$region_datasets)[regionsFeatures_forward@metadata$region_datasets$size>maxsize]
regionsFeatures_forward_complete=regionsFeatures_forward
# 10 times
for(i in 1:10){
  index=lapply(regions_maxsize,function(id) sort(sample(regionsFeatures_forward_complete@metadata$region_datasets[id,'size'],maxsize)))
  IPD_filtered=mapply(function(IPD,index) as.matrix(IPD[,index]),
                      regionsFeatures_forward_complete@features$IPD_Forward[regions_maxsize],index,SIMPLIFY=FALSE)
  regions_filtered=mapply(function(regions,index) regions[index,],
                          regionsFeatures_forward_complete@regions[regions_maxsize],index,SIMPLIFY=FALSE)
  length_features_filtered=mapply(function(length,index) length[index],
                                  regionsFeatures_forward_complete@length_features$IPD_Forward[regions_maxsize],index,SIMPLIFY=FALSE)
  regionsFeatures_forward=regionsFeatures_forward_complete
  regionsFeatures_forward@metadata$region_datasets[regions_maxsize,'size']=maxsize
  regionsFeatures_forward@regions[regions_maxsize]=GRangesList(regions_filtered)
  regionsFeatures_forward@features$IPD_Forward[regions_maxsize]=IPD_filtered
  regionsFeatures_forward@length_features$IPD_Forward[regions_maxsize]=length_features_filtered
  validObject(regionsFeatures_forward)
  
  # Add microsatellite sequences
  mcols(regionsFeatures_forward@regions) <- mcols(regionsFeatures_forward_complete@regions)
  
  save(regionsFeatures_forward,file=paste0('IPD_forward_An_different_lengths_max10000_',i,'.RData'))
}

# Test for difference using IWT
# 10 random subsample of maximum 10.000 regions for each region dataset
for(i in 1:10){
  load(paste0('IPD_forward_An_different_lengths_max10000_',i,'.RData'))
  
  # Test all features
  result_multi_quantile=IWTomicsTest(regionsFeatures_forward,
                                     setdiff(idRegions(regionsFeatures_forward),'Control'),rep('Control',nRegions(regionsFeatures_forward)-1),
                                     statistics='quantile',probs=c(0.05,0.25,0.5,0.75,0.95),B=10000)
  save(result_multi_quantile,file=paste0('IPD_forward_An_different_lengths_max10000_',i,'_results_quantiles_5_25_50_75_95.RData'))

  # Plot test results
  pdf(paste0('IPD_forward_An_different_lengths_max10000_',i,'_results_quantiles_5_25_50_75_95.pdf'),width=8.5,height=10)
  plotTest_same_sample_size(result_multi_quantile,ask=FALSE,col=c(rep('blue',nRegions(result_multi_quantile)-1),'red'),xlab='Coordinates',
                            average=FALSE,size=TRUE,size_perc=TRUE,position=TRUE,cex.axis=1.5,cex.lab=1.5,cex.main=2,probs=c(0.05,0.25,0.5,0.75,0.95))
  dev.off()
  
  # Test all features
  # number of control regions equal to feature regions
  result_multi_quantile=IWTomicsTest_same_sample_size(regionsFeatures_forward,
                                                      setdiff(idRegions(regionsFeatures_forward),'Control'),rep('Control',nRegions(regionsFeatures_forward)-1),
                                                      statistics='quantile',probs=c(0.05,0.25,0.5,0.75,0.95),B=10000)
  save(result_multi_quantile,file=paste0('IPD_forward_An_different_lengths_max10000_',i,'_results_same_sample_size_quantiles_5_25_50_75_95.RData'))
  
  # Plot test results
  pdf(paste0('IPD_forward_An_different_lengths_max10000_',i,'_results_same_sample_size_quantiles_5_25_50_75_95.pdf'),width=8.5,height=10)
  plotTest_same_sample_size(result_multi_quantile,ask=FALSE,col=c(rep('blue',nRegions(result_multi_quantile)-1),'red'),xlab='Coordinates',
                            average=FALSE,size=TRUE,size_perc=TRUE,position=TRUE,cex.axis=1.5,cex.lab=1.5,cex.main=2,probs=c(0.05,0.25,0.5,0.75,0.95))
  dev.off()
}

# Plot only reproducible results
alpha=0.05
min_significant=1
load('IPD_forward_An_different_lengths_max10000_1_results_quantiles_5_25_50_75_95.RData')
reproducible_pval_matrix=lapply(seq_len(nTests(result_multi_quantile)),function(name) matrix(TRUE,100,100))
for(i in 1:10){
  load(paste0('IPD_forward_An_different_lengths_max10000_',i,'_results_quantiles_5_25_50_75_95.RData'))
  reproducible_pval_matrix=mapply(function(reproducible_pval,result){
    reproducible_pval=(reproducible_pval & (result$IPD_Forward$adjusted_pval_matrix<=alpha))
  },reproducible_pval_matrix,result_multi_quantile@test$result,SIMPLIFY=FALSE)
  load(paste0('IPD_forward_An_different_lengths_max10000_',i,'_results_same_sample_size_quantiles_5_25_50_75_95.RData'))
  reproducible_pval_matrix=mapply(function(reproducible_pval,result){
    reproducible_pval=(reproducible_pval & (result$IPD_Forward$adjusted_pval_matrix<=alpha))
  },reproducible_pval_matrix,result_multi_quantile@test$result,SIMPLIFY=FALSE)
}
reproducible_pval_matrix=mapply(function(reproducible_pval){
  reproducible_pval[rowSums(reproducible_pval)<min_significant,]=matrix(FALSE,nrow=sum(rowSums(reproducible_pval)<min_significant),ncol=100)
  return(reproducible_pval)
},reproducible_pval_matrix,SIMPLIFY=FALSE)
save(reproducible_pval_matrix,file='IPD_forward_An_different_lengths_max10000_reproducible_results_quantiles_5_25_50_75_95.RData')

load('IPD_forward_An_different_lengths.RData')
pdf('boxplot_vs_control_forward_An_different_lengths_reproducible_results_quantiles_5_25_50_75_95.pdf',width=9,height=6)
for(i in seq_len(nrow(regionsFeatures_forward)-1)){
  plot_no_size_control(regionsFeatures_forward,reproducible_pval=reproducible_pval_matrix[[i]],
                       type='boxplot',probs=c(0.05,0.25,0.5,0.75,0.95),size_perc=TRUE,position=TRUE,average=FALSE,cex=1.4,cex.axis=1.4,cex.lab=1.4,cex.main=2.5,
                       id_regions_subset=c('Control',idRegions(regionsFeatures_forward)[i+1]),col=c('red','blue'),ask=FALSE,xlab='Coordinates',
                       ylim=c(0,11))
}
dev.off()




### Tn datasets, grouped by feature lengths ###

# Create dataset
region='Tn'
lengths=length_100_range[[region]]
x=regionsFeatures_forward[c('Control',region),]
new_ids=paste0(region,lengths)
reg=x@regions[[region]]
feat=x@features[["IPD_Forward"]][[region]]
len=x@length_features[["IPD_Forward"]][[region]]
index=lapply(lengths,function(len) which(reg$length==len))
if(100 %in% lengths)
  index[[which(lengths==100)]]=which(reg$length>=100)
names(index)=new_ids
size_new=unlist(lapply(index,length))
feat_new=lapply(index,function(index) as.matrix(feat[,index]))
reg_new=lapply(index,function(index) reg[index])
len_new=lapply(index,function(index) len[index])
reg_datasets_new=as.data.frame(matrix(x@metadata$region_datasets[region,],nrow=length(lengths),ncol=3,byrow=TRUE),row.names=new_ids)
names(reg_datasets_new)=names(x@metadata$region_datasets)
reg_datasets_new$name=paste0(reg_datasets_new$name,' - Length ',lengths)
reg_datasets_new$size=size_new
x@metadata$region_datasets=rbind(x@metadata$region_datasets,reg_datasets_new)
x@regions=c(x@regions,GRangesList(reg_new))
mcols(x@regions[new_ids])=mcols(x@regions[region])
for(new_id in new_ids){
  NAbefore=pmax(unique(mcols(x@regions[[new_id]])$NAbefore),0)
  len=unique(pmin(unique(mcols(x@regions[[new_id]])$length),100))
  if(length(len)>1)
    stop('Something wrong with the length.')
  mcols(x@regions[new_id])[seq_len(min(NAbefore))]=NA
  mcols(x@regions[new_id])[100-seq_len(100-max(NAbefore)-len)+1]=NA
}
x@metadata$feature_datasets=x@metadata$feature_datasets[,c('name',paste0('file_',c(c('Control',region),rep(region,length(lengths)))),'resolution')]
names(x@metadata$feature_datasets)=c('name',paste0('file_',c('Control',region,new_ids)),'resolution')
x@features[["IPD_Forward"]]=c(x@features[["IPD_Forward"]],feat_new)
x@length_features[["IPD_Forward"]]=c(x@length_features[["IPD_Forward"]],len_new)
regionsFeatures_forward=x[setdiff(idRegions(x),region),]
save(regionsFeatures_forward,file='IPD_forward_Tn_different_lengths.RData')

# Subsample regions, maximum 10.000 for each region dataset
# 10 random subsamples
# to test for differences using IWT
maxsize=10000
regions_maxsize=rownames(regionsFeatures_forward@metadata$region_datasets)[regionsFeatures_forward@metadata$region_datasets$size>maxsize]
regionsFeatures_forward_complete=regionsFeatures_forward
# 10 times
for(i in 1:10){
  index=lapply(regions_maxsize,function(id) sort(sample(regionsFeatures_forward_complete@metadata$region_datasets[id,'size'],maxsize)))
  IPD_filtered=mapply(function(IPD,index) as.matrix(IPD[,index]),
                      regionsFeatures_forward_complete@features$IPD_Forward[regions_maxsize],index,SIMPLIFY=FALSE)
  regions_filtered=mapply(function(regions,index) regions[index,],
                          regionsFeatures_forward_complete@regions[regions_maxsize],index,SIMPLIFY=FALSE)
  length_features_filtered=mapply(function(length,index) length[index],
                                  regionsFeatures_forward_complete@length_features$IPD_Forward[regions_maxsize],index,SIMPLIFY=FALSE)
  regionsFeatures_forward=regionsFeatures_forward_complete
  regionsFeatures_forward@metadata$region_datasets[regions_maxsize,'size']=maxsize
  regionsFeatures_forward@regions[regions_maxsize]=GRangesList(regions_filtered)
  regionsFeatures_forward@features$IPD_Forward[regions_maxsize]=IPD_filtered
  regionsFeatures_forward@length_features$IPD_Forward[regions_maxsize]=length_features_filtered
  validObject(regionsFeatures_forward)
  
  # Add microsatellite sequences
  mcols(regionsFeatures_forward@regions) <- mcols(regionsFeatures_forward_complete@regions)
  
  save(regionsFeatures_forward,file=paste0('IPD_forward_Tn_different_lengths_max10000_',i,'.RData'))
}

# Test for difference using IWT
# 10 random subsample of maximum 10.000 regions for each region dataset
for(i in 1:10){
  load(paste0('IPD_forward_Tn_different_lengths_max10000_',i,'.RData'))
  
  # Test all features
  result_multi_quantile=IWTomicsTest(regionsFeatures_forward,
                                     setdiff(idRegions(regionsFeatures_forward),'Control'),rep('Control',nRegions(regionsFeatures_forward)-1),
                                     statistics='quantile',probs=c(0.05,0.25,0.5,0.75,0.95),B=10000)
  save(result_multi_quantile,file=paste0('IPD_forward_Tn_different_lengths_max10000_',i,'_results_quantiles_5_25_50_75_95.RData'))
  
  # Plot test results
  pdf(paste0('IPD_forward_Tn_different_lengths_max10000_',i,'_results_quantiles_5_25_50_75_95.pdf'),width=8.5,height=10)
  plotTest_same_sample_size(result_multi_quantile,ask=FALSE,col=c(rep('blue',nRegions(result_multi_quantile)-1),'red'),xlab='Coordinates',
                            average=FALSE,size=TRUE,size_perc=TRUE,position=TRUE,cex.axis=1.5,cex.lab=1.5,cex.main=2,probs=c(0.05,0.25,0.5,0.75,0.95))
  dev.off()
  
  # Test all features
  # number of control regions equal to feature regions
  result_multi_quantile=IWTomicsTest_same_sample_size(regionsFeatures_forward,
                                                      setdiff(idRegions(regionsFeatures_forward),'Control'),rep('Control',nRegions(regionsFeatures_forward)-1),
                                                      statistics='quantile',probs=c(0.05,0.25,0.5,0.75,0.95),B=10000)
  save(result_multi_quantile,file=paste0('IPD_forward_Tn_different_lengths_max10000_',i,'_results_same_sample_size_quantiles_5_25_50_75_95.RData'))
  
  # Plot test results
  pdf(paste0('IPD_forward_Tn_different_lengths_max10000_',i,'_results_same_sample_size_quantiles_5_25_50_75_95.pdf'),width=8.5,height=10)
  plotTest_same_sample_size(result_multi_quantile,ask=FALSE,col=c(rep('blue',nRegions(result_multi_quantile)-1),'red'),xlab='Coordinates',
                            average=FALSE,size=TRUE,size_perc=TRUE,position=TRUE,cex.axis=1.5,cex.lab=1.5,cex.main=2,probs=c(0.05,0.25,0.5,0.75,0.95))
  dev.off()
}

# Plot only reproducible results
alpha=0.05
min_significant=1
load('IPD_forward_Tn_different_lengths_max10000_1_results_quantiles_5_25_50_75_95.RData')
reproducible_pval_matrix=lapply(seq_len(nTests(result_multi_quantile)),function(name) matrix(TRUE,100,100))
for(i in 1:10){
  load(paste0('IPD_forward_Tn_different_lengths_max10000_',i,'_results_quantiles_5_25_50_75_95.RData'))
  reproducible_pval_matrix=mapply(function(reproducible_pval,result){
    reproducible_pval=(reproducible_pval & (result$IPD_Forward$adjusted_pval_matrix<=alpha))
  },reproducible_pval_matrix,result_multi_quantile@test$result,SIMPLIFY=FALSE)
  load(paste0('IPD_forward_Tn_different_lengths_max10000_',i,'_results_same_sample_size_quantiles_5_25_50_75_95.RData'))
  reproducible_pval_matrix=mapply(function(reproducible_pval,result){
    reproducible_pval=(reproducible_pval & (result$IPD_Forward$adjusted_pval_matrix<=alpha))
  },reproducible_pval_matrix,result_multi_quantile@test$result,SIMPLIFY=FALSE)
}
reproducible_pval_matrix=mapply(function(reproducible_pval){
  reproducible_pval[rowSums(reproducible_pval)<min_significant,]=matrix(FALSE,nrow=sum(rowSums(reproducible_pval)<min_significant),ncol=100)
  return(reproducible_pval)
},reproducible_pval_matrix,SIMPLIFY=FALSE)
save(reproducible_pval_matrix,file='IPD_forward_Tn_different_lengths_max10000_reproducible_results_quantiles_5_25_50_75_95.RData')

load('IPD_forward_Tn_different_lengths.RData')
pdf('boxplot_vs_control_forward_Tn_different_lengths_reproducible_results_quantiles_5_25_50_75_95.pdf',width=9,height=6)
for(i in seq_len(nrow(regionsFeatures_forward)-1)){
  plot_no_size_control(regionsFeatures_forward,reproducible_pval=reproducible_pval_matrix[[i]],
                       type='boxplot',probs=c(0.05,0.25,0.5,0.75,0.95),size_perc=TRUE,position=TRUE,average=FALSE,cex=1.4,cex.axis=1.4,cex.lab=1.4,cex.main=2.5,
                       id_regions_subset=c('Control',idRegions(regionsFeatures_forward)[i+1]),col=c('red','blue'),ask=FALSE,xlab='Coordinates',
                       ylim=c(0,11))
}
dev.off()















#######################
##### MINUS STRAND ####
#######################

load('IPD_reverse.RData')

# Feature length histograms
length=lapply(regionsFeatures_reverse@regions[-84],function(region) region$length)
length100=lapply(length,pmin,100) # Set maximum feature length to 100 (window length)
index=lapply(length100,order)
length_100_range=length100
pdf("length_distribution_reverse.pdf",width=7,height=5)
for(i in seq_len(nrow(regionsFeatures_reverse)-1)){
  counts=hist(length100[[i]],breaks=0:100+0.5,xlim=c(0.5,100.5),col="blue",xlab="Length",main=nameRegions(regionsFeatures_reverse)[i])$counts
  length_100_range[[i]]=which(counts>=10)
}
dev.off()

# Boxplot of log mean IPD in windows, grouped by feature lengths
meanIPD=lapply(regionsFeatures_reverse@features$IPD_Reverse[-84],function(feature) colMeans(feature,na.rm=TRUE))
pdf("length_vs_meanIPD_log_reverse_boxplot.pdf",width=7,height=7)
for(i in seq_len(nrow(regionsFeatures_reverse)-1)){
  box=boxplot(log(meanIPD[[i]][index[[i]]]+0.01)~factor(length100[[i]][index[[i]]]),plot=FALSE)
  box$stats=as.matrix(Reduce(cbind,tapply(log(meanIPD[[i]][index[[i]]]+0.01),factor(length100[[i]][index[[i]]]),quantile,probs=c(0.05,0.25,0.5,0.75,0.95),na.rm=TRUE,SIMPLIFY=FALSE)))
  perc=round(box$n/max(box$n)*100000,0)
  col_perc=rev(heat.colors(100000))[perc]
  bxp(box,at=sort(unique(length100[[i]])),xlim=c(0,100),outline=FALSE,pars=list(whisklty=3,staplewex=0.8,boxfill=col_perc),
      xlab="Length",ylab="Log Mean IPD",main=nameRegions(regionsFeatures_reverse)[i])
}
dev.off()

# Boxplot of log mean IPD in features, grouped by feature lengths
meanIPD_feature=mapply(function(feature,region){
  pos=apply(as.matrix(mcols(region)),1,
            function(info){
              info=pmax(pmin(info[1:2],100),0)
              info[2]=min((100-info[1]),info[2])
              return(c(rep(NA,info[2]),rep(1,info[1]),rep(NA,100-sum(info))))})
  feature=feature*pos
  meanIPD=colMeans(feature,na.rm=TRUE)
  meanIPD[is.nan(meanIPD)]=NA
  return(meanIPD)
},
regionsFeatures_reverse@features$IPD_Reverse[-84],regionsFeatures_reverse@regions[-84],SIMPLIFY=FALSE)
pdf("length_vs_meanIPD_log_feature_reverse_boxplot.pdf",width=7,height=7)
for(i in seq_len(nrow(regionsFeatures_reverse)-1)){
  box=boxplot(log(meanIPD_feature[[i]][index[[i]]]+0.01)~factor(length100[[i]][index[[i]]]),plot=FALSE)
  box$stats=as.matrix(Reduce(cbind,tapply(log(meanIPD_feature[[i]][index[[i]]]+0.01),factor(length100[[i]][index[[i]]]),quantile,probs=c(0.05,0.25,0.5,0.75,0.95),na.rm=TRUE,SIMPLIFY=FALSE)))
  perc=round(box$n/max(box$n)*100000,0)
  col_perc=rev(heat.colors(100000))[perc]
  bxp(box,at=sort(unique(length100[[i]])),xlim=c(0,100),outline=FALSE,pars=list(whisklty=3,staplewex=0.8,boxfill=col_perc),
      xlab="Length",ylab="Log Mean IPD",main=nameRegions(regionsFeatures_reverse)[i])
}
dev.off()

# Pointwise boxplots, grouped by feature lengths
for(region in idRegions(regionsFeatures_reverse)){
  if(length(length_100_range[[region]])>0){
    pdf(paste0('boxplot_vs_control_reverse_',region,'.pdf'),width=9,height=6)
    plot_no_size_control(regionsFeatures_reverse,lengths=length_100_range[[region]],
                         type='boxplot',probs=c(0.05,0.25,0.5,0.75,0.95),size_perc=TRUE,position=TRUE,average=FALSE,cex=1.4,cex.axis=1.4,cex.lab=1.4,cex.main=2.5,
                         id_regions_subset=c('Control',region),col=c('red','blue'),ask=FALSE,xlab='Coordinates')
    dev.off()
  }
}




### An datasets, grouped by feature lengths ###

# Create dataset
region='An'
lengths=length_100_range[[region]]
x=regionsFeatures_reverse[c('Control',region),]
new_ids=paste0(region,lengths)
reg=x@regions[[region]]
feat=x@features[["IPD_Reverse"]][[region]]
len=x@length_features[["IPD_Reverse"]][[region]]
index=lapply(lengths,function(len) which(reg$length==len))
if(100 %in% lengths)
  index[[which(lengths==100)]]=which(reg$length>=100)
names(index)=new_ids
size_new=unlist(lapply(index,length))
feat_new=lapply(index,function(index) as.matrix(feat[,index]))
reg_new=lapply(index,function(index) reg[index])
len_new=lapply(index,function(index) len[index])
reg_datasets_new=as.data.frame(matrix(x@metadata$region_datasets[region,],nrow=length(lengths),ncol=3,byrow=TRUE),row.names=new_ids)
names(reg_datasets_new)=names(x@metadata$region_datasets)
reg_datasets_new$name=paste0(reg_datasets_new$name,' - Length ',lengths)
reg_datasets_new$size=size_new
x@metadata$region_datasets=rbind(x@metadata$region_datasets,reg_datasets_new)
x@regions=c(x@regions,GRangesList(reg_new))
mcols(x@regions[new_ids])=mcols(x@regions[region])
for(new_id in new_ids){
  NAbefore=pmax(unique(mcols(x@regions[[new_id]])$NAbefore),0)
  len=unique(pmin(unique(mcols(x@regions[[new_id]])$length),100))
  if(length(len)>1)
    stop('Something wrong with the length.')
  mcols(x@regions[new_id])[seq_len(min(NAbefore))]=NA
  mcols(x@regions[new_id])[100-seq_len(100-max(NAbefore)-len)+1]=NA
}
x@metadata$feature_datasets=x@metadata$feature_datasets[,c('name',paste0('file_',c(c('Control',region),rep(region,length(lengths)))),'resolution')]
names(x@metadata$feature_datasets)=c('name',paste0('file_',c('Control',region,new_ids)),'resolution')
x@features[["IPD_Reverse"]]=c(x@features[["IPD_Reverse"]],feat_new)
x@length_features[["IPD_Reverse"]]=c(x@length_features[["IPD_Reverse"]],len_new)
regionsFeatures_reverse=x[setdiff(idRegions(x),region),]
save(regionsFeatures_reverse,file='IPD_reverse_An_different_lengths.RData')

# Subsample regions, maximum 10.000 for each region dataset
# 10 random subsamples
# to test for differences using IWT
maxsize=10000
regions_maxsize=rownames(regionsFeatures_reverse@metadata$region_datasets)[regionsFeatures_reverse@metadata$region_datasets$size>maxsize]
regionsFeatures_reverse_complete=regionsFeatures_reverse
# 10 times
for(i in 1:10){
  index=lapply(regions_maxsize,function(id) sort(sample(regionsFeatures_reverse_complete@metadata$region_datasets[id,'size'],maxsize)))
  IPD_filtered=mapply(function(IPD,index) as.matrix(IPD[,index]),
                      regionsFeatures_reverse_complete@features$IPD_Reverse[regions_maxsize],index,SIMPLIFY=FALSE)
  regions_filtered=mapply(function(regions,index) regions[index,],
                          regionsFeatures_reverse_complete@regions[regions_maxsize],index,SIMPLIFY=FALSE)
  length_features_filtered=mapply(function(length,index) length[index],
                                  regionsFeatures_reverse_complete@length_features$IPD_Reverse[regions_maxsize],index,SIMPLIFY=FALSE)
  regionsFeatures_reverse=regionsFeatures_reverse_complete
  regionsFeatures_reverse@metadata$region_datasets[regions_maxsize,'size']=maxsize
  regionsFeatures_reverse@regions[regions_maxsize]=GRangesList(regions_filtered)
  regionsFeatures_reverse@features$IPD_Reverse[regions_maxsize]=IPD_filtered
  regionsFeatures_reverse@length_features$IPD_Reverse[regions_maxsize]=length_features_filtered
  validObject(regionsFeatures_reverse)
  
  # Add microsatellite sequences
  mcols(regionsFeatures_reverse@regions) <- mcols(regionsFeatures_reverse_complete@regions)
  
  save(regionsFeatures_reverse,file=paste0('IPD_reverse_An_different_lengths_max10000_',i,'.RData'))
}

# Test for difference using IWT
# 10 random subsample of maximum 10.000 regions for each region dataset
for(i in 1:10){
  load(paste0('IPD_reverse_An_different_lengths_max10000_',i,'.RData'))
  
  # Test all features
  result_multi_quantile=IWTomicsTest(regionsFeatures_reverse,
                                     setdiff(idRegions(regionsFeatures_reverse),'Control'),rep('Control',nRegions(regionsFeatures_reverse)-1),
                                     statistics='quantile',probs=c(0.05,0.25,0.5,0.75,0.95),B=10000)
  save(result_multi_quantile,file=paste0('IPD_reverse_An_different_lengths_max10000_',i,'_results_quantiles_5_25_50_75_95.RData'))
  
  # Plot test results
  pdf(paste0('IPD_reverse_An_different_lengths_max10000_',i,'_results_quantiles_5_25_50_75_95.pdf'),width=8.5,height=10)
  plotTest_same_sample_size(result_multi_quantile,ask=FALSE,col=c(rep('blue',nRegions(result_multi_quantile)-1),'red'),xlab='Coordinates',
                            average=FALSE,size=TRUE,size_perc=TRUE,position=TRUE,cex.axis=1.5,cex.lab=1.5,cex.main=2,probs=c(0.05,0.25,0.5,0.75,0.95))
  dev.off()
  
  # Test all features
  # number of control regions equal to feature regions
  result_multi_quantile=IWTomicsTest_same_sample_size(regionsFeatures_reverse,
                                                      setdiff(idRegions(regionsFeatures_reverse),'Control'),rep('Control',nRegions(regionsFeatures_reverse)-1),
                                                      statistics='quantile',probs=c(0.05,0.25,0.5,0.75,0.95),B=10000)
  save(result_multi_quantile,file=paste0('IPD_reverse_An_different_lengths_max10000_',i,'_results_same_sample_size_quantiles_5_25_50_75_95.RData'))
  
  # Plot test results
  pdf(paste0('IPD_reverse_An_different_lengths_max10000_',i,'_results_same_sample_size_quantiles_5_25_50_75_95.pdf'),width=8.5,height=10)
  plotTest_same_sample_size(result_multi_quantile,ask=FALSE,col=c(rep('blue',nRegions(result_multi_quantile)-1),'red'),xlab='Coordinates',
                            average=FALSE,size=TRUE,size_perc=TRUE,position=TRUE,cex.axis=1.5,cex.lab=1.5,cex.main=2,probs=c(0.05,0.25,0.5,0.75,0.95))
  dev.off()
}

# Plot only reproducible results
alpha=0.05
min_significant=1
load('IPD_reverse_An_different_lengths_max10000_1_results_quantiles_5_25_50_75_95.RData')
reproducible_pval_matrix=lapply(seq_len(nTests(result_multi_quantile)),function(name) matrix(TRUE,100,100))
for(i in 1:10){
  load(paste0('IPD_reverse_An_different_lengths_max10000_',i,'_results_quantiles_5_25_50_75_95.RData'))
  reproducible_pval_matrix=mapply(function(reproducible_pval,result){
    reproducible_pval=(reproducible_pval & (result$IPD_Reverse$adjusted_pval_matrix<=alpha))
  },reproducible_pval_matrix,result_multi_quantile@test$result,SIMPLIFY=FALSE)
  load(paste0('IPD_reverse_An_different_lengths_max10000_',i,'_results_same_sample_size_quantiles_5_25_50_75_95.RData'))
  reproducible_pval_matrix=mapply(function(reproducible_pval,result){
    reproducible_pval=(reproducible_pval & (result$IPD_Reverse$adjusted_pval_matrix<=alpha))
  },reproducible_pval_matrix,result_multi_quantile@test$result,SIMPLIFY=FALSE)
}
reproducible_pval_matrix=mapply(function(reproducible_pval){
  reproducible_pval[rowSums(reproducible_pval)<min_significant,]=matrix(FALSE,nrow=sum(rowSums(reproducible_pval)<min_significant),ncol=100)
  return(reproducible_pval)
},reproducible_pval_matrix,SIMPLIFY=FALSE)
save(reproducible_pval_matrix,file='IPD_reverse_An_different_lengths_max10000_reproducible_results_quantiles_5_25_50_75_95.RData')

load('IPD_reverse_An_different_lengths.RData')
pdf('boxplot_vs_control_reverse_An_different_lengths_reproducible_results_quantiles_5_25_50_75_95.pdf',width=9,height=6)
for(i in seq_len(nrow(regionsFeatures_reverse)-1)){
  plot_no_size_control(regionsFeatures_reverse,reproducible_pval=reproducible_pval_matrix[[i]],
                       type='boxplot',probs=c(0.05,0.25,0.5,0.75,0.95),size_perc=TRUE,position=TRUE,average=FALSE,cex=1.4,cex.axis=1.4,cex.lab=1.4,cex.main=2.5,
                       id_regions_subset=c('Control',idRegions(regionsFeatures_reverse)[i+1]),col=c('red','blue'),ask=FALSE,xlab='Coordinates',
                       ylim=c(0,11))
}
dev.off()




### Tn datasets, grouped by feature lengths ###

# Create dataset
region='Tn'
lengths=length_100_range[[region]]
x=regionsFeatures_reverse[c('Control',region),]
new_ids=paste0(region,lengths)
reg=x@regions[[region]]
feat=x@features[["IPD_Reverse"]][[region]]
len=x@length_features[["IPD_Reverse"]][[region]]
index=lapply(lengths,function(len) which(reg$length==len))
if(100 %in% lengths)
  index[[which(lengths==100)]]=which(reg$length>=100)
names(index)=new_ids
size_new=unlist(lapply(index,length))
feat_new=lapply(index,function(index) as.matrix(feat[,index]))
reg_new=lapply(index,function(index) reg[index])
len_new=lapply(index,function(index) len[index])
reg_datasets_new=as.data.frame(matrix(x@metadata$region_datasets[region,],nrow=length(lengths),ncol=3,byrow=TRUE),row.names=new_ids)
names(reg_datasets_new)=names(x@metadata$region_datasets)
reg_datasets_new$name=paste0(reg_datasets_new$name,' - Length ',lengths)
reg_datasets_new$size=size_new
x@metadata$region_datasets=rbind(x@metadata$region_datasets,reg_datasets_new)
x@regions=c(x@regions,GRangesList(reg_new))
mcols(x@regions[new_ids])=mcols(x@regions[region])
for(new_id in new_ids){
  NAbefore=pmax(unique(mcols(x@regions[[new_id]])$NAbefore),0)
  len=unique(pmin(unique(mcols(x@regions[[new_id]])$length),100))
  if(length(len)>1)
    stop('Something wrong with the length.')
  mcols(x@regions[new_id])[seq_len(min(NAbefore))]=NA
  mcols(x@regions[new_id])[100-seq_len(100-max(NAbefore)-len)+1]=NA
}
x@metadata$feature_datasets=x@metadata$feature_datasets[,c('name',paste0('file_',c(c('Control',region),rep(region,length(lengths)))),'resolution')]
names(x@metadata$feature_datasets)=c('name',paste0('file_',c('Control',region,new_ids)),'resolution')
x@features[["IPD_Reverse"]]=c(x@features[["IPD_Reverse"]],feat_new)
x@length_features[["IPD_Reverse"]]=c(x@length_features[["IPD_Reverse"]],len_new)
regionsFeatures_reverse=x[setdiff(idRegions(x),region),]
save(regionsFeatures_reverse,file='IPD_reverse_Tn_different_lengths.RData')

# Subsample regions, maximum 10.000 for each region dataset
# 10 random subsamples
# to test for differences using IWT
maxsize=10000
regions_maxsize=rownames(regionsFeatures_reverse@metadata$region_datasets)[regionsFeatures_reverse@metadata$region_datasets$size>maxsize]
regionsFeatures_reverse_complete=regionsFeatures_reverse
# 10 times
for(i in 1:10){
  index=lapply(regions_maxsize,function(id) sort(sample(regionsFeatures_reverse_complete@metadata$region_datasets[id,'size'],maxsize)))
  IPD_filtered=mapply(function(IPD,index) as.matrix(IPD[,index]),
                      regionsFeatures_reverse_complete@features$IPD_Reverse[regions_maxsize],index,SIMPLIFY=FALSE)
  regions_filtered=mapply(function(regions,index) regions[index,],
                          regionsFeatures_reverse_complete@regions[regions_maxsize],index,SIMPLIFY=FALSE)
  length_features_filtered=mapply(function(length,index) length[index],
                                  regionsFeatures_reverse_complete@length_features$IPD_Reverse[regions_maxsize],index,SIMPLIFY=FALSE)
  regionsFeatures_reverse=regionsFeatures_reverse_complete
  regionsFeatures_reverse@metadata$region_datasets[regions_maxsize,'size']=maxsize
  regionsFeatures_reverse@regions[regions_maxsize]=GRangesList(regions_filtered)
  regionsFeatures_reverse@features$IPD_Reverse[regions_maxsize]=IPD_filtered
  regionsFeatures_reverse@length_features$IPD_Reverse[regions_maxsize]=length_features_filtered
  validObject(regionsFeatures_reverse)
  
  # Add microsatellite sequences
  mcols(regionsFeatures_reverse@regions) <- mcols(regionsFeatures_reverse_complete@regions)
  
  save(regionsFeatures_reverse,file=paste0('IPD_reverse_Tn_different_lengths_max10000_',i,'.RData'))
}

# Test for difference using IWT
# 10 random subsample of maximum 10.000 regions for each region dataset
for(i in 1:10){
  load(paste0('IPD_reverse_Tn_different_lengths_max10000_',i,'.RData'))
  
  # Test all features
  result_multi_quantile=IWTomicsTest(regionsFeatures_reverse,
                                     setdiff(idRegions(regionsFeatures_reverse),'Control'),rep('Control',nRegions(regionsFeatures_reverse)-1),
                                     statistics='quantile',probs=c(0.05,0.25,0.5,0.75,0.95),B=10000)
  save(result_multi_quantile,file=paste0('IPD_reverse_Tn_different_lengths_max10000_',i,'_results_quantiles_5_25_50_75_95.RData'))
  
  # Plot test results
  pdf(paste0('IPD_reverse_Tn_different_lengths_max10000_',i,'_results_quantiles_5_25_50_75_95.pdf'),width=8.5,height=10)
  plotTest_same_sample_size(result_multi_quantile,ask=FALSE,col=c(rep('blue',nRegions(result_multi_quantile)-1),'red'),xlab='Coordinates',
                            average=FALSE,size=TRUE,size_perc=TRUE,position=TRUE,cex.axis=1.5,cex.lab=1.5,cex.main=2,probs=c(0.05,0.25,0.5,0.75,0.95))
  dev.off()
  
  # Test all features
  # number of control regions equal to feature regions
  result_multi_quantile=IWTomicsTest_same_sample_size(regionsFeatures_reverse,
                                                      setdiff(idRegions(regionsFeatures_reverse),'Control'),rep('Control',nRegions(regionsFeatures_reverse)-1),
                                                      statistics='quantile',probs=c(0.05,0.25,0.5,0.75,0.95),B=10000)
  save(result_multi_quantile,file=paste0('IPD_reverse_Tn_different_lengths_max10000_',i,'_results_same_sample_size_quantiles_5_25_50_75_95.RData'))
  
  # Plot test results
  pdf(paste0('IPD_reverse_Tn_different_lengths_max10000_',i,'_results_same_sample_size_quantiles_5_25_50_75_95.pdf'),width=8.5,height=10)
  plotTest_same_sample_size(result_multi_quantile,ask=FALSE,col=c(rep('blue',nRegions(result_multi_quantile)-1),'red'),xlab='Coordinates',
                            average=FALSE,size=TRUE,size_perc=TRUE,position=TRUE,cex.axis=1.5,cex.lab=1.5,cex.main=2,probs=c(0.05,0.25,0.5,0.75,0.95))
  dev.off()
}

# Plot only reproducible results
alpha=0.05
min_significant=1
load('IPD_reverse_Tn_different_lengths_max10000_1_results_quantiles_5_25_50_75_95.RData')
reproducible_pval_matrix=lapply(seq_len(nTests(result_multi_quantile)),function(name) matrix(TRUE,100,100))
for(i in 1:10){
  load(paste0('IPD_reverse_Tn_different_lengths_max10000_',i,'_results_quantiles_5_25_50_75_95.RData'))
  reproducible_pval_matrix=mapply(function(reproducible_pval,result){
    reproducible_pval=(reproducible_pval & (result$IPD_Reverse$adjusted_pval_matrix<=alpha))
  },reproducible_pval_matrix,result_multi_quantile@test$result,SIMPLIFY=FALSE)
  load(paste0('IPD_reverse_Tn_different_lengths_max10000_',i,'_results_same_sample_size_quantiles_5_25_50_75_95.RData'))
  reproducible_pval_matrix=mapply(function(reproducible_pval,result){
    reproducible_pval=(reproducible_pval & (result$IPD_Reverse$adjusted_pval_matrix<=alpha))
  },reproducible_pval_matrix,result_multi_quantile@test$result,SIMPLIFY=FALSE)
}
reproducible_pval_matrix=mapply(function(reproducible_pval){
  reproducible_pval[rowSums(reproducible_pval)<min_significant,]=matrix(FALSE,nrow=sum(rowSums(reproducible_pval)<min_significant),ncol=100)
  return(reproducible_pval)
},reproducible_pval_matrix,SIMPLIFY=FALSE)
save(reproducible_pval_matrix,file='IPD_reverse_Tn_different_lengths_max10000_reproducible_results_quantiles_5_25_50_75_95.RData')

load('IPD_reverse_Tn_different_lengths.RData')
pdf('boxplot_vs_control_reverse_Tn_different_lengths_reproducible_results_quantiles_5_25_50_75_95.pdf',width=9,height=6)
for(i in seq_len(nrow(regionsFeatures_reverse)-1)){
  plot_no_size_control(regionsFeatures_reverse,reproducible_pval=reproducible_pval_matrix[[i]],
                       type='boxplot',probs=c(0.05,0.25,0.5,0.75,0.95),size_perc=TRUE,position=TRUE,average=FALSE,cex=1.4,cex.axis=1.4,cex.lab=1.4,cex.main=2.5,
                       id_regions_subset=c('Control',idRegions(regionsFeatures_reverse)[i+1]),col=c('red','blue'),ask=FALSE,xlab='Coordinates',
                       ylim=c(0,11))
}
dev.off()



