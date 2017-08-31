require(IWTomics)
source('IWTomics_modified_functions.r')



# INPUT FROM FILE load_IPD_data.r


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
              only_significant=FALSE,gaps_tests=c(7,11,17,37),align_lab='Motif center',ask=FALSE,xlab='Coordinates',cellwidth=10,cellheight=10)
  
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
              only_significant=FALSE,gaps_tests=c(7,11,17,37),align_lab='Motif center',ask=FALSE,xlab='Coordinates',cellwidth=10,cellheight=10)
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

# Load mean depth in each position
depth=as.matrix(read.delim('forward/depth_mean.txt',
                           sep='\t',header=FALSE,row.names=1,col.names=0:100))
# depth_mean.txt is a text file (tab separated) in the following format: id   depth_1   ...   depth_100
# with the average number of reads used to compute the IPD in wach of the 100 nucleotides
# Direct_Repeats	24.07356475	24.05188188	24.0728078	24.04015518	24.05769788	24.01441598	24.01398419	24.06720134	24.04421089	24.04415299	24.02012393	24.03552197	24.02312858	24.03398564	24.01835394	24.01609448	24.01464841	23.98292033	24.0125351	24.00602131	23.96355689	23.99108254	23.97287517	24.0145603	23.97828667	24.01836085	23.9884788	23.94800081	23.98552319	23.97848809	23.98905522	24.00344568	23.97955816	24.01818708	23.97524251	24.01016478	24.02145339	23.99088146	24.02183419	24.05116845	24.03075052	24.02232388	24.02232776	23.97417188	23.95456651	23.91122118	23.94982484	23.90875923	23.87848248	23.83410606	23.8482355	23.87390599	23.88136232	23.86868804	23.83164759	23.78518025	23.77242178	23.70282335	23.68041147	23.64443221	23.57346974	23.58616093	23.63882128	23.6174238	23.67478921	23.67214255	23.69950468	23.73646241	23.69235226	23.76979022	23.74214019	23.76077761	23.81020249	23.80385852	23.81526337	23.83096266	23.85785449	23.8654654	23.83893492	23.88461093	23.87777681	23.88404832	23.89246813	23.90569153	23.89341793	23.90741867	23.89201089	23.90731354	23.93448656	23.90126854	23.93518116	23.8789213	23.93776066	23.9201228	23.90386955	23.93003853	23.88798147	23.91288734	23.94520389	23.94453455
# Inverted_Repeats	24.11423662	24.1033751	24.08757856	24.09426607	24.09771062	24.09319811	24.10347168	24.10347389	24.11615649	24.09304843	24.11192682	24.09687067	24.09054268	24.09301938	24.10900626	24.09893244	24.09243581	24.08326324	24.10912867	24.12452023	24.07547323	24.11131391	24.0735086	24.13078516	24.12736048	24.09181887	24.12241918	24.11018559	24.0844525	24.13411778	24.11763425	24.08641551	24.11596229	24.08126051	24.09855754	24.13140553	24.07512056	24.10440207	24.11115831	24.11827996	24.14429383	24.09330777	24.11619318	24.1538606	24.2019341	24.04360124	24.11433768	24.07683489	24.02232531	23.99456056	24.14336431	24.19017915	24.19660236	24.15415159	24.14748226	24.07154248	23.92504392	24.05486734	24.12291169	24.050141	24.11593006	24.08429841	24.13672847	24.10229929	24.05566087	24.0587949	24.12579465	24.10219199	24.08706446	24.06503942	24.08796977	24.10474149	24.04729945	24.10933214	24.10988267	24.08807657	24.10424312	24.08209484	24.10918581	24.12383843	24.06126052	24.0807112	24.12684415	24.06240567	24.12555316	24.06995575	24.10269208	24.10519912	24.08216795	24.09287417	24.10870651	24.09801921	24.08031782	24.10199229	24.06836331	24.10121591	24.10292174	24.09275483	24.10635477	24.10779188
# Mirror_Repeats	23.95236663	23.94548006	23.95271625	23.94765985	23.94673281	23.95740523	23.96028564	23.92284948	23.93351877	23.91640983	23.92974368	23.9270896	23.94450152	23.95539853	23.96068457	23.93622992	23.93275037	23.90510675	23.9342036	23.93517753	23.95451926	23.96276529	23.96103343	23.9530467	23.94775615	23.95725489	23.96685041	23.9500877	23.94703082	23.96158954	23.92687633	23.98621796	23.97481834	23.96993235	23.96141511	23.99338512	23.97647206	23.98398938	23.98128382	23.99566536	24.03605703	24.05291907	24.04733268	24.01971937	23.95727781	23.77094244	23.83559687	23.93340349	23.98877474	24.02680965	24.025057	23.95126901	23.84942114	23.93048664	24.01696104	24.00684108	23.94643081	23.89115783	23.85646269	23.8215646	23.8681718	23.89428715	23.88696502	23.9012454	23.91001804	23.90135304	23.9119998	23.89944372	23.88489966	23.90956927	23.90378592	23.91679194	23.91613841	23.93618461	23.92640573	23.9387121	23.94025811	23.93648049	23.93264677	23.94264308	23.95171515	23.93760649	23.94184706	23.95486781	23.96130811	23.92064049	23.907232	23.90789045	23.94259872	23.91460787	23.92439354	23.95176406	23.93818592	23.94777466	23.87908096	23.92071566	23.93635521	23.93565522	23.96393936	23.94001052
# ...
depth=t((depth/matrix(depth[nrow(depth),],nrow=nrow(depth),ncol=100,byrow=TRUE))[-nrow(depth),])

load('IPD_forward_max10000_1_results_mean.RData')
plotSummary_reproducible(result_mean,groupby='feature',filenames=paste0('IPD_forward_max10000_reproducible_results_mean_summary.pdf'),
                         reproducible_pval_matrix=reproducible_pval_matrix,only_significant=TRUE,
                         gaps_tests=c(7,11,17,37),align_lab='Motif center',xlab='',ylab='',ask=FALSE,cellwidth=5,cellheight=12,fontsize=12)
load('IPD_forward.RData')
pdf('boxplot_vs_control_forward_reproducible_results_mean.pdf',width=9,height=6)
for(i in seq_len(nrow(regionsFeatures_forward)-1)){
  plot_no_size_control(regionsFeatures_forward,reproducible_pval=reproducible_pval_matrix[[i]],depth=depth,
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
              only_significant=FALSE,gaps_tests=c(7,11,17,37),align_lab='Motif center',ask=FALSE,xlab='Coordinates',cellwidth=10,cellheight=10)
  
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
              only_significant=FALSE,gaps_tests=c(7,11,17,37),align_lab='Motif center',ask=FALSE,xlab='Coordinates',cellwidth=10,cellheight=10)
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

# Load mean depth in each position
depth=as.matrix(read.delim('forward/depth_mean.txt',
                           sep='\t',header=FALSE,row.names=1,col.names=0:100))
# depth_mean.txt is a text file (tab separated) in the following format: id   depth_1   ...   depth_100
# with the average number of reads used to compute the IPD in wach of the 100 nucleotides
# Direct_Repeats	24.07356475	24.05188188	24.0728078	24.04015518	24.05769788	24.01441598	24.01398419	24.06720134	24.04421089	24.04415299	24.02012393	24.03552197	24.02312858	24.03398564	24.01835394	24.01609448	24.01464841	23.98292033	24.0125351	24.00602131	23.96355689	23.99108254	23.97287517	24.0145603	23.97828667	24.01836085	23.9884788	23.94800081	23.98552319	23.97848809	23.98905522	24.00344568	23.97955816	24.01818708	23.97524251	24.01016478	24.02145339	23.99088146	24.02183419	24.05116845	24.03075052	24.02232388	24.02232776	23.97417188	23.95456651	23.91122118	23.94982484	23.90875923	23.87848248	23.83410606	23.8482355	23.87390599	23.88136232	23.86868804	23.83164759	23.78518025	23.77242178	23.70282335	23.68041147	23.64443221	23.57346974	23.58616093	23.63882128	23.6174238	23.67478921	23.67214255	23.69950468	23.73646241	23.69235226	23.76979022	23.74214019	23.76077761	23.81020249	23.80385852	23.81526337	23.83096266	23.85785449	23.8654654	23.83893492	23.88461093	23.87777681	23.88404832	23.89246813	23.90569153	23.89341793	23.90741867	23.89201089	23.90731354	23.93448656	23.90126854	23.93518116	23.8789213	23.93776066	23.9201228	23.90386955	23.93003853	23.88798147	23.91288734	23.94520389	23.94453455
# Inverted_Repeats	24.11423662	24.1033751	24.08757856	24.09426607	24.09771062	24.09319811	24.10347168	24.10347389	24.11615649	24.09304843	24.11192682	24.09687067	24.09054268	24.09301938	24.10900626	24.09893244	24.09243581	24.08326324	24.10912867	24.12452023	24.07547323	24.11131391	24.0735086	24.13078516	24.12736048	24.09181887	24.12241918	24.11018559	24.0844525	24.13411778	24.11763425	24.08641551	24.11596229	24.08126051	24.09855754	24.13140553	24.07512056	24.10440207	24.11115831	24.11827996	24.14429383	24.09330777	24.11619318	24.1538606	24.2019341	24.04360124	24.11433768	24.07683489	24.02232531	23.99456056	24.14336431	24.19017915	24.19660236	24.15415159	24.14748226	24.07154248	23.92504392	24.05486734	24.12291169	24.050141	24.11593006	24.08429841	24.13672847	24.10229929	24.05566087	24.0587949	24.12579465	24.10219199	24.08706446	24.06503942	24.08796977	24.10474149	24.04729945	24.10933214	24.10988267	24.08807657	24.10424312	24.08209484	24.10918581	24.12383843	24.06126052	24.0807112	24.12684415	24.06240567	24.12555316	24.06995575	24.10269208	24.10519912	24.08216795	24.09287417	24.10870651	24.09801921	24.08031782	24.10199229	24.06836331	24.10121591	24.10292174	24.09275483	24.10635477	24.10779188
# Mirror_Repeats	23.95236663	23.94548006	23.95271625	23.94765985	23.94673281	23.95740523	23.96028564	23.92284948	23.93351877	23.91640983	23.92974368	23.9270896	23.94450152	23.95539853	23.96068457	23.93622992	23.93275037	23.90510675	23.9342036	23.93517753	23.95451926	23.96276529	23.96103343	23.9530467	23.94775615	23.95725489	23.96685041	23.9500877	23.94703082	23.96158954	23.92687633	23.98621796	23.97481834	23.96993235	23.96141511	23.99338512	23.97647206	23.98398938	23.98128382	23.99566536	24.03605703	24.05291907	24.04733268	24.01971937	23.95727781	23.77094244	23.83559687	23.93340349	23.98877474	24.02680965	24.025057	23.95126901	23.84942114	23.93048664	24.01696104	24.00684108	23.94643081	23.89115783	23.85646269	23.8215646	23.8681718	23.89428715	23.88696502	23.9012454	23.91001804	23.90135304	23.9119998	23.89944372	23.88489966	23.90956927	23.90378592	23.91679194	23.91613841	23.93618461	23.92640573	23.9387121	23.94025811	23.93648049	23.93264677	23.94264308	23.95171515	23.93760649	23.94184706	23.95486781	23.96130811	23.92064049	23.907232	23.90789045	23.94259872	23.91460787	23.92439354	23.95176406	23.93818592	23.94777466	23.87908096	23.92071566	23.93635521	23.93565522	23.96393936	23.94001052
# ...
depth=t((depth/matrix(depth[nrow(depth),],nrow=nrow(depth),ncol=100,byrow=TRUE))[-nrow(depth),])

load('IPD_forward_max10000_1_results_median.RData')
plotSummary_reproducible(result_median,groupby='feature',filenames=paste0('IPD_forward_max10000_reproducible_results_median_summary.pdf'),
                         reproducible_pval_matrix=reproducible_pval_matrix,only_significant=TRUE,
                         gaps_tests=c(7,11,17,37),align_lab='Motif center',xlab='',ylab='',ask=FALSE,cellwidth=5,cellheight=12,fontsize=12)
load('IPD_forward.RData')
pdf('boxplot_vs_control_forward_reproducible_results_median.pdf',width=9,height=6)
for(i in seq_len(nrow(regionsFeatures_forward)-1)){
  plot_no_size_control(regionsFeatures_forward,reproducible_pval=reproducible_pval_matrix[[i]],depth=depth,
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
              only_significant=FALSE,gaps_tests=c(7,11,17,37),align_lab='Motif center',ask=FALSE,xlab='Coordinates',cellwidth=10,cellheight=10)
  
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
              only_significant=FALSE,gaps_tests=c(7,11,17,37),align_lab='Motif center',ask=FALSE,xlab='Coordinates',cellwidth=10,cellheight=10)
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

# Load mean depth in each position
depth=as.matrix(read.delim('forward/depth_mean.txt',
                           sep='\t',header=FALSE,row.names=1,col.names=0:100))
# depth_mean.txt is a text file (tab separated) in the following format: id   depth_1   ...   depth_100
# with the average number of reads used to compute the IPD in wach of the 100 nucleotides
# Direct_Repeats	24.07356475	24.05188188	24.0728078	24.04015518	24.05769788	24.01441598	24.01398419	24.06720134	24.04421089	24.04415299	24.02012393	24.03552197	24.02312858	24.03398564	24.01835394	24.01609448	24.01464841	23.98292033	24.0125351	24.00602131	23.96355689	23.99108254	23.97287517	24.0145603	23.97828667	24.01836085	23.9884788	23.94800081	23.98552319	23.97848809	23.98905522	24.00344568	23.97955816	24.01818708	23.97524251	24.01016478	24.02145339	23.99088146	24.02183419	24.05116845	24.03075052	24.02232388	24.02232776	23.97417188	23.95456651	23.91122118	23.94982484	23.90875923	23.87848248	23.83410606	23.8482355	23.87390599	23.88136232	23.86868804	23.83164759	23.78518025	23.77242178	23.70282335	23.68041147	23.64443221	23.57346974	23.58616093	23.63882128	23.6174238	23.67478921	23.67214255	23.69950468	23.73646241	23.69235226	23.76979022	23.74214019	23.76077761	23.81020249	23.80385852	23.81526337	23.83096266	23.85785449	23.8654654	23.83893492	23.88461093	23.87777681	23.88404832	23.89246813	23.90569153	23.89341793	23.90741867	23.89201089	23.90731354	23.93448656	23.90126854	23.93518116	23.8789213	23.93776066	23.9201228	23.90386955	23.93003853	23.88798147	23.91288734	23.94520389	23.94453455
# Inverted_Repeats	24.11423662	24.1033751	24.08757856	24.09426607	24.09771062	24.09319811	24.10347168	24.10347389	24.11615649	24.09304843	24.11192682	24.09687067	24.09054268	24.09301938	24.10900626	24.09893244	24.09243581	24.08326324	24.10912867	24.12452023	24.07547323	24.11131391	24.0735086	24.13078516	24.12736048	24.09181887	24.12241918	24.11018559	24.0844525	24.13411778	24.11763425	24.08641551	24.11596229	24.08126051	24.09855754	24.13140553	24.07512056	24.10440207	24.11115831	24.11827996	24.14429383	24.09330777	24.11619318	24.1538606	24.2019341	24.04360124	24.11433768	24.07683489	24.02232531	23.99456056	24.14336431	24.19017915	24.19660236	24.15415159	24.14748226	24.07154248	23.92504392	24.05486734	24.12291169	24.050141	24.11593006	24.08429841	24.13672847	24.10229929	24.05566087	24.0587949	24.12579465	24.10219199	24.08706446	24.06503942	24.08796977	24.10474149	24.04729945	24.10933214	24.10988267	24.08807657	24.10424312	24.08209484	24.10918581	24.12383843	24.06126052	24.0807112	24.12684415	24.06240567	24.12555316	24.06995575	24.10269208	24.10519912	24.08216795	24.09287417	24.10870651	24.09801921	24.08031782	24.10199229	24.06836331	24.10121591	24.10292174	24.09275483	24.10635477	24.10779188
# Mirror_Repeats	23.95236663	23.94548006	23.95271625	23.94765985	23.94673281	23.95740523	23.96028564	23.92284948	23.93351877	23.91640983	23.92974368	23.9270896	23.94450152	23.95539853	23.96068457	23.93622992	23.93275037	23.90510675	23.9342036	23.93517753	23.95451926	23.96276529	23.96103343	23.9530467	23.94775615	23.95725489	23.96685041	23.9500877	23.94703082	23.96158954	23.92687633	23.98621796	23.97481834	23.96993235	23.96141511	23.99338512	23.97647206	23.98398938	23.98128382	23.99566536	24.03605703	24.05291907	24.04733268	24.01971937	23.95727781	23.77094244	23.83559687	23.93340349	23.98877474	24.02680965	24.025057	23.95126901	23.84942114	23.93048664	24.01696104	24.00684108	23.94643081	23.89115783	23.85646269	23.8215646	23.8681718	23.89428715	23.88696502	23.9012454	23.91001804	23.90135304	23.9119998	23.89944372	23.88489966	23.90956927	23.90378592	23.91679194	23.91613841	23.93618461	23.92640573	23.9387121	23.94025811	23.93648049	23.93264677	23.94264308	23.95171515	23.93760649	23.94184706	23.95486781	23.96130811	23.92064049	23.907232	23.90789045	23.94259872	23.91460787	23.92439354	23.95176406	23.93818592	23.94777466	23.87908096	23.92071566	23.93635521	23.93565522	23.96393936	23.94001052
# ...
depth=t((depth/matrix(depth[nrow(depth),],nrow=nrow(depth),ncol=100,byrow=TRUE))[-nrow(depth),])

load('IPD_forward_max10000_1_results_quantiles_5_25_50_75_95.RData')
plotSummary_reproducible(result_multi_quantile,groupby='feature',filenames=paste0('IPD_forward_max10000_reproducible_results_quantiles_5_25_50_75_95_summary.pdf'),
                         reproducible_pval_matrix=reproducible_pval_matrix,only_significant=TRUE,
                         gaps_tests=c(7,11,17,37),align_lab='Motif center',xlab='',ylab='',ask=FALSE,cellwidth=5,cellheight=12,fontsize=12)
load('IPD_forward.RData')
pdf('boxplot_vs_control_forward_reproducible_results_quantiles_5_25_50_75_95.pdf',width=9,height=6)
for(i in seq_len(nrow(regionsFeatures_forward)-1)){
  plot_no_size_control(regionsFeatures_forward,reproducible_pval=reproducible_pval_matrix[[i]],depth=depth,
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

# Load mean depth in each position
depth=as.matrix(read.delim('depth_mean_Gquad.txt',
                           sep='\t',header=FALSE,row.names=1,col.names=0:100))
# depth_mean_Gquad.txt is a text file (tab separated) in the following format: id   depth_1   ...   depth_100
# with the average number of reads used to compute the IPD in wach of the 100 nucleotides
# GGGAGGGAGGGAGGG	20.70568928	20.89059081	21.08871851	21.04918033	20.78446389	20.91356674	20.95519126	20.99234973	20.74535519	20.93333333	20.9321663	20.92997812	20.74835886	20.85355191	20.91794311	20.90710383	20.53005464	20.80765027	20.8512035	20.8284153	20.55191257	20.72021858	20.81291028	20.88852459	20.39781421	20.67978142	20.80437158	20.91693989	20.36284153	20.70601093	20.76805252	20.95846995	20.30163934	20.69726776	20.67103825	20.90054645	20.33989071	21	20.42185792	20.80853392	20.24590164	21.01420765	20.21600877	22.0010929	20.93770492	16.46280088	21.26120219	21.96174863	20.9715847	16.28008753	21.25027322	21.87540984	20.86229508	16.03169399	21.29180328	21.65027322	20.45901639	15.90350877	21.75628415	20.84918033	20.37923497	20.84371585	21.12021858	21.24153005	20.87759563	21.16502732	21.24043716	21.41530055	21.08415301	20.86338798	21.35956284	21.26338798	21.15409836	20.93442623	21.26885246	21.36174863	21.12349727	20.94535519	21.32932166	21.25027322	21.08087432	21.12910284	21.32677596	21.25464481	21.13442623	21.0788609	21.21530055	21.3715847	21.12472648	21.20787746	21.21772429	21.21639344	21.12240437	21.23632385	21.22319475	21.21202186	21.23986857	21.23304158	21.2702407	21.32459016
# GGGAGGGAGGGAGGGAGGG	20.46470588	20.77286136	20.64411765	20.55325444	20.60294118	20.56764706	20.58823529	20.47647059	20.54572271	20.64117647	20.50294118	20.39411765	20.55588235	20.63823529	20.35294118	20.62241888	20.55294118	20.63235294	20.35882353	20.50882353	20.33823529	20.66764706	20.10294118	20.59705882	20.58284024	20.70294118	20.18823529	20.55294118	20.37758112	20.62352941	19.85588235	20.58112094	20.38529412	20.67352941	19.70294118	20.80825959	19.93235294	20.44247788	19.86764706	20.61061947	19.78761062	21.50294118	20.53235294	16.30678466	20.85294118	21.50882353	20.56176471	16.13017751	20.81764706	21.47352941	20.44705882	16.33727811	20.82890855	21.36176471	20.46470588	15.93823529	20.93529412	21.07647059	20.07058824	15.72941176	21.33235294	20.38235294	19.79411765	20.52352941	20.63823529	20.74117647	20.44411765	20.79056047	20.99410029	21.02647059	20.61764706	20.58529412	20.86470588	21.02941176	20.62941176	20.75294118	20.93529412	20.94411765	20.69117647	20.61470588	20.86470588	20.93235294	20.63529412	20.80588235	20.89117647	21.09144543	20.9	20.7079646	20.88823529	20.99705015	20.85250737	20.77941176	20.85	21.07352941	20.96165192	20.84411765	21.08529412	21.10588235	20.75882353	21.03529412
# GGGAGGGAGGGAGGGAGGGAGGG	21.31404959	21.37190083	21.33884298	21.51239669	21.09917355	20.81818182	21.12396694	21.40495868	20.91735537	21.00826446	21.24793388	21.52066116	21.11570248	21.3553719	21.09917355	21.17355372	20.88429752	21.30578512	21.20661157	21.09917355	20.20661157	21.42975207	21.07438017	20.96694215	20.71900826	21.23966942	21.17355372	21.23140496	20.50413223	21.52066116	21.00826446	21.39669421	20.38016529	21.33057851	20.98347107	21.38842975	20.42975207	21.25619835	20.34710744	22.24793388	20.95041322	17.15702479	21.52066116	22.04132231	21.15702479	16.88429752	21.51239669	22.04132231	21.05785124	16.78512397	21.52892562	22.05785124	21.00826446	16.76859504	21.50413223	22	20.98347107	16.4214876	21.55371901	21.32231405	20	15.5785124	21.47107438	20.60330579	19.90909091	21.21487603	20.98347107	21.38016529	20.68595041	21.19834711	21.36363636	21.55371901	20.90082645	21.37190083	21.17355372	21.28099174	21.23140496	21.41322314	21.27272727	21.37190083	21.11570248	21.23966942	21.2892562	21.38016529	20.76859504	21.10743802	21.14876033	21.3553719	21.05785124	21.00826446	21.32231405	21.4214876	20.90082645	20.68595041	21.01652893	21.45454545	21.14049587	20.71900826	21.01652893	21.23966942
depth=t((depth/matrix(depth[nrow(depth),],nrow=nrow(depth),ncol=100,byrow=TRUE))[-nrow(depth),])

load('IPD_forward_Gquad.RData')
pdf('boxplot_vs_control_forward_Gquad_reproducible_results_quantiles_5_25_50_75_95.pdf',width=9,height=6)
for(i in seq_len(nrow(regionsFeatures_forward_GQuadruplexMotifs_subset)-1)){
  plot_no_size_control(regionsFeatures_forward_GQuadruplexMotifs_subset,reproducible_pval=reproducible_pval_matrix[[i]],depth=depth,
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
              only_significant=FALSE,gaps_tests=c(7,11,17,37),align_lab='Motif center',ask=FALSE,xlab='Coordinates',cellwidth=10,cellheight=10)
  
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
              only_significant=FALSE,gaps_tests=c(7,11,17,37),align_lab='Motif center',ask=FALSE,xlab='Coordinates',cellwidth=10,cellheight=10)
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

# Load mean depth in each position
depth=as.matrix(read.delim('reverse/depth_mean.txt',
                           sep='\t',header=FALSE,row.names=1,col.names=0:100))
# depth_mean.txt is a text file (tab separated) in the following format: id   depth_1   ...   depth_100
# with the average number of reads used to compute the IPD in wach of the 100 nucleotides
# Direct_Repeats	24.07356475	24.05188188	24.0728078	24.04015518	24.05769788	24.01441598	24.01398419	24.06720134	24.04421089	24.04415299	24.02012393	24.03552197	24.02312858	24.03398564	24.01835394	24.01609448	24.01464841	23.98292033	24.0125351	24.00602131	23.96355689	23.99108254	23.97287517	24.0145603	23.97828667	24.01836085	23.9884788	23.94800081	23.98552319	23.97848809	23.98905522	24.00344568	23.97955816	24.01818708	23.97524251	24.01016478	24.02145339	23.99088146	24.02183419	24.05116845	24.03075052	24.02232388	24.02232776	23.97417188	23.95456651	23.91122118	23.94982484	23.90875923	23.87848248	23.83410606	23.8482355	23.87390599	23.88136232	23.86868804	23.83164759	23.78518025	23.77242178	23.70282335	23.68041147	23.64443221	23.57346974	23.58616093	23.63882128	23.6174238	23.67478921	23.67214255	23.69950468	23.73646241	23.69235226	23.76979022	23.74214019	23.76077761	23.81020249	23.80385852	23.81526337	23.83096266	23.85785449	23.8654654	23.83893492	23.88461093	23.87777681	23.88404832	23.89246813	23.90569153	23.89341793	23.90741867	23.89201089	23.90731354	23.93448656	23.90126854	23.93518116	23.8789213	23.93776066	23.9201228	23.90386955	23.93003853	23.88798147	23.91288734	23.94520389	23.94453455
# Inverted_Repeats	24.11423662	24.1033751	24.08757856	24.09426607	24.09771062	24.09319811	24.10347168	24.10347389	24.11615649	24.09304843	24.11192682	24.09687067	24.09054268	24.09301938	24.10900626	24.09893244	24.09243581	24.08326324	24.10912867	24.12452023	24.07547323	24.11131391	24.0735086	24.13078516	24.12736048	24.09181887	24.12241918	24.11018559	24.0844525	24.13411778	24.11763425	24.08641551	24.11596229	24.08126051	24.09855754	24.13140553	24.07512056	24.10440207	24.11115831	24.11827996	24.14429383	24.09330777	24.11619318	24.1538606	24.2019341	24.04360124	24.11433768	24.07683489	24.02232531	23.99456056	24.14336431	24.19017915	24.19660236	24.15415159	24.14748226	24.07154248	23.92504392	24.05486734	24.12291169	24.050141	24.11593006	24.08429841	24.13672847	24.10229929	24.05566087	24.0587949	24.12579465	24.10219199	24.08706446	24.06503942	24.08796977	24.10474149	24.04729945	24.10933214	24.10988267	24.08807657	24.10424312	24.08209484	24.10918581	24.12383843	24.06126052	24.0807112	24.12684415	24.06240567	24.12555316	24.06995575	24.10269208	24.10519912	24.08216795	24.09287417	24.10870651	24.09801921	24.08031782	24.10199229	24.06836331	24.10121591	24.10292174	24.09275483	24.10635477	24.10779188
# Mirror_Repeats	23.95236663	23.94548006	23.95271625	23.94765985	23.94673281	23.95740523	23.96028564	23.92284948	23.93351877	23.91640983	23.92974368	23.9270896	23.94450152	23.95539853	23.96068457	23.93622992	23.93275037	23.90510675	23.9342036	23.93517753	23.95451926	23.96276529	23.96103343	23.9530467	23.94775615	23.95725489	23.96685041	23.9500877	23.94703082	23.96158954	23.92687633	23.98621796	23.97481834	23.96993235	23.96141511	23.99338512	23.97647206	23.98398938	23.98128382	23.99566536	24.03605703	24.05291907	24.04733268	24.01971937	23.95727781	23.77094244	23.83559687	23.93340349	23.98877474	24.02680965	24.025057	23.95126901	23.84942114	23.93048664	24.01696104	24.00684108	23.94643081	23.89115783	23.85646269	23.8215646	23.8681718	23.89428715	23.88696502	23.9012454	23.91001804	23.90135304	23.9119998	23.89944372	23.88489966	23.90956927	23.90378592	23.91679194	23.91613841	23.93618461	23.92640573	23.9387121	23.94025811	23.93648049	23.93264677	23.94264308	23.95171515	23.93760649	23.94184706	23.95486781	23.96130811	23.92064049	23.907232	23.90789045	23.94259872	23.91460787	23.92439354	23.95176406	23.93818592	23.94777466	23.87908096	23.92071566	23.93635521	23.93565522	23.96393936	23.94001052
# ...
depth=t((depth/matrix(depth[nrow(depth),],nrow=nrow(depth),ncol=100,byrow=TRUE))[-nrow(depth),])

load('IPD_reverse_max10000_1_results_mean.RData')
plotSummary_reproducible(result_mean,groupby='feature',filenames=paste0('IPD_reverse_max10000_reproducible_results_mean_summary.pdf'),
                         reproducible_pval_matrix=reproducible_pval_matrix,only_significant=TRUE,
                         gaps_tests=c(7,11,17,37),align_lab='Motif center',xlab='',ylab='',ask=FALSE,cellwidth=5,cellheight=12,fontsize=12)
load('IPD_reverse.RData')
pdf('boxplot_vs_control_reverse_reproducible_results_mean.pdf',width=9,height=6)
for(i in seq_len(nrow(regionsFeatures_reverse)-1)){
  plot_no_size_control(regionsFeatures_reverse,reproducible_pval=reproducible_pval_matrix[[i]],depth=depth,
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
              only_significant=FALSE,gaps_tests=c(7,11,17,37),align_lab='Motif center',ask=FALSE,xlab='Coordinates',cellwidth=10,cellheight=10)
  
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
              only_significant=FALSE,gaps_tests=c(7,11,17,37),align_lab='Motif center',ask=FALSE,xlab='Coordinates',cellwidth=10,cellheight=10)
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

# Load mean depth in each position
depth=as.matrix(read.delim('reverse/depth_mean.txt',
                           sep='\t',header=FALSE,row.names=1,col.names=0:100))
# depth_mean.txt is a text file (tab separated) in the following format: id   depth_1   ...   depth_100
# with the average number of reads used to compute the IPD in wach of the 100 nucleotides
# Direct_Repeats	24.07356475	24.05188188	24.0728078	24.04015518	24.05769788	24.01441598	24.01398419	24.06720134	24.04421089	24.04415299	24.02012393	24.03552197	24.02312858	24.03398564	24.01835394	24.01609448	24.01464841	23.98292033	24.0125351	24.00602131	23.96355689	23.99108254	23.97287517	24.0145603	23.97828667	24.01836085	23.9884788	23.94800081	23.98552319	23.97848809	23.98905522	24.00344568	23.97955816	24.01818708	23.97524251	24.01016478	24.02145339	23.99088146	24.02183419	24.05116845	24.03075052	24.02232388	24.02232776	23.97417188	23.95456651	23.91122118	23.94982484	23.90875923	23.87848248	23.83410606	23.8482355	23.87390599	23.88136232	23.86868804	23.83164759	23.78518025	23.77242178	23.70282335	23.68041147	23.64443221	23.57346974	23.58616093	23.63882128	23.6174238	23.67478921	23.67214255	23.69950468	23.73646241	23.69235226	23.76979022	23.74214019	23.76077761	23.81020249	23.80385852	23.81526337	23.83096266	23.85785449	23.8654654	23.83893492	23.88461093	23.87777681	23.88404832	23.89246813	23.90569153	23.89341793	23.90741867	23.89201089	23.90731354	23.93448656	23.90126854	23.93518116	23.8789213	23.93776066	23.9201228	23.90386955	23.93003853	23.88798147	23.91288734	23.94520389	23.94453455
# Inverted_Repeats	24.11423662	24.1033751	24.08757856	24.09426607	24.09771062	24.09319811	24.10347168	24.10347389	24.11615649	24.09304843	24.11192682	24.09687067	24.09054268	24.09301938	24.10900626	24.09893244	24.09243581	24.08326324	24.10912867	24.12452023	24.07547323	24.11131391	24.0735086	24.13078516	24.12736048	24.09181887	24.12241918	24.11018559	24.0844525	24.13411778	24.11763425	24.08641551	24.11596229	24.08126051	24.09855754	24.13140553	24.07512056	24.10440207	24.11115831	24.11827996	24.14429383	24.09330777	24.11619318	24.1538606	24.2019341	24.04360124	24.11433768	24.07683489	24.02232531	23.99456056	24.14336431	24.19017915	24.19660236	24.15415159	24.14748226	24.07154248	23.92504392	24.05486734	24.12291169	24.050141	24.11593006	24.08429841	24.13672847	24.10229929	24.05566087	24.0587949	24.12579465	24.10219199	24.08706446	24.06503942	24.08796977	24.10474149	24.04729945	24.10933214	24.10988267	24.08807657	24.10424312	24.08209484	24.10918581	24.12383843	24.06126052	24.0807112	24.12684415	24.06240567	24.12555316	24.06995575	24.10269208	24.10519912	24.08216795	24.09287417	24.10870651	24.09801921	24.08031782	24.10199229	24.06836331	24.10121591	24.10292174	24.09275483	24.10635477	24.10779188
# Mirror_Repeats	23.95236663	23.94548006	23.95271625	23.94765985	23.94673281	23.95740523	23.96028564	23.92284948	23.93351877	23.91640983	23.92974368	23.9270896	23.94450152	23.95539853	23.96068457	23.93622992	23.93275037	23.90510675	23.9342036	23.93517753	23.95451926	23.96276529	23.96103343	23.9530467	23.94775615	23.95725489	23.96685041	23.9500877	23.94703082	23.96158954	23.92687633	23.98621796	23.97481834	23.96993235	23.96141511	23.99338512	23.97647206	23.98398938	23.98128382	23.99566536	24.03605703	24.05291907	24.04733268	24.01971937	23.95727781	23.77094244	23.83559687	23.93340349	23.98877474	24.02680965	24.025057	23.95126901	23.84942114	23.93048664	24.01696104	24.00684108	23.94643081	23.89115783	23.85646269	23.8215646	23.8681718	23.89428715	23.88696502	23.9012454	23.91001804	23.90135304	23.9119998	23.89944372	23.88489966	23.90956927	23.90378592	23.91679194	23.91613841	23.93618461	23.92640573	23.9387121	23.94025811	23.93648049	23.93264677	23.94264308	23.95171515	23.93760649	23.94184706	23.95486781	23.96130811	23.92064049	23.907232	23.90789045	23.94259872	23.91460787	23.92439354	23.95176406	23.93818592	23.94777466	23.87908096	23.92071566	23.93635521	23.93565522	23.96393936	23.94001052
# ...
depth=t((depth/matrix(depth[nrow(depth),],nrow=nrow(depth),ncol=100,byrow=TRUE))[-nrow(depth),])

load('IPD_reverse_max10000_1_results_median.RData')
plotSummary_reproducible(result_median,groupby='feature',filenames=paste0('IPD_reverse_max10000_reproducible_results_median_summary.pdf'),
                         reproducible_pval_matrix=reproducible_pval_matrix,only_significant=TRUE,
                         gaps_tests=c(7,11,17,37),align_lab='Motif center',xlab='',ylab='',ask=FALSE,cellwidth=5,cellheight=12,fontsize=12)
load('IPD_reverse.RData')
pdf('boxplot_vs_control_reverse_reproducible_results_median.pdf',width=9,height=6)
for(i in seq_len(nrow(regionsFeatures_reverse)-1)){
  plot_no_size_control(regionsFeatures_reverse,reproducible_pval=reproducible_pval_matrix[[i]],depth=depth,
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
              only_significant=FALSE,gaps_tests=c(7,11,17,37),align_lab='Motif center',ask=FALSE,xlab='Coordinates',cellwidth=10,cellheight=10)
  
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
              only_significant=FALSE,gaps_tests=c(7,11,17,37),align_lab='Motif center',ask=FALSE,xlab='Coordinates',cellwidth=10,cellheight=10)
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

# Load mean depth in each position
depth=as.matrix(read.delim('reverse/depth_mean.txt',
                           sep='\t',header=FALSE,row.names=1,col.names=0:100))
# depth_mean.txt is a text file (tab separated) in the following format: id   depth_1   ...   depth_100
# with the average number of reads used to compute the IPD in wach of the 100 nucleotides
# Direct_Repeats	24.07356475	24.05188188	24.0728078	24.04015518	24.05769788	24.01441598	24.01398419	24.06720134	24.04421089	24.04415299	24.02012393	24.03552197	24.02312858	24.03398564	24.01835394	24.01609448	24.01464841	23.98292033	24.0125351	24.00602131	23.96355689	23.99108254	23.97287517	24.0145603	23.97828667	24.01836085	23.9884788	23.94800081	23.98552319	23.97848809	23.98905522	24.00344568	23.97955816	24.01818708	23.97524251	24.01016478	24.02145339	23.99088146	24.02183419	24.05116845	24.03075052	24.02232388	24.02232776	23.97417188	23.95456651	23.91122118	23.94982484	23.90875923	23.87848248	23.83410606	23.8482355	23.87390599	23.88136232	23.86868804	23.83164759	23.78518025	23.77242178	23.70282335	23.68041147	23.64443221	23.57346974	23.58616093	23.63882128	23.6174238	23.67478921	23.67214255	23.69950468	23.73646241	23.69235226	23.76979022	23.74214019	23.76077761	23.81020249	23.80385852	23.81526337	23.83096266	23.85785449	23.8654654	23.83893492	23.88461093	23.87777681	23.88404832	23.89246813	23.90569153	23.89341793	23.90741867	23.89201089	23.90731354	23.93448656	23.90126854	23.93518116	23.8789213	23.93776066	23.9201228	23.90386955	23.93003853	23.88798147	23.91288734	23.94520389	23.94453455
# Inverted_Repeats	24.11423662	24.1033751	24.08757856	24.09426607	24.09771062	24.09319811	24.10347168	24.10347389	24.11615649	24.09304843	24.11192682	24.09687067	24.09054268	24.09301938	24.10900626	24.09893244	24.09243581	24.08326324	24.10912867	24.12452023	24.07547323	24.11131391	24.0735086	24.13078516	24.12736048	24.09181887	24.12241918	24.11018559	24.0844525	24.13411778	24.11763425	24.08641551	24.11596229	24.08126051	24.09855754	24.13140553	24.07512056	24.10440207	24.11115831	24.11827996	24.14429383	24.09330777	24.11619318	24.1538606	24.2019341	24.04360124	24.11433768	24.07683489	24.02232531	23.99456056	24.14336431	24.19017915	24.19660236	24.15415159	24.14748226	24.07154248	23.92504392	24.05486734	24.12291169	24.050141	24.11593006	24.08429841	24.13672847	24.10229929	24.05566087	24.0587949	24.12579465	24.10219199	24.08706446	24.06503942	24.08796977	24.10474149	24.04729945	24.10933214	24.10988267	24.08807657	24.10424312	24.08209484	24.10918581	24.12383843	24.06126052	24.0807112	24.12684415	24.06240567	24.12555316	24.06995575	24.10269208	24.10519912	24.08216795	24.09287417	24.10870651	24.09801921	24.08031782	24.10199229	24.06836331	24.10121591	24.10292174	24.09275483	24.10635477	24.10779188
# Mirror_Repeats	23.95236663	23.94548006	23.95271625	23.94765985	23.94673281	23.95740523	23.96028564	23.92284948	23.93351877	23.91640983	23.92974368	23.9270896	23.94450152	23.95539853	23.96068457	23.93622992	23.93275037	23.90510675	23.9342036	23.93517753	23.95451926	23.96276529	23.96103343	23.9530467	23.94775615	23.95725489	23.96685041	23.9500877	23.94703082	23.96158954	23.92687633	23.98621796	23.97481834	23.96993235	23.96141511	23.99338512	23.97647206	23.98398938	23.98128382	23.99566536	24.03605703	24.05291907	24.04733268	24.01971937	23.95727781	23.77094244	23.83559687	23.93340349	23.98877474	24.02680965	24.025057	23.95126901	23.84942114	23.93048664	24.01696104	24.00684108	23.94643081	23.89115783	23.85646269	23.8215646	23.8681718	23.89428715	23.88696502	23.9012454	23.91001804	23.90135304	23.9119998	23.89944372	23.88489966	23.90956927	23.90378592	23.91679194	23.91613841	23.93618461	23.92640573	23.9387121	23.94025811	23.93648049	23.93264677	23.94264308	23.95171515	23.93760649	23.94184706	23.95486781	23.96130811	23.92064049	23.907232	23.90789045	23.94259872	23.91460787	23.92439354	23.95176406	23.93818592	23.94777466	23.87908096	23.92071566	23.93635521	23.93565522	23.96393936	23.94001052
# ...
depth=t((depth/matrix(depth[nrow(depth),],nrow=nrow(depth),ncol=100,byrow=TRUE))[-nrow(depth),])

load('IPD_reverse_max10000_1_results_quantiles_5_25_50_75_95.RData')
plotSummary_reproducible(result_multi_quantile,groupby='feature',filenames=paste0('IPD_reverse_max10000_reproducible_results_quantiles_5_25_50_75_95_summary.pdf'),
                         reproducible_pval_matrix=reproducible_pval_matrix,only_significant=TRUE,
                         gaps_tests=c(7,11,17,37),align_lab='Motif center',xlab='',ylab='',ask=FALSE,cellwidth=5,cellheight=12,fontsize=12)
load('IPD_reverse.RData')
pdf('boxplot_vs_control_reverse_reproducible_results_quantiles_5_25_50_75_95.pdf',width=9,height=6)
for(i in seq_len(nrow(regionsFeatures_reverse)-1)){
  plot_no_size_control(regionsFeatures_reverse,reproducible_pval=reproducible_pval_matrix[[i]],depth=depth,
                       type='boxplot',probs=c(0.05,0.25,0.5,0.75,0.95),size_perc=TRUE,position=TRUE,average=FALSE,cex=1.4,cex.axis=1.4,cex.lab=1.4,cex.main=2.5,
                       id_regions_subset=c('Control',idRegions(regionsFeatures_reverse)[i]),col=c('red','blue'),ask=FALSE,xlab='Coordinates')
}
dev.off()




