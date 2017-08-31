require(IWTomics)
source('IWTomics_modified_functions.r')



#######################
##### PLUS STRAND #####
#######################

### Complete datasets ###

datasets=read.table("datasets.txt",sep="\t",header=TRUE,stringsAsFactors=FALSE)
# datasets.txt is a text files (tab separated) in the following format:
# id	              name	              regionFile
# Direct_Repeats	  Direct  Repeats	    DirectRepeats_filtered
# Inverted_Repeats	Inverted Repeats	  InvertedRepeats_filtered
# Mirror_Repeats	  Mirror Repeats	    MirrorRepeats_filtered
# GQuadPlus	        G Quadruplex plus	  GQuadPlus_filtered
# GQuadMinus	      G Quadruplex minus	GQuadMinus_filtered
# ...

features_datasets=read.table("features_datasets_forward.txt",sep="\t",header=TRUE,stringsAsFactors=FALSE)
# features_datasets_forward.txt is a text files (tab separated) in the following format:
# id	        name	          Direct_Repeats	        Inverted_Repeats	        Mirror_Repeats	        GQuadPlus	          GQuadMinus          ...
# IPD_Forward	IPD plus strand	DirectRepeats_filtered	InvertedRepeats_filtered	MirrorRepeats_filtered	GQuadPlus_filtered	GQuadMinus_filtered ...

# The files should be in the following format (tab separated): chr   start   end   IPD_1   IPD_2   ...   IPD_100
# chr1	264949	265048	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	2.102	0.513	0.294	NA	0.557	1.19	0.242	3.062	1.303	0.215	0.496	1.106	0.731	0.032	0.921	1.065	NA	NA	NA	0.475	0.849	0.622	NA
# chr1	344598	344697	1.058	1.329	1.55	0.617	0.72	0.747	0.85	0.246	1.25	1.51	0.802	1.167	0.119	1.87	0.505	1.766	0.76	1.222	2.462	12.283	1.308	1.552	1.945	1.22	0.608	0.973	0.647	0.684	0.243	0.879	0.599	0.254	0.236	1.455	0.238	0.629	0.438	3.804	1.493	0.731	1.066	0.344	3.036	1.859	0.606	0.681	0.236	1.389	0.646	1.121	1.049	1.863	2.862	0.409	1.954	0.655	0.821	1.133	0.449	0.99	0.712	0.942	2.117	1.987	0.728	0.078	0.528	2.171	3.097	0.763	0.501	0.432	0.672	0.446	0.833	0.338	0.824	0.796	6.577	1.039	0.762	0.346	0.396	1.009	3.847	0.838	0.552	0.779	0.345	1.017	4.794	0.433	1.058	0.448	4.361	0.463	0.141	0.555	0.314	0.299
# ...

# Read regions and IPD
regionsFeatures_forward=IWTomicsData(datasets$regionFile,features_datasets[,datasets$id],'center',
                                     datasets$id,datasets$name,
                                     features_datasets$id,features_datasets$name,
                                     path='forward/',
                                     start.are.0based=FALSE)

# Add feature length, number of nucleotides before the start of the feature and composition for every region
regionsFeatures_forward@regions=GRangesList(lapply(regionsFeatures_forward@regions, 
                                                   function(region){
                                                     mcols(region)=data.frame(length=rep(NA,length(region)),NAbefore=rep(NA,length(region)),
                                                                              A=rep(NA,length(region)),T=rep(NA,length(region)),G=rep(NA,length(region)),C=rep(NA,length(region)),
                                                                              AA=rep(NA,length(region)),AT=rep(NA,length(region)),AG=rep(NA,length(region)),AC=rep(NA,length(region)),
                                                                              TA=rep(NA,length(region)),TT=rep(NA,length(region)),TG=rep(NA,length(region)),TC=rep(NA,length(region)),
                                                                              GA=rep(NA,length(region)),GT=rep(NA,length(region)),GG=rep(NA,length(region)),GC=rep(NA,length(region)),
                                                                              CA=rep(NA,length(region)),CT=rep(NA,length(region)),CG=rep(NA,length(region)),CC=rep(NA,length(region)))
                                                     return(region)
                                                   }))
for(id_region in idRegions(regionsFeatures_forward)[1:7]){ # non-microsatellite region datasets
  tmp=read.delim(paste0('forward/',regionsFeatures_forward@metadata$region_datasets[id_region,'file'],'_length'),
                 sep='\t',header=FALSE)
  # The files should have the same name as above, but with "_length" at the end
  # The files should be in the following format (tab separated): chr   start   end   length
  # chr1	264949	265048	25
  # chr1	344598	344697	25
  # ...
  names(tmp)=c('chr','start','end','length')
  tmp$NAbefore=round((100.1-tmp$length)/2)
  tmp=makeGRangesFromDataFrame(tmp,keep.extra.columns=TRUE)
  index=match(regions(regionsFeatures_forward)[[id_region]],tmp)
  if(sum(is.na(index)))
    stop('there are NA...')
  tmp2=read.delim(paste0('forward/',regionsFeatures_forward@metadata$region_datasets[id_region,'file'],'_composition'),
                  sep='\t',header=FALSE)
  # The files should have the same name as above, but with "_composition" at the end
  # The files should be in the following format (tab separated): chr   start   end   A   T   G   C
  # chr1	264949	265048	35	35	18	12
  # chr1	344598	344697	35	36	13	16
  # ...
  names(tmp2)=c('chr','start','end','A','T','G','C')
  tmp2=makeGRangesFromDataFrame(tmp2,keep.extra.columns=TRUE)
  index2=match(regions(regionsFeatures_forward)[[id_region]],tmp2)
  if(sum(is.na(index2)))
    stop('there are NA...')
  tmp3=read.delim(paste0('forward/',regionsFeatures_forward@metadata$region_datasets[id_region,'file'],'_composition_di'),
                  sep='\t',header=FALSE)
  # The files should have the same name as above, but with "_composition" at the end
  # The files should be in the following format (tab separated): chr   start   end   AA   AT   AG   AC ...
  # chr1	264949	265048	5 7 11	1  ...
  # chr1	344598	344697	3 10  3	13  ...
  # ...
  names(tmp3)=c('chr','start','end','AA','AT','AG','AC','TA','TT','TG','TC','GA','GT','GG','GC','CA','CT','CG','CC')
  tmp3=makeGRangesFromDataFrame(tmp3,keep.extra.columns=TRUE)
  index3=match(regions(regionsFeatures_forward)[[id_region]],tmp3)
  if(sum(is.na(index3)))
    stop('there are NA...')
  mcols(regionsFeatures_forward@regions[[id_region]])<-cbind(mcols(tmp[index]),mcols(tmp2[index2]),mcols(tmp3[index3]))
}
for(id_region in idRegions(regionsFeatures_forward)[8:83]){ # microsatellite region datasets
  tmp=read.delim(paste0('forward/',regionsFeatures_forward@metadata$region_datasets[id_region,'file'],'_length'),
                 sep='\t',header=FALSE)
  # The files should have the same name as above, but with "_length" at the end
  # The files should be in the following format (tab separated): chr   start   end   length   NAbefore
  # chr10	1628182	1628281	12	46
  # chr10	2286394	2286493	12	44
  # ...
  names(tmp)=c('chr','start','end','length','NAbefore')
  tmp=makeGRangesFromDataFrame(tmp,keep.extra.columns=TRUE)
  index=match(regions(regionsFeatures_forward)[[id_region]],tmp)
  if(sum(is.na(index)))
    stop('there are NA...')
  tmp2=read.delim(paste0('forward/',regionsFeatures_forward@metadata$region_datasets[id_region,'file'],'_composition'),
                  sep='\t',header=FALSE)
  # The files should have the same name as above, but with "_composition" at the end
  # The files should be in the following format (tab separated): chr   start   end   A   T   G   C
  # chr1	3838648	3838747	46	20	18	16
  # chr1	4003529	4003628	37	26	15	22
  # ...
  names(tmp2)=c('chr','start','end','A','T','G','C')
  tmp2=makeGRangesFromDataFrame(tmp2,keep.extra.columns=TRUE)
  index2=match(regions(regionsFeatures_forward)[[id_region]],tmp2)
  if(sum(is.na(index2)))
    stop('there are NA...')
  tmp3=read.delim(paste0('forward/',regionsFeatures_forward@metadata$region_datasets[id_region,'file'],'_composition_di'),
                  sep='\t',header=FALSE)
  # The files should have the same name as above, but with "_composition" at the end
  # The files should be in the following format (tab separated): chr   start   end   AA   AT   AG   AC ...
  # chr1	264949	265048	5 7 11	1  ...
  # chr1	344598	344697	3 10  3	13  ...
  # ...
  names(tmp3)=c('chr','start','end','AA','AT','AG','AC','TA','TT','TG','TC','GA','GT','GG','GC','CA','CT','CG','CC')
  tmp3=makeGRangesFromDataFrame(tmp3,keep.extra.columns=TRUE)
  index3=match(regions(regionsFeatures_forward)[[id_region]],tmp3)
  if(sum(is.na(index3)))
    stop('there are NA...')
  mcols(regionsFeatures_forward@regions[[id_region]])<-cbind(mcols(tmp[index]),mcols(tmp2[index2]),mcols(tmp3[index3]))
}
id_region=idRegions(regionsFeatures_forward)[84] # controls
tmp2=read.delim(paste0('forward/',regionsFeatures_forward@metadata$region_datasets[id_region,'file'],'_composition'),
                sep='\t',header=FALSE)
tmp2=cbind(tmp2[,1:3],NA,NA,tmp2[,4:7])
names(tmp2)=c('chr','start','end','length','NAbefore','A','T','G','C')
tmp2=makeGRangesFromDataFrame(tmp2,keep.extra.columns=TRUE)
index2=match(regions(regionsFeatures_forward)[[id_region]],tmp2)
if(sum(is.na(index2)))
  stop('there are NA...')
tmp3=read.delim(paste0('forward/',regionsFeatures_forward@metadata$region_datasets[id_region,'file'],'_composition_di'),
                sep='\t',header=FALSE)
names(tmp3)=c('chr','start','end','AA','AT','AG','AC','TA','TT','TG','TC','GA','GT','GG','GC','CA','CT','CG','CC')
tmp3=makeGRangesFromDataFrame(tmp3,keep.extra.columns=TRUE)
index3=match(regions(regionsFeatures_forward)[[id_region]],tmp3)
if(sum(is.na(index3)))
  stop('there are NA...')
mcols(regionsFeatures_forward@regions[[id_region]])<-cbind(mcols(tmp2[index2]),mcols(tmp3[index3]))

# Add microsatellite sequences
tmp=read.delim('forward/microsat_letters.txt',
               sep='\t',header=FALSE,row.names=1,col.names=0:100)
# microsat_letters.txt is a text file (tab separated) in the following format: id   letter_1   ...   letter_100
# with the aligned microsatellite sequences in the 100 bp windows (NA in the position where there is no microsatellite)
# ACTn	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	T	A	C	T	A	C	T	A	C	T	A	C	T	A	C	T	A	C	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA
# AGCn	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	C	A	G	C	A	G	C	A	G	C	A	G	C	A	G	C	A	G	C	A	G	C	A	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA
# ...
tmp=rbind(matrix(NA,ncol=100,nrow=7),as.matrix(tmp),matrix(NA,ncol=100,nrow=1))
row.names(tmp)=idRegions(regionsFeatures_forward)
tmp=as.data.frame(tmp,stringsAsFactors=FALSE)
mcols(regionsFeatures_forward@regions) <- tmp

save(regionsFeatures_forward,file='IPD_forward.RData')
write.table(regionsFeatures_forward@metadata$region_datasets,file='size_forward.txt',quote=FALSE,sep='\t')

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

# Pointwise boxplots
pdf('boxplot_vs_control_forward_position.pdf',width=9,height=6)
for(i in seq_len(nrow(regionsFeatures_forward)-1)){
  plot_no_size_control(regionsFeatures_forward,depth=depth,type='boxplot',probs=c(0.05,0.25,0.5,0.75,0.95),size_perc=TRUE,position=TRUE,average=FALSE,cex=1.4,cex.axis=1.4,cex.lab=1.4,cex.main=2.5,
                       id_regions_subset=c('Control',idRegions(regionsFeatures_forward)[i]),col=c('red','blue'),ask=FALSE,xlab='Coordinates')
}
dev.off()


# Subsample regions, maximum 10.000 for each region dataset
# 10 random subsamples
# to test for differences using IWT
# see IWT_IPD.r
maxsize=10000
regions_maxsize=rownames(regionsFeatures_forward_complete@metadata$region_datasets)[regionsFeatures_forward_complete@metadata$region_datasets$size>maxsize]
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
  
  save(regionsFeatures_forward,file=paste0('IPD_forward_max10000_',i,'.RData'))
}




### G quadruplex datasets ###

datasets=read.table("datasets_Gquad.txt",sep="\t",header=TRUE,stringsAsFactors=FALSE)
# datasets_Gquad.txt is a text files (tab separated) in the following format:
# id	                    name	            regionFile
# GGGAGGGAGGGAGGG	        G Quadruplex plus	GGGAGGGAGGGAGGG_filtered
# GGGAGGGAGGGAGGGAGGG	    G Quadruplex plus	GGGAGGGAGGGAGGGAGGG_filtered
# GGGAGGGAGGGAGGGAGGGAGGG	G Quadruplex plus	GGGAGGGAGGGAGGGAGGGAGGG_filtered
# GGGAGGGAGGTGGGGGGG	    G Quadruplex plus	GGGAGGGAGGTGGGGGGG_filtered
# GGGAGGGAGGTGGGGGGGG	    G Quadruplex plus	GGGAGGGAGGTGGGGGGGG_filtered
# ...

features_datasets=read.table("features_datasets_forward.txt",sep="\t",header=TRUE,stringsAsFactors=FALSE)
# features_datasets_forward.txt is a text files (tab separated) in the following format:
# id          name            GGGAGGGAGGGAGGG           GGGAGGGAGGGAGGGAGGG           GGGAGGGAGGGAGGGAGGGAGGG           GGGAGGGAGGTGGGGGGG          GGGAGGGAGGTGGGGGGGG           ...
# IPD_Forward IPD plus strand GGGAGGGAGGGAGGG_filtered  GGGAGGGAGGGAGGGAGGG_filtered  GGGAGGGAGGGAGGGAGGGAGGG_filtered  GGGAGGGAGGTGGGGGGG_filtered GGGAGGGAGGTGGGGGGGG_filtered  ...

# The files should be in the following format (tab separated): chr   start   end   IPD_1   IPD_2   ...   IPD_100
# chr1	2406011	2406110	0.744	0.985	0.508	0.857	0.252	0.989	0.266	0.337	0.275	1.037	0.983	1.483	0.548	0.385	0.625	4.175	2.082	0.83	0.36	1.258	1.118	0.215	0.585	1.344	1.428	1.244	0.669	0.744	0.51	1.426	1.969	0.23	0.795	0.815	1.923	1.055	0.399	1.474	0.438	0.747	1.162	0.879	1.261	1.149	0.738	0.956	0.618	1.806	0.707	0.245	0.562	0.343	1.236	0.705	2.086	0.743	1.03	1.438	0.318	0.905	0.883	1.179	1.237	0.658	0.717	0.974	0.197	0.547	0.151	0.756	0.613	1.184	1.409	0.908	0.938	1.386	1.612	0.317	3.36	0.427	1.072	0.799	0.698	0.423	0.961	0.904	0.641	0.463	2.669	0.371	0.182	1.754	3.926	0.405	0.454	1.894	3.026	0.707	0.725	1.135
# chr1	4997697	4997796	0.254	1.651	1.039	0.701	0.45	0.832	2.124	0.329	0.329	1.743	0.586	2.662	0.805	0.702	1.043	0.54	0.451	1.087	3.136	1.006	0.39	1.744	0.576	0.381	0.783	0.254	0.717	0.125	0.911	1.002	1.152	0.433	1.139	2.235	2.026	0.692	0.592	0.598	0.676	0.272	0.6	0.638	0.686	1.715	1.765	0.889	0.647	0.606	0.814	0.267	0.699	0.766	0.818	1.019	0.57	1.184	1.286	NA	2.505	0.759	0.596	0.502	0.286	0.651	0.336	0.526	1.673	0.616	0.66	0.335	1.041	1.594	0.886	0.559	0.688	0.555	0.354	0.446	0.507	0.445	0.772	0.607	0.758	0.669	1.631	0.616	1.357	0.926	0.692	1.72	0.849	0.468	0.734	0.459	0.8	0.235	0.887	0.774	1.63	1.093
# ...

# Read regions and IPD
regionsFeatures_forward_GQuadruplexMotifs_subset=IWTomicsData(datasets$regionFile,features_datasets[,datasets$id],'center',
                                                              datasets$id,datasets$name,
                                                              features_datasets$id,features_datasets$name,
                                                              path='forward/',
                                                              start.are.0based=FALSE)

# Add feature length and number of nucleotides before the start of the feature
regionsFeatures_forward_GQuadruplexMotifs_subset@regions=GRangesList(lapply(regionsFeatures_forward_GQuadruplexMotifs_subset@regions,
                                                                            function(region){
                                                                              mcols(region)=data.frame(length=rep(NA,length(region)),NAbefore=rep(NA,length(region)))
                                                                              return(region)
                                                                            }))
for(id_region in idRegions(regionsFeatures_forward_GQuadruplexMotifs_subset)[1:20]){
  tmp=read.delim# The files should have the same name as above, but with "_length" at the end
  # The files should be in the following format (tab separated): chr   start   end   length
  # chr1	2406011	2406110	15
  # chr1	4997697	4997796	15
  # ...
  names(tmp)=c('chr','start','end','length')
  tmp$NAbefore=round((100.1-tmp$length)/2)
  tmp=makeGRangesFromDataFrame(tmp,keep.extra.columns=TRUE)
  index=match(regions(regionsFeatures_forward_GQuadruplexMotifs_subset)[[id_region]],tmp)
  if(sum(is.na(index)))
    stop('there are NA...')
  mcols(regionsFeatures_forward_GQuadruplexMotifs_subset@regions[[id_region]])<-mcols(tmp[index])
}

# Add G quadruplex sequences
tmp=read.delim('forward/Gquad_letters.txt',
               sep='\t',header=FALSE,row.names=1,col.names=0:100)
# Gquad_letters.txt is a text file (tab separated) in the following format: id   letter_1   ...   letter_100
# GGGAGGGAGGGAGGG	    NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	G	G	G	A	G	G	G	A	G	G	G	A	G	G	G	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA
# GGGAGGGAGGGAGGGAGGG	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	G	G	G	A	G	G	G	A	G	G	G	A	G	G	G	A	G	G	G	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA
# ...
tmp=rbind(as.matrix(tmp),matrix(NA,ncol=100,nrow=1))
row.names(tmp)=idRegions(regionsFeatures_forward_GQuadruplexMotifs_subset)
tmp=as.data.frame(tmp,stringsAsFactors=FALSE)
mcols(regionsFeatures_forward_GQuadruplexMotifs_subset@regions) <- tmp

save(regionsFeatures_forward_GQuadruplexMotifs_subset,file='IPD_forward_Gquad.RData')

# Load mean depth in each position
depth=as.matrix(read.delim('depth_mean_Gquad.txt',
                           sep='\t',header=FALSE,row.names=1,col.names=0:100))
# depth_mean_Gquad.txt is a text file (tab separated) in the following format: id   depth_1   ...   depth_100
# with the average number of reads used to compute the IPD in wach of the 100 nucleotides
# GGGAGGGAGGGAGGG	20.70568928	20.89059081	21.08871851	21.04918033	20.78446389	20.91356674	20.95519126	20.99234973	20.74535519	20.93333333	20.9321663	20.92997812	20.74835886	20.85355191	20.91794311	20.90710383	20.53005464	20.80765027	20.8512035	20.8284153	20.55191257	20.72021858	20.81291028	20.88852459	20.39781421	20.67978142	20.80437158	20.91693989	20.36284153	20.70601093	20.76805252	20.95846995	20.30163934	20.69726776	20.67103825	20.90054645	20.33989071	21	20.42185792	20.80853392	20.24590164	21.01420765	20.21600877	22.0010929	20.93770492	16.46280088	21.26120219	21.96174863	20.9715847	16.28008753	21.25027322	21.87540984	20.86229508	16.03169399	21.29180328	21.65027322	20.45901639	15.90350877	21.75628415	20.84918033	20.37923497	20.84371585	21.12021858	21.24153005	20.87759563	21.16502732	21.24043716	21.41530055	21.08415301	20.86338798	21.35956284	21.26338798	21.15409836	20.93442623	21.26885246	21.36174863	21.12349727	20.94535519	21.32932166	21.25027322	21.08087432	21.12910284	21.32677596	21.25464481	21.13442623	21.0788609	21.21530055	21.3715847	21.12472648	21.20787746	21.21772429	21.21639344	21.12240437	21.23632385	21.22319475	21.21202186	21.23986857	21.23304158	21.2702407	21.32459016
# GGGAGGGAGGGAGGGAGGG	20.46470588	20.77286136	20.64411765	20.55325444	20.60294118	20.56764706	20.58823529	20.47647059	20.54572271	20.64117647	20.50294118	20.39411765	20.55588235	20.63823529	20.35294118	20.62241888	20.55294118	20.63235294	20.35882353	20.50882353	20.33823529	20.66764706	20.10294118	20.59705882	20.58284024	20.70294118	20.18823529	20.55294118	20.37758112	20.62352941	19.85588235	20.58112094	20.38529412	20.67352941	19.70294118	20.80825959	19.93235294	20.44247788	19.86764706	20.61061947	19.78761062	21.50294118	20.53235294	16.30678466	20.85294118	21.50882353	20.56176471	16.13017751	20.81764706	21.47352941	20.44705882	16.33727811	20.82890855	21.36176471	20.46470588	15.93823529	20.93529412	21.07647059	20.07058824	15.72941176	21.33235294	20.38235294	19.79411765	20.52352941	20.63823529	20.74117647	20.44411765	20.79056047	20.99410029	21.02647059	20.61764706	20.58529412	20.86470588	21.02941176	20.62941176	20.75294118	20.93529412	20.94411765	20.69117647	20.61470588	20.86470588	20.93235294	20.63529412	20.80588235	20.89117647	21.09144543	20.9	20.7079646	20.88823529	20.99705015	20.85250737	20.77941176	20.85	21.07352941	20.96165192	20.84411765	21.08529412	21.10588235	20.75882353	21.03529412
# GGGAGGGAGGGAGGGAGGGAGGG	21.31404959	21.37190083	21.33884298	21.51239669	21.09917355	20.81818182	21.12396694	21.40495868	20.91735537	21.00826446	21.24793388	21.52066116	21.11570248	21.3553719	21.09917355	21.17355372	20.88429752	21.30578512	21.20661157	21.09917355	20.20661157	21.42975207	21.07438017	20.96694215	20.71900826	21.23966942	21.17355372	21.23140496	20.50413223	21.52066116	21.00826446	21.39669421	20.38016529	21.33057851	20.98347107	21.38842975	20.42975207	21.25619835	20.34710744	22.24793388	20.95041322	17.15702479	21.52066116	22.04132231	21.15702479	16.88429752	21.51239669	22.04132231	21.05785124	16.78512397	21.52892562	22.05785124	21.00826446	16.76859504	21.50413223	22	20.98347107	16.4214876	21.55371901	21.32231405	20	15.5785124	21.47107438	20.60330579	19.90909091	21.21487603	20.98347107	21.38016529	20.68595041	21.19834711	21.36363636	21.55371901	20.90082645	21.37190083	21.17355372	21.28099174	21.23140496	21.41322314	21.27272727	21.37190083	21.11570248	21.23966942	21.2892562	21.38016529	20.76859504	21.10743802	21.14876033	21.3553719	21.05785124	21.00826446	21.32231405	21.4214876	20.90082645	20.68595041	21.01652893	21.45454545	21.14049587	20.71900826	21.01652893	21.23966942
depth=t((depth/matrix(depth[nrow(depth),],nrow=nrow(depth),ncol=100,byrow=TRUE))[-nrow(depth),])

# Pointwise boxplots
pdf('boxplot_vs_control_Gquad_forward_position.pdf',width=9,height=6)
for(i in seq_len(nrow(regionsFeatures_forward_GQuadruplexMotifs_subset)-1)){
  plot_no_size_control(regionsFeatures_forward_GQuadruplexMotifs_subset,depth=depth,type='boxplot',probs=c(0.05,0.25,0.5,0.75,0.95),size_perc=TRUE,position=TRUE,average=FALSE,cex=1.4,cex.axis=1.4,cex.lab=1.4,cex.main=2.5,
                       id_regions_subset=c('Control',idRegions(regionsFeatures_forward_GQuadruplexMotifs_subset)[i]),col=c('red','blue'),ask=FALSE,xlab='Coordinates')
}
dev.off()


# Subsample regions, maximum 10.000 for each region dataset
# 10 random subsamples
# to test for differences using IWT
# see IWT_IPD.r
regionsFeatures_forward_GQuadruplexMotifs_subset@metadata$region_datasets['Control','size']=10000
regionsFeatures_forward_GQuadruplexMotifs_subset@length_features$IPD_Forward$Control=rep(100,10000)
for(i in 1:10){
  load(paste0('IPD_forward_max10000_',i,'.RData'))
  regionsFeatures_forward_GQuadruplexMotifs_subset@features$IPD_Forward$Control=
    regionsFeatures_forward@features$IPD_Forward$Control
  validObject(regionsFeatures_forward_GQuadruplexMotifs_subset)
  
  save(regionsFeatures_forward_GQuadruplexMotifs_subset,file=paste0('IPD_forward_Gquad_max10000_',i,'.RData'))
}




### A-phased repeats datasets ###

datasets=read.table("datasets_Aphased.txt",sep="\t",header=TRUE,stringsAsFactors=FALSE)
# datasets_Aphased.txt is a text files (tab separated) in the following format:
# id	                      name	            regionFile
# AAACCTTGGACAATACCTGGCTTT	A-phased repeats	AAACCTTGGACAATACCTGGCTTT_filtered
# AAAGAAGCTAAAAACCTTGAAAAAA	A-phased repeats	AAAGAAGCTAAAAACCTTGAAAAAA_filtered
# AAAGAGAATAAAATACCTAGGAAT	A-phased repeats	AAAGAGAATAAAATACCTAGGAAT_filtered
# AAAGAGAATAAAATACCTAGGAATA	A-phased repeats	AAAGAGAATAAAATACCTAGGAATA_filtered
# AAAGAGAATAAAATACTTAGGAAT	A-phased repeats	AAAGAGAATAAAATACTTAGGAAT_filtered
# ...

features_datasets=read.table("features_datasets_forward.txt",sep="\t",header=TRUE,stringsAsFactors=FALSE)
# features_datasets_forward.txt is a text files (tab separated) in the following format:
# id          name            AAACCTTGGACAATACCTGGCTTT           AAAGAAGCTAAAAACCTTGAAAAAA           AAAGAGAATAAAATACCTAGGAAT           AAAGAGAATAAAATACCTAGGAATA          AAAGAGAATAAAATACTTAGGAAT           ...
# IPD_Forward IPD plus strand AAACCTTGGACAATACCTGGCTTT_filtered  AAAGAAGCTAAAAACCTTGAAAAAA_filtered  AAAGAGAATAAAATACCTAGGAAT_filtered  AAAGAGAATAAAATACCTAGGAATA_filtered AAAGAGAATAAAATACTTAGGAAT_filtered  ...

# The files should be in the following format (tab separated): chr   start   end   IPD_1   IPD_2   ...   IPD_100
# chr1	5634277	  5634376 	1.554	0.441	0.375	0.569	1.198	0.495	2.259	0.297	1.481	0.627	1.069	0.93	1.114	1.802	0.519	1.146	1.968	0.324	0.589	0.291	0.615	2.821	0.683	2.058	1.738	0.532	0.8	0.422	1.107	4.071	0.495	1.191	2.022	1.66	1.908	1.777	2.077	0.651	1.091	0.929	0.428	2.213	2.49	1.367	0.247	1.357	0.989	0.848	2.427	0.807	1.117	1.079	0.787	1.089	0.601	0.297	1.355	2.033	0.57	0.767	0.884	0.722	0.616	1.105	0.767	3.065	1.197	0.932	0.425	1.265	0.854	5.791	0.199	0.575	1.026	0.764	0.632	0.525	0.584	0.373	1.835	0.88	1.811	0.791	3.424	0.505	0.537	0.873	0.658	0.78	0.55	0.557	2.537	0.696	0.628	0.359	1.025	0.662	0.482	2.141
# chr1	15462409	15462508	0.594	0.291	0.552	0.592	NA	0.801	1.852	0.215	0.654	0.589	0.273	0.442	0.509	0.881	0.164	0.722	0.28	0.344	0.137	0.13	0.396	1.162	0.225	2.067	1.544	0.458	0.311	0.278	0.293	1.272	0.838	0.661	1.976	0.113	0.961	0.944	0.894	0.435	0.389	0.251	0.39	0.22	0.918	0.914	NA	0.886	0.286	0.39	1.203	1.16	0.476	0.323	0.273	0.747	0.909	0.229	0.452	1.651	0.715	0.292	0.186	NA	0.676	NA	0.57	0.584	0.514	0.397	0.322	0.557	0.277	3.123	0.474	0.492	0.286	0.604	0.632	NA	1.396	0.545	0.862	0.417	0.462	3.567	1.701	0.515	0.161	0.405	0.679	0.435	0.376	0.5	0.303	0.292	0.531	1.251	0.123	NA	NA	0.295
# ...

# Read regions and IPD
regionsFeatures_forward_APhasedRepeats_subset=IWTomicsData(datasets$regionFile,features_datasets[,datasets$id],'center',
                                                           datasets$id,datasets$name,
                                                           features_datasets$id,features_datasets$name,
                                                           path='forward/',
                                                           start.are.0based=FALSE)

# Add feature length and number of nucleotides before the start of the feature
regionsFeatures_forward_APhasedRepeats_subset@regions=GRangesList(lapply(regionsFeatures_forward_APhasedRepeats_subset@regions,
                                                                         function(region){
                                                                           mcols(region)=data.frame(length=rep(NA,length(region)),NAbefore=rep(NA,length(region)))
                                                                           return(region)
                                                                         }))
for(id_region in idRegions(regionsFeatures_forward_APhasedRepeats_subset)[1:20]){
  tmp=read.delim(paste0('forward/',regionsFeatures_forward_APhasedRepeats_subset@metadata$region_datasets[id_region,'file'],'_length'),
                 sep='\t',header=FALSE)
  # The files should have the same name as above, but with "_length" at the end
  # The files should be in the following format (tab separated): chr   start   end   length
  # chr1	5634277	  5634376 	24
  # chr1	15462409	15462508	24
  # ...
  names(tmp)=c('chr','start','end','length')
  tmp$NAbefore=round((100.1-tmp$length)/2)
  tmp=makeGRangesFromDataFrame(tmp,keep.extra.columns=TRUE)
  index=match(regions(regionsFeatures_forward_APhasedRepeats_subset)[[id_region]],tmp)
  if(sum(is.na(index)))
    stop('there are NA...')
  mcols(regionsFeatures_forward_APhasedRepeats_subset@regions[[id_region]])<-mcols(tmp[index])
}

# Add A-phased repeat sequences
tmp=read.delim('forward/Aphased_letters.txt',
               sep='\t',header=FALSE,row.names=1,col.names=0:100,colClasses='factor')
# Aphased_letters.txt is a text file (tab separated) in the following format: id   letter_1   ...   letter_100
# AAACCTTGGACAATACCTGGCTTT	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	A	A	A	C	C	T	T	G	G	A	C	A	A	T	A	C	C	T	G	G	C	T	T	T	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA
# AAAGAAGCTAAAAACCTTGAAAAAA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	A	A	A	G	A	A	G	C	T	A	A	A	A	A	C	C	T	T	G	A	A	A	A	A	A	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA
# ...
tmp=rbind(as.matrix(tmp),matrix(NA,ncol=100,nrow=1))
row.names(tmp)=idRegions(regionsFeatures_forward_APhasedRepeats_subset)
tmp=as.data.frame(tmp,stringsAsFactors=FALSE)
mcols(regionsFeatures_forward_APhasedRepeats_subset@regions) <- tmp

save(regionsFeatures_forward_APhasedRepeats_subset,file='IPD_forward_Aphased.RData')

# Pointwise boxplots
pdf('boxplot_vs_control_Aphased_forward_position.pdf',width=9,height=6)
for(i in seq_len(nrow(regionsFeatures_forward_APhasedRepeats_subset)-1)){
  plot_no_size_control(regionsFeatures_forward_APhasedRepeats_subset,type='boxplot',probs=c(0.05,0.25,0.5,0.75,0.95),size_perc=TRUE,position=TRUE,average=FALSE,cex=1.4,cex.axis=1.4,cex.lab=1.4,cex.main=2.5,
                       id_regions_subset=c('Control',idRegions(regionsFeatures_forward_APhasedRepeats_subset)[i]),col=c('red','blue'),ask=FALSE,xlab='Coordinates')
}
dev.off()


# Subsample regions, maximum 10.000 for each region dataset
# 10 random subsamples
# to test for differences using IWT
# see IWT_IPD.r
regionsFeatures_forward_APhasedRepeats_subset@metadata$region_datasets['Control','size']=10000
regionsFeatures_forward_APhasedRepeats_subset@length_features$IPD_Forward$Control=rep(100,10000)
for(i in 1:10){
  load(paste0('IPD_forward_max10000_',i,'.RData'))
  regionsFeatures_forward_APhasedRepeats_subset@features$IPD_Forward$Control=
    regionsFeatures_forward@features$IPD_Forward$Control
  validObject(regionsFeatures_forward_APhasedRepeats_subset)
  
  save(regionsFeatures_forward_APhasedRepeats_subset,file=paste0('IPD_forward_Aphased_max10000_',i,'.RData'))
}















#######################
##### MINUS STRAND ####
#######################

### Complete datasets ###

datasets=read.table("datasets.txt",sep="\t",header=TRUE,stringsAsFactors=FALSE)
# datasets.txt is a text files (tab separated) in the following format:
# id	              name	              regionFile
# Direct_Repeats	  Direct  Repeats	    DirectRepeats_filtered
# Inverted_Repeats	Inverted Repeats	  InvertedRepeats_filtered
# Mirror_Repeats	  Mirror Repeats	    MirrorRepeats_filtered
# GQuadPlus	        G Quadruplex plus	  GQuadPlus_filtered
# GQuadMinus	      G Quadruplex minus	GQuadMinus_filtered
# ...

features_datasets=read.table("features_datasets_reverse.txt",sep="\t",header=TRUE,stringsAsFactors=FALSE)
# features_datasets_reverse.txt is a text files (tab separated) in the following format:
# id	        name	          Direct_Repeats	        Inverted_Repeats	        Mirror_Repeats	        GQuadPlus	          GQuadMinus          ...
# IPD_Reverse	IPD minus strand	DirectRepeats_filtered	InvertedRepeats_filtered	MirrorRepeats_filtered	GQuadPlus_filtered	GQuadMinus_filtered ...

# The files should be in the following format (tab separated): chr   start   end   IPD_1   IPD_2   ...   IPD_100
# chr1	264949	265048	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	2.102	0.513	0.294	NA	0.557	1.19	0.242	3.062	1.303	0.215	0.496	1.106	0.731	0.032	0.921	1.065	NA	NA	NA	0.475	0.849	0.622	NA
# chr1	344598	344697	1.058	1.329	1.55	0.617	0.72	0.747	0.85	0.246	1.25	1.51	0.802	1.167	0.119	1.87	0.505	1.766	0.76	1.222	2.462	12.283	1.308	1.552	1.945	1.22	0.608	0.973	0.647	0.684	0.243	0.879	0.599	0.254	0.236	1.455	0.238	0.629	0.438	3.804	1.493	0.731	1.066	0.344	3.036	1.859	0.606	0.681	0.236	1.389	0.646	1.121	1.049	1.863	2.862	0.409	1.954	0.655	0.821	1.133	0.449	0.99	0.712	0.942	2.117	1.987	0.728	0.078	0.528	2.171	3.097	0.763	0.501	0.432	0.672	0.446	0.833	0.338	0.824	0.796	6.577	1.039	0.762	0.346	0.396	1.009	3.847	0.838	0.552	0.779	0.345	1.017	4.794	0.433	1.058	0.448	4.361	0.463	0.141	0.555	0.314	0.299
# ...

# Read regions and IPD
regionsFeatures_reverse=IWTomicsData(datasets$regionFile,features_datasets[,datasets$id],'center',
                                     datasets$id,datasets$name,
                                     features_datasets$id,features_datasets$name,
                                     path='reverse/',
                                     start.are.0based=FALSE)

# Add feature length, number of nucleotides before the start of the feature and composition for every region
regionsFeatures_reverse@regions=GRangesList(lapply(regionsFeatures_reverse@regions, 
                                                   function(region){
                                                     mcols(region)=data.frame(length=rep(NA,length(region)),NAbefore=rep(NA,length(region)),
                                                                              A=rep(NA,length(region)),T=rep(NA,length(region)),G=rep(NA,length(region)),C=rep(NA,length(region)),
                                                                              AA=rep(NA,length(region)),AT=rep(NA,length(region)),AG=rep(NA,length(region)),AC=rep(NA,length(region)),
                                                                              TA=rep(NA,length(region)),TT=rep(NA,length(region)),TG=rep(NA,length(region)),TC=rep(NA,length(region)),
                                                                              GA=rep(NA,length(region)),GT=rep(NA,length(region)),GG=rep(NA,length(region)),GC=rep(NA,length(region)),
                                                                              CA=rep(NA,length(region)),CT=rep(NA,length(region)),CG=rep(NA,length(region)),CC=rep(NA,length(region)))
                                                     return(region)
                                                   }))
for(id_region in idRegions(regionsFeatures_reverse)[1:7]){ # non-microsatellite region datasets
  tmp=read.delim(paste0('reverse/',regionsFeatures_reverse@metadata$region_datasets[id_region,'file'],'_length'),
                 sep='\t',header=FALSE)
  # The files should have the same name as above, but with "_length" at the end
  # The files should be in the following format (tab separated): chr   start   end   length
  # chr1	264949	265048	25
  # chr1	344598	344697	25
  # ...
  names(tmp)=c('chr','start','end','length')
  tmp$NAbefore=round((100.1-tmp$length)/2)
  tmp=makeGRangesFromDataFrame(tmp,keep.extra.columns=TRUE)
  index=match(regions(regionsFeatures_reverse)[[id_region]],tmp)
  if(sum(is.na(index)))
    stop('there are NA...')
  tmp2=read.delim(paste0('reverse/',regionsFeatures_reverse@metadata$region_datasets[id_region,'file'],'_composition'),
                  sep='\t',header=FALSE)
  # The files should have the same name as above, but with "_composition" at the end
  # The files should be in the following format (tab separated): chr   start   end   A   T   G   C
  # chr1	264949	265048	35	35	18	12
  # chr1	344598	344697	35	36	13	16
  # ...
  names(tmp2)=c('chr','start','end','A','T','G','C')
  tmp2=makeGRangesFromDataFrame(tmp2,keep.extra.columns=TRUE)
  index2=match(regions(regionsFeatures_reverse)[[id_region]],tmp2)
  if(sum(is.na(index2)))
    stop('there are NA...')
  tmp3=read.delim(paste0('reverse/',regionsFeatures_reverse@metadata$region_datasets[id_region,'file'],'_composition_di'),
                  sep='\t',header=FALSE)
  # The files should have the same name as above, but with "_composition" at the end
  # The files should be in the following format (tab separated): chr   start   end   AA   AT   AG   AC ...
  # chr1	264949	265048	5 7 11	1  ...
  # chr1	344598	344697	3 10  3	13  ...
  # ...
  names(tmp3)=c('chr','start','end','AA','AT','AG','AC','TA','TT','TG','TC','GA','GT','GG','GC','CA','CT','CG','CC')
  tmp3=makeGRangesFromDataFrame(tmp3,keep.extra.columns=TRUE)
  index3=match(regions(regionsFeatures_reverse)[[id_region]],tmp3)
  if(sum(is.na(index3)))
    stop('there are NA...')
  mcols(regionsFeatures_reverse@regions[[id_region]])<-cbind(mcols(tmp[index]),mcols(tmp2[index2]),mcols(tmp3[index3]))
}
for(id_region in idRegions(regionsFeatures_reverse)[8:83]){ # microsatellite region datasets
  tmp=read.delim(paste0('reverse/',regionsFeatures_reverse@metadata$region_datasets[id_region,'file'],'_length'),
                 sep='\t',header=FALSE)
  # The files should have the same name as above, but with "_length" at the end
  # The files should be in the following format (tab separated): chr   start   end   length   NAbefore
  # chr10	1628182	1628281	12	46
  # chr10	2286394	2286493	12	44
  # ...
  names(tmp)=c('chr','start','end','length','NAbefore')
  tmp=makeGRangesFromDataFrame(tmp,keep.extra.columns=TRUE)
  index=match(regions(regionsFeatures_reverse)[[id_region]],tmp)
  if(sum(is.na(index)))
    stop('there are NA...')
  tmp2=read.delim(paste0('reverse/',regionsFeatures_reverse@metadata$region_datasets[id_region,'file'],'_composition'),
                  sep='\t',header=FALSE)
  # The files should have the same name as above, but with "_composition" at the end
  # The files should be in the following format (tab separated): chr   start   end   A   T   G   C
  # chr1	3838648	3838747	46	20	18	16
  # chr1	4003529	4003628	37	26	15	22
  # ...
  names(tmp2)=c('chr','start','end','A','T','G','C')
  tmp2=makeGRangesFromDataFrame(tmp2,keep.extra.columns=TRUE)
  index2=match(regions(regionsFeatures_reverse)[[id_region]],tmp2)
  if(sum(is.na(index2)))
    stop('there are NA...')
  tmp3=read.delim(paste0('reverse/',regionsFeatures_reverse@metadata$region_datasets[id_region,'file'],'_composition_di'),
                  sep='\t',header=FALSE)
  # The files should have the same name as above, but with "_composition" at the end
  # The files should be in the following format (tab separated): chr   start   end   AA   AT   AG   AC ...
  # chr1	264949	265048	5 7 11	1  ...
  # chr1	344598	344697	3 10  3	13  ...
  # ...
  names(tmp3)=c('chr','start','end','AA','AT','AG','AC','TA','TT','TG','TC','GA','GT','GG','GC','CA','CT','CG','CC')
  tmp3=makeGRangesFromDataFrame(tmp3,keep.extra.columns=TRUE)
  index3=match(regions(regionsFeatures_reverse)[[id_region]],tmp3)
  if(sum(is.na(index3)))
    stop('there are NA...')
  mcols(regionsFeatures_reverse@regions[[id_region]])<-cbind(mcols(tmp[index]),mcols(tmp2[index2]),mcols(tmp3[index3]))
}
id_region=idRegions(regionsFeatures_reverse)[84] # controls
tmp2=read.delim(paste0('reverse/',regionsFeatures_reverse@metadata$region_datasets[id_region,'file'],'_composition'),
                sep='\t',header=FALSE)
tmp2=cbind(tmp2[,1:3],NA,NA,tmp2[,4:7])
names(tmp2)=c('chr','start','end','length','NAbefore','A','T','G','C')
tmp2=makeGRangesFromDataFrame(tmp2,keep.extra.columns=TRUE)
index2=match(regions(regionsFeatures_reverse)[[id_region]],tmp2)
if(sum(is.na(index2)))
  stop('there are NA...')
tmp3=read.delim(paste0('reverse/',regionsFeatures_reverse@metadata$region_datasets[id_region,'file'],'_composition_di'),
                sep='\t',header=FALSE)
names(tmp3)=c('chr','start','end','AA','AT','AG','AC','TA','TT','TG','TC','GA','GT','GG','GC','CA','CT','CG','CC')
tmp3=makeGRangesFromDataFrame(tmp3,keep.extra.columns=TRUE)
index3=match(regions(regionsFeatures_reverse)[[id_region]],tmp3)
if(sum(is.na(index3)))
  stop('there are NA...')
mcols(regionsFeatures_reverse@regions[[id_region]])<-cbind(mcols(tmp2[index2]),mcols(tmp3[index3]))

# Add microsatellite sequences
tmp=read.delim('reverse/microsat_letters.txt',
               sep='\t',header=FALSE,row.names=1,col.names=0:100)
# microsat_letters.txt is a text file (tab separated) in the following format: id   letter_1   ...   letter_100
# with the aligned microsatellite sequences in the 100 bp windows (NA in the position where there is no microsatellite)
# ACTn	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	T	A	C	T	A	C	T	A	C	T	A	C	T	A	C	T	A	C	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA
# AGCn	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	C	A	G	C	A	G	C	A	G	C	A	G	C	A	G	C	A	G	C	A	G	C	A	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA
# ...
tmp=rbind(matrix(NA,ncol=100,nrow=7),as.matrix(tmp),matrix(NA,ncol=100,nrow=1))
row.names(tmp)=idRegions(regionsFeatures_reverse)
tmp=as.data.frame(tmp,stringsAsFactors=FALSE)
mcols(regionsFeatures_reverse@regions) <- tmp

save(regionsFeatures_reverse,file='IPD_reverse.RData')
write.table(regionsFeatures_reverse@metadata$region_datasets,file='size_reverse.txt',quote=FALSE,sep='\t')

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

# Pointwise boxplots
pdf('boxplot_vs_control_reverse_position.pdf',width=9,height=6)
for(i in seq_len(nrow(regionsFeatures_reverse)-1)){
  plot_no_size_control(regionsFeatures_reverse,depth=depth,type='boxplot',probs=c(0.05,0.25,0.5,0.75,0.95),size_perc=TRUE,position=TRUE,average=FALSE,cex=1.4,cex.axis=1.4,cex.lab=1.4,cex.main=2.5,
                       id_regions_subset=c('Control',idRegions(regionsFeatures_reverse)[i]),col=c('red','blue'),ask=FALSE,xlab='Coordinates')
}
dev.off()


# Subsample regions, maximum 10.000 for each region dataset
# 10 random subsamples
# to test for differences using IWT
# see IWT_IPD.r
maxsize=10000
regions_maxsize=rownames(regionsFeatures_reverse_complete@metadata$region_datasets)[regionsFeatures_reverse_complete@metadata$region_datasets$size>maxsize]
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
  
  save(regionsFeatures_reverse,file=paste0('IPD_reverse_max10000_',i,'.RData'))
}
