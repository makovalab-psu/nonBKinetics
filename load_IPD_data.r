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
                                                                              A=rep(NA,length(region)),T=rep(NA,length(region)),G=rep(NA,length(region)),C=rep(NA,length(region)))
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
  mcols(regionsFeatures_forward@regions[[id_region]])<-cbind(mcols(tmp[index]),mcols(tmp2[index2]))
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
  mcols(regionsFeatures_forward@regions[[id_region]])<-cbind(mcols(tmp[index]),mcols(tmp2[index2]))
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
mcols(regionsFeatures_forward@regions[[id_region]])<-mcols(tmp2[index2])

# Add microsatellite sequences
tmp=read.delim('forward/microsat_letters.txt',
               sep='\t',header=FALSE,row.names=1,col.names=0:100)
# microsat_letters.txt is a text file (tab separated) in the following format: id   letter_1   ...   letter_100
# with the aligned microsatellite sequences in the 100 bp windows (NA in the position where there is no microsatellite)
# ACTn	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	T	A	C	T	A	C	T	A	C	T	A	C	T	A	C	T	A	C	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA
# AGCn	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	C	A	G	C	A	G	C	A	G	C	A	G	C	A	G	C	A	G	C	A	G	C	A	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA
tmp=rbind(matrix(NA,ncol=100,nrow=7),as.matrix(tmp),matrix(NA,ncol=100,nrow=1))
row.names(tmp)=idRegions(regionsFeatures_forward)
tmp=as.data.frame(tmp,stringsAsFactors=FALSE)
mcols(regionsFeatures_forward@regions) <- tmp

save(regionsFeatures_forward,file='IPD_forward.RData')
write.table(regionsFeatures_forward@metadata$region_datasets,file='size_forward.txt',quote=FALSE,sep='\t')

# Pointwise boxplots
pdf('boxplot_vs_control_forward_position.pdf',width=9,height=6)
for(i in seq_len(nrow(regionsFeatures_forward)-1)){
  plot_no_size_control(regionsFeatures_forward,type='boxplot',probs=c(0.05,0.25,0.5,0.75,0.95),size_perc=TRUE,position=TRUE,average=FALSE,cex=1.4,cex.axis=1.4,cex.lab=1.4,cex.main=2.5,
                       id_regions_subset=c('Control',idRegions(regionsFeatures_forward)[i]),col=c('red','blue'),ask=FALSE,xlab='Coordinates')
}
dev.off()


# Filter out regions with NA
# to analyze the effect of sequence composition on IPD
# see sequence_composition_IPD.r
regionsFeatures_forward_complete=regionsFeatures_forward
tot_NA=lapply(regionsFeatures_forward_complete@features$IPD_Forward,function(IPD) colSums(is.na(IPD)))
index=lapply(tot_NA,function(tot_NA) which(tot_NA==0))
size_filtered=unlist(lapply(index,length))
IPD_filtered=mapply(function(IPD,index) as.matrix(IPD[,index]),regionsFeatures_forward_complete@features$IPD_Forward,index)
regions_filtered=mapply(function(regions,index) regions[index],regionsFeatures_forward_complete@regions,index)
length_features_filtered=mapply(function(length,index) length[index],regionsFeatures_forward_complete@length_features$IPD_Forward,index)
regions_10=row.names(regionsFeatures_forward_complete@metadata$region_datasets)[size_filtered>=10]
regionsFeatures_forward@metadata$region_datasets=regionsFeatures_forward_complete@metadata$region_datasets[regions_10,]
regionsFeatures_forward@metadata$region_datasets$size=size_filtered[regions_10]
regionsFeatures_forward@regions=GRangesList(regions_filtered)[regions_10]
regionsFeatures_forward@metadata$feature_datasets=regionsFeatures_forward_complete@metadata$feature_datasets[,c('name',paste0('file_',regions_10),'resolution')]
regionsFeatures_forward@features$IPD_Forward=IPD_filtered[regions_10]
regionsFeatures_forward@length_features$IPD_Forward=length_features_filtered[regions_10]
validObject(regionsFeatures_forward)

# Add microsatellite sequences
mcols(regionsFeatures_forward@regions) <- tmp[idRegions(regionsFeatures_forward_complete),]

save(regionsFeatures_forward,file='IPD_forward_NA0.RData')
write.table(regionsFeatures_forward@metadata$region_datasets,file='size_forward_NA0.txt',quote=FALSE,sep='\t')

# Pointwise boxplots
pdf('boxplot_vs_control_forward_NA0_position.pdf',width=9,height=6)
for(i in seq_len(nrow(regionsFeatures_forward)-1)){
  plot_no_size_control(regionsFeatures_forward,type='boxplot',probs=c(0.05,0.25,0.5,0.75,0.95),size_perc=TRUE,position=TRUE,average=FALSE,cex=1.4,cex.axis=1.4,cex.lab=1.4,cex.main=2.5,
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
tmp=rbind(as.matrix(tmp),matrix(NA,ncol=100,nrow=1))
row.names(tmp)=idRegions(regionsFeatures_forward_GQuadruplexMotifs_subset)
tmp=as.data.frame(tmp,stringsAsFactors=FALSE)
mcols(regionsFeatures_forward_GQuadruplexMotifs_subset@regions) <- tmp

save(regionsFeatures_forward_GQuadruplexMotifs_subset,file='IPD_forward_Gquad.RData')

# Pointwise boxplots
pdf('boxplot_vs_control_Gquad_forward_position.pdf',width=9,height=6)
for(i in seq_len(nrow(regionsFeatures_forward_GQuadruplexMotifs_subset)-1)){
  plot_no_size_control(regionsFeatures_forward_GQuadruplexMotifs_subset,type='boxplot',probs=c(0.05,0.25,0.5,0.75,0.95),size_perc=TRUE,position=TRUE,average=FALSE,cex=1.4,cex.axis=1.4,cex.lab=1.4,cex.main=2.5,
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
                                                                              A=rep(NA,length(region)),T=rep(NA,length(region)),G=rep(NA,length(region)),C=rep(NA,length(region)))
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
  mcols(regionsFeatures_reverse@regions[[id_region]])<-cbind(mcols(tmp[index]),mcols(tmp2[index2]))
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
  mcols(regionsFeatures_reverse@regions[[id_region]])<-cbind(mcols(tmp[index]),mcols(tmp2[index2]))
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
mcols(regionsFeatures_reverse@regions[[id_region]])<-mcols(tmp2[index2])

# Add microsatellite sequences
tmp=read.delim('reverse/microsat_letters.txt',
               sep='\t',header=FALSE,row.names=1,col.names=0:100)
# microsat_letters.txt is a text file (tab separated) in the following format: id   letter_1   ...   letter_100
# with the aligned microsatellite sequences in the 100 bp windows (NA in the position where there is no microsatellite)
# ACTn	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	T	A	C	T	A	C	T	A	C	T	A	C	T	A	C	T	A	C	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA
# AGCn	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	C	A	G	C	A	G	C	A	G	C	A	G	C	A	G	C	A	G	C	A	G	C	A	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA
tmp=rbind(matrix(NA,ncol=100,nrow=7),as.matrix(tmp),matrix(NA,ncol=100,nrow=1))
row.names(tmp)=idRegions(regionsFeatures_reverse)
tmp=as.data.frame(tmp,stringsAsFactors=FALSE)
mcols(regionsFeatures_reverse@regions) <- tmp

save(regionsFeatures_reverse,file='IPD_reverse.RData')
write.table(regionsFeatures_reverse@metadata$region_datasets,file='size_reverse.txt',quote=FALSE,sep='\t')

# Pointwise boxplots
pdf('boxplot_vs_control_reverse_position.pdf',width=9,height=6)
for(i in seq_len(nrow(regionsFeatures_reverse)-1)){
  plot_no_size_control(regionsFeatures_reverse,type='boxplot',probs=c(0.05,0.25,0.5,0.75,0.95),size_perc=TRUE,position=TRUE,average=FALSE,cex=1.4,cex.axis=1.4,cex.lab=1.4,cex.main=2.5,
                       id_regions_subset=c('Control',idRegions(regionsFeatures_reverse)[i]),col=c('red','blue'),ask=FALSE,xlab='Coordinates')
}
dev.off()


# Filter out regions with NA
# to analyze the effect of sequence composition on IPD
# see sequence_composition_IPD.r
regionsFeatures_reverse_complete=regionsFeatures_reverse
tot_NA=lapply(regionsFeatures_reverse_complete@features$IPD_Reverse,function(IPD) colSums(is.na(IPD)))
index=lapply(tot_NA,function(tot_NA) which(tot_NA==0))
size_filtered=unlist(lapply(index,length))
IPD_filtered=mapply(function(IPD,index) as.matrix(IPD[,index]),regionsFeatures_reverse_complete@features$IPD_Reverse,index)
regions_filtered=mapply(function(regions,index) regions[index],regionsFeatures_reverse_complete@regions,index)
length_features_filtered=mapply(function(length,index) length[index],regionsFeatures_reverse_complete@length_features$IPD_Reverse,index)
regions_10=row.names(regionsFeatures_reverse_complete@metadata$region_datasets)[size_filtered>=10]
regionsFeatures_reverse@metadata$region_datasets=regionsFeatures_reverse_complete@metadata$region_datasets[regions_10,]
regionsFeatures_reverse@metadata$region_datasets$size=size_filtered[regions_10]
regionsFeatures_reverse@regions=GRangesList(regions_filtered)[regions_10]
regionsFeatures_reverse@metadata$feature_datasets=regionsFeatures_reverse_complete@metadata$feature_datasets[,c('name',paste0('file_',regions_10),'resolution')]
regionsFeatures_reverse@features$IPD_Reverse=IPD_filtered[regions_10]
regionsFeatures_reverse@length_features$IPD_Reverse=length_features_filtered[regions_10]
validObject(regionsFeatures_reverse)

# Add microsatellite sequences
mcols(regionsFeatures_reverse@regions) <- tmp[idRegions(regionsFeatures_reverse_complete),]

save(regionsFeatures_reverse,file='IPD_reverse_NA0.RData')
write.table(regionsFeatures_reverse@metadata$region_datasets,file='size_reverse_NA0.txt',quote=FALSE,sep='\t')

# Pointwise boxplot
pdf('boxplot_vs_control_reverse_NA0_position.pdf',width=9,height=6)
for(i in seq_len(nrow(regionsFeatures_reverse)-1)){
  plot_no_size_control(regionsFeatures_reverse,type='boxplot',probs=c(0.05,0.25,0.5,0.75,0.95),size_perc=TRUE,position=TRUE,average=FALSE,cex=1.4,cex.axis=1.4,cex.lab=1.4,cex.main=2.5,
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


