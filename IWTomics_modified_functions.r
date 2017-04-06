plot_no_size_control <- function(x,type='boxplot',method='pearson',
                                 N_regions=pmin(lengthRegions(x),ifelse(type=='curves',10,ifelse(type=='pairs',1000,+Inf))),
                                 probs=c(0.25,0.5,0.75),average=TRUE,size=TRUE,
                                 size_perc=FALSE, # NEW plot the sample size in each position as percentage of the total sample size
                                 position=FALSE, # NEW plot the feature positions
                                 reproducible_pval=NULL, # NEW matrix of reproducible p-values (TRUE or FALSE). The TRUE are plotted as black squares under the pointwise boxplot
                                 scale_threshold=NULL, # NEW scale_threshold for reproducible_pval
                                 lengths=NULL, # NEW list with features length to be plotted separetely, for analysis of different lengths
                                 zero_line=TRUE, # NEK plot a horizontal line at y=0
                                 id_regions_subset=idRegions(x),id_features_subset=idFeatures(x),
                                 log_scale=FALSE,log_shift=0,col=1+seq_along(id_regions_subset),
                                 plot=TRUE,ask=TRUE,xlab='Windows',ylim=NULL,...){
  # MODIFIED FROM BIOCONDUCTOR PACKAGE IWTomics
  # FUNCTION plot.IWTomicsData
  
  if(sum(!(id_regions_subset %in% idRegions(x))))
    stop('invalid id_regions_subset. The region datasets provided are not listed in x.')
  if(sum(!(id_features_subset %in% idFeatures(x))))
    stop('invalid id_features_subset. The features provided are not listed in x.')
  if(!(type %in% c('curves','boxplot','pairs','pairsSmooth')))
    stop('invalid plot type \'',type,'\'. Available types are \'curves\', \'boxplot\', \'pairs\' and \'pairsSmooth\'.')
  if(!(method %in% c('pearson','kendall','spearman')))
    stop('invalid method type \'',method,'\'. Available methods are \'pearson\', \'kendall\' and \'spearman\'.')
  
  ### start NEW ###
  if(!is.null(reproducible_pval)){
    if(is.null(scale_threshold))
      scale_threshold=nrow(reproducible_pval)
    reproducible_pval=reproducible_pval[ncol(reproducible_pval)-scale_threshold+1,]
  }
  ### end NEW ###
  x=x[id_regions_subset,id_features_subset]
  ### start NEW ###
  if(!is.null(lengths)){
    new_ids=paste0(setdiff(id_regions_subset,'Control'),lengths)
    reg=x@regions[[setdiff(id_regions_subset,'Control')]]
    feat=x@features[[id_features_subset]][[setdiff(id_regions_subset,'Control')]]
    len=x@length_features[[id_features_subset]][[setdiff(id_regions_subset,'Control')]]
    index=lapply(lengths,function(len) which(reg$length==len))
    if(100 %in% lengths)
      index[[which(lengths==100)]]=which(reg$length>=100)
    names(index)=new_ids
    size_new=unlist(lapply(index,length))
    feat_new=lapply(index,function(index) as.matrix(feat[,index]))
    reg_new=lapply(index,function(index) reg[index])
    len_new=lapply(index,function(index) len[index])
    reg_datasets_new=as.data.frame(matrix(x@metadata$region_datasets[setdiff(id_regions_subset,'Control'),],nrow=length(lengths),ncol=3,byrow=TRUE),row.names=new_ids)
    names(reg_datasets_new)=names(x@metadata$region_datasets)
    reg_datasets_new$size=size_new
    x@metadata$region_datasets=rbind(x@metadata$region_datasets,reg_datasets_new)
    x@regions=c(x@regions,GRangesList(reg_new))
    mcols(x@regions[new_ids])=mcols(x@regions[setdiff(id_regions_subset,'Control')])
    for(new_id in new_ids){
      NAbefore=pmax(unique(mcols(x@regions[[new_id]])$NAbefore),0)
      len=unique(pmin(unique(mcols(x@regions[[new_id]])$length),100))
      if(length(len)>1)
        stop('Something wrong with the length.')
      mcols(x@regions[new_id])[seq_len(min(NAbefore))]=NA
      mcols(x@regions[new_id])[100-seq_len(100-max(NAbefore)-len)+1]=NA
    }
    x@metadata$feature_datasets=x@metadata$feature_datasets[,c('name',paste0('file_',c(id_regions_subset,rep(setdiff(id_regions_subset,'Control'),length(lengths)))),'resolution')]
    names(x@metadata$feature_datasets)=c('name',paste0('file_',c(id_regions_subset,new_ids)),'resolution')
    x@features[[id_features_subset]]=c(x@features[[id_features_subset]],feat_new)
    x@length_features[[id_features_subset]]=c(x@length_features[[id_features_subset]],len_new)
  }else{
    new_ids=setdiff(id_regions_subset,'Control')
  }
  ### end NEW ###
  
  features_plot=features(x)
  if(log_scale){
    features_plot=lapply(features_plot,lapply,'+',log_shift)
    features_plot=lapply(features_plot,lapply,log)
    if(sum(unlist(lapply(features_plot,lapply,is.infinite))))
      stop('logarithm of 0.')
  }
  
  if(type %in% c('pairs','pairsSmooth')){
    if(length(unique(resolution(x)))!=1)
      stop('type\'',type,'\' but selected features with different resolution. Smooth data first to have the same resolution.')
    N_regions=rep(N_regions,length.out=length(id_regions_subset))
    N_regions=pmin(lengthRegions(x),N_regions)
    index_plot=mapply(sample,lengthRegions(x),N_regions,SIMPLIFY=FALSE)
    shuffle=sample(sum(N_regions))
    features_plot=lapply(features_plot,function(feature_plot) Reduce(cbind,mapply(function(feature,index_plot) feature[,index_plot],feature_plot,index_plot,SIMPLIFY=FALSE))[,shuffle])
    col=rep(col,length.out=length(N_regions))
    col_plot=rep(rep(col,N_regions)[shuffle],each=nrow(features_plot[[1]]))
    features_plot=do.call(cbind,lapply(features_plot,as.vector))
    features_cor=cor(features_plot[!is.na(rowSums(features_plot)),],method=method)
    z=list(features_plot=features_plot,features_cor=features_cor,type=type)
    if(plot){
      panel.cor <- function(x,y,digits=2,method_cor=method,prefix="",cex.cor, ...){
        usr <- par("usr"); on.exit(par(usr))
        par(usr = c(0, 1, 0, 1))
        index=(!is.na(x))&(!is.na(y))
        r <- cor(x[index], y[index],method=method_cor)
        txt <- format(c(r, 0.123456789), digits = digits)[1]
        txt <- paste0(prefix, txt)
        if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
        text(0.5, 0.5, txt, cex = cex.cor*(0.4+0.8*abs(r)))
      }
      devAskNewPage(ask)
      if(type=='pairs'){
        upper.panel=points
      }else{
        upper.panel=function(...) {par(new=TRUE);smoothScatter(...,nrpoints=0,colramp=colorRampPalette(c("white",col_plot[1])),add=TRUE)}
      }
      pairs(features_plot,main='Features correlation',labels=nameFeatures(x),col=col_plot,pch=3,lower.panel=panel.cor,upper.panel=upper.panel,...)
      invisible(z)
    }else{
      z
    }
  }else{
    if(alignment(x)=='left')
      x_plot=lapply(features_plot,function(feature_plot) 0.5:(nrow(feature_plot[[1]])-0.5))
    if(alignment(x)=='right')
      x_plot=lapply(features_plot,function(feature_plot) -(nrow(feature_plot[[1]])-0.5):(-0.5))
    if(alignment(x)=='center')
      x_plot=lapply(features_plot,function(feature_plot) seq_len(nrow(feature_plot[[1]]))-nrow(feature_plot[[1]])/2-0.5)
    if(alignment(x)=='scale'){
      length_features_plot=mapply(function(feature_plot,length_feature) mapply(function(feature,length) pmax(apply(feature,2,function(feature) rev(which(!is.na(feature)))[1]),length),
                                                                               feature_plot,length_feature,SIMPLIFY=FALSE),
                                  features_plot,x@length_features,SIMPLIFY=FALSE)
      if(sum(unlist(lapply(length_features_plot,function(length_feature) length(unique(unlist(length_feature)))!=1)))){
        if(type=='boxplot')
          stop('type \'boxplot\' is incompatible with \'scale\' alignment and regions of different length. Smooth data first.')
        if(average){
          warning('average=TRUE is incompatible with \'scale\' alignment and regions of different length. Setting average=FALSE.')
          average=FALSE
        }
        if(size){
          warning('size=TRUE is incompatible with \'scale\' alignment and regions of different length. Setting size=FALSE.')
          size=FALSE
        }
        x_plot=lapply(length_features_plot,
                      function(length_feature){
                        length_max=max(unlist(length_feature))
                        x_plot=lapply(length_feature,function(length) do.call(cbind,lapply(length,function(length) c(seq(0,1,length.out=length),rep(NA,length_max-length)))))
                        return(x_plot)
                      })
      }else{
        x_plot=lapply(features_plot,function(feature_plot) seq(0,1,length.out=nrow(feature_plot[[1]])))
      }
    }
    if(length(col)!=length(id_regions_subset)){
      warning('number of colors in \'col\' different from the number of region datasets considered.')
      col=rep(col,length.out=length(id_regions_subset))
    }
    if(average)
      features_average=lapply(features_plot,function(feature_plot) Reduce(cbind,lapply(feature_plot,rowMeans,na.rm=TRUE)))
    if(size){
      if(size_perc){
        ### NEW ###
        features_position_size=lapply(features_plot,function(feature_plot) do.call(cbind,lapply(rev(feature_plot[!(idRegions(x) %in% c('Control','Control_new'))]),function(feature) rowSums(!is.na(feature))/ncol(feature))))
      }else{
        features_position_size=lapply(features_plot,function(feature_plot) do.call(cbind,lapply(rev(feature_plot[!(idRegions(x) %in% c('Control','Control_new'))]),function(feature) rowSums(!is.na(feature)))))
      }
    }
    ### start NEW ###
    if(position){
      features_position=do.call(cbind,
                                lapply(rev(regions(x)[!(idRegions(x) %in% c('Control','Control_new'))]),
                                       function(region) 
                                         rowSums(apply(as.matrix(mcols(region)),1,
                                                       function(info){
                                                         info=pmax(pmin(info[1:2],100),0)
                                                         info[2]=min((100-info[1]),info[2])
                                                         return(c(rep(0,info[2]),rep(1,info[1]),rep(0,100-sum(info)))/length(region))
                                                       }))))
    }else{
      features_position=NULL
    }
    ### end NEW ###
    if(type=='curves'){
      N_regions=rep(N_regions,length.out=length(id_regions_subset))
      N_regions=pmin(lengthRegions(x)[id_regions_subset],N_regions)
      index_plot=mapply(sample,lengthRegions(x)[id_regions_subset],N_regions,SIMPLIFY=FALSE)
      shuffle=sample(sum(N_regions))
      features_plot=lapply(features_plot,function(feature_plot) Reduce(cbind,mapply(function(feature,index_plot) as.matrix(feature[,index_plot]),feature_plot,index_plot,SIMPLIFY=FALSE))[,shuffle])
      if(is.list(x_plot[[1]]))
        x_plot=lapply(x_plot,function(x_plot) Reduce(cbind,mapply(function(x,index_plot) x[,index_plot],x_plot,index_plot,SIMPLIFY=FALSE))[,shuffle])
      col_plot=rep(col,N_regions)[shuffle]
    }
    if(type=='boxplot'){
      features_plot=lapply(features_plot,function(feature_plot) lapply(feature_plot,function(feature) t(apply(feature,1,quantile,na.rm=TRUE,probs=probs))))
      col_plot=lapply(col,function(col) c(rgb(colorRamp(c('white',col))((1:4)/4)[-1,],alpha=50,max=255),col)) ### CHANGED ###
      names(col_plot)=id_regions_subset
      ### start NEW ###
      if(!is.null(lengths)){
        col_plot=c(col_plot,rep(col_plot[2],length(lengths)))
        names(col_plot)=c(id_regions_subset,new_ids)
      }
      ### end NEW ###
    }
    if(plot){
      devAskNewPage(ask)
      if(size){
        layout(matrix(1:3,nrow=3),heights=c(5,1,1))
        mar.left=6
      }else{
        mar.left=4
      }
      par(oma=c(0,0,0,8))
      if(is.null(ylim)){
        if(average){
          ylim=mapply(function(feature,average) pmin(range(c(unlist(feature),unlist(average)),na.rm=TRUE),c(0,+Inf)),features_plot,features_average,SIMPLIFY=FALSE)
        }else{
          ylim=lapply(features_plot,function(feature) pmin(range(unlist(feature),na.rm=TRUE),c(0,+Inf)))
        }
      }else{
        ylim=lapply(features_plot,function(feature) ylim)
      }
      ### start NEW ###
      id_regions_subset_old=id_regions_subset
      ### end NEW ###
      for(id_feature in id_features_subset){
        for(id in new_ids){ ### NEW ###
          ### start NEW ###
          if(!is.null(lengths)){
            id_regions_subset=c('Control',id)
          }
          ### end NEW ###
          par(mar=c(5,mar.left,4,6))
          if(type=='curves'){
            matplot(x_plot[[id_feature]],features_plot[[id_feature]],type='l',col=col_plot,ylim=ylim[[id_feature]],
                    main=paste(nameRegions(x)[id_regions_subset[2]],'vs',nameRegions(x)[id_regions_subset[1]]),xlab=xlab,
                    ylab=paste0(ifelse(log_scale,'log ',''),nameFeatures(x)[id_feature]),...)
          }
          if(type=='boxplot'){
            plot(1,type="n",xlim=range(x_plot[[id_feature]]),ylim=ylim[[id_feature]],
                 main=paste(nameRegions(x)[id_regions_subset[2]],'vs',nameRegions(x)[id_regions_subset[1]]),xlab=xlab,
                 ylab=paste0(ifelse(log_scale,'log ',''),nameFeatures(x)[id_feature]),...)
            for(id_region in id_regions_subset){
              ### CHANGED ###
              polygon(c(x_plot[[id_feature]],rev(x_plot[[id_feature]])),c(features_plot[[id_feature]][[id_region]][,1],rev(features_plot[[id_feature]][[id_region]][,length(probs)])),
                      col=col_plot[[id_region]][1],border=col_plot[[id_region]][3])
              polygon(c(x_plot[[id_feature]],rev(x_plot[[id_feature]])),c(features_plot[[id_feature]][[id_region]][,2],rev(features_plot[[id_feature]][[id_region]][,length(probs)-1])),
                      col=col_plot[[id_region]][2],border=col_plot[[id_region]][3])
              lines(x_plot[[id_feature]],features_plot[[id_feature]][[id_region]][,3],col=col_plot[[id_region]][4],lty=1,lwd=2)
            }
          }
          if(average)
            matplot(x_plot[[id_feature]],features_average[[id_feature]],type='l',col=col,lty=1,lwd=2,add=TRUE)
          ### start NEW ###
          if(sum(!is.na(as.data.frame(mcols(x@regions[setdiff(id_regions_subset,c('Control','Control_new'))])))))
            text(x_plot[[id_feature]],rep(0,100),as.data.frame(mcols(x@regions[setdiff(id_regions_subset,c('Control','Control_new'))])),cex=0.8,adj=c(0.5,-0.2))
          if(zero_line)
            abline(h=0)
          if(!is.null(reproducible_pval)){
            x_rect=x_plot[[id_feature]][reproducible_pval]
            points(x_rect,rep(par('usr')[3]*0.5,length(x_rect)),col='black',pch=22,bg='black')
          }
          ### end NEW ###
          args=as.list(match.call())
          if(is.null(args$cex)){
            cex=ifelse(is.null(args$cex.lab),1,args$cex.lab)
          }else{
            cex=args$cex
          }
          legend(par('usr')[2],mean(par('usr')[3:4]),legend=nameRegions(x)[rev(id_regions_subset)],xpd=NA,bty='n',lty=1,lwd=2,col=rev(col),yjust=0.5,seg.len=1,cex=cex)
          if(size){
            par(mar=c(1,mar.left,2,6))
            if(size_perc){ 
              ### NEW ###
              #col_perc=heat.colors(1000000)
              #col_perc=rev(c(col_perc[c(seq(1,700001,70),round(seq(700001,1000000,length.out=90000))[-1])]))
              #col_perc=rainbow(1000000,start=0,end=0.25)
              #col_perc=rev(c(col_perc[c(seq(1,600001,60),round(seq(600001,1000000,length.out=90000))[-1])]))
              col_perc=heat.colors(1000000)
              col_perc=rev(c(col_perc[seq(1,700001,70)],colorRampPalette(c(col_perc[700001],"#8080FFFF"))(60001)[-1],colorRampPalette(c("#8080FFFF","#FFFFFFFF"))(35001)[1:30000][-1]))
              col_perc=c("white",rep(col_perc,each=10))
              
              if(position){
                if(!is.null(lengths)){
                  image(x_plot[[id_feature]],seq_len(2*(length(id_regions_subset)-1)),
                        cbind(features_position[,id],features_position_size[[id_feature]][,id]),
                        col=col_perc,xlim=par('usr')[1:2],ylim=range(seq_len(2*(length(id_regions_subset)-1)))+c(-0.5,0.5),
                        zlim=c(0,1),axes=FALSE,xlab='',ylab='',...)
                }else{
                  image(x_plot[[id_feature]],seq_len(2*(length(id_regions_subset)-1)),
                        cbind(features_position[,!(colnames(features_position) %in% c('Control','Control_new'))],
                              features_position_size[[id_feature]][,!(colnames(features_position_size[[id_feature]]) %in% c('Control','Control_new'))]),
                        col=col_perc,xlim=par('usr')[1:2],ylim=range(seq_len(2*(length(id_regions_subset)-1)))+c(-0.5,0.5),
                        zlim=c(0,1),axes=FALSE,xlab='',ylab='',...)
                }
              }else{
                image(x_plot[[id_feature]],seq_along(setdiff(id_regions_subset,c('Control','Control_new'))),
                      as.matrix(features_position_size[[id_feature]][,!(colnames(features_position_size[[id_feature]]) %in% c('Control','Control_new'))]),
                      col=col_perc,xlim=par('usr')[1:2],ylim=range(seq_along(setdiff(id_regions_subset,c('Control','Control_new'))))+c(-0.5,0.5),
                      zlim=c(0,1),axes=FALSE,xlab='',ylab='',...)
              }
            }else{
              image(x_plot[[id_feature]],seq_along(setdiff(id_regions_subset,c('Control','Control_new'))),
                    as.matrix(features_position_size[[id_feature]][,!(colnames(features_position_size[[id_feature]]) %in% c('Control','Control_new'))]),
                    col=cm.colors(101),xlim=par('usr')[1:2],ylim=range(seq_along(setdiff(id_regions_subset,c('Control','Control_new'))))+c(-0.5,0.5),
                    axes=FALSE,xlab='',ylab='',...)
            }
            axis(side=3,...)
            if(position){
              ### NEW ###
              if(!is.null(lengths)){
                index_pos=which(features_position[,id]>0)
              }else{
                index_pos=which(features_position>0)
              }
              index_pos=c(min(index_pos),max(index_pos))
              #lines(x_plot[[id_feature]][rep(index_pos[1],each=2)]-0.5,c(0.5,1.5),col='black',lty=2)
              #lines(x_plot[[id_feature]][rep(index_pos[2],each=2)]+0.5,c(0.5,1.5),col='black',lty=2)
              axis(side=2,at=seq_len(2*(length(id_regions_subset)-1)),labels=c('Feature','IPD values'),tick=FALSE,las=1,line=-2,...)
              axis(side=4,at=mean(seq_len(2*(length(id_regions_subset)-1))),labels=paste0('Total:\n',lengthRegions(x)[setdiff(id_regions_subset,c('Control','Control_new'))],' windows'),tick=FALSE,las=1,line=-2,...)
              x.rect=range(x_plot[[id_feature]])+c(-1,1)*diff(x_plot[[id_feature]])[1]/2
              rect(x.rect[1],seq_len(2*(length(id_regions_subset)-1))-0.5,x.rect[2],seq_len(2*(length(id_regions_subset)-1))+0.5,border='black',...)
            }else{
              axis(side=2,at=seq_along(setdiff(id_regions_subset,c('Control','Control_new'))),labels=c('Sample size'),tick=FALSE,las=1,line=-2,...)
              axis(side=4,at=seq_along(setdiff(id_regions_subset,c('Control','Control_new'))),labels=paste0('Total:\n',lengthRegions(x)[setdiff(id_regions_subset,c('Control','Control_new'))],' windows'),tick=FALSE,las=1,line=-2,...)
              x.rect=range(x_plot[[id_feature]])+c(-1,1)*diff(x_plot[[id_feature]])[1]/2
              rect(x.rect[1],seq_along(setdiff(id_regions_subset,c('Control','Control_new')))-0.5,x.rect[2],seq_along(setdiff(id_regions_subset,c('Control','Control_new')))+0.5,border='black',...)
            }
            par(mar=c(5,mar.left,0.5,6))
            if(size_perc){
              ### NEW ###
              image(seq(-50,50,length.out=1001),1,as.matrix(seq.int(1001)),xlim=c(-100,100),axes=FALSE,xlab='Percentage',ylab='',col=col_perc,...)
              labels=seq(0,1,length.out=5)
            }else{
              image(-50:50,1,as.matrix(seq.int(101)),xlim=c(-100,100),axes=FALSE,xlab='Sample size',ylab='',col=cm.colors(101),...)
              labels=seq(min(features_position_size[[id_feature]][,!(colnames(features_position_size[[id_feature]]) %in% c('Control','Control_new'))]),max(features_position_size[[id_feature]][,!(colnames(features_position_size[[id_feature]]) %in% c('Control','Control_new'))]),length.out=5)
            }
            if(length(unique(labels))==1)
              labels=-2:2+labels
            if(size_perc){
              ### NEW ###
              axis(side=1,at=c(-50,-25,0,25,50),labels=paste0(labels*100,'%'),...)
              rect(-50.1,par('usr')[3],50.1,par('usr')[4],border='black',...)
            }else{
              axis(side=1,at=c(-50,-25,0,25,50),labels=labels,...)
              rect(-50.5,par('usr')[3],50.5,par('usr')[4],border='black',...)
            }
          }
        }
      }
      if(!average)
        features_average=NULL
      if(!size)
        features_position_size=NULL
      z=list(x_plot=x_plot,features_plot=features_plot,features_average=features_average,features_position_size=features_position_size,features_position=features_position,type=type,col=col,col_plot=col_plot)
      invisible(z)
    }else{
      if(!average)
        features_average=NULL
      if(!size)
        features_position_size=NULL
      z=list(x_plot=x_plot,features_plot=features_plot,features_average=features_average,features_position_size=features_position_size,features_position=features_position,type=type,col=col,col_plot=col_plot)
      z
    }
  }
}


plotTest_same_sample_size <- function(regionsFeatures,alpha=0.05,scale_threshold=NULL,nlevel=100,type='boxplot',
                                      N_regions=pmin(lengthRegions(regionsFeatures),10),
                                      probs=c(0.25,0.5,0.75),average=TRUE,size=TRUE,
                                      size_perc=FALSE, # NEW plot the sample size in each position as percentage of the total sample size
                                      position=FALSE, # NEW plot the feature positions
                                      id_features_subset=idFeatures(regionsFeatures),
                                      col=1+seq_len(nRegions(regionsFeatures)),
                                      ask=TRUE,xlab='Windows',ylim=NULL,...){
  # MODIFIED FROM BIOCONDUCTOR PACKAGE IWTomics
  # FUNCTION plotTest
  
  if(class(regionsFeatures)!='IWTomicsData')
    stop('invalid regionsFeatures. IWTomicsData object expected.')
  if(!validObject(regionsFeatures))
    stop('invalid regionsFeatures.')
  if(is.null(regionsFeatures@test))
    stop('No test results present in regionFeatures.')
  if(sum(!(id_features_subset %in% idFeatures(regionsFeatures))))
    stop('invalid id_features_subset. The features provided are not listed in regionFeatures.')
  if(!(type %in% c('curves','boxplot')))
    stop('invalid plot type \'',type,'\'. Available types are \'curves\' and \'boxplot\'.')
  ylim_ext=ylim
  col=rep(col,length.out=nRegions(regionsFeatures))
  names(col)=idRegions(regionsFeatures)
  N_regions=rep(N_regions,length.out=nRegions(regionsFeatures))
  names(N_regions)=idRegions(regionsFeatures)
  
  devAskNewPage(ask)
  for(i in seq_len(nTests(regionsFeatures))){
    if(is.null(scale_threshold)){
      scale_threshold_i=lapply(.testResults(regionsFeatures)[[i]][id_features_subset],function(feature) feature$max_scale)
    }else{
      if(is.list(scale_threshold)){
        scale_threshold_i=as.list(rep(scale_threshold[[i]],length.out=length(id_features_subset)))
      }else{
        scale_threshold_i=as.list(rep(scale_threshold,length.out=length(id_features_subset)))
      }
      names(scale_threshold_i)=id_features_subset
    }
    id_region1_i=idRegionsTest(regionsFeatures,i)[[1]][1]
    id_region2_i=idRegionsTest(regionsFeatures,i)[[1]][2]
    if((!is.null(id_region2_i))&&(id_region2_i==''))
      id_region2_i=NULL
    for(id_feature in id_features_subset){
      if(.testResults(regionsFeatures)[[i]][[id_feature]]$test=='1pop'){
        plot_data=plot_no_size_control(regionsFeatures,type=type,
                                       N_regions=N_regions[id_region1_i],probs=probs,average=average,size=size,
                                       size_perc=size_perc,position=position, ### MODIFIED ###
                                       id_regions_subset=id_region1_i,id_features_subset=id_feature,
                                       col=col[id_region1_i],plot=FALSE)
      }else{
        ### start NEW ###
        regionsFeatures_sub=regionsFeatures
        if(!is.null(regionsFeatures@test$result[[i]][[id_feature]]$index.data2)){
          if('Control' %in% idRegions(regionsFeatures)){
            regionsFeatures_sub@metadata$region_datasets['Control','size']=length(regionsFeatures@test$result[[i]][[id_feature]]$index.data2)
            regionsFeatures_sub@features[[id_feature]]$Control=regionsFeatures@features[[id_feature]]$Control[,regionsFeatures@test$result[[i]][[id_feature]]$index.data2]
            regionsFeatures_sub@length_features[[id_feature]]$Control=regionsFeatures@length_features[[id_feature]]$Control[regionsFeatures@test$result[[i]][[id_feature]]$index.data2]
          }else{
            regionsFeatures_sub@metadata$region_datasets['Control_new','size']=length(regionsFeatures@test$result[[i]][[id_feature]]$index.data2)
            regionsFeatures_sub@features[[id_feature]]$Control_new=regionsFeatures@features[[id_feature]]$Control_new[,regionsFeatures@test$result[[i]][[id_feature]]$index.data2]
            regionsFeatures_sub@length_features[[id_feature]]$Control_new=regionsFeatures@length_features[[id_feature]]$Control_new[regionsFeatures@test$result[[i]][[id_feature]]$index.data2]
            
          }
        }
        ### end NEW ###
        plot_data=plot_no_size_control(regionsFeatures_sub,type=type,
                                       N_regions=N_regions[c(id_region1_i,id_region2_i)],probs=probs,average=average,
                                       size=size,size_perc=size_perc,position=position, ### MODIFIED ###
                                       id_regions_subset=c(id_region1_i,id_region2_i),id_features_subset=id_feature,
                                       col=col[c(id_region1_i,id_region2_i)],plot=FALSE)
      }
      if(size){
        layout(rbind(1:2,cbind(3:6,0)),widths=c(8,2),heights=c(4,2,2,1,1))
        mar.left=8
      }else{
        layout(rbind(1:2,cbind(3:4,0)),widths=c(8,2),heights=c(4,2,2))
        mar.left=8
      }
      p=ncol(.testResults(regionsFeatures)[[i]][[id_feature]]$adjusted_pval_matrix)
      max_scale=.testResults(regionsFeatures)[[i]][[id_feature]]$max_scale
      if((scale_threshold_i[[id_feature]]<1)||(scale_threshold_i[[id_feature]]>p)){
        warning('invalid scale_threshold. Setting it to the default value.',call.=FALSE,immediate.=TRUE)
        scale_threshold_i[[id_feature]]=max_scale
      }
      
      # Adjusted p-value heatmap
      par(mar=c(4,mar.left,6.5,3.5),mgp=c(2.5,1,0),xpd=FALSE,las=1)
      col_heatmap=rainbow(nlevel,start=0.15,end=0.67)[nlevel:1]
      x_plot=plot_data$x_plot[[id_feature]]
      image(x_plot,1:max_scale-0.5,t(.testResults(regionsFeatures)[[i]][[id_feature]]$adjusted_pval_matrix[p:(p-max_scale+1),]), 
            col=col_heatmap,ylab="Maximum interval length",main=paste(paste(nameRegions(regionsFeatures)[id_region1_i],"vs",nameRegions(regionsFeatures)[id_region2_i]),"\nAdjusted p-value heatmap",sep='\n'), 
            xlab=xlab,zlim=c(0,1),mgp=c(3,0.8,0),...)
      abline(h=scale_threshold_i[[id_feature]])
      abline(a=-2*max_scale/(max(x_plot)-min(x_plot)+diff(x_plot[1:2]))*(x_plot[1]-diff(x_plot[1:2])/2),
             b=2*max_scale/(max(x_plot)-min(x_plot)+diff(x_plot[1:2])))
      abline(a=2*max_scale/(max(x_plot)-min(x_plot)+diff(x_plot[1:2]))*(rev(x_plot)[1]+diff(x_plot[1:2])/2),
             b=-2*max_scale/(max(x_plot)-min(x_plot)+diff(x_plot[1:2])))
      box()
      par(mar=c(4,0.5,6.5,4))
      image(1,seq.int(nlevel)-(nlevel+1)/2,t(as.matrix(seq.int(nlevel))),axes=FALSE,xlab='',ylab='',col=col_heatmap)
      axis(side=4,at=seq(0,nlevel,length.out=6)-nlevel/2,labels=format(seq(0,1,length.out=6),2),padj=0.4,...)
      box()
      
      # Adjusted p-value plot at scale_threshold
      par(mar=c(4,mar.left,2.6,3.5))
      pval_scale_threshold=adjusted_pval(regionsFeatures,i,id_feature,scale_threshold_i)[[1]][[1]]
      plot(1,type="n",xlim=c(x_plot[1]-diff(x_plot[1:2])/2,rev(x_plot)[1]+diff(x_plot[1:2])/2),ylim=c(0,1),ylab="p-value",xlab=xlab,xaxs="i",
           main=paste0("Adjusted p-values - Threshold ",scale_threshold_i[[id_feature]]),mgp=c(3,0.8,0),...)
      low.p.value=which((pval_scale_threshold<0.05)&(!is.infinite(pval_scale_threshold)))
      if(length(low.p.value)){
        for(j in seq_along(low.p.value))
          rect(x_plot[low.p.value[j]]-diff(x_plot[1:2])/2,par("usr")[3],
               x_plot[low.p.value[j]]+diff(x_plot[1:2])/2,par("usr")[4],col="gray90",density=-2,border=NA)
      }
      box()
      for (j in 0:10)
        abline(h=j/10,col="lightgray",lty="dotted")
      points(x_plot,pval_scale_threshold,pch=16)
      
      # Data plot
      if(is.null(ylim_ext)){
        if(average){
          ylim=pmin(range(c(unlist(plot_data$features_plot[[1]]),unlist(plot_data$features_average[[1]])),na.rm=TRUE),c(0,+Inf))
        }else{
          ylim=pmin(range(unlist(plot_data$features_plot[[1]]),na.rm=TRUE),c(0,+Inf))
        }
      }else{
        ylim=ylim_ext
      }
      par(mar=c(4,mar.left,1,3.5),las=0)
      if(type=='curves'){
        matplot(x_plot,plot_data$features_plot[[id_feature]],type='l',col=plot_data$col_plot,xlim=par('usr')[1:2],ylim=ylim,xaxs="i",
                main=nameFeatures(regionsFeatures)[id_feature],xlab=xlab,
                ylab=nameFeatures(regionsFeatures)[id_feature],...)
      }
      if(type=='boxplot'){
        plot(1,type="n",xlim=par('usr')[1:2],ylim=ylim,xaxs="i",xlab=xlab,
             ylab=nameFeatures(regionsFeatures)[id_feature],mgp=c(3,0.8,0),...) ### MODIFIED ###
        for(id_region in c(id_region2_i,id_region1_i)){
          ### MODIFIED ###
          polygon(c(x_plot,rev(x_plot)),c(plot_data$features_plot[[id_feature]][[id_region]][,1],rev(plot_data$features_plot[[id_feature]][[id_region]][,length(probs)])),
                  col=plot_data$col_plot[[id_region]][1],border=plot_data$col_plot[[id_region]][3])
          polygon(c(x_plot,rev(x_plot)),c(plot_data$features_plot[[id_feature]][[id_region]][,2],rev(plot_data$features_plot[[id_feature]][[id_region]][,length(probs)-1])),
                  col=plot_data$col_plot[[id_region]][2],border=plot_data$col_plot[[id_region]][3])
          lines(x_plot,plot_data$features_plot[[id_feature]][[id_region]][,3],col=plot_data$col_plot[[id_region]][4],lty=1,lwd=2)
        }
      }
      if(average)
        matplot(x_plot,plot_data$features_average[[id_feature]],type='l',col=col[c(id_region1_i,id_region2_i)],lty=1,lwd=2,add=TRUE)
      ### start NEW ###
      if(sum(!is.na(as.data.frame(mcols(regionsFeatures@regions[id_region1_i])))))
        text(x_plot,rep(0,100),as.data.frame(mcols(regionsFeatures@regions[id_region1_i])),cex=0.8,adj=c(0.5,-0.2))
      if(.testResults(regionsFeatures)[[i]][[id_feature]]$test=='1pop')
        lines(x_plot,rep(.testResults(regionsFeatures)[[i]][[id_feature]]$mu,p),col='black',type='l',lwd=2)
      ### end NEW ###
      args=as.list(match.call())
      if(is.null(args$cex)){
        cex=ifelse(is.null(args$cex.lab),1,args$cex.lab)
      }else{
        cex=args$cex
      }
      legend(par('usr')[2],mean(par('usr')[3:4]),legend=nameRegions(regionsFeatures)[c(id_region1_i,id_region2_i)],xpd=NA,bty='n',lty=1,lwd=2,col=col[c(id_region1_i,id_region2_i)],yjust=0.5,cex=cex)
      if(size){
        par(mar=c(1,mar.left,2.5,3.5))
        if(size_perc){
          ### NEW ###
          #col_perc=heat.colors(1000000)
          #col_perc=rev(c(col_perc[c(seq(1,700001,70),round(seq(700001,1000000,length.out=90000))[-1])]))
          #col_perc=rainbow(1000000,start=0,end=0.25)
          #col_perc=rev(c(col_perc[c(seq(1,600001,60),round(seq(600001,1000000,length.out=90000))[-1])]))
          col_perc=heat.colors(1000000)
          col_perc=rev(c(col_perc[seq(1,700001,70)],colorRampPalette(c(col_perc[700001],"#8080FFFF"))(60001)[-1],colorRampPalette(c("#8080FFFF","#FFFFFFFF"))(35001)[1:30000][-1]))
          col_perc=c("white",rep(col_perc,each=10))
          if(position){
            image(x_plot,seq_len(2*length(id_region1_i)),
                  cbind(plot_data$features_position,
                        plot_data$features_position_size[[id_feature]]),
                  col=col_perc,xlim=par('usr')[1:2],ylim=range(seq_len(2*length(id_region1_i)))+c(-0.5,0.5),
                  zlim=c(0,1),axes=FALSE,xlab='',ylab='',mgp=c(3,0.8,0),...)
          }else{
            image(x_plot,seq_len(ncol(plot_data$features_position_size[[id_feature]])),plot_data$features_position_size[[id_feature]],
                  col=col_perc,xlim=par('usr')[1:2],ylim=range(seq_len(ncol(plot_data$features_position_size[[id_feature]])))+c(-0.5,0.5),
                  zlim=c(0,1),axes=FALSE,xlab='',ylab='',mgp=c(3,0.8,0),...)
          }
        }else{
          image(x_plot,seq_along(id_region1_i),plot_data$features_position_size[[id_feature]],
                col=cm.colors(101),xlim=par('usr')[1:2],ylim=range(seq_along(id_region1_i))+c(-0.5,0.5),
                axes=FALSE,xlab='',ylab='',...)
        }
        axis(side=3,...)
        if(position){
          ### NEW ###
          axis(side=2,at=seq_len(2*length(id_region1_i)),labels=c('Feature','Sample size'),mgp=c(3,1.5,0),tick=FALSE,las=1,line=-1,...)
          axis(side=4,at=mean(seq_len(2*length(id_region1_i))),labels=paste0('Total:\n',lengthRegions(regionsFeatures)[id_region1_i],' windows'),mgp=c(3,1.5,0),tick=FALSE,las=1,line=-1,...)
          x.rect=range(x_plot)+c(-1,1)*diff(x_plot)[1]/2
          rect(x.rect[1],seq_len(2*length(id_region1_i))-0.5,x.rect[2],seq_len(2*length(id_region1_i))+0.5,border='black')
        }else{
          axis(side=2,at=seq_len(ncol(plot_data$features_position_size[[id_feature]])),labels=rep('Sample size',ncol(plot_data$features_position_size[[id_feature]])),mgp=c(3,1.5,0),tick=FALSE,las=1,line=-1,...)
          axis(side=4,at=seq_along(id_region1_i),labels=paste0('Total:\n',lengthRegions(regionsFeatures)[id_region1_i],' windows'),mgp=c(3,1.5,0),tick=FALSE,las=1,line=-1,...)
          x.rect=range(x_plot)+c(-1,1)*diff(x_plot)[1]/2
          rect(x.rect[1],seq_len(ncol(plot_data$features_position_size[[id_feature]]))-0.5,x.rect[2],seq_len(ncol(plot_data$features_position_size[[id_feature]]))+0.5,border='black')
        }
        par(mar=c(5,mar.left,0.5,2))
        if(size_perc){
          ### NEW ###
          image(seq(-50,50,length.out=1001),1,as.matrix(seq.int(1001)),xlim=c(-100,100),axes=FALSE,xlab='Percentage',ylab='',mgp=c(3,1.5,0),col=col_perc,...)
          labels=seq(0,1,length.out=5)
        }else{
          image(-50:50,1,as.matrix(seq.int(101)),xlim=c(-100,100),axes=FALSE,xlab='Sample size',ylab='',col=cm.colors(101),...)
          labels=seq(min(plot_data$features_position_size[[id_feature]]),max(plot_data$features_position_size[[id_feature]]),length.out=5)
        }
        ### start NEW ###
        if(length(unique(labels))==1)
          labels=-2:2+labels
        ### end NEW ###
        if(size_perc){
          ### NEW ###
          axis(side=1,at=c(-50,-25,0,25,50),labels=paste0(labels*100,'%'),...)
          rect(-50.1,par('usr')[3],50.1,par('usr')[4],border='black')
        }else{
          axis(side=1,at=c(-50,-25,0,25,50),labels=labels,...)
          rect(-50.5,par('usr')[3],50.5,par('usr')[4],border='black')
        }
      }
    }
  }
}

IWTomicsTest_same_sample_size <- function(regionsFeatures,
                                          id_region1=idRegions(regionsFeatures)[1],id_region2=NULL, 
                                          id_features_subset=idFeatures(regionsFeatures),
                                          mu=0,statistics="mean",probs=0.5,max_scale=NULL,paired=FALSE,B=1000){
  # MODIFIED FROM BIOCONDUCTOR PACKAGE IWTomics
  # FUNCTION IWTomicsTest
  
  if(class(regionsFeatures)!='IWTomicsData')
    stop('invalid regionsFeatures. IWTomicsData object expected.')
  if(!validObject(regionsFeatures))
    stop('invalid regionsFeatures.')
  if((length(id_region2)==1)&&(id_region2==''))
    id_region2=NULL
  regionsFeatures_tot=initialize(regionsFeatures,
                                 test=list(input=list(id_region1=id_region1,id_region2=id_region2,
                                                      id_features_subset=id_features_subset,mu=mu,
                                                      statistics=statistics,probs=probs,max_scale=max_scale,
                                                      paired=paired,B=B),
                                           result=list()))
  regionsFeatures=regionsFeatures[setdiff(union(id_region1,id_region2),''),id_features_subset]
  regionsFeatures=initialize(regionsFeatures,
                             test=list(input=list(id_region1=id_region1,id_region2=id_region2,
                                                  id_features_subset=id_features_subset,mu=mu,
                                                  statistics=statistics,probs=probs,max_scale=max_scale,
                                                  paired=paired,B=B),
                                       result=list()))
  if(!is.null(id_region2)){
    if(paired){
      if(sum(lengthRegions(regionsFeatures)[id_region1]!=lengthRegions(regionsFeatures)[id_region2],na.rm=TRUE))
        stop('cannot perform a paired test with different sample size in the two region datasets. ')
    }
  }
  
  regionsFeatures@test$result=vector("list",length(id_region1))
  for(i in seq_along(id_region1)){
    id_region1_i=id_region1[i]
    id_region2_i=id_region2[i]
    if((!is.null(id_region2_i))&&(id_region2_i==''))
      id_region2_i=NULL
    if(is.null(id_region2_i)){
      write2='...'
    }else{
      write2=paste0(' vs. \'',nameRegions(regionsFeatures[id_region2_i,]),'\'...')
    }
    message('Performing IWT for \'',nameRegions(regionsFeatures)[id_region1_i],'\'',write2)
    features_test_i=lapply(features(regionsFeatures),function(feature) feature[c(id_region1_i,id_region2_i)])
    features_test_i_not_NA=lapply(features_test_i,function(feature_test) rowSums(!is.na(Reduce(cbind,feature_test)))!=0)
    features_test_i=mapply(function(feature_test,features_test_i_not_NA) lapply(feature_test,function(feature) as.matrix(feature[features_test_i_not_NA,])),
                           features_test_i,features_test_i_not_NA,SIMPLIFY=FALSE)
    if(alignment(regionsFeatures)=='scale'){
      if(sum(unlist(lapply(features_test_i,function(feature_test) lapply(feature_test,function(feature) sum(is.na(feature[nrow(feature),])))))))
        stop('IWTomicsTest is incompatible with \'scale\' alignment and regions of different length. Smooth data first.')
    }
    if(sum((unlist(lapply(mu,length))!=1)&(unlist(lapply(mu,length))!=unlist(lapply(features_test_i,function(feature) nrow(feature[[1]]))))))
      stop('invalid mu. Grids different from feature grids.')
    mu_i=as.list(rep(mu,length.out=length(id_features_subset)))
    names(mu_i)=id_features_subset
    if(is.null(max_scale)){
      max_scale_i=lapply(features_test_i,function(feature_test) nrow(feature_test[[id_region1_i]]))
    }else{
      if(is.list(max_scale)){
        max_scale_i=as.list(rep(max_scale[[i]],length.out=length(id_features_subset)))
      }else{
        max_scale_i=as.list(rep(max_scale,length.out=length(id_features_subset)))
      }
      names(max_scale_i)=id_features_subset
    }
    regionsFeatures@test$result[[i]]=list()
    for(id_feature in id_features_subset){
      message('  Performing IWT for feature \'',nameFeatures(regionsFeatures)[id_feature],'\'...')
      p=nrow(features_test_i[[id_feature]][[id_region1_i]])
      if((max_scale_i[[id_feature]]<1)||(max_scale_i[[id_feature]]>p))
        warning('invalid max_scale. Setting it to the default value.',call.=FALSE,immediate.=TRUE)
      max_scale_i[[id_feature]]=ifelse(max_scale_i[[id_feature]]<1,p,min(max_scale_i[[id_feature]],p))
      
      if(is.null(id_region2_i)){
        regionsFeatures@test$result[[i]][[id_feature]]=.IWTomics.1pop(features_test_i[[id_feature]][[id_region1_i]],
                                                                      mu=mu_i[[id_feature]],statistics=statistics,probs=probs,
                                                                      max_scale=max_scale_i[[id_feature]],B=B)
        if(regionsFeatures@test$result[[i]][[id_feature]]$exact)
          warning('number of iteration B greater than number of possible permutations. Exact p-values computed.',call.=FALSE,immediate.=TRUE)
      }else{
        ### MODIFIED ###
        tmp=.IWTomics.2pop_same_sample_size(features_test_i[[id_feature]][[id_region1_i]],features_test_i[[id_feature]][[id_region2_i]],
                                            mu=mu_i[[id_feature]],statistics=statistics,probs=probs,max_scale=max_scale_i[[id_feature]],paired=paired,B=B)
        if(TRUE %in% tmp$no.pval)
          warning('p-value not fully computable in some points, because of too many NAs present.',call.=FALSE,immediate.=TRUE)
        regionsFeatures@test$result[[i]][[id_feature]]=tmp$result
        if(regionsFeatures@test$result[[i]][[id_feature]]$exact)
          warning('number of iteration B greater than number of possible permutations. Exact p-values computed.',call.=FALSE,immediate.=TRUE)
      }
      regionsFeatures@test$result[[i]][[id_feature]]$notNA=features_test_i_not_NA[[id_feature]]
    }
  }
  
  return(regionsFeatures)
}

.IWTomics.2pop_same_sample_size <- function(data1,data2,mu=0,statistics='mean',probs=0.5,max_scale=nrow(data1),paired=FALSE,B=1000,recycle=FALSE){
  # MODIFIED FROM BIOCONDUCTOR PACKAGE IWTomics
  # FUNCTION .IWTomics.2pop
  
  # recycle     if TRUE, edges are recycled.
  # max_scale   max interval length for the adjustment.
  
  if(statistics=='median'){
    statistics='quantile'
    probs=0.5
  }
  n1=ncol(data1)
  n2=ncol(data2)
  ### start NEW ###
  if(n1<n2){
    index.data2=sample(n2,n1)
    data2=data2[,index.data2]
  }else{
    index.data2=seq_len(n2)
  }
  ### end NEW ###
  n2=ncol(data2)
  n=n1+n2
  data1=data1-mu
  if(paired){
    exact=(B>=(2^n1))
  }else{
    exact=(B>=choose(n,n1))
  }
  # Non computable p-values (all NAs in one group)
  allNA=(rowSums(!is.na(data1))==0)|(rowSums(!is.na(data2))==0)
  data1=data1[!allNA,]
  data2=data2[!allNA,]
  max_scale=min(max_scale,nrow(data1))
  
  p=nrow(data1)
  result=list(test='2pop',mu=mu,max_scale=max_scale)
  data=cbind(data1,data2)
  
  # Univariate permutations
  message('    Point-wise tests...')
  if(statistics=='mean'){
    T0_plot=rowMeans(data1,na.rm=TRUE)-rowMeans(data2,na.rm=TRUE)
    T0=(T0_plot)^2
  }
  if(statistics=='quantile'){
    T0_plot=apply(data1,1,quantile,probs=probs,na.rm=TRUE)-apply(data2,1,quantile,probs=probs,na.rm=TRUE)
    T0=(T0_plot)^2
    if(is.matrix(T0_plot)){
      T0_plot=colSums(T0_plot)
      T0=colSums(T0)
    }
  }
  if(statistics=='variance'){
    T0_plot=apply(data1,1,var,na.rm=TRUE)/apply(data2,1,var,na.rm=TRUE)
    T0=(T0_plot)^2
  }
  if(exact){
    if(paired){
      T_perm=do.call(cbind,lapply(seq_len(n1+1)-1,
                                  function(m){
                                    group_change=combn(n1,m)
                                    T_perm=apply(group_change,2,
                                                 function(change){
                                                   data_perm=data
                                                   data_perm[,c(change,n1+change)]=data_perm[,c(n1+change,change)]
                                                   if(statistics=='mean')
                                                     return((rowMeans(as.matrix(data_perm[,seq.int(n1)]),na.rm=TRUE)-rowMeans(as.matrix(data_perm[,n1+seq.int(n2)]),na.rm=TRUE))^2)
                                                   if(statistics=='quantile'){
                                                     T_perm=(apply(as.matrix(data_perm[,seq.int(n1)]),1,quantile,probs=probs,na.rm=TRUE)-apply(as.matrix(data_perm[,n1+seq.int(n2)]),1,quantile,probs=probs,na.rm=TRUE))^2
                                                     if(is.matrix(T_perm))
                                                       return(colSums(T_perm))
                                                     return(T_perm)
                                                   }
                                                   if(statistics=='variance')
                                                     return((apply(as.matrix(data_perm[,seq.int(n1)]),1,var,na.rm=TRUE)/apply(as.matrix(data_perm[,n1+seq.int(n2)]),1,var,na.rm=TRUE))^2)
                                                 })
                                    return(T_perm)
                                  }))
    }else{
      first_group=combn(n,n1)
      T_perm=apply(first_group,2,
                   function(group){
                     data_perm=data[,c(group,setdiff(seq_len(n1+n2),group))]
                     if(statistics=='mean')
                       return((rowMeans(as.matrix(data_perm[,seq.int(n1)]),na.rm=TRUE)-rowMeans(as.matrix(data_perm[,n1+seq.int(n2)]),na.rm=TRUE))^2)
                     if(statistics=='quantile'){
                       T_perm=(apply(as.matrix(data_perm[,seq.int(n1)]),1,quantile,probs=probs,na.rm=TRUE)-apply(as.matrix(data_perm[,n1+seq.int(n2)]),1,quantile,probs=probs,na.rm=TRUE))^2
                       if(is.matrix(T_perm))
                         return(colSums(T_perm))
                       return(T_perm)
                     }
                     if(statistics=='variance')
                       return((apply(as.matrix(data_perm[,seq.int(n1)]),1,var,na.rm=TRUE)/apply(as.matrix(data_perm[,n1+seq.int(n2)]),1,var,na.rm=TRUE))^2)
                   })
    }
  }else{
    T_perm=do.call(cbind,lapply(seq.int(B-1),
                                function(perm){
                                  if(paired){
                                    couple.perm=rbinom(n1,1,0.5)
                                    data_perm=data[,c(n1*couple.perm,-n1*couple.perm)+seq.int(2*n1)]
                                  }else{
                                    permutation=sample(n,n1)
                                    data_perm=data[,c(permutation,setdiff(seq_len(n),permutation))]
                                  }
                                  if(statistics=='mean')
                                    return((rowMeans(as.matrix(data_perm[,seq.int(n1)]),na.rm=TRUE)-rowMeans(as.matrix(data_perm[,n1+seq.int(n2)]),na.rm=TRUE))^2)
                                  if(statistics=='quantile'){
                                    T_perm=(apply(as.matrix(data_perm[,seq.int(n1)]),1,quantile,probs=probs,na.rm=TRUE)-apply(as.matrix(data_perm[,n1+seq.int(n2)]),1,quantile,probs=probs,na.rm=TRUE))^2
                                    if(is.matrix(T_perm))
                                      return(colSums(T_perm))
                                    return(T_perm)
                                  }
                                  if(statistics=='variance')
                                    return((apply(as.matrix(data_perm[,seq.int(n1)]),1,var,na.rm=TRUE)/apply(as.matrix(data_perm[,n1+seq.int(n2)]),1,var,na.rm=TRUE))^2)
                                }))
    T_perm=cbind(T_perm,T0)
  }
  
  # Not fully computable p-values (some permutations do not produce any test statistics because of the NAs)
  if(statistics=='variance'){
    no.pval=rowSums(is.na(data))>=(min(n1,n2)-1)
  }else{
    no.pval=rowSums(is.na(data))>=min(n1,n2)
  }
  #T_perm[no.pval,]=NaN # do not compute any p-value when it is not fully computable
  # do not compute any p-value if NaN is in T_perm
  #if(statistics!='variance'){
  #  pval=rowSums(T_perm>=T0)/B
  #}else{
  #  pval=pmin(2*rowSums(T_perm>=T0)/B,2*rowSums(T_perm<=T0)/B)
  #}
  # compute p-value omitting NaN
  if(statistics!='variance'){
    pval=rowSums(T_perm>=T0,na.rm=TRUE)/rowSums(!is.nan(T_perm))
  }else{
    pval=pmin(2*rowSums(T_perm>=T0,na.rm=TRUE)/rowSums(!is.nan(T_perm)),2*rowSums(T_perm<=T0,na.rm=TRUE)/rowSums(!is.nan(T_perm)))
  }
  
  # Combination
  message('    Interval-wise tests...')
  # Asymmetric combination matrix:
  matrix_pval_asymm=matrix(nrow=p,ncol=p)
  matrix_pval_asymm[p,]=pval
  T0_2x=c(T0,T0)
  T_perm_2x=rbind(T_perm,T_perm)
  maxrow=p-max_scale+1
  #message('      Creating the p-value matrix:')
  if(recycle==TRUE){
    for(i in (p-1):maxrow){ # rows
      for(j in seq.int(p)){ # columns
        inf=j
        sup=(p-i)+j
        T0_temp=sum(T0_2x[inf:sup])
        T_temp=colSums(T_perm_2x[inf:sup,])
        # do not compute any p-value if NaN is in T_temp
        #matrix_pval_asymm[i,j]=ifelse(statistics!='variance',sum(T_temp>=T0_temp)/B,min(2*sum(T_temp>=T0_temp)/B,2*sum(T_temp<=T0_temp)/B))
        # compute p-value omitting NaN
        matrix_pval_asymm[i,j]=ifelse(statistics!='variance',sum(T_temp>=T0_temp,na.rm=TRUE)/sum(!is.nan(T_temp)),min(2*sum(T_temp>=T0_temp,na.rm=TRUE)/sum(!is.nan(T_temp)),2*sum(T_temp<=T0_temp,na.rm=TRUE)/sum(!is.nan(T_temp))))
      }
      #message('               end of row ',p-i+1,' out of ',p,'...')
    }
  }else{
    for(i in (p-1):maxrow){ # rows
      for(j in seq.int(i)){ # columns
        inf=j
        sup=(p-i)+j
        T0_temp=sum(T0_2x[inf:sup])
        T_temp=colSums(T_perm_2x[inf:sup,])
        # do not compute any p-value if NaN is in T_temp
        #matrix_pval_asymm[i,j]=ifelse(statistics!='variance',sum(T_temp>=T0_temp)/B,min(2*sum(T_temp>=T0_temp)/B,2*sum(T_temp<=T0_temp)/B))
        # compute p-value omitting NaN
        matrix_pval_asymm[i,j]=ifelse(statistics!='variance',sum(T_temp>=T0_temp,na.rm=TRUE)/sum(!is.nan(T_temp)),min(2*sum(T_temp>=T0_temp,na.rm=TRUE)/sum(!is.nan(T_temp)),2*sum(T_temp<=T0_temp,na.rm=TRUE)/sum(!is.nan(T_temp))))
      }
      #message('               end of row ',p-i+1,' out of ',p,'...')
    }
  }
  corrected.pval.matrix=.pval.correct(matrix_pval_asymm,maxrow)
  corrected.pval=corrected.pval.matrix[maxrow,]
  
  result$T0_plot=rep(NA,length(allNA))
  result$T0_plot[!allNA]=T0_plot
  result$adjusted_pval=rep(NA,length(allNA))
  result$adjusted_pval[!allNA]=corrected.pval
  result$adjusted_pval_matrix=matrix(NA,nrow=length(allNA),ncol=length(allNA))
  result$adjusted_pval_matrix[(sum(allNA)+1):length(allNA),!allNA]=corrected.pval.matrix
  result$unadjusted_pval=rep(NA,length(allNA))
  result$unadjusted_pval[!allNA]=pval
  result$pval_matrix=matrix(NA,nrow=length(allNA),ncol=length(allNA))
  result$pval_matrix[(sum(allNA)+1):length(allNA),!allNA]=matrix_pval_asymm
  result$exact=exact
  result$index.data2=index.data2
  class(result)='ITWomics.2pop'
  return(list(result=result,no.pval=no.pval))
}

plotSummary_reproducible <- function(regionsFeatures,alpha=0.05,only_significant=FALSE,scale_threshold=NULL,nlevel=100,
                                     groupby='test',test=1:nTests(regionsFeatures),gaps_tests=NULL,
                                     reproducible_pval_matrix=NULL, # NEW matrix of reproducible p-values (TRUE or FALSE).
                                     id_features_subset=idFeaturesTest(regionsFeatures),gaps_features=NULL,
                                     ask=TRUE,filenames=NA,align_lab=NA,cellwidth=NA,cellheight=NA,
                                     xlab='Windows',ylab=ifelse(groupby=='test','Features','Tests'),...){
  # MODIFIED FROM BIOCONDUCTOR PACKAGE IWTomics
  # FUNCTION plotSummary
  
  if(class(regionsFeatures)!='IWTomicsData')
    stop('invalid regionsFeatures. IWTomicsData object expected.')
  if(!validObject(regionsFeatures))
    stop('invalid regionsFeatures.')
  if(is.null(regionsFeatures@test))
    stop('No test results present in regionFeatures.')
  if(!(groupby %in% c('test','feature')))
    stop('invalid groupby \'',groupby,'\'. Available groupby are \'test\' and \'feature\'.')
  if(sum(!(test %in% 1:nTests(regionsFeatures))))
    stop('invalid test number.')
  if(sum(!(id_features_subset %in% idFeaturesTest(regionsFeatures))))
    stop('invalid id_features_subset. The features provided are not listed in regionsFeatures.')
  if(is.null(scale_threshold)){
    scale_threshold=lapply(.testResults(regionsFeatures)[test],
                           function(test){
                             unlist(lapply(test[id_features_subset],function(feature) feature$max_scale))
                           })
  }else{
    if(is.list(scale_threshold)){
      scale_threshold=as.list(rep(scale_threshold,length.out=length(test)))
      scale_threshold=lapply(scale_threshold,function(scale_threshold) rep(scale_threshold,length.out=length(id_features_subset)))
    }else{
      scale_threshold=lapply(seq_along(test),function(test) return(rep(scale_threshold,length.out=length(id_features_subset))))
    }
    scale_threshold=lapply(scale_threshold,
                           function(scale_threshold){
                             names(scale_threshold)=id_features_subset
                             return(scale_threshold)})
  }
  if(alignment(regionsFeatures)=='scale')
    align_lab=NA
  
  devAskNewPage(ask)
  if(groupby=='test'){
    if(length(unique(resolution(regionsFeatures)[id_features_subset]))!=1)
      stop('groupby \'test\' but selected features with different resolution. Smooth data first to have the same resolution.')
    if(!is.na(filenames)[1])
      if(length(filenames)<length(test))
        stop('filenames should contain one file path for each test as defined by groupby.')
    if(!is.null(gaps_features))
      if(FALSE %in% (gaps_features %in% seq_along(id_features_subset)))
        stop('invalid gap_features.')
    for(i in test){
      filename_i=filenames[i]
      scale_threshold_i=scale_threshold[[i]]
      id_region1_i=idRegionsTest(regionsFeatures,i)[[1]][1]
      id_region2_i=idRegionsTest(regionsFeatures,i)[[1]][2]
      if((!is.null(id_region2_i))&&(id_region2_i==''))
        id_region2_i=NULL
      id_features_subset_i=id_features_subset
      gaps_features_i=gaps_features
      result_i=.testResults(regionsFeatures,i,id_features_subset_i)[[1]]
      pval_scale_threshold=Reduce(cbind,adjusted_pval(regionsFeatures,i,id_features_subset_i,scale_threshold_i)[[1]])
      # Significant p-value
      significant=pval_scale_threshold
      significant[significant<1e-5]=1e-5
      significant[significant>alpha]=1
      if(only_significant){
        if(!is.null(gaps_features_i)){
          with_gaps=rep(NA,length(id_features_subset_i)+length(gaps_features_i))
          with_gaps[-(gaps_features_i+seq_along(gaps_features_i))]=id_features_subset_i
        }
        id_features_subset_i=id_features_subset_i[colSums(significant<=alpha,na.rm=TRUE)>0]
        significant=significant[,id_features_subset_i]
        result_i=result_i[id_features_subset_i]
        scale_threshold_i=scale_threshold_i[id_features_subset_i]
        if(!is.null(gaps_features_i)){
          with_gaps=with_gaps[with_gaps %in% c(id_features_subset_i,NA)]
          gaps_features_i=which(is.na(with_gaps))-seq_len(sum(is.na(with_gaps)))
          gaps_features_i=gaps_features_i[gaps_features_i>0]
          if(length(gaps_features_i)==0)
            gaps_features_i=NULL
        }
      }
      log_significant=as.matrix(-log10(significant))
      colnames(log_significant)=nameFeatures(regionsFeatures)[id_features_subset_i]
      T0_plot=as.matrix(Reduce(cbind,lapply(result_i,function(feature) feature$T0_plot)))
      colnames(T0_plot)=id_features_subset_i
      # Recode the T0_plot matrix in 1 and -1
      if(testInput(regionsFeatures)$statistics=='variance'){
        T0_plot[T0_plot<1]=-1
        T0_plot[T0_plot>1]=1
      }else{
        T0_plot[T0_plot<0]=-1
        T0_plot[T0_plot>=0]=1
      }
      # Put a sign to pvalues based on T0_plot
      log_significant=log_significant*T0_plot
      
      if(alignment(regionsFeatures)=='left')
        row.names(log_significant)=seq_len(nrow(log_significant))
      if(alignment(regionsFeatures)=='right')
        row.names(log_significant)=-(nrow(log_significant):1)
      if(alignment(regionsFeatures)=='center')
        if(nrow(log_significant)%%2!=0){
          row.names(log_significant)=seq_len(nrow(log_significant))-(nrow(log_significant)+1)/2
        }else{
          row.names(log_significant)=setdiff(seq(-nrow(log_significant)/2,nrow(log_significant)/2),0)
        }
      if(alignment(regionsFeatures)=='scale')
        row.names(log_significant)=round(seq(0,1,length.out=nrow(log_significant)),2)
      
      if(result_i[[1]]$test=="1pop"){
        main=paste0(nameRegions(regionsFeatures)[id_region1_i])
      }else{
        main=paste0(nameRegions(regionsFeatures)[id_region1_i],' vs ',nameRegions(regionsFeatures)[id_region2_i])
      }
      .pheatmap(t(log_significant),gaps_row=gaps_features_i,border_color="grey60",filename=filename_i,
                breaks=seq(-5,5,length.out=nlevel),color=colorRampPalette(c("navy","white","red"))(n=nlevel-1),
                legend.main='-log10(p-value)',main=main,cex.main=2,xlab=xlab,ylab=ylab,thresholds=scale_threshold_i,
                zero_lab=align_lab,zero_lab_pos=alignment(regionsFeatures),cellwidth=cellwidth,cellheight=cellheight,...)
    }
  } else {
    if(!is.na(filenames[1]))
      if(length(filenames)<length(id_features_subset))
        stop('filenames should contain one file path for each test as defined by groupby.')
    if(!is.na(filenames)[1])
      if(length(filenames)<length(id_features_subset))
        stop('filenames should contain one file path for each test as defined by groupby.')
    if(!is.null(gaps_tests))
      if(FALSE %in% (gaps_tests %in% 1:nTests(regionsFeatures)))
        stop('invalid gap_tests.')
    for(i in seq_along(id_features_subset)){
      filename_i=filenames[i]
      scale_threshold_i=lapply(scale_threshold,function(scale_threshold) scale_threshold[id_features_subset[i]])
      id_region1_i=testInput(regionsFeatures)$id_region1[test]
      id_region2_i=testInput(regionsFeatures)$id_region2[test]
      if(is.null(id_region2_i))
        id_region2_i=rep('',length(id_region1_i))
      id_features_subset_i=id_features_subset[i]
      gaps_tests_i=gaps_tests
      gaps_tests_i=match(gaps_tests_i,test)
      gaps_tests_i=gaps_tests_i[!is.na(gaps_tests_i)]
      result_i=.testResults(regionsFeatures,test,id_features_subset_i)
      ### start NEW ###
      if(!is.null(reproducible_pval_matrix)){
        pval_scale_threshold=mapply(function(test,pval){
          pval_scale_threshold=rep(NA,length(test[[1]]$notNA))
          pval_scale_threshold[test[[1]]$notNA]=pval[[1]]
          return(pval_scale_threshold)
        },result_i,adjusted_pval_reproducible(regionsFeatures,test,id_features_subset_i,scale_threshold_i,reproducible_pval_matrix))
      }else{
        pval_scale_threshold=mapply(function(test,pval){
          pval_scale_threshold=rep(NA,length(test[[1]]$notNA))
          pval_scale_threshold[test[[1]]$notNA]=pval[[1]]
          return(pval_scale_threshold)
        },result_i,adjusted_pval(regionsFeatures,test,id_features_subset_i,scale_threshold_i))
      }
      ### end NEW ###
      
      # Significant p-value
      significant=pval_scale_threshold
      significant[significant<1e-5]=1e-5
      if(is.null(reproducible_pval_matrix))
        significant[significant>alpha]=1
      if(only_significant){
        if(!is.null(gaps_tests_i)){
          with_gaps=rep(NA,length(id_region1_i)+length(gaps_tests_i))
          with_gaps[-(gaps_tests_i+seq_along(gaps_tests_i))]=seq_along(id_region1_i)
        }
        index=which(colSums(significant<=alpha,na.rm=TRUE)>0)
        id_region1_i=id_region1_i[index]
        id_region2_i=id_region2_i[index]
        significant=significant[,index]
        result_i=result_i[index]
        scale_threshold_i=scale_threshold_i[index]
        if(!is.null(gaps_tests_i)){
          with_gaps=with_gaps[with_gaps %in% c(index,NA)]
          gaps_tests_i=which(is.na(with_gaps))-seq_len(sum(is.na(with_gaps)))
          gaps_tests_i=gaps_tests_i[gaps_tests_i>0]
          if(length(gaps_tests_i)==0)
            gaps_tests_i=NULL
        }
      }
      log_significant=as.matrix(-log10(significant))
      colnames(log_significant)[id_region2_i!='']=nameRegions(regionsFeatures)[id_region1_i[id_region2_i!='']]
      colnames(log_significant)[id_region2_i=='']=nameRegions(regionsFeatures)[id_region1_i[id_region2_i=='']]
      T0_plot=as.matrix(Reduce(cbind,lapply(result_i,
                                            function(test){
                                              T0_plot=rep(NA,length(test[[1]]$notNA))
                                              T0_plot[test[[1]]$notNA]=test[[1]]$T0_plot
                                              return(T0_plot)
                                            })))
      colnames(T0_plot)[id_region2_i!='']=nameRegions(regionsFeatures)[id_region1_i[id_region2_i!='']]
      colnames(T0_plot)[id_region2_i=='']=nameRegions(regionsFeatures)[id_region1_i[id_region2_i=='']]
      # Recode the T0_plot matrix in 1 and -1
      if(testInput(regionsFeatures)$statistics=='variance'){
        T0_plot[T0_plot<1]=-1
        T0_plot[T0_plot>1]=1
      }else{
        T0_plot[T0_plot<0]=-1
        T0_plot[T0_plot>=0]=1
      }
      # Put a sign to pvalues based on T0_plot
      log_significant=log_significant*T0_plot
      
      #if(alignment(regionsFeatures)=='left')
      #  row.names(log_significant)=seq_len(nrow(log_significant))
      #if(alignment(regionsFeatures)=='right')
      #  row.names(log_significant)=-(nrow(log_significant):1)
      #if(alignment(regionsFeatures)=='center')
      #  if(nrow(log_significant)%%2!=0){
      #    row.names(log_significant)=seq_len(nrow(log_significant))-(nrow(log_significant)+1)/2
      #  }else{
      #    row.names(log_significant)=setdiff(seq(-nrow(log_significant)/2,nrow(log_significant)/2),0)
      #  }
      #if(alignment(regionsFeatures)=='scale')
      #  row.names(log_significant)=round(seq(0,1,length.out=nrow(log_significant)),2)
      
      main=paste0(nameFeatures(regionsFeatures)[id_features_subset_i])
      .pheatmap(t(log_significant),gaps_row=gaps_tests_i,border_color="grey60",filename=filename_i,
                breaks=seq(-5,5,length.out=nlevel),color=colorRampPalette(c("navy","white","red"))(n=nlevel-1),
                legend.main='-log10(p-value)',main=main,cex.main=2,xlab=xlab,ylab=ylab,
                zero_lab=align_lab,zero_lab_pos=alignment(regionsFeatures),cellwidth=cellwidth,cellheight=cellheight,...)
    }
  }
}

adjusted_pval_reproducible <- function(x,test,id_features_subset,scale_threshold,reproducible_pval_matrix){
  # MODIFIED FROM BIOCONDUCTOR PACKAGE IWTomics
  # FUNCTION adjusted_pval
  
  if(nTests(x)==0)
    return(NULL)
  if(sum(!(test %in% 1:nTests(x))))
    stop('invalid test number.')
  if(sum(!(id_features_subset %in% idFeaturesTest(x))))
    stop('invalid id_features_subset.')
  if(is.list(scale_threshold)){
    scale_threshold=lapply(scale_threshold,
                           function(scale){
                             scale=as.list(rep(scale,length.out=length(id_features_subset)))
                             names(scale)=id_features_subset
                             return(scale)
                           })
  }else{
    scale_threshold=lapply(test,
                           function(i){
                             scale=as.list(rep(scale_threshold,length.out=length(id_features_subset)))
                             names(scale)=id_features_subset
                             return(scale)
                           })
  }
  pval=mapply(function(results,scale,reproducible_pval) mapply(function(result,scale){
    pval=result$adjusted_pval_matrix
    if((scale<1)||(result$max_scale<scale)){
      warning('invalid scale_threshold. Setting it to the default value.',call.=FALSE,immediate.=TRUE)
      scale=result$max_scale
    }
    pval=pval[ncol(pval)-scale+1,]
    pval[!reproducible_pval[ncol(reproducible_pval)-scale+1,]]=1
    return(pval)
  },results,scale,SIMPLIFY=FALSE),
  .testResults(x,test,id_features_subset),scale_threshold,reproducible_pval_matrix,SIMPLIFY=FALSE)
  names(pval)=paste0('test',test)
  return(pval)
}

# INTERNAL FUNCTIONS FROM BIOCONDUCTOR PACKAGE IWTomics
setGeneric(".testResults",function(x,test,id_features_subset,...) standardGeneric(".testResults"))
setMethod(".testResults",c("IWTomicsData","vector","character"),
          function(x,test,id_features_subset){
            if(nTests(x)==0)
              return(NULL)
            if(sum(!(test %in% 1:nTests(x))))
              stop('invalid test number.')
            if(sum(!(id_features_subset %in% idFeaturesTest(x))))
              stop('invalid id_features_subset.')
            lapply(x@test$result[test],function(results) results[id_features_subset])
          })
setMethod(".testResults",c("IWTomicsData","missing","character"),
          function(x,test,id_features_subset){
            if(nTests(x)==0)
              return(NULL)
            .testResults(x,1:nTests(x),id_features_subset)
          })
setMethod(".testResults",c("IWTomicsData","vector","missing"),
          function(x,test,id_features_subset){
            if(nTests(x)==0)
              return(NULL)
            .testResults(x,test,idFeaturesTest(x))
          })
setMethod(".testResults",c("IWTomicsData","missing","missing"),
          function(x){
            if(nTests(x)==0)
              return(NULL)
            .testResults(x,1:nTests(x),idFeaturesTest(x))
          })

setGeneric(".testResults",function(x,test,id_features_subset,...) standardGeneric(".testResults"))
setMethod(".testResults",c("IWTomicsData","vector","character"),
          function(x,test,id_features_subset){
            if(nTests(x)==0)
              return(NULL)
            if(sum(!(test %in% 1:nTests(x))))
              stop('invalid test number.')
            if(sum(!(id_features_subset %in% idFeaturesTest(x))))
              stop('invalid id_features_subset.')
            lapply(x@test$result[test],function(results) results[id_features_subset])
          })
setMethod(".testResults",c("IWTomicsData","missing","character"),
          function(x,test,id_features_subset){
            if(nTests(x)==0)
              return(NULL)
            .testResults(x,1:nTests(x),id_features_subset)
          })
setMethod(".testResults",c("IWTomicsData","vector","missing"),
          function(x,test,id_features_subset){
            if(nTests(x)==0)
              return(NULL)
            .testResults(x,test,idFeaturesTest(x))
          })
setMethod(".testResults",c("IWTomicsData","missing","missing"),
          function(x){
            if(nTests(x)==0)
              return(NULL)
            .testResults(x,1:nTests(x),idFeaturesTest(x))
          })

.pheatmap <- function (mat, color = colorRampPalette(c("navy","white","red"))(100), 
                       breaks = NA, border_color = "grey60", 
                       cellwidth = NA, cellheight = NA, scale = "none", 
                       legend = TRUE, legend.main = '',legend_breaks = NA, legend_labels = NA, 
                       show_rownames = TRUE, show_colnames = TRUE, main = NA, fontsize = 10, 
                       fontsize_row = fontsize, fontsize_col = fontsize, display_numbers = FALSE, 
                       number_format = "%.2f", number_color = "grey30", fontsize_number = 0.8 * fontsize, 
                       gaps_row = NULL, gaps_col = NULL, labels_row = NULL, 
                       labels_col = NULL, filename = NA, width = NA, height = NA, 
                       silent = FALSE, xlab = NA, ylab = NA, cex.main = 1.3, thresholds = NA, zero_lab = NA, zero_lab_pos = "center", ...) {
  # modified from pheatmap package
  require(grid)
  require(gtable)
  
  # Set labels
  if (is.null(labels_row)) {
    labels_row = rownames(mat)
  }
  if (is.null(labels_col)) {
    labels_col = colnames(mat)
  }
  
  # Preprocess matrix
  mat = as.matrix(mat)
  if (scale != "none") {
    mat = .scale_mat(mat, scale)
    if (.is.na2(breaks)) {
      breaks = .generate_breaks(mat, length(color), center = TRUE)
    }
  }
  
  # Format numbers to be displayed in cells
  if (is.matrix(display_numbers) | is.data.frame(display_numbers)) {
    if (nrow(display_numbers) != nrow(mat) | ncol(display_numbers) != ncol(mat)) {
      stop("If display_numbers provided as matrix, its dimensions have to match with mat")
    }
    display_numbers = as.matrix(display_numbers)
    fmat = matrix(as.character(display_numbers), nrow = nrow(display_numbers), ncol = ncol(display_numbers))
    fmat_draw = TRUE
  } else {
    if (display_numbers) {
      fmat = matrix(sprintf(number_format, mat), nrow = nrow(mat), ncol = ncol(mat))
      fmat_draw = TRUE
    } else {
      fmat = matrix(NA, nrow = nrow(mat), ncol = ncol(mat))
      fmat_draw = FALSE
    }
  }
  
  attr(fmat, "draw") = fmat_draw
  
  # Colors and scales
  if (!.is.na2(legend_breaks) & !.is.na2(legend_labels)) {
    if (length(legend_breaks) != length(legend_labels)) {
      stop("Lengths of legend_breaks and legend_labels must be the same")
    }
  }
  if (.is.na2(breaks)) {
    breaks = .generate_breaks(as.vector(mat), length(color))
  }
  if (legend & .is.na2(legend_breaks)) {
    legend = grid.pretty(range(as.vector(breaks)))
    names(legend) = legend
  } else if (legend & !.is.na2(legend_breaks)) {
    legend = legend_breaks[legend_breaks >= min(breaks) & legend_breaks <= max(breaks)]
    if (!.is.na2(legend_labels)) {
      legend_labels = legend_labels[legend_breaks >= min(breaks) & legend_breaks <= max(breaks)]
      names(legend) = legend_labels
    } else {
      names(legend) = legend
    }
  } else {
    legend = NA
  }
  mat = .scale_colours(mat, col = color, breaks = breaks)
  
  if (!show_rownames) {
    labels_row = NULL
  }
  
  if (!show_colnames) {
    labels_col = NULL
  }
  
  # Draw heatmap
  gt = .heatmap_motor(mat, border_color = border_color, cellwidth = cellwidth, 
                      cellheight = cellheight, filename = filename, width = width, 
                      height = height, breaks = breaks, color = color, legend = legend, legend.main = legend.main,
                      main = main, fontsize = fontsize, fontsize_row = fontsize_row, 
                      fontsize_col = fontsize_col, fmat = fmat, fontsize_number = fontsize_number, 
                      number_color = number_color, gaps_row = gaps_row, gaps_col = gaps_col, thresholds = thresholds, zero_lab = zero_lab, zero_lab_pos = zero_lab_pos, 
                      labels_row = labels_row, labels_col = labels_col, xlab = xlab, ylab = ylab, cex.main = cex.main,...)
  
  if (is.na(filename) & !silent) {
    grid.newpage()
    grid.draw(gt)
  }
  
  invisible(list(gtable = gt))
}

.lo <- function (rown, coln, nrow, ncol, cellheight = NA, cellwidth = NA, 
                 legend, legend.main, main, fontsize, fontsize_row, thresholds, zero_lab, 
                 fontsize_col, gaps_row, gaps_col, xlab, ylab, cex.main = cex.main, ...) {
  # modified from pheatmap package
  
  # Get height of colnames and length of rownames
  if (!is.null(coln[1])) {
    longest_coln = which.max(strwidth(coln, units = "in"))
    gp = list(fontsize = fontsize_col, ...)
    coln_height = unit(1, "grobheight", textGrob(coln[longest_coln], rot = 90, gp = do.call(gpar, gp))) + unit(10, "bigpts")
  } else {
    coln_height = unit(5, "bigpts")
  }
  
  if (!is.null(rown[1])) {
    longest_rown = which.max(strwidth(rown, units = "in"))
    gp = list(fontsize = fontsize_row, ...)
    rown_width = unit(1, "grobwidth", textGrob(rown[longest_rown], gp = do.call(gpar, gp))) + unit(5, "bigpts")
  } else {
    rown_width = unit(5, "bigpts")
  }
  
  gp = list(fontsize = fontsize, ...)
  # Axis label positions
  if (!.is.na2(xlab)) {
    xlab_height = unit(1, "grobheight", textGrob(xlab, gp = do.call(gpar, gp))) + unit(5, "bigpts")
  } else {
    xlab_height = unit(1, "bigpts")
  }
  if (!.is.na2(ylab)) {
    ylab_width = unit(1, "grobheight", textGrob(ylab, rot=270, gp = do.call(gpar, gp))) + unit(1, "bigpts")
  } else {
    ylab_width = unit(5, "bigpts")
  }
  
  # Threshold position
  if (!.is.na2(thresholds)){
    thresholds_width = unit(1.1, "grobwidth", textGrob("Threshold", gp = gpar(...)))+unit(2, "bigpts")
    thresholds_lab_height = unit(1, "grobheight", textGrob("Threshold", just = "bottom", gp = do.call(gpar, gp))) + unit(5, "bigpts")
  } else {
    thresholds_width = unit(15, "bigpts")
    thresholds_lab_height = unit(0, "bigpts")
  }
  
  # Zero label position
  if(!.is.na2(zero_lab)){
    zero_lab_height = unit(1, "grobheight", textGrob("Threshold", just = "bottom", gp = do.call(gpar, gp))) + unit(5, "bigpts")
  } 
  
  # Legend position
  if (!.is.na2(legend)) {
    longest_break = which.max(nchar(names(legend)))
    longest_break = unit(1.1, "grobwidth", textGrob(as.character(names(legend))[longest_break], gp = do.call(gpar, gp)))
    legend_width = unit(25, "bigpts") + longest_break * 1.2
  } else {
    legend_width = unit(0, "bigpts")
  }
  if (!.is.na2(legend.main)) {
    legend_main_width = unit(1.1, "grobwidth", textGrob(legend.main, rot=90, gp = gpar(...)))+unit(6, "bigpts")
  } else {
    legend_main_width = unit(0, "bigpts")
  }
  
  # Set main title height
  if (is.na(main)) {
    main_height = unit(0, "npc")
  } else {
    main_height = unit(2.2, "grobheight", textGrob(main, gp = gpar(fontsize = cex.main * fontsize, ...)))
  }
  textheight = unit(fontsize, "bigpts")
  
  # Set cell sizes
  if (is.na(cellwidth)) {
    mat_width = unit(1, "npc") - thresholds_width - rown_width - ylab_width - legend_width - legend_main_width
  } else {
    mat_width = unit(cellwidth * ncol, "bigpts") + length(gaps_col) * unit(10, "bigpts")
  }
  
  if (is.na(cellheight)) {
    mat_height = unit(1, "npc") - main_height - thresholds_lab_height - coln_height - xlab_height 
  } else {
    mat_height = unit(cellheight * nrow, "bigpts") + length(gaps_row) * unit(10, "bigpts")
  }
  
  # Produce gtable
  gt = gtable(widths = unit.c(unit(0, "bigpts"), unit(0, "bigpts"), thresholds_width, mat_width, rown_width, ylab_width, legend_width, legend_main_width, unit(0, "bigpts")), 
              heights = unit.c(main_height, unit(0, "bigpts"), unit(0, "bigpts"), thresholds_lab_height, mat_height, coln_height, zero_lab_height, xlab_height), vp = viewport(gp = do.call(gpar,  gp)))
  cw = convertWidth(mat_width - (length(gaps_col) * unit(10, "bigpts")), "bigpts", valueOnly = TRUE)/ncol
  ch = convertHeight(mat_height - (length(gaps_row) * unit(10, "bigpts")), "bigpts", valueOnly = TRUE)/nrow
  
  # Return minimal cell dimension in bigpts to decide if borders are drawn
  mindim = min(cw, ch)
  
  res = list(gt = gt, mindim = mindim)
  
  return(res)
}

.find_coordinates <- function (n, gaps, m = 1:n) {
  if (length(gaps) == 0) {
    return(list(coord = unit(m/n, "npc"), size = unit(1/n, "npc")))
  }
  
  if (max(gaps) > n) {
    stop("Gaps do not match with matrix size")
  }
  
  size = (1/n) * (unit(1, "npc") - length(gaps) * unit("10", "bigpts"))
  gaps2 = apply(sapply(gaps, function(gap, x) { x > gap }, m), 1, sum)
  
  coord = m * size + (gaps2 * unit("10", "bigpts"))
  return(list(coord = coord, size = size))
}

.draw_matrix=function (matrix, border_color, gaps_rows, gaps_cols, fmat, fontsize_number, number_color, zero_lab, zero_lab_pos) {
  # modified from pheatmap package
  
  n = nrow(matrix)
  m = ncol(matrix)
  
  coord_x = .find_coordinates(m, gaps_cols)
  coord_y = .find_coordinates(n, gaps_rows)
  
  x = coord_x$coord - 0.5 * coord_x$size
  y = unit(1, "npc") - (coord_y$coord - 0.5 * coord_y$size)
  
  res = gList()
  for(i in seq_len(nrow(matrix))){
    index_i=which(!is.na(matrix[i,]))
    x_i = (coord_x$coord - 0.5 * coord_x$size)[index_i]
    y_i = y[i]
    coord = expand.grid(y = y_i, x = x_i)
    res[[paste0("rect",i)]] = rectGrob(x = coord$x, y = coord$y, width = coord_x$size, 
                                       height = coord_y$size, gp = gpar(fill = matrix[i,][index_i], col = border_color))
  }
  
  if (attr(fmat, "draw")) {
    res[["text"]] = textGrob(x = coord$x, y = coord$y, label = fmat, gp = gpar(col = number_color, fontsize = fontsize_number))
  }
  
  if(!.is.na2(zero_lab)) {
    if(zero_lab_pos=='center')
      x = rep(unit(0.5, "npc"),2)
    if(zero_lab_pos=='right')
      x = rep(unit(1, "npc"),2)
    if(zero_lab_pos=='left')
      x = rep(unit(0, "npc"),2)
    
    #res[["zero_line"]] = linesGrob(x = x, y = unit(c(0,1), "npc")+unit(c(0,3), "bigpts"),gp=gpar(lwd = 2))
    res[["zero_line"]] = linesGrob(x = x, y = y[c(length(y),1)]+c(-0.5,0.5)*coord_y$size+unit(c(-3,0), "bigpts"),gp=gpar(lwd = 2.2))
  }
  
  res = gTree(children = res)
  
  return(res)
}

.draw_colnames <- function (coln, gaps, ...) {
  # modified from pheatmap package
  
  coord = .find_coordinates(length(coln), gaps)
  x = coord$coord - 0.5 * coord$size
  
  res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3, "bigpts"), just='right', rot = 90, gp = gpar(...))
  
  return(res)
}

.draw_rownames <- function (rown, gaps, ...) {
  # modified from pheatmap package
  
  coord = .find_coordinates(length(rown), gaps)
  y = unit(1, "npc") - (coord$coord - 0.5 * coord$size)
  
  res = textGrob(rown, x = unit(3, "bigpts"), y = y, just='left', gp = gpar(...))
  
  return(res)
}

.draw_thresholds <- function (thresholds, gaps, ...) {
  
  coord = .find_coordinates(length(thresholds), gaps)
  y = unit(1, "npc") - (coord$coord - 0.5 * coord$size)
  
  res = textGrob(thresholds, x =unit(0.9, "npc") ,y = y, just='right', gp = gpar(...))
  
  return(res)
}

.draw_thresholds_lab <- function(thresholds_lab, ...) {
  
  res = textGrob(thresholds_lab, x = unit(0.9, "npc"), y = unit(15, "bigpts"), just=c('right','top'), gp = gpar(...))
  
  return(res)
}

.draw_zero_lab <- function(zero_lab, zero_lab_pos, ...) {
  
  if(zero_lab_pos=='center')
    x = unit(0.5, "npc")
  if(zero_lab_pos=='right')
    x = unit(1, "npc")
  if(zero_lab_pos=='left')
    x = unit(0, "npc")
  
  res = textGrob(zero_lab, x = x, y = unit(0, "bigpts"), just=c(zero_lab_pos,'bottom'), gp = gpar(...))
  
  return(res)
}

.draw_legend <- function (color, breaks, legend,...) {
  # modified from pheatmap package
  
  height = min(unit(1, "npc"), unit(150, "bigpts"))
  
  legend_pos = (legend - min(breaks)) / (max(breaks) - min(breaks))
  legend_pos = height * legend_pos + (unit(1, "npc") - height)
  
  breaks = (breaks - min(breaks))/(max(breaks) - min(breaks))
  breaks = height * breaks + (unit(1, "npc") - height)
  
  h = breaks[-1] - breaks[-length(breaks)]
  
  rect = rectGrob(x = 0, y = breaks[-length(breaks)], width = unit(10, "bigpts"), height = h, 
                  hjust = 0, vjust = 0, gp = gpar(fill = color, col = "#FFFFFF00"))
  text = textGrob(names(legend), x = unit(14, "bigpts"), y = legend_pos, hjust = 0, gp = gpar(...))
  
  res = grobTree(rect, text)
  return(res)
}

.draw_legend_main <- function (breaks, legend, legend.main,...) {
  height = min(unit(1, "npc"), unit(150, "bigpts"))
  
  legend_main_pos=(sum(range(legend))/2-min(breaks))/(max(breaks) - min(breaks))
  legend_main_pos=height * legend_main_pos + (unit(1, "npc") - height)
  
  label = textGrob(legend.main, x = 0, y = legend_main_pos, just='centre', rot=270, gp = gpar(...))
  
  res = grobTree(label)
  return(res)
}

.draw_main = function(text, ...){
  # modified from pheatmap package
  
  res = textGrob(text, gp = gpar(fontface = "bold", ...))
  
  return(res)
}

.draw_xlab = function(text, ...){
  
  res = textGrob(text, just = "bottom", gp = gpar(...))
  
  return(res)
}

.draw_ylab = function(text, ...){
  
  res = textGrob(text, just = "center", rot = 270, gp = gpar(...))
  
  return(res)
}

.heatmap_motor <- function (matrix, border_color, cellwidth, cellheight, filename, width, height, 
                            breaks, color, legend, legend.main, main, fontsize, fontsize_row, 
                            fontsize_col, fmat, fontsize_number, number_color, gaps_col, thresholds, zero_lab, zero_lab_pos,
                            gaps_row, labels_row, labels_col, xlab, ylab, cex.main,...) {
  # modified from pheatmap package
  
  # Set layout
  lo = .lo(coln = labels_col, rown = labels_row, nrow = nrow(matrix), 
           ncol = ncol(matrix), cellwidth = cellwidth, cellheight = cellheight, 
           legend = legend, legend.main=legend.main, 
           main = main, fontsize = fontsize, fontsize_row = fontsize_row, 
           thresholds = thresholds, zero_lab = zero_lab,
           fontsize_col = fontsize_col, gaps_row = gaps_row, gaps_col = gaps_col, 
           xlab = xlab, ylab = ylab, cex.main = cex.main, ...)
  
  res = lo$gt
  mindim = lo$mindim
  
  if (!is.na(filename)) {
    if (is.na(height)) {
      height = convertHeight(gtable_height(res), "inches", valueOnly = TRUE)
    }
    if (is.na(width)) {
      width = convertWidth(gtable_width(res), "inches", valueOnly = TRUE)
    }
    # Get file type
    r = regexpr("\\.[a-zA-Z]*$", filename)
    if (r == -1) 
      stop("Improper filename")
    ending = substr(filename, r + 1, r + attr(r, "match.length"))
    
    f = switch(ending, pdf = function(x, ...) pdf(x, ...), 
               png = function(x, ...) png(x, units = "in", res = 300, ...), 
               jpeg = function(x, ...) jpeg(x, units = "in", res = 300, ...), 
               jpg = function(x, ...) jpeg(x, units = "in", res = 300, ...), 
               tiff = function(x, ...) tiff(x, units = "in", res = 300, compression = "lzw", ...), 
               bmp = function(x, ...) bmp(x, units = "in", res = 300, ...), 
               stop("File type should be: pdf, png, bmp, jpg, tiff"))
    
    f(filename, height = height, width = width)
    gt = .heatmap_motor(matrix, cellwidth = cellwidth, cellheight = cellheight, 
                        border_color = border_color, breaks = breaks, 
                        color = color, legend = legend, legend.main=legend.main, filename = NA, 
                        main = main, fontsize = fontsize, fontsize_row = fontsize_row, 
                        fontsize_col = fontsize_col, fmat = fmat, fontsize_number = fontsize_number, thresholds = thresholds, zero_lab = zero_lab, zero_lab_pos = zero_lab_pos,
                        number_color = number_color, labels_row = labels_row, xlab = xlab, ylab = ylab, cex.main = cex.main, 
                        labels_col = labels_col, gaps_col = gaps_col, gaps_row = gaps_row, ...)
    grid.draw(gt)
    dev.off()
    
    return(gt)
  }
  
  # Omit border color if cell size is too small 
  if (mindim < 3) 
    border_color = NA
  
  # Draw title
  if (!is.na(main)) {
    elem = .draw_main(main, fontsize = cex.main * fontsize, ...)
    res = gtable_add_grob(res, elem, t = 1, l = 4, name = "main")
  }
  
  # Draw thresholds
  if (!.is.na2(thresholds)) {
    elem = .draw_thresholds(thresholds, gaps_row, fontsize = fontsize)
    res = gtable_add_grob(res, elem, t = 5, l = 3, clip = "off", name = "thresholds")
    elem = .draw_thresholds_lab("Threshold", fontsize = fontsize)
    res = gtable_add_grob(res, elem, t = 4, l = 3, clip = "off", name = "thresholds_lab")
  } else {
    elem = .draw_thresholds(rep('',nrow(matrix)), gaps_row, fontsize = fontsize)
    res = gtable_add_grob(res, elem, t = 5, l = 3, clip = "off", name = "thresholds")
  }
  
  # Draw matrix
  elem = .draw_matrix(matrix, border_color, gaps_row, gaps_col, fmat, fontsize_number, number_color, zero_lab, zero_lab_pos)
  res = gtable_add_grob(res, elem, t = 5, l = 4, clip = "off", name = "matrix")
  
  # Draw zero label
  if (!.is.na2(zero_lab)) {
    elem = .draw_zero_lab(zero_lab, zero_lab_pos, fontsize = fontsize)
    res = gtable_add_grob(res, elem, t = 7, l = 4, clip = "off", name = "zero_lab")
  }
  
  # Draw colnames
  if (length(labels_col) != 0) {
    pars = list(labels_col, gaps = gaps_col, fontsize = fontsize_col, ...)
    elem = do.call(.draw_colnames, pars)
    res = gtable_add_grob(res, elem, t = 6, l = 4, clip = "off", name = "col_names")
  }
  
  # Draw rownames
  if (length(labels_row) != 0) {
    pars = list(labels_row, gaps = gaps_row, fontsize = fontsize_row,  ...)
    elem = do.call(.draw_rownames, pars)
    res = gtable_add_grob(res, elem, t = 5, l = 5, clip = "off", name = "row_names")
  }
  
  # Draw xlab
  if (!.is.na2(xlab)){
    elem = .draw_xlab(xlab,fontsize = fontsize)
    res = gtable_add_grob(res, elem, t = 7, l = 4, clip = "off", name = "xlab")
  }
  
  # Draw ylab
  if (!.is.na2(ylab)){
    elem = .draw_ylab(ylab,fontsize = fontsize)
    res = gtable_add_grob(res, elem, t = 5, l = 6, clip = "off", name = "ylab")
  }
  
  # Draw legend
  if (!.is.na2(legend)) {
    elem = .draw_legend(color, breaks, legend, fontsize = fontsize, ...)
    #t = ifelse(is.null(labels_row), 5, 4)
    res = gtable_add_grob(res, elem, t = 5, l = 7, b = 5, clip = "off", name = "legend")
  }
  if (!.is.na2(legend.main)) {
    elem = .draw_legend_main(breaks, legend, legend.main, fontsize = fontsize, ...)
    #t = ifelse(is.null(labels_row), 4, 3)
    res = gtable_add_grob(res, elem, t = 5, l = 8, b = 5, clip = "off", name = "legend_main")
  }
  
  return(res)
}

.generate_breaks = function(x, n, center = FALSE){
  # modified from pheatmap package
  
  if(center){
    m = max(abs(c(min(x, na.rm = TRUE), max(x, na.rm = TRUE))))
    res = seq(-m, m, length.out = n + 1)
  }
  else{
    res = seq(min(x, na.rm = TRUE), max(x, na.rm = TRUE), length.out = n + 1)
  }
  
  return(res)
}

.scale_vec_colours <- function (x, col, breaks = NA){
  # modified from pheatmap package
  
  return(col[as.numeric(cut(x, breaks = breaks, include.lowest = TRUE))])
}

.scale_colours <- function (mat, col, breaks = NA) {
  # modified from pheatmap package
  
  mat = as.matrix(mat)
  return(matrix(.scale_vec_colours(as.vector(mat), col = col, breaks = breaks), nrow(mat), 
                ncol(mat), dimnames = list(rownames(mat), colnames(mat))))
}

.scale_rows = function(x){
  # modified from pheatmap package
  
  m = apply(x, 1, mean, na.rm = TRUE)
  s = apply(x, 1, sd, na.rm = TRUE)
  return((x - m) / s)
}

.scale_mat = function(mat, scale){
  # modified from pheatmap package
  
  if(!(scale %in% c("none", "row", "column"))){
    stop("scale argument shoud take values: 'none', 'row' or 'column'")
  }
  mat = switch(scale, none = mat, row = .scale_rows(mat), column = t(.scale_rows(t(mat))))
  return(mat)
}

.is.na2 = function(x){
  # modified from pheatmap package
  
  if(is.list(x) | length(x) > 1){
    return(FALSE)
  }
  if(length(x) == 0){
    return(TRUE)
  }
  
  return(is.na(x))
}

.identity2 = function(x, ...){
  # modified from pheatmap package
  
  return(x)
}

.pval.correct <- function(pval.matrix,maxrow){
  p=ncol(pval.matrix)
  pval.matrix.2x=cbind(pval.matrix[,p:1],pval.matrix[,p:1])
  corrected.pval.matrix=matrix(nrow=p,ncol=p)
  corrected.pval.matrix[p,]=rev(pval.matrix[p,])
  for(var in 1:p){
    pval_var=pval.matrix.2x[p,var]
    start=var
    end=var
    for(row in (p-1):maxrow){
      end=end+1
      pval_cone=pval.matrix.2x[row,start:end]
      pval_var=suppressWarnings(max(pval_var,pval_cone,na.rm=TRUE))
      corrected.pval.matrix[row,var]=pval_var
    }
  }
  corrected.pval.matrix=corrected.pval.matrix[,p:1]
  return(corrected.pval.matrix)
}
