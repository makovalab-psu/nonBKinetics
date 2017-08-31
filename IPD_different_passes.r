# Scalar version of the test used in IWTomics package
test_scalar <- function(data1,data2,mu=0,statistics='mean',probs=0.5,paired=FALSE,B=1000){
  # data1 and data2 matrices with 1 row
  
  if(statistics=='median'){
    statistics='quantile'
    probs=0.5
  }
  n1=length(data1)
  n2=length(data2)
  n=n1+n2
  data1=data1-mu
  if(paired){
    exact=(B>=(2^n1))
  }else{
    exact=(B>=choose(n,n1))
  }
  
  p=1
  result=list(test='2pop',mu=mu)
  data=c(data1,data2)
  
  # Univariate permutations
  message('    Point-wise tests...')
  if(statistics=='mean'){
    T0_plot=mean(data1,na.rm=TRUE)-mean(data2,na.rm=TRUE)
    T0=(T0_plot)^2
  }
  if(statistics=='quantile'){
    T0_plot=quantile(data1,probs=probs,na.rm=TRUE)-quantile(data2,probs=probs,na.rm=TRUE)
    T0=(T0_plot)^2
    T0_plot=sum(T0_plot)
    T0=sum(T0)
  }
  if(statistics=='variance'){
    T0_plot=var(data1,na.rm=TRUE)/var(data2,na.rm=TRUE)
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
                                                   data_perm[c(change,n1+change)]=data_perm[c(n1+change,change)]
                                                   if(statistics=='mean')
                                                     return((mean(data_perm[seq.int(n1)],na.rm=TRUE)-mean(data_perm[n1+seq.int(n2)],na.rm=TRUE))^2)
                                                   if(statistics=='quantile'){
                                                     T_perm=(quantile(data_perm[seq.int(n1)],probs=probs,na.rm=TRUE)-quantile(data_perm[n1+seq.int(n2)],probs=probs,na.rm=TRUE))^2
                                                     return(sum(T_perm))
                                                   }
                                                   if(statistics=='variance')
                                                     return((var(data_perm[seq.int(n1)],na.rm=TRUE)/var(data_perm[n1+seq.int(n2)],na.rm=TRUE))^2)
                                                 })
                                    return(T_perm)
                                  }))
    }else{
      first_group=combn(n,n1)
      T_perm=apply(first_group,2,
                   function(group){
                     data_perm=data[c(group,setdiff(seq_len(n1+n2),group))]
                     if(statistics=='mean')
                       return((mean(data_perm[seq.int(n1)],na.rm=TRUE)-mean(data_perm[n1+seq.int(n2)],na.rm=TRUE))^2)
                     if(statistics=='quantile'){
                       T_perm=(quantile(data_perm[seq.int(n1)],probs=probs,na.rm=TRUE)-quantile(data_perm[n1+seq.int(n2)],probs=probs,na.rm=TRUE))^2
                       return(sum(T_perm))
                     }
                     if(statistics=='variance')
                       return((var(data_perm[seq.int(n1)],na.rm=TRUE)/var(data_perm[n1+seq.int(n2)],na.rm=TRUE))^2)
                   })
    }
  }else{
    T_perm=do.call(cbind,lapply(seq.int(B-1),
                                function(perm){
                                  if(paired){
                                    couple.perm=rbinom(n1,1,0.5)
                                    data_perm=data[c(n1*couple.perm,-n1*couple.perm)+seq.int(2*n1)]
                                  }else{
                                    permutation=sample(n,n1)
                                    data_perm=data[c(permutation,setdiff(seq_len(n),permutation))]
                                  }
                                  if(statistics=='mean')
                                    return((mean(data_perm[seq.int(n1)],na.rm=TRUE)-mean(data_perm[n1+seq.int(n2)],na.rm=TRUE))^2)
                                  if(statistics=='quantile'){
                                    T_perm=(quantile(data_perm[seq.int(n1)],probs=probs,na.rm=TRUE)-quantile(data_perm[n1+seq.int(n2)],probs=probs,na.rm=TRUE))^2
                                    return(sum(T_perm))
                                  }
                                  if(statistics=='variance')
                                    return((var(data_perm[seq.int(n1)],na.rm=TRUE)/var(data_perm[n1+seq.int(n2)],na.rm=TRUE))^2)
                                }))
    T_perm=c(T_perm,T0)
  }
  
  # Not fully computable p-values (some permutations do not produce any test statistics because of the NAs)
  #if(statistics=='variance'){
  #  no.pval=rowSums(is.na(data))>=(min(n1,n2)-1)
  #}else{
  #  no.pval=rowSums(is.na(data))>=min(n1,n2)
  #}
  #T_perm[no.pval,]=NaN # do not compute any p-value when it is not fully computable
  # do not compute any p-value if NaN is in T_perm
  #if(statistics!='variance'){
  #  pval=rowSums(T_perm>=T0)/B
  #}else{
  #  pval=pmin(2*rowSums(T_perm>=T0)/B,2*rowSums(T_perm<=T0)/B)
  #}
  # compute p-value omitting NaN
  if(statistics!='variance'){
    pval=sum(T_perm>=T0,na.rm=TRUE)/sum(!is.nan(T_perm))
  }else{
    pval=pmin(2*sum(T_perm>=T0,na.rm=TRUE)/sum(!is.nan(T_perm)),2*sum(T_perm<=T0,na.rm=TRUE)/sum(!is.nan(T_perm)))
  }
  
  result$T0_plot=T0_plot
  result$unadjusted_pval=pval
  result$exact=exact
  class(result)='ITWomics.2pop'
  return(list(result=result))
}


################################
# G-quadruplex vs Control plus #
################################


# Load G-quadruplex data
load('IPDlist_passes.GQuadPlusFeatureOnly.mf.gff.intersect.IPDs.txt.rda')

# Molecules that start from strand 0
IPD_start_0=list(pass1=IPDlistPASS1_0,
                 pass2=IPDlistPASS2_1,
                 pass3=IPDlistPASS3_0,
                 pass4=IPDlistPASS4_1)

# Molecules that start from strand 1
IPD_start_1=list(pass1=IPDlistPASS1_1,
                 pass2=IPDlistPASS2_0,
                 pass3=IPDlistPASS3_1,
                 pass4=IPDlistPASS4_0)

# Compute log meanIPD
IPD_start_0_mean_Gquad=lapply(IPD_start_0,
                              function(IPD) unlist(lapply(IPD,function(IPD) log(mean(IPD,na.rm=TRUE)+0.01))))
IPD_start_1_mean_Gquad=lapply(IPD_start_1,
                              function(IPD) unlist(lapply(IPD,function(IPD) log(mean(IPD,na.rm=TRUE)+0.01))))



# Load matching motif-free regions
load('IPDlist_passes.GQuadPlusFeatureOnly.mfEmptyTmp.mf.gff.intersect.IPDs.txt.rda')

# Molecules that start from strand 0
IPD_start_0=list(pass1=IPDlistPASS1_0,
                 pass2=IPDlistPASS2_1,
                 pass3=IPDlistPASS3_0,
                 pass4=IPDlistPASS4_1)

# Molecules that start from strand 1
IPD_start_1=list(pass1=IPDlistPASS1_1,
                 pass2=IPDlistPASS2_0,
                 pass3=IPDlistPASS3_1,
                 pass4=IPDlistPASS4_0)

# Compute log meanIPD
IPD_start_0_mean_Control=lapply(IPD_start_0,
                                function(IPD) unlist(lapply(IPD,function(IPD) log(mean(IPD,na.rm=TRUE)+0.01))))
IPD_start_1_mean_Control=lapply(IPD_start_1,
                                function(IPD) unlist(lapply(IPD,function(IPD) log(mean(IPD,na.rm=TRUE)+0.01))))


IPD_start_0_mean=unlist(mapply(function(Gquad,Control) list(Gquad=Gquad,Control=Control),IPD_start_0_mean_Gquad,IPD_start_0_mean_Control,SIMPLIFY=FALSE),recursive=FALSE)
IPD_start_1_mean=unlist(mapply(function(Gquad,Control) list(Gquad=Gquad,Control=Control),IPD_start_1_mean_Gquad,IPD_start_1_mean_Control,SIMPLIFY=FALSE),recursive=FALSE)



# Test motifs vs motif-free regions
pval_1=test_scalar(IPD_start_0_mean[[1]],IPD_start_0_mean[[2]],
                   statistics='quantile',probs=c(0.05,0.25,0.50,0.75,0.95),paired=FALSE,B=10000)$result$unadjusted_pval
pval_2=test_scalar(IPD_start_0_mean[[3]],IPD_start_0_mean[[4]],
                   statistics='quantile',probs=c(0.05,0.25,0.50,0.75,0.95),paired=FALSE,B=10000)$result$unadjusted_pval
pval_3=test_scalar(IPD_start_0_mean[[5]],IPD_start_0_mean[[6]],
                   statistics='quantile',probs=c(0.05,0.25,0.50,0.75,0.95),paired=FALSE,B=10000)$result$unadjusted_pval
pval_4=test_scalar(IPD_start_0_mean[[7]],IPD_start_0_mean[[8]],
                   statistics='quantile',probs=c(0.05,0.25,0.50,0.75,0.95),paired=FALSE,B=10000)$result$unadjusted_pval
box=boxplot(IPD_start_0_mean,plot=FALSE)
box$stats[c(1,5),]=Reduce(cbind,lapply(IPD_start_0_mean,quantile,probs=c(0.05,0.95),na.rm=TRUE))
par(mar=c(7.5,4,4,2)+0.1)
bxp(box,at=c(1:2,4:5,7:8,10:11),outline=FALSE,ylim=c(-3,1),show.names=FALSE,
    main='Start from G4- as template',ylab='Log mean IPD',las=3)
text(c(1.5,4.5,7.5,10.5),-2.5,labels=c(pval_1,pval_2,pval_3,pval_4),col=c('black','red','black','red'))
axis(1,at=c(1:2,4:5,7:8,10:11),labels=paste0(c('G4+','Motif-free'),' pass',rep(1:4,each=2)),las=2)

# Test motifs vs motif-free regions
pval_1=test_scalar(IPD_start_1_mean[[1]],IPD_start_1_mean[[2]],
                   statistics='quantile',probs=c(0.05,0.25,0.50,0.75,0.95),paired=FALSE,B=10000)$result$unadjusted_pval
pval_2=test_scalar(IPD_start_1_mean[[3]],IPD_start_1_mean[[4]],
                   statistics='quantile',probs=c(0.05,0.25,0.50,0.75,0.95),paired=FALSE,B=10000)$result$unadjusted_pval
pval_3=test_scalar(IPD_start_1_mean[[5]],IPD_start_1_mean[[6]],
                   statistics='quantile',probs=c(0.05,0.25,0.50,0.75,0.95),paired=FALSE,B=10000)$result$unadjusted_pval
pval_4=test_scalar(IPD_start_1_mean[[7]],IPD_start_1_mean[[8]],
                   statistics='quantile',probs=c(0.05,0.25,0.50,0.75,0.95),paired=FALSE,B=10000)$result$unadjusted_pval
box=boxplot(IPD_start_1_mean,plot=FALSE)
box$stats[c(1,5),]=Reduce(cbind,lapply(IPD_start_1_mean,quantile,probs=c(0.05,0.95),na.rm=TRUE))
par(mar=c(7.5,4,4,2)+0.1)
bxp(box,at=c(1:2,4:5,7:8,10:11),outline=FALSE,ylim=c(-3,1),show.names=FALSE,
    main='Start from G4+ as template',ylab='Log mean IPD',las=3)
text(c(1.5,4.5,7.5,10.5),-2.5,labels=c(pval_1,pval_2,pval_3,pval_4),col=c('red','black','red','black'))
axis(1,at=c(1:2,4:5,7:8,10:11),labels=paste0(c('G4+','Motif-free'),' pass',rep(1:4,each=2)),las=2)










################################
# G-quadruplex vs Control plus #
################################


# Load G-quadruplex data
load('IPDlist_passes.GQuadMinusFeatureOnly.mf.gff.intersect.IPDs.txt.rda')

# Molecules that start from strand 0
IPD_start_0=list(pass1=IPDlistPASS1_0,
                 pass2=IPDlistPASS2_1,
                 pass3=IPDlistPASS3_0,
                 pass4=IPDlistPASS4_1)

# Molecules that start from strand 1
IPD_start_1=list(pass1=IPDlistPASS1_1,
                 pass2=IPDlistPASS2_0,
                 pass3=IPDlistPASS3_1,
                 pass4=IPDlistPASS4_0)

# Compute log meanIPD
IPD_start_0_mean_Gquad=lapply(IPD_start_0,
                              function(IPD) unlist(lapply(IPD,function(IPD) log(mean(IPD,na.rm=TRUE)+0.01))))
IPD_start_1_mean_Gquad=lapply(IPD_start_1,
                              function(IPD) unlist(lapply(IPD,function(IPD) log(mean(IPD,na.rm=TRUE)+0.01))))



# Load matching motif-free regions
load('IPDlist_passes.GQuadMinusFeatureOnly.mfEmptyTmp.mf.gff.intersect.IPDs.txt.rda')

# Molecules that start from strand 0
IPD_start_0=list(pass1=IPDlistPASS1_0,
                 pass2=IPDlistPASS2_1,
                 pass3=IPDlistPASS3_0,
                 pass4=IPDlistPASS4_1)

# Molecules that start from strand 1
IPD_start_1=list(pass1=IPDlistPASS1_1,
                 pass2=IPDlistPASS2_0,
                 pass3=IPDlistPASS3_1,
                 pass4=IPDlistPASS4_0)

# Compute log meanIPD
IPD_start_0_mean_Control=lapply(IPD_start_0,
                                function(IPD) unlist(lapply(IPD,function(IPD) log(mean(IPD,na.rm=TRUE)+0.01))))
IPD_start_1_mean_Control=lapply(IPD_start_1,
                                function(IPD) unlist(lapply(IPD,function(IPD) log(mean(IPD,na.rm=TRUE)+0.01))))

IPD_start_0_mean=unlist(mapply(function(Gquad,Control) list(Gquad=Gquad,Control=Control),IPD_start_0_mean_Gquad,IPD_start_0_mean_Control,SIMPLIFY=FALSE),recursive=FALSE)
IPD_start_1_mean=unlist(mapply(function(Gquad,Control) list(Gquad=Gquad,Control=Control),IPD_start_1_mean_Gquad,IPD_start_1_mean_Control,SIMPLIFY=FALSE),recursive=FALSE)



# Test motifs vs motif-free regions
pval_1=test_scalar(IPD_start_0_mean[[1]],IPD_start_0_mean[[2]],
                   statistics='quantile',probs=c(0.05,0.25,0.50,0.75,0.95),paired=FALSE,B=10000)$result$unadjusted_pval
pval_2=test_scalar(IPD_start_0_mean[[3]],IPD_start_0_mean[[4]],
                   statistics='quantile',probs=c(0.05,0.25,0.50,0.75,0.95),paired=FALSE,B=10000)$result$unadjusted_pval
pval_3=test_scalar(IPD_start_0_mean[[5]],IPD_start_0_mean[[6]],
                   statistics='quantile',probs=c(0.05,0.25,0.50,0.75,0.95),paired=FALSE,B=10000)$result$unadjusted_pval
pval_4=test_scalar(IPD_start_0_mean[[7]],IPD_start_0_mean[[8]],
                   statistics='quantile',probs=c(0.05,0.25,0.50,0.75,0.95),paired=FALSE,B=10000)$result$unadjusted_pval
box=boxplot(IPD_start_0_mean,plot=FALSE)
box$stats[c(1,5),]=Reduce(cbind,lapply(IPD_start_0_mean,quantile,probs=c(0.05,0.95),na.rm=TRUE))
par(mar=c(7,4,4,2)+0.1)
bxp(box,at=c(1:2,4:5,7:8,10:11),outline=FALSE,ylim=c(-3,1.8),
    main='G-quadruplex rev vs Control - Start from 0 strand',ylab='Log mean IPD',las=3)
text(c(1.5,4.5,7.5,10.5),-2.5,labels=c(pval_1,pval_2,pval_3,pval_4),col=c('red','black','red','black'))

# Test motifs vs motif-free regions
pval_1=test_scalar(IPD_start_1_mean[[1]],IPD_start_1_mean[[2]],
                   statistics='quantile',probs=c(0.05,0.25,0.50,0.75,0.95),paired=FALSE,B=10000)$result$unadjusted_pval
pval_2=test_scalar(IPD_start_1_mean[[3]],IPD_start_1_mean[[4]],
                   statistics='quantile',probs=c(0.05,0.25,0.50,0.75,0.95),paired=FALSE,B=10000)$result$unadjusted_pval
pval_3=test_scalar(IPD_start_1_mean[[5]],IPD_start_1_mean[[6]],
                   statistics='quantile',probs=c(0.05,0.25,0.50,0.75,0.95),paired=FALSE,B=10000)$result$unadjusted_pval
pval_4=test_scalar(IPD_start_1_mean[[7]],IPD_start_1_mean[[8]],
                   statistics='quantile',probs=c(0.05,0.25,0.50,0.75,0.95),paired=FALSE,B=10000)$result$unadjusted_pval
box=boxplot(IPD_start_1_mean,plot=FALSE)
box$stats[c(1,5),]=Reduce(cbind,lapply(IPD_start_1_mean,quantile,probs=c(0.05,0.95),na.rm=TRUE))
par(mar=c(7,4,4,2)+0.1)
bxp(box,at=c(1:2,4:5,7:8,10:11),outline=FALSE,ylim=c(-3,1.8),
    main='G-quadruplex rev vs Control - Start from 1 strand',ylab='Log mean IPD',las=3)
text(c(1.5,4.5,7.5,10.5),-2.5,labels=c(pval_1,pval_2,pval_3,pval_4),col=c('black','red','black','red'))

