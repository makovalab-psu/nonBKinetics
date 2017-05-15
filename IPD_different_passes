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


setwd("C:/Users/MarziaAngela/Documents/Marzia/Progetti/PacBio times/Monika")



# G-quadruplex (chr7 missing)
load('IPDlist_passes_Gquad.rda')

IPD_start_0=list(pass1=IPDlistPASS1_0,
                 pass2=IPDlistPASS2_1,
                 pass3=IPDlistPASS3_0,
                 pass4=IPDlistPASS4_1,
                 pass5=IPDlistPASS5_0)

IPD_start_1=list(pass1=IPDlistPASS1_1,
                 pass2=IPDlistPASS2_0,
                 pass3=IPDlistPASS3_1,
                 pass4=IPDlistPASS4_0,
                 pass5=IPDlistPASS5_1)




IPD_start_0_mean=lapply(IPD_start_0,
                        function(IPD) unlist(lapply(IPD,function(IPD) mean(log(IPD+0.01),na.rm=TRUE))))
IPD_start_1_mean=lapply(IPD_start_1,
                        function(IPD) unlist(lapply(IPD,function(IPD) mean(log(IPD+0.01),na.rm=TRUE))))

box=boxplot(IPD_start_0_mean,plot=FALSE)
box$stats[c(1,5),]=Reduce(cbind,lapply(IPD_start_0_mean,quantile,probs=c(0.05,0.95)))
bxp(box,outline=FALSE,main=paste0('Start from 0 strand (',length(IPD_start_0_mean$pass1),' molecules)'),ylab='Log mean IPD')

box=boxplot(IPD_start_0_mean[c(1,3,5)],plot=FALSE)
box$stats[c(1,5),]=Reduce(cbind,lapply(IPD_start_0_mean[c(1,3,5)],quantile,probs=c(0.05,0.95)))
bxp(box,outline=FALSE,main='Strand without G-quadruplex',ylab='Log mean IPD')
#matplot(Reduce(rbind,IPD_start_0_mean[c(1,3,5)])[,1:20],type='l',lty=1,add=TRUE)
box=boxplot(cbind(IPD_start_0_mean[[3]]-IPD_start_0_mean[[1]],
                  IPD_start_0_mean[[5]]-IPD_start_0_mean[[1]]),plot=FALSE)
box$stats[c(1,5),]=apply(cbind(IPD_start_0_mean[[3]]-IPD_start_0_mean[[1]],
                               IPD_start_0_mean[[5]]-IPD_start_0_mean[[1]]),2,quantile,probs=c(0.05,0.95))
bxp(box,outline=FALSE,ylim=c(-1.2,1.2),main='Strand without G-quadruplex',ylab='Log mean IPD',show.names=FALSE)
axis(1,1:2,c('pass 3-1','pass 5-1'))
abline(h=0,col='red')
#pval_3_1=signif(t.test(IPD_start_0_mean[[3]],IPD_start_0_mean[[1]],paired=TRUE)$p.val,2)
#pval_5_1=signif(t.test(IPD_start_0_mean[[5]],IPD_start_0_mean[[1]],paired=TRUE)$p.val,2)
pval_3_1=test_scalar(IPD_start_0_mean[[3]],IPD_start_0_mean[[1]],
                     statistics='quantile',probs=c(0.05,0.25,0.50,0.75,0.95),paired=TRUE,B=10000)$result$unadjusted_pval
pval_5_1=test_scalar(IPD_start_0_mean[[5]],IPD_start_0_mean[[1]],
                     statistics='quantile',probs=c(0.05,0.25,0.50,0.75,0.95),paired=TRUE,B=10000)$result$unadjusted_pval
text(1:2,-1.1,labels=c(pval_3_1,pval_5_1))

box=boxplot(IPD_start_0_mean[c(2,4)],plot=FALSE)
box$stats[c(1,5),]=Reduce(cbind,lapply(IPD_start_0_mean[c(2,4)],quantile,probs=c(0.05,0.95)))
bxp(box,outline=FALSE,main='Strand with G-quadruplex',ylab='Log mean IPD')
#matplot(Reduce(rbind,IPD_start_0_mean[c(2,4)])[,1:20],type='l',lty=1,add=TRUE)
box=boxplot(IPD_start_0_mean[[4]]-IPD_start_0_mean[[2]],plot=FALSE)
box$stats[c(1,5),]=quantile(IPD_start_0_mean[[4]]-IPD_start_0_mean[[2]],probs=c(0.05,0.95))
bxp(box,outline=FALSE,ylim=c(-1.2,1.2),main='Strand with G-quadruplex',ylab='Log mean IPD',show.names=FALSE)
axis(1,1,'pass 4-2')
abline(h=0,col='red')
#pval_4_2=signif(t.test(IPD_start_0_mean[[4]],IPD_start_0_mean[[2]],paired=TRUE)$p.val,2)
pval_4_2=test_scalar(IPD_start_0_mean[[4]],IPD_start_0_mean[[2]],
                     statistics='quantile',probs=c(0.05,0.25,0.50,0.75,0.95),paired=TRUE,B=10000)$result$unadjusted_pval
text(1,-1.1,labels=pval_4_2)

box=boxplot(IPD_start_0_mean[c(3,5)],plot=FALSE)
box$stats[c(1,5),]=Reduce(cbind,lapply(IPD_start_0_mean[c(3,5)],quantile,probs=c(0.05,0.95)))
bxp(box,outline=FALSE,main='Strand without G-quadruplex',ylab='Log mean IPD')
#matplot(Reduce(rbind,IPD_start_0_mean[c(3,5)])[,1:20],type='l',lty=1,add=TRUE)
box=boxplot(IPD_start_0_mean[[5]]-IPD_start_0_mean[[3]],plot=FALSE)
box$stats[c(1,5),]=quantile(IPD_start_0_mean[[5]]-IPD_start_0_mean[[3]],probs=c(0.05,0.95))
bxp(box,outline=FALSE,ylim=c(-1.2,1.2),main='Strand without G-quadruplex',ylab='Log mean IPD',show.names=FALSE)
axis(1,1,'pass 5-3')
abline(h=0,col='red')
#pval_5_3=signif(t.test(IPD_start_0_mean[[5]],IPD_start_0_mean[[3]],paired=TRUE)$p.val,2)
pval_5_3=test_scalar(IPD_start_0_mean[[5]],IPD_start_0_mean[[3]],
                     statistics='quantile',probs=c(0.05,0.25,0.50,0.75,0.95),paired=TRUE,B=10000)$result$unadjusted_pval
text(1,-1.1,labels=pval_5_3)



box=boxplot(IPD_start_1_mean,plot=FALSE)
box$stats[c(1,5),]=Reduce(cbind,lapply(IPD_start_1_mean,quantile,probs=c(0.05,0.95)))
bxp(box,outline=FALSE,main=paste0('Start from 1 strand (',length(IPD_start_1_mean$pass1),' molecules)'),ylab='Log mean IPD')

box=boxplot(IPD_start_1_mean[c(1,3,5)],plot=FALSE)
box$stats[c(1,5),]=Reduce(cbind,lapply(IPD_start_1_mean[c(1,3,5)],quantile,probs=c(0.05,0.95)))
bxp(box,outline=FALSE,main='Strand with G-quadruplex',ylab='Log mean IPD')
#matplot(Reduce(rbind,IPD_start_1_mean[c(1,3,5)])[,1:20],type='l',lty=1,add=TRUE)
box=boxplot(cbind(IPD_start_1_mean[[3]]-IPD_start_1_mean[[1]],
                  IPD_start_1_mean[[5]]-IPD_start_1_mean[[1]]),plot=FALSE)
box$stats[c(1,5),]=apply(cbind(IPD_start_1_mean[[3]]-IPD_start_1_mean[[1]],
                               IPD_start_1_mean[[5]]-IPD_start_1_mean[[1]]),2,quantile,probs=c(0.05,0.95))
bxp(box,outline=FALSE,ylim=c(-1.2,1.2),main='Strand with G-quadruplex',ylab='Log mean IPD',show.names=FALSE)
axis(1,1:2,c('pass 3-1','pass 5-1'))
abline(h=0,col='red')
#pval_3_1=signif(t.test(IPD_start_1_mean[[3]],IPD_start_1_mean[[1]],paired=TRUE)$p.val,2)
#pval_5_1=signif(t.test(IPD_start_1_mean[[5]],IPD_start_1_mean[[1]],paired=TRUE)$p.val,2)
pval_3_1=test_scalar(IPD_start_1_mean[[3]],IPD_start_1_mean[[1]],
                     statistics='quantile',probs=c(0.05,0.25,0.50,0.75,0.95),paired=TRUE,B=10000)$result$unadjusted_pval
pval_5_1=test_scalar(IPD_start_1_mean[[5]],IPD_start_1_mean[[1]],
                     statistics='quantile',probs=c(0.05,0.25,0.50,0.75,0.95),paired=TRUE,B=10000)$result$unadjusted_pval
text(1:2,-1.1,labels=c(pval_3_1,pval_5_1))

box=boxplot(IPD_start_1_mean[c(2,4)],plot=FALSE)
box$stats[c(1,5),]=Reduce(cbind,lapply(IPD_start_1_mean[c(2,4)],quantile,probs=c(0.05,0.95)))
bxp(box,outline=FALSE,main='Strand without G-quadruplex',ylab='Log mean IPD')
#matplot(Reduce(rbind,IPD_start_1_mean[c(2,4)])[,1:20],type='l',lty=1,add=TRUE)
box=boxplot(IPD_start_1_mean[[4]]-IPD_start_1_mean[[2]],plot=FALSE)
box$stats[c(1,5),]=quantile(IPD_start_1_mean[[4]]-IPD_start_1_mean[[2]],probs=c(0.05,0.95))
bxp(box,outline=FALSE,ylim=c(-1.2,1.2),main='Strand without G-quadruplex',ylab='Log mean IPD',show.names=FALSE)
axis(1,1,'pass 4-2')
abline(h=0,col='red')
#pval_4_2=signif(t.test(IPD_start_1_mean[[4]],IPD_start_1_mean[[2]],paired=TRUE)$p.val,2)
pval_4_2=test_scalar(IPD_start_1_mean[[4]],IPD_start_1_mean[[2]],
                     statistics='quantile',probs=c(0.05,0.25,0.50,0.75,0.95),paired=TRUE,B=10000)$result$unadjusted_pval
text(1,-1.1,labels=pval_4_2)

box=boxplot(IPD_start_1_mean[c(3,5)],plot=FALSE)
box$stats[c(1,5),]=Reduce(cbind,lapply(IPD_start_1_mean[c(3,5)],quantile,probs=c(0.05,0.95)))
bxp(box,outline=FALSE,main='Strand with G-quadruplex',ylab='Log mean IPD')
#matplot(Reduce(rbind,IPD_start_1_mean[c(3,5)])[,1:20],type='l',lty=1,add=TRUE)
box=boxplot(IPD_start_1_mean[[5]]-IPD_start_1_mean[[3]],plot=FALSE)
box$stats[c(1,5),]=quantile(IPD_start_1_mean[[5]]-IPD_start_1_mean[[3]],probs=c(0.05,0.95))
bxp(box,outline=FALSE,ylim=c(-1.2,1.2),main='Strand with G-quadruplex',ylab='Log mean IPD',show.names=FALSE)
axis(1,1,'pass 5-3')
abline(h=0,col='red')
#pval_5_3=signif(t.test(IPD_start_1_mean[[5]],IPD_start_1_mean[[3]],paired=TRUE)$p.val,2)
pval_5_3=test_scalar(IPD_start_1_mean[[5]],IPD_start_1_mean[[3]],
                     statistics='quantile',probs=c(0.05,0.25,0.50,0.75,0.95),paired=TRUE,B=10000)$result$unadjusted_pval
text(1,-1.1,labels=pval_5_3)










# Control
load('IPDlist_passes_control.rda')

IPD_start_0=list(pass1=IPDlistPASS1_0,
                 pass2=IPDlistPASS2_1,
                 pass3=IPDlistPASS3_0,
                 pass4=IPDlistPASS4_1,
                 pass5=IPDlistPASS5_0)

IPD_start_1=list(pass1=IPDlistPASS1_1,
                 pass2=IPDlistPASS2_0,
                 pass3=IPDlistPASS3_1,
                 pass4=IPDlistPASS4_0,
                 pass5=IPDlistPASS5_1)




IPD_start_0_mean=lapply(IPD_start_0,
                        function(IPD) unlist(lapply(IPD,function(IPD) mean(log(IPD+0.01),na.rm=TRUE))))
IPD_start_1_mean=lapply(IPD_start_1,
                        function(IPD) unlist(lapply(IPD,function(IPD) mean(log(IPD+0.01),na.rm=TRUE))))

box=boxplot(IPD_start_0_mean,plot=FALSE)
box$stats[c(1,5),]=Reduce(cbind,lapply(IPD_start_0_mean,quantile,probs=c(0.05,0.95)))
bxp(box,outline=FALSE,main=paste0('Start from 0 strand (',length(IPD_start_0_mean$pass1),' molecules)'),ylab='Log mean IPD')
#matplot(Reduce(rbind,IPD_start_0_mean)[,1:20],type='l',lty=1,add=TRUE)
matplot(Reduce(rbind,IPD_start_0_mean)[,1:20],type='l',lty=1,
        main=paste0('Start from 0 strand (',length(IPD_start_0_mean$pass1),' molecules)'),ylab='Log mean IPD')
box=boxplot(cbind(IPD_start_0_mean[[2]]-IPD_start_0_mean[[1]],IPD_start_0_mean[[3]]-IPD_start_0_mean[[1]],
                  IPD_start_0_mean[[4]]-IPD_start_0_mean[[1]],IPD_start_0_mean[[5]]-IPD_start_0_mean[[1]]),plot=FALSE)
box$stats[c(1,5),]=apply(cbind(IPD_start_0_mean[[2]]-IPD_start_0_mean[[1]],IPD_start_0_mean[[3]]-IPD_start_0_mean[[1]],
                               IPD_start_0_mean[[4]]-IPD_start_0_mean[[1]],IPD_start_0_mean[[5]]-IPD_start_0_mean[[1]]),2,quantile,probs=c(0.05,0.95))
bxp(box,outline=FALSE,ylim=c(-1.2,1.2),main=paste0('Start from 0 strand (',length(IPD_start_0_mean$pass1),' molecules)'),ylab='Log mean IPD',show.names=FALSE)
axis(1,1:4,c('pass 2-1','pass 3-1','pass 4-1','pass 5-1'))
abline(h=0,col='red')
#pval_2_1=signif(t.test(IPD_start_0_mean[[2]],IPD_start_0_mean[[1]],paired=TRUE)$p.val,2)
#pval_3_1=signif(t.test(IPD_start_0_mean[[3]],IPD_start_0_mean[[1]],paired=TRUE)$p.val,2)
#pval_4_1=signif(t.test(IPD_start_0_mean[[4]],IPD_start_0_mean[[1]],paired=TRUE)$p.val,2)
#pval_5_1=signif(t.test(IPD_start_0_mean[[5]],IPD_start_0_mean[[1]],paired=TRUE)$p.val,2)
pval_2_1=test_scalar(IPD_start_0_mean[[2]],IPD_start_0_mean[[1]],
                     statistics='quantile',probs=c(0.05,0.25,0.50,0.75,0.95),paired=TRUE,B=10000)$result$unadjusted_pval
pval_3_1=test_scalar(IPD_start_0_mean[[3]],IPD_start_0_mean[[1]],
                     statistics='quantile',probs=c(0.05,0.25,0.50,0.75,0.95),paired=TRUE,B=10000)$result$unadjusted_pval
pval_4_1=test_scalar(IPD_start_0_mean[[4]],IPD_start_0_mean[[1]],
                     statistics='quantile',probs=c(0.05,0.25,0.50,0.75,0.95),paired=TRUE,B=10000)$result$unadjusted_pval
pval_5_1=test_scalar(IPD_start_0_mean[[5]],IPD_start_0_mean[[1]],
                     statistics='quantile',probs=c(0.05,0.25,0.50,0.75,0.95),paired=TRUE,B=10000)$result$unadjusted_pval
text(1:4,-1.1,labels=c(pval_2_1,pval_3_1,pval_4_1,pval_5_1))




box=boxplot(IPD_start_1_mean,plot=FALSE)
box$stats[c(1,5),]=Reduce(cbind,lapply(IPD_start_1_mean,quantile,probs=c(0.05,0.95)))
bxp(box,outline=FALSE,main=paste0('Start from 1 strand (',length(IPD_start_1_mean$pass1),' molecules)'),ylab='Log mean IPD')
#matplot(Reduce(rbind,IPD_start_1_mean)[,1:20],type='l',lty=1,add=TRUE)
matplot(Reduce(rbind,IPD_start_1_mean)[,1:20],type='l',lty=1,
        main=paste0('Start from 1 strand (',length(IPD_start_1_mean$pass1),' molecules)'),ylab='Log mean IPD')

box=boxplot(cbind(IPD_start_1_mean[[2]]-IPD_start_1_mean[[1]],IPD_start_1_mean[[3]]-IPD_start_1_mean[[1]],
                  IPD_start_1_mean[[4]]-IPD_start_1_mean[[1]],IPD_start_1_mean[[5]]-IPD_start_1_mean[[1]]),plot=FALSE)
box$stats[c(1,5),]=apply(cbind(IPD_start_1_mean[[2]]-IPD_start_1_mean[[1]],IPD_start_1_mean[[3]]-IPD_start_1_mean[[1]],
                               IPD_start_1_mean[[4]]-IPD_start_1_mean[[1]],IPD_start_1_mean[[5]]-IPD_start_1_mean[[1]]),2,quantile,probs=c(0.05,0.95))
bxp(box,outline=FALSE,ylim=c(-1.2,1.2),main=paste0('Start from 1 strand (',length(IPD_start_1_mean$pass1),' molecules)'),ylab='Log mean IPD',show.names=FALSE)
axis(1,1:4,c('pass 2-1','pass 3-1','pass 4-1','pass 5-1'))
abline(h=0,col='red')
#pval_2_1=signif(t.test(IPD_start_1_mean[[2]],IPD_start_1_mean[[1]],paired=TRUE)$p.val,2)
#pval_3_1=signif(t.test(IPD_start_1_mean[[3]],IPD_start_1_mean[[1]],paired=TRUE)$p.val,2)
#pval_4_1=signif(t.test(IPD_start_1_mean[[4]],IPD_start_1_mean[[1]],paired=TRUE)$p.val,2)
#pval_5_1=signif(t.test(IPD_start_1_mean[[5]],IPD_start_1_mean[[1]],paired=TRUE)$p.val,2)
pval_2_1=test_scalar(IPD_start_1_mean[[2]],IPD_start_1_mean[[1]],
                     statistics='quantile',probs=c(0.05,0.25,0.50,0.75,0.95),paired=TRUE,B=10000)$result$unadjusted_pval
pval_3_1=test_scalar(IPD_start_1_mean[[3]],IPD_start_1_mean[[1]],
                     statistics='quantile',probs=c(0.05,0.25,0.50,0.75,0.95),paired=TRUE,B=10000)$result$unadjusted_pval
pval_4_1=test_scalar(IPD_start_1_mean[[4]],IPD_start_1_mean[[1]],
                     statistics='quantile',probs=c(0.05,0.25,0.50,0.75,0.95),paired=TRUE,B=10000)$result$unadjusted_pval
pval_5_1=test_scalar(IPD_start_1_mean[[5]],IPD_start_1_mean[[1]],
                     statistics='quantile',probs=c(0.05,0.25,0.50,0.75,0.95),paired=TRUE,B=10000)$result$unadjusted_pval
text(1:4,-1.1,labels=c(pval_2_1,pval_3_1,pval_4_1,pval_5_1))




IPD_start_01_mean=mapply(function(IPD_0,IPD_1) unlist(c(lapply(IPD_0,function(IPD) mean(log(IPD+0.01),na.rm=TRUE)),
                                                        lapply(IPD_1,function(IPD) mean(log(IPD+0.01),na.rm=TRUE)))),
                         IPD_start_0,IPD_start_1,SIMPLIFY=FALSE)

box=boxplot(IPD_start_01_mean,plot=FALSE)
box$stats[c(1,5),]=Reduce(cbind,lapply(IPD_start_01_mean,quantile,probs=c(0.05,0.95)))
bxp(box,outline=FALSE,main=paste0('Start from both strands (',length(IPD_start_01_mean$pass1),' molecules)'),ylab='Log mean IPD')
#matplot(Reduce(rbind,IPD_start_01_mean)[,1:20],type='l',lty=1,add=TRUE)
matplot(Reduce(rbind,IPD_start_01_mean)[,sample(length(IPD_start_01_mean$pass1),20)],type='l',lty=1,
        main=paste0('Start from both strands (',length(IPD_start_01_mean$pass1),' molecules)'),ylab='Log mean IPD')

box=boxplot(cbind(IPD_start_01_mean[[2]]-IPD_start_01_mean[[1]],IPD_start_01_mean[[3]]-IPD_start_01_mean[[1]],
                  IPD_start_01_mean[[4]]-IPD_start_01_mean[[1]],IPD_start_01_mean[[5]]-IPD_start_01_mean[[1]]),plot=FALSE)
box$stats[c(1,5),]=apply(cbind(IPD_start_01_mean[[2]]-IPD_start_01_mean[[1]],IPD_start_01_mean[[3]]-IPD_start_01_mean[[1]],
                               IPD_start_01_mean[[4]]-IPD_start_01_mean[[1]],IPD_start_01_mean[[5]]-IPD_start_01_mean[[1]]),2,quantile,probs=c(0.05,0.95))
bxp(box,outline=FALSE,ylim=c(-1.2,1.2),main=paste0('Start from both strands (',length(IPD_start_01_mean$pass1),' molecules)'),ylab='Log mean IPD',show.names=FALSE)
axis(1,1:4,c('pass 2-1','pass 3-1','pass 4-1','pass 5-1'))
abline(h=0,col='red')
#pval_2_1=signif(t.test(IPD_start_01_mean[[2]],IPD_start_01_mean[[1]],paired=TRUE)$p.val,2)
#pval_3_1=signif(t.test(IPD_start_01_mean[[3]],IPD_start_01_mean[[1]],paired=TRUE)$p.val,2)
#pval_4_1=signif(t.test(IPD_start_01_mean[[4]],IPD_start_01_mean[[1]],paired=TRUE)$p.val,2)
#pval_5_1=signif(t.test(IPD_start_01_mean[[5]],IPD_start_01_mean[[1]],paired=TRUE)$p.val,2)
pval_2_1=test_scalar(IPD_start_01_mean[[2]],IPD_start_01_mean[[1]],
                     statistics='quantile',probs=c(0.05,0.25,0.50,0.75,0.95),paired=TRUE,B=10000)$result$unadjusted_pval
pval_3_1=test_scalar(IPD_start_01_mean[[3]],IPD_start_01_mean[[1]],
                     statistics='quantile',probs=c(0.05,0.25,0.50,0.75,0.95),paired=TRUE,B=10000)$result$unadjusted_pval
pval_4_1=test_scalar(IPD_start_01_mean[[4]],IPD_start_01_mean[[1]],
                     statistics='quantile',probs=c(0.05,0.25,0.50,0.75,0.95),paired=TRUE,B=10000)$result$unadjusted_pval
pval_5_1=test_scalar(IPD_start_01_mean[[5]],IPD_start_01_mean[[1]],
                     statistics='quantile',probs=c(0.05,0.25,0.50,0.75,0.95),paired=TRUE,B=10000)$result$unadjusted_pval
text(1:4,-1.1,labels=c(pval_2_1,pval_3_1,pval_4_1,pval_5_1))




box=boxplot(cbind(IPD_start_01_mean[[3]]-IPD_start_01_mean[[2]],
                  IPD_start_01_mean[[4]]-IPD_start_01_mean[[2]],IPD_start_01_mean[[5]]-IPD_start_01_mean[[2]]),plot=FALSE)
box$stats[c(1,5),]=apply(cbind(IPD_start_01_mean[[3]]-IPD_start_01_mean[[2]],
                               IPD_start_01_mean[[4]]-IPD_start_01_mean[[2]],IPD_start_01_mean[[5]]-IPD_start_01_mean[[2]]),2,quantile,probs=c(0.05,0.95))
bxp(box,outline=FALSE,ylim=c(-1.2,1.2),main=paste0('Start from both strands (',length(IPD_start_01_mean$pass1),' molecules)'),ylab='Log mean IPD',show.names=FALSE)
axis(1,1:3,c('pass 3-2','pass 4-2','pass 5-2'))
abline(h=0,col='red')
#pval_3_2=signif(t.test(IPD_start_01_mean[[3]],IPD_start_01_mean[[2]],paired=TRUE)$p.val,2)
#pval_4_2=signif(t.test(IPD_start_01_mean[[4]],IPD_start_01_mean[[2]],paired=TRUE)$p.val,2)
#pval_5_2=signif(t.test(IPD_start_01_mean[[5]],IPD_start_01_mean[[2]],paired=TRUE)$p.val,2)
pval_3_2=test_scalar(IPD_start_01_mean[[3]],IPD_start_01_mean[[2]],
                     statistics='quantile',probs=c(0.05,0.25,0.50,0.75,0.95),paired=TRUE,B=10000)$result$unadjusted_pval
pval_4_2=test_scalar(IPD_start_01_mean[[4]],IPD_start_01_mean[[2]],
                     statistics='quantile',probs=c(0.05,0.25,0.50,0.75,0.95),paired=TRUE,B=10000)$result$unadjusted_pval
pval_5_2=test_scalar(IPD_start_01_mean[[5]],IPD_start_01_mean[[2]],
                     statistics='quantile',probs=c(0.05,0.25,0.50,0.75,0.95),paired=TRUE,B=10000)$result$unadjusted_pval
text(1:3,-1.1,labels=c(pval_3_2,pval_4_2,pval_5_2))




box=boxplot(cbind(IPD_start_01_mean[[4]]-IPD_start_01_mean[[3]],IPD_start_01_mean[[5]]-IPD_start_01_mean[[3]],
                  IPD_start_01_mean[[5]]-IPD_start_01_mean[[4]]),plot=FALSE)
box$stats[c(1,5),]=apply(cbind(IPD_start_01_mean[[4]]-IPD_start_01_mean[[3]],IPD_start_01_mean[[5]]-IPD_start_01_mean[[3]],
                               IPD_start_01_mean[[5]]-IPD_start_01_mean[[4]]),2,quantile,probs=c(0.05,0.95))
bxp(box,outline=FALSE,ylim=c(-1.2,1.2),main=paste0('Start from both strands (',length(IPD_start_01_mean$pass1),' molecules)'),ylab='Log mean IPD',show.names=FALSE)
axis(1,1:3,c('pass 4-3','pass 5-3','pass 5-4'))
abline(h=0,col='red')
#pval_4_3=signif(t.test(IPD_start_01_mean[[4]],IPD_start_01_mean[[3]],paired=TRUE)$p.val,2)
#pval_5_3=signif(t.test(IPD_start_01_mean[[5]],IPD_start_01_mean[[3]],paired=TRUE)$p.val,2)
#pval_5_4=signif(t.test(IPD_start_01_mean[[5]],IPD_start_01_mean[[4]],paired=TRUE)$p.val,2)
pval_4_3=test_scalar(IPD_start_01_mean[[4]],IPD_start_01_mean[[3]],
                     statistics='quantile',probs=c(0.05,0.25,0.50,0.75,0.95),paired=TRUE,B=10000)$result$unadjusted_pval
pval_5_3=test_scalar(IPD_start_01_mean[[5]],IPD_start_01_mean[[3]],
                     statistics='quantile',probs=c(0.05,0.25,0.50,0.75,0.95),paired=TRUE,B=10000)$result$unadjusted_pval
pval_5_4=test_scalar(IPD_start_01_mean[[5]],IPD_start_01_mean[[4]],
                     statistics='quantile',probs=c(0.05,0.25,0.50,0.75,0.95),paired=TRUE,B=10000)$result$unadjusted_pval
text(1:3,-1.1,labels=c(pval_4_3,pval_5_3,pval_5_4))
