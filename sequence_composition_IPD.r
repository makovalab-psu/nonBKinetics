require(IWTomics)
source('IWTomics_modified_functions.r')
library(compositions)



############################
############################
###                      ###
###  WINDOWS WITHOUT NA  ###
###                      ###
############################
############################

load('IPD_forward_NA0.RData')
load('IPD_reverse_NA0.RData')

meanIPD_forward=lapply(regionsFeatures_forward@features$IPD_Forward,function(feature) colMeans(feature,na.rm=TRUE))
meanIPD_reverse=lapply(regionsFeatures_reverse@features$IPD_Reverse,function(feature) colMeans(feature,na.rm=TRUE))

meanIPD_composition_forward=mapply(function(meanIPD,region){
                                     cbind(meanIPD,as.data.frame(mcols(region)[,3:6]))
                                   },meanIPD_forward,regions(regionsFeatures_forward),SIMPLIFY=FALSE)
meanIPD_composition_reverse=mapply(function(meanIPD,region){
                                     cbind(meanIPD,as.data.frame(mcols(region)[,3:6]))
                                   },meanIPD_reverse,regions(regionsFeatures_reverse),SIMPLIFY=FALSE)
meanIPD_composition=mapply(function(forward,reverse) rbind(forward,reverse),
                           meanIPD_composition_forward[intersect(names(meanIPD_composition_forward),names(meanIPD_composition_reverse))],
                           meanIPD_composition_reverse[intersect(names(meanIPD_composition_forward),names(meanIPD_composition_reverse))],
                           SIMPLIFY=FALSE)

save(meanIPD_composition_forward,meanIPD_composition_reverse,meanIPD_composition,file='meanIPD_composition.RData')


load('meanIPD_composition.RData')
# scatterplot controls
pdf('meanIPD_composition_control.pdf',8,8)
par(mfrow=c(2,2))
ymax=max(meanIPD_composition_forward$Control$meanIPD,meanIPD_composition_reverse$Control$meanIPD)
# forward
smoothScatter(meanIPD_composition_forward$Control$A,meanIPD_composition_forward$Control$meanIPD,
              nrpoints=0,colramp=colorRampPalette(c("white",'red')),xlim=c(0,100),ylim=c(0,ymax),
              main='Mean IPD vs A',ylab='Mean IPD plus strand',xlab='A')
legend('topright',legend=paste0('r = ',round(cor(meanIPD_composition_forward$Control$A,meanIPD_composition_forward$Control$meanIPD),4)),bty='n')
smoothScatter(meanIPD_composition_forward$Control$T,meanIPD_composition_forward$Control$meanIPD,
              nrpoints=0,colramp=colorRampPalette(c("white",'red')),xlim=c(0,100),ylim=c(0,ymax),
              main='Mean IPD vs T',ylab='Mean IPD plus strand',xlab='T')
legend('topright',legend=paste0('r = ',round(cor(meanIPD_composition_forward$Control$T,meanIPD_composition_forward$Control$meanIPD),4)),bty='n')
smoothScatter(meanIPD_composition_forward$Control$G,meanIPD_composition_forward$Control$meanIPD,
              nrpoints=0,colramp=colorRampPalette(c("white",'red')),xlim=c(0,100),ylim=c(0,ymax),
              main='Mean IPD vs G',ylab='Mean IPD plus strand',xlab='G')
legend('topright',legend=paste0('r = ',round(cor(meanIPD_composition_forward$Control$G,meanIPD_composition_forward$Control$meanIPD),4)),bty='n')
smoothScatter(meanIPD_composition_forward$Control$C,meanIPD_composition_forward$Control$meanIPD,
              nrpoints=0,colramp=colorRampPalette(c("white",'red')),xlim=c(0,100),ylim=c(0,ymax),
              main='Mean IPD vs C',ylab='Mean IPD plus strand',xlab='C')
legend('topright',legend=paste0('r = ',round(cor(meanIPD_composition_forward$Control$C,meanIPD_composition_forward$Control$meanIPD),4)),bty='n')
# reverse
smoothScatter(meanIPD_composition_reverse$Control$A,meanIPD_composition_reverse$Control$meanIPD,
              nrpoints=0,colramp=colorRampPalette(c("white",'red')),xlim=c(0,100),ylim=c(0,ymax),
              main='Mean IPD vs A',ylab='Mean IPD minus strand',xlab='A')
legend('topright',legend=paste0('r = ',round(cor(meanIPD_composition_reverse$Control$A,meanIPD_composition_reverse$Control$meanIPD),4)),bty='n')
smoothScatter(meanIPD_composition_reverse$Control$T,meanIPD_composition_reverse$Control$meanIPD,
              nrpoints=0,colramp=colorRampPalette(c("white",'red')),xlim=c(0,100),ylim=c(0,ymax),
              main='Mean IPD vs T',ylab='Mean IPD minus strand',xlab='T')
legend('topright',legend=paste0('r = ',round(cor(meanIPD_composition_reverse$Control$T,meanIPD_composition_reverse$Control$meanIPD),4)),bty='n')
smoothScatter(meanIPD_composition_reverse$Control$G,meanIPD_composition_reverse$Control$meanIPD,
              nrpoints=0,colramp=colorRampPalette(c("white",'red')),xlim=c(0,100),ylim=c(0,ymax),
              main='Mean IPD vs G',ylab='Mean IPD minus strand',xlab='G')
legend('topright',legend=paste0('r = ',round(cor(meanIPD_composition_reverse$Control$G,meanIPD_composition_reverse$Control$meanIPD),4)),bty='n')
smoothScatter(meanIPD_composition_reverse$Control$C,meanIPD_composition_reverse$Control$meanIPD,
              nrpoints=0,colramp=colorRampPalette(c("white",'red')),xlim=c(0,100),ylim=c(0,ymax),
              main='Mean IPD vs C',ylab='Mean IPD minus strand',xlab='C')
legend('topright',legend=paste0('r = ',round(cor(meanIPD_composition_reverse$Control$C,meanIPD_composition_reverse$Control$meanIPD),4)),bty='n')
# together
smoothScatter(meanIPD_composition$Control$A,meanIPD_composition$Control$meanIPD,
              nrpoints=0,colramp=colorRampPalette(c("white",'red')),xlim=c(0,100),ylim=c(0,ymax),
              main='Mean IPD vs A',ylab='Mean IPD',xlab='A')
legend('topright',legend=paste0('r = ',round(cor(meanIPD_composition$Control$A,meanIPD_composition$Control$meanIPD),4)),bty='n')
smoothScatter(meanIPD_composition$Control$T,meanIPD_composition$Control$meanIPD,
              nrpoints=0,colramp=colorRampPalette(c("white",'red')),xlim=c(0,100),ylim=c(0,ymax),
              main='Mean IPD vs T',ylab='Mean IPD',xlab='T')
legend('topright',legend=paste0('r = ',round(cor(meanIPD_composition$Control$T,meanIPD_composition$Control$meanIPD),4)),bty='n')
smoothScatter(meanIPD_composition$Control$G,meanIPD_composition$Control$meanIPD,
              nrpoints=0,colramp=colorRampPalette(c("white",'red')),xlim=c(0,100),ylim=c(0,ymax),
              main='Mean IPD vs G',ylab='Mean IPD',xlab='G')
legend('topright',legend=paste0('r = ',round(cor(meanIPD_composition$Control$G,meanIPD_composition$Control$meanIPD),4)),bty='n')
smoothScatter(meanIPD_composition$Control$C,meanIPD_composition$Control$meanIPD,
              nrpoints=0,colramp=colorRampPalette(c("white",'red')),xlim=c(0,100),ylim=c(0,ymax),
              main='Mean IPD vs C',ylab='Mean IPD',xlab='C')
legend('topright',legend=paste0('r = ',round(cor(meanIPD_composition$Control$C,meanIPD_composition$Control$meanIPD),4)),bty='n')
dev.off()

pdf('meanIPD_composition_log_control.pdf',8,8)
par(mfrow=c(2,2))
ymin=log(min(meanIPD_composition_forward$Control$meanIPD,meanIPD_composition_reverse$Control$meanIPD)+0.01)
ymax=log(max(meanIPD_composition_forward$Control$meanIPD,meanIPD_composition_reverse$Control$meanIPD)+0.01)
# forward
smoothScatter(meanIPD_composition_forward$Control$A,log(meanIPD_composition_forward$Control$meanIPD+0.01),
              nrpoints=0,colramp=colorRampPalette(c("white",'red')),xlim=c(0,100),ylim=c(ymin,ymax),
              main='Mean IPD vs A',ylab='Log mean IPD plus strand',xlab='A')
legend('topright',legend=paste0('r = ',round(cor(meanIPD_composition_forward$Control$A,log(meanIPD_composition_forward$Control$meanIPD+0.01)),4)),bty='n')
smoothScatter(meanIPD_composition_forward$Control$T,log(meanIPD_composition_forward$Control$meanIPD+0.01),
              nrpoints=0,colramp=colorRampPalette(c("white",'red')),xlim=c(0,100),ylim=c(ymin,ymax),
              main='Mean IPD vs T',ylab='Log mean IPD plus strand',xlab='T')
legend('topright',legend=paste0('r = ',round(cor(meanIPD_composition_forward$Control$T,log(meanIPD_composition_forward$Control$meanIPD+0.01)),4)),bty='n')
smoothScatter(meanIPD_composition_forward$Control$G,log(meanIPD_composition_forward$Control$meanIPD+0.01),
              nrpoints=0,colramp=colorRampPalette(c("white",'red')),xlim=c(0,100),ylim=c(ymin,ymax),
              main='Mean IPD vs G',ylab='Log mean IPD plus strand',xlab='G')
legend('topright',legend=paste0('r = ',round(cor(meanIPD_composition_forward$Control$G,log(meanIPD_composition_forward$Control$meanIPD+0.01)),4)),bty='n')
smoothScatter(meanIPD_composition_forward$Control$C,log(meanIPD_composition_forward$Control$meanIPD+0.01),
              nrpoints=0,colramp=colorRampPalette(c("white",'red')),xlim=c(0,100),ylim=c(ymin,ymax),
              main='Mean IPD vs C',ylab='Log mean IPD plus strand',xlab='C')
legend('topright',legend=paste0('r = ',round(cor(meanIPD_composition_forward$Control$C,log(meanIPD_composition_forward$Control$meanIPD+0.01)),4)),bty='n')
# reverse
smoothScatter(meanIPD_composition_reverse$Control$A,log(meanIPD_composition_reverse$Control$meanIPD+0.01),
              nrpoints=0,colramp=colorRampPalette(c("white",'red')),xlim=c(0,100),ylim=c(ymin,ymax),
              main='Mean IPD vs A',ylab='Log mean IPD minus strand',xlab='A')
legend('topright',legend=paste0('r = ',round(cor(meanIPD_composition_reverse$Control$A,log(meanIPD_composition_reverse$Control$meanIPD+0.01)),4)),bty='n')
smoothScatter(meanIPD_composition_reverse$Control$T,log(meanIPD_composition_reverse$Control$meanIPD+0.01),
              nrpoints=0,colramp=colorRampPalette(c("white",'red')),xlim=c(0,100),ylim=c(ymin,ymax),
              main='Mean IPD vs T',ylab='Log mean IPD minus strand',xlab='T')
legend('topright',legend=paste0('r = ',round(cor(meanIPD_composition_reverse$Control$T,log(meanIPD_composition_reverse$Control$meanIPD+0.01)),4)),bty='n')
smoothScatter(meanIPD_composition_reverse$Control$G,log(meanIPD_composition_reverse$Control$meanIPD+0.01),
              nrpoints=0,colramp=colorRampPalette(c("white",'red')),xlim=c(0,100),ylim=c(ymin,ymax),
              main='Mean IPD vs G',ylab='Log mean IPD minus strand',xlab='G')
legend('topright',legend=paste0('r = ',round(cor(meanIPD_composition_reverse$Control$G,log(meanIPD_composition_reverse$Control$meanIPD+0.01)),4)),bty='n')
smoothScatter(meanIPD_composition_reverse$Control$C,log(meanIPD_composition_reverse$Control$meanIPD+0.01),
              nrpoints=0,colramp=colorRampPalette(c("white",'red')),xlim=c(0,100),ylim=c(ymin,ymax),
              main='Mean IPD vs C',ylab='Log mean IPD minus strand',xlab='C')
legend('topright',legend=paste0('r = ',round(cor(meanIPD_composition_reverse$Control$C,log(meanIPD_composition_reverse$Control$meanIPD+0.01)),4)),bty='n')
# together
smoothScatter(meanIPD_composition$Control$A,log(meanIPD_composition$Control$meanIPD+0.01),
              nrpoints=0,colramp=colorRampPalette(c("white",'red')),xlim=c(0,100),ylim=c(ymin,ymax),
              main='Mean IPD vs A',ylab='Log mean IPD',xlab='A')
legend('topright',legend=paste0('r = ',round(cor(meanIPD_composition$Control$A,log(meanIPD_composition$Control$meanIPD+0.01)),4)),bty='n')
smoothScatter(meanIPD_composition$Control$T,log(meanIPD_composition$Control$meanIPD+0.01),
              nrpoints=0,colramp=colorRampPalette(c("white",'red')),xlim=c(0,100),ylim=c(ymin,ymax),
              main='Mean IPD vs T',ylab='Log mean IPD',xlab='T')
legend('topright',legend=paste0('r = ',round(cor(meanIPD_composition$Control$T,log(meanIPD_composition$Control$meanIPD+0.01)),4)),bty='n')
smoothScatter(meanIPD_composition$Control$G,log(meanIPD_composition$Control$meanIPD+0.01),
              nrpoints=0,colramp=colorRampPalette(c("white",'red')),xlim=c(0,100),ylim=c(ymin,ymax),
              main='Mean IPD vs G',ylab='Log mean IPD',xlab='G')
legend('topright',legend=paste0('r = ',round(cor(meanIPD_composition$Control$G,log(meanIPD_composition$Control$meanIPD+0.01)),4)),bty='n')
smoothScatter(meanIPD_composition$Control$C,log(meanIPD_composition$Control$meanIPD+0.01),
              nrpoints=0,colramp=colorRampPalette(c("white",'red')),xlim=c(0,100),ylim=c(ymin,ymax),
              main='Mean IPD vs C',ylab='Log mean IPD',xlab='C')
legend('topright',legend=paste0('r = ',round(cor(meanIPD_composition$Control$C,log(meanIPD_composition$Control$meanIPD+0.01)),4)),bty='n')
dev.off()


pdf('meanIPD_composition_control_correlation.pdf',8,8)
panel.cor <- function(x,y,digits=2,method_cor='pearson',prefix="",cex.cor, ...){
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  index=(!is.na(x))&(!is.na(y))
  r <- cor(x[index], y[index],method=method_cor)
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor*(0.4+0.8*abs(r)))
}
upper.panel=function(...) {par(new=TRUE);smoothScatter(...,nrpoints=0,colramp=colorRampPalette(c("white","red")),add=TRUE)}
par(mfrow=c(2,2))
# forward
pairs(meanIPD_composition_forward$Control[,2:5],main='Composition correlation plus strand',labels=c('A','T','G','C'),pch=3,lower.panel=panel.cor,upper.panel=upper.panel)
# reverse
pairs(meanIPD_composition_reverse$Control[,2:5],main='Composition correlation minus strand',labels=c('A','T','G','C'),pch=3,lower.panel=panel.cor,upper.panel=upper.panel)
# together
pairs(meanIPD_composition$Control[,2:5],main='Composition correlation',labels=c('A','T','G','C'),pch=3,lower.panel=panel.cor,upper.panel=upper.panel)
dev.off()











##############################
# linear regression controls #
#  compositional regressors  #
##############################
# Isometric log ratio transform
meanIPD_composition_forward=lapply(meanIPD_composition_forward,function(region) cbind(region,ilr(region[,-1])))
meanIPD_composition_reverse=lapply(meanIPD_composition_reverse,function(region) cbind(region,ilr(region[,-1])))
meanIPD_composition=lapply(meanIPD_composition,function(region) cbind(region,ilr(region[,-1])))
save(meanIPD_composition_forward,meanIPD_composition_reverse,meanIPD_composition,file='meanIPD_composition.RData')


load('meanIPD_composition.RData')

# forward
lm_forward_control=lm(meanIPD~V1+V2+V3,meanIPD_composition_forward$Control)
summary(lm_forward_control)
# reverse
lm_reverse_control=lm(meanIPD~V1+V2+V3,meanIPD_composition_reverse$Control)
summary(lm_reverse_control)
# together
dummy=c(rep(0,nrow(meanIPD_composition_forward$Control)),rep(1,nrow(meanIPD_composition_reverse$Control)))
lm_control=lm(meanIPD~V1+V2+V3+dummy+V1*dummy+V2*dummy+V3*dummy,meanIPD_composition$Control)
summary(lm_control)
lm_control=lm(meanIPD~V1+V2+V3,meanIPD_composition$Control)
summary(lm_control)
ilrInv(coef(lm_control)[-1],orig=meanIPD_composition$Control[,2:5])


pdf('meanIPD_composition_lm_ilr.pdf',width=7,height=7)
smoothScatter(meanIPD_composition$Control$meanIPD,lm_control$fitted.values,
              nrpoints=0,colramp=colorRampPalette(c("white",'red')),
              main='Fitted mean IPD vs observed',ylab='Fitted',xlab='Observed')
abline(0,1)
smoothScatter(lm_control$residuals,
              nrpoints=0,colramp=colorRampPalette(c("white",'red')),
              main='Residuals',ylab='Residuals',xlab='Index')
abline(h=0)
par(mfrow=c(2,2))
smoothScatter(meanIPD_composition$Control$V1,lm_control$residuals,
              nrpoints=0,colramp=colorRampPalette(c("white",'red')),
              main='Residuals vs V1',ylab='Residuals',xlab='V1')
abline(h=0)
smoothScatter(meanIPD_composition$Control$V1,lm_control$residuals,
              nrpoints=0,colramp=colorRampPalette(c("white",'red')),
              main='Residuals vs V2',ylab='Residuals',xlab='V2')
abline(h=0)
smoothScatter(meanIPD_composition$Control$V3,lm_control$residuals,
              nrpoints=0,colramp=colorRampPalette(c("white",'red')),
              main='Residuals vs V3',ylab='Residuals',xlab='V3')
abline(h=0)
dev.off()


# linear regression log controls
# forward
lm_forward_control_log=lm(log(meanIPD+0.01)~V1+V2+V3,meanIPD_composition_forward$Control)
summary(lm_forward_control_log)
# reverse
lm_reverse_control_log=lm(log(meanIPD+0.01)~V1+V2+V3,meanIPD_composition_reverse$Control)
summary(lm_reverse_control_log)
# together
dummy=c(rep(0,nrow(meanIPD_composition_forward$Control)),rep(1,nrow(meanIPD_composition_reverse$Control)))
lm_control_log=lm(log(meanIPD+0.01)~V1+V2+V3+dummy+V1*dummy+V2*dummy+V3*dummy,meanIPD_composition$Control)
summary(lm_control_log)
lm_control_log=lm(log(meanIPD+0.01)~V1+V2+V3,meanIPD_composition$Control)
summary(lm_control_log)
ilrInv(coef(lm_control_log)[-1],orig=meanIPD_composition$Control[,2:5])


pdf('meanIPD_composition_lm_log_ilr.pdf',width=7,height=7)
smoothScatter(log(meanIPD_composition$Control$meanIPD+0.01),lm_control_log$fitted.values,
              nrpoints=0,colramp=colorRampPalette(c("white",'red')),
              main='Fitted mean IPD vs observed',ylab='Fitted',xlab='Observed')
abline(0,1)
smoothScatter(lm_control_log$residuals,
              nrpoints=0,colramp=colorRampPalette(c("white",'red')),
              main='Residuals',ylab='Residuals',xlab='Index')
abline(h=0)
par(mfrow=c(2,2))
smoothScatter(meanIPD_composition$Control$V1,lm_control_log$residuals,
              nrpoints=0,colramp=colorRampPalette(c("white",'red')),
              main='Residuals vs V1',ylab='Residuals',xlab='V1')
abline(h=0)
smoothScatter(meanIPD_composition$Control$V2,lm_control_log$residuals,
              nrpoints=0,colramp=colorRampPalette(c("white",'red')),
              main='Residuals vs V2',ylab='Residuals',xlab='V2')
abline(h=0)
smoothScatter(meanIPD_composition$Control$V3,lm_control_log$residuals,
              nrpoints=0,colramp=colorRampPalette(c("white",'red')),
              main='Residuals vs V3',ylab='Residuals',xlab='V3')
abline(h=0)
dev.off()




###########################
# prediction IPD features #
###########################
meanIPD_forward_predict=lapply(meanIPD_composition_forward,
                               function(meanIPD){
                                 meanIPD_predicted=predict(lm_control,meanIPD)
                                 cbind(meanIPD,meanIPD_predicted)
                               })
meanIPD_reverse_predict=lapply(meanIPD_composition_reverse,
                               function(meanIPD){
                                 meanIPD_predicted=predict(lm_control,meanIPD)
                                 cbind(meanIPD,meanIPD_predicted)
                               })
meanIPD_predict=lapply(meanIPD_composition,
                       function(meanIPD){
                         meanIPD_predicted=predict(lm_control,meanIPD)
                         cbind(meanIPD,meanIPD_predicted)
                       })

# scatterplot controls
pdf('meanIPD_composition_control_ilr.pdf',8,8)
par(mfrow=c(2,2))
ymax=max(meanIPD_predict$Control[,c(1,9)])
# together
smoothScatter(meanIPD_predict$Control$A,meanIPD_predict$Control$meanIPD,
              nrpoints=0,colramp=colorRampPalette(c("white",'red')),xlim=c(0,100),ylim=c(0,ymax),
              main='Mean IPD vs A',ylab='Mean IPD',xlab='A')
smoothScatter(meanIPD_predict$Control$A,meanIPD_predict$Control$meanIPD_predicted,
              nrpoints=0,colramp=colorRampPalette(c(rgb(1,1,1,0),rgb(0,0,1,0.5)),alpha=TRUE),xlim=c(0,100),ylim=c(0,ymax),
              add=TRUE)
legend('topright',legend=paste0('r = ',round(cor(meanIPD_composition$Control$A,meanIPD_composition$Control$meanIPD),4)),bty='n')
smoothScatter(meanIPD_predict$Control$T,meanIPD_predict$Control$meanIPD,
              nrpoints=0,colramp=colorRampPalette(c("white",'red')),xlim=c(0,100),ylim=c(0,ymax),
              main='Mean IPD vs T',ylab='Mean IPD',xlab='T')
smoothScatter(meanIPD_predict$Control$T,meanIPD_predict$Control$meanIPD_predicted,
              nrpoints=0,colramp=colorRampPalette(c(rgb(1,1,1,0),rgb(0,0,1,0.5)),alpha=TRUE),xlim=c(0,100),ylim=c(0,ymax),
              add=TRUE)
legend('topright',legend=paste0('r = ',round(cor(meanIPD_composition$Control$T,meanIPD_composition$Control$meanIPD),4)),bty='n')
smoothScatter(meanIPD_predict$Control$G,meanIPD_predict$Control$meanIPD,
              nrpoints=0,colramp=colorRampPalette(c("white",'red')),xlim=c(0,100),ylim=c(0,ymax),
              main='Mean IPD vs G',ylab='Mean IPD',xlab='G')
smoothScatter(meanIPD_predict$Control$G,meanIPD_predict$Control$meanIPD_predicted,
              nrpoints=0,colramp=colorRampPalette(c(rgb(1,1,1,0),rgb(0,0,1,0.5)),alpha=TRUE),xlim=c(0,100),ylim=c(0,ymax),
              add=TRUE)
legend('topright',legend=paste0('r = ',round(cor(meanIPD_composition$Control$G,meanIPD_composition$Control$meanIPD),4)),bty='n')
smoothScatter(meanIPD_predict$Control$C,meanIPD_predict$Control$meanIPD,
              nrpoints=0,colramp=colorRampPalette(c("white",'red')),xlim=c(0,100),ylim=c(0,ymax),
              main='Mean IPD vs C',ylab='Mean IPD',xlab='C')
smoothScatter(meanIPD_predict$Control$C,meanIPD_predict$Control$meanIPD_predicted,
              nrpoints=0,colramp=colorRampPalette(c(rgb(1,1,1,0),rgb(0,0,1,0.5)),alpha=TRUE),xlim=c(0,100),ylim=c(0,ymax),
              add=TRUE)
legend('topright',legend=paste0('r = ',round(cor(meanIPD_composition$Control$C,meanIPD_composition$Control$meanIPD),4)),bty='n')
dev.off()


meanIPD_forward_residual=lapply(meanIPD_forward_predict,
                                function(meanIPD) meanIPD$meanIPD-meanIPD$meanIPD_predicted)
meanIPD_reverse_residual=lapply(meanIPD_reverse_predict,
                                function(meanIPD) meanIPD$meanIPD-meanIPD$meanIPD_predicted)
datasets=read.table("C:/Users/MarziaAngela/Google Drive/PacBio/datasets.txt",sep="\t",header=TRUE,stringsAsFactors=FALSE,row.names=1)

pdf('meanIPD_composition_boxplot_residual_ilr.pdf',width=7,height=15)
par(mar=c(5,9,4,2)+0.1)
pval_forward=unlist(lapply(meanIPD_forward_residual,function(res) t.test(res)$p.value))
pval_forward=p.adjust(pval_forward,method='bonferroni')
pval_forward_code=rep('',length(pval_forward))
pval_forward_code[pval_forward<=0.1]='.'
pval_forward_code[pval_forward<=0.05]='*'
pval_forward_code[pval_forward<=0.01]='**'
pval_forward_code[pval_forward<=0.001]='***'
pval_forward_code[pval_forward<=0.0001]='****'
statistic_forward=unlist(lapply(meanIPD_forward_residual,function(res) t.test(res)$statistic))
col_pval_forward=ifelse(statistic_forward>0,'red','blue')
box=boxplot(rev(meanIPD_forward_residual),plot=FALSE)
box$stats[c(1,5),]=Reduce(cbind,lapply(rev(meanIPD_forward_residual),quantile,probs=c(0.05,0.95),na.rm=TRUE,SIMPLIFY=FALSE))
ylim=range(box$stats)+c(-0.2,0)
box$names=rev(datasets[names(meanIPD_forward_residual),'name'])
bxp(box,at=c(-4,-1:42,44:63,65:70,72:75,77:83),border=c('red',rep('blue',82)),
    horizontal=TRUE,cex.main=1.8,ylim=ylim,xlim=c(-4,83),outline=FALSE,pars=list(whisklty=3,staplewex=0.8),
    xlab='Residual mean IPD',main='Residual IPD plus strand',las=1)
abline(v=box$stats[,1],col='red')
text(x=ylim[1]+0.1,y=c(-4,-1:42,44:63,65:70,72:75,77:83),labels=rev(pval_forward_code),col=rev(col_pval_forward),cex=1.8)
pval_reverse=unlist(lapply(meanIPD_reverse_residual,function(res) t.test(res)$p.value))
pval_reverse=p.adjust(pval_reverse,method='bonferroni')
pval_reverse_code=rep('',length(pval_reverse))
pval_reverse_code[pval_reverse<=0.1]='.'
pval_reverse_code[pval_reverse<=0.05]='*'
pval_reverse_code[pval_reverse<=0.01]='**'
pval_reverse_code[pval_reverse<=0.001]='***'
pval_reverse_code[pval_reverse<=0.0001]='****'
statistic_reverse=unlist(lapply(meanIPD_reverse_residual,function(res) t.test(res)$statistic))
col_pval_reverse=ifelse(statistic_reverse>0,'red','blue')
box=boxplot(rev(meanIPD_reverse_residual),plot=FALSE)
box$stats[c(1,5),]=Reduce(cbind,lapply(rev(meanIPD_reverse_residual),quantile,probs=c(0.05,0.95),na.rm=TRUE,SIMPLIFY=FALSE))
box$names=rev(datasets[names(meanIPD_reverse_residual),'name'])
bxp(box,at=c(-4,-1:42,44:63,65:70,72:75,77:83),border=c('red',rep('blue',82)),
    horizontal=TRUE,cex.main=1.8,ylim=ylim,xlim=c(-4,83),outline=FALSE,pars=list(whisklty=3,staplewex=0.8),
    xlab='Residual mean IPD',main='Residual IPD minus strand',las=1)
abline(v=box$stats[,1],col='red')
text(x=ylim[1]+0.1,y=c(-4,-1:42,44:63,65:70,72:75,77:83),labels=rev(pval_reverse_code),col=rev(col_pval_reverse),cex=1.8)
dev.off()


# log
meanIPD_forward_log_predict=lapply(meanIPD_composition_forward,
                                   function(meanIPD){
                                     meanIPD_predicted=predict(lm_control_log,meanIPD)
                                     meanIPD$meanIPD=log(meanIPD$meanIPD+0.01)
                                     cbind(meanIPD,meanIPD_predicted)
                                   })
meanIPD_reverse_log_predict=lapply(meanIPD_composition_reverse,
                                   function(meanIPD){
                                     meanIPD_predicted=predict(lm_control_log,meanIPD)
                                     meanIPD$meanIPD=log(meanIPD$meanIPD+0.01)
                                     cbind(meanIPD,meanIPD_predicted)
                                   })
meanIPD_log_predict=lapply(meanIPD_composition,
                           function(meanIPD){
                             meanIPD_predicted=predict(lm_control_log,meanIPD)
                             meanIPD$meanIPD=log(meanIPD$meanIPD+0.01)
                             cbind(meanIPD,meanIPD_predicted)
                           })

# scatterplot controls
pdf('meanIPD_composition_log_control_ilr.pdf',8,8)
par(mfrow=c(2,2))
ymin=min(meanIPD_log_predict$Control[,c(1,9)])
ymax=max(meanIPD_log_predict$Control[,c(1,9)])
# together
smoothScatter(meanIPD_log_predict$Control$A,meanIPD_log_predict$Control$meanIPD,
              nrpoints=0,colramp=colorRampPalette(c("white",'red')),xlim=c(0,100),ylim=c(ymin,ymax),
              main='Mean IPD vs A',ylab='Log mean IPD',xlab='A')
smoothScatter(meanIPD_log_predict$Control$A,meanIPD_log_predict$Control$meanIPD_predicted,
              nrpoints=0,colramp=colorRampPalette(c(rgb(1,1,1,0),rgb(0,0,1,0.5)),alpha=TRUE),xlim=c(0,100),ylim=c(ymin,ymax),
              add=TRUE)
legend('topright',legend=paste0('r = ',round(cor(meanIPD_composition_forward$Control$A,log(meanIPD_composition_forward$Control$meanIPD+0.01)),4)),bty='n')
smoothScatter(meanIPD_log_predict$Control$T,meanIPD_log_predict$Control$meanIPD,
              nrpoints=0,colramp=colorRampPalette(c("white",'red')),xlim=c(0,100),ylim=c(ymin,ymax),
              main='Mean IPD vs T',ylab='Log mean IPD',xlab='T')
smoothScatter(meanIPD_log_predict$Control$T,meanIPD_log_predict$Control$meanIPD_predicted,
              nrpoints=0,colramp=colorRampPalette(c(rgb(1,1,1,0),rgb(0,0,1,0.5)),alpha=TRUE),xlim=c(0,100),ylim=c(ymin,ymax),
              add=TRUE)
legend('topright',legend=paste0('r = ',round(cor(meanIPD_composition_forward$Control$T,log(meanIPD_composition_forward$Control$meanIPD+0.01)),4)),bty='n')
smoothScatter(meanIPD_log_predict$Control$G,meanIPD_log_predict$Control$meanIPD,
              nrpoints=0,colramp=colorRampPalette(c("white",'red')),xlim=c(0,100),ylim=c(ymin,ymax),
              main='Mean IPD vs G',ylab='Log mean IPD',xlab='G')
smoothScatter(meanIPD_log_predict$Control$G,meanIPD_log_predict$Control$meanIPD_predicted,
              nrpoints=0,colramp=colorRampPalette(c(rgb(1,1,1,0),rgb(0,0,1,0.5)),alpha=TRUE),xlim=c(0,100),ylim=c(ymin,ymax),
              add=TRUE)
legend('topright',legend=paste0('r = ',round(cor(meanIPD_composition_forward$Control$G,log(meanIPD_composition_forward$Control$meanIPD+0.01)),4)),bty='n')
smoothScatter(meanIPD_log_predict$Control$C,meanIPD_log_predict$Control$meanIPD,
              nrpoints=0,colramp=colorRampPalette(c("white",'red')),xlim=c(0,100),ylim=c(ymin,ymax),
              main='Mean IPD vs C',ylab='Log mean IPD',xlab='C')
smoothScatter(meanIPD_log_predict$Control$C,meanIPD_log_predict$Control$meanIPD_predicted,
              nrpoints=0,colramp=colorRampPalette(c(rgb(1,1,1,0),rgb(0,0,1,0.5)),alpha=TRUE),xlim=c(0,100),ylim=c(ymin,ymax),
              add=TRUE)
legend('topright',legend=paste0('r = ',round(cor(meanIPD_composition_forward$Control$C,log(meanIPD_composition_forward$Control$meanIPD+0.01)),4)),bty='n')
dev.off()


meanIPD_forward_log_residual=lapply(meanIPD_forward_log_predict,
                                    function(meanIPD) meanIPD$meanIPD-meanIPD$meanIPD_predicted)
meanIPD_reverse_log_residual=lapply(meanIPD_reverse_log_predict,
                                    function(meanIPD) meanIPD$meanIPD-meanIPD$meanIPD_predicted)
datasets=read.table("C:/Users/MarziaAngela/Google Drive/PacBio/datasets.txt",sep="\t",header=TRUE,stringsAsFactors=FALSE,row.names=1)

pdf('meanIPD_composition_boxplot_residual_log_ilr.pdf',width=7,height=15)
par(mar=c(5,9,4,2)+0.1)
pval_forward=unlist(lapply(meanIPD_forward_log_residual,function(res) t.test(res)$p.value))
pval_forward=p.adjust(pval_forward,method='bonferroni')
pval_forward_code=rep('',length(pval_forward))
pval_forward_code[pval_forward<=0.1]='.'
pval_forward_code[pval_forward<=0.05]='*'
pval_forward_code[pval_forward<=0.01]='**'
pval_forward_code[pval_forward<=0.001]='***'
pval_forward_code[pval_forward<=0.0001]='****'
statistic_forward=unlist(lapply(meanIPD_forward_log_residual,function(res) t.test(res)$statistic))
col_pval_forward=ifelse(statistic_forward>0,'red','blue')
box=boxplot(rev(meanIPD_forward_log_residual),plot=FALSE)
box$stats[c(1,5),]=Reduce(cbind,lapply(rev(meanIPD_forward_log_residual),quantile,probs=c(0.05,0.95),na.rm=TRUE,SIMPLIFY=FALSE))
ylim=range(box$stats)+c(-0.2,0)
box$names=rev(datasets[names(meanIPD_forward_log_residual),'name'])
bxp(box,at=c(-4,-1:42,44:63,65:70,72:75,77:83),border=c('red',rep('blue',82)),
    horizontal=TRUE,cex.main=1.8,ylim=ylim,xlim=c(-4,83),outline=FALSE,pars=list(whisklty=3,staplewex=0.8),
    xlab='Residual log mean IPD',main='Residual IPD plus strand',las=1)
abline(v=box$stats[,1],col='red')
text(x=ylim[1]+0.1,y=c(-4,-1:42,44:63,65:70,72:75,77:83),labels=rev(pval_forward_code),col=rev(col_pval_forward),cex=1.8)
pval_reverse=unlist(lapply(meanIPD_reverse_log_residual,function(res) t.test(res)$p.value))
pval_reverse=p.adjust(pval_reverse,method='bonferroni')
pval_reverse_code=rep('',length(pval_reverse))
pval_reverse_code[pval_reverse<=0.1]='.'
pval_reverse_code[pval_reverse<=0.05]='*'
pval_reverse_code[pval_reverse<=0.01]='**'
pval_reverse_code[pval_reverse<=0.001]='***'
pval_reverse_code[pval_reverse<=0.0001]='****'
statistic_reverse=unlist(lapply(meanIPD_reverse_log_residual,function(res) t.test(res)$statistic))
col_pval_reverse=ifelse(statistic_reverse>0,'red','blue')
box=boxplot(rev(meanIPD_reverse_log_residual),plot=FALSE)
box$stats[c(1,5),]=Reduce(cbind,lapply(rev(meanIPD_reverse_log_residual),quantile,probs=c(0.05,0.95),na.rm=TRUE,SIMPLIFY=FALSE))
box$names=rev(datasets[names(meanIPD_reverse_log_residual),'name'])
bxp(box,at=c(-4,-1:42,44:63,65:70,72:75,77:83),border=c('red',rep('blue',82)),
    horizontal=TRUE,cex.main=1.8,ylim=ylim,xlim=c(-4,83),outline=FALSE,pars=list(whisklty=3,staplewex=0.8),
    xlab='Residual log mean IPD',main='Residual IPD minus strand',las=1)
abline(v=box$stats[,1],col='red')
text(x=ylim[1]+0.1,y=c(-4,-1:42,44:63,65:70,72:75,77:83),labels=rev(pval_reverse_code),col=rev(col_pval_reverse),cex=1.8)
dev.off()
