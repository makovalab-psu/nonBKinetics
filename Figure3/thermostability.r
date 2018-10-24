# Input files intramolecular.out and intermolecular.out 
# are text file (tab delimited) containing the IPDs for all occurrences 
# of the 10 most common G-quadruplex motifs (divided in intra and intermolecular)
# with the format: motif  IPD_1   IPD_2   ...   IPD_l
# GGGGTGGGGGGAGGGGGGAGGG	0.518	0.705	0.649	0.768	1.462	1.007	0.623	1.01	1.322	2.933	6.575	4.915	1.026	0.92	2.137	0.809	0.788	0.729	0.888	0.538	0.535	0.753
# GGGGTGGGGGGAGGGGGGAGGG	0.953	1.16	1.084	1.16	0.929	1.026	0.782	0.877	1.213	2.082	3.357	2.446	1.051	0.842	0.819	1.026	1.005	0.737	1.122	0.605	0.66	0.895
# ...
intramolecular <- read.table('intramolecular.out', fill=TRUE)
G1 <- intramolecular[intramolecular$V1=='GGGGTGGGGGGAGGGGGGAGGG',]
G3 <- intramolecular[intramolecular$V1=='GGGGTCGGGGGAGGGGGGAGGG',]
G4 <- intramolecular[intramolecular$V1=='GGGGTGGGGGGAGTGGGGAGGG',]
G8 <- intramolecular[intramolecular$V1=='GGGGTTGGGGGAGGGGGGAGGG',]
G9 <- intramolecular[intramolecular$V1=='GGGGTGGGGGGAGGGGGAGGG',]
G10 <- intramolecular[intramolecular$V1=='GGGGTGGGGGGAGCGGGGAGGG',]

intermolecular <- read.table('intermolecular.out', fill=TRUE)
G2 <- intermolecular[intermolecular$V1=='GGGAGGGAGGTGGGGGGG',]
G5 <- intermolecular[intermolecular$V1=='GGGAGGGAGGGAGGGAGGG',]
G6 <- intermolecular[intermolecular$V1=='GGGAGGGAGGTGGGGGGGG',]
G7 <- intermolecular[intermolecular$V1=='GGGTGGAGGGTGGGAGGAGGG',]


# Compute meanIPD
mean_G1 <- vector()
for (i in 1:length(G1[,1])) {mean_G1 <- append(mean_G1, mean(as.numeric(G1[i,2:length(G1[1,])]),na.rm=TRUE))}
mean_G1=mean_G1[!is.nan(mean_G1)]
mean_G2 <- vector()
for (i in 1:length(G2[,1])) {mean_G2 <- append(mean_G2, mean(as.numeric(G2[i,2:length(G2[1,])]),na.rm=TRUE))}
mean_G2=mean_G2[!is.nan(mean_G2)]
mean_G3 <- vector()
for (i in 1:length(G3[,1])) {mean_G3 <- append(mean_G3, mean(as.numeric(G3[i,2:length(G3[1,])]),na.rm=TRUE))}
mean_G3=mean_G3[!is.nan(mean_G3)]
mean_G4 <- vector()
for (i in 1:length(G4[,1])) {mean_G4 <- append(mean_G4, mean(as.numeric(G4[i,2:length(G4[1,])]),na.rm=TRUE))}
mean_G4=mean_G4[!is.nan(mean_G4)]
mean_G5 <- vector()
for (i in 1:length(G5[,1])) {mean_G5 <- append(mean_G5, mean(as.numeric(G5[i,2:length(G5[1,])]),na.rm=TRUE))}
mean_G5=mean_G5[!is.nan(mean_G5)]
mean_G6 <- vector()
for (i in 1:length(G6[,1])) {mean_G6 <- append(mean_G6, mean(as.numeric(G6[i,2:length(G6[1,])]),na.rm=TRUE))}
mean_G6=mean_G6[!is.nan(mean_G6)]
mean_G7 <- vector()
for (i in 1:length(G7[,1])) {mean_G7 <- append(mean_G7, mean(as.numeric(G7[i,2:length(G7[1,])]),na.rm=TRUE))}
mean_G7=mean_G7[!is.nan(mean_G7)]
mean_G8 <- vector()
for (i in 1:length(G8[,1])) {mean_G8 <- append(mean_G8, mean(as.numeric(G8[i,2:length(G8[1,])]),na.rm=TRUE))}
mean_G8=mean_G8[!is.nan(mean_G8)]
mean_G9 <- vector()
for (i in 1:length(G9[,1])) {mean_G9 <- append(mean_G9, mean(as.numeric(G9[i,2:length(G9[1,])]),na.rm=TRUE))}
mean_G9=mean_G9[!is.nan(mean_G9)]
mean_G10 <- vector()
for (i in 1:length(G10[,1])) {mean_G10 <- append(mean_G10, mean(as.numeric(G10[i,2:length(G10[1,])]),na.rm=TRUE))}
mean_G10=mean_G10[!is.nan(mean_G10)]

intra_G <- c(mean_G1,mean_G3,mean_G4,mean_G8,mean_G9,mean_G10)
inter_G <- c(mean_G2,mean_G5,mean_G6,mean_G7)


# Compute log meanIPD
log_mean_G=list(G1=log(mean_G1+0.01),G2=log(mean_G2+0.01),G3=log(mean_G3+0.01),G4=log(mean_G4+0.01),G5=log(mean_G5+0.01),
                G6=log(mean_G6+0.01),G7=log(mean_G7+0.01),G8=log(mean_G8+0.01),G9=log(mean_G9+0.01),G10=log(mean_G10+0.01))








# Melting temperature TM from experiments
TMs <- c(74.3,64.8,74.8,69,69,68,61.5,73.2,71.9,68.5)
intra_TMs <- rep(TMs[c(1,3,4,8,9,10)],c(length(mean_G1),length(mean_G3),length(mean_G4),length(mean_G8),length(mean_G9),length(mean_G10)))
inter_TMs <- rep(TMs[c(2,5,6,7)],c(length(mean_G2),length(mean_G5),length(mean_G6),length(mean_G7)))


# Delta epsilon from experiments
epsilons <- c(248,298,216,209,300,300,282,211,281,216)
intra_epsilons <- rep(epsilons[c(1,3,4,8,9,10)],c(length(mean_G1),length(mean_G3),length(mean_G4),length(mean_G8),length(mean_G9),length(mean_G10)))
inter_epsilons <- rep(epsilons[c(2,5,6,7)],c(length(mean_G2),length(mean_G5),length(mean_G6),length(mean_G7)))



###############################
# REGRESSION MODEL WITH BOTH  #
# INTRA AND INTERMOLECULAR G4 #
# AND MOLECULARITY AS BINARY  #
#         PREDICTOR           #
###############################

# Fit final regression model for Tm (after selecting relevant predictors)
group = rep(0:1,c(length(c(mean_G1,mean_G3,mean_G4,mean_G8,mean_G9,mean_G10)),length(c(mean_G2,mean_G5,mean_G6,mean_G7))))
res=lm(c(log(intra_G+0.01),log(inter_G+0.01))~c(intra_TMs,inter_TMs)+c(intra_TMs,inter_TMs)*group+group)
summary(res)
plot(res)

pdf('termostability_Tm_intra_inter.pdf',width=6,height=6)
par(mar=c(5,4,2,2)+0.1)
box=boxplot(log_mean_G,plot=FALSE)
box$stats[c(1,5),]=Reduce(cbind,lapply(log_mean_G,quantile,probs=c(0.05,0.95),na.rm=TRUE))
bxp(box,at=TMs,pars=list(xaxt='n',boxfill=c("cyan","gold","cyan","cyan","gold","gold","gold","cyan","cyan","cyan"),boxwex=diff(range(TMs))/50+0.03,whisklty=3,staplewex=0.8),
        main="",ylab="Log mean IPD",xlab=expression(T[m]*" [°C]"),outline=FALSE,ylim=c(-0.55,0.85),
        cex.main=2,cex.axis=1,cex.lab=1)
text(TMs,c(box$stats[5,1:3]+0.05,box$stats[1,4]-0.05,box$stats[5,5:10]+0.05),
     c("G1","G2","G3","G4","G5","G6","G7","G8","G9","G10"), 
     cex=1,col="black")
axis(1)

coef=res$coefficients
abline(coef[1],coef[2],col="cyan",lwd=3)
abline(coef[1]+coef[3],coef[2]+coef[4],col="gold",lwd=3)
bxp(box,at=TMs,pars=list(xaxt='n',boxfill=c("cyan","gold","cyan","cyan","gold","gold","gold","cyan","cyan","cyan"),boxwex=diff(range(TMs))/50+0.03,whisklty=3,staplewex=0.8),
    outline=FALSE,ylim=c(-0.55,0.85),
    cex.main=2,cex.axis=1,cex.lab=1,add=TRUE)

legend('topleft',legend=c('Inter','Intra'),col=c('gold','cyan'),lwd=3,border=FALSE)
dev.off()




# Fit final regression model for delta-epsilon (after selecting relevant predictors)
group = rep(0:1,c(length(c(mean_G1,mean_G3,mean_G4,mean_G8,mean_G9,mean_G10)),length(c(mean_G2,mean_G5,mean_G6,mean_G7))))
res=lm(c(log(intra_G+0.01),log(inter_G+0.01))~c(intra_epsilons,inter_epsilons)+c(intra_epsilons,inter_epsilons)*group-group)
summary(res)
plot(res)

pdf('termostability_deltaepsilon_intra_inter.pdf',width=6,height=6)
par(mar=c(5,4,2,2)+0.1)
box=boxplot(log_mean_G,plot=FALSE)
box$stats[c(1,5),]=Reduce(cbind,lapply(log_mean_G,quantile,probs=c(0.05,0.95),na.rm=TRUE))
bxp(box,at=epsilons,pars=list(xaxt='n',boxfill=c("cyan","gold","cyan","cyan","gold","gold","gold","cyan","cyan","cyan"),boxwex=diff(range(epsilons))/50,whisklty=3,staplewex=0.8),
    main="",ylab="Log mean IPD",xlab=expression(Delta*epsilon*" ["*M^{-1}*cm^{-1}*"]"),outline=FALSE,ylim=c(-0.55,0.85),
    cex.main=2,cex.axis=1,cex.lab=1)
text(epsilons,c(box$stats[5,1:2]+0.05,box$stats[1,3:4]-0.05,box$stats[5,5]+0.05,box$stats[1,6:7]-0.05,box$stats[5,8:10]+0.05),
     c("G1","G2","G3","G4","G5","G6","G7","G8","G9","G10"), 
     cex=1,col="black")
axis(1)

coef=res$coefficients
abline(coef[1],coef[2],col="cyan",lwd=3)
abline(coef[1],coef[2]+coef[3],col="gold",lwd=3)
bxp(box,at=epsilons,pars=list(xaxt='n',boxfill=c("cyan","gold","cyan","cyan","gold","gold","gold","cyan","cyan","cyan"),boxwex=diff(range(epsilons))/50,whisklty=3,staplewex=0.8),
    outline=FALSE,ylim=c(-0.55,0.85),
    cex.main=2,cex.axis=1,cex.lab=1,add=TRUE)

legend('topleft',legend=c('Inter','Intra'),col=c('gold','cyan'),lwd=3,border=FALSE)
dev.off()















###############################
# REGRESSION MODEL WITH ONLY  #
#      INTRAMOLECULAR G4      #
###############################

# ONLY INTRAMOLECULAR
select=c(1,3,4,8,9,10)

# Fit final regression model for Tm (after selecting relevant predictors)
res=lm(log(intra_G+0.01)~intra_TMs)
summary(res)
plot(res)

pdf('termostability_Tm_intra.pdf',width=6,height=6)
par(mar=c(5,4,2,2)+0.1)
box=boxplot(log_mean_G[select],plot=FALSE)
box$stats[c(1,5),]=Reduce(cbind,lapply(log_mean_G[select],quantile,probs=c(0.05,0.95),na.rm=TRUE))
bxp(box,at=TMs[select],pars=list(xaxt='n',boxfill=c("cyan","gold","cyan","cyan","gold","gold","gold","cyan","cyan","cyan")[select],boxwex=diff(range(TMs[select]))/50+0.025,whisklty=3,staplewex=0.8),
    main="",ylab="Log mean IPD",xlab=expression(T[m]*" [°C]"),outline=FALSE,ylim=c(-0.45,0.85),
    cex.main=2,cex.axis=1,cex.lab=1)
text(TMs[select],c(box$stats[5,1:2]+0.05,box$stats[1,3]-0.05,box$stats[5,4:6]+0.05),
     c("G1","G2","G3","G4","G5","G6","G7","G8","G9","G10")[select], 
     cex=1,col="black")
axis(1)

coef=res$coefficients
abline(coef[1],coef[2],col="cyan",lwd=3)
bxp(box,at=TMs[select],pars=list(xaxt='n',boxfill=c("cyan","gold","cyan","cyan","gold","gold","gold","cyan","cyan","cyan")[select],boxwex=diff(range(TMs[select]))/50+0.025,whisklty=3,staplewex=0.8),
    outline=FALSE,ylim=c(-0.45,0.85),axes=FALSE,
    cex.main=2,cex.axis=1,cex.lab=1,add=TRUE)
dev.off()




# Fit final regression model for delta-epsilon (after selecting relevant predictors)
res=lm(log(intra_G+0.01)~intra_epsilons)
summary(res)
plot(res)

pdf('termostability_deltaepsilon_intra.pdf',width=6,height=6)
par(mar=c(5,4,2,2)+0.1)
box=boxplot(log_mean_G[select],plot=FALSE)
box$stats[c(1,5),]=Reduce(cbind,lapply(log_mean_G[select],quantile,probs=c(0.05,0.95),na.rm=TRUE))
bxp(box,at=epsilons[select],pars=list(xaxt='n',boxfill=c("cyan","gold","cyan","cyan","gold","gold","gold","cyan","cyan","cyan")[select],boxwex=diff(range(epsilons[select]))/50,whisklty=3,staplewex=0.8),
    main="",ylab="Log mean IPD",xlab=expression(Delta*epsilon*" ["*M^{-1}*cm^{-1}*"]"),outline=FALSE,ylim=c(-0.45,0.85),
    cex.main=2,cex.axis=1,cex.lab=1)
text(epsilons[select],c(box$stats[5,1]+0.05,box$stats[1,2:3]-0.05,box$stats[5,4:6]+0.05),
     c("G1","G2","G3","G4","G5","G6","G7","G8","G9","G10")[select], 
     cex=1,col="black")
axis(1)

coef=res$coefficients
abline(coef[1],coef[2],col="cyan",lwd=3)
bxp(box,at=epsilons[select],pars=list(xaxt='n',boxfill=c("cyan","gold","cyan","cyan","gold","gold","gold","cyan","cyan","cyan")[select],boxwex=diff(range(epsilons[select]))/50,whisklty=3,staplewex=0.8),
    outline=FALSE,ylim=c(-0.45,0.85),axes=FALSE,
    cex.main=2,cex.axis=1,cex.lab=1,add=TRUE)
dev.off()


