PDvsPbErrorsvsDiv <- read.table('GQuadPlus_Compo_IPD_PbError_Div_1kG_residuals_FILTERED_CONTROLS_phyloP', header = TRUE)

q1 <- quantile(IPDvsPbErrorsvsDiv$DIVERGENCESNP,probs = c(0.05))
q2 <- quantile(IPDvsPbErrorsvsDiv$DIVERGENCESNP,probs = c(0.25))
q3 <- quantile(IPDvsPbErrorsvsDiv$DIVERGENCESNP,probs = c(0.75))
q4 <- quantile(IPDvsPbErrorsvsDiv$DIVERGENCESNP,probs = c(0.93))
q5 <- quantile(IPDvsPbErrorsvsDiv$DIVERGENCESNP,probs = c(0.95))
q6 <- quantile(IPDvsPbErrorsvsDiv$DIVERGENCESNP,probs = c(0.96))
q7 <- quantile(IPDvsPbErrorsvsDiv$DIVERGENCESNP,probs = c(0.97))
q8 <- quantile(IPDvsPbErrorsvsDiv$DIVERGENCESNP,probs = c(0.98))


q1_Ctrl<- quantile(IPDvsPbErrorsvsDiv_Ctrl$DIVERGENCESNP,probs = c(0.05,0.95))[1]
q2_Ctrl <- quantile(IPDvsPbErrorsvsDiv_Ctrl$DIVERGENCESNP,probs = c(0.05,0.95))[2]

df1 <- IPDvsPbErrorsvsDiv[which(IPDvsPbErrorsvsDiv$DIVERGENCESNP <= q1),]
df2 <- IPDvsPbErrorsvsDiv[which(IPDvsPbErrorsvsDiv$DIVERGENCESNP <= q2),]
df3 <- IPDvsPbErrorsvsDiv[which(IPDvsPbErrorsvsDiv$DIVERGENCESNP >= q3),]
df4 <- IPDvsPbErrorsvsDiv[which(IPDvsPbErrorsvsDiv$DIVERGENCESNP >= q4),]
df5 <- IPDvsPbErrorsvsDiv[which(IPDvsPbErrorsvsDiv$DIVERGENCESNP >= q5),]
df6 <- IPDvsPbErrorsvsDiv[which(IPDvsPbErrorsvsDiv$DIVERGENCESNP >= q6),]
df7 <- IPDvsPbErrorsvsDiv[which(IPDvsPbErrorsvsDiv$DIVERGENCESNP >= q7),]
df8 <- IPDvsPbErrorsvsDiv[which(IPDvsPbErrorsvsDiv$DIVERGENCESNP >= q8),]


df1_Ctrl <- IPDvsPbErrorsvsDiv_Ctrl[which(IPDvsPbErrorsvsDiv_Ctrl$DIVERGENCESNP <= q1),]
df2_Ctrl <- IPDvsPbErrorsvsDiv_Ctrl[which(IPDvsPbErrorsvsDiv_Ctrl$DIVERGENCESNP >= q3),]


median(df1$PACBIOSNP)
median(df2$PACBIOSNP)
median(df3$PACBIOSNP)
median(df4$PACBIOSNP)
median(df5$PACBIOSNP)
median(df6$PACBIOSNP)
median(df7$PACBIOSNP)
median(df8$PACBIOSNP)

median(df1$LOGMEANIPD_RES, na.rm = TRUE)
median(df2$LOGMEANIPD_RES, na.rm = TRUE)
median(df3$LOGMEANIPD_RES, na.rm = TRUE)
median(df4$LOGMEANIPD_RES, na.rm = TRUE)
median(df5$LOGMEANIPD_RES, na.rm = TRUE)
median(df6$LOGMEANIPD_RES, na.rm = TRUE)
median(df7$LOGMEANIPD_RES, na.rm = TRUE)
median(df8$LOGMEANIPD_RES, na.rm = TRUE)

median(df1$MeanPhyloP, na.rm = TRUE)
median(df2$MeanPhyloP, na.rm = TRUE)

boxplot(df1$MeanPhyloP,df2$MeanPhyloP, outline = FALSE)
boxplot(df1_Ctrl$MeanPhyloP,df2_Ctrl$MeanPhyloP, outline = FALSE)

boxplot(df1$MeanPhyloP,df1_Ctrl$MeanPhyloP,df2$MeanPhyloP,df2_Ctrl$MeanPhyloP, outline = FALSE, names = c('G4+', 'Ctrl', 'G4+', 'Ctrl'), col = c('red', 'red', 'blue', 'blue'), main = 'Mean PhyloP in Motifs' )
legend(3.5, 1.5, legend=c("Low Divergence", "High Divergence"),
       col=c("red", "blue"), lty=1:1, cex=0.8)

p_pacbio=rep(NA,1000)
p_ipd=rep(NA,1000)


for(i in 1:1000){
  df1_sub=df1[sample(1:nrow(df1),nrow(df3)),]
  p_pacbio[i]=t.test(log(df1_sub$PACBIOSNP+0.01),log(df3$PACBIOSNP+0.01))$p.value
  p_ipd[i]=t.test(df1_sub$LOGMEANIPD_RES,df3$LOGMEANIPD_RES)$p.value
}
hist(p_pacbio, breaks = 50)
hist(p_ipd, breaks = 50)
median(p_pacbio)
median(p_ipd)

p_pacbio=rep(NA,1000)
p_ipd=rep(NA,1000)

for(i in 1:1000){
  df1_sub=df1[sample(1:nrow(df1),nrow(df4)),]
  p_pacbio[i]=t.test(log(df1_sub$PACBIOSNP+0.01),log(df4$PACBIOSNP+0.01))$p.value
  p_ipd[i]=t.test(df1_sub$LOGMEANIPD_RES,df4$LOGMEANIPD_RES)$p.value
}

median(p_pacbio)
median(p_ipd)

p_pacbio=rep(NA,1000)
p_ipd=rep(NA,1000)

for(i in 1:1000){
  df1_sub=df1[sample(1:nrow(df1),nrow(df5)),]
  p_pacbio[i]=t.test(log(df1_sub$PACBIOSNP+0.01),log(df5$PACBIOSNP+0.01))$p.value
  p_ipd[i]=t.test(df1_sub$LOGMEANIPD_RES,df5$LOGMEANIPD_RES)$p.value
}
median(p_pacbio)
median(p_ipd)

p_pacbio=rep(NA,1000)
p_ipd=rep(NA,1000)

for(i in 1:1000){
  df1_sub=df1[sample(1:nrow(df1),nrow(df6)),]
  p_pacbio[i]=t.test(log(df1_sub$PACBIOSNP+0.01),log(df6$PACBIOSNP+0.01))$p.value
  p_ipd[i]=t.test(df1_sub$LOGMEANIPD_RES,df6$LOGMEANIPD_RES)$p.value
}
median(p_pacbio)
median(p_ipd)

p_pacbio=rep(NA,1000)
p_ipd=rep(NA,1000)

for(i in 1:1000){
  df1_sub=df1[sample(1:nrow(df1),nrow(df7)),]
  p_pacbio[i]=t.test(log(df1_sub$PACBIOSNP+0.01),log(df7$PACBIOSNP+0.01))$p.value
  p_ipd[i]=t.test(df1_sub$LOGMEANIPD_RES,df7$LOGMEANIPD_RES)$p.value
}
median(p_pacbio)
median(p_ipd)

p_pacbio=rep(NA,1000)
p_ipd=rep(NA,1000)

for(i in 1:1000){
  df1_sub=df1[sample(1:nrow(df1),nrow(df8)),]
  p_pacbio[i]=t.test(log(df1_sub$PACBIOSNP+0.01),log(df8$PACBIOSNP+0.01))$p.value
  p_ipd[i]=t.test(df1_sub$LOGMEANIPD_RES,df8$LOGMEANIPD_RES)$p.value
}
median(p_pacbio)
median(p_ipd)

plot(density(log(df1$PACBIOSNP+0.01),na.rm=TRUE), xlim=c(-4.5,-2), ylim=c(0,2), main = 'Divergence', xlab = "log error rate", lwd=3, axes=FALSE, cex.lab=1.4)
axis(side = 1, lwd = 3,)
axis(side = 2, lwd = 3)
lines(density(log(df3$PACBIOSNP+0.01),na.rm=TRUE),col='purple', lwd=3)
lines(density(log(df4$PACBIOSNP+0.01),na.rm=TRUE),col='red', lwd=3)
lines(density(log(df5$PACBIOSNP+0.01),na.rm=TRUE),col='pink', lwd=3)
lines(density(log(df6$PACBIOSNP+0.01),na.rm=TRUE),col='yellow', lwd=3)
abline(v=mean(log(df1$PACBIOSNP+0.01),na.rm = TRUE), lwd=2)
abline(v=mean(log(df3$PACBIOSNP+0.01),na.rm = TRUE),col='purple', lwd=2)
abline(v=mean(log(df4$PACBIOSNP+0.01),na.rm = TRUE),col='red', lwd=2)
abline(v=mean(log(df5$PACBIOSNP+0.01),na.rm = TRUE),col='pink', lwd=2)
abline(v=mean(log(df6$PACBIOSNP+0.01),na.rm = TRUE),col='yellow', lwd=2)


legend(-3, 1.5, legend=c("Top 3%","Top 5%", "Top 7%", "Top 25%", "Bottom 5%"),
       col=c("yellow","pink","red","purple","black"), lty=1, cex=0.8, lwd = 3)

plot(density(df1$LOGMEANIPD_RES,na.rm=TRUE), xlim=c(-1,2), ylim=c(0,2), main = 'Divergence', xlab = "log residuals IPD", lwd=3, axes=FALSE, cex.lab=1.4)
axis(side = 1, lwd = 3,)
axis(side = 2, lwd = 3)
lines(density(df3$LOGMEANIPD_RES,na.rm=TRUE),col='purple', lwd=3)
lines(density(df4$LOGMEANIPD_RES,na.rm=TRUE),col='red', lwd=3)
lines(density(df5$LOGMEANIPD_RES,na.rm=TRUE),col='pink', lwd=3)
lines(density(df6$LOGMEANIPD_RES,na.rm=TRUE),col='yellow', lwd=3)
abline(v=mean(df1$LOGMEANIPD_RES,na.rm = TRUE), lwd=2)
abline(v=mean(df3$LOGMEANIPD_RES,na.rm = TRUE),col='purple', lwd=2)
abline(v=mean(df4$LOGMEANIPD_RES,na.rm = TRUE),col='red', lwd=2)
abline(v=mean(df5$LOGMEANIPD_RES,na.rm = TRUE),col='pink', lwd=2)
abline(v=mean(df6$LOGMEANIPD_RES,na.rm = TRUE),col='yellow', lwd=2)


legend(1, 1.5, legend=c("Top 3%","Top 5%", "Top 7%", "Top 25%", "Bottom 5%"),
       col=c("yellow","pink","red","purple","black"), lty=1, cex=0.8, lwd = 3)


q1 <- quantile(IPDvsPbErrorsvsDiv$DIVERSITYSNP,probs = c(0.05))
q2 <- quantile(IPDvsPbErrorsvsDiv$DIVERSITYSNP,probs = c(0.25))
q3 <- quantile(IPDvsPbErrorsvsDiv$DIVERSITYSNP,probs = c(0.75))
q4 <- quantile(IPDvsPbErrorsvsDiv$DIVERSITYSNP,probs = c(0.93))
q5 <- quantile(IPDvsPbErrorsvsDiv$DIVERSITYSNP,probs = c(0.95))
q6 <- quantile(IPDvsPbErrorsvsDiv$DIVERSITYSNP,probs = c(0.97))


df1 <- IPDvsPbErrorsvsDiv[which(IPDvsPbErrorsvsDiv$DIVERSITYSNP <= q1),]
df2 <- IPDvsPbErrorsvsDiv[which(IPDvsPbErrorsvsDiv$DIVERSITYSNP <= q2),]
df3 <- IPDvsPbErrorsvsDiv[which(IPDvsPbErrorsvsDiv$DIVERSITYSNP >= q3),]
df4 <- IPDvsPbErrorsvsDiv[which(IPDvsPbErrorsvsDiv$DIVERSITYSNP >= q4),]
df5 <- IPDvsPbErrorsvsDiv[which(IPDvsPbErrorsvsDiv$DIVERSITYSNP >= q5),]
df6 <- IPDvsPbErrorsvsDiv[which(IPDvsPbErrorsvsDiv$DIVERSITYSNP >= q6),]


median(df1$PACBIOSNP)
median(df2$PACBIOSNP)
median(df3$PACBIOSNP)
median(df4$PACBIOSNP)
median(df5$PACBIOSNP)
median(df6$PACBIOSNP)
median(df7$PACBIOSNP)
median(df8$PACBIOSNP)

median(df1$LOGMEANIPD_RES, na.rm = TRUE)
median(df2$LOGMEANIPD_RES, na.rm = TRUE)
median(df3$LOGMEANIPD_RES, na.rm = TRUE)
median(df4$LOGMEANIPD_RES, na.rm = TRUE)
median(df5$LOGMEANIPD_RES, na.rm = TRUE)
median(df6$LOGMEANIPD_RES, na.rm = TRUE)
median(df7$LOGMEANIPD_RES, na.rm = TRUE)
median(df8$LOGMEANIPD_RES, na.rm = TRUE)

median(df1$MeanPhyloP, na.rm = TRUE)
median(df2$MeanPhyloP, na.rm = TRUE)

p_pacbio=rep(NA,1000)
p_ipd=rep(NA,1000)

#for(i in 1:1000){
#  df1_sub=df1[sample(1:nrow(df1),nrow(df3)),]
#p_pacbio[i]=t.test(log(df1_sub$PACBIOSNP+0.01),log(df3$PACBIOSNP+0.01))$p.value
#p_ipd[i]=t.test(df1_sub$LOGMEANIPD_RES,df3$LOGMEANIPD_RES)$p.value
#}

for(i in 1:1000){
  df3_sub=df3[sample(1:nrow(df3),nrow(df1)),]
  p_pacbio[i]=t.test(log(df3_sub$PACBIOSNP+0.01),log(df1$PACBIOSNP+0.01))$p.value
  p_ipd[i]=t.test(df3_sub$LOGMEANIPD_RES,df1$LOGMEANIPD_RES)$p.value
}

hist(p_pacbio, breaks = 50)
hist(p_ipd, breaks = 50)
median(p_pacbio)
median(p_ipd)

p_pacbio=rep(NA,1000)
p_ipd=rep(NA,1000)

for(i in 1:1000){
  df4_sub=df4[sample(1:nrow(df4),nrow(df1)),]
  p_pacbio[i]=t.test(log(df4_sub$PACBIOSNP+0.01),log(df1$PACBIOSNP+0.01))$p.value
  p_ipd[i]=t.test(df4_sub$LOGMEANIPD_RES,df1$LOGMEANIPD_RES)$p.value
}

median(p_pacbio)
median(p_ipd)

p_pacbio=rep(NA,1000)
p_ipd=rep(NA,1000)

for(i in 1:1000){
  df5_sub=df5[sample(1:nrow(df5),nrow(df1)),]
  p_pacbio[i]=t.test(log(df5_sub$PACBIOSNP+0.01),log(df1$PACBIOSNP+0.01))$p.value
  p_ipd[i]=t.test(df5_sub$LOGMEANIPD_RES,df1$LOGMEANIPD_RES)$p.value
}

median(p_pacbio)
median(p_ipd)

p_pacbio=rep(NA,1000)
p_ipd=rep(NA,1000)

for(i in 1:1000){
  df1_sub=df1[sample(1:nrow(df1),nrow(df6)),]
  p_pacbio[i]=t.test(log(df1_sub$PACBIOSNP+0.01),log(df6$PACBIOSNP+0.01))$p.value
  p_ipd[i]=t.test(df1_sub$LOGMEANIPD_RES,df6$LOGMEANIPD_RES)$p.value
}
median(p_pacbio)
median(p_ipd)

p_pacbio=rep(NA,1000)
p_ipd=rep(NA,1000)

for(i in 1:1000){
  df1_sub=df1[sample(1:nrow(df1),nrow(df7)),]
  p_pacbio[i]=t.test(log(df1_sub$PACBIOSNP+0.01),log(df7$PACBIOSNP+0.01))$p.value
  p_ipd[i]=t.test(df1_sub$LOGMEANIPD_RES,df7$LOGMEANIPD_RES)$p.value
}
median(p_pacbio)
median(p_ipd)

p_pacbio=rep(NA,1000)
p_ipd=rep(NA,1000)

for(i in 1:1000){
  df1_sub=df1[sample(1:nrow(df1),nrow(df8)),]
  p_pacbio[i]=t.test(log(df1_sub$PACBIOSNP+0.01),log(df8$PACBIOSNP+0.01))$p.value
  p_ipd[i]=t.test(df1_sub$LOGMEANIPD_RES,df8$LOGMEANIPD_RES)$p.value
}
median(p_pacbio)
median(p_ipd)


plot(density(log(df1$PACBIOSNP+0.01),na.rm=TRUE), xlim=c(-4.5,-2), ylim=c(0,2), main = 'Diversity', xlab = "log error rate", lwd=3, axes=FALSE, cex.lab=1.4)
axis(side = 1, lwd = 3,)
axis(side = 2, lwd = 3)
lines(density(log(df3$PACBIOSNP+0.01),na.rm=TRUE),col='purple', lwd=3)
lines(density(log(df4$PACBIOSNP+0.01),na.rm=TRUE),col='red', lwd=3)
lines(density(log(df5$PACBIOSNP+0.01),na.rm=TRUE),col='pink', lwd=3)
lines(density(log(df6$PACBIOSNP+0.01),na.rm=TRUE),col='yellow', lwd=3)
abline(v=mean(log(df1$PACBIOSNP+0.01),na.rm = TRUE), lwd=2)
abline(v=mean(log(df3$PACBIOSNP+0.01),na.rm = TRUE),col='purple', lwd=2)
abline(v=mean(log(df4$PACBIOSNP+0.01),na.rm = TRUE),col='red', lwd=2)
abline(v=mean(log(df5$PACBIOSNP+0.01),na.rm = TRUE),col='pink', lwd=2)
abline(v=mean(log(df6$PACBIOSNP+0.01),na.rm = TRUE),col='yellow', lwd=2)



legend(-3, 1.5, legend=c("Top 3%","Top 5%", "Top 7%", "Top 25%", "Bottom 5%"),
       col=c("yellow","pink","red","purple","black"), lty=1, cex=0.8, lwd = 3)

plot(density(df1$LOGMEANIPD_RES,na.rm=TRUE), xlim=c(-1,2), ylim=c(0,2), main = 'Diversity', xlab = "log residuals IPD", lwd=3, axes=FALSE, cex.lab=1.4)
axis(side = 1, lwd = 3,)
axis(side = 2, lwd = 3)
lines(density(df3$LOGMEANIPD_RES,na.rm=TRUE),col='purple', lwd=3)
lines(density(df4$LOGMEANIPD_RES,na.rm=TRUE),col='red', lwd=3)
lines(density(df5$LOGMEANIPD_RES,na.rm=TRUE),col='pink', lwd=3)
lines(density(df6$LOGMEANIPD_RES,na.rm=TRUE),col='yellow', lwd=3)
abline(v=mean(df1$LOGMEANIPD_RES,na.rm = TRUE), lwd=2)
abline(v=mean(df3$LOGMEANIPD_RES,na.rm = TRUE),col='purple', lwd=2)
abline(v=mean(df4$LOGMEANIPD_RES,na.rm = TRUE),col='red', lwd=2)
abline(v=mean(df5$LOGMEANIPD_RES,na.rm = TRUE),col='pink', lwd=2)
abline(v=mean(df6$LOGMEANIPD_RES,na.rm = TRUE),col='yellow', lwd=2)


legend(1, 1.5, legend=c("Top 3%","Top 5%", "Top 7%", "Top 25%", "Bottom 5%"),
       col=c("yellow","pink","red","purple","black"), lty=1, cex=0.8, lwd = 3)


hist(IPDvsPbErrorsvsDiv$DIVERSITYSNP, breaks = 100)
temp <- hist(IPDvsPbErrorsvsDiv$DIVERSITYSNP, plot = FALSE, breaks = 20) #get histogram data
plot(x = temp$mids, y = log(temp$counts), type = "h", lwd = 10, main = "Diversity", xlab = "Diversity in motifs", ylab = 'log(counts)', axes=FALSE,cex.lab=1.4, xlim = c(0,0.12), ylim=c(0,10))
axis(side=1,lwd=3)
axis(side=2,lwd=3)



hist(IPDvsPbErrorsvsDiv$DIVERGENCESNP, breaks = 100)
temp <- hist(IPDvsPbErrorsvsDiv$DIVERGENCESNP, plot = FALSE, breaks = 20) #get histogram data
plot(x = temp$mids, y = log(temp$counts), type = "h", lwd =10, main = "Divergence", xlab = "Divergence in motifs", ylab = "log(counts)", axes=FALSE,cex.lab=1.4, xlim = c(0,0.35), ylim=c(0,10))
axis(side=1,lwd=3)
axis(side=2,lwd=3)
