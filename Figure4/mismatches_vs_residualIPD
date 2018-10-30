IPDvsPbErrorsvsDivG4Plus <- read.table('GQuadPlus_Compo_IPD_PbError_Div_1kG_residuals_FILTERED_CONTROLS', header = TRUE)
IPDisNotNA = which(!is.na(IPDvsPbErrorsvsDivG4Plus$MEANIPD))
IPDvsPbErrorsvsDivG4Plus <- IPDvsPbErrorsvsDivG4Plus[IPDisNotNA,]


IPDvsPbErrorsvsDivG4Minus <- read.table('GQuadMinus_Compo_IPD_PbError_Div_1kG_residuals_FILTERED_CONTROLS', header = TRUE)
IPDisNotNA = which(!is.na(IPDvsPbErrorsvsDivG4Minus$MEANIPD))
IPDvsPbErrorsvsDivG4Minus <- IPDvsPbErrorsvsDivG4Minus[IPDisNotNA,]


IPDvsErrorsControl <- read.table('Empty_Compo_IPD_PbError_Div_1kG_residuals_FILTERED_CONTROLS', header = TRUE)
IPDisNotNA = which(!is.na(IPDvsErrorsControl$MEANIPD))
IPDvsErrorsControl <- IPDvsErrorsControl[IPDisNotNA,]

set.seed(44)

IPDvsErrorsControlSample <- IPDvsErrorsControl[sample(nrow(IPDvsErrorsControl), nrow(IPDvsPbErrorsvsDivG4Plus) + nrow(IPDvsPbErrorsvsDivG4Minus)),]

IPDvsErrors <- rbind(IPDvsPbErrorsvsDivG4Plus,IPDvsPbErrorsvsDivG4Minus,IPDvsErrorsControlSample)

IPDvsErrors <- cbind(IPDvsErrors,rep(1:0,c(nrow(IPDvsPbErrorsvsDivG4Plus),nrow(IPDvsPbErrorsvsDivG4Minus)+nrow(IPDvsErrorsControlSample))))

IPDvsErrors <- cbind(IPDvsErrors,
                    c(rep(0:1,c(nrow(IPDvsPbErrorsvsDivG4Plus),nrow(IPDvsPbErrorsvsDivG4Minus))),rep(0:1,c(nrow(IPDvsErrorsControlSample),0)))
)

colnames(IPDvsErrors) <- c('CHROM','START','END','MEANIPD','PACBIOSNP','DIVERGENCESNP','DIVERSITYSNP','LOGMEANIPD','LOGMEANIPD_PRED','LOGMEANIPD_RES','TestG4P','TestG4M')
IPDvsErrors=IPDvsErrors[sample(nrow(IPDvsErrors)),]

fit = lm(log(PACBIOSNP+0.01)~LOGMEANIPD_RES+TestG4P+TestG4P*LOGMEANIPD_RES+TestG4M+TestG4M*LOGMEANIPD_RES, data = IPDvsErrors)
summary(fit)


### Figure 4
x11(15,7)
par(cex = 1.2, cex.axis= 1.2, cex.lab = 1.2)
alpha=0.15
colors=c('gray50','dodgerblue','tomato')
palette(adjustcolor(colors, alpha.f = alpha))
plot(log(PACBIOSNP+0.01)~LOGMEANIPD_RES,  pch = 20, data = IPDvsErrors,col=IPDvsErrors$TestG4M*1+IPDvsErrors$TestG4P*2+1,main='', xlim = c(-0.5,2.5), ylim = c(-4.6,-2.5), xlab="Residual log(mean IPD)", ylab="Log SMRT Mismatches", cex=1.4, axes=FALSE)
abline(fit$coefficients[1],fit$coefficients[2],col=colors[1],lwd=3)
abline(fit$coefficients[1]+fit$coefficients[3],fit$coefficients[2]+fit$coefficients[5],col=colors[3],lwd=3)
abline(fit$coefficients[1]+fit$coefficients[4],fit$coefficients[2]+fit$coefficients[6],col=colors[2],lwd=3)
axis(side = 1, lwd = 3)
axis(side = 2, lwd = 3)

legend('bottomright', legend=c("G4+", "G4-", "Motif-free"),
       col=colors[3:1], lty=1, cex=1, lwd = 5, bty = "n", y.intersp = 0.8)




### Figure 4 plus three separate plots for the three groups
x11(15,10)
par(cex = 1.2, cex.axis= 1.2, cex.lab = 1.2, mfrow=c(2,2))
plot(log(PACBIOSNP+0.01)~LOGMEANIPD_RES,  pch = 20, data = IPDvsErrors,col=IPDvsErrors$TestG4M*1+IPDvsErrors$TestG4P*2+1,main='', xlim = c(-0.5,2.5), ylim = c(-4.6,-2.5), xlab="Residual log(mean IPD)", ylab="Log SMRT Mismatches", cex=1.4, axes=FALSE)
abline(fit$coefficients[1],fit$coefficients[2],col=colors[1],lwd=3)
abline(fit$coefficients[1]+fit$coefficients[3],fit$coefficients[2]+fit$coefficients[5],col=colors[3],lwd=3)
abline(fit$coefficients[1]+fit$coefficients[4],fit$coefficients[2]+fit$coefficients[6],col=colors[2],lwd=3)
axis(side = 1, lwd = 3)
axis(side = 2, lwd = 3)
legend('bottomright', legend=c("G4+", "G4-", "Motif-free"),
       col=colors[3:1], lty=1, cex=1, lwd = 5, bty = "n", y.intersp = 0.8)

plot(log(PACBIOSNP+0.01)[(IPDvsErrors$TestG4M+IPDvsErrors$TestG4P)==0]~LOGMEANIPD_RES[(IPDvsErrors$TestG4M+IPDvsErrors$TestG4P)==0],  
     pch = 20, data = IPDvsErrors,col=1,main='', xlim = c(-0.5,2.5), ylim = c(-4.6,-2.5), xlab="Residual log(mean IPD)", ylab="Log SMRT Mismatches", cex=1.4, axes=FALSE)
abline(fit$coefficients[1],fit$coefficients[2],col=colors[1],lwd=3)
abline(fit$coefficients[1]+fit$coefficients[3],fit$coefficients[2]+fit$coefficients[5],col=colors[3],lwd=3)
abline(fit$coefficients[1]+fit$coefficients[4],fit$coefficients[2]+fit$coefficients[6],col=colors[2],lwd=3)
axis(side = 1, lwd = 3)
axis(side = 2, lwd = 3)
legend('bottomright', legend=c("G4+", "G4-", "Motif-free"),
       col=colors[3:1], lty=1, cex=1, lwd = 5, bty = "n", y.intersp = 0.8)

plot(log(PACBIOSNP+0.01)[IPDvsErrors$TestG4M==1]~LOGMEANIPD_RES[IPDvsErrors$TestG4M==1],  
     pch = 20, data = IPDvsErrors,col=2,main='', xlim = c(-0.5,2.5), ylim = c(-4.6,-2.5), xlab="Residual log(mean IPD)", ylab="Log SMRT Mismatches", cex=1.4, axes=FALSE)
abline(fit$coefficients[1],fit$coefficients[2],col=colors[1],lwd=3)
abline(fit$coefficients[1]+fit$coefficients[3],fit$coefficients[2]+fit$coefficients[5],col=colors[3],lwd=3)
abline(fit$coefficients[1]+fit$coefficients[4],fit$coefficients[2]+fit$coefficients[6],col=colors[2],lwd=3)
axis(side = 1, lwd = 3)
axis(side = 2, lwd = 3)
legend('bottomright', legend=c("G4+", "G4-", "Motif-free"),
       col=colors[3:1], lty=1, cex=1, lwd = 5, bty = "n", y.intersp = 0.8)

plot(log(PACBIOSNP+0.01)[IPDvsErrors$TestG4P==1]~LOGMEANIPD_RES[IPDvsErrors$TestG4P==1],  
     pch = 20, data = IPDvsErrors,col=3,main='', xlim = c(-0.5,2.5), ylim = c(-4.6,-2.5), xlab="Residual log(mean IPD)", ylab="Log SMRT Mismatches", cex=1.4, axes=FALSE)
abline(fit$coefficients[1],fit$coefficients[2],col=colors[1],lwd=3)
abline(fit$coefficients[1]+fit$coefficients[3],fit$coefficients[2]+fit$coefficients[5],col=colors[3],lwd=3)
abline(fit$coefficients[1]+fit$coefficients[4],fit$coefficients[2]+fit$coefficients[6],col=colors[2],lwd=3)
axis(side = 1, lwd = 3)
axis(side = 2, lwd = 3)
legend('bottomright', legend=c("G4+", "G4-", "Motif-free"),
       col=colors[3:1], lty=1, cex=1, lwd = 5, bty = "n", y.intersp = 0.8)

