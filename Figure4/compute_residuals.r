require(IWTomics)
library(compositions)




##################################
### mononucleotide composition ###
##################################

files=c('GQuadPlus_Compo_IPD_PbError_Div_1kG','GQuadMinus_Compo_IPD_PbError_Div_1kG','APhased_Compo_IPD_PbError_Div_1kG',
        'Direct_Compo_IPD_PbError_Div_1kG','Inverted_Compo_IPD_PbError_Div_1kG','Mirror_Compo_IPD_PbError_Div_1kG',
        'ZDNA_Compo_IPD_PbError_Div_1kG','Empty_Compo_IPD_PbError_Div_1kG')
IDs=c('G4plus','G4minus','A_phased','Direct','Inverted','Mirror','Z_DNA','Control')
names=c('G4+','G4-','A phased rep','Direct rep','Inverted rep','Mirror rep','Z-DNA','Motif-free')


folder='compositions_residuals/'
composition=vector('list',length(files))
for(i in seq_along(files)){
  composition[[i]]=read.table(paste0(folder,files[i],'_composition'),sep='\t',header=FALSE)
  colnames(composition[[i]])=c('CHROM','START','END','MEANIPD','PACBIOSNP','DIVERGENCESNP','DIVERSITYSNP','N','A','T','G','C')
}
names(composition)=IDs
save(composition,file='composition_residuals.RData')


load('composition_residuals.RData')

# linear regression controls
# compositional regressors

load('regression_composition_results_log.RData')

# Isometric log ratio transform
meanIPD_composition=lapply(composition,function(region) cbind(data.frame(meanIPD=region$MEANIPD),ilr(region[,9:12])))

# regression with only filtered (non-repetitive regions) forward controls
lm_control_log_new=lm(log(meanIPD+0.01)~V1+V2+V3,meanIPD_composition$Control)
summary(lm_control_log_new)

# prediction features
meanIPD_log_predict_new=lapply(meanIPD_composition,
                               function(meanIPD){
                                 meanIPD_predicted=predict(lm_control_log_new,meanIPD)
                                 meanIPD$meanIPD=log(meanIPD$meanIPD+0.01)
                                 cbind(meanIPD,meanIPD_predicted)
                               })

# compute residuals
meanIPD_log_residual_new=lapply(meanIPD_log_predict_new,
                                function(meanIPD) meanIPD$meanIPD-meanIPD$meanIPD_predicted)

# write residuals in files
composition_residuals=vector('list',length(files))
for(i in seq_along(files)){
  composition_residuals[[i]]=cbind(composition[[i]][1:7],data.frame(LOGMEANIPD=meanIPD_log_predict_new[[i]]$meanIPD,
                                                                    LOGMEANIPD_PRED=meanIPD_log_predict_new[[i]]$meanIPD_predicted,
                                                                    LOGMEANIPD_RES=meanIPD_log_residual_new[[i]]))
  write.table(composition_residuals[[i]],paste0(files[i],'_residuals_FILTERED_CONTROLS'),sep='\t',row.names=FALSE,quote=FALSE)
}
names(composition_residuals)=IDs
