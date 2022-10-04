rm(list=ls())
library(scales) # For percent function
library(caret)
library(data.table)

#setwd("L:/Lab/NCCT_ExpoCast/ExpoCast2022/Dawson_PFAS_HALFLIFE/PFAS_HL_QSAR_2021")
readsuff="JFW100322-noLogD" 

load(file=paste("RData/ClassificationModel_3Bin_EndoSimDisc_HLH_", readsuff, ".RData",sep=""))
load(file=paste("RData/ClassificationModel_4Bin_EndoSimDisc_HLH_", readsuff, ".RData",sep=""))
load(file=paste("RData/ClassificationModel_5Bin_EndoSimDisc_HLH_", readsuff, ".RData",sep=""))
load(file=paste("RData/ClassificationModel_PostRFE_HLHBin4_", readsuff, ".RData", sep=""))
load(file=paste("RData/ClassificationModel_YRand_HLHBin4_", readsuff, ".RData", sep=""))
load(file=paste0("RData/Tox21_AllMods_ADindicate_", readsuff,".RData"))
           
# Abstract
print(paste("The classification model had an average accuracy of, ",
  percent(classmod4$results[["Accuracy"]],accuracy=0.1),
  " Â± ",
  percent(classmod4$results[["AccuracySD"]],accuracy=0.1),
  " across species and chemicals.",
  sep=""))

print(paste("In contrast, y-randomized training data had an average accuracy of ",
  percent(classmod4YROverall$results[["Accuracy"]],accuracy=0.1),
  " ± ",
  percent(classmod4YROverall$results[["AccuracySD"]],accuracy=0.1),
  ".",
  sep="")) 


# Y-Randomization:

print(paste0("A model using t½ values randomized across all species-by-PFAS combinations had low predictive value (accuracy of ",
  percent(classmod4YROverall$results[["Accuracy"]],accuracy=0.1),
  " ± ",
  percent(classmod4YROverall$results[["AccuracySD"]],accuracy=0.1),
  "). However, the models fit with t½ values randomized within species (accuracy of ",
  percent(classmod4YRspecies$results[["Accuracy"]],accuracy=0.1),
  " ± ",
  percent(classmod4YRspecies$results[["AccuracySD"]],accuracy=0.1),
  ") and within chemicals (accuracy of ",
  percent(classmod4YRchem$results[["Accuracy"]],accuracy=0.1),
  " ± ",
  percent(classmod4YRchem$results[["AccuracySD"]],accuracy=0.1),
  ") did show that some variation t½ is accounted for by differences at the species and chemical level."))






# In thalf model domain:
in.domain <- subset(DPFASwAD,ClassModDomain==1&Species%in%c("Mouse","Rat","Monkey","Human"))
humanAD.M <- DPFASwAD[DPFASwAD$Species=="Human" & 
  DPFASwAD$Sex=="Male" & 
  DPFASwAD$DosingAdj=="Other" &
  DPFASwAD$ClassModDomain==1,]
humanAD.M <-table(humanAD.M$ClassPredFull)/sum(table(humanAD.M$ClassModDomain))
humanAD.F <- DPFASwAD[DPFASwAD$Species=="Human" & 
  DPFASwAD$Sex=="Female" & 
  DPFASwAD$DosingAdj=="Other" &
  DPFASwAD$ClassModDomain==1,]
humanAD.F <-table(humanAD.F$ClassPredFull)/sum(table(humanAD.F$ClassModDomain))
print(paste0("For these ",
  length(unique(in.domain$DTXSID)),
  " chemicals, human t½ was predicted to be distributed such that ",
  percent(humanAD.F[["4"]]),
  " were classified in Bin 4, ",
  percent(humanAD.F[["3"]]),
  " were classified in Bin 3, and ",
  percent(humanAD.F[["2"]]),
  " were classified in Bin 2."))
   
# In AM domain:
in.AM.domain <- subset(DPFASwAD,AMAD==1&Species%in%c("Mouse","Rat","Monkey","Human"))

humanAMAD.M <- DPFASwAD[DPFASwAD$Species=="Human" & 
  DPFASwAD$Sex=="Male" & 
  DPFASwAD$DosingAdj=="Other" &
  DPFASwAD$AMAD==1,]
humanAMAD.M <-table(humanAMAD.M$ClassPredFull)/sum(table(humanAMAD.M$ClassModDomain))

humanAMAD.F <- DPFASwAD[DPFASwAD$Species=="Human" & 
  DPFASwAD$Sex=="Female" & 
  DPFASwAD$DosingAdj=="Other" &
  DPFASwAD$AMAD==1,]
humanAMAD.F <-table(humanAMAD.F$ClassPredFull)/sum(table(humanAMAD.F$ClassModDomain))

print(paste0("For humans a majority (",
  percent(humanAMAD.F[["2"]]),
  ") of this subset of chemicals were predicted to fall into Bin 2, followed by ",
  percent(humanAMAD.F[["4"]]),
  " in Bin 4 and ",
  percent(humanAMAD.F[["3"]]),
  " in Bin 3."))




# Results 3.1
print(paste("The models had cross-validated accuracies of ",
       percent(classmod3$results[["Accuracy"]],accuracy=0.1), ", ",
       percent(classmod4$results[["Accuracy"]],accuracy=0.1), ", and ",
       percent(classmod5$results[["Accuracy"]],accuracy=0.1), ", respectively.",sep=""))
print(paste("Cohens's Kappa{Tarald, 1989} was ",
       signif(classmod3$results[["Kappa"]],3), ", ",
       signif(classmod4$results[["Kappa"]],3), ", and ",
       signif(classmod5$results[["Kappa"]],3), ", respectively.",sep=""))       

load(paste0("RData/ClassMod_RFE_FullSet_", readsuff,".RData"))

# The accuracies and kapps below all are not as good as the 16 variable model:
rfeClass$results

# It looks like any of the models with >2 variables has accuracy > 0.8:
varImp(rfeClass)

# Importance:
varImp(classmod4,scale=FALSE)
varImp(classmod4)

# Occurence of different bins for humans:

print(paste0(length(unique(subset(DPFASwAD,Species=="Human" & ClassModDomain==1)$DTXSID)),
  " of ",
  length(unique(subset(DPFASwAD,Species=="Human")$DTXSID)),
  " chemicals are within the AD."))

 print(paste0("For humans a majority (",
  percent(humanAD.F[["4"]]),
  ") of this subset of chemicals were predicted to fall into Bin 4, followed by ",
  percent(humanAD.F[["2"]]),
  " in Bin 2, and ",
  percent(humanAD.F[["3"]]),
  " in Bin 1."))

print(paste0("we found that the majority (",
  percent(length(unique(subset(DPFASwAD,ClassModDomain==1)$DTXSID)) / length(unique(DPFASwAD$DTXSID))),
  ") of these chemicals fall into the domain of the model."))


print(paste0("Further restricting chemicals to those also within the ADs of the OPERA models serving as model predictors (All Model ADs / AMAD) potentially further reduces the list to ",
  length(unique(subset(DPFASwAD,AMAD==1)$DTXSID)),
  " of the ",
  length(unique(DPFASwAD$DTXSID)),
  "."))
  
