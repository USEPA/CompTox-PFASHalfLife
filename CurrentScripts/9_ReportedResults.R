rm(list=ls())
library(scales) # For percent function
library(caret)
library(data.table)

#setwd("L:/Lab/NCCT_ExpoCast/ExpoCast2022/Dawson_PFAS_HALFLIFE/PFAS_HL_QSAR_2021")
readsuff <- "JFW100322-noLogD" 
readsuff.data <- "JFW100322"

load(file=paste("RData/ClassificationModel_3Bin_EndoSimDisc_HLH_", readsuff, ".RData",sep=""))
load(file=paste("RData/ClassificationModel_4Bin_EndoSimDisc_HLH_", readsuff, ".RData",sep=""))
load(file=paste("RData/ClassificationModel_5Bin_EndoSimDisc_HLH_", readsuff, ".RData",sep=""))
load(file=paste("RData/ClassificationModel_PostRFE_HLHBin4_", readsuff, ".RData", sep=""))
load(file=paste("RData/ClassificationModel_YRand_HLHBin4_", readsuff, ".RData", sep=""))
load(file=paste0("RData/Tox21_AllMods_ADindicate_", readsuff,".RData"))
load(file=paste0("RData/PFAS_11Chemicals_QSARdataset_", readsuff.data,".RData"))
load(paste0("RData/ClassMod_RFE_FullSet_", readsuff,".RData"))
load(paste("RData/ClassyfierSubSets_",readsuff,".RData",sep=""))
  
# Calculate chemicals in domain for humans:
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

           
# Abstract
print("ABSTRACT")
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


print(paste0("When applied to USEPA’s largest list of ",
  length(unique(subset(DPFASwAD,Species=="Human")$DTXSID)),
  " PFAS, ",
  length(unique(subset(DPFASwAD,Species=="Human" & ClassModDomain==1)$DTXSID)),
  " compounds were estimated to be within domain of the model."))



  print(paste0("For these ",
               length(unique(subset(DPFASwAD,Species=="Human" & ClassModDomain==1)$DTXSID)),
               " chemicals, human t½ was predicted to be distributed such that ",
               percent(humanAD.F[["4"]]),
               " were classified in Bin 4, ",
               percent(humanAD.F[["3"]]),
               " were classified in Bin 3, and ",
               percent(humanAD.F[["2"]]),
               " were classified in Bin 2."))

# Results 3.1
print("RESULTS 3.1")
print("FIRST PARAGRAPH")
  
fit.all <- lm(log(HLH)~log(BW),data=pfasds) 
fit.pfoa <- lm(log(HLH)~log(BW),data=subset(pfasds,CASRN=="335-67-1"))
 
print(paste0("For PFOA we found the TK t½ scales only weakly across species with bodyweight (R2 = ",
      signif(summary(fit.pfoa)$adj.r.squared,2),
      "). This scaling was on average even less for the other chemicals in our data set (R2 = ",
      signif(summary(fit.all)$adj.r.squared,2),
      " overall)."))
         
print(paste("The models had cross-validated accuracies of ",
       percent(classmod3$results[["Accuracy"]],accuracy=0.1), ", ",
       percent(classmod4$results[["Accuracy"]],accuracy=0.1), ", and ",
       percent(classmod5$results[["Accuracy"]],accuracy=0.1), ", respectively.",sep=""))
print(paste("Cohens's Kappa {Tarald, 1989} was ",
       signif(classmod3$results[["Kappa"]],3), ", ",
       signif(classmod4$results[["Kappa"]],3), ", and ",
       signif(classmod5$results[["Kappa"]],3), ", respectively.",sep=""))                  
conf <- classmod4$finalModel$confusion
conf <- conf[,-ncol(conf)]
oob <- 1 - (sum(diag(conf))/sum(conf))
print(paste0("The four-bin model has an error rate of ",
       percent(oob),
       " and a no information rate of ",
             percent(max(apply(conf,1,sum))/100),
             ". The no information rate is an effective \"null hypothesis\" in which all chemicals were predicted to be in the most common bin."))     
          
print("SECOND PARAGRAPH")
print(paste0("A model using t½ values randomized across all species-by-PFAS combinations had low predictive value (accuracy of ",
  percent(classmod4YROverall$results[["Accuracy"]],accuracy=0.1),
  " ± ",
  percent(classmod4YROverall$results[["AccuracySD"]],accuracy=0.1),
  "). Y-randomization showed that some variation in t� is accounted for by differences at the species and chemical level. The models for with t� with training data randomized within species but not chemicals (that is, the chemicals were correct) had an accuracy of ",
  percent(classmod4YRspecies$results[["Accuracy"]],accuracy=0.1),
  " ± ",
  percent(classmod4YRspecies$results[["AccuracySD"]],accuracy=0.1),
  "). The models where training data chemical identities were randomized, but not species, had an accuracy of ",
  percent(classmod4YRchem$results[["Accuracy"]],accuracy=0.1),
  " ± ",
  percent(classmod4YRchem$results[["AccuracySD"]],accuracy=0.1),
  ")."))

# The accuracies and kappas below all are not as good as the 15 variable model:
rfeClass$results

print("TABLE THREE")
# Importance:
varImp(classmod4,scale=FALSE)
varImp(classmod4)
imp <- varImp(classmod4)
imp <- data.frame(imp$importance)
imp$Descriptor <- rownames(imp)
imp$Class <- "Other"
imp[imp$Descriptor %in% c(
  "GlomTotSA_KW_ratio"
  "ProxTubDiam"), "Class"] <- "Kidney"
imp[imp$Descriptor %in% c(
  "AVERAGE_MASS",
  "LogP_pred",
  "LogVP_pred",
  "LogWS_pred",
  "LogKOA_pred"), "Class"] <- "OPERA"
imp[imp$Descriptor %in% c(
  "TSPC_107.92.6",
  "TSPC_142.62.1",
  "TSPC_111.16.0"), "Class"] <- "Endogenous"
gabstract <- ggplot(imp, aes(x=Descriptor, y=Overall, fill=Class)) +
  geom_bar(stat="identity")+
# Horizontal bar plot
  coord_flip()

print(gabstract)
# In thalf model domain:


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

print("RESULTS 3.2.1")
print("FIRST PARAGRAPH")

print(paste0("we found that the majority (",
  percent(length(unique(subset(DPFASwAD,ClassModDomain==1)$DTXSID)) / length(unique(DPFASwAD$DTXSID))),
  ") of these chemicals fall into the domain of the model."))
  
print(paste0("Across the four species ",
             length(unique(subset(DPFASwAD, ClassModDomain==1)$DTXSID)),
             " PFAS were within the AD (Fig. 3A)."))
             
print(paste0("For humans (over both sexes and dosing methods), ",
             length(unique(subset(DPFASwAD,Species=="Human" & ClassModDomain==1)$DTXSID)),
             " chemicals were estimated to be within AD. Of these, ",
             percent(humanAD.F[["4"]]),
             " were classified in t½ Bin 4, ",
             percent(humanAD.F[["3"]]),
             " were classified in Bin 3, and ",
             percent(humanAD.F[["2"]]),
             " were classified in Bin 2."))

print(paste0("The AM domain further reduces the list to ",
             length(unique(subset(DPFASwAD,AMAD==1)$DTXSID)),
             " of the ",
             length(unique(DPFASwAD$DTXSID)),
             " chemicals."))

print(paste0("For humans, a majority (",
             percent(humanAMAD.F[["4"]]),
             ") of this subset of chemicals were predicted to fall into Bin 4, followed by ",
             percent(humanAMAD.F[["2"]]),
             " in Bin 2 and ",
             percent(humanAMAD.F[["3"]]),
             " in Bin 3."))

print(paste0(length(unique(TCSub$DTXSID)),
             " of the PFAS were in these three classes (Fig. 3B)."))

HAMAD=TCSub[TCSub$Species=="Human" & TCSub$DosingAdj=="Other" & TCSub$Sex=="Female",]
classy.hum <- table(HAMAD$ClassPredFull)/sum(table(HAMAD$ClassPredFull))

print(paste0("For humans, a majority (",
             percent(classy.hum[["4"]]),
             ") of this subset of chemicals were predicted to fall into Bin 4, followed by ",
             percent(classy.hum[["2"]]),
             " in Bin 2 and ",
             percent(classy.hum[["3"]]),
             " in Bin 3."))

print("FIGURE 3 CAPTION")
print(paste0("Distributions of predicted t� for A) ",
             length(unique(subset(DPFASwAD, ClassModDomain==1)$DTXSID)),
             " PFAS within the AD of the model, and B)",
             length(unique(TCSub$DTXSID)),
             " PFAS classified in the same 3 classes as the 11 training set chemicals via ClassyFire."))

             
             
print("RESULTS 3.2.2")
print("SECOND PARAGRAPH")




print("RESULTS 3.2.2")
print("SECOND PARAGRAPH")
print(paste0("For humans (over both sexes and dosing methods), a majority (",
             percent(humanAMAD.F[["4"]]),
             ") of this subset of chemicals were predicted to fall into Bin 4, followed by ",
             percent(humanAMAD.F[["2"]]),
             " in Bin 2 and ",
             percent(humanAMAD.F[["3"]]),
             " in Bin 3."))

print("RESULTS 3.2.3")
print("FIRST PARAGRAPH")

print(paste0("we found that ",
             length(unique(subset(DPFASwAD, ClassModDomain==1)$DTXSID)),
             " PFAS were within the AD. Restricting predictions to only those chemicals whose properties were within the ADs of the OPERA predictors reduced this to ",
             length(unique(subset(DPFASwAD,AMAD==1)$DTXSID)),
             " PFAS. Alternatively, using the ClassyFire chemical structure ontology [96] restricted predictions to ",
             length(unique(TCSub$DTXSID)),
             " PFAS."))

print("SECOND PARAGRAPH")

load(file=paste0("RData/Tox21_AllMods_ADindicate_", readsuff.data,".RData"))
logdindomain <- length(unique(subset(DPFASwAD,ClassModDomain==1)$DTXSID))
load(file=paste0("RData/Tox21_AllMods_ADindicate_", readsuff,".RData"))

print(paste0("increased the number of chemicals for which predictions could be made (from ",
             logdindomain,
             " to ",
             length(unique(subset(DPFASwAD,ClassModDomain==1)$DTXSID)),
             ")."))


print("RESULTS 3.3")
print("SECOND PARAGRAPH")

print(paste0("The majority (",
  percent(humanAD.F[["4"]]),
  ") of PFAS were predicted to be in the longest t� category in humans. "))







