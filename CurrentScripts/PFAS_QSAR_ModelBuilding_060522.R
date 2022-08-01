# Clear the workspace:
rm(list=ls())

# Change to the shared drive:
#setwd("L:/Lab/NCCT_ExpoCast/ExpoCast2022/Dawson_PFAS_HALFLIFE/PFAS_HL_QSAR_2021")

# Specify which RData files we are working with:
readsuff="JFW060522" #Note, this is the suffix for previous work on this page
writesuff="060622" #This is the suffix for ongoing work. 

####Total number predictors prior to pruning:
15:81
85
90
94:110
17
185
186
length(c(15:81, 85, 90,94:110, 1:17,185,186))
####105

# Get rid of logD to increase applicability:
NOLOGD <- FALSE

#Load packages
packages=c("readxl", "openxlsx", "httk", "corrplot", "caret", "randomForest", "OneR", "gdata", "scales")
sapply(packages, require, character.only=TRUE)
##

#Load chemical set assembled from dataset script
load(paste0("PFAS_11Chemicals_QSARdataset_",readsuff,".RData", sep="")) #loads pfasds
dim(pfasds)
length(unique(pfasds$DTXSID)) #11 pfas chemicals
length(which(is.na(pfasds$HLH))) #No missing HLH data



#Substitute opera values from latest Opera run
#3. Merge with Opera predictors run by Wambaugh usig OPERA 2.7
Opera <- read.csv(paste0("PFAS_Catalog/DSSToxPFAS_DTXSID_SMILES_Oct21-1-smi_OPERA2.7Pred.csv"))
which(names(Opera)%in%names(pfasds))
substitutelist=which(names(pfasds)%in%names(Opera))
importlist=which(names(Opera)%in%names(pfasds)[substitutelist])
Operaimport=data.frame("DTXSID"=Opera$MoleculeID, Opera[importlist])
pfasds=pfasds[,-c(substitutelist)]
names(pfasds)[substitutelist]
pfasds=merge(pfasds, Operaimport,by="DTXSID", all.x=TRUE, all.y=FALSE)
pfasdsreorder=pfasds$TrainOrder
save(pfasdsreorder, file=paste0("Training_Data_Reorder_For_Prediction_",writesuff,".RData"))

#Remove missing HLH Data rows
#pfasds=pfasds[-c(which(is.na(pfasds$HLH))),]
#dim(pfasds)

#pfasds1=pfasds
#Data check and get rid of columns with NA's
#NAcheck=apply(pfasds1, 2, FUN=function(x){which(is.na(x))})
#NAlist=as.vector(unlist(lapply(NAcheck, length)))
#pfasds1=pfasds1[,-c(which(NAlist>0))]
#names(pfasds1)

#Variables only
unique(pfasds$PREFERRED_NAME)
#pfasds=pfasds1
#dim(pfasds)

#Screen out opera models in which < 50% are in decided domain for this project (both not in the global domain & with local domain < 0.5) 
##rename henry's law and LogOH model to prevent confusion later
names(pfasds)[which(names(pfasds)=="LogHL_pred")]="LogHenLaw_pred"
names(pfasds)[which(names(pfasds)=="AD_HL")]="AD_HenLaw"
names(pfasds)[which(names(pfasds)=="AD_index_HL")]="AD_index_HenLaw"
names(pfasds)[which(names(pfasds)=="AD_QSAR_HL")]="AD_QSAR_HenLaw"
names(pfasds)[which(names(pfasds)=="LogOH_pred")]="LogAOH_pred"
names(pfasds)[which(names(pfasds)=="LogOH_pred")]="LogAOH_pred"

models=c("LogP", "MP", "BP", "VP", "WS", "HenLaw", "RT", "KOA","pKa", "LogD","AOH", "BCF", "BioDeg", "ReadyBiodeg", "KM", "LogKoc" )
##Pull out single records of each chemical 
uniquechems=pfasds[which(duplicated(pfasds$DTXSID)==FALSE),]
uniquechems 


#Comeup with list of model predictors to remove
removevec=NULL
for (i in models) {
  modname=paste0("AD_QSAR_", i)
  sub=uniquechems[,which(names(uniquechems)==modname)]
  sumsub=sum(sub)
  if (sumsub<=nrow(uniquechems)*0.5) {
    removevec=c(removevec,which(names(pfasds)==modname))
  } else {
    next
  }
}

removenames=names(pfasds)[removevec]
removenames=do.call("rbind",strsplit(removenames, "_"))[,3] #this splits on "_", rbinds them into matrix,and selects out the third column 
removecols=unlist(lapply(removenames, function(x) {grep(x,names(pfasds))} ))
names(pfasds[,removecols])
#Inadvertently grabs Kd_hL_FABP; take this out. Check to make sure that this is correct column 
removecols=removecols[-c(1:3)]
names(pfasds[,removecols]) #Good
pfasds=pfasds[,-c(removecols)] #take out columns for models not meeting criteria. 
names(pfasds)


#Select out predictor columns
idcols=which(names(pfasds)%in%c("DTXSID", "CASRN", "PREFERRED_NAME"))
HLHcol=which(names(pfasds)%in%c("HLH"))
HLHscol=which(names(pfasds)%in%c("HLH_LSsc"))
catcols=which(names(pfasds)%in%c("Species", "Sex", "Type", "DosingAdj"))
#I am disincluding Maximum correlation here because of its difficulty in application to new chemicals 
#Endocols=c(which(names(pfasds)%in%c("MaxEndoPC", "MaxEndoM")), which(grepl('\\TS', names(pfasds))==TRUE) )
Endocols=which(grepl('\\TS', names(pfasds))==TRUE)

KDcols=which(names(pfasds)%in%c("Kd_hL_FABP"))
CMCcols=which(names(pfasds)%in%c("CMC_Pred"))
KAcols=which(names(pfasds)%in%c("Ka_perM_SerAlb_Han"))
Kidneycols=which(names(pfasds)%in%c("BW", "KW", "Neph_Num","KW_BW_ratio","Neph_BW_ratio", "GlomSA",                       
"GlomTotSA", "GlomTotSA_BW_ratio", "GlomTotSA_KW_ratio", "ProxTubLen", "ProxTubDiam", "ProxTubVol",                   
"ProxTubSA" , "ProxTubTotalVol", "ProxTubTotSA", "GlomTotSA_ProxTubTotVol_ratio", "ProxTubTotSA_ProxTotVol_ratio"))
#operacols=c(which(names(pfasds)=="AVERAGE_MASS"), which(grepl('\\AD_' , names(pfasds)) |  grepl('\\_pred$' , names(pfasds)))) #Note, the grepl command here is looking for "_pred", you have to preceed this with // to look for this literal. Also, the $ is an end of line anchor
operacols=c(which(names(pfasds)=="AVERAGE_MASS"), which(grepl('\\_pred$' , names(pfasds)))) #Note, the grepl command here is looking for "_pred", you have to preceed this with // to look for this literal. Also, the $ is an end of line anchor
COCcols=which(names(pfasds)=="COC_aliphatic")
removekidneypreds=which(names(pfasds)%in%c("Log10Neph_Num_pred", "GlomSA_pred", "ProxTubLen_pred", "ProxTubDiam_pred"))
operacols=operacols[-c(which(operacols%in%removekidneypreds))]



#Excluding chain length, since it's not a clear-cut predictor 
#clcols=which(names(pfasds)%in%c("TCChainLength")) #Note; only including total chain length here, after hearing back from Mark Stryner and James McCord


############Endo Variables discretized########
#NOTE: Conversation with Richard Judson suggested that the endogenous similarity variables
#should be discretized based on some high threshold, with the idea that a compound
#has to be very similar to mimic endogenous activity(like, 90% or so). So, similaritities
#are discretized for each PFAS first. Then variables with <5% variance zero variance are removed. Then,
#highly correlated variables are screened out. 

simthreshold=0.9
Endodisc=pfasds[, Endocols]
range(apply(Endodisc, 2, FUN="max"))
Endodisc=ifelse(Endodisc>=simthreshold, 1,0)

#####Pulling out numeric variables
##Removing CMC because relies on dragon predictor
#allvariablemat=data.frame("HLH"=pfasds[,HLHcol],pfasds[,catcols], Endodisc, "Kd_hL_FABP"=pfasds[,KDcols], "CMC_Pred"= pfasds[,CMCcols], "Ka_perM_SerAlb_Han"=pfasds[,KAcols], pfasds[,Kidneycols], pfasds[,operacols],"COC_aliphatic"= pfasds[, COCcols])
allvariablemat=data.frame("HLH"=pfasds[,HLHcol],pfasds[,catcols], Endodisc, "Kd_hL_FABP"=pfasds[,KDcols], "Ka_perM_SerAlb_Han"=pfasds[,KAcols], pfasds[,Kidneycols], pfasds[,operacols],"COC_aliphatic"= pfasds[, COCcols])
numvariablemat=data.frame("HLH"=pfasds[,HLHcol], Endodisc, "Kd_hL_FABP"=pfasds[,KDcols], "Ka_perM_SerAlb_Han" = pfasds[,KAcols], pfasds[,Kidneycols], pfasds[,operacols],"COC_aliphatic"= pfasds[, COCcols])

#Cleaning out low variance variables 
deletevarlist=NULL
varthresh <- 0.05
for(i in 1:length(numvariablemat)){
  col=numvariablemat[,i]
  if(abs(sd(col)/mean(col)) < varthresh | is.na(sd(col)) | sd(col)==0) #Note, used abs here
    deletevarlist=c(deletevarlist, i)
}
write.csv(colnames(numvariablemat)[deletevarlist], 
  paste0("DroppedPredictors_LowVariance_",
  varthresh,"_", writesuff,".csv"))
numvariablemat=numvariablemat[,-deletevarlist] #This reduces the list to 40 numeric variables 


#Remove LogD
if (NOLOGD) {
  names(numvariablemat)
  numvariablemat=numvariablemat[,-c(which(names(numvariablemat)=="LogD55_pred" | 
    names(numvariablemat)=="LogD74_pred")) ]
}

corrmat=cor(numvariablemat, method="spearman")
corrplot(corrmat, method="square")
sort(corrmat[,1])

#Note, most of kidney physiological parameters are correlated with each other, but are not correlated with the other factors
#Similarity values of a number of endogenous chemicals are highly correlated with other 
#Note: corrmat[,which(colnames(corrmat)=="HL-H")] shows that there is no correlation with HL-H that 
#rises above the CorrelationThreshold of 0.9
#Note some Endochemcials have pretty strong negative correlations: 
#Average mass, RT and LogKOA are highly correlated.
plot(log(pfasds$AVERAGE_MASS),log(pfasds$HLH))
plot((pfasds$`TSPC_142-62-1`),log(pfasds$HLH))
plot((pfasds$ProxTubLen),log(pfasds$HLH))

###Reduce total number of predictors by identifying highly(>0.9) correlated predictors
#Note; moving over to the findCorrelation tool, which uses the average overall correlation as the metric instead of the
#the correlation with HLH. This ends up with slightly more predictors 


CorrelationThreshold=0.9
Losers=findCorrelation(corrmat,cutoff = 0.9, exact=TRUE, names=TRUE, verbose=TRUE)
write.csv(Losers, paste0(
  "DroppedPredictors_Correlation_",
  CorrelationThreshold,"_", writesuff,".csv"))

#Now, assemble reduced dataset 
pfasdsred=numvariablemat
pfasdsred=pfasdsred[, -c(which(colnames(pfasdsred)%in%Losers))] #take out losers
colnames(pfasdsred)
pfasdsred=pfasdsred[, -c(which(colnames(pfasdsred)=="HLH"))] #take out HL-H
colnames(pfasdsred)
dim(pfasdsred)
names(pfasdsred) #this reduces it down to a list of 16 variables 

corrmat2=cor(pfasdsred, method="spearman")
check=which(corrmat2>0.9)
corrmat2[check] #Good; all measures are less than 0.9 correlated. 
#Center and Scale based on recorded values of variables
pfasdsredsc=scale(pfasdsred) #scale and center numeric varaibles
pfasdsredsc<-cbind(pfasdsredsc, pfasds[,catcols]) #add categorical variables
pfasdsredsc$Sex<-as.factor(pfasdsredsc$Sex)
pfasdsredsc$Species<-as.factor(pfasdsredsc$Species)
pfasdsredsc$Type<-as.factor(pfasdsredsc$Type)
pfasdsredsc$DosingAdj<-as.factor(pfasdsredsc$DosingAdj)
dim(pfasdsredsc) #Note, this produces a variable matrix with 20 variables 
names(pfasdsredsc)
pfasdsredsc_ds=pfasdsredsc

save(pfasdsred, file=paste("PFAS_Training_Set_Endo_Disc_HLH_unscaled_", writesuff,".RData", sep=""))
load(paste("PFAS_Training_Set_Endo_Disc_HLH_unscaled_", writesuff,".RData", sep=""))
save(pfasdsredsc_ds, file=paste("PFAS_Training_Set_Scaled_Centered_Endo_Disc_HLH_", writesuff,".RData", sep=""))
load(paste0("PFAS_Training_Set_Scaled_Centered_Endo_Disc_HLH_", writesuff,".RData", sep=""))
####For scaling other sets
train_msd=pfasdsred[1:2,]
for(i in 1:dim(train_msd)[2]){
  train_msd[1,i] <- (mean(pfasdsred[,i],na.rm=T))
  train_msd[2,i] <- (sd(pfasdsred[,i],na.rm=T))
}
save(train_msd, file=paste0("Mean_SD_Trainingset_",writesuff,".RData"))
load(paste0("Mean_SD_Trainingset_",writesuff,".RData"))

load(file=paste("PFAS_Training_Set_Scaled_Centered_Endo_Disc_HLH_",writesuff,".RData", sep=""))

#1
############Full MODEL BUILDING#####################
#Needs package "OneR" for bins
##Assign factorial condition to factorial variables,namely gender and species
##Build model using pfasdsred
#First Try Random Forest classification
HLH=pfasds$HLH
save(HLH, file=paste0("HLH_continous_",writesuff,".RData"))
hist(log(HLH))
load(paste0("HLH_continous_",writesuff,".RData"))

#5 Bin determination
HLHBin5=bin(pfasds$HLH, nbins=5,labels=c(1,2,3,4,5), method="content")
HLHBin=cbind(pfasds$HLH, HLHBin5)
HLHBin=HLHBin[order(HLHBin[,1]),] #=
#Shifting it to 0-0.333 day, 0.33-3 days, 3 days-1 month, 1 month-3 months, > 3 months
HLHBin5<-cut(HLH, breaks=c(0,8, 72, 720, 2160, Inf), labels=c(1,2,3,4,5))
table(HLHBin5)/sum(as.numeric(HLHBin5))
hist(as.numeric(HLHBin5))
save(HLHBin5, file=paste0("HLH_Bins_5_", writesuff,".RData"))

#4 Bin determination 
HLHBin4=bin(pfasds$HLH, nbins=4,labels=c(1,2,3,4), method="content")
HLHBin4=cbind(pfasds$HLH, HLHBin4)
HLHBin4=HLHBin4[order(HLHBin4[,1]),] #=


#Going to set bins at <12 hours, 12 hours-1 week, 1 week-2 months, and >2 months
HLHBin4<-cut(HLH, breaks=c(0,12, 168, 1440,Inf), labels=c(1,2,3,4)) #this breaks down to half day, half-day to 1 week, 1 week and 2 months, and >2 months 
hist(as.numeric(HLHBin4))
table(HLHBin4)/length(HLH)
save(HLHBin4, file=paste0("HLH_Bins_4_", writesuff,".RData"))
load(paste0("HLH_Bins_4_", writesuff,".RData"))

#subset(train, HLHBin4==4)

#3 Bin determination 
HLHBin3=bin(pfasds$HLH, nbins=3,labels=c(1,2,3), method="content")
HLHBin3=cbind(pfasds$HLH, HLHBin3)
HLHBin3=HLHBin3[order(HLHBin3[,1]),] #=
#less than 1.5 days, 1.5days-36 days, >36 days)
HLHBin3<-cut(HLH, breaks=c(0,36,864,Inf), labels=c(1,2,3))  #changed this less than 1.5 days, between 1.5 days and 36 days, and > 36 days 
table(HLHBin3)/sum(as.numeric(HLHBin3))
save(HLHBin3, file=paste0("HLH_Bins_3_", writesuff,".RData"))

train=data.frame(pfasdsredsc_ds)
train=train[,-c(which(names(train)=="Species" | names(train)=="Type"))]
save(train, file=paste0("FullTrainingSet_",writesuff,".RData"))
load(paste0("FullTrainingSet_",writesuff,".RData"))
######Conduct 10 fold cross validation on overall dataset using caret
control=trainControl(method="repeatedcv", number = 10, repeats = 10, savePredictions = "all")
seed = 12345
metric <- "Accuracy"

#3Bin
set.seed(seed)
out <- tuneRF(x=train, 
              y=HLHBin3,
              mtryStart=floor(length(train[1,])^(1/2)), #j is sizes 
              ntreeTry=500,
              stepFactor=1.25,
              improve=0.005,
              trace=FALSE, 
              plot=FALSE, 
              dBest=FALSE)
bestmtry <- out[which(out[,2]==min(out[,2])),1]
if(length(bestmtry)>1){bestmtry= max(bestmtry)}
tunegrid <- expand.grid(.mtry=bestmtry)
set.seed(seed)
classmod3=train(x=train,y=HLHBin3, method="rf", metric=metric, trControl=control, tuneGrid = tunegrid)
classmod3$bestTune
classmod3$finalModel
classmod3$trainingData
classmod3$results # the average across all cross-validations and reps
hist(classmod3$resample$Accuracy)
save(classmod3, file=paste("ClassificationModel_3Bin_EndoSimDisc_HLH_", writesuff, ".RData",sep=""))
load(paste("ClassificationModel_3Bin_EndoSimDisc_HLH_", writesuff, ".RData",sep=""))
#84.4 accuracy

#4Bin
set.seed(seed)
out <- tuneRF(x=train, 
              y=HLHBin4,
              mtryStart=floor(length(train[1,])^(1/2)), #j is sizes 
              ntreeTry=500,
              stepFactor=1.25,
              improve=0.005,
              trace=FALSE, 
              plot=FALSE, 
              dBest=FALSE)
bestmtry <- out[which(out[,2]==min(out[,2])),1]
if(length(bestmtry)>1){bestmtry= max(bestmtry)}
tunegrid <- expand.grid(.mtry=bestmtry)
set.seed(seed)
classmod4=train(x=train,y=HLHBin4, method="rf", metric=metric, trControl=control, tuneGrid = tunegrid)
classmod4$bestTune
classmod4$finalModel
classmod4$trainingData
classmod4$results # the average across all cross-validations and reps
#86.6% accuracy  
hist(classmod4$resample$Accuracy)
table(HLHBin4)/sum(table(HLHBin4))
varImp(classmod4)

save(classmod4, file=paste("ClassificationModel_4Bin_EndoSimDisc_HLH_", writesuff, ".RData",sep=""))
load(paste("ClassificationModel_4Bin_EndoSimDisc_HLH_", writesuff, ".RData",sep=""))



#5Bin 
set.seed(seed)
out <- tuneRF(x=train, 
              y=HLHBin5,
              mtryStart=floor(length(train[1,])^(1/2)), #j is sizes 
              ntreeTry=500,
              stepFactor=1.25,
              improve=0.005,
              trace=FALSE, 
              plot=FALSE, 
              dBest=FALSE)
bestmtry <- out[which(out[,2]==min(out[,2])),1]
if(length(bestmtry)>1){bestmtry= max(bestmtry)}
tunegrid <- expand.grid(.mtry=bestmtry)
set.seed(seed)
classmod5=train(x=train,y=HLHBin5, method="rf", metric=metric, trControl=control, tuneGrid = tunegrid)
classmod5$bestTune
classmod5$finalModel
classmod5$trainingData
classmod5$results # the average across all cross-validations and reps
hist(classmod5$resample$Accuracy)
save(classmod5, file=paste("ClassificationModel_5Bin_EndoSimDisc_HLH_", writesuff, ".RData",sep=""))
#76.2% accuracy
#load(paste("ClassificationModel_5Bin_EndoSimDisc_HLH_", suff, ".RData",sep=""))
varImp(classmod5)

# ########RF Regression 
# metric <- "RMSE"
# set.seed(seed)
# out <- tuneRF(x=train, 
#               y=HLH,
#               mtryStart=floor(length(train[1,])^(1/2)), #j is sizes 
#               ntreeTry=500,
#               stepFactor=1.25,
#               improve=0.005,
#               trace=TRUE, 
#               plot=TRUE, 
#               dBest=FALSE)
# bestmtry <- out[which(out[,2]==min(out[,2])),1]
# set.seed(seed)
# regmod=train(x=train,y=HLH, method="rf", importance=TRUE, metric=metric, trControl=control, tuneGrid = tunegrid)
# regmod$results
# hist(regmod$resample$Rsquared)
# regmod$finalModel #R2=0.68+- 0.32
#                   #78.89% of residuals
#                   #12333.17 +- 8898.399
# varImp(regmod)
# save(regmod, file=paste("RegressionModel_EndoSimDisc_HLH_", writesuff, ".RData", sep=""))
# load(paste("RegressionModel_EndoSimDisc_HLH_", writesuff, ".RData", sep=""))

####
load(paste("ClassificationModel_3Bin_EndoSimDisc_HLH_", writesuff, ".RData",sep=""))
load(paste("ClassificationModel_4Bin_EndoSimDisc_HLH_", writesuff, ".RData",sep=""))
load(paste("ClassificationModel_5Bin_EndoSimDisc_HLH_", writesuff, ".RData",sep=""))
#load(paste("RegressionModel_EndoSimDisc_HLH_", writesuff, ".RData", sep=""))

#Compare to model with and out Chain Length and a parameter
##2:Recursive feature elimination 
#################

#Note: rfe kind of works in a hierchical way. The main function rfe, require specifications for things in the rfeControl function. And, the rfeControl
#function requires specification of a set of functions. There are convenience functions for RF already set up with rfFuncs below:
#rfFuncs$summary
#rfFuncs$fit
#rfFuncs$pred
#rfFuncs$rank
#rfFuncs$selectSize
#rfFuncs$selectVar

#However, you could specify your own helper function 
####

#You need the "ctrl" object to pass into rfe, but the functions argument of that function could reference the list of other functions, like rfRFE above. 
#So, if you wanted to use the pickSizeTolerance instead of pickSizeBest, you could simply specify it as the argument for SelectSize above. 

#These are kidney descriptors retaied in the training set 
kc=c("BW", "KW", "Neph_Num","KW_BW_ratio","Neph_BW_ratio", "GlomSA",                       
   "GlomTotSA", "GlomTotSA_BW_ratio", "GlomTotSA_KW_ratio", "ProxTubLen", "ProxTubDiam", "ProxTubVol",                   
   "ProxTubSA" , "ProxTubTotalVol", "ProxTubTotSA", "GlomTotSA_ProxTubTotVol_ratio", "ProxTubTotSA_ProxTotVol_ratio")
kd=names(train)[which(names(train)%in%kc)]

save(kd, file=paste0("Kidney_Descriptors_in_Training_Set_", writesuff,".RData"))

#Classification 
ctrl <- rfeControl(functions = rfFuncs,
                   method = "repeatedcv",
                   repeats = 10,
                   verbose = TRUE,
                   rerank=TRUE,
                   number=10)

# #RFE process for regression model 
# subsets=c(1:16)
# set.seed(12345)
# rfeReg <- rfe(train, HLH,
#                 rfeControl = ctrl,
#                 sizes=subsets,
#                 metric="RMSE",
#                 maximize=FALSE)
# save(rfeReg, file=paste0("RFE_Reg_", writesuff, ".RData"))
# rfeReg
# rfeReg$bestSubset
# rfeReg$optVariables
# RMresults=rfeReg$results
# regvars=rfeReg$optVariables



#Ok, so based on 
# #Now, fit regression model with regvars
# metric <- "RMSE"
# set.seed(seed)
# out <- tuneRF(x=train[,which(names(train)%in%regvars)], 
#               y=HLH,
#               mtryStart=floor(length(train[1,which(names(train)%in%regvars)])^(1/2)), #j is sizes 
#               ntreeTry=500,
#               stepFactor=1.25,
#               improve=0.005,
#               trace=TRUE, 
#               plot=TRUE, 
#               dBest=FALSE)
# bestmtry <- out[which(out[,2]==min(out[,2])),1]
# if(length(bestmtry)>1){bestmtry= max(bestmtry)}
# tunegrid <- expand.grid(.mtry=bestmtry)
# 
# set.seed(seed)
# regmodred=train(x=train[,which(names(train)%in%regvars)],y=HLH, method="rf", importance=TRUE, metric=metric, trControl=control, tuneGrid = tunegrid)
# regmodred$results #R2 of 0.6850606 for cross-validated model using random seed 12345 
# 
# hist(regmodred$resample$Rsquared)
# regmodred$finalModel #Final model %Var explained 84.02
# #save(regmodred, file=paste("RegressionModel_PostRFE_HLH_", writesuff, ".RData", sep=""))
# load(file=paste("RegressionModel_PostRFE_HLH_", writesuff, ".RData", sep=""))

##Classification model RFE process
subsets=c(2:16)
set.seed(12345)
rfeClass <- rfe(train, HLHBin4,
                rfeControl = ctrl,
                sizes=subsets,
                metric="Accuracy",
                maximize=TRUE)

save(rfeClass, file=paste0("ClassMod_RFE_FullSet_", writesuff,".RData"))
load(paste0("ClassMod_RFE_FullSet_", writesuff,".RData"))

rfeClass$optVariables
rfeClass$results
classvars=rfeClass$optVariables

#Now, fit classification model with classvars
metric <- "Accuracy"
set.seed(seed)
out <- tuneRF(x=train[,which(names(train)%in%classvars)], 
              y=HLHBin4,
              mtryStart=floor(length(train[1,which(names(train)%in%classvars)])^(1/2)), #j is sizes 
              ntreeTry=1000,
              stepFactor=1.25,
              improve=0.005,
              trace=TRUE, 
              plot=TRUE, 
              dBest=FALSE)
bestmtry <- out[which(out[,2]==min(out[,2])),1]
if(length(bestmtry)>1){bestmtry= max(bestmtry)}
tunegrid <- expand.grid(.mtry=bestmtry)
set.seed(seed)
classmodred=train(x=train[,which(names(train)%in%classvars)],y=HLHBin4, method="rf",  metric=metric, trControl=control, tuneGrid = tunegrid)
classmodred$results #R2 of 0.866 for cross-validated model 
classmodred
hist(classmodred$resample$Accuracy)
classmodred$finalModel #OOB 12.09%
save(classmodred, file=paste("ClassificationModel_PostRFE_HLHBin4_", writesuff, ".RData", sep=""))
load(file=paste("ClassificationModel_PostRFE_HLHBin4_", writesuff, ".RData", sep=""))


#3. Variable Importance and Histograms of average Rsquared and classication error rate 


#full models and reduced models
##Full models
load(paste("ClassificationModel_4Bin_EndoSimDisc_HLH_", writesuff, ".RData",sep=""))
# load(paste("RegressionModel_EndoSimDisc_HLH_", writesuff, ".RData", sep=""))
##Reduced Models# Note: the reduced classification model is the same as the full model 
# load(file=paste("RegressionModel_PostRFE_HLH_", writesuff, ".RData", sep=""))


par(mfrow=c(1,1),
    mar=c(5,4.0,4.1,1))


##Variable importance 
varimportlist=list() #data list

#Classification model 
#save data
ClassFullVarImp=varImp(classmod4, estimate="Accuracy", scale=TRUE)
classmodimp=as.data.frame(varImp(classmod4, estimate="Accuracy", scale=FALSE)$importance)
names(classmodimp)[1]="Change_RawAccuracy"
classmodimpsc=as.data.frame(varImp(classmod4,estimate="Accuracy", scale=TRUE)$importance)
names(classmodimpsc)[1]="Change_ScaledAccuracy"
classmodimp=cbind(classmodimp, classmodimpsc)
classmodimp=classmodimp[order(classmodimp$Change_RawAccuracy, decreasing=TRUE),]
varimportlist=c(varimportlist, "ClassMod_Full&Red"=list( classmodimp ))

#Save figure
jpeg(filename=paste0("Figures/ClassModFull_and_Reduced_Importance", writesuff,".png", sep=""), width=600, height=600)
plot(ClassFullVarImp, main="Classification Model:Full and Reduced")
dev.off()


#Regression Model 

# #save data
# regmodfullimp=as.data.frame(varImp(regmod, estimate="RMSE", scale=FALSE)$importance)
# names(regmodfullimp)[1]="Change_RawRSME"
# regmodfullimpsc=as.data.frame(varImp(regmod,estimate="RSME", scale=TRUE)$importance)
# names(regmodfullimpsc)[1]="Change_ScaledRSME"
# regmodredimp=as.data.frame(varImp(regmodred, estimate="RMSE", scale=FALSE)$importance)
# names(regmodredimp)[1]="Change_RawRSME"
# regmodredimpsc=as.data.frame(varImp(regmodred,estimate="RSME", scale=TRUE)$importance)
# names(regmodredimpsc)[1]="Change_ScaledRSME"
# 
# regmodfullimp=cbind(regmodfullimp, regmodfullimpsc)
# regmodfullimp=regmodfullimp[order(regmodfullimp$Change_RawRSME, decreasing=TRUE),]
# varimportlist=c(varimportlist, "RegMod_Full"=list( regmodfullimp ))
# 
# regmodredimp=cbind(regmodredimp, regmodredimpsc)
# regmodredimp=regmodredimp[order(regmodredimp$Change_RawRSME, decreasing=TRUE),]
# varimportlist=c(varimportlist, "RegMod_Red"=list( regmodredimp ))
# 
# save(varimportlist, file=paste0("VariableImportance_Class_Reg_Mods_", writesuff, ".RData"))


#Save figure
jpeg(filename=paste0("Figures/ClassModFull_and_Reduced_Importance", writesuff,".png", sep=""), width=600, height=600)
plot(ClassFullVarImp, main="Classification Model:Full and Reduced")
dev.off()

# jpeg(filename=paste0("Figures/RegModFull_Importance", writesuff,".png", sep=""), width=600, height=600)
# plot(varImp(regmod,estimate="RSME", scale=TRUE), main="Regression Model:Full")
# dev.off()
# 
# jpeg(filename=paste0("Figures/RegMod_Reduced_Importance", writesuff,".png", sep=""), width=600, height=600)
# plot(varImp(regmodred,estimate="RSME", scale=TRUE), main="Regression Model: Reduced")
# dev.off()
# 
# plot(RegFullVarImp)
# plot(RegRedVarImp)


#Figures of the Rsquares and Accuracies of the models
# jpeg(filename=paste0("Figures/RegModFull_R2Hist_HLH", writesuff,".png", sep=""), width=600, height=600)
# hist(regmod$resample$Rsquared, main="RF-Full Regression of HL: R2 distribution", xlab="R2", cex.main=0.9)
# dev.off()
# 
# jpeg(paste("Figures/RegModRed_R2Hist_HLH", writesuff,".jpeg", sep=""),width=600, height=600)
# hist(regmodred$resample$Rsquared, main="RF-Reduced Regression of HL: R2 distribution", xlab="R2", cex.main=0.9)
# dev.off()


jpeg(paste("Figures/ClassModFull_R2Hist_HLH", writesuff,".jpeg", sep=""),width=600, height=600)
hist(classmod4$resample$Accuracy, main="RF Classification of HL: Avg OOB Error", xlab="% OOB Error", cex.main=0.9)
dev.off()


#3. Predicted vs Observed Half-life 

#####
set.seed(seed)
classpredfull=predict(classmod4)
# set.seed(seed)
# regpredfull=predict(regmod)
# set.seed(seed)
# regpredred=predict(regmodred)

#Save data and prediction set
pfasdsredsc_ds$ClassPredFull=classpredfull
# pfasdsredsc_ds$RegPredFull=regpredfull
# pfasdsredsc_ds$RegPredRed=regpredred
pfasdsredsc_ds$HLH=pfasds[,HLHcol]
pfasdsredsc_ds=data.frame(pfasds[,idcols], pfasdsredsc_ds)
save(pfasdsredsc_ds, file=paste0("PFAS_QSAR_Predictors_Predictions_", writesuff,".RData"))
#This loads pfasdsredsc_ds with HLH
load(file=paste0("PFAS_QSAR_Predictors_Predictions_", writesuff,".RData"))

#Continue on with analysis of fit
pfasdsfitted=pfasds
# pfasdsfitted$RegPredFull=regpredfull
# pfasdsfitted$RegPredRed=regpredred
pfasdsfitted$ClassPredFull=classpredfull

# ###Now, scale observed and predicted by lifespan
# pfasdsfitted$HLH_PerLF=pfasdsfitted$HLH/pfasdsfitted$HrLifeSpan
# pfasdsfitted$RegPredFull_PerLF=pfasdsfitted$RegPredFull/pfasdsfitted$HrLifeSpan
# pfasdsfitted$RegPredRed_PerLF=pfasdsfitted$RegPredRed/pfasdsfitted$HrLifeSpan

##Calcuate means for each bin for empirical data
HLHBin4<-cut(HLH, breaks=c(0,12, 168, 1440,Inf), labels=c(1,2,3,4)) #this breaks down to half day, half-day to 1 week, 1 week and 2 months, and >2 months 
agtabHLH=aggregate(HLH~HLHBin4,data=pfasdsfitted, FUN="median")
pfasdsfitted$HLHBin4=HLHBin4
pfasdsfitted$Bin4ClassFullMean=ifelse(pfasdsfitted$ClassPredFull==1, agtabHLH[1,2], ifelse(pfasdsfitted$ClassPredFull==2, agtabHLH[2,2],ifelse(pfasdsfitted$ClassPredFull==3, agtabHLH[3,2],agtabHLH[4,2])))
pfasdsfitted$Bin4ClassFullMean_perLS=pfasdsfitted$Bin4ClassFullMean/pfasdsfitted$HrLifeSpan
save(agtabHLH, file=paste0("Median_HLH_per_Bin_4_", writesuff,".RData"))


#Filling in abbreviations
scientific_10 <- function(x) {                                  
  out <- gsub("1e", "10^", scientific_format()(x))              
  out <- gsub("//+","",out)                                     
  out <- gsub("10//^01","10",out)                               
  out <- parse(text=gsub("10//^00","1",out))                    
}  


pfasdsfitted$Sex=as.factor(pfasdsfitted$Sex)
pfasdsfitted$Species=as.factor(pfasdsfitted$Species)
pfasdsfitted$Type=as.factor(pfasdsfitted$Type)
pfasdsfitted$DosingAdj=as.factor(pfasdsfitted$DosingAdj)

# #Regression: Obs vs predicted
# #Full
# ObsPredPlotRFRegFull_Chem_Species <- ggplot(pfasdsfitted, aes(y=HLH, x=RegPredFull)) +
#   geom_point(size=7,aes(shape=PREFERRED_NAME,color=Species)) +
#   scale_x_log10(label=scientific_10)+ scale_y_log10(label=scientific_10) +
#   scale_shape_manual(values=c(1,2,3,4,5,6,7,8,9,10,11))+
#   xlab("RF Regression Predicted Serum Half-Life (hr)") +
#   ylab(expression(paste(italic("In Vivo"), " Serum Half-Life (hr)",sep=""))) +
#   ggtitle("11 PFAS Chemicals in 4 Species", subtitle = "RF Regression Model: Full")+
#   theme(axis.text=element_text(size=15), 
#         plot.title=element_text(size=20),
#         plot.subtitle = element_text(size=17),
#         axis.title=element_text(size=17),
#         legend.title=element_text(size=20),
#         legend.text=element_text(size=15)) +
#   labs(shape="Chemical")
# 
# ObsPredPlotRFRegFull_Chem_Species
# png(paste0("Figures/ObsvPred_RFRegFull_Pred_by_Species_Chem",writesuff,".png",sep=""), width=1200, height=744)
# ObsPredPlotRFRegFull_Chem_Species
# dev.off()
# 
# ###
# #Reduced
# ObsPredPlotRFRegRed_Chem_Species <- ggplot(pfasdsfitted, aes(y=HLH, x=RegPredRed)) +
#   geom_point(size=7,aes(shape=PREFERRED_NAME,color=Species)) +
#   scale_x_log10(label=scientific_10)+ scale_y_log10(label=scientific_10) +
#   scale_shape_manual(values=c(1,2,3,4,5,6,7,8,9,10,11))+
#   xlab("RF Regression Predicted Serum Half-Life (hr)") +
#   ylab(expression(paste(italic("In Vivo"), " Serum Half-Life (hr)",sep=""))) +
#   ggtitle("11 PFAS Chemicals in 4 Species",subtitle = "RF Regression Model: Reduced")+
#   theme(axis.text=element_text(size=15), 
#         plot.title=element_text(size=20),
#         plot.subtitle = element_text(size=17),
#         axis.title=element_text(size=17),
#         legend.title=element_text(size=20),
#         legend.text=element_text(size=15)) +
#   labs(shape="Chemical")
# 
# ObsPredPlotRFRegRed_Chem_Species
# png(paste("Figures/ObsvPred_RFRegRed_Pred_by_Species_Chem",writesuff,".png",sep=""), width=1200, height=744)
# ObsPredPlotRFRegRed_Chem_Species
# dev.off()


#Classification Model:Obs vs predict using mean values for classes
#Full
jitter <- position_jitter(width = 0.3, height = 0)
ObsPredPlotRFClassFull_Chem_Species <- ggplot(pfasdsfitted, aes(y=HLH, x=Bin4ClassFullMean)) +
  geom_point(size=7,position = jitter, aes(shape=PREFERRED_NAME,color=Species)) +
  scale_x_log10(label=scientific_10)+ scale_y_log10(label=scientific_10) +
  scale_shape_manual(values=c(1,2,3,4,5,6,7,8,9,10,11))+
  xlab("Classification Model Bin Means of Predicted Serum Half-Life (hr)") +
  ylab(expression(paste(italic("In vivo"), " Serum Half-Life (hr)",sep=""))) +
  ggtitle("11 PFAS Chemicals in 4 Species",subtitle = "Random Forest Classification Model")+
  
  theme(axis.text=element_text(size=15), 
        plot.title=element_text(size=20),
        plot.subtitle = element_text(size=17),
        axis.title=element_text(size=17),
        legend.title=element_text(size=20),
        legend.text=element_text(size=15)) +
  labs(shape="Chemical")


obsplot=ObsPredPlotRFClassFull_Chem_Species +
  geom_hline(yintercept=24,linetype="dashed", color = "red") +
  geom_hline(yintercept=168, linetype="dashed", color = "blue") + 
  geom_hline(yintercept=1440,linetype="dashed", color = "black") 

png(paste("Figures/ObsvPred_RFClassPredFull_andRed_by_Species_Chem_",writesuff,".png",sep=""), width=1200, height=744)
obsplot
dev.off()




#######Y-randomization#############
load(file=paste0("PFAS_QSAR_Predictors_Predictions_", writesuff,".RData"))
pfaseval=pfasdsredsc_ds #with HLH
pfaseval$TrainOrder=seq(1,length(pfaseval[,1])) #This is so we can re-order the y-randomization lines by species and chem so everything lines up except for the y-values
species=unique(pfaseval$Species)
sex=unique(pfaseval$Sex)
dosing=unique(pfaseval$DosingAdj)
chems=unique(pfaseval$DTXSID)
set.seed(seed)

##Y-rand approaches
###Overall Y-rand
YROverall=sample(1:91,91,replace=FALSE)
HLHYRandOverall=pfaseval[YROverall, "HLH"]
HLHYROverall_Bin<-cut(HLHYRandOverall, breaks=c(0,12, 168, 1440,Inf), labels=c(1,2,3,4)) #this breaks down to half day, half-day to 1 week, 1 week and 2 months, and >2 months 


###Y-chemicals within species
HLHYRspecies=NULL
for (i in species){
  sub=pfaseval[pfaseval$Species==i,]
  sub1=sample(1:dim(sub)[1],dim(sub)[1],replace=FALSE)
  sub$HLHYRspecies=sub[sub1, "HLH"]
  HLHYRspecies=rbind(HLHYRspecies, sub)
}
HLHYRspecies=HLHYRspecies[order(HLHYRspecies$TrainOrder),]
HLHYRspecies_Bin<-cut(HLHYRspecies$HLHYRspecies, breaks=c(0,12, 168, 1440,Inf), labels=c(1,2,3,4)) #this breaks down to half day, half-day to 1 week, 1 week and 2 months, and >2 months 


###Y-rand within chemicals
HLHYRchem=NULL
for (i in chems){
  sub=pfaseval[pfaseval$DTXSID==i,]
  sub1=sample(1:dim(sub)[1],dim(sub)[1],replace=FALSE)
  sub$HLHYRchem=sub[sub1, "HLH"]
  HLHYRchem=rbind(HLHYRchem, sub)
}
HLHYRchem=HLHYRchem[order(HLHYRchem$TrainOrder),]
HLHYRchem_Bin<-cut(HLHYRchem$HLHYRchem, breaks=c(0,12, 168, 1440,Inf), labels=c(1,2,3,4)) #this breaks down to half day, half-day to 1 week, 1 week and 2 months, and >2 months 

#load training set
load(paste0("FullTrainingSet_",writesuff,".RData"))
names(train)

#Classification Model 
#Overall y-randomization
set.seed(seed)
out <- tuneRF(x=train, 
              y=HLHYROverall_Bin,
              mtryStart=floor(length(train[1,])^(1/2)), #j is sizes 
              ntreeTry=500,
              stepFactor=1.25,
              improve=0.005,
              trace=FALSE, 
              plot=FALSE, 
              dBest=FALSE)
bestmtry <- out[which(out[,2]==min(out[,2])),1]
if(length(bestmtry)>1){bestmtry= max(bestmtry)}
tunegrid <- expand.grid(.mtry=bestmtry)
set.seed(seed)
classmod4YROverall=train(x=train,y=HLHYROverall_Bin, method="rf", metric=metric, trControl=control, tuneGrid = tunegrid)
classmod4YROverall$bestTune
classmod4YROverall$finalModel
classmod4YROverall$trainingData
classmod4YROverall$results # the average across all cross-validations and reps
#16.0+-11.3% accuracy  
hist(classmod4YROverall$resample$Accuracy)
mean(classmod4YROverall$resample$Accuracy)
sd(classmod4YROverall$resample$Accuracy)
table(HLHYROverall_Bin)/sum(table(HLHYROverall_Bin))
varImp(classmod4YROverall)

#Y-randomization across chemicals within species 
set.seed(seed)
out <- tuneRF(x=train, 
              y=HLHYRspecies_Bin,
              mtryStart=floor(length(train[1,])^(1/2)), #j is sizes 
              ntreeTry=500,
              stepFactor=1.25,
              improve=0.005,
              trace=FALSE, 
              plot=FALSE, 
              dBest=FALSE)
bestmtry <- out[which(out[,2]==min(out[,2])),1]
if(length(bestmtry)>1){bestmtry= max(bestmtry)}
tunegrid <- expand.grid(.mtry=bestmtry)
set.seed(seed)
classmod4YRspecies=train(x=train,y=HLHYRspecies_Bin, method="rf", metric=metric, trControl=control, tuneGrid = tunegrid)
classmod4YRspecies$bestTune
classmod4YRspecies$finalModel
classmod4YRspecies$trainingData
classmod4YRspecies$results 
#41.5+-14.46% accuracy  
hist(classmod4YRspecies$resample$Accuracy)
mean(classmod4YRspecies$resample$Accuracy)
sd(classmod4YRspecies$resample$Accuracy)
table(HLHYRspecies_Bin)/sum(table(HLHYRspecies_Bin))
varImp(classmod4YRspecies)


#Y-randomization across species within chemicals
set.seed(seed)
out <- tuneRF(x=train, 
              y=HLHYRchem_Bin,
              mtryStart=floor(length(train[1,])^(1/2)), #j is sizes 
              ntreeTry=500,
              stepFactor=1.25,
              improve=0.005,
              trace=FALSE, 
              plot=FALSE, 
              dBest=FALSE)
bestmtry <- out[which(out[,2]==min(out[,2])),1]
if(length(bestmtry)>1){bestmtry= max(bestmtry)}
tunegrid <- expand.grid(.mtry=bestmtry)
set.seed(seed)
classmod4YRchem=train(x=train,y=HLHYRchem_Bin, method="rf", metric=metric, trControl=control, tuneGrid = tunegrid)
classmod4YRchem$bestTune
classmod4YRchem$finalModel
classmod4YRchem$trainingData
classmod4YRchem$results 
#52.9+-14.9% accuracy  
hist(classmod4YRchem$resample$Accuracy)
mean(classmod4YRchem$resample$Accuracy)
sd(classmod4YRchem$resample$Accuracy)
table(HLHYRchem_Bin)/sum(table(HLHYRchem_Bin))
varImp(classmod4YRchem)

#Make y-randomization list for export
YRData=data.frame(pfaseval,HLHYRandOverall, HLHYROverall_Bin, "HLHYRspecies"=HLHYRspecies$HLHYRspecies, HLHYRspecies_Bin,"HLHYRchem"=HLHYRchem$HLHYRchem, HLHYRchem_Bin  )
write.csv(YRData, file="All_Training_Data_ " )
YRresults=list(YRData,classmod4YROverall,classmod4YRspecies,classmod4YRchem)
names(YRresults)=c("YR_Data","YR_Overall", "YR_by_Species", "YR_by_Chem")
save(YRresults, file=paste0("YR_Analysis_",writesuff,".RData"))
save(classmod4YROverall,classmod4YRspecies,classmod4YRchem,
  file=paste("ClassificationModel_YRand_HLHBin4_", writesuff, ".RData", sep=""))


####Make final export set for training chemicals

CTD=data.frame(YRData,"Bin4ClassMean"=pfasdsfitted$Bin4ClassFullMean)
write.csv(CTD, file=paste0("Complete_Training_Data_", writesuff,".csv"))

##########Calculate concordance and RSME###########
#Classification Model
#Final models used to predict the training data
confusionMatrix(pfasdsfitted$HLHBin4, pfasdsfitted$ClassPredFull) #0.978
#accuracy=0.978

#Regression model
# mod1=lm(pfasdsfitted$HLH~pfasdsfitted$RegPredFull) 
# plot(pfasdsfitted$HLH, pfasdsfitted$RegPredFull)
# summary(mod1) #R2=0.976
# mod2=lm(pfasdsfitted$HLH~pfasdsfitted$RegPredRed) 
# summary(mod2)#R2=0.967
# 
# RMSE(pfasdsfitted$HLH, pfasdsfitted$RegPredFull) #6975.011
# RMSE(pfasdsfitted$HLH, pfasdsfitted$RegPredRed) #7629.578
# 
# pfasdsfitted$RegFullRedDiff=pfasdsfitted$RegPredFull-pfasdsfitted$RegPredRed
# hist(pfasdsfitted$RegFullRedDiff)
# RMSE(pfasdsfitted$RegPredFull, pfasdsfitted$RegPredRed)#1293.801
# MAE(pfasdsfitted$RegPredFull, pfasdsfitted$RegPredRed)#575.96
# #hist(log10(abs(pfasdsfitted$RegFullRedDiff)))
# median(pfasdsfitted$RegFullRedDiff)-33.23
# #RMSE=1356 ; this might be a good choice, since although most differences are less than 1000 hours, they range up to 
# #MAE=819






###Training Set Summary
#chemical coverage
suffbuilding="DED060121"
load(paste0("PredictorSet_pre-interpolation_", suffbuilding))
df5$count=1
df5$Ka_perM_SerAlb_Han=ifelse(is.na(df5$Ka_perM_SerAlb_Han),NA, 1) 
aggregate(df5$count~df5$Ka_perM_SerAlb_Han*df5$DTXSID*df5$Sex, FUN='sum')
df5$Kd_hL_FABP=ifelse(is.na(df5$Kd_hL_FABP),NA, 1) 
aggregate(df5$count~df5$Kd_hL_FABP*df5$DTXSID*df5$Sex, FUN='sum')

load(paste("PFAS_Training_Set_Endo_Disc_HLH_unscaled_", writesuff,".RData", sep=""))
tdd=t(data.frame(apply(pfasdsred, 2, "summary")))
saveRDS(tdd,paste0("Training_Data_Description_",writesuff,".rds"))

load(paste0("PFAS_11Chemicals_QSARdataset_",readsuff,".RData", sep="")) #loads pfasds
a=aggregate(ProxTubDiam~Species, data=pfasds, FUN="median")
b=aggregate(ProxTubLen~Species, data=pfasds, FUN="median")
c=aggregate(KW_BW_ratio~Species, data=pfasds, FUN="median")
d=aggregate(GlomTotSA_KW_ratio ~Species, data=pfasds, FUN="median")
data.frame(a,b[,2], c[,2], d[,2])
e=Reduce(function(x, y) merge(x, y, all=TRUE, by="Species"), list(a, b,c,d))



#Look at VP and LogD7.4
###LogD74_pred
ggplot(pfasds, aes(x=LogP_pred, y=log(HLH+1))) +
  geom_point(size=3, shape=16)+
  xlab(expression(LogD[7.4])) +
  ylab(expression(paste("Log ",italic("t")[1/2]))) +
  ggtitle(expression(paste("Log ",italic("t")[1/2], " as a function of ", LogP)))

###LogVP_pred
ggplot(pfasds, aes(x=LogVP_pred, y=log(HLH+1))) +
  geom_point(size=3, shape=16)+
  xlab("Log Vapor Pressure (mmHg)") +
  ylab(expression(paste("Log ",italic("t")[1/2]))) +
ggtitle(expression(paste("Log ",italic("t")[1/2], " as a function of Log Vapor Pressure")))

#Look at correlations with top 5 predictors 
d=pfasdsredsc_ds
plot(pfasds)
plot(pfasds$AVERAGE_MASS, log(pfasds$HLH))
cor.test(pfasds$AVERAGE_MASS, log(pfasds$HLH), method="spearman") #0.71
plot(pfasds$ProxTubLen,log(pfasds$HLH))
cor.test(pfasds$ProxTubLen, log(pfasds$HLH), method="spearman") #0.303
plot((pfasds$ProxTubSA),log(pfasds$HLH))
cor.test(pfasds$ProxTubSA, log(pfasds$HLH), method="spearman") #0.45
plot((pfasds$LogVP_pred),log(pfasds$HLH))
cor.test((pfasds$LogVP_pred),log(pfasds$HLH), method="spearman") #-40
#plot(d$Kd_hL_FABP, log(pfasds$HLH))
#cor.test(d$Kd_hL_FABP, log(pfasds$HLH), method="spearman")#-0.18
plot(d$LogP_pred, log(pfasds$HLH))
cor.test(d$LogP_pred, log(pfasds$HLH), method="spearman")#-0.18
plot(d$COC_aliphatic, log(pfasds$HLH))
cor.test(d$COC_aliphatic, log(pfasds$HLH), method="spearman")#-0.18

plot((pfasds$`TSPC_142-62-1`),log(pfasds$HLH))
