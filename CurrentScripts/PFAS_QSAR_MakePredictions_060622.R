###############Use of Model to Predict DSSTox PFAS Set#########################
#Order of operations: First, in discrete models and continuous models and determine the 
#variables we need. Next, bring in the DSSTox PFAS chemicals, and merge with the Opera
#predictors for each chemical. Next, merge in the similarity data. Next, merge in 
#the kidney data since it has species specific info, and lastly, merge in data from Han, which 
#includes a species column, and includes humans, rats, and cattle. I will also extend QSAR predications to 
#dog, which has kidney measurements in Oliver et al. 1968.  

# Clear the workspace:
rm(list=ls())

# Change to the shared drive:
#setwd("L:/Lab/NCCT_ExpoCast/ExpoCast2022/Dawson_PFAS_HALFLIFE/PFAS_HL_QSAR_2021")

# Specify which RData files we are working with:
modelbuildsuff="060622" 
trainingsetsuff="JFW060522"
writesuff="JFW060622" #This is the suffix for ongoing work. 

# Random number geenerator seed:
seed <- "12345"

#Want to extend the model to pigs which we have some ability to do. 

packages=c("readxl","MLmetrics","stringr","scales","ggplot2")
sapply(packages, require,character.only=TRUE) #Note, the "character.only" argument is necessary her


#NOte, this loads in the the DSSTox PFAS list PFASStructV3(on the CompTox Dashboard), with 8163 chemicals. As of 06/10/21, we
#only  have COC and similarity data for the 6648 chemicals loaded in PFASSTructV2. In order to get
#predictions for the 2000 or so additional chemials, I need to get ToxPrint predictions for COC, and 
#the Python script for similarity values. 

#Pre-processing
PFASdata=read_excel("PFAS_Catalog/DSSTox-PFAS-MasterList.xlsx",sheet=2)
whitepaperchems <- read_excel("PFAS_Catalog/pfas_catalog 2021-08-23.xlsx")
qsar.smiles <- read.csv("PFAS_Catalog/final_list_OPPT_QSAR-ready_smi.smi",sep="\t",header=FALSE)
colnames(qsar.smiles) <- c("qsar-smiles","dtxsid")
whitepaperchems <- merge(whitepaperchems,qsar.smiles,all.x=TRUE)
whitepaperchems$INCHIKEY <- NA
whitepaperchems$INCHI_STRING <- NA

whitepaperchems <- whitepaperchems[,c(
  "dtxsid",
  "name",
  "casrn",
  "qsar-smiles",
  "INCHI_STRING",
  "INCHIKEY",
  "name",
  "formula",
  "mw")]
PFASdata <- PFASdata[,3:11]
colnames(whitepaperchems) <- colnames(PFASdata)

# The SMILES from "whitepaperchems" are more likely to be QSAR-ready:
dim(PFASdata)
PFASdata <- subset(PFASdata, !(DTXSID %in%whitepaperchems$DTXSID))
dim(PFASdata)
# Merge the two listS:
PFASdata <- rbind(PFASdata, whitepaperchems)
dim(PFASdata)
# Remove any entries with NA SMILES:
PFASdata <- subset(PFASdata, !is.na(QSAR_READY_SMILES))
dim(PFASdata)
# Remove dupicated ID's:
PFASdata <- subset(PFASdata, !duplicated(DTXSID))
dim(PFASdata)

saveRDS(PFASdata, file="DSSToxList_For_Example.rds")
write.table(PFASdata[,c("QSAR_READY_SMILES","DTXSID")], file=
  paste("DSSToxPFAS_DTXSID_SMILES_",
  writesuff,
  ".smi", sep=""), 
  row.names = FALSE, 
  quote=FALSE, 
  col.names = FALSE, 
  sep="\t"  )

 #Output to Opera
write.table(PFASdata[,c("QSAR_READY_SMILES","DTXSID")], paste("DSSToxPFAS_DTXSID_SMILES_",
  writesuff,
  ".smi", sep=""), row.names = FALSE, quote=FALSE, col.names = FALSE, sep="\t"  )


#1. Bring in variables from models 
load(paste0("PFAS_Training_Set_Endo_Disc_HLH_unscaled_",modelbuildsuff,".RData", sep=""))
# Load pfasdsrd for name index:
pfasdsredHLH=pfasdsred

#Lists of non-categorical variables that ended up in the training set; not categorical variables, including species, sex, and get added
AllVars <- names(pfasdsredHLH)

#2. Bring in DSSTox chemicals with Average Mass from dashboard
PFASdata <- subset(PFASdata, select=c("DTXSID", "PREFERRED_NAME", "CASRN"))

#3. Merge with Opera predictors
Opera <- read.csv(paste0("PFAS_Catalog/DSSToxPFAS_DTXSID_SMILES_Oct21-1-smi_OPERA2.7Pred.csv"))

##Reduce to only variable columns of interest, and also pull AD info 
names(Opera)

#PUll out molecular weight and name the same as the training data
AVERAGE_MASS=Opera$MolWeight

Operapreds=names(Opera)[unlist(lapply(AllVars, function(x) {grep(x,names(Opera))} ))] #Find Opera columns
#[1] "LogP_pred"   "LogVP_pred"  "LogWS_pred"  "LogD74_pred"


#Can't pull AD's automatically because they aren't listed in a straight forwrd manner; just have to code in manualy
Operacols=c("MoleculeID","MolWeight", Operapreds, "AD_LogP", "AD_VP", "AD_WS","AD_KOA")

##Remove logp range column 
Operacols=Operacols[-c(which(Operacols=="LogP_predRange"))]

Opera=Opera[,c(which(names(Opera)=="RowID"),which(names(Opera)%in%Operacols))]

#Set name of average mass to same as model 
names(Opera)[which(names(Opera)=="MolWeight")]="AVERAGE_MASS"


##Merge with PFASdata
PFASdata1=merge(PFASdata, Opera, by.x="DTXSID", by.y="MoleculeID", all.x=TRUE)



####Check against training data
# Load pfads:
load(paste0("PFAS_11Chemicals_QSARdataset_",trainingsetsuff,".RData", sep="")) #loads pfasds
load(paste0("Training_Data_Reorder_For_Prediction_",modelbuildsuff,".RData")) #pfasdsorder
pfasds=pfasds[pfasdsreorder,]
head(pfasds)
#Select out predictor columns
idcols=which(names(pfasds)%in%c("DTXSID", "CASRN", "PREFERRED_NAME"))
pfasdsred1=data.frame(pfasds[,idcols], pfasdsred)




ps=PFASdata1[which(PFASdata1$DTXSID%in%pfasdsred1$DTXSID),]
pfasdsred1=pfasdsred1[pfasdsred1$DTXSID=="DTXSID3031862",]
ps=ps[ps$DTXSID=="DTXSID3031862",]
pfasdsred1
ps 

pfasdsred1=data.frame(pfasds[,idcols], pfasdsred)
ps=PFASdata1[which(PFASdata1$DTXSID%in%pfasdsred1$DTXSID),]
pfasdsred1=pfasdsred1[pfasdsred1$DTXSID=="DTXSID1037303",]
ps=ps[ps$DTXSID=="DTXSID1037303",]
pfasdsred1
ps 

pfasdsred1=data.frame(pfasds[,idcols], pfasdsred)
ps=PFASdata1[which(PFASdata1$DTXSID%in%pfasdsred1$DTXSID),]
pfasdsred1=pfasdsred1[pfasdsred1$DTXSID=="DTXSID7040150",]
ps=ps[ps$DTXSID=="DTXSID7040150",]
pfasdsred1
ps 

"DTXSID7040150"

#Good, have very close values for unscaled opera predictors and average mass to training set. 

#####
#Now, merge with similarity data from Prachi 
Endoinfo<-read.csv("Predictors/EndoSubset-EPA_PFAS_DSSTox-Similarity.csv")

# Extra chemicals for Richard Judson:
endo <- read.csv("Predictors/EndoSubset-EPA_PFAS_DSSTox-Similarity_final_list_OPPT.csv")
#endo.melt <- melt(endo,id=c("CASRN","ENDO.DTXSID..EPA.PFAS."))
#endo.cast <- dcast(subset(endo.melt,variable=="Tanimoto.Score..PubChem."),ENDO.DTXSID..EPA.PFAS.~CASRN)
dim(endo)
endo <- subset(endo, !(ENDO.DTXSID..EPA.PFAS. %in% unique(Endoinfo$ENDO.DTXSID..EPA.PFAS.)))
dim(endo)

dim(Endoinfo)
Endoinfo <- rbind(Endoinfo,endo) 
dim(Endoinfo)
Endoinfo <- subset(Endoinfo,!duplicated(Endoinfo))
dim(Endoinfo)



###Grab the similarity values of the chemicals that are in the model; back calculate from abbreviations 
SimVarsPC=AllVars[which(grepl('TSPC', AllVars))]
####strip-off TS_ and replace "." with "-"
SimVarsPC=substring(SimVarsPC, 6)
SimVarsPC=str_replace_all(SimVarsPC, "[.]", "-")

###Next, pull similarity data from EndoinfoTotal1
Endosubset=Endoinfo[which(Endoinfo$CASRN%in%SimVarsPC),]
length(unique(Endosubset$ENDO.DTXSID..EPA.PFAS.)) #8384


###Fill in similarity values for each chemical in the model 
PFASdata2=PFASdata1
threshold=0.9                
for(i in SimVarsPC){
  sub=subset(Endosubset, CASRN==i, select=c("CASRN","ENDO.DTXSID..EPA.PFAS.", "Tanimoto.Score..PubChem."))
  names(sub)[2:3]=c("DTXSID", paste("TSPC_",i,sep=""))                         
  sub=sub[,-1]
  sub[,2]=ifelse(is.na(sub[,2]),NA,ifelse(sub[,2]>=threshold,1,0))
  PFASdata2=merge(PFASdata2, sub, all.x=TRUE,all.y=FALSE)
}
names(PFASdata2)
length(unique(PFASdata2$DTXSID))
length(which(is.na(PFASdata2$`TSPC_142-62-1`)==TRUE))

###Next add max similarity: Note, this would only be a max within the subset of available chemicals, so only 30 or so; this could be different if all chemicals were added
##This takes a bit of time because there are so many chemicals
#for(i in 1:length(PFASdata2$CASRN)){
#  a=subset(Endoinfo, ENDO.DTXSID..EPA.PFAS. == PFASdata2$DTXSID[i])
#  maxsim=max(a$Tanimoto.Score..PubChem.)
#  PFASdata2$MaxEndoPC[i]=maxsim  
#}

#PFASdata2
#simthreshold=0.9
#Endodisc=pfasds[, Endocols]
#range(apply(Endodisc, 2, FUN="max"))
#Endodisc=ifelse(Endodisc>=simthreshold, 1,0)

#Note that for some PFAS, there are no similar chemicals 
# could do this with an apply function too; for future reference
#sapply(PFASdata2$CASRN, function(x){
#  a=subset(Endoinfo, ENDO.DTXSID..EPA.PFAS. == x)
#  maxsim=max(a$Tanimoto.Score..PubChem.)
#  return(maxsim)
#})

###When there are no similar chemicals, this is set to zero. 
#PFASdata2$MaxEndoPC=ifelse(is.infinite(PFASdata2$MaxEndoPC), 0, PFASdata2$MaxEndoPC)

#########Check against training data
pfasdsred1=data.frame(pfasds[,idcols], pfasdsred)
ps=PFASdata2[which(PFASdata2$DTXSID%in%pfasdsred1$DTXSID),]
pfasdsred1=pfasdsred1[c(pfasdsred1$DTXSID=="DTXSID3031860" | pfasdsred1$DTXSID=="DTXSID3031864"),]
ps=ps[c(ps$DTXSID=="DTXSID3031860" | ps$DTXSID=="DTXSID3031864"),]
pfasdsred1[1:2,]
ps #Good, maintaining similar descretized similarity values  


#########
#Next, merge in kidney data to create slots for different species
###Adding in Sex to Kidney dataset 
load(file="Predictors/Kidney_Predictions_MultSpecies_DED032921.RData")
load(paste0("Kidney_Descriptors_in_Training_Set_",trainingsetsuff,".RData"))


kidney=dplyr::rename(kidney, "Species"="Mammal")
#"KW_BW_ratio"        "GlomTotSA_KW_ratio"
#"ProxTubLen"         "ProxTubVol"
kidney=subset(kidney,select=c("Species", "Type",kd))
kidney=kidney[-c(which(kidney$Species=="Whale" | kidney$Species=="Elephant")),]
#Parse down to the four parameters of interest
#Adding in categories for species, sex, and dosing approach here 
#Merging Kidney dataset into PFAS dataset
Sex=c("Male", "Female")
#Species=c("Human", "Monkey", "Rat", "Mouse")
DosingAdj=c("IV", "Oral", "Other")
Species=c("Human", "Monkey", "Rat", "Mouse", "Cattle", "Dog", "Rabbit", "Horse", "Chicken")
CatTab=expand.grid("Sex"=Sex, "Species"=Species, "DosingAdj"=DosingAdj)
kidney=merge(kidney, CatTab, by="Species", all=TRUE)

PFASdata3=PFASdata2
#Adding in Species and Sex to PFAS list
PFASdata3=merge(PFASdata3, Species)
PFASdata3=dplyr::rename(PFASdata3, "Species"="y")

PFASdata3=merge(PFASdata3, kidney, by="Species", all=TRUE)
dim(PFASdata3)
length(unique(PFASdata3$DTXSID))#8163 #9898



#########Check against training data
catcols=which(names(pfasds)%in%c("Species", "Sex", "Type", "DosingAdj"))
pfasdsred1=data.frame(pfasds[,idcols], pfasds[,catcols],pfasdsred)
ps=PFASdata3[which(PFASdata3$DTXSID%in%pfasdsred1$DTXSID),]
pfasdsred1=pfasdsred1[c(pfasdsred1$Species=="Rat" & pfasdsred1$Sex=="Male" & pfasdsred1$DosingAdj=="Other"),]
ps=ps[c(ps$Species=="Rat" & ps$Sex=="Male" & ps$DosingAdj=="Other"),]
ps
pfasdsred1[1:2,]

#Good, the rat and human data is the same  


#########

##Bring in data from Zhang)
#Z=read.xlsx("Predictors/Zhang et al.2013_PFAS_Protein_Dissociation Constants.xlsx")
#Z$Species="Human"
#Z=Z[-c(which(names(Z)=="ChemicalName"))]
#PFASdata3[which(PFASdata3$CASRN%in%Z$CASRN==FALSE),1:4] #chemicals in dataset not in Zhang
#dim(PFASdata3)
##De-salting some chemicals
##PFOS
##Note: Perfluoropentane sulfonic acid is not available
##    : Perfuoroheptansulfonic acid is not available
#fillin=Z[which(Z$ChemicalAbbreviation=="PFOS"),which(names(Z)=="CASRN")]
#Z[which(Z$CASRN==fillin),which(names(Z)=="CASRN")]="1763-23-1"
#fillin=Z[which(Z$ChemicalAbbreviation=="PFHxS"),which(names(Z)=="CASRN")]
#Z[which(Z$CASRN==fillin),which(names(Z)=="CASRN")]<-"355-46-4" 
#fillin=Z[which(Z$ChemicalAbbreviation=="PFBS"),which(names(Z)=="CASRN")]
#Z[which(Z$CASRN==fillin),which(names(Z)=="CASRN")]<-"375-73-5" 
##Z=subset(Z, select=c("CASRN", "Species", "Kd_hL_FABP"))
##Changing non-detects to 0
##Z$Kd_hL_FABP=as.numeric(ifelse(Z$Kd_hL_FABP=="ND", 0, Z$Kd_hL_FABP))
#
##PFASdata4=merge(PFASdata3, Z[,c("CASRN", "Kd_hL_FABP")],by="CASRN", all.x=TRUE, all.y=FALSE) #only merge by DTXSID; the data is only in humans here. 
#length(unique(PFASdata4$DTXSID))#9898

PFASdata4 <- PFASdata3
names(PFASdata4)

#Replace "-" with "."
names(PFASdata4)<-str_replace_all(names(PFASdata4), "-", ".") #change - to .

#Read in add in COC
COC=read.csv("Predictors/DSSTox_PFAS_ToxPrints_COC_aliphatic.csv") 
dim(COC)
ToxPrints <- read.csv("PFAS_Catalog/ToxPrints-COC.csv")
ToxPrints <- ToxPrints[,c("DTXSID","bond.COC_ether_aliphatic")]
dim(ToxPrints)
colnames(ToxPrints) <- colnames(COC)
COC <- rbind(COC,ToxPrints)
dim(COC)
COC <- subset(COC,!is.na(COC_aliphatic))
dim(COC)
COC <- subset(COC,!duplicated(DTXSID))
dim(COC)
save(COC, file="PFAS_COC_DSSTox.RData")


PFASdata5=merge(PFASdata4, COC, by="DTXSID")

length(unique(PFASdata5$DTXSID))#6509 #8209

names(PFASdata5)

#Read in Data from Han
H=read_excel("Predictors/Han_et al.2012_AlbuminProteinBinding_Data.xlsx")
H$Reference_Ka_SerAlb="Han 2011"
H[which(H$PREFERRED_NAME=="Perfluorobutanesulfonate"), which(names(H)=="PREFERRED_NAME")]<-"Perfluorobutanesulfonic acid"

PFASdata6=merge(PFASdata5, H[,c("PREFERRED_NAME", "Species", "Ka_perM_SerAlb_Han")], by=c("PREFERRED_NAME", "Species"), all.x=TRUE, all.y=FALSE)
dim(PFASdata6)
length(unique(PFASdata6$DTXSID))#6509 #8297

####Check against training data
pfasdsred1=data.frame(pfasds[,idcols], pfasds[,catcols],pfasdsred)
ps=PFASdata6[which(PFASdata6$DTXSID%in%pfasdsred1$DTXSID),]
pfasdsred1=pfasdsred1[c(pfasdsred1$Species=="Human" & pfasdsred1$Sex=="Male" & pfasdsred1$DosingAdj=="Other" & pfasdsred1$DTXSID=="DTXSID3031862" ),]
ps=ps[c(ps$Species=="Human" & ps$Sex=="Male" & ps$DosingAdj=="Other" & ps$DTXSID=="DTXSID3031862" ),]
ps
pfasdsred1


######


AllVars%in%names(PFASdata6) #All continous variables are present 
IDcols=PFASdata6[,c("CASRN", "DTXSID", "Species")] #
OperaADcols=PFASdata6[,grep("AD", names(PFASdata6))]
length(unique(IDcols$DTXSID))#6509 #8297

PFASdata7=PFASdata6[,which(names(PFASdata6)%in%AllVars)] #only continous variables 
names(PFASdata7)



######
#Grab categorical sets and adjust the numerics to be on the same scale as the training data 
PFASdata7cat=subset(PFASdata6, select=c("Species","Sex",  "Type", "DosingAdj"))

#PFASdata6 is all numeric
sapply(PFASdata7, class)


######Fill in Average Values: will do at end for all numeric values 
#Find extent of missing data
missdat=apply(PFASdata7, 1, function(x){length(which(is.na(x)==TRUE))})
PFASdata7[343,]

hist(missdat)#Looks like the vast majority are missing at least 1 parameter,which is the FABP parameter. Probably should
#determine a better way to group those then. Instituting an adhoc decision to remove all chemicals in which >2 parameters
#have to be imputed via an average. 

removemissing=which(missdat>2) #Will need to remove from IDcols and categorical cols as well 
PFASdata7=PFASdata7[-removemissing,]

#Load averages of training set 
load(paste0("Mean_SD_Trainingset_",modelbuildsuff,".RData"))

#Input averages of trainset for the rest for rest
for(i in 1:ncol(PFASdata7)){
  meancol=train_msd[1,which(names(train_msd)==names(PFASdata7)[i])]
  sub=which(is.na(PFASdata7[,i])==TRUE)  
  if(length(sub)>0){
    PFASdata7[sub,i]=meancol
  }}

which(is.na(PFASdata7)) #Should be zero 
#Make sure each cat is a factor 


  #Make sure categorical variables are factors,and remove chemicals with missing data
  PFASdata7cat$Sex<-as.factor(PFASdata7cat$Sex)
  PFASdata7cat$Type<-as.factor(PFASdata7cat$Type)
  PFASdata7cat$DosingAdj<-as.factor(PFASdata7cat$DosingAdj)
  PFASdata7cat=PFASdata7cat[-removemissing,]
  IDcols=IDcols[-removemissing,] 
  OperaADcols=OperaADcols[-removemissing,]
  length(unique(IDcols$DTXSID)) #6603; so this mainly results from the ad hoc decision to remove chemicals with >2 missing columns
  

  ###Scaling and centering on the training data mean and sd
  PFASdata7sc=PFASdata7
  for(i in 1:ncol(PFASdata7sc)){
    meancol=train_msd[1,which(names(train_msd)==names(PFASdata7sc)[i])]
    sdcol=train_msd[2,which(names(train_msd)==names(PFASdata7sc)[i])]
    PFASdata7sc[,i]=(PFASdata7sc[,i]-meancol)/sdcol  
  }


#Total Dataset to Predict
pfastopred=cbind(PFASdata7cat, PFASdata7sc)

#######Training data check
load(file=paste0("PFAS_QSAR_Predictors_Predictions_", modelbuildsuff,".RData"))
pfasdsred1=data.frame(pfasds[,idcols], pfasds[,catcols],pfasdsredsc_ds)
ps=cbind(IDcols,PFASdata7cat,PFASdata7sc)
ps=ps[which(ps$DTXSID%in%pfasdsred1$DTXSID),]
pfasdsred1=pfasdsred1[c(pfasdsred1$Species=="Human" & pfasdsred1$Sex=="Male" & pfasdsred1$DosingAdj=="Other"),]
ps=ps[c(ps$Species=="Human" & ps$Sex=="Male" & ps$DosingAdj=="Other"),]
#Note; at least estimate for WS is very different; DTXSID3031860. So, I need to get new predictions
#from another source the model. 

set=NULL
for(i in unique(ps$DTXSID)){
  sub=ps[ps$DTXSID==i,which(names(ps)%in%c("DTXSID","TSPC_107.92.6",  "AVERAGE_MASS", "LogWS_pred", "LogP_pred", "LogVP_pred"))]
  sub2=pfasdsred1[pfasdsred1$DTXSID==i,which(names(pfasdsred1)%in%c("DTXSID","TSPC_107.92.6", "AVERAGE_MASS", "LogWS_pred", "LogP_pred", "LogVP_pred"))]
  order1=order(names(sub))
  sub=sub[,order1]
  order2=order(names(sub2))
  sub2=sub2[,order2]
  set=rbind(set,sub,sub2)
  }

set=rbind(ps[,c("DTXSID","TSPC_107.92.6", "AVERAGE_MASS", "LogWS_pred", "LogP_pred")],pfasdsred1[,c("DTXSID","TSPC_107.92.6", "AVERAGE_MASS", "LogWS_pred","LogP_pred")])
#set; very values for average mass and LogWS_pred


#Make model predictions 
#Regression model 
#load("RegressionModel_EndoSimDisc_HLH_DED061621.RData")
#set.seed(seed)
#pfastopred_regFull=predict(regmod, newdata = pfastopred)

#load("RegressionModel_PostRFE_HLH_DED061621.RData")
#set.seed(seed)
#pfastopred_regRed=predict(regmodred, newdata = pfastopred)

#Classification model predictions 
load(paste0("ClassificationModel_4Bin_EndoSimDisc_HLH_",modelbuildsuff,".RData"))
set.seed(seed)
pfastopred_classFull=predict(classmod4, newdata = pfastopred) #Need to re-arrange the order


length(unique(IDcols$CASRN))
length(unique(IDcols$DTXSID))

#Making complete set for Tox21 
completedataset=data.frame(IDcols[,c("CASRN", "DTXSID")], OperaADcols, pfastopred,  "ClassPredFull"=pfastopred_classFull)

#Comaring RegPredFull with RegPredRed
#MAE(completedataset$RegPredFull, completedataset$RegPredRed)
#5442
#RMSE(completedataset$RegPredFull, completedataset$RegPredRed)
#7060 #pretty close agreement between the complete and reduced model 


save(completedataset, file=paste("PFAS_HLpreds_HLH_HLHLSsc_DSSToxchemicals_coarse_imputation_",writesuff,".RData", sep="")) 
#load(paste("PFAS_HLpreds_HLH_HLHLSsc_DSSToxchemicals_coarse_imputation_",writesuff,".RData", sep=""))
length(unique(completedataset$CASRN))
# Load medians of the four bins:
load(paste("Median_HLH_per_Bin_4_",
  modelbuildsuff,
  ".RData",sep=""))
agtabHLH
# Convert to elimination rate (1/h):
agtabHLH$kelim.h <- signif(log(2)/agtabHLH$HLH,3)
# Median Vd from data set is 0.202 L/kg BW
agtabHLH$CLtot.Lpkgbwpday <- signif(agtabHLH$kelim.h*0.202*24,3)
# PFOA health advisory CL 0.00014 L/kg BW/day #predicted by model here to be 0.000115
# PFOS health advisory CL 0.000081 L/kg BW/day
agtabHLH$CLtot.Lpkgbwpday
# 1 mg/kg BW/day dose rate -> Css (mg/L)
agtabHLH$Css.mgpL <- signif(1/agtabHLH$CLtot.Lpkgbwpday,3)

completedataset$CLtot.Lpkgbwpday <- NA
completedataset$Css.mgpL <- NA
for (this.bin in seq(1,4))
  for (this.col in c("CLtot.Lpkgbwpday","Css.mgpL"))
    completedataset[completedataset[,"ClassPredFull"]==this.bin, this.col] <-
      agtabHLH[agtabHLH[,"HLHBin4"]==this.bin, this.col]
      
  completedataset[completedataset$CASRN=="335-67-1" & completedataset$Species=="Human",]
  completedataset[completedataset$CASRN=="13252-13-6" & completedataset$Species=="Human",]
  

  

####Applicability Domain####
#To be consistent with Mansouri et al, if something is in the domain, Domain=1, if outside=0. 
#If if its in the domain=1, if outside the doma
#This is only using the method by Roy et al. 2015, using the descriptors that are scaled by those of the training set. So, 

#May need to find some other idea of AD because apparently everythign is out of the AD now
ADfunction=function(model, data){
  vars=names(model$trainingData)
  vars=vars[-c(length(vars))] #removes last column which is predictions
  cats=which(vars%in%c("Species", "Sex", "DosingAdj"))
  Vars1=vars
  if(length(cats>0)){
    Vars1=vars[-c(which(vars%in%c("Species", "Sex",  "DosingAdj")))]}
  a=data
  b=a[,c("CASRN", "DTXSID","Species", "Sex",  "DosingAdj")]
  a=a[,which(names(a)%in%Vars1)]
  a=abs(a)
  SImin=apply(a,1,min) #Find the min, max, sd, and mean
  SImax=apply(a,1,max)
  SImean=apply(a,1,mean)
  SIsd=apply(a,1,sd)
  SIdf=data.frame(data,SImin, SImax, SImean, SIsd)
  SIdf$SI90=SIdf$SImean + (1.28 * SIdf$SIsd)
  SIdf$Domain=1 #something is assumed to be inside the AD unless it it found to not be 
  
  #Based Algorithm specificed in Roy et al. 2005)
  SIdf$Domain=ifelse(SIdf$SImin>3, 0, ifelse((SIdf$SImean+(1.28*SIsd))>3,0,1))
  
  return(SIdf)}

CMFAD=ADfunction(model=classmod4,data=completedataset)
#RMFAD=ADfunction(model=regmod,data=completedataset)
#RMRAD=ADfunction(model=regmodred,data=completedataset)

names(CMFAD)[which(names(CMFAD)=="Domain")]="ClassModDomain"
#names(RMFAD)[which(names(RMFAD)=="Domain")]="RegFullModDomain"
#names(RMRAD)[which(names(RMRAD)=="Domain")]="RegRedModDomain"

names(CMFAD)[which(names(CMFAD)=="SI90")]="ClassModSI90"
#names(RMFAD)[which(names(RMFAD)=="SI90")]="RegFullModSI90"
#names(RMRAD)[which(names(RMRAD)=="SI90")]="RegRedModSI90"

#DPFASwAD=data.frame(CMFAD[,-c(which(names(CMFAD)%in%c("SImin",  "SImax", "SImean","SIsd")))], RMFAD[,c("RegFullModSI90", "RegFullModDomain")], RMRAD[,c("RegRedModSI90", "RegRedModDomain")]) 
DPFASwAD=data.frame(CMFAD[,-c(which(names(CMFAD)%in%c("SImin",  "SImax", "SImean","SIsd")))]) 

#Pull in life span information prior to saving
#library(openxlsx)
#LF=as.data.frame(read_excel("Predictors/LifeSpan_from_Risa_05_25_20.xlsx"))
# 
# unique(completedataset$Species)
# 
# PerLS=ifelse(DPFASwAD$Species=="Human", LF[which(LF$Species=="Human"), "HrLifeSpan"],  
#              ifelse(DPFASwAD$Species=="Monkey", LF[which(LF$Species=="Monkey"), "HrLifeSpan"],
#                     ifelse(DPFASwAD$Species=="Rat", LF[which(LF$Species=="Rat"), "HrLifeSpan"],
#                            ifelse(DPFASwAD$Species=="Mouse", LF[which(LF$Species=="Mouse"), "HrLifeSpan"],
#                                   ifelse(DPFASwAD$Species=="Dog", LF[which(LF$Species=="Dog"), "HrLifeSpan"], NA)))))
# 
# 
# DPFASwAD$PerLS=PerLS

DPFASwAD$AMAD <- apply(DPFASwAD[,Operacols[regexpr("AD",Operacols)>-1]],1,prod)

save(DPFASwAD, file=paste0("Tox21_AllMods_ADindicate_", writesuff,".RData"))
write.csv(DPFASwAD, file=paste0("Tox21_AllMods_ADindicate_", writesuff,".csv"))
#load(file=paste0("Tox21_AllMods_ADindicate_", writesuff,".RData"))
#DPFASwAD=read.csv(paste0("Tox21_AllMods_ADindicate_", writesuff,".csv"))
DPFASwAD[DPFASwAD$DTXSID=="DTXSID8047553" & DPFASwAD$Species=="Human",]

  DPFASwAD[DPFASwAD$CASRN=="335-67-1",]

#  PFOS=DTXSID3031865
#  PFOA=DTXSID40892486
#Gex=DTXSID70880215 

subset1=DPFASwAD[DPFASwAD$Sex=="Male" & DPFASwAD$DosingAdj=="Oral",]
subset1$Count=ifelse(is.na(subset1$Css.mgpL), NA,1)
aggregate(Count~Species, data=subset1, FUN="sum")  
  

#Note, only 1917 chemicals are listed as being in the domain of the reduced regresson model 
length(unique(DPFASwAD[DPFASwAD$ClassModDomain==1, "CASRN"])) #6117
#length(unique(DPFASwAD[DPFASwAD$RegFullModDomain==1, "CASRN"])) #3368
#length(unique(DPFASwAD[DPFASwAD$RegRedModDomain==1, "CASRN"])) #3368

###Domain tables
 dom=DPFASwAD
 dom=subset(dom, Species%in%c("Human","Monkey", "Rat", "Mouse", "Dog"))
 dom=subset(dom, Sex=="Female" & DosingAdj=="Other")
 CDtab=aggregate(ClassModDomain~DosingAdj+Sex + Species,data=dom, FUN = "sum")
# RRtab=aggregate(RegRedModDomain~Species,data=dom, FUN = "sum")
# RFtab=aggregate(RegFullModDomain~Species,data=dom, FUN = "sum")
 dom2=subset(dom, AD_LogP==1 & AD_VP==1 & AD_WS==1)
 CDtab2=aggregate(ClassModDomain~DosingAdj+Sex + Species,data=dom2, FUN = "sum")
# RRtab2=aggregate(RegRedModDomain~Species,data=dom2, FUN = "sum")
# 
# 
 domtab1=data.frame("Model"="Classification",CDtab,"All_Domains"=CDtab2[,"ClassModDomain"])
# domtab2=data.frame("Model"="Regression_Reduced", RRtab, "All_Domains"=RRtab2[,"RegRedModDomain"])
# names(domtab1)[3]="ModelDomain"
# names(domtab2)[3]="ModelDomain"
# domtab=rbind(domtab1, domtab2)
# domtab
# save()
####

 
 

#Modellist=list("ClassModFull"=classmod4,  "RegModFull"=regmod, "RegModRed"=regmodred)
#save(Modellist, file=paste0("Class_Reg_Full_Reduced_MOdels_",writesuff,".RData"))
#load(paste0("Class_Reg_Full_Reduced_MOdels_",writesuff,".RData"))
#This will serve as an interesting discussion of in versus out of the AD, and how it differs by species and by model. 

###Plots ######
#Extract all modvars, including catetorical (sex and DosingAdj)

Modvars=row.names(varImp(classmod4)$importance)
###Plot out some predictions: for numbers of chemicals per bin per species
#Plan for plots:
#1. Make 1 plot apiece for each species for of all of the models(3) for the raw distribution of predicted half-lives
#For regression models, this should amount to a box and whisker plot. For the classification models, this should
#amount to the stacked bar chart. Either way, we should have specise on the X axis and distribution of half-lives on the y, probably by log10.  
#So, for 5 species(human,monkey, rat, mouse, and cow(to showcase extension))* 3 models)
#2. Make 1 plot apiece for each of the three models of the distributions scaled by species lifespans. These will all be box and whisker plots
#ranging from 0 to 1. 

DPFASwAD$Count=1

#Gen X
DPFASwAD[which(DPFASwAD$CASRN=="13252-13-6" & DPFASwAD$Species=="Human"),] 
#Classification Model putting it in Bin 2 for humans; good
#Regression model putting it at 7615 hours for the full model and 12019 hours for the reduced model 

#Proportions of predictions for humans
human=DPFASwAD[DPFASwAD$Species=="Human" & DPFASwAD$Sex=="Male" & DPFASwAD$DosingAdj=="Other",]
table(human$ClassPredFull)/sum(table(human$ClassModDomain))

#########Compare with training set chemicals
ps=DPFASwAD[which(DPFASwAD$DTXSID%in%pfasdsredsc_ds$DTXSID),]
ps=ps[ps$Species=="Human" & ps$Sex=="Male" & ps$DosingAdj=="Other",]
dshuman=pfasdsredsc_ds
dshuman=dshuman[c(dshuman$Species=="Human" & dshuman$Sex=="Male"),]

merge(ps[c("DTXSID", "ClassPredFull")], dshuman[c("DTXSID", "ClassPredFull")], by="DTXSID")

#Note: while the classification model predictiosn for the training chemicals is the same for the 
#the training model and when applied here, three of the chemicals have different regression model results. 
#the predictors are the same, so this must result from model instability issues. 
chems=c("DTXSID1037303", "DTXSID4059916")

sub1=subset(dshuman, c(DTXSID%in%chems))
sub2=subset(ps, c(DTXSID%in%chems))
merge(sub1, sub2,  by="DTXSID")

#####

seq=seq(1,1000,1)

b=3*seq^(2/3)
plot(seq,b)


# Chemicals in domain:
length(unique(subset(DPFASwAD,ClassModDomain==1)$DTXSID))

# PFUnDA:
subset(DPFASwAD,DTXSID=="DTXSID8047553" & Species=="Human")


in.domain <- subset(DPFASwAD,ClassModDomain==1&Species%in%c("Mouse","Rat","Monkey","Human"))
in.domain$SpeciesSex=paste0(in.domain$Species, "_", in.domain$Sex)
in.domain$SpeciesDose=paste0(in.domain$Species, "_", in.domain$DosingAdj)
in.domain$SpeciesDose=paste0(in.domain$Species, "_", in.domain$DosingAdj)

library(ggplot2)

#Now, plot out the distribution for the alkyl halides, carboxylic acids and sulfonic acid derivatives
#Insert male and female symbols
######Prepare plotting features
#Labels
library(showtext)
female = intToUtf8(9792)
male = intToUtf8(9794)
IV="Iv"
Oral="Or"
Other="Ot"
####Numbers of chemicals by bin
namelabelsSex=c(paste0("Dog:",female), paste0("Dog:",male), paste0("Human:",female), paste0("Human:",male), paste0("Monkey:",female), paste0("Monkey:",male), paste0("Mouse:",female), paste0("Mouse:",male), paste0("Rat:",female), paste0("Rat:",male))
namelabelsDose=c(paste0("D:", IV), paste0("D:", Oral), paste0("D:", Other), paste0("H:",IV), paste0("H:",Oral),paste0("H:",Other), paste0("Mn:",IV), paste0("Mn:",Oral),paste0("Mn:",Other), paste0("Mo:",IV), paste0("Mo:",Oral),paste0("Mo:",Other), paste0("Rat:",IV), paste0("Rat:",Oral),paste0("Rat:",Other))
#########




pfascountCFSS=aggregate(Count~ ClassPredFull +DosingAdj+SpeciesSex , data=in.domain , FUN="sum")
Plot_Species_by_Sex=ggplot(data=in.domain[in.domain$DosingAdj=="IV",], aes(x=SpeciesSex, y=Count, fill=ClassPredFull)) +
  geom_bar(stat="identity")+
  xlab("Species/Sex")+
  ylab("Number of PFAS chemicals")+
  ggtitle(paste("A) Toxicokinetic Half-Life of ",length(unique(in.domain$DTXSID)), " PFAS Chemicals in Model Domain",sep=""))+
  scale_fill_discrete(name="Serum Half Life", labels = c("<0.5 Day", "<Week", "< 2 Months", "> 2 Months"))+   
  theme(
    plot.title=element_text(size =20), 
    # plot.subtitle = element_text(size=17),
    axis.text.x =element_text(size=15, vjust=0.3, hjust=0.5),
    axis.text.y=element_text(size=15),
    axis.title=element_text(size=20),
    legend.position="bottom",
    legend.text=element_text(size=15),
    legend.title=element_text(size=20))
Plot_Species_by_Sex=
  Plot_Species_by_Sex+  
  scale_x_discrete(labels= c("Dog_Female" = namelabelsSex[1], "Dog_Male"=namelabelsSex[2], "Human_Female"=namelabelsSex[3],
                             "Human_Male"=namelabelsSex[4],"Monkey_Female"=namelabelsSex[5], "Monkey_Male"=namelabelsSex[6], 
                             "Mouse_Female"=namelabelsSex[7], "Mouse_Male"=namelabelsSex[8], "Rat_Female"=namelabelsSex[9], 
                             "Rat_Male"=namelabelsSex[10]))

Plot_Species_by_Sex
png(paste("Figures/THalfPred_Species_by_Sex_",writesuff,".png",sep=""), width=800, height=600)
Plot_Species_by_Sex
dev.off()

head(in.domain[,Operacols[regexpr("AD",Operacols)>-1]])

