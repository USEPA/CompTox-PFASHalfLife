###############Use of Model to Predict DSSTox PFAS Set#########################
#Order of operations: First, in discrete models and continuous models and determine the 
#variables we need. Next, bring in the DSSTox PFAS chemicals, and merge with the Opera
#predictors for each chemical. Next, merge in the similarity data. Next, merge in 
#the kidney data since it has species specific info, and lastly, merge in data from Han, which 
#includes a species column, and includes humans, rats, and cattle. I will also extend QSAR predications to 
#dog, which has kidney measurements in Oliver et al. 1968.  


writesuff="DED061621_2"

#Want to extend the model to pigs which we have some ability to do. 

packages=c("readxl","openxlsx", "gdata", "stringr","scales")
sapply(packages, require,character.only=TRUE) #Note, the "character.only" argument is necessary her


#NOte, this loads in the the DSSTox PFAS list PFASStructV3(on the CompTox Dashboard), with 8163 chemicals. As of 06/10/21, we
#only have COC and similarity data for the 6648 chemicals loaded in PFASSTructV2. In order to get
#predictions for the 2000 or so additional chemials, I need to get ToxPrint predictions for COC, and 
#the Python script for similarity values. 

#Pre-processing
PFASdata=read_xls("DSSTox_ComptoxDashboard_060921.xls")
#Output to Opera
PFASdata_for_Opera=subset(PFASdata, select=c("SMILES"))
write.table(PFASdata_for_Opera, "DSSToxPFAS_DTXSID_SMILES_060921.txt", row.names = FALSE, quote=FALSE, col.names = FALSE  )


#1. Bring in variables from models 
##LOad name index:
load(paste("PFAS_Training_Set_Endo_Disc_HLH_unscaledDED061621.RData", sep=""))

pfasdsredHLH=pfasdsred

#Lists of non-categorical variables that ended up in the training set; not categorical variables, including species, sex, and get added
AllVars=names(pfasdsredHLH)

#2. Bring in DSSTox chemicals with Average Mass from dashboard
#NOte: May need to de-salt or use only the acid/base conjugate compound
#can probably figure out later with Wambaugh#
PFASdata=subset(PFASdata, select=c("DTXSID", "PREFERRED_NAME", "CASRN"))
PFASdata$RowID=seq(1,nrow(PFASdata),1)

#Salts
head(PFASdata[which(grepl('oate', PFASdata$PREFERRED_NAME)),])
#Note, 6 salts here



#3. Merge with Opera predictors
Opera=read.csv(paste0("DSSTox_PFAS_061021/DSSToxPFAS_DTXSID_SMILES_060921-smi_OPERA2.6Pred.csv"))
MolRowID=strsplit(Opera$MoleculeID, "_")
MolRowID=as.numeric(sapply(MolRowID, "[[", 2))
Opera$RowID=MolRowID


##Reduce to only variable columns of interest, and also pull AD info 
which(names(Opera))

Operapreds=names(Opera)[unlist(lapply(AllVars, function(x) {grep(x,names(Opera))} ))] #Find Opera columns
[1] "LogP_pred"   "LogVP_pred"  "LogWS_pred"  "LogD74_pred"
#Can't pull AD's automatically because they aren't listed in a straight forwrd manner; just have to code in manualy
Operacols=c(Operapreds, "AD_LogP", "AD_VP", "AD_WS", "AD_LogD")
Opera=Opera[,c(which(names(Opera)=="RowID"),which(names(Opera)%in%Operacols))]

##Merge with PFASdata
PFASdata1=merge(PFASdata, Opera, by.x="RowID", all.x=TRUE)


####Check against training data

pfasdsred1=data.frame(pfasds[,idcols], pfasdsred)
ps=PFASdata1[which(PFASdata1$DTXSID%in%pfasdsred1$DTXSID),]
pfasdsred1=pfasdsred1[pfasdsred1$DTXSID=="DTXSID3031862",]
ps=ps[ps$DTXSID=="DTXSID3031862",]
pfasdsred1
ps #Good, have identical values for unscaled opera predictors and average mass to training set. 

#####


#Now, merge with similarity data from Prachi 
library(openxlsx)
library(stringr)
Endoinfo=read.csv("EndoSubset-EPA_PFAS_DSSTox-Similarity.csv")


###Grab the similarity values of the chemicals that are in the model; back calculate from abbreviations 
SimVarsPC=AllVars[which(grepl('TSPC', AllVars))]
####strip-off TS_ and replace "." with "-"
SimVarsPC=substring(SimVarsPC, 6)
SimVarsPC=str_replace_all(SimVarsPC, "[.]", "-")

###Next, pull similarity data from EndoinfoTotal1
Endosubset=Endoinfo[which(Endoinfo$CASRN%in%SimVarsPC),]
length(unique(Endosubset$ENDO.DTXSID..EPA.PFAS.)) #6433


###Fill in similarity values for each chemical in the model 
PFASdata2=PFASdata1
threshold=0.9
for(i in SimVarsPC){
  sub=subset(Endosubset, CASRN==i, select=c("CASRN","ENDO.DTXSID..EPA.PFAS.", "Tanimoto.Score..PubChem."))
  names(sub)[2:3]=c("DTXSID", paste("TSPC_",i,sep=""))
  sub=sub[,-1]
  sub[,2]=ifelse(sub[,2]>=threshold,1,0)
  PFASdata2=merge(PFASdata2, sub, all.x=TRUE,all.y=FALSE)
}
names(PFASdata2)
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
ps #Good, maintaining similar descritized similarity values  


#########
#Next, merge in kidney data to create slots for different species
###Adding in Sex to Kidney dataset 
load(file="Kidney_Predictions_MultSpecies_DED032921.RData")
load(paste0("Kidney_Descriptors_in_Training_Set_DED061621.RData"))


kidney=dplyr::rename(kidney, "Species"="Mammal")
#"KW_BW_ratio"        "GlomTotSA_KW_ratio"
"ProxTubLen"         "ProxTubVol"
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
length(unique(PFASdata3$DTXSID))#8163

#########Check against training data
pfasdsred1=data.frame(pfasds[,idcols], pfasds[,catcols],pfasdsred)
ps=PFASdata3[which(PFASdata3$DTXSID%in%pfasdsred1$DTXSID),]
pfasdsred1=pfasdsred1[c(pfasdsred1$Species=="Rat" & pfasdsred1$Sex=="Male" & pfasdsred1$DosingAdj=="Other"),]
ps=ps[c(ps$Species=="Rat" & ps$Sex=="Male" & ps$DosingAdj=="Other"),]
ps
pfasdsred1[1:2,]

#Good, the rat and human data is the same  


#########

#Bring in data from Zhang)
Z=read.xlsx("Zhang et al.2013_PFAS_Protein_Dissociation Constants.xlsx")
Z$Species="Human"
Z=Z[-c(which(names(Z)=="ChemicalName"))]
df3[which(df3$CASRN%in%Z$CASRN==FALSE),1:4] #chemicals in dataset not in Zhang
dim(df3)
#De-salting some chemicals
#PFOS
#Note: Perfluoropentane sulfonic acid is not available
#    : Perfuoroheptansulfonic acid is not available
fillin=Z[which(Z$ChemicalAbbreviation=="PFOS"),which(names(Z)=="CASRN")]
Z[which(Z$CASRN==fillin),which(names(Z)=="CASRN")]="1763-23-1"
fillin=Z[which(Z$ChemicalAbbreviation=="PFHxS"),which(names(Z)=="CASRN")]
Z[which(Z$CASRN==fillin),which(names(Z)=="CASRN")]<-"355-46-4" 
fillin=Z[which(Z$ChemicalAbbreviation=="PFBS"),which(names(Z)=="CASRN")]
Z[which(Z$CASRN==fillin),which(names(Z)=="CASRN")]<-"375-73-5" 
Z=subset(Z, select=c("CASRN", "Species", "Kd_hL_FABP"))
#Changing non-detects to 0
Z$Kd_hL_FABP=as.numeric(ifelse(Z$Kd_hL_FABP=="ND", 0, Z$Kd_hL_FABP))

PFASdata4=merge(PFASdata3, Z[,c("CASRN", "Kd_hL_FABP")],by="CASRN", all.x=TRUE, all.y=FALSE) #only merge by DTXSID; the data is only in humans here. 
length(unique(PFASdata4$DTXSID))#8163

names(PFASdata4)

#Replace "-" with "."
names(PFASdata4)<-str_replace_all(names(PFASdata4), "-", ".") #change - to .

#Read in add in COC
COC=read.csv("DSSTox_PFAS_ToxPrints_COC_aliphatic.csv") 
PFASdata5=merge(PFASdata4, COC, by="DTXSID")

length(unique(PFASdata5$DTXSID))#6509

names(PFASdata5)

#Read in Data from Han
H=read.xlsx("Han_et al.2012_AlbuminProteinBinding_Data.xlsx")
H$Reference_Ka_SerAlb="Han 2011"
H[which(H$PREFERRED_NAME=="Perfluorobutanesulfonate"), which(names(H)=="PREFERRED_NAME")]<-"Perfluorobutanesulfonic acid"

PFASdata6=merge(PFASdata5, H[,c("PREFERRED_NAME", "Species", "Ka_perM_SerAlb_Han")], by=c("PREFERRED_NAME", "Species"), all.x=TRUE, all.y=FALSE)
dim(PFASdata6)
length(unique(PFASdata6$DTXSID))#6509

####Check against training data
pfasdsred1=data.frame(pfasds[,idcols], pfasds[,catcols],pfasdsred)
ps=PFASdata6[which(PFASdata6$DTXSID%in%pfasdsred1$DTXSID),]
pfasdsred1=pfasdsred1[c(pfasdsred1$Species=="Human" & pfasdsred1$Sex=="Male" & pfasdsred1$DosingAdj=="Other" & pfasdsred1$DTXSID=="DTXSID3031862" ),]
ps=ps[c(ps$Species=="Human" & ps$Sex=="Male" & ps$DosingAdj=="Other" & ps$DTXSID=="DTXSID3031862" ),]
ps
pfasdsred1
set=rbind(ps[,c("DTXSID","TSPC_107.92.6", "Kd_hL_FABP", "LogD74_pred", "LogP_pred")],pfasdsred1[,c("DTXSID","TSPC_107.92.6", "Kd_hL_FABP", "LogD74_pred", "LogP_pred")])
set

######


AllVars%in%names(PFASdata6) #All continous variables are present 
IDcols=PFASdata6[,c("CASRN", "DTXSID", "Species")] #
OperaADcols=PFASdata6[,grep("AD", names(PFASdata6))]
length(unique(IDcols$DTXSID))#6509

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
load(paste0("Mean_SD_Trainingset_DED061621.RData"))

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
length(unique(IDcols$DTXSID)) #6181
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

pfasdsred1=data.frame(pfasds[,idcols], pfasds[,catcols],pfasdsredsc_ds)
ps=cbind(IDcols,PFASdata7cat,PFASdata7sc)
ps=ps[which(ps$DTXSID%in%pfasdsred1$DTXSID),]
pfasdsred1=pfasdsred1[c(pfasdsred1$Species=="Human" & pfasdsred1$Sex=="Male" & pfasdsred1$DosingAdj=="Other" & pfasdsred1$DTXSID=="DTXSID70880215" ),]
ps=ps[c(ps$Species=="Human" & ps$Sex=="Male" & ps$DosingAdj=="Other" & ps$DTXSID=="DTXSID70880215" ),]
set=rbind(ps[,c("DTXSID","TSPC_107.92.6", "Kd_hL_FABP", "LogD74_pred", "LogWS_pred")],pfasdsred1[,c("DTXSID","TSPC_107.92.6", "Kd_hL_FABP", "LogD74_pred", "LogWS_pred")])
set

######


#Make model predictions 
#Regression model 
set.seed(seed)
pfastopred_regFull=predict(regmod, newdata = pfastopred)
set.seed(seed)
pfastopred_regRed=predict(regmodred, newdata = pfastopred)

#Classification model 
set.seed(seed)
pfastopred_classFull=predict(classmod4, newdata = pfastopred) #Need to re-arrange the order


length(unique(IDcols$CASRN))
length(unique(IDcols$DTXSID))

#Making complete set for Tox21 
completedataset=data.frame(IDcols[,c("CASRN", "DTXSID")], OperaADcols, pfastopred,  "ClassPredFull"=pfastopred_classFull, "RegPredFull"=pfastopred_regFull, "RegPredRed"=pfastopred_regRed)

#Comaring RegPredFull with RegPredRed
MAE(completedataset$RegPredFull, completedataset$RegPredRed)
#5579.881
RMSE(completedataset$RegPredFull, completedataset$RegPredRed)
#7215 #pretty close agreement between the complete and reduced model 


save(completedataset, file=paste("PFAS_HLpreds_HLH_HLHLSsc_DSSToxchemicals_coarse_imputation_",writesuff,".RData", sep="")) 
load(paste("PFAS_HLpreds_HLH_HLHLSsc_DSSToxchemicals_coarse_imputation_",writesuff,".RData", sep=""))


####Applicability Domain####
#To be consistent with Mansouri et al, if something is in the domain, Domain=1, if outside=0. 
If if its in the domain=1, if outside the doma
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
RMFAD=ADfunction(model=regmod,data=completedataset)
RMRAD=ADfunction(model=regmodred,data=completedataset)

names(CMFAD)[which(names(CMFAD)=="Domain")]="ClassModDomain"
names(RMFAD)[which(names(RMFAD)=="Domain")]="RegFullModDomain"
names(RMRAD)[which(names(RMRAD)=="Domain")]="RegRedModDomain"

names(CMFAD)[which(names(CMFAD)=="SI90")]="ClassModSI90"
names(RMFAD)[which(names(RMFAD)=="SI90")]="RegFullModSI90"
names(RMRAD)[which(names(RMRAD)=="SI90")]="RegRedModSI90"

DPFASwAD=data.frame(CMFAD[,-c(which(names(CMFAD)%in%c("SImin",  "SImax", "SImean","SIsd")))], RMFAD[,c("RegFullModSI90", "RegFullModDomain")], RMRAD[,c("RegRedModSI90", "RegRedModDomain")]) 

#Pull in life span information prior to saving
library(openxlsx)
LF=read.xlsx("LifeSpan_from_Risa_05_25_20.xlsx")

unique(completedataset$Species)

PerLS=ifelse(DPFASwAD$Species=="Human", LF[which(LF$Species=="Human"), "HrLifeSpan"],  
             ifelse(DPFASwAD$Species=="Monkey", LF[which(LF$Species=="Monkey"), "HrLifeSpan"],
                    ifelse(DPFASwAD$Species=="Rat", LF[which(LF$Species=="Rat"), "HrLifeSpan"],
                           ifelse(DPFASwAD$Species=="Mouse", LF[which(LF$Species=="Mouse"), "HrLifeSpan"],
                                  ifelse(DPFASwAD$Species=="Dog", LF[which(LF$Species=="Dog"), "HrLifeSpan"], NA)))))


DPFASwAD$PerLS=PerLS

save(DPFASwAD, file=paste0("Tox21_AllMods_ADindicate_", writesuff,".RData"))
write.csv(DPFASwAD, file=paste0("Tox21_AllMods_ADindicate_", writesuff,".csv"))
load(file=paste0("Tox21_AllMods_ADindicate_", writesuff,".RData"))

DPFASwAD[DPFASwAD$DTXSID=="DTXSID8047553" & DPFASwAD$Species=="Human",]

#Note, only 1917 chemicals are listed as being in the domain of the reduced regresson model 
length(unique(DPFASwAD[DPFASwAD$ClassModDomain==1, "CASRN"])) #3128
length(unique(DPFASwAD[DPFASwAD$RegFullModDomain==1, "CASRN"])) #3128
length(unique(DPFASwAD[DPFASwAD$RegRedModDomain==1, "CASRN"])) #3128

###Domain tables
dom=DPFASwAD
dom=subset(dom, Species%in%c("Human","Monkey", "Rat", "Mouse", "Dog") & Sex=="Female" & DosingAdj=="Other")
CDtab=aggregate(ClassModDomain~Species,data=dom, FUN = "sum")
RRtab=aggregate(RegRedModDomain~Species,data=dom, FUN = "sum")
RFtab=aggregate(RegFullModDomain~Species,data=dom, FUN = "sum")
dom2=subset(dom, AD_LogP==1 & AD_VP==1 & AD_WS==1 & AD_LogD==1)
CDtab2=aggregate(ClassModDomain~Species,data=dom2, FUN = "sum")
RRtab2=aggregate(RegRedModDomain~Species,data=dom2, FUN = "sum")


domtab1=data.frame("Model"="Classification",CDtab,"All_Domains"=CDtab2[,"ClassModDomain"])
domtab2=data.frame("Model"="Regression_Reduced", RRtab, "All_Domains"=RRtab2[,"RegRedModDomain"])
names(domtab1)[3]="ModelDomain"
names(domtab2)[3]="ModelDomain"
domtab=rbind(domtab1, domtab2)
domtab
####


Modellist=list("ClassModFull"=classmod4,  "RegModFull"=regmod, "RegModRed"=regmodred)
save(Modellist, file=paste0("Class_Reg_Full_Reduced_MOdels_",writesuff,".RData"))

#This will serve as an interesting discussion of in versus out of the AD, and how it differs by species and by model. 
#I'll have 


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

#########Compare with training set chemicals
load(paste0("PFAS_QSAR_Predictors_Predictions_", writesuff,".RData"))
ps=DPFASwAD[which(DPFASwAD$DTXSID%in%pfasdsredsc_ds$DTXSID),]
ps=ps[ps$Species=="Human" & ps$Sex=="Male" & ps$DosingAdj=="Other",]
dshuman=pfasdsredsc_ds
dshuman=dshuman[c(dshuman$Species=="Human" & dshuman$Sex=="Male"),]

merge(ps[c("DTXSID", "ClassPredFull")], dshuman[c("DTXSID", "ClassPredFull")], by="DTXSID")
merge(ps[c("DTXSID", "RegPredFull")], dshuman[c("DTXSID", "RegPredFull")], by="DTXSID")
merge(ps[c("DTXSID", "RegPredRed")], dshuman[c("DTXSID", "RegPredRed")], by="DTXSID")

#Note: while the classification model predictiosn for the training chemicals is the same for the 
#the training model and when applied here, three of the chemicals have different regression model results. 
#the predictors are the same, so this must result from model instability issues. 
chems=c("DTXSID1037303", "DTXSID4059916")

sub1=subset(dshuman, c(DTXSID%in%chems))
sub2=subset(ps, c(DTXSID%in%chems))
merge(sub1, sub2,  by="DTXSID")

#####


#Limit the complete dataset to only the species in the list 
ds=subset(DPFASwAD, Species%in%c("Human", "Monkey", "Mouse", "Rat", "Dog"), Sex="Male")
ds$RegpredFull_byLS=ds$RegPredFull/ds$PerLS
ds$RegpredRed_byLS=ds$RegPredRed/ds$PerLS
ds$SpeciesSex=paste0(ds$Species, "_", ds$Sex)
ds$SpeciesDose=paste0(ds$Species, "_", ds$DosingAdj)
pfascount=aggregate(Count~ ClassPredFull * SpeciesDose  , data=ds, FUN="sum")
pfascount=aggregate(Count~ ClassPredFull * SpeciesSex  , data=ds, FUN="sum")
#Identify specific combinations that we want to try
#ds$SpeciesSexDose=paste0(ds$Species, "_", ds$Sex, "_", ds$DosingAdj)
#pfascountCFSD=aggregate(Count~ ClassPredFull * SpeciesDose * Sex , data=ds, FUN="sum")
#pfascountCFSS=aggregate(Count~ ClassPredFull * Species*Sex, data=ds, FUN="sum")
#pfascountCFSS=aggregate(Count~ ClassPredFull * Species, data=ds, FUN="sum")

#Note, I will need get the bin medians for each bin 

##CLassification Model 
####Numbers of chemicals by bin
namelabels=c("Dog:Female", "Dog:Male", "Human:Female", "Human:Male", "Monkey:Female", "Monkey:Male", "Mouse:Female", "Mouse:Male", "Rat:Female", "Rat:Male")
PredictionPlotCF_Raw=ggplot(data=pfascount, aes(x=SpeciesDose, y=Count, fill=ClassPredFull)) +
  geom_bar(stat="identity")+
  xlab("Species/Sex")+
  ylab("Number of PFAS chemicals")+
  
  ggtitle(paste("RF Classification Model: Serum Half-Life of ",length(unique(completedataset$CASRN)), " PFAS Chemicals",sep=""))+
  scale_fill_discrete(name="Serum Half Life", labels = c("<1 Day", "<Week", "< 2 Months", "> 2 Months"))+   
  theme(
    plot.title=element_text(size =15), 
    axis.text.x =element_text(angle=90, size=15, vjust=0.3, hjust=1),
    axis.text.y=element_text(size=15),
    axis.title=element_text(size=17),
    legend.position="bottom",
    legend.text=element_text(size=15),
    legend.title=element_text(size=20))#+
#scale_x_discrete(labels= namelabels)

PredictionPlotCF_Raw

png(paste("Figures/PFAS_DSSTox_",length(unique(completedataset$CASRN)),"chemicals_HLHBinPred_coarse_imputation_ClassPred_",writesuff,".png",sep=""), width=800, height = 668)
PredictionPlotCF_Raw
dev.off()

png(paste("Figures/PFAS_DSSTox_",length(unique(completedataset$CASRN)),"chemicals_HLHBinPred_coarse_imputation_ClassPred_ByDose_",writesuff,".png",sep=""), width=800, height = 668)
PredictionPlotCF_Raw
dev.off()



#Load median bin values and scale chemicals by life span of species 
load(paste0("Median_HLH_per_Bin_4_", writesuff,".RData"))
ds$ClassHLHPred=ifelse(ds$ClassPredFull==1, agtabHLH[1,2], ifelse(ds$ClassPredFull==2, agtabHLH[2,2],ifelse(ds$ClassPredFull==3, agtabHLH[3,2],agtabHLH[4,2])))
ds$ClassHLHPred_PerLS=ds$ClassHLHPred/ds$PerLS
human=ds[c(ds$Species=="Human" & ds$Sex=="Male"),]
human
quantile(human$ClassHLHPred_PerLS )
table(human$ClassHLHPred_PerLS )

rat=ds[c(ds$Species=="Rat"),]
table(rat$ClassHLHPred_PerLS )

dog=ds[c(ds$Species=="Dog"),]
table(dog$ClassHLHPred_PerLS )

PredictionPlotCF_LFScaled=ggplot(data = ds, aes(x = Species,  y = ClassHLHPred_PerLS)) +
  labs(x = "Species/Sex", 
       y = "Serum Half-Lives as a Proportion of Average Lifespan")+ 
  ggtitle("RF Classification Model: Distribution of Predicted Serum Half-Lives\n as Proportion of Average Lifespan", 
          subtitle = paste0(length(unique(ds$CASRN)), " PFAS Chemicals of DSSTox List" )) +
  theme(plot.title = element_text(size = 15, face = "bold"), 
        axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.text=element_text(size=15))+
  #scale_y_log10(limits=c(1e-6, 1)) +
  ylim(0,1)+
  geom_boxplot()

PredictionPlotCF_LFScaled
png(paste("Figures/PFAS_DSSTox_",length(unique(completedataset$CASRN)),"chemicals_ClassPred_byLifeSpan",writesuff,".png",sep=""), width=800, height = 668)
PredictionPlotCF_LFScaled
dev.off()


######Regression Models 
###Prediction by Full Regression MOdel 
PredictionPlotRF=ggplot(data = ds, aes(x = Species,  y = RegPredFull)) +
  labs(x = "Species", 
       y = "Serum Half-Lives in Hours")+ 
  ggtitle("RF Full Regression Model: Distribution of Predicted Serum Half-Lives (Hrs)", 
          subtitle = paste0(length(unique(ds$CASRN)), " PFAS Chemicals of DSSTox List" )) +
  theme(plot.title = element_text(size = 15, face = "bold"), 
        axis.title.y = element_text(size = 12),
        axis.title.x = element_text(size = 12))+
  #scale_y_log10(label=comma)+
  scale_y_continuous(labels = comma)  +
  #  scale_y_log10() +
  geom_boxplot()

PredictionPlotRF
png(paste("Figures/PFAS_DSSTox_",length(unique(completedataset$CASRN)),"chemicals_RegPredFull_",writesuff,".png",sep=""), width=800, height = 668)
PredictionPlotRF
dev.off()

###Scaled by halflife 
PredictionPlotRF_LS=ggplot(data = ds, aes(x = Species,  y = RegpredFull_byLS)) +
  labs(x = "Species", 
       y = "Serum Half-Lives as a Proportion of Average Lifespan")+ 
  ggtitle("RF Full Regression Model: Distribution of Predicted Serum Half-Lives (Hrs)", 
          subtitle = paste0(length(unique(ds$CASRN)), " PFAS Chemicals of DSSTox List" )) +
  theme(plot.title = element_text(size = 15, face = "bold"), 
        axis.title.y = element_text(size = 12),
        axis.title.x = element_text(size = 12))+
  
  scale_y_log10() +
  geom_boxplot()

PredictionPlotRF_LS
png(paste("Figures/PFAS_DSSTox_",length(unique(completedataset$CASRN)),"chemicals_RegPredFull_byLifeSpan_",writesuff,".png",sep=""), width=800, height = 668)
PredictionPlotRF_LS
dev.off()



###Prediction by Reduced Regression MOdel 
library(scales)
namelabels1=c("Dog:IV", "Dog:Oral","Dog:Other","Human:IV", "Human:Oral","Human:Other","Monkey:IV", "Monkey:Oral","Monkey:Other","Mouse:IV", "Mouse:Oral","Mouse:Other","Rat:IV", "Rat:Oral","Rat:Other")
PredictionPlotRR_Raw=ggplot(data = ds, aes(x = SpeciesDose,  y = RegPredRed)) +
  labs(x = "Species", 
       y = "Serum Half-Lives in Hours")+ 
  ggtitle("RF Reduced Regression Model: Distribution of Predicted Serum Half-Lives (Hrs)", 
          subtitle = paste0(length(unique(ds$CASRN)), " PFAS Chemicals of DSSTox List" )) +
  theme(
    plot.title=element_text(size =15), 
    axis.text.x =element_text(angle=90, size=15, vjust=0.3, hjust=1),
    axis.text.y=element_text(size=15),
    axis.title=element_text(size=17),
    legend.position="bottom",
    legend.text=element_text(size=15),
    legend.title=element_text(size=20))+
  scale_x_discrete(labels= namelabels1)+
  scale_y_continuous(labels = comma) +
  geom_boxplot()

PredictionPlotRR_Raw
png(paste("Figures/PFAS_DSSTox_",length(unique(completedataset$CASRN)),"chemicals_RegPredRed",writesuff,".png",sep=""), width=800, height = 668)
PredictionPlotRR_Raw
dev.off()

###Scaled by halflife 
PredictionPlotRR_LS=ggplot(data = ds, aes(x = Species,  y = RegpredRed_byLS)) +
  labs(x = "Species", 
       y = "Serum Half-Lives as a Proportion of Average Lifespan")+ 
  ggtitle("RF Reduced Regression Model: Distribution of Predicted Serum Half-Lives (Hrs)", 
          subtitle = paste0(length(unique(ds$CASRN)), " PFAS Chemicals of DSSTox List" )) +
  theme(plot.title = element_text(size = 15, face = "bold"), 
        axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.text = element_text(size=15))+
  scale_y_log10() +
  geom_boxplot()

PredictionPlotRR_LS
png(paste("Figures/PFAS_DSSTox_",length(unique(completedataset$CASRN)),"chemicals_RegPredRed_byLifeSpan_",writesuff,".png",sep=""), width=800, height = 668)
PredictionPlotRR_LS
dev.off()

#Quantile of mouse 1/2 as function of reproductive half-life. 
quantile(ds[ds$Species=="Mouse", "RegpredRed_byLS"])


###Reduce Models by inclusion in all domains###############################################
#Only include chemicals in domain of all opera model
sub1=subset(ds, AD_LogP==1 & AD_VP==1 & AD_WS==1 & AD_LogD==1)
length(unique(sub1$DTXSID)) #420
#Only Include chemicals in domain of Classification model
sub2=subset(sub1, ClassModDomain==1) 
length(unique(sub2$DTXSID)) #252 chemicals 

#Only Include chemicals in domain of reduced Regression model
sub3=subset(sub1, RegRedModDomain==1) 
length(unique(sub3$DTXSID)) #252 chemicals 


##Classification Model 
####Numbers of chemicals by bin
pfascount=aggregate(Count~ ClassPredFull * Species  , data=sub2, FUN="sum")

PredictionPlotCF_Raw_InDomain=ggplot(data=pfascount, aes(x=Species, y=Count, fill=ClassPredFull)) +
  geom_bar(stat="identity")+
  xlab("Species")+
  ylab("Number of PFAS chemicals")+
  ggtitle(paste("RF Classification Model: Serum Half-Life of ",length(unique(completedataset$CASRN)), " PFAS Chemicals",sep=""),
          subtitle = paste0(length(unique(sub2$CASRN)), " PFAS Chemicals in Domain" ))+
  scale_fill_discrete(name="Serum Half Life", labels = c("<1 Day", "<Week", "< 2 Months", "> 2 Months"))+   
  theme(
    plot.title=element_text(size =15), 
    axis.text.x =element_text(angle=90, size=15),
    axis.text.y=element_text(size=15),
    axis.title=element_text(size=17),
    legend.position="bottom",
    legend.text=element_text(size=15),
    legend.title=element_text(size=20))
PredictionPlotCF_Raw_InDomain

png(paste("Figures/PFAS_DSSTox_",length(unique(completedataset$CASRN)),"chemicals_HLHBinPred_coarse_imputation_ClassPred_InDomain_",writesuff,".png",sep=""), width=800, height = 668)
PredictionPlotCF_Raw_InDomain
dev.off()

#Load median bin values and scale chemicals by life span of species 
PredictionPlotCF_LFScaled_InDomain=ggplot(data = sub2, aes(x = Species,  y = ClassHLHPred_PerLS)) +
  labs(x = "Species", 
       y = "Serum Half-Lives as a Proportion of Average Lifespan")+ 
  ggtitle("RF Classification Model: Distribution of Predicted Serum Half-Lives\n as Proportion of Average Lifespan",
          subtitle = paste0("Of ", length(unique(ds$CASRN)), " PFAS Chemicals of DSSTox List, ", length(unique(sub2$CASRN)), " in Model Domain"))+ 
  theme(plot.title = element_text(size = 15, face = "bold"), 
        axis.title.y = element_text(size = 12),
        axis.title.x = element_text(size = 12))+
  #  scale_y_log10(limits=c(1e-6, 1)) +
  geom_boxplot()

PredictionPlotCF_LFScaled_InDomain
png(paste("Figures/PFAS_DSSTox_",length(unique(completedataset$CASRN)),"chemicals_ClassPred_byLifeSpan_InDomain_",writesuff,".png",sep=""), width=800, height = 668)
PredictionPlotCF_LFScaled_InDomain
dev.off()


######Regression Models 
###Prediction by Full Regression MOdel 
PredictionPlotRF_Raw_InDomain=ggplot(data = sub1, aes(x = Species,  y = RegPredFull)) +
  labs(x = "Species", 
       y = "Serum Half-Lives in Hours")+ 
  ggtitle("RF Classification Model: Distribution of Predicted Serum Half-Lives\n as Proportion of Average Lifespan",
          subtitle = paste0("Of ", length(unique(ds$CASRN)), " PFAS Chemicals of DSSTox List, ", length(unique(sub2$CASRN)), " in Model Domain"))+ 
  theme(plot.title = element_text(size = 15, face = "bold"), 
        axis.title.y = element_text(size = 12),
        axis.title.x = element_text(size = 12))+
  scale_y_log10() +
  geom_boxplot()

PredictionPlotRF_Raw_InDomain
png(paste("Figures/PFAS_DSSTox_",length(unique(completedataset$CASRN)),"chemicals_RegPredFull_InDomain_",writesuff,".png",sep=""), width=800, height = 668)
PredictionPlotRF_Raw_InDomain
dev.off()

###Scaled by halflife 
PredictionPlotRF_LS_InDomain=ggplot(data = sub1, aes(x = Species,  y = RegpredFull_byLS)) +
  labs(x = "Species", 
       y = "Serum Half-Lives as a Proportion of Average Lifespan")+ 
  ggtitle("RF Full Regression Model: Distribution of Predicted Serum Half-Lives (Hrs)", 
          subtitle = paste0("Of ", length(unique(ds$CASRN)), " PFAS Chemicals of DSSTox List, ", length(unique(sub2$CASRN)), " in Model Domain"))+ 
  theme(plot.title = element_text(size = 15, face = "bold"), 
        axis.title.y = element_text(size = 12),
        axis.title.x = element_text(size = 12))+
  scale_y_log10() +
  geom_boxplot()

PredictionPlotRF_LS_InDomain
png(paste("Figures/PFAS_DSSTox_",length(unique(completedataset$CASRN)),"chemicals_RegPredFull_byLifeSpan_InDomain",writesuff,".png",sep=""), width=800, height = 668)
PredictionPlotRF_LS_InDomain
dev.off()



###Prediction by Reduced Regression MOdel 
PredictionPlotRR_Raw_InDomain=ggplot(data = sub3, aes(x = Species,  y = RegPredRed)) +
  labs(x = "Species", 
       y = "Serum Half-Lives in Hours")+ 
  ggtitle("RF Reduced Regression Model: Distribution of Predicted Serum Half-Lives (Hrs)", 
          subtitle = paste0("Of ", length(unique(ds$CASRN)), " PFAS Chemicals of DSSTox List, ", length(unique(sub3$CASRN)), " in Model Domain"))+ 
  theme(plot.title = element_text(size = 15, face = "bold"), 
        axis.title.y = element_text(size = 12),
        axis.title.x = element_text(size = 12))+
  scale_y_log10() +
  geom_boxplot()

PredictionPlotRR_Raw_InDomain
png(paste("Figures/PFAS_DSSTox_",length(unique(completedataset$CASRN)),"chemicals_RegPredRed_InDomain_",writesuff,".png",sep=""), width=800, height = 668)
PredictionPlotRR_Raw_InDomain
dev.off()

###Scaled by halflife 
PredictionPlotRR_LS_InDomain=ggplot(data = sub3, aes(x = Species,  y = RegpredRed_byLS)) +
  labs(x = "Species", 
       y = "Serum Half-Lives as a Proportion of Average Lifespan")+ 
  ggtitle("RF Reduced Regression Model: Distribution of Predicted Serum Half-Lives (Hrs)", 
          subtitle = paste0("Of ", length(unique(ds$CASRN)), " PFAS Chemicals of DSSTox List, ", length(unique(sub3$CASRN)), " in Model Domain"))+ 
  theme(plot.title = element_text(size = 15, face = "bold"), 
        axis.title.y = element_text(size = 12),
        axis.title.x = element_text(size = 12))+
  scale_y_log10() +
  geom_boxplot()

PredictionPlotRR_LS_InDomain
png(paste("Figures/PFAS_DSSTox_",length(unique(completedataset$CASRN)),"chemicals_RegPredRed_byLifeSpan_InDomain_",writesuff,".png",sep=""), width=800, height = 668)
PredictionPlotRR_LS_InDomain
dev.off()


range(ds[ds$Species=="Human", "RegpredRed_byLS"]) 
(ds[ds$Species=="Mouse", "RegpredRed_byLS"]) 

#perfluoroundecnoic acid for Jeff Municci
DPFASwAD[DPFASwAD$DTXSID=="DTXSID8047553" & DPFASwAD$Species=="Human", c("CASRN", "DTXSID", "ClassPredFull", "RegPredFull")]

