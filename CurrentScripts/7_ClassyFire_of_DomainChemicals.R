# Clear the workspace:
rm(list=ls())
try(dev.off())

#remotes::install_github('aberHRML/classyfireR')
library(classyfireR)
library(readxl)
library(data.table) 
library(purrr)
library(tidyr)
suff <- "JFW100322-noLogD" #Note, this is the suffix for previous work on this page
seed <- "12345"
trainingsetsuff="JFW100322"
writesuff<-"JFW100322-noLogD"


#IN this script, we are going to read in all of the chemicals from the DSSTox dataset to see how
#chemicals in the model domain versus chemicals not in the domain differ. 

# Read in data
#setwd('L:/Lab/NCCT_ExpoCast/ExpoCast2021/Dawson_PFAS_HALFLIFE/PFAS_HL_QSAR_2021')
#PFASChems <- read_excel("All_Chemicals_CatModDomain_CompToxDashboard_",writesuff,".xlsx", sheet=2)
PFASChems <- read.csv(paste("RData/AllChemicals_in_CatMod_Domain_",writesuff,".csv",sep=""))

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
# Select relevant columns

#Writeout PFAS data to get inchikeys on comptox dashboard
#write.csv(PFASdata, file=paste0("All_PFAS_DSSTox_", writesuff, ".csv"))
PFASdata=read.csv("PFAS_Catalog/DSSTox_AllPFAS_with_INCHI_030722.csv")
length(which(is.na(PFASdata$INCHIKEY))) #All have INCHIKeys

#Get classifications for all of the chemicals on the PFAS data list
inchiloc=which(names(PFASdata)=="INCHIKEY")
classifylist=sapply(PFASdata[,inchiloc], get_classification)
saveRDS(classifylist, file=paste0("RData/ClassyFire_output_PFAS_Chems_DSSTox_HLHClassMod_Domain_",writesuff,".rds"))
classifylist=readRDS(paste0("RData/ClassyFire_output_PFAS_Chems_DSSTox_HLHClassMod_Domain_",writesuff,".rds"))
classifylist[[1]]

#Parse out the columns
# Add in columns to the data.table to reflect the classification of each chemical
which(is.null(classifylist))
which(is.na(classifylist))

PFASdata1=as.data.table(PFASdata)

for (i in c(1:length(classifylist))){
  #print(i)
  if(is.null(classifylist[[i]])==TRUE){next}
  temp <- classification((classifylist[[i]]))$Classification
  classifiers <- c(temp, rep('', 11 - length(temp)))
  rowloc=which(paste0("InChIKey=",PFASdata$INCHIKEY)==classyfireR::meta(classifylist[[i]])$inchikey)
  #print(classifiers)
  PFASdata1[rowloc, c("kingdom", "superclass", "class", "subclass", "level5", "level6", "level7", "level8", "level9", "level10", "level11") := as.list(classifiers)]
}


#Want to do three main questions: analyses:
#1. How does the DSSTox set out of my domain differ from the one that's in it?
#2. How does the DSSTox set differ from the training data?
#3. How does the in-domain chemicals compare to the training data? 



#Read in list of domains and extract out only chemical for humans in the AMAD
load(file=paste0("RData/Tox21_AllMods_ADindicate_", writesuff,".RData"))
ds=subset(DPFASwAD, Species%in%c("Human", "Monkey", "Mouse", "Rat", "Dog"))
ds$SpeciesSex=paste0(ds$Species, "_", ds$Sex)
ds$SpeciesDose=paste0(ds$Species, "_", ds$DosingAdj)
ds$SpeciesDose=paste0(ds$Species, "_", ds$DosingAdj)
#ds=ds[ds$ClassModDomain==1,]
#dsAMAD=ds[ds$AD_LogP==1 & ds$AD_LogP==1 & ds$AD_VP==1 & ds$Species=="Human",]
#PFASdata=PFASdata[which(PFASdata$DTXSID%in%ds$DTXSID==TRUE),]
#PFASdata=as.data.table(PFASdata)

TotalDataset=merge(PFASdata1,ds, by="DTXSID")
saveRDS(TotalDataset, paste0("RData/ClassyFire_PFAS_DSSTox_HLHClassMod_AllData_", writesuff,".rds")) 
TD=readRDS(paste0("RData/ClassyFire_PFAS_DSSTox_HLHClassMod_AllData_", writesuff,".rds"))
#TD=TotalDataset
TD=as.data.frame(TD)

#Refining down to only Human, Female, Other dosing,AMAD; previous check found that IV/Oral dosing and male have the same chemicals

load(paste0("RData/PFAS_11Chemicals_QSARdataset_",trainingsetsuff,".RData", sep="")) #loads pfasds
Trainchem=unique(pfasds$DTXSID)
Trainchemsub=TD[which(TD$DTXSID%in%Trainchem),]
Trainchemsub=Trainchemsub[Trainchemsub$Species=="Human" & Trainchemsub$DosingAdj=="Other" & Trainchemsub$Sex=="Female",]
Trainchemsub
table(Trainchemsub$kingdom)
table(Trainchemsub$superclass)
table(Trainchemsub$class)
table(Trainchemsub$subclass)
table(Trainchemsub$level5) #; this level shows that most of the training chems are per either perflouro carboxylic acid directives or pf sulfonic acid derivatives, with slight differences for GexN and the other one. 
table(Trainchemsub$level6)
Trainchemsub$class


#Level 5 or subclass

##Restrict TD and to AMAD with only one case to eliminate duplicates
TD=as.data.frame(TD)
length(unique(TD$DTXSID))
TD1=TD[TD$Species=="Human" & TD$Sex=="Female" & TD$DosingAdj=="Other",]
df=data.frame(table(TD1$class))
df=df[order(df$Freq, decreasing = TRUE),]
df$Per=df$Freq/sum(df$Freq) #8% alkyl fluorides; 10% are phospate esters,32% are not classified at that level
df$CumPer=cumsum(df$Per)  #50% are made up of unclassified, alkyl flouorides, or phosphate esters
dim(df)

dim(table(TD1$class))
table(TD1$class)/sum(table(TD1$class))

OutDomain=TD1[TD1$ClassModDomain==0,] 
dim(OutDomain) #591 chemicals
df1=data.frame(table(OutDomain$class))
df1=df1[order(df1$Freq, decreasing = TRUE),]
df1$Per=df1$Freq/sum(df1$Freq) #8% alkyl fluorides; 10% are phospate esters,32% are not classified at that level
df1$CumPer=cumsum(df1$Per)  #50% are made up of unclassified, alkyl flouorides, or phosphate esters
dim(df1)


InDomain=TD1[TD1$ClassModDomain==1,]
df2=data.frame(table(InDomain$class))
df2=df2[order(df2$Freq, decreasing = TRUE),]
df2$Per=df2$Freq/sum(df2$Freq) #11% are alkyl fluorides, 24% are not classified, phosphate esters are only 1.3%; fatty acids 
df2$CumPer=cumsum(df2$Per)#; unclassified, alkyl flourides, fatty acids and conjugates
dim(df2)

jaccard <- function(a, b) {
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection
  return (intersection/union)
}

jaccard(df2$Var1, df1$Var1)



AMAD=TD[TD$AD_LogP==1 & TD$AD_VP==1 & TD$AD_WS==1 & TD$ClassModDomain==1 & TD$Species=="Human" & TD$DosingAdj=="Other" & TD$Sex=="Female",]
dim(table(AMAD$class))
df3=data.frame(table(AMAD$class))
df3=df3[order(df3$Freq, decreasing = TRUE),]
df3$Per=df3$Freq/sum(df3$Freq) #11% are alkyl fluorides, 24% are not classified, phosphate esters are only 1.3%; fatty acids 
df3$CumPer=cumsum(df3$Per)#; unclassified, alkyl flourides, fatty acids and conjugates
dim(df3)

table(AMAD$subclass)
table(Trainchemsub$subclass)
table(Trainchemsub$level5)





#Proportion of chemicals from training class making up domains
length(which(AMAD$class=="Alkyl halides" | AMAD$class== "Carboxylic acids and derivatives" | AMAD$class=="Organic sulfonic acids and derivatives"))/dim(AMAD)[1]  


am=data.frame(table(AMAD$subclass)/sum(table(AMAD$subclass)))
dim(table(Trainchemsub$class))




#Now, can parse by the two domains(HLHDomain, AMAD):
AMAD=as.data.frame(TotalDataset)
OUT=AMAD[c(AMAD$AD_LogP==0 | AMAD$AD_VP==0 | AMAD$AD_WS==0) & AMAD$ClassModDomain==1 & AMAD$Species=="Human" & AMAD$Sex=="Female" & AMAD$DosingAdj=="Other",]
dim(OUT)

LP=OUT[OUT$AD_LogP==0 & OUT$AD_VP==1 & OUT$AD_WS==1 & OUT$ClassModDomain==1 & OUT$Species=="Human" & OUT$Sex=="Female" & OUT$DosingAdj=="Other",]
VP=OUT[OUT$AD_LogP==1 & OUT$AD_VP==0 & OUT$AD_WS==1 & OUT$ClassModDomain==1 & OUT$Species=="Human" & OUT$Sex=="Female" & OUT$DosingAdj=="Other",]
WS=OUT[OUT$AD_LogP==1 & OUT$AD_VP==1 & OUT$AD_WS==0 & OUT$ClassModDomain==1 & OUT$Species=="Human" & OUT$Sex=="Female" & OUT$DosingAdj=="Other",]
WSLP=OUT[OUT$AD_LogP==0 & OUT$AD_VP==1 & OUT$AD_WS==0 & OUT$ClassModDomain==1 & OUT$Species=="Human" & OUT$Sex=="Female" & OUT$DosingAdj=="Other",]
WSVP=OUT[OUT$AD_LogP==1 & OUT$AD_VP==0 & OUT$AD_WS==0 & OUT$ClassModDomain==1 & OUT$Species=="Human" & OUT$Sex=="Female" & OUT$DosingAdj=="Other",]
LPVP=OUT[OUT$AD_LogP==0 & OUT$AD_VP==0 & OUT$AD_WS==1 & OUT$ClassModDomain==1 & OUT$Species=="Human" & OUT$Sex=="Female" & OUT$DosingAdj=="Other",]
WSLPVP=OUT[OUT$AD_LogP==0 & OUT$AD_VP==0 & OUT$AD_WS==0 & OUT$ClassModDomain==1 & OUT$Species=="Human" & OUT$Sex=="Female" & OUT$DosingAdj=="Other",]

dim(LP)[1]/3248
dim(VP)[1]/3248
dim(WS)[1]/3248
dim(WSLP)[1]/3248
dim(WSVP)[1]/3248
dim(LPVP)[1]/3248
dim(WSLPVP)[1]/3248


#Parse which opera domains have trouble, and need some additional research
table(AMAD$AD_LogP)/sum(table(AMAD$AD_LogP))
table(AMAD$AD_VP)/sum(table(AMAD$AD_VP))
table(AMAD$AD_WS)/sum(table(AMAD$AD_WS))


###TK outputs for PFOA and PFOS
Trainchemsub[Trainchemsub$DTXSID=="DTXSID3031864" |  Trainchemsub$DTXSID=="DTXSID8031865", ]


##Now, look at predictions for only chemicals in teh AMAD that fall into the same 3 classes as the training set 
library(ggplot2)
AMAD=TD
AMAD=AMAD[c(AMAD$AD_LogP==1 | AMAD$AD_VP==1 | AMAD$AD_WS==1) & AMAD$ClassModDomain==1 ,]
TCSub=AMAD[AMAD$class=="Alkyl halides" | AMAD$class== "Carboxylic acids and derivatives" | AMAD$class=="Organic sulfonic acids and derivatives",]
NotTCSub=AMAD[AMAD$class!="Alkyl halides" | AMAD$class!= "Carboxylic acids and derivatives" | AMAD$class!="Organic sulfonic acids and derivatives",]

TCSub$Count=1


hist(as.numeric(NotTCSub$AVERAGE_MASS.x))
hist(as.numeric(TCSub$AVERAGE_MASS.x), add=TRUE, col="blue")

save(NotTCSub,TCSub,AMAD,file=paste("RData/ClassyfierSubSets_",writesuff,".RData",sep=""))

library(ggplot2)

load(paste("RData/ClassyfierSubSets_",writesuff,".RData",sep=""))

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



 
pfascountCFSS=aggregate(Count~ ClassPredFull +DosingAdj+SpeciesSex , 
  data=subset(TCSub,Species%in%c("Mouse","Rat","Monkey","Human")), 
  FUN="sum")

Plot_Species_by_Sex=ggplot(data=pfascountCFSS[pfascountCFSS$DosingAdj=="Other",], aes(x=SpeciesSex, y=Count, fill=ClassPredFull)) +
  geom_bar(stat="identity")+
  xlab("Species/Sex")+
  ylab("Number of PFAS chemicals")+
  ggtitle(paste("B) ",length(unique(TCSub$DTXSID)), " PFAS Chemicals Matching Training Set Class",sep=""))+
  scale_fill_discrete(name="Serum Half Life", labels = c("<0.5 Day", "<Week", "< 2 Months", "> 2 Months"))+   
  theme(
    plot.title=element_text(size =20), 
    # plot.subtitle = element_text(size=17),
    axis.text.x =element_text(size=15, vjust=0.3, hjust=0.5),
    axis.text.y=element_text(size=15),
    axis.title=element_text(size=20),
    legend.position="none",
    legend.text=element_text(size=15),
    legend.title=element_text(size=20))
Plot_Species_by_Sex=
  Plot_Species_by_Sex+  
  scale_x_discrete(labels= c("Dog_Female" = namelabelsSex[1], "Dog_Male"=namelabelsSex[2], "Human_Female"=namelabelsSex[3],
                             "Human_Male"=namelabelsSex[4],"Monkey_Female"=namelabelsSex[5], "Monkey_Male"=namelabelsSex[6], 
                             "Mouse_Female"=namelabelsSex[7], "Mouse_Male"=namelabelsSex[8], "Rat_Female"=namelabelsSex[9], 
                             "Rat_Male"=namelabelsSex[10]))


Plot_Species_by_Sex
png(paste("Figures/THalfPred_Species_by_Sex_B_",writesuff,".png",sep=""), width=800, height=550)
Plot_Species_by_Sex
dev.off()

#Human stats
HAMAD=TCSub[TCSub$Species=="Human" & TCSub$DosingAdj=="Other" & TCSub$Sex=="Female",]
table(HAMAD$ClassPredFull)/sum(table(HAMAD$ClassPredFull))


#Rats stats
RAMADF=TCSub[TCSub$Species=="Rat" & TCSub$DosingAdj=="Other" & TCSub$Sex=="Female",]
RAMADM=TCSub[TCSub$Species=="Rat" & TCSub$DosingAdj=="Other" & TCSub$Sex=="Male",]

table(RAMADF$ClassPredFull)/sum(table(RAMADF$ClassPredFull))
table(RAMADM$ClassPredFull)/sum(table(RAMADF$ClassPredFull))



##Export table for SI
#Join DPFASwAD with classification table 
ExportTab=merge(DPFASwAD, PFASdata1[,c("DTXSID", "PREFERRED_NAME", "INCHIKEY","kingdom", "superclass", "class","subclass")], by="DTXSID", all.x=TRUE )
write.csv(ExportTab, file=paste0("RData/MasterExport_ModPred_AD_TK_CF_", writesuff,".csv"))
TDADex=TDAD[,c(1,4:11,26:29,38:41,66,65,43,45:60,63,64)]

unique(DPFASwAD$Species)

head(DPFASwAD)
save(DPFASwAD,file="RDAta/")


