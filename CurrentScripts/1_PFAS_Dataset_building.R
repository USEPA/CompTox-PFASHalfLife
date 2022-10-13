# Clear the workspace:

#load necessary packages:
rm(list=ls())
packages=c("readxl","openxlsx", "httk")
sapply(packages, require,character.only=TRUE) #Note, the "character.only" argument is necessary her

# Identify which run we're creating:
suff="JFW100322"

#In this script, we are using PFAS 1/2 life data curated by Chris Lau recently(as of 2/2020) from various species and PFAS chemicals. 
#We are then combinining it with information from Oliver 1968, endogenous compound similarity from Prachi, Critical Micelle Concentration
#data from Bhhatarai, liver protein disassociation constants from Zhang, Serum Albumin binding constants from Han 


#In this script, we are using PFAS 1/2 life data curated by Chris Lau recently(as of 2/2020) from various species and PFAS chemicals. 
#We are then combinining it with information from Oliver 1968, endogenous compound similarity from Prachi, Critical Micelle Concentration
#data from Bhhatarai, liver protein disassociation constants from Zhang, Serum Albumin binding constants from Han 

#setwd("L:/Lab/NCCT_ExpoCast/ExpoCast2020/Dawson_TK_QSAR_PFAS/PFAS_HL_QSAR/ScriptsDataOutput")
#Merging Oliver 1968 information on Kidney physiological parameters


#New data
df=read.xlsx("PFAS_HalfLifeData/PFAS_HLH_Revised_DED_060121_matched_names.xlsx", sheet =1)

#Note: there are 89 lines


oliver=read.xlsx("Predictors/Oliver 1968 Table IX Single Mammal kidney_DED122019.xlsx", sheet =2)
oliver$Reference_KidneyPhys="Oliver1968"
oliver1=oliver[oliver$Mammal%in%unique(df$Species),] #only bring in species that we have PFAS data for. We'll add in  kidney info a bit later for mouse and monkey based on the whole list of species here.  
df1=merge(df, oliver1, by.x=c("Species"), by.y=c("Mammal"), all.x=TRUE)


#Note: Using similarity (Jaccard) of both PubChem fingerprints(which simply count whether a compound has certain features)

#Importing endogenous chemical similarity calculated by Prachi Pradeep and Richard Judson
Endoinfo=read.xlsx("Predictors/DawsonList-HMDBSimilarity.xlsx", sheet = 2)
names(Endoinfo)=c("CAS_pfas", "DTXSID_pfas", "CAS_endo", "DTXSID_endo","Tanimoto_PC", "Tanimoto_M")

#Find the EndoMax and Min, and and fill in tanimoto scores for the max score for each chemical, and then the list of chemicals
#in the max and min set 
maxlistPC=NULL
minlistPC=NULL
maxlistM=NULL
minlistM=NULL

EndoinfoPC=Endoinfo[-c(which(is.na(Endoinfo$Tanimoto_PC))),]
EndoinfoM=Endoinfo[-c(which(is.na(Endoinfo$Tanimoto_M))),]

for(i in unique(df1$DTXSID)){
  sub1=subset(EndoinfoPC, DTXSID_pfas == i)
  sub2=subset(EndoinfoM, DTXSID_pfas == i)
  maxlistPC=rbind(maxlistPC, sub1[sub1$Tanimoto_PC == max(sub1$Tanimoto_PC),])
  minlistPC=rbind(minlistPC, sub1[sub1$Tanimoto_PC == min(sub1$Tanimoto_PC),])
  maxlistM=rbind(maxlistM,sub2[sub2$Tanimoto_M == max(sub2$Tanimoto_M),])
  minlistM=rbind(minlistM, sub2[sub2$Tanimoto_M == min(sub2$Tanimoto_M),])
}

#Max chemicals 
maxlist1PC=aggregate(Tanimoto_PC~DTXSID_pfas, data=maxlistPC, FUN="mean") #This is in case multiple chemicals had the exact same similarity
maxlist1M=aggregate(Tanimoto_M~DTXSID_pfas, data=maxlistM, FUN="mean")

df2.1=merge(df1, maxlist1PC, by.x="DTXSID", by.y = "DTXSID_pfas", all = TRUE)
names(df2.1)[length(df2.1)]="MaxEndoPC"
df2.2=merge(df2.1, maxlist1M, by.x="DTXSID", by.y = "DTXSID_pfas", all = TRUE)
names(df2.2)[length(df2.2)]="MaxEndoM"

#Endoinfo=na.omit(Endoinfo)
#Max and Min Chemicals
df2.3=df2.2
maxminPC=unique(c(maxlistPC$CAS_endo, minlistPC$CAS_endo))
for(i in maxminPC){
  sub=Endoinfo[Endoinfo$CAS_endo==i,]  
  sub=subset(sub, select = c("CAS_pfas", "CAS_endo", "Tanimoto_PC"))
  df2.3=merge(df2.3, sub[, c("CAS_pfas", "Tanimoto_PC")], by.x="CASRN", by.y="CAS_pfas", all.x=TRUE)
  names(df2.3)[length(df2.3[1,])]=paste("TSPC_",sub$CAS_endo[1],sep="")  
}


maxminM=unique(c(maxlistM$CAS_endo, minlistM$CAS_endo))
for(i in maxminM){
  sub=Endoinfo[Endoinfo$CAS_endo==i,]  
  sub=subset(sub, select = c("CAS_pfas", "CAS_endo", "Tanimoto_M"))
  df2.3=merge(df2.3, sub[, c("CAS_pfas", "Tanimoto_M")], by.x="CASRN", by.y="CAS_pfas", all.x=TRUE)
  names(df2.3)[length(df2.3[1,])]=paste("TSM_",sub$CAS_endo[1],sep="")  
}

names(df2.3)

df2=df2.3
dim(df2) #should be 89 lines 

#Importing Bhhatrai 2010 Vapor Pressure(VP), Aquatic Solubility(AqS), and Critical Micelle Concentration(CMC) data
B=read.xlsx("Predictors/Bhhatarai_2010_PFAS_VP_AqS_CMC_data.xlsx")
#Replace "-" with NA's 
B=apply(B,2,function(x){replace(x, which(x=="-"), NA)})
B=apply(B,2,function(x){sub("\\*\\*", "", x)})
B=apply(B,2,function(x){sub("\\*", "", x)})

#Strip leading zeros off of CAS numbers
B[,1]=sub("^[0]+", "", B[,1])
AD=which(colnames(B)==c("AqS_AD") | colnames(B)==c("VP_AD") |colnames(B)==c("CMC_AD") )
numseq=seq(1,length(B[1,]),1)[-c(1,AD)]
B=as.data.frame(B)
for( i in numseq){
  B[,i] <- as.numeric(as.character(B[,i]))
}
#Cis and Trans 306-94-5 have identical values in the table; taking out cis
B=B[-c(which(B$CASRN=="306-94-5 cis")),]
B$Species = "Human"
B$Reference_CMC="Bhhatarai 2010"

#Subset by only CMC info
B=subset(B, select=c("CASRN", "Species", "CMC_Exp", "CMC_Pred",  "Reference_CMC"))
B$CASRN%in%df2$CASRN
df2$CASRN%in%B$CASRN #Note: Doesn't include 375-92-8, or Perfluoroheptanesulfonic acid
df3=merge(df2, B, by=c("CASRN", "Species"), all.x=TRUE)
dim(df3)

#Importing Protein Binding Affinity(Dissociation constants) from Zhang et al. 2013
Z=read.xlsx("Predictors/Zhang et al.2013_PFAS_Protein_Dissociation Constants.xlsx")
Z$Species="Human"
Z$Reference_Kd_hl_FABP="Zhang 2013"
#Z=Z[-c(which(names(Z)=="ChemicalName"))]
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
Z=Z[,-c(which(names(Z)=="DTXSID")),]
Z=Z[,-c(which(names(Z)=="PREFERRED_NAME"))]
#Changing non-detects to 0
#Z$Kd_hL_FABP=as.numeric(ifelse(Z$Kd_hL_FABP=="ND", 0, Z$Kd_hL_FABP))
df3[which(df3$CASRN%in%Z$CASRN==FALSE),1:4] #chemicals in dataset not in Zhang

df4=merge(df3, Z, by=c("CASRN", "Species"), all.x=TRUE)
cbind(df4$CASRN, df4$PREFERRED_NAME, df4$Kd_hL_FABP)
dim(df4)

#Importing Han et al. 2012
H=read.xlsx("Predictors/Han_et al.2012_AlbuminProteinBinding_Data.xlsx")
H$Reference_Ka_SerAlb="Han 2011"
H[which(H$PREFERRED_NAME=="Perfluorobutanesulfonate"), which(names(H)=="PREFERRED_NAME")]<-"Perfluorobutanesulfonic acid"


df5=merge(df4, H, by=c("PREFERRED_NAME", "Species"), all.x=TRUE)
dim(df5)
save(df5,file=paste0("RData/PredictorSet_pre-interpolation_", suff))

##Input missing data within chemical 

##Columns to try to fill in#
#Go chemical at a time; one data for one species exists fill in from the rest; either fill in or average existing values
uniquecasrn=unique(df5$CASRN)
checkcols=c("MaxEndoPC", "MaxEndoM", "CMC_Exp", "CMC_Pred", "Kd_hL_FABP", "Ka_perM_SerAlb_Han")
subset(df5, select=checkcols)
df6=NULL
for( i in 1:length(uniquecasrn)){
  sub=df5[df5$CASRN==uniquecasrn[i],]
  for(j in 1:length(checkcols)){
    sub1=sub[,which(names(sub)==checkcols[j])]
    if(length(which(is.na(sub1)==F))==1){
      sub[,which(names(sub)==checkcols[j])]=sub1[which(is.na(sub1)==F)]
    } else if (length(which(is.na(sub1)==F)) > 1 & length(which(is.na(sub1)==F)) < length(sub1))
    {sub[which(is.na(sub[,names(sub)==checkcols[j]])==T),which(names(sub)==checkcols[j])]=mean(sub1[which(is.na(sub1)==F)])} else {next}
  }
  df6=rbind(df6, sub)
}
dim(df6)
subset(df6,select=c("CASRN","PREFERRED_NAME",checkcols))

####Input missing data across chemicals;average missing values across species 
df7=df6

for(j in 1:length(checkcols)){
  sub=df6[,which(names(df6)==checkcols[j])]
  if (length(which(is.na(sub)==F)) > 1 & length(which(is.na(sub1)==F)) < length(sub))
  {
    df7[which(is.na(df7[,which(names(df7)==checkcols[j])])),which(names(df7)==checkcols[j])]=mean(sub[which(is.na(sub)==F)])
  } else { next  
  }
}
subset(df7,select=checkcols)
dim(df7)
#Filling in Kidney data: Use linear regressions using logKW or BW 
kidney=read.xlsx("Predictors/Oliver 1968_SupplementalData_JFW100322.xlsx", sheet =2)[1:11,]
kidney$Neph_Num <- as.numeric(kidney$Neph_Num)
kidney$BW <- as.numeric(kidney$BW)
kidney$KW <- as.numeric(kidney$KW)
kidney$GlomSA <- as.numeric(kidney$GlomSA)

kidney[!is.na(kidney$Neph_Num),"Log10Neph_Num"]=log10(kidney[!is.na(kidney$Neph_Num),"Neph_Num"])
kidney$Log10BW=log10(kidney$BW)
kidney$Log10KW=log10(kidney$KW)
k1=na.omit(kidney)

#Need to do regressions with either 1:body weight or 2: kidney weight for four variables, including: Neph_Num,GlomSA, ProxTubLen,ProxTubDiam.  
#Number nephrons
#mod1=lm(Log10Neph_Num ~ Log10BW, data=k1)
#summary(mod1)
#plot(k1$Log10BW, k1$Log10Neph_Num, log="xy")
#text(k1$Log10BW, k1$Log10Neph_Num,  labels=k1$Mammal) 
#abline(a=mod1$coefficients[1], b=mod1$coefficients[2], untf=TRUE)
mod2=lm(Log10Neph_Num ~ Log10KW, data=k1)
summary(mod2)
plot(k1$Log10KW, k1$Log10Neph_Num, log="xy")
text(k1$Log10KW, k1$Log10Neph_Num,  labels=k1$Mammal) 
abline(a=mod2$coefficients[1], b=mod2$coefficients[2], untf=TRUE)
Log10Neph_Num_pred=predict(mod2,data.frame(Log10KW=log10(kidney$KW))) 
plot(kidney$Log10KW, Log10Neph_Num_pred, log="xy")
text(kidney$Log10KW, Log10Neph_Num_pred,  labels=kidney$Mammal) 
kidney$Log10Neph_Num_pred=Log10Neph_Num_pred

#Glomerular surface area 
#mod1=lm(GlomSA ~ Log10BW, data=k1[-c(which(k1$Mammal=="Elephant" | k1$Mammal=="Whale")), ])
#summary(mod1)
#plot(k1$Log10BW, k1$GlomSA, log="x")
#text(k1$Log10BW, k1$GlomSA,  labels=k1$Mammal) 
#abline(a=mod1$coefficients[1], b=mod1$coefficients[2], untf=TRUE)
k1$Log10GlomSA <- log10(k1$GlomSA)
mod2 <- lm(Log10GlomSA ~ Log10KW, data=k1[-c(which(k1$Mammal=="Elephant" | k1$Mammal=="Whale")),])
# Glomerular surface area model accuracy:
summary(mod2)
plot(k1$Log10KW, k1$Log10GlomSA)
text(k1$Log10KW, k1$Log10GlomSA,  labels=k1$Mammal) 
abline(a=mod2$coefficients[1], b=mod2$coefficients[2], untf=TRUE)

Log10GlomSA_pred <- predict(mod2, data.frame(Log10KW=log10(kidney$KW))) 
plot(kidney$Log10BW, Log10GlomSA_pred, log="x")
text(kidney$Log10BW, Log10GlomSA_pred,  labels=kidney$Mammal) 
kidney$GlomSA_pred <- 10^Log10GlomSA_pred

plot(kidney$GlomSA, kidney$GlomSA_pred, log="xy")
text(kidney$GlomSA, kidney$GlomSA_pred,  labels=kidney$Mammal) 



#Note concerned about mouse and chicken. Mouse looks very small; Chicken lines up with monkey's and rabbits

#ProxTubLen
mod2=lm(ProxTubLen ~ Log10KW, data=k1[-c(which(k1$Mammal=="Whale")),])
summary(mod2)
plot(k1$Log10KW, k1$ProxTubLen, log="xy")
text(k1$Log10KW, k1$ProxTubLen,  labels=k1$Mammal) 
abline(a=mod2$coefficients[1], b=mod2$coefficients[2], untf=TRUE)

ProxTubLen_pred=predict(mod2,data.frame(Log10KW=log10(kidney$KW))) 
plot(kidney$Log10KW, ProxTubLen_pred, log="x")
text(kidney$Log10KW, ProxTubLen_pred,  labels=kidney$Mammal) 
kidney$ProxTubLen_pred=ProxTubLen_pred
points(k1$Log10KW, k1$ProxTubLen, col="red")
text(k1$Log10KW, k1$ProxTubLen,  labels=k1$Mammal, col = "red") 


#ProxTubDiam
plot(k1$Log10KW, k1$ProxTubDiam)
text(k1$Log10KW, k1$ProxTubDiam,  labels=k1$Mammal) 
mod2=lm(ProxTubDiam ~ Log10KW, data=k1)
summary(mod2)
plot(k1$Log10KW, k1$ProxTubDiam, log="xy")
text(k1$Log10KW, k1$ProxTubDiam,  labels=k1$Mammal) 
abline(a=mod2$coefficients[1], b=mod2$coefficients[2], untf=TRUE)

ProxTubDiam_pred=predict(mod2,data.frame(Log10KW=log10(kidney$KW))) 
plot(kidney$Log10KW, ProxTubDiam_pred, log="x")
text(kidney$Log10KW, ProxTubDiam_pred,  labels=kidney$Mammal) 
kidney$ProxTubDiam_pred=ProxTubDiam_pred
points(k1$Log10KW, k1$ProxTubDiam, col="red")
text(k1$Log10KW, k1$ProxTubDiam,  labels=k1$Mammal, col = "red") 
#For Proximal tubule length and diameter, the small and medium sized animals (Rat, Dog, and Rabbit) are pretty close predictions based on 
#the model, even with not a great fit. Chicken, Monkey, and Dog are in the same ball-park size wise, so hopefully
#that means there will be a similar correspondence between kidney size. The fit is worse
#for humans and large animals.  

plot(kidney$KW, kidney$BW)

#Note: Kidney type was altered based on a bit of reading about kidneys so that types were more consistent; essentially changed everything to either unipapillary or multirenculated. 

kidney$Neph_Num<-ifelse(is.na(kidney$Neph_Num), 10^(kidney$Log10Neph_Num_pred), kidney$Neph_Num)
kidney$GlomSA<-ifelse(is.na(kidney$GlomSA), kidney$GlomSA_pred, kidney$GlomSA)
kidney$ProxTubLen<-ifelse(is.na(kidney$ProxTubLen), kidney$ProxTubLen_pred, kidney$ProxTubLen)
kidney$ProxTubDiam<-ifelse(is.na(kidney$ProxTubDiam), kidney$ProxTubDiam_pred, kidney$ProxTubDiam)
kidney$KW_BW_ratio<-kidney$KW/kidney$BW
kidney$Neph_BW_ratio<-kidney$Neph_Num/kidney$BW
kidney$GlomTotSA<-kidney$GlomSA*kidney$Neph_Num
kidney$GlomTotSA_BW_ratio<-kidney$GlomTotSA/kidney$BW
kidney$GlomTotSA_KW_ratio<-kidney$GlomTotSA/kidney$KW
kidney$ProxTubVol<-kidney$ProxTubLen * (kidney$ProxTubDiam/2)^2 * pi 
kidney$ProxTubSA <- kidney$ProxTubLen * kidney$ProxTubDiam * pi
kidney$ProxTubTotalVol<-kidney$Neph_Num * kidney$ProxTubVol
kidney$ProxTubTotSA<-kidney$Neph_Num * kidney$ProxTubSA
kidney$GlomTotSA_ProxTubTotVol_ratio<-kidney$GlomTotSA/kidney$ProxTubTotalVol
kidney$ProxTubTotSA_ProxTotVol_ratio<-kidney$ProxTubTotSA/kidney$ProxTubTotalVol

###
#Join updated kidney stuff with df7
kidneynames=names(kidney)
kidneynames=kidneynames[which(kidneynames!="Mammal")]
df8=df7[,-c(which(names(df7)%in%kidneynames))]
df8=merge(df8, kidney, by.x="Species", by.y="Mammal", all.x=TRUE)
dim(df8)

save(kidney, file=paste0("RData/Kidney_Predictions_MultSpecies_", suff,".RData"))


#Lastly, add opera predictors;
#NOte, it looks like I'll have to grab average mass from outside this dataset
#Also, these are notes from the OPERA application documentation
#AD_Var(0/1):Global applicability domain considering the whole chemical space of the model
#AD_index[0-1]:Local applicability domain based on the similarity to the 5 nearest neighbors
#Conf_index[0]:Accuracy estimate based on the predictions of the 5 nearest neighbors

#Notes from Mansouri supplemental documentation about AD of models: 
#- If a chemical is considered outside the global AD with a low
#local AD-index, the prediction can be unreliable
#- If a chemical is considered outside the global AD but the local
#AD-index is average or relatively high, this means the query chemical is
#on the boundaries of the training set but has quite similar neighbors.
#The prediction can be trusted.
#- If a chemical is considered inside the global AD but the local
#AD-index is average or relatively low, this means the query chemical
#fell in a "gap" of the chemical space of the model but still within the
#boudaries of the training set and surrounded with training chemicals.
#The prediction should be considered with caution.
#- If a chemical is considered inside the global AD with a high
#local AD-index, the prediction should be considered reliable.

#From this, I'll use the criteria of either in the global domain or have a local index of at least 50%


#Add opera predictions enerated by OpERA program(stand alone application v2.6), not pulled from dashboard
opera=read.csv("Predictors/DSSToxPFAS_DTXSID_SMILES_Oct21-1-smi_OPERA2.7Pred.csv")
idvec=which(names(opera)%in%c("MoleculeID"))
op=which(grepl('\\AD_' , names(opera)) |  grepl('\\_pred$' , names(opera))) #Note, the grepl command here is looking for "_pred", you have to preceed this with \\ to look for this literal. Also, the $ is an end of line anchor
#notop=which(grepl('\\AD_index' , names(opera))) #Note, the grepl command here is looking for "_pred", you have to preceed this with \\ to look for this literal. Also, the $ is an end of line anchor
opvec=c(idvec, op)
#opvec=opvec[-c(which(opvec%in%notop))]
opera=opera[,opvec] #removing superfluous columns 
names(opera)[1]= "DTXSID"
#Make numeric
opera=data.frame("DTXSID"=opera[,1],apply(opera[,2:length(opera[1,])], 2, as.numeric))
##Fill in missing info for F53B with means of the rest
#opera[1,c(2:13)]<-colMeans(na.omit(opera[,2:13]))
#Take out an extra line of PFOS in there for some reason
opera=opera[which(opera$DTXSID%in%df8$DTXSID),]
#uknown column:SpAD_Dzm removed
opera=opera[,-c(which(names(opera)=="SpAD_Dzm"))]


#Determine applicability for modeling here based on AD information
models=c("LogP", "MP", "BP", "VP", "WS", "HL", "RT", "KOA","pKa", "LogD")
for(i in 1:length(models)){
  glob=opera[,which(names(opera)==paste0("AD_", models[i]))]
  local=opera[,which(names(opera)==paste0("AD_index_", models[i]))]
  opera=data.frame(opera,ifelse(c(glob==1 | local >=0.5),1,0))
  names(opera)[length(opera[1,])]=paste0("AD_QSAR_",models[i])}

df9=merge(df8, opera, by="DTXSID", keep="all")
names(df9)
dim(df9)
#Add in AverageMass
Mass=read.csv("Predictors/AverageMass_FromDashBoard_11PFAS_DED_04232020.csv")
names(Mass)[1]="DTXSID"

df10=merge(df9, Mass[,c(1,5)], by="DTXSID", keep="all")
dim(df10)

#Impute missing similarity values for some chemicals
endocols=which(grepl('\\TS' , names(df10)))
for(i in endocols){
  dfsub=df10[,i]
  dfsub[which(is.na(dfsub))]=mean(na.omit(dfsub))
  df10[,i]=dfsub
}
dim(df10)

#To account for ether bond, load in aliphatic COC toxprint
COCbond=read.xlsx("Predictors/toxprint_V2_vs_PFAS_list_Smiles_05_26_20.xlsx", sheet=2)
df11=merge(df10, COCbond, by="DTXSID", all= TRUE)
dim(df11)

#Lastly, add a column where halflife is scaled by lifetime
Lifespan=read.xlsx("Predictors/LifeSpan_from_Risa_05_25_20.xlsx")
Lifespan=subset(Lifespan, select=c("Species", "HrLifeSpan"))

df12=merge(df11, Lifespan, by="Species", all.x=TRUE, all.y=FALSE)
dim(df12)


#Add chain length as a predictor variable
##For PFAS chemicals with ether bond, including two chain lengths, TCChainLength, in which the oxygen
#is counted as a carbon(per Chris Lau), and LCChainLength, in which only the longest continously carbon
#chain is counted.

chainlist=as.data.frame(matrix(c("Perfluoroheptanoic acid",   7,7,
                                 "Perfluorooctanesulfonic acid",   8,8,
                                 "Perfluorohexanoic acid",   6,6,
                                 "Perfluorodecanoic acid",   10,10,
                                 "Perfluorononanoic acid",   9,9,
                                 "Perfluorohexanesulfonic acid",   6,6,
                                 "Perfluorobutanoic acid",   4,4,
                                 "F53B",   9,6,
                                 "Perfluorobutanesulfonic acid",   4,4,
                                 "Perfluorooctanoic acid",  8,8,
                                 "GenX" ,  3, 6), nrow=11, ncol=3, byrow = TRUE))
names(chainlist)=c("Chemical", "TCChainLength","LCChainLength" )
for(i in 1:length(chainlist[,1])){
  df12[which(df12$PREFERRED_NAME==chainlist[i,1]), "TCChainLength"]=chainlist[i,2]
  df12[which(df12$PREFERRED_NAME==chainlist[i,1]), "LCChainLength"]=chainlist[i,3]
}

df12$TCChainLength=as.numeric(df12$TCChainLength)
df12$LCChainLength=as.numeric(df12$LCChainLength)

pfasds=df12
#NOte: Check suffix before saving
suff
pfasds$TrainOrder=seq(1:length(pfasds[,1]))
save(pfasds, file=paste("RData/PFAS_11Chemicals_QSARdataset_", suff, ".RData      ", sep=""))
unique(pfasds$DTXSID)



#####Generate PFAS List if not present
unpf=unique(pfasds$DTXSID)
unpflist=NULL
for(i in unpf){
  sub=subset(pfasds, DTXSID == i, select=c("PREFERRED_NAME", "DTXSID"))
  unpflist=rbind(unpflist, sub[1,])}

#Generate Smiles codes for each PFAS
print(unique(df$DTXSID))

pfasds.unique <- subset(pfasds,!duplicated(DTXSID))
num.FABP <- 2
print(paste("Two FABP descriptors were available for",
  sum(!is.na(as.numeric(pfasds.unique$Kd_hL_FABP))),
  "structures."))              

num.physio <- sum((regexpr("Neph",colnames(pfasds.unique))!=-1 |
  regexpr("Prox",colnames(pfasds.unique))!=-1 |
  regexpr("Glom",colnames(pfasds.unique))!=-1 |
  regexpr("KW",colnames(pfasds.unique))!=-1) &
  regexpr("Log",colnames(pfasds.unique))==-1
  ) + 1 # Add BW
print(paste(num.physio,
  "physiological descriptors"))
  
num.CMC <- sum((regexpr("CMC",colnames(pfasds.unique))!=-1) &
  regexpr("Reference",colnames(pfasds.unique))==-1
  )
print(paste(num.CMC,
  "critical micellular concentration descriptors"))

num.endo <- sum(regexpr("TSM",colnames(pfasds.unique))!=-1 |
  regexpr("TSPC",colnames(pfasds.unique))!=-1 
  ) + 2 # Add max vals
print(paste(num.endo,
  "endongenous ligand descriptors"))

num.album <- sum((regexpr("SerAlb",colnames(pfasds.unique))!=-1) &
  regexpr("Reference",colnames(pfasds.unique))==-1
  )
print(paste(num.album,
  "serum albumin binding descriptors"))

num.physico <- sum(regexpr("AD",colnames(opera))==-1) -1 +4# -DTXSID +mass, +2chainlength, +ethane
print(paste(num.physico,
  "physico-chemical descriptors"))

num.cat <- 2
print(paste(num.cat, "categorical descriptors"))

print(paste(num.FABP + num.physio + num.CMC + num.endo + num.album + 
  num.physico + num.cat,
  "total descriptors considered"))
  
print(paste(num.FABP +num.physio + num.endo + num.album + num.physico,
  "examined with RFE"))
  

