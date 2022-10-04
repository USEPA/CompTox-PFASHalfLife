library(openxlsx)
library(fitdistrplus)
#setwd("L:/Lab/NCCT_ExpoCast/ExpoCast2020/Dawson_TK_QSAR_PFAS/PFAS_HL_QSAR/Manuscript")
suff="DED042722"
setwd("L:/Lab/NCCT_ExpoCast/ExpoCast2021/Dawson_PFAS_HALFLIFE/PFAS_HL_QSAR_2021")
dat=read.xlsx("PFAS_Tables_For_Manuscripts_DED062121.xlsx", sheet=3)

#Eliminate missing species/chemical combos from the dataset
dat=dat[-c(which(is.na(dat$Mean))),]

chem=unique(dat$Chemical)
sp=unique(dat$Species)
sex=c("Male", "Female")
dat$DosingAdj=ifelse(dat$Dosing!="Oral" & dat$Dosing!="IV", "Other", dat$Dosing )
dosing=unique(dat$DosingAdj)

#Translate SD's and CI95's into SE's. Note, things presented as ranges are treated as standard errors
dat$SE= ifelse(dat$Boudary=="SE"| dat$Boudary=="Range", dat$Mean-dat$Lower, 
               ifelse(dat$Boudary=="SD", (dat$Mean-dat$Lower)/sqrt(dat$Study_N),
                      ifelse(dat$Boudary=="CI95", (dat$Mean-dat$Lower)/1.96, NA)))



#Translate mean and SE everything into hours
dat$HLH<-ifelse(dat$Study.Unit=="Days", dat$Mean*24, ifelse(dat$Study.Unit=="Years", dat$Mean*24*365,dat$Mean))
dat$HLH_SE<-ifelse(dat$Study.Unit=="Days", dat$SE*24, ifelse(dat$Study.Unit=="Years", dat$SE*24*365,dat$SE))




#Find all combos of interest
paramtab=expand.grid(chem, sp, sex, dosing)
names(paramtab)=c("Chemical", "Species", "Sex", "Dosing")



modlist=NULL
meantab=NULL
#First, Loop through entire list to determine which ones need which kind of processing rule, and
#make means for other cases
for(i in 1:length(paramtab[,1])){
  sub=subset(dat, dat$Chemical==paramtab[i, "Chemical"] & dat$Species==paramtab[i, "Species"] & dat$Sex==paramtab[i, "Sex"] & dat$DosingAdj==paramtab[i, "Dosing"])
  if(dim(sub)[1]==0) {next} 
  else if (dim(sub)[1]==1){
    sub$NObs=1
    meantab=rbind(meantab, sub)}
    else if (dim(sub)[1] > 1 &  any(which(sub$SE>0))==FALSE){
      sub1=sub[1,]
      sub1$NObs=dim(sub)[1]
      sub1$Mean=mean(sub$Mean)
      meantab=rbind(meantab, sub1)
    }else{
    sub1=paramtab[i,]
    sub1$NObs=dim(sub)[1]
    modlist=rbind(modlist, sub1)
    }}


#Now,take develop distributions for the remaining classes
#Have a couple of situations. First, everything has a HLH_SE, so we can simply sample
#across distributions, combine the distributions, find the a mean. 


#Next, fill in values for ranges of values where low and high values are missing
moddat=NULL
for (i in 1:length(modlist[,1])){
  sub=subset(dat, dat$Chemical==modlist[i, "Chemical"] & dat$Species==modlist[i, "Species"] & dat$Sex==modlist[i, "Sex"] & dat$DosingAdj==modlist[i, "Dosing"])
  nabound=which(is.na(sub$Lower))
  if(length(nabound)==0){
  moddat=rbind(moddat, sub)} else {
  sub1=sub[which(!is.na(sub$Lower)),]
  Nof1=which(sub1$Study_N==1)
  if(length(Nof1)>0){
  Nof1tab=sub1[Nof1,]
  sub1=sub1[-Nof1,]
  Nof1tab$HLH_SE=0
  moddat=rbind(moddat, Nof1tab)
  }
  sub1$SD=sub1$HLH_SE * sqrt(sub1$Study_N)  
  sub1SD=mean(sub1$SD)
  sub$HLH_SE[nabound]=sub1SD/sqrt(sub$Study_N[nabound])
  moddat=rbind(moddat, sub)
    }
  }
  

#For each of the chemical with multiple points, sample across distributions
#Combine together 
MeanDat=subset(meantab, select=c("Chemical", "Sex", "Species", "DosingAdj","Phase", "OriginalValue", "Orignal.Unit", "HLH", "HLH_SE","Source", "NObs")) 
names(MeanDat)[which(names(MeanDat)=="HLH_SE")]="HLH_SE_Range"
MeanDat$FittedDist=NA
MeanDat$Multimodal=NA
MeanDat$HLH_ME=NA
MeanDat$HLHSE_ME=NA
MeanDat$HLHSS_ME=NA
for (i in 1:length(modlist[,1])){
sub=subset(moddat, moddat$Chemical==modlist[i, "Chemical"] & moddat$Species==modlist[i, "Species"] & moddat$Sex==modlist[i, "Sex"] & moddat$DosingAdj==modlist[i, "Dosing"])
sub$HLH_SE_Range = sub$HLH_SE

UB=sub$HLH+sub$HLH_SE_Range
LB=sub$HLH-sub$HLH_SE_Range
Multimodal=ifelse(min(UB)<max(LB), "Yes", "No") #This notates whether the distributions
#don't over lap. This is mainly a concern when you've got only two distributions. 

distvec=NULL
set.seed(12345)
for(k in 1:100){ #This runs through 100 times, sampling within distributions in proportion to study sample size
for(j in 1:dim(sub)[1]){
  LBadjust=ifelse(sub$HLH[j]-sub$HLH_SE_Range[j]<0, 0,sub$HLH[j]-sub$HLH_SE_Range[j])
  distvec=c(distvec,runif(sub$Study_N[j], LBadjust, sub$HLH[j]+sub$HLH_SE_Range[j]))
}}


 norm=fitdist(distvec, distr = "norm") 
 lnorm=fitdist(distvec, distr = "lnorm")

 tab=gofstat(list(norm, lnorm))
 stat=which(tab$aic==min(tab$aic))
 sub$HLH=ifelse(stat==1, norm$estimate[1], exp(lnorm$estimate[1]))
 sub$HLH_ME=ifelse(stat==1, norm$estimate[1], lnorm$estimate[1])
 sub$HLHSE_ME=ifelse(stat==1, norm$estimate[2], lnorm$estimate[2])
 sub$HLHSS_ME=length(distvec)
 
 
 sub$Source=paste("Multiple")
 sub$FittedDist=ifelse(stat==1, "norm","lnorm")
 sub$NObs=dim(sub)[1]
 sub$Multimodal=Multimodal
 sub=sub[1,]
 sub=subset(sub, select=c("Chemical", "Sex", "Species", "DosingAdj","Phase", "OriginalValue", "Orignal.Unit", "HLH", "HLH_SE_Range", "Source", "FittedDist", "NObs", "Multimodal", "HLH_ME", "HLHSE_ME", "HLHSS_ME")) 
 MeanDat=rbind(MeanDat, sub)
}

MeanDat=MeanDat[,-c(which(names(MeanDat)=="HLH_SE_Range"))]
head(MeanDat)
tail(MeanDat)

unique(MeanDat$Chemical)


#####
#Old dataset
load("L:/Lab/NCCT_ExpoCast/ExpoCast2020/Dawson_TK_QSAR_PFAS/PFAS_HL_QSAR/ScriptsDataOutput/PFAS_HL_QSAR/PFAS_11Chemicals_QSARdataset_DED04222020_TBG.RData")
#Use old PFAS set to associate CASRN's and DTXID's
unpfas=unique(pfasds$DTXSID)
unpfaslist=NULL
for(i in unpfas){
  sub=subset(pfasds, DTXSID==i, select=c("DTXSID", "CASRN", "PREFERRED_NAME"))
  sub=sub[1,]
  unpfaslist=rbind(unpfaslist, sub)
}

#
unique(MeanDat$Chemical)
unpfaslist
mergenames=c("PFHpA (C7)", "PFOS (C8)", "GenX", 
             "PFOA (C8)", "PFNA (C9)", "PFDA (C10)", 
             "PFHxS (C6)", "PFHxA (C6)", "F-53B", "PFBS (C4)", "PFBA (C4)" )
unpfaslist$mergenames=mergenames

MeanDatMerge=merge(MeanDat, unpfaslist, by.x="Chemical", by.y="mergenames", all.y=TRUE, all.x=FALSE)

names(MeanDatMerge)

#Reduce and re-organize MeanDatMerge columns

MeanDatMerge1=subset(MeanDatMerge, select=c("CASRN", "DTXSID", "PREFERRED_NAME", "Species", "Sex", "DosingAdj","HLH", "HLH_ME", "HLHSE_ME", "HLHSS_ME", "Source","NObs", "Multimodal", "FittedDist"))
write.xlsx(MeanDatMerge1, file=paste0("PFAS_HLH_Revised_", suff,".xlsx"))
options(scipen=999)


##to find arithmetic mean from log transformed values

a=distvec
mean(distvec)
exp(mean(a)+0.5*var(a))


(exp(var(a))-1) * exp(2*mean(a)+var(a))
var(distvec)


a=rlnorm(900, 1, 0.4)
hist(a)
exp(mean(a))
mean(exp(a))
hist(exp(a))



##Note, that the distribution we've got there samples from individual studies via the uniform distribution, and then puts 
##all of the studies together, fits a distribution to it, and then samples from that. 
