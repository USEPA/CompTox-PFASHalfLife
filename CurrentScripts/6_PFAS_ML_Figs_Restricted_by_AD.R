# Clear the workspace:
rm(list=ls())
try(dev.off())

#setwd("L:/Lab/NCCT_ExpoCast/ExpoCast2021/Dawson_PFAS_HALFLIFE/PFAS_HL_QSAR_2021")
writesuff <- "JFW100322-noLogD"
suff <- "JFW100322-noLogD" #Note, this is the suffix for previous work on this page
seed <- "12345"
packages=c("readxl","MLmetrics","stringr","scales","ggplot2")
sapply(packages, require,character.only=TRUE) #Note, the "character.only" argument is necessary her

#In this script, I am going to restrict the domain of the predicted chemicals to only those that fall under
#all of the domain. Then, I'm going to make new figures and describe the chemical characteristics
#of the chemicals in domain.

load(file=paste0("RData/Tox21_AllMods_ADindicate_", suff,".RData"))

#Limit the complete dataset to only the species in the list 
ds=subset(DPFASwAD, Species%in%c("Human", "Monkey", "Mouse", "Rat", "Dog"))
ds$SpeciesSex=paste0(ds$Species, "_", ds$Sex)
ds$SpeciesDose=paste0(ds$Species, "_", ds$DosingAdj)
ds$SpeciesDose=paste0(ds$Species, "_", ds$DosingAdj)
ds=ds[ds$ClassModDomain==1,]
ds$Count=1
dsAMAD=ds[ds$AD_LogP==1 & ds$AD_WS==1 & ds$AD_VP==1 & ds$ClassModDomain==1,]
write.csv(DPFASwAD, file=paste0("RData/All_DSSTox_Predictions_With_Domains_",writesuff,".csv"))


#Human chemical subset
human=dsAMAD[dsAMAD$Species=="Human",]
data.frame(unique(human$DTXSID), row.names = NULL)

#All chemicals
allchems=data.frame(unique(ds$DTXSID), row.names = NULL)
write.csv(allchems, file=paste0("RData/AllChemicals_in_CatMod_Domain_",writesuff,".csv"))

########
#Insert male and female symbols
library(showtext)
#Labels
female = intToUtf8(9792)
male = intToUtf8(9794)
IV="Iv"
Oral="Or"
Other="Ot"
#######

####Numbers of chemicals by bin
namelabelsSex=c(paste0("Dog:",female), paste0("Dog:",male), paste0("Human:",female), paste0("Human:",male), paste0("Monkey:",female), paste0("Monkey:",male), paste0("Mouse:",female), paste0("Mouse:",male), paste0("Rat:",female), paste0("Rat:",male))
namelabelsDose=c(paste0("D:", IV), paste0("D:", Oral), paste0("D:", Other), paste0("H:",IV), paste0("H:",Oral),paste0("H:",Other), paste0("Mn:",IV), paste0("Mn:",Oral),paste0("Mn:",Other), paste0("Mo:",IV), paste0("Mo:",Oral),paste0("Mo:",Other), paste0("Rat:",IV), paste0("Rat:",Oral),paste0("Rat:",Other))



##RF Classification Model 
###First, figures of chemicals in domain of 1/2 life model without regard to the underlying OPERA models
pfascountCFSS=aggregate(Count~ ClassPredFull +DosingAdj+SpeciesSex , data=ds, FUN="sum")

pfascountCFSD=aggregate(Count~ ClassPredFull + SpeciesDose + Sex , data=ds, FUN="sum")

Plot_Species_by_Sex=ggplot(data=pfascountCFSS[pfascountCFSS$DosingAdj=="Other",], aes(x=SpeciesSex, y=Count, fill=ClassPredFull)) +
  geom_bar(stat="identity")+
  xlab("Species/Sex")+
  ylab("Number of PFAS chemicals")+
  
  ggtitle(paste("RF Classification Model: Serum Half-Life of ",length(unique(ds$CASRN)), " PFAS Chemicals in Half-life Model Domain",sep=""))+
  scale_fill_discrete(name="Serum Half Life", labels = c("<1 Day", "<Week", "< 2 Months", "> 2 Months"))+   
  theme(
    plot.title=element_text(size =15), 
    # plot.subtitle = element_text(size=17),
    axis.text.x =element_text(size=15, vjust=0.3, hjust=0.5),
    axis.text.y=element_text(size=15),
    axis.title=element_text(size=17),
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
png(paste("Figures/PFAS_DSSTox_",length(unique(ds$CASRN)),"chemicals_ClassPred_ClassMod_Sp_Sex_",writesuff,".png",sep=""), width=800, height = 668)
Plot_Species_by_Sex
dev.off()

####Species By Dose
#Note: reducing the plot to only include male sex to avoid replication
Plot_Species_by_Dose=ggplot(data=pfascountCFSD[pfascountCFSD$Sex=="Male",], aes(x=SpeciesDose, y=Count, fill=ClassPredFull)) +
  geom_bar(stat="identity")+
  xlab("Species & Dose")+
  ylab("Number of PFAS chemicals")+
  
  ggtitle(paste("RF Classification Model: Serum Half-Life of ",length(unique(ds$CASRN)), " PFAS Chemicals in Half-Life Model Domain",sep=""))+
  scale_fill_discrete(name="Serum Half Life", labels = c("<1 Day", "<Week", "< 2 Months", "> 2 Months"))+   
  theme(
    plot.title=element_text(size =15), 
    # plot.subtitle = element_text(size=17),
    axis.text.x =element_text(size=12, vjust=0.3, hjust=0.5),
    axis.text.y=element_text(size=15),
    axis.title=element_text(size=17),
    legend.position="bottom",
    legend.text=element_text(size=15),
    legend.title=element_text(size=20))
Plot_Species_by_Dose=
  Plot_Species_by_Dose+  
  scale_x_discrete(labels= c("Dog_IV" = namelabelsDose[1], "Dog_Oral"=namelabelsDose[2], "Dog_Other"= namelabelsDose[3],
                             "Human_IV"=namelabelsDose[4], "Human_Oral"=namelabelsDose[5],"Human_Other"=namelabelsDose[5],
                             "Monkey_IV"=namelabelsDose[7], "Monkey_Oral"=namelabelsDose[8], "Monkey_Other"=namelabelsDose[9],
                             "Mouse_IV"=namelabelsDose[10], "Mouse_Oral"=namelabelsDose[11], "Mouse_Other"=namelabelsDose[12],
                             "Rat_IV"=namelabelsDose[13], "Rat_Oral"=namelabelsDose[14],"Rat_Other"=namelabelsDose[15]))




Plot_Species_by_Dose

png(paste("Figures/PFAS_DSSTox_",length(unique(ds$CASRN)),"chemicals_ClassPred_ClassMod_Sp_Dose_",writesuff,".png",sep=""), width=800, height = 668)
Plot_Species_by_Dose
dev.off()


#######Next inclusion only chemicals in all model domains(AMAD)
pfascountCFSS=aggregate(Count~ ClassPredFull +DosingAdj+SpeciesSex , data=dsAMAD, FUN="sum")
pfascountCFSD=aggregate(Count~ ClassPredFull + SpeciesDose + Sex , data=dsAMAD, FUN="sum")
pfascount=aggregate(Count~ ClassPredFull +Species +Sex +DosingAdj , data=dsAMAD, FUN="sum")
pfascount1=aggregate(Count~ ClassPredFull + Sex+ DosingAdj+Species, data=dsAMAD, FUN="sum")

#Human stats
HAMAD=dsAMAD[dsAMAD$AD_LogP==1 & dsAMAD$AD_VP==1 & dsAMAD$AD_WS==1 & dsAMAD$ClassModDomain==1 & dsAMAD$Species=="Human" & dsAMAD$DosingAdj=="Other" & dsAMAD$Sex=="Female",]
table(HAMAD$ClassPredFull)/sum(table(HAMAD$ClassPredFull))

#Rats stats
RAMADM=dsAMAD[dsAMAD$AD_LogP==1 & dsAMAD$AD_VP==1 & dsAMAD$AD_WS==1 & dsAMAD$ClassModDomain==1 & dsAMAD$Species=="Rat" & dsAMAD$DosingAdj=="Other" & dsAMAD$Sex=="Male",]
RAMADF=dsAMAD[dsAMAD$AD_LogP==1 & dsAMAD$AD_VP==1 & dsAMAD$AD_WS==1 & dsAMAD$ClassModDomain==1 & dsAMAD$Species=="Rat" & dsAMAD$DosingAdj=="Other" & dsAMAD$Sex=="Female",]
table(RAMADF$ClassPredFull)/sum(table(RAMADF$ClassPredFull))
table(RAMADM$ClassPredFull)/sum(table(RAMADF$ClassPredFull))


#Dose Stats
dsAMAD$DoseBin=paste0(dsAMAD$DosingAdj,"_", dsAMAD$ClassPredFull)
HAMAD=dsAMAD[dsAMAD$AD_LogP==1 & dsAMAD$AD_VP==1 & dsAMAD$AD_WS==1 & dsAMAD$ClassModDomain==1 & dsAMAD$Species=="Human" & dsAMAD$Sex=="Female",]
table(HAMAD$DoseBin)/2764




#Plots
#Note: reducing the plot to only include the "Other" dosing methodology to avoid replication
Plot_Species_by_Sex=ggplot(data=pfascountCFSS[pfascountCFSS$DosingAdj=="Other",], aes(x=SpeciesSex, y=Count, fill=ClassPredFull)) +
  geom_bar(stat="identity")+
  xlab("Species/Sex")+
  ylab("Number of PFAS chemicals")+
  
  ggtitle(paste("RF Classification Model: Serum Half-Life of ",length(unique(dsAMAD$CASRN)), " PFAS Chemicals in All Model Domains",sep=""))+
  scale_fill_discrete(name="Serum Half Life", labels = c("< 12 hrs", "12 hrs - 1 wk", "1 wk - 2 mth", "> 2 mth"))+   
  theme(
    plot.title=element_text(size =15), 
    # plot.subtitle = element_text(size=17),
    axis.text.x =element_text(size=15, vjust=0.3, hjust=0.5),
    axis.text.y=element_text(size=15),
    axis.title=element_text(size=17),
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
png(paste("Figures/PFAS_DSSTox_",length(unique(dsAMAD$CASRN)),"chemicals_ClassPred_AMAD_Sp_Sex_",writesuff,".png",sep=""), width=800, height = 668)
Plot_Species_by_Sex
dev.off()

####Species By Dose
#Note: reducing the plot to only include male sex to avoid replication
Plot_Species_by_Dose=ggplot(data=pfascountCFSD[pfascountCFSD$Sex=="Male",], aes(x=SpeciesDose, y=Count, fill=ClassPredFull)) +
  geom_bar(stat="identity")+
  xlab("Species & Dose")+
  ylab("Number of PFAS chemicals")+
  
  ggtitle(paste("RF Classification Model: Serum Half-Life of ",length(unique(dsAMAD$CASRN)), " PFAS Chemicals in All Model Domains",sep=""))+
  scale_fill_discrete(name="Serum Half Life", labels = c("<1 Day", "<Week", "< 2 Months", "> 2 Months"))+   
  theme(
    plot.title=element_text(size =15), 
    # plot.subtitle = element_text(size=17),
    axis.text.x =element_text(size=12, vjust=0.3, hjust=0.5),
    axis.text.y=element_text(size=15),
    axis.title=element_text(size=17),
    legend.position="bottom",
    legend.text=element_text(size=15),
    legend.title=element_text(size=20))
Plot_Species_by_Dose=
  Plot_Species_by_Dose+  
  scale_x_discrete(labels= c("Dog_IV" = namelabelsDose[1], "Dog_Oral"=namelabelsDose[2], "Dog_Other"= namelabelsDose[3],
                              "Human_IV"=namelabelsDose[4], "Human_Oral"=namelabelsDose[5],"Human_Other"=namelabelsDose[6],
                             "Monkey_IV"=namelabelsDose[7], "Monkey_Oral"=namelabelsDose[8], "Monkey_Other"=namelabelsDose[9],
                             "Mouse_IV"=namelabelsDose[10], "Mouse_Oral"=namelabelsDose[11], "Mouse_Other"=namelabelsDose[12],
                             "Rat_IV"=namelabelsDose[13], "Rat_Oral"=namelabelsDose[14],"Rat_Other"=namelabelsDose[15]))




Plot_Species_by_Dose

png(paste("Figures/PFAS_DSSTox_",length(unique(dsAMAD$CASRN)),"chemicals_ClassPred_AMAD_Sp_Dose_",writesuff,".png",sep=""), width=800, height = 668)
Plot_Species_by_Dose
dev.off()




############Older Code#####Regression, Life-span scaled approaches
# ##Descriptive statistics of half-lives per lifespan
# #Load median bin values and scale chemicals by life span of species 
# load(paste0("Median_HLH_per_Bin_4_DED061621.RData"))
# ds$ClassHLHPred=ifelse(ds$ClassPredFull==1, agtabHLH[1,2], ifelse(ds$ClassPredFull==2, agtabHLH[2,2],ifelse(ds$ClassPredFull==3, agtabHLH[3,2],agtabHLH[4,2])))
# ds$ClassHLHPred_PerLS=ds$ClassHLHPred/ds$PerLS
# 
# #Quantiles of chemicals
# #Classification Model 
# quantile(ds[ds$Species=="Human", which(names(ds)=="ClassHLHPred_PerLS")] )
# quantile(ds[ds$Species=="Monkey", which(names(ds)=="ClassHLHPred_PerLS")] )
# quantile(ds[ds$Species=="Rat", which(names(ds)=="ClassHLHPred_PerLS")] )
# quantile(ds[ds$Species=="Mouse", which(names(ds)=="ClassHLHPred_PerLS")] )
# quantile(ds[ds$Species=="Dog", which(names(ds)=="ClassHLHPred_PerLS")],probs = seq(0, 1, 0.1) )
# 
# #Regression models
# quantile(ds[ds$Species=="Human", which(names(ds)=="RegpredRed_byLS")],probs = seq(0, 1, 0.02) )
# quantile(ds[ds$Species=="Monkey", which(names(ds)=="RegpredRed_byLS")],probs = seq(0, 1, 0.02) )
# quantile(ds[ds$Species=="Rat", which(names(ds)=="RegpredRed_byLS")],probs = seq(0, 1, 0.02) )
# quantile(ds[ds$Species=="Mouse", which(names(ds)=="RegpredRed_byLS")],probs = seq(0, 1, 0.02) )
# quantile(ds[ds$Species=="Dog", which(names(ds)=="RegpredRed_byLS")], probs = seq(0, 1, 0.1) )
# 
# PredictionPlotCF_LFScaled=ggplot(data = ds, aes(x = Species,  y = ClassHLHPred_PerLS)) +
#   labs(x = "Species/Sex", 
#        y = "Serum Half-Lives as a Proportion of Average Lifespan")+ 
#   ggtitle("RF Classification Model: Distribution of Predicted Serum Half-Lives\n as Proportion of Average Lifespan", 
#           subtitle = paste0(length(unique(ds$CASRN)), " PFAS Chemicals of DSSTox List" )) +
#   theme(plot.title = element_text(size = 15, face = "bold"), 
#         plot.subtitle = element_text(size=17)
#         axis.title.y = element_text(size = 20),
#         axis.title.x = element_text(size = 20),
#         axis.text=element_text(size=15))+
#   #scale_y_log10(limits=c(1e-6, 1)) +
#   ylim(0,1)+
#   geom_boxplot()
# 
# PredictionPlotCF_LFScaled
# png(paste("Figures/PFAS_DSSTox_",length(unique(completedataset$CASRN)),"chemicals_ClassPred_byLifeSpan",writesuff,".png",sep=""), width=800, height = 668)
# PredictionPlotCF_LFScaled
# dev.off()


# ######Regression Models 
# ###Prediction by Full Regression MOdel 
# PredictionPlotRF=ggplot(data = ds, aes(x = Species,  y = RegPredFull)) +
#   labs(x = "Species", 
#        y = "Serum Half-Lives in Hours")+ 
#   ggtitle("RF Full Regression Model: Distribution of Predicted Serum Half-Lives (Hrs)", 
#           subtitle = paste0(length(unique(ds$CASRN)), " PFAS Chemicals of DSSTox List" )) +
#   theme(plot.title = element_text(size = 15, face = "bold"), 
#         axis.title.y = element_text(size = 12),
#         axis.title.x = element_text(size = 12))+
#   #scale_y_log10(label=comma)+
#   scale_y_continuous(labels = comma)  +
#   #  scale_y_log10() +
#   geom_boxplot()
# 
# PredictionPlotRF
# png(paste("Figures/PFAS_DSSTox_",length(unique(completedataset$CASRN)),"chemicals_RegPredFull_",writesuff,".png",sep=""), width=800, height = 668)
# PredictionPlotRF
# dev.off()
# 
# ###Scaled by halflife 
# PredictionPlotRF_LS=ggplot(data = ds, aes(x = Species,  y = RegpredFull_byLS)) +
#   labs(x = "Species", 
#        y = "Serum Half-Lives as a Proportion of Average Lifespan")+ 
#   ggtitle("RF Full Regression Model: Distribution of Predicted Serum Half-Lives (Hrs)", 
#           subtitle = paste0(length(unique(ds$CASRN)), " PFAS Chemicals of DSSTox List" )) +
#   theme(plot.title = element_text(size = 15, face = "bold"), 
#         axis.title.y = element_text(size = 12),
#         axis.title.x = element_text(size = 12))+
#   
#   scale_y_log10() +
#   geom_boxplot()
# 
# PredictionPlotRF_LS
# png(paste("Figures/PFAS_DSSTox_",length(unique(completedataset$CASRN)),"chemicals_RegPredFull_byLifeSpan_",writesuff,".png",sep=""), width=800, height = 668)
# PredictionPlotRF_LS
# dev.off()



# ###Prediction by Reduced Regression MOdel 
# library(scales)
# namelabels1=c("Dog:IV", "Dog:Oral","Dog:Other","Human:IV", "Human:Oral","Human:Other","Monkey:IV", "Monkey:Oral","Monkey:Other","Mouse:IV", "Mouse:Oral","Mouse:Other","Rat:IV", "Rat:Oral","Rat:Other")
# PredictionPlotRR_Raw=ggplot(data = ds, aes(x = SpeciesDose,  y = RegPredRed)) +
#   labs(x = "Species and Dosing Method", 
#        y = "Serum Half-Lives in Hours")+ 
#   ggtitle("RF Reduced Regression Model: Distribution of Predicted Serum Half-Lives (Hrs)", 
#           subtitle = paste0(length(unique(ds$CASRN)), " PFAS Chemicals of DSSTox List" )) +
#   theme(
#     plot.title=element_text(size =20), 
#     plot.subtitle = element_text(size=17),
#     axis.text.x =element_text(angle=90, size=5, vjust=0.3, hjust=1),
#     axis.text.y=element_text(size=15),
#     axis.title=element_text(size=17),
#     legend.position="bottom",
#     legend.text=element_text(size=15),
#     legend.title=element_text(size=20))+
#   scale_y_log10(label=scientific_10,name="Serum Half-life (Hrs)", 
#                 sec.axis = sec_axis(trans=~./8760, name="Serum Half-life (Yrs)")) +
#   #  scale_x_discrete(labels= namelabels1)
#   #  scale_y_continuous(labels = comma, name="Serum Half-lives in Hours", 
#   #                    sec.axis = sec_axis(trans=~./8760, name="Serum Half-lives in Years"))+ 
#   geom_boxplot()
# PredictionPlotRR_Raw
# 
# png(paste("Figures/PFAS_DSSTox_",length(unique(completedataset$CASRN)),"chemicals_RegPredRed",writesuff,".png",sep=""), width=800, height = 668)
# PredictionPlotRR_Raw
# dev.off()
# 
# ###Scaled by halflife 
# PredictionPlotRR_LS=ggplot(data = ds, aes(x = Species,  y = RegpredRed_byLS)) +
#   labs(x = "Species", 
#        y = "Serum Half-Lives as a Proportion of Average Lifespan")+ 
#   ggtitle("RF Reduced Regression Model: Distribution of Predicted Serum Half-Lives (Hrs)", 
#           subtitle = paste0(length(unique(ds$CASRN)), " PFAS Chemicals of DSSTox List" )) +
#   theme(plot.title = element_text(size = 15, face = "bold"), 
#         axis.title.y = element_text(size = 20),
#         axis.title.x = element_text(size = 20),
#         axis.text = element_text(size=15))+
#   scale_y_log10() +
#   geom_boxplot()
# 
# PredictionPlotRR_LS
# png(paste("Figures/PFAS_DSSTox_",length(unique(completedataset$CASRN)),"chemicals_RegPredRed_byLifeSpan_",writesuff,".png",sep=""), width=800, height = 668)
# PredictionPlotRR_LS
# dev.off()
# 
# #Quantile of mouse 1/2 as function of reproductive half-life. 
# quantile(ds[ds$Species=="Mouse", "RegpredRed_byLS"])



# #Load median bin values and scale chemicals by life span of species 
# PredictionPlotCF_LFScaled_InDomain=ggplot(data = sub2, aes(x = Species,  y = ClassHLHPred_PerLS)) +
#   labs(x = "Species", 
#        y = "Serum Half-Lives as a Proportion of Average Lifespan")+ 
#   ggtitle("RF Classification Model: Distribution of Predicted Serum Half-Lives\n as Proportion of Average Lifespan",
#           subtitle = paste0("Of ", length(unique(ds$CASRN)), " PFAS Chemicals of DSSTox List, ", length(unique(sub2$CASRN)), " in Model Domain"))+ 
#   theme(plot.title = element_text(size = 15, face = "bold"), 
#         axis.title.y = element_text(size = 12),
#         axis.title.x = element_text(size = 12))+
#   #  scale_y_log10(limits=c(1e-6, 1)) +
#   geom_boxplot()
# 
# PredictionPlotCF_LFScaled_InDomain
# png(paste("Figures/PFAS_DSSTox_",length(unique(completedataset$CASRN)),"chemicals_ClassPred_byLifeSpan_InDomain_",writesuff,".png",sep=""), width=800, height = 668)
# PredictionPlotCF_LFScaled_InDomain
# dev.off()


# ######Regression Models 
# ###Prediction by Full Regression MOdel 
# PredictionPlotRF_Raw_InDomain=ggplot(data = sub1, aes(x = Species,  y = RegPredFull)) +
#   labs(x = "Species", 
#        y = "Serum Half-Lives in Hours")+ 
#   ggtitle("RF Classification Model: Distribution of Predicted Serum Half-Lives\n as Proportion of Average Lifespan",
#           subtitle = paste0("Of ", length(unique(ds$CASRN)), " PFAS Chemicals of DSSTox List, ", length(unique(sub2$CASRN)), " in Model Domain"))+ 
#   theme(plot.title = element_text(size = 15, face = "bold"), 
#         axis.title.y = element_text(size = 12),
#         axis.title.x = element_text(size = 12))+
#   scale_y_log10() +
#   geom_boxplot()
# 
# PredictionPlotRF_Raw_InDomain
# png(paste("Figures/PFAS_DSSTox_",length(unique(completedataset$CASRN)),"chemicals_RegPredFull_InDomain_",writesuff,".png",sep=""), width=800, height = 668)
# PredictionPlotRF_Raw_InDomain
# dev.off()

# ###Scaled by halflife 
# PredictionPlotRF_LS_InDomain=ggplot(data = sub1, aes(x = Species,  y = RegpredFull_byLS)) +
#   labs(x = "Species", 
#        y = "Serum Half-Lives as a Proportion of Average Lifespan")+ 
#   ggtitle("RF Full Regression Model: Distribution of Predicted Serum Half-Lives (Hrs)", 
#           subtitle = paste0("Of ", length(unique(ds$CASRN)), " PFAS Chemicals of DSSTox List, ", length(unique(sub2$CASRN)), " in Model Domain"))+ 
#   theme(plot.title = element_text(size = 15, face = "bold"), 
#         axis.title.y = element_text(size = 12),
#         axis.title.x = element_text(size = 12))+
#   scale_y_log10() +
#   geom_boxplot()
# 
# PredictionPlotRF_LS_InDomain
# png(paste("Figures/PFAS_DSSTox_",length(unique(completedataset$CASRN)),"chemicals_RegPredFull_byLifeSpan_InDomain",writesuff,".png",sep=""), width=800, height = 668)
# PredictionPlotRF_LS_InDomain
# dev.off()



# ###Prediction by Reduced Regression MOdel 
# PredictionPlotRR_Raw_InDomain=ggplot(data = sub3, aes(x = Species,  y = RegPredRed)) +
#   labs(x = "Species", 
#        y = "Serum Half-Lives in Hours")+ 
#   ggtitle("RF Reduced Regression Model: Distribution of Predicted Serum Half-Lives (Hrs)", 
#           subtitle = paste0("Of ", length(unique(ds$CASRN)), " PFAS Chemicals of DSSTox List, ", length(unique(sub3$CASRN)), " in Model Domain"))+ 
#   theme(plot.title = element_text(size = 15, face = "bold"), 
#         axis.title.y = element_text(size = 12),
#         axis.title.x = element_text(size = 12))+
#   scale_y_log10() +
#   geom_boxplot()
# 
# PredictionPlotRR_Raw_InDomain
# png(paste("Figures/PFAS_DSSTox_",length(unique(completedataset$CASRN)),"chemicals_RegPredRed_InDomain_",writesuff,".png",sep=""), width=800, height = 668)
# PredictionPlotRR_Raw_InDomain
# dev.off()

# ###Scaled by halflife 
# PredictionPlotRR_LS_InDomain=ggplot(data = sub3, aes(x = Species,  y = RegpredRed_byLS)) +
#   labs(x = "Species", 
#        y = "Serum Half-Lives as a Proportion of Average Lifespan")+ 
#   ggtitle("RF Reduced Regression Model: Distribution of Predicted Serum Half-Lives (Hrs)", 
#           subtitle = paste0("Of ", length(unique(ds$CASRN)), " PFAS Chemicals of DSSTox List, ", length(unique(sub3$CASRN)), " in Model Domain"))+ 
#   theme(plot.title = element_text(size = 15, face = "bold"), 
#         axis.title.y = element_text(size = 12),
#         axis.title.x = element_text(size = 12))+
#   scale_y_log10() +
#   geom_boxplot()
# 
# PredictionPlotRR_LS_InDomain
# png(paste("Figures/PFAS_DSSTox_",length(unique(completedataset$CASRN)),"chemicals_RegPredRed_byLifeSpan_InDomain_",writesuff,".png",sep=""), width=800, height = 668)
# PredictionPlotRR_LS_InDomain
# dev.off()


#range(ds[ds$Species=="Human", "RegpredRed_byLS"]) 
#(ds[ds$Species=="Mouse", "RegpredRed_byLS"]) 

#perfluoroundecnoic acid for Jeff Municci
#DPFASwAD[DPFASwAD$DTXSID=="DTXSID8047553" & DPFASwAD$Species=="Human", c("CASRN", "DTXSID", "ClassPredFull", "RegPredFull")]
