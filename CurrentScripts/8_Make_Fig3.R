# Clear the workspace:
rm(list=ls())
try(dev.off())

readsuff <- "JFW100322-noLogD"
writesuff <- "JFW100422-noLogD"

library(ggplot2)

load(paste("RData/ClassyfierSubSets_",readsuff,".RData",sep=""))

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




pfascountCFSS=aggregate(Count~ ClassPredFull +DosingAdj+SpeciesSex , data=TCSub, FUN="sum")
Plot_Species_by_Sex=ggplot(data=pfascountCFSS[pfascountCFSS$DosingAdj=="Other",], aes(x=SpeciesSex, y=Count, fill=ClassPredFull)) +
  geom_bar(stat="identity")+
  xlab("Species/Sex")+
  ylab("Number of PFAS chemicals")+
  
  ggtitle(paste("RF Classification Model: Serum Half-Life of ",length(unique(TCSub$DTXSID)), " PFAS Chemicals in Half-life Model Domain",sep=""))+
  scale_fill_discrete(name="Serum Half Life", labels = c("<0.5 Day", "<Week", "< 2 Months", "> 2 Months"))+   
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

 
print(Plot_Species_by_Sex)
jpeg(paste("Figures/PFAS_DSSTox_",length(unique(TCSub$DTXSID)),"_TrainingClassSubset_ClassPred_ClassMod_Sp_Sex_",writesuff,".jpg",sep=""), width=800, height = 668)
print(Plot_Species_by_Sex)
try(dev.off())

#Human stats
HAMAD=TCSub[TCSub$Species=="Human" & TCSub$DosingAdj=="Other" & TCSub$Sex=="Female",]
table(HAMAD$ClassPredFull)/sum(table(HAMAD$ClassPredFull))


#Rats stats
RAMADF=TCSub[TCSub$Species=="Rat" & TCSub$DosingAdj=="Other" & TCSub$Sex=="Female",]
RAMADM=TCSub[TCSub$Species=="Rat" & TCSub$DosingAdj=="Other" & TCSub$Sex=="Male",]

table(RAMADF$ClassPredFull)/sum(table(RAMADF$ClassPredFull))
table(RAMADM$ClassPredFull)/sum(table(RAMADF$ClassPredFull))



