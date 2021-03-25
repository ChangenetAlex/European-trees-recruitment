library(reshape2)
rm(list = ls())
gc()

Dir <- c("~/SynologyDrive/FUNDIV - NFI - Europe/")
setwd(Dir)
source(paste0(Dir,"Myscripts/Fundiv.project/RecrutProject.function4.Modelsummary.R"))

## Apply the function:
Modname <- "rW.yO.sumO.mRY1"
#Modname <- "SCALEDRrate.yr.QualiF.TF.noY.nbinom1.V3" WARNING
Model.All(Modname)

## Load files (csv)
tableresults <- read.csv("~/SynologyDrive/FUNDIV - NFI - Europe/our-data/species/Recrut.All/Model_All.csv")
tableresults <- read.csv("~/SynologyDrive/FUNDIV - NFI - Europe/our-data/species/Recrut.All/Model_All.M.csv")
tableresults <- read.csv("~/SynologyDrive/FUNDIV - NFI - Europe/our-data/species/Recrut.All/Model_All_R.csv")


## Load files (rds)
tableresults <- readRDS(paste0("~/SynologyDrive/FUNDIV - NFI - Europe/our-data/species/Recrut.All/Model_All_",Modname,".rds"))
summary(as.factor(tableresults$variable))
## Simple effects 

#tableresults <- tableresults[tableresults$signif<=0.05,] # Keep only significant values

tableresults.zi <- tableresults[tableresults$model=="zi",] # Keep only zi model
tableresults.cond <- tableresults[tableresults$model=="cond",] # Keep only zi model
## Summary 
summary(as.factor(tableresults.zi[tableresults.zi$estimate<0,"variable"]))
summary(as.factor(tableresults.zi[tableresults.zi$estimate>0,"variable"]))
summary(as.factor(tableresults.zi[tableresults.zi$estimate==0,"variable"]))

summary(as.factor(tableresults.cond[tableresults.cond$estimate<0,"variable"]))
summary(as.factor(tableresults.cond[tableresults.cond$estimate>0,"variable"]))
summary(as.factor(tableresults.cond[tableresults.cond$estimate==0,"variable"]))

# Change names 
tableresults$variable <- gsub("logdbhJ.IMall.plot.mean", "DBH", tableresults$variable, ignore.case = FALSE, perl = FALSE,fixed = FALSE, useBytes = FALSE)
tableresults$variable <- gsub("logtreeNbrJ", "D", tableresults$variable, ignore.case = FALSE, perl = FALSE,fixed = FALSE, useBytes = FALSE)
tableresults$variable <- gsub("mean_spei12", "SPEI", tableresults$variable, ignore.case = FALSE, perl = FALSE,fixed = FALSE, useBytes = FALSE)
tableresults$variable <- gsub("Plotcat1", "TE", tableresults$variable, ignore.case = FALSE, perl = FALSE,fixed = FALSE, useBytes = FALSE)
tableresults$variable <- gsub("Plotcat2", "LE", tableresults$variable, ignore.case = FALSE, perl = FALSE,fixed = FALSE, useBytes = FALSE)
tableresults$variable <- gsub("sqrtBA.O.plot.1", "Inter", tableresults$variable, ignore.case = FALSE, perl = FALSE,fixed = FALSE, useBytes = FALSE)
tableresults$variable <- gsub("logBAj.plot.1", "Intra", tableresults$variable, ignore.case = FALSE, perl = FALSE,fixed = FALSE, useBytes = FALSE)
tableresults$variable <- gsub("sqrtBAIj.plot.1.mean.J.I", "G", tableresults$variable, ignore.case = FALSE, perl = FALSE,fixed = FALSE, useBytes = FALSE)
tableresults$variable <- gsub("logsp.mortality.plot.rate.yr", "M", tableresults$variable, ignore.case = FALSE, perl = FALSE,fixed = FALSE, useBytes = FALSE)
tableresults$variable <- gsub("(Intercept)", "Intercept", tableresults$variable, ignore.case = FALSE, perl = FALSE,fixed = FALSE, useBytes = FALSE)
tableresults$variable <- gsub(":", " * ", tableresults$variable, ignore.case = FALSE, perl = FALSE,fixed = FALSE, useBytes = FALSE)
unique(tableresults$variable)


tableresults.zi <- tableresults[tableresults$model=="zi",] # Keep only zi model
tableresults.cond <- tableresults[tableresults$model=="cond",] # Keep only zi model

# Remove intercept and Z-values  
#tableresults.zi <- tableresults.zi[tableresults.zi$variable!="(Intercept)",c(1:5,7)] # Keep only zi model
#tableresults.cond <- tableresults.cond[tableresults.cond$variable!="(Intercept)",c(1:5,7)] # Keep only cond model

# Remove z-values 
tableresults.zi <- tableresults.zi[,c(1:5,7)] # Keep only zi model
tableresults.cond <- tableresults.cond[,c(1:5,7)] # Keep only cond model

# Round values 
tableresults.zi[,c(4:6)] <- round(tableresults.zi[,c(4:6)],2)
tableresults.cond[,c(4:6)] <- round(tableresults.cond[,c(4:6)],2)

# Replace significiance by *
tableresults.zi[tableresults.zi$signif<=0.05,"signif"] <- "*"
tableresults.zi[tableresults.zi$signif>0.05,"signif"] <- ""
tableresults.cond[tableresults.cond$signif<=0.05,"signif"] <- "*"
tableresults.cond[tableresults.cond$signif>0.05,"signif"] <- ""


# Paste estimate and signif
tableresults.zi$final <- paste0(tableresults.zi$estimates,"\n(",tableresults.zi$SE,")\n",tableresults.zi$signif)
tableresults.cond$final <- paste0(tableresults.cond$estimates,"\n(",tableresults.cond$SE,")\n",tableresults.cond$signif)


# Delete useless column 
tableresults.zi <- tableresults.zi[,c(1,3,7)]
tableresults.cond <- tableresults.cond[,c(1,3,7)]

# Reorder the table by species or variable
tableresults.zi <- dcast(data=tableresults.zi, variable ~ species, drop = FALSE)
tableresults.zi[,"signifsum"] <- apply(tableresults.zi, 1, function(x) sum(grepl("*",x,fixed=T)==T))

tableresults.cond <- dcast(data=tableresults.cond, variable ~ species, drop = FALSE)
tableresults.cond[,"signifsum"] <- apply(tableresults.cond, 1, function(x) sum(grepl("*",x,fixed=T)==T))

# tableresults.zi1 <- dcast(data=tableresults.zi, species ~ variable, drop = FALSE)
# tableresults.zi1["signifsum",] <- apply(tableresults.zi1, 2, function(x) sum(grepl("*",x,fixed=T)==T))
# tableresults.cond1 <- dcast(data=tableresults.cond, species ~ variable, drop = FALSE)
# tableresults.cond1["signifsum",] <- apply(tableresults.cond1, 2, function(x) sum(grepl("*",x,fixed=T)==T))


#### Ajout fréquence et nombre de values negative or positive 
tableresults.cond[,"signifNEG"] <- apply(tableresults.cond, 1, function(x) sum(grepl("^-.*\\*$",x,fixed=F)==T))
tableresults.cond[,"signifPOS"] <- apply(tableresults.cond, 1, function(x) sum(grepl("^[0-9].*\\*$",x,fixed=F)==T))
tableresults.zi[,"signifPOS"] <- apply(tableresults.zi, 1, function(x) sum(grepl("^-.*\\*$",x,fixed=F)==T)) ### INVERSED 
tableresults.zi[,"signifNEG"] <- apply(tableresults.zi, 1, function(x) sum(grepl("^[0-9].*\\*$",x,fixed=F)==T)) ## IDEM inverse compared to conditional part


rownames(tableresults.cond) <-tableresults.cond$variable
rownames(tableresults.zi) <-tableresults.zi$variable

tableresults.cond[,"FREQsignif"] <- apply(tableresults.cond, 1, function(x) round(as.numeric(x["signifsum"])/length(which(!is.na(x[2:19]))),2))
tableresults.cond[,"FREQsignifPOS"] <- apply(tableresults.cond, 1, function(x) round(as.numeric(x["signifPOS"])/length(which(!is.na(x[2:19]))),2))
tableresults.cond[,"FREQsignifNEG"] <- apply(tableresults.cond, 1, function(x) round(as.numeric(x["signifNEG"])/length(which(!is.na(x[2:19]))),2))
tableresults.zi[,"FREQsignif"] <- apply(tableresults.zi, 1, function(x) round(as.numeric(x["signifsum"])/length(which(!is.na(x[2:19]))),2))
tableresults.zi[,"FREQsignifPOS"] <- apply(tableresults.zi, 1, function(x) round(as.numeric(x["signifPOS"])/length(which(!is.na(x[2:19]))),2))
tableresults.zi[,"FREQsignifNEG"] <- apply(tableresults.zi, 1, function(x) round(as.numeric(x["signifNEG"])/length(which(!is.na(x[2:19]))),2))

      

# Delete intercept ! 
tableresults.zi <- tableresults.zi[!tableresults.zi$variable=="(Intercept)",]
tableresults.cond <- tableresults.cond[!tableresults.cond$variable=="(Intercept)",]


### Here obtain barplot of frequency for CON models 
library(dplyr)
library(ggplot2)
library(cowplot)

COND <- melt(tableresults.cond, id=c("variable","FREQsignif"),measure.vars = c("FREQsignifPOS", "FREQsignifNEG"),variable.name = "Signe")
COND_sorted <- arrange(COND, variable, desc(Signe)) 
COND_cumsum <- ddply(COND_sorted, "variable",
                   transform, label_ypos=cumsum(value)-0.5*value)
p <- ggplot(data=COND_cumsum, aes(x=reorder(variable,-value), y=value,fill=Signe)) +
  theme(panel.background = element_rect(fill="white", colour="black", size=1,
                                        linetype=1, color="black")) +
  theme(panel.grid.major = element_line(size=0, linetype=0,lineend="square", color="red"),
        panel.grid.minor = element_line(size=0, linetype=0,lineend="round", color="blue"))+
  geom_bar(stat="identity",color="black",size=0.5)+
  geom_text(aes(y=label_ypos, label=ifelse(value==0,NA,value*100)),color="black", vjust=-0.01, size=4,fontface='bold')+
  labs(y=paste0("Significiance frequency"), x="Variable",title=paste0('b) Main effects and interaction importance across species (Conditionnal model)'))+
  scale_colour_manual(values = c("blue","red"),name="Effect sign",labels=c("Positive","Negative")) +
  scale_fill_manual(values = c("FREQsignifPOS"="steelblue","FREQsignifNEG"="indianred"),name="Effect sign",labels=c("Positive","Negative"))+
  scale_y_continuous(labels = scales::percent,limits = c(0,1), expand = c(0,0.01)) +
  theme(text = element_text(face="bold"),
        legend.direction ="vertical",
        legend.position = c(.95, .95),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.box.background = element_rect(color="black", size=0.5),
        legend.box.margin = margin(6, 6, 6, 6),
        axis.text.x = element_text(size=11,color="black",angle = 45,hjust = 1,vjust = 1),axis.text.y = element_text(size=18,color="black"),
        axis.title.x = element_text(size=18),axis.title.y = element_text(size=18),
        legend.key = element_rect(fill = "white", colour = "black",size = 1),
        legend.key.heigh = unit(3,"line"),
        legend.key.width = unit(1.5,"line"),
        legend.text=element_text(size=14),
        legend.title=element_text(size=16),
        #legend.justification = "center",
        legend.margin = margin(0,0,0,0),
        #legend.background=element_rect(fill="white",colour="black",size=0.5,linetype="solid"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.3),
        axis.line = element_line(colour="black"),
        plot.title = element_text(size=17,hjust = 0,vjust = 0),
        plot.margin = margin(0.2,0.2,0.2,0.2, "cm"),
        plot.caption = element_text(face="bold.italic"))
p
save_plot(filename = paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/paper2/Freq.Barplot.COND.pdf"),plot = p,base_width = 10, base_height = 7, dpi = 500 ,units = "in",nrow=1,ncol=1)



## Here the same barplots for Zero inflated part of the model ! 
ZI <- melt(tableresults.zi, id=c("variable","FREQsignif"),measure.vars = c("FREQsignifPOS", "FREQsignifNEG"),variable.name = "Signe")
ZI_sorted <- arrange(ZI, variable, desc(Signe)) 
ZI_cumsum <- ddply(ZI_sorted, "variable",
                     transform, label_ypos=cumsum(value)-0.5*value)
p1 <- ggplot(data=ZI_cumsum, aes(x=reorder(variable,-value), y=value,fill=Signe)) +
  theme(panel.background = element_rect(fill="white", colour="black", size=1,
                                        linetype=1, color="black")) +
  theme(panel.grid.major = element_line(size=0, linetype=0,lineend="square", color="red"),
        panel.grid.minor = element_line(size=0, linetype=0,lineend="round", color="blue"))+
  geom_bar(stat="identity",color="black",size=0.5)+
  geom_text(aes(y=label_ypos, label=ifelse(value==0,NA,value*100)),color="black", vjust=-0.01, size=4,fontface='bold')+
  labs(y=paste0("Significiance frequency"), x="Variable",title=paste0('a) Main effects and interaction importance across species (ZI model)'))+
  scale_colour_manual(values = c("blue","red"),name="Effect sign",labels=c("Positive","Negative")) +
  scale_fill_manual(values = c("FREQsignifPOS"="steelblue","FREQsignifNEG"="indianred"),name="Effect sign",labels=c("Positive","Negative"))+
  scale_y_continuous(labels = scales::percent,limits = c(0,1), expand = c(0,0.01)) +
  theme(text = element_text(face="bold"),
        legend.direction ="vertical",
        legend.position = c(.95, .95),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.box.background = element_rect(color="black", size=0.5),
        legend.box.margin = margin(6, 6, 6, 6),
        axis.text.x = element_text(size=11,color="black",angle = 45,hjust = 1,vjust = 1),axis.text.y = element_text(size=18,color="black"),
        axis.title.x = element_text(size=18),axis.title.y = element_text(size=18),
        legend.key = element_rect(fill = "white", colour = "black",size = 1),
        legend.key.heigh = unit(3,"line"),
        legend.key.width = unit(1.5,"line"),
        legend.text=element_text(size=14),
        legend.title=element_text(size=16),
        #legend.justification = "center",
        legend.margin = margin(0,0,0,0),
        #legend.background=element_rect(fill="white",colour="black",size=0.5,linetype="solid"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.3),
        axis.line = element_line(colour="black"),
        plot.title = element_text(size=17,hjust = 0,vjust = 0),
        plot.margin = margin(0.2,0.2,0.2,0.2, "cm"),
        plot.caption = element_text(face="bold.italic"))
p1
save_plot(filename = paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/paper2/Freq.Barplot.ZI.pdf"),plot = p1,base_width = 10, base_height = 7, dpi = 500 ,units = "in",nrow=1,ncol=1)


library(patchwork)
pall <- p1+theme(plot.margin = unit(c(0,0,0,0), "cm"))+
  p+theme(plot.margin = unit(c(0,0,0,0), "cm"))+
  plot_layout(ncol=1,nrow=2)

save_plot(filename = paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/paper2/Figure2.pdf"),
          plot = pall,base_width = 10, dpi = 1000, base_height = 7,units = "in",nrow=2,ncol=1)
# save_plot(filename = paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/Paper1/Figure3/Ag.Prediction.Mort.bin.Part1.pdf"),
#           plot = pall,base_width = 5, base_height = 5, dpi = 300 ,units = "in",nrow=3,ncol=3)









## Order column 
colnames(tableresults.cond)
rownames(tableresults.cond)
tableresults.zi$variable
tableresults.zi <- tableresults.zi[c(1,26,3,16,2,8,6,17,10,20,11,21,14,24,15,25,13,23,12,22,5,4,9,7,19,18),]
tableresults.cond <- tableresults.cond[c(1,26,3,16,2,8,6,17,10,20,11,21,14,24,15,25,13,23,12,22,5,4,9,7,19,18),]


tableresults.zi <- tableresults.zi[c(7,1,2,5,17,9,6,18,11,21,15,25,12,22,16,26,14,24,13,23,4,3,10,8,20,19),]
tableresults.cond <- tableresults.cond[c(7,1,2,5,17,9,6,18,11,21,15,25,12,22,16,26,14,24,13,23,4,3,10,8,20,19),]


write.csv2(tableresults.zi,file = paste0("~/SynologyDrive/FUNDIV - NFI - Europe/our-data/species/Recrut.All/ziMod_",Modname,".csv"))
write.csv2(tableresults.cond,file = paste0("~/SynologyDrive/FUNDIV - NFI - Europe/our-data/species/Recrut.All/condMod_",Modname,".csv"))




#### Obtain the figure we want (first we need to obtain the myDFsim ALL data for all species script 22B.intercation.recruitment.R)
### Then load the script 24
Dir <- c("~/SynologyDrive/FUNDIV - NFI - Europe/")
setwd(Dir)
source(paste0(Dir,"Myscripts/Fundiv.project/24B.Recruitment.fig.R"))


## Check Plotcat BAj (mu p and responses)
## Check Plotcat Mortality (mu p and responses)
# BAJ*BAI
# BAJ*dbh 
# BAO*dbh 
# BAO*BAI 

nameVAR
PredVAR

Interaction(Modname="rW.yO.sumO.mRY1",name3 = nameVAR[4], nameSeq = nameVAR[2],Pred = PredVAR[2],save = F)


Interaction(Modname="rW.yO.sumO.mRY1",name3 = nameVAR[4], nameSeq = nameVAR[1],Pred = PredVAR[3],save = T)
Interaction(Modname="rW.yO.sumO.mRY1",name3 = nameVAR[3], nameSeq = nameVAR[6],Pred = PredVAR[2],save = T)
Interaction(Modname="rW.yO.sumO.mRY1",name3 = nameVAR[8], nameSeq = nameVAR[7],Pred = PredVAR[3],save = T)
Interaction(Modname="rW.yO.sumO.mRY1",name3 = nameVAR[8], nameSeq = nameVAR[4],Pred = PredVAR[2],save = T)

Interaction(Modname="rW.yO.sumO.mRY1",name3 = nameVAR[4], nameSeq = nameVAR[1],Pred = PredVAR[2],save = T)
Interaction(Modname="rW.yO.sumO.mRY1",name3 = nameVAR[3], nameSeq = nameVAR[1],Pred = PredVAR[3],save = T)
Interaction(Modname="rW.yO.sumO.mRY1",name3 = nameVAR[8], nameSeq = nameVAR[7],Pred = PredVAR[2],save = T)

Interaction(Modname="rW.yO.sumO.mRY1",name3 = nameVAR[4], nameSeq = nameVAR[1],Pred = PredVAR[3],save = T)

Interaction(Modname="rW.yO.sumO.mRY1",name3 = nameVAR[8], nameSeq = nameVAR[1],Pred = PredVAR[3],save = T)
Interaction(Modname="rW.yO.sumO.mRY1",name3 = nameVAR[8], nameSeq = nameVAR[2],Pred = PredVAR[3],save = T)
Interaction(Modname="rW.yO.sumO.mRY1",name3 = nameVAR[8], nameSeq = nameVAR[4],Pred = PredVAR[3],save = T)



Interaction(name3 = nameVAR[8], nameSeq = nameVAR[7],Pred = PredVAR[2],save = T)
Interaction(name3 = nameVAR[8], nameSeq = nameVAR[7],Pred = PredVAR[3],save =T)

Interaction(name3 = nameVAR[1], nameSeq = nameVAR[4],Pred = PredVAR[2],save = T)
Interaction(name3 = nameVAR[1], nameSeq = nameVAR[4],Pred = PredVAR[3],save =T)
Interaction(name3 = nameVAR[4], nameSeq = nameVAR[1],Pred = PredVAR[2],save = T)
Interaction(name3 = nameVAR[4], nameSeq = nameVAR[1],Pred = PredVAR[3],save =T)


Interaction(name3 = nameVAR[6], nameSeq = nameVAR[4],Pred = PredVAR[2],save = T)
Interaction(name3 = nameVAR[6], nameSeq = nameVAR[4],Pred = PredVAR[3],save =T)
Interaction(name3 = nameVAR[4], nameSeq = nameVAR[6],Pred = PredVAR[2],save = T)
Interaction(name3 = nameVAR[4], nameSeq = nameVAR[6],Pred = PredVAR[3],save =T)


Interaction(name3 = nameVAR[1], nameSeq = nameVAR[3],Pred = PredVAR[2],save = T)
Interaction(name3 = nameVAR[1], nameSeq = nameVAR[3],Pred = PredVAR[3],save =T)
Interaction(name3 = nameVAR[3], nameSeq = nameVAR[1],Pred = PredVAR[2],save = T)
Interaction(name3 = nameVAR[3], nameSeq = nameVAR[1],Pred = PredVAR[3],save =T)


Interaction(name3 = nameVAR[6], nameSeq = nameVAR[3],Pred = PredVAR[2],save = T)
Interaction(name3 = nameVAR[6], nameSeq = nameVAR[3],Pred = PredVAR[3],save =T)
Interaction(name3 = nameVAR[3], nameSeq = nameVAR[6],Pred = PredVAR[2],save = T)
Interaction(name3 = nameVAR[3], nameSeq = nameVAR[6],Pred = PredVAR[3],save =T)

# And check SPEI in the case of. (SPEI*Plotcat and SPEI*BAj)
Interaction(name3 = nameVAR[2], nameSeq = nameVAR[4],Pred = PredVAR[2],save = T)
Interaction(name3 = nameVAR[4], nameSeq = nameVAR[2],Pred = PredVAR[3],save =T)

#### Premodel : Use script number 28. 
### VIF 
Dir <- c("/Users/alexandrechangenet/SynologyDrive/FUNDIV - NFI - Europe/")
setwd(Dir)
source(paste0(Dir,"Myscripts/Fundiv.project/27.VIF.Recrut.R"))
VIF.R(Modname)


###### Diagnostic of the model: (script 26)














