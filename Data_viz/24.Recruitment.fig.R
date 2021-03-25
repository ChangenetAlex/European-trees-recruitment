rm(list = ls())
gc()
library(parallel)
library(forecast)
library(glmmTMB)
library(bbmle)
library(ggplot2)
library(reshape2)
library(stringr)

# Dir <- c("/home/achangenet/SynologyDrive/FUNDIV - NFI - Europe/")
Dir <- c("/Users/alexandrechangenet/SynologyDrive/FUNDIV - NFI - Europe/")
setwd(Dir)
#source(paste0(Dir,"Myscripts/Fundiv.project/RecrutProject.function1.Saving.R"))
Allcode <- c("PINPINA","ABIALB","ACEPSE","ALNGLU","BETPEN","CASSAT","FAGSYL",
             "FRAEXC","PICABI","PINHAL","PINNIG","PINPIN","PINSYL","POPNIG",
             "POPTRE","QUEILE","QUEPET","QUEPYR","QUEROB","QUESUB")
Allseuil <- c(0.8,0.7,0.55,0.6,0.8,0.7,0.7,0.7,0.8,0.7,0.7,0.7,0.8,0.7,0.7,0.8,0.7,0.7,0.8,0.7)
AllseuilC <- c(0.6,0.6,0.5,0.5,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6)


####################################################
Lvl <- 30                                         ## Number of values for the interaction graphs
Modname <- "SCALEDRrate.yr.QualiF.TF.noY.nbinom1.V2" ## Name of the model 
nameVAR <- c(                                     ##
  "sqrtBAIj.plot.1.mean.J.I",                     ##
  "mean_spei12",                                  ##
  "sqrtBA.O.plot.1",                              ##
  "logBAj.plot.1",                                ##
  "logtreeNbrJ",                                     ## treeNbrJ
  "logdbhJ.IMall.plot.mean",                      ##
  "logsp.mortality.plot.rate.yr")                ##
PredVAR <- c(                                     ## Predictions to be on what part of the model ? 
  "Predictedresp",                                ##
  "Predictedzprob",                               ##
  "Predictedmu")                                  ##
Mylabels <- c(                                    ##
  "Core","Trailing Edge","Leading Edge")          ## Normale 
Mycol <- c("black", "red", "blue")                ## 1 trailing and 2 leading 
####################################################  

for (j in c(1:13,16:length(Allcode))){
CODE <- Allcode[j]
seuilC <- AllseuilC[j]
seuil <- Allseuil[j]
# Dir = paste0("/home/achangenet/SynologyDrive/FUNDIV - NFI - Europe/our-data/species/",CODE,"/CLIMAP/Recrut.Mortality.2020/")
Dir = paste0("/Users/alexandrechangenet/SynologyDrive/FUNDIV - NFI - Europe/our-data/species/",CODE,"/Recrut.Mortality.2020/")

setwd(Dir)
for (i in 1:length(nameVAR)){
  name <- nameVAR[i]
# try(assign(paste0("myDF",CODE,name),readRDS(paste0("/home/achangenet/SynologyDrive/FUNDIV - NFI - Europe/our-data/species/Recrut.All/",Modname,"/myDF.SIM_30_", CODE,".rds"))[((Lvl*3)*(i-1)+1):(Lvl*3*i),c("Plotcat",name,"Predictedresp","SEresp","Predictedmu","SEmu","Predictedzprob","SEzprob")]),silent=T)
try(assign(paste0("myDF",CODE,name),readRDS(paste0("/Users/alexandrechangenet/SynologyDrive/FUNDIV - NFI - Europe/our-data/species/Recrut.All/",Modname,"/myDF.SIM_30_", CODE,".rds"))[((Lvl*3)*(i-1)+1):(Lvl*3*i),c("Plotcat",name,"Predictedresp","SEresp","Predictedmu","SEmu","Predictedzprob","SEzprob")]),silent=T)
}}

for (i in 1:length(nameVAR)){
  name <- nameVAR[i]
assign(paste0("myDF_",name),do.call("rbind",sapply(ls(pattern = name), function(x) cbind(get(x)), simplify = FALSE)))
assign(paste0("species"),str_sub(str_extract(rownames(get(paste0("myDF_",name))), "DF.(.*?)[lmst]"),3,-2L))
assign(paste0("myDF_",name),cbind(get(paste0("myDF_",name)),species))
}

## If the species is one of the following we need to exchange the values for plotcat 1 and plotcat 2

REV <- c("PINPIN","PINNIG","QUEPET","QUEROB","PINPINA","POPTRE","ALNGLU","QUEPYR","QUESUB")
myDF_logBAj.plot.1$Plotcat2 <- myDF_logBAj.plot.1$Plotcat
myDF_logBAj.plot.1[myDF_logBAj.plot.1$species%in%REV & myDF_logBAj.plot.1$Plotcat2==1,"Plotcat"] <- 2
myDF_logBAj.plot.1[myDF_logBAj.plot.1$species%in%REV & myDF_logBAj.plot.1$Plotcat2==2,"Plotcat"] <- 1

myDF_sqrtBAIj.plot.1.mean.J.I$Plotcat2 <- myDF_sqrtBAIj.plot.1.mean.J.I$Plotcat
myDF_sqrtBAIj.plot.1.mean.J.I[myDF_sqrtBAIj.plot.1.mean.J.I$species%in%REV & myDF_sqrtBAIj.plot.1.mean.J.I$Plotcat2==1,"Plotcat"] <- 2
myDF_sqrtBAIj.plot.1.mean.J.I[myDF_sqrtBAIj.plot.1.mean.J.I$species%in%REV & myDF_sqrtBAIj.plot.1.mean.J.I$Plotcat2==2,"Plotcat"] <- 1

myDF_mean_spei12$Plotcat2 <- myDF_mean_spei12$Plotcat
myDF_mean_spei12[myDF_mean_spei12$species%in%REV & myDF_mean_spei12$Plotcat2==1,"Plotcat"] <- 2
myDF_mean_spei12[myDF_mean_spei12$species%in%REV & myDF_mean_spei12$Plotcat2==2,"Plotcat"] <- 1

myDF_logtreeNbrJ$Plotcat2 <- myDF_logtreeNbrJ$Plotcat
myDF_logtreeNbrJ[myDF_logtreeNbrJ$species%in%REV & myDF_logtreeNbrJ$Plotcat2==1,"Plotcat"] <- 2
myDF_logtreeNbrJ[myDF_logtreeNbrJ$species%in%REV & myDF_logtreeNbrJ$Plotcat2==2,"Plotcat"] <- 1

myDF_sqrtBA.O.plot.1$Plotcat2 <- myDF_sqrtBA.O.plot.1$Plotcat
myDF_sqrtBA.O.plot.1[myDF_sqrtBA.O.plot.1$species%in%REV & myDF_sqrtBA.O.plot.1$Plotcat2==1,"Plotcat"] <- 2
myDF_sqrtBA.O.plot.1[myDF_sqrtBA.O.plot.1$species%in%REV & myDF_sqrtBA.O.plot.1$Plotcat2==2,"Plotcat"] <- 1

myDF_logdbhJ.IMall.plot.mean$Plotcat2 <- myDF_logdbhJ.IMall.plot.mean$Plotcat
myDF_logdbhJ.IMall.plot.mean[myDF_logdbhJ.IMall.plot.mean$species%in%REV & myDF_logdbhJ.IMall.plot.mean$Plotcat2==1,"Plotcat"] <- 2
myDF_logdbhJ.IMall.plot.mean[myDF_logdbhJ.IMall.plot.mean$species%in%REV & myDF_logdbhJ.IMall.plot.mean$Plotcat2==2,"Plotcat"] <- 1

myDF_logsp.mortality.plot.rate.yr$Plotcat2 <- myDF_logsp.mortality.plot.rate.yr$Plotcat
myDF_logsp.mortality.plot.rate.yr[myDF_logsp.mortality.plot.rate.yr$species%in%REV & myDF_logsp.mortality.plot.rate.yr$Plotcat2==1,"Plotcat"] <- 2
myDF_logsp.mortality.plot.rate.yr[myDF_logsp.mortality.plot.rate.yr$species%in%REV & myDF_logsp.mortality.plot.rate.yr$Plotcat2==2,"Plotcat"] <- 1

myDF_sp.recruitment.plot.rate.yr.2$Plotcat2 <- myDF_sp.recruitment.plot.rate.yr.2$Plotcat
myDF_sp.recruitment.plot.rate.yr.2[myDF_sp.recruitment.plot.rate.yr.2$species%in%REV & myDF_sp.recruitment.plot.rate.yr.2$Plotcat2==1,"Plotcat"] <- 2
myDF_sp.recruitment.plot.rate.yr.2[myDF_sp.recruitment.plot.rate.yr.2$species%in%REV & myDF_sp.recruitment.plot.rate.yr.2$Plotcat2==2,"Plotcat"] <- 1


## Plot 

nameVAR
PredVAR
name <- nameVAR[5]
Pred <- PredVAR[3] # For predvar 2 -> need to use 1-get(pred)

for (i in 1:length(nameVAR)){
name <- nameVAR[i]
p1<-ggplot(get(paste0("myDF_",name)), # All
           aes(x=get(name), y=get(Pred), color=Plotcat,shape=Plotcat,group=Plotcat))+
  geom_line(size=1)+ #no size
  geom_point(size=1.5,stroke=1.5)+
  facet_wrap(vars(species),scales = "free")+
  #geom_linerange(size=1.5)+ # no size
  scale_alpha_manual(values = c(0.3,1),guide=F)+
  scale_color_manual("",values = Mycol,labels=Mylabels)+
  scale_shape_manual("", values=c(1,2,3),labels=Mylabels)+
  labs(y=Pred, x=name)+
  theme_light(base_size = 15)+
  theme(text = element_text(face="bold"),legend.position = c(0.5,0.95),legend.direction ="horizontal",
        axis.text.x = element_text(size=13,color="black"),axis.text.y = element_text(size=13,color="black"),
        legend.background=element_rect(fill="white",colour="black",size=0.2),
        panel.border = element_rect(colour = "black", fill=NA, size=0.8),
        axis.line = element_line(colour="black"),
        plot.title = element_text(size=18,hjust = 0.5),
        plot.caption = element_text(face="bold.italic"))
# ggsave(filename = paste0("/home/achangenet/SynologyDrive/FUNDIV - NFI - Europe/our-data/species/Recrut.All/",Pred,"_",name,"_",Modname,".png"),plot = p1, width = 6, height = 6, dpi=300) # To modify 
# ggsave(filename = paste0("/home/achangenet/SynologyDrive/FUNDIV - NFI - Europe/our-data/species/Recrut.All/",Pred,"_",name,"_",Modname,".png"),plot = p1, scale=0.8, width = 20,height = 16, dpi=400) # To modify 
ggsave(filename = paste0("/Users/alexandrechangenet/SynologyDrive/FUNDIV - NFI - Europe/our-data/species/Recrut.All/",Pred,"_",name,"_",Modname,".png"),plot = p1, scale=0.8, width = 20,height = 16, dpi=400) # To modify 
print(p1)

}

### Ajout partie grisé pour chaque espèces + regarder les trois réponses !!! (ou seulement 2)








