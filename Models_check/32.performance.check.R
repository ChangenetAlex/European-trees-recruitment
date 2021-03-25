rm(list = ls())
gc()
library(glmmTMB)
library(car)
library(emmeans)
library(effects)
library(multcomp)
library(MuMIn)
library(DHARMa)
library(broom)
library(broom.mixed)
library(dotwhisker)
library(ggplot2); theme_set(theme_bw())
library(texreg)
library(xtable)
library(huxtable)
library(plyr)
library(beanplot)
library(parallel)

Modname <- "SCALEDRrate.yr.QualiF.TF.noY.nbinom1.V3" ## Name of the model 
Modname <- "rW.yO.sumO.mRY1"

Allcode <- c("ABIALB","ACEPSE","ALNGLU","BETPEN","CASSAT","FAGSYL",
             "FRAEXC","PICABI","PINHAL","PINNIG","PINPIN","PINPINA","PINSYL","POPNIG",
             "POPTRE","QUEILE","QUEPET","QUEPYR","QUEROB","QUESUB")
Allseuil <- c(0.7,0.55,0.6,0.8,0.7,0.7,0.7,0.8,0.7,0.7,0.7,0.8,0.8,0.7,0.7,0.8,0.7,0.7,0.8,0.7)
AllseuilC <- c(0.6,0.5,0.5,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6)
mod_ALL <- list()
A <- list()
dfplot <- list()
simulationOutput <- list()
sizeN = 1800
resol = 15
for (i in c(1:13,16:length(Allcode))){
  # for (i in c(1:3)){
  CODE <- Allcode[i]
  seuilC <- AllseuilC[i]
  seuil <- Allseuil[i]
  Dir = paste0("~/SynologyDrive/FUNDIV - NFI - Europe/our-data/species/",CODE,"/CLIMAP/Recrut.Mortality.2020/")
  setwd(Dir)
  # Charge the dfplot for simulation and residuals test dispersion etc... 
  dfplot[[CODE]] <- readRDS(paste0("dfplotFinalV3",CODE,seuil,"R.M.rds"))
  
  #Load the desired model of the desired species 
  Dir <- c(paste0(Dir,Modname,"/")) # Directory 
  setwd(Dir)
  mod_ALL[[CODE]] <- get(load(file = paste0(Modname,".rda")))
  rm(x)
  
  A[[CODE]] <- mod_ALL[[CODE]][[5]] ## Keep just the dataframe
  colnames(A[[CODE]])[which(colnames(A[[CODE]])=="offset(log(yearsbetweensurveys))")] <- "yearsbetweensurveys"
  colnames(A[[CODE]])[which(colnames(A[[CODE]])=="offset(log(sp.SUM.ha.weight2))")] <- "sp.SUM.ha.weight2"
  A[[CODE]][,"yearsbetweensurveys"] <- exp(A[[CODE]][,"yearsbetweensurveys"])
  A[[CODE]][,"sp.SUM.ha.weight2"] <- exp(A[[CODE]][,"sp.SUM.ha.weight2"])
  # 
  mod_ALL[[CODE]][["call"]][["data"]] <- substitute(A[[CODE]])
  #mod_ALL[[CODE]][["call"]][["data"]] <- substitute(dfplot[[CODE]])
  
  ### Test #### 
  #we define the data as the part of the dataset on which the model has been fitted. (for every loop)
  #dfplot <- mod_ALL[[CODE]]$frame
  ## Test ###
  
  #Simulate the residuals based on the right dfplot base (above) and the right species
  #simulationOutput[[CODE]] <- simulateResiduals(mod_ALL[[CODE]], n = 3000,re.form=NA) 
  png(file=paste0("~/SynologyDrive/FUNDIV - NFI - Europe/our-data/species/Recrut.All/",Modname,"/diagnostic.test/",CODE,"_Checkmodel_all.png"),width=sizeN,height=sizeN*0.9,res=sizeN/resol) # Save it 
  p1 <- check_model(mod_ALL[[CODE]], dot_size = 2, line_size = 0.8, panel = TRUE, check = c("outliers","normality","homogeneity","qq","reqq","ncv")) #all in one + the plot !
  print(p1)
  graphics.off()
  
  print("OK")
  
  png(file=paste0("~/SynologyDrive/FUNDIV - NFI - Europe/our-data/species/Recrut.All/",Modname,"/diagnostic.test/",CODE,"_Checkmodel_VIF.png"),width=sizeN,height=sizeN*0.9,res=sizeN/resol) # Save it 
  p2 <- check_model(mod_ALL[[CODE]], dot_size = 2, line_size = 0.8, panel = TRUE, check = c("vif")) #all in one + the plot !
  print(p2)
  graphics.off()
  
  print("OK both")
  
}


