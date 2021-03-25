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
library(performance)
i <- 1
Allcode <- c("ABIALB","ACEPSE","ALNGLU","BETPEN","CASSAT","FAGSYL",
             "FRAEXC","PICABI","PINHAL","PINNIG","PINPIN","PINPINA","PINSYL","POPNIG",
             "POPTRE","QUEILE","QUEPET","QUEPYR","QUEROB","QUESUB")
Allseuil <- c(0.7,0.55,0.6,0.8,0.7,0.7,0.7,0.8,0.7,0.7,0.7,0.8,0.8,0.7,0.7,0.8,0.7,0.7,0.8,0.7)
AllseuilC <- c(0.6,0.5,0.5,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6)

NAME <- c("SCALEDRrate.yr.QualiF.TF.noY.nbinom1.V3","SCALEDRrate.yr.QualiF.TF.YO.nbinom1","SCALEDRrate.yr.QualiF.TF.YF.nbinom1","rW.yF.sumF.mC", "rW.yF.sumF.mR", "rW.yF.sumF.mRY", "rW.yF.sumO.mC", 
  "rW.yF.sumO.mR", "rW.yF.sumO.mRY", "rW.yO.sumF.mC", "rW.yO.sumF.mR", 
  "rW.yO.sumF.mRY", "rW.yO.sumO.mC", "rW.yO.sumO.mR", "rW.yO.sumO.mRY")


NAME <- c("SCALEDRrate.yr.QualiF.TF.noY.nbinom1.V3",
          "rW.yO.sumO.mRY",
          "rW.yO.sumO.mRY1",
          "rW.yO.sumO.mRY1.nT")

ModAIC <- list()
for (i in c(1:13,16:length(Allcode))){
  # for (i in c(1:3)){
  CODE <- Allcode[i]
  seuilC <- AllseuilC[i]
  seuil <- Allseuil[i]
  Dir = paste0("~/SynologyDrive/FUNDIV - NFI - Europe/our-data/species/",CODE,"/CLIMAP/Recrut.Mortality.2020/")
  setwd(Dir)
  # Charge the dfplot for simulation and residuals test dispersion etc... 
  dfplot <- readRDS(paste0("dfplotFinalV3",CODE,seuil,"R.M.rds"))
  ModAIC[[CODE]] <- as.data.frame(matrix(data=NA,nrow=length(NAME),ncol=14))
  colnames(ModAIC[[CODE]]) <- c("aic","cor","zi.mort.rate.yr","zi.mort.rate.yr.W","zi.mort.rate.W","zi.mort.W","ziyears","zisumWeight",
                                "cond.mort.rate.yr","cond.mort.rate.yr.W","cond.mort.rate.W","cond.mort.W","condyears","condsumWeight")
  rownames(ModAIC[[CODE]]) <- NAME
  
  for (j in c(1:length(NAME))){ 
  Modname <- NAME[j]
  #Load the desired model of the desired species 
  Dir <- c(paste0("~/SynologyDrive/FUNDIV - NFI - Europe/our-data/species/",CODE,"/CLIMAP/Recrut.Mortality.2020/",Modname,"/")) # Directory 
  setwd(Dir)
  mod_ALL <- try(get(load(file = paste0(Modname,".rda"))),silent=T)
  if (inherits(mod_ALL,"try-error")==T) {next}
  rm(x)
  
  ModAIC[[CODE]][Modname,"aic"] <- AIC(mod_ALL)
  ModAIC[[CODE]][Modname,"cor"] <- cor(mod_ALL$frame[,1],predict(mod_ALL,type = "response"))^2
  ModAIC[[CODE]][Modname,"zi.mort.rate.yr"] <- as.numeric(try(summary(mod_ALL)[["coefficients"]][["zi"]]["logsp.mortality.plot.rate.yr",4],silent=T))
  ModAIC[[CODE]][Modname,"zi.mort.rate.yr.W"] <- as.numeric(try(summary(mod_ALL)[["coefficients"]][["zi"]]["sp.mortality.ha.weight.rate.yr",4],silent=T))
  ModAIC[[CODE]][Modname,"zi.mort.rate.W"] <- as.numeric(try(summary(mod_ALL)[["coefficients"]][["zi"]]["sp.mortality.ha.weight.rate",4],silent=T))
  ModAIC[[CODE]][Modname,"zi.mort.W"] <- as.numeric(try(summary(mod_ALL)[["coefficients"]][["zi"]]["sp.mortality.ha.weight",4],silent=T))
  ModAIC[[CODE]][Modname,"ziyears"] <- as.numeric(try(summary(mod_ALL)[["coefficients"]][["zi"]]["yearsbetweensurveys",4],silent=T))
  ModAIC[[CODE]][Modname,"zisumWeight"] <- as.numeric(try(summary(mod_ALL)[["coefficients"]][["zi"]]["sp.SUM.ha.weight2",4],silent=T))
  ModAIC[[CODE]][Modname,"cond.mort.rate.yr"] <- as.numeric(try(summary(mod_ALL)[["coefficients"]][["cond"]]["logsp.mortality.plot.rate.yr",4],silent=T))
  ModAIC[[CODE]][Modname,"cond.mort.rate.yr.W"] <- as.numeric(try(summary(mod_ALL)[["coefficients"]][["cond"]]["sp.mortality.ha.weight.rate.yr",4],silent=T))
  ModAIC[[CODE]][Modname,"cond.mort.rate.W"] <- as.numeric(try(summary(mod_ALL)[["coefficients"]][["cond"]]["sp.mortality.ha.weight.rate",4],silent=T))
  ModAIC[[CODE]][Modname,"cond.mort.W"] <- as.numeric(try(summary(mod_ALL)[["coefficients"]][["cond"]]["sp.mortality.ha.weight",4],silent=T))
  ModAIC[[CODE]][Modname,"condyears"] <- as.numeric(try(summary(mod_ALL)[["coefficients"]][["cond"]]["yearsbetweensurveys",4],silent=T))
  ModAIC[[CODE]][Modname,"condsumWeight"] <- as.numeric(try(summary(mod_ALL)[["coefficients"]][["cond"]]["sp.SUM.ha.weight2",4],silent=T))
  rm(mod_ALL)
  ## effet 
  }
  ### Test #### 
  #we define the data as the part of the dataset on which the model has been fitted. (for every loop)
  #dfplot <- mod_ALL$frame
  #remove dfplot of the current species to ensure it won't be used for next species 
  rm(dfplot) 
  rm(mod_ALL)
}

saveRDS(ModAIC,file="~/SynologyDrive/FUNDIV - NFI - Europe/our-data/species/Recrut.All/MultiModAIC2.rds")
A <- readRDS(file="~/SynologyDrive/FUNDIV - NFI - Europe/our-data/species/Recrut.All/MultiModAIC.rds")
