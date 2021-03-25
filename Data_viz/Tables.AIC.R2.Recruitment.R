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

Modname <- "rW.yO.sumO.mRY1" ## Name of the model 
Modname <- "SCALEDRrate.yr.QualiF.TF.noY.nbinom1.V3"
Allcode <- c("ABIALB","ACEPSE","ALNGLU","BETPEN","CASSAT","FAGSYL",
             "FRAEXC","PICABI","PINHAL","PINNIG","PINPIN","PINPINA","PINSYL","POPNIG",
             "POPTRE","QUEILE","QUEPET","QUEPYR","QUEROB","QUESUB")
Allseuil <- c(0.7,0.55,0.6,0.8,0.7,0.7,0.7,0.8,0.7,0.7,0.7,0.8,0.8,0.7,0.7,0.8,0.7,0.7,0.8,0.7)
AllseuilC <- c(0.6,0.5,0.5,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6)
ModAIC <- list()
RMSE <- list()
for (i in c(1:13,16:length(Allcode))){
  # for (i in c(1:3)){
  CODE <- Allcode[i]
  seuilC <- AllseuilC[i]
  seuil <- Allseuil[i]
  Dir = paste0("~/SynologyDrive/FUNDIV - NFI - Europe/our-data/species/",CODE,"/CLIMAP/Recrut.Mortality.2020/")
  setwd(Dir)
  # Charge the dfplot for simulation and residuals test dispersion etc... 
  #dfplot <- readRDS(paste0("dfplotFinalV2",CODE,seuil,"R.M.rds"))
  #dfplot <- readRDS(paste0("dfplotFinalV3",CODE,seuil,"R.M.rds"))
  
  #Load the desired model of the desired species 
  Dir <- c(paste0(Dir,Modname,"/")) # Directory 
  setwd(Dir)
  mod_ALL <- get(load(file = paste0(Modname,".rda")))
  rm(x)
  dfplot <- mod_ALL$frame
  colnames(dfplot)[which(colnames(dfplot)=="offset(log(yearsbetweensurveys))")] <- "yearsbetweensurveys"
  colnames(dfplot)[which(colnames(dfplot)=="offset(log(sp.SUM.ha.weight2))")] <- "sp.SUM.ha.weight2"
  dfplot[,"yearsbetweensurveys"] <- exp(dfplot[,"yearsbetweensurveys"])
  dfplot[,"sp.SUM.ha.weight2"] <- exp(dfplot[,"sp.SUM.ha.weight2"])
  
  #colnames(dfplot)[3:4] <- c("yearsbetweensurveys","sp.SUM.ha.weight2") 
  RMSE[[CODE]] <- performance_rmse(mod_ALL, normalized = T)
  ModAIC[[CODE]] <- performance::performance(mod_ALL)
  ModAIC[[CODE]][1,"cor"] <- cor(mod_ALL$frame[,1],predict(mod_ALL,type = "response"))^2
  
  ### Test #### 
  #we define the data as the part of the dataset on which the model has been fitted. (for every loop)
  #remove dfplot of the current species to ensure it won't be used for next species 
  rm(dfplot) 
  rm(mod_ALL)
}


ModAIC
# ModAIC <- mapply(function(x,y) {
#   x[,"species"] <- y
#   return(x)},
#   x=ModAIC,y=names(ModAIC),SIMPLIFY = F)

ModAIC2 <- lapply(ModAIC, function(x) {x <- x[,c("AIC", "BIC", "R2_conditional", "R2_marginal", "RMSE", "cor")]})
ModAIC2 <- do.call(rbind,ModAIC2)
ModAIC2 <- round(as.data.frame(ModAIC2),2)
RMSE <- do.call(rbind,RMSE)
RMSE <- round(as.data.frame(RMSE),2)
ModAIC2 <- ModAIC2[,-c(5)]
ModAIC3 <- cbind(ModAIC2,RMSE)
write.table(ModAIC3,
           file = paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/paper2/AIC.R2.",Modname,".csv"),
           sep = ",",
           dec = ".")




