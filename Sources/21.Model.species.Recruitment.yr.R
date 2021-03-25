rm(list = ls())
gc()
library(parallel)
library(forecast)
library(glmmTMB)
library(bbmle)

eval_fork <- function(..., timeout=600){
  starttime <- Sys.time(); # Depart 
  #this limit must always be higher than the timeout on the fork!
  setTimeLimit(timeout+5);      
  myfork <- parallel::mcparallel({
    eval(...)
  }, silent=FALSE);
  #wait max n seconds for a result.
  myresult <- parallel::mccollect(myfork, wait=FALSE, timeout=timeout);
  enddtime <- Sys.time(); # Ajout 
  totaltime <- as.numeric(enddtime - starttime, units="secs")
  
  #waits for max another 2 seconds if proc looks dead 
  while(is.null(myresult) && totaltime < timeout && totaltime < 4) {
    Sys.sleep(.1)
    enddtime <- Sys.time();
    totaltime <- as.numeric(enddtime - starttime, units="secs")
    myresult <- parallel::mccollect(myfork, wait = FALSE, timeout = timeout);
  }
  #kill fork after collect has returned
  tools::pskill(myfork$pid, tools::SIGKILL);    
  tools::pskill(-1 * myfork$pid, tools::SIGKILL);  
  #clean up:
  parallel::mccollect(myfork, wait=FALSE);
  if(is.null(myresult)){
    message("R call did not return within ", timeout, " seconds. Terminating process. => ",Model," and ",CODE," , Next iteration")
    stop("R call did not return within ", timeout, " seconds. Terminating process.", call.=FALSE);      
  }
  
  #move this to distinguish between timeout and NULL returns
  myresult <- myresult[[1]];
  #reset timer
  setTimeLimit();     
  #forks don't throw errors themselves
  
  if(inherits(myresult,"try-error")){
    message("R call did not return within ", timeout, " seconds. Terminating process. => ",Model," and ",CODE," , Next iteration")
    stop(attr(myresult, "condition"));
  } else return(myresult) #send the buffered response if not an error
}

Dir <- c("/Users/alexandrechangenet/SynologyDrive/FUNDIV - NFI - Europe/")
setwd(Dir)
source(paste0(Dir,"Myscripts/Fundiv.project/RecrutProject.function1.Saving.R"))
Allcode <- c("PINPINA","ABIALB","ACEPSE","ALNGLU","BETPEN","CASSAT","FAGSYL","FRAEXC","PICABI","PINHAL","PINNIG","PINPIN","PINSYL","POPNIG","POPTRE","QUEILE","QUEPET",
             "QUEPYR","QUEROB","QUESUB")
Allseuil <- c(0.8,0.7,0.55,0.6,0.8,0.7,0.7,0.7,0.8,0.7,0.7,0.7,0.8,0.7,0.7,0.8,0.7,0.7,0.8,0.7)
AllseuilC <- c(0.6,0.6,0.5,0.5,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6)

# Errors.files <- file(paste0(Dir,"our-data/species/Errors.allmodels.Recrut.Rout"), open="wt")
# sink(Errors.files, type="message")
# do the same with 0.7

for (i in c(5,8,11,17)){
# for (i in c(1:13,16:length(Allcode))){
  CODE <- Allcode[i]
  seuilC <- AllseuilC[i]
  seuil <- Allseuil[i]
  Dir = (paste0("/Users/alexandrechangenet/SynologyDrive/FUNDIV - NFI - Europe/our-data/species/",CODE,"/Recrut.Mortality.2020/")) # Directory 
  setwd(Dir)
  dfplot <- try(readRDS(paste0("dfplotFinalV2",CODE,seuil,"R.M.rds")),silent = T) #Base de données plot

## Add interaction with competition 
  rm(SCALEDRrate.yr.QualiF.TF.noY.nbinom1.V3)
  try(SCALEDRrate.yr.QualiF.TF.noY.nbinom1.V3 <- glmmTMB(sp.recruitment.plot.rate.yr.2 ~ logdbhJ.IMall.plot.mean + (1|country) + logtreeNbrJ + 
                                                             Plotcat*mean_spei12 + Plotcat*sqrtBA.O.plot.1 + Plotcat*logBAj.plot.1 + Plotcat*sqrtBAIj.plot.1.mean.J.I + logsp.mortality.plot.rate.yr +
                                                           mean_spei12*logBAj.plot.1 + logdbhJ.IMall.plot.mean*logBAj.plot.1 + sqrtBAIj.plot.1.mean.J.I*logBAj.plot.1 + 
                                                           mean_spei12*sqrtBA.O.plot.1 + logdbhJ.IMall.plot.mean*sqrtBA.O.plot.1 + sqrtBAIj.plot.1.mean.J.I*sqrtBA.O.plot.1,
                                                           zi=~.,family=nbinom1, data=dfplot),silent=T)
  try(Saving(SCALEDRrate.yr.QualiF.TF.noY.nbinom1.V3),silent=T)
  try(print(summary(SCALEDRrate.yr.QualiF.TF.noY.nbinom1.V3)),silent=T)
  if (length(summary(as.factor(SCALEDRrate.yr.QualiF.TF.noY.nbinom1.V3$frame$country)))==1){
    rm(SCALEDRrate.yr.QualiF.TF.noY.nbinom1.V3)
    try(SCALEDRrate.yr.QualiF.TF.noY.nbinom1.V3 <- glmmTMB(sp.recruitment.plot.rate.yr.2 ~ logdbhJ.IMall.plot.mean + logtreeNbrJ + 
                                                             Plotcat*mean_spei12 + Plotcat*sqrtBA.O.plot.1 + Plotcat*logBAj.plot.1 + Plotcat*sqrtBAIj.plot.1.mean.J.I + logsp.mortality.plot.rate.yr +
                                                             mean_spei12*logBAj.plot.1 + logdbhJ.IMall.plot.mean*logBAj.plot.1 + sqrtBAIj.plot.1.mean.J.I*logBAj.plot.1 + 
                                                             mean_spei12*sqrtBA.O.plot.1 + logdbhJ.IMall.plot.mean*sqrtBA.O.plot.1 + sqrtBAIj.plot.1.mean.J.I*sqrtBA.O.plot.1,
                                                           zi=~.,family=nbinom1, data=dfplot),silent=T)
    try(Saving(SCALEDRrate.yr.QualiF.TF.noY.nbinom1.V3),silent=T)
    try(print(summary(SCALEDRrate.yr.QualiF.TF.noY.nbinom1.V3)),silent=T)
    
  }
  
}

## Parallelization
glmmTMBControl(parallel = 10)
##

myDF$sp.recruitment.plot.rate.2
mHbinom2$frame$sp.recruitment.plot.rate.2
mHbinom2

summary(mHbinom1)
summary(mbinom1)
  
  

?write.table        
        
        
### Same with annual rate 
### Same by removing 




dfplot$sp.recruitment.plot.rate
## Verif ###
by(dfplot2$sp.mortality.plot.rate,dfplot2$country,summary)
boxplot(dfplot2$sp.mortality.plot.rate~dfplot2$country)
by(dfplot2$sp.mortality.plot.rate.yr,dfplot2$country,summary)
boxplot(dfplot2$sp.mortality.plot.rate.yr~dfplot2$country)

by(dfplot2$sp.mortality.plot.rate,dfplot2$yearsbetweensurveys,summary)
by(dfplot2$sp.mortality.plot.rate.yr,dfplot2$yearsbetweensurveys,summary)


by(dfplot2$BAIj.plot.bis,dfplot2$country,summary)
boxplot(dfplot2$BAIj.plot.bis~dfplot2$country)

by(dfplot2$BAIj.plot.1,dfplot2$country,summary)
boxplot(dfplot2$BAIj.plot.1~dfplot2$country)

by(dfplot2$BAIj.plot,dfplot2$country,summary)
boxplot(dfplot2$BAIj.plot~dfplot2$country)

by(dfplot2$BAIj.plot.mean,dfplot2$country,summary)
boxplot(dfplot2$BAIj.plot.mean~dfplot2$country)

by(dfplot2$BAIj.plot.1.mean,dfplot2$country,summary)
boxplot(dfplot2$BAIj.plot.1.mean~dfplot2$country)

hist(dfplot2$sp.recruitment.plot.rate)
colnames(dfplot2)
Mbin1A <- glmer(sp.recruitment.plot.rate ~ treeNbr + yearsbetweensurveys + 
                  sqrtBAIj.plot.1.mean + logdbh.plot.mean + mean_spei12 + logBAj.plot.1 + sqrtBA.O.plot.1 + Plotcat + 
                  mean_spei12:logBAj.plot.1 + mean_spei12:sqrtBA.O.plot.1 + mean_spei12:Plotcat + logBAj.plot.1:Plotcat + sqrtBA.O.plot.1:Plotcat + (1|country), 
                data=dfplot2, family = poisson)

hist(dfplot2$sp.recruitment.plot.rate)
hist(dfplot2$sp.mortality.plot.rate)
hist(dfplot2$sp.mortality.plot.rate.yr)
boxplot(dfplot2$sp.mortality.plot.rate.yr)
summary(dfplot2$sp.mortality.plot.rate.yr)
summary(dfplot2$sp.mortality.plot.rate)
summary(dfplot2$sp.recruitment.plot.rate)
test <- dfplot2[is.na(dfplot2$sp.recruitment.plot.rate),]
test
dfplot2[1:50,]
colnames(dfplot2)
dfplot2$sp.recruitment.plot.rate.yr <-round((dfplot2$sp.recruitment.plot.rate * 1000 / dfplot2$yearsbetweensurveys))
hist(dfplot2$sp.recruitment.plot.rate.yr)
boxplot(dfplot2$sp.recruitment.plot.rate.yr)
summary(dfplot2$sp.recruitment.plot.rate.yr)
