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

Dir <- c("~/SynologyDrive/FUNDIV - NFI - Europe/")
setwd(Dir)
source(paste0(Dir,"Myscripts/Fundiv.project/RecrutProject.function1.Saving.R"))
Allcode <- c("PINPINA","ABIALB","ACEPSE","ALNGLU","BETPEN","CASSAT","FAGSYL","FRAEXC","PICABI","PINHAL","PINNIG","PINPIN","PINSYL","POPNIG","POPTRE","QUEILE","QUEPET",
             "QUEPYR","QUEROB","QUESUB")
Allseuil <- c(0.8,0.7,0.55,0.6,0.8,0.7,0.7,0.7,0.8,0.7,0.7,0.7,0.8,0.7,0.7,0.8,0.7,0.7,0.8,0.7)
AllseuilC <- c(0.6,0.6,0.5,0.5,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6)

DFPLOT <- list()
for (i in c(1:length(Allcode))){
  CODE <- Allcode[i]
  seuilC <- AllseuilC[i]
  seuil <- Allseuil[i]
  Dir = paste0("~/SynologyDrive/FUNDIV - NFI - Europe/our-data/species/",CODE,"/CLIMAP/Recrut.Mortality.2020/")
  setwd(Dir)
  DFPLOT[[CODE]] <- try(readRDS(paste0("dfplotFinalV3",CODE,seuil,"R.M.rds")),silent = T) #Base de données plot
}
rm("CODE") 

mcmapply(function(dfplot,CODE){
  # CODE <- Allcode[i]
  # seuilC <- AllseuilC[i]
  # seuil <- Allseuil[i]
  Dir = (paste0("~/SynologyDrive/FUNDIV - NFI - Europe/our-data/species/",CODE,"/CLIMAP/Recrut.Mortality.2020/")) # Directory 
  setwd(Dir)
  #dfplot <- try(readRDS(paste0("dfplotFinalV2",CODE,seuil,"R.M.rds")),silent = T) #Base de données plot
  #dfplot <- try(readRDS(paste0("dfplotFinalV3",CODE,seuil,"R.M.rds")),silent = T) #Base de données plot
  
  # Add recruitment as a count per year and as a total count
  dfplot[,"sp.recruitment.plot.rate.yr.2"] <- round((dfplot[,"sp.recruitment.plot.rate"]*1000/dfplot[,"yearsbetweensurveys"])) # taux de recruitment.ha/year = response 1 first model serie
  dfplot[,"sp.recruitment.plot.rate.2"] <- round(dfplot[,"sp.recruitment.plot.rate"] * 1000)
  
  # Add recruitment type 2
  dfplot[,"sp.recruitment.ha.weight.2"] <- round(dfplot[,"sp.recruitment.ha.weight"]) # mortality count / ha (weight). NON ramené au plot mais ramené à l'hectare = SOMME
  dfplot[,"sp.recruitment.ha.weight.rate"] <- round(dfplot[,"sp.recruitment.ha.weight"]/dfplot[,"sp.SUM.ha.weight2"]*1000) # mortality count / ha (weight) ramené au plot et a l'hectare
  dfplot[,"sp.recruitment.ha.weight.rate.yr"] <- round((dfplot[,"sp.recruitment.ha.weight"]/dfplot[,"sp.SUM.ha.weight2"]*1000)/dfplot[,"yearsbetweensurveys"]) # ramené au plot a l'hectare et à l'année
  
  # Recruitment type 3
  #dfplot[,"treeNbrJ.R"] # number of event on the plot only ! 
  dfplot[,"sp.recruitment.plot.number.rate"] <- round(dfplot[,"treeNbrJ.R"]/dfplot[,"treeNbrJ.IR"]*1000) # number of event ramené à la taille du plot (nbr arbre)
  dfplot[,"sp.recruitment.plot.number.rate.yr"] <- round((dfplot[,"treeNbrJ.R"]/dfplot[,"treeNbrJ.IR"]*1000)/dfplot[,"yearsbetweensurveys"]) # ramené aussi à l'année 
  
  
  # Mortality rate type 1 => 
  #dfplot[,"sp.mortality.plot.rate"]     # mortality rate / ha (weight) = ramené au poids total => To compare with weight dead/weight total and it should be the same !
  #dfplot[,"sp.mortality.plot.rate.yr"]  # mortality rate / ha (weight)/year (ramené à l'année)
  
  # Mortality type 2
  #dfplot[,"sp.mortality.ha.weight"] # mortality count / ha (weight). NON ramené au plot mais ramené à l'hectare = SOMME
  dfplot[,"sp.mortality.ha.weight.rate"] <- dfplot[,"sp.mortality.ha.weight"]/dfplot[,"sp.SUM.ha.weight1"] # mortality count / ha (weight) ramené au plot et a l'hectare
  dfplot[,"sp.mortality.ha.weight.rate.yr"] <- dfplot[,"sp.mortality.ha.weight.rate"]/dfplot[,"yearsbetweensurveys"] # ramené au plot a l'hectare et à l'année

  # Mortality type 3
  #dfplot[,"treeNbrJ.M"] # number of event on the plot only ! 
  dfplot[,"sp.mortality.plot.number.rate"] <- dfplot[,"treeNbrJ.M"]/dfplot[,"treeNbrJ.IMall"] # number of event ramené à la taille du plot (nbr arbre)
  dfplot[,"sp.mortality.plot.number.rate.yr"] <- dfplot[,"sp.mortality.plot.number.rate"]/dfplot[,"yearsbetweensurveys"] # ramené aussi à l'année 

  # Correlation recruitment: 
  # cor.test(dfplot$treeNbrJ.R,dfplot$sp.recruitment.ha.weight)
  # cor.test(dfplot$sp.recruitment.plot.rate,dfplot$sp.recruitment.ha.weight.rate,na.rm=T)
  # cor.test(dfplot$sp.recruitment.plot.rate,dfplot$sp.recruitment.plot.number.rate,na.rm=T)
  # cor.test(dfplot$sp.recruitment.ha.weight.rate,dfplot$sp.recruitment.plot.number.rate,na.rm=T)
  # 
  # # Test correlation entre les différents types de mortalité:
  # cor.test(dfplot$treeNbrJ.M,dfplot$sp.mortality.ha.weight)
  # cor.test(dfplot$sp.mortality.plot.rate,dfplot$sp.mortality.ha.weight.rate,na.rm=T)
  # cor.test(dfplot$sp.mortality.plot.rate,dfplot$sp.mortality.plot.number.rate,na.rm=T)
  # cor.test(dfplot$sp.mortality.ha.weight.rate,dfplot$sp.mortality.plot.number.rate,na.rm=T)
  # 
  
  if (CODE%in%Allcode[c(2:4,6,7,9,13,19)]){
    # years offset both part, sum offset both, mortality weight rate year
    try(rW.yO.sumO.mRY1 <- glmmTMB(sp.recruitment.ha.weight.2 ~ logdbhJ.IMall.plot.mean + 
                                    offset(log(yearsbetweensurveys)) + 
                                    offset(log(sp.SUM.ha.weight2)) + 
                                    (1 | country) + logtreeNbrJ + Plotcat * mean_spei12 + Plotcat *  sqrtBA.O.plot.1 + Plotcat * logBAj.plot.1 + Plotcat * sqrtBAIj.plot.1.mean.J.I + 
                                    Plotcat * logsp.mortality.plot.rate.yr +
                                    mean_spei12 * logBAj.plot.1 + logdbhJ.IMall.plot.mean * logBAj.plot.1 + sqrtBAIj.plot.1.mean.J.I * logBAj.plot.1 + mean_spei12 * sqrtBA.O.plot.1 + logdbhJ.IMall.plot.mean *sqrtBA.O.plot.1 + sqrtBAIj.plot.1.mean.J.I * sqrtBA.O.plot.1,
                                  zi=~logdbhJ.IMall.plot.mean + 
                                    offset(log(yearsbetweensurveys)) + 
                                    offset(log(sp.SUM.ha.weight2)) + 
                                    (1 | country) + logtreeNbrJ + Plotcat * mean_spei12 + Plotcat *  sqrtBA.O.plot.1 + Plotcat * logBAj.plot.1 + Plotcat * sqrtBAIj.plot.1.mean.J.I + 
                                    Plotcat * logsp.mortality.plot.rate.yr +
                                    mean_spei12 * logBAj.plot.1 + logdbhJ.IMall.plot.mean * logBAj.plot.1 + sqrtBAIj.plot.1.mean.J.I * logBAj.plot.1 + mean_spei12 * sqrtBA.O.plot.1 + logdbhJ.IMall.plot.mean *sqrtBA.O.plot.1 + sqrtBAIj.plot.1.mean.J.I * sqrtBA.O.plot.1,family=nbinom1, data=dfplot),silent=T)
    
    try(Saving(rW.yO.sumO.mRY1,CODE),silent=T)
    try(print(summary(rW.yO.sumO.mRY1)),silent=T)
    rm(rW.yO.sumO.mRY1)
    
    # idem minus treenbr
    try(rW.yO.sumO.mRY1.nT <- glmmTMB(sp.recruitment.ha.weight.2 ~ logdbhJ.IMall.plot.mean + 
                                      offset(log(yearsbetweensurveys)) + 
                                      offset(log(sp.SUM.ha.weight2)) + 
                                      (1 | country) + Plotcat * mean_spei12 + Plotcat *  sqrtBA.O.plot.1 + Plotcat * logBAj.plot.1 + Plotcat * sqrtBAIj.plot.1.mean.J.I + 
                                      Plotcat * logsp.mortality.plot.rate.yr +
                                      mean_spei12 * logBAj.plot.1 + logdbhJ.IMall.plot.mean * logBAj.plot.1 + sqrtBAIj.plot.1.mean.J.I * logBAj.plot.1 + mean_spei12 * sqrtBA.O.plot.1 + logdbhJ.IMall.plot.mean *sqrtBA.O.plot.1 + sqrtBAIj.plot.1.mean.J.I * sqrtBA.O.plot.1,
                                    zi=~logdbhJ.IMall.plot.mean + 
                                      offset(log(yearsbetweensurveys)) + 
                                      offset(log(sp.SUM.ha.weight2)) + 
                                      (1 | country) + Plotcat * mean_spei12 + Plotcat *  sqrtBA.O.plot.1 + Plotcat * logBAj.plot.1 + Plotcat * sqrtBAIj.plot.1.mean.J.I + 
                                      Plotcat * logsp.mortality.plot.rate.yr +
                                      mean_spei12 * logBAj.plot.1 + logdbhJ.IMall.plot.mean * logBAj.plot.1 + sqrtBAIj.plot.1.mean.J.I * logBAj.plot.1 + mean_spei12 * sqrtBA.O.plot.1 + logdbhJ.IMall.plot.mean *sqrtBA.O.plot.1 + sqrtBAIj.plot.1.mean.J.I * sqrtBA.O.plot.1,family=nbinom1, data=dfplot),silent=T)
    
    try(Saving(rW.yO.sumO.mRY1.nT,CODE),silent=T)
    try(print(summary(rW.yO.sumO.mRY1.nT)),silent=T)
    rm(rW.yO.sumO.mRY1.nT)
    
  } 
  
  if (CODE%in%Allcode[c(1,10,12,16,18,20)]){
      
      # years offset both part, sum offset both, mortality weight rate year
      try(rW.yO.sumO.mRY1 <- glmmTMB(sp.recruitment.ha.weight.2 ~ logdbhJ.IMall.plot.mean + 
                                      offset(log(yearsbetweensurveys)) + 
                                      offset(log(sp.SUM.ha.weight2)) + 
                                      logtreeNbrJ + Plotcat * mean_spei12 + Plotcat *  sqrtBA.O.plot.1 + Plotcat * logBAj.plot.1 + Plotcat * sqrtBAIj.plot.1.mean.J.I + 
                                      Plotcat * logsp.mortality.plot.rate.yr +
                                      mean_spei12 * logBAj.plot.1 + logdbhJ.IMall.plot.mean * logBAj.plot.1 + sqrtBAIj.plot.1.mean.J.I * logBAj.plot.1 + mean_spei12 * sqrtBA.O.plot.1 + logdbhJ.IMall.plot.mean *sqrtBA.O.plot.1 + sqrtBAIj.plot.1.mean.J.I * sqrtBA.O.plot.1,
                                    zi=~logdbhJ.IMall.plot.mean + 
                                      offset(log(yearsbetweensurveys)) + 
                                      offset(log(sp.SUM.ha.weight2)) + 
                                      logtreeNbrJ + Plotcat * mean_spei12 + Plotcat *  sqrtBA.O.plot.1 + Plotcat * logBAj.plot.1 + Plotcat * sqrtBAIj.plot.1.mean.J.I + 
                                      Plotcat * logsp.mortality.plot.rate.yr +
                                      mean_spei12 * logBAj.plot.1 + logdbhJ.IMall.plot.mean * logBAj.plot.1 + sqrtBAIj.plot.1.mean.J.I * logBAj.plot.1 + mean_spei12 * sqrtBA.O.plot.1 + logdbhJ.IMall.plot.mean *sqrtBA.O.plot.1 + sqrtBAIj.plot.1.mean.J.I * sqrtBA.O.plot.1,family=nbinom1, data=dfplot),silent=T)
      
      try(Saving(rW.yO.sumO.mRY1,CODE),silent=T)
      try(print(summary(rW.yO.sumO.mRY1)),silent=T)
      rm(rW.yO.sumO.mRY1)
      
      try(rW.yO.sumO.mRY1.nT <- glmmTMB(sp.recruitment.ha.weight.2 ~ logdbhJ.IMall.plot.mean + 
                                      offset(log(yearsbetweensurveys)) + 
                                      offset(log(sp.SUM.ha.weight2)) + 
                                      Plotcat * mean_spei12 + Plotcat *  sqrtBA.O.plot.1 + Plotcat * logBAj.plot.1 + Plotcat * sqrtBAIj.plot.1.mean.J.I + 
                                      Plotcat * logsp.mortality.plot.rate.yr +
                                      mean_spei12 * logBAj.plot.1 + logdbhJ.IMall.plot.mean * logBAj.plot.1 + sqrtBAIj.plot.1.mean.J.I * logBAj.plot.1 + mean_spei12 * sqrtBA.O.plot.1 + logdbhJ.IMall.plot.mean *sqrtBA.O.plot.1 + sqrtBAIj.plot.1.mean.J.I * sqrtBA.O.plot.1,
                                    zi=~logdbhJ.IMall.plot.mean + 
                                      offset(log(yearsbetweensurveys)) + 
                                      offset(log(sp.SUM.ha.weight2)) + 
                                      Plotcat * mean_spei12 + Plotcat *  sqrtBA.O.plot.1 + Plotcat * logBAj.plot.1 + Plotcat * sqrtBAIj.plot.1.mean.J.I + 
                                      Plotcat * logsp.mortality.plot.rate.yr +
                                      mean_spei12 * logBAj.plot.1 + logdbhJ.IMall.plot.mean * logBAj.plot.1 + sqrtBAIj.plot.1.mean.J.I * logBAj.plot.1 + mean_spei12 * sqrtBA.O.plot.1 + logdbhJ.IMall.plot.mean *sqrtBA.O.plot.1 + sqrtBAIj.plot.1.mean.J.I * sqrtBA.O.plot.1,family=nbinom1, data=dfplot),silent=T)
      
      try(Saving(rW.yO.sumO.mRY1.nT,CODE),silent=T)
      try(print(summary(rW.yO.sumO.mRY1.nT)),silent=T)
      rm(rW.yO.sumO.mRY1.nT)
      
      
      }
  
  if (CODE%in%Allcode[c(5,8,17)]){
        # years offset both part, sum offset both, mortality weight rate year
        try(rW.yO.sumO.mRY1 <- glmmTMB(sp.recruitment.ha.weight.2 ~ logdbhJ.IMall.plot.mean + 
                                        offset(log(yearsbetweensurveys)) + 
                                        offset(log(sp.SUM.ha.weight2)) + 
                                        (1 | country) + logtreeNbrJ + Plotcat * mean_spei12 + Plotcat *  sqrtBA.O.plot.1 + Plotcat * logBAj.plot.1 + Plotcat * sqrtBAIj.plot.1.mean.J.I + 
                                        logsp.mortality.plot.rate.yr +
                                        mean_spei12 * logBAj.plot.1 + logdbhJ.IMall.plot.mean * logBAj.plot.1 + sqrtBAIj.plot.1.mean.J.I * logBAj.plot.1 + mean_spei12 * sqrtBA.O.plot.1 + logdbhJ.IMall.plot.mean *sqrtBA.O.plot.1 + sqrtBAIj.plot.1.mean.J.I * sqrtBA.O.plot.1,
                                      zi=~logdbhJ.IMall.plot.mean + 
                                        offset(log(yearsbetweensurveys)) + 
                                        offset(log(sp.SUM.ha.weight2)) + 
                                        (1 | country) + logtreeNbrJ + Plotcat * mean_spei12 + Plotcat *  sqrtBA.O.plot.1 + Plotcat * logBAj.plot.1 + Plotcat * sqrtBAIj.plot.1.mean.J.I + 
                                        logsp.mortality.plot.rate.yr +
                                        mean_spei12 * logBAj.plot.1 + logdbhJ.IMall.plot.mean * logBAj.plot.1 + sqrtBAIj.plot.1.mean.J.I * logBAj.plot.1 + mean_spei12 * sqrtBA.O.plot.1 + logdbhJ.IMall.plot.mean *sqrtBA.O.plot.1 + sqrtBAIj.plot.1.mean.J.I * sqrtBA.O.plot.1,family=nbinom1, data=dfplot),silent=T)
        
        try(Saving(rW.yO.sumO.mRY1,CODE),silent=T)
        try(print(summary(rW.yO.sumO.mRY1)),silent=T)
        rm(rW.yO.sumO.mRY1)
        
        try(rW.yO.sumO.mRY1.nT <- glmmTMB(sp.recruitment.ha.weight.2 ~ logdbhJ.IMall.plot.mean + 
                                        offset(log(yearsbetweensurveys)) + 
                                        offset(log(sp.SUM.ha.weight2)) + 
                                        (1 | country) + Plotcat * mean_spei12 + Plotcat *  sqrtBA.O.plot.1 + Plotcat * logBAj.plot.1 + Plotcat * sqrtBAIj.plot.1.mean.J.I + 
                                        logsp.mortality.plot.rate.yr +
                                        mean_spei12 * logBAj.plot.1 + logdbhJ.IMall.plot.mean * logBAj.plot.1 + sqrtBAIj.plot.1.mean.J.I * logBAj.plot.1 + mean_spei12 * sqrtBA.O.plot.1 + logdbhJ.IMall.plot.mean *sqrtBA.O.plot.1 + sqrtBAIj.plot.1.mean.J.I * sqrtBA.O.plot.1,
                                      zi=~logdbhJ.IMall.plot.mean + 
                                        offset(log(yearsbetweensurveys)) + 
                                        offset(log(sp.SUM.ha.weight2)) + 
                                        (1 | country) + Plotcat * mean_spei12 + Plotcat *  sqrtBA.O.plot.1 + Plotcat * logBAj.plot.1 + Plotcat * sqrtBAIj.plot.1.mean.J.I + 
                                        logsp.mortality.plot.rate.yr +
                                        mean_spei12 * logBAj.plot.1 + logdbhJ.IMall.plot.mean * logBAj.plot.1 + sqrtBAIj.plot.1.mean.J.I * logBAj.plot.1 + mean_spei12 * sqrtBA.O.plot.1 + logdbhJ.IMall.plot.mean *sqrtBA.O.plot.1 + sqrtBAIj.plot.1.mean.J.I * sqrtBA.O.plot.1,family=nbinom1, data=dfplot),silent=T)
        
        try(Saving(rW.yO.sumO.mRY1.nT,CODE),silent=T)
        try(print(summary(rW.yO.sumO.mRY1.nT)),silent=T)
        rm(rW.yO.sumO.mRY1.nT)
      
        }
  
  
  if (CODE%in%Allcode[c(11)]){
          # years offset both part, sum offset both, mortality weight rate year
          try(rW.yO.sumO.mRY1 <- glmmTMB(sp.recruitment.ha.weight.2 ~ logdbhJ.IMall.plot.mean + 
                                          offset(log(yearsbetweensurveys)) + 
                                          offset(log(sp.SUM.ha.weight2)) + 
                                          logtreeNbrJ + Plotcat * mean_spei12 + Plotcat *  sqrtBA.O.plot.1 + Plotcat * logBAj.plot.1 + Plotcat * sqrtBAIj.plot.1.mean.J.I + 
                                          logsp.mortality.plot.rate.yr +
                                          mean_spei12 * logBAj.plot.1 + logdbhJ.IMall.plot.mean * logBAj.plot.1 + sqrtBAIj.plot.1.mean.J.I * logBAj.plot.1 + mean_spei12 * sqrtBA.O.plot.1 + logdbhJ.IMall.plot.mean *sqrtBA.O.plot.1 + sqrtBAIj.plot.1.mean.J.I * sqrtBA.O.plot.1,
                                        zi=~logdbhJ.IMall.plot.mean + 
                                          offset(log(yearsbetweensurveys)) + 
                                          offset(log(sp.SUM.ha.weight2)) + 
                                          logtreeNbrJ + Plotcat * mean_spei12 + Plotcat *  sqrtBA.O.plot.1 + Plotcat * logBAj.plot.1 + Plotcat * sqrtBAIj.plot.1.mean.J.I + 
                                          logsp.mortality.plot.rate.yr +
                                          mean_spei12 * logBAj.plot.1 + logdbhJ.IMall.plot.mean * logBAj.plot.1 + sqrtBAIj.plot.1.mean.J.I * logBAj.plot.1 + mean_spei12 * sqrtBA.O.plot.1 + logdbhJ.IMall.plot.mean *sqrtBA.O.plot.1 + sqrtBAIj.plot.1.mean.J.I * sqrtBA.O.plot.1,family=nbinom1, data=dfplot),silent=T)
          
          try(Saving(rW.yO.sumO.mRY1,CODE),silent=T)
          try(print(summary(rW.yO.sumO.mRY1)),silent=T)
          rm(rW.yO.sumO.mRY1)
          
          try(rW.yO.sumO.mRY1.nT <- glmmTMB(sp.recruitment.ha.weight.2 ~ logdbhJ.IMall.plot.mean + 
                                          offset(log(yearsbetweensurveys)) + 
                                          offset(log(sp.SUM.ha.weight2)) + 
                                          Plotcat * mean_spei12 + Plotcat *  sqrtBA.O.plot.1 + Plotcat * logBAj.plot.1 + Plotcat * sqrtBAIj.plot.1.mean.J.I + 
                                          logsp.mortality.plot.rate.yr +
                                          mean_spei12 * logBAj.plot.1 + logdbhJ.IMall.plot.mean * logBAj.plot.1 + sqrtBAIj.plot.1.mean.J.I * logBAj.plot.1 + mean_spei12 * sqrtBA.O.plot.1 + logdbhJ.IMall.plot.mean *sqrtBA.O.plot.1 + sqrtBAIj.plot.1.mean.J.I * sqrtBA.O.plot.1,
                                        zi=~logdbhJ.IMall.plot.mean + 
                                          offset(log(yearsbetweensurveys)) + 
                                          offset(log(sp.SUM.ha.weight2)) + 
                                          Plotcat * mean_spei12 + Plotcat *  sqrtBA.O.plot.1 + Plotcat * logBAj.plot.1 + Plotcat * sqrtBAIj.plot.1.mean.J.I + 
                                          logsp.mortality.plot.rate.yr +
                                          mean_spei12 * logBAj.plot.1 + logdbhJ.IMall.plot.mean * logBAj.plot.1 + sqrtBAIj.plot.1.mean.J.I * logBAj.plot.1 + mean_spei12 * sqrtBA.O.plot.1 + logdbhJ.IMall.plot.mean *sqrtBA.O.plot.1 + sqrtBAIj.plot.1.mean.J.I * sqrtBA.O.plot.1,family=nbinom1, data=dfplot),silent=T)
          
          try(Saving(rW.yO.sumO.mRY1.nT,CODE),silent=T)
          try(print(summary(rW.yO.sumO.mRY1.nT)),silent=T)
          rm(rW.yO.sumO.mRY1.nT)
          
  }
  
  },dfplot=DFPLOT,CODE=Allcode,SIMPLIFY=F,mc.cores=20)
        
        
      
      


      
    
    
# years fixed both part, sum fixed both, mortality weight 
try(rW.yF.sumF.mC <- glmmTMB(sp.recruitment.ha.weight.2 ~ logdbhJ.IMall.plot.mean + 
                               yearsbetweensurveys + 
                               sp.SUM.ha.weight2 + 
                               (1 | country) + logtreeNbrJ + Plotcat * mean_spei12 + Plotcat *  sqrtBA.O.plot.1 + Plotcat * logBAj.plot.1 + Plotcat * sqrtBAIj.plot.1.mean.J.I + 
                               Plotcat * sp.mortality.ha.weight +
                               mean_spei12 * logBAj.plot.1 + logdbhJ.IMall.plot.mean * logBAj.plot.1 + sqrtBAIj.plot.1.mean.J.I * logBAj.plot.1 + mean_spei12 * sqrtBA.O.plot.1 + logdbhJ.IMall.plot.mean *sqrtBA.O.plot.1 + sqrtBAIj.plot.1.mean.J.I * sqrtBA.O.plot.1,
                             zi=~.,family=nbinom1, data=dfplot),silent=T)

try(Saving(rW.yF.sumF.mC,CODE),silent=T)
try(print(summary(rW.yF.sumF.mC)),silent=T)
rm(rW.yF.sumF.mC)

    # years fixed both part, sum fixed both, mortality weight rate
    try(rW.yF.sumF.mR <- glmmTMB(sp.recruitment.ha.weight.2 ~ logdbhJ.IMall.plot.mean + 
                                   yearsbetweensurveys + 
                                   sp.SUM.ha.weight2 + 
                                   (1 | country) + logtreeNbrJ + Plotcat * mean_spei12 + Plotcat *  sqrtBA.O.plot.1 + Plotcat * logBAj.plot.1 + Plotcat * sqrtBAIj.plot.1.mean.J.I + 
                                   Plotcat * sp.mortality.ha.weight.rate +
                                   mean_spei12 * logBAj.plot.1 + logdbhJ.IMall.plot.mean * logBAj.plot.1 + sqrtBAIj.plot.1.mean.J.I * logBAj.plot.1 + mean_spei12 * sqrtBA.O.plot.1 + logdbhJ.IMall.plot.mean *sqrtBA.O.plot.1 + sqrtBAIj.plot.1.mean.J.I * sqrtBA.O.plot.1,
                                 zi=~.,family=nbinom1, data=dfplot),silent=T)
    try(Saving(rW.yF.sumF.mR,CODE),silent=T)
    try(print(summary(rW.yF.sumF.mR)),silent=T)
    rm(rW.yF.sumF.mR)
    
    # years fixed both part, sum fixed both, mortality weight rate year
    try(rW.yF.sumF.mRY <- glmmTMB(sp.recruitment.ha.weight.2 ~ logdbhJ.IMall.plot.mean + 
                                    yearsbetweensurveys + 
                                    sp.SUM.ha.weight2 + 
                                    (1 | country) + logtreeNbrJ + Plotcat * mean_spei12 + Plotcat *  sqrtBA.O.plot.1 + Plotcat * logBAj.plot.1 + Plotcat * sqrtBAIj.plot.1.mean.J.I + 
                                    Plotcat * sp.mortality.ha.weight.rate.yr +
                                    mean_spei12 * logBAj.plot.1 + logdbhJ.IMall.plot.mean * logBAj.plot.1 + sqrtBAIj.plot.1.mean.J.I * logBAj.plot.1 + mean_spei12 * sqrtBA.O.plot.1 + logdbhJ.IMall.plot.mean *sqrtBA.O.plot.1 + sqrtBAIj.plot.1.mean.J.I * sqrtBA.O.plot.1,
                                  zi=~.,family=nbinom1, data=dfplot),silent=T)
    try(Saving(rW.yF.sumF.mRY,CODE),silent=T)
    try(print(summary(rW.yF.sumF.mRY)),silent=T)
    rm(rW.yF.sumF.mRY)
    
    # years offset both part, sum fixed both, mortality weight
    try(rW.yO.sumF.mC <- glmmTMB(sp.recruitment.ha.weight.2 ~ logdbhJ.IMall.plot.mean + 
                                   offset(log(yearsbetweensurveys)) + 
                                   sp.SUM.ha.weight2 + 
                                   (1 | country) + logtreeNbrJ + Plotcat * mean_spei12 + Plotcat *  sqrtBA.O.plot.1 + Plotcat * logBAj.plot.1 + Plotcat * sqrtBAIj.plot.1.mean.J.I + 
                                   Plotcat * sp.mortality.ha.weight +
                                   mean_spei12 * logBAj.plot.1 + logdbhJ.IMall.plot.mean * logBAj.plot.1 + sqrtBAIj.plot.1.mean.J.I * logBAj.plot.1 + mean_spei12 * sqrtBA.O.plot.1 + logdbhJ.IMall.plot.mean *sqrtBA.O.plot.1 + sqrtBAIj.plot.1.mean.J.I * sqrtBA.O.plot.1,
                                 zi=~logdbhJ.IMall.plot.mean + 
                                   offset(log(yearsbetweensurveys)) + 
                                   sp.SUM.ha.weight2 + 
                                   (1 | country) + logtreeNbrJ + Plotcat * mean_spei12 + Plotcat *  sqrtBA.O.plot.1 + Plotcat * logBAj.plot.1 + Plotcat * sqrtBAIj.plot.1.mean.J.I + 
                                   Plotcat * sp.mortality.ha.weight +
                                   mean_spei12 * logBAj.plot.1 + logdbhJ.IMall.plot.mean * logBAj.plot.1 + sqrtBAIj.plot.1.mean.J.I * logBAj.plot.1 + mean_spei12 * sqrtBA.O.plot.1 + logdbhJ.IMall.plot.mean *sqrtBA.O.plot.1 + sqrtBAIj.plot.1.mean.J.I * sqrtBA.O.plot.1,family=nbinom1, data=dfplot),silent=T)
    
    try(Saving(rW.yO.sumF.mC,CODE),silent=T)
    try(print(summary(rW.yO.sumF.mC)),silent=T)
    rm(rW.yO.sumF.mC)
    
    # years offset both part, sum fixed both, mortality weight rate
    try(rW.yO.sumF.mR <- glmmTMB(sp.recruitment.ha.weight.2 ~ logdbhJ.IMall.plot.mean + 
                                   offset(log(yearsbetweensurveys)) + 
                                   sp.SUM.ha.weight2 + 
                                   (1 | country) + logtreeNbrJ + Plotcat * mean_spei12 + Plotcat *  sqrtBA.O.plot.1 + Plotcat * logBAj.plot.1 + Plotcat * sqrtBAIj.plot.1.mean.J.I + 
                                   Plotcat * sp.mortality.ha.weight.rate +
                                   mean_spei12 * logBAj.plot.1 + logdbhJ.IMall.plot.mean * logBAj.plot.1 + sqrtBAIj.plot.1.mean.J.I * logBAj.plot.1 + mean_spei12 * sqrtBA.O.plot.1 + logdbhJ.IMall.plot.mean *sqrtBA.O.plot.1 + sqrtBAIj.plot.1.mean.J.I * sqrtBA.O.plot.1,
                                 zi=~logdbhJ.IMall.plot.mean + 
                                   offset(log(yearsbetweensurveys)) + 
                                   sp.SUM.ha.weight2 + 
                                   (1 | country) + logtreeNbrJ + Plotcat * mean_spei12 + Plotcat *  sqrtBA.O.plot.1 + Plotcat * logBAj.plot.1 + Plotcat * sqrtBAIj.plot.1.mean.J.I + 
                                   Plotcat * sp.mortality.ha.weight.rate +
                                   mean_spei12 * logBAj.plot.1 + logdbhJ.IMall.plot.mean * logBAj.plot.1 + sqrtBAIj.plot.1.mean.J.I * logBAj.plot.1 + mean_spei12 * sqrtBA.O.plot.1 + logdbhJ.IMall.plot.mean *sqrtBA.O.plot.1 + sqrtBAIj.plot.1.mean.J.I * sqrtBA.O.plot.1,family=nbinom1, data=dfplot),silent=T)
    
    try(Saving(rW.yO.sumF.mR,CODE),silent=T)
    try(print(summary(rW.yO.sumF.mR)),silent=T)
    rm(rW.yO.sumF.mR)
    
    # years fixed both part, sum fixed both, mortality weight rate year
    try(rW.yO.sumF.mRY <- glmmTMB(sp.recruitment.ha.weight.2 ~ logdbhJ.IMall.plot.mean + 
                                   offset(log(yearsbetweensurveys)) + 
                                   sp.SUM.ha.weight2 + 
                                   (1 | country) + logtreeNbrJ + Plotcat * mean_spei12 + Plotcat *  sqrtBA.O.plot.1 + Plotcat * logBAj.plot.1 + Plotcat * sqrtBAIj.plot.1.mean.J.I + 
                                   Plotcat * sp.mortality.ha.weight.rate.yr +
                                   mean_spei12 * logBAj.plot.1 + logdbhJ.IMall.plot.mean * logBAj.plot.1 + sqrtBAIj.plot.1.mean.J.I * logBAj.plot.1 + mean_spei12 * sqrtBA.O.plot.1 + logdbhJ.IMall.plot.mean *sqrtBA.O.plot.1 + sqrtBAIj.plot.1.mean.J.I * sqrtBA.O.plot.1,
                                 zi=~logdbhJ.IMall.plot.mean + 
                                   offset(log(yearsbetweensurveys)) + 
                                   sp.SUM.ha.weight2 + 
                                   (1 | country) + logtreeNbrJ + Plotcat * mean_spei12 + Plotcat *  sqrtBA.O.plot.1 + Plotcat * logBAj.plot.1 + Plotcat * sqrtBAIj.plot.1.mean.J.I + 
                                   Plotcat * sp.mortality.ha.weight.rate.yr +
                                   mean_spei12 * logBAj.plot.1 + logdbhJ.IMall.plot.mean * logBAj.plot.1 + sqrtBAIj.plot.1.mean.J.I * logBAj.plot.1 + mean_spei12 * sqrtBA.O.plot.1 + logdbhJ.IMall.plot.mean *sqrtBA.O.plot.1 + sqrtBAIj.plot.1.mean.J.I * sqrtBA.O.plot.1,family=nbinom1, data=dfplot),silent=T)
    
    try(Saving(rW.yO.sumF.mRY,CODE),silent=T)
    try(print(summary(rW.yO.sumF.mRY)),silent=T)
    rm(rW.yO.sumF.mRY)
    
    # years fixed both part, sum offset both, mortality weight
    try(rW.yF.sumO.mC <- glmmTMB(sp.recruitment.ha.weight.2 ~ logdbhJ.IMall.plot.mean + 
                                   yearsbetweensurveys + 
                                   offset(log(sp.SUM.ha.weight2)) + 
                                   (1 | country) + logtreeNbrJ + Plotcat * mean_spei12 + Plotcat *  sqrtBA.O.plot.1 + Plotcat * logBAj.plot.1 + Plotcat * sqrtBAIj.plot.1.mean.J.I + 
                                   Plotcat * sp.mortality.ha.weight +
                                   mean_spei12 * logBAj.plot.1 + logdbhJ.IMall.plot.mean * logBAj.plot.1 + sqrtBAIj.plot.1.mean.J.I * logBAj.plot.1 + mean_spei12 * sqrtBA.O.plot.1 + logdbhJ.IMall.plot.mean *sqrtBA.O.plot.1 + sqrtBAIj.plot.1.mean.J.I * sqrtBA.O.plot.1,
                                 zi=~logdbhJ.IMall.plot.mean + 
                                   yearsbetweensurveys + 
                                   offset(log(sp.SUM.ha.weight2)) + 
                                   (1 | country) + logtreeNbrJ + Plotcat * mean_spei12 + Plotcat *  sqrtBA.O.plot.1 + Plotcat * logBAj.plot.1 + Plotcat * sqrtBAIj.plot.1.mean.J.I + 
                                   Plotcat * sp.mortality.ha.weight +
                                   mean_spei12 * logBAj.plot.1 + logdbhJ.IMall.plot.mean * logBAj.plot.1 + sqrtBAIj.plot.1.mean.J.I * logBAj.plot.1 + mean_spei12 * sqrtBA.O.plot.1 + logdbhJ.IMall.plot.mean *sqrtBA.O.plot.1 + sqrtBAIj.plot.1.mean.J.I * sqrtBA.O.plot.1,family=nbinom1, data=dfplot),silent=T)
    
    try(Saving(rW.yF.sumO.mC,CODE),silent=T)
    try(print(summary(rW.yF.sumO.mC)),silent=T)
    rm(rW.yF.sumO.mC)
    
    # years fixed both part, sum offset both, mortality weight rate
    try(rW.yF.sumO.mR <- glmmTMB(sp.recruitment.ha.weight.2 ~ logdbhJ.IMall.plot.mean + 
                                   yearsbetweensurveys + 
                                   offset(log(sp.SUM.ha.weight2)) +
                                   (1 | country) + logtreeNbrJ + Plotcat * mean_spei12 + Plotcat *  sqrtBA.O.plot.1 + Plotcat * logBAj.plot.1 + Plotcat * sqrtBAIj.plot.1.mean.J.I + 
                                   Plotcat * sp.mortality.ha.weight.rate +
                                   mean_spei12 * logBAj.plot.1 + logdbhJ.IMall.plot.mean * logBAj.plot.1 + sqrtBAIj.plot.1.mean.J.I * logBAj.plot.1 + mean_spei12 * sqrtBA.O.plot.1 + logdbhJ.IMall.plot.mean *sqrtBA.O.plot.1 + sqrtBAIj.plot.1.mean.J.I * sqrtBA.O.plot.1,
                                 zi=~logdbhJ.IMall.plot.mean + 
                                   yearsbetweensurveys + 
                                   offset(log(sp.SUM.ha.weight2)) + 
                                   (1 | country) + logtreeNbrJ + Plotcat * mean_spei12 + Plotcat *  sqrtBA.O.plot.1 + Plotcat * logBAj.plot.1 + Plotcat * sqrtBAIj.plot.1.mean.J.I + 
                                   Plotcat * sp.mortality.ha.weight.rate +
                                   mean_spei12 * logBAj.plot.1 + logdbhJ.IMall.plot.mean * logBAj.plot.1 + sqrtBAIj.plot.1.mean.J.I * logBAj.plot.1 + mean_spei12 * sqrtBA.O.plot.1 + logdbhJ.IMall.plot.mean *sqrtBA.O.plot.1 + sqrtBAIj.plot.1.mean.J.I * sqrtBA.O.plot.1,family=nbinom1, data=dfplot),silent=T)
    
    try(Saving(rW.yF.sumO.mR,CODE),silent=T)
    try(print(summary(rW.yF.sumO.mR)),silent=T)
    rm(rW.yF.sumO.mR)
    
    # years fixed both part, sum offset both, mortality weight rate year
    try(rW.yF.sumO.mRY <- glmmTMB(sp.recruitment.ha.weight.2 ~ logdbhJ.IMall.plot.mean + 
                                    yearsbetweensurveys + 
                                    offset(log(sp.SUM.ha.weight2)) + 
                                    (1 | country) + logtreeNbrJ + Plotcat * mean_spei12 + Plotcat *  sqrtBA.O.plot.1 + Plotcat * logBAj.plot.1 + Plotcat * sqrtBAIj.plot.1.mean.J.I + 
                                    Plotcat * sp.mortality.ha.weight.rate.yr +
                                    mean_spei12 * logBAj.plot.1 + logdbhJ.IMall.plot.mean * logBAj.plot.1 + sqrtBAIj.plot.1.mean.J.I * logBAj.plot.1 + mean_spei12 * sqrtBA.O.plot.1 + logdbhJ.IMall.plot.mean *sqrtBA.O.plot.1 + sqrtBAIj.plot.1.mean.J.I * sqrtBA.O.plot.1,
                                  zi=~logdbhJ.IMall.plot.mean + 
                                    yearsbetweensurveys + 
                                    offset(log(sp.SUM.ha.weight2)) + 
                                    (1 | country) + logtreeNbrJ + Plotcat * mean_spei12 + Plotcat *  sqrtBA.O.plot.1 + Plotcat * logBAj.plot.1 + Plotcat * sqrtBAIj.plot.1.mean.J.I + 
                                    Plotcat * sp.mortality.ha.weight.rate.yr +
                                    mean_spei12 * logBAj.plot.1 + logdbhJ.IMall.plot.mean * logBAj.plot.1 + sqrtBAIj.plot.1.mean.J.I * logBAj.plot.1 + mean_spei12 * sqrtBA.O.plot.1 + logdbhJ.IMall.plot.mean *sqrtBA.O.plot.1 + sqrtBAIj.plot.1.mean.J.I * sqrtBA.O.plot.1,family=nbinom1, data=dfplot),silent=T)
    
    try(Saving(rW.yF.sumO.mRY,CODE),silent=T)
    try(print(summary(rW.yF.sumO.mRY)),silent=T)
    rm(rW.yF.sumO.mRY)
    
    
    # years offset both part, sum offset both, mortality weight
    try(rW.yO.sumO.mC <- glmmTMB(sp.recruitment.ha.weight.2 ~ logdbhJ.IMall.plot.mean + 
                                   offset(log(yearsbetweensurveys)) + 
                                   offset(log(sp.SUM.ha.weight2)) + 
                                   (1 | country) + logtreeNbrJ + Plotcat * mean_spei12 + Plotcat *  sqrtBA.O.plot.1 + Plotcat * logBAj.plot.1 + Plotcat * sqrtBAIj.plot.1.mean.J.I + 
                                   Plotcat * sp.mortality.ha.weight +
                                   mean_spei12 * logBAj.plot.1 + logdbhJ.IMall.plot.mean * logBAj.plot.1 + sqrtBAIj.plot.1.mean.J.I * logBAj.plot.1 + mean_spei12 * sqrtBA.O.plot.1 + logdbhJ.IMall.plot.mean *sqrtBA.O.plot.1 + sqrtBAIj.plot.1.mean.J.I * sqrtBA.O.plot.1,
                                 zi=~logdbhJ.IMall.plot.mean + 
                                   offset(log(yearsbetweensurveys)) + 
                                   offset(log(sp.SUM.ha.weight2)) + 
                                   (1 | country) + logtreeNbrJ + Plotcat * mean_spei12 + Plotcat *  sqrtBA.O.plot.1 + Plotcat * logBAj.plot.1 + Plotcat * sqrtBAIj.plot.1.mean.J.I + 
                                   Plotcat * sp.mortality.ha.weight +
                                   mean_spei12 * logBAj.plot.1 + logdbhJ.IMall.plot.mean * logBAj.plot.1 + sqrtBAIj.plot.1.mean.J.I * logBAj.plot.1 + mean_spei12 * sqrtBA.O.plot.1 + logdbhJ.IMall.plot.mean *sqrtBA.O.plot.1 + sqrtBAIj.plot.1.mean.J.I * sqrtBA.O.plot.1,family=nbinom1, data=dfplot),silent=T)
    
    try(Saving(rW.yO.sumO.mC,CODE),silent=T)
    try(print(summary(rW.yO.sumO.mC)),silent=T)
    rm(rW.yO.sumO.mC)
    
    # years offset both part, sum offset both, mortality weight rate
    try(rW.yO.sumO.mR <- glmmTMB(sp.recruitment.ha.weight.2 ~ logdbhJ.IMall.plot.mean + 
                                   offset(log(yearsbetweensurveys)) + 
                                   offset(log(sp.SUM.ha.weight2)) +
                                   (1 | country) + logtreeNbrJ + Plotcat * mean_spei12 + Plotcat *  sqrtBA.O.plot.1 + Plotcat * logBAj.plot.1 + Plotcat * sqrtBAIj.plot.1.mean.J.I + 
                                   Plotcat * sp.mortality.ha.weight.rate +
                                   mean_spei12 * logBAj.plot.1 + logdbhJ.IMall.plot.mean * logBAj.plot.1 + sqrtBAIj.plot.1.mean.J.I * logBAj.plot.1 + mean_spei12 * sqrtBA.O.plot.1 + logdbhJ.IMall.plot.mean *sqrtBA.O.plot.1 + sqrtBAIj.plot.1.mean.J.I * sqrtBA.O.plot.1,
                                 zi=~logdbhJ.IMall.plot.mean + 
                                   offset(log(yearsbetweensurveys)) + 
                                   offset(log(sp.SUM.ha.weight2)) + 
                                   (1 | country) + logtreeNbrJ + Plotcat * mean_spei12 + Plotcat *  sqrtBA.O.plot.1 + Plotcat * logBAj.plot.1 + Plotcat * sqrtBAIj.plot.1.mean.J.I + 
                                   Plotcat * sp.mortality.ha.weight.rate +
                                   mean_spei12 * logBAj.plot.1 + logdbhJ.IMall.plot.mean * logBAj.plot.1 + sqrtBAIj.plot.1.mean.J.I * logBAj.plot.1 + mean_spei12 * sqrtBA.O.plot.1 + logdbhJ.IMall.plot.mean *sqrtBA.O.plot.1 + sqrtBAIj.plot.1.mean.J.I * sqrtBA.O.plot.1,family=nbinom1, data=dfplot),silent=T)
    
    try(Saving(rW.yO.sumO.mR,CODE),silent=T)
    try(print(summary(rW.yO.sumO.mR)),silent=T)
    rm(rW.yO.sumO.mR)
    

  }
  
  if (CODE%in%Allcode[c(1,10,12,16,18,20)]){
  #if (i%in%c(1,10,12,16,18,20)){
    
    # years fixed both part, sum fixed both, mortality weight 
    try(rW.yF.sumF.mC <- glmmTMB(sp.recruitment.ha.weight.2 ~ logdbhJ.IMall.plot.mean + 
                                   yearsbetweensurveys + 
                                   sp.SUM.ha.weight2 + 
                                   logtreeNbrJ + Plotcat * mean_spei12 + Plotcat *  sqrtBA.O.plot.1 + Plotcat * logBAj.plot.1 + Plotcat * sqrtBAIj.plot.1.mean.J.I + 
                                   Plotcat * sp.mortality.ha.weight +
                                   mean_spei12 * logBAj.plot.1 + logdbhJ.IMall.plot.mean * logBAj.plot.1 + sqrtBAIj.plot.1.mean.J.I * logBAj.plot.1 + mean_spei12 * sqrtBA.O.plot.1 + logdbhJ.IMall.plot.mean *sqrtBA.O.plot.1 + sqrtBAIj.plot.1.mean.J.I * sqrtBA.O.plot.1,
                                 zi=~.,family=nbinom1, data=dfplot),silent=T)
    
    try(Saving(rW.yF.sumF.mC,CODE),silent=T)
    try(print(summary(rW.yF.sumF.mC)),silent=T)
    rm(rW.yF.sumF.mC)
    
    # years fixed both part, sum fixed both, mortality weight rate
    try(rW.yF.sumF.mR <- glmmTMB(sp.recruitment.ha.weight.2 ~ logdbhJ.IMall.plot.mean + 
                                   yearsbetweensurveys + 
                                   sp.SUM.ha.weight2 + 
                                   logtreeNbrJ + Plotcat * mean_spei12 + Plotcat *  sqrtBA.O.plot.1 + Plotcat * logBAj.plot.1 + Plotcat * sqrtBAIj.plot.1.mean.J.I + 
                                   Plotcat * sp.mortality.ha.weight.rate +
                                   mean_spei12 * logBAj.plot.1 + logdbhJ.IMall.plot.mean * logBAj.plot.1 + sqrtBAIj.plot.1.mean.J.I * logBAj.plot.1 + mean_spei12 * sqrtBA.O.plot.1 + logdbhJ.IMall.plot.mean *sqrtBA.O.plot.1 + sqrtBAIj.plot.1.mean.J.I * sqrtBA.O.plot.1,
                                 zi=~.,family=nbinom1, data=dfplot),silent=T)
    try(Saving(rW.yF.sumF.mR,CODE),silent=T)
    try(print(summary(rW.yF.sumF.mR)),silent=T)
    rm(rW.yF.sumF.mR)
    
    # years fixed both part, sum fixed both, mortality weight rate year
    try(rW.yF.sumF.mRY <- glmmTMB(sp.recruitment.ha.weight.2 ~ logdbhJ.IMall.plot.mean + 
                                    yearsbetweensurveys + 
                                    sp.SUM.ha.weight2 + 
                                    logtreeNbrJ + Plotcat * mean_spei12 + Plotcat *  sqrtBA.O.plot.1 + Plotcat * logBAj.plot.1 + Plotcat * sqrtBAIj.plot.1.mean.J.I + 
                                    Plotcat * sp.mortality.ha.weight.rate.yr +
                                    mean_spei12 * logBAj.plot.1 + logdbhJ.IMall.plot.mean * logBAj.plot.1 + sqrtBAIj.plot.1.mean.J.I * logBAj.plot.1 + mean_spei12 * sqrtBA.O.plot.1 + logdbhJ.IMall.plot.mean *sqrtBA.O.plot.1 + sqrtBAIj.plot.1.mean.J.I * sqrtBA.O.plot.1,
                                  zi=~.,family=nbinom1, data=dfplot),silent=T)
    try(Saving(rW.yF.sumF.mRY,CODE),silent=T)
    try(print(summary(rW.yF.sumF.mRY)),silent=T)
    rm(rW.yF.sumF.mRY)
    
    # years offset both part, sum fixed both, mortality weight
    try(rW.yO.sumF.mC <- glmmTMB(sp.recruitment.ha.weight.2 ~ logdbhJ.IMall.plot.mean + 
                                   offset(log(yearsbetweensurveys)) + 
                                   sp.SUM.ha.weight2 + 
                                   logtreeNbrJ + Plotcat * mean_spei12 + Plotcat *  sqrtBA.O.plot.1 + Plotcat * logBAj.plot.1 + Plotcat * sqrtBAIj.plot.1.mean.J.I + 
                                   Plotcat * sp.mortality.ha.weight +
                                   mean_spei12 * logBAj.plot.1 + logdbhJ.IMall.plot.mean * logBAj.plot.1 + sqrtBAIj.plot.1.mean.J.I * logBAj.plot.1 + mean_spei12 * sqrtBA.O.plot.1 + logdbhJ.IMall.plot.mean *sqrtBA.O.plot.1 + sqrtBAIj.plot.1.mean.J.I * sqrtBA.O.plot.1,
                                 zi=~logdbhJ.IMall.plot.mean + 
                                   offset(log(yearsbetweensurveys)) + 
                                   sp.SUM.ha.weight2 + 
                                   logtreeNbrJ + Plotcat * mean_spei12 + Plotcat *  sqrtBA.O.plot.1 + Plotcat * logBAj.plot.1 + Plotcat * sqrtBAIj.plot.1.mean.J.I + 
                                   Plotcat * sp.mortality.ha.weight +
                                   mean_spei12 * logBAj.plot.1 + logdbhJ.IMall.plot.mean * logBAj.plot.1 + sqrtBAIj.plot.1.mean.J.I * logBAj.plot.1 + mean_spei12 * sqrtBA.O.plot.1 + logdbhJ.IMall.plot.mean *sqrtBA.O.plot.1 + sqrtBAIj.plot.1.mean.J.I * sqrtBA.O.plot.1,family=nbinom1, data=dfplot),silent=T)
    
    try(Saving(rW.yO.sumF.mC,CODE),silent=T)
    try(print(summary(rW.yO.sumF.mC)),silent=T)
    rm(rW.yO.sumF.mC)
    
    # years offset both part, sum fixed both, mortality weight rate
    try(rW.yO.sumF.mR <- glmmTMB(sp.recruitment.ha.weight.2 ~ logdbhJ.IMall.plot.mean + 
                                   offset(log(yearsbetweensurveys)) + 
                                   sp.SUM.ha.weight2 + 
                                   logtreeNbrJ + Plotcat * mean_spei12 + Plotcat *  sqrtBA.O.plot.1 + Plotcat * logBAj.plot.1 + Plotcat * sqrtBAIj.plot.1.mean.J.I + 
                                   Plotcat * sp.mortality.ha.weight.rate +
                                   mean_spei12 * logBAj.plot.1 + logdbhJ.IMall.plot.mean * logBAj.plot.1 + sqrtBAIj.plot.1.mean.J.I * logBAj.plot.1 + mean_spei12 * sqrtBA.O.plot.1 + logdbhJ.IMall.plot.mean *sqrtBA.O.plot.1 + sqrtBAIj.plot.1.mean.J.I * sqrtBA.O.plot.1,
                                 zi=~logdbhJ.IMall.plot.mean + 
                                   offset(log(yearsbetweensurveys)) + 
                                   sp.SUM.ha.weight2 + 
                                   logtreeNbrJ + Plotcat * mean_spei12 + Plotcat *  sqrtBA.O.plot.1 + Plotcat * logBAj.plot.1 + Plotcat * sqrtBAIj.plot.1.mean.J.I + 
                                   Plotcat * sp.mortality.ha.weight.rate +
                                   mean_spei12 * logBAj.plot.1 + logdbhJ.IMall.plot.mean * logBAj.plot.1 + sqrtBAIj.plot.1.mean.J.I * logBAj.plot.1 + mean_spei12 * sqrtBA.O.plot.1 + logdbhJ.IMall.plot.mean *sqrtBA.O.plot.1 + sqrtBAIj.plot.1.mean.J.I * sqrtBA.O.plot.1,family=nbinom1, data=dfplot),silent=T)
    
    try(Saving(rW.yO.sumF.mR,CODE),silent=T)
    try(print(summary(rW.yO.sumF.mR)),silent=T)
    rm(rW.yO.sumF.mR)
    
    # years fixed both part, sum fixed both, mortality weight rate year
    try(rW.yO.sumF.mRY <- glmmTMB(sp.recruitment.ha.weight.2 ~ logdbhJ.IMall.plot.mean + 
                                    offset(log(yearsbetweensurveys)) + 
                                    sp.SUM.ha.weight2 + 
                                    logtreeNbrJ + Plotcat * mean_spei12 + Plotcat *  sqrtBA.O.plot.1 + Plotcat * logBAj.plot.1 + Plotcat * sqrtBAIj.plot.1.mean.J.I + 
                                    Plotcat * sp.mortality.ha.weight.rate.yr +
                                    mean_spei12 * logBAj.plot.1 + logdbhJ.IMall.plot.mean * logBAj.plot.1 + sqrtBAIj.plot.1.mean.J.I * logBAj.plot.1 + mean_spei12 * sqrtBA.O.plot.1 + logdbhJ.IMall.plot.mean *sqrtBA.O.plot.1 + sqrtBAIj.plot.1.mean.J.I * sqrtBA.O.plot.1,
                                  zi=~logdbhJ.IMall.plot.mean + 
                                    offset(log(yearsbetweensurveys)) + 
                                    sp.SUM.ha.weight2 + 
                                    logtreeNbrJ + Plotcat * mean_spei12 + Plotcat *  sqrtBA.O.plot.1 + Plotcat * logBAj.plot.1 + Plotcat * sqrtBAIj.plot.1.mean.J.I + 
                                    Plotcat * sp.mortality.ha.weight.rate.yr +
                                    mean_spei12 * logBAj.plot.1 + logdbhJ.IMall.plot.mean * logBAj.plot.1 + sqrtBAIj.plot.1.mean.J.I * logBAj.plot.1 + mean_spei12 * sqrtBA.O.plot.1 + logdbhJ.IMall.plot.mean *sqrtBA.O.plot.1 + sqrtBAIj.plot.1.mean.J.I * sqrtBA.O.plot.1,family=nbinom1, data=dfplot),silent=T)
    
    try(Saving(rW.yO.sumF.mRY,CODE),silent=T)
    try(print(summary(rW.yO.sumF.mRY)),silent=T)
    rm(rW.yO.sumF.mRY)
    
    # years fixed both part, sum offset both, mortality weight
    try(rW.yF.sumO.mC <- glmmTMB(sp.recruitment.ha.weight.2 ~ logdbhJ.IMall.plot.mean + 
                                   yearsbetweensurveys + 
                                   offset(log(sp.SUM.ha.weight2)) + 
                                   logtreeNbrJ + Plotcat * mean_spei12 + Plotcat *  sqrtBA.O.plot.1 + Plotcat * logBAj.plot.1 + Plotcat * sqrtBAIj.plot.1.mean.J.I + 
                                   Plotcat * sp.mortality.ha.weight +
                                   mean_spei12 * logBAj.plot.1 + logdbhJ.IMall.plot.mean * logBAj.plot.1 + sqrtBAIj.plot.1.mean.J.I * logBAj.plot.1 + mean_spei12 * sqrtBA.O.plot.1 + logdbhJ.IMall.plot.mean *sqrtBA.O.plot.1 + sqrtBAIj.plot.1.mean.J.I * sqrtBA.O.plot.1,
                                 zi=~logdbhJ.IMall.plot.mean + 
                                   yearsbetweensurveys + 
                                   offset(log(sp.SUM.ha.weight2)) + 
                                   logtreeNbrJ + Plotcat * mean_spei12 + Plotcat *  sqrtBA.O.plot.1 + Plotcat * logBAj.plot.1 + Plotcat * sqrtBAIj.plot.1.mean.J.I + 
                                   Plotcat * sp.mortality.ha.weight +
                                   mean_spei12 * logBAj.plot.1 + logdbhJ.IMall.plot.mean * logBAj.plot.1 + sqrtBAIj.plot.1.mean.J.I * logBAj.plot.1 + mean_spei12 * sqrtBA.O.plot.1 + logdbhJ.IMall.plot.mean *sqrtBA.O.plot.1 + sqrtBAIj.plot.1.mean.J.I * sqrtBA.O.plot.1,family=nbinom1, data=dfplot),silent=T)
    
    try(Saving(rW.yF.sumO.mC,CODE),silent=T)
    try(print(summary(rW.yF.sumO.mC)),silent=T)
    rm(rW.yF.sumO.mC)
    
    # years fixed both part, sum offset both, mortality weight rate
    try(rW.yF.sumO.mR <- glmmTMB(sp.recruitment.ha.weight.2 ~ logdbhJ.IMall.plot.mean + 
                                   yearsbetweensurveys + 
                                   offset(log(sp.SUM.ha.weight2)) +
                                   logtreeNbrJ + Plotcat * mean_spei12 + Plotcat *  sqrtBA.O.plot.1 + Plotcat * logBAj.plot.1 + Plotcat * sqrtBAIj.plot.1.mean.J.I + 
                                   Plotcat * sp.mortality.ha.weight.rate +
                                   mean_spei12 * logBAj.plot.1 + logdbhJ.IMall.plot.mean * logBAj.plot.1 + sqrtBAIj.plot.1.mean.J.I * logBAj.plot.1 + mean_spei12 * sqrtBA.O.plot.1 + logdbhJ.IMall.plot.mean *sqrtBA.O.plot.1 + sqrtBAIj.plot.1.mean.J.I * sqrtBA.O.plot.1,
                                 zi=~logdbhJ.IMall.plot.mean + 
                                   yearsbetweensurveys + 
                                   offset(log(sp.SUM.ha.weight2)) + 
                                   logtreeNbrJ + Plotcat * mean_spei12 + Plotcat *  sqrtBA.O.plot.1 + Plotcat * logBAj.plot.1 + Plotcat * sqrtBAIj.plot.1.mean.J.I + 
                                   Plotcat * sp.mortality.ha.weight.rate +
                                   mean_spei12 * logBAj.plot.1 + logdbhJ.IMall.plot.mean * logBAj.plot.1 + sqrtBAIj.plot.1.mean.J.I * logBAj.plot.1 + mean_spei12 * sqrtBA.O.plot.1 + logdbhJ.IMall.plot.mean *sqrtBA.O.plot.1 + sqrtBAIj.plot.1.mean.J.I * sqrtBA.O.plot.1,family=nbinom1, data=dfplot),silent=T)
    
    try(Saving(rW.yF.sumO.mR,CODE),silent=T)
    try(print(summary(rW.yF.sumO.mR)),silent=T)
    rm(rW.yF.sumO.mR)
    
    # years fixed both part, sum offset both, mortality weight rate year
    try(rW.yF.sumO.mRY <- glmmTMB(sp.recruitment.ha.weight.2 ~ logdbhJ.IMall.plot.mean + 
                                    yearsbetweensurveys + 
                                    offset(log(sp.SUM.ha.weight2)) + 
                                    logtreeNbrJ + Plotcat * mean_spei12 + Plotcat *  sqrtBA.O.plot.1 + Plotcat * logBAj.plot.1 + Plotcat * sqrtBAIj.plot.1.mean.J.I + 
                                    Plotcat * sp.mortality.ha.weight.rate.yr +
                                    mean_spei12 * logBAj.plot.1 + logdbhJ.IMall.plot.mean * logBAj.plot.1 + sqrtBAIj.plot.1.mean.J.I * logBAj.plot.1 + mean_spei12 * sqrtBA.O.plot.1 + logdbhJ.IMall.plot.mean *sqrtBA.O.plot.1 + sqrtBAIj.plot.1.mean.J.I * sqrtBA.O.plot.1,
                                  zi=~logdbhJ.IMall.plot.mean + 
                                    yearsbetweensurveys + 
                                    offset(log(sp.SUM.ha.weight2)) + 
                                    logtreeNbrJ + Plotcat * mean_spei12 + Plotcat *  sqrtBA.O.plot.1 + Plotcat * logBAj.plot.1 + Plotcat * sqrtBAIj.plot.1.mean.J.I + 
                                    Plotcat * sp.mortality.ha.weight.rate.yr +
                                    mean_spei12 * logBAj.plot.1 + logdbhJ.IMall.plot.mean * logBAj.plot.1 + sqrtBAIj.plot.1.mean.J.I * logBAj.plot.1 + mean_spei12 * sqrtBA.O.plot.1 + logdbhJ.IMall.plot.mean *sqrtBA.O.plot.1 + sqrtBAIj.plot.1.mean.J.I * sqrtBA.O.plot.1,family=nbinom1, data=dfplot),silent=T)
    
    try(Saving(rW.yF.sumO.mRY,CODE),silent=T)
    try(print(summary(rW.yF.sumO.mRY)),silent=T)
    rm(rW.yF.sumO.mRY)
    
    
    # years offset both part, sum offset both, mortality weight
    try(rW.yO.sumO.mC <- glmmTMB(sp.recruitment.ha.weight.2 ~ logdbhJ.IMall.plot.mean + 
                                   offset(log(yearsbetweensurveys)) + 
                                   offset(log(sp.SUM.ha.weight2)) + 
                                   logtreeNbrJ + Plotcat * mean_spei12 + Plotcat *  sqrtBA.O.plot.1 + Plotcat * logBAj.plot.1 + Plotcat * sqrtBAIj.plot.1.mean.J.I + 
                                   Plotcat * sp.mortality.ha.weight +
                                   mean_spei12 * logBAj.plot.1 + logdbhJ.IMall.plot.mean * logBAj.plot.1 + sqrtBAIj.plot.1.mean.J.I * logBAj.plot.1 + mean_spei12 * sqrtBA.O.plot.1 + logdbhJ.IMall.plot.mean *sqrtBA.O.plot.1 + sqrtBAIj.plot.1.mean.J.I * sqrtBA.O.plot.1,
                                 zi=~logdbhJ.IMall.plot.mean + 
                                   offset(log(yearsbetweensurveys)) + 
                                   offset(log(sp.SUM.ha.weight2)) + 
                                   logtreeNbrJ + Plotcat * mean_spei12 + Plotcat *  sqrtBA.O.plot.1 + Plotcat * logBAj.plot.1 + Plotcat * sqrtBAIj.plot.1.mean.J.I + 
                                   Plotcat * sp.mortality.ha.weight +
                                   mean_spei12 * logBAj.plot.1 + logdbhJ.IMall.plot.mean * logBAj.plot.1 + sqrtBAIj.plot.1.mean.J.I * logBAj.plot.1 + mean_spei12 * sqrtBA.O.plot.1 + logdbhJ.IMall.plot.mean *sqrtBA.O.plot.1 + sqrtBAIj.plot.1.mean.J.I * sqrtBA.O.plot.1,family=nbinom1, data=dfplot),silent=T)
    
    try(Saving(rW.yO.sumO.mC,CODE),silent=T)
    try(print(summary(rW.yO.sumO.mC)),silent=T)
    rm(rW.yO.sumO.mC)
    
    # years offset both part, sum offset both, mortality weight rate
    try(rW.yO.sumO.mR <- glmmTMB(sp.recruitment.ha.weight.2 ~ logdbhJ.IMall.plot.mean + 
                                   offset(log(yearsbetweensurveys)) + 
                                   offset(log(sp.SUM.ha.weight2)) +
                                   logtreeNbrJ + Plotcat * mean_spei12 + Plotcat *  sqrtBA.O.plot.1 + Plotcat * logBAj.plot.1 + Plotcat * sqrtBAIj.plot.1.mean.J.I + 
                                   Plotcat * sp.mortality.ha.weight.rate +
                                   mean_spei12 * logBAj.plot.1 + logdbhJ.IMall.plot.mean * logBAj.plot.1 + sqrtBAIj.plot.1.mean.J.I * logBAj.plot.1 + mean_spei12 * sqrtBA.O.plot.1 + logdbhJ.IMall.plot.mean *sqrtBA.O.plot.1 + sqrtBAIj.plot.1.mean.J.I * sqrtBA.O.plot.1,
                                 zi=~logdbhJ.IMall.plot.mean + 
                                   offset(log(yearsbetweensurveys)) + 
                                   offset(log(sp.SUM.ha.weight2)) + 
                                   logtreeNbrJ + Plotcat * mean_spei12 + Plotcat *  sqrtBA.O.plot.1 + Plotcat * logBAj.plot.1 + Plotcat * sqrtBAIj.plot.1.mean.J.I + 
                                   Plotcat * sp.mortality.ha.weight.rate +
                                   mean_spei12 * logBAj.plot.1 + logdbhJ.IMall.plot.mean * logBAj.plot.1 + sqrtBAIj.plot.1.mean.J.I * logBAj.plot.1 + mean_spei12 * sqrtBA.O.plot.1 + logdbhJ.IMall.plot.mean *sqrtBA.O.plot.1 + sqrtBAIj.plot.1.mean.J.I * sqrtBA.O.plot.1,family=nbinom1, data=dfplot),silent=T)
    
    try(Saving(rW.yO.sumO.mR,CODE),silent=T)
    try(print(summary(rW.yO.sumO.mR)),silent=T)
    rm(rW.yO.sumO.mR)
    

    
  }
    
  if (CODE%in%Allcode[c(5,8,17)]){
    
      # years fixed both part, sum fixed both, mortality weight 
      try(rW.yF.sumF.mC <- glmmTMB(sp.recruitment.ha.weight.2 ~ logdbhJ.IMall.plot.mean + 
                                     yearsbetweensurveys + 
                                     sp.SUM.ha.weight2 + 
                                     (1 | country) + logtreeNbrJ + Plotcat * mean_spei12 + Plotcat *  sqrtBA.O.plot.1 + Plotcat * logBAj.plot.1 + Plotcat * sqrtBAIj.plot.1.mean.J.I + 
                                     sp.mortality.ha.weight +
                                     mean_spei12 * logBAj.plot.1 + logdbhJ.IMall.plot.mean * logBAj.plot.1 + sqrtBAIj.plot.1.mean.J.I * logBAj.plot.1 + mean_spei12 * sqrtBA.O.plot.1 + logdbhJ.IMall.plot.mean *sqrtBA.O.plot.1 + sqrtBAIj.plot.1.mean.J.I * sqrtBA.O.plot.1,
                                   zi=~.,family=nbinom1, data=dfplot),silent=T)
      
      try(Saving(rW.yF.sumF.mC,CODE),silent=T)
      try(print(summary(rW.yF.sumF.mC)),silent=T)
      rm(rW.yF.sumF.mC)
      
      # years fixed both part, sum fixed both, mortality weight rate
      try(rW.yF.sumF.mR <- glmmTMB(sp.recruitment.ha.weight.2 ~ logdbhJ.IMall.plot.mean + 
                                     yearsbetweensurveys + 
                                     sp.SUM.ha.weight2 + 
                                     (1 | country) + logtreeNbrJ + Plotcat * mean_spei12 + Plotcat *  sqrtBA.O.plot.1 + Plotcat * logBAj.plot.1 + Plotcat * sqrtBAIj.plot.1.mean.J.I + 
                                     sp.mortality.ha.weight.rate +
                                     mean_spei12 * logBAj.plot.1 + logdbhJ.IMall.plot.mean * logBAj.plot.1 + sqrtBAIj.plot.1.mean.J.I * logBAj.plot.1 + mean_spei12 * sqrtBA.O.plot.1 + logdbhJ.IMall.plot.mean *sqrtBA.O.plot.1 + sqrtBAIj.plot.1.mean.J.I * sqrtBA.O.plot.1,
                                   zi=~.,family=nbinom1, data=dfplot),silent=T)
      try(Saving(rW.yF.sumF.mR,CODE),silent=T)
      try(print(summary(rW.yF.sumF.mR)),silent=T)
      rm(rW.yF.sumF.mR)
      
      # years fixed both part, sum fixed both, mortality weight rate year
      try(rW.yF.sumF.mRY <- glmmTMB(sp.recruitment.ha.weight.2 ~ logdbhJ.IMall.plot.mean + 
                                      yearsbetweensurveys + 
                                      sp.SUM.ha.weight2 + 
                                      (1 | country) + logtreeNbrJ + Plotcat * mean_spei12 + Plotcat *  sqrtBA.O.plot.1 + Plotcat * logBAj.plot.1 + Plotcat * sqrtBAIj.plot.1.mean.J.I + 
                                      sp.mortality.ha.weight.rate.yr +
                                      mean_spei12 * logBAj.plot.1 + logdbhJ.IMall.plot.mean * logBAj.plot.1 + sqrtBAIj.plot.1.mean.J.I * logBAj.plot.1 + mean_spei12 * sqrtBA.O.plot.1 + logdbhJ.IMall.plot.mean *sqrtBA.O.plot.1 + sqrtBAIj.plot.1.mean.J.I * sqrtBA.O.plot.1,
                                    zi=~.,family=nbinom1, data=dfplot),silent=T)
      try(Saving(rW.yF.sumF.mRY,CODE),silent=T)
      try(print(summary(rW.yF.sumF.mRY)),silent=T)
      rm(rW.yF.sumF.mRY)
      
      # years offset both part, sum fixed both, mortality weight
      try(rW.yO.sumF.mC <- glmmTMB(sp.recruitment.ha.weight.2 ~ logdbhJ.IMall.plot.mean + 
                                     offset(log(yearsbetweensurveys)) + 
                                     sp.SUM.ha.weight2 + 
                                     (1 | country) + logtreeNbrJ + Plotcat * mean_spei12 + Plotcat *  sqrtBA.O.plot.1 + Plotcat * logBAj.plot.1 + Plotcat * sqrtBAIj.plot.1.mean.J.I + 
                                     sp.mortality.ha.weight +
                                     mean_spei12 * logBAj.plot.1 + logdbhJ.IMall.plot.mean * logBAj.plot.1 + sqrtBAIj.plot.1.mean.J.I * logBAj.plot.1 + mean_spei12 * sqrtBA.O.plot.1 + logdbhJ.IMall.plot.mean *sqrtBA.O.plot.1 + sqrtBAIj.plot.1.mean.J.I * sqrtBA.O.plot.1,
                                   zi=~logdbhJ.IMall.plot.mean + 
                                     offset(log(yearsbetweensurveys)) + 
                                     sp.SUM.ha.weight2 + 
                                     (1 | country) + logtreeNbrJ + Plotcat * mean_spei12 + Plotcat *  sqrtBA.O.plot.1 + Plotcat * logBAj.plot.1 + Plotcat * sqrtBAIj.plot.1.mean.J.I + 
                                     sp.mortality.ha.weight +
                                     mean_spei12 * logBAj.plot.1 + logdbhJ.IMall.plot.mean * logBAj.plot.1 + sqrtBAIj.plot.1.mean.J.I * logBAj.plot.1 + mean_spei12 * sqrtBA.O.plot.1 + logdbhJ.IMall.plot.mean *sqrtBA.O.plot.1 + sqrtBAIj.plot.1.mean.J.I * sqrtBA.O.plot.1,family=nbinom1, data=dfplot),silent=T)
      
      try(Saving(rW.yO.sumF.mC,CODE),silent=T)
      try(print(summary(rW.yO.sumF.mC)),silent=T)
      rm(rW.yO.sumF.mC)
      
      # years offset both part, sum fixed both, mortality weight rate
      try(rW.yO.sumF.mR <- glmmTMB(sp.recruitment.ha.weight.2 ~ logdbhJ.IMall.plot.mean + 
                                     offset(log(yearsbetweensurveys)) + 
                                     sp.SUM.ha.weight2 + 
                                     (1 | country) + logtreeNbrJ + Plotcat * mean_spei12 + Plotcat *  sqrtBA.O.plot.1 + Plotcat * logBAj.plot.1 + Plotcat * sqrtBAIj.plot.1.mean.J.I + 
                                     sp.mortality.ha.weight.rate +
                                     mean_spei12 * logBAj.plot.1 + logdbhJ.IMall.plot.mean * logBAj.plot.1 + sqrtBAIj.plot.1.mean.J.I * logBAj.plot.1 + mean_spei12 * sqrtBA.O.plot.1 + logdbhJ.IMall.plot.mean *sqrtBA.O.plot.1 + sqrtBAIj.plot.1.mean.J.I * sqrtBA.O.plot.1,
                                   zi=~logdbhJ.IMall.plot.mean + 
                                     offset(log(yearsbetweensurveys)) + 
                                     sp.SUM.ha.weight2 + 
                                     (1 | country) + logtreeNbrJ + Plotcat * mean_spei12 + Plotcat *  sqrtBA.O.plot.1 + Plotcat * logBAj.plot.1 + Plotcat * sqrtBAIj.plot.1.mean.J.I + 
                                     sp.mortality.ha.weight.rate +
                                     mean_spei12 * logBAj.plot.1 + logdbhJ.IMall.plot.mean * logBAj.plot.1 + sqrtBAIj.plot.1.mean.J.I * logBAj.plot.1 + mean_spei12 * sqrtBA.O.plot.1 + logdbhJ.IMall.plot.mean *sqrtBA.O.plot.1 + sqrtBAIj.plot.1.mean.J.I * sqrtBA.O.plot.1,family=nbinom1, data=dfplot),silent=T)
      
      try(Saving(rW.yO.sumF.mR,CODE),silent=T)
      try(print(summary(rW.yO.sumF.mR)),silent=T)
      rm(rW.yO.sumF.mR)
      
      # years fixed both part, sum fixed both, mortality weight rate year
      try(rW.yO.sumF.mRY <- glmmTMB(sp.recruitment.ha.weight.2 ~ logdbhJ.IMall.plot.mean + 
                                      offset(log(yearsbetweensurveys)) + 
                                      sp.SUM.ha.weight2 + 
                                      (1 | country) + logtreeNbrJ + Plotcat * mean_spei12 + Plotcat *  sqrtBA.O.plot.1 + Plotcat * logBAj.plot.1 + Plotcat * sqrtBAIj.plot.1.mean.J.I + 
                                      sp.mortality.ha.weight.rate.yr +
                                      mean_spei12 * logBAj.plot.1 + logdbhJ.IMall.plot.mean * logBAj.plot.1 + sqrtBAIj.plot.1.mean.J.I * logBAj.plot.1 + mean_spei12 * sqrtBA.O.plot.1 + logdbhJ.IMall.plot.mean *sqrtBA.O.plot.1 + sqrtBAIj.plot.1.mean.J.I * sqrtBA.O.plot.1,
                                    zi=~logdbhJ.IMall.plot.mean + 
                                      offset(log(yearsbetweensurveys)) + 
                                      sp.SUM.ha.weight2 + 
                                      (1 | country) + logtreeNbrJ + Plotcat * mean_spei12 + Plotcat *  sqrtBA.O.plot.1 + Plotcat * logBAj.plot.1 + Plotcat * sqrtBAIj.plot.1.mean.J.I + 
                                      sp.mortality.ha.weight.rate.yr +
                                      mean_spei12 * logBAj.plot.1 + logdbhJ.IMall.plot.mean * logBAj.plot.1 + sqrtBAIj.plot.1.mean.J.I * logBAj.plot.1 + mean_spei12 * sqrtBA.O.plot.1 + logdbhJ.IMall.plot.mean *sqrtBA.O.plot.1 + sqrtBAIj.plot.1.mean.J.I * sqrtBA.O.plot.1,family=nbinom1, data=dfplot),silent=T)
      
      try(Saving(rW.yO.sumF.mRY,CODE),silent=T)
      try(print(summary(rW.yO.sumF.mRY)),silent=T)
      rm(rW.yO.sumF.mRY)
      
      # years fixed both part, sum offset both, mortality weight
      try(rW.yF.sumO.mC <- glmmTMB(sp.recruitment.ha.weight.2 ~ logdbhJ.IMall.plot.mean + 
                                     yearsbetweensurveys + 
                                     offset(log(sp.SUM.ha.weight2)) + 
                                     (1 | country) + logtreeNbrJ + Plotcat * mean_spei12 + Plotcat *  sqrtBA.O.plot.1 + Plotcat * logBAj.plot.1 + Plotcat * sqrtBAIj.plot.1.mean.J.I + 
                                     sp.mortality.ha.weight +
                                     mean_spei12 * logBAj.plot.1 + logdbhJ.IMall.plot.mean * logBAj.plot.1 + sqrtBAIj.plot.1.mean.J.I * logBAj.plot.1 + mean_spei12 * sqrtBA.O.plot.1 + logdbhJ.IMall.plot.mean *sqrtBA.O.plot.1 + sqrtBAIj.plot.1.mean.J.I * sqrtBA.O.plot.1,
                                   zi=~logdbhJ.IMall.plot.mean + 
                                     yearsbetweensurveys + 
                                     offset(log(sp.SUM.ha.weight2)) + 
                                     (1 | country) + logtreeNbrJ + Plotcat * mean_spei12 + Plotcat *  sqrtBA.O.plot.1 + Plotcat * logBAj.plot.1 + Plotcat * sqrtBAIj.plot.1.mean.J.I + 
                                     sp.mortality.ha.weight +
                                     mean_spei12 * logBAj.plot.1 + logdbhJ.IMall.plot.mean * logBAj.plot.1 + sqrtBAIj.plot.1.mean.J.I * logBAj.plot.1 + mean_spei12 * sqrtBA.O.plot.1 + logdbhJ.IMall.plot.mean *sqrtBA.O.plot.1 + sqrtBAIj.plot.1.mean.J.I * sqrtBA.O.plot.1,family=nbinom1, data=dfplot),silent=T)
      
      try(Saving(rW.yF.sumO.mC,CODE),silent=T)
      try(print(summary(rW.yF.sumO.mC)),silent=T)
      rm(rW.yF.sumO.mC)
      
      # years fixed both part, sum offset both, mortality weight rate
      try(rW.yF.sumO.mR <- glmmTMB(sp.recruitment.ha.weight.2 ~ logdbhJ.IMall.plot.mean + 
                                     yearsbetweensurveys + 
                                     offset(log(sp.SUM.ha.weight2)) +
                                     (1 | country) + logtreeNbrJ + Plotcat * mean_spei12 + Plotcat *  sqrtBA.O.plot.1 + Plotcat * logBAj.plot.1 + Plotcat * sqrtBAIj.plot.1.mean.J.I + 
                                     sp.mortality.ha.weight.rate +
                                     mean_spei12 * logBAj.plot.1 + logdbhJ.IMall.plot.mean * logBAj.plot.1 + sqrtBAIj.plot.1.mean.J.I * logBAj.plot.1 + mean_spei12 * sqrtBA.O.plot.1 + logdbhJ.IMall.plot.mean *sqrtBA.O.plot.1 + sqrtBAIj.plot.1.mean.J.I * sqrtBA.O.plot.1,
                                   zi=~logdbhJ.IMall.plot.mean + 
                                     yearsbetweensurveys + 
                                     offset(log(sp.SUM.ha.weight2)) + 
                                     (1 | country) + logtreeNbrJ + Plotcat * mean_spei12 + Plotcat *  sqrtBA.O.plot.1 + Plotcat * logBAj.plot.1 + Plotcat * sqrtBAIj.plot.1.mean.J.I + 
                                     sp.mortality.ha.weight.rate +
                                     mean_spei12 * logBAj.plot.1 + logdbhJ.IMall.plot.mean * logBAj.plot.1 + sqrtBAIj.plot.1.mean.J.I * logBAj.plot.1 + mean_spei12 * sqrtBA.O.plot.1 + logdbhJ.IMall.plot.mean *sqrtBA.O.plot.1 + sqrtBAIj.plot.1.mean.J.I * sqrtBA.O.plot.1,family=nbinom1, data=dfplot),silent=T)
      
      try(Saving(rW.yF.sumO.mR,CODE),silent=T)
      try(print(summary(rW.yF.sumO.mR)),silent=T)
      rm(rW.yF.sumO.mR)
      
      # years fixed both part, sum offset both, mortality weight rate year
      try(rW.yF.sumO.mRY <- glmmTMB(sp.recruitment.ha.weight.2 ~ logdbhJ.IMall.plot.mean + 
                                      yearsbetweensurveys + 
                                      offset(log(sp.SUM.ha.weight2)) + 
                                      (1 | country) + logtreeNbrJ + Plotcat * mean_spei12 + Plotcat *  sqrtBA.O.plot.1 + Plotcat * logBAj.plot.1 + Plotcat * sqrtBAIj.plot.1.mean.J.I + 
                                      sp.mortality.ha.weight.rate.yr +
                                      mean_spei12 * logBAj.plot.1 + logdbhJ.IMall.plot.mean * logBAj.plot.1 + sqrtBAIj.plot.1.mean.J.I * logBAj.plot.1 + mean_spei12 * sqrtBA.O.plot.1 + logdbhJ.IMall.plot.mean *sqrtBA.O.plot.1 + sqrtBAIj.plot.1.mean.J.I * sqrtBA.O.plot.1,
                                    zi=~logdbhJ.IMall.plot.mean + 
                                      yearsbetweensurveys + 
                                      offset(log(sp.SUM.ha.weight2)) + 
                                      (1 | country) + logtreeNbrJ + Plotcat * mean_spei12 + Plotcat *  sqrtBA.O.plot.1 + Plotcat * logBAj.plot.1 + Plotcat * sqrtBAIj.plot.1.mean.J.I + 
                                      sp.mortality.ha.weight.rate.yr +
                                      mean_spei12 * logBAj.plot.1 + logdbhJ.IMall.plot.mean * logBAj.plot.1 + sqrtBAIj.plot.1.mean.J.I * logBAj.plot.1 + mean_spei12 * sqrtBA.O.plot.1 + logdbhJ.IMall.plot.mean *sqrtBA.O.plot.1 + sqrtBAIj.plot.1.mean.J.I * sqrtBA.O.plot.1,family=nbinom1, data=dfplot),silent=T)
      
      try(Saving(rW.yF.sumO.mRY,CODE),silent=T)
      try(print(summary(rW.yF.sumO.mRY)),silent=T)
      rm(rW.yF.sumO.mRY)
      
      
      # years offset both part, sum offset both, mortality weight
      try(rW.yO.sumO.mC <- glmmTMB(sp.recruitment.ha.weight.2 ~ logdbhJ.IMall.plot.mean + 
                                     offset(log(yearsbetweensurveys)) + 
                                     offset(log(sp.SUM.ha.weight2)) + 
                                     (1 | country) + logtreeNbrJ + Plotcat * mean_spei12 + Plotcat *  sqrtBA.O.plot.1 + Plotcat * logBAj.plot.1 + Plotcat * sqrtBAIj.plot.1.mean.J.I + 
                                     sp.mortality.ha.weight +
                                     mean_spei12 * logBAj.plot.1 + logdbhJ.IMall.plot.mean * logBAj.plot.1 + sqrtBAIj.plot.1.mean.J.I * logBAj.plot.1 + mean_spei12 * sqrtBA.O.plot.1 + logdbhJ.IMall.plot.mean *sqrtBA.O.plot.1 + sqrtBAIj.plot.1.mean.J.I * sqrtBA.O.plot.1,
                                   zi=~logdbhJ.IMall.plot.mean + 
                                     offset(log(yearsbetweensurveys)) + 
                                     offset(log(sp.SUM.ha.weight2)) + 
                                     (1 | country) + logtreeNbrJ + Plotcat * mean_spei12 + Plotcat *  sqrtBA.O.plot.1 + Plotcat * logBAj.plot.1 + Plotcat * sqrtBAIj.plot.1.mean.J.I + 
                                     sp.mortality.ha.weight +
                                     mean_spei12 * logBAj.plot.1 + logdbhJ.IMall.plot.mean * logBAj.plot.1 + sqrtBAIj.plot.1.mean.J.I * logBAj.plot.1 + mean_spei12 * sqrtBA.O.plot.1 + logdbhJ.IMall.plot.mean *sqrtBA.O.plot.1 + sqrtBAIj.plot.1.mean.J.I * sqrtBA.O.plot.1,family=nbinom1, data=dfplot),silent=T)
      
      try(Saving(rW.yO.sumO.mC,CODE),silent=T)
      try(print(summary(rW.yO.sumO.mC)),silent=T)
      rm(rW.yO.sumO.mC)
      
      # years offset both part, sum offset both, mortality weight rate
      try(rW.yO.sumO.mR <- glmmTMB(sp.recruitment.ha.weight.2 ~ logdbhJ.IMall.plot.mean + 
                                     offset(log(yearsbetweensurveys)) + 
                                     offset(log(sp.SUM.ha.weight2)) +
                                     (1 | country) + logtreeNbrJ + Plotcat * mean_spei12 + Plotcat *  sqrtBA.O.plot.1 + Plotcat * logBAj.plot.1 + Plotcat * sqrtBAIj.plot.1.mean.J.I + 
                                     sp.mortality.ha.weight.rate +
                                     mean_spei12 * logBAj.plot.1 + logdbhJ.IMall.plot.mean * logBAj.plot.1 + sqrtBAIj.plot.1.mean.J.I * logBAj.plot.1 + mean_spei12 * sqrtBA.O.plot.1 + logdbhJ.IMall.plot.mean *sqrtBA.O.plot.1 + sqrtBAIj.plot.1.mean.J.I * sqrtBA.O.plot.1,
                                   zi=~logdbhJ.IMall.plot.mean + 
                                     offset(log(yearsbetweensurveys)) + 
                                     offset(log(sp.SUM.ha.weight2)) + 
                                     (1 | country) + logtreeNbrJ + Plotcat * mean_spei12 + Plotcat *  sqrtBA.O.plot.1 + Plotcat * logBAj.plot.1 + Plotcat * sqrtBAIj.plot.1.mean.J.I + 
                                     sp.mortality.ha.weight.rate +
                                     mean_spei12 * logBAj.plot.1 + logdbhJ.IMall.plot.mean * logBAj.plot.1 + sqrtBAIj.plot.1.mean.J.I * logBAj.plot.1 + mean_spei12 * sqrtBA.O.plot.1 + logdbhJ.IMall.plot.mean *sqrtBA.O.plot.1 + sqrtBAIj.plot.1.mean.J.I * sqrtBA.O.plot.1,family=nbinom1, data=dfplot),silent=T)
      
      try(Saving(rW.yO.sumO.mR,CODE),silent=T)
      try(print(summary(rW.yO.sumO.mR)),silent=T)
      rm(rW.yO.sumO.mR)
      
      
    
      
      
  }
  if (CODE%in%Allcode[c(11)]){
          #if (i%in%c(11)){
        
        # years fixed both part, sum fixed both, mortality weight 
        try(rW.yF.sumF.mC <- glmmTMB(sp.recruitment.ha.weight.2 ~ logdbhJ.IMall.plot.mean + 
                                       yearsbetweensurveys + 
                                       sp.SUM.ha.weight2 + 
                                       logtreeNbrJ + Plotcat * mean_spei12 + Plotcat *  sqrtBA.O.plot.1 + Plotcat * logBAj.plot.1 + Plotcat * sqrtBAIj.plot.1.mean.J.I + 
                                       sp.mortality.ha.weight +
                                       mean_spei12 * logBAj.plot.1 + logdbhJ.IMall.plot.mean * logBAj.plot.1 + sqrtBAIj.plot.1.mean.J.I * logBAj.plot.1 + mean_spei12 * sqrtBA.O.plot.1 + logdbhJ.IMall.plot.mean *sqrtBA.O.plot.1 + sqrtBAIj.plot.1.mean.J.I * sqrtBA.O.plot.1,
                                     zi=~.,family=nbinom1, data=dfplot),silent=T)
        
        try(Saving(rW.yF.sumF.mC,CODE),silent=T)
        try(print(summary(rW.yF.sumF.mC)),silent=T)
        rm(rW.yF.sumF.mC)
        
        # years fixed both part, sum fixed both, mortality weight rate
        try(rW.yF.sumF.mR <- glmmTMB(sp.recruitment.ha.weight.2 ~ logdbhJ.IMall.plot.mean + 
                                       yearsbetweensurveys + 
                                       sp.SUM.ha.weight2 + 
                                       logtreeNbrJ + Plotcat * mean_spei12 + Plotcat *  sqrtBA.O.plot.1 + Plotcat * logBAj.plot.1 + Plotcat * sqrtBAIj.plot.1.mean.J.I + 
                                       sp.mortality.ha.weight.rate +
                                       mean_spei12 * logBAj.plot.1 + logdbhJ.IMall.plot.mean * logBAj.plot.1 + sqrtBAIj.plot.1.mean.J.I * logBAj.plot.1 + mean_spei12 * sqrtBA.O.plot.1 + logdbhJ.IMall.plot.mean *sqrtBA.O.plot.1 + sqrtBAIj.plot.1.mean.J.I * sqrtBA.O.plot.1,
                                     zi=~.,family=nbinom1, data=dfplot),silent=T)
        try(Saving(rW.yF.sumF.mR,CODE),silent=T)
        try(print(summary(rW.yF.sumF.mR)),silent=T)
        rm(rW.yF.sumF.mR)
        
        # years fixed both part, sum fixed both, mortality weight rate year
        try(rW.yF.sumF.mRY <- glmmTMB(sp.recruitment.ha.weight.2 ~ logdbhJ.IMall.plot.mean + 
                                        yearsbetweensurveys + 
                                        sp.SUM.ha.weight2 + 
                                        logtreeNbrJ + Plotcat * mean_spei12 + Plotcat *  sqrtBA.O.plot.1 + Plotcat * logBAj.plot.1 + Plotcat * sqrtBAIj.plot.1.mean.J.I + 
                                        sp.mortality.ha.weight.rate.yr +
                                        mean_spei12 * logBAj.plot.1 + logdbhJ.IMall.plot.mean * logBAj.plot.1 + sqrtBAIj.plot.1.mean.J.I * logBAj.plot.1 + mean_spei12 * sqrtBA.O.plot.1 + logdbhJ.IMall.plot.mean *sqrtBA.O.plot.1 + sqrtBAIj.plot.1.mean.J.I * sqrtBA.O.plot.1,
                                      zi=~.,family=nbinom1, data=dfplot),silent=T)
        try(Saving(rW.yF.sumF.mRY,CODE),silent=T)
        try(print(summary(rW.yF.sumF.mRY)),silent=T)
        rm(rW.yF.sumF.mRY)
        
        # years offset both part, sum fixed both, mortality weight
        try(rW.yO.sumF.mC <- glmmTMB(sp.recruitment.ha.weight.2 ~ logdbhJ.IMall.plot.mean + 
                                       offset(log(yearsbetweensurveys)) + 
                                       sp.SUM.ha.weight2 + 
                                       logtreeNbrJ + Plotcat * mean_spei12 + Plotcat *  sqrtBA.O.plot.1 + Plotcat * logBAj.plot.1 + Plotcat * sqrtBAIj.plot.1.mean.J.I + 
                                       sp.mortality.ha.weight +
                                       mean_spei12 * logBAj.plot.1 + logdbhJ.IMall.plot.mean * logBAj.plot.1 + sqrtBAIj.plot.1.mean.J.I * logBAj.plot.1 + mean_spei12 * sqrtBA.O.plot.1 + logdbhJ.IMall.plot.mean *sqrtBA.O.plot.1 + sqrtBAIj.plot.1.mean.J.I * sqrtBA.O.plot.1,
                                     zi=~logdbhJ.IMall.plot.mean + 
                                       offset(log(yearsbetweensurveys)) + 
                                       sp.SUM.ha.weight2 + 
                                       logtreeNbrJ + Plotcat * mean_spei12 + Plotcat *  sqrtBA.O.plot.1 + Plotcat * logBAj.plot.1 + Plotcat * sqrtBAIj.plot.1.mean.J.I + 
                                       sp.mortality.ha.weight +
                                       mean_spei12 * logBAj.plot.1 + logdbhJ.IMall.plot.mean * logBAj.plot.1 + sqrtBAIj.plot.1.mean.J.I * logBAj.plot.1 + mean_spei12 * sqrtBA.O.plot.1 + logdbhJ.IMall.plot.mean *sqrtBA.O.plot.1 + sqrtBAIj.plot.1.mean.J.I * sqrtBA.O.plot.1,family=nbinom1, data=dfplot),silent=T)
        
        try(Saving(rW.yO.sumF.mC,CODE),silent=T)
        try(print(summary(rW.yO.sumF.mC)),silent=T)
        rm(rW.yO.sumF.mC)
        
        # years offset both part, sum fixed both, mortality weight rate
        try(rW.yO.sumF.mR <- glmmTMB(sp.recruitment.ha.weight.2 ~ logdbhJ.IMall.plot.mean + 
                                       offset(log(yearsbetweensurveys)) + 
                                       sp.SUM.ha.weight2 + 
                                       logtreeNbrJ + Plotcat * mean_spei12 + Plotcat *  sqrtBA.O.plot.1 + Plotcat * logBAj.plot.1 + Plotcat * sqrtBAIj.plot.1.mean.J.I + 
                                       sp.mortality.ha.weight.rate +
                                       mean_spei12 * logBAj.plot.1 + logdbhJ.IMall.plot.mean * logBAj.plot.1 + sqrtBAIj.plot.1.mean.J.I * logBAj.plot.1 + mean_spei12 * sqrtBA.O.plot.1 + logdbhJ.IMall.plot.mean *sqrtBA.O.plot.1 + sqrtBAIj.plot.1.mean.J.I * sqrtBA.O.plot.1,
                                     zi=~logdbhJ.IMall.plot.mean + 
                                       offset(log(yearsbetweensurveys)) + 
                                       sp.SUM.ha.weight2 + 
                                       logtreeNbrJ + Plotcat * mean_spei12 + Plotcat *  sqrtBA.O.plot.1 + Plotcat * logBAj.plot.1 + Plotcat * sqrtBAIj.plot.1.mean.J.I + 
                                       sp.mortality.ha.weight.rate +
                                       mean_spei12 * logBAj.plot.1 + logdbhJ.IMall.plot.mean * logBAj.plot.1 + sqrtBAIj.plot.1.mean.J.I * logBAj.plot.1 + mean_spei12 * sqrtBA.O.plot.1 + logdbhJ.IMall.plot.mean *sqrtBA.O.plot.1 + sqrtBAIj.plot.1.mean.J.I * sqrtBA.O.plot.1,family=nbinom1, data=dfplot),silent=T)
        
        try(Saving(rW.yO.sumF.mR,CODE),silent=T)
        try(print(summary(rW.yO.sumF.mR)),silent=T)
        rm(rW.yO.sumF.mR)
        
        # years fixed both part, sum fixed both, mortality weight rate year
        try(rW.yO.sumF.mRY <- glmmTMB(sp.recruitment.ha.weight.2 ~ logdbhJ.IMall.plot.mean + 
                                        offset(log(yearsbetweensurveys)) + 
                                        sp.SUM.ha.weight2 + 
                                        logtreeNbrJ + Plotcat * mean_spei12 + Plotcat *  sqrtBA.O.plot.1 + Plotcat * logBAj.plot.1 + Plotcat * sqrtBAIj.plot.1.mean.J.I + 
                                        sp.mortality.ha.weight.rate.yr +
                                        mean_spei12 * logBAj.plot.1 + logdbhJ.IMall.plot.mean * logBAj.plot.1 + sqrtBAIj.plot.1.mean.J.I * logBAj.plot.1 + mean_spei12 * sqrtBA.O.plot.1 + logdbhJ.IMall.plot.mean *sqrtBA.O.plot.1 + sqrtBAIj.plot.1.mean.J.I * sqrtBA.O.plot.1,
                                      zi=~logdbhJ.IMall.plot.mean + 
                                        offset(log(yearsbetweensurveys)) + 
                                        sp.SUM.ha.weight2 + 
                                        logtreeNbrJ + Plotcat * mean_spei12 + Plotcat *  sqrtBA.O.plot.1 + Plotcat * logBAj.plot.1 + Plotcat * sqrtBAIj.plot.1.mean.J.I + 
                                        sp.mortality.ha.weight.rate.yr +
                                        mean_spei12 * logBAj.plot.1 + logdbhJ.IMall.plot.mean * logBAj.plot.1 + sqrtBAIj.plot.1.mean.J.I * logBAj.plot.1 + mean_spei12 * sqrtBA.O.plot.1 + logdbhJ.IMall.plot.mean *sqrtBA.O.plot.1 + sqrtBAIj.plot.1.mean.J.I * sqrtBA.O.plot.1,family=nbinom1, data=dfplot),silent=T)
        
        try(Saving(rW.yO.sumF.mRY,CODE),silent=T)
        try(print(summary(rW.yO.sumF.mRY)),silent=T)
        rm(rW.yO.sumF.mRY)
        
        # years fixed both part, sum offset both, mortality weight
        try(rW.yF.sumO.mC <- glmmTMB(sp.recruitment.ha.weight.2 ~ logdbhJ.IMall.plot.mean + 
                                       yearsbetweensurveys + 
                                       offset(log(sp.SUM.ha.weight2)) + 
                                       logtreeNbrJ + Plotcat * mean_spei12 + Plotcat *  sqrtBA.O.plot.1 + Plotcat * logBAj.plot.1 + Plotcat * sqrtBAIj.plot.1.mean.J.I + 
                                       sp.mortality.ha.weight +
                                       mean_spei12 * logBAj.plot.1 + logdbhJ.IMall.plot.mean * logBAj.plot.1 + sqrtBAIj.plot.1.mean.J.I * logBAj.plot.1 + mean_spei12 * sqrtBA.O.plot.1 + logdbhJ.IMall.plot.mean *sqrtBA.O.plot.1 + sqrtBAIj.plot.1.mean.J.I * sqrtBA.O.plot.1,
                                     zi=~logdbhJ.IMall.plot.mean + 
                                       yearsbetweensurveys + 
                                       offset(log(sp.SUM.ha.weight2)) + 
                                       logtreeNbrJ + Plotcat * mean_spei12 + Plotcat *  sqrtBA.O.plot.1 + Plotcat * logBAj.plot.1 + Plotcat * sqrtBAIj.plot.1.mean.J.I + 
                                       sp.mortality.ha.weight +
                                       mean_spei12 * logBAj.plot.1 + logdbhJ.IMall.plot.mean * logBAj.plot.1 + sqrtBAIj.plot.1.mean.J.I * logBAj.plot.1 + mean_spei12 * sqrtBA.O.plot.1 + logdbhJ.IMall.plot.mean *sqrtBA.O.plot.1 + sqrtBAIj.plot.1.mean.J.I * sqrtBA.O.plot.1,family=nbinom1, data=dfplot),silent=T)
        
        try(Saving(rW.yF.sumO.mC,CODE),silent=T)
        try(print(summary(rW.yF.sumO.mC)),silent=T)
        rm(rW.yF.sumO.mC)
        
        # years fixed both part, sum offset both, mortality weight rate
        try(rW.yF.sumO.mR <- glmmTMB(sp.recruitment.ha.weight.2 ~ logdbhJ.IMall.plot.mean + 
                                       yearsbetweensurveys + 
                                       offset(log(sp.SUM.ha.weight2)) +
                                       logtreeNbrJ + Plotcat * mean_spei12 + Plotcat *  sqrtBA.O.plot.1 + Plotcat * logBAj.plot.1 + Plotcat * sqrtBAIj.plot.1.mean.J.I + 
                                       sp.mortality.ha.weight.rate +
                                       mean_spei12 * logBAj.plot.1 + logdbhJ.IMall.plot.mean * logBAj.plot.1 + sqrtBAIj.plot.1.mean.J.I * logBAj.plot.1 + mean_spei12 * sqrtBA.O.plot.1 + logdbhJ.IMall.plot.mean *sqrtBA.O.plot.1 + sqrtBAIj.plot.1.mean.J.I * sqrtBA.O.plot.1,
                                     zi=~logdbhJ.IMall.plot.mean + 
                                       yearsbetweensurveys + 
                                       offset(log(sp.SUM.ha.weight2)) + 
                                       logtreeNbrJ + Plotcat * mean_spei12 + Plotcat *  sqrtBA.O.plot.1 + Plotcat * logBAj.plot.1 + Plotcat * sqrtBAIj.plot.1.mean.J.I + 
                                       sp.mortality.ha.weight.rate +
                                       mean_spei12 * logBAj.plot.1 + logdbhJ.IMall.plot.mean * logBAj.plot.1 + sqrtBAIj.plot.1.mean.J.I * logBAj.plot.1 + mean_spei12 * sqrtBA.O.plot.1 + logdbhJ.IMall.plot.mean *sqrtBA.O.plot.1 + sqrtBAIj.plot.1.mean.J.I * sqrtBA.O.plot.1,family=nbinom1, data=dfplot),silent=T)
        
        try(Saving(rW.yF.sumO.mR,CODE),silent=T)
        try(print(summary(rW.yF.sumO.mR)),silent=T)
        rm(rW.yF.sumO.mR)
        
        # years fixed both part, sum offset both, mortality weight rate year
        try(rW.yF.sumO.mRY <- glmmTMB(sp.recruitment.ha.weight.2 ~ logdbhJ.IMall.plot.mean + 
                                        yearsbetweensurveys + 
                                        offset(log(sp.SUM.ha.weight2)) + 
                                        logtreeNbrJ + Plotcat * mean_spei12 + Plotcat *  sqrtBA.O.plot.1 + Plotcat * logBAj.plot.1 + Plotcat * sqrtBAIj.plot.1.mean.J.I + 
                                        sp.mortality.ha.weight.rate.yr +
                                        mean_spei12 * logBAj.plot.1 + logdbhJ.IMall.plot.mean * logBAj.plot.1 + sqrtBAIj.plot.1.mean.J.I * logBAj.plot.1 + mean_spei12 * sqrtBA.O.plot.1 + logdbhJ.IMall.plot.mean *sqrtBA.O.plot.1 + sqrtBAIj.plot.1.mean.J.I * sqrtBA.O.plot.1,
                                      zi=~logdbhJ.IMall.plot.mean + 
                                        yearsbetweensurveys + 
                                        offset(log(sp.SUM.ha.weight2)) + 
                                        logtreeNbrJ + Plotcat * mean_spei12 + Plotcat *  sqrtBA.O.plot.1 + Plotcat * logBAj.plot.1 + Plotcat * sqrtBAIj.plot.1.mean.J.I + 
                                        sp.mortality.ha.weight.rate.yr +
                                        mean_spei12 * logBAj.plot.1 + logdbhJ.IMall.plot.mean * logBAj.plot.1 + sqrtBAIj.plot.1.mean.J.I * logBAj.plot.1 + mean_spei12 * sqrtBA.O.plot.1 + logdbhJ.IMall.plot.mean *sqrtBA.O.plot.1 + sqrtBAIj.plot.1.mean.J.I * sqrtBA.O.plot.1,family=nbinom1, data=dfplot),silent=T)
        
        try(Saving(rW.yF.sumO.mRY,CODE),silent=T)
        try(print(summary(rW.yF.sumO.mRY)),silent=T)
        rm(rW.yF.sumO.mRY)
        
        
        # years offset both part, sum offset both, mortality weight
        try(rW.yO.sumO.mC <- glmmTMB(sp.recruitment.ha.weight.2 ~ logdbhJ.IMall.plot.mean + 
                                       offset(log(yearsbetweensurveys)) + 
                                       offset(log(sp.SUM.ha.weight2)) + 
                                       logtreeNbrJ + Plotcat * mean_spei12 + Plotcat *  sqrtBA.O.plot.1 + Plotcat * logBAj.plot.1 + Plotcat * sqrtBAIj.plot.1.mean.J.I + 
                                       sp.mortality.ha.weight +
                                       mean_spei12 * logBAj.plot.1 + logdbhJ.IMall.plot.mean * logBAj.plot.1 + sqrtBAIj.plot.1.mean.J.I * logBAj.plot.1 + mean_spei12 * sqrtBA.O.plot.1 + logdbhJ.IMall.plot.mean *sqrtBA.O.plot.1 + sqrtBAIj.plot.1.mean.J.I * sqrtBA.O.plot.1,
                                     zi=~logdbhJ.IMall.plot.mean + 
                                       offset(log(yearsbetweensurveys)) + 
                                       offset(log(sp.SUM.ha.weight2)) + 
                                       logtreeNbrJ + Plotcat * mean_spei12 + Plotcat *  sqrtBA.O.plot.1 + Plotcat * logBAj.plot.1 + Plotcat * sqrtBAIj.plot.1.mean.J.I + 
                                       sp.mortality.ha.weight +
                                       mean_spei12 * logBAj.plot.1 + logdbhJ.IMall.plot.mean * logBAj.plot.1 + sqrtBAIj.plot.1.mean.J.I * logBAj.plot.1 + mean_spei12 * sqrtBA.O.plot.1 + logdbhJ.IMall.plot.mean *sqrtBA.O.plot.1 + sqrtBAIj.plot.1.mean.J.I * sqrtBA.O.plot.1,family=nbinom1, data=dfplot),silent=T)
        
        try(Saving(rW.yO.sumO.mC,CODE),silent=T)
        try(print(summary(rW.yO.sumO.mC)),silent=T)
        rm(rW.yO.sumO.mC)
        
        # years offset both part, sum offset both, mortality weight rate
        try(rW.yO.sumO.mR <- glmmTMB(sp.recruitment.ha.weight.2 ~ logdbhJ.IMall.plot.mean + 
                                       offset(log(yearsbetweensurveys)) + 
                                       offset(log(sp.SUM.ha.weight2)) +
                                       logtreeNbrJ + Plotcat * mean_spei12 + Plotcat *  sqrtBA.O.plot.1 + Plotcat * logBAj.plot.1 + Plotcat * sqrtBAIj.plot.1.mean.J.I + 
                                       sp.mortality.ha.weight.rate +
                                       mean_spei12 * logBAj.plot.1 + logdbhJ.IMall.plot.mean * logBAj.plot.1 + sqrtBAIj.plot.1.mean.J.I * logBAj.plot.1 + mean_spei12 * sqrtBA.O.plot.1 + logdbhJ.IMall.plot.mean *sqrtBA.O.plot.1 + sqrtBAIj.plot.1.mean.J.I * sqrtBA.O.plot.1,
                                     zi=~logdbhJ.IMall.plot.mean + 
                                       offset(log(yearsbetweensurveys)) + 
                                       offset(log(sp.SUM.ha.weight2)) + 
                                       logtreeNbrJ + Plotcat * mean_spei12 + Plotcat *  sqrtBA.O.plot.1 + Plotcat * logBAj.plot.1 + Plotcat * sqrtBAIj.plot.1.mean.J.I + 
                                       sp.mortality.ha.weight.rate +
                                       mean_spei12 * logBAj.plot.1 + logdbhJ.IMall.plot.mean * logBAj.plot.1 + sqrtBAIj.plot.1.mean.J.I * logBAj.plot.1 + mean_spei12 * sqrtBA.O.plot.1 + logdbhJ.IMall.plot.mean *sqrtBA.O.plot.1 + sqrtBAIj.plot.1.mean.J.I * sqrtBA.O.plot.1,family=nbinom1, data=dfplot),silent=T)
        
        try(Saving(rW.yO.sumO.mR,CODE),silent=T)
        try(print(summary(rW.yO.sumO.mR)),silent=T)
        rm(rW.yO.sumO.mR)
        



    
    
    
    
    # 2 combinaisons pour le years betweenb surveyx + 2* sum + 3 options de mortalité
  
  
  
  
  # Nouvelle series avec offset VS effet fixe sur les deux parties et une en modélisant un countage 
  # offset sur le nombre d'arbre total du plot et offset sur le temps d'observation (nbr years) 
  
  
#   
#   try(SCALEDRrate.yr.QualiF.TF.YF.nbinom1 <- glmmTMB(sp.recruitment.plot.rate.yr.2 ~ logdbhJ.IMall.plot.mean + yearsbetweensurveys +  (1 | country) + logtreeNbrJ + Plotcat * mean_spei12 + Plotcat *  
#                                                        sqrtBA.O.plot.1 + Plotcat * logBAj.plot.1 + Plotcat * sqrtBAIj.plot.1.mean.J.I + Plotcat * logsp.mortality.plot.rate.yr + mean_spei12 * logBAj.plot.1 +  
#                                                        logdbhJ.IMall.plot.mean * logBAj.plot.1 + sqrtBAIj.plot.1.mean.J.I * logBAj.plot.1 + mean_spei12 * sqrtBA.O.plot.1 + logdbhJ.IMall.plot.mean *  
#                                                        sqrtBA.O.plot.1 + sqrtBAIj.plot.1.mean.J.I * sqrtBA.O.plot.1, 
#                                                      zi=~.,family=nbinom1, data=dfplot),silent=T)
# try(Saving(SCALEDRrate.yr.QualiF.TF.YF.nbinom1),silent=T) 
#   try(print(summary(SCALEDRrate.yr.QualiF.TF.YF.nbinom1)),silent=T)
#   rm(SCALEDRrate.yr.QualiF.TF.YF.nbinom1)
#   
#   
#   
#   #### Years offset both part !!!!!!! ####### 
#   try(SCALEDRrate.yr.QualiF.TF.YO.nbinom1 <- glmmTMB(sp.recruitment.plot.rate.yr.2 ~ logdbhJ.IMall.plot.mean + offset(log(yearsbetweensurveys)) +  (1 | country) + logtreeNbrJ + Plotcat * mean_spei12 + Plotcat *  
#                                                                     sqrtBA.O.plot.1 + Plotcat * logBAj.plot.1 + Plotcat * sqrtBAIj.plot.1.mean.J.I + Plotcat * logsp.mortality.plot.rate.yr + mean_spei12 * logBAj.plot.1 +  
#                                                                     logdbhJ.IMall.plot.mean * logBAj.plot.1 + sqrtBAIj.plot.1.mean.J.I * logBAj.plot.1 + mean_spei12 * sqrtBA.O.plot.1 + logdbhJ.IMall.plot.mean *  
#                                                                     sqrtBA.O.plot.1 + sqrtBAIj.plot.1.mean.J.I * sqrtBA.O.plot.1, 
#                                                                   zi=~logdbhJ.IMall.plot.mean + offset(log(yearsbetweensurveys)) +  (1 | country) + logtreeNbrJ + Plotcat * mean_spei12 + Plotcat *  
#                                                                     sqrtBA.O.plot.1 + Plotcat * logBAj.plot.1 + Plotcat * sqrtBAIj.plot.1.mean.J.I + Plotcat * logsp.mortality.plot.rate.yr + mean_spei12 * logBAj.plot.1 +  
#                                                                     logdbhJ.IMall.plot.mean * logBAj.plot.1 + sqrtBAIj.plot.1.mean.J.I * logBAj.plot.1 + mean_spei12 * sqrtBA.O.plot.1 + logdbhJ.IMall.plot.mean *  
#                                                                     sqrtBA.O.plot.1 + sqrtBAIj.plot.1.mean.J.I * sqrtBA.O.plot.1
#                                                                   ,family=nbinom1, data=dfplot),silent=T)
#   try(Saving(SCALEDRrate.yr.QualiF.TF.YO.nbinom1),silent=T) 
#   try(print(summary(SCALEDRrate.yr.QualiF.TF.YO.nbinom1)),silent=T)
#   rm(SCALEDRrate.yr.QualiF.TF.YO.nbinom1)
#   
#   #### Years fixed second part only
#   try(SCALEDRrate.yr.QualiF.TF.YF2.nbinom1 <- glmmTMB(sp.recruitment.plot.rate.yr.2 ~ logdbhJ.IMall.plot.mean +  (1 | country) + logtreeNbrJ + Plotcat * mean_spei12 + Plotcat *  
#                                                                     sqrtBA.O.plot.1 + Plotcat * logBAj.plot.1 + Plotcat * sqrtBAIj.plot.1.mean.J.I + Plotcat * logsp.mortality.plot.rate.yr + mean_spei12 * logBAj.plot.1 +  
#                                                                     logdbhJ.IMall.plot.mean * logBAj.plot.1 + sqrtBAIj.plot.1.mean.J.I * logBAj.plot.1 + mean_spei12 * sqrtBA.O.plot.1 + logdbhJ.IMall.plot.mean *  
#                                                                     sqrtBA.O.plot.1 + sqrtBAIj.plot.1.mean.J.I * sqrtBA.O.plot.1, 
#                                                                   zi=~.+yearsbetweensurveys,family=nbinom1, data=dfplot),silent=T)
#   try(Saving(SCALEDRrate.yr.QualiF.TF.YF2.nbinom1),silent=T) 
#   try(print(summary(SCALEDRrate.yr.QualiF.TF.YF2.nbinom1)),silent=T)
#   rm(SCALEDRrate.yr.QualiF.TF.YF2.nbinom1)
#   
#   #### Years offset second part only
#   try(SCALEDRrate.yr.QualiF.TF.YO2.nbinom1 <- glmmTMB(sp.recruitment.plot.rate.yr.2 ~ logdbhJ.IMall.plot.mean  +  (1 | country) + logtreeNbrJ + Plotcat * mean_spei12 + Plotcat *  
#                                                                     sqrtBA.O.plot.1 + Plotcat * logBAj.plot.1 + Plotcat * sqrtBAIj.plot.1.mean.J.I + Plotcat * logsp.mortality.plot.rate.yr + mean_spei12 * logBAj.plot.1 +  
#                                                                     logdbhJ.IMall.plot.mean * logBAj.plot.1 + sqrtBAIj.plot.1.mean.J.I * logBAj.plot.1 + mean_spei12 * sqrtBA.O.plot.1 + logdbhJ.IMall.plot.mean *  
#                                                                     sqrtBA.O.plot.1 + sqrtBAIj.plot.1.mean.J.I * sqrtBA.O.plot.1, 
#                                                                   zi=~logdbhJ.IMall.plot.mean + offset(log(yearsbetweensurveys)) +  (1 | country) + logtreeNbrJ + Plotcat * mean_spei12 + Plotcat *  
#                                                                     sqrtBA.O.plot.1 + Plotcat * logBAj.plot.1 + Plotcat * sqrtBAIj.plot.1.mean.J.I + Plotcat * logsp.mortality.plot.rate.yr + mean_spei12 * logBAj.plot.1 +  
#                                                                     logdbhJ.IMall.plot.mean * logBAj.plot.1 + sqrtBAIj.plot.1.mean.J.I * logBAj.plot.1 + mean_spei12 * sqrtBA.O.plot.1 + logdbhJ.IMall.plot.mean *  
#                                                                     sqrtBA.O.plot.1 + sqrtBAIj.plot.1.mean.J.I * sqrtBA.O.plot.1
#                                                                   ,family=nbinom1, data=dfplot),silent=T)
#   try(Saving(SCALEDRrate.yr.QualiF.TF.YO2.nbinom1),silent=T) 
#   try(print(summary(SCALEDRrate.yr.QualiF.TF.YO2.nbinom1)),silent=T)
#   rm(SCALEDRrate.yr.QualiF.TF.YO2.nbinom1)
#   

  
  # 
  # 
  # ## Model de base sans effet country
  # #### Years Fixed in both part!!!!!!! #######
  # 
  # try(SCALEDRrate.yr.QualiF.TF.YF.nbinom1 <- glmmTMB(sp.recruitment.plot.rate.yr.2 ~ logdbhJ.IMall.plot.mean + yearsbetweensurveys + logtreeNbrJ + Plotcat * mean_spei12 + Plotcat *
  #                                                                   sqrtBA.O.plot.1 + Plotcat * logBAj.plot.1 + Plotcat * sqrtBAIj.plot.1.mean.J.I + Plotcat * logsp.mortality.plot.rate.yr + mean_spei12 * logBAj.plot.1 +
  #                                                                   logdbhJ.IMall.plot.mean * logBAj.plot.1 + sqrtBAIj.plot.1.mean.J.I * logBAj.plot.1 + mean_spei12 * sqrtBA.O.plot.1 + logdbhJ.IMall.plot.mean *
  #                                                                   sqrtBA.O.plot.1 + sqrtBAIj.plot.1.mean.J.I * sqrtBA.O.plot.1,
  #                                                                 zi=~.,family=nbinom1, data=dfplot),silent=T)
  # try(Saving(SCALEDRrate.yr.QualiF.TF.YF.nbinom1),silent=T)
  # try(print(summary(SCALEDRrate.yr.QualiF.TF.YF.nbinom1)),silent=T)
  # rm(SCALEDRrate.yr.QualiF.TF.YF.nbinom1)
  # 
  # 
  # #### Years offset both part !!!!!!! #######
  # try(SCALEDRrate.yr.QualiF.TF.YO.nbinom1 <- glmmTMB(sp.recruitment.plot.rate.yr.2 ~ logdbhJ.IMall.plot.mean + offset(log(yearsbetweensurveys)) + logtreeNbrJ + Plotcat * mean_spei12 + Plotcat *
  #                                                                   sqrtBA.O.plot.1 + Plotcat * logBAj.plot.1 + Plotcat * sqrtBAIj.plot.1.mean.J.I + Plotcat * logsp.mortality.plot.rate.yr + mean_spei12 * logBAj.plot.1 +
  #                                                                   logdbhJ.IMall.plot.mean * logBAj.plot.1 + sqrtBAIj.plot.1.mean.J.I * logBAj.plot.1 + mean_spei12 * sqrtBA.O.plot.1 + logdbhJ.IMall.plot.mean *
  #                                                                   sqrtBA.O.plot.1 + sqrtBAIj.plot.1.mean.J.I * sqrtBA.O.plot.1,
  #                                                                 zi=~logdbhJ.IMall.plot.mean + offset(log(yearsbetweensurveys)) + logtreeNbrJ + Plotcat * mean_spei12 + Plotcat *
  #                                                                   sqrtBA.O.plot.1 + Plotcat * logBAj.plot.1 + Plotcat * sqrtBAIj.plot.1.mean.J.I + Plotcat * logsp.mortality.plot.rate.yr + mean_spei12 * logBAj.plot.1 +
  #                                                                   logdbhJ.IMall.plot.mean * logBAj.plot.1 + sqrtBAIj.plot.1.mean.J.I * logBAj.plot.1 + mean_spei12 * sqrtBA.O.plot.1 + logdbhJ.IMall.plot.mean *
  #                                                                   sqrtBA.O.plot.1 + sqrtBAIj.plot.1.mean.J.I * sqrtBA.O.plot.1
  #                                                                 ,family=nbinom1, data=dfplot),silent=T)
  # try(Saving(SCALEDRrate.yr.QualiF.TF.YO.nbinom1),silent=T)
  # try(print(summary(SCALEDRrate.yr.QualiF.TF.YO.nbinom1)),silent=T)
  # rm(SCALEDRrate.yr.QualiF.TF.YO.nbinom1)
  # 
  # #### Years fixed second part only
  # try(SCALEDRrate.yr.QualiF.TF.YF2.nbinom1 <- glmmTMB(sp.recruitment.plot.rate.yr.2 ~ logdbhJ.IMall.plot.mean + logtreeNbrJ + Plotcat * mean_spei12 + Plotcat *
  #                                                                    sqrtBA.O.plot.1 + Plotcat * logBAj.plot.1 + Plotcat * sqrtBAIj.plot.1.mean.J.I + Plotcat * logsp.mortality.plot.rate.yr + mean_spei12 * logBAj.plot.1 +
  #                                                                    logdbhJ.IMall.plot.mean * logBAj.plot.1 + sqrtBAIj.plot.1.mean.J.I * logBAj.plot.1 + mean_spei12 * sqrtBA.O.plot.1 + logdbhJ.IMall.plot.mean *
  #                                                                    sqrtBA.O.plot.1 + sqrtBAIj.plot.1.mean.J.I * sqrtBA.O.plot.1,
  #                                                                  zi=~.+yearsbetweensurveys,family=nbinom1, data=dfplot),silent=T)
  # try(Saving(SCALEDRrate.yr.QualiF.TF.YF2.nbinom1),silent=T)
  # try(print(summary(SCALEDRrate.yr.QualiF.TF.YF2.nbinom1)),silent=T)
  # rm(SCALEDRrate.yr.QualiF.TF.YF2.nbinom1)
  # 
  # #### Years offset second part only
  # try(SCALEDRrate.yr.QualiF.TF.YO2.nbinom1 <- glmmTMB(sp.recruitment.plot.rate.yr.2 ~ logdbhJ.IMall.plot.mean + logtreeNbrJ + Plotcat * mean_spei12 + Plotcat *
  #                                                                    sqrtBA.O.plot.1 + Plotcat * logBAj.plot.1 + Plotcat * sqrtBAIj.plot.1.mean.J.I + Plotcat * logsp.mortality.plot.rate.yr + mean_spei12 * logBAj.plot.1 +
  #                                                                    logdbhJ.IMall.plot.mean * logBAj.plot.1 + sqrtBAIj.plot.1.mean.J.I * logBAj.plot.1 + mean_spei12 * sqrtBA.O.plot.1 + logdbhJ.IMall.plot.mean *
  #                                                                    sqrtBA.O.plot.1 + sqrtBAIj.plot.1.mean.J.I * sqrtBA.O.plot.1,
  #                                                                  zi=~logdbhJ.IMall.plot.mean + offset(log(yearsbetweensurveys)) + logtreeNbrJ + Plotcat * mean_spei12 + Plotcat *
  #                                                                    sqrtBA.O.plot.1 + Plotcat * logBAj.plot.1 + Plotcat * sqrtBAIj.plot.1.mean.J.I + Plotcat * logsp.mortality.plot.rate.yr + mean_spei12 * logBAj.plot.1 +
  #                                                                    logdbhJ.IMall.plot.mean * logBAj.plot.1 + sqrtBAIj.plot.1.mean.J.I * logBAj.plot.1 + mean_spei12 * sqrtBA.O.plot.1 + logdbhJ.IMall.plot.mean *
  #                                                                    sqrtBA.O.plot.1 + sqrtBAIj.plot.1.mean.J.I * sqrtBA.O.plot.1
  #                                                                  ,family=nbinom1, data=dfplot),silent=T)
  # try(Saving(SCALEDRrate.yr.QualiF.TF.YO2.nbinom1),silent=T)
  # try(print(summary(SCALEDRrate.yr.QualiF.TF.YO2.nbinom1)),silent=T)
  # rm(SCALEDRrate.yr.QualiF.TF.YO2.nbinom1)


  # 

  ## Model sans interaction mortality
  #### Years Fixed in both part!!!!!!! #######

  # try(SCALEDRrate.yr.QualiF.TF.YF.nbinom1 <- glmmTMB(sp.recruitment.plot.rate.yr.2 ~ logdbhJ.IMall.plot.mean + yearsbetweensurveys +  (1 | country) + logtreeNbrJ + Plotcat * mean_spei12 + Plotcat *
  #                                                                   sqrtBA.O.plot.1 + Plotcat * logBAj.plot.1 + Plotcat * sqrtBAIj.plot.1.mean.J.I + logsp.mortality.plot.rate.yr + mean_spei12 * logBAj.plot.1 +
  #                                                                   logdbhJ.IMall.plot.mean * logBAj.plot.1 + sqrtBAIj.plot.1.mean.J.I * logBAj.plot.1 + mean_spei12 * sqrtBA.O.plot.1 + logdbhJ.IMall.plot.mean *
  #                                                                   sqrtBA.O.plot.1 + sqrtBAIj.plot.1.mean.J.I * sqrtBA.O.plot.1,
  #                                                                 zi=~.,family=nbinom1, data=dfplot),silent=T)
  # try(Saving(SCALEDRrate.yr.QualiF.TF.YF.nbinom1),silent=T)
  # try(print(summary(SCALEDRrate.yr.QualiF.TF.YF.nbinom1)),silent=T)
  # rm(SCALEDRrate.yr.QualiF.TF.YF.nbinom1)
  # 
  # #### Years offset both part !!!!!!! #######
  # try(SCALEDRrate.yr.QualiF.TF.YO.nbinom1 <- glmmTMB(sp.recruitment.plot.rate.yr.2 ~ logdbhJ.IMall.plot.mean + offset(log(yearsbetweensurveys)) +  (1 | country) + logtreeNbrJ + Plotcat * mean_spei12 + Plotcat *
  #                                                                   sqrtBA.O.plot.1 + Plotcat * logBAj.plot.1 + Plotcat * sqrtBAIj.plot.1.mean.J.I + logsp.mortality.plot.rate.yr + mean_spei12 * logBAj.plot.1 +
  #                                                                   logdbhJ.IMall.plot.mean * logBAj.plot.1 + sqrtBAIj.plot.1.mean.J.I * logBAj.plot.1 + mean_spei12 * sqrtBA.O.plot.1 + logdbhJ.IMall.plot.mean *
  #                                                                   sqrtBA.O.plot.1 + sqrtBAIj.plot.1.mean.J.I * sqrtBA.O.plot.1,
  #                                                                 zi=~logdbhJ.IMall.plot.mean + offset(log(yearsbetweensurveys)) +  (1 | country) + logtreeNbrJ + Plotcat * mean_spei12 + Plotcat *
  #                                                                   sqrtBA.O.plot.1 + Plotcat * logBAj.plot.1 + Plotcat * sqrtBAIj.plot.1.mean.J.I + logsp.mortality.plot.rate.yr + mean_spei12 * logBAj.plot.1 +
  #                                                                   logdbhJ.IMall.plot.mean * logBAj.plot.1 + sqrtBAIj.plot.1.mean.J.I * logBAj.plot.1 + mean_spei12 * sqrtBA.O.plot.1 + logdbhJ.IMall.plot.mean *
  #                                                                   sqrtBA.O.plot.1 + sqrtBAIj.plot.1.mean.J.I * sqrtBA.O.plot.1
  #                                                                 ,family=nbinom1, data=dfplot),silent=T)
  # try(Saving(SCALEDRrate.yr.QualiF.TF.YO.nbinom1),silent=T)
  # try(print(summary(SCALEDRrate.yr.QualiF.TF.YO.nbinom1)),silent=T)
  # rm(SCALEDRrate.yr.QualiF.TF.YO.nbinom1)
  # 
  # #### Years fixed second part only
  # try(SCALEDRrate.yr.QualiF.TF.YF2.nbinom1 <- glmmTMB(sp.recruitment.plot.rate.yr.2 ~ logdbhJ.IMall.plot.mean +  (1 | country) + logtreeNbrJ + Plotcat * mean_spei12 + Plotcat *
  #                                                                    sqrtBA.O.plot.1 + Plotcat * logBAj.plot.1 + Plotcat * sqrtBAIj.plot.1.mean.J.I + logsp.mortality.plot.rate.yr + mean_spei12 * logBAj.plot.1 +
  #                                                                    logdbhJ.IMall.plot.mean * logBAj.plot.1 + sqrtBAIj.plot.1.mean.J.I * logBAj.plot.1 + mean_spei12 * sqrtBA.O.plot.1 + logdbhJ.IMall.plot.mean *
  #                                                                    sqrtBA.O.plot.1 + sqrtBAIj.plot.1.mean.J.I * sqrtBA.O.plot.1,
  #                                                                  zi=~.+yearsbetweensurveys,family=nbinom1, data=dfplot),silent=T)
  # try(Saving(SCALEDRrate.yr.QualiF.TF.YF2.nbinom1),silent=T)
  # try(print(summary(SCALEDRrate.yr.QualiF.TF.YF2.nbinom1)),silent=T)
  # rm(SCALEDRrate.yr.QualiF.TF.YF2.nbinom1)
  # 
  # 
  # #### Years offset second part only
  # try(SCALEDRrate.yr.QualiF.TF.YO2.nbinom1 <- glmmTMB(sp.recruitment.plot.rate.yr.2 ~ logdbhJ.IMall.plot.mean  +  (1 | country) + logtreeNbrJ + Plotcat * mean_spei12 + Plotcat *
  #                                                                    sqrtBA.O.plot.1 + Plotcat * logBAj.plot.1 + Plotcat * sqrtBAIj.plot.1.mean.J.I + logsp.mortality.plot.rate.yr + mean_spei12 * logBAj.plot.1 +
  #                                                                    logdbhJ.IMall.plot.mean * logBAj.plot.1 + sqrtBAIj.plot.1.mean.J.I * logBAj.plot.1 + mean_spei12 * sqrtBA.O.plot.1 + logdbhJ.IMall.plot.mean *
  #                                                                    sqrtBA.O.plot.1 + sqrtBAIj.plot.1.mean.J.I * sqrtBA.O.plot.1,
  #                                                                  zi=~logdbhJ.IMall.plot.mean + offset(log(yearsbetweensurveys)) +  (1 | country) + logtreeNbrJ + Plotcat * mean_spei12 + Plotcat *
  #                                                                    sqrtBA.O.plot.1 + Plotcat * logBAj.plot.1 + Plotcat * sqrtBAIj.plot.1.mean.J.I + logsp.mortality.plot.rate.yr + mean_spei12 * logBAj.plot.1 +
  #                                                                    logdbhJ.IMall.plot.mean * logBAj.plot.1 + sqrtBAIj.plot.1.mean.J.I * logBAj.plot.1 + mean_spei12 * sqrtBA.O.plot.1 + logdbhJ.IMall.plot.mean *
  #                                                                    sqrtBA.O.plot.1 + sqrtBAIj.plot.1.mean.J.I * sqrtBA.O.plot.1
  #                                                                  ,family=nbinom1, data=dfplot),silent=T)
  # try(Saving(SCALEDRrate.yr.QualiF.TF.YO2.nbinom1),silent=T)
  # try(print(summary(SCALEDRrate.yr.QualiF.TF.YO2.nbinom1)),silent=T)
  # rm(SCALEDRrate.yr.QualiF.TF.YO2.nbinom1)

  #
  # Model without both
  #### Years Fixed in both part!!!!!!! #######
# 
#   try(SCALEDRrate.yr.QualiF.TF.YF.nbinom1 <- glmmTMB(sp.recruitment.plot.rate.yr.2 ~ logdbhJ.IMall.plot.mean + yearsbetweensurveys + logtreeNbrJ + Plotcat * mean_spei12 + Plotcat *
#                                                                     sqrtBA.O.plot.1 + Plotcat * logBAj.plot.1 + Plotcat * sqrtBAIj.plot.1.mean.J.I + logsp.mortality.plot.rate.yr + mean_spei12 * logBAj.plot.1 +
#                                                                     logdbhJ.IMall.plot.mean * logBAj.plot.1 + sqrtBAIj.plot.1.mean.J.I * logBAj.plot.1 + mean_spei12 * sqrtBA.O.plot.1 + logdbhJ.IMall.plot.mean *
#                                                                     sqrtBA.O.plot.1 + sqrtBAIj.plot.1.mean.J.I * sqrtBA.O.plot.1,
#                                                                   zi=~.,family=nbinom1, data=dfplot),silent=T)
#   try(Saving(SCALEDRrate.yr.QualiF.TF.YF.nbinom1),silent=T)
#   try(print(summary(SCALEDRrate.yr.QualiF.TF.YF.nbinom1)),silent=T)
#   rm(SCALEDRrate.yr.QualiF.TF.YF.nbinom1)
# 
#   #### Years offset both part !!!!!!! #######
#   try(SCALEDRrate.yr.QualiF.TF.YO.nbinom1 <- glmmTMB(sp.recruitment.plot.rate.yr.2 ~ logdbhJ.IMall.plot.mean + offset(log(yearsbetweensurveys)) + logtreeNbrJ + Plotcat * mean_spei12 + Plotcat *
#                                                                     sqrtBA.O.plot.1 + Plotcat * logBAj.plot.1 + Plotcat * sqrtBAIj.plot.1.mean.J.I + logsp.mortality.plot.rate.yr + mean_spei12 * logBAj.plot.1 +
#                                                                     logdbhJ.IMall.plot.mean * logBAj.plot.1 + sqrtBAIj.plot.1.mean.J.I * logBAj.plot.1 + mean_spei12 * sqrtBA.O.plot.1 + logdbhJ.IMall.plot.mean *
#                                                                     sqrtBA.O.plot.1 + sqrtBAIj.plot.1.mean.J.I * sqrtBA.O.plot.1,
#                                                                   zi=~logdbhJ.IMall.plot.mean + offset(log(yearsbetweensurveys)) + logtreeNbrJ + Plotcat * mean_spei12 + Plotcat *
#                                                                     sqrtBA.O.plot.1 + Plotcat * logBAj.plot.1 + Plotcat * sqrtBAIj.plot.1.mean.J.I + logsp.mortality.plot.rate.yr + mean_spei12 * logBAj.plot.1 +
#                                                                     logdbhJ.IMall.plot.mean * logBAj.plot.1 + sqrtBAIj.plot.1.mean.J.I * logBAj.plot.1 + mean_spei12 * sqrtBA.O.plot.1 + logdbhJ.IMall.plot.mean *
#                                                                     sqrtBA.O.plot.1 + sqrtBAIj.plot.1.mean.J.I * sqrtBA.O.plot.1
#                                                                   ,family=nbinom1, data=dfplot),silent=T)
#   try(Saving(SCALEDRrate.yr.QualiF.TF.YO.nbinom1),silent=T)
#   try(print(summary(SCALEDRrate.yr.QualiF.TF.YO.nbinom1)),silent=T)
#   rm(SCALEDRrate.yr.QualiF.TF.YO.nbinom1)
# 
#   #### Years fixed second part only
#   try(SCALEDRrate.yr.QualiF.TF.YF2.nbinom1 <- glmmTMB(sp.recruitment.plot.rate.yr.2 ~ logdbhJ.IMall.plot.mean + logtreeNbrJ + Plotcat * mean_spei12 + Plotcat *
#                                                                      sqrtBA.O.plot.1 + Plotcat * logBAj.plot.1 + Plotcat * sqrtBAIj.plot.1.mean.J.I + logsp.mortality.plot.rate.yr + mean_spei12 * logBAj.plot.1 +
#                                                                      logdbhJ.IMall.plot.mean * logBAj.plot.1 + sqrtBAIj.plot.1.mean.J.I * logBAj.plot.1 + mean_spei12 * sqrtBA.O.plot.1 + logdbhJ.IMall.plot.mean *
#                                                                      sqrtBA.O.plot.1 + sqrtBAIj.plot.1.mean.J.I * sqrtBA.O.plot.1,
#                                                                    zi=~.+yearsbetweensurveys,family=nbinom1, data=dfplot),silent=T)
#   try(Saving(SCALEDRrate.yr.QualiF.TF.YF2.nbinom1),silent=T)
#   try(print(summary(SCALEDRrate.yr.QualiF.TF.YF2.nbinom1)),silent=T)
#   rm(SCALEDRrate.yr.QualiF.TF.YF2.nbinom1)
# 
#   #### Years offset second part only
#   try(SCALEDRrate.yr.QualiF.TF.YO2.nbinom1 <- glmmTMB(sp.recruitment.plot.rate.yr.2 ~ logdbhJ.IMall.plot.mean + logtreeNbrJ + Plotcat * mean_spei12 + Plotcat *
#                                                                      sqrtBA.O.plot.1 + Plotcat * logBAj.plot.1 + Plotcat * sqrtBAIj.plot.1.mean.J.I + logsp.mortality.plot.rate.yr + mean_spei12 * logBAj.plot.1 +
#                                                                      logdbhJ.IMall.plot.mean * logBAj.plot.1 + sqrtBAIj.plot.1.mean.J.I * logBAj.plot.1 + mean_spei12 * sqrtBA.O.plot.1 + logdbhJ.IMall.plot.mean *
#                                                                      sqrtBA.O.plot.1 + sqrtBAIj.plot.1.mean.J.I * sqrtBA.O.plot.1,
#                                                                    zi=~logdbhJ.IMall.plot.mean + offset(log(yearsbetweensurveys)) + logtreeNbrJ + Plotcat * mean_spei12 + Plotcat *
#                                                                      sqrtBA.O.plot.1 + Plotcat * logBAj.plot.1 + Plotcat * sqrtBAIj.plot.1.mean.J.I + logsp.mortality.plot.rate.yr + mean_spei12 * logBAj.plot.1 +
#                                                                      logdbhJ.IMall.plot.mean * logBAj.plot.1 + sqrtBAIj.plot.1.mean.J.I * logBAj.plot.1 + mean_spei12 * sqrtBA.O.plot.1 + logdbhJ.IMall.plot.mean *
#                                                                      sqrtBA.O.plot.1 + sqrtBAIj.plot.1.mean.J.I * sqrtBA.O.plot.1
#                                                                    ,family=nbinom1, data=dfplot),silent=T)
#   try(Saving(SCALEDRrate.yr.QualiF.TF.YO2.nbinom1),silent=T)
#   try(print(summary(SCALEDRrate.yr.QualiF.TF.YO2.nbinom1)),silent=T)
#   rm(SCALEDRrate.yr.QualiF.TF.YO2.nbinom1)
}
  
  
  
  
  