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

Dir <- c("/home/achangenet/SynologyDrive/FUNDIV - NFI - Europe/")
setwd(Dir)
source(paste0(Dir,"Myscripts/Fundiv.project/RecrutProject.function1.Saving.R"))
Allcode <- c("PINPINA","ABIALB","ACEPSE","ALNGLU","BETPEN","CASSAT","FAGSYL","FRAEXC","PICABI","PINHAL","PINNIG","PINPIN","PINSYL","POPNIG","POPTRE","QUEILE","QUEPET",
             "QUEPYR","QUEROB","QUESUB")
Allseuil <- c(0.8,0.7,0.55,0.6,0.8,0.7,0.7,0.7,0.8,0.7,0.7,0.7,0.8,0.7,0.7,0.8,0.7,0.7,0.8,0.7)
AllseuilC <- c(0.6,0.6,0.5,0.5,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6)

Errors.files <- file(paste0(Dir,"our-data/species/Errors.allmodels.Recrut.Rout"), open="wt")
sink(Errors.files, type="message")
# do the same with 0.7
for (i in 1:length(Allcode)){
  CODE <- Allcode[i]
  seuilC <- AllseuilC[i]
  seuil <- Allseuil[i]
  Dir = (paste0("/home/achangenet/SynologyDrive/FUNDIV - NFI - Europe/our-data/species/",CODE,"/CLIMAP/Recrut.Mortality.2020/")) # Directory 
  setwd(Dir)
  dfplot <- try(readRDS(paste0("dfplotFinal",CODE,seuil,"R.M.rds")),silent = T) #Base de données plot
  
  # Add recruitment as a count per year and as a total count
  dfplot$sp.recruitment.plot.rate.yr.2 <- round((dfplot$sp.recruitment.plot.rate * 1000 / dfplot$yearsbetweensurveys))
  dfplot$sp.recruitment.plot.rate.2 <- round(dfplot$sp.recruitment.plot.rate * 1000)
  dfplot$sp.mortality.plot.rate.2 <- round(dfplot$sp.mortality.plot.rate * 1000)
  

### Recrut rate Marginality Quali Fixed Recrut Treenbr fixe yearsbetweenstudy fixed
  try(Rrate.QualiF.TF.YF.nbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys + treeNbrJ +
                Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.mortality.plot.rate,
              zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QualiF.TF.YF.nbinom1),silent=T) 
  ### Recrut rate Marginality Quanti Fixed Recrut Treenbr fixe yearsbetweenstudy fixed
  try(Rrate.QuantiF.TF.YF.nbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys + treeNbrJ +
                                                        M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.mortality.plot.rate,
                                                      zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QuantiF.TF.YF.nbinom1),silent=T) 

  ### Recrut rate Marginality Quali Aleatoire Recrut Treenbr fixe yearsbetweenstudy fixed 
  try(Rrate.QualiA.TF.YF.nbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys + treeNbrJ +
                                                        (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate,
                                                      zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QualiA.TF.YF.nbinom1),silent=T) 
  
  ### Recrut rate Marginality Quanti Aleatoire Recrut Treenbr fixe yearsbetweenstudy fixed
  try(Rrate.QuantiA.TF.YF.nbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys + treeNbrJ +
                                                         (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate,
                                                       zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QuantiA.TF.YF.nbinom1),silent=T) 
  
  
  
  ## Offset tree 
  ### Recrut rate Marginality Quali Fixed Recrut Treenbr offset yearsbetweenstudy fixed
  try(Rrate.QualiF.TO.YF.nbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys + offset(log(treeNbrJ)) +
                                                        Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.mortality.plot.rate,
                                                      zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QualiF.TO.YF.nbinom1),silent=T) 
  
  ### Recrut rate Marginality Quanti Fixed Recrut Treenbr offset yearsbetweenstudy fixed
  try(Rrate.QuantiF.TO.YF.nbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys + offset(log(treeNbrJ)) +
                                                         M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.mortality.plot.rate,
                                                       zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QuantiF.TO.YF.nbinom1),silent=T) 
  
  ### Recrut rate Marginality Quali Aleatoire Recrut Treenbr offset yearsbetweenstudy fixed 
  try(Rrate.QualiA.TO.YF.nbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys + offset(log(treeNbrJ)) +
                                                        (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate,
                                                      zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QualiA.TO.YF.nbinom1),silent=T) 
  
  ### Recrut rate Marginality Quanti Aleatoire Recrut Treenbr offset yearsbetweenstudy fixed
  try(Rrate.QuantiA.TO.YF.nbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys + offset(log(treeNbrJ)) +
                                                         (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate,
                                                       zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QuantiA.TO.YF.nbinom1),silent=T) 
  
  
  ## No tree parameter 
  ### Recrut rate Marginality Quali Fixed Recrut no treenbr yearsbetweenstudy fixed
  try(Rrate.QualiF.noT.YF.nbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys  +
                                                        Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.mortality.plot.rate,
                                                      zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QualiF.noT.YF.nbinom1),silent=T) 
  
  ### Recrut rate Marginality Quanti Fixed Recrut no treenbr yearsbetweenstudy fixed
  try(Rrate.QuantiF.noT.YF.nbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys  +
                                                         M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.mortality.plot.rate,
                                                       zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QuantiF.noT.YF.nbinom1),silent=T) 
  
  ### Recrut rate Marginality Quali AleanoTire Recrut no treenbr yearsbetweenstudy fixed 
  try(Rrate.QualiA.noT.YF.nbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys  +
                                                        (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate,
                                                      zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QualiA.noT.YF.nbinom1),silent=T) 
  
  ### Recrut rate Marginality Quanti AleanoTire Recrut no treenbr yearsbetweenstudy fixed
  try(Rrate.QuantiA.noT.YF.nbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys  +
                                                         (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate,
                                                       zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QuantiA.noT.YF.nbinom1),silent=T) 


  
  
  #### Years between study as an offset 
  ### Recrut rate Marginality Quali Fixed Recrut Treenbr fixe yearsbetween study offset
  try(Rrate.QualiF.TF.YO.nbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys)) + treeNbrJ +
                                                        Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.mortality.plot.rate,
                                                      zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QualiF.TF.YO.nbinom1),silent=T) 
  
  ### Recrut rate Marginality Quanti Fixed Recrut Treenbr fixe yearsbetween study offset
  try(Rrate.QuantiF.TF.YO.nbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys)) + treeNbrJ +
                                                         M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.mortality.plot.rate,
                                                       zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QuantiF.TF.YO.nbinom1),silent=T) 
  
  ### Recrut rate Marginality Quali Aleatoire Recrut Treenbr fixe yearsbetween study offset 
  try(Rrate.QualiA.TF.YO.nbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys)) + treeNbrJ +
                                                        (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate,
                                                      zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QualiA.TF.YO.nbinom1),silent=T) 
  
  ### Recrut rate Marginality Quanti Aleatoire Recrut Treenbr fixe yearsbetween study offset
  try(Rrate.QuantiA.TF.YO.nbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys)) + treeNbrJ +
                                                         (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate,
                                                       zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QuantiA.TF.YO.nbinom1),silent=T) 
  
  
  
  ## Offset tree 
  ### Recrut rate Marginality Quali Fixed Recrut Treenbr offset yearsbetween study offset
  try(Rrate.QualiF.TO.YO.nbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys)) + offset(log(treeNbrJ)) +
                                                        Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.mortality.plot.rate,
                                                      zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QualiF.TO.YO.nbinom1),silent=T) 
  
  ### Recrut rate Marginality Quanti Fixed Recrut Treenbr offset yearsbetween study offset
  try(Rrate.QuantiF.TO.YO.nbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys)) + offset(log(treeNbrJ)) +
                                                         M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.mortality.plot.rate,
                                                       zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QuantiF.TO.YO.nbinom1),silent=T) 
  
  ### Recrut rate Marginality Quali Aleatoire Recrut Treenbr offset yearsbetween study offset 
  try(Rrate.QualiA.TO.YO.nbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys)) + offset(log(treeNbrJ)) +
                                                        (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate,
                                                      zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QualiA.TO.YO.nbinom1),silent=T) 
  
  ### Recrut rate Marginality Quanti Aleatoire Recrut Treenbr offset yearsbetween study offset
  try(Rrate.QuantiA.TO.YO.nbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys)) + offset(log(treeNbrJ)) +
                                                         (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate,
                                                       zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QuantiA.TO.YO.nbinom1),silent=T) 
  
  
  ## No tree parameter 
  ### Recrut rate Marginality Quali Fixed Recrut no treenbr yearsbetween study offset
  try(Rrate.QualiF.noT.YO.nbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys))  +
                                                         Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.mortality.plot.rate,
                                                       zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QualiF.noT.YO.nbinom1),silent=T) 
  
  ### Recrut rate Marginality Quanti Fixed Recrut no treenbr yearsbetween study offset
  try(Rrate.QuantiF.noT.YO.nbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys))  +
                                                          M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.mortality.plot.rate,
                                                        zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QuantiF.noT.YO.nbinom1),silent=T) 
  
  ### Recrut rate Marginality Quali AleanoTire Recrut no treenbr yearsbetween study offset 
  try(Rrate.QualiA.noT.YO.nbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys))  +
                                                         (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate,
                                                       zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QualiA.noT.YO.nbinom1),silent=T) 
  
  ### Recrut rate Marginality Quanti AleanoTire Recrut no treenbr yearsbetween study offset
  try(Rrate.QuantiA.noT.YO.nbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys))  +
                                                          (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate,
                                                        zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QuantiA.noT.YO.nbinom1),silent=T) 
  
  ## Yearsbetweenstudy not as a predictor
  ### Recrut rate Marginality Quali Fixed Recrut Treenbr fixe no years
  try(Rrate.QualiF.TF.noY.nbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + treeNbrJ +
                                                        Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.mortality.plot.rate,
                                                      zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QualiF.TF.noY.nbinom1),silent=T) 
  
  ### Recrut rate Marginality Quanti Fixed Recrut Treenbr fixe no years
  try(Rrate.QuantiF.TF.noY.nbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) +  treeNbrJ +
                                                         M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.mortality.plot.rate,
                                                       zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QuantiF.TF.noY.nbinom1),silent=T) 
  
  ### Recrut rate Marginality Quali Aleatoire Recrut Treenbr fixe no years 
  try(Rrate.QualiA.TF.noY.nbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) +  treeNbrJ +
                                                        (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate,
                                                      zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QualiA.TF.noY.nbinom1),silent=T) 
  
  ### Recrut rate Marginality Quanti Aleatoire Recrut Treenbr fixe no years
  try(Rrate.QuantiA.TF.noY.nbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) +  treeNbrJ +
                                                         (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate,
                                                       zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QuantiA.TF.noY.nbinom1),silent=T) 
  
  
  
  ## Offset tree 
  ### Recrut rate Marginality Quali Fixed Recrut Treenbr offset no years
  try(Rrate.QualiF.TO.noY.nbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) +  offset(log(treeNbrJ)) +
                                                        Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.mortality.plot.rate,
                                                      zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QualiF.TO.noY.nbinom1),silent=T) 
  
  ### Recrut rate Marginality Quanti Fixed Recrut Treenbr offset no years
  try(Rrate.QuantiF.TO.noY.nbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) +  offset(log(treeNbrJ)) +
                                                         M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.mortality.plot.rate,
                                                       zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QuantiF.TO.noY.nbinom1),silent=T) 
  
  ### Recrut rate Marginality Quali Aleatoire Recrut Treenbr offset no years 
  try(Rrate.QualiA.TO.noY.nbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) +  offset(log(treeNbrJ)) +
                                                        (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate,
                                                      zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QualiA.TO.noY.nbinom1),silent=T) 
  
  ### Recrut rate Marginality Quanti Aleatoire Recrut Treenbr offset no years
  try(Rrate.QuantiA.TO.noY.nbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) +  offset(log(treeNbrJ)) +
                                                         (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate,
                                                       zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QuantiA.TO.noY.nbinom1),silent=T) 
  
  
  ## No tree parameter 
  ### Recrut rate Marginality Quali Fixed Recrut no treenbr no years
  try(Rrate.QualiF.noT.noY.nbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + 
                                                         Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.mortality.plot.rate,
                                                       zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QualiF.noT.noY.nbinom1),silent=T) 
  
  ### Recrut rate Marginality Quanti Fixed Recrut no treenbr no years
  try(Rrate.QuantiF.noT.noY.nbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + 
                                                          M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.mortality.plot.rate,
                                                        zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QuantiF.noT.noY.nbinom1),silent=T) 
  
  ### Recrut rate Marginality Quali AleanoTire Recrut no treenbr no years 
  try(Rrate.QualiA.noT.noY.nbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + 
                                                         (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate,
                                                       zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QualiA.noT.noY.nbinom1),silent=T) 
  
  ### Recrut rate Marginality Quanti AleanoTire Recrut no treenbr no years
  try(Rrate.QuantiA.noT.noY.nbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + 
                                                          (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate,
                                                        zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QuantiA.noT.noY.nbinom1),silent=T) 
  
  
  ##################
  # Idem distrib 2 ##
  ##################
  
  ### Recrut rate Marginality Quali Fixed Recrut Treenbr fixe yearsbetweenstudy fixed
  try(Rrate.QualiF.TF.YF.nbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys + treeNbrJ +
                                                        Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.mortality.plot.rate,
                                                      zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QualiF.TF.YF.nbinom2),silent=T) 
  
  ### Recrut rate Marginality Quanti Fixed Recrut Treenbr fixe yearsbetweenstudy fixed
  try(Rrate.QuantiF.TF.YF.nbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys + treeNbrJ +
                                                         M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.mortality.plot.rate,
                                                       zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QuantiF.TF.YF.nbinom2),silent=T) 
  
  ### Recrut rate Marginality Quali Aleatoire Recrut Treenbr fixe yearsbetweenstudy fixed 
  try(Rrate.QualiA.TF.YF.nbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys + treeNbrJ +
                                                        (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate,
                                                      zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QualiA.TF.YF.nbinom2),silent=T) 
  
  ### Recrut rate Marginality Quanti Aleatoire Recrut Treenbr fixe yearsbetweenstudy fixed
  try(Rrate.QuantiA.TF.YF.nbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys + treeNbrJ +
                                                         (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate,
                                                       zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QuantiA.TF.YF.nbinom2),silent=T) 
  
  
  
  ## Offset tree 
  ### Recrut rate Marginality Quali Fixed Recrut Treenbr offset yearsbetweenstudy fixed
  try(Rrate.QualiF.TO.YF.nbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys + offset(log(treeNbrJ)) +
                                                        Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.mortality.plot.rate,
                                                      zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QualiF.TO.YF.nbinom2),silent=T) 
  
  ### Recrut rate Marginality Quanti Fixed Recrut Treenbr offset yearsbetweenstudy fixed
  try(Rrate.QuantiF.TO.YF.nbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys + offset(log(treeNbrJ)) +
                                                         M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.mortality.plot.rate,
                                                       zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QuantiF.TO.YF.nbinom2),silent=T) 
  
  ### Recrut rate Marginality Quali Aleatoire Recrut Treenbr offset yearsbetweenstudy fixed 
  try(Rrate.QualiA.TO.YF.nbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys + offset(log(treeNbrJ)) +
                                                        (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate,
                                                      zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QualiA.TO.YF.nbinom2),silent=T) 
  
  ### Recrut rate Marginality Quanti Aleatoire Recrut Treenbr offset yearsbetweenstudy fixed
  try(Rrate.QuantiA.TO.YF.nbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys + offset(log(treeNbrJ)) +
                                                         (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate,
                                                       zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QuantiA.TO.YF.nbinom2),silent=T) 
  
  
  ## No tree parameter 
  ### Recrut rate Marginality Quali Fixed Recrut no treenbr yearsbetweenstudy fixed
  try(Rrate.QualiF.noT.YF.nbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys  +
                                                         Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.mortality.plot.rate,
                                                       zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QualiF.noT.YF.nbinom2),silent=T) 
  
  ### Recrut rate Marginality Quanti Fixed Recrut no treenbr yearsbetweenstudy fixed
  try(Rrate.QuantiF.noT.YF.nbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys  +
                                                          M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.mortality.plot.rate,
                                                        zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QuantiF.noT.YF.nbinom2),silent=T) 
  
  ### Recrut rate Marginality Quali AleanoTire Recrut no treenbr yearsbetweenstudy fixed 
  try(Rrate.QualiA.noT.YF.nbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys  +
                                                         (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate,
                                                       zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QualiA.noT.YF.nbinom2),silent=T) 
  
  ### Recrut rate Marginality Quanti AleanoTire Recrut no treenbr yearsbetweenstudy fixed
  try(Rrate.QuantiA.noT.YF.nbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys  +
                                                          (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate,
                                                        zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QuantiA.noT.YF.nbinom2),silent=T) 
  
  
  
  
  #### Years between study as an offset 
  ### Recrut rate Marginality Quali Fixed Recrut Treenbr fixe yearsbetween study offset
  try(Rrate.QualiF.TF.YO.nbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys)) + treeNbrJ +
                                                        Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.mortality.plot.rate,
                                                      zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QualiF.TF.YO.nbinom2),silent=T) 
  
  ### Recrut rate Marginality Quanti Fixed Recrut Treenbr fixe yearsbetween study offset
  try(Rrate.QuantiF.TF.YO.nbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys)) + treeNbrJ +
                                                         M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.mortality.plot.rate,
                                                       zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QuantiF.TF.YO.nbinom2),silent=T) 
  
  ### Recrut rate Marginality Quali Aleatoire Recrut Treenbr fixe yearsbetween study offset 
  try(Rrate.QualiA.TF.YO.nbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys)) + treeNbrJ +
                                                        (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate,
                                                      zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QualiA.TF.YO.nbinom2),silent=T) 
  
  ### Recrut rate Marginality Quanti Aleatoire Recrut Treenbr fixe yearsbetween study offset
  try(Rrate.QuantiA.TF.YO.nbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys)) + treeNbrJ +
                                                         (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate,
                                                       zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QuantiA.TF.YO.nbinom2),silent=T) 
  
  
  
  ## Offset tree 
  ### Recrut rate Marginality Quali Fixed Recrut Treenbr offset yearsbetween study offset
  try(Rrate.QualiF.TO.YO.nbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys)) + offset(log(treeNbrJ)) +
                                                        Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.mortality.plot.rate,
                                                      zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QualiF.TO.YO.nbinom2),silent=T) 
  
  ### Recrut rate Marginality Quanti Fixed Recrut Treenbr offset yearsbetween study offset
  try(Rrate.QuantiF.TO.YO.nbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys)) + offset(log(treeNbrJ)) +
                                                         M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.mortality.plot.rate,
                                                       zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QuantiF.TO.YO.nbinom2),silent=T) 
  
  ### Recrut rate Marginality Quali Aleatoire Recrut Treenbr offset yearsbetween study offset 
  try(Rrate.QualiA.TO.YO.nbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys)) + offset(log(treeNbrJ)) +
                                                        (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate,
                                                      zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QualiA.TO.YO.nbinom2),silent=T) 
  
  ### Recrut rate Marginality Quanti Aleatoire Recrut Treenbr offset yearsbetween study offset
  try(Rrate.QuantiA.TO.YO.nbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys)) + offset(log(treeNbrJ)) +
                                                         (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate,
                                                       zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QuantiA.TO.YO.nbinom2),silent=T) 
  
  
  ## No tree parameter 
  ### Recrut rate Marginality Quali Fixed Recrut no treenbr yearsbetween study offset
  try(Rrate.QualiF.noT.YO.nbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys))  +
                                                         Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.mortality.plot.rate,
                                                       zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QualiF.noT.YO.nbinom2),silent=T) 
  
  ### Recrut rate Marginality Quanti Fixed Recrut no treenbr yearsbetween study offset
  try(Rrate.QuantiF.noT.YO.nbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys))  +
                                                          M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.mortality.plot.rate,
                                                        zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QuantiF.noT.YO.nbinom2),silent=T) 
  
  ### Recrut rate Marginality Quali AleanoTire Recrut no treenbr yearsbetween study offset 
  try(Rrate.QualiA.noT.YO.nbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys))  +
                                                         (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate,
                                                       zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QualiA.noT.YO.nbinom2),silent=T) 
  
  ### Recrut rate Marginality Quanti AleanoTire Recrut no treenbr yearsbetween study offset
  try(Rrate.QuantiA.noT.YO.nbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys))  +
                                                          (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate,
                                                        zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QuantiA.noT.YO.nbinom2),silent=T) 
  
  ## Yearsbetweenstudy not as a predictor
  ### Recrut rate Marginality Quali Fixed Recrut Treenbr fixe no years
  try(Rrate.QualiF.TF.noY.nbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + treeNbrJ +
                                                         Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.mortality.plot.rate,
                                                       zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QualiF.TF.noY.nbinom2),silent=T) 
  
  ### Recrut rate Marginality Quanti Fixed Recrut Treenbr fixe no years
  try(Rrate.QuantiF.TF.noY.nbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) +  treeNbrJ +
                                                          M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.mortality.plot.rate,
                                                        zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QuantiF.TF.noY.nbinom2),silent=T) 
  
  ### Recrut rate Marginality Quali Aleatoire Recrut Treenbr fixe no years 
  try(Rrate.QualiA.TF.noY.nbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) +  treeNbrJ +
                                                         (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate,
                                                       zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QualiA.TF.noY.nbinom2),silent=T) 
  
  ### Recrut rate Marginality Quanti Aleatoire Recrut Treenbr fixe no years
  try(Rrate.QuantiA.TF.noY.nbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) +  treeNbrJ +
                                                          (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate,
                                                        zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QuantiA.TF.noY.nbinom2),silent=T) 
  
  
  
  ## Offset tree 
  ### Recrut rate Marginality Quali Fixed Recrut Treenbr offset no years
  try(Rrate.QualiF.TO.noY.nbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) +  offset(log(treeNbrJ)) +
                                                         Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.mortality.plot.rate,
                                                       zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QualiF.TO.noY.nbinom2),silent=T) 
  
  ### Recrut rate Marginality Quanti Fixed Recrut Treenbr offset no years
  try(Rrate.QuantiF.TO.noY.nbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) +  offset(log(treeNbrJ)) +
                                                          M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.mortality.plot.rate,
                                                        zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QuantiF.TO.noY.nbinom2),silent=T) 
  
  ### Recrut rate Marginality Quali Aleatoire Recrut Treenbr offset no years 
  try(Rrate.QualiA.TO.noY.nbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) +  offset(log(treeNbrJ)) +
                                                         (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate,
                                                       zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QualiA.TO.noY.nbinom2),silent=T) 
  
  ### Recrut rate Marginality Quanti Aleatoire Recrut Treenbr offset no years
  try(Rrate.QuantiA.TO.noY.nbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) +  offset(log(treeNbrJ)) +
                                                          (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate,
                                                        zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QuantiA.TO.noY.nbinom2),silent=T) 
  
  
  ## No tree parameter 
  ### Recrut rate Marginality Quali Fixed Recrut no treenbr no years
  try(Rrate.QualiF.noT.noY.nbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + 
                                                          Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.mortality.plot.rate,
                                                        zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QualiF.noT.noY.nbinom2),silent=T) 
  
  ### Recrut rate Marginality Quanti Fixed Recrut no treenbr no years
  try(Rrate.QuantiF.noT.noY.nbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + 
                                                           M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.mortality.plot.rate,
                                                         zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QuantiF.noT.noY.nbinom2),silent=T) 
  
  ### Recrut rate Marginality Quali AleanoTire Recrut no treenbr no years 
  try(Rrate.QualiA.noT.noY.nbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + 
                                                          (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate,
                                                        zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QualiA.noT.noY.nbinom2),silent=T) 
  
  ### Recrut rate Marginality Quanti AleanoTire Recrut no treenbr no years
  try(Rrate.QuantiA.noT.noY.nbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + 
                                                           (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate,
                                                         zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QuantiA.noT.noY.nbinom2),silent=T) 
  
  #####################
  ### Idem distrib 3 ##
  #####################
  
  ### Recrut rate Marginality Quali Fixed Recrut Treenbr fixe yearsbetweenstudy fixed
  try(Rrate.QualiF.TF.YF.Hnbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys + treeNbrJ +
                                                        Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.mortality.plot.rate,
                                                      zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QualiF.TF.YF.Hnbinom1),silent=T) 
  
  ### Recrut rate Marginality Quanti Fixed Recrut Treenbr fixe yearsbetweenstudy fixed
  try(Rrate.QuantiF.TF.YF.Hnbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys + treeNbrJ +
                                                         M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.mortality.plot.rate,
                                                       zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QuantiF.TF.YF.Hnbinom1),silent=T) 
  
  ### Recrut rate Marginality Quali Aleatoire Recrut Treenbr fixe yearsbetweenstudy fixed 
  try(Rrate.QualiA.TF.YF.Hnbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys + treeNbrJ +
                                                        (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate,
                                                      zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QualiA.TF.YF.Hnbinom1),silent=T) 
  
  ### Recrut rate Marginality Quanti Aleatoire Recrut Treenbr fixe yearsbetweenstudy fixed
  try(Rrate.QuantiA.TF.YF.Hnbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys + treeNbrJ +
                                                         (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate,
                                                       zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QuantiA.TF.YF.Hnbinom1),silent=T) 
  
  
  
  ## Offset tree 
  ### Recrut rate Marginality Quali Fixed Recrut Treenbr offset yearsbetweenstudy fixed
  try(Rrate.QualiF.TO.YF.Hnbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys + offset(log(treeNbrJ)) +
                                                        Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.mortality.plot.rate,
                                                      zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QualiF.TO.YF.Hnbinom1),silent=T) 
  
  ### Recrut rate Marginality Quanti Fixed Recrut Treenbr offset yearsbetweenstudy fixed
  try(Rrate.QuantiF.TO.YF.Hnbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys + offset(log(treeNbrJ)) +
                                                         M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.mortality.plot.rate,
                                                       zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QuantiF.TO.YF.Hnbinom1),silent=T) 
  
  ### Recrut rate Marginality Quali Aleatoire Recrut Treenbr offset yearsbetweenstudy fixed 
  try(Rrate.QualiA.TO.YF.Hnbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys + offset(log(treeNbrJ)) +
                                                        (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate,
                                                      zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QualiA.TO.YF.Hnbinom1),silent=T) 
  
  ### Recrut rate Marginality Quanti Aleatoire Recrut Treenbr offset yearsbetweenstudy fixed
  try(Rrate.QuantiA.TO.YF.Hnbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys + offset(log(treeNbrJ)) +
                                                         (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate,
                                                       zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QuantiA.TO.YF.Hnbinom1),silent=T) 
  
  
  ## No tree parameter 
  ### Recrut rate Marginality Quali Fixed Recrut no treenbr yearsbetweenstudy fixed
  try(Rrate.QualiF.noT.YF.Hnbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys  +
                                                         Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.mortality.plot.rate,
                                                       zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QualiF.noT.YF.Hnbinom1),silent=T) 
  
  ### Recrut rate Marginality Quanti Fixed Recrut no treenbr yearsbetweenstudy fixed
  try(Rrate.QuantiF.noT.YF.Hnbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys  +
                                                          M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.mortality.plot.rate,
                                                        zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QuantiF.noT.YF.Hnbinom1),silent=T) 
  
  ### Recrut rate Marginality Quali AleanoTire Recrut no treenbr yearsbetweenstudy fixed 
  try(Rrate.QualiA.noT.YF.Hnbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys  +
                                                         (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate,
                                                       zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QualiA.noT.YF.Hnbinom1),silent=T) 
  
  ### Recrut rate Marginality Quanti AleanoTire Recrut no treenbr yearsbetweenstudy fixed
  try(Rrate.QuantiA.noT.YF.Hnbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys  +
                                                          (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate,
                                                        zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QuantiA.noT.YF.Hnbinom1),silent=T) 
  
  
  
  
  #### Years between study as an offset 
  ### Recrut rate Marginality Quali Fixed Recrut Treenbr fixe yearsbetween study offset
  try(Rrate.QualiF.TF.YO.Hnbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys)) + treeNbrJ +
                                                        Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.mortality.plot.rate,
                                                      zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QualiF.TF.YO.Hnbinom1),silent=T) 
  
  ### Recrut rate Marginality Quanti Fixed Recrut Treenbr fixe yearsbetween study offset
  try(Rrate.QuantiF.TF.YO.Hnbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys)) + treeNbrJ +
                                                         M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.mortality.plot.rate,
                                                       zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QuantiF.TF.YO.Hnbinom1),silent=T) 
  
  ### Recrut rate Marginality Quali Aleatoire Recrut Treenbr fixe yearsbetween study offset 
  try(Rrate.QualiA.TF.YO.Hnbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys)) + treeNbrJ +
                                                        (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate,
                                                      zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QualiA.TF.YO.Hnbinom1),silent=T) 
  
  ### Recrut rate Marginality Quanti Aleatoire Recrut Treenbr fixe yearsbetween study offset
  try(Rrate.QuantiA.TF.YO.Hnbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys)) + treeNbrJ +
                                                         (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate,
                                                       zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QuantiA.TF.YO.Hnbinom1),silent=T) 
  
  
  
  ## Offset tree 
  ### Recrut rate Marginality Quali Fixed Recrut Treenbr offset yearsbetween study offset
  try(Rrate.QualiF.TO.YO.Hnbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys)) + offset(log(treeNbrJ)) +
                                                        Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.mortality.plot.rate,
                                                      zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QualiF.TO.YO.Hnbinom1),silent=T) 
  
  ### Recrut rate Marginality Quanti Fixed Recrut Treenbr offset yearsbetween study offset
  try(Rrate.QuantiF.TO.YO.Hnbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys)) + offset(log(treeNbrJ)) +
                                                         M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.mortality.plot.rate,
                                                       zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QuantiF.TO.YO.Hnbinom1),silent=T) 
  
  ### Recrut rate Marginality Quali Aleatoire Recrut Treenbr offset yearsbetween study offset 
  try(Rrate.QualiA.TO.YO.Hnbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys)) + offset(log(treeNbrJ)) +
                                                        (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate,
                                                      zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QualiA.TO.YO.Hnbinom1),silent=T) 
  
  ### Recrut rate Marginality Quanti Aleatoire Recrut Treenbr offset yearsbetween study offset
  try(Rrate.QuantiA.TO.YO.Hnbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys)) + offset(log(treeNbrJ)) +
                                                         (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate,
                                                       zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QuantiA.TO.YO.Hnbinom1),silent=T) 
  
  
  ## No tree parameter 
  ### Recrut rate Marginality Quali Fixed Recrut no treenbr yearsbetween study offset
  try(Rrate.QualiF.noT.YO.Hnbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys))  +
                                                         Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.mortality.plot.rate,
                                                       zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QualiF.noT.YO.Hnbinom1),silent=T) 
  
  ### Recrut rate Marginality Quanti Fixed Recrut no treenbr yearsbetween study offset
  try(Rrate.QuantiF.noT.YO.Hnbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys))  +
                                                          M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.mortality.plot.rate,
                                                        zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QuantiF.noT.YO.Hnbinom1),silent=T) 
  
  ### Recrut rate Marginality Quali AleanoTire Recrut no treenbr yearsbetween study offset 
  try(Rrate.QualiA.noT.YO.Hnbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys))  +
                                                         (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate,
                                                       zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QualiA.noT.YO.Hnbinom1),silent=T) 
  
  ### Recrut rate Marginality Quanti AleanoTire Recrut no treenbr yearsbetween study offset
  try(Rrate.QuantiA.noT.YO.Hnbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys))  +
                                                          (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate,
                                                        zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QuantiA.noT.YO.Hnbinom1),silent=T) 
  
  ## Yearsbetweenstudy not as a predictor
  ### Recrut rate Marginality Quali Fixed Recrut Treenbr fixe no years
  try(Rrate.QualiF.TF.noY.Hnbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + treeNbrJ +
                                                         Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.mortality.plot.rate,
                                                       zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QualiF.TF.noY.Hnbinom1),silent=T) 
  
  ### Recrut rate Marginality Quanti Fixed Recrut Treenbr fixe no years
  try(Rrate.QuantiF.TF.noY.Hnbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) +  treeNbrJ +
                                                          M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.mortality.plot.rate,
                                                        zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QuantiF.TF.noY.Hnbinom1),silent=T) 
  
  ### Recrut rate Marginality Quali Aleatoire Recrut Treenbr fixe no years 
  try(Rrate.QualiA.TF.noY.Hnbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) +  treeNbrJ +
                                                         (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate,
                                                       zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QualiA.TF.noY.Hnbinom1),silent=T) 
  
  ### Recrut rate Marginality Quanti Aleatoire Recrut Treenbr fixe no years
  try(Rrate.QuantiA.TF.noY.Hnbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) +  treeNbrJ +
                                                          (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate,
                                                        zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QuantiA.TF.noY.Hnbinom1),silent=T) 
  
  
  
  ## Offset tree 
  ### Recrut rate Marginality Quali Fixed Recrut Treenbr offset no years
  try(Rrate.QualiF.TO.noY.Hnbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) +  offset(log(treeNbrJ)) +
                                                         Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.mortality.plot.rate,
                                                       zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QualiF.TO.noY.Hnbinom1),silent=T) 
  
  ### Recrut rate Marginality Quanti Fixed Recrut Treenbr offset no years
  try(Rrate.QuantiF.TO.noY.Hnbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) +  offset(log(treeNbrJ)) +
                                                          M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.mortality.plot.rate,
                                                        zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QuantiF.TO.noY.Hnbinom1),silent=T) 
  
  ### Recrut rate Marginality Quali Aleatoire Recrut Treenbr offset no years 
  try(Rrate.QualiA.TO.noY.Hnbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) +  offset(log(treeNbrJ)) +
                                                         (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate,
                                                       zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QualiA.TO.noY.Hnbinom1),silent=T) 
  
  ### Recrut rate Marginality Quanti Aleatoire Recrut Treenbr offset no years
  try(Rrate.QuantiA.TO.noY.Hnbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) +  offset(log(treeNbrJ)) +
                                                          (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate,
                                                        zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QuantiA.TO.noY.Hnbinom1),silent=T) 
  
  
  ## No tree parameter 
  ### Recrut rate Marginality Quali Fixed Recrut no treenbr no years
  try(Rrate.QualiF.noT.noY.Hnbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + 
                                                          Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.mortality.plot.rate,
                                                        zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QualiF.noT.noY.Hnbinom1),silent=T) 
  
  ### Recrut rate Marginality Quanti Fixed Recrut no treenbr no years
  try(Rrate.QuantiF.noT.noY.Hnbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + 
                                                           M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.mortality.plot.rate,
                                                         zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QuantiF.noT.noY.Hnbinom1),silent=T) 
  
  ### Recrut rate Marginality Quali AleanoTire Recrut no treenbr no years 
  try(Rrate.QualiA.noT.noY.Hnbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + 
                                                          (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate,
                                                        zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QualiA.noT.noY.Hnbinom1),silent=T) 
  
  ### Recrut rate Marginality Quanti AleanoTire Recrut no treenbr no years
  try(Rrate.QuantiA.noT.noY.Hnbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + 
                                                           (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate,
                                                         zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QuantiA.noT.noY.Hnbinom1),silent=T) 
  
  ##############
  ## Distrib4 ##
  ##############
  
  ### Recrut rate Marginality Quali Fixed Recrut Treenbr fixe yearsbetweenstudy fixed
  try(Rrate.QualiF.TF.YF.Hnbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys + treeNbrJ +
                                                        Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.mortality.plot.rate,
                                                      zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QualiF.TF.YF.Hnbinom2),silent=T) 
  
  ### Recrut rate Marginality Quanti Fixed Recrut Treenbr fixe yearsbetweenstudy fixed
  try(Rrate.QuantiF.TF.YF.Hnbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys + treeNbrJ +
                                                         M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.mortality.plot.rate,
                                                       zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QuantiF.TF.YF.Hnbinom2),silent=T) 
  
  ### Recrut rate Marginality Quali Aleatoire Recrut Treenbr fixe yearsbetweenstudy fixed 
  try(Rrate.QualiA.TF.YF.Hnbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys + treeNbrJ +
                                                        (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate,
                                                      zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QualiA.TF.YF.Hnbinom2),silent=T) 
  
  ### Recrut rate Marginality Quanti Aleatoire Recrut Treenbr fixe yearsbetweenstudy fixed
  try(Rrate.QuantiA.TF.YF.Hnbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys + treeNbrJ +
                                                         (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate,
                                                       zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QuantiA.TF.YF.Hnbinom2),silent=T) 
  
  
  
  ## Offset tree 
  ### Recrut rate Marginality Quali Fixed Recrut Treenbr offset yearsbetweenstudy fixed
  try(Rrate.QualiF.TO.YF.Hnbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys + offset(log(treeNbrJ)) +
                                                        Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.mortality.plot.rate,
                                                      zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QualiF.TO.YF.Hnbinom2),silent=T) 
  
  ### Recrut rate Marginality Quanti Fixed Recrut Treenbr offset yearsbetweenstudy fixed
  try(Rrate.QuantiF.TO.YF.Hnbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys + offset(log(treeNbrJ)) +
                                                         M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.mortality.plot.rate,
                                                       zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QuantiF.TO.YF.Hnbinom2),silent=T) 
  
  ### Recrut rate Marginality Quali Aleatoire Recrut Treenbr offset yearsbetweenstudy fixed 
  try(Rrate.QualiA.TO.YF.Hnbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys + offset(log(treeNbrJ)) +
                                                        (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate,
                                                      zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QualiA.TO.YF.Hnbinom2),silent=T) 
  
  ### Recrut rate Marginality Quanti Aleatoire Recrut Treenbr offset yearsbetweenstudy fixed
  try(Rrate.QuantiA.TO.YF.Hnbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys + offset(log(treeNbrJ)) +
                                                         (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate,
                                                       zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QuantiA.TO.YF.Hnbinom2),silent=T) 
  
  
  ## No tree parameter 
  ### Recrut rate Marginality Quali Fixed Recrut no treenbr yearsbetweenstudy fixed
  try(Rrate.QualiF.noT.YF.Hnbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys  +
                                                         Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.mortality.plot.rate,
                                                       zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QualiF.noT.YF.Hnbinom2),silent=T) 
  
  ### Recrut rate Marginality Quanti Fixed Recrut no treenbr yearsbetweenstudy fixed
  try(Rrate.QuantiF.noT.YF.Hnbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys  +
                                                          M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.mortality.plot.rate,
                                                        zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QuantiF.noT.YF.Hnbinom2),silent=T) 
  
  ### Recrut rate Marginality Quali AleanoTire Recrut no treenbr yearsbetweenstudy fixed 
  try(Rrate.QualiA.noT.YF.Hnbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys  +
                                                         (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate,
                                                       zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QualiA.noT.YF.Hnbinom2),silent=T) 
  
  ### Recrut rate Marginality Quanti AleanoTire Recrut no treenbr yearsbetweenstudy fixed
  try(Rrate.QuantiA.noT.YF.Hnbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys  +
                                                          (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate,
                                                        zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QuantiA.noT.YF.Hnbinom2),silent=T) 
  
  
  
  
  #### Years between study as an offset 
  ### Recrut rate Marginality Quali Fixed Recrut Treenbr fixe yearsbetween study offset
  try(Rrate.QualiF.TF.YO.Hnbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys)) + treeNbrJ +
                                                        Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.mortality.plot.rate,
                                                      zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QualiF.TF.YO.Hnbinom2),silent=T) 
  
  ### Recrut rate Marginality Quanti Fixed Recrut Treenbr fixe yearsbetween study offset
  try(Rrate.QuantiF.TF.YO.Hnbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys)) + treeNbrJ +
                                                         M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.mortality.plot.rate,
                                                       zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QuantiF.TF.YO.Hnbinom2),silent=T) 
  
  ### Recrut rate Marginality Quali Aleatoire Recrut Treenbr fixe yearsbetween study offset 
  try(Rrate.QualiA.TF.YO.Hnbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys)) + treeNbrJ +
                                                        (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate,
                                                      zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QualiA.TF.YO.Hnbinom2),silent=T) 
  
  ### Recrut rate Marginality Quanti Aleatoire Recrut Treenbr fixe yearsbetween study offset
  try(Rrate.QuantiA.TF.YO.Hnbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys)) + treeNbrJ +
                                                         (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate,
                                                       zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QuantiA.TF.YO.Hnbinom2),silent=T) 
  
  
  
  ## Offset tree 
  ### Recrut rate Marginality Quali Fixed Recrut Treenbr offset yearsbetween study offset
  try(Rrate.QualiF.TO.YO.Hnbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys)) + offset(log(treeNbrJ)) +
                                                        Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.mortality.plot.rate,
                                                      zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QualiF.TO.YO.Hnbinom2),silent=T) 
  
  ### Recrut rate Marginality Quanti Fixed Recrut Treenbr offset yearsbetween study offset
  try(Rrate.QuantiF.TO.YO.Hnbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys)) + offset(log(treeNbrJ)) +
                                                         M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.mortality.plot.rate,
                                                       zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QuantiF.TO.YO.Hnbinom2),silent=T) 
  
  ### Recrut rate Marginality Quali Aleatoire Recrut Treenbr offset yearsbetween study offset 
  try(Rrate.QualiA.TO.YO.Hnbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys)) + offset(log(treeNbrJ)) +
                                                        (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate,
                                                      zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QualiA.TO.YO.Hnbinom2),silent=T) 
  
  ### Recrut rate Marginality Quanti Aleatoire Recrut Treenbr offset yearsbetween study offset
  try(Rrate.QuantiA.TO.YO.Hnbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys)) + offset(log(treeNbrJ)) +
                                                         (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate,
                                                       zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QuantiA.TO.YO.Hnbinom2),silent=T) 
  
  
  ## No tree parameter 
  ### Recrut rate Marginality Quali Fixed Recrut no treenbr yearsbetween study offset
  try(Rrate.QualiF.noT.YO.Hnbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys))  +
                                                         Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.mortality.plot.rate,
                                                       zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QualiF.noT.YO.Hnbinom2),silent=T) 
  
  ### Recrut rate Marginality Quanti Fixed Recrut no treenbr yearsbetween study offset
  try(Rrate.QuantiF.noT.YO.Hnbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys))  +
                                                          M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.mortality.plot.rate,
                                                        zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QuantiF.noT.YO.Hnbinom2),silent=T) 
  
  ### Recrut rate Marginality Quali AleanoTire Recrut no treenbr yearsbetween study offset 
  try(Rrate.QualiA.noT.YO.Hnbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys))  +
                                                         (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate,
                                                       zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QualiA.noT.YO.Hnbinom2),silent=T) 
  
  ### Recrut rate Marginality Quanti AleanoTire Recrut no treenbr yearsbetween study offset
  try(Rrate.QuantiA.noT.YO.Hnbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys))  +
                                                          (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate,
                                                        zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QuantiA.noT.YO.Hnbinom2),silent=T) 
  
  ## Yearsbetweenstudy not as a predictor
  ### Recrut rate Marginality Quali Fixed Recrut Treenbr fixe no years
  try(Rrate.QualiF.TF.noY.Hnbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + treeNbrJ +
                                                         Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.mortality.plot.rate,
                                                       zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QualiF.TF.noY.Hnbinom2),silent=T) 
  
  ### Recrut rate Marginality Quanti Fixed Recrut Treenbr fixe no years
  try(Rrate.QuantiF.TF.noY.Hnbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) +  treeNbrJ +
                                                          M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.mortality.plot.rate,
                                                        zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QuantiF.TF.noY.Hnbinom2),silent=T) 
  
  ### Recrut rate Marginality Quali Aleatoire Recrut Treenbr fixe no years 
  try(Rrate.QualiA.TF.noY.Hnbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) +  treeNbrJ +
                                                         (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate,
                                                       zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QualiA.TF.noY.Hnbinom2),silent=T) 
  
  ### Recrut rate Marginality Quanti Aleatoire Recrut Treenbr fixe no years
  try(Rrate.QuantiA.TF.noY.Hnbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) +  treeNbrJ +
                                                          (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate,
                                                        zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QuantiA.TF.noY.Hnbinom2),silent=T) 
  
  
  
  ## Offset tree 
  ### Recrut rate Marginality Quali Fixed Recrut Treenbr offset no years
  try(Rrate.QualiF.TO.noY.Hnbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) +  offset(log(treeNbrJ)) +
                                                         Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.mortality.plot.rate,
                                                       zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QualiF.TO.noY.Hnbinom2),silent=T) 
  
  ### Recrut rate Marginality Quanti Fixed Recrut Treenbr offset no years
  try(Rrate.QuantiF.TO.noY.Hnbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) +  offset(log(treeNbrJ)) +
                                                          M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.mortality.plot.rate,
                                                        zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QuantiF.TO.noY.Hnbinom2),silent=T) 
  
  ### Recrut rate Marginality Quali Aleatoire Recrut Treenbr offset no years 
  try(Rrate.QualiA.TO.noY.Hnbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) +  offset(log(treeNbrJ)) +
                                                         (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate,
                                                       zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QualiA.TO.noY.Hnbinom2),silent=T) 
  
  ### Recrut rate Marginality Quanti Aleatoire Recrut Treenbr offset no years
  try(Rrate.QuantiA.TO.noY.Hnbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) +  offset(log(treeNbrJ)) +
                                                          (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate,
                                                        zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QuantiA.TO.noY.Hnbinom2),silent=T) 
  
  
  ## No tree parameter 
  ### Recrut rate Marginality Quali Fixed Recrut no treenbr no years
  try(Rrate.QualiF.noT.noY.Hnbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + 
                                                          Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.mortality.plot.rate,
                                                        zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QualiF.noT.noY.Hnbinom2),silent=T) 
  
  ### Recrut rate Marginality Quanti Fixed Recrut no treenbr no years
  try(Rrate.QuantiF.noT.noY.Hnbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + 
                                                           M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.mortality.plot.rate,
                                                         zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QuantiF.noT.noY.Hnbinom2),silent=T) 
  
  ### Recrut rate Marginality Quali AleanoTire Recrut no treenbr no years 
  try(Rrate.QualiA.noT.noY.Hnbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + 
                                                          (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate,
                                                        zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QualiA.noT.noY.Hnbinom2),silent=T) 
  
  ### Recrut rate Marginality Quanti AleanoTire Recrut no treenbr no years
  try(Rrate.QuantiA.noT.noY.Hnbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + 
                                                           (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate,
                                                         zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.QuantiA.noT.noY.Hnbinom2),silent=T) 
  
  
  
  
  #########################
  # IDEM but with rate.yr #
  #########################
  
  
  ### Recrut rate.yr.yr Marginality Quali Fixed Recrut Treenbr fixe yearsbetweenstudy fixed
  try(Rrate.yr.yr.QualiF.TF.YF.nbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys + treeNbrJ +
                                                        Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.mortality.plot.rate.yr.yr,
                                                      zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.yr.QualiF.TF.YF.nbinom1),silent=T) 
  
  ### Recrut rate.yr Marginality Quanti Fixed Recrut Treenbr fixe yearsbetweenstudy fixed
  try(Rrate.yr.QuantiF.TF.YF.nbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys + treeNbrJ +
                                                         M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.mortality.plot.rate.yr,
                                                       zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QuantiF.TF.YF.nbinom1),silent=T) 
  
  ### Recrut rate.yr Marginality Quali Aleatoire Recrut Treenbr fixe yearsbetweenstudy fixed 
  try(Rrate.yr.QualiA.TF.YF.nbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys + treeNbrJ +
                                                        (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate.yr,
                                                      zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QualiA.TF.YF.nbinom1),silent=T) 
  
  ### Recrut rate.yr Marginality Quanti Aleatoire Recrut Treenbr fixe yearsbetweenstudy fixed
  try(Rrate.yr.QuantiA.TF.YF.nbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys + treeNbrJ +
                                                         (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate.yr,
                                                       zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QuantiA.TF.YF.nbinom1),silent=T) 
  
  
  
  ## Offset tree 
  ### Recrut rate.yr Marginality Quali Fixed Recrut Treenbr offset yearsbetweenstudy fixed
  try(Rrate.yr.QualiF.TO.YF.nbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys + offset(log(treeNbrJ)) +
                                                        Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.mortality.plot.rate.yr,
                                                      zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QualiF.TO.YF.nbinom1),silent=T) 
  
  ### Recrut rate.yr Marginality Quanti Fixed Recrut Treenbr offset yearsbetweenstudy fixed
  try(Rrate.yr.QuantiF.TO.YF.nbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys + offset(log(treeNbrJ)) +
                                                         M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.mortality.plot.rate.yr,
                                                       zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QuantiF.TO.YF.nbinom1),silent=T) 
  
  ### Recrut rate.yr Marginality Quali Aleatoire Recrut Treenbr offset yearsbetweenstudy fixed 
  try(Rrate.yr.QualiA.TO.YF.nbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys + offset(log(treeNbrJ)) +
                                                        (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate.yr,
                                                      zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QualiA.TO.YF.nbinom1),silent=T) 
  
  ### Recrut rate.yr Marginality Quanti Aleatoire Recrut Treenbr offset yearsbetweenstudy fixed
  try(Rrate.yr.QuantiA.TO.YF.nbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys + offset(log(treeNbrJ)) +
                                                         (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate.yr,
                                                       zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QuantiA.TO.YF.nbinom1),silent=T) 
  
  
  ## No tree parameter 
  ### Recrut rate.yr Marginality Quali Fixed Recrut no treenbr yearsbetweenstudy fixed
  try(Rrate.yr.QualiF.noT.YF.nbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys  +
                                                         Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.mortality.plot.rate.yr,
                                                       zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QualiF.noT.YF.nbinom1),silent=T) 
  
  ### Recrut rate.yr Marginality Quanti Fixed Recrut no treenbr yearsbetweenstudy fixed
  try(Rrate.yr.QuantiF.noT.YF.nbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys  +
                                                          M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.mortality.plot.rate.yr,
                                                        zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QuantiF.noT.YF.nbinom1),silent=T) 
  
  ### Recrut rate.yr Marginality Quali AleanoTire Recrut no treenbr yearsbetweenstudy fixed 
  try(Rrate.yr.QualiA.noT.YF.nbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys  +
                                                         (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate.yr,
                                                       zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QualiA.noT.YF.nbinom1),silent=T) 
  
  ### Recrut rate.yr Marginality Quanti AleanoTire Recrut no treenbr yearsbetweenstudy fixed
  try(Rrate.yr.QuantiA.noT.YF.nbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys  +
                                                          (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate.yr,
                                                        zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QuantiA.noT.YF.nbinom1),silent=T) 
  
  
  
  
  #### Years between study as an offset 
  ### Recrut rate.yr Marginality Quali Fixed Recrut Treenbr fixe yearsbetween study offset
  try(Rrate.yr.QualiF.TF.YO.nbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys)) + treeNbrJ +
                                                        Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.mortality.plot.rate.yr,
                                                      zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QualiF.TF.YO.nbinom1),silent=T) 
  
  ### Recrut rate.yr Marginality Quanti Fixed Recrut Treenbr fixe yearsbetween study offset
  try(Rrate.yr.QuantiF.TF.YO.nbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys)) + treeNbrJ +
                                                         M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.mortality.plot.rate.yr,
                                                       zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QuantiF.TF.YO.nbinom1),silent=T) 
  
  ### Recrut rate.yr Marginality Quali Aleatoire Recrut Treenbr fixe yearsbetween study offset 
  try(Rrate.yr.QualiA.TF.YO.nbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys)) + treeNbrJ +
                                                        (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate.yr,
                                                      zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QualiA.TF.YO.nbinom1),silent=T) 
  
  ### Recrut rate.yr Marginality Quanti Aleatoire Recrut Treenbr fixe yearsbetween study offset
  try(Rrate.yr.QuantiA.TF.YO.nbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys)) + treeNbrJ +
                                                         (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate.yr,
                                                       zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QuantiA.TF.YO.nbinom1),silent=T) 
  
  
  
  ## Offset tree 
  ### Recrut rate.yr Marginality Quali Fixed Recrut Treenbr offset yearsbetween study offset
  try(Rrate.yr.QualiF.TO.YO.nbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys)) + offset(log(treeNbrJ)) +
                                                        Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.mortality.plot.rate.yr,
                                                      zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QualiF.TO.YO.nbinom1),silent=T) 
  
  ### Recrut rate.yr Marginality Quanti Fixed Recrut Treenbr offset yearsbetween study offset
  try(Rrate.yr.QuantiF.TO.YO.nbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys)) + offset(log(treeNbrJ)) +
                                                         M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.mortality.plot.rate.yr,
                                                       zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QuantiF.TO.YO.nbinom1),silent=T) 
  
  ### Recrut rate.yr Marginality Quali Aleatoire Recrut Treenbr offset yearsbetween study offset 
  try(Rrate.yr.QualiA.TO.YO.nbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys)) + offset(log(treeNbrJ)) +
                                                        (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate.yr,
                                                      zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QualiA.TO.YO.nbinom1),silent=T) 
  
  ### Recrut rate.yr Marginality Quanti Aleatoire Recrut Treenbr offset yearsbetween study offset
  try(Rrate.yr.QuantiA.TO.YO.nbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys)) + offset(log(treeNbrJ)) +
                                                         (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate.yr,
                                                       zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QuantiA.TO.YO.nbinom1),silent=T) 
  
  
  ## No tree parameter 
  ### Recrut rate.yr Marginality Quali Fixed Recrut no treenbr yearsbetween study offset
  try(Rrate.yr.QualiF.noT.YO.nbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys))  +
                                                         Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.mortality.plot.rate.yr,
                                                       zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QualiF.noT.YO.nbinom1),silent=T) 
  
  ### Recrut rate.yr Marginality Quanti Fixed Recrut no treenbr yearsbetween study offset
  try(Rrate.yr.QuantiF.noT.YO.nbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys))  +
                                                          M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.mortality.plot.rate.yr,
                                                        zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QuantiF.noT.YO.nbinom1),silent=T) 
  
  ### Recrut rate.yr Marginality Quali AleanoTire Recrut no treenbr yearsbetween study offset 
  try(Rrate.yr.QualiA.noT.YO.nbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys))  +
                                                         (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate.yr,
                                                       zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QualiA.noT.YO.nbinom1),silent=T) 
  
  ### Recrut rate.yr Marginality Quanti AleanoTire Recrut no treenbr yearsbetween study offset
  try(Rrate.yr.QuantiA.noT.YO.nbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys))  +
                                                          (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate.yr,
                                                        zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QuantiA.noT.YO.nbinom1),silent=T) 
  
  ## Yearsbetweenstudy not as a predictor
  ### Recrut rate.yr Marginality Quali Fixed Recrut Treenbr fixe no years
  try(Rrate.yr.QualiF.TF.noY.nbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) + treeNbrJ +
                                                         Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.mortality.plot.rate.yr,
                                                       zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QualiF.TF.noY.nbinom1),silent=T) 
  
  ### Recrut rate.yr Marginality Quanti Fixed Recrut Treenbr fixe no years
  try(Rrate.yr.QuantiF.TF.noY.nbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) +  treeNbrJ +
                                                          M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.mortality.plot.rate.yr,
                                                        zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QuantiF.TF.noY.nbinom1),silent=T) 
  
  ### Recrut rate.yr Marginality Quali Aleatoire Recrut Treenbr fixe no years 
  try(Rrate.yr.QualiA.TF.noY.nbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) +  treeNbrJ +
                                                         (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate.yr,
                                                       zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QualiA.TF.noY.nbinom1),silent=T) 
  
  ### Recrut rate.yr Marginality Quanti Aleatoire Recrut Treenbr fixe no years
  try(Rrate.yr.QuantiA.TF.noY.nbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) +  treeNbrJ +
                                                          (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate.yr,
                                                        zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QuantiA.TF.noY.nbinom1),silent=T) 
  
  
  
  ## Offset tree 
  ### Recrut rate.yr Marginality Quali Fixed Recrut Treenbr offset no years
  try(Rrate.yr.QualiF.TO.noY.nbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) +  offset(log(treeNbrJ)) +
                                                         Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.mortality.plot.rate.yr,
                                                       zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QualiF.TO.noY.nbinom1),silent=T) 
  
  ### Recrut rate.yr Marginality Quanti Fixed Recrut Treenbr offset no years
  try(Rrate.yr.QuantiF.TO.noY.nbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) +  offset(log(treeNbrJ)) +
                                                          M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.mortality.plot.rate.yr,
                                                        zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QuantiF.TO.noY.nbinom1),silent=T) 
  
  ### Recrut rate.yr Marginality Quali Aleatoire Recrut Treenbr offset no years 
  try(Rrate.yr.QualiA.TO.noY.nbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) +  offset(log(treeNbrJ)) +
                                                         (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate.yr,
                                                       zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QualiA.TO.noY.nbinom1),silent=T) 
  
  ### Recrut rate.yr Marginality Quanti Aleatoire Recrut Treenbr offset no years
  try(Rrate.yr.QuantiA.TO.noY.nbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) +  offset(log(treeNbrJ)) +
                                                          (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate.yr,
                                                        zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QuantiA.TO.noY.nbinom1),silent=T) 
  
  
  ## No tree parameter 
  ### Recrut rate.yr Marginality Quali Fixed Recrut no treenbr no years
  try(Rrate.yr.QualiF.noT.noY.nbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) + 
                                                          Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.mortality.plot.rate.yr,
                                                        zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QualiF.noT.noY.nbinom1),silent=T) 
  
  ### Recrut rate.yr Marginality Quanti Fixed Recrut no treenbr no years
  try(Rrate.yr.QuantiF.noT.noY.nbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) + 
                                                           M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.mortality.plot.rate.yr,
                                                         zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QuantiF.noT.noY.nbinom1),silent=T) 
  
  ### Recrut rate.yr Marginality Quali AleanoTire Recrut no treenbr no years 
  try(Rrate.yr.QualiA.noT.noY.nbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) + 
                                                          (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate.yr,
                                                        zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QualiA.noT.noY.nbinom1),silent=T) 
  
  ### Recrut rate.yr Marginality Quanti AleanoTire Recrut no treenbr no years
  try(Rrate.yr.QuantiA.noT.noY.nbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) + 
                                                           (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate.yr,
                                                         zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QuantiA.noT.noY.nbinom1),silent=T) 
  
  
  ##################
  # Idem distrib 2 ##
  ##################
  
  ### Recrut rate.yr Marginality Quali Fixed Recrut Treenbr fixe yearsbetweenstudy fixed
  try(Rrate.yr.QualiF.TF.YF.nbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys + treeNbrJ +
                                                        Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.mortality.plot.rate.yr,
                                                      zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QualiF.TF.YF.nbinom2),silent=T) 
  
  ### Recrut rate.yr Marginality Quanti Fixed Recrut Treenbr fixe yearsbetweenstudy fixed
  try(Rrate.yr.QuantiF.TF.YF.nbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys + treeNbrJ +
                                                         M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.mortality.plot.rate.yr,
                                                       zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QuantiF.TF.YF.nbinom2),silent=T) 
  
  ### Recrut rate.yr Marginality Quali Aleatoire Recrut Treenbr fixe yearsbetweenstudy fixed 
  try(Rrate.yr.QualiA.TF.YF.nbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys + treeNbrJ +
                                                        (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate.yr,
                                                      zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QualiA.TF.YF.nbinom2),silent=T) 
  
  ### Recrut rate.yr Marginality Quanti Aleatoire Recrut Treenbr fixe yearsbetweenstudy fixed
  try(Rrate.yr.QuantiA.TF.YF.nbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys + treeNbrJ +
                                                         (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate.yr,
                                                       zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QuantiA.TF.YF.nbinom2),silent=T) 
  
  
  
  ## Offset tree 
  ### Recrut rate.yr Marginality Quali Fixed Recrut Treenbr offset yearsbetweenstudy fixed
  try(Rrate.yr.QualiF.TO.YF.nbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys + offset(log(treeNbrJ)) +
                                                        Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.mortality.plot.rate.yr,
                                                      zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QualiF.TO.YF.nbinom2),silent=T) 
  
  ### Recrut rate.yr Marginality Quanti Fixed Recrut Treenbr offset yearsbetweenstudy fixed
  try(Rrate.yr.QuantiF.TO.YF.nbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys + offset(log(treeNbrJ)) +
                                                         M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.mortality.plot.rate.yr,
                                                       zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QuantiF.TO.YF.nbinom2),silent=T) 
  
  ### Recrut rate.yr Marginality Quali Aleatoire Recrut Treenbr offset yearsbetweenstudy fixed 
  try(Rrate.yr.QualiA.TO.YF.nbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys + offset(log(treeNbrJ)) +
                                                        (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate.yr,
                                                      zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QualiA.TO.YF.nbinom2),silent=T) 
  
  ### Recrut rate.yr Marginality Quanti Aleatoire Recrut Treenbr offset yearsbetweenstudy fixed
  try(Rrate.yr.QuantiA.TO.YF.nbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys + offset(log(treeNbrJ)) +
                                                         (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate.yr,
                                                       zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QuantiA.TO.YF.nbinom2),silent=T) 
  
  
  ## No tree parameter 
  ### Recrut rate.yr Marginality Quali Fixed Recrut no treenbr yearsbetweenstudy fixed
  try(Rrate.yr.QualiF.noT.YF.nbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys  +
                                                         Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.mortality.plot.rate.yr,
                                                       zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QualiF.noT.YF.nbinom2),silent=T) 
  
  ### Recrut rate.yr Marginality Quanti Fixed Recrut no treenbr yearsbetweenstudy fixed
  try(Rrate.yr.QuantiF.noT.YF.nbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys  +
                                                          M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.mortality.plot.rate.yr,
                                                        zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QuantiF.noT.YF.nbinom2),silent=T) 
  
  ### Recrut rate.yr Marginality Quali AleanoTire Recrut no treenbr yearsbetweenstudy fixed 
  try(Rrate.yr.QualiA.noT.YF.nbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys  +
                                                         (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate.yr,
                                                       zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QualiA.noT.YF.nbinom2),silent=T) 
  
  ### Recrut rate.yr Marginality Quanti AleanoTire Recrut no treenbr yearsbetweenstudy fixed
  try(Rrate.yr.QuantiA.noT.YF.nbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys  +
                                                          (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate.yr,
                                                        zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QuantiA.noT.YF.nbinom2),silent=T) 
  
  
  
  
  #### Years between study as an offset 
  ### Recrut rate.yr Marginality Quali Fixed Recrut Treenbr fixe yearsbetween study offset
  try(Rrate.yr.QualiF.TF.YO.nbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys)) + treeNbrJ +
                                                        Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.mortality.plot.rate.yr,
                                                      zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QualiF.TF.YO.nbinom2),silent=T) 
  
  ### Recrut rate.yr Marginality Quanti Fixed Recrut Treenbr fixe yearsbetween study offset
  try(Rrate.yr.QuantiF.TF.YO.nbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys)) + treeNbrJ +
                                                         M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.mortality.plot.rate.yr,
                                                       zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QuantiF.TF.YO.nbinom2),silent=T) 
  
  ### Recrut rate.yr Marginality Quali Aleatoire Recrut Treenbr fixe yearsbetween study offset 
  try(Rrate.yr.QualiA.TF.YO.nbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys)) + treeNbrJ +
                                                        (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate.yr,
                                                      zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QualiA.TF.YO.nbinom2),silent=T) 
  
  ### Recrut rate.yr Marginality Quanti Aleatoire Recrut Treenbr fixe yearsbetween study offset
  try(Rrate.yr.QuantiA.TF.YO.nbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys)) + treeNbrJ +
                                                         (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate.yr,
                                                       zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QuantiA.TF.YO.nbinom2),silent=T) 
  
  
  
  ## Offset tree 
  ### Recrut rate.yr Marginality Quali Fixed Recrut Treenbr offset yearsbetween study offset
  try(Rrate.yr.QualiF.TO.YO.nbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys)) + offset(log(treeNbrJ)) +
                                                        Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.mortality.plot.rate.yr,
                                                      zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QualiF.TO.YO.nbinom2),silent=T) 
  
  ### Recrut rate.yr Marginality Quanti Fixed Recrut Treenbr offset yearsbetween study offset
  try(Rrate.yr.QuantiF.TO.YO.nbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys)) + offset(log(treeNbrJ)) +
                                                         M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.mortality.plot.rate.yr,
                                                       zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QuantiF.TO.YO.nbinom2),silent=T) 
  
  ### Recrut rate.yr Marginality Quali Aleatoire Recrut Treenbr offset yearsbetween study offset 
  try(Rrate.yr.QualiA.TO.YO.nbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys)) + offset(log(treeNbrJ)) +
                                                        (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate.yr,
                                                      zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QualiA.TO.YO.nbinom2),silent=T) 
  
  ### Recrut rate.yr Marginality Quanti Aleatoire Recrut Treenbr offset yearsbetween study offset
  try(Rrate.yr.QuantiA.TO.YO.nbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys)) + offset(log(treeNbrJ)) +
                                                         (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate.yr,
                                                       zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QuantiA.TO.YO.nbinom2),silent=T) 
  
  
  ## No tree parameter 
  ### Recrut rate.yr Marginality Quali Fixed Recrut no treenbr yearsbetween study offset
  try(Rrate.yr.QualiF.noT.YO.nbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys))  +
                                                         Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.mortality.plot.rate.yr,
                                                       zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QualiF.noT.YO.nbinom2),silent=T) 
  
  ### Recrut rate.yr Marginality Quanti Fixed Recrut no treenbr yearsbetween study offset
  try(Rrate.yr.QuantiF.noT.YO.nbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys))  +
                                                          M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.mortality.plot.rate.yr,
                                                        zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QuantiF.noT.YO.nbinom2),silent=T) 
  
  ### Recrut rate.yr Marginality Quali AleanoTire Recrut no treenbr yearsbetween study offset 
  try(Rrate.yr.QualiA.noT.YO.nbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys))  +
                                                         (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate.yr,
                                                       zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QualiA.noT.YO.nbinom2),silent=T) 
  
  ### Recrut rate.yr Marginality Quanti AleanoTire Recrut no treenbr yearsbetween study offset
  try(Rrate.yr.QuantiA.noT.YO.nbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys))  +
                                                          (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate.yr,
                                                        zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QuantiA.noT.YO.nbinom2),silent=T) 
  
  ## Yearsbetweenstudy not as a predictor
  ### Recrut rate.yr Marginality Quali Fixed Recrut Treenbr fixe no years
  try(Rrate.yr.QualiF.TF.noY.nbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) + treeNbrJ +
                                                         Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.mortality.plot.rate.yr,
                                                       zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QualiF.TF.noY.nbinom2),silent=T) 
  
  ### Recrut rate.yr Marginality Quanti Fixed Recrut Treenbr fixe no years
  try(Rrate.yr.QuantiF.TF.noY.nbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) +  treeNbrJ +
                                                          M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.mortality.plot.rate.yr,
                                                        zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QuantiF.TF.noY.nbinom2),silent=T) 
  
  ### Recrut rate.yr Marginality Quali Aleatoire Recrut Treenbr fixe no years 
  try(Rrate.yr.QualiA.TF.noY.nbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) +  treeNbrJ +
                                                         (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate.yr,
                                                       zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QualiA.TF.noY.nbinom2),silent=T) 
  
  ### Recrut rate.yr Marginality Quanti Aleatoire Recrut Treenbr fixe no years
  try(Rrate.yr.QuantiA.TF.noY.nbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) +  treeNbrJ +
                                                          (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate.yr,
                                                        zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QuantiA.TF.noY.nbinom2),silent=T) 
  
  
  
  ## Offset tree 
  ### Recrut rate.yr Marginality Quali Fixed Recrut Treenbr offset no years
  try(Rrate.yr.QualiF.TO.noY.nbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) +  offset(log(treeNbrJ)) +
                                                         Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.mortality.plot.rate.yr,
                                                       zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QualiF.TO.noY.nbinom2),silent=T) 
  
  ### Recrut rate.yr Marginality Quanti Fixed Recrut Treenbr offset no years
  try(Rrate.yr.QuantiF.TO.noY.nbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) +  offset(log(treeNbrJ)) +
                                                          M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.mortality.plot.rate.yr,
                                                        zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QuantiF.TO.noY.nbinom2),silent=T) 
  
  ### Recrut rate.yr Marginality Quali Aleatoire Recrut Treenbr offset no years 
  try(Rrate.yr.QualiA.TO.noY.nbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) +  offset(log(treeNbrJ)) +
                                                         (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate.yr,
                                                       zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QualiA.TO.noY.nbinom2),silent=T) 
  
  ### Recrut rate.yr Marginality Quanti Aleatoire Recrut Treenbr offset no years
  try(Rrate.yr.QuantiA.TO.noY.nbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) +  offset(log(treeNbrJ)) +
                                                          (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate.yr,
                                                        zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QuantiA.TO.noY.nbinom2),silent=T) 
  
  
  ## No tree parameter 
  ### Recrut rate.yr Marginality Quali Fixed Recrut no treenbr no years
  try(Rrate.yr.QualiF.noT.noY.nbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) + 
                                                          Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.mortality.plot.rate.yr,
                                                        zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QualiF.noT.noY.nbinom2),silent=T) 
  
  ### Recrut rate.yr Marginality Quanti Fixed Recrut no treenbr no years
  try(Rrate.yr.QuantiF.noT.noY.nbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) + 
                                                           M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.mortality.plot.rate.yr,
                                                         zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QuantiF.noT.noY.nbinom2),silent=T) 
  
  ### Recrut rate.yr Marginality Quali AleanoTire Recrut no treenbr no years 
  try(Rrate.yr.QualiA.noT.noY.nbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) + 
                                                          (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate.yr,
                                                        zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QualiA.noT.noY.nbinom2),silent=T) 
  
  ### Recrut rate.yr Marginality Quanti AleanoTire Recrut no treenbr no years
  try(Rrate.yr.QuantiA.noT.noY.nbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) + 
                                                           (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate.yr,
                                                         zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QuantiA.noT.noY.nbinom2),silent=T) 
  
  #####################
  ### Idem distrib 3 ##
  #####################
  
  ### Recrut rate.yr Marginality Quali Fixed Recrut Treenbr fixe yearsbetweenstudy fixed
  try(Rrate.yr.QualiF.TF.YF.Hnbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys + treeNbrJ +
                                                         Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.mortality.plot.rate.yr,
                                                       zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QualiF.TF.YF.Hnbinom1),silent=T) 
  
  ### Recrut rate.yr Marginality Quanti Fixed Recrut Treenbr fixe yearsbetweenstudy fixed
  try(Rrate.yr.QuantiF.TF.YF.Hnbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys + treeNbrJ +
                                                          M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.mortality.plot.rate.yr,
                                                        zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QuantiF.TF.YF.Hnbinom1),silent=T) 
  
  ### Recrut rate.yr Marginality Quali Aleatoire Recrut Treenbr fixe yearsbetweenstudy fixed 
  try(Rrate.yr.QualiA.TF.YF.Hnbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys + treeNbrJ +
                                                         (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate.yr,
                                                       zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QualiA.TF.YF.Hnbinom1),silent=T) 
  
  ### Recrut rate.yr Marginality Quanti Aleatoire Recrut Treenbr fixe yearsbetweenstudy fixed
  try(Rrate.yr.QuantiA.TF.YF.Hnbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys + treeNbrJ +
                                                          (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate.yr,
                                                        zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QuantiA.TF.YF.Hnbinom1),silent=T) 
  
  
  
  ## Offset tree 
  ### Recrut rate.yr Marginality Quali Fixed Recrut Treenbr offset yearsbetweenstudy fixed
  try(Rrate.yr.QualiF.TO.YF.Hnbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys + offset(log(treeNbrJ)) +
                                                         Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.mortality.plot.rate.yr,
                                                       zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QualiF.TO.YF.Hnbinom1),silent=T) 
  
  ### Recrut rate.yr Marginality Quanti Fixed Recrut Treenbr offset yearsbetweenstudy fixed
  try(Rrate.yr.QuantiF.TO.YF.Hnbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys + offset(log(treeNbrJ)) +
                                                          M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.mortality.plot.rate.yr,
                                                        zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QuantiF.TO.YF.Hnbinom1),silent=T) 
  
  ### Recrut rate.yr Marginality Quali Aleatoire Recrut Treenbr offset yearsbetweenstudy fixed 
  try(Rrate.yr.QualiA.TO.YF.Hnbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys + offset(log(treeNbrJ)) +
                                                         (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate.yr,
                                                       zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QualiA.TO.YF.Hnbinom1),silent=T) 
  
  ### Recrut rate.yr Marginality Quanti Aleatoire Recrut Treenbr offset yearsbetweenstudy fixed
  try(Rrate.yr.QuantiA.TO.YF.Hnbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys + offset(log(treeNbrJ)) +
                                                          (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate.yr,
                                                        zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QuantiA.TO.YF.Hnbinom1),silent=T) 
  
  
  ## No tree parameter 
  ### Recrut rate.yr Marginality Quali Fixed Recrut no treenbr yearsbetweenstudy fixed
  try(Rrate.yr.QualiF.noT.YF.Hnbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys  +
                                                          Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.mortality.plot.rate.yr,
                                                        zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QualiF.noT.YF.Hnbinom1),silent=T) 
  
  ### Recrut rate.yr Marginality Quanti Fixed Recrut no treenbr yearsbetweenstudy fixed
  try(Rrate.yr.QuantiF.noT.YF.Hnbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys  +
                                                           M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.mortality.plot.rate.yr,
                                                         zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QuantiF.noT.YF.Hnbinom1),silent=T) 
  
  ### Recrut rate.yr Marginality Quali AleanoTire Recrut no treenbr yearsbetweenstudy fixed 
  try(Rrate.yr.QualiA.noT.YF.Hnbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys  +
                                                          (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate.yr,
                                                        zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QualiA.noT.YF.Hnbinom1),silent=T) 
  
  ### Recrut rate.yr Marginality Quanti AleanoTire Recrut no treenbr yearsbetweenstudy fixed
  try(Rrate.yr.QuantiA.noT.YF.Hnbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys  +
                                                           (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate.yr,
                                                         zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QuantiA.noT.YF.Hnbinom1),silent=T) 
  
  
  
  
  #### Years between study as an offset 
  ### Recrut rate.yr Marginality Quali Fixed Recrut Treenbr fixe yearsbetween study offset
  try(Rrate.yr.QualiF.TF.YO.Hnbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys)) + treeNbrJ +
                                                         Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.mortality.plot.rate.yr,
                                                       zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QualiF.TF.YO.Hnbinom1),silent=T) 
  
  ### Recrut rate.yr Marginality Quanti Fixed Recrut Treenbr fixe yearsbetween study offset
  try(Rrate.yr.QuantiF.TF.YO.Hnbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys)) + treeNbrJ +
                                                          M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.mortality.plot.rate.yr,
                                                        zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QuantiF.TF.YO.Hnbinom1),silent=T) 
  
  ### Recrut rate.yr Marginality Quali Aleatoire Recrut Treenbr fixe yearsbetween study offset 
  try(Rrate.yr.QualiA.TF.YO.Hnbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys)) + treeNbrJ +
                                                         (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate.yr,
                                                       zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QualiA.TF.YO.Hnbinom1),silent=T) 
  
  ### Recrut rate.yr Marginality Quanti Aleatoire Recrut Treenbr fixe yearsbetween study offset
  try(Rrate.yr.QuantiA.TF.YO.Hnbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys)) + treeNbrJ +
                                                          (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate.yr,
                                                        zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QuantiA.TF.YO.Hnbinom1),silent=T) 
  
  
  
  ## Offset tree 
  ### Recrut rate.yr Marginality Quali Fixed Recrut Treenbr offset yearsbetween study offset
  try(Rrate.yr.QualiF.TO.YO.Hnbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys)) + offset(log(treeNbrJ)) +
                                                         Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.mortality.plot.rate.yr,
                                                       zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QualiF.TO.YO.Hnbinom1),silent=T) 
  
  ### Recrut rate.yr Marginality Quanti Fixed Recrut Treenbr offset yearsbetween study offset
  try(Rrate.yr.QuantiF.TO.YO.Hnbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys)) + offset(log(treeNbrJ)) +
                                                          M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.mortality.plot.rate.yr,
                                                        zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QuantiF.TO.YO.Hnbinom1),silent=T) 
  
  ### Recrut rate.yr Marginality Quali Aleatoire Recrut Treenbr offset yearsbetween study offset 
  try(Rrate.yr.QualiA.TO.YO.Hnbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys)) + offset(log(treeNbrJ)) +
                                                         (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate.yr,
                                                       zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QualiA.TO.YO.Hnbinom1),silent=T) 
  
  ### Recrut rate.yr Marginality Quanti Aleatoire Recrut Treenbr offset yearsbetween study offset
  try(Rrate.yr.QuantiA.TO.YO.Hnbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys)) + offset(log(treeNbrJ)) +
                                                          (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate.yr,
                                                        zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QuantiA.TO.YO.Hnbinom1),silent=T) 
  
  
  ## No tree parameter 
  ### Recrut rate.yr Marginality Quali Fixed Recrut no treenbr yearsbetween study offset
  try(Rrate.yr.QualiF.noT.YO.Hnbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys))  +
                                                          Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.mortality.plot.rate.yr,
                                                        zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QualiF.noT.YO.Hnbinom1),silent=T) 
  
  ### Recrut rate.yr Marginality Quanti Fixed Recrut no treenbr yearsbetween study offset
  try(Rrate.yr.QuantiF.noT.YO.Hnbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys))  +
                                                           M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.mortality.plot.rate.yr,
                                                         zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QuantiF.noT.YO.Hnbinom1),silent=T) 
  
  ### Recrut rate.yr Marginality Quali AleanoTire Recrut no treenbr yearsbetween study offset 
  try(Rrate.yr.QualiA.noT.YO.Hnbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys))  +
                                                          (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate.yr,
                                                        zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QualiA.noT.YO.Hnbinom1),silent=T) 
  
  ### Recrut rate.yr Marginality Quanti AleanoTire Recrut no treenbr yearsbetween study offset
  try(Rrate.yr.QuantiA.noT.YO.Hnbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys))  +
                                                           (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate.yr,
                                                         zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QuantiA.noT.YO.Hnbinom1),silent=T) 
  
  ## Yearsbetweenstudy not as a predictor
  ### Recrut rate.yr Marginality Quali Fixed Recrut Treenbr fixe no years
  try(Rrate.yr.QualiF.TF.noY.Hnbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) + treeNbrJ +
                                                          Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.mortality.plot.rate.yr,
                                                        zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QualiF.TF.noY.Hnbinom1),silent=T) 
  
  ### Recrut rate.yr Marginality Quanti Fixed Recrut Treenbr fixe no years
  try(Rrate.yr.QuantiF.TF.noY.Hnbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) +  treeNbrJ +
                                                           M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.mortality.plot.rate.yr,
                                                         zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QuantiF.TF.noY.Hnbinom1),silent=T) 
  
  ### Recrut rate.yr Marginality Quali Aleatoire Recrut Treenbr fixe no years 
  try(Rrate.yr.QualiA.TF.noY.Hnbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) +  treeNbrJ +
                                                          (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate.yr,
                                                        zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QualiA.TF.noY.Hnbinom1),silent=T) 
  
  ### Recrut rate.yr Marginality Quanti Aleatoire Recrut Treenbr fixe no years
  try(Rrate.yr.QuantiA.TF.noY.Hnbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) +  treeNbrJ +
                                                           (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate.yr,
                                                         zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QuantiA.TF.noY.Hnbinom1),silent=T) 
  
  
  
  ## Offset tree 
  ### Recrut rate.yr Marginality Quali Fixed Recrut Treenbr offset no years
  try(Rrate.yr.QualiF.TO.noY.Hnbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) +  offset(log(treeNbrJ)) +
                                                          Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.mortality.plot.rate.yr,
                                                        zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QualiF.TO.noY.Hnbinom1),silent=T) 
  
  ### Recrut rate.yr Marginality Quanti Fixed Recrut Treenbr offset no years
  try(Rrate.yr.QuantiF.TO.noY.Hnbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) +  offset(log(treeNbrJ)) +
                                                           M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.mortality.plot.rate.yr,
                                                         zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QuantiF.TO.noY.Hnbinom1),silent=T) 
  
  ### Recrut rate.yr Marginality Quali Aleatoire Recrut Treenbr offset no years 
  try(Rrate.yr.QualiA.TO.noY.Hnbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) +  offset(log(treeNbrJ)) +
                                                          (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate.yr,
                                                        zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QualiA.TO.noY.Hnbinom1),silent=T) 
  
  ### Recrut rate.yr Marginality Quanti Aleatoire Recrut Treenbr offset no years
  try(Rrate.yr.QuantiA.TO.noY.Hnbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) +  offset(log(treeNbrJ)) +
                                                           (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate.yr,
                                                         zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QuantiA.TO.noY.Hnbinom1),silent=T) 
  
  
  ## No tree parameter 
  ### Recrut rate.yr Marginality Quali Fixed Recrut no treenbr no years
  try(Rrate.yr.QualiF.noT.noY.Hnbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) + 
                                                           Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.mortality.plot.rate.yr,
                                                         zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QualiF.noT.noY.Hnbinom1),silent=T) 
  
  ### Recrut rate.yr Marginality Quanti Fixed Recrut no treenbr no years
  try(Rrate.yr.QuantiF.noT.noY.Hnbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) + 
                                                            M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.mortality.plot.rate.yr,
                                                          zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QuantiF.noT.noY.Hnbinom1),silent=T) 
  
  ### Recrut rate.yr Marginality Quali AleanoTire Recrut no treenbr no years 
  try(Rrate.yr.QualiA.noT.noY.Hnbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) + 
                                                           (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate.yr,
                                                         zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QualiA.noT.noY.Hnbinom1),silent=T) 
  
  ### Recrut rate.yr Marginality Quanti AleanoTire Recrut no treenbr no years
  try(Rrate.yr.QuantiA.noT.noY.Hnbinom1 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) + 
                                                            (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate.yr,
                                                          zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QuantiA.noT.noY.Hnbinom1),silent=T) 
  
  ##############
  ## Distrib4 ##
  ##############
  
  ### Recrut rate.yr Marginality Quali Fixed Recrut Treenbr fixe yearsbetweenstudy fixed
  try(Rrate.yr.QualiF.TF.YF.Hnbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys + treeNbrJ +
                                                         Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.mortality.plot.rate.yr,
                                                       zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QualiF.TF.YF.Hnbinom2),silent=T) 
  
  ### Recrut rate.yr Marginality Quanti Fixed Recrut Treenbr fixe yearsbetweenstudy fixed
  try(Rrate.yr.QuantiF.TF.YF.Hnbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys + treeNbrJ +
                                                          M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.mortality.plot.rate.yr,
                                                        zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QuantiF.TF.YF.Hnbinom2),silent=T) 
  
  ### Recrut rate.yr Marginality Quali Aleatoire Recrut Treenbr fixe yearsbetweenstudy fixed 
  try(Rrate.yr.QualiA.TF.YF.Hnbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys + treeNbrJ +
                                                         (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate.yr,
                                                       zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QualiA.TF.YF.Hnbinom2),silent=T) 
  
  ### Recrut rate.yr Marginality Quanti Aleatoire Recrut Treenbr fixe yearsbetweenstudy fixed
  try(Rrate.yr.QuantiA.TF.YF.Hnbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys + treeNbrJ +
                                                          (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate.yr,
                                                        zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QuantiA.TF.YF.Hnbinom2),silent=T) 
  
  
  
  ## Offset tree 
  ### Recrut rate.yr Marginality Quali Fixed Recrut Treenbr offset yearsbetweenstudy fixed
  try(Rrate.yr.QualiF.TO.YF.Hnbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys + offset(log(treeNbrJ)) +
                                                         Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.mortality.plot.rate.yr,
                                                       zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QualiF.TO.YF.Hnbinom2),silent=T) 
  
  ### Recrut rate.yr Marginality Quanti Fixed Recrut Treenbr offset yearsbetweenstudy fixed
  try(Rrate.yr.QuantiF.TO.YF.Hnbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys + offset(log(treeNbrJ)) +
                                                          M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.mortality.plot.rate.yr,
                                                        zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QuantiF.TO.YF.Hnbinom2),silent=T) 
  
  ### Recrut rate.yr Marginality Quali Aleatoire Recrut Treenbr offset yearsbetweenstudy fixed 
  try(Rrate.yr.QualiA.TO.YF.Hnbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys + offset(log(treeNbrJ)) +
                                                         (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate.yr,
                                                       zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QualiA.TO.YF.Hnbinom2),silent=T) 
  
  ### Recrut rate.yr Marginality Quanti Aleatoire Recrut Treenbr offset yearsbetweenstudy fixed
  try(Rrate.yr.QuantiA.TO.YF.Hnbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys + offset(log(treeNbrJ)) +
                                                          (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate.yr,
                                                        zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QuantiA.TO.YF.Hnbinom2),silent=T) 
  
  
  ## No tree parameter 
  ### Recrut rate.yr Marginality Quali Fixed Recrut no treenbr yearsbetweenstudy fixed
  try(Rrate.yr.QualiF.noT.YF.Hnbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys  +
                                                          Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.mortality.plot.rate.yr,
                                                        zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QualiF.noT.YF.Hnbinom2),silent=T) 
  
  ### Recrut rate.yr Marginality Quanti Fixed Recrut no treenbr yearsbetweenstudy fixed
  try(Rrate.yr.QuantiF.noT.YF.Hnbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys  +
                                                           M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.mortality.plot.rate.yr,
                                                         zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QuantiF.noT.YF.Hnbinom2),silent=T) 
  
  ### Recrut rate.yr Marginality Quali AleanoTire Recrut no treenbr yearsbetweenstudy fixed 
  try(Rrate.yr.QualiA.noT.YF.Hnbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys  +
                                                          (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate.yr,
                                                        zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QualiA.noT.YF.Hnbinom2),silent=T) 
  
  ### Recrut rate.yr Marginality Quanti AleanoTire Recrut no treenbr yearsbetweenstudy fixed
  try(Rrate.yr.QuantiA.noT.YF.Hnbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys  +
                                                           (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate.yr,
                                                         zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QuantiA.noT.YF.Hnbinom2),silent=T) 
  
  
  
  
  #### Years between study as an offset 
  ### Recrut rate.yr Marginality Quali Fixed Recrut Treenbr fixe yearsbetween study offset
  try(Rrate.yr.QualiF.TF.YO.Hnbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys)) + treeNbrJ +
                                                         Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.mortality.plot.rate.yr,
                                                       zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QualiF.TF.YO.Hnbinom2),silent=T) 
  
  ### Recrut rate.yr Marginality Quanti Fixed Recrut Treenbr fixe yearsbetween study offset
  try(Rrate.yr.QuantiF.TF.YO.Hnbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys)) + treeNbrJ +
                                                          M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.mortality.plot.rate.yr,
                                                        zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QuantiF.TF.YO.Hnbinom2),silent=T) 
  
  ### Recrut rate.yr Marginality Quali Aleatoire Recrut Treenbr fixe yearsbetween study offset 
  try(Rrate.yr.QualiA.TF.YO.Hnbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys)) + treeNbrJ +
                                                         (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate.yr,
                                                       zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QualiA.TF.YO.Hnbinom2),silent=T) 
  
  ### Recrut rate.yr Marginality Quanti Aleatoire Recrut Treenbr fixe yearsbetween study offset
  try(Rrate.yr.QuantiA.TF.YO.Hnbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys)) + treeNbrJ +
                                                          (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate.yr,
                                                        zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QuantiA.TF.YO.Hnbinom2),silent=T) 
  
  
  
  ## Offset tree 
  ### Recrut rate.yr Marginality Quali Fixed Recrut Treenbr offset yearsbetween study offset
  try(Rrate.yr.QualiF.TO.YO.Hnbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys)) + offset(log(treeNbrJ)) +
                                                         Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.mortality.plot.rate.yr,
                                                       zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QualiF.TO.YO.Hnbinom2),silent=T) 
  
  ### Recrut rate.yr Marginality Quanti Fixed Recrut Treenbr offset yearsbetween study offset
  try(Rrate.yr.QuantiF.TO.YO.Hnbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys)) + offset(log(treeNbrJ)) +
                                                          M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.mortality.plot.rate.yr,
                                                        zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QuantiF.TO.YO.Hnbinom2),silent=T) 
  
  ### Recrut rate.yr Marginality Quali Aleatoire Recrut Treenbr offset yearsbetween study offset 
  try(Rrate.yr.QualiA.TO.YO.Hnbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys)) + offset(log(treeNbrJ)) +
                                                         (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate.yr,
                                                       zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QualiA.TO.YO.Hnbinom2),silent=T) 
  
  ### Recrut rate.yr Marginality Quanti Aleatoire Recrut Treenbr offset yearsbetween study offset
  try(Rrate.yr.QuantiA.TO.YO.Hnbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys)) + offset(log(treeNbrJ)) +
                                                          (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate.yr,
                                                        zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QuantiA.TO.YO.Hnbinom2),silent=T) 
  
  
  ## No tree parameter 
  ### Recrut rate.yr Marginality Quali Fixed Recrut no treenbr yearsbetween study offset
  try(Rrate.yr.QualiF.noT.YO.Hnbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys))  +
                                                          Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.mortality.plot.rate.yr,
                                                        zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QualiF.noT.YO.Hnbinom2),silent=T) 
  
  ### Recrut rate.yr Marginality Quanti Fixed Recrut no treenbr yearsbetween study offset
  try(Rrate.yr.QuantiF.noT.YO.Hnbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys))  +
                                                           M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.mortality.plot.rate.yr,
                                                         zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QuantiF.noT.YO.Hnbinom2),silent=T) 
  
  ### Recrut rate.yr Marginality Quali AleanoTire Recrut no treenbr yearsbetween study offset 
  try(Rrate.yr.QualiA.noT.YO.Hnbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys))  +
                                                          (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate.yr,
                                                        zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QualiA.noT.YO.Hnbinom2),silent=T) 
  
  ### Recrut rate.yr Marginality Quanti AleanoTire Recrut no treenbr yearsbetween study offset
  try(Rrate.yr.QuantiA.noT.YO.Hnbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys))  +
                                                           (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate.yr,
                                                         zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QuantiA.noT.YO.Hnbinom2),silent=T) 
  
  ## Yearsbetweenstudy not as a predictor
  ### Recrut rate.yr Marginality Quali Fixed Recrut Treenbr fixe no years
  try(Rrate.yr.QualiF.TF.noY.Hnbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) + treeNbrJ +
                                                          Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.mortality.plot.rate.yr,
                                                        zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QualiF.TF.noY.Hnbinom2),silent=T) 
  
  ### Recrut rate.yr Marginality Quanti Fixed Recrut Treenbr fixe no years
  try(Rrate.yr.QuantiF.TF.noY.Hnbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) +  treeNbrJ +
                                                           M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.mortality.plot.rate.yr,
                                                         zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QuantiF.TF.noY.Hnbinom2),silent=T) 
  
  ### Recrut rate.yr Marginality Quali Aleatoire Recrut Treenbr fixe no years 
  try(Rrate.yr.QualiA.TF.noY.Hnbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) +  treeNbrJ +
                                                          (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate.yr,
                                                        zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QualiA.TF.noY.Hnbinom2),silent=T) 
  
  ### Recrut rate.yr Marginality Quanti Aleatoire Recrut Treenbr fixe no years
  try(Rrate.yr.QuantiA.TF.noY.Hnbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) +  treeNbrJ +
                                                           (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate.yr,
                                                         zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QuantiA.TF.noY.Hnbinom2),silent=T) 
  
  
  
  ## Offset tree 
  ### Recrut rate.yr Marginality Quali Fixed Recrut Treenbr offset no years
  try(Rrate.yr.QualiF.TO.noY.Hnbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) +  offset(log(treeNbrJ)) +
                                                          Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.mortality.plot.rate.yr,
                                                        zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QualiF.TO.noY.Hnbinom2),silent=T) 
  
  ### Recrut rate.yr Marginality Quanti Fixed Recrut Treenbr offset no years
  try(Rrate.yr.QuantiF.TO.noY.Hnbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) +  offset(log(treeNbrJ)) +
                                                           M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.mortality.plot.rate.yr,
                                                         zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QuantiF.TO.noY.Hnbinom2),silent=T) 
  
  ### Recrut rate.yr Marginality Quali Aleatoire Recrut Treenbr offset no years 
  try(Rrate.yr.QualiA.TO.noY.Hnbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) +  offset(log(treeNbrJ)) +
                                                          (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate.yr,
                                                        zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QualiA.TO.noY.Hnbinom2),silent=T) 
  
  ### Recrut rate.yr Marginality Quanti Aleatoire Recrut Treenbr offset no years
  try(Rrate.yr.QuantiA.TO.noY.Hnbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) +  offset(log(treeNbrJ)) +
                                                           (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate.yr,
                                                         zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QuantiA.TO.noY.Hnbinom2),silent=T) 
  
  
  ## No tree parameter 
  ### Recrut rate.yr Marginality Quali Fixed Recrut no treenbr no years
  try(Rrate.yr.QualiF.noT.noY.Hnbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) + 
                                                           Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.mortality.plot.rate.yr,
                                                         zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QualiF.noT.noY.Hnbinom2),silent=T) 
  
  ### Recrut rate.yr Marginality Quanti Fixed Recrut no treenbr no years
  try(Rrate.yr.QuantiF.noT.noY.Hnbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) + 
                                                            M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.mortality.plot.rate.yr,
                                                          zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QuantiF.noT.noY.Hnbinom2),silent=T) 
  
  ### Recrut rate.yr Marginality Quali AleanoTire Recrut no treenbr no years 
  try(Rrate.yr.QualiA.noT.noY.Hnbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) + 
                                                           (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate.yr,
                                                         zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QualiA.noT.noY.Hnbinom2),silent=T) 
  
  ### Recrut rate.yr Marginality Quanti AleanoTire Recrut no treenbr no years
  try(Rrate.yr.QuantiA.noT.noY.Hnbinom2 <- eval_fork(glmmTMB(sp.recruitment.plot.rate.yr.2 ~ dbhJ.IMall.plot.mean + (1|country) + 
                                                            (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.mortality.plot.rate.yr,
                                                          zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Rrate.yr.QuantiA.noT.noY.Hnbinom2),silent=T) 
  
  #########################
  ##### Mortality tot #####
  ########################
  
  ### Recrut rate Marginality Quali Fixed Recrut Treenbr fixe yearsbetweenstudy fixed
  try(Mrate.QualiF.TF.YF.nbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys + treeNbrJ +
                                                        Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.recruitment.plot.rate.2,
                                                      zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QualiF.TF.YF.nbinom1),silent=T) 
  
  ### Recrut rate Marginality Quanti Fixed Recrut Treenbr fixe yearsbetweenstudy fixed
  try(Mrate.QuantiF.TF.YF.nbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys + treeNbrJ +
                                                         M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.recruitment.plot.rate.2,
                                                       zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QuantiF.TF.YF.nbinom1),silent=T) 
  
  ### Recrut rate Marginality Quali Aleatoire Recrut Treenbr fixe yearsbetweenstudy fixed 
  try(Mrate.QualiA.TF.YF.nbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys + treeNbrJ +
                                                        (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.2,
                                                      zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QualiA.TF.YF.nbinom1),silent=T) 
  
  ### Recrut rate Marginality Quanti Aleatoire Recrut Treenbr fixe yearsbetweenstudy fixed
  try(Mrate.QuantiA.TF.YF.nbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys + treeNbrJ +
                                                         (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.2,
                                                       zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QuantiA.TF.YF.nbinom1),silent=T) 
  
  
  
  ## Offset tree 
  ### Recrut rate Marginality Quali Fixed Recrut Treenbr offset yearsbetweenstudy fixed
  try(Mrate.QualiF.TO.YF.nbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys + offset(log(treeNbrJ)) +
                                                        Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.recruitment.plot.rate.2,
                                                      zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QualiF.TO.YF.nbinom1),silent=T) 
  
  ### Recrut rate Marginality Quanti Fixed Recrut Treenbr offset yearsbetweenstudy fixed
  try(Mrate.QuantiF.TO.YF.nbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys + offset(log(treeNbrJ)) +
                                                         M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.recruitment.plot.rate.2,
                                                       zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QuantiF.TO.YF.nbinom1),silent=T) 
  
  ### Recrut rate Marginality Quali Aleatoire Recrut Treenbr offset yearsbetweenstudy fixed 
  try(Mrate.QualiA.TO.YF.nbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys + offset(log(treeNbrJ)) +
                                                        (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.2,
                                                      zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QualiA.TO.YF.nbinom1),silent=T) 
  
  ### Recrut rate Marginality Quanti Aleatoire Recrut Treenbr offset yearsbetweenstudy fixed
  try(Mrate.QuantiA.TO.YF.nbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys + offset(log(treeNbrJ)) +
                                                         (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.2,
                                                       zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QuantiA.TO.YF.nbinom1),silent=T) 
  
  
  ## No tree parameter 
  ### Recrut rate Marginality Quali Fixed Recrut no treenbr yearsbetweenstudy fixed
  try(Mrate.QualiF.noT.YF.nbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys  +
                                                         Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.recruitment.plot.rate.2,
                                                       zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QualiF.noT.YF.nbinom1),silent=T) 
  
  ### Recrut rate Marginality Quanti Fixed Recrut no treenbr yearsbetweenstudy fixed
  try(Mrate.QuantiF.noT.YF.nbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys  +
                                                          M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.recruitment.plot.rate.2,
                                                        zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QuantiF.noT.YF.nbinom1),silent=T) 
  
  ### Recrut rate Marginality Quali AleanoTire Recrut no treenbr yearsbetweenstudy fixed 
  try(Mrate.QualiA.noT.YF.nbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys  +
                                                         (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.2,
                                                       zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QualiA.noT.YF.nbinom1),silent=T) 
  
  ### Recrut rate Marginality Quanti AleanoTire Recrut no treenbr yearsbetweenstudy fixed
  try(Mrate.QuantiA.noT.YF.nbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys  +
                                                          (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.2,
                                                        zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QuantiA.noT.YF.nbinom1),silent=T) 
  
  
  
  
  #### Years between study as an offset 
  ### Recrut rate Marginality Quali Fixed Recrut Treenbr fixe yearsbetween study offset
  try(Mrate.QualiF.TF.YO.nbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys)) + treeNbrJ +
                                                        Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.recruitment.plot.rate.2,
                                                      zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QualiF.TF.YO.nbinom1),silent=T) 
  
  ### Recrut rate Marginality Quanti Fixed Recrut Treenbr fixe yearsbetween study offset
  try(Mrate.QuantiF.TF.YO.nbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys)) + treeNbrJ +
                                                         M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.recruitment.plot.rate.2,
                                                       zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QuantiF.TF.YO.nbinom1),silent=T) 
  
  ### Recrut rate Marginality Quali Aleatoire Recrut Treenbr fixe yearsbetween study offset 
  try(Mrate.QualiA.TF.YO.nbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys)) + treeNbrJ +
                                                        (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.2,
                                                      zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QualiA.TF.YO.nbinom1),silent=T) 
  
  ### Recrut rate Marginality Quanti Aleatoire Recrut Treenbr fixe yearsbetween study offset
  try(Mrate.QuantiA.TF.YO.nbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys)) + treeNbrJ +
                                                         (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.2,
                                                       zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QuantiA.TF.YO.nbinom1),silent=T) 
  
  
  
  ## Offset tree 
  ### Recrut rate Marginality Quali Fixed Recrut Treenbr offset yearsbetween study offset
  try(Mrate.QualiF.TO.YO.nbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys)) + offset(log(treeNbrJ)) +
                                                        Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.recruitment.plot.rate.2,
                                                      zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QualiF.TO.YO.nbinom1),silent=T) 
  
  ### Recrut rate Marginality Quanti Fixed Recrut Treenbr offset yearsbetween study offset
  try(Mrate.QuantiF.TO.YO.nbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys)) + offset(log(treeNbrJ)) +
                                                         M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.recruitment.plot.rate.2,
                                                       zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QuantiF.TO.YO.nbinom1),silent=T) 
  
  ### Recrut rate Marginality Quali Aleatoire Recrut Treenbr offset yearsbetween study offset 
  try(Mrate.QualiA.TO.YO.nbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys)) + offset(log(treeNbrJ)) +
                                                        (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.2,
                                                      zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QualiA.TO.YO.nbinom1),silent=T) 
  
  ### Recrut rate Marginality Quanti Aleatoire Recrut Treenbr offset yearsbetween study offset
  try(Mrate.QuantiA.TO.YO.nbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys)) + offset(log(treeNbrJ)) +
                                                         (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.2,
                                                       zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QuantiA.TO.YO.nbinom1),silent=T) 
  
  
  ## No tree parameter 
  ### Recrut rate Marginality Quali Fixed Recrut no treenbr yearsbetween study offset
  try(Mrate.QualiF.noT.YO.nbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys))  +
                                                         Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.recruitment.plot.rate.2,
                                                       zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QualiF.noT.YO.nbinom1),silent=T) 
  
  ### Recrut rate Marginality Quanti Fixed Recrut no treenbr yearsbetween study offset
  try(Mrate.QuantiF.noT.YO.nbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys))  +
                                                          M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.recruitment.plot.rate.2,
                                                        zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QuantiF.noT.YO.nbinom1),silent=T) 
  
  ### Recrut rate Marginality Quali AleanoTire Recrut no treenbr yearsbetween study offset 
  try(Mrate.QualiA.noT.YO.nbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys))  +
                                                         (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.2,
                                                       zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QualiA.noT.YO.nbinom1),silent=T) 
  
  ### Recrut rate Marginality Quanti AleanoTire Recrut no treenbr yearsbetween study offset
  try(Mrate.QuantiA.noT.YO.nbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys))  +
                                                          (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.2,
                                                        zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QuantiA.noT.YO.nbinom1),silent=T) 
  
  ## Yearsbetweenstudy not as a predictor
  ### Recrut rate Marginality Quali Fixed Recrut Treenbr fixe no years
  try(Mrate.QualiF.TF.noY.nbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + treeNbrJ +
                                                         Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.recruitment.plot.rate.2,
                                                       zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QualiF.TF.noY.nbinom1),silent=T) 
  
  ### Recrut rate Marginality Quanti Fixed Recrut Treenbr fixe no years
  try(Mrate.QuantiF.TF.noY.nbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) +  treeNbrJ +
                                                          M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.recruitment.plot.rate.2,
                                                        zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QuantiF.TF.noY.nbinom1),silent=T) 
  
  ### Recrut rate Marginality Quali Aleatoire Recrut Treenbr fixe no years 
  try(Mrate.QualiA.TF.noY.nbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) +  treeNbrJ +
                                                         (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.2,
                                                       zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QualiA.TF.noY.nbinom1),silent=T) 
  
  ### Recrut rate Marginality Quanti Aleatoire Recrut Treenbr fixe no years
  try(Mrate.QuantiA.TF.noY.nbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) +  treeNbrJ +
                                                          (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.2,
                                                        zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QuantiA.TF.noY.nbinom1),silent=T) 
  
  
  
  ## Offset tree 
  ### Recrut rate Marginality Quali Fixed Recrut Treenbr offset no years
  try(Mrate.QualiF.TO.noY.nbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) +  offset(log(treeNbrJ)) +
                                                         Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.recruitment.plot.rate.2,
                                                       zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QualiF.TO.noY.nbinom1),silent=T) 
  
  ### Recrut rate Marginality Quanti Fixed Recrut Treenbr offset no years
  try(Mrate.QuantiF.TO.noY.nbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) +  offset(log(treeNbrJ)) +
                                                          M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.recruitment.plot.rate.2,
                                                        zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QuantiF.TO.noY.nbinom1),silent=T) 
  
  ### Recrut rate Marginality Quali Aleatoire Recrut Treenbr offset no years 
  try(Mrate.QualiA.TO.noY.nbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) +  offset(log(treeNbrJ)) +
                                                         (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.2,
                                                       zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QualiA.TO.noY.nbinom1),silent=T) 
  
  ### Recrut rate Marginality Quanti Aleatoire Recrut Treenbr offset no years
  try(Mrate.QuantiA.TO.noY.nbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) +  offset(log(treeNbrJ)) +
                                                          (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.2,
                                                        zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QuantiA.TO.noY.nbinom1),silent=T) 
  
  
  ## No tree parameter 
  ### Recrut rate Marginality Quali Fixed Recrut no treenbr no years
  try(Mrate.QualiF.noT.noY.nbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + 
                                                          Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.recruitment.plot.rate.2,
                                                        zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QualiF.noT.noY.nbinom1),silent=T) 
  
  ### Recrut rate Marginality Quanti Fixed Recrut no treenbr no years
  try(Mrate.QuantiF.noT.noY.nbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + 
                                                           M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.recruitment.plot.rate.2,
                                                         zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QuantiF.noT.noY.nbinom1),silent=T) 
  
  ### Recrut rate Marginality Quali AleanoTire Recrut no treenbr no years 
  try(Mrate.QualiA.noT.noY.nbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + 
                                                          (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.2,
                                                        zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QualiA.noT.noY.nbinom1),silent=T) 
  
  ### Recrut rate Marginality Quanti AleanoTire Recrut no treenbr no years
  try(Mrate.QuantiA.noT.noY.nbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + 
                                                           (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.2,
                                                         zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QuantiA.noT.noY.nbinom1),silent=T) 
  
  
  ##################
  # Idem distrib 2 ##
  ##################
  
  ### Recrut rate Marginality Quali Fixed Recrut Treenbr fixe yearsbetweenstudy fixed
  try(Mrate.QualiF.TF.YF.nbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys + treeNbrJ +
                                                        Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.recruitment.plot.rate.2,
                                                      zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QualiF.TF.YF.nbinom2),silent=T) 
  
  ### Recrut rate Marginality Quanti Fixed Recrut Treenbr fixe yearsbetweenstudy fixed
  try(Mrate.QuantiF.TF.YF.nbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys + treeNbrJ +
                                                         M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.recruitment.plot.rate.2,
                                                       zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QuantiF.TF.YF.nbinom2),silent=T) 
  
  ### Recrut rate Marginality Quali Aleatoire Recrut Treenbr fixe yearsbetweenstudy fixed 
  try(Mrate.QualiA.TF.YF.nbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys + treeNbrJ +
                                                        (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.2,
                                                      zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QualiA.TF.YF.nbinom2),silent=T) 
  
  ### Recrut rate Marginality Quanti Aleatoire Recrut Treenbr fixe yearsbetweenstudy fixed
  try(Mrate.QuantiA.TF.YF.nbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys + treeNbrJ +
                                                         (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.2,
                                                       zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QuantiA.TF.YF.nbinom2),silent=T) 
  
  
  
  ## Offset tree 
  ### Recrut rate Marginality Quali Fixed Recrut Treenbr offset yearsbetweenstudy fixed
  try(Mrate.QualiF.TO.YF.nbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys + offset(log(treeNbrJ)) +
                                                        Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.recruitment.plot.rate.2,
                                                      zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QualiF.TO.YF.nbinom2),silent=T) 
  
  ### Recrut rate Marginality Quanti Fixed Recrut Treenbr offset yearsbetweenstudy fixed
  try(Mrate.QuantiF.TO.YF.nbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys + offset(log(treeNbrJ)) +
                                                         M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.recruitment.plot.rate.2,
                                                       zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QuantiF.TO.YF.nbinom2),silent=T) 
  
  ### Recrut rate Marginality Quali Aleatoire Recrut Treenbr offset yearsbetweenstudy fixed 
  try(Mrate.QualiA.TO.YF.nbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys + offset(log(treeNbrJ)) +
                                                        (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.2,
                                                      zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QualiA.TO.YF.nbinom2),silent=T) 
  
  ### Recrut rate Marginality Quanti Aleatoire Recrut Treenbr offset yearsbetweenstudy fixed
  try(Mrate.QuantiA.TO.YF.nbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys + offset(log(treeNbrJ)) +
                                                         (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.2,
                                                       zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QuantiA.TO.YF.nbinom2),silent=T) 
  
  
  ## No tree parameter 
  ### Recrut rate Marginality Quali Fixed Recrut no treenbr yearsbetweenstudy fixed
  try(Mrate.QualiF.noT.YF.nbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys  +
                                                         Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.recruitment.plot.rate.2,
                                                       zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QualiF.noT.YF.nbinom2),silent=T) 
  
  ### Recrut rate Marginality Quanti Fixed Recrut no treenbr yearsbetweenstudy fixed
  try(Mrate.QuantiF.noT.YF.nbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys  +
                                                          M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.recruitment.plot.rate.2,
                                                        zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QuantiF.noT.YF.nbinom2),silent=T) 
  
  ### Recrut rate Marginality Quali AleanoTire Recrut no treenbr yearsbetweenstudy fixed 
  try(Mrate.QualiA.noT.YF.nbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys  +
                                                         (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.2,
                                                       zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QualiA.noT.YF.nbinom2),silent=T) 
  
  ### Recrut rate Marginality Quanti AleanoTire Recrut no treenbr yearsbetweenstudy fixed
  try(Mrate.QuantiA.noT.YF.nbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys  +
                                                          (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.2,
                                                        zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QuantiA.noT.YF.nbinom2),silent=T) 
  
  
  
  
  #### Years between study as an offset 
  ### Recrut rate Marginality Quali Fixed Recrut Treenbr fixe yearsbetween study offset
  try(Mrate.QualiF.TF.YO.nbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys)) + treeNbrJ +
                                                        Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.recruitment.plot.rate.2,
                                                      zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QualiF.TF.YO.nbinom2),silent=T) 
  
  ### Recrut rate Marginality Quanti Fixed Recrut Treenbr fixe yearsbetween study offset
  try(Mrate.QuantiF.TF.YO.nbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys)) + treeNbrJ +
                                                         M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.recruitment.plot.rate.2,
                                                       zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QuantiF.TF.YO.nbinom2),silent=T) 
  
  ### Recrut rate Marginality Quali Aleatoire Recrut Treenbr fixe yearsbetween study offset 
  try(Mrate.QualiA.TF.YO.nbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys)) + treeNbrJ +
                                                        (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.2,
                                                      zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QualiA.TF.YO.nbinom2),silent=T) 
  
  ### Recrut rate Marginality Quanti Aleatoire Recrut Treenbr fixe yearsbetween study offset
  try(Mrate.QuantiA.TF.YO.nbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys)) + treeNbrJ +
                                                         (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.2,
                                                       zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QuantiA.TF.YO.nbinom2),silent=T) 
  
  
  
  ## Offset tree 
  ### Recrut rate Marginality Quali Fixed Recrut Treenbr offset yearsbetween study offset
  try(Mrate.QualiF.TO.YO.nbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys)) + offset(log(treeNbrJ)) +
                                                        Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.recruitment.plot.rate.2,
                                                      zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QualiF.TO.YO.nbinom2),silent=T) 
  
  ### Recrut rate Marginality Quanti Fixed Recrut Treenbr offset yearsbetween study offset
  try(Mrate.QuantiF.TO.YO.nbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys)) + offset(log(treeNbrJ)) +
                                                         M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.recruitment.plot.rate.2,
                                                       zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QuantiF.TO.YO.nbinom2),silent=T) 
  
  ### Recrut rate Marginality Quali Aleatoire Recrut Treenbr offset yearsbetween study offset 
  try(Mrate.QualiA.TO.YO.nbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys)) + offset(log(treeNbrJ)) +
                                                        (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.2,
                                                      zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QualiA.TO.YO.nbinom2),silent=T) 
  
  ### Recrut rate Marginality Quanti Aleatoire Recrut Treenbr offset yearsbetween study offset
  try(Mrate.QuantiA.TO.YO.nbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys)) + offset(log(treeNbrJ)) +
                                                         (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.2,
                                                       zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QuantiA.TO.YO.nbinom2),silent=T) 
  
  
  ## No tree parameter 
  ### Recrut rate Marginality Quali Fixed Recrut no treenbr yearsbetween study offset
  try(Mrate.QualiF.noT.YO.nbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys))  +
                                                         Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.recruitment.plot.rate.2,
                                                       zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QualiF.noT.YO.nbinom2),silent=T) 
  
  ### Recrut rate Marginality Quanti Fixed Recrut no treenbr yearsbetween study offset
  try(Mrate.QuantiF.noT.YO.nbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys))  +
                                                          M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.recruitment.plot.rate.2,
                                                        zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QuantiF.noT.YO.nbinom2),silent=T) 
  
  ### Recrut rate Marginality Quali AleanoTire Recrut no treenbr yearsbetween study offset 
  try(Mrate.QualiA.noT.YO.nbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys))  +
                                                         (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.2,
                                                       zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QualiA.noT.YO.nbinom2),silent=T) 
  
  ### Recrut rate Marginality Quanti AleanoTire Recrut no treenbr yearsbetween study offset
  try(Mrate.QuantiA.noT.YO.nbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys))  +
                                                          (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.2,
                                                        zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QuantiA.noT.YO.nbinom2),silent=T) 
  
  ## Yearsbetweenstudy not as a predictor
  ### Recrut rate Marginality Quali Fixed Recrut Treenbr fixe no years
  try(Mrate.QualiF.TF.noY.nbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + treeNbrJ +
                                                         Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.recruitment.plot.rate.2,
                                                       zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QualiF.TF.noY.nbinom2),silent=T) 
  
  ### Recrut rate Marginality Quanti Fixed Recrut Treenbr fixe no years
  try(Mrate.QuantiF.TF.noY.nbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) +  treeNbrJ +
                                                          M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.recruitment.plot.rate.2,
                                                        zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QuantiF.TF.noY.nbinom2),silent=T) 
  
  ### Recrut rate Marginality Quali Aleatoire Recrut Treenbr fixe no years 
  try(Mrate.QualiA.TF.noY.nbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) +  treeNbrJ +
                                                         (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.2,
                                                       zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QualiA.TF.noY.nbinom2),silent=T) 
  
  ### Recrut rate Marginality Quanti Aleatoire Recrut Treenbr fixe no years
  try(Mrate.QuantiA.TF.noY.nbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) +  treeNbrJ +
                                                          (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.2,
                                                        zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QuantiA.TF.noY.nbinom2),silent=T) 
  
  
  
  ## Offset tree 
  ### Recrut rate Marginality Quali Fixed Recrut Treenbr offset no years
  try(Mrate.QualiF.TO.noY.nbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) +  offset(log(treeNbrJ)) +
                                                         Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.recruitment.plot.rate.2,
                                                       zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QualiF.TO.noY.nbinom2),silent=T) 
  
  ### Recrut rate Marginality Quanti Fixed Recrut Treenbr offset no years
  try(Mrate.QuantiF.TO.noY.nbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) +  offset(log(treeNbrJ)) +
                                                          M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.recruitment.plot.rate.2,
                                                        zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QuantiF.TO.noY.nbinom2),silent=T) 
  
  ### Recrut rate Marginality Quali Aleatoire Recrut Treenbr offset no years 
  try(Mrate.QualiA.TO.noY.nbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) +  offset(log(treeNbrJ)) +
                                                         (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.2,
                                                       zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QualiA.TO.noY.nbinom2),silent=T) 
  
  ### Recrut rate Marginality Quanti Aleatoire Recrut Treenbr offset no years
  try(Mrate.QuantiA.TO.noY.nbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) +  offset(log(treeNbrJ)) +
                                                          (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.2,
                                                        zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QuantiA.TO.noY.nbinom2),silent=T) 
  
  
  ## No tree parameter 
  ### Recrut rate Marginality Quali Fixed Recrut no treenbr no years
  try(Mrate.QualiF.noT.noY.nbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + 
                                                          Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.recruitment.plot.rate.2,
                                                        zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QualiF.noT.noY.nbinom2),silent=T) 
  
  ### Recrut rate Marginality Quanti Fixed Recrut no treenbr no years
  try(Mrate.QuantiF.noT.noY.nbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + 
                                                           M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.recruitment.plot.rate.2,
                                                         zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QuantiF.noT.noY.nbinom2),silent=T) 
  
  ### Recrut rate Marginality Quali AleanoTire Recrut no treenbr no years 
  try(Mrate.QualiA.noT.noY.nbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + 
                                                          (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.2,
                                                        zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QualiA.noT.noY.nbinom2),silent=T) 
  
  ### Recrut rate Marginality Quanti AleanoTire Recrut no treenbr no years
  try(Mrate.QuantiA.noT.noY.nbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + 
                                                           (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.2,
                                                         zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QuantiA.noT.noY.nbinom2),silent=T) 
  
  #####################
  ### Idem distrib 3 ##
  #####################
  
  ### Recrut rate Marginality Quali Fixed Recrut Treenbr fixe yearsbetweenstudy fixed
  try(Mrate.QualiF.TF.YF.Hnbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys + treeNbrJ +
                                                         Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.recruitment.plot.rate.2,
                                                       zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QualiF.TF.YF.Hnbinom1),silent=T) 
  
  ### Recrut rate Marginality Quanti Fixed Recrut Treenbr fixe yearsbetweenstudy fixed
  try(Mrate.QuantiF.TF.YF.Hnbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys + treeNbrJ +
                                                          M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.recruitment.plot.rate.2,
                                                        zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QuantiF.TF.YF.Hnbinom1),silent=T) 
  
  ### Recrut rate Marginality Quali Aleatoire Recrut Treenbr fixe yearsbetweenstudy fixed 
  try(Mrate.QualiA.TF.YF.Hnbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys + treeNbrJ +
                                                         (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.2,
                                                       zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QualiA.TF.YF.Hnbinom1),silent=T) 
  
  ### Recrut rate Marginality Quanti Aleatoire Recrut Treenbr fixe yearsbetweenstudy fixed
  try(Mrate.QuantiA.TF.YF.Hnbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys + treeNbrJ +
                                                          (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.2,
                                                        zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QuantiA.TF.YF.Hnbinom1),silent=T) 
  
  
  
  ## Offset tree 
  ### Recrut rate Marginality Quali Fixed Recrut Treenbr offset yearsbetweenstudy fixed
  try(Mrate.QualiF.TO.YF.Hnbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys + offset(log(treeNbrJ)) +
                                                         Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.recruitment.plot.rate.2,
                                                       zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QualiF.TO.YF.Hnbinom1),silent=T) 
  
  ### Recrut rate Marginality Quanti Fixed Recrut Treenbr offset yearsbetweenstudy fixed
  try(Mrate.QuantiF.TO.YF.Hnbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys + offset(log(treeNbrJ)) +
                                                          M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.recruitment.plot.rate.2,
                                                        zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QuantiF.TO.YF.Hnbinom1),silent=T) 
  
  ### Recrut rate Marginality Quali Aleatoire Recrut Treenbr offset yearsbetweenstudy fixed 
  try(Mrate.QualiA.TO.YF.Hnbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys + offset(log(treeNbrJ)) +
                                                         (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.2,
                                                       zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QualiA.TO.YF.Hnbinom1),silent=T) 
  
  ### Recrut rate Marginality Quanti Aleatoire Recrut Treenbr offset yearsbetweenstudy fixed
  try(Mrate.QuantiA.TO.YF.Hnbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys + offset(log(treeNbrJ)) +
                                                          (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.2,
                                                        zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QuantiA.TO.YF.Hnbinom1),silent=T) 
  
  
  ## No tree parameter 
  ### Recrut rate Marginality Quali Fixed Recrut no treenbr yearsbetweenstudy fixed
  try(Mrate.QualiF.noT.YF.Hnbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys  +
                                                          Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.recruitment.plot.rate.2,
                                                        zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QualiF.noT.YF.Hnbinom1),silent=T) 
  
  ### Recrut rate Marginality Quanti Fixed Recrut no treenbr yearsbetweenstudy fixed
  try(Mrate.QuantiF.noT.YF.Hnbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys  +
                                                           M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.recruitment.plot.rate.2,
                                                         zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QuantiF.noT.YF.Hnbinom1),silent=T) 
  
  ### Recrut rate Marginality Quali AleanoTire Recrut no treenbr yearsbetweenstudy fixed 
  try(Mrate.QualiA.noT.YF.Hnbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys  +
                                                          (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.2,
                                                        zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QualiA.noT.YF.Hnbinom1),silent=T) 
  
  ### Recrut rate Marginality Quanti AleanoTire Recrut no treenbr yearsbetweenstudy fixed
  try(Mrate.QuantiA.noT.YF.Hnbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys  +
                                                           (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.2,
                                                         zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QuantiA.noT.YF.Hnbinom1),silent=T) 
  
  
  
  
  #### Years between study as an offset 
  ### Recrut rate Marginality Quali Fixed Recrut Treenbr fixe yearsbetween study offset
  try(Mrate.QualiF.TF.YO.Hnbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys)) + treeNbrJ +
                                                         Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.recruitment.plot.rate.2,
                                                       zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QualiF.TF.YO.Hnbinom1),silent=T) 
  
  ### Recrut rate Marginality Quanti Fixed Recrut Treenbr fixe yearsbetween study offset
  try(Mrate.QuantiF.TF.YO.Hnbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys)) + treeNbrJ +
                                                          M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.recruitment.plot.rate.2,
                                                        zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QuantiF.TF.YO.Hnbinom1),silent=T) 
  
  ### Recrut rate Marginality Quali Aleatoire Recrut Treenbr fixe yearsbetween study offset 
  try(Mrate.QualiA.TF.YO.Hnbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys)) + treeNbrJ +
                                                         (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.2,
                                                       zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QualiA.TF.YO.Hnbinom1),silent=T) 
  
  ### Recrut rate Marginality Quanti Aleatoire Recrut Treenbr fixe yearsbetween study offset
  try(Mrate.QuantiA.TF.YO.Hnbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys)) + treeNbrJ +
                                                          (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.2,
                                                        zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QuantiA.TF.YO.Hnbinom1),silent=T) 
  
  
  
  ## Offset tree 
  ### Recrut rate Marginality Quali Fixed Recrut Treenbr offset yearsbetween study offset
  try(Mrate.QualiF.TO.YO.Hnbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys)) + offset(log(treeNbrJ)) +
                                                         Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.recruitment.plot.rate.2,
                                                       zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QualiF.TO.YO.Hnbinom1),silent=T) 
  
  ### Recrut rate Marginality Quanti Fixed Recrut Treenbr offset yearsbetween study offset
  try(Mrate.QuantiF.TO.YO.Hnbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys)) + offset(log(treeNbrJ)) +
                                                          M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.recruitment.plot.rate.2,
                                                        zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QuantiF.TO.YO.Hnbinom1),silent=T) 
  
  ### Recrut rate Marginality Quali Aleatoire Recrut Treenbr offset yearsbetween study offset 
  try(Mrate.QualiA.TO.YO.Hnbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys)) + offset(log(treeNbrJ)) +
                                                         (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.2,
                                                       zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QualiA.TO.YO.Hnbinom1),silent=T) 
  
  ### Recrut rate Marginality Quanti Aleatoire Recrut Treenbr offset yearsbetween study offset
  try(Mrate.QuantiA.TO.YO.Hnbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys)) + offset(log(treeNbrJ)) +
                                                          (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.2,
                                                        zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QuantiA.TO.YO.Hnbinom1),silent=T) 
  
  
  ## No tree parameter 
  ### Recrut rate Marginality Quali Fixed Recrut no treenbr yearsbetween study offset
  try(Mrate.QualiF.noT.YO.Hnbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys))  +
                                                          Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.recruitment.plot.rate.2,
                                                        zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QualiF.noT.YO.Hnbinom1),silent=T) 
  
  ### Recrut rate Marginality Quanti Fixed Recrut no treenbr yearsbetween study offset
  try(Mrate.QuantiF.noT.YO.Hnbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys))  +
                                                           M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.recruitment.plot.rate.2,
                                                         zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QuantiF.noT.YO.Hnbinom1),silent=T) 
  
  ### Recrut rate Marginality Quali AleanoTire Recrut no treenbr yearsbetween study offset 
  try(Mrate.QualiA.noT.YO.Hnbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys))  +
                                                          (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.2,
                                                        zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QualiA.noT.YO.Hnbinom1),silent=T) 
  
  ### Recrut rate Marginality Quanti AleanoTire Recrut no treenbr yearsbetween study offset
  try(Mrate.QuantiA.noT.YO.Hnbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys))  +
                                                           (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.2,
                                                         zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QuantiA.noT.YO.Hnbinom1),silent=T) 
  
  ## Yearsbetweenstudy not as a predictor
  ### Recrut rate Marginality Quali Fixed Recrut Treenbr fixe no years
  try(Mrate.QualiF.TF.noY.Hnbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + treeNbrJ +
                                                          Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.recruitment.plot.rate.2,
                                                        zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QualiF.TF.noY.Hnbinom1),silent=T) 
  
  ### Recrut rate Marginality Quanti Fixed Recrut Treenbr fixe no years
  try(Mrate.QuantiF.TF.noY.Hnbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) +  treeNbrJ +
                                                           M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.recruitment.plot.rate.2,
                                                         zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QuantiF.TF.noY.Hnbinom1),silent=T) 
  
  ### Recrut rate Marginality Quali Aleatoire Recrut Treenbr fixe no years 
  try(Mrate.QualiA.TF.noY.Hnbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) +  treeNbrJ +
                                                          (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.2,
                                                        zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QualiA.TF.noY.Hnbinom1),silent=T) 
  
  ### Recrut rate Marginality Quanti Aleatoire Recrut Treenbr fixe no years
  try(Mrate.QuantiA.TF.noY.Hnbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) +  treeNbrJ +
                                                           (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.2,
                                                         zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QuantiA.TF.noY.Hnbinom1),silent=T) 
  
  
  
  ## Offset tree 
  ### Recrut rate Marginality Quali Fixed Recrut Treenbr offset no years
  try(Mrate.QualiF.TO.noY.Hnbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) +  offset(log(treeNbrJ)) +
                                                          Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.recruitment.plot.rate.2,
                                                        zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QualiF.TO.noY.Hnbinom1),silent=T) 
  
  ### Recrut rate Marginality Quanti Fixed Recrut Treenbr offset no years
  try(Mrate.QuantiF.TO.noY.Hnbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) +  offset(log(treeNbrJ)) +
                                                           M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.recruitment.plot.rate.2,
                                                         zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QuantiF.TO.noY.Hnbinom1),silent=T) 
  
  ### Recrut rate Marginality Quali Aleatoire Recrut Treenbr offset no years 
  try(Mrate.QualiA.TO.noY.Hnbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) +  offset(log(treeNbrJ)) +
                                                          (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.2,
                                                        zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QualiA.TO.noY.Hnbinom1),silent=T) 
  
  ### Recrut rate Marginality Quanti Aleatoire Recrut Treenbr offset no years
  try(Mrate.QuantiA.TO.noY.Hnbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) +  offset(log(treeNbrJ)) +
                                                           (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.2,
                                                         zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QuantiA.TO.noY.Hnbinom1),silent=T) 
  
  
  ## No tree parameter 
  ### Recrut rate Marginality Quali Fixed Recrut no treenbr no years
  try(Mrate.QualiF.noT.noY.Hnbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + 
                                                           Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.recruitment.plot.rate.2,
                                                         zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QualiF.noT.noY.Hnbinom1),silent=T) 
  
  ### Recrut rate Marginality Quanti Fixed Recrut no treenbr no years
  try(Mrate.QuantiF.noT.noY.Hnbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + 
                                                            M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.recruitment.plot.rate.2,
                                                          zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QuantiF.noT.noY.Hnbinom1),silent=T) 
  
  ### Recrut rate Marginality Quali AleanoTire Recrut no treenbr no years 
  try(Mrate.QualiA.noT.noY.Hnbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + 
                                                           (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.2,
                                                         zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QualiA.noT.noY.Hnbinom1),silent=T) 
  
  ### Recrut rate Marginality Quanti AleanoTire Recrut no treenbr no years
  try(Mrate.QuantiA.noT.noY.Hnbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + 
                                                            (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.2,
                                                          zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QuantiA.noT.noY.Hnbinom1),silent=T) 
  
  ##############
  ## Distrib4 ##
  ##############
  
  ### Recrut rate Marginality Quali Fixed Recrut Treenbr fixe yearsbetweenstudy fixed
  try(Mrate.QualiF.TF.YF.Hnbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys + treeNbrJ +
                                                         Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.recruitment.plot.rate.2,
                                                       zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QualiF.TF.YF.Hnbinom2),silent=T) 
  
  ### Recrut rate Marginality Quanti Fixed Recrut Treenbr fixe yearsbetweenstudy fixed
  try(Mrate.QuantiF.TF.YF.Hnbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys + treeNbrJ +
                                                          M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.recruitment.plot.rate.2,
                                                        zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QuantiF.TF.YF.Hnbinom2),silent=T) 
  
  ### Recrut rate Marginality Quali Aleatoire Recrut Treenbr fixe yearsbetweenstudy fixed 
  try(Mrate.QualiA.TF.YF.Hnbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys + treeNbrJ +
                                                         (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.2,
                                                       zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QualiA.TF.YF.Hnbinom2),silent=T) 
  
  ### Recrut rate Marginality Quanti Aleatoire Recrut Treenbr fixe yearsbetweenstudy fixed
  try(Mrate.QuantiA.TF.YF.Hnbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys + treeNbrJ +
                                                          (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.2,
                                                        zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QuantiA.TF.YF.Hnbinom2),silent=T) 
  
  
  
  ## Offset tree 
  ### Recrut rate Marginality Quali Fixed Recrut Treenbr offset yearsbetweenstudy fixed
  try(Mrate.QualiF.TO.YF.Hnbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys + offset(log(treeNbrJ)) +
                                                         Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.recruitment.plot.rate.2,
                                                       zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QualiF.TO.YF.Hnbinom2),silent=T) 
  
  ### Recrut rate Marginality Quanti Fixed Recrut Treenbr offset yearsbetweenstudy fixed
  try(Mrate.QuantiF.TO.YF.Hnbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys + offset(log(treeNbrJ)) +
                                                          M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.recruitment.plot.rate.2,
                                                        zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QuantiF.TO.YF.Hnbinom2),silent=T) 
  
  ### Recrut rate Marginality Quali Aleatoire Recrut Treenbr offset yearsbetweenstudy fixed 
  try(Mrate.QualiA.TO.YF.Hnbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys + offset(log(treeNbrJ)) +
                                                         (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.2,
                                                       zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QualiA.TO.YF.Hnbinom2),silent=T) 
  
  ### Recrut rate Marginality Quanti Aleatoire Recrut Treenbr offset yearsbetweenstudy fixed
  try(Mrate.QuantiA.TO.YF.Hnbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys + offset(log(treeNbrJ)) +
                                                          (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.2,
                                                        zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QuantiA.TO.YF.Hnbinom2),silent=T) 
  
  
  ## No tree parameter 
  ### Recrut rate Marginality Quali Fixed Recrut no treenbr yearsbetweenstudy fixed
  try(Mrate.QualiF.noT.YF.Hnbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys  +
                                                          Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.recruitment.plot.rate.2,
                                                        zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QualiF.noT.YF.Hnbinom2),silent=T) 
  
  ### Recrut rate Marginality Quanti Fixed Recrut no treenbr yearsbetweenstudy fixed
  try(Mrate.QuantiF.noT.YF.Hnbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys  +
                                                           M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.recruitment.plot.rate.2,
                                                         zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QuantiF.noT.YF.Hnbinom2),silent=T) 
  
  ### Recrut rate Marginality Quali AleanoTire Recrut no treenbr yearsbetweenstudy fixed 
  try(Mrate.QualiA.noT.YF.Hnbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys  +
                                                          (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.2,
                                                        zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QualiA.noT.YF.Hnbinom2),silent=T) 
  
  ### Recrut rate Marginality Quanti AleanoTire Recrut no treenbr yearsbetweenstudy fixed
  try(Mrate.QuantiA.noT.YF.Hnbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys  +
                                                           (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.2,
                                                         zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QuantiA.noT.YF.Hnbinom2),silent=T) 
  
  
  
  
  #### Years between study as an offset 
  ### Recrut rate Marginality Quali Fixed Recrut Treenbr fixe yearsbetween study offset
  try(Mrate.QualiF.TF.YO.Hnbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys)) + treeNbrJ +
                                                         Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.recruitment.plot.rate.2,
                                                       zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QualiF.TF.YO.Hnbinom2),silent=T) 
  
  ### Recrut rate Marginality Quanti Fixed Recrut Treenbr fixe yearsbetween study offset
  try(Mrate.QuantiF.TF.YO.Hnbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys)) + treeNbrJ +
                                                          M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.recruitment.plot.rate.2,
                                                        zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QuantiF.TF.YO.Hnbinom2),silent=T) 
  
  ### Recrut rate Marginality Quali Aleatoire Recrut Treenbr fixe yearsbetween study offset 
  try(Mrate.QualiA.TF.YO.Hnbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys)) + treeNbrJ +
                                                         (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.2,
                                                       zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QualiA.TF.YO.Hnbinom2),silent=T) 
  
  ### Recrut rate Marginality Quanti Aleatoire Recrut Treenbr fixe yearsbetween study offset
  try(Mrate.QuantiA.TF.YO.Hnbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys)) + treeNbrJ +
                                                          (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.2,
                                                        zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QuantiA.TF.YO.Hnbinom2),silent=T) 
  
  
  
  ## Offset tree 
  ### Recrut rate Marginality Quali Fixed Recrut Treenbr offset yearsbetween study offset
  try(Mrate.QualiF.TO.YO.Hnbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys)) + offset(log(treeNbrJ)) +
                                                         Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.recruitment.plot.rate.2,
                                                       zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QualiF.TO.YO.Hnbinom2),silent=T) 
  
  ### Recrut rate Marginality Quanti Fixed Recrut Treenbr offset yearsbetween study offset
  try(Mrate.QuantiF.TO.YO.Hnbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys)) + offset(log(treeNbrJ)) +
                                                          M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.recruitment.plot.rate.2,
                                                        zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QuantiF.TO.YO.Hnbinom2),silent=T) 
  
  ### Recrut rate Marginality Quali Aleatoire Recrut Treenbr offset yearsbetween study offset 
  try(Mrate.QualiA.TO.YO.Hnbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys)) + offset(log(treeNbrJ)) +
                                                         (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.2,
                                                       zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QualiA.TO.YO.Hnbinom2),silent=T) 
  
  ### Recrut rate Marginality Quanti Aleatoire Recrut Treenbr offset yearsbetween study offset
  try(Mrate.QuantiA.TO.YO.Hnbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys)) + offset(log(treeNbrJ)) +
                                                          (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.2,
                                                        zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QuantiA.TO.YO.Hnbinom2),silent=T) 
  
  
  ## No tree parameter 
  ### Recrut rate Marginality Quali Fixed Recrut no treenbr yearsbetween study offset
  try(Mrate.QualiF.noT.YO.Hnbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys))  +
                                                          Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.recruitment.plot.rate.2,
                                                        zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QualiF.noT.YO.Hnbinom2),silent=T) 
  
  ### Recrut rate Marginality Quanti Fixed Recrut no treenbr yearsbetween study offset
  try(Mrate.QuantiF.noT.YO.Hnbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys))  +
                                                           M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.recruitment.plot.rate.2,
                                                         zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QuantiF.noT.YO.Hnbinom2),silent=T) 
  
  ### Recrut rate Marginality Quali AleanoTire Recrut no treenbr yearsbetween study offset 
  try(Mrate.QualiA.noT.YO.Hnbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys))  +
                                                          (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.2,
                                                        zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QualiA.noT.YO.Hnbinom2),silent=T) 
  
  ### Recrut rate Marginality Quanti AleanoTire Recrut no treenbr yearsbetween study offset
  try(Mrate.QuantiA.noT.YO.Hnbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys))  +
                                                           (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.2,
                                                         zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QuantiA.noT.YO.Hnbinom2),silent=T) 
  
  ## Yearsbetweenstudy not as a predictor
  ### Recrut rate Marginality Quali Fixed Recrut Treenbr fixe no years
  try(Mrate.QualiF.TF.noY.Hnbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + treeNbrJ +
                                                          Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.recruitment.plot.rate.2,
                                                        zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QualiF.TF.noY.Hnbinom2),silent=T) 
  
  ### Recrut rate Marginality Quanti Fixed Recrut Treenbr fixe no years
  try(Mrate.QuantiF.TF.noY.Hnbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) +  treeNbrJ +
                                                           M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.recruitment.plot.rate.2,
                                                         zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QuantiF.TF.noY.Hnbinom2),silent=T) 
  
  ### Recrut rate Marginality Quali Aleatoire Recrut Treenbr fixe no years 
  try(Mrate.QualiA.TF.noY.Hnbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) +  treeNbrJ +
                                                          (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.2,
                                                        zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QualiA.TF.noY.Hnbinom2),silent=T) 
  
  ### Recrut rate Marginality Quanti Aleatoire Recrut Treenbr fixe no years
  try(Mrate.QuantiA.TF.noY.Hnbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) +  treeNbrJ +
                                                           (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.2,
                                                         zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QuantiA.TF.noY.Hnbinom2),silent=T) 
  
  
  
  ## Offset tree 
  ### Recrut rate Marginality Quali Fixed Recrut Treenbr offset no years
  try(Mrate.QualiF.TO.noY.Hnbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) +  offset(log(treeNbrJ)) +
                                                          Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.recruitment.plot.rate.2,
                                                        zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QualiF.TO.noY.Hnbinom2),silent=T) 
  
  ### Recrut rate Marginality Quanti Fixed Recrut Treenbr offset no years
  try(Mrate.QuantiF.TO.noY.Hnbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) +  offset(log(treeNbrJ)) +
                                                           M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.recruitment.plot.rate.2,
                                                         zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QuantiF.TO.noY.Hnbinom2),silent=T) 
  
  ### Recrut rate Marginality Quali Aleatoire Recrut Treenbr offset no years 
  try(Mrate.QualiA.TO.noY.Hnbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) +  offset(log(treeNbrJ)) +
                                                          (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.2,
                                                        zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QualiA.TO.noY.Hnbinom2),silent=T) 
  
  ### Recrut rate Marginality Quanti Aleatoire Recrut Treenbr offset no years
  try(Mrate.QuantiA.TO.noY.Hnbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) +  offset(log(treeNbrJ)) +
                                                           (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.2,
                                                         zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QuantiA.TO.noY.Hnbinom2),silent=T) 
  
  
  ## No tree parameter 
  ### Recrut rate Marginality Quali Fixed Recrut no treenbr no years
  try(Mrate.QualiF.noT.noY.Hnbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + 
                                                           Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.recruitment.plot.rate.2,
                                                         zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QualiF.noT.noY.Hnbinom2),silent=T) 
  
  ### Recrut rate Marginality Quanti Fixed Recrut no treenbr no years
  try(Mrate.QuantiF.noT.noY.Hnbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + 
                                                            M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.recruitment.plot.rate.2,
                                                          zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QuantiF.noT.noY.Hnbinom2),silent=T) 
  
  ### Recrut rate Marginality Quali AleanoTire Recrut no treenbr no years 
  try(Mrate.QualiA.noT.noY.Hnbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + 
                                                           (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.2,
                                                         zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QualiA.noT.noY.Hnbinom2),silent=T) 
  
  ### Recrut rate Marginality Quanti AleanoTire Recrut no treenbr no years
  try(Mrate.QuantiA.noT.noY.Hnbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + 
                                                            (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.2,
                                                          zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.QuantiA.noT.noY.Hnbinom2),silent=T) 
  
  
  
  #########################
  ##### Mortality tot yr #####
  ########################
  
  ### Recrut rate Marginality Quali Fixed Recrut Treenbr fixe yearsbetweenstudy fixed
  try(Mrate.yr.QualiF.TF.YF.nbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys + treeNbrJ +
                                                        Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.recruitment.plot.rate.yr.2,
                                                      zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QualiF.TF.YF.nbinom1),silent=T) 
  
  ### Recrut rate Marginality Quanti Fixed Recrut Treenbr fixe yearsbetweenstudy fixed
  try(Mrate.yr.QuantiF.TF.YF.nbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys + treeNbrJ +
                                                         M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.recruitment.plot.rate.yr.2,
                                                       zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QuantiF.TF.YF.nbinom1),silent=T) 
  
  ### Recrut rate Marginality Quali Aleatoire Recrut Treenbr fixe yearsbetweenstudy fixed 
  try(Mrate.yr.QualiA.TF.YF.nbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys + treeNbrJ +
                                                        (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.yr.2,
                                                      zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QualiA.TF.YF.nbinom1),silent=T) 
  
  ### Recrut rate Marginality Quanti Aleatoire Recrut Treenbr fixe yearsbetweenstudy fixed
  try(Mrate.yr.QuantiA.TF.YF.nbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys + treeNbrJ +
                                                         (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.yr.2,
                                                       zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QuantiA.TF.YF.nbinom1),silent=T) 
  
  
  
  ## Offset tree 
  ### Recrut rate Marginality Quali Fixed Recrut Treenbr offset yearsbetweenstudy fixed
  try(Mrate.yr.QualiF.TO.YF.nbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys + offset(log(treeNbrJ)) +
                                                        Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.recruitment.plot.rate.yr.2,
                                                      zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QualiF.TO.YF.nbinom1),silent=T) 
  
  ### Recrut rate Marginality Quanti Fixed Recrut Treenbr offset yearsbetweenstudy fixed
  try(Mrate.yr.QuantiF.TO.YF.nbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys + offset(log(treeNbrJ)) +
                                                         M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.recruitment.plot.rate.yr.2,
                                                       zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QuantiF.TO.YF.nbinom1),silent=T) 
  
  ### Recrut rate Marginality Quali Aleatoire Recrut Treenbr offset yearsbetweenstudy fixed 
  try(Mrate.yr.QualiA.TO.YF.nbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys + offset(log(treeNbrJ)) +
                                                        (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.yr.2,
                                                      zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QualiA.TO.YF.nbinom1),silent=T) 
  
  ### Recrut rate Marginality Quanti Aleatoire Recrut Treenbr offset yearsbetweenstudy fixed
  try(Mrate.yr.QuantiA.TO.YF.nbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys + offset(log(treeNbrJ)) +
                                                         (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.yr.2,
                                                       zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QuantiA.TO.YF.nbinom1),silent=T) 
  
  
  ## No tree parameter 
  ### Recrut rate Marginality Quali Fixed Recrut no treenbr yearsbetweenstudy fixed
  try(Mrate.yr.QualiF.noT.YF.nbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys  +
                                                         Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.recruitment.plot.rate.yr.2,
                                                       zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QualiF.noT.YF.nbinom1),silent=T) 
  
  ### Recrut rate Marginality Quanti Fixed Recrut no treenbr yearsbetweenstudy fixed
  try(Mrate.yr.QuantiF.noT.YF.nbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys  +
                                                          M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.recruitment.plot.rate.yr.2,
                                                        zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QuantiF.noT.YF.nbinom1),silent=T) 
  
  ### Recrut rate Marginality Quali AleanoTire Recrut no treenbr yearsbetweenstudy fixed 
  try(Mrate.yr.QualiA.noT.YF.nbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys  +
                                                         (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.yr.2,
                                                       zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QualiA.noT.YF.nbinom1),silent=T) 
  
  ### Recrut rate Marginality Quanti AleanoTire Recrut no treenbr yearsbetweenstudy fixed
  try(Mrate.yr.QuantiA.noT.YF.nbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys  +
                                                          (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.yr.2,
                                                        zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QuantiA.noT.YF.nbinom1),silent=T) 
  
  
  
  
  #### Years between study as an offset 
  ### Recrut rate Marginality Quali Fixed Recrut Treenbr fixe yearsbetween study offset
  try(Mrate.yr.QualiF.TF.YO.nbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys)) + treeNbrJ +
                                                        Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.recruitment.plot.rate.yr.2,
                                                      zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QualiF.TF.YO.nbinom1),silent=T) 
  
  ### Recrut rate Marginality Quanti Fixed Recrut Treenbr fixe yearsbetween study offset
  try(Mrate.yr.QuantiF.TF.YO.nbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys)) + treeNbrJ +
                                                         M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.recruitment.plot.rate.yr.2,
                                                       zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QuantiF.TF.YO.nbinom1),silent=T) 
  
  ### Recrut rate Marginality Quali Aleatoire Recrut Treenbr fixe yearsbetween study offset 
  try(Mrate.yr.QualiA.TF.YO.nbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys)) + treeNbrJ +
                                                        (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.yr.2,
                                                      zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QualiA.TF.YO.nbinom1),silent=T) 
  
  ### Recrut rate Marginality Quanti Aleatoire Recrut Treenbr fixe yearsbetween study offset
  try(Mrate.yr.QuantiA.TF.YO.nbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys)) + treeNbrJ +
                                                         (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.yr.2,
                                                       zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QuantiA.TF.YO.nbinom1),silent=T) 
  
  
  
  ## Offset tree 
  ### Recrut rate Marginality Quali Fixed Recrut Treenbr offset yearsbetween study offset
  try(Mrate.yr.QualiF.TO.YO.nbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys)) + offset(log(treeNbrJ)) +
                                                        Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.recruitment.plot.rate.yr.2,
                                                      zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QualiF.TO.YO.nbinom1),silent=T) 
  
  ### Recrut rate Marginality Quanti Fixed Recrut Treenbr offset yearsbetween study offset
  try(Mrate.yr.QuantiF.TO.YO.nbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys)) + offset(log(treeNbrJ)) +
                                                         M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.recruitment.plot.rate.yr.2,
                                                       zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QuantiF.TO.YO.nbinom1),silent=T) 
  
  ### Recrut rate Marginality Quali Aleatoire Recrut Treenbr offset yearsbetween study offset 
  try(Mrate.yr.QualiA.TO.YO.nbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys)) + offset(log(treeNbrJ)) +
                                                        (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.yr.2,
                                                      zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QualiA.TO.YO.nbinom1),silent=T) 
  
  ### Recrut rate Marginality Quanti Aleatoire Recrut Treenbr offset yearsbetween study offset
  try(Mrate.yr.QuantiA.TO.YO.nbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys)) + offset(log(treeNbrJ)) +
                                                         (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.yr.2,
                                                       zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QuantiA.TO.YO.nbinom1),silent=T) 
  
  
  ## No tree parameter 
  ### Recrut rate Marginality Quali Fixed Recrut no treenbr yearsbetween study offset
  try(Mrate.yr.QualiF.noT.YO.nbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys))  +
                                                         Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.recruitment.plot.rate.yr.2,
                                                       zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QualiF.noT.YO.nbinom1),silent=T) 
  
  ### Recrut rate Marginality Quanti Fixed Recrut no treenbr yearsbetween study offset
  try(Mrate.yr.QuantiF.noT.YO.nbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys))  +
                                                          M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.recruitment.plot.rate.yr.2,
                                                        zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QuantiF.noT.YO.nbinom1),silent=T) 
  
  ### Recrut rate Marginality Quali AleanoTire Recrut no treenbr yearsbetween study offset 
  try(Mrate.yr.QualiA.noT.YO.nbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys))  +
                                                         (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.yr.2,
                                                       zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QualiA.noT.YO.nbinom1),silent=T) 
  
  ### Recrut rate Marginality Quanti AleanoTire Recrut no treenbr yearsbetween study offset
  try(Mrate.yr.QuantiA.noT.YO.nbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys))  +
                                                          (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.yr.2,
                                                        zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QuantiA.noT.YO.nbinom1),silent=T) 
  
  ## Yearsbetweenstudy not as a predictor
  ### Recrut rate Marginality Quali Fixed Recrut Treenbr fixe no years
  try(Mrate.yr.QualiF.TF.noY.nbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) + treeNbrJ +
                                                         Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.recruitment.plot.rate.yr.2,
                                                       zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QualiF.TF.noY.nbinom1),silent=T) 
  
  ### Recrut rate Marginality Quanti Fixed Recrut Treenbr fixe no years
  try(Mrate.yr.QuantiF.TF.noY.nbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) +  treeNbrJ +
                                                          M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.recruitment.plot.rate.yr.2,
                                                        zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QuantiF.TF.noY.nbinom1),silent=T) 
  
  ### Recrut rate Marginality Quali Aleatoire Recrut Treenbr fixe no years 
  try(Mrate.yr.QualiA.TF.noY.nbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) +  treeNbrJ +
                                                         (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.yr.2,
                                                       zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QualiA.TF.noY.nbinom1),silent=T) 
  
  ### Recrut rate Marginality Quanti Aleatoire Recrut Treenbr fixe no years
  try(Mrate.yr.QuantiA.TF.noY.nbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) +  treeNbrJ +
                                                          (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.yr.2,
                                                        zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QuantiA.TF.noY.nbinom1),silent=T) 
  
  
  
  ## Offset tree 
  ### Recrut rate Marginality Quali Fixed Recrut Treenbr offset no years
  try(Mrate.yr.QualiF.TO.noY.nbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) +  offset(log(treeNbrJ)) +
                                                         Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.recruitment.plot.rate.yr.2,
                                                       zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QualiF.TO.noY.nbinom1),silent=T) 
  
  ### Recrut rate Marginality Quanti Fixed Recrut Treenbr offset no years
  try(Mrate.yr.QuantiF.TO.noY.nbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) +  offset(log(treeNbrJ)) +
                                                          M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.recruitment.plot.rate.yr.2,
                                                        zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QuantiF.TO.noY.nbinom1),silent=T) 
  
  ### Recrut rate Marginality Quali Aleatoire Recrut Treenbr offset no years 
  try(Mrate.yr.QualiA.TO.noY.nbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) +  offset(log(treeNbrJ)) +
                                                         (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.yr.2,
                                                       zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QualiA.TO.noY.nbinom1),silent=T) 
  
  ### Recrut rate Marginality Quanti Aleatoire Recrut Treenbr offset no years
  try(Mrate.yr.QuantiA.TO.noY.nbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) +  offset(log(treeNbrJ)) +
                                                          (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.yr.2,
                                                        zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QuantiA.TO.noY.nbinom1),silent=T) 
  
  
  ## No tree parameter 
  ### Recrut rate Marginality Quali Fixed Recrut no treenbr no years
  try(Mrate.yr.QualiF.noT.noY.nbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) + 
                                                          Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.recruitment.plot.rate.yr.2,
                                                        zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QualiF.noT.noY.nbinom1),silent=T) 
  
  ### Recrut rate Marginality Quanti Fixed Recrut no treenbr no years
  try(Mrate.yr.QuantiF.noT.noY.nbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) + 
                                                           M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.recruitment.plot.rate.yr.2,
                                                         zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QuantiF.noT.noY.nbinom1),silent=T) 
  
  ### Recrut rate Marginality Quali AleanoTire Recrut no treenbr no years 
  try(Mrate.yr.QualiA.noT.noY.nbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) + 
                                                          (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.yr.2,
                                                        zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QualiA.noT.noY.nbinom1),silent=T) 
  
  ### Recrut rate Marginality Quanti AleanoTire Recrut no treenbr no years
  try(Mrate.yr.QuantiA.noT.noY.nbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) + 
                                                           (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.yr.2,
                                                         zi=~.,family=nbinom1, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QuantiA.noT.noY.nbinom1),silent=T) 
  
  
  ##################
  # Idem distrib 2 ##
  ##################
  
  ### Recrut rate Marginality Quali Fixed Recrut Treenbr fixe yearsbetweenstudy fixed
  try(Mrate.yr.QualiF.TF.YF.nbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys + treeNbrJ +
                                                        Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.recruitment.plot.rate.yr.2,
                                                      zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QualiF.TF.YF.nbinom2),silent=T) 
  
  ### Recrut rate Marginality Quanti Fixed Recrut Treenbr fixe yearsbetweenstudy fixed
  try(Mrate.yr.QuantiF.TF.YF.nbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys + treeNbrJ +
                                                         M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.recruitment.plot.rate.yr.2,
                                                       zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QuantiF.TF.YF.nbinom2),silent=T) 
  
  ### Recrut rate Marginality Quali Aleatoire Recrut Treenbr fixe yearsbetweenstudy fixed 
  try(Mrate.yr.QualiA.TF.YF.nbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys + treeNbrJ +
                                                        (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.yr.2,
                                                      zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QualiA.TF.YF.nbinom2),silent=T) 
  
  ### Recrut rate Marginality Quanti Aleatoire Recrut Treenbr fixe yearsbetweenstudy fixed
  try(Mrate.yr.QuantiA.TF.YF.nbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys + treeNbrJ +
                                                         (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.yr.2,
                                                       zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QuantiA.TF.YF.nbinom2),silent=T) 
  
  
  
  ## Offset tree 
  ### Recrut rate Marginality Quali Fixed Recrut Treenbr offset yearsbetweenstudy fixed
  try(Mrate.yr.QualiF.TO.YF.nbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys + offset(log(treeNbrJ)) +
                                                        Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.recruitment.plot.rate.yr.2,
                                                      zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QualiF.TO.YF.nbinom2),silent=T) 
  
  ### Recrut rate Marginality Quanti Fixed Recrut Treenbr offset yearsbetweenstudy fixed
  try(Mrate.yr.QuantiF.TO.YF.nbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys + offset(log(treeNbrJ)) +
                                                         M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.recruitment.plot.rate.yr.2,
                                                       zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QuantiF.TO.YF.nbinom2),silent=T) 
  
  ### Recrut rate Marginality Quali Aleatoire Recrut Treenbr offset yearsbetweenstudy fixed 
  try(Mrate.yr.QualiA.TO.YF.nbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys + offset(log(treeNbrJ)) +
                                                        (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.yr.2,
                                                      zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QualiA.TO.YF.nbinom2),silent=T) 
  
  ### Recrut rate Marginality Quanti Aleatoire Recrut Treenbr offset yearsbetweenstudy fixed
  try(Mrate.yr.QuantiA.TO.YF.nbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys + offset(log(treeNbrJ)) +
                                                         (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.yr.2,
                                                       zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QuantiA.TO.YF.nbinom2),silent=T) 
  
  
  ## No tree parameter 
  ### Recrut rate Marginality Quali Fixed Recrut no treenbr yearsbetweenstudy fixed
  try(Mrate.yr.QualiF.noT.YF.nbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys  +
                                                         Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.recruitment.plot.rate.yr.2,
                                                       zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QualiF.noT.YF.nbinom2),silent=T) 
  
  ### Recrut rate Marginality Quanti Fixed Recrut no treenbr yearsbetweenstudy fixed
  try(Mrate.yr.QuantiF.noT.YF.nbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys  +
                                                          M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.recruitment.plot.rate.yr.2,
                                                        zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QuantiF.noT.YF.nbinom2),silent=T) 
  
  ### Recrut rate Marginality Quali AleanoTire Recrut no treenbr yearsbetweenstudy fixed 
  try(Mrate.yr.QualiA.noT.YF.nbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys  +
                                                         (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.yr.2,
                                                       zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QualiA.noT.YF.nbinom2),silent=T) 
  
  ### Recrut rate Marginality Quanti AleanoTire Recrut no treenbr yearsbetweenstudy fixed
  try(Mrate.yr.QuantiA.noT.YF.nbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys  +
                                                          (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.yr.2,
                                                        zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QuantiA.noT.YF.nbinom2),silent=T) 
  
  
  
  
  #### Years between study as an offset 
  ### Recrut rate Marginality Quali Fixed Recrut Treenbr fixe yearsbetween study offset
  try(Mrate.yr.QualiF.TF.YO.nbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys)) + treeNbrJ +
                                                        Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.recruitment.plot.rate.yr.2,
                                                      zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QualiF.TF.YO.nbinom2),silent=T) 
  
  ### Recrut rate Marginality Quanti Fixed Recrut Treenbr fixe yearsbetween study offset
  try(Mrate.yr.QuantiF.TF.YO.nbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys)) + treeNbrJ +
                                                         M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.recruitment.plot.rate.yr.2,
                                                       zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QuantiF.TF.YO.nbinom2),silent=T) 
  
  ### Recrut rate Marginality Quali Aleatoire Recrut Treenbr fixe yearsbetween study offset 
  try(Mrate.yr.QualiA.TF.YO.nbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys)) + treeNbrJ +
                                                        (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.yr.2,
                                                      zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QualiA.TF.YO.nbinom2),silent=T) 
  
  ### Recrut rate Marginality Quanti Aleatoire Recrut Treenbr fixe yearsbetween study offset
  try(Mrate.yr.QuantiA.TF.YO.nbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys)) + treeNbrJ +
                                                         (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.yr.2,
                                                       zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QuantiA.TF.YO.nbinom2),silent=T) 
  
  
  
  ## Offset tree 
  ### Recrut rate Marginality Quali Fixed Recrut Treenbr offset yearsbetween study offset
  try(Mrate.yr.QualiF.TO.YO.nbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys)) + offset(log(treeNbrJ)) +
                                                        Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.recruitment.plot.rate.yr.2,
                                                      zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QualiF.TO.YO.nbinom2),silent=T) 
  
  ### Recrut rate Marginality Quanti Fixed Recrut Treenbr offset yearsbetween study offset
  try(Mrate.yr.QuantiF.TO.YO.nbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys)) + offset(log(treeNbrJ)) +
                                                         M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.recruitment.plot.rate.yr.2,
                                                       zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QuantiF.TO.YO.nbinom2),silent=T) 
  
  ### Recrut rate Marginality Quali Aleatoire Recrut Treenbr offset yearsbetween study offset 
  try(Mrate.yr.QualiA.TO.YO.nbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys)) + offset(log(treeNbrJ)) +
                                                        (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.yr.2,
                                                      zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QualiA.TO.YO.nbinom2),silent=T) 
  
  ### Recrut rate Marginality Quanti Aleatoire Recrut Treenbr offset yearsbetween study offset
  try(Mrate.yr.QuantiA.TO.YO.nbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys)) + offset(log(treeNbrJ)) +
                                                         (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.yr.2,
                                                       zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QuantiA.TO.YO.nbinom2),silent=T) 
  
  
  ## No tree parameter 
  ### Recrut rate Marginality Quali Fixed Recrut no treenbr yearsbetween study offset
  try(Mrate.yr.QualiF.noT.YO.nbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys))  +
                                                         Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.recruitment.plot.rate.yr.2,
                                                       zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QualiF.noT.YO.nbinom2),silent=T) 
  
  ### Recrut rate Marginality Quanti Fixed Recrut no treenbr yearsbetween study offset
  try(Mrate.yr.QuantiF.noT.YO.nbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys))  +
                                                          M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.recruitment.plot.rate.yr.2,
                                                        zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QuantiF.noT.YO.nbinom2),silent=T) 
  
  ### Recrut rate Marginality Quali AleanoTire Recrut no treenbr yearsbetween study offset 
  try(Mrate.yr.QualiA.noT.YO.nbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys))  +
                                                         (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.yr.2,
                                                       zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QualiA.noT.YO.nbinom2),silent=T) 
  
  ### Recrut rate Marginality Quanti AleanoTire Recrut no treenbr yearsbetween study offset
  try(Mrate.yr.QuantiA.noT.YO.nbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys))  +
                                                          (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.yr.2,
                                                        zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QuantiA.noT.YO.nbinom2),silent=T) 
  
  ## Yearsbetweenstudy not as a predictor
  ### Recrut rate Marginality Quali Fixed Recrut Treenbr fixe no years
  try(Mrate.yr.QualiF.TF.noY.nbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) + treeNbrJ +
                                                         Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.recruitment.plot.rate.yr.2,
                                                       zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QualiF.TF.noY.nbinom2),silent=T) 
  
  ### Recrut rate Marginality Quanti Fixed Recrut Treenbr fixe no years
  try(Mrate.yr.QuantiF.TF.noY.nbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) +  treeNbrJ +
                                                          M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.recruitment.plot.rate.yr.2,
                                                        zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QuantiF.TF.noY.nbinom2),silent=T) 
  
  ### Recrut rate Marginality Quali Aleatoire Recrut Treenbr fixe no years 
  try(Mrate.yr.QualiA.TF.noY.nbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) +  treeNbrJ +
                                                         (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.yr.2,
                                                       zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QualiA.TF.noY.nbinom2),silent=T) 
  
  ### Recrut rate Marginality Quanti Aleatoire Recrut Treenbr fixe no years
  try(Mrate.yr.QuantiA.TF.noY.nbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) +  treeNbrJ +
                                                          (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.yr.2,
                                                        zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QuantiA.TF.noY.nbinom2),silent=T) 
  
  
  
  ## Offset tree 
  ### Recrut rate Marginality Quali Fixed Recrut Treenbr offset no years
  try(Mrate.yr.QualiF.TO.noY.nbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) +  offset(log(treeNbrJ)) +
                                                         Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.recruitment.plot.rate.yr.2,
                                                       zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QualiF.TO.noY.nbinom2),silent=T) 
  
  ### Recrut rate Marginality Quanti Fixed Recrut Treenbr offset no years
  try(Mrate.yr.QuantiF.TO.noY.nbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) +  offset(log(treeNbrJ)) +
                                                          M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.recruitment.plot.rate.yr.2,
                                                        zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QuantiF.TO.noY.nbinom2),silent=T) 
  
  ### Recrut rate Marginality Quali Aleatoire Recrut Treenbr offset no years 
  try(Mrate.yr.QualiA.TO.noY.nbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) +  offset(log(treeNbrJ)) +
                                                         (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.yr.2,
                                                       zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QualiA.TO.noY.nbinom2),silent=T) 
  
  ### Recrut rate Marginality Quanti Aleatoire Recrut Treenbr offset no years
  try(Mrate.yr.QuantiA.TO.noY.nbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) +  offset(log(treeNbrJ)) +
                                                          (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.yr.2,
                                                        zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QuantiA.TO.noY.nbinom2),silent=T) 
  
  
  ## No tree parameter 
  ### Recrut rate Marginality Quali Fixed Recrut no treenbr no years
  try(Mrate.yr.QualiF.noT.noY.nbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) + 
                                                          Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.recruitment.plot.rate.yr.2,
                                                        zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QualiF.noT.noY.nbinom2),silent=T) 
  
  ### Recrut rate Marginality Quanti Fixed Recrut no treenbr no years
  try(Mrate.yr.QuantiF.noT.noY.nbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) + 
                                                           M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.recruitment.plot.rate.yr.2,
                                                         zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QuantiF.noT.noY.nbinom2),silent=T) 
  
  ### Recrut rate Marginality Quali AleanoTire Recrut no treenbr no years 
  try(Mrate.yr.QualiA.noT.noY.nbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) + 
                                                          (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.yr.2,
                                                        zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QualiA.noT.noY.nbinom2),silent=T) 
  
  ### Recrut rate Marginality Quanti AleanoTire Recrut no treenbr no years
  try(Mrate.yr.QuantiA.noT.noY.nbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) + 
                                                           (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.yr.2,
                                                         zi=~.,family=nbinom2, data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QuantiA.noT.noY.nbinom2),silent=T) 
  
  #####################
  ### Idem distrib 3 ##
  #####################
  
  ### Recrut rate Marginality Quali Fixed Recrut Treenbr fixe yearsbetweenstudy fixed
  try(Mrate.yr.QualiF.TF.YF.Hnbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys + treeNbrJ +
                                                         Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.recruitment.plot.rate.yr.2,
                                                       zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QualiF.TF.YF.Hnbinom1),silent=T) 
  
  ### Recrut rate Marginality Quanti Fixed Recrut Treenbr fixe yearsbetweenstudy fixed
  try(Mrate.yr.QuantiF.TF.YF.Hnbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys + treeNbrJ +
                                                          M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.recruitment.plot.rate.yr.2,
                                                        zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QuantiF.TF.YF.Hnbinom1),silent=T) 
  
  ### Recrut rate Marginality Quali Aleatoire Recrut Treenbr fixe yearsbetweenstudy fixed 
  try(Mrate.yr.QualiA.TF.YF.Hnbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys + treeNbrJ +
                                                         (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.yr.2,
                                                       zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QualiA.TF.YF.Hnbinom1),silent=T) 
  
  ### Recrut rate Marginality Quanti Aleatoire Recrut Treenbr fixe yearsbetweenstudy fixed
  try(Mrate.yr.QuantiA.TF.YF.Hnbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys + treeNbrJ +
                                                          (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.yr.2,
                                                        zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QuantiA.TF.YF.Hnbinom1),silent=T) 
  
  
  
  ## Offset tree 
  ### Recrut rate Marginality Quali Fixed Recrut Treenbr offset yearsbetweenstudy fixed
  try(Mrate.yr.QualiF.TO.YF.Hnbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys + offset(log(treeNbrJ)) +
                                                         Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.recruitment.plot.rate.yr.2,
                                                       zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QualiF.TO.YF.Hnbinom1),silent=T) 
  
  ### Recrut rate Marginality Quanti Fixed Recrut Treenbr offset yearsbetweenstudy fixed
  try(Mrate.yr.QuantiF.TO.YF.Hnbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys + offset(log(treeNbrJ)) +
                                                          M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.recruitment.plot.rate.yr.2,
                                                        zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QuantiF.TO.YF.Hnbinom1),silent=T) 
  
  ### Recrut rate Marginality Quali Aleatoire Recrut Treenbr offset yearsbetweenstudy fixed 
  try(Mrate.yr.QualiA.TO.YF.Hnbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys + offset(log(treeNbrJ)) +
                                                         (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.yr.2,
                                                       zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QualiA.TO.YF.Hnbinom1),silent=T) 
  
  ### Recrut rate Marginality Quanti Aleatoire Recrut Treenbr offset yearsbetweenstudy fixed
  try(Mrate.yr.QuantiA.TO.YF.Hnbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys + offset(log(treeNbrJ)) +
                                                          (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.yr.2,
                                                        zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QuantiA.TO.YF.Hnbinom1),silent=T) 
  
  
  ## No tree parameter 
  ### Recrut rate Marginality Quali Fixed Recrut no treenbr yearsbetweenstudy fixed
  try(Mrate.yr.QualiF.noT.YF.Hnbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys  +
                                                          Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.recruitment.plot.rate.yr.2,
                                                        zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QualiF.noT.YF.Hnbinom1),silent=T) 
  
  ### Recrut rate Marginality Quanti Fixed Recrut no treenbr yearsbetweenstudy fixed
  try(Mrate.yr.QuantiF.noT.YF.Hnbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys  +
                                                           M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.recruitment.plot.rate.yr.2,
                                                         zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QuantiF.noT.YF.Hnbinom1),silent=T) 
  
  ### Recrut rate Marginality Quali AleanoTire Recrut no treenbr yearsbetweenstudy fixed 
  try(Mrate.yr.QualiA.noT.YF.Hnbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys  +
                                                          (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.yr.2,
                                                        zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QualiA.noT.YF.Hnbinom1),silent=T) 
  
  ### Recrut rate Marginality Quanti AleanoTire Recrut no treenbr yearsbetweenstudy fixed
  try(Mrate.yr.QuantiA.noT.YF.Hnbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys  +
                                                           (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.yr.2,
                                                         zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QuantiA.noT.YF.Hnbinom1),silent=T) 
  
  
  
  
  #### Years between study as an offset 
  ### Recrut rate Marginality Quali Fixed Recrut Treenbr fixe yearsbetween study offset
  try(Mrate.yr.QualiF.TF.YO.Hnbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys)) + treeNbrJ +
                                                         Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.recruitment.plot.rate.yr.2,
                                                       zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QualiF.TF.YO.Hnbinom1),silent=T) 
  
  ### Recrut rate Marginality Quanti Fixed Recrut Treenbr fixe yearsbetween study offset
  try(Mrate.yr.QuantiF.TF.YO.Hnbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys)) + treeNbrJ +
                                                          M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.recruitment.plot.rate.yr.2,
                                                        zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QuantiF.TF.YO.Hnbinom1),silent=T) 
  
  ### Recrut rate Marginality Quali Aleatoire Recrut Treenbr fixe yearsbetween study offset 
  try(Mrate.yr.QualiA.TF.YO.Hnbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys)) + treeNbrJ +
                                                         (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.yr.2,
                                                       zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QualiA.TF.YO.Hnbinom1),silent=T) 
  
  ### Recrut rate Marginality Quanti Aleatoire Recrut Treenbr fixe yearsbetween study offset
  try(Mrate.yr.QuantiA.TF.YO.Hnbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys)) + treeNbrJ +
                                                          (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.yr.2,
                                                        zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QuantiA.TF.YO.Hnbinom1),silent=T) 
  
  
  
  ## Offset tree 
  ### Recrut rate Marginality Quali Fixed Recrut Treenbr offset yearsbetween study offset
  try(Mrate.yr.QualiF.TO.YO.Hnbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys)) + offset(log(treeNbrJ)) +
                                                         Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.recruitment.plot.rate.yr.2,
                                                       zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QualiF.TO.YO.Hnbinom1),silent=T) 
  
  ### Recrut rate Marginality Quanti Fixed Recrut Treenbr offset yearsbetween study offset
  try(Mrate.yr.QuantiF.TO.YO.Hnbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys)) + offset(log(treeNbrJ)) +
                                                          M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.recruitment.plot.rate.yr.2,
                                                        zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QuantiF.TO.YO.Hnbinom1),silent=T) 
  
  ### Recrut rate Marginality Quali Aleatoire Recrut Treenbr offset yearsbetween study offset 
  try(Mrate.yr.QualiA.TO.YO.Hnbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys)) + offset(log(treeNbrJ)) +
                                                         (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.yr.2,
                                                       zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QualiA.TO.YO.Hnbinom1),silent=T) 
  
  ### Recrut rate Marginality Quanti Aleatoire Recrut Treenbr offset yearsbetween study offset
  try(Mrate.yr.QuantiA.TO.YO.Hnbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys)) + offset(log(treeNbrJ)) +
                                                          (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.yr.2,
                                                        zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QuantiA.TO.YO.Hnbinom1),silent=T) 
  
  
  ## No tree parameter 
  ### Recrut rate Marginality Quali Fixed Recrut no treenbr yearsbetween study offset
  try(Mrate.yr.QualiF.noT.YO.Hnbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys))  +
                                                          Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.recruitment.plot.rate.yr.2,
                                                        zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QualiF.noT.YO.Hnbinom1),silent=T) 
  
  ### Recrut rate Marginality Quanti Fixed Recrut no treenbr yearsbetween study offset
  try(Mrate.yr.QuantiF.noT.YO.Hnbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys))  +
                                                           M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.recruitment.plot.rate.yr.2,
                                                         zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QuantiF.noT.YO.Hnbinom1),silent=T) 
  
  ### Recrut rate Marginality Quali AleanoTire Recrut no treenbr yearsbetween study offset 
  try(Mrate.yr.QualiA.noT.YO.Hnbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys))  +
                                                          (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.yr.2,
                                                        zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QualiA.noT.YO.Hnbinom1),silent=T) 
  
  ### Recrut rate Marginality Quanti AleanoTire Recrut no treenbr yearsbetween study offset
  try(Mrate.yr.QuantiA.noT.YO.Hnbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys))  +
                                                           (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.yr.2,
                                                         zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QuantiA.noT.YO.Hnbinom1),silent=T) 
  
  ## Yearsbetweenstudy not as a predictor
  ### Recrut rate Marginality Quali Fixed Recrut Treenbr fixe no years
  try(Mrate.yr.QualiF.TF.noY.Hnbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) + treeNbrJ +
                                                          Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.recruitment.plot.rate.yr.2,
                                                        zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QualiF.TF.noY.Hnbinom1),silent=T) 
  
  ### Recrut rate Marginality Quanti Fixed Recrut Treenbr fixe no years
  try(Mrate.yr.QuantiF.TF.noY.Hnbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) +  treeNbrJ +
                                                           M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.recruitment.plot.rate.yr.2,
                                                         zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QuantiF.TF.noY.Hnbinom1),silent=T) 
  
  ### Recrut rate Marginality Quali Aleatoire Recrut Treenbr fixe no years 
  try(Mrate.yr.QualiA.TF.noY.Hnbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) +  treeNbrJ +
                                                          (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.yr.2,
                                                        zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QualiA.TF.noY.Hnbinom1),silent=T) 
  
  ### Recrut rate Marginality Quanti Aleatoire Recrut Treenbr fixe no years
  try(Mrate.yr.QuantiA.TF.noY.Hnbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) +  treeNbrJ +
                                                           (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.yr.2,
                                                         zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QuantiA.TF.noY.Hnbinom1),silent=T) 
  
  
  
  ## Offset tree 
  ### Recrut rate Marginality Quali Fixed Recrut Treenbr offset no years
  try(Mrate.yr.QualiF.TO.noY.Hnbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) +  offset(log(treeNbrJ)) +
                                                          Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.recruitment.plot.rate.yr.2,
                                                        zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QualiF.TO.noY.Hnbinom1),silent=T) 
  
  ### Recrut rate Marginality Quanti Fixed Recrut Treenbr offset no years
  try(Mrate.yr.QuantiF.TO.noY.Hnbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) +  offset(log(treeNbrJ)) +
                                                           M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.recruitment.plot.rate.yr.2,
                                                         zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QuantiF.TO.noY.Hnbinom1),silent=T) 
  
  ### Recrut rate Marginality Quali Aleatoire Recrut Treenbr offset no years 
  try(Mrate.yr.QualiA.TO.noY.Hnbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) +  offset(log(treeNbrJ)) +
                                                          (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.yr.2,
                                                        zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QualiA.TO.noY.Hnbinom1),silent=T) 
  
  ### Recrut rate Marginality Quanti Aleatoire Recrut Treenbr offset no years
  try(Mrate.yr.QuantiA.TO.noY.Hnbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) +  offset(log(treeNbrJ)) +
                                                           (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.yr.2,
                                                         zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QuantiA.TO.noY.Hnbinom1),silent=T) 
  
  
  ## No tree parameter 
  ### Recrut rate Marginality Quali Fixed Recrut no treenbr no years
  try(Mrate.yr.QualiF.noT.noY.Hnbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) + 
                                                           Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.recruitment.plot.rate.yr.2,
                                                         zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QualiF.noT.noY.Hnbinom1),silent=T) 
  
  ### Recrut rate Marginality Quanti Fixed Recrut no treenbr no years
  try(Mrate.yr.QuantiF.noT.noY.Hnbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) + 
                                                            M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.recruitment.plot.rate.yr.2,
                                                          zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QuantiF.noT.noY.Hnbinom1),silent=T) 
  
  ### Recrut rate Marginality Quali AleanoTire Recrut no treenbr no years 
  try(Mrate.yr.QualiA.noT.noY.Hnbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) + 
                                                           (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.yr.2,
                                                         zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QualiA.noT.noY.Hnbinom1),silent=T) 
  
  ### Recrut rate Marginality Quanti AleanoTire Recrut no treenbr no years
  try(Mrate.yr.QuantiA.noT.noY.Hnbinom1 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) + 
                                                            (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.yr.2,
                                                          zi=~.,family="truncated_nbinom1", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QuantiA.noT.noY.Hnbinom1),silent=T) 
  
  ##############
  ## Distrib4 ##
  ##############
  
  ### Recrut rate Marginality Quali Fixed Recrut Treenbr fixe yearsbetweenstudy fixed
  try(Mrate.yr.QualiF.TF.YF.Hnbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys + treeNbrJ +
                                                         Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.recruitment.plot.rate.yr.2,
                                                       zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QualiF.TF.YF.Hnbinom2),silent=T) 
  
  ### Recrut rate Marginality Quanti Fixed Recrut Treenbr fixe yearsbetweenstudy fixed
  try(Mrate.yr.QuantiF.TF.YF.Hnbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys + treeNbrJ +
                                                          M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.recruitment.plot.rate.yr.2,
                                                        zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QuantiF.TF.YF.Hnbinom2),silent=T) 
  
  ### Recrut rate Marginality Quali Aleatoire Recrut Treenbr fixe yearsbetweenstudy fixed 
  try(Mrate.yr.QualiA.TF.YF.Hnbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys + treeNbrJ +
                                                         (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.yr.2,
                                                       zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QualiA.TF.YF.Hnbinom2),silent=T) 
  
  ### Recrut rate Marginality Quanti Aleatoire Recrut Treenbr fixe yearsbetweenstudy fixed
  try(Mrate.yr.QuantiA.TF.YF.Hnbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys + treeNbrJ +
                                                          (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.yr.2,
                                                        zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QuantiA.TF.YF.Hnbinom2),silent=T) 
  
  
  
  ## Offset tree 
  ### Recrut rate Marginality Quali Fixed Recrut Treenbr offset yearsbetweenstudy fixed
  try(Mrate.yr.QualiF.TO.YF.Hnbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys + offset(log(treeNbrJ)) +
                                                         Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.recruitment.plot.rate.yr.2,
                                                       zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QualiF.TO.YF.Hnbinom2),silent=T) 
  
  ### Recrut rate Marginality Quanti Fixed Recrut Treenbr offset yearsbetweenstudy fixed
  try(Mrate.yr.QuantiF.TO.YF.Hnbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys + offset(log(treeNbrJ)) +
                                                          M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.recruitment.plot.rate.yr.2,
                                                        zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QuantiF.TO.YF.Hnbinom2),silent=T) 
  
  ### Recrut rate Marginality Quali Aleatoire Recrut Treenbr offset yearsbetweenstudy fixed 
  try(Mrate.yr.QualiA.TO.YF.Hnbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys + offset(log(treeNbrJ)) +
                                                         (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.yr.2,
                                                       zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QualiA.TO.YF.Hnbinom2),silent=T) 
  
  ### Recrut rate Marginality Quanti Aleatoire Recrut Treenbr offset yearsbetweenstudy fixed
  try(Mrate.yr.QuantiA.TO.YF.Hnbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys + offset(log(treeNbrJ)) +
                                                          (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.yr.2,
                                                        zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QuantiA.TO.YF.Hnbinom2),silent=T) 
  
  
  ## No tree parameter 
  ### Recrut rate Marginality Quali Fixed Recrut no treenbr yearsbetweenstudy fixed
  try(Mrate.yr.QualiF.noT.YF.Hnbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys  +
                                                          Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.recruitment.plot.rate.yr.2,
                                                        zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QualiF.noT.YF.Hnbinom2),silent=T) 
  
  ### Recrut rate Marginality Quanti Fixed Recrut no treenbr yearsbetweenstudy fixed
  try(Mrate.yr.QuantiF.noT.YF.Hnbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys  +
                                                           M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.recruitment.plot.rate.yr.2,
                                                         zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QuantiF.noT.YF.Hnbinom2),silent=T) 
  
  ### Recrut rate Marginality Quali AleanoTire Recrut no treenbr yearsbetweenstudy fixed 
  try(Mrate.yr.QualiA.noT.YF.Hnbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys  +
                                                          (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.yr.2,
                                                        zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QualiA.noT.YF.Hnbinom2),silent=T) 
  
  ### Recrut rate Marginality Quanti AleanoTire Recrut no treenbr yearsbetweenstudy fixed
  try(Mrate.yr.QuantiA.noT.YF.Hnbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys  +
                                                           (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.yr.2,
                                                         zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QuantiA.noT.YF.Hnbinom2),silent=T) 
  
  
  
  
  #### Years between study as an offset 
  ### Recrut rate Marginality Quali Fixed Recrut Treenbr fixe yearsbetween study offset
  try(Mrate.yr.QualiF.TF.YO.Hnbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys)) + treeNbrJ +
                                                         Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.recruitment.plot.rate.yr.2,
                                                       zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QualiF.TF.YO.Hnbinom2),silent=T) 
  
  ### Recrut rate Marginality Quanti Fixed Recrut Treenbr fixe yearsbetween study offset
  try(Mrate.yr.QuantiF.TF.YO.Hnbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys)) + treeNbrJ +
                                                          M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.recruitment.plot.rate.yr.2,
                                                        zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QuantiF.TF.YO.Hnbinom2),silent=T) 
  
  ### Recrut rate Marginality Quali Aleatoire Recrut Treenbr fixe yearsbetween study offset 
  try(Mrate.yr.QualiA.TF.YO.Hnbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys)) + treeNbrJ +
                                                         (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.yr.2,
                                                       zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QualiA.TF.YO.Hnbinom2),silent=T) 
  
  ### Recrut rate Marginality Quanti Aleatoire Recrut Treenbr fixe yearsbetween study offset
  try(Mrate.yr.QuantiA.TF.YO.Hnbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys)) + treeNbrJ +
                                                          (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.yr.2,
                                                        zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QuantiA.TF.YO.Hnbinom2),silent=T) 
  
  
  
  ## Offset tree 
  ### Recrut rate Marginality Quali Fixed Recrut Treenbr offset yearsbetween study offset
  try(Mrate.yr.QualiF.TO.YO.Hnbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys)) + offset(log(treeNbrJ)) +
                                                         Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.recruitment.plot.rate.yr.2,
                                                       zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QualiF.TO.YO.Hnbinom2),silent=T) 
  
  ### Recrut rate Marginality Quanti Fixed Recrut Treenbr offset yearsbetween study offset
  try(Mrate.yr.QuantiF.TO.YO.Hnbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys)) + offset(log(treeNbrJ)) +
                                                          M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.recruitment.plot.rate.yr.2,
                                                        zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QuantiF.TO.YO.Hnbinom2),silent=T) 
  
  ### Recrut rate Marginality Quali Aleatoire Recrut Treenbr offset yearsbetween study offset 
  try(Mrate.yr.QualiA.TO.YO.Hnbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys)) + offset(log(treeNbrJ)) +
                                                         (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.yr.2,
                                                       zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QualiA.TO.YO.Hnbinom2),silent=T) 
  
  ### Recrut rate Marginality Quanti Aleatoire Recrut Treenbr offset yearsbetween study offset
  try(Mrate.yr.QuantiA.TO.YO.Hnbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys)) + offset(log(treeNbrJ)) +
                                                          (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.yr.2,
                                                        zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QuantiA.TO.YO.Hnbinom2),silent=T) 
  
  
  ## No tree parameter 
  ### Recrut rate Marginality Quali Fixed Recrut no treenbr yearsbetween study offset
  try(Mrate.yr.QualiF.noT.YO.Hnbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys))  +
                                                          Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.recruitment.plot.rate.yr.2,
                                                        zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QualiF.noT.YO.Hnbinom2),silent=T) 
  
  ### Recrut rate Marginality Quanti Fixed Recrut no treenbr yearsbetween study offset
  try(Mrate.yr.QuantiF.noT.YO.Hnbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys))  +
                                                           M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.recruitment.plot.rate.yr.2,
                                                         zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QuantiF.noT.YO.Hnbinom2),silent=T) 
  
  ### Recrut rate Marginality Quali AleanoTire Recrut no treenbr yearsbetween study offset 
  try(Mrate.yr.QualiA.noT.YO.Hnbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys))  +
                                                          (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.yr.2,
                                                        zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QualiA.noT.YO.Hnbinom2),silent=T) 
  
  ### Recrut rate Marginality Quanti AleanoTire Recrut no treenbr yearsbetween study offset
  try(Mrate.yr.QuantiA.noT.YO.Hnbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) + offset(log(yearsbetweensurveys))  +
                                                           (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.yr.2,
                                                         zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QuantiA.noT.YO.Hnbinom2),silent=T) 
  
  ## Yearsbetweenstudy not as a predictor
  ### Recrut rate Marginality Quali Fixed Recrut Treenbr fixe no years
  try(Mrate.yr.QualiF.TF.noY.Hnbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) + treeNbrJ +
                                                          Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.recruitment.plot.rate.yr.2,
                                                        zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QualiF.TF.noY.Hnbinom2),silent=T) 
  
  ### Recrut rate Marginality Quanti Fixed Recrut Treenbr fixe no years
  try(Mrate.yr.QuantiF.TF.noY.Hnbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) +  treeNbrJ +
                                                           M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.recruitment.plot.rate.yr.2,
                                                         zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QuantiF.TF.noY.Hnbinom2),silent=T) 
  
  ### Recrut rate Marginality Quali Aleatoire Recrut Treenbr fixe no years 
  try(Mrate.yr.QualiA.TF.noY.Hnbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) +  treeNbrJ +
                                                          (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.yr.2,
                                                        zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QualiA.TF.noY.Hnbinom2),silent=T) 
  
  ### Recrut rate Marginality Quanti Aleatoire Recrut Treenbr fixe no years
  try(Mrate.yr.QuantiA.TF.noY.Hnbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) +  treeNbrJ +
                                                           (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.yr.2,
                                                         zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QuantiA.TF.noY.Hnbinom2),silent=T) 
  
  
  
  ## Offset tree 
  ### Recrut rate Marginality Quali Fixed Recrut Treenbr offset no years
  try(Mrate.yr.QualiF.TO.noY.Hnbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) +  offset(log(treeNbrJ)) +
                                                          Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.recruitment.plot.rate.yr.2,
                                                        zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QualiF.TO.noY.Hnbinom2),silent=T) 
  
  ### Recrut rate Marginality Quanti Fixed Recrut Treenbr offset no years
  try(Mrate.yr.QuantiF.TO.noY.Hnbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) +  offset(log(treeNbrJ)) +
                                                           M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.recruitment.plot.rate.yr.2,
                                                         zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QuantiF.TO.noY.Hnbinom2),silent=T) 
  
  ### Recrut rate Marginality Quali Aleatoire Recrut Treenbr offset no years 
  try(Mrate.yr.QualiA.TO.noY.Hnbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) +  offset(log(treeNbrJ)) +
                                                          (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.yr.2,
                                                        zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QualiA.TO.noY.Hnbinom2),silent=T) 
  
  ### Recrut rate Marginality Quanti Aleatoire Recrut Treenbr offset no years
  try(Mrate.yr.QuantiA.TO.noY.Hnbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) +  offset(log(treeNbrJ)) +
                                                           (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.yr.2,
                                                         zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QuantiA.TO.noY.Hnbinom2),silent=T) 
  
  
  ## No tree parameter 
  ### Recrut rate Marginality Quali Fixed Recrut no treenbr no years
  try(Mrate.yr.QualiF.noT.noY.Hnbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) + 
                                                           Plotcat*mean_spei12 + Plotcat*BA.O.plot.1 + Plotcat*BAj.plot.1 + Plotcat*BAIj.plot.1.mean.J.I + Plotcat*sp.recruitment.plot.rate.yr.2,
                                                         zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QualiF.noT.noY.Hnbinom2),silent=T) 
  
  ### Recrut rate Marginality Quanti Fixed Recrut no treenbr no years
  try(Mrate.yr.QuantiF.noT.noY.Hnbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) + 
                                                            M.Plotcat*mean_spei12 + M.Plotcat*BA.O.plot.1 + M.Plotcat*BAj.plot.1 + M.Plotcat*BAIj.plot.1.mean.J.I + M.Plotcat*sp.recruitment.plot.rate.yr.2,
                                                          zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QuantiF.noT.noY.Hnbinom2),silent=T) 
  
  ### Recrut rate Marginality Quali AleanoTire Recrut no treenbr no years 
  try(Mrate.yr.QualiA.noT.noY.Hnbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) + 
                                                           (1|Plotcat) + mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.yr.2,
                                                         zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QualiA.noT.noY.Hnbinom2),silent=T) 
  
  ### Recrut rate Marginality Quanti AleanoTire Recrut no treenbr no years
  try(Mrate.yr.QuantiA.noT.noY.Hnbinom2 <- eval_fork(glmmTMB(sp.mortality.plot.rate.yr ~ dbhJ.IMall.plot.mean + (1|country) + 
                                                            (1|M.plotcat)+mean_spei12 + BA.O.plot.1 + BAj.plot.1 + BAIj.plot.1.mean.J.I + sp.recruitment.plot.rate.yr.2,
                                                          zi=~.,family="truncated_nbinom2", data=dfplot),timeout = 660),silent=T)
  try(Saving(Mrate.yr.QuantiA.noT.noY.Hnbinom2),silent=T) 
  A <- ls(pattern = "binom")
  
}

  
  #Load all + check all
  
  
  #### to be continued
  A <- dput(ls(pattern = test))
  B <- paste0("get")
  dput(AICtab(mbinom1,mbinom2,mHbinom1,mHbinom2))
  AICtab(dput(A))
  
  
  
  
  
  
  
  

  test <- glmmTMB(sp.recruitment.plot.rate.2 ~ dbhJ.IMall.plot.mean + (1|country) + yearsbetweensurveys + treeNbrJ +
            (1|Plotcat)*mean_spei12 + (1|Plotcat)*BA.O.plot.1 + (1|Plotcat)*BAj.plot.1 + (1|Plotcat)*BAIj.plot.1.mean.J.I + (1|Plotcat)*sp.mortality.plot.rate,
          zi=~.,family=nbinom1, data=dfplot)
try(Saving(Rrate.QualiF.TF.YF.nbinom1),silent=T) 































### For each single one except one give average values
data <- mHbinom1$frame
## Fixer tous les paramètres à leur moyennes 
## Sauf on veut 1/3 Plocat chaque cat
## Et on veut aussi faire varier un paramètre particulier (20 niveau)
myDF = as.data.frame(matrix(ncol = length(colnames(data)), nrow = 20*3, NA)) # Df vierge
colnames(myDF) <- colnames(data) # Bon noms de colonnes

#myDF[,1] <- as.numeric(rep(mean(data[,1]),20))
myDF[,2] <- as.numeric(rep(mean(data[,2]),20))
myDF[,3] <- as.numeric(rep(mean(data[,3]),20)) # country
myDF[,4] <- as.numeric(rep(mean(data[,4]),20))
myDF[,5] <- as.numeric(rep(mean(data[,5]),20))
myDF[,6] <- as.factor(c(rep("0",20),rep("1",20),rep("2",20)))
myDF[,7] <- as.numeric(rep(mean(data[,7]),20))
myDF[,8] <- as.numeric(rep(mean(data[,8]),20))
myDF[,9] <- as.numeric(rep(mean(data[,9]),20))
myDF[,10] <- as.numeric(rep(mean(data[,10]),20))
myDF[,11] <- as.numeric(rep(mean(data[,11]),20))

## Celui qu'on veut faire varier: 
name <- "BAIj.plot.1.mean.J.I"
myDF[,name] <- seq(min(data[,name]),max(data[,name]),length=20)


myDF$Predicted <- predict(
  mbinom2,
  newdata = myDF,
  newparams = NULL,
  se.fit = F,
  allow.new.levels = FALSE,
  type = c("response"),
  na.action = na.pass,
  debug = FALSE
)

p1<-ggplot(data=myDF, # All
           aes(x=get(name), y=Predicted, col=Plotcat,shape=Plotcat))+
  geom_line()+
  geom_point(size=2)

p1
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

Mbin1A <- fitme(sp.recruitment.plot.rate.yr ~ treeNbr + yearsbetweensurveys + 
                  sqrtBAIj.plot.1.mean + logdbh.plot.mean + mean_spei12 + logBAj.plot.1 + sqrtBA.O.plot.1 + Plotcat + 
                  mean_spei12:logBAj.plot.1 + mean_spei12:sqrtBA.O.plot.1 + mean_spei12:Plotcat + logBAj.plot.1:Plotcat + sqrtBA.O.plot.1:Plotcat + (1|country), 
                data=dfplot2, family=negbin(),method='REML')

Mbin2A <- fitme(sp.recruitment.plot.rate.yr ~ offset(treeNbr) + offset(yearsbetweensurveys) + 
                  sqrtBAIj.plot.1.mean + logdbh.plot.mean + mean_spei12 + logBAj.plot.1 + sqrtBA.O.plot.1 + Plotcat + 
                  mean_spei12:logBAj.plot.1 + mean_spei12:sqrtBA.O.plot.1 + mean_spei12:Plotcat + logBAj.plot.1:Plotcat + sqrtBA.O.plot.1:Plotcat + (1|country), 
                data=dfplot2, family=negbin(),method='REML')



# offset(sampl_effort)


##### Models #####
    
boxplot(dfplot2$BAIj.plot.1.mean~dfplot2$Plotcat)
plot(log(dfplot2[dfplot2$sp.mortality.plot.rate.yr!=0,"sp.mortality.plot.rate.yr"]+1)~log(dfplot2[dfplot2$sp.mortality.plot.rate.yr!=0,"BAIj.plot.1.mean"]+1))
coef <- lm(log(dfplot2[dfplot2$sp.mortality.plot.rate.yr!=0,"sp.mortality.plot.rate.yr"]+1)~log(dfplot2[dfplot2$sp.mortality.plot.rate.yr!=0,"BAIj.plot.1.mean"]+1))
abline(reg=coef)


coef <- lm(log(dfplot2[dfplot2$sp.mortality.plot.rate.yr!=0,"sp.mortality.plot.rate.yr"]+1)~log(dfplot2[dfplot2$sp.mortality.plot.rate.yr!=0,"BAIj.plot.1.mean"]+1) + dfplot2[dfplot2$sp.mortality.plot.rate.yr!=0,"Plotcat"])
coef
anova(coef)

RE <- dfplot2[dfplot2$Plotcat=="1",]
LE <- dfplot2[dfplot2$Plotcat=="2",]
C <- dfplot2[dfplot2$Plotcat=="0",]
plot(log(RE[RE$sp.mortality.plot.rate.yr!=0,"sp.mortality.plot.rate.yr"]+1)~log(RE[RE$sp.mortality.plot.rate.yr!=0,"BAIj.plot.1.mean"]+1))
coefC <- lm(log(RE[RE$sp.mortality.plot.rate.yr!=0,"sp.mortality.plot.rate.yr"]+1)~log(RE[RE$sp.mortality.plot.rate.yr!=0,"BAIj.plot.1.mean"]+1))
abline(reg=coefC,col="red")

plot(log(C[C$sp.mortality.plot.rate.yr!=0,"sp.mortality.plot.rate.yr"]+1)~log(C[C$sp.mortality.plot.rate.yr!=0,"BAIj.plot.1.mean"]+1))
coefRE <- lm(log(C[C$sp.mortality.plot.rate.yr!=0,"sp.mortality.plot.rate.yr"]+1)~log(C[C$sp.mortality.plot.rate.yr!=0,"BAIj.plot.1.mean"]+1))
abline(reg=coefRE,col="green")

plot(log(LE[LE$sp.mortality.plot.rate.yr!=0,"sp.mortality.plot.rate.yr"]+1)~log(LE[LE$sp.mortality.plot.rate.yr!=0,"BAIj.plot.1.mean"]+1))
coefLE <- lm(log(LE[LE$sp.mortality.plot.rate.yr!=0,"sp.mortality.plot.rate.yr"]+1)~log(LE[LE$sp.mortality.plot.rate.yr!=0,"BAIj.plot.1.mean"]+1))
abline(reg=coefLE,col="blue")


### Do the same for the recruitment once it is calculated


plot(log(dfplot2$)~log(dfplot2$BAIj.plot.1.mean+1))

plot(log(dfplot2$sp.mortality.plot.rate.yr+1)~log(dfplot2$recruitment.plot.rate+1))

coef <- lm(log(dfplot2$sp.mortality.plot.rate.yr+1)~log(dfplot2$BAIj.plot.1.mean+1))


anova(test)
summary(test)
test$coefficients
?abline
coef$coefficients
length(which(dfplot2$BAIj.plot.1.mean<0)) # Problem how can it be inferior to zero so often

RE <- dfplot2[dfplot2$Plotcat=="1",]
LE <- dfplot2[dfplot2$Plotcat=="2",]
C <- dfplot2[dfplot2$Plotcat=="0",]

plot(log(RE$sp.mortality.plot.rate.yr+1)~log(RE$BAIj.plot.1.mean+1))
plot(LE$sp.mortality.plot.rate.yr~LE$BAIj.plot.1.mean)
plot(log(C$sp.mortality.plot.rate.yr+1)~log(C$BAIj.plot.1.mean+1))
cor.test(RE$sp.mortality.plot.rate.yr,RE$BAIj.plot.1.mean)
cor.test(C$sp.mortality.plot.rate.yr,C$BAIj.plot.1.mean)
cor.test(LE$sp.mortality.plot.rate.yr,LE$BAIj.plot.1.mean)


head(dfplot2)

    
    #Mbin1A = bio1_climate_mean.30 ~ bio12_climate_mean.30 + sqrtBA.ha.plot.1 + sqrtBA.O.plot.1
    Model <- deparse(substitute(Mbin1A))
    try(Mbin1A <- eval_fork(fitme(sp.mort.bin ~ sqrtBAIj.plot.1 + sqrtBAIj.plot.1.mean + logdbh.plot.mean + treeNbr + yearsbetweensurveys + bio1_climate_mean.30 + bio12_climate_mean.30 + min_spei12 + mean_spei12 + sqrtBA.ha.plot.1 + sqrtBA.O.plot.1 + Plotcat + I(bio1_climate_mean.30^2) + I(bio12_climate_mean.30^2) + I(min_spei12^2) + I(mean_spei12^2) + bio1_climate_mean.30:bio12_climate_mean.30 + bio1_climate_mean.30:min_spei12 + bio1_climate_mean.30:mean_spei12 + bio1_climate_mean.30:sqrtBA.ha.plot.1 + bio1_climate_mean.30:sqrtBA.O.plot.1 + bio12_climate_mean.30:min_spei12 + bio12_climate_mean.30:mean_spei12 + bio12_climate_mean.30:sqrtBA.ha.plot.1 + bio12_climate_mean.30:sqrtBA.O.plot.1 + min_spei12:mean_spei12 + min_spei12:sqrtBA.ha.plot.1 + min_spei12:sqrtBA.O.plot.1 + mean_spei12:sqrtBA.ha.plot.1 + mean_spei12:sqrtBA.O.plot.1 + sqrtBA.ha.plot.1:sqrtBA.O.plot.1 +  + bio1_climate_mean.30:Plotcat + bio12_climate_mean.30:Plotcat + min_spei12:Plotcat + mean_spei12:Plotcat + sqrtBA.ha.plot.1:Plotcat + sqrtBA.O.plot.1:Plotcat + (1|country), data=dfplot2, family = binomial,method='REML'),timeout = 660),silent=T)
    try(Saving(Mbin1A),silent=T) 
    rm(list = c('Mbin1A')) 
    gc() 