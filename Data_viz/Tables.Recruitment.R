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

Modname <- "SCALEDRrate.yr.QualiF.TF.noY.nbinom1.V3" ## Name of the model 
Modname <- "rW.yO.sumO.mRY1"
Allcode <- c("ABIALB","ACEPSE","ALNGLU","BETPEN","CASSAT","FAGSYL",
             "FRAEXC","PICABI","PINHAL","PINNIG","PINPIN","PINPINA","PINSYL","POPNIG",
             "POPTRE","QUEILE","QUEPET","QUEPYR","QUEROB","QUESUB")
Allseuil <- c(0.7,0.55,0.6,0.8,0.7,0.7,0.7,0.8,0.7,0.7,0.7,0.8,0.8,0.7,0.7,0.8,0.7,0.7,0.8,0.7)
AllseuilC <- c(0.6,0.5,0.5,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6)
mod_ALL <- list()
dfplot <- list()
dftree <- list()
for (i in c(1:13,16:length(Allcode))){
  # for (i in c(1:3)){
  CODE <- Allcode[i]
  seuilC <- AllseuilC[i]
  seuil <- Allseuil[i]
  Dir = paste0("~/SynologyDrive/FUNDIV - NFI - Europe/our-data/species/",CODE,"/CLIMAP/Recrut.Mortality.2020/")
  setwd(Dir)
  # Charge the dfplot for simulation and residuals test dispersion etc... 
  #dfplot[[CODE]] <- readRDS(paste0("dfplotFinalV2",CODE,seuil,"R.M.rds")) # Final database scaled and trasnformed (cf script 10second recruit and 10second recruit edited)
  dfplot[[CODE]] <- readRDS(paste0("dfplotV2", CODE, seuil, "R.M.rds")) #Base de données plot full No scale !!!
  dftree[[CODE]] <- readRDS(paste0("Mydf3_",CODE,"_",seuil,"_",seuilC,"R.M.rds"))
  #Load the desired model of the desired species 
  Dir <- c(paste0(Dir,Modname,"/")) # Directory 
  setwd(Dir)
  mod_ALL[[CODE]] <- get(load(file = paste0(Modname,".rda")))$frame
  rm(x)
}

# Here we obtain a table with the two qualitative variables and the total number of plot
A <- mapply(function(x,xx) {
  if (!"country"%in%colnames(x)){
    # Si le country est vide <- mettre le niveau de celui ci (ni = 0 dans df plot ni = FR)
    x[,"country"] <- levels(xx$country)[table(xx$country)!=0 & levels(xx$country)!="FR"]
  }
  x[,"country"] <- factor(x[,"country"],levels=levels(xx[,"country"]))
  return(x)
},x=mod_ALL,xx=dfplot,SIMPLIFY = F)

# Here table of the number of years (mean, mini and maxi)
D <- lapply(dfplot,function(x) {
  c("Census interval"=paste0(round(mean(x[x$country!="FR","yearsbetweensurveys"]),2),
           "\n(",min(x[x$country!="FR","yearsbetweensurveys"]),
           "-",max(x[x$country!="FR","yearsbetweensurveys"]),')'))
})

D <- do.call(rbind,D)

# Here table of country and marginality 
B <- lapply(A,function(x) {
  c("Total" = nrow(x),table(x$country),table(x$Plotcat))})

B <- do.call(rbind,B)

# Need number total of trees (minus france) + nuèmber of recruit over the number total of trees 
dftreeTot <- lapply(dftree,function(x) {
  c("NTrees"=nrow(x[x$country!="FR",]),"NRecrut (and %)"=paste0(nrow(x[x$treestatus_th=="1" & x$country!="FR",]),"\n(",
    round((nrow(x[x$treestatus_th=="1" & x$country!="FR",])/nrow(x[x$country!="FR",]))*100,2),"%)"))
})
dftreeTot <- do.call(rbind,dftreeTot)
B <- cbind(dftreeTot,B,D)
B <- as.data.frame(B)
unlist(B)
write.csv(B,file = paste0("/home/achangenet/SynologyDrive/FUNDIV - NFI - Europe/Redaction/paper2/Qualitable",Modname,".csv"))



## 07/08/2020 => Add a table with the average dbh of the recruited trees by countrys !!!!!! To save for later
dftreebis <- lapply(dftree,function(x){
  x <- x[x$country!="FR" & x$treestatus_th==1,c("country","dbh2","yearsbetweensurveys")]
})
dftreebis <- do.call(rbind,dftreebis)
dftreebis <- split(dftreebis,dftreebis$country)


dftreethird <- lapply(dftreebis,function(x) {
  c("dbh recruited"=paste0(round(mean(x[,"dbh2"]),2),
                             "\n(",round(min(x[,"dbh2"]),2),
                             "-",round(max(x[,"dbh2"]),2),')'),
    "census interval"=paste0(round(mean(x[,"yearsbetweensurveys"]),2),
                           "\n(",round(min(x[,"yearsbetweensurveys"]),2),
                           "-",round(max(x[,"yearsbetweensurveys"]),2),')'))
})
dftreethird <- do.call(rbind,dftreethird)


## Compare with the average of all trees: 
dftreebis <- lapply(dftree,function(x){
  x <- x[x$country!="FR",c("country","dbh1","dbh2","yearsbetweensurveys")]
})
dftreebis <- do.call(rbind,dftreebis)
dftreebis <- split(dftreebis,dftreebis$country)


dftreefourth <- lapply(dftreebis,function(x) {
  c("dbh1 all"=paste0(round(mean(x[,"dbh1"],na.rm = T),2),
                           "\n(",round(min(x[,"dbh1"],na.rm = T),2),
                           "-",round(max(x[,"dbh1"],na.rm = T),2),')'),
    "dbh2 all"=paste0(round(mean(x[,"dbh2"],na.rm = T),2),
                      "\n(",round(min(x[,"dbh2"],na.rm = T),2),
                      "-",round(max(x[,"dbh2"],na.rm = T),2),')'))
})
dftreefourth <- do.call(rbind,dftreefourth)

DBH <- cbind(dftreefourth,dftreethird)
write.csv2(DBH,file = paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/paper2/DBH.",Modname,".csv"))


# Need to record in a table 
# Need to do that for the individuals that correpond in the dfplot database
# Need to check if 
C <- mapply(function (x,xx){
  apply(x[which(rownames(x)%in%rownames(xx)),c("sp.recruitment.plot.rate.yr.2", "dbhJ.IMall.plot.mean", 
             "treeNbrJ","mean_spei12", "BA.O.plot.1", 
             "BAj.plot.1", "BAIj.plot.1.mean.J.I",
             "sp.mortality.plot.rate.yr")],2,function(y){paste0(round(mean(y),2),"\n(",round(sd(y),2),")\n",round(min(y),2),"-",round(max(y),2))})
  },xx=mod_ALL,x=dfplot,SIMPLIFY = F)

C <- do.call(rbind,C)

Rweight <- lapply(mod_ALL,function(x){
  paste0(round(mean(x[,1]),2),"\n(",round(sd(x[,1]),2),")\n",round(min(x[,1]),2),"-",round(max(x[,1]),2))})
Rweight <- do.call(rbind,Rweight)
D <- cbind(C,Rweight)
colnames(D)[9] <- "sp.recrut.weight"
write.csv2(D,file = paste0("/home/achangenet/SynologyDrive/FUNDIV - NFI - Europe/Redaction/paper2/QuantitableV2",Modname,".csv"))

