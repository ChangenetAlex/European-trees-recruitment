### In this script we investigate the effects of simple effects on predicted count and predicted occurrence of recruitment
# One objective is to pinpoint threshold effects (on density for instance and also DBH for which main effects is oppposite in both part of the model)

rm(list = ls())
gc()
library(parallel)
library(forecast)
library(glmmTMB)
library(bbmle)
library(ggplot2)
library(reshape2)

Dir <- c("~/SynologyDrive/FUNDIV - NFI - Europe/")
setwd(Dir)
Allcode <- c("PINPINA","ABIALB","ACEPSE","ALNGLU","BETPEN","CASSAT","FAGSYL",
             "FRAEXC","PICABI","PINHAL","PINNIG","PINPIN","PINSYL","POPNIG",
             "POPTRE","QUEILE","QUEPET","QUEPYR","QUEROB","QUESUB")
Allseuil <- c(0.8,0.7,0.55,0.6,0.8,0.7,0.7,0.7,0.8,0.7,0.7,0.7,0.8,0.7,0.7,0.8,0.7,0.7,0.8,0.7)
AllseuilC <- c(0.6,0.6,0.5,0.5,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6)


####################################################
Lvl <- 30                                         ## Number of values for the interaction graphs
Modname <- "rW.yO.sumO.mRY1" ## Name of the model 
nameVAR <- c(                                     ##
  "sqrtBAIj.plot.1.mean.J.I",                     ##
  "mean_spei12",                                  ##
  "sqrtBA.O.plot.1",                              ##
  "logBAj.plot.1",                                ##
  "logtreeNbrJ",                                  ##
  "logdbhJ.IMall.plot.mean",                      ##
  "logsp.mortality.plot.rate.yr",                 ##
  "Plotcat")                                      ##
####################################################

### for all these variables we will simulate the response in CONd or ZI part along X values. ACCOUNTING NLY FOR VALUES THAT ARE IN THE RANGE 
mod_ALL <- list()

## First we load the model and dfplot 

for (i in c(1:13,16:length(Allcode))){
#for (i in c(1:2)){
    
  # for (i in c(1:3)){
  CODE <- Allcode[i]
  seuilC <- AllseuilC[i]
  seuil <- Allseuil[i]
  Dir = paste0("~/SynologyDrive/FUNDIV - NFI - Europe/our-data/species/",CODE,"/CLIMAP/Recrut.Mortality.2020/")
  setwd(Dir)

  #Load the desired model of the desired species 
  Dir <- c(paste0(Dir,Modname,"/")) # Directory 
  setwd(Dir)
  mod_ALL[[CODE]] <- get(load(file = paste0(Modname,".rda")))
  rm(x)
  
  ### Test #### 
  #we define the data as the part of the dataset on which the model has been fitted. (for every loop)
  #dfplot <- mod_ALL[[CODE]]$frame
  ## Test ##
}

TEST <- mclapply(mod_ALL, function(x){

### For each single one except one give average values
data <- x$frame

## We create a matrix of 30 values that we will predict by variable + 30 by level of plotcat. 
myDF = as.data.frame(matrix(ncol = length(colnames(data)), nrow = Lvl*length(nameVAR[nameVAR%in%colnames(data)]),NA)) 
colnames(myDF) <- colnames(data) # Bon noms de colonnes

#####
# 1 # Remplir tout à la moyenne et tous les Plotcat à 0 par défaut. 
#####
myDF[,nameVAR[nameVAR%in%colnames(data) & nameVAR!="Plotcat"]] <- apply(data[,nameVAR[nameVAR%in%colnames(data) & nameVAR!="Plotcat"]],2,
                                                                        function(y) as.numeric(rep(mean(y),nrow(myDF))))
myDF[,"Plotcat"] <- as.numeric(rep("0",nrow(myDF)))
if (length(unlist(ranef(x)))!=0){myDF[,"country"] <- NA} # Si random effect, put all values to NA. 


#####
# 2 # Tous les LvL*rows => Faire varier sur le range une variable. Sauf les LvL*2*dernières rows ou on fixe a 1 et 2 pour plotcat. 
#####

#Boucle imbriquée: pour chaque variable (qui n'est pas celle selectionnée en 2), 
# la variable est divisée en 30 (nameSeq) a répétée 3 fois par variable qui sera mis en face des 3 niveau de la variable name3
COUNT <- 0  
for (i in 1:(length(nameVAR[nameVAR%in%colnames(data)]))){  # Il y aura autant de tour de boucle que de variable. Et 1 next  
  nameSeq <- nameVAR[nameVAR%in%colnames(data)][i]
  if (nameSeq=="Plotcat"){ # sila variable est plotcat on passe a la suite
    COUNT <- COUNT+1
    next} #idem 
  else {myDF[(Lvl*(i-1)+1-(Lvl*COUNT)):(Lvl*(i)-(Lvl*COUNT)),nameSeq] <- seq(min(data[,nameSeq]),max(data[,nameSeq]),length=Lvl)}
}


#####
# 3 # Pour les dernière lignes répetitionn of 0,1,2
#####
myDF[(1+nrow(myDF)-Lvl):nrow(myDF),"Plotcat"] <- as.numeric(c(rep("0",Lvl/3),rep("1",Lvl/3),rep("2",Lvl/3)))

## Added on the 15 august : change the name of the offset column by the real name of the variable and put it to the exponential (delog)
myDF[,"offset(log(yearsbetweensurveys))"] <- mean(exp(data[,"offset(log(yearsbetweensurveys))"]))
myDF[,"offset(log(sp.SUM.ha.weight2))"] <- mean(exp(data[,"offset(log(sp.SUM.ha.weight2))"]))

colnames(myDF)[colnames(myDF)=="offset(log(yearsbetweensurveys))"] <- "yearsbetweensurveys" 
colnames(myDF)[colnames(myDF)=="offset(log(sp.SUM.ha.weight2))"] <- "sp.SUM.ha.weight2"


## Predicted response # Maybe add an arument re.form = NA
Predictedresp <- predict(x,newdata = myDF,se.fit = T,type = c("response"))
## Predicted mu (count part)
Predictedmu <- predict(x,newdata = myDF,se.fit = T,type = c("conditional"))
## Predicted probability to observe a zero event (= proba of a non recruit event) 
Predictedzprob <- predict(x,newdata = myDF,se.fit = T,type = c("zprob"))
## Predicted link (count part)
Predictedmulink <- predict(x,newdata = myDF,se.fit = T,type = c("link"))
## Predicted probability to observe a zero event (= proba of a non recruit event) #link scaled 
Predictedzlink <- predict(x,newdata = myDF,se.fit = T,type = c("zlink"))

myDF$Predictedresp <- Predictedresp$fit
myDF$SEresp <- Predictedresp$se.fit  
myDF$Predictedmu <- Predictedmu$fit
myDF$SEmu <- Predictedmu$se.fit
myDF$Predictedzprob <- Predictedzprob$fit
myDF$SEzprob <- Predictedzprob$se.fit
myDF$Predictedzlink <- Predictedzlink$fit
myDF$SEzlink <- Predictedzlink$se.fit
myDF$Predictedmulink <- Predictedmulink$fit
myDF$SEmulink <- Predictedmulink$se.fit


### Split all data in one dataset par variable 
# COUNT <- 0  
# for (i in 1:(length(nameVAR[nameVAR%in%colnames(data)]))){  # Il y aura autant de tour de boucle que de variable. Et 1 next  
#   nameSeq <- nameVAR[nameVAR%in%colnames(data)][i]
#   if (nameSeq!="Plotcat"){ # sila variable est plotcat on passe a la suite
#     assign(paste0(nameSeq,".df"),myDF[(Lvl*(i-1)+1-(Lvl*COUNT)):(Lvl*(i)-(Lvl*COUNT)),])
#   }else if (nameSeq=="Plotcat"){assign(paste0(nameSeq,".df"),myDF[(1+nrow(myDF)-Lvl):nrow(myDF),])}
# }
;return(myDF)},mc.cores = 18)


TEST2 <- mclapply(TEST, function(x) {
if (length(which(colnames(x)=="country")!=0)){
    x <- x[,-which(colnames(x)=="country")]
  }
  ;x},mc.cores = 18)

## Convert it to a df
TEST2 <- mapply(function(x,y) {
  x[,"species"] <- y
  return(x)
},
x=TEST2,y=names(TEST2),SIMPLIFY = F)

#Save these mean and SD predictions as well as the data that were used to generate it
saveRDS(get("TEST2"), paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/paper2/myDF.SIMPLE.EFFECT.SIM_",Lvl,".rds"))










