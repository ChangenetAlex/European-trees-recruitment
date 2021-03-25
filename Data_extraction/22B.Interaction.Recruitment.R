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
source(paste0(Dir,"Myscripts/Fundiv.project/RecrutProject.function1.Saving.R"))
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
for (i in c(1:13,16:length(Allcode))){
#for (i in c(1)){
    
  CODE <- Allcode[i]
  seuilC <- AllseuilC[i]
  seuil <- Allseuil[i]
  # Dir = paste0("/home/achangenet/SynologyDrive/FUNDIV - NFI - Europe/our-data/species/",CODE,"/CLIMAP/Recrut.Mortality.2020/")
  Dir = paste0("~/SynologyDrive/FUNDIV - NFI - Europe/our-data/species/",CODE,"/CLIMAP/Recrut.Mortality.2020/")
  
  setwd(Dir)
  ## Charge the dfplot that is located in the species recrut file 
  dfplot <- try(readRDS(paste0("dfplotFinalV3",CODE,seuil,"R.M.rds")),silent = T) #Base de données plot
  
  # Also load dfplo2 or three
  
  # Change directory the the model we want to investigate 
  Dir <- c(paste0(Dir,Modname,"/")) # Directory 
  setwd(Dir)
  
  # Add recruitment as a count per year and as a total count
  # dfplot$sp.recruitment.plot.rate.yr.2 <- round((dfplot$sp.recruitment.plot.rate * 1000 / dfplot$yearsbetweensurveys))
  # dfplot$sp.recruitment.plot.rate.2 <- round(dfplot$sp.recruitment.plot.rate * 1000)
  # dfplot$sp.mortality.plot.rate.2 <- round(dfplot$sp.mortality.plot.rate * 1000)
  
  # Charge the model 
  try(assign(Modname,get(load(file = paste0(Modname,".rda")))),silent = F)
  rm(x)
  Dir <- c("~/SynologyDrive/FUNDIV - NFI - Europe/our-data/species/Recrut.All/")
  dir.create(paste0(Dir,Modname,"/")) # Creation of a new file for this model
  Dir <- paste0(Dir,Modname,"/")
  
  ### For each single one except one give average values
  data <- get(Modname)$frame
  myDF = as.data.frame(matrix(ncol = length(colnames(data)), nrow = Lvl*3*(length(nameVAR[nameVAR%in%colnames(data)])-1)*(length(nameVAR[nameVAR%in%colnames(data)])-1), NA)) # Df vierge
  colnames(myDF) <- colnames(data) # Bon noms de colonnes
  
  #####
  # 1 # Remplir tout à la moyenne et tous les Plotcat à 0 par défaut. 
  #####
  
  myDF[,nameVAR[nameVAR%in%colnames(data) & nameVAR!="Plotcat"]] <- apply(data[,nameVAR[nameVAR%in%colnames(data) & nameVAR!="Plotcat"]],2,
                                                                          function(x) as.numeric(rep(mean(x),nrow(myDF))))
  myDF[,"Plotcat"] <- as.numeric(rep("0",nrow(myDF)))
  if (length(unlist(ranef(get(Modname))))!=0){myDF[,"country"] <- NA} # Si random effect, put all values to NA. 
  
  
  #####
  # 2 # Variable que l'on veut regarder (fixer et faire varier)
  #####
  
  # Si la variable est Plotcat alors répétition des 3 niveau X nombre d'interaction possibles = mon nombre de variable -1
  # Si la variable est une autre alors répétition en 3 niveaux (low, high, moyenne) X nombre d'interaction possibles = mon nombre de variable -1
  
  for (j in 1:length(nameVAR[nameVAR%in%colnames(data)])){  # Il y aura autant de tour de boucle que de variable. 
    name3 <- nameVAR[nameVAR%in%colnames(data)][j]
    if (name3=="Plotcat"){myDF[((Lvl*3)*(j-1)*(length(nameVAR[nameVAR%in%colnames(data)])-2)+1):(Lvl*3*j*(length(nameVAR[nameVAR%in%colnames(data)])-2))
                               ,name3] <- as.numeric(c(rep("0",Lvl),rep("1",Lvl),rep("2",Lvl)))
    } else {myDF[((Lvl*3)*(j-1)*(length(nameVAR[nameVAR%in%colnames(data)])-2)+1):(Lvl*3*j*(length(nameVAR[nameVAR%in%colnames(data)])-2))
                 ,name3] <- c(as.numeric(rep(mean(data[,name3]),Lvl)),
                              rep(as.numeric(quantile(data[,name3],0.995)),length=Lvl),
                              rep(as.numeric(quantile(data[,name3],0.005)), length=Lvl))}
    
    #####
    # 3 # Toutes les interactions. (Autant de tour de boucle que d'interaction possible par variable == nbr de variable -1)
    #####
    
    #Boucle imbriquée: pour chaque variable (qui n'est pas celle selectionnée en 2), 
    # la variable est divisée en 30 (nameSeq) a répétée 3 fois par variable qui sera mis en face des 3 niveau de la variable name3
    COUNT <- 0  
    for (i in 1:(length(nameVAR[nameVAR%in%colnames(data)]))){  # Il y aura autant de tour de boucle que de variable. Et 1 next  
      nameSeq <- nameVAR[nameVAR%in%colnames(data)][i]
      if (nameSeq==name3){
        COUNT <- COUNT+1
        next} # lui dire de sauter Lvl ligne 
      else if (nameSeq=="Plotcat"){
        COUNT <- COUNT+1
        next} #idem 
      else {myDF[((Lvl*3)*(j-1)*(length(nameVAR[nameVAR%in%colnames(data)])-2)+((Lvl*3)*(i-1)+1)-(Lvl*3*COUNT)):
                   ((Lvl*3)*(j-1)*(length(nameVAR[nameVAR%in%colnames(data)])-2)+(Lvl*3*i)-(Lvl*3*COUNT)),
                 nameSeq] <- seq(min(data[,nameSeq]),max(data[,nameSeq]),length=Lvl)}
    }
  }
  
  #####
  # 4 # Pour le dernier set d'interactions : ajouter 
  #####
  myDF[(1+nrow(myDF)-Lvl*3):nrow(myDF),"Plotcat"] <- as.numeric(c(rep("0",Lvl),rep("1",Lvl),rep("2",Lvl)))
  
  ## Added on the 15 august : change the name of the offset column by the real name of the variable and put it to the exponential (delog)
  myDF[,"offset(log(yearsbetweensurveys))"] <- mean(exp(data[,"offset(log(yearsbetweensurveys))"]))
  myDF[,"offset(log(sp.SUM.ha.weight2))"] <- mean(exp(data[,"offset(log(sp.SUM.ha.weight2))"]))
  
  colnames(myDF)[colnames(myDF)=="offset(log(yearsbetweensurveys))"] <- "yearsbetweensurveys" 
  colnames(myDF)[colnames(myDF)=="offset(log(sp.SUM.ha.weight2))"] <- "sp.SUM.ha.weight2"
  
  
  ## Predicted response # Maybe add an arument re.form = NA
  Predictedresp <- predict(get(Modname),newdata = myDF,se.fit = T,type = c("response"))
  ## Predicted mu (count part)
  Predictedmu <- predict(get(Modname),newdata = myDF,se.fit = T,type = c("conditional"))
  ## Predicted probability to observe a zero event (= proba of a non recruit event) 
  Predictedzprob <- predict(get(Modname),newdata = myDF,se.fit = T,type = c("zprob"))
  
  myDF$Predictedresp <- Predictedresp$fit
  myDF$SEresp <- Predictedresp$se.fit  
  myDF$Predictedmu <- Predictedmu$fit
  myDF$SEmu <- Predictedmu$se.fit
  myDF$Predictedzprob <- Predictedzprob$fit
  myDF$SEzprob <- Predictedzprob$se.fit
  
  #Save these mean and SD predictions as well as the data that were used to generate it
  # saveRDS(get("myDF"), paste0("/home/achangenet/SynologyDrive/FUNDIV - NFI - Europe/our-data/species/Recrut.All/",Modname,"/myDF.SIM_",Lvl,"_", CODE,".rds"))
  saveRDS(get("myDF"), paste0("~/SynologyDrive/FUNDIV - NFI - Europe/our-data/species/Recrut.All/",Modname,"/myDF.SIM_All_",Lvl,"_", CODE,".rds"))
  
  yMAX.resp <- max(myDF[,"Predictedresp"])+1.96*max(myDF[,"SEresp"]) # Here is the maximum values that will be the limit for the plot
  yMAX.mu <- max(myDF[,"Predictedmu"])+1.96*max(myDF[,"SEmu"]) # Here is the maximum values that will be the limit for the plot
  yMAX.resp <- max(myDF[,"Predictedzprob"])+1.96*max(myDF[,"SEzprob"]) # Here is the maximum values that will be the limit for the plot
  
}







