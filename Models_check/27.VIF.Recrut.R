rm(list = ls())
gc()
library(lattice)
library(tidyr)
library(usdm)

#### 1 #### Load our data 
Dir <- c("~/SynologyDrive/FUNDIV - NFI - Europe/")
setwd(Dir)
source(paste0(Dir,"Myscripts/Fundiv.project/RecrutProject.function2.Extract.R"))

##########################################################
Modname <- c("SCALEDRrate.yr.QualiF.TF.noY.nbinom1.V3") ## Name of the model 
##########################################################

Allcode <- c("PINPINA","ABIALB","ACEPSE","ALNGLU","BETPEN","CASSAT","FAGSYL",
             "FRAEXC","PICABI","PINHAL","PINNIG","PINPIN","PINSYL","POPNIG",
             "POPTRE","QUEILE","QUEPET","QUEPYR","QUEROB","QUESUB")
Allseuil <- c(0.8,0.7,0.55,0.6,0.8,0.7,0.7,0.7,0.8,0.7,0.7,0.7,0.8,0.7,0.7,0.8,0.7,0.7,0.8,0.7)
AllseuilC <- c(0.6,0.6,0.5,0.5,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6)


### Function ta calculate VIF on all species and record the results directly in the file where the model is recorded
VIF.R <- function(Modname){
  VIFall <- list()
  for (i in c(1:13,16:length(Allcode))){
    CODE <- Allcode[i]
    seuilC <- AllseuilC[i]
    seuil <- Allseuil[i]
    Dir = paste0("~/SynologyDrive/FUNDIV - NFI - Europe/our-data/species/",CODE,"/Recrut.Mortality.2020/")
    setwd(Dir)
    dfplot <- try(readRDS(paste0("dfplotFinalV2",CODE,seuil,"R.M.rds")),silent = T) #Base de données plot
    # Change directory the the model we want to investigate 
    Dir <- c(paste0(Dir,Modname,"/")) # Directory 
    setwd(Dir)
    
    # Charge the model 
    try(assign(Modname,get(load(file = paste0(Modname,".rda")))),silent = F)
    rm(x)
    data <- get(Modname)$frame
    data$M.Plotcat <- dfplot$M.Plotcat[match(rownames(data), rownames(dfplot),incomparables = NA)] # Obtain the values of continous marginality from dfplot
    
    Explain <- Extraction.R(get(Modname)) # Use the new extraction function to obtain the variables we need 
    Explain[Explain=="Plotcat"] <- "M.Plotcat" # Replace the plotcat categorical variable by the continous one
    Resp <- sub(" ","", unlist(strsplit(paste0(getCall(get(Modname))[2]), "~",fixed = T))[1], ignore.case = FALSE,fixed = T) # take the response
    
    VIF <- usdm::vif(data[,Explain])
    VIFSTEP <- vifstep(data[,Explain],th=8)
    colnames(VIF)[1] <- "Variables"
    colnames(VIF)[2] <- paste0(CODE)
    # Create the directory 
    
    Dir <- paste0("~/SynologyDrive/FUNDIV - NFI - Europe/our-data/species/Recrut.All/",Modname,"/")
    setwd(Dir)
    VIFall[[CODE]] <- VIF[2] 
    capture.output(print(VIF),file = paste0(Modname,"_All.Vif.txt"),append=T)
    capture.output(c(CODE,VIFSTEP),file = paste0(Modname,"_All.Vifstep.txt"),append=T)
  }
  VIFall <- as.data.frame(VIFall)
  rownames(VIFall) <- VIF[,1]
  VIFall <- t(VIFall)
  write.table(round(VIFall,2),file = paste0("VIFall.csv"),row.names = T,sep = ";",col.names = NA)
}

