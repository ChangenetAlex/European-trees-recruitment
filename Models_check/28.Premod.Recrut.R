rm(list = ls())
gc()
dev.off()

Dir <- c("~/SynologyDrive/FUNDIV - NFI - Europe/")
setwd(Dir)
source(paste0(Dir,"Myscripts/Fundiv.project/RecrutProject.function2.Extract.R"))
source(paste0(Dir,"Myscripts/Fundiv.project/RecrutProject.function3.Premodel.R"))

#######################################################
Modname <- "SCALEDRrate.yr.QualiF.TF.noY.nbinom1.V3" ## Name of the model 
#######################################################

Allcode <- c("PINPINA","ABIALB","ACEPSE","ALNGLU","BETPEN","CASSAT","FAGSYL",
             "FRAEXC","PICABI","PINHAL","PINNIG","PINPIN","PINSYL","POPNIG",
             "POPTRE","QUEILE","QUEPET","QUEPYR","QUEROB","QUESUB")
Allseuil <- c(0.8,0.7,0.55,0.6,0.8,0.7,0.7,0.7,0.8,0.7,0.7,0.7,0.8,0.7,0.7,0.8,0.7,0.7,0.8,0.7)
AllseuilC <- c(0.6,0.6,0.5,0.5,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6)
# Load the function 


for (k in c(1:13,16:length(Allcode))){
# for (k in c(length(Allcode))){
    
  CODE <- Allcode[k]
  seuilC <- AllseuilC[k]
  seuil <- Allseuil[k]
  Dir = paste0("/Users/alexandrechangenet/SynologyDrive/FUNDIV - NFI - Europe/our-data/species/",CODE,"/Recrut.Mortality.2020/")
  setwd(Dir)
  dfplot <- try(readRDS(paste0("dfplotFinalV2",CODE,seuil,"R.M.rds")),silent = T) #Base de données plot
  # Change directory the the model we want to investigate 
  Dir <- c(paste0(Dir,Modname,"/")) # Directory 
  setwd(Dir)
  
  # Charge the model 
  try(assign(Modname,get(load(file = paste0(Modname,".rda")))),silent = F)
  rm(x)
  data <- get(Modname)$frame
  Explain <- Extraction.R(get(Modname)) # Use the new extraction function to obtain the variables we need 
  Explain <- c(Explain[-which(Explain=="Plotcat")],"Plotcat") # Get Plotcat at the last position 
  Resp <- sub(" ","", unlist(strsplit(paste0(getCall(get(Modname))[2]), "~",fixed = T))[1], ignore.case = FALSE,fixed = T) # take the response
  
  ###########################
  ## Premodel analysis ###### 
  ###########################
  
  # Set directory to a common file for all models 
  Dir <- paste0("/Users/alexandrechangenet/SynologyDrive/FUNDIV - NFI - Europe/our-data/species/Recrut.All/",Modname,"/")
  setwd(Dir)
  
  # Perform the anaysis on all variable 4 and 3 with Plotcat at last
    Explain2 <- Explain[c(1:4,length(Explain))] # Apply the function on the first 4 variables 
    print(Explain2)
    Premodel(z=data,Resp=Resp,Explain=Explain2,size=3,save=T)
    Explain2 <- Explain[c(5:7,length(Explain))] # Apply it on the three remaininbg variables 
    print(Explain2)
    Premodel(z=data,Resp=Resp,Explain=Explain2,size=2,save=T)
}
  
