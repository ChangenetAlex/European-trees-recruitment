## Edited version of the 18/03/2020 to add all calculated variables 
# + marginality continous + scaled variables without centering and do 
# it also with mortality and treenbr to normalize them
rm(list = ls())
gc()
Dir <- c("/Users/alexandrechangenet/SynologyDrive/FUNDIV - NFI - Europe/")
setwd(Dir)
# Load the function 
source(paste0(Dir,"Myscripts/Fundiv.project/10second.Recruit.Calculation.EDITED.R"))
Allcode <- c("PINPINA","ABIALB","ACEPSE","ALNGLU","BETPEN","CASSAT","FAGSYL","FRAEXC","PICABI","PINHAL","PINNIG","PINPIN","PINSYL","POPNIG","POPTRE","QUEILE","QUEPET",
             "QUEPYR","QUEROB","QUESUB")
Allseuil <- c(0.8,0.7,0.55,0.6,0.8,0.7,0.7,0.7,0.8,0.7,0.7,0.7,0.8,0.7,0.7,0.8,0.7,0.7,0.8,0.7)
AllseuilC <- c(0.6,0.6,0.5,0.5,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6)

Dir <- c("/Users/alexandrechangenet/SynologyDrive/FUNDIV - NFI - Europe/")
setwd(Dir)
RecrutCalculations <- file(paste0(Dir,"our-data/species/RecrutCalculations.Rout"), open="wt")
sink(RecrutCalculations, type="output")

for (i in 1:length(Allcode)){
  CODE <- Allcode[i]
  seuilC <- AllseuilC[i]
  seuil <- Allseuil[i]
  try(MyDFs(CODE,seuil,seuilC),silent=T)
}
close(RecrutCalculations) 


###########################################################################################################


#### First version of the script (not all variable, 
# marginality to add below, centered and scaled on 
# all variable except treeNbrJ and mortality)

rm(list = ls())
gc()
Dir <- c("/home/achangenet/SynologyDrive/FUNDIV - NFI - Europe/")
setwd(Dir)
# Load the function 
source(paste0(Dir,"Myscripts/Fundiv.project/10second.Recruit.Calculation.R"))
Allcode <- c("PINPINA","ABIALB","ACEPSE","ALNGLU","BETPEN","CASSAT","FAGSYL","FRAEXC","PICABI","PINHAL","PINNIG","PINPIN","PINSYL","POPNIG","POPTRE","QUEILE","QUEPET",
             "QUEPYR","QUEROB","QUESUB")
Allseuil <- c(0.8,0.7,0.55,0.6,0.8,0.7,0.7,0.7,0.8,0.7,0.7,0.7,0.8,0.7,0.7,0.8,0.7,0.7,0.8,0.7)
AllseuilC <- c(0.6,0.6,0.5,0.5,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6)

Dir <- c("/home/achangenet/SynologyDrive/FUNDIV - NFI - Europe/")
setwd(Dir)
RecrutCalculations <- file(paste0(Dir,"our-data/species/RecrutCalculations.Rout"), open="wt")
sink(RecrutCalculations, type="output")

for (i in 1:length(Allcode)){
  CODE <- Allcode[i]
  seuilC <- AllseuilC[i]
  seuil <- Allseuil[i]
  try(MyDFs(CODE,seuil,seuilC),silent=T)
}
close(RecrutCalculations) 


#### Add continuous marginality index in the dfplot 

rm(list = ls())
gc()

#### 1 #### Load our data 
Dir <- c("/home/achangenet/SynologyDrive/FUNDIV - NFI - Europe/")
setwd(Dir)


i = 1
Allcode <- c("PINPINA","ABIALB","ACEPSE","ALNGLU","BETPEN","CASSAT","FAGSYL","FRAEXC","PICABI","PINHAL","PINNIG","PINPIN","PINSYL","POPNIG","POPTRE","QUEILE","QUEPET",
             "QUEPYR","QUEROB","QUESUB")
Allseuil <- c(0.8,0.7,0.55,0.6,0.8,0.7,0.7,0.7,0.8,0.7,0.7,0.7,0.8,0.7,0.7,0.8,0.7,0.7,0.8,0.7)
AllseuilC <- c(0.6,0.6,0.5,0.5,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6)

for (i in 1:length(Allcode)){
  Dir = (paste0("/home/achangenet/SynologyDrive/FUNDIV - NFI - Europe/our-data/species/",Allcode[i])) # Directory 
  setwd(Dir)
  assign(paste0("dfplot"),readRDS(paste0("CLIMAP/Recrut.Mortality.2020/dfplot2",Allcode[i],Allseuil[i],"R.M.rds"))) #Base de données plot
  load(paste0(Dir,"/CLIMAP/CLIMAP2019/Climap_",Allcode[i],"_",Allseuil[i],"_",AllseuilC[i],".RData"),envir = globalenv())
  ras <- rasterFromXYZ(test6[,c(1:2,26)],crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0 ") #New raster with our three categories 
  M.Plotcat <- raster::extract(ras,dfplot[,2:3])
  Myplot = cbind(dfplot,M.Plotcat)
  saveRDS(get("Myplot"), paste0(Dir, "/CLIMAP/Recrut.Mortality.2020/dfplotFinal", Allcode[i], Allseuil[i], "R.M.rds"))
  
}

