### Function to put together all summary of all species together

SavingInfo = " This function just require an argument Modname:
Model.All( 
Modname = 'the model name comon for all species' \n
Example : 
Model.All(x = 'SCALEDRrate.yr.QualiF.TF.noY.nbinom1.V3')\n"
bannerBreak = "\n*********************************************************************\n"
cat(paste0(bannerBreak,SavingInfo,bannerBreak,"\n"))


rm(list = ls())
gc()
Dir <- c("~/SynologyDrive/FUNDIV - NFI - Europe/")
setwd(Dir)
Allcode <- c("PINPINA","ABIALB","ACEPSE","ALNGLU","BETPEN","CASSAT","FAGSYL",
             "FRAEXC","PICABI","PINHAL","PINNIG","PINPIN","PINSYL","POPNIG",
             "POPTRE","QUEILE","QUEPET","QUEPYR","QUEROB","QUESUB")
Allseuil <- c(0.8,0.7,0.55,0.6,0.8,0.7,0.7,0.7,0.8,0.7,0.7,0.7,0.8,0.7,0.7,0.8,0.7,0.7,0.8,0.7)
AllseuilC <- c(0.6,0.6,0.5,0.5,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6)
Modname <- "SCALEDRrate.yr.QualiF.TF.noY.nbinom1.V3" ## Name of 
REV <- c("PINPIN","PINNIG","QUEPET","QUEROB","PINPINA","POPTRE","ALNGLU","QUEPYR","QUESUB")

Model.All <- function(Modname)
{
  for (j in c(1:13,16:length(Allcode))){
    CODE <- Allcode[j]
    seuilC <- AllseuilC[j]
    seuil <- Allseuil[j]
    # Dir = paste0("/home/achangenet/SynologyDrive/FUNDIV - NFI - Europe/our-data/species/",CODE,"/CLIMAP/Recrut.Mortality.2020/")
    Dir = paste0("~/SynologyDrive/FUNDIV - NFI - Europe/our-data/species/",CODE,"/CLIMAP/Recrut.Mortality.2020/")
    setwd(Dir)
    assign(paste0("df_",CODE),read.csv2(paste0(Dir,Modname,"/",Modname,"_Table.csv")))
  }
  
  test <- sapply(ls(pattern = "df_"), function(x) cbind(get(x)), simplify = FALSE)
  CondMod <- lapply(test,function(x) x[c(which(x$Estimate=="condmodel")+1):(which(x$Estimate=="zimodel")-1),])
  CondMod <- do.call("rbind",CondMod)
  CondMod$species <- gsub("df_(.*)\\.[0-9]*","\\1",rownames(CondMod)) ### modif 
  CondMod$model <- "cond"
  
  CondMod$X2 <- CondMod$X
  CondMod[CondMod$species%in%REV & grepl("Plotcat1",CondMod$X2),"X"] <- gsub("(Plotcat1)(.*)",paste0("Plotcat2","\\2"),CondMod[CondMod$species%in%REV & grepl("Plotcat1",CondMod$X2),"X2"])
  CondMod[CondMod$species%in%REV & grepl("Plotcat2",CondMod$X2),"X"] <- gsub("(Plotcat2)(.*)",paste0("Plotcat1","\\2"),CondMod[CondMod$species%in%REV & grepl("Plotcat2",CondMod$X2),"X2"])
  CondMod[CondMod$species%in%REV & grepl("Plotcat",CondMod$X),c("X","X2")]
  
  
  ZiMod <- lapply(test,function(x) x[c(which(x$Estimate=="zimodel")+1):(which(x$X=="AIC")-1),])
  ZiMod <- do.call("rbind",ZiMod)
  ZiMod$species <- gsub("df_(.*)\\.[0-9]*","\\1",rownames(ZiMod)) ### modif 
  ZiMod$model <- "zi"
  
  ZiMod$X2 <- ZiMod$X
  ZiMod[ZiMod$species%in%REV & grepl("Plotcat1",ZiMod$X2),"X"] <- gsub("(Plotcat1)(.*)",paste0("Plotcat2","\\2"),ZiMod[ZiMod$species%in%REV & grepl("Plotcat1",ZiMod$X2),"X2"])
  ZiMod[ZiMod$species%in%REV & grepl("Plotcat2",ZiMod$X2),"X"] <- gsub("(Plotcat2)(.*)",paste0("Plotcat1","\\2"),ZiMod[ZiMod$species%in%REV & grepl("Plotcat2",ZiMod$X2),"X2"])
  ZiMod[ZiMod$species%in%REV & grepl("Plotcat",ZiMod$X),c("X","X2")]
  
  rm(list=ls(pattern="df_"))
  DF <- rbind(CondMod[,c(6:7,1:5)],ZiMod[,c(6:7,1:5)])
  colnames(DF) <- c("species","model","variable","estimates","SE","Z-value","signif")
  
  DF$variable <- as.character(DF$variable)
  DF$estimates <- as.numeric(as.character(DF$estimates))
  DF$SE <- as.numeric(as.character(DF$SE))
  DF$`Z-value` <- as.numeric(as.character(DF$`Z-value`))
  DF$signif <- as.numeric(as.character(DF$signif))
  
  saveRDS(DF, file = paste0("~/SynologyDrive/FUNDIV - NFI - Europe/our-data/species/Recrut.All/Model_All_",Modname,".rds"))
}

