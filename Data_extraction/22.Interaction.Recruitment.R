rm(list = ls())
gc()
library(parallel)
library(forecast)
library(glmmTMB)
library(bbmle)
library(ggplot2)
library(reshape2)

# Dir <- c("/home/achangenet/SynologyDrive/FUNDIV - NFI - Europe/")
Dir <- c("/Users/alexandrechangenet/SynologyDrive/FUNDIV - NFI - Europe/")
setwd(Dir)
source(paste0(Dir,"Myscripts/Fundiv.project/RecrutProject.function1.Saving.R"))
Allcode <- c("PINPINA","ABIALB","ACEPSE","ALNGLU","BETPEN","CASSAT","FAGSYL",
             "FRAEXC","PICABI","PINHAL","PINNIG","PINPIN","PINSYL","POPNIG",
             "POPTRE","QUEILE","QUEPET","QUEPYR","QUEROB","QUESUB")
Allseuil <- c(0.8,0.7,0.55,0.6,0.8,0.7,0.7,0.7,0.8,0.7,0.7,0.7,0.8,0.7,0.7,0.8,0.7,0.7,0.8,0.7)
AllseuilC <- c(0.6,0.6,0.5,0.5,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6)


####################################################
Lvl <- 30                                         ## Number of values for the interaction graphs
Modname <- "SCALEDRrate.yr.QualiF.TF.noY.nbinom1.V2" ## Name of the model 
nameVAR <- c(                                     ##
  "sqrtBAIj.plot.1.mean.J.I",                     ##
  "mean_spei12",                                  ##
  "sqrtBA.O.plot.1",                              ##
  "logBAj.plot.1",                                ##
  "logtreeNbrJ",                                  ##
  "logdbhJ.IMall.plot.mean",                      ##
  "logsp.mortality.plot.rate.yr")                 ##
####################################################  
for (i in c(1:13,16:length(Allcode))){
  CODE <- Allcode[i]
  seuilC <- AllseuilC[i]
  seuil <- Allseuil[i]
  # Dir = paste0("/home/achangenet/SynologyDrive/FUNDIV - NFI - Europe/our-data/species/",CODE,"/CLIMAP/Recrut.Mortality.2020/")
  Dir = paste0("/Users/alexandrechangenet/SynologyDrive/FUNDIV - NFI - Europe/our-data/species/",CODE,"/Recrut.Mortality.2020/")
  
  setwd(Dir)
  ## Charge the dfplot that is located in the species recrut file 
  dfplot <- try(readRDS(paste0("dfplotFinalV2",CODE,seuil,"R.M.rds")),silent = T) #Base de données plot
  
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
  ### For each single one except one give average values
  data <- get(Modname)$frame
  ## Fixer tous les paramètres à leur moyennes 
  ## Sauf on veut 1/3 Plocat chaque cat
  myDF = as.data.frame(matrix(ncol = length(colnames(data)), nrow = Lvl*3, NA)) # Df vierge
  colnames(myDF) <- colnames(data) # Bon noms de colonnes
  # Dir <- c("/home/achangenet/SynologyDrive/FUNDIV - NFI - Europe/our-data/species/Recrut.All/")
  Dir <- c("/Users/alexandrechangenet/SynologyDrive/FUNDIV - NFI - Europe/our-data/species/Recrut.All/")
  dir.create(paste0(Dir,Modname,"/")) # Creation of a new file for this model
  Dir <- paste0(Dir,Modname,"/")
  ## Variable que l'on veut regarder (fixer et faire varier)
  for (i in 1:length(nameVAR[nameVAR%in%colnames(data)])){
    name <- nameVAR[nameVAR%in%colnames(data)][i]
    myDF[((Lvl*3)*(i-1)+1):(Lvl*3*i),nameVAR[nameVAR%in%colnames(data)]] <- apply(data[,nameVAR[nameVAR%in%colnames(data)]],2,
                                                                                  function(x) as.numeric(rep(mean(x),Lvl*3)))
    myDF[((Lvl*3)*(i-1)+1):(Lvl*3*i),name] <- seq(min(data[,name]),max(data[,name]),length=Lvl)
  }
  if (length(unlist(ranef(get(Modname))))!=0){myDF[((Lvl*3)*(i-1)+1):(Lvl*3*i),"country"] <- NA} # country
  myDF[,"Plotcat"] <- as.factor(c(rep("0",Lvl),rep("1",Lvl),rep("2",Lvl))) # Here our data are ready to be used for bootstrapping 
  
  ## Predicted response 
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
  saveRDS(get("myDF"), paste0("/Users/alexandrechangenet/SynologyDrive/FUNDIV - NFI - Europe/our-data/species/Recrut.All/",Modname,"/myDF.SIM_",Lvl,"_", CODE,".rds"))
  
  yMAX.resp <- max(myDF[,"Predictedresp"])+1.96*max(myDF[,"SEresp"]) # Here is the maximum values that will be the limit for the plot
  yMAX.mu <- max(myDF[,"Predictedmu"])+1.96*max(myDF[,"SEmu"]) # Here is the maximum values that will be the limit for the plot
  yMAX.resp <- max(myDF[,"Predictedzprob"])+1.96*max(myDF[,"SEzprob"]) # Here is the maximum values that will be the limit for the plot
  
  
  ## GGPLOT 
  
  # for (i in 1:length(nameVAR[nameVAR%in%colnames(data)])){
  #   name <- nameVAR[nameVAR%in%colnames(data)][i]
  #   dir.create(paste0(Dir,nameVAR[nameVAR%in%colnames(data)][i],"/"))
  #   setwd(paste0(Dir,nameVAR[nameVAR%in%colnames(data)][i],"/"))
  #   
  #   ## Extrapolated value in lighter colors
  #   Extrapolated0 <- seq(min(data[,name]),max(data[,name]),length=Lvl) >= min(data[data[,"Plotcat"]==0,name]) & seq(min(data[,name]),max(data[,name]),length=Lvl) <= max(data[data[,"Plotcat"]==0,name])
  #   Extrapolated1 <- seq(min(data[,name]),max(data[,name]),length=Lvl) >= min(data[data[,"Plotcat"]==1,name]) & seq(min(data[,name]),max(data[,name]),length=Lvl) <= max(data[data[,"Plotcat"]==1,name])
  #   Extrapolated2 <- seq(min(data[,name]),max(data[,name]),length=Lvl) >= min(data[data[,"Plotcat"]==2,name]) & seq(min(data[,name]),max(data[,name]),length=Lvl) <= max(data[data[,"Plotcat"]==2,name])
  #   Extrapolated <- as.factor(c(Extrapolated0,Extrapolated1,Extrapolated2))
  #   
  #   ## Colors and labels according the the reverse species or not (Trailing edge and leading edge reversed in some species)
  #   if (CODE%in%c("PINPIN","PINNIG","QUEPET","QUEROB","PINPINA","POPTRE","ALNGLU","QUEPYR","QUESUB")){ 
  #     Mycol <- c("black", "blue", "red") # Pas normale
  #     Mylabels <- c("Core","Leading Edge","Trailing Edge") # pas normale
  #   }else {
  #     Mylabels <- c("Core","Trailing Edge","Leading Edge") #Normale 
  #     Mycol <- c("black", "red", "blue")}
  #   
  #   
  #   p1<-ggplot(data=myDF[((Lvl*3)*(i-1)+1):(Lvl*3*i),], # All
  #              aes(x=get(name), y=Predictedresp, col=Plotcat,shape=Plotcat,group=Plotcat,alpha=Extrapolated,
  #                  ymin=Predictedresp-SEresp, ymax=Predictedresp+SEresp))+
  #     geom_line(size=1.5)+ #no size
  #     geom_point(size=2,stroke=2)+
  #     geom_linerange(size=1.5)+ # no size
  #     scale_alpha_manual(values = c(0.3,1),guide=F)+
  #     scale_color_manual("",values = Mycol,labels=Mylabels)+
  #     scale_shape_manual("", values=c(1,2,3),labels=Mylabels)+
  #     labs(title=paste0(CODE),y=paste0("Predicted response"), x=paste0(name))+
  #     theme_light(base_size = 15)+
  #     theme(text = element_text(face="bold"),legend.position = c(0.5,0.9),legend.direction ="horizontal",
  #           axis.text.x = element_text(size=13,color="black"),axis.text.y = element_text(size=13,color="black"),
  #           legend.background=element_rect(fill="white",colour="black",size=0.2),
  #           panel.border = element_rect(colour = "black", fill=NA, size=0.8),
  #           axis.line = element_line(colour="black"),
  #           plot.title = element_text(size=18,hjust = 0.5),
  #           plot.caption = element_text(face="bold.italic"))
  #   assign(paste0("presp_",CODE,"_",name),p1,envir = .GlobalEnv) # added to save the figure object because it is really long
  #   saveRDS(p1,file = paste0("presp_",CODE,"_",name,".rds"))
  #   ggsave(filename = paste0(CODE,"RESPONSE_",name,"_",Modname,".png"),plot = p1, width = 6, height = 6, dpi=300) # To modify 
  #   
  #   
  #   p2<-ggplot(data=myDF[((Lvl*3)*(i-1)+1):(Lvl*3*i),], # All
  #              aes(x=get(name), y=Predictedmu, col=Plotcat,shape=Plotcat,group=Plotcat,alpha=Extrapolated,
  #                  ymin=Predictedmu-SEmu, ymax=Predictedmu+SEmu))+
  #     geom_line(size=1.5)+ #no size
  #     geom_point(size=2,stroke=2)+
  #     geom_linerange(size=1.5)+ # no size
  #     scale_alpha_manual(values = c(0.3,1),guide=F)+
  #     scale_color_manual("",values = Mycol,labels=Mylabels)+
  #     scale_shape_manual("", values=c(1,2,3),labels=Mylabels)+
  #     labs(title=paste0(CODE),y=paste0("Predicted Mu"), x=paste0(name))+
  #     theme_light(base_size = 15)+
  #     theme(text = element_text(face="bold"),legend.position = c(0.5,0.9),legend.direction ="horizontal",
  #           axis.text.x = element_text(size=13,color="black"),axis.text.y = element_text(size=13,color="black"),
  #           legend.background=element_rect(fill="white",colour="black",size=0.2),
  #           panel.border = element_rect(colour = "black", fill=NA, size=0.8),
  #           axis.line = element_line(colour="black"),
  #           plot.title = element_text(size=18,hjust = 0.5),
  #           plot.caption = element_text(face="bold.italic"))
  #   assign(paste0("pmu_",CODE,"_",name),p2,envir = .GlobalEnv) # added to save the figure object because it is really long
  #   saveRDS(p2,file = paste0("pmu_",CODE,"_",name,".rds"))
  #   ggsave(filename = paste0(CODE,"MU_",name,"_",Modname,".png"),plot = p2, width = 6, height = 6, dpi=300) # To modify 
  #   
  #   p3<-ggplot(data=myDF[((Lvl*3)*(i-1)+1):(Lvl*3*i),], # All
  #              aes(x=get(name), y=1-Predictedzprob, col=Plotcat,shape=Plotcat,group=Plotcat,alpha=Extrapolated,
  #                  ymin=1-Predictedzprob-SEzprob, ymax=1-Predictedzprob+SEzprob))+
  #     geom_line(size=1.5)+ #no size
  #     geom_point(size=2,stroke=2)+
  #     geom_linerange(size=1.5)+ # no size
  #     scale_alpha_manual(values = c(0.3,1),guide=F)+
  #     scale_color_manual("",values = Mycol,labels=Mylabels)+
  #     scale_shape_manual("", values=c(1,2,3),labels=Mylabels)+
  #     labs(title=paste0(CODE),y=paste0("1-Predicted zprob"), x=paste0(name))+
  #     theme_light(base_size = 15)+
  #     theme(text = element_text(face="bold"),legend.position = c(0.5,0.9),legend.direction ="horizontal",
  #           axis.text.x = element_text(size=13,color="black"),axis.text.y = element_text(size=13,color="black"),
  #           legend.background=element_rect(fill="white",colour="black",size=0.2),
  #           panel.border = element_rect(colour = "black", fill=NA, size=0.8),
  #           axis.line = element_line(colour="black"),
  #           plot.title = element_text(size=18,hjust = 0.5),
  #           plot.caption = element_text(face="bold.italic"))
  #   assign(paste0("pzprob_",CODE,"_",name),p3,envir = .GlobalEnv) # added to save the figure object because it is really long
  #   saveRDS(p3,file = paste0("pzprob_",CODE,"_",name,".rds"))
  #   ggsave(filename = paste0(CODE,"1-ZPROB_",name,"_",Modname,".png"),plot = p3, width = 6, height = 6, dpi=300) # To modify 
  # }
}
