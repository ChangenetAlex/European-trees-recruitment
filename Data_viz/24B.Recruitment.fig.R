# Alex on the 24/03/2020
# Plot any interaction for all species at the same time with ggplot

SavingInfo = "interaction( 
name3 = 'the variable with three levels',
nameSeq= 'the variable to vary',
Pred='one response or one part of the model',
Modname='SCALEDRrate.yr.QualiF.TF.noY.nbinom1.V2,' for instance
save=False or True \n
Example : 
interaction(name3 = 'Plotcat',nameSeq = 'mean_spei12',Pred='Predictedmu',Modname='SCALEDRrate.yr.QualiF.TF.noY.nbinom1.V2',save=F)\n"
bannerBreak = "\n*********************************************************************\n"
cat(paste0(bannerBreak,SavingInfo,bannerBreak,"\n"))


rm(list = ls())
gc()
library(parallel)
library(forecast)
library(glmmTMB)
library(bbmle)
library(ggplot2)
library(reshape2)
library(stringr)

# Dir <- c("/home/achangenet/SynologyDrive/FUNDIV - NFI - Europe/")
Dir <- c("~/SynologyDrive/FUNDIV - NFI - Europe/")
setwd(Dir)
#source(paste0(Dir,"Myscripts/Fundiv.project/RecrutProject.function1.Saving.R"))
Allcode <- c("PINPINA","ABIALB","ACEPSE","ALNGLU","BETPEN","CASSAT","FAGSYL",
             "FRAEXC","PICABI","PINHAL","PINNIG","PINPIN","PINSYL","POPNIG",
             "POPTRE","QUEILE","QUEPET","QUEPYR","QUEROB","QUESUB")
Allseuil <- c(0.8,0.7,0.55,0.6,0.8,0.7,0.7,0.7,0.8,0.7,0.7,0.7,0.8,0.7,0.7,0.8,0.7,0.7,0.8,0.7)
AllseuilC <- c(0.6,0.6,0.5,0.5,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6)


######################################################
Lvl <- 30                                           ## Number of values for the interaction graphs
nameVAR <- c(                                       ##
  "sqrtBAIj.plot.1.mean.J.I",                       ##
  "mean_spei12",                                    ##
  "sqrtBA.O.plot.1",                                ##
  "logBAj.plot.1",                                  ##
  "logtreeNbrJ",                                    ##  All variables
  "logdbhJ.IMall.plot.mean",                        ##
  "logsp.mortality.plot.rate.yr",                   ##
  "Plotcat")                                        ##
                                                    ##
nameReal <- c(                                      ##  Real names for the graphs 
  "Species mean\ngrowth rate",                       ##
  "Relative drought index",                         ##
  "Interspecific \ncompetition",                      ##
  "Intraspecific \ncompetition",                      ##
  "Species tree density (ha)",                      ## 
  "Species average age\n(mean dbh)",                 ##
  "Annual \nmortality rate",                            ##
  "Climatic \nmarginality")                           ##
                                                    ##
A <- expand.grid(nameVAR,nameVAR)                   ## All combination of variable 
A <- A[,c(2,1)]                                     ##
A <- A[-c(1,8,10,16,19,24,28,32,                    ##
          37,40,46,48,55,56,64),]                   ##
COMB <- paste0(A[,1],A[,2],collapse=NULL)           ##
                                                    ##
PredVAR <- c(                                       ## Predictions to be on what part of the model ? 
  "Predictedresp",                                  ##
  "Predictedzprob",                                 ##
  "Predictedmu")                                    ##
                                                    ##
Predname <- c(                                      ##
  "Predicted number of saplings given the occurence probability (mu x p)",   # Real name on the graph
  "Occurence probability of recruitment in the census interval (p)",
  "Predicted number of saplings in the census interval (mu)")
                                                    ##
######################################################  

name3 = "sqrtBA.O.plot.1"
nameSeq <- "logdbhJ.IMall.plot.mean"
Pred="Predictedzprob"
Modname="SCALEDRrate.yr.QualiF.TF.noY.nbinom1.V3"
Modname <- "rW.yO.sumO.mRY1" ## Name of the model 


Interaction <- function(name3 = "Plotcat",nameSeq,Pred="Predictedmu",Modname="SCALEDRrate.yr.QualiF.TF.noY.nbinom1.V3",save=F){
  
  ## 1 ## Load all species data and keep the lines corresponding to our variables
  ## and load only species for which the interaction we want is significant
  rm(list=ls(pattern = "myDF"))
  tableresults <- readRDS(paste0("~/SynologyDrive/FUNDIV - NFI - Europe/our-data/species/Recrut.All/Model_All_",Modname,".rds"))
  if (Pred=="Predictedmu") (tableresults <- tableresults[tableresults$model=="cond",])
  if (Pred=="Predictedzprob") (tableresults <- tableresults[tableresults$model=="zi",])
  tableresults$variable <- gsub("Plotcat1","Plotcat",tableresults$variable)
  tableresults$variable <- gsub("Plotcat2","Plotcat",tableresults$variable)
  tableresults <- unique(tableresults[tableresults$variable%in%c(paste0(name3,":",nameSeq),paste0(nameSeq,":",name3))
                               & tableresults$signif<=0.05,
                               "species"]) # Keep only significant values
  
  ## Tableresults is the vector containing the species for which there is a signif effect. 
  
  ## Load the data values for only these species:(data used for predictions and also real values !!!)
  myDF <- list()
  mod_ALL <- list()
  for (j in which(Allcode%in%tableresults)){
  #for (j in c(1:13,16:length(Allcode))){
    CODE <- Allcode[j]
    seuilC <- AllseuilC[j]
    seuil <- Allseuil[j]
    Dir = paste0("~/SynologyDrive/FUNDIV - NFI - Europe/our-data/species/",CODE,"/CLIMAP/Recrut.Mortality.2020/")
    #Dir = paste0("~/SynologyDrive/FUNDIV - NFI - Europe/our-data/species/",CODE,"/Recrut.Mortality.2020/")
    setwd(Dir)
    myDF[[CODE]] <- readRDS(paste0("~/SynologyDrive/FUNDIV - NFI - Europe/our-data/species/Recrut.All/",Modname,"/myDF.SIM_All_30_", CODE,".rds"))[(1+(which(paste0(name3,nameSeq)==COMB)-1)*Lvl*3):(which(paste0(name3,nameSeq)==COMB)*Lvl*3),c(nameSeq,name3,"Predictedresp","SEresp","Predictedmu","SEmu","Predictedzprob","SEzprob")]
    # Real value"s 
    Dir = paste0("~/SynologyDrive/FUNDIV - NFI - Europe/our-data/species/",CODE,"/CLIMAP/Recrut.Mortality.2020/",Modname,"/")
    setwd(Dir)
    # Charge the model and the data
    mod_ALL[[CODE]] <- get(load(file = paste0(Modname,".rda")))
    rm(x)
    }
  
  
  ## 2 ## Put together all species with the right names and Convert the lists into dfs
  
  mod_ALL <- lapply(mod_ALL, function(x) {
    if (length(which(colnames(x[[5]])=="country")!=0)){
      x[[5]] <- x[[5]][,-which(colnames(x[[5]])=="country")]
    }
    ;x[[5]]})
  mod_ALL <- mapply(function(x,y) {
    x[,"species"] <- y
    return(x)
  },
  x=mod_ALL,y=names(mod_ALL),SIMPLIFY = F)


  myDF <- mapply(function(x,y) {
    x[,"species"] <- y
    return(x)
  },
  x=myDF,y=names(myDF),SIMPLIFY = F)
  myDF <- do.call(rbind,myDF)
  
  ## 2bis ## If the species is one of the following we need to exchange the values for plotcat 1 and plotcat 2
  ## If plotcat is not the factor variable, we need to replace the values of the factor variable by a real factor
  ## Chose label and colors according to this factor (plotcat or other variables)
  ## Finally we need to find out for each category of the factor (0,1,2 or low high, mean) wether the values for predictions are observed or not in reality. (Extrapolation)
  
  REV <- c("PINPIN","PINNIG","QUEPET","QUEROB","PINPINA","POPTRE","ALNGLU","QUEPYR","QUESUB")
  if (name3=="Plotcat"){
    myDF$Plotcat2 <- myDF$Plotcat
    myDF[myDF$species%in%REV & myDF$Plotcat2==1,"Plotcat"] <- 2
    myDF[myDF$species%in%REV & myDF$Plotcat2==2,"Plotcat"] <- 1
    myDF$Plotcat <- as.factor(myDF$Plotcat)
    
    Mylabels <- c("0"="Core",                       
                  "1"="Trailing Edge",               
                  "2"="Leading Edge")
    Mycol <- c("0"="black",
               "1"="red",
               "2"="blue")
    
    myDF <- split(myDF,myDF$species)
    myDF <- mapply(function(x,y) {
      x[x[,name3]=="0","extra"] <- seq(min(y[,nameSeq]),max(y[,nameSeq]),length=Lvl) >= min(y[y[,name3]=="0",nameSeq]) & seq(min(y[,nameSeq]),max(y[,nameSeq]),length=Lvl) <= max(y[y[,name3]=="0",nameSeq])
      x[x[,name3]=="1","extra"] <- seq(min(y[,nameSeq]),max(y[,nameSeq]),length=Lvl) >= min(y[y[,name3]=="1",nameSeq]) & seq(min(y[,nameSeq]),max(y[,nameSeq]),length=Lvl) <= max(y[y[,name3]=="1",nameSeq])
      x[x[,name3]=="2","extra"] <- seq(min(y[,nameSeq]),max(y[,nameSeq]),length=Lvl) >= min(y[y[,name3]=="2",nameSeq]) & seq(min(y[,nameSeq]),max(y[,nameSeq]),length=Lvl) <= max(y[y[,name3]=="2",nameSeq])
      return(x)
    },
    x=myDF,y=mod_ALL,SIMPLIFY = F)
    
  }else{
    colnames(myDF)[which(colnames(myDF)==name3)] <- paste0(name3,"2")
    myDF[,name3] <- as.factor(c(rep("mean",Lvl),rep("high",Lvl),rep("low",Lvl)))
    Mylabels <- c("mean"="Average level",
                  "high"="High level",   
                  "low"="Low level")    
    Mycol <- c("mean"="black",
               "high"="red",
               "low"="blue")
    myDF <- split(myDF,myDF$species)
    myDF <- mapply(function(x,y) {
      x[x[,name3]=="mean","extra"] <- seq(min(y[,nameSeq]),max(y[,nameSeq]),length=Lvl) >= min(y[y[,name3]>quantile(y[,name3],round(ecdf(y[,name3])(mean(y[,name3])),3)-0.05) & y[,name3]<quantile(y[,name3],round(ecdf(y[,name3])(mean(y[,name3])),3)+0.05),nameSeq]) & seq(min(y[,nameSeq]),max(y[,nameSeq]),length=Lvl) <= max(y[y[,name3]>quantile(y[,name3],round(ecdf(y[,name3])(mean(y[,name3])),3)-0.05) & y[,name3]<quantile(y[,name3],round(ecdf(y[,name3])(mean(y[,name3])),3)+0.05),nameSeq])
      x[x[,name3]=="high","extra"] <- seq(min(y[,nameSeq]),max(y[,nameSeq]),length=Lvl) >= min(y[y[,name3]>=quantile(y[,name3],0.95),nameSeq]) & seq(min(y[,nameSeq]),max(y[,nameSeq]),length=Lvl) <= max(y[y[,name3]>=quantile(y[,name3],0.95),nameSeq])
      x[x[,name3]=="low","extra"] <- seq(min(y[,nameSeq]),max(y[,nameSeq]),length=Lvl) >= min(y[y[,name3]<=quantile(y[,name3],0.05),nameSeq]) & seq(min(y[,nameSeq]),max(y[,nameSeq]),length=Lvl) <= max(y[y[,name3]<=quantile(y[,name3],0.05),nameSeq])
      return(x)
    },
    x=myDF,y=mod_ALL,SIMPLIFY = F)}
  
  
  
#mod_ALL[,name3]>quantile(mod_ALL[,name3],0.005) & mod_ALL[,name3]<quantile(mod_ALL[,name3],0.995)
#round(ecdf(mod_ALL[,name3])(mean(mod_ALL[,name3])),3)-0.025
#round(ecdf(mod_ALL[,name3])(mean(mod_ALL[,name3])),3)+0.025


  ## Finally put back in a df (not a list)
  
myDF <- do.call(rbind,myDF)
mod_ALL <- do.call(rbind,mod_ALL)
  

  ## 3 ### Plot
  try(rm(p1),silent = T)
  if (Pred==PredVAR[2]){p1<-ggplot(get(paste0("myDF")),aes(x=get(nameSeq), y=1-get(Pred), color=get(name3),shape=get(name3),group=get(name3),alpha=extra))
  }else {p1<-ggplot(get(paste0("myDF")),aes(x=get(nameSeq), y=get(Pred), color=get(name3),shape=get(name3),group=get(name3),alpha=extra))}
  p1 <- p1+geom_line(size=1)+
    geom_point(size=1.5,stroke=1.5)+
    facet_wrap(vars(species),scales = "free")+
    #geom_linerange(size=1.5)+ # no size
    scale_alpha_manual(name=paste0("",nameReal[nameVAR==name3]),values = c(0.3,1),guide=F)+
    scale_color_manual(name=paste0("",nameReal[nameVAR==name3]),values = Mycol,labels=Mylabels)+
    scale_shape_manual(name=paste0("",nameReal[nameVAR==name3]), values=rep(c(1,2,3)),labels=Mylabels)+
    labs(y=Predname[PredVAR==Pred], x=nameReal[nameVAR==nameSeq],title=paste0(Predname[PredVAR==Pred]," VS ",nameReal[nameVAR==nameSeq]))+
    theme_light(base_size = 15)+
    theme(text = element_text(face="bold"),legend.position = c(0.8,0.87),legend.direction ="vertical",
          axis.text.x = element_text(size=13,color="black"),axis.text.y = element_text(size=13,color="black"),
          legend.text=element_text(size=15),
          legend.key.size = unit(2.2,"line"),
          legend.background=element_rect(fill="white",colour="black",size=0.2),
          panel.border = element_rect(colour = "black", fill=NA, size=0.8),
          axis.line = element_line(colour="black"),
          plot.title = element_text(size=18,hjust = 0.5),
          plot.caption = element_text(face="bold.italic"))
  
  print(p1)
  if (save==T){ggsave(filename = paste0("~/SynologyDrive/FUNDIV - NFI - Europe/our-data/species/Recrut.All/",Modname,"/",Pred,"_VS_",nameSeq,"_for_3_level_of_",name3,"_",Modname,".png"),plot = p1, scale=0.8, width = 20,height = 16, dpi=400)} # To modify 
}


Interaction(name3 = nameVAR[6],nameSeq=nameVAR[6],Pred=PredVAR[2],Modname="SCALEDRrate.yr.QualiF.TF.noY.nbinom1.V3",save=F)

  


