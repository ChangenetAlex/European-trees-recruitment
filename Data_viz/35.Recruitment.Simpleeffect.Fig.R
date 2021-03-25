# Alex on the 24/03/2020
# Plot any effect for all species at the same time with ggplot

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
nameReal <- c(                                      ##  Real names for the graphs ### To change if we want the plot with transformed variables !
  "Species mean\ngrowth rate (square rooted)",                       ##
  "Relative drought index",                         ##
  "Interspecific \ncompetition (square rooted)",                      ##
  "Intraspecific \ncompetition (log)",                      ##
  "Tree density (log)",                      ## 
  "Species mean DBH\n(log)",                 ##
  "Annual \nmortality rate (log)",                            ##
  "Climatic \nmarginality")                           ##

##
PredVAR <- c(                                       ## Predictions to be on what part of the model ? 
  "Predictedresp",                                  ##
  "Predictedzprob",                                 ##
  "Predictedmu",
  "Predictedmulink",
  "Predictedzlink")                                    ##
##
Predname <- c(                                      ##
  "Predicted number of saplings given the occurence probability (mu x p)",   # Real name on the graph
  "Recruitment occurrence (p)",
  "Predicted recruitment count (mu)",
  "Predicted recruitment count (mu) on the link scale",
  "Recruitment occurrence (p) on the link scale")
##
######################################################  

VARIABLE = "sqrtBA.O.plot.1"
Pred="Predictedmu"
Modname <- "rW.yO.sumO.mRY1" ## Name of the model 
Lvl <- 30

Interaction <- function(VARIABLE = "Plotcat",Pred="Predictedmu",FACET=T,save=T){
myDF <- readRDS(paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/paper2/myDF.SIMPLE.EFFECT.SIM_",Lvl,".rds"))
Species18 <- names(myDF) # All species name

#############################################################
## Here need to detransform (here we just delog or desqraed)#
## Need to load the dfplot with orginal values              #
#############################################################

# dfplot <- mapply(function(x,y){
#   Dir = paste0("~/SynologyDrive/FUNDIV - NFI - Europe/our-data/species/",x,"/CLIMAP/Recrut.Mortality.2020/")
#   readRDS(paste0(Dir, "/dfplotV2", x, y, "R.M.rds"))
# },x=Allcode,y=Allseuil,SIMPLIFY = F)
# dfplot <- dfplot[which(names(dfplot)%in%names(myDF))] #keep inly species for which we fitted a model
# 
# # Then need to descaled the target (based on the log variable )
# myDF <- mapply(function(x,y){
#   if (VARIABLE%in%c("sqrtBAIj.plot.1.mean.J.I","sqrtBA.O.plot.1")){
#     x[,VARIABLE] <- x[,VARIABLE]*sqrt(mean(y[,VARIABLE]^2,na.rm=T)) # scaled on the loged !!
#     x[,VARIABLE] <- x[,VARIABLE]^2
#   }else if (VARIABLE%in%c("logBAj.plot.1",
#                            "logtreeNbrJ",
#                            "logdbhJ.IMall.plot.mean",
#                            "logsp.mortality.plot.rate.yr")){
#     x[,VARIABLE] <- x[,VARIABLE]*sqrt(mean(y[,VARIABLE]^2,na.rm=T)) # scaled on the loged !!
#     x[,VARIABLE] <- exp(x[,VARIABLE])-1
#   }
#   return(x)
# },
#x=myDF,y=dfplot,SIMPLIFY = F)
#############################################################


Myshape18 <- c(0:17)
Mycol18 <- c(          # All colors 
  "#F8766D",
  "#E88526",
  "green",
  "red3",
  "gray70",
  "#5EB300",
  "gray50",
  "#00BF74",
  "#00C19F",
  "black",
  "#00B9E3",
  "gray30",
  "lightskyblue",
  "blue3",
  "#DB72FB",
  "gray15",
  "deeppink3",
  "orange")
  
  ## 1 ## Load all species data and keep the lines corresponding to our variables
  ## and load only species for which the interaction we want is significant
  tableresults <- readRDS(paste0("~/SynologyDrive/FUNDIV - NFI - Europe/our-data/species/Recrut.All/Model_All_",Modname,".rds"))
  if (Pred=="Predictedmu") (tableresults <- tableresults[tableresults$model=="cond",])
  if (Pred=="Predictedzprob") (tableresults <- tableresults[tableresults$model=="zi",])
  tableresults$variable <- gsub("Plotcat1","Plotcat",tableresults$variable)
  tableresults$variable <- gsub("Plotcat2","Plotcat",tableresults$variable)
  tableresults <- unique(tableresults[tableresults$variable%in%c(VARIABLE) & tableresults$signif<=0.05,
                                      "species"]) # Keep only significant values
  
  ## Tableresults is the vector containing the species for which there is a signif effect. 
  
  ## Load the data values for only these species:(data used for predictions and also real values !!!)

  i <- which(VARIABLE==nameVAR)
  monDF <- myDF[which(names(myDF)%in%tableresults)] #keep inly species for which this effect is significant 
  SpeciesSignif <- names(monDF) #Here take only the significant species names
  MycolSignif <- Mycol18[which(Species18%in%SpeciesSignif)] # Here association of the signif species names with the associated colors (according to number)
  MyshapeSignif <- Myshape18[which(Species18%in%SpeciesSignif)] # Here association of the signif species names with the associated colors (according to number)
  
  
  
  monDF <- lapply(monDF,function(x) x[(Lvl*(i-1)+1):(Lvl*(i)),
                                     c(VARIABLE,"Predictedresp","SEresp",
                                       "Predictedmu","SEmu",
                                       "Predictedzprob","SEzprob",
                                       "Predictedmulink","SEmulink",
                                       "Predictedzlink","SEzlink",
                                       "species")]) ## keep the lines we want
  
  
  myDF <- do.call(rbind,monDF)
  
  
  
  
  
  ## 2bis ## If the species is one of the following we need to exchange the values for plotcat 1 and plotcat 2
  ## If plotcat is not the factor variable, we need to replace the values of the factor variable by a real factor
  ## Chose label and colors according to this factor (plotcat or other variables)
  ## Finally we need to find out for each category of the factor (0,1,2 or low high, mean) wether the values for predictions are observed or not in reality. (Extrapolation)
  
  REV <- c("PINPIN","PINNIG","QUEPET","QUEROB","PINPINA","POPTRE","ALNGLU","QUEPYR","QUESUB")
  if (VARIABLE=="Plotcat"){
    myDF$Plotcat2 <- myDF$Plotcat
    myDF[myDF$species%in%REV & myDF$Plotcat2==1,"Plotcat"] <- 2
    myDF[myDF$species%in%REV & myDF$Plotcat2==2,"Plotcat"] <- 1
    myDF$Plotcat <- as.factor(myDF$Plotcat)}
  
  
  ## 3 ### Plot
  try(rm(p1),silent = T)
  if (Pred==PredVAR[2]|Pred==PredVAR[5]){p1<-ggplot(get(paste0("myDF")),aes(x=get(VARIABLE), y=1-get(Pred),color=species,group=species,shape=species))
  #+geom_ribbon(aes(ymin = 1-get(Pred)-SEzprob , ymax = 1-get(Pred)+SEzprob),alpha=0.1)
  }else {p1<-ggplot(get(paste0("myDF")),aes(x=get(VARIABLE), y=get(Pred),color=species,group=species,shape=species))}
    #+geom_ribbon(aes(ymin = get(Pred)-SEmu, ymax = get(Pred)+SEmu,color=species,group=species))}
  p1 <- p1+geom_line(size=1)+geom_point(size=3,stroke=1)
  if(FACET==T){p1 <- p1+facet_wrap(vars(species),scales = "free")+scale_color_manual(values = rep("black",length(tableresults)))
  }else{p1 <- p1}
  #geom_linerange(size=1.5)+ # no size
  #scale_alpha_manual(name=paste0("",nameReal[nameVAR==name3]),values = c(0.3,1),guide=F)+
  #scale_shape_manual(name=paste0("",nameReal[nameVAR==name3]), values=rep(c(1,2,3)),labels=Mylabels)+
  p1 <- p1+labs(y=Predname[PredVAR==Pred], x=nameReal[nameVAR==VARIABLE])+
    # title=paste0(Predname[PredVAR==Pred]," VS ",nameReal[nameVAR==VARIABLE])
    scale_color_manual(name="",values = MycolSignif,labels=SpeciesSignif)+
    scale_shape_manual(name="",values = MyshapeSignif,labels=SpeciesSignif)+
    #ylim(0,max(myDF[,Pred])*1.1)+
    #scale_y_continuous(expand = c(0.1,0)) + # for the loop (most of the figure)
    scale_y_continuous(expand = c(0.05,0)) + # after the loop
    scale_x_continuous(expand = c(0.03,0)) +
    theme_light(base_size = 15)+
    theme(text = element_text(face="bold"),legend.position = c(.01, .99),
          legend.justification = c("left", "top"),legend.direction ="horizontal",
          axis.text.x = element_text(size=13,color="black"),axis.text.y = element_text(size=13,color="black"),
          legend.text=element_text(size=12),
          legend.key.size = unit(1.5,"line"),
          legend.margin = margin(1,1,1,1),
          legend.box.margin = margin(1,1,1,1),
          legend.title = element_blank(),
          legend.background=element_rect(fill="white",colour="black",size=0.2),
          panel.border = element_rect(colour = "black", fill=NA, size=0.8),
          axis.line = element_line(colour="black"),
          plot.title = element_text(size=18,hjust = 0.5),
          plot.caption = element_text(face="bold.italic"))
  
  print(p1)
  if (save==T){
    save_plot(filename = paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/paper2/",Pred,"_",VARIABLE,".png"),
                         plot = p1,base_width = 12, base_height = 7, dpi = 400 ,units = "in",1)
    saveRDS(p1,file = paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/paper2/",Pred,"_",VARIABLE,".rds"))}
}

# This function save all plots one by one. 
# Then need to load them two by two to assemble each effect with each otehrs. 


#Interaction(VARIABLE = nameVAR[5],Pred=PredVAR[3],FACET=F,save=T)

### Location of the legend not good in few: 

for (i in 1:7){
  for (j in 2:3){
    Interaction(VARIABLE = nameVAR[i],Pred=PredVAR[j],FACET=F,save=T)
  }
}  


# ## these line to change the size of axis and adapt legend
# Interaction(VARIABLE = nameVAR[7],Pred=PredVAR[2],FACET=F,save=T)
# Interaction(VARIABLE = nameVAR[6],Pred=PredVAR[2],FACET=F,save=T)
# Interaction(VARIABLE = nameVAR[5],Pred=PredVAR[2],FACET=F,save=T)
# Interaction(VARIABLE = nameVAR[4],Pred=PredVAR[3],FACET=F,save=T)
# Interaction(VARIABLE = nameVAR[4],Pred=PredVAR[2],FACET=F,save=T)
# Interaction(VARIABLE = nameVAR[3],Pred=PredVAR[2],FACET=F,save=T)


library(patchwork)

MyVAR <- nameVAR[1]
A <- list.files(paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/paper2/"),pattern = paste0(MyVAR,"*.rds$"),full.names = T)
# Load only the first half 
mapply(function(x,y){
  assign(paste0("p.",y),readRDS(x),envir = .GlobalEnv)
},x=A,y=PredVAR[c(2:3)]) ## be carreful with the order 
pall <- p.Predictedmu+theme(legend.position = c(0.99,0.85),legend.justification = c("right"))+guides(colour = guide_legend(ncol = 1))+
  p.Predictedzprob+theme(legend.position = c(0.01,0.99),legend.justification = c("left","top"))+guides(colour = guide_legend(ncol = 2))+
  plot_layout(ncol=2)
pall
save_plot(filename = paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/paper2/",MyVAR,"Effect_bothmodel.png"),
          plot = pall,base_width = 8, base_height = 8, dpi = 300 ,units = "in",nrow=1,ncol=2)


MyVAR <- nameVAR[2]
A <- list.files(paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/paper2/"),pattern = paste0(MyVAR,"*.rds$"),full.names = T)
# Load only the first half 
mapply(function(x,y){
  assign(paste0("p.",y),readRDS(x),envir = .GlobalEnv)
},x=A,y=PredVAR[c(2:3)]) ## be carreful with the order 
pall <- p.Predictedmu+
  p.Predictedzprob+
  plot_layout(ncol=2)
pall
save_plot(filename = paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/paper2/",MyVAR,"Effect_bothmodel.png"),
          plot = pall,base_width = 8, base_height = 8, dpi = 300 ,units = "in",nrow=1,ncol=2)



MyVAR <- nameVAR[3]
A <- list.files(paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/paper2/Figure.Simple.effect.pred/"),pattern = paste0(MyVAR,"*.rds$"),full.names = T)
# Load only the first half 
mapply(function(x,y){
  assign(paste0("p2.",y),readRDS(x),envir = .GlobalEnv)
},x=A,y=PredVAR[c(2:3)]) ## be carreful with the order 
pall <- p.Predictedmu+
  theme(legend.position = c(0.99,0.35),legend.justification = c("right"))+guides(colour = guide_legend(ncol = 3))+
  p.Predictedzprob+
  theme(legend.position = c(0.99,0.35),legend.justification = c("right"))+guides(colour = guide_legend(ncol = 3))+
  plot_layout(ncol=2)
pall
save_plot(filename = paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/paper2/",MyVAR,"Effect_bothmodel.pdf"),
          plot = pall,base_width = 8, base_height = 8, dpi = 300 ,units = "in",nrow=1,ncol=2)

MyVAR <- nameVAR[4]
A <- list.files(paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/paper2/Figure.Simple.effect.pred/"),pattern = paste0(MyVAR,"*.rds$"),full.names = T)
# Load only the first half 
mapply(function(x,y){
  assign(paste0("p1.",y),readRDS(x),envir = .GlobalEnv)
},x=A,y=PredVAR[c(2:3)]) ## be carreful with the order 
pall2 <- p.Predictedmu+theme(legend.position = c(.99, .99),legend.justification = c("right", "top"))+guides(colour = guide_legend(ncol = 2))+
  p.Predictedzprob+theme(legend.position = c(.99, .99),legend.justification = c("right", "top"))+guides(colour = guide_legend(ncol = 2))+
  plot_layout(ncol=2)
pall2
save_plot(filename = paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/paper2/",MyVAR,"Effect_bothmodel.pdf"),
          plot = pall2,base_width = 8, base_height = 8, dpi = 300 ,units = "in",nrow=1,ncol=2)


p <- 
  p1.Predictedmu+theme(legend.position = c(.99, .99),legend.justification = c("right", "top"))+guides(colour = guide_legend(ncol = 2))+
  ggtitle('a) Intraspecific competition effect on p')+theme(plot.title = element_text(size = 16))+
  
  p1.Predictedzprob+theme(legend.position = c(.99, .99),legend.justification = c("right", "top"))+guides(colour = guide_legend(ncol = 2))+
  ggtitle('b) Intraspecific competition effect on mu')+theme(plot.title = element_text(size = 16))+
  
  p2.Predictedmu+
  theme(legend.position = c(0.99,0.35),legend.justification = c("right"))+guides(colour = guide_legend(ncol = 3))+
  ggtitle('c) Interspecific competition effect on p')+theme(plot.title = element_text(size = 16))+
  
  p2.Predictedzprob+
  theme(legend.position = c(0.99,0.35),legend.justification = c("right"))+guides(colour = guide_legend(ncol = 3))+
  ggtitle('d) Interspecific competition effect on mu')+theme(plot.title = element_text(size = 16))+
  
  plot_layout(ncol=2,nrow=2)

save_plot(filename = paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/paper2/Competition_simple.pdf"),
          plot = p,base_width = 8, dpi = 300, base_height = 8,units = "in",nrow=2,ncol=2)
# save_plot(filename = paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/Paper1/Figure3/Ag.Prediction.Mort.bin.Part1.pdf"),
#           plot = pall,base_width = 5, base_height = 5, dpi = 300 ,units = "in",nrow=3,ncol=3)




MyVAR <- nameVAR[5]
A <- list.files(paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/paper2/"),pattern = paste0(MyVAR,"*.rds$"),full.names = T)
# Load only the first half 
mapply(function(x,y){
  assign(paste0("p.",y),readRDS(x),envir = .GlobalEnv)
},x=A,y=PredVAR[c(2:3)]) ## be carreful with the order 
pall <- p.Predictedmu+theme(legend.position = c(.99, .5),legend.justification = c("right", "center"))+guides(colour = guide_legend(ncol = 3))+
  p.Predictedzprob+theme(legend.position = c(.01, .99),legend.justification = c("left", "top"))+guides(colour = guide_legend(ncol = 3))+
  plot_layout(ncol=2)
pall
save_plot(filename = paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/paper2/",MyVAR,"Effect_bothmodel.png"),
          plot = pall,base_width = 8, base_height = 8, dpi = 300 ,units = "in",nrow=1,ncol=2)



MyVAR <- nameVAR[6]
A <- list.files(paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/paper2/"),pattern = paste0(MyVAR,"*.rds$"),full.names = T)
# Load only the first half 
mapply(function(x,y){
  assign(paste0("p.",y),readRDS(x),envir = .GlobalEnv)
},x=A,y=PredVAR[c(2:3)]) ## be carreful with the order 
pall <- p.Predictedmu+theme(legend.position = c(.99, .2),legend.justification = c("right", "center"))+guides(colour = guide_legend(ncol = 2))+
  p.Predictedzprob+theme(legend.position = c(.01, .99),legend.justification = c("left", "top"))+guides(colour = guide_legend(ncol = 3))+
  plot_layout(ncol=2)
pall
save_plot(filename = paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/paper2/",MyVAR,"Effect_bothmodel.png"),
          plot = pall,base_width = 8, base_height = 8, dpi = 300 ,units = "in",nrow=1,ncol=2)









