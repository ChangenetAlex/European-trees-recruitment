rm(list = ls())
gc()
require(rgdal)
require(adegenet)
require(ade4)
require(parallel)
require(fields)
library(rworldmap)
library(lattice)
require(spatial.tools)
library(maptools)
require(rworldxtra)
library(rgeos)
library(RStoolbox)
require(stringr)
library(data.table)
library(ggplot2) 
library(rangeBuilder)
library(knitr)
library(spaMM)
library(fields)
library(piecewiseSEM)
library(raster)
library(rasterVis)
library(cowplot)
library(sf)
library(rnaturalearth)
library(grid)
library("ggplotify")


### Two maps per species: Observed recruit continous
### Observed recruit quartile


QUANT.all <- data.frame(matrix(nrow=18,ncol=4,NA))
colnames(QUANT.all) <- c("Q1 fit","Q3 fit","Q1 obs","Q3 obs")
Modname <- "SCALEDRrate.yr.QualiF.TF.noY.nbinom1.V3" ## Name of the model 
Modname <- "rW.yO.sumO.mRY1"
i <- 1
Allcode <- c("ABIALB","ACEPSE","ALNGLU","BETPEN","CASSAT","FAGSYL",
             "FRAEXC","PICABI","PINHAL","PINNIG","PINPIN","PINPINA","PINSYL"
             ,"QUEILE","QUEPET","QUEPYR","QUEROB","QUESUB")
Allseuil <- c(0.7,0.55,0.6,0.8,0.7,0.7,0.7,0.8,0.7,0.7,0.7,0.8,0.8,0.8,0.7,0.7,0.8,0.7)
AllseuilC <- c(0.6,0.5,0.5,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6)
rownames(QUANT.all) <-Allcode


#ici besoin du même nom pour tous 

## Faire la même avec 4 couleurs 

#### 3 maps can be made: 
## Observed annual recruitment 
## Fitted annual recruitment 
## 

#for (i in c(8:12,14,16,18)){ #8
#for (i in c(4:7,13,15)){  #6
#for (i in c(1:3,17)){ #4

for (i in c(1:18)){ #4

  if (Allcode[i]%in%c("PINPIN","PINNIG","QUEPET","QUEROB","PINPINA","POPTRE","ALNGLU","QUEPYR","QUESUB")){
    Mycol <- c("0"="green","2"="red","1"="blue","10"="gray")
    Mylab <- c("0"="Core","2"="Trailing edge","1"="Leading edge","10"="Transition area")
  }else{
    Mycol <- c("0"="green","1"="red","2"="blue","10"="gray")
    Mylab <- c("0"="Core","1"="Trailing edge","2"="Leading edge","10"="Transition area")
    } 
  
  Dir = c("~/SynologyDrive/FUNDIV - NFI - Europe/our-data/") # Directory 
  test6 <- get(load(file=paste0(Dir,"species/",Allcode[i],"/CLIMAP/CLIMAP2019/Climap_",Allcode[i],"_",Allseuil[i],"_",AllseuilC[i],".RData")))
  Dir = paste0("~/SynologyDrive/FUNDIV - NFI - Europe/our-data/species/",Allcode[i],"/CLIMAP/Recrut.Mortality.2020/")
  setwd(Dir)
  # Charge the dfplot for simulation and residuals test dispersion etc... 
  dfplot <- readRDS(paste0("dfplotFinalV3",Allcode[i],Allseuil[i],"R.M.rds"))
  #Load the desired model of the desired species 
  Dir <- c(paste0(Dir,Modname,"/")) # Directory 
  setwd(Dir)
  x <- get(load(file = paste0(Modname,".rda")))
  A <- x[[5]]
  colnames(A)[which(colnames(A)=="offset(log(yearsbetweensurveys))")] <- "yearsbetweensurveys"
  colnames(A)[which(colnames(A)=="offset(log(sp.SUM.ha.weight2))")] <- "sp.SUM.ha.weight2"
  A[,"yearsbetweensurveys"] <- exp(A[,"yearsbetweensurveys"])
  A[,"sp.SUM.ha.weight2"] <- exp(A[,"sp.SUM.ha.weight2"])
  x$frame[,"fit"] <- predict(x,A,type=c("response"))
  ### Test #### 
  #we define the data as the part of the dataset on which the model has been fitted. (for every loop)
  x <- x$frame
  # Recup latitude and longitude 
  r1 <- data.frame(cbind(dfplot[which(rownames(dfplot)%in%rownames(x)),c("longitude","latitude")],x$fit,x$sp.recruitment.ha.weight.2,x$Plotcat))
  colnames(r1) <- c("X","Y","Z","Zbis","Plotcat")
  QUANT.all[which(rownames(QUANT.all)==Allcode[i]),1] <- signif(quantile(r1$Z,0.25,type=7),digits=2)
  QUANT.all[which(rownames(QUANT.all)==Allcode[i]),2] <- signif(quantile(r1$Z,0.75,type=7),digits=2)
  QUANT.all[which(rownames(QUANT.all)==Allcode[i]),3] <- signif(quantile(r1$Zbis,0.25,type=7),digits=2)
  QUANT.all[which(rownames(QUANT.all)==Allcode[i]),4] <- signif(quantile(r1$Zbis,0.75,type=7),digits=2)
  print(QUANT.all)
  r1$bin <- ifelse(r1$Z>quantile(r1$Z,0.75,type=7),"0",ifelse(r1$Z<quantile(r1$Z,0.25,type=7),"2","1")) ## Percentage of mortality as classes 
  r1$binbis <- ifelse(r1$Zbis>quantile(r1$Zbis,0.75,type=7),"0",ifelse(r1$Zbis<quantile(r1$Zbis,0.25,type=7),"2","1")) ## Percentage of mortality as classes 
  # print(summary(r1$Z))
  # print(summary(r1$Zbis))
  
  r1$Z <- log(1+r1$Z)
  r1$Zbis <- log(1+r1$Zbis)
  
  ## Top obtain a map that is square: 
  #extract the largest range between X and Y
  Yrange <- range(r1$Y)[2]-range(r1$Y)[1]
  Xrange <- range(r1$X)[2]-range(r1$X)[1]
  rangeDiff <- abs(Yrange-Xrange)/2
  
  ## add the difference of range divieded by two at bioth end of the smallest range
  if (Yrange>Xrange){
    MyYrange <- range(r1$Y)
    MyXrange <- c(min(r1$X)-rangeDiff,max(r1$X)+rangeDiff)
  }else{
    MyXrange <- range(r1$X)
    MyYrange <- c(min(r1$Y)-rangeDiff,max(r1$Y)+rangeDiff)}
  dfrange <- as.data.frame(cbind(MyXrange,MyYrange))
  ### First plot the inset map
  par(mfrow=c(1,1))
  worldmap <- getMap(resolution = "high")
  #europe <- worldmap[which(worldmap$REGION=="Europe"),]
  #europe <- worldmap[which(grepl("Europe", worldmap$GEO3) & as.character(worldmap$NAME) != "Russia"),] 
  europe <- worldmap[which(grepl("Europe", worldmap$GEO3)|grepl("North Africa", worldmap$GEO3)),] 
  europe <-spTransform(europe,"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0 ") # Convert to the right system
  #europe1 = crop(europe,c(-15, 45, 35, 70))
  europe1 = crop(europe,c(-10, 32, 36, 70))
  #europe = crop(europe,c(min(r1$X), max(r1$X), min(r1$Y), max(r1$Y))) ## cut the europe with the minimum rectangle to obtain the data 
  europe = crop(europe,c(MyXrange,MyYrange))
  # Here is the inset map in which we will zoom
  insetMap <- ggplot() + 
    labs(x = NULL, y = NULL)+
    geom_raster(data=test6[,c(1:2,29)],aes(x=x, y=y,fill=as.factor(groupes)),na.rm=T)+#fill="gray80",color="gray80",interpolate=T) +
    xlim(-10,32)+
    ylim(36,70)+
    coord_fixed(1.3) + 
    geom_polygon(data = europe1, aes(x=long, y = lat,group=group),fill=NA,col="black") + #gray 60 & gray 60
    #geom_point(data = r1, aes(x = X, y = Y,col=Z,fill=Z),shape=21,size=1)+
    scale_fill_manual(values = Mycol,name="Climatic areas: ",labels=Mylab)+
    #geom_rect(data = data.frame(),aes(xmin=min(r1$X), xmax=max(r1$X), ymin=min(r1$Y), ymax=max(r1$Y)),colour = "red", fill = NA,size=1)+ ## Rectangle around the area
    geom_rect(data=dfrange,
              aes(xmin=ifelse(MyXrange[1]<-10,min(r1$X),MyXrange[1]), xmax=ifelse(MyXrange[2]>32,max(r1$X),MyXrange[2]),
                  ymin=ifelse(MyYrange[1]<36,min(r1$Y),MyYrange[1]), ymax=ifelse(MyYrange[2]>70,max(r1$Y),MyYrange[2])),
              colour = "red",fill=NA,size=1)+ ## Square !!! 
    theme(text = element_text(face="bold"),
          legend.direction ="horizontal",legend.position = "none",
          legend.key = element_rect(fill = "white", colour = "black",size = 1.5),
          legend.key.heigh = unit(1.5,"line"),
          legend.key.width = unit(4,"line"),
          legend.box.just = "right",
          legend.box.background = element_rect(color="black", size=0.5),
          legend.box.margin = margin(6, 6, 6, 6),
          legend.text=element_text(size=20),
          legend.title=element_text(size=22),
          legend.justification = "center",
          legend.margin = margin(0,0,0,0),
          legend.background=element_rect(fill="white",colour="black",size=0,linetype="solid"),
          panel.border = element_rect(colour = "black", fill=NA, size=0.3),
          # axis.text.x = element_text(size=18,color="black"),axis.text.y = element_text(size=18,color="black"),
          # axis.title.x = element_text(size=18),axis.title.y = element_text(size=18),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.grid = element_blank(),
          plot.margin = margin(0,0,0,0, "cm"),#margin(0,0,0,0,unit="pt"),
          panel.background = element_rect(fill="white", colour="black", size=0.4,linetype=1))
  #legend <- get_legend(insetMap)
  #saveRDS(legend,paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/paper2/Legend.Marginality.rds"))
  # 
  Inset<-plot_grid(insetMap,align="hv") #Need to trasnform it to a grid (for the inset)
  myinset <- ggplotGrob(Inset) #transform to a grob 

  
  test7 <- test6[test6$x>=MyXrange[1]&test6$x<=MyXrange[2]&test6$y>=MyYrange[1]&test6$y<=MyYrange[2],c(1:2)] # define new range for the raster object 
  
  # Now plot the map with observed recruit as a contuinous variable 
  p <- ggplot() + 
    theme(panel.background = element_rect(fill="white", colour="black", size=3,
                                          linetype=1, color="black")) +
    theme(panel.grid.major = element_line(size=0, linetype=0,lineend="square", color="red"),
          panel.grid.minor = element_line(size=0, linetype=0,lineend="round", color="blue"))+
    scale_y_continuous(expand = c(0.02,0)) +
    scale_x_continuous(expand = c(0.02,0)) +
    geom_raster(data=test7[,c(1:2)],aes(x=x, y=y),fill="gray80",color="gray80",interpolate=F,na.rm=T) +
    geom_polygon(data = europe, aes(x=long, y = lat,group=group),fill=NA,col="black") + #gray 60 & gray 60
    geom_point(data = r1, aes(x = X, y = Y,col=Zbis,fill=Zbis),shape=21,size=1)+
    coord_fixed(1.3,xlim=MyXrange,ylim = MyYrange,expand=T) +
    # xlim(min(r1$X),max(r1$X))+
    # ylim(min(r1$Y),max(r1$Y))+
    scale_fill_gradient(name=NULL,low = "yellow", high = "red4", na.value = NA,limits=c(min(r1$Zbis,r1$Z),max(r1$Zbis,r1$Z)))+
    scale_color_gradient(name=NULL, low = "yellow", high = "red4", na.value = NA,limits=c(min(r1$Zbis,r1$Z),max(r1$Zbis,r1$Z)))+
    guides(shape = guide_legend(override.aes = list(size = 5)))+
    labs(y=paste0("Latitude"), x="Longitude",title=paste0('Recruitment rates of ',Allcode[i]," (log)"))+
    theme(text = element_text(face="bold"),
          legend.position = "right",
          legend.justification = c(1, 1),
          legend.direction ="vertical",
          # legend.box.just = "center",
          # legend.box.background = element_rect(color="white", size=0.5),
          # legend.box.margin = margin(0, 0, 0, -5),
          axis.text.x = element_text(size=18,color="black"),axis.text.y = element_text(size=18,color="black"),
          axis.title.x = element_text(size=18),axis.title.y = element_text(size=18),
          legend.key = element_rect(fill = "gray60", colour = "black",size = 1),
          legend.key.heigh = unit(4,"line"),
          legend.key.width = unit(0.8,"line"),
          legend.text=element_text(size=14),
          #legend.justification = "center",
          legend.margin = margin(0,0,0,0),
          legend.background=element_rect(fill="white",colour="black",size=0,linetype="solid"),
          panel.border = element_rect(colour = "black", fill=NA, size=0.3),
          axis.line = element_line(colour="black"),
          plot.title = element_text(size=17,hjust = 0,vjust = 0),
          plot.margin = margin(0,0,0,0, "cm"),
          plot.caption = element_text(face="bold.italic"))
  p
  # To check how it looks
  # Here save the object 
  # save_plot(filename = paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/paper2/Pobs.",Allcode[i],".png"),plot = pobs,base_width = 12.34, base_height = 7.04, dpi = 300 ,units = "in",nrow=1,ncol=1)
  # saveRDS(pall,paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/paper2/Pobs.",Allcode[i],".rds"))
  # #assign(paste0("pall",Allcode[i]),pall2,envir = .GlobalEnv)


  # Here us the zoomed map with predicted recruitment by quantile. 
  #   p <- ggplot() + 
  #   theme(panel.background = element_rect(fill="white", colour="black", size=3,
  #                                         linetype=1, color="black")) +
  #   theme(panel.grid.major = element_line(size=0, linetype=0,lineend="square", color="red"),
  #         panel.grid.minor = element_line(size=0, linetype=0,lineend="round", color="blue"))+
  #   #geom_raster(data=test6[,c(1:2)],aes(x=x, y=y),fill="gray80",color="gray80",interpolate=T,na.rm=T) +
  #   geom_polygon(data = europe, aes(x=long, y = lat,group=group),fill=NA,col="black") + #gray 60 & gray 60
  #   coord_fixed(1.3,xlim=c(min(r1$X),max(r1$X)),ylim=c(min(r1$Y),max(r1$Y))) + 
  #   #coord_cartesian(1.3,xlim=c(min(r1$X),max(r1$X)),ylim=c(min(r1$Y),max(r1$Y)),clip = 'off') + 
  #   #geom_point(data = r1, aes(x = X, y = Y,col=Zbis,fill=Zbis),shape=21,size=1)+
  #   geom_point(data = r1, aes(x = X, y = Y, group=bin, col = bin, shape = bin, size=bin),alpha=0.55)+
  #   scale_colour_manual(values = c("red","green","blue"),name="Recruitment events by plot",labels = c("High (>Q3)", "Medium \n(Q1 - Q3)","Low (<Q1)")) +
  #   scale_shape_manual(values = c(20,20,20),name="Recruitment events by plot",labels = c("High (>Q3)", "Medium \n(Q1 - Q3)","Low (<Q1)")) +
  #   scale_size_manual(values= c(2.5,2.5,2.5),name="Recruitment events by plot",labels = c("High (>Q3)", "Medium \n(Q1 - Q3)","Low (<Q1)")) +
  #   guides(shape = guide_legend(override.aes = list(size = 15)))+
  #   #coord_cartesian(xlim=c(min(r1$X),max(r1$X)),ylim=c(min(r1$Y),max(r1$Y)),clip = 'off')+
  #   labs(y=paste0("Latitude"), x="Longitude",title=paste0('Predicted recruitment rates of ',Allcode[i]))+
  #   theme(text = element_text(face="bold"),
  #         legend.direction ="vertical",#legend.position = "none",
  #         axis.text.x = element_text(size=18,color="black"),axis.text.y = element_text(size=18,color="black"),
  #         axis.title.x = element_text(size=18),axis.title.y = element_text(size=18),
  #         legend.key = element_rect(fill = "white", colour = "black",size = 1),
  #         legend.key.heigh = unit(4,"line"),
  #         legend.key.width = unit(1.5,"line"),
  #         legend.text=element_text(size=16),
  #         legend.title=element_text(size=18),
  #         legend.justification = "center",
  #         legend.margin = margin(0,0,0,0),
  #         legend.background=element_rect(fill="white",colour="black",size=0,linetype="solid"),
  #         panel.border = element_rect(colour = "black", fill=NA, size=0.3),
  #         axis.line = element_line(colour="black"),
  #         plot.title = element_text(size=17,hjust = 0,vjust = 0),
  #         plot.margin = margin(0,0,0,0, "cm"),
  #         plot.caption = element_text(face="bold.italic"))
  # 
  # p check how it looks
  
  # Here save the object 
  #save_plot(filename = paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/paper2/TESTTTTTTT",Allcode[i],".png"),plot = ptest,base_width = 12, base_height = 7.04, dpi = 300 ,units = "in",nrow=1,ncol=1)
  # saveRDS(pall2,paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/paper2/Ppred.",Allcode[i],".rds"))
  # #assign(paste0("pall",Allcode[i]),pall2,envir = .GlobalEnv)
  
  #### On recup les coordonées sur l'un des deux plot 
  xrange <- unlist(ggplot_build(p)$layout$panel_params[[1]][1])
  yrange <- unlist(ggplot_build(p)$layout$panel_params[[1]][8])
  xMin = min(xrange)
  xMax = max(xrange)
  xdelta = xMax-xMin
  yMin = min(yrange)
  yMax = max(yrange)
  ydelta = yMax-yMin
  
  ### Now plot the figures with the inset according to the locations (two groups, bottom right corner and topleft)
  if (Allcode[i]%in%c("PICABI","QUEPYR","PINNIG","PINHAL","PINPINA","PINPIN","QUEILE","QUESUB","PINSYL","CASSAT")){
    pallV2 <- p + annotation_custom(myinset, xmin=xMax-xdelta*0.38-(xdelta*0.38)*0.03, xmax=xMax-(xdelta*0.38)*0.03, 
                                    ymin = yMin-ydelta*0.02, ymax=yMin+ydelta*0.38-ydelta*0.02)  
    save_plot(filename = paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/paper2/SINSET.Pobs.",Allcode[i],".png"),plot = pallV2,base_width = 10, base_height = 7.04, dpi = 300 ,units = "in",nrow=1,ncol=1)
    saveRDS(pallV2,paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/paper2/INSET.Pobs.",Allcode[i],".rds"))

  }else if (Allcode[i]%in%c("FRAEXC","QUEROB","ACEPSE","ABIALB","QUEPET","BETPEN","FAGSYL","ALNGLU")){
    pallV2 <- p + annotation_custom(myinset, xMin+(xdelta*0.38)*0.03, xMin+(xdelta*0.38)+(xdelta*0.38)*0.03, 
                                    ymin = yMax-ydelta*0.38+ydelta*0.02, ymax=yMax+ydelta*0.02)
    save_plot(filename = paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/paper2/SINSET.Pobs.",Allcode[i],".png"),plot = pallV2,base_width = 10.04, base_height = 7.04, dpi = 300 ,units = "in",nrow=1,ncol=1)
    saveRDS(pallV2,paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/paper2/INSET.Pobs.",Allcode[i],".rds"))
  }
  
  # if (Allcode[i]%in%c("PICABI","QUEPYR")){
  #   pallV2 <- p + annotation_custom(myinset, xmin=xMax-xdelta*0.38-(xdelta*0.38)*0.01, xmax=xMax-(xdelta*0.38)*0.01, 
  #                                              ymin = yMin+(ydelta*0.38)*0.02, ymax=yMin+ydelta*0.38+(ydelta*0.38)*0.02)  
  #   save_plot(filename = paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/paper2/SINSET.Pobs.",Allcode[i],".png"),plot = pallV2,base_width = 10, base_height = 7.04, dpi = 300 ,units = "in",nrow=1,ncol=1)
  #   saveRDS(pallV2,paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/paper2/INSET.Pobs.",Allcode[i],".rds"))
  #   
  #   
  #   }else if (Allcode[i]%in%c("PINNIG","PINHAL")){
  #     pallV2 <- p + annotation_custom(myinset, xmin=xMax-xdelta*0.38+(xdelta*0.38)*0.06, xmax=xMax+(xdelta*0.38)*0.06, 
  #                                     ymin = yMin+(ydelta*0.38)*0.02, ymax=yMin+ydelta*0.38+(ydelta*0.38)*0.02)   
  #     save_plot(filename = paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/paper2/SINSET.Pobs.",Allcode[i],".png"),plot = pallV2,base_width = 10, base_height = 7.04, dpi = 300 ,units = "in",nrow=1,ncol=1)
  #     saveRDS(pallV2,paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/paper2/INSET.Pobs.",Allcode[i],".rds"))
  #     
  #   }else if (Allcode[i]%in%c("PINPINA")){
  #     pallV2 <- p + annotation_custom(myinset, xmin=xMax-xdelta*0.38+(xdelta*0.38)*0.09, xmax=xMax+(xdelta*0.38)*0.09, 
  #                                     ymin = yMin+(ydelta*0.38)*0.02, ymax=yMin+ydelta*0.38+(ydelta*0.38)*0.02)   
  #     save_plot(filename = paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/paper2/SINSET.Pobs.",Allcode[i],".png"),plot = pallV2,base_width = 10, base_height = 7.04, dpi = 300 ,units = "in",nrow=1,ncol=1)
  #     saveRDS(pallV2,paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/paper2/INSET.Pobs.",Allcode[i],".rds"))
  #   
  #   }else if (Allcode[i]%in%c("PINPIN","QUEILE","QUESUB")){ ##OK 
  #     pallV2 <- p + annotation_custom(myinset, xmin=xMax-xdelta*0.38+(xdelta*0.38)*0.08, xmax=xMax+(xdelta*0.38)*0.08, 
  #                                     ymin = yMin+(ydelta*0.38)*0.02, ymax=yMin+ydelta*0.38+(ydelta*0.38)*0.02)   
  #     save_plot(filename = paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/paper2/SINSET.Pobs.",Allcode[i],".png"),plot = pallV2,base_width = 10, base_height = 7.04, dpi = 300 ,units = "in",nrow=1,ncol=1)
  #     saveRDS(pallV2,paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/paper2/INSET.Pobs.",Allcode[i],".rds"))
  #     
  #   
  #   }else if (Allcode[i]%in%c("FRAEXC")){
  #     pallV2 <- p + annotation_custom(myinset, xMin-(xdelta*0.38)*0.09, xMin+(xdelta*0.38)-+(xdelta*0.38)*0.09, 
  #                                     ymin = yMax-(ydelta*0.38)*0.02-ydelta*0.38, ymax=yMax-(ydelta*0.38)*0.02) 
  #     save_plot(filename = paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/paper2/SINSET.Pobs.",Allcode[i],".png"),plot = pallV2,base_width = 10, base_height = 7.04, dpi = 300 ,units = "in",nrow=1,ncol=1)
  #     saveRDS(pallV2,paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/paper2/INSET.Pobs.",Allcode[i],".rds"))
  #     
  #   }else if (Allcode[i]%in%c("QUEROB","ACEPSE","ABIALB")){
  #     pallV2 <- p + annotation_custom(myinset, xMin-(xdelta*0.38)*0.08, xMin+(xdelta*0.38)-+(xdelta*0.38)*0.08, 
  #                                     ymin = yMax-(ydelta*0.38)*0.03-ydelta*0.38, ymax=yMax-(ydelta*0.38)*0.03) 
  #     save_plot(filename = paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/paper2/SINSET.Pobs.",Allcode[i],".png"),plot = pallV2,base_width = 10, base_height = 7.04, dpi = 300 ,units = "in",nrow=1,ncol=1)
  #     saveRDS(pallV2,paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/paper2/INSET.Pobs.",Allcode[i],".rds"))
  #     
  #   }else if (Allcode[i]%in%c("QUEPET")){
  #     pallV2 <- p + annotation_custom(myinset, xMin-(xdelta*0.38)*0.09, xMin+(xdelta*0.38)-+(xdelta*0.38)*0.09, 
  #                                     ymin = yMax-(ydelta*0.38)*0.03-ydelta*0.38, ymax=yMax-(ydelta*0.38)*0.03) 
  #     save_plot(filename = paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/paper2/SINSET.Pobs.",Allcode[i],".png"),plot = pallV2,base_width = 10, base_height = 7.04, dpi = 300 ,units = "in",nrow=1,ncol=1)
  #     saveRDS(pallV2,paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/paper2/INSET.Pobs.",Allcode[i],".rds"))
  #     
  #     
  # }else if (Allcode[i]%in%c("PINSYL","CASSAT")){ # Normalement OK
  #   pallV2 <- p + annotation_custom(myinset, xMin+(xdelta*0.38)*0.02, xMin+(xdelta*0.38)+(xdelta*0.38)*0.02, 
  #                                   ymin = yMax-(ydelta*0.38)*0.02-ydelta*0.38, ymax=yMax-(ydelta*0.38)*0.02)
  #   save_plot(filename = paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/paper2/SINSET.Pobs.",Allcode[i],".png"),plot = pallV2,base_width = 10.04, base_height = 7.04, dpi = 300 ,units = "in",nrow=1,ncol=1)
  #   saveRDS(pallV2,paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/paper2/INSET.Pobs.",Allcode[i],".rds"))
  #   
  # }else if (Allcode[i]%in%c("BETPEN","FAGSYL")){ # Normalement OK
  #   pallV2 <- p + annotation_custom(myinset, xMin+(xdelta*0.38)*0.03, xMin+(xdelta*0.38)+(xdelta*0.38)*0.03, 
  #                                   ymin = yMax-(ydelta*0.38)*0.02-ydelta*0.38, ymax=yMax-(ydelta*0.38)*0.02)
  #   save_plot(filename = paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/paper2/SINSET.Pobs.",Allcode[i],".png"),plot = pallV2,base_width = 10.04, base_height = 7.04, dpi = 300 ,units = "in",nrow=1,ncol=1)
  #   saveRDS(pallV2,paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/paper2/INSET.Pobs.",Allcode[i],".rds"))
  #   
  # }else if (Allcode[i]%in%c("ALNGLU")){ # Normalement OK
  #   pallV2 <- p + annotation_custom(myinset, xMin+(xdelta*0.38)*0.03, xMin+(xdelta*0.38)+(xdelta*0.38)*0.03, 
  #                                   ymin = yMax-(ydelta*0.38)*0.02-ydelta*0.38, ymax=yMax-(ydelta*0.38)*0.02)
  #   save_plot(filename = paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/paper2/SINSET.Pobs.",Allcode[i],".png"),plot = pallV2,base_width = 10.04, base_height = 7.04, dpi = 300 ,units = "in",nrow=1,ncol=1)
  #   saveRDS(pallV2,paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/paper2/INSET.Pobs.",Allcode[i],".rds"))
  #   } 
  rm("x")
}
 

## Check resolution and cowplot
# legend <- get_legend(insetMap)
# legend2 <- get_legend(p)
# saveRDS(legend,paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/paper2/Legend.Marginality.rds"))
# saveRDS(legend2,paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/paper2/Legend.Quantile.rds"))


# raster version 
#A <- list.files(paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/paper2/",Modname,"/Figure.maps"),pattern = ".rds",full.names = T)
#version map only 
A <- list.files(paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/paper2"),pattern = ".rds",full.names = T)

# Load only the first half 
mapply(function(x,y){
  assign(paste0("pall",y),readRDS(x),envir = .GlobalEnv)
},x=A[1:9],y=Allcode[1:9])


pall <- plot_grid(
  pallABIALB+theme(axis.title.x=element_blank(),
                   #axis.title.y=element_blank(),
                   plot.margin = unit(c(0,0,0.2,0), "cm")),
  pallACEPSE+theme(axis.title.x=element_blank(),
                  axis.title.y=element_blank(),
                  plot.margin = unit(c(0,0,0.2,0), "cm")),
  pallALNGLU+theme(axis.title.x=element_blank(),
                   axis.title.y=element_blank(),
                   plot.margin = unit(c(0,0,0.2,0), "cm")),
  pallBETPEN+theme(axis.title.x=element_blank(),
                   #axis.title.y=element_blank(),
                   plot.margin = unit(c(0,0,0.2,0), "cm")),
  pallCASSAT+theme(axis.title.x=element_blank(),
                   axis.title.y=element_blank(),
                   plot.margin = unit(c(0,0,0.2,0), "cm")),
  pallFAGSYL+theme(axis.title.x=element_blank(),
                   axis.title.y=element_blank(),
                   plot.margin = unit(c(0,0,0.2,0), "cm")),
  pallFRAEXC+theme(#axis.title.x=element_blank(),
                   #axis.title.y=element_blank(),
                   plot.margin = unit(c(0,0,0,0), "cm")),
  pallPICABI+theme(#axis.title.x=element_blank(),
                   axis.title.y=element_blank(),
                   plot.margin = unit(c(0,0,0,0), "cm")),
  pallPINHAL+theme(#axis.title.x=element_blank(),
                    axis.title.y=element_blank(),
                    plot.margin = unit(c(0,0,0,0), "cm")),ncol=3,nrow=3)

save_plot(filename = paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/paper2/",Modname,"/Figure.maps/MapsAllPART1.Recrut2.png"),plot = pall,base_width = 5, base_height = 6, dpi = 300 ,units = "in",nrow=3,ncol=3)
save_plot(filename = paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/paper2/",Modname,"/Figure.maps/MapsAllPART1.Recrut.pdf"),plot = pall,base_width = 5, base_height = 6, dpi = 300 ,units = "in",nrow=3,ncol=3)


## Free the space used !!
rm(list = ls())
gc()

Modname <- "rW.yO.sumO.mRY1"
Allcode <- c("ABIALB","ACEPSE","ALNGLU","BETPEN","CASSAT","FAGSYL",
             "FRAEXC","PICABI","PINHAL","PINNIG","PINPIN","PINPINA","PINSYL"
             ,"QUEILE","QUEPET","QUEPYR","QUEROB","QUESUB")
## Load the second half and the legend:
#### load the legend
Marginality <- readRDS(paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/paper2/Legend.Marginality.rds"))

A <- list.files(paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/paper2"),pattern = ".rds",full.names = T)
# Load the second half 
mapply(function(x,y){
  assign(paste0("pall",y),readRDS(x),envir = .GlobalEnv)
},x=A[10:18],y=Allcode[10:18])



pallbis <- plot_grid(
  pallPINNIG+theme(axis.title.x=element_blank(),
                   #axis.title.y=element_blank(),
                   plot.margin = unit(c(0,0,0.2,0), "cm")),
  pallPINPIN+theme(axis.title.x=element_blank(),
                   axis.title.y=element_blank(),
                   plot.margin = unit(c(0,0,0.2,0), "cm")),
  pallPINPINA+theme(axis.title.x=element_blank(),
                   axis.title.y=element_blank(),
                   plot.margin = unit(c(0,0,0.2,0), "cm")),
  pallPINSYL+theme(axis.title.x=element_blank(),
                   #axis.title.y=element_blank(),
                   plot.margin = unit(c(0,0,0.2,0), "cm")),
  pallQUEILE+theme(axis.title.x=element_blank(),
                   axis.title.y=element_blank(),
                   plot.margin = unit(c(0,0,0.2,0), "cm")),
  pallQUEPET+theme(axis.title.x=element_blank(),
                   axis.title.y=element_blank(),
                   plot.margin = unit(c(0,0,0.2,0), "cm")),
  pallQUEPYR+theme(#axis.title.x=element_blank(),
    #axis.title.y=element_blank(),
    plot.margin = unit(c(0,0,0,0), "cm")),
  pallQUEROB+theme(#axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    plot.margin = unit(c(0,0,0,0), "cm")),
  pallQUESUB+theme(#axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    plot.margin = unit(c(0,0,0,0), "cm")),ncol=3,nrow=3)
# save the 9 species without legend
save_plot(filename = paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/paper2/",Modname,"/Figure.maps/MapsAllPART1.Recrut2.part2.png"),plot = pall,base_width = 5, base_height = 6, dpi = 300 ,units = "in",nrow=3,ncol=3)


pall.2 <- plot_grid(pallbis,Marginality,rel_heights = c(1.1,0.05),nrow=2,ncol=1)
# save with the legend
save_plot(filename = paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/paper2/",Modname,"/Figure.maps/MapsAllPART1.Recrut2.part2.legend.png"),plot = pall.2,base_width = 5, base_height = 6, dpi = 300 ,units = "in",nrow=3,ncol=3)
save_plot(filename = paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/paper2/",Modname,"/Figure.maps/MapsAllPART1.Recrut2.part2.legend.pdf"),plot = pall.2,base_width = 5, base_height = 6, dpi = 300 ,units = "in",nrow=3,ncol=3)

