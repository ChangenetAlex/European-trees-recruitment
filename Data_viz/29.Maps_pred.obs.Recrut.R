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
library(svglite)
library("ggplotify")


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
for (i in c(1:3,17)){ #4
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
  
  print(summary(r1$Z))
  print(summary(r1$Zbis))
  
  r1$Z <- log(1+r1$Z)
  r1$Zbis <- log(1+r1$Zbis)
  
  par(mfrow=c(1,1))
  worldmap <- getMap(resolution = "high")
  #europe <- worldmap[which(worldmap$REGION=="Europe"),]
  #europe <- worldmap[which(grepl("Europe", worldmap$GEO3) & as.character(worldmap$NAME) != "Russia"),] 
  europe <- worldmap[which(grepl("Europe", worldmap$GEO3)|grepl("North Africa", worldmap$GEO3)),] 
  europe <-spTransform(europe,"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0 ") # Convert to the right system
  #europe1 = crop(europe,c(-15, 45, 35, 70))
  europe1 = crop(europe,c(-10, 32, 36, 70))
  europe = crop(europe,c(min(r1$X), max(r1$X), min(r1$Y), max(r1$Y)))
  plot(europe1)
  # Here is the inset map in which we will zoom
  insetMap <-   ggplot() + 
    labs(x = NULL, y = NULL)+
    geom_raster(data=test6,aes(x=x, y=y),fill="gray80",color="gray80",interpolate=T) +
    xlim(-10,32)+
    ylim(36,70)+
    coord_fixed(1.3) + geom_polygon(data = europe1, aes(x=long, y = lat,group=group),fill=NA,col="black") + #gray 60 & gray 60
    #geom_point(data = r1, aes(x = X, y = Y,col=Z,fill=Z),shape=21,size=1)+
    scale_fill_gradient(low = "yellow", high = "red4", na.value = NA)+
    scale_color_gradient(low = "yellow", high = "red4", na.value = NA)+
    geom_rect(data = data.frame(),aes(xmin=min(r1$X), xmax=max(r1$X), ymin=min(r1$Y), ymax=max(r1$Y)),colour = "red", fill = NA,size=1)+
    theme(legend.position = "none",
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            panel.grid = element_blank(),
            plot.margin = margin(0,0,0,0, "cm"),#margin(0,0,0,0,unit="pt"),
            panel.border = element_rect(fill = NA, colour = "black"),
            panel.background = element_rect(fill="white", colour="black", size=1.2,linetype=1))
  
  # Here us the zoomed map
  
  ppred <- ggplot() + 
    theme(panel.background = element_rect(fill="white", colour="black", size=3,
                                                        linetype=1, color="black"),legend.key.size = unit(1.5, "cm"),legend.position = c(0.5,-0.2)) +
     theme(panel.grid.major = element_line(size=0, linetype=0,lineend="square", color="red"),
           panel.grid.minor = element_line(size=0, linetype=0,lineend="round", color="blue"))+
    geom_raster(data=test6,aes(x=x, y=y),fill="gray80",color="gray80",interpolate=T) +
    geom_polygon(data = europe, aes(x=long, y = lat,group=group),fill=NA,col="black") + #gray 60 & gray 60
    coord_fixed(1.3) + geom_point(data = r1, aes(x = X, y = Y,col=Z,fill=Z),shape=21,size=1)+
    scale_fill_gradient(name=NULL,low = "yellow", high = "red4", na.value = NA,limits=c(min(r1$Zbis,r1$Z),max(r1$Zbis,r1$Z)))+
    scale_color_gradient(name=NULL,low = "yellow", high = "red4", na.value = NA,limits=c(min(r1$Zbis,r1$Z),max(r1$Zbis,r1$Z)))
    
  ppred <- ppred + guides(shape = guide_legend(override.aes = list(size = 5)))+
    xlim(min(r1$X),max(r1$X))+
    ylim(min(r1$Y),max(r1$Y))+
    labs(y=paste0("Latitude"), x="Longitude",title=paste0('b) Fitted recruitment rates of ',Allcode[i]," (log)"))+
    theme(text = element_text(face="bold"),legend.direction ="horizontal",legend.position = c(0.5,-0.2),
          axis.text.x = element_text(size=18,color="black"),axis.text.y = element_text(size=0,color="white"),
          axis.title.x = element_text(size=18,color="black"),axis.title.y = element_text(size=0,color="white"),
          legend.key = element_rect(fill = "gray60", colour = "black",size = 0.1),
          legend.key.heigh = unit(0.5,"line"),
          legend.key.width = unit(2,"line"),
          legend.text=element_text(size=12),
          legend.background=element_rect(fill="white",colour="black",size=0,linetype="solid"),
          #panel.border = element_rect(colour = "black", fill=NA, size=0.3),
          axis.line = element_line(colour="black"),
          plot.title = element_text(size=17,hjust = 0,vjust = 0),
          plot.margin = margin(0,0,1,0, "cm"),
          plot.caption = element_text(face="bold.italic"))

  
  pobs <- ggplot() + 
    theme(panel.background = element_rect(fill="white", colour="black", size=3,
                                          linetype=1, color="black")) +
    theme(panel.grid.major = element_line(size=0, linetype=0,lineend="square", color="red"),
          panel.grid.minor = element_line(size=0, linetype=0,lineend="round", color="blue"))+
    geom_raster(data=test6,aes(x=x, y=y),fill="gray80",color="gray80",interpolate=T) +
    geom_polygon(data = europe, aes(x=long, y = lat,group=group),fill=NA,col="black") + #gray 60 & gray 60
    coord_fixed(1.3) + geom_point(data = r1, aes(x = X, y = Y,col=Zbis,fill=Zbis),shape=21,size=1)+
    scale_fill_gradient(name=NULL,low = "yellow", high = "red4", na.value = NA,limits=c(min(r1$Zbis,r1$Z),max(r1$Zbis,r1$Z)))+
    scale_color_gradient(name=NULL, low = "yellow", high = "red4", na.value = NA,limits=c(min(r1$Zbis,r1$Z),max(r1$Zbis,r1$Z)))
  
  pobs <- pobs + guides(shape = guide_legend(override.aes = list(size = 5)))+
    xlim(min(r1$X),max(r1$X))+
    ylim(min(r1$Y),max(r1$Y))+
    labs(y=paste0("Latitude"), x="Longitude",title=paste0('a) Observed recruitment rates of ',Allcode[i]," (log)"))+
    theme(text = element_text(face="bold"),legend.direction ="horizontal",legend.position = c(0.48,-0.2),
          axis.text.x = element_text(size=18,color="black"),axis.text.y = element_text(size=18,color="black"),
          axis.title.x = element_text(size=18),axis.title.y = element_text(size=18),
          legend.key = element_rect(fill = "gray60", colour = "black",size = 0.1),
          legend.key.heigh = unit(0.5,"line"),
          legend.key.width = unit(2,"line"),
          legend.text=element_text(size=12),
          legend.background=element_rect(fill="white",colour="black",size=0,linetype="solid"),
          #panel.border = element_rect(colour = "black", fill=NA, size=0.3),
          axis.line = element_line(colour="black"),
          plot.title = element_text(size=17,hjust = 0,vjust = 0),
          plot.margin = margin(0,0.1,1,0, "cm"),
          plot.caption = element_text(face="bold.italic"))
  
  
  #pall<-plot_grid(pobs + theme(legend.position = "none"),
  
pall<-plot_grid(pobs + theme(legend.position = "none"),ppred + theme(legend.direction ="horizontal",legend.position = c(-0.07,-0.13),
                                                                    legend.margin = margin(0,0,0,0),
                                                                    legend.text.align=0,
                                                                    legend.justification = "center",
                                                                    legend.key.heigh = unit(1,"line"),
                                                                    legend.key.width = unit(4,"line"),
                                                                    legend.text=element_text(size=13)),align="hv")


  
  myinset <- ggplotGrob(insetMap)
  # pall2 <- pall+annotation_custom(myinset, xmin = 0.33, xmax = 0.63, # PICABI, PINHAL, PINNIG, PINPIN, PINPINA, QUEILE, QUEPYR, QUESUB
  #                                 ymin = 0.26, ymax = 0.56)
  # 
  #pall2 <- pall+annotation_custom(myinset, xmin = 0, xmax = 0.30, # BETPEN,CASSAT,FAGSYL,FRAEXC,PINSYL,QUEPET
  #                               ymin = 0.58, ymax = 0.88)
  # 
  pall2 <- pall+annotation_custom(myinset, xmin = 0.02, xmax = 0.32, # ABIALB,ACEPSE,ALNGLU, QUEROB
                                  ymin = 0.47, ymax = 0.77)

  ggsave(paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/paper2/",Modname,"/Figure.maps/",Allcode[i],".mapW.png"),plot=pall2,width = 12.37, height = 7.04,units = "in",dpi=400)
  assign(paste0("pall",Allcode[i]),pall2,envir = .GlobalEnv)
  saveRDS(pall2,paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/paper2/",Modname,"/Figure.maps/",Allcode[i],".mapW.rds"))
  graphics.off()
}

A <- list.files(paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/paper2/",Modname,"/Figure.maps/"),pattern = ".rds",full.names = T)

mapply(function(x,y){
  assign(paste0("pall",y),readRDS(x),envir = .GlobalEnv)
},x=A,y=Allcode)

pall <- plot_grid(
  pallABIALB,pallACEPSE,pallALNGLU,pallBETPEN,pallCASSAT,pallFAGSYL, 
  pallFRAEXC, pallPICABI, pallPINHAL, pallPINNIG, pallPINPIN, pallPINPINA, 
  pallPINSYL, pallQUEILE, pallQUEPET, pallQUEPYR, pallQUEROB,pallQUESUB,align="hv",ncol=3,nrow=6)
save_plot(filename = paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/paper2/",Modname,"/Figure.maps/MapsAllV2.Recrut.pdf"),plot = pall,base_width = 12.34, base_height = 7.04, dpi = 800 ,units = "in",nrow=6,ncol=3)


pall <- plot_grid(
  pallABIALB,pallACEPSE,pallALNGLU,pallBETPEN,pallCASSAT,pallFAGSYL,align="hv",ncol=2,nrow=3)
save_plot(filename = paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/paper2/",Modname,"/Figure.maps/MapsAllPART1.Recrut.svg"),plot = pall,base_width = 12.34, base_height = 7.04, dpi = 800 ,units = "in",nrow=3,ncol=2)

pall <- plot_grid(
  pallFRAEXC, pallPICABI, pallPINHAL, pallPINNIG, pallPINPIN, pallPINPINA,align="hv",ncol=2,nrow=3)
save_plot(filename = paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/paper2/",Modname,"/Figure.maps/MapsAllPART2.Recrut.svg"),plot = pall,base_width = 12.34, base_height = 7.04, dpi = 800 ,units = "in",nrow=3,ncol=2)


pall <- plot_grid(
  pallPINSYL, pallQUEILE, pallQUEPET, pallQUEPYR, pallQUEROB,pallQUESUB,align="hv",ncol=2,nrow=3)
save_plot(filename = paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/paper2/",Modname,"/Figure.maps/MapsAllPART3.Recrut.svg"),plot = pall,base_width = 12.34, base_height = 7.04, dpi = 800 ,units = "in",nrow=3,ncol=2)


 



vp_b <- viewport(width = 1, height = 1, x = 0.5, y = 0.5)  # the larger map
vp_a <- viewport(width = 0.3, height = 0.3, x = 0.48, y = 0.41)  # PICABI, PINHAL, PINNIG, PINPIN, PINPINA, QUEILE, QUEPYR, QUESUB
#vp_a <- viewport(width = 0.3, height = 0.3, x = 0.15, y = 0.73)  # BETPEN,CASSAT,FAGSYL,FRAEXC,PINSYL,QUEPET
#vp_a <- viewport(width = 0.3, height = 0.3, x = 0.17, y = 0.62)  # ABIALB,ACEPSE,ALNGLU, QUEROB
print(pall, vp = vp_b)
print(insetMap, vp = vp_a)