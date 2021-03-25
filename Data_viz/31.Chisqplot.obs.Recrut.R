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
library(grid)            #
library(cowplot)          #

Modname <- "SCALEDRrate.yr.QualiF.TF.noY.nbinom1.V3" ## Name of the model 
Modname <- "rW.yO.sumO.mRY1"
Allcode <- c("ABIALB","ACEPSE","ALNGLU","BETPEN","CASSAT","FAGSYL",
             "FRAEXC","PICABI","PINHAL","PINNIG","PINPIN","PINPINA","PINSYL"
             ,"QUEILE","QUEPET","QUEPYR","QUEROB","QUESUB")
Allseuil <- c(0.7,0.55,0.6,0.8,0.7,0.7,0.7,0.8,0.7,0.7,0.7,0.8,0.8,0.8,0.7,0.7,0.8,0.7)
AllseuilC <- c(0.6,0.5,0.5,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6)
pdf(file=paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/paper2/",Modname,"/Figure.maps/","_","Intersection.Recrut.pdf"), width = 12.37, height = 7.04) # Save it 
#for (i in c(1:18)){ #4
for (i in c(1:18)){ #4
    
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
  
  # Make sure Plotcat are in right order
  REV <- c("PINPIN","PINNIG","QUEPET","QUEROB","PINPINA","POPTRE","ALNGLU","QUEPYR","QUESUB")
  if (Allcode[i]%in%REV){
    x$Plotcat2 <- x$Plotcat
    x[x$Plotcat2==1,"Plotcat"] <- 2
    x[x$Plotcat2==2,"Plotcat"] <- 1
  } 
  
  # Recup latitude and longitude 
  r1 <- data.frame(cbind(dfplot[which(rownames(dfplot)%in%rownames(x)),c("longitude","latitude")],x$fit,x$sp.recruitment.ha.weight.2,x$Plotcat))
  colnames(r1) <- c("X","Y","Z","Zbis","Cat")
  r1$Cat <- as.character(r1$Cat)
  r1[r1$Cat=="0","Cat"] <- "Core"
  r1[r1$Cat=="2","Cat"] <- "LE" # To change when reverse. 2=LE in normal. 2=RE in reverse
  r1[r1$Cat=="1","Cat"] <- "TE"
  r1$bin <- cut(r1$Zbis,3)
  A <- table(r1$bin,r1$Cat)
  B <- chisq.test(A)
  #A <- melt(prop.table(table(r1$bin,r1$Cat),2))
  A <- melt(B$residuals*-1)
  A[,"x2"] <- B[1]
  A[,"pvalue"] <- B[3]
  A[,c(3:5)] <- round(A[,c(3:5)],2)
  colnames(A) <- c("Range breaks","Zone","Proportion","x2","pvalue")
  
  p <- ggplot(data=A,aes(x=Zone,y=Proportion,fill=factor(`Range breaks`))) + 
    theme(panel.background = element_rect(fill="white", colour="black", size=3,linetype=1, color="black"),
          legend.key.size = unit(0.5, "cm"),legend.key.width = unit(1,"cm"),legend.position = c(0.65,1.05)) + # 0.955,0.85 verticale.
    theme(panel.grid.major = element_line(colour="blue", size=4, linetype=0,lineend="square", color="red"),
          panel.grid.minor = element_line(colour="blue", size=1, linetype=0,lineend="round", color="blue"))+
    theme(legend.background = element_rect(fill="white",size=1, linetype="solid",
                                           colour ="black"))+
    #theme(legend.position = "top")+ ## Allow us to disable the full legend
    geom_bar(stat="identity",position="dodge",width = 0.65,colour="black")+
    geom_text(aes(label=Proportion),vjust=-0.2,colour="black",position=position_dodge(0.7),size=5,fontface = "bold")+
    #geom_text(aes(label=`Range breaks`),vjust=-0.1,colour="black",position=position_dodge(0.7),size=3,fontface = "bold")+
    #annotate("label", x = 2, y = 0, label = paste0("X-squared = ",A$x2[1],"; p-value = ",A$pvalue[1]),size=4.7,fontface = "bold")+
    scale_fill_manual(values = c("blue", "green", "red"),name="Range breaks",
                       labels=levels(A$`Range breaks`))+
          
    labs(title=paste0(Allcode[i],": X-squared = ",A$x2[1],";\np-value = ",A$pvalue[1]),
         y=paste0("Repartition of residuals across range breaks"), x="Zone")+
     theme(text = element_text(face="bold"),legend.direction ="horizontal",
          axis.text.x = element_text(size=16,color="black"),axis.text.y = element_text(size=16,color="black"),
          axis.title.x = element_text(size=16),axis.title.y = element_text(size=16),
          legend.title=element_text(size=13), 
          legend.text=element_text(size=13),
          legend.key = element_rect(fill = "gray60", colour = "black",size = 0.5),
          legend.background=element_rect(fill="white",colour="black",size=0.5),
          panel.border = element_rect(colour = "black", fill=NA, size=0.5),
          axis.line = element_line(colour="black"),
          plot.title = element_text(size=18,hjust = 0.01,vjust = 0),
          plot.caption = element_text(face="bold.italic"))
print(p)
assign(paste0("p",Allcode[i]),p,envir = .GlobalEnv)
  #ggsave(filename = paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/paper2/",Modname,"/Figure.maps/",Allcode[i],"_","Intersection.Recrut.pdf"),plot = p, width = 12.37, height = 7.04, dpi=300,units = "in")
}
dev.off()

pall <- plot_grid(
  pABIALB,pACEPSE,pALNGLU,pBETPEN,pCASSAT,pFAGSYL, 
    pFRAEXC, pPICABI, pPINHAL, pPINNIG, pPINPIN, pPINPINA, 
    pPINSYL, pQUEILE, pQUEPET, pQUEPYR, pQUEROB,pQUESUB,align="hv",ncol=3,nrow=6)
save_plot(filename = paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/paper2/",Modname,"/Figure.maps/Intersection.Recrut.pdf"),plot = pall,base_width = 12.34, base_height = 7.04, dpi = 800 ,units = "in",nrow=7,ncol=3)
  
  
  