rm(list = ls())
gc()
library(glmmTMB)
library(car)
library(emmeans)
library(effects)
library(multcomp)
library(MuMIn)
library(DHARMa)
library(broom)
library(broom.mixed)
library(dotwhisker)
library(ggplot2); theme_set(theme_bw())
library(texreg)
library(xtable)
library(huxtable)
library(plyr)
library(beanplot)
library(parallel)

Modname <- "SCALEDRrate.yr.QualiF.TF.noY.nbinom1.V3" ## Name of the model 
Modname <- "rW.yO.sumO.mRY1"

Allcode <- c("ABIALB","ACEPSE","ALNGLU","BETPEN","CASSAT","FAGSYL",
             "FRAEXC","PICABI","PINHAL","PINNIG","PINPIN","PINPINA","PINSYL","POPNIG",
             "POPTRE","QUEILE","QUEPET","QUEPYR","QUEROB","QUESUB")
Allseuil <- c(0.7,0.55,0.6,0.8,0.7,0.7,0.7,0.8,0.7,0.7,0.7,0.8,0.8,0.7,0.7,0.8,0.7,0.7,0.8,0.7)
AllseuilC <- c(0.6,0.5,0.5,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6)
mod_ALL <- list()
simulationOutput <- list()
for (i in c(1:13,16:length(Allcode))){
  # for (i in c(1:3)){
  CODE <- Allcode[i]
  seuilC <- AllseuilC[i]
  seuil <- Allseuil[i]
  Dir = paste0("~/SynologyDrive/FUNDIV - NFI - Europe/our-data/species/",CODE,"/CLIMAP/Recrut.Mortality.2020/")
  setwd(Dir)
  # Charge the dfplot for simulation and residuals test dispersion etc... 
  dfplot <- readRDS(paste0("dfplotFinalV3",CODE,seuil,"R.M.rds"))
  
  #Load the desired model of the desired species 
  Dir <- c(paste0(Dir,Modname,"/")) # Directory 
  setwd(Dir)
  mod_ALL[[CODE]] <- get(load(file = paste0(Modname,".rda")))
  rm(x)
  
  ### Test #### 
  #we define the data as the part of the dataset on which the model has been fitted. (for every loop)
  #dfplot <- mod_ALL[[CODE]]$frame
  ## Test ###
  
  #Simulate the residuals based on the right dfplot base (above) and the right species
  #simulationOutput[[CODE]] <- simulateResiduals(mod_ALL[[CODE]], n = 3000,re.form=NA) 
  
  #remove dfplot of the current species to ensure it won't be used for next species 
  rm(dfplot) 
}

VAR <- c("sp.recruitment.ha.weight.2", "logdbhJ.IMall.plot.mean", 
         "logtreeNbrJ", "mean_spei12", "sqrtBA.O.plot.1", "logBAj.plot.1", 
         "sqrtBAIj.plot.1.mean.J.I", "logsp.mortality.plot.rate.yr")

mod_ALL <- mclapply(mod_ALL, function(x) {
  A <- x[[5]]
  colnames(A)[which(colnames(A)=="offset(log(yearsbetweensurveys))")] <- "yearsbetweensurveys"
  colnames(A)[which(colnames(A)=="offset(log(sp.SUM.ha.weight2))")] <- "sp.SUM.ha.weight2"
  A[,"yearsbetweensurveys"] <- exp(A[,"yearsbetweensurveys"])
  A[,"sp.SUM.ha.weight2"] <- exp(A[,"sp.SUM.ha.weight2"])
  x[[5]][,"fit"] <- predict(x,A,type=c("response"))
  x[[5]][,"diff"] <- x[[5]][,"sp.recruitment.ha.weight.2"]-x[[5]][,"fit"]
  if (length(which(colnames(x[[5]])=="country")!=0)){
    x[[5]] <- x[[5]][,-which(colnames(x[[5]])=="country")]
  }
  ;x[[5]]},mc.cores = 18)

## Convert it to a df
mod_ALL <- mapply(function(x,y) {
  x[,"species"] <- y
  return(x)
},
x=mod_ALL,y=names(mod_ALL),SIMPLIFY = F)

mod_ALL <- do.call(rbind,mod_ALL) # From list to df

# Make sure Plotcat are in right order
REV <- c("PINPIN","PINNIG","QUEPET","QUEROB","PINPINA","POPTRE","ALNGLU","QUEPYR","QUESUB")
mod_ALL$Plotcat2 <- mod_ALL$Plotcat
mod_ALL[mod_ALL$species%in%REV & mod_ALL$Plotcat2==1,"Plotcat"] <- 2
mod_ALL[mod_ALL$species%in%REV & mod_ALL$Plotcat2==2,"Plotcat"] <- 1

######################
## Here dharma check #
######################
Dir <- paste0("~/SynologyDrive/FUNDIV - NFI - Europe/our-data/species/Recrut.All/",Modname,"/")
dir.create(paste0(Dir,"diagnostic.test/")) # Creation of a new file for this model
setwd(paste0(Dir,"diagnostic.test/"))
# Split it again into a list
mod_ALL <- split(mod_ALL,mod_ALL$species)

# Chisq2
Mytest <- mclapply(mod_ALL, function(x) chisq.test(x[,"sp.recruitment.ha.weight.2"],x[,"fit"]),mc.cores = 18)

# !!! Need to add the column names acording to the test that is performed !!!!

#Dharma
sizeN = 1800
resol = 15
# First test dispersion of the sclaed residuals + test 
png(file=paste0(Dir,"diagnostic.test/","Interger_Dispersion.test.",sizeN,".png"),width=sizeN,height=sizeN*0.9,res=sizeN/resol) # Save it 
par(mfrow=c(4,5))
Dispersion <- mapply(function(x,y){
  testDispersion(x)
  mtext(y)},
  x=simulationOutput,y=names(simulationOutput)) # Okay to save all in a single plot
dev.off()
#Plots is saved but values of tests are not recorded anywhere !!
Dispersion <- lapply(simulationOutput,function(x){testDispersion(x,plot=F)})
# Values are recorded now !

# Second test QQplot = unifiormity distrib of the scaled residuals + test
# Te"st over and underdisp here
png(file=paste0(Dir,"diagnostic.test/","Integer.Uniformity.test.",sizeN,".png"),width=sizeN,height=sizeN*0.9,res=sizeN/resol) # Save it 
par(mfrow=c(4,5))
Uniformity <- mapply(function(x,y){
  testUniformity(x)
  mtext(y)},
  x=simulationOutput, y=names(simulationOutput),SIMPLIFY = F) # Okay to save all in a single plot
dev.off()
#Plot is saved but values of the tests are nowhere
Uniformity <- lapply(simulationOutput,function(x){testUniformity(x,plot=F)})


# Zero inf figure + tests 
png(file=paste0(Dir,"diagnostic.test/","ZI.test.",sizeN,".png"),width=sizeN,height=sizeN*0.9,res=sizeN/resol) # Save it 
par(mfrow=c(4,5))
ZeroInf <- lapply(simulationOutput,testZeroInflation)
dev.off()
# Plots is saved and values are recorded 
#Need to save the plots 


# Just plot heteroscedasticite (no test ??)
Residuals <- lapply(simulationOutput, function(x) plotResiduals(x)) #Just the plot mais pas standardized
Quantiles <- mclapply(simulationOutput, function(x) testQuantiles(x, predictor = NULL, quantiles = c(0.25, 0.5,
                                                                                                   0.75), plot = T),mc.cores = 18)
# Same thing as two previous ones with smooth lines and quantile regressions (longer but more accurate)
png(file=paste0(Dir,"diagnostic.test/","quantileslocation.test.",sizeN,".png"),width=sizeN,height=sizeN*0.9,res=sizeN/resol) # Save it 
par(mfrow=c(4,5))
Residuals2 <- mapply(function(x,y){
  plotResiduals(x, rank = TRUE, quantreg = T, smoothScatter = TRUE,
                main=paste0("Residual vs. predicted for ",y))
},
x=simulationOutput,y=names(simulationOutput),SIMPLIFY = F)
dev.off()

png(file=paste0(Dir,"diagnostic.test/","quantileslocation.test3.",sizeN,".png"),width=sizeN,height=sizeN*0.9,res=sizeN/resol) # Save it 
par(mfrow=c(4,5))
Residuals2 <- lapply(simulationOutput, function(x){
  plotResiduals(x, rank = TRUE, quantreg = T, smoothScatter = TRUE)
})
dev.off()

#Need to save the plots and to redo it with names of the species 

# Plot and tests
png(file=paste0(Dir,"diagnostic.test/","Outliers.test3.",sizeN,".png"),width=sizeN,height=sizeN*0.9,res=sizeN/resol) # Save it 
par(mfrow=c(4,5))
Outliers <- mapply(function(x,y){
  testOutliers(x)
  mtext(y)},
  x=simulationOutput, y=names(simulationOutput)) # Okay to save all in a single plot
dev.off()
#Need to save the plots 
Outliers <- lapply(simulationOutput,function(x){testOutliers(x,plot=F)})

#Saving these plots by species 
# On these we have two figures and the three tests are 
# KS test = uniformity test cf uniformiy,
# disp test = dispersion test (cf Dispersion),
# and last = outlier test (cf outliers)
All <- mapply(function(x,y){
  png(file=paste0(Dir,"diagnostic.test/Species.allplot/",y,"Allplot.test.",sizeN,".png"),width=sizeN,height=sizeN*0.9,res=sizeN/resol) # Save it 
  plot(x)
  dev.off()
},
x=simulationOutput,y=names(simulationOutput),SIMPLIFY = F)


# Add the colnames 
mapply(function(x,y,z,u,ZI,R,o) {
  cat(y,
      paste0(x[["statistic"]]," "," (",x[["p.value"]],")"),
      paste0(z[["statistic"]]," (",z[["p.value"]],")"),
      paste0(u[["statistic"]]," (",u[["p.value"]],")"),
      paste0(ZI[["statistic"]]," (",ZI[["p.value"]],")"),R[["p.value"]],o[["p.value"]],"\n",
      file = paste0(Dir,"diagnostic.test/All.TEST_Table.csv"),sep=";",append = T)
},
x=Mytest,y=names(Mytest),z=Dispersion,u=Uniformity,ZI=ZeroInf,R=Residuals2,o=Outliers,SIMPLIFY=F)


A <- mapply(function(x,y,z,u,ZI,R,o) {
  data.frame(y,
             paste0(round(x[["statistic"]],2),"\n(",round(x[["p.value"]],2),")"),
             paste0(round(z[["statistic"]],2),"\n(",round(z[["p.value"]],2),")"),
             paste0(round(u[["statistic"]],2),"\n(",round(u[["p.value"]],2),")"),
             paste0(round(ZI[["statistic"]],2),"\n(",round(ZI[["p.value"]],2),")"),round(R[["p.value"]],2),round(o[["p.value"]],2))
},
x=Mytest,y=names(Mytest),z=Dispersion,u=Uniformity,ZI=ZeroInf,R=Residuals2,o=Outliers,SIMPLIFY=F)
A <- do.call(rbind,A) # From list to df
colnames(A) <- c("Species","Chi-square\n(p-value)",
                 "Dispersion\n(p-value)",
                 "Uniformity KS\n(p-value)",
                 "Zero-inflation\n(p-value)",
                 "Location of quantiles\np-value",
                 "Outlier\np-value")
write.table(A, sep=",",
            file=paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/paper2/Alltest.",Modname,".csv"),
            row.names = T)


# Mytest
# Dispersion
# Uniformity
# ZeroInf
# Residuals2
# Outliers 


# And test against all variables !!!!!! to detect overdisp 
for (i in c(1:length(VAR))){
  PRED <- VAR[i]
  png(file=paste0(Dir,"diagnostic.test/predictors/",PRED,"_",sizeN,".png"),width=sizeN,height=sizeN*0.9,res=sizeN/resol) # Save it 
  par(mfrow=c(4,5))
  mapply(function(x,y) {
    plotResiduals(x,y[,PRED])
  },x=simulationOutput,y=mod_ALL,SIMPLIFY = F)
  dev.off()
}

## Three test in once + three visual 
TestResid <- lapply(simulationOutput, function(x) testResiduals(x)) # uniformity, disp, and outliers tests
par(mfrow=c(1,1))
plot(simulationOutput$PINPINA) #2 uniformity and manque plot VS resid
lapply(simulationOutput, plot)
plot(simulationOutput$PINPINA$scaledResiduals)

# test with factors (country or plotcat)

######################
## End dharma check ##
######################

## Convert all quantitative variable to qualitative variables. Usefull ????
mod_ALL <-lapply(mod_ALL, function(x) {
  for (i in 1:length(VAR)) {
    x[,paste0(VAR[i],"_cut")] <- cut(x[,VAR[i]],
                                     breaks=seq(range(x[,VAR[i]])[1],
                                                range(x[,VAR[i]])[2],
                                                length.out=10),
                                     labels=FALSE,
                                     include.lowest = T)
  }
  return(x)
})

mod_ALL <- do.call(rbind,mod_ALL)

##############################
## Distribution obs and fit ##
##############################
dplot <- ggplot(mod_ALL) + 
  geom_line(aes(x=log(1+sp.recruitment.ha.weight.2),color="observed"),size=1.5,key_glyph = "abline",stat = "density") +
  geom_line(aes(x=log(1+fit),color="fit"),size=1.5,key_glyph = "abline", stat = "density") +
  #key_glyph = "rect"
  facet_wrap(vars(species),scales = "free")+
  scale_color_manual("Counts",values = c("observed"="black","fit"="blue"),labels=c("observed"="Empirical distribution","fit"="Fitted distribution"))+
  labs(y="Frequency", x="Recruitment counts",las=1)+
  theme_light(base_size = 15)+
  theme(text = element_text(face="bold"),legend.position = c(0.8,0.1),legend.direction ="vertical",
        axis.text.x = element_text(size=11,color="black"),axis.text.y = element_text(size=13,color="black"),
        #axis.text.x = element_blank(),axis.text.y = element_text(size=13,color="black"),
        legend.background=element_rect(fill="white",colour="black",size=0.2),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 18),
        legend.key.size = unit(2, "cm"),
        panel.border = element_rect(colour = "black",fill=NA, size=0.8),
        axis.line = element_line(colour="black"),
        plot.title = element_text(size=18,hjust = 0.5),
        plot.caption = element_text(face="bold.italic"))
dplot
ggsave(filename = paste0("~/SynologyDrive/FUNDIV - NFI - Europe/our-data/species/Recrut.All/",Modname,"/LOG.Distributions_",Modname,".png"),plot = dplot, scale=0.8, width = 20,height = 16, dpi=400) # To modify 

##############################
##     Fit VS residuals     ##
##############################

p <- ggplot(mod_ALL, aes(x=fit, y=diff)) + 
  geom_point(size=0.3,alpha=0.1)+
  geom_abline(slope = 0,intercept = 0,size=0.3)+
  #geom_violin()+
  facet_wrap(vars(species),scales="free")+
  #scale_color_manual("",values = Mycol,labels=Mylabels)+
  geom_jitter(size=0.3,alpha=0.05,position=position_jitter(0.4))+
  #scale_y_continuous(breaks = NULL)+
  #geom_boxplot(width=0.1,outlier.size = 0)+
  labs(y="Residuals", x="Predicted values")+
  theme_light(base_size = 15)+
  theme(text = element_text(face="bold"),legend.position = c(0.8,0.05),legend.direction ="horizontal",
        axis.text.x = element_text(size=11,color="black"),axis.text.y = element_text(size=13,color="black"),
        #axis.text.x = element_blank(),axis.text.y = element_text(size=13,color="black"),
        legend.background=element_rect(fill="white",colour="black",size=0.2),
        panel.border = element_rect(colour = "black", fill=NA, size=0.8),
        axis.line = element_line(colour="black"),
        plot.title = element_text(size=18,hjust = 0.5),
        plot.caption = element_text(face="bold.italic"))
p
ggsave(filename = paste0("~/SynologyDrive/FUNDIV - NFI - Europe/our-data/species/Recrut.All/",Modname,"/Resid.VS.Fit_",Modname,".png"),plot = p, scale=0.8, width = 20,height = 16, dpi=400) # To modify 

##############################
##     Obs VS residuals     ##
############################## 

p2 <- ggplot(mod_ALL, aes(x=sp.recruitment.ha.weight.2, y=diff)) + 
  geom_point(size=0.3,alpha=0.1)+
  geom_abline(slope = 0,intercept = 0,size=0.3)+
  #geom_violin()+
  facet_wrap(vars(species),scales="free")+
  #scale_color_manual("",values = Mycol,labels=Mylabels)+
  geom_jitter(size=0.3,alpha=0.05,position=position_jitter(0.4))+
  #geom_boxplot(width=0.1,outlier.size = 0)+
  labs(y="Residuals", x="Observed values")+
  theme_light(base_size = 15)+
  theme(text = element_text(face="bold"),legend.position = c(0.8,0.05),legend.direction ="horizontal",
        axis.text.x = element_text(size=11,color="black"),axis.text.y = element_text(size=13,color="black"),
        #axis.text.x = element_blank(),axis.text.y = element_text(size=13,color="black"),
        legend.background=element_rect(fill="white",colour="black",size=0.2),
        panel.border = element_rect(colour = "black", fill=NA, size=0.8),
        axis.line = element_line(colour="black"),
        plot.title = element_text(size=18,hjust = 0.5),
        plot.caption = element_text(face="bold.italic"))
p2
ggsave(filename = paste0("~/SynologyDrive/FUNDIV - NFI - Europe/our-data/species/Recrut.All/",Modname,"/Resid.VS.Obs_",Modname,".png"),plot = p2, scale=0.8, width = 20,height = 16, dpi=400) # To modify 

##############################
##     Observed VS Fit      ##
##############################

p3 <- ggplot(mod_ALL, aes(x=sp.recruitment.ha.weight.2, y=fit)) + 
  geom_point(size=0.3,alpha=0.1)+
  geom_abline(slope = 1,intercept = 0,size=0.3)+
  #geom_violin()+
  facet_wrap(vars(species),scales="free")+
  #scale_color_manual("",values = Mycol,labels=Mylabels)+
  geom_jitter(size=0.3,alpha=0.05,position=position_jitter(0.4))+
  #geom_boxplot(width=0.1,outlier.size = 0)+
  labs(y="Fitted values", x="Observed values")+
  theme_light(base_size = 15)+
  theme(text = element_text(face="bold"),legend.position = c(0.8,0.05),legend.direction ="horizontal",
        axis.text.x = element_text(size=11,color="black"),axis.text.y = element_text(size=13,color="black"),
        #axis.text.x = element_blank(),axis.text.y = element_text(size=13,color="black"),
        legend.background=element_rect(fill="white",colour="black",size=0.2),
        panel.border = element_rect(colour = "black", fill=NA, size=0.8),
        axis.line = element_line(colour="black"),
        plot.title = element_text(size=18,hjust = 0.5),
        plot.caption = element_text(face="bold.italic"))
p3
ggsave(filename = paste0("~/SynologyDrive/FUNDIV - NFI - Europe/our-data/species/Recrut.All/",Modname,"/Fitted.VS.Obs_",Modname,".png"),plot = p3, scale=0.8, width = 20,height = 16, dpi=400) # To modify 



######################################################################################################
###                                                                                                ###
### Here are the plots of fit vs pred as a function as another variable (qualitatively transformed)### ## Not usefull anylonger 
###                                                                                                ###
######################################################################################################

colnames(mod_ALL)[colnames(mod_ALL)=="Plotcat"] <- "Plotcat_cut"
## Remplacer Plotcat par la variable que l'on veut + checkl Plotcat 1 et 2

Variation <- c("Plotcat_cut","sp.recruitment.plot.rate.yr.2_cut", "logdbhJ.IMall.plot.mean_cut", 
               "logtreeNbrJ_cut", "mean_spei12_cut", "sqrtBA.O.plot.1_cut", 
               "logBAj.plot.1_cut", "sqrtBAIj.plot.1.mean.J.I_cut", "logsp.mortality.plot.rate.yr_cut")

for (i in 1:length(Variation)){
  mod_ALL[,paste0("SPE_",Variation[i])] <- paste0(mod_ALL$species,mod_ALL[,Variation[i]])
}

# Variation <- Variation[4]
# FitVSobs <- ddply(mod_ALL, ~get(Variation)+species+Plotcat, summarize, mobs=mean(sp.recruitment.plot.rate.yr.2),sdobs=sd(sp.recruitment.plot.rate.yr.2),mfit=mean(fit),sdfit=sd(fit))
# colnames(FitVSobs)[1] <- Variation

Variation
for (i in 1:length(Variation)){
  VarCat <- Variation[i]
  VarSPECat <- paste0("SPE_",VarCat)
  
  if (VarCat!="Plotcat_cut"){
    Mylabels <- seq(1:9)    
    Mycol <- brewer.pal(n = 9, name = "YlOrBr")
  }else  {
    Mylabels <- c("0"="Core",                       
                  "1"="Trailing Edge",               
                  "2"="Leading Edge")
    Mycol <- c("0"="black",
               "1"="red",
               "2"="blue")}
  
  p <- ggplot(mod_ALL, aes(x=get(VarSPECat), y=diff,colour=as.factor(get(VarCat)))) + 
    geom_violin()+
    facet_wrap(vars(species),scales="free")+
    scale_color_manual("",values = Mycol,labels=Mylabels)+
    geom_jitter(aes(colour=as.factor(get(VarCat))),size=0.3,alpha=0.05,position=position_jitter(0.4))+
    geom_boxplot(aes(colour=as.factor(get(VarCat))),width=0.1,outlier.size = 0)+
    labs(y="Obs-fit", x="Species")+
    theme_light(base_size = 15)+
    theme(text = element_text(face="bold"),legend.position = c(0.8,0.05),legend.direction ="horizontal",
          # axis.text.x = element_text(size=11,color="black",angle=90),axis.text.y = element_text(size=13,color="black"),
          axis.text.x = element_blank(),axis.text.y = element_text(size=13,color="black"),
          legend.background=element_rect(fill="white",colour="black",size=0.2),
          panel.border = element_rect(colour = "black", fill=NA, size=0.8),
          axis.line = element_line(colour="black"),
          plot.title = element_text(size=18,hjust = 0.5),
          plot.caption = element_text(face="bold.italic"))
  p
  ggsave(filename = paste0("~/SynologyDrive/FUNDIV - NFI - Europe/our-data/species/Recrut.All/",Modname,"/Obs_Fit_",VarCat,"_",Modname,".png"),plot = p, scale=0.8, width = 20,height = 16, dpi=400) # To modify 
}

p1<-ggplot(myDF_logBAj.plot.1, # All
           aes(x=logBAj.plot.1, y=Predictedresp, color=Plotcat,shape=Plotcat,group=Plotcat))+
  geom_line(size=1)+ #no size
  geom_point(size=1.5,stroke=1.5)+
  facet_grid(vars(Plotcat), vars(species))+
  #geom_linerange(size=1.5)+ # no size
  scale_alpha_manual(values = c(0.3,1),guide=F)+
  scale_color_manual("",values = Mycol,labels=Mylabels)+
  scale_shape_manual("", values=c(1,2,3),labels=Mylabels)+
  labs(title=paste0(CODE),y=paste0("Predicted response"), x=paste0(name))+
  theme_light(base_size = 15)+
  theme(text = element_text(face="bold"),legend.position = c(0.5,0.9),legend.direction ="horizontal",
        axis.text.x = element_text(size=13,color="black"),axis.text.y = element_text(size=13,color="black"),
        legend.background=element_rect(fill="white",colour="black",size=0.2),
        panel.border = element_rect(colour = "black", fill=NA, size=0.8),
        axis.line = element_line(colour="black"),
        plot.title = element_text(size=18,hjust = 0.5),
        plot.caption = element_text(face="bold.italic"))

p1



