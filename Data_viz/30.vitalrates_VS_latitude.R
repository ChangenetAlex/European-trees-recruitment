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

Modname <- "SCALEDRrate.yr.QualiF.TF.noY.nbinom1.V3" ## Name of the model 
Modname <- "rW.yO.sumO.mRY1"

Allcode <- c("ABIALB","ACEPSE","ALNGLU","BETPEN","CASSAT","FAGSYL",
             "FRAEXC","PICABI","PINHAL","PINNIG","PINPIN","PINPINA","PINSYL","POPNIG",
             "POPTRE","QUEILE","QUEPET","QUEPYR","QUEROB","QUESUB")
Allseuil <- c(0.7,0.55,0.6,0.8,0.7,0.7,0.7,0.8,0.7,0.7,0.7,0.8,0.8,0.7,0.7,0.8,0.7,0.7,0.8,0.7)
AllseuilC <- c(0.6,0.5,0.5,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6)
mod_ALL <- list()
dfplot <- list()

for (i in c(1:13,16:length(Allcode))){
  # for (i in c(1:3)){
  CODE <- Allcode[i]
  seuilC <- AllseuilC[i]
  seuil <- Allseuil[i]
  Dir = paste0("~/SynologyDrive/FUNDIV - NFI - Europe/our-data/species/",CODE,"/CLIMAP/Recrut.Mortality.2020/")
  setwd(Dir)
  # Charge the dfplot for simulation and residuals test dispersion etc... 
  dfplot[[CODE]] <- readRDS(paste0("dfplotFinalV2",CODE,seuil,"R.M.rds"))
  
  #Load the desired model of the desired species 
  Dir <- c(paste0(Dir,Modname,"/")) # Directory 
  setwd(Dir)
  mod_ALL[[CODE]] <- get(load(file = paste0(Modname,".rda")))
  rm(x)
}


mod_ALL <- lapply(mod_ALL, function(x) {
  x[[5]][,"fit"] <- predict(x,x[[5]],type=c("response"))
  x[[5]][,"diff"] <- x[[5]][,"sp.recruitment.plot.rate.yr.2"]-x[[5]][,"fit"]
  if (length(which(colnames(x[[5]])=="country")!=0)){
    x[[5]] <- x[[5]][,-which(colnames(x[[5]])=="country")]
  }
  ;x[[5]]})

## Convert it to a df
mod_ALL <- mapply(function(x,y) {
  x[,"species"] <- y
  return(x)
},
x=mod_ALL,y=names(mod_ALL),SIMPLIFY = F)

mod_ALL <- mapply(function (x,xx){
  xx[,"latitude"] <- x[which(rownames(x)%in%rownames(xx)),c("latitude")]
  return(xx)
},xx=mod_ALL,x=dfplot,SIMPLIFY = F)

unit <- 0.5 #

mod_ALL.sum <- lapply(mod_ALL, function(x){
  S <- seq(floor(min(x[,"latitude"])),ceiling(max(x[,"latitude"])),unit)
  print(S[length(S)])
  S[length(S)] <- ceiling(max(x[,"latitude"]))
  print(S[length(S)])
  x[,"latitude_c"] <- cut(x[,"latitude"],S,right=F)
  y <- ddply(x, ~latitude_c, .drop = F,summarize,  #drop = F keep all categorie. 
             mR=median(sp.recruitment.plot.rate.yr.2),
             muR=mean(sp.recruitment.plot.rate.yr.2),
             sdR=sd(sp.recruitment.plot.rate.yr.2),
             mFit=median(fit),
             muFit=mean(fit),
             sdFit=sd(fit),
             mM=median(exp(logsp.mortality.plot.rate.yr)-1),
             muM=mean(exp(logsp.mortality.plot.rate.yr)-1),
             sdM=sd(exp(logsp.mortality.plot.rate.yr)-1),
             mG=median(sqrtBAIj.plot.1.mean.J.I*sqrtBAIj.plot.1.mean.J.I),
             muG=mean(sqrtBAIj.plot.1.mean.J.I*sqrtBAIj.plot.1.mean.J.I),
             sdG=sd(sqrtBAIj.plot.1.mean.J.I*sqrtBAIj.plot.1.mean.J.I))
  y[,"latitude"] <- seq(floor(min(x[,"latitude"])),ceiling(max(x[,"latitude"])),unit)[-1]
  y[,"species"] <- x[1,"species"]
  ;y})


mod_ALL.2 <- do.call(rbind,mod_ALL) # F
mod_ALL.sum.2 <- do.call(rbind,mod_ALL.sum)
#mod_ALL <- split(x = mod_ALL,f = mod_ALL$species)
#mod_ALL.sum.2[mod_ALL.sum.2$muR=="0","muR"] <- NA

# Mortality and growth
dplot <- ggplot(mod_ALL.sum.2)+
  geom_smooth(aes(x=latitude, y=log(muR+1),color="Recrut"),size=1.5,se=F) +
  geom_smooth(aes(x=latitude, y=log(muG+1),color="Growth"),size=1.5,se=F) +
  geom_smooth(aes(x=latitude, y=log(muM+1),color="Mortality"),size=1.5,se=F) +
  
  geom_point(aes(x=latitude, y=log(muG+1),color="Growth"),shape=19, size=2)+
  geom_point(aes(x=latitude, y=log(muM+1),color="Mortality"),shape=19, size=2)+
  geom_point(aes(x=latitude, y=log(muR+1),color="Recrut"),shape=19, size=2)+
  
  # geom_ribbon(aes(x=latitude, y=log(1+muR), ymin=log(1+muR)-log(sdR), ymax=log(1+muR)+log(sdR), color="Recrut"),alpha=0.1, linetype=2)+
  # geom_ribbon(aes(x=latitude, y=log(1+muM), ymin=log(1+muM)-log(sdM), ymax=log(1+muM)+log(sdM), color="Mortality"),alpha=0.1, linetype=2)+
  # geom_ribbon(aes(x=latitude, y=log(1+muG), ymin=log(1+muG)-log(sdG), ymax=log(1+muG)+log(sdG), color="Growth"),alpha=0.1, linetype=2)+
  
  facet_wrap(vars(species),scales = "free")+
  scale_color_manual("Demographic traits",values = c("Mortality"="black","Growth"="blue","Recrut"="red"),labels=c("Mortality"="Mortality","Growth"="Growth","Recrut"="Recruitment"))+
  scale_fill_manual("Demographic traits",values = c("Mortality"="black","Growth"="blue","Recrut"="red"),labels=c("Mortality"="Mortality","Growth"="Growth","Recrut"="Recruitment"))+
  labs(y="annual rate (log)", x="Latitude")+
  theme_light(base_size = 15)+
  theme(text = element_text(face="bold"),legend.position = c(0.8,0.1),legend.direction ="vertical",
        axis.text.x = element_text(size=11,color="black"),axis.text.y = element_text(size=13,color="black"),
        #axis.text.x = element_blank(),axis.text.y = element_text(size=13,color="black"),
        legend.background=element_rect(fill="white",colour="black",size=0.5),
        legend.key.size = unit(3.5, "lines"),
        legend.key.width = unit(8,"lines"),
        legend.text = element_text(size=14,color="black"),
        panel.border = element_rect(colour = "black",fill=NA, size=0.8),
        axis.line = element_line(colour="black"),
        plot.title = element_text(size=18,hjust = 0.5),
        plot.caption = element_text(face="bold.italic"))
dplot
ggsave(filename = paste0("~/SynologyDrive/FUNDIV - NFI - Europe/our-data/species/Recrut.All/",Modname,"/G.M.R_VS_latitude_0.5.log",Modname,".png"),plot = dplot, scale=0.8, width = 20,height = 16, dpi=400) # To modify 



unit <- 0.2 #


##### Idem SPEI !!!!!
mod_ALL.sum.SPEI <- lapply(mod_ALL, function(x){
  S <- seq(floor(min(x[,"mean_spei12"])),ceiling(max(x[,"mean_spei12"])),unit)
  print(S[length(S)])
  S[length(S)] <- ceiling(max(x[,"mean_spei12"]))
  print(S[length(S)])
  x[,"mean_spei12_c"] <- cut(x[,"mean_spei12"],S,right=F)
  y <- ddply(x, ~mean_spei12_c, .drop = F,summarize,  #drop = F keep all categorie. 
             mR=median(sp.recruitment.plot.rate.yr.2),
             muR=mean(sp.recruitment.plot.rate.yr.2),
             sdR=sd(sp.recruitment.plot.rate.yr.2),
             mFit=median(fit),
             muFit=mean(fit),
             sdFit=sd(fit),
             mM=median(exp(logsp.mortality.plot.rate.yr)-1),
             muM=mean(exp(logsp.mortality.plot.rate.yr)-1),
             sdM=sd(exp(logsp.mortality.plot.rate.yr)-1),
             mG=median(sqrtBAIj.plot.1.mean.J.I*sqrtBAIj.plot.1.mean.J.I),
             muG=mean(sqrtBAIj.plot.1.mean.J.I*sqrtBAIj.plot.1.mean.J.I),
             sdG=sd(sqrtBAIj.plot.1.mean.J.I*sqrtBAIj.plot.1.mean.J.I))
  y[,"mean_spei12"] <- seq(floor(min(x[,"mean_spei12"])),ceiling(max(x[,"mean_spei12"])),unit)[-1]
  y[,"species"] <- x[1,"species"]
  ;y})


mod_ALL.sum.SPEI.2 <- do.call(rbind,mod_ALL.sum.SPEI) # F

# Figure
# Mortality and growth
dplot <- ggplot(mod_ALL.sum.SPEI.2)+
  geom_smooth(aes(x=mean_spei12, y=log(muR+1),color="Recrut"),size=1.5,se=F) +
  geom_smooth(aes(x=mean_spei12, y=log(muG+1),color="Growth"),size=1.5,se=F) +
  geom_smooth(aes(x=mean_spei12, y=log(muM+1),color="Mortality"),size=1.5,se=F) +
  
  geom_point(aes(x=mean_spei12, y=log(muG+1),color="Growth"),shape=19, size=2)+
  geom_point(aes(x=mean_spei12, y=log(muM+1),color="Mortality"),shape=19, size=2)+
  geom_point(aes(x=mean_spei12, y=log(muR+1),color="Recrut"),shape=19, size=2)+
  
  # geom_ribbon(aes(x=latitude, y=log(1+muR), ymin=log(1+muR)-log(sdR), ymax=log(1+muR)+log(sdR), color="Recrut"),alpha=0.1, linetype=2)+
  # geom_ribbon(aes(x=latitude, y=log(1+muM), ymin=log(1+muM)-log(sdM), ymax=log(1+muM)+log(sdM), color="Mortality"),alpha=0.1, linetype=2)+
  # geom_ribbon(aes(x=latitude, y=log(1+muG), ymin=log(1+muG)-log(sdG), ymax=log(1+muG)+log(sdG), color="Growth"),alpha=0.1, linetype=2)+
  
  facet_wrap(vars(species),scales = "free")+
  scale_color_manual("Demographic traits",values = c("Mortality"="black","Growth"="blue","Recrut"="red"),labels=c("Mortality"="Mortality","Growth"="Growth","Recrut"="Recruitment"))+
  scale_fill_manual("Demographic traits",values = c("Mortality"="black","Growth"="blue","Recrut"="red"),labels=c("Mortality"="Mortality","Growth"="Growth","Recrut"="Recruitment"))+
  labs(y="annual rate (log)", x="mean_spei12")+
  theme_light(base_size = 15)+
  theme(text = element_text(face="bold"),legend.position = c(0.8,0.1),legend.direction ="vertical",
        axis.text.x = element_text(size=11,color="black"),axis.text.y = element_text(size=13,color="black"),
        #axis.text.x = element_blank(),axis.text.y = element_text(size=13,color="black"),
        legend.background=element_rect(fill="white",colour="black",size=0.5),
        legend.key.size = unit(3.5, "lines"),
        legend.key.width = unit(8,"lines"),
        legend.text = element_text(size=14,color="black"),
        panel.border = element_rect(colour = "black",fill=NA, size=0.8),
        axis.line = element_line(colour="black"),
        plot.title = element_text(size=18,hjust = 0.5),
        plot.caption = element_text(face="bold.italic"))
dplot





# This figure to obtain the density of fit and predict. 

dplot <- ggplot(mod_ALL) + 
  geom_line(aes(x=sp.recruitment.plot.rate.yr.2,color="observed"),size=1.5,key_glyph = "abline",stat = "density") +
  geom_line(aes(x=fit,color="fit"),size=1.5,key_glyph = "abline",stat = "density") +
  #key_glyph = "rect"
  facet_wrap(vars(species),scales = "free")+
  scale_color_manual("Counts",values = c("observed"="black","fit"="blue"),labels=c("observed"="Empirical distribution","fit"="Fitted distribution"))+
  labs(y="Probability", x="Recruitment counts")+
  theme_light(base_size = 15)+
  theme(text = element_text(face="bold"),legend.position = c(0.8,0.05),legend.direction ="vertical",
        #axis.text.x = element_text(size=11,color="black",angle=90),axis.text.y = element_text(size=13,color="black"),
        axis.text.x = element_blank(),axis.text.y = element_text(size=13,color="black"),
        legend.background=element_rect(fill="white",colour="black",size=0.2),
        panel.border = element_rect(colour = "black",fill=NA, size=0.8),
        axis.line = element_line(colour="black"),
        plot.title = element_text(size=18,hjust = 0.5),
        plot.caption = element_text(face="bold.italic"))
dplot
ggsave(filename = paste0("~/SynologyDrive/FUNDIV - NFI - Europe/our-data/species/Recrut.All/",Modname,"/Distributions_",Modname,".png"),plot = dplot, scale=0.8, width = 20,height = 16, dpi=400) # To modify 



# Fitted recruit VS obs recruit
dplot  <- ggplot(mod_ALL.sum)+
  geom_line(aes(x=latitude, y=muFit,color="fit"),size=1.5) +
  #geom_ribbon(aes(x=latitude, y=muFit, ymin=muFit-sdFit, ymax=muFit+sdFit, color="fit", fill="fit"),alpha=0.1, linetype=2)+
  #geom_point(aes(x=latitude, y=muFit,color="fit"),shape=15, size=2)+
  geom_line(aes(x=latitude, y=muR,color="observed"),size=1.5) +
  #geom_ribbon(aes(x=latitude, y=muR, ymin=muR-sdR, ymax=muR+sdR, color="observed", fill="observed"),alpha=0.1, linetype=2)+
  #geom_point(aes(x=latitude, y=muR,color="observed"),shape=15, size=2)+
  facet_wrap(vars(species),scales = "free")+
  scale_color_manual("Counts",values = c("observed"="black","fit"="blue"),labels=c("observed"="Empirical distribution","fit"="Fitted distribution"))+
  scale_fill_manual("Counts",values = c("observed"="black","fit"="blue"),labels=c("observed"="Empirical distribution","fit"="Fitted distribution"))+
  labs(y="mu", x="latitude")+
  theme_light(base_size = 15)+
  theme(text = element_text(face="bold"),legend.position = c(0.8,0.05),legend.direction ="vertical",
        axis.text.x = element_text(size=11,color="black"),axis.text.y = element_text(size=13,color="black"),
        #axis.text.x = element_blank(),axis.text.y = element_text(size=13,color="black"),
        legend.background=element_rect(fill="white",colour="black",size=0.2),
        panel.border = element_rect(colour = "black",fill=NA, size=0.8),
        axis.line = element_line(colour="black"),
        plot.title = element_text(size=18,hjust = 0.5),
        plot.caption = element_text(face="bold.italic"))
dplot
ggsave(filename = paste0("~/SynologyDrive/FUNDIV - NFI - Europe/our-data/species/Recrut.All/",Modname,"/Recrut_VS_latitude_0.5",Modname,".png"),plot = dplot, scale=0.8, width = 20,height = 16, dpi=400) # To modify 


# Mortality and growth
dplot <- ggplot(mod_ALL.sum)+
  geom_line(aes(x=latitude, y=log(muR+1),color="Recrut"),size=1.5) +
  geom_line(aes(x=latitude, y=muG,color="Growth"),size=1.5) +
  #geom_ribbon(aes(x=latitude, y=muG, ymin=muG-sdG, ymax=muG+sdG, color="Growth", fill="Growth"),alpha=0.1, linetype=2)+
  geom_line(aes(x=latitude, y=muM,color="Mortality"),size=1.5) +
  #geom_ribbon(aes(x=latitude, y=muM, ymin=muM-sdM, ymax=muM+sdM, color="Mortality", fill="Mortality"),alpha=0.1, linetype=2)+
  #geom_line(aes(x=latitude, y=sqrt(muR),color="Recruitment"),size=1.5) +
  #geom_point(aes(x=latitude, y=sqrt(muR),color="Recruitment"),shape=15, size=2)+
  facet_wrap(vars(species),scales = "free_x")+
  scale_color_manual("Demographic traits",values = c("Mortality"="black","Growth"="blue","Recrut"="red"),labels=c("Mortality"="Mortality","Growth"="Growth","Recrut"="Recruitment"))+
  scale_fill_manual("Demographic traits",values = c("Mortality"="black","Growth"="blue","Recrut"="red"),labels=c("Mortality"="Mortality","Growth"="Growth","Recrut"="Recruitment"))+
  labs(y="annual rate (log)", x="Latitude")+
  theme_light(base_size = 15)+
  theme(text = element_text(face="bold"),legend.position = c(0.8,0.1),legend.direction ="vertical",
        axis.text.x = element_text(size=11,color="black"),axis.text.y = element_text(size=13,color="black"),
        #axis.text.x = element_blank(),axis.text.y = element_text(size=13,color="black"),
        legend.background=element_rect(fill="white",colour="black",size=0.5),
        legend.key.size = unit(3.5, "lines"),
        legend.key.width = unit(8,"lines"),
        legend.text = element_text(size=14,color="black"),
        panel.border = element_rect(colour = "black",fill=NA, size=0.8),
        axis.line = element_line(colour="black"),
        plot.title = element_text(size=18,hjust = 0.5),
        plot.caption = element_text(face="bold.italic"))
dplot
ggsave(filename = paste0("~/SynologyDrive/FUNDIV - NFI - Europe/our-data/species/Recrut.All/",Modname,"/G.M.R_VS_latitude_0.5",Modname,".png"),plot = dplot, scale=0.8, width = 20,height = 16, dpi=400) # To modify 



# Mortality growth and recruit: 
dplot <- ggplot(mod_ALL.sum)+
  geom_line(aes(x=latitude, y=muG,color="Growth"),size=1.5) +
  #geom_ribbon(aes(x=latitude, y=muG, ymin=muG-sdG, ymax=muG+sdG, color="Growth", fill="Growth"),alpha=0.1, linetype=2)+
  geom_line(aes(x=latitude, y=muM,color="Mortality"),size=1.5) +
  #geom_ribbon(aes(x=latitude, y=muM, ymin=muM-sdM, ymax=muM+sdM, color="Mortality", fill="Mortality"),alpha=0.1, linetype=2)+
  geom_line(aes(x=latitude, y=log(muR),color="Recruitment"),size=1.5) +
  #geom_point(aes(x=latitude, y=sqrt(muR),color="Recruitment"),shape=15, size=2)+
  facet_wrap(vars(species),scales = "free")+
  scale_color_manual("Demographic traits",values = c("Mortality"="black","Growth"="blue","Recruitment"="red"),labels=c("Mortality"="Mortality","Growth"="Growth","Recruitment"="Recruitment"))+
  scale_fill_manual("Demographic traits",values = c("Mortality"="black","Growth"="blue","Recruitment"="red"),labels=c("Mortality"="Mortality","Growth"="Growth","Recruitment"="Recruitment"))+
  labs(y="mu", x="latitude")+
  theme_light(base_size = 15)+
  theme(text = element_text(face="bold"),legend.position = c(0.8,0.05),legend.direction ="vertical",
        axis.text.x = element_text(size=11,color="black"),axis.text.y = element_text(size=13,color="black"),
        #axis.text.x = element_blank(),axis.text.y = element_text(size=13,color="black"),
        legend.background=element_rect(fill="white",colour="black",size=0.2),
        legend.key.size = unit(3, "lines"),
        panel.border = element_rect(colour = "black",fill=NA, size=0.8),
        axis.line = element_line(colour="black"),
        plot.title = element_text(size=18,hjust = 0.5),
        plot.caption = element_text(face="bold.italic"))
dplot
ggsave(filename = paste0("~/SynologyDrive/FUNDIV - NFI - Europe/our-data/species/Recrut.All/",Modname,"/G.M.R_VS_latitude_0.5",Modname,".png"),plot = dplot, scale=0.8, width = 20,height = 16, dpi=400) # To modify 


### Ajout des ecarts types








## shape my mean data for ggplots (abs and relative)
ind.rel.imp.clim.bio.long <- melt(t.mean.rel, id = "latitude", measure = c(EffectAll))
names(ind.rel.imp.clim.bio.long) <- c("latitude", "variable", "mean")
ind.abs.imp.clim.bio.long <- melt(t.mean.abs, id = "latitude", measure = c(EffectAll))
names(ind.abs.imp.clim.bio.long) <- c("latitude", "variable", "mean")
## error bands for ggplots (se)
t.se.rel <- t.sd.rel[,(ncol(t.sd.rel)-length(EffectAll)):ncol(t.sd.rel)]
colnames(t.se.rel) <- sub("_se","", colnames(t.se.rel), ignore.case = FALSE,fixed = T) # remove the SE in the name to match all the variable
ind.rel.imp.clim.bio.long.se <- melt(t.se.rel, id = "latitude", measure = c(EffectAll))
names(ind.rel.imp.clim.bio.long.se) <- c("latitude", "variable", "se")
t.se.abs <- t.sd.abs[,(ncol(t.sd.abs)-length(EffectAll)):ncol(t.sd.abs)] #same thing for relative 
colnames(t.se.abs) <- sub("_se","", colnames(t.se.abs), ignore.case = FALSE,fixed = T)
ind.abs.imp.clim.bio.long.se <- melt(t.se.abs, id = "latitude", measure = c(EffectAll))
names(ind.abs.imp.clim.bio.long.se) <- c("latitude", "variable", "se")

# Merge the mean and the se dataframes
clim.bio.rel <- merge(ind.rel.imp.clim.bio.long, ind.rel.imp.clim.bio.long.se, by=c("latitude", "variable"))
clim.bio.rel$lwr <- clim.bio.rel$mean-(1.96*clim.bio.rel$se) # lower intervalle
clim.bio.rel$upr <- clim.bio.rel$mean+(1.96*clim.bio.rel$se) # higher intervalle 
clim.bio.rel$lwr <- ifelse(clim.bio.rel$lwr<0,0,clim.bio.rel$lwr) 
clim.bio.rel$upr <- ifelse(clim.bio.rel$upr>1,1,clim.bio.rel$upr)

clim.bio.abs <- merge(ind.abs.imp.clim.bio.long, ind.abs.imp.clim.bio.long.se, by=c("latitude", "variable"))
clim.bio.abs$lwr <- clim.bio.abs$mean-(1.96*clim.bio.abs$se)
clim.bio.abs$upr <- clim.bio.abs$mean+(1.96*clim.bio.abs$se) 
#clim.bio.abs$lwr <- ifelse(clim.bio.abs$lwr<0,0,clim.bio.abs$lwr) # NON SENSE TO SCALE IN ABSOLUTE
#clim.bio.abs$upr <- ifelse(clim.bio.abs$upr>1,1,clim.bio.abs$upr) # NON SENSE TO SCALE IN ABSOLUTE

# Remove infinite values by NA in order to plot correctly 
l.Inf1 <- nrow(clim.bio.abs[sapply(clim.bio.abs[,3:6], function(x) is.infinite(x)),3:6])
l.Inf2 <- nrow(clim.bio.rel[sapply(clim.bio.rel[,3:6], function(x) is.infinite(x)),3:6])
if (l.Inf1!=0){
  message(paste0("There is ",l.Inf1," Infinite values in the clim.bio.abs database, it will be replaced by NA values in order to be plotted"))
  clim.bio.abs[,3:6] <- sapply(clim.bio.abs[,3:6], function(x) {x[is.infinite(x)] <- NA; return(x)})}
if (l.Inf2!=0){
  message(paste0("There is ",l.Inf2," Infinite values in the clim.bio.rel database, it will be replaced by NA values in order to be plotted"))
  clim.bio.rel[,3:6] <- sapply(clim.bio.rel[,3:6], function(x) {x[is.infinite(x)] <- NA; return(x)})}

## Save it all 
save(clim.bio.rel, file=paste0("clim_bio_rel_",deparse(substitute(x)),".RData"))
save(clim.bio.abs, file=paste0("clim_bio_abs_",deparse(substitute(x)),".RData"))

