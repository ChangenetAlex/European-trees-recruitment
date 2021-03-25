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

Allcode <- c("PINPINA","ABIALB","ACEPSE","ALNGLU","BETPEN","CASSAT","FAGSYL",
             "FRAEXC","PICABI","PINHAL","PINNIG","PINPIN","PINSYL","POPNIG",
             "POPTRE","QUEILE","QUEPET","QUEPYR","QUEROB","QUESUB")
Allseuil <- c(0.8,0.7,0.55,0.6,0.8,0.7,0.7,0.7,0.8,0.7,0.7,0.7,0.8,0.7,0.7,0.8,0.7,0.7,0.8,0.7)
AllseuilC <- c(0.6,0.6,0.5,0.5,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6)


for (i in c(1:13,16:length(Allcode))){
CODE <- Allcode[i]
seuilC <- AllseuilC[i]
seuil <- Allseuil[i]
Dir = paste0("/Users/alexandrechangenet/SynologyDrive/FUNDIV - NFI - Europe/our-data/species/",CODE,"/Recrut.Mortality.2020/")
# Dir = paste0("/home/achangenet/SynologyDrive/FUNDIV - NFI - Europe/our-data/species/",CODE,"/CLIMAP/Recrut.Mortality.2020/")
setwd(Dir)
## Charge the dfplot that is located in the species recrut file 
dfplot <- try(readRDS(paste0("dfplotFinalV2",CODE,seuil,"R.M.rds")),silent = T) #Base de données plot

# Change directory the the model we want to investigate 
Dir <- c(paste0(Dir,Modname,"/")) # Directory 
setwd(Dir)

# Charge the model and the data 
try(assign(Modname,get(load(file = paste0(Modname,".rda")))),silent = F)
rm(x)
data <- get(Modname)$frame

## Response predites in a new dataframe 
data$fit <- as.data.frame(predict(get(Modname),type=c("response")))
data$diff <- data$sp.recruitment.plot.rate.yr.2-data$fit
data$species <- CODE

# Obtain the mean observation for each category and the mean fit (with SD)
realobs=ddply(data, ~Plotcat, summarize, m=mean(sp.recruitment.plot.rate.yr.2),sd=sd(sp.recruitment.plot.rate.yr.2))
realobs$species <- CODE
realfit=ddply(datadf, ~Plotcat, summarize, m=mean(fit),sd=sd(fit))
realfit$species <- CODE
assign(paste0(CODE,"realfit"),realfit)
assign(paste0(CODE,"realobs"),realobs)
assign(paste0(CODE,"datadf"),datadf)
rm("realobs")
rm("realfit")
rm("datadf")
}


realfit <- do.call("rbind",sapply(ls(pattern = "realfit"), function(x) rbind(get(x)), simplify = FALSE))
realobs <- do.call("rbind",sapply(ls(pattern = "realobs"), function(x) rbind(get(x)), simplify = FALSE))
datadf <- do.call("rbind",sapply(ls(pattern = "datadf"), function(x) rbind(get(x)[,c("species","diff","Plotcat")]), simplify = FALSE))

datadf$specat <- paste0(datadf$species,datadf$Plotcat)
ggplot(realfit, aes(species, m, colour=Plotcat))+
  geom_point(data=realfit, aes(x=species, y=m, colour=Plotcat),size=3)+
  # geom_errorbar(data=realfit, aes(ymin=m-sd, ymax=m+sd))+
  geom_point(data=realobs, aes(x=species, y=m, colour=Plotcat),position = position_nudge(x=0.1),alpha=0.5,size=3)+
  #geom_errorbar(data=realobs, aes(ymin=m-sd, ymax=m+sd),position = position_nudge(x=0.1))+
  ylab("Average abundance \n including presences and absences")+
  xlab("Species")
  geom_line(size=1)+ #no size
  geom_point(size=1.5,stroke=1.5)+
  facet_wrap(vars(species),scales = "free")+
  #geom_linerange(size=1.5)+ # no size
  scale_alpha_manual(values = c(0.3,1),guide=F)+
  scale_color_manual("",values = Mycol,labels=Mylabels)+
  scale_shape_manual("", values=c(1,2,3),labels=Mylabels)+
  labs(y=Pred, x=name)+
  theme_light(base_size = 15)+
  theme(text = element_text(face="bold"),legend.position = c(0.5,0.95),legend.direction ="horizontal",
        axis.text.x = element_text(size=13,color="black"),axis.text.y = element_text(size=13,color="black"),
        legend.background=element_rect(fill="white",colour="black",size=0.2),
        panel.border = element_rect(colour = "black", fill=NA, size=0.8),
        axis.line = element_line(colour="black"),
        plot.title = element_text(size=18,hjust = 0.5),
        plot.caption = element_text(face="bold.italic"))
# ggsave(filename = paste0("/home/achangenet/SynologyDrive/FUNDIV - NFI - Europe/our-data/species/Recrut.All/",Pred,"_",name,"_",Modname,".png"),plot = p1, width = 6, height = 6, dpi=300) # To modify 
# ggsave(filename = paste0("/home/achangenet/SynologyDrive/FUNDIV - NFI - Europe/our-data/species/Recrut.All/",Pred,"_",name,"_",Modname,".png"),plot = p1, scale=0.8, width = 20,height = 16, dpi=400) # To modify 
ggsave(filename = paste0("/Users/alexandrechangenet/SynologyDrive/FUNDIV - NFI - Europe/our-data/species/Recrut.All/",Pred,"_",name,"_",Modname,".png"),plot = p1, scale=0.8, width = 20,height = 16, dpi=400) # To modify 
print(p1)

REV <- c("PINPIN","PINNIG","QUEPET","QUEROB","PINPINA","POPTRE","ALNGLU","QUEPYR","QUESUB")
datadf$Plotcat2 <- datadf$Plotcat
datadf[datadf$species%in%REV & datadf$Plotcat2==1,"Plotcat"] <- 2
datadf[datadf$species%in%REV & datadf$Plotcat2==2,"Plotcat"] <- 1
Mycol <- c("black", "red", "blue")
Mylabels <- c("Core","Trailing Edge","Leading Edge")

p <- ggplot(datadf, aes(x=specat, y=diff,colour=Plotcat)) + 
  geom_violin()+
  facet_wrap(vars(species),scales="free")+
scale_color_manual("",values = Mycol,labels=Mylabels)+
geom_jitter(aes(colour=Plotcat),size=0.3,alpha=0.05,position=position_jitter(0.4))+
  geom_boxplot(aes(colour=Plotcat),width=0.1,outlier.size = 0)+
  labs(y="Obs-fit", x="Species")+
theme_light(base_size = 15)+
  theme(text = element_text(face="bold"),legend.position = c(0.8,0.05),legend.direction ="horizontal",
        axis.text.x = element_text(size=13,color="black"),axis.text.y = element_text(size=13,color="black"),
        legend.background=element_rect(fill="white",colour="black",size=0.2),
        panel.border = element_rect(colour = "black", fill=NA, size=0.8),
        axis.line = element_line(colour="black"),
        plot.title = element_text(size=18,hjust = 0.5),
        plot.caption = element_text(face="bold.italic"))
p
ggsave(filename = paste0("/Users/alexandrechangenet/SynologyDrive/FUNDIV - NFI - Europe/our-data/species/Recrut.All/",Modname,"/Obs_Fit",Modname,".png"),plot = p, scale=0.8, width = 20,height = 16, dpi=400) # To modify 





## Difference entre chaque observation et les fit. 

pt.jitter <- function(x, y, bg="gray", col="darkgray", alpha=1, data=NULL){
  points(y~I(as.numeric(as.factor(x))+runif(length(x),-0.15, 0.15)),
         pch=21, bg=adjustcolor(bg, alpha.f = alpha), col=col,cex=0.5)
}

beanplot(diff~specat,
         data=datadf,
         log="",
         xlab="", col.axis="darkgray",
         what=c(0,1,1,0),
         cex.lab=1.5,
         # cutmin=0,
         cex.axis=1,
         # ylim=c(0,300),
         border="gray",
         #border=colors[c("Bp", "Pp", "Qi", "Qr")],
         #names=c("BRAN", "CDS1","CDS2","CDS3", "MTF1", "MTF2","MTF3"),
         ylab = "Chao1 index",
         col=as.list("white"),
         main="Fungi")
pt.jitter(datadf$specat, datadf$diff,col =c("0"="yellow","1"="blue","2"="red")[datadf$Plotcat],
          bg = c("0"="yellow","1"="blue","2"="red")[datadf$Plotcat], alpha = 0.6)





newdata = unique(data[,c("logtreeNbrJ","Plotcat")])
temp = predict(get(Modname), newdata, se.fit=TRUE, type="response")



# Ici on regroupe le nombre moyen de recruited par plotcat et niveau de logtreeNbrJ 
real=ddply(datadf, ~Plotcat+logtreeNbrJ, summarize, m=mean(sp.recruitment.plot.rate.yr.2))
ggplot(datadf, aes(logtreeNbrJ, fit, colour=Plotcat))+geom_point()+
  geom_errorbar(aes(ymin=fit-se.fit, ymax=fit+se.fit))+
  geom_point(data=real, aes(x=logtreeNbrJ, y=m) )+
  ylab("Average abundance \n including presences and absences")+
  xlab("Species")

real=ddply(datadf, ~Plotcat+logBAj.plot.1, summarize, m=mean(sp.recruitment.plot.rate.yr.2))
ggplot(datadf, aes(logBAj.plot.1, fit, colour=Plotcat))+geom_point()+
  geom_errorbar(aes(ymin=fit-se.fit, ymax=fit+se.fit))+
  geom_point(data=real, aes(x=logBAj.plot.1, y=m) )+
  ylab("Average abundance \n including presences and absences")+
  xlab("Species")









?ddply

chisq.test(Test,data$sp.recruitment.plot.rate.yr.2)
Test <- predict(get(Modname),type=c("response"))
Test2 <- predict(SCALEDRrate.yr.QualiF.TF.noY.nbinom1)
data$sp.recruitment.plot.rate.yr.2
plot(Test~data$sp.recruitment.plot.rate.yr.2)



plot(get(Modname)$fitted,data$sp.recruitment.plot.rate.yr.2)
SCALEDRrate.yr.QualiF.TF.noY.nbinom1.V2$obj

owls_nb1_simres <- simulateResiduals(get(Modname))
plot(owls_nb1_simres)
Anova(get(Modname),type="III")



simulationOutput <- simulateResiduals(get(Modname), n = 1000)
simulationOutput$scaledResiduals
plot(simulationOutput)
plot(simulationOutput$scaledResiduals)


plotResiduals(get(Modname)$frame$logdbhJ.IMall.plot.mean, simulationOutput$scaledResiduals)
plotResiduals(get(Modname)$frame$treeNbrJ, simulationOutput$scaledResiduals)
plotResiduals(get(Modname)$frame$mean_spei12, simulationOutput$scaledResiduals)
plotResiduals(get(Modname)$frame$sqrtBA.O.plot.1, simulationOutput$scaledResiduals)
plotResiduals(get(Modname)$frame$sqrtBAIj.plot.1.mean.J.I, simulationOutput$scaledResiduals)
plotResiduals(get(Modname)$frame$sp.mortality.plot.rate.yr, simulationOutput$scaledResiduals)


testResiduals(simulationOutput)
testUniformity(simulationOutput)
testZeroinflation(simulationOutput)
testDispersion(simulationOutput)


simulationOutput = recalculateResiduals(simulationOutput, group = testData$group)



library(effects)
plot(allEffects(get(Modname)))

system.time(owls_nb1_d1 <- drop1(get(Modname),test="Chisq"))
print(owls_nb1_d1)

owls_nb1_dredge <- MuMIn::dredge(get(Modname))
op <- par(mar=c(2,5,14,3))
plot(owls_nb1_dredge)
par(op)
# Model averaging
model.avg(owls_nb1_dredge)


### Broom 
library(dplyr)
t1 <- broom.mixed::tidy(get(Modname), conf.int = TRUE)
  #if (packageVersion("dotwhisker")>"0.4.1") {
    ## to get this version (which fixes various dotwhisker problems)
    ## use devtools::install_github("bbolker/broom.mixed") or
    ## wait for pull request acceptance/submission to CRAN/etc.
dwplot(get(Modname))+geom_vline(xintercept=0,lty=2)
# depends on the version
SCALEDRrate.yr.QualiF.TF.noY.nbinom1$coefficients <- TRUE ## hack!
dwplot(SCALEDRrate.yr.QualiF.TF.noY.nbinom1,by_2sd=FALSE)+geom_vline(xintercept=0,lty=2)












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










ggplot(myDF_logBAj.plot.1, aes(x=logBAj.plot.1, y=Predictedresp, color=Plotcat))+
  geom_line

ggplot(myDF_logBAj.plot.1, aes(x=logBAj.plot.1, y=Predictedresp, color=Plotcat))+
  facet_wrap(vars(species))+
  geom_line

ggplot(myDF_logBAj.plot.1, aes(x=logBAj.plot.1, y=Predictedresp, color=Plotcat))+
  facet_wrap(vars(species))+
  geom_line







rownames(myDF_logBAj.plot.1)

Esp <- str_sub(Esp, 3)



## Charge the dfplot that is located in the species recrut file 
library(stringr)# Also load dfplo2 or three

# Change directory the the model we want to investigate 
Dir <- c(paste0(Dir,Modname,"/")) # Directory 
setwd(Dir)



ggplot(ssd, aes(x=absence, color=Plotcat))+
  geom_density(adjust=4)+
  geom_point(data=real, aes(x=absence, y=1, color=Plotcat), size=2)+
  xlab("Probability that salamanders are not observed")+ylab(NULL)


## Representation 












## Rep like in the paper with the code given

library(plyr)
data("Salamanders")
Salamanders
data <- get(Modname)$frame
sims=simulate(get(Modname), seed = 1, nsim = 1000)

colnames(data)

simdatlist=lapply(sims, function(count){cbind(count, data)
})
simdatsums=lapply(simdatlist, function(x){
  ddply(x, ~Plotcat+treeNbrJ, summarize,
        absence=mean(count==0),
        mu=mean(count))
})
ssd=do.call(rbind, simdatsums)
real = ddply(data, ~Plotcat+treeNbrJ, summarize,
             absence=mean(sp.recruitment.plot.rate.yr.2==0),
             mu=mean(sp.recruitment.plot.rate.yr.2))
ggplot(ssd, aes(x=absence, color=Plotcat))+
  geom_density(adjust=4)+
  geom_point(data=real, aes(x=absence, y=1, color=Plotcat), size=2)+
  xlab("Probability that salamanders are not observed")+ylab(NULL)

ggplot(ssd, aes(x=mu, color=Plotcat))+
  geom_density(adjust=4)+
  geom_point(data=real, aes(x=mu, y=.5, color=Plotcat), size=2)+
  xlab("Abundance including presences and absences")+ylab(NULL)



###
Test <- as.data.frame(predict(get(Modname),type=c("response"),se=T))
data$sp.recruitment.plot.rate.yr.2
datadf <- cbind(data,Test)

datadf$sp.recruitment.plot.rate.yr.2
data$fit




cor.test(datadf$sp.recruitment.plot.rate.yr.2,datadf$fit)
datadf$sqrtBA.O.plot.1

real.count = ddply(datadf, ~Plotcat+treeNbrJ, summarize, m=median(sp.recruitment.plot.rate.yr.2), mu=mean(sp.recruitment.plot.rate.yr.2))
pred.count = ddply(datadf, ~Plotcat+treeNbrJ, summarize, m=median(fit), mu=mean(fit))

ggplot(datadf, aes(x=Plotcat, y=sp.recruitment.plot.rate.yr.2))+
  geom_point(shape=15, size=2)+
  #geom_errorbar(aes(ymin=ucount.low, ymax=ucount.high))+
  geom_point(data=datadf, aes(x=Plotcat, y=fit), shape=0, size=2)
#geom_point(data=real.count, aes(x=spp, y=mu, colour=mined), shape=5, size=2)+
#ylab("Abundance \n including presences and absences")+
#xlab("Species")








datadf$d <- datadf$sp.recruitment.plot.rate.yr.2-datadf$fit

plot(datadf$d~datadf$sp.recruitment.plot.rate.yr.2)
test=ddply(datadf, ~sp.recruitment.plot.rate.yr.2, summarize, m=mean(d))
plot(test$m~test$sp.recruitment.plot.rate.yr.2)

test=ddply(datadf, ~mean_spei12, summarize, m=mean(fit))
plot(test$m~test$mean_spei12)





