### Comparaison méthodes pour le bootstrap ###


###Use first the standard error method within the predict function 


Predictedresp <- predict(
  get(Modname),
  newdata = myDF,
  newparams = NULL,
  se.fit = T,
  allow.new.levels = FALSE,
  type = c("response"),
  na.action = na.pass,
  debug = FALSE
)
## Predicted mu (count part)

Predictedmu <- predict(
  get(Modname),
  newdata = myDF,
  newparams = NULL,
  se.fit = T,
  allow.new.levels = FALSE,
  type = c("conditional"),
  na.action = na.pass,
  debug = FALSE
)

## Predicted probability to observe a zero event (= proba of a non recruit event) 
Predictedzprob <- predict(
  get(Modname),
  newdata = myDF,
  newparams = NULL,
  se.fit = T,
  allow.new.levels = FALSE,
  type = c("zprob"),
  na.action = na.pass,
  debug = FALSE
)


## Predicted mu (count part)

Predictedlink <- predict(
  get(Modname),
  newdata = myDF,
  newparams = NULL,
  se.fit = T,
  allow.new.levels = FALSE,
  type = c("link"),
  na.action = na.pass,
  debug = FALSE
)

## Predicted probability to observe a zero event (= proba of a non recruit event) 
Predictedzlink <- predict(
  get(Modname),
  newdata = myDF,
  newparams = NULL,
  se.fit = T,
  allow.new.levels = FALSE,
  type = c("zlink"),
  na.action = na.pass,
  debug = FALSE
)



myDF$Predictedresp <- Predictedresp$fit
myDF$SEresp <- Predictedresp$se.fit  
myDF$Predictedmu <- Predictedmu$fit
myDF$SEmu <- Predictedmu$se.fit
myDF$Predictedzprob <- Predictedzprob$fit
myDF$SEzprob <- Predictedzprob$se.fit
myDF$Predictedlink <- Predictedlink$fit
myDF$SElink <- Predictedlink$se.fit
myDF$Predictedzlink <- Predictedzlink$fit
myDF$SEzlink <- Predictedzlink$se.fit


## Extrapolated value in lighter colors
Extrapolated0 <- seq(min(data[,name]),max(data[,name]),length=Lvl) >= min(data[data[,"Plotcat"]==0,name]) & seq(min(data[,name]),max(data[,name]),length=Lvl) <= max(data[data[,"Plotcat"]==0,name])
Extrapolated1 <- seq(min(data[,name]),max(data[,name]),length=Lvl) >= min(data[data[,"Plotcat"]==1,name]) & seq(min(data[,name]),max(data[,name]),length=Lvl) <= max(data[data[,"Plotcat"]==1,name])
Extrapolated2 <- seq(min(data[,name]),max(data[,name]),length=Lvl) >= min(data[data[,"Plotcat"]==2,name]) & seq(min(data[,name]),max(data[,name]),length=Lvl) <= max(data[data[,"Plotcat"]==2,name])
Extrapolated <- as.factor(c(Extrapolated0,Extrapolated1,Extrapolated2))


## Colors and labels according the the reverse species or not (Trailing edge and leading edge reversed in some species)
if (CODE%in%c("PINPIN","PINNIG","QUEPET","QUEROB","PINPINA","POPTRE","ALNGLU","QUEPYR","QUESUB")){ 
  Mycol <- c("black", "blue", "red") # Pas normale
  Mylabels <- c("Core","Leading Edge","Trailing Edge") # pas normale
}else {
  Mylabels <- c("Core","Trailing Edge","Leading Edge") #Normale 
  Mycol <- c("black", "red", "blue")}


p1<-ggplot(data=myDF, # All
           aes(x=get(name), y=Predictedresp, col=Plotcat,shape=Plotcat,group=Plotcat,alpha=Extrapolated,
               ymin=Predictedresp-SEresp, ymax=Predictedresp+SEresp))+
  geom_line(size=1.5)+ #no size
  geom_point(size=2,stroke=2)+
  geom_linerange(size=1.5)+ # no size
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
# ggsave(filename = paste0(CODE,"RESPONSE_",name,"_",Modname,".png"),plot = p1, width = 6, height = 6, dpi=300) # To modify 
p1

p2<-ggplot(data=myDF, # All
           aes(x=get(name), y=Predictedmu, col=Plotcat,shape=Plotcat,group=Plotcat,alpha=Extrapolated,
               ymin=Predictedmu-SEmu, ymax=Predictedmu+SEmu))+
  geom_line(size=1.5)+ #no size
  geom_point(size=2,stroke=2)+
  geom_linerange(size=1.5)+ # no size
  scale_alpha_manual(values = c(0.3,1),guide=F)+
  scale_color_manual("",values = Mycol,labels=Mylabels)+
  scale_shape_manual("", values=c(1,2,3),labels=Mylabels)+
  labs(title=paste0(CODE),y=paste0("Predicted Mu"), x=paste0(name))+
  theme_light(base_size = 15)+
  theme(text = element_text(face="bold"),legend.position = c(0.5,0.9),legend.direction ="horizontal",
        axis.text.x = element_text(size=13,color="black"),axis.text.y = element_text(size=13,color="black"),
        legend.background=element_rect(fill="white",colour="black",size=0.2),
        panel.border = element_rect(colour = "black", fill=NA, size=0.8),
        axis.line = element_line(colour="black"),
        plot.title = element_text(size=18,hjust = 0.5),
        plot.caption = element_text(face="bold.italic"))
p2
# ggsave(filename = paste0(CODE,"MU_",name,"_",Modname,".png"),plot = p2, width = 6, height = 6, dpi=300) # To modify 

p3<-ggplot(data=myDF, # All
           aes(x=get(name), y=1-Predictedzprob, col=Plotcat,shape=Plotcat,group=Plotcat,alpha=Extrapolated,
               ymin=1-Predictedzprob-SEzprob, ymax=1-Predictedzprob+SEzprob))+
  geom_line(size=1.5)+ #no size
  geom_point(size=2,stroke=2)+
  geom_linerange(size=1.5)+ # no size
  scale_alpha_manual(values = c(0.3,1),guide=F)+
  scale_color_manual("",values = Mycol,labels=Mylabels)+
  scale_shape_manual("", values=c(1,2,3),labels=Mylabels)+
  labs(title=paste0(CODE),y=paste0("1-Predicted zprob"), x=paste0(name))+
  theme_light(base_size = 15)+
  theme(text = element_text(face="bold"),legend.position = c(0.5,0.9),legend.direction ="horizontal",
        axis.text.x = element_text(size=13,color="black"),axis.text.y = element_text(size=13,color="black"),
        legend.background=element_rect(fill="white",colour="black",size=0.2),
        panel.border = element_rect(colour = "black", fill=NA, size=0.8),
        axis.line = element_line(colour="black"),
        plot.title = element_text(size=18,hjust = 0.5),
        plot.caption = element_text(face="bold.italic"))
p3
# ggsave(filename = paste0(CODE,"1-ZPROB_",name,"_",Modname,".png"),plot = p3, width = 6, height = 6, dpi=300) # To modify 


p4<-ggplot(data=myDF, # All
           aes(x=get(name), y=1-Predictedzlink, col=Plotcat,shape=Plotcat,group=Plotcat,alpha=Extrapolated,
               ymin=1-Predictedzlink-SEzlink, ymax=1-Predictedzlink+SEzlink))+
  geom_line(size=1.5)+ #no size
  geom_point(size=2,stroke=2)+
  geom_linerange(size=1.5)+ # no size
  scale_alpha_manual(values = c(0.3,1),guide=F)+
  scale_color_manual("",values = Mycol,labels=Mylabels)+
  scale_shape_manual("", values=c(1,2,3),labels=Mylabels)+
  labs(title=paste0(CODE),y=paste0("1-Predicted zlink"), x=paste0(name))+
  theme_light(base_size = 15)+
  theme(text = element_text(face="bold"),legend.position = c(0.5,0.9),legend.direction ="horizontal",
        axis.text.x = element_text(size=13,color="black"),axis.text.y = element_text(size=13,color="black"),
        legend.background=element_rect(fill="white",colour="black",size=0.2),
        panel.border = element_rect(colour = "black", fill=NA, size=0.8),
        axis.line = element_line(colour="black"),
        plot.title = element_text(size=18,hjust = 0.5),
        plot.caption = element_text(face="bold.italic"))
p4
# ggsave(filename = paste0(CODE,"1-ZPROB_",name,"_",Modname,".png"),plot = p3, width = 6, height = 6, dpi=300) # To modify 


p5<-ggplot(data=myDF, # All
           aes(x=get(name), y=Predictedlink, col=Plotcat,shape=Plotcat,group=Plotcat,alpha=Extrapolated,
               ymin=Predictedlink-SElink, ymax=Predictedlink+SElink))+
  geom_line(size=1.5)+ #no size
  geom_point(size=2,stroke=2)+
  geom_linerange(size=1.5)+ # no size
  scale_alpha_manual(values = c(0.3,1),guide=F)+
  scale_color_manual("",values = Mycol,labels=Mylabels)+
  scale_shape_manual("", values=c(1,2,3),labels=Mylabels)+
  labs(title=paste0(CODE),y=paste0("1-Predicted zprob"), x=paste0(name))+
  theme_light(base_size = 15)+
  theme(text = element_text(face="bold"),legend.position = c(0.5,0.9),legend.direction ="horizontal",
        axis.text.x = element_text(size=13,color="black"),axis.text.y = element_text(size=13,color="black"),
        legend.background=element_rect(fill="white",colour="black",size=0.2),
        panel.border = element_rect(colour = "black", fill=NA, size=0.8),
        axis.line = element_line(colour="black"),
        plot.title = element_text(size=18,hjust = 0.5),
        plot.caption = element_text(face="bold.italic"))
# ggsave(filename = paste0(CODE,"1-ZPROB_",name,"_",Modname,".png"),plot = p3, width = 6, height = 6, dpi=300) # To modify 
p5



## Bootstrap
X.cond = model.matrix(lme4::nobars(formula(get(Modname))[-2]), myDF)
beta.cond = fixef(get(Modname))$cond
pred.cond = X.cond %*% beta.cond
ziformula = get(Modname)$modelInfo$allForm$ziformula #$
X.zi = model.matrix(lme4::nobars(ziformula), myDF)
beta.zi = fixef(get(Modname))$zi
pred.zi = X.zi %*% beta.zi

pred.ucount = exp(pred.cond)*(1-plogis(pred.zi))
cor.test(myDF$Predictedresp,pred.ucount)
# library(MASS)
set.seed(101)

pred.condpar.psim = mvrnorm(1000000,mu=beta.cond,Sigma=vcov(get(Modname))$cond)
pred.cond.psim = X.cond %*% t(pred.condpar.psim)
pred.zipar.psim = mvrnorm(1000000,mu=beta.zi,Sigma=vcov(get(Modname))$zi)
pred.zi.psim = X.zi %*% t(pred.zipar.psim)
pred.ucount.psim = exp(pred.cond.psim)*(1-plogis(pred.zi.psim))


ci.cond = t(apply(pred.cond.psim,1,quantile,c(0.025,0.975)))
ci.cond = data.frame(ci.cond)
names(ci.cond) = c("cond.low","cond.high")
pred.cond = data.frame(myDF, pred.cond, ci.cond)


ci.zi = t(apply(pred.zi.psim,1,quantile,c(0.025,0.975)))
ci.zi = data.frame(ci.zi)
names(ci.zi) = c("zi.low","zi.high")
pred.zi = data.frame(myDF, pred.zi, ci.zi)


ci.ucount = t(apply(pred.ucount.psim,1,quantile,c(0.025,0.975)))
ci.ucount = data.frame(ci.ucount)
names(ci.ucount) = c("ucount.low","ucount.high")
pred.ucount = data.frame(myDF, pred.ucount, ci.ucount)


p6<-ggplot(data=pred.ucount, # All
           aes(x=get(name), y=pred.ucount, col=Plotcat,shape=Plotcat,group=Plotcat,alpha=Extrapolated,
               ymin=ucount.low, ymax=ucount.high))+
  geom_line(size=1.5)+ #no size
  geom_point(size=2,stroke=2)+
  geom_linerange(size=1.5)+ # no size
  scale_alpha_manual(values = c(0.3,1),guide=F)+
  scale_color_manual("",values = Mycol,labels=Mylabels)+
  scale_shape_manual("", values=c(1,2,3),labels=Mylabels)+
  labs(title=paste0(CODE),y=paste0("Predicted Resp sim"), x=paste0(name))+
  theme_light(base_size = 15)+
  theme(text = element_text(face="bold"),legend.position = c(0.5,0.9),legend.direction ="horizontal",
        axis.text.x = element_text(size=13,color="black"),axis.text.y = element_text(size=13,color="black"),
        legend.background=element_rect(fill="white",colour="black",size=0.2),
        panel.border = element_rect(colour = "black", fill=NA, size=0.8),
        axis.line = element_line(colour="black"),
        plot.title = element_text(size=18,hjust = 0.5),
        plot.caption = element_text(face="bold.italic"))
p6
# ggsave(filename = paste0(CODE,"1-ZPROB_",name,"_",Modname,".png"),plot = p3, width = 6, height = 6, dpi=300) # To modify 


p7<-ggplot(data=pred.cond, # All
           aes(x=get(name), y=pred.cond, col=Plotcat,shape=Plotcat,group=Plotcat,alpha=Extrapolated,
               ymin=cond.low, ymax=cond.high))+
  geom_line(size=1.5)+ #no size
  geom_point(size=2,stroke=2)+
  geom_linerange(size=1.5)+ # no size
  scale_alpha_manual(values = c(0.3,1),guide=F)+
  scale_color_manual("",values = Mycol,labels=Mylabels)+
  scale_shape_manual("", values=c(1,2,3),labels=Mylabels)+
  labs(title=paste0(CODE),y=paste0("Predicted cond sim"), x=paste0(name))+
  theme_light(base_size = 15)+
  theme(text = element_text(face="bold"),legend.position = c(0.5,0.9),legend.direction ="horizontal",
        axis.text.x = element_text(size=13,color="black"),axis.text.y = element_text(size=13,color="black"),
        legend.background=element_rect(fill="white",colour="black",size=0.2),
        panel.border = element_rect(colour = "black", fill=NA, size=0.8),
        axis.line = element_line(colour="black"),
        plot.title = element_text(size=18,hjust = 0.5),
        plot.caption = element_text(face="bold.italic"))
# ggsave(filename = paste0(CODE,"1-ZPROB_",name,"_",Modname,".png"),plot = p3, width = 6, height = 6, dpi=300) # To modify 
p7




p8<-ggplot(data=pred.zi, # All
           aes(x=get(name), y=1-pred.zi, col=Plotcat,shape=Plotcat,group=Plotcat,alpha=Extrapolated,
               ymin=1-zi.low, ymax=1-zi.high))+
  geom_line(size=1.5)+ #no size
  geom_point(size=2,stroke=2)+
  geom_linerange(size=1.5)+ # no size
  scale_alpha_manual(values = c(0.3,1),guide=F)+
  scale_color_manual("",values = Mycol,labels=Mylabels)+
  scale_shape_manual("", values=c(1,2,3),labels=Mylabels)+
  labs(title=paste0(CODE),y=paste0("1-Predicted zlink sim"), x=paste0(name))+
  theme_light(base_size = 15)+
  theme(text = element_text(face="bold"),legend.position = c(0.5,0.9),legend.direction ="horizontal",
        axis.text.x = element_text(size=13,color="black"),axis.text.y = element_text(size=13,color="black"),
        legend.background=element_rect(fill="white",colour="black",size=0.2),
        panel.border = element_rect(colour = "black", fill=NA, size=0.8),
        axis.line = element_line(colour="black"),
        plot.title = element_text(size=18,hjust = 0.5),
        plot.caption = element_text(face="bold.italic"))
# ggsave(filename = paste0(CODE,"RESPONSE_",name,"_",Modname,".png"),plot = p1, width = 6, height = 6, dpi=300) # To modify 
p8




p1 # response # Mieux 
p6 # resp simul

p2 # mu
p5 # mulink
p7 # mulink simul # Pire

p3 # 1-zprob
p4 # 1-zlink
p8 # 1-zlink simul # Pire










