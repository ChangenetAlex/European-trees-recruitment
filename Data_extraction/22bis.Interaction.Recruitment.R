rm(list = ls())
gc()
library(parallel)
library(forecast)
library(glmmTMB)
library(bbmle)
library(ggplot2)
library(reshape2)

Dir <- c("/home/achangenet/SynologyDrive/FUNDIV - NFI - Europe/")
setwd(Dir)
source(paste0(Dir,"Myscripts/Fundiv.project/RecrutProject.function1.Saving.R"))
Allcode <- c("PINPINA","ABIALB","ACEPSE","ALNGLU","BETPEN","CASSAT","FAGSYL",
             "FRAEXC","PICABI","PINHAL","PINNIG","PINPIN","PINSYL","POPNIG",
             "POPTRE","QUEILE","QUEPET","QUEPYR","QUEROB","QUESUB")
Allseuil <- c(0.8,0.7,0.55,0.6,0.8,0.7,0.7,0.7,0.8,0.7,0.7,0.7,0.8,0.7,0.7,0.8,0.7,0.7,0.8,0.7)
AllseuilC <- c(0.6,0.6,0.5,0.5,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6)


####################################################
Lvl <- 30                                         ## Number of values for the interaction graphs
nBoot <- 1000                                     ## 
nCoeur <- 20                                      ## 
Yportion <- 0.75                                  ## 
Modname <- "SCALEDRrate.yr.QualiF.TF.noY.nbinom1" ## Name of the model 
nameVAR <- c(                                     ##
  "BAIj.plot.1.mean.J.I",                         ##
  "mean_spei12",                                  ##
  "BA.O.plot.1",                                  ##
  "BAj.plot.1",                                   ##
  "treeNbrJ",                                     ##
  "dbhJ.IMall.plot.mean",                         ##
  "sp.mortality.plot.rate.yr")                    ##
####################################################  


for (i in 2:length(Allcode)){
# for (i in 1:2){
  CODE <- Allcode[i]
  seuilC <- AllseuilC[i]
  seuil <- Allseuil[i]
  Dir = paste0("/home/achangenet/SynologyDrive/FUNDIV - NFI - Europe/our-data/species/",CODE,"/CLIMAP/Recrut.Mortality.2020/")
  setwd(Dir)
  ## Charge the dfplot that is located in the species recrut file 
  dfplot <- try(readRDS(paste0("dfplotFinal",CODE,seuil,"R.M.rds")),silent = T) #Base de données plot
  
  # Also load dfplo2 or three
  
  # Change directory the the model we want to investigate 
  Dir <- c(paste0(Dir,Modname,"/")) # Directory 
  setwd(Dir)
  
  # Add recruitment as a count per year and as a total count
  dfplot$sp.recruitment.plot.rate.yr.2 <- round((dfplot$sp.recruitment.plot.rate * 1000 / dfplot$yearsbetweensurveys))
  dfplot$sp.recruitment.plot.rate.2 <- round(dfplot$sp.recruitment.plot.rate * 1000)
  dfplot$sp.mortality.plot.rate.2 <- round(dfplot$sp.mortality.plot.rate * 1000)
  
  # Charge the model 
  try(assign(Modname,get(load(file = paste0(Modname,".rda")))),silent = F)
  rm(x)
  ### For each single one except one give average values
  data <- get(Modname)$frame
  ## Fixer tous les paramètres à leur moyennes 
  ## Sauf on veut 1/3 Plocat chaque cat
  myDF = as.data.frame(matrix(ncol = length(colnames(data)), nrow = Lvl*3, NA)) # Df vierge
  colnames(myDF) <- colnames(data) # Bon noms de colonnes
  Dir <- c("/home/achangenet/SynologyDrive/FUNDIV - NFI - Europe/our-data/species/Recrut.All/")
  dir.create(paste0(Dir,Modname,"/")) # Creation of a new file for this model
  Dir <- paste0(Dir,Modname,"/")
  
  ## Variable que l'on veut regarder (fixer et faire varier)
  for (i in 1:length(nameVAR)){
    name <- nameVAR[i]
    
    if(try(is.null(dim(data[,paste0("sqrt",name)])),silent=T)==T){name <- paste0("sqrt",name)}
    if(try(is.null(dim(data[,paste0("log",name)])),silent=T)==T){name <- paste0("log",name)}
    
     if (length(unlist(ranef(get(Modname))))!=0){
       
      myDF[((Lvl*3)*(i-1)+1):(Lvl*3*i),2] <- as.numeric(rep(mean(data[,2]),Lvl*3))
      myDF[((Lvl*3)*(i-1)+1):(Lvl*3*i),3] <- NA # country
      myDF[((Lvl*3)*(i-1)+1):(Lvl*3*i),4] <- as.numeric(rep(mean(data[,4]),Lvl*3))
      myDF[((Lvl*3)*(i-1)+1):(Lvl*3*i),6] <- as.numeric(rep(mean(data[,6]),Lvl*3))
      myDF[((Lvl*3)*(i-1)+1):(Lvl*3*i),7] <- as.numeric(rep(mean(data[,7]),Lvl*3))
      myDF[((Lvl*3)*(i-1)+1):(Lvl*3*i),8] <- as.numeric(rep(mean(data[,8]),Lvl*3))
      myDF[((Lvl*3)*(i-1)+1):(Lvl*3*i),9] <- as.numeric(rep(mean(data[,9]),Lvl*3))
      myDF[((Lvl*3)*(i-1)+1):(Lvl*3*i),10] <- as.numeric(rep(mean(data[,10]),Lvl*3))
      myDF[((Lvl*3)*(i-1)+1):(Lvl*3*i),name] <- seq(min(data[,name]),max(data[,name]),length=Lvl)

      } else {
        
        myDF[((Lvl*3)*(i-1)+1):(Lvl*3*i),2] <- as.numeric(rep(mean(data[,2]),Lvl*3))
        myDF[((Lvl*3)*(i-1)+1):(Lvl*3*i),3] <- as.numeric(rep(mean(data[,3]),Lvl*3))
        myDF[((Lvl*3)*(i-1)+1):(Lvl*3*i),5] <- as.numeric(rep(mean(data[,5]),Lvl*3))
        myDF[((Lvl*3)*(i-1)+1):(Lvl*3*i),6] <- as.numeric(rep(mean(data[,6]),Lvl*3))
        myDF[((Lvl*3)*(i-1)+1):(Lvl*3*i),7] <- as.numeric(rep(mean(data[,7]),Lvl*3))
        myDF[((Lvl*3)*(i-1)+1):(Lvl*3*i),8] <- as.numeric(rep(mean(data[,8]),Lvl*3))
        myDF[((Lvl*3)*(i-1)+1):(Lvl*3*i),9] <- as.numeric(rep(mean(data[,9]),Lvl*3))
        myDF[((Lvl*3)*(i-1)+1):(Lvl*3*i),name] <- seq(min(data[,name]),max(data[,name]),length=Lvl)}

  }
  
  myDF[,"Plotcat"] <- as.factor(c(rep("0",Lvl),rep("1",Lvl),rep("2",Lvl))) # Here our data are ready to be used for bootstrapping 

  # Here we fit the same model on 80% of our data (we do it X time) and do prediction

  pred.boot =rep(NA,nBoot)
  pred.boot.resp =rep(NA,nBoot)
  pred.boot.mu =rep(NA,nBoot)
  pred.boot.zi =rep(NA,nBoot)
  Bootpred <- mclapply(pred.boot,function(train){
      calibrate.data = sample(1:nrow(data), Yportion*nrow(data)) # 66% of our data
      train = data[calibrate.data,]
      sub_binomial = update(get(Modname), data=train) # Fit my model on my sample 
      
    ## Predicted response (count * proba)
    pred.boot.resp <- predict(sub_binomial,newdata = myDF,type = c("response"))
    
    ## Predicted mu (count part)
    pred.boot.mu <- predict(sub_binomial,newdata = myDF,type = c("conditional"))
    
    ## Predicted probability to observe a zero event (= proba of a non recruit event) 
    pred.boot.zi <- predict(sub_binomial,newdata = myDF,type = c("zprob"))
    
    pred.boot <- list(pred.boot.resp,pred.boot.mu,pred.boot.zi)
    },mc.cores=nCoeur,mc.silent=T)
  
  # Keep the predictions we want 
  Bootpred <- as.data.frame(Bootpred) # obtain a dataframe with nBoot * 3 columns MAYBE BUG BECAUSE OF ERROR
  
  
  pred.boot.resp <- Bootpred[,seq((ncol(Bootpred)/nBoot)-2,ncol(Bootpred)-2,3),] # First columns being predicted response
  pred.boot.mu <- Bootpred[,seq((ncol(Bootpred)/nBoot)-1,ncol(Bootpred)-1,3),] # Second columns being pred mu
  pred.boot.zi <- Bootpred[,seq((ncol(Bootpred)/nBoot),ncol(Bootpred),3),] # Thirs being predicted zi

  message(paste0("Species is ",CODE," with model ",Modname, " and variable ",name))
  Myerrors <- unlist(Bootpred[sapply(Bootpred, function(x) inherits(x, "try-error"))]) # Take the errors out 
    if (length(Myerrors)!=0) {print(Myerrors) # If any, print it on the screen 
      Bootpred <- Bootpred[sapply(Bootpred, function(x) !inherits(x, "try-error"))] # Keep the non-errors
      message(paste0(nBoot-length(Bootpred)," Bootstraps on ",nBoot," failed because of an error")) 
      pred.boot.resp <- pred.boot.resp[sapply(pred.boot.resp, function(x) !inherits(x, "try-error"))] # Keep the non-errors
      message(paste0("among which ",nBoot-length(pred.boot.resp)," Bootstraps on ",nBoot," failed because of an error in response estimation")) 
      pred.boot.mu <- pred.boot.mu[sapply(pred.boot.mu, function(x) !inherits(x, "try-error"))] # Keep the non-errors
      message(paste0("among which ",nBoot-length(pred.boot.mu)," Bootstraps on ",nBoot," failed because of an error in mu estimation")) 
      pred.boot.zi <- pred.boot.zi[sapply(pred.boot.zi, function(x) !inherits(x, "try-error"))] # Keep the non-errors
      message(paste0("among which ",nBoot-length(pred.boot.zi)," Bootstraps on ",nBoot," failed because of an error in zi estimation")) 
    }
  
  message(paste0("Species is ",CODE," with model ",Modname, " and variable ",name))
    message(paste0("There is ",length(Bootpred[Bootpred>2000])," predicted values above 2000, the maximum being ",max(Bootpred),". These are not going to be replaced by NA"))
    # Above <- as.data.frame(table(unlist(sapply(as.list(Bootpred),function(x) which(x >2000)))))
    # message(paste0(Above$Var1," * ",Above$Freq,sep="\n"),"Are the lines for which a huge values as been estimated and the number of times \n \n")
    #Bootpred[Bootpred>2000] <- NA
    
    
    # Calculate the mean and SD for the three predictions 
    myDF[,"pred.boot.resp.means"] = rowMeans(pred.boot.resp,na.rm=T) # mean of each row for the range
    myDF[,"pred.boot.resp.SD"] = apply(pred.boot.resp, 1, function(x) sd(x,na.rm=T)) # sd of each row for the range 
    yMAX.resp <- max(myDF[,"pred.boot.resp.means"])+1.96*max(myDF[,"pred.boot.resp.SD"]) # Here is the maximum values that will be the limit for the plot
    
    myDF[,"pred.boot.mu.means"] = rowMeans(pred.boot.mu,na.rm=T) # mean of each row for the range
    myDF[,"pred.boot.mu.SD"] = apply(pred.boot.mu, 1, function(x) sd(x,na.rm=T)) # sd of each row for the range 
    yMAX.mu <- max(myDF[,"pred.boot.mu.means"])+1.96*max(myDF[,"pred.boot.mu.SD"]) # Here is the maximum values that will be the limit for the plot
    
    myDF[,"pred.boot.zi.means"] = rowMeans(pred.boot.zi,na.rm=T) # mean of each row for the range
    myDF[,"pred.boot.zi.SD"] = apply(pred.boot.zi, 1, function(x) sd(x,na.rm=T)) # sd of each row for the range 
    yMAX.zi <- max(myDF[,"pred.boot.zi.means"])+1.96*max(myDF[,"pred.boot.zi.SD"]) # Here is the maximum values that will be the limit for the plot
    # yMAX <- max(c(Means_Bootpred0,Means_Bootpred1))+max(c(SD_Bootpred0,SD_Bootpred1)) # Here is the maximum values when removing one category (need to modify the number)
    
    #Save these mean and SD predictions as well as the data that were used to generate it
    saveRDS(get("myDF"), paste0("/home/achangenet/SynologyDrive/FUNDIV - NFI - Europe/our-data/species/Recrut.All/SCALEDRrate.yr.QualiF.TF.noY.nbinom1/myDF.SIM_",Lvl,"_", CODE,".rds"))


    ## GGPLOT 
    
    for (i in 1:length(nameVAR)){
      name <- nameVAR[i]
      dir.create(paste0(Dir,nameVAR[i],"/"))
      setwd(paste0(Dir,nameVAR[i],"/"))
      
      if(try(is.null(dim(data[,paste0("sqrt",name)])),silent=T)==T){name <- paste0("sqrt",name)}
      if(try(is.null(dim(data[,paste0("log",name)])),silent=T)==T){name <- paste0("log",name)}
    
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
    
    p1<-ggplot(data=myDF[((Lvl*3)*(i-1)+1):(Lvl*3*i),], # All
               aes(x=get(name), y=pred.boot.resp.means, col=Plotcat,shape=Plotcat,group=Plotcat,alpha=Extrapolated,
                   ymin=pred.boot.resp.means-pred.boot.resp.SD, ymax=pred.boot.resp.means+pred.boot.resp.SD))+
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
    assign(paste0("presp_",CODE,"_",name),p1,envir = .GlobalEnv) # added to save the figure object because it is really long
    saveRDS(p1,file = paste0("presp_",CODE,"_",name,".rds"))
    ggsave(filename = paste0(CODE,"RESPONSE_",name,"_",Modname,".png"),plot = p1, width = 6, height = 6, dpi=300) # To modify 
    
    
    p2<-ggplot(data=myDF[((Lvl*3)*(i-1)+1):(Lvl*3*i),], # All
               aes(x=get(name), y=pred.boot.mu.means, col=Plotcat,shape=Plotcat,group=Plotcat,alpha=Extrapolated,
                   ymin=pred.boot.mu.means-pred.boot.mu.SD, ymax=pred.boot.mu.means+pred.boot.mu.SD))+
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
    assign(paste0("pmu_",CODE,"_",name),p2,envir = .GlobalEnv) # added to save the figure object because it is really long
    saveRDS(p2,file = paste0("pmu_",CODE,"_",name,".rds"))
    ggsave(filename = paste0(CODE,"MU_",name,"_",Modname,".png"),plot = p2, width = 6, height = 6, dpi=300) # To modify 
    
    p3<-ggplot(data=myDF[((Lvl*3)*(i-1)+1):(Lvl*3*i),], # All
               aes(x=get(name), y=1-pred.boot.zi.means, col=Plotcat,shape=Plotcat,group=Plotcat,alpha=Extrapolated,
                   ymin=1-pred.boot.zi.means-pred.boot.zi.SD, ymax=1-pred.boot.zi.means+pred.boot.zi.SD))+
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
    assign(paste0("pzprob_",CODE,"_",name),p3,envir = .GlobalEnv) # added to save the figure object because it is really long
    saveRDS(p3,file = paste0("pzprob_",CODE,"_",name,".rds"))
    ggsave(filename = paste0(CODE,"1-ZPROB_",name,"_",Modname,".png"),plot = p3, width = 6, height = 6, dpi=300) # To modify 
    }
}

  
  