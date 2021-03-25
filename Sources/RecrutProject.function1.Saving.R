# Alex on the 17/02/2020
# Script to save the output of my models with package GLMMTMB

SavingInfo = "Saving 
A R function to save your model and the summary in a friendly csv format  

The package GLMMTMB is needed
Input data: Your model (i.e. Hnbinom1)

create and save three elements: 
(1) A directory named after your model to be found in the recrut.mortality.2020 directory 
(2) A csv file with the estimated coefs, standards deviations and significance levels 
(3) Your model "
bannerBreak = "\n*********************************************************************\n"
cat(paste0(bannerBreak,SavingInfo,bannerBreak,"\n"))

library(glmmTMB)

Saving <- function(x,y){
  #Dir <- paste0("/home/achangenet/SynologyDrive/FUNDIV - NFI - Europe/our-data/species/",CODE,"/CLIMAP/Recrut.Mortality.2020")
  Dir <- paste0("~/SynologyDrive/FUNDIV - NFI - Europe/our-data/species/",y,"/CLIMAP/Recrut.Mortality.2020")
  if (length(grep(substitute(x),pattern = "get",fixed=T,value=F,invert=F))!=0){ # Added to simplify the model selection process 
    z <- paste0(Mymod,num)
    dir.create(paste0(Dir,"/",z,"/"))
    save(x, file = paste0(Dir,"/",z,"/",z,".rda")) # save the model as an RDA file 
  }else {dir.create(paste0(Dir,"/",deparse(substitute(x)),"/"))
    save(x, file = paste0(Dir,"/",deparse(substitute(x)),"/",deparse(substitute(x)),".rda"))
    z <- deparse(substitute(x))} # save the model as an RDA file 
  A <- summary(x)
  B <- as.data.frame(rbind(A[["AICtab"]]))
  B <- rbind(colnames(B),B)
  C <- c(A[["sigma"]],rep(NA,4))
  D <- rbind(B,C)
  
  if(length(try(as.data.frame(ranef(x))))!=1){
    E <- as.data.frame(ranef(x))
    rownames(E) <- paste0(E[,1],"_",E[,2],"_",E[,3],"_",E[,4])
    E <- round(E[,c(5,6)],3)
    E <- rbind(colnames(E),E)
    E$. <- NA
    E$.. <- NA
    E$... <- NA
    F <- rbind(rep("condmodel",4),round(A[["coefficients"]][["cond"]],3),rep("zimodel",4),round(A[["coefficients"]][["zi"]],3))
    F <- as.data.frame(cbind(rownames(F),F))
    colnames(E) <- colnames(F)
    colnames(D) <- colnames(F)
    G <- rbind(F,D,E)
  }else
  {F <- rbind(rep("condmodel",4),round(A[["coefficients"]][["cond"]],3),rep("zimodel",4),round(A[["coefficients"]][["zi"]],3))
  F <- as.data.frame(cbind(rownames(F),F))
  colnames(D) <- colnames(F)
  G <- rbind(F,D)}
  write.table(G[,-c(1)],paste0(Dir,"/",z,"/",z,"_Table.csv"),sep=";",col.names = NA, row.names = G$V1) # Variable outputs in a table.csv
}
