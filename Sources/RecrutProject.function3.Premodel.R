# Alex on the 17/03/2020
# Script to visualize my data before modeling and saving (or not) all the plots in one file 
# Function extraction + function Bootstrap

SavingInfo = "Premodel(
z = dfplot2, which is my dataframe with all variables I want to analyse
Resp=Resp, Which is my response variable, by defaut this is sp.mortality.plot.rate.yr
Explain=Explain, Which are the variable I want to study, by defaut a few are included 
size=4, For the cleveland plots, this is the par(mfrow=c(size,size)) parameter
save = T, to save in the new directory premodel or just display my plots)

Require the tidyr, lattice and usdm packages.

!!! WARNING !!! 
Before running this function you need to be in the directory corresponding to the SPECIES you want to study. 
!!! WARNING !!!

Example : Premodel(dfplot2,Resp=sp.mortality.ba,Explain=c(x,y...),size=3,save=T)

"
bannerBreak = "\n*********************************************************************\n"
cat(paste0(bannerBreak,SavingInfo,bannerBreak,"\n"))

## Packages ##

#library(usdm)
library(lattice)
library(tidyr)
library(usdm)
library(dplyr)

sizeN = 1800
resol = 15

panel.hist <- function(x, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col = "cyan", ...)
}

panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}

# Dotchart
Premodel <- function(z,Resp,Explain,size,save=F){
  
  Name <- if(all(list.dirs(Dir,full.names=F)!="Premodel")) dir.create(paste0(Dir,"/Premodel/"))
  #A <- cbind(Resp,Explain[-c(length(Explain))]) # Erreur 
  A <- c(Resp,Explain[-c(length(Explain))])
  c1 <- rep(z[,c(Resp)],times=length(Explain))
  z$Plotcat <- as.numeric(as.character(z$Plotcat))
  c2 <- gather(z[,c(Explain)])
  test <- cbind(c1,c2[,c(2,1)])

  if (save==T) png(file=paste0(Dir,"Premodel/",CODE,Resp,"_",Explain[1],"_Scatterplots.",sizeN,".png"),width=sizeN,height=sizeN*0.9,res=sizeN/resol)
    par(mfrow=c(1,1))
    print(xyplot(test[,1]~test[,2] | test[,3],col=1,
           type = c("p", "smooth"),
           scales=list(alternating=TRUE,
                       x=list(relation="free"),
                       y=list(relation="same")),
           xlab="Explanatory variables",
           ylab="Mortality",
           panel=function(x,y){
             panel.grid(h=-1,v=2)
             panel.points(x,y,col=1,cex=0.5,pch=19)
             #panel.loess(x,y,lwd=2) draw the lines
           }))
    
    dev.off()
    if (save==T) png(file=paste0(Dir,"Premodel/",CODE,Resp,"_",Explain[1],"_ClevelandPlot.",sizeN,".png"),width=sizeN,height=sizeN*0.9,res=sizeN/resol) # Save it 
    par(mfrow=c(size,2))
    for (i in 1:length(A)){
      dotchart(z[,c(A[i])],
               groups=factor(z[,c(last(Explain))]),
               ylab = last(Explain),
               xlab = A[i],
               main = "Cleveland dotplot",
               pch = c(1:3))
    }
    graphics.off()
    if (save==T) png(file=paste0(Dir,"Premodel/",CODE,Resp,"_",Explain[1],"_Hist.",sizeN,".png"),width=sizeN,height=sizeN*0.9,res=sizeN/resol) # Save it 
    par(mfrow=c(2,2))
    for (i in 2:length(A)){
      hist(z[,c(A[i])],breaks=40,xlab=c(paste0("Centered & normed Scale ",A[i])),ylab="Frequency",main="")
    }
    dev.off()
    if (save==T) png(file=paste0(Dir,"Premodel/",CODE,Resp,"_",Explain[1],"_PairsPlot.",sizeN,".png"),width=sizeN,height=sizeN*0.9,res=sizeN/resol) # Save it 
    par(mfrow=c(1,1))
    pairs(z[,c(A)],panel = panel.smooth,upper.panel=panel.cor,diag.panel = panel.hist)
    dev.off()
}
