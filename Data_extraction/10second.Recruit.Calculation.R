# Alex on the 17/01/2020
# Script to calcultae the many additional variables at the plot scales 
# Trees number etc ... ... 

SavingInfo = "MyDFs(dir='bureau' or 'home',
CODE = 'BETPEN',
seuil = 0.8
seuilC = 0.6)
This function extract plot data in three tables (one with management, one scaled and one non scaled. Mortality binary data are added, as well as SPEI indexes that were calculated previously"
bannerBreak = "\n*********************************************************************\n"
cat(paste0(bannerBreak,SavingInfo,bannerBreak,"\n"))


library(parallel)
MyDFs <- function(CODE,seuil,seuilC) {
  Dir <- c("/home/achangenet/SynologyDrive/FUNDIV - NFI - Europe/")
  setwd(Dir)
  Dir = c(paste0(Dir, "our-data/species/", CODE, "/CLIMAP"))
  setwd(Dir)
  list.files(Dir, pattern = paste0(".rds"))
  df <- readRDS(paste0("Mydf2_", CODE, "_", seuil, "_", seuilC, ".rds")) #Base de données plot full
  dir.create(path=paste0(Dir,"/Recrut.Mortality.2020"))
  Dir <- paste0(Dir,"/Recrut.Mortality.2020")
  setwd(Dir)
  
  ###################################################################
  ###                                                             ###
  ###   Here is the df (scaled plots but some are managed)        ###
  ###                                                             ###
  ###################################################################
  Mydf <- df
  Mydf[Mydf$country=="FR" & Mydf$treestatus_th==0,"treestatus_th"] <- 4  ### Dead french trees noted 3 instead of 0. (Like other inventories) (for calculation of treeNBR)
  Mydf[Mydf$country=="FR" & Mydf$treestatus_th==1,"treestatus_th"] <- 2  ### Ingrowth french trees noted 2 instead of 1. (Like other inventories) (fo calculation of treeNBR)
  Mydf[Mydf$treestatus_th==6,"treestatus_th"] <- 5  ### For PICABI and PINSYL some trees (42 and 142 were noted as 6 and are effectively dead which false the calculations later on)
  
  ### First add treenumber calculation for the species only ! ! !
  dat.fundiv=split(Mydf[],as.character(Mydf$plotcode))
  data <- do.call(rbind, mclapply(dat.fundiv,function(df){
    i=df[1,"plotcode"]
    dat=matrix(NA,nrow(Mydf[Mydf$plotcode==i,]),7)
    dat[,7]=c(1:nrow(Mydf))[Mydf$plotcode==i]
    if(nrow(df)>0) {
      dat[,1] <- as.numeric(nrow(df)) # Number of lines = total number of trees of the species
      dat[,2] <- as.numeric(nrow(df[df$treestatus_th=="2" & !is.na(df$bachange_ha_yr.1),])) #ingrowth trees that have a values of growth (excluding some french trees without ir5 values for the calculation of mean BAI)
      dat[,3] <- as.numeric(nrow(df[df$treestatus_th=="2" & !is.na(df$dbh1),])) #ingrowth trees that have a values of dbh (to include trees used to calculate the dbh sum of the plot). +6000 compare to the previous one
      dat[,4] <- as.numeric(nrow(df[df$treestatus_th=="1",])) #recruited
      dat[,5] <- as.numeric(nrow(df[df$treestatus_th%in%c("3","4","5")  & !is.na(df$dbh1),])) #dead trees with values of dbh (excluding some french trees for the caclulation of meandbh)
      dat[,6] <- as.numeric(nrow(df[df$treestatus_th=="4" & !is.na(df$dbh1),])) # Real dead (used in calculation of mortality rate)
      # The previous line should include only real dead trees (4) but in the calculation of the total dbh of the plot (script 10bis.meanDBH),
      # all trees with a value of dbh1 were included. (Except just recruited trees)
      }
    as.data.frame(dat)},mc.cores=10,mc.silent=T))
  
  treeNbrJ=as.numeric(data[order(data[,7]),1]) # total lines of the species
  treeNbrJ.I.BAI=as.numeric(data[order(data[,7]),2]) #ingrowth with a value of bai
  treeNbrJ.I.dbh=as.numeric(data[order(data[,7]),3]) #ingrowth with a value of dbh1
  treeNbrJ.R=as.numeric(data[order(data[,7]),4]) # recruited
  treeNbrJ.M=as.numeric(data[order(data[,7]),6]) # dead (three statuts)
  treeNbrJ.Mall=as.numeric(data[order(data[,7]),5]) # dead (just the fourth statuts)
  
  Mydf[,"treeNbrJ"]=treeNbrJ # number of lines = number of individuals (including the mistakes, see above)
  Mydf[,"treeNbrJ.I.BAI"]=treeNbrJ.I.BAI # number of individuals ingrowth
  Mydf[,"treeNbrJ.I.dbh"]=treeNbrJ.I.dbh # number of individuals ingrowth
  Mydf[,"treeNbrJ.R"]=treeNbrJ.R # number of recruted
  Mydf[,"treeNbrJ.M"]=treeNbrJ.M # number of real dead (statut 4)
  Mydf[,"treeNbrJ.Mall"]=treeNbrJ.Mall # number of dead ( + statut 3 and statut 5 also). For check at the end  of the script
  
  Mydf[,"treeNbrJ.IMall"]=treeNbrJ.Mall+treeNbrJ.I.dbh   ## Ingrowth + Dead and 3 and 5  = Trees accounted for in mortality rate and for dbhmean
  Mydf[,"treeNbrJ.IR"]=treeNbrJ.R+treeNbrJ.I.BAI   ## Ingrowth + Recruted = Trees accounted for in rectuitment rate

  # Change name of the variable dbh.plot by dbh.plot.J
  colnames(Mydf)[333] <- "dbh.plot.J"
  ### Second : Fonction that calculate the sum dbh for each plot done already. 
  Mydf[,"BAI.plot.mean"]=Mydf$BAI.plot/Mydf$treeNbr # Mean Basal area increment / ha / plot. 
  #It is a sum at the plot/ha level divided by the total number of trees of the plot (all species)
  #(maybe should be the number of ingrowth trees of the plot but too long to calculate again)
  # This is sometimes negative and is therefore put to zero below
  
  Mydf[,"BAIj.plot.1.mean.J.I"]=Mydf$BAIj.plot.1/Mydf$treeNbrJ.I.BAI # Idem but only for ingrowth trees (this is the right measure) but check problem.
  # Check how we obtained the BAIj.plot.1: Calculated on all trees with a value of bachange.ha.yr.1
  # Check those that have a value of bachange.ha.yr.1
  summary(Mydf$BAIj.plot.1) ## Some values are NA 6371
  summary(Mydf$BAIj.plot.1.mean.J.I) # OK no infinite values 
  Mydf[,"dbhJ.IMall.plot.mean"]=Mydf$dbh.plot.J/Mydf$treeNbrJ.IMall # Here is the average dbh of the species at the plot (proxy of the age based on dead + false dead and ingrowth trees). 
  summary(Mydf$dbhJ.IMall.plot.mean)
  summary(Mydf$dbh.plot.J)
  
  #test <- Mydf[Mydf$treestatus_th=="2" & is.na(Mydf$bachange_ha_yr.1) & !is.na(Mydf$dbh1),c("BAIj.plot.1","bachange_ha_yr.1","dbhJ.IMall.plot.mean","dbh.plot.J","dbh1","treeNbrJ","treeNbrJ.I","treeNbrJ.M","treeNbrJ.R","treeNbrJ.Mall","treeNbrJ.IMall","treestatus_th","plotcode")]
  ## Check if there is probleme 
  test <- Mydf[!is.na(Mydf$bachange_ha_yr.1),] # Those that were used to calculate BAIj.plot ! Supposed to be only ingrowth trees
  test2 <- Mydf[!is.na(Mydf$dbh1),] # Those used to calculate sum dbh (only ingrowth but also dead ans false dead (4,3,5))
  test3 <- Mydf[Mydf$treestatus_th=="1" & Mydf$weight1!=0 & !is.na(Mydf$weight1),] # Those used to calculate mortality rate that were not suppose to (recruit with a weight at t1)
  test4 <- Mydf[Mydf$treestatus_th%in%c("3","4","5") & Mydf$weight2!=0 & !is.na(Mydf$weight2),] # Those used to calculate recruit rate that were not supposed to (dead or false dead with a weight at t2)
  test5 <- Mydf[!is.na(Mydf$BAIj.plot.1) & Mydf$treeNbrJ.I==0,] # Check possibility to be divided by zero = BAI value but no trees (noarmally no)
  test6 <- Mydf[!is.na(Mydf$dbh.plot.J) & Mydf$treeNbrJ.IMall==0,] # Idem with dbh 
  
  #test6[,c("BAIj.plot.1","bachange_ha_yr.1","dbh.plot.J","dbh1","treeNbrJ","treeNbrJ.I","treeNbrJ.M","treeNbrJ.R","treeNbrJ.Mall","treeNbrJ.IMall","treestatus_th","plotcode")] # 5630 infinite values              
  # Check with this lin is there is only ingrowth individuals (present at both time)
  capture.output(paste0("Individuals used to calculate BAIj.plot (Supposed to be only ingrowth trees): "),
                 file = paste0(Dir, "/df.CheckIndividualsForcalculation.txt")) 
  capture.output(summary(as.factor(test$treestatus_th)),file = paste0(Dir, "/df.CheckIndividualsForcalculation.txt"),append = T) 
  
  capture.output(paste0("Individuals used to calculate dbh.plot (supposed to be dead (3,4,5) and ingrowth individuals: "),
                 file = paste0(Dir, "/df.CheckIndividualsForcalculation.txt"),append = T)  
  capture.output(summary(as.factor(test2$treestatus_th)),file = paste0(Dir, "/df.CheckIndividualsForcalculation.txt"),append = T)
  
  capture.output(paste0("Recrut whose weight 1 is not 0 by country. These have hence been accounted for in the calculation of the density of dead: "),
                 file = paste0(Dir, "/df.CheckIndividualsForcalculation.txt"),append = T)  
  capture.output(summary(test3$country),file = paste0(Dir, "/df.CheckIndividualsForcalculation.txt"),append = T)
  
  capture.output(paste0("Dead whose weight 2 is not 0 by country. These have been accounted for in the calculation of the density of recruited: "),
                 file = paste0(Dir, "/df.CheckIndividualsForcalculation.txt"),append = T)  
  capture.output(summary(as.factor(test4$treestatus_th)),file = paste0(Dir, "/df.CheckIndividualsForcalculation.txt"),append = T)
  
  capture.output(paste0("Trees with BAIj at the plot scale but with no statut 2. Supposed to be 0 "),
                       file = paste0(Dir, "/df.CheckIndividualsForcalculation.txt"),append = T)
  capture.output(summary(as.factor(test5$treestatus_th)),file = paste0(Dir, "/df.CheckIndividualsForcalculation.txt"),append = T)
                              
  capture.output(paste0("Trees with dbhj at the plot scale but with no statut 2 3 4 or 5. Supposed to be 0 "),
                       file = paste0(Dir, "/df.CheckIndividualsForcalculation.txt"),append = T)
  capture.output(summary(as.factor(test6$treestatus_th)),file = paste0(Dir, "/df.CheckIndividualsForcalculation.txt"),append = T)
  
  capture.output(paste0("Differences of number of trees between nbrJ.Idbh and nbrTreeJ.I.bai :",  nrow(Mydf[Mydf$treestatus_th=="2" & !is.na(Mydf$bachange_ha_yr.1),])-nrow(Mydf[Mydf$treestatus_th=="2" & !is.na(Mydf$dbh1),]))
  ,file = paste0(Dir, "/df.CheckIndividualsForcalculation.txt"),append = T)
  
  capture.output(paste0("Sum at the plot/ha level divided by the total number of trees of the plot (all species). Not based only on the ingrowth trees but should be. This number is the number of individuals with a value < 0 and number of plots"),
                 file = paste0(Dir, "/df.CheckIndividualsForcalculation.txt"),append = T)
  capture.output(paste0(nrow(Mydf[Mydf$BAI.plot.mean<0 & !is.na(Mydf$BAI.plot.mean),])," & ",length(unique(Mydf[Mydf$BAI.plot.mean<0 & !is.na(Mydf$BAI.plot.mean),"plotcode"]))),
                 file = paste0(Dir, "/df.CheckIndividualsForcalculation.txt"),append = T)
  
  
  Mydf[Mydf$BAI.plot.mean<0 & !is.na(Mydf$BAI.plot.mean),"BAI.plot.mean"] <- 0 # Mean Basal area increment / ha / plot. Put to zero the negative values !
  
  saveRDS(get("Mydf"), paste0(Dir, "/Mydf3_", CODE, "_", seuil, "_", seuilC, "R.M.rds")) # To save on a new file Recrut.Mort
  ### Go to the plot scale
  
  dfplot <-Mydf[!duplicated(Mydf$plotcode), -c(5:7, 12:27, 334:336)] # Keep unique plot codes and remove colonnes 334 and 335 (BAIj.plot.mean and dbh.plot.mean)

  # Divide by the number of years (perthousand unit)
  dfplot[dfplot$country == "FR", "yearsbetweensurveys"] <- 5 # Add this line really matters
  
  
  ### Here add the calculation of the recruitment / years + sp.mort.bin + mortality rate 
  
  #  Sp Mortality rate as a binomial # Added on the 17/05/2018
  dfplot$sp.mort.bin <- dfplot[, "sp.mortality.plot.rate"]
  dfplot[!is.na(dfplot$sp.mortality.plot.rate) & dfplot$sp.mortality.plot.rate != 0, 326] = 1
  dfplot[!is.na(dfplot$sp.mortality.plot.rate) & dfplot$sp.mortality.plot.rate == 0, 326] = 0
  dfplot[, 326] = as.numeric(dfplot[, 326])
  
  ## Add recrut and mort as a annual rate 
  dfplot$sp.mortality.plot.rate.yr <-round((dfplot$sp.mortality.plot.rate * 1000 / dfplot$yearsbetweensurveys))
  dfplot$sp.recruitment.plot.rate.yr <-dfplot$sp.recruitment.plot.rate / dfplot$yearsbetweensurveys
  
  ### Add recrut as a binomial response
  dfplot$sp.recrut.bin <- dfplot[, "sp.recruitment.plot.rate"]
  dfplot[!is.na(dfplot$sp.recruitment.plot.rate) & dfplot$sp.recruitment.plot.rate != 0, "sp.recrut.bin"] = 1
  dfplot[!is.na(dfplot$sp.recruitment.plot.rate) & dfplot$sp.recruitment.plot.rate == 0, "sp.recrut.bin"] = 0
  # Summary as a latex doc
  capture.output(print(summary(as.factor(dfplot$sp.recruitment.plot.rate.yr))),file = paste0(Dir, "/dfplotR.M.sp.recruitment.plot.rate.yr.txt")) 
  capture.output(print(summary(as.factor(dfplot$sp.mortality.plot.rate.yr))),file = paste0(Dir, "/dfplotR.M.sp.mortality.plot.rate.yr.txt"))  
  
  ## Added on the 22th june Add the minimum and mean SPEI # Remove SPEI min
  Years <- c("spei01","spei03","spei06","spei12","spei18","spei24","spei36","spei48")
  for (i in Years){
    load(file = paste0("/home/achangenet/Documents/FUNDIV - NFI - Europe/our-data/climate/SPEI/plots_",i,"_1981_2015.RData"))
    dfplot[,paste0("mean_",i)] <- spei_plots_all$mean_spei_survey_years[match(dfplot$plotcode, spei_plots_all$plotcode,incomparables = NA)]
  }
  
  saveRDS(get("dfplot"), paste0(Dir, "/dfplot", CODE, seuil, "R.M.rds")) # Work at the plot scale NO SCALE
  dfplot <- readRDS(paste0("dfplot", CODE, seuil, "R.M.rds")) #Base de données plot full
  
  # Scale (based on all available data)
  ####### Added 05/07 # No scaled data but without NA and trasnfo ! ! ############
  #                                                                             ##
  dfplot <-dfplot[!is.na(dfplot$Plotcat) & dfplot$Plotcat != 10, ]              ## Unique plots that falls within the species distribution area remove NA and transition zone
  dfplot$Plotcat <- as.factor(dfplot$Plotcat)                                   ##
  print(nrow(dfplot[is.na(dfplot$sp.mortality.plot.rate.yr), ]))
  print(nrow(dfplot[is.na(dfplot$sp.recruitment.plot.rate.yr), ]))              ### The french plot will be excluded from the analysis 
  #dfplottest <-dfplot[!is.na(dfplot$sp.mortality.plot.count.yr), ]             ## Remove 2 NA or 14 (sp.count) To remove later

  Transfo <- c("BA.ha.plot.1","BA.O.plot.1","BAj.plot.1",                       ##
               "dbhJ.IMall.plot.mean","BAI.plot.mean",                          ## Problem negative values with BAU.plot.mean + pbm hist dbhJ.IMall...
               "BAIj.plot.1.mean.J.I")                                          ## Variable to transform 
  for (i in 1:length(Transfo)){                                                 ##
    dfplot[,paste0("log",Transfo[i])] <- log(dfplot[,Transfo[i]]+1)             ## Log transfo
    dfplot[,paste0("sqrt",Transfo[i])] <- sqrt(dfplot[,Transfo[i]])}            ## Sqrt transfo ## Probleme de NAN et inf OK
  saveRDS(get("dfplot"), paste0(Dir, "/dfplotbis", CODE, seuil, "R.M.rds"))     ## Save this database as Dfplot normal 
  #                                                                             ##
  ################################################################################
  

  dfplot[, c(10:17, 38:43, 47:50, 60:314,323:325,338:349)] <-scale(dfplot[, c(10:17, 38:43, 47:50, 60:314,323:325,338:349)],
                                                                   center = TRUE, scale =TRUE) #Center and reduced ALL except SPEI and number of trees
  dfplot2 <- dfplot
  
  capture.output(print(summary(as.factor(dfplot2[dfplot2$sp.mortality.plot.rate.yr ==0, "Plotcat"]))),file = paste0(Dir, "/dfplot2R.M.sp.mortality.rate.yr.txt")) # Output as a latex wrapped in a txt file
  capture.output(print(summary(as.factor(dfplot2[dfplot2$sp.mortality.plot.rate.yr > 0, "Plotcat"]))),file = paste0(Dir, "/dfplot2R.M.sp.mortality.rate.yr.txt"),append = T) # Output as a latex wrapped in a txt file
  
  # Idem with recruitment
  capture.output(print(summary(as.factor(dfplot2[dfplot2$sp.recruitment.plot.rate.yr ==0, "Plotcat"]))),file = paste0(Dir, "/dfplot2R.M.sp.recruitment.rate.yr.txt")) # Output as a latex wrapped in a txt file
  capture.output(print(summary(as.factor(dfplot2[dfplot2$sp.recruitment.plot.rate.yr > 0, "Plotcat"]))),file = paste0(Dir, "/dfplot2R.M.sp.recruitment.rate.yr.txt"),append = T) # Output as a latex wrapped in a txt file
  
  # At that point; there is not much plots left
  
  ###################################################################
  ###                                                             ###
  ###   Here is the dfplot2 (scaled plots but some are managed)   ###
  ###                                                             ###
  ###################################################################
  saveRDS(get("dfplot2"), paste0(Dir, "/dfplot2", CODE, seuil, "R.M.rds"))
  dfplot2 <- readRDS(paste0("dfplot2", CODE, seuil, "R.M.rds"))
  
  # Remove gestion(management 2 >= 1) & gest ???
  dfplot2 <-dfplot2[is.na(dfplot2$management2) | dfplot2$management2 == 0, ] # 255 NA (591,627 & 51)
  dfplot2 <- dfplot2[!is.na(dfplot2$management2) & dfplot2$management2 == 0, ] #  (581,539,51)
  #summary(as.factor(dfplot2[dfplot2$sp.mortality.plot.count.yr > 0, "Plotcat"])) # Check my categories
  # If any, remove it
  dfplot2$gest <- as.factor(dfplot2$gest)
  dfplot3 <- dfplot2[-which(dfplot2$gest == 2 | dfplot2$gest == 1), ] # Remove 459 + 537
  
  ##############################################################
  ###                                                        ###
  ###   Here is the dfplot3 (scaled plots and not managed)   ###
  ###                                                        ###
  ##############################################################
  capture.output(print(summary(as.factor(dfplot3[dfplot3$sp.mortality.plot.rate.yr ==
                                                   0, "Plotcat"]))),
                 file = paste0(Dir, "/dfplot3R.M.sp.mortality.rate.yr.txt")) # Output as a latex wrapped in a txt file
  capture.output(print(summary(as.factor(dfplot3[dfplot3$sp.mortality.plot.rate.yr > 0, "Plotcat"]))),
    file = paste0(Dir, "/dfplot3R.M.sp.mortality.rate.yr.txt"),
    append = T
  ) # Output as a latex wrapped in a txt file
  saveRDS(get("dfplot3"), paste0(Dir, "/dfplot3", CODE, seuil, "R.M.rds"))
  dfplot3 <- readRDS(paste0("dfplot3", CODE, seuil, "R.M.rds"))

  
###############################################################################################
##### Plot rate mortality and recruitment against number of trees | year and save the plot ####
###############################################################################################

  
  ## Test to investigate minimum number of trees at the plot scale 
  
  png(file=paste0(CODE,"_","Years&treenbrVSTotalRATE.png"),width = 12.37, height = 7.04,units = "in",res=600) # Save it 
  par(mfrow=c(2,2))
  plot(dfplot$sp.recruitment.plot.rate.yr~dfplot$yearsbetweensurveys)
  plot(dfplot$sp.recruitment.plot.rate.yr~dfplot$treeNbrJ.IR)
  plot(dfplot$sp.recruitment.plot.rate.yr~dfplot$treeNbrJ)
  plot(dfplot$sp.recruitment.plot.rate.yr~dfplot$treeNbrJ.R)
  dev.off()
  
  png(file=paste0(CODE,"_","Years&treenbrVSTotal.png"),width = 12.37, height = 7.04,units = "in",res=600) # Save it 
  par(mfrow=c(2,2))
  plot(dfplot$sp.recruitment.plot.rate~dfplot$yearsbetweensurveys)
  plot(dfplot$sp.recruitment.plot.rate~dfplot$treeNbrJ.IR)
  plot(dfplot$sp.recruitment.plot.rate~dfplot$treeNbrJ)
  plot(dfplot$sp.recruitment.plot.rate~dfplot$treeNbrJ.R)
  dev.off()
  
  png(file=paste0(CODE,"_","TreeIRnbrVSrate*4.png"),width = 12.37, height = 7.04,units = "in",res=600) # Save it 
  par(mfrow=c(2,2))
  test <- dfplot[dfplot$treeNbrJ.IR>3, ] # Investigate this number to have the right number. 
  plot(test$sp.recruitment.plot.rate.yr~test$treeNbrJ.IR,main=paste0(nrow(test)," plot / ",nrow(dfplot)))
  test <- dfplot[dfplot$treeNbrJ.IR>5, ] # Investigate this number to have the right number. 
  plot(test$sp.recruitment.plot.rate.yr~test$treeNbrJ.IR,main=paste0(nrow(test)," plot / ",nrow(dfplot)))
  test <- dfplot[dfplot$treeNbrJ.IR>7, ] # Investigate this number to have the right number. 
  plot(test$sp.recruitment.plot.rate.yr~test$treeNbrJ.IR,main=paste0(nrow(test)," plot / ",nrow(dfplot)))
  test <- dfplot[dfplot$treeNbrJ.IR>10, ] # Investigate this number to have the right number. 
  plot(test$sp.recruitment.plot.rate.yr~test$treeNbrJ.IR,main=paste0(nrow(test)," plot / ",nrow(dfplot)))
  dev.off()
  
  png(file=paste0(CODE,"_","TreeIRnbrVSrate*4.png"),width = 12.37, height = 7.04,units = "in",res=600) # Save it 
  par(mfrow=c(2,2))
  test <- dfplot[dfplot$treeNbrJ>3, ] # Investigate this number to have the right number. 
  plot(test$sp.recruitment.plot.rate.yr~test$treeNbrJ,main=paste0(nrow(test)," plot / ",nrow(dfplot)))
  test <- dfplot[dfplot$treeNbrJ>5, ] # Investigate this number to have the right number. 
  plot(test$sp.recruitment.plot.rate.yr~test$treeNbrJ,main=paste0(nrow(test)," plot / ",nrow(dfplot)))
  test <- dfplot[dfplot$treeNbrJ>7, ] # Investigate this number to have the right number. 
  plot(test$sp.recruitment.plot.rate.yr~test$treeNbrJ,main=paste0(nrow(test)," plot / ",nrow(dfplot)))
  test <- dfplot[dfplot$treeNbrJ>10, ] # Investigate this number to have the right number. 
  plot(test$sp.recruitment.plot.rate.yr~test$treeNbrJ,main=paste0(nrow(test)," plot / ",nrow(dfplot)))
  dev.off()
  
}

  




  