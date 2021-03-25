MyDFs <- function(CODE,seuil,seuilC) {
  Dir <- c("/Users/alexandrechangenet/SynologyDrive/FUNDIV - NFI - Europe/")
  setwd(Dir)
  Dir = c(paste0(Dir, "our-data/species/", CODE,"/Recrut.Mortality.2020"))
  setwd(Dir)
  list.files(Dir, pattern = paste0(".rds"))
  dfplot <- readRDS(paste0("dfplot", CODE, seuil, "R.M.rds")) #Base de données plot full
  dfplot$sp.recruitment.plot.rate.yr.2 <- round((dfplot$sp.recruitment.plot.rate * 1000 / dfplot$yearsbetweensurveys))

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
               "BAIj.plot.1.mean.J.I","treeNbrJ","sp.mortality.plot.rate.yr")   ## Variable to transform 
  for (i in 1:length(Transfo)){                                                 ##
    dfplot[,paste0("log",Transfo[i])] <- log(dfplot[,Transfo[i]]+1)             ## Log transfo
    dfplot[,paste0("sqrt",Transfo[i])] <- sqrt(dfplot[,Transfo[i]])}            ## Sqrt transfo ## Probleme de NAN et inf OK
  saveRDS(get("dfplot"), paste0(Dir, "/dfplotV2", CODE, seuil, "R.M.rds"))     ## Save this database as Dfplot normal 
  #                                                                             ##
  ################################################################################
  
  #Center and reduced ALL except SPEI. 
  # Edit 18/03/2020 => Just scaled not centered. And also for mortality rates/year and number of trees
  # Check lines 
  dfplot[, c(10:17, 38:43, 47:50, 60:315,323:325,327,339:354)] <-scale(dfplot[, c(10:17, 38:43, 47:50, 60:315,323:325,327,339:354)],
                                                                   center = F, scale =TRUE) 
  dfplot2 <- dfplot
  
  capture.output(print(summary(as.factor(dfplot2[dfplot2$sp.mortality.plot.rate.yr ==0, "Plotcat"]))),file = paste0(Dir, "/dfplot2V2R.M.sp.mortality.rate.yr.txt")) # Output as a latex wrapped in a txt file
  capture.output(print(summary(as.factor(dfplot2[dfplot2$sp.mortality.plot.rate.yr > 0, "Plotcat"]))),file = paste0(Dir, "/dfplot2V2R.M.sp.mortality.rate.yr.txt"),append = T) # Output as a latex wrapped in a txt file
  
  # Idem with recruitment
  capture.output(print(summary(as.factor(dfplot2[dfplot2$sp.recruitment.plot.rate.yr ==0, "Plotcat"]))),file = paste0(Dir, "/dfplot2V2R.M.sp.recruitment.rate.yr.txt")) # Output as a latex wrapped in a txt file
  capture.output(print(summary(as.factor(dfplot2[dfplot2$sp.recruitment.plot.rate.yr > 0, "Plotcat"]))),file = paste0(Dir, "/dfplot2V2R.M.sp.recruitment.rate.yr.txt"),append = T) # Output as a latex wrapped in a txt file
  
  # At that point; there is not much plots left
  
  ###################################################################
  ###                                                             ###
  ###   Here is the dfplot2 (scaled plots but some are managed)   ###
  ###                                                             ###
  ###################################################################
  
  ## Add line to add marginality as a continous variable (from the previous version)
  dfplotFinalV1 <- readRDS(paste0("dfplotFinal", CODE, seuil, "R.M.rds"))
  dfplot2$M.Plotcat <- dfplotFinalV1$M.Plotcat[match(rownames(dfplot2), rownames(dfplotFinalV1),incomparables = NA)]
  saveRDS(get("dfplot2"), paste0(Dir, "/dfplotFinalV2", CODE, seuil, "R.M.rds"))
  
  
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
                 file = paste0(Dir, "/dfplot3V2R.M.sp.mortality.rate.yr.txt")) # Output as a latex wrapped in a txt file
  capture.output(print(summary(as.factor(dfplot3[dfplot3$sp.mortality.plot.rate.yr > 0, "Plotcat"]))),
                 file = paste0(Dir, "/dfplot3V2R.M.sp.mortality.rate.yr.txt"),
                 append = T
  ) # Output as a latex wrapped in a txt file
  saveRDS(get("dfplot3"), paste0(Dir, "/dfplot3V2", CODE, seuil, "R.M.rds"))
  dfplot3 <- readRDS(paste0("dfplot3", CODE, seuil, "R.M.rds"))
  
  
  ###############################################################################################
  ##### Plot rate mortality and recruitment against number of trees | year and save the plot ####
  ###############################################################################################
  
  
  ## Test to investigate minimum number of trees at the plot scale 
  
  # png(file=paste0(CODE,"_","Years&treenbrVSTotalRATE.png"),width = 12.37, height = 7.04,units = "in",res=600) # Save it 
  # par(mfrow=c(2,2))
  # plot(dfplot$sp.recruitment.plot.rate.yr~dfplot$yearsbetweensurveys)
  # plot(dfplot$sp.recruitment.plot.rate.yr~dfplot$treeNbrJ.IR)
  # plot(dfplot$sp.recruitment.plot.rate.yr~dfplot$treeNbrJ)
  # plot(dfplot$sp.recruitment.plot.rate.yr~dfplot$treeNbrJ.R)
  # dev.off()
  # 
  # png(file=paste0(CODE,"_","Years&treenbrVSTotal.png"),width = 12.37, height = 7.04,units = "in",res=600) # Save it 
  # par(mfrow=c(2,2))
  # plot(dfplot$sp.recruitment.plot.rate~dfplot$yearsbetweensurveys)
  # plot(dfplot$sp.recruitment.plot.rate~dfplot$treeNbrJ.IR)
  # plot(dfplot$sp.recruitment.plot.rate~dfplot$treeNbrJ)
  # plot(dfplot$sp.recruitment.plot.rate~dfplot$treeNbrJ.R)
  # dev.off()
  # 
  # png(file=paste0(CODE,"_","TreeIRnbrVSrate*4.png"),width = 12.37, height = 7.04,units = "in",res=600) # Save it 
  # par(mfrow=c(2,2))
  # test <- dfplot[dfplot$treeNbrJ.IR>3, ] # Investigate this number to have the right number. 
  # plot(test$sp.recruitment.plot.rate.yr~test$treeNbrJ.IR,main=paste0(nrow(test)," plot / ",nrow(dfplot)))
  # test <- dfplot[dfplot$treeNbrJ.IR>5, ] # Investigate this number to have the right number. 
  # plot(test$sp.recruitment.plot.rate.yr~test$treeNbrJ.IR,main=paste0(nrow(test)," plot / ",nrow(dfplot)))
  # test <- dfplot[dfplot$treeNbrJ.IR>7, ] # Investigate this number to have the right number. 
  # plot(test$sp.recruitment.plot.rate.yr~test$treeNbrJ.IR,main=paste0(nrow(test)," plot / ",nrow(dfplot)))
  # test <- dfplot[dfplot$treeNbrJ.IR>10, ] # Investigate this number to have the right number. 
  # plot(test$sp.recruitment.plot.rate.yr~test$treeNbrJ.IR,main=paste0(nrow(test)," plot / ",nrow(dfplot)))
  # dev.off()
  # 
  # png(file=paste0(CODE,"_","TreeIRnbrVSrate*4.png"),width = 12.37, height = 7.04,units = "in",res=600) # Save it 
  # par(mfrow=c(2,2))
  # test <- dfplot[dfplot$treeNbrJ>3, ] # Investigate this number to have the right number. 
  # plot(test$sp.recruitment.plot.rate.yr~test$treeNbrJ,main=paste0(nrow(test)," plot / ",nrow(dfplot)))
  # test <- dfplot[dfplot$treeNbrJ>5, ] # Investigate this number to have the right number. 
  # plot(test$sp.recruitment.plot.rate.yr~test$treeNbrJ,main=paste0(nrow(test)," plot / ",nrow(dfplot)))
  # test <- dfplot[dfplot$treeNbrJ>7, ] # Investigate this number to have the right number. 
  # plot(test$sp.recruitment.plot.rate.yr~test$treeNbrJ,main=paste0(nrow(test)," plot / ",nrow(dfplot)))
  # test <- dfplot[dfplot$treeNbrJ>10, ] # Investigate this number to have the right number. 
  # plot(test$sp.recruitment.plot.rate.yr~test$treeNbrJ,main=paste0(nrow(test)," plot / ",nrow(dfplot)))
  # dev.off()
  
}






