library(reshape2)
library(plyr)
library(ggplot2)
library(cowplot)
library(forecat)
rm(list = ls())
gc()
Modname <- "rW.yO.sumO.mRY1"

tableresults <- readRDS(paste0("~/SynologyDrive/FUNDIV - NFI - Europe/our-data/species/Recrut.All/Model_All_",Modname,".rds"))


# Change names 
tableresults$variable <- gsub("logdbhJ.IMall.plot.mean", "DBH", tableresults$variable, ignore.case = FALSE, perl = FALSE,fixed = FALSE, useBytes = FALSE)
tableresults$variable <- gsub("logtreeNbrJ", "D", tableresults$variable, ignore.case = FALSE, perl = FALSE,fixed = FALSE, useBytes = FALSE)
tableresults$variable <- gsub("mean_spei12", "SPEI", tableresults$variable, ignore.case = FALSE, perl = FALSE,fixed = FALSE, useBytes = FALSE)
tableresults$variable <- gsub("Plotcat1", "TE", tableresults$variable, ignore.case = FALSE, perl = FALSE,fixed = FALSE, useBytes = FALSE)
tableresults$variable <- gsub("Plotcat2", "LE", tableresults$variable, ignore.case = FALSE, perl = FALSE,fixed = FALSE, useBytes = FALSE)
tableresults$variable <- gsub("sqrtBA.O.plot.1", "Inter", tableresults$variable, ignore.case = FALSE, perl = FALSE,fixed = FALSE, useBytes = FALSE)
tableresults$variable <- gsub("logBAj.plot.1", "Intra", tableresults$variable, ignore.case = FALSE, perl = FALSE,fixed = FALSE, useBytes = FALSE)
tableresults$variable <- gsub("sqrtBAIj.plot.1.mean.J.I", "G", tableresults$variable, ignore.case = FALSE, perl = FALSE,fixed = FALSE, useBytes = FALSE)
tableresults$variable <- gsub("logsp.mortality.plot.rate.yr", "M", tableresults$variable, ignore.case = FALSE, perl = FALSE,fixed = FALSE, useBytes = FALSE)
tableresults$variable <- gsub("(Intercept)", "Intercept", tableresults$variable, ignore.case = FALSE, perl = FALSE,fixed = FALSE, useBytes = FALSE)
tableresults$variable <- gsub(":", " * ", tableresults$variable, ignore.case = FALSE, perl = FALSE,fixed = FALSE, useBytes = FALSE)
unique(tableresults$variable)


tableresults.zi <- tableresults[tableresults$model=="zi",] # Keep only zi model
tableresults.cond <- tableresults[tableresults$model=="cond",] # Keep only zi model

# Remove intercept and Z-values  
#tableresults.zi <- tableresults.zi[tableresults.zi$variable!="(Intercept)",c(1:5,7)] # Keep only zi model
#tableresults.cond <- tableresults.cond[tableresults.cond$variable!="(Intercept)",c(1:5,7)] # Keep only cond model

# Remove z-values 
tableresults.zi <- tableresults.zi[,c(1:5,7)] # Keep only zi model
tableresults.cond <- tableresults.cond[,c(1:5,7)] # Keep only cond model

# Round values 
tableresults.zi[,c(4:6)] <- round(tableresults.zi[,c(4:6)],2)
tableresults.cond[,c(4:6)] <- round(tableresults.cond[,c(4:6)],2)

# Replace significiance by *
tableresults.zi[tableresults.zi$signif<=0.05,"signif"] <- "*"
tableresults.zi[tableresults.zi$signif>0.05,"signif"] <- ""
tableresults.cond[tableresults.cond$signif<=0.05,"signif"] <- "*"
tableresults.cond[tableresults.cond$signif>0.05,"signif"] <- ""

## Reorder the factor in the order we want 
tableresults.cond2 <- tableresults.cond %>%
  mutate(name = fct_relevel(variable,
                            "(Intercept)", "DBH", "D", "TE", "LE", "SPEI", "Inter", "Intra", 
                              "G", "M", "TE * SPEI", "LE * SPEI", "TE * Inter", "LE * Inter", 
                              "TE * Intra", "LE * Intra", "TE * G", "LE * G", "TE * M", "LE * M", 
                              "SPEI * Intra", "SPEI * Inter","DBH * Intra","DBH * Inter","Intra * G", 
                              "Inter * G"))

tableresults.zi2 <- tableresults.zi %>%
  mutate(name = fct_relevel(variable,
                            "(Intercept)", "DBH", "D", "TE", "LE", "SPEI", "Inter", "Intra", 
                            "G", "M", "TE * SPEI", "LE * SPEI", "TE * Inter", "LE * Inter", 
                            "TE * Intra", "LE * Intra", "TE * G", "LE * G", "TE * M", "LE * M", 
                            "SPEI * Intra", "SPEI * Inter","DBH * Intra","DBH * Inter","Intra * G", 
                            "Inter * G"))

### Save cond part 

p <- ggplot(tableresults.cond2[tableresults.cond2$variable!="(Intercept)",],aes(estimates, name,color=signif)) + 
  #geom_point(size=0.5) +
  geom_pointrange(aes(xmin=estimates-SE, xmax=estimates+SE),fatten = 0.5,)+
  facet_wrap(~ species,scales = "free_x")+
  scale_y_discrete(limits=rev)+
  geom_vline(xintercept = 0, colour = "black", linetype = 2)+
  scale_color_manual(values=c("blue","red"),name="Significance",labels=c("Not significant","Significant"))+
  theme_bw() + xlab("Coefficient Estimate") + ylab("") +
  ggtitle("b) Coefficient estimates by species in conditionnal model")+
  theme(panel.background = element_rect(fill="white", colour="black", size=1,
                                        linetype=1, color="black"))+
  theme(text = element_text(face="bold"),
        legend.direction ="vertical",
        legend.position = c(0.90, 0.17),
        legend.justification = c(1,1),
        legend.box.just = "right",
        legend.box.background = element_rect(color="black", size=0.5),
        legend.box.margin = margin(6, 12, 6, 12),
        axis.text.x = element_text(size=11,color="black"),axis.text.y = element_text(size=11,color="black",hjust = 1),
        axis.title.x = element_text(size=14),axis.title.y = element_text(size=12),
        legend.key = element_rect(fill = "white", colour = "black",size = 1),
        legend.key.heigh = unit(3,"line"),
        legend.key.width = unit(1.5,"line"),
        legend.text=element_text(size=14),
        legend.title=element_text(size=16),
        #legend.justification = "center",
        legend.margin = margin(0,0,0,0),
        #legend.background=element_rect(fill="white",colour="black",size=0.5,linetype="solid"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.3),
        axis.line = element_line(colour="black"),
        plot.title = element_text(size=17,hjust = 0,vjust = 0),
        plot.margin = margin(0.2,0.2,0.2,0.2, "cm"),
        plot.caption = element_text(face="bold.italic"))
save_plot(filename = paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/paper2/Coef.Dwplot.COND.png"),plot = p,base_width = 10, base_height = 15, dpi = 500 ,units = "in",nrow=1,ncol=1)


### Save ZI part

p2 <- ggplot(tableresults.zi2[tableresults.zi2$variable!="(Intercept)",],aes(estimates, name,color=signif)) + 
  #geom_point(size=0.5) +
  geom_pointrange(aes(xmin=estimates-SE, xmax=estimates+SE),fatten = 0.5,)+
  facet_wrap(~ species,scales = "free_x")+
  scale_y_discrete(limits=rev)+
  geom_vline(xintercept = 0, colour = "black", linetype = 2)+
  scale_color_manual(values=c("blue","red"),name="Significance",labels=c("Not significant","Significant"))+
  theme_bw() + xlab("Coefficient Estimate") + ylab("") +
  ggtitle("a) Coefficient estimates by species in zero-inflated model")+
  theme(panel.background = element_rect(fill="white", colour="black", size=1,
                                        linetype=1, color="black"))+
  theme(text = element_text(face="bold"),
        legend.direction ="vertical",
        legend.position = c(0.90, 0.17),
        legend.justification = c(1,1),
        legend.box.just = "right",
        legend.box.background = element_rect(color="black", size=0.5),
        legend.box.margin = margin(6, 12, 6, 12),
        axis.text.x = element_text(size=11,color="black"),axis.text.y = element_text(size=11,color="black",hjust = 1),
        axis.title.x = element_text(size=14),axis.title.y = element_text(size=12),
        legend.key = element_rect(fill = "white", colour = "black",size = 1),
        legend.key.heigh = unit(3,"line"),
        legend.key.width = unit(1.5,"line"),
        legend.text=element_text(size=14),
        legend.title=element_text(size=16),
        #legend.justification = "center",
        legend.margin = margin(0,0,0,0),
        #legend.background=element_rect(fill="white",colour="black",size=0.5,linetype="solid"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.3),
        axis.line = element_line(colour="black"),
        plot.title = element_text(size=17,hjust = 0,vjust = 0),
        plot.margin = margin(0.2,0.2,0.2,0.2, "cm"),
        plot.caption = element_text(face="bold.italic"))
save_plot(filename = paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/paper2/Coef.Dwplot.ZI.png"),plot = p2,base_width = 10, base_height = 15, dpi = 500 ,units = "in",nrow=1,ncol=1)

## Assemble both plot together
library(patchwork)
pall <- p2+theme(plot.margin = unit(c(0,0,0,0), "cm"))+
  p+theme(plot.margin = unit(c(0,0,0,0), "cm"))+
  plot_layout(ncol=2)


## Save this plot and add a line
png(paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/paper2/Coef.Dwplot.both.png"),width = 20, height = 25, res = 500 ,units = "in")
print(pall)
grid.draw(linesGrob(x = unit(c(0.505, 0.505), "npc"), y = unit(c(0, 01), "npc"),gp = gpar(col="black", lwd=4,lty=2)))
dev.off()

