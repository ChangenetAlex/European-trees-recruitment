#### Obtain the figure we want (first we need to obtain the myDFsim ALL data for all species script 22B.intercation.recruitment.R)
### Then load the script 24
rm(list = ls())
gc()

Dir <- c("~/SynologyDrive/FUNDIV - NFI - Europe/")
setwd(Dir)
source(paste0(Dir,"Myscripts/Fundiv.project/25bis.Stacked_interactions.R"))


## Check Plotcat BAj (mu p and responses)
## Check Plotcat Mortality (mu p and responses)
# BAJ*BAI
# BAJ*dbh 
# BAO*dbh 
# BAO*BAI 

nameVAR
PredVAR

Interaction(Modname="rW.yO.sumO.mRY1",name3 = nameVAR[8], nameSeq = nameVAR[2],Pred = PredVAR[3],save = T,free=F)
Interaction(Modname="rW.yO.sumO.mRY1",name3 = nameVAR[8], nameSeq = nameVAR[3],Pred = PredVAR[3],save = T,free=F)
Interaction(Modname="rW.yO.sumO.mRY1",name3 = nameVAR[8], nameSeq = nameVAR[4],Pred = PredVAR[3],save = T,free=F)
Interaction(Modname="rW.yO.sumO.mRY1",name3 = nameVAR[8], nameSeq = nameVAR[7],Pred = PredVAR[3],save = T,free=F)

Interaction(Modname="rW.yO.sumO.mRY1",name3 = nameVAR[3], nameSeq = nameVAR[2],Pred = PredVAR[2],save = T,free=F)
Interaction(Modname="rW.yO.sumO.mRY1",name3 = nameVAR[3], nameSeq = nameVAR[6],Pred = PredVAR[2],save = T,free=F)
Interaction(Modname="rW.yO.sumO.mRY1",name3 = nameVAR[4], nameSeq = nameVAR[2],Pred = PredVAR[2],save = T,free=F)
Interaction(Modname="rW.yO.sumO.mRY1",name3 = nameVAR[4], nameSeq = nameVAR[1],Pred = PredVAR[2],save = T,free=F)


library(patchwork)

MyVAR <- nameVAR[8]

A <- list.files(paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/paper2/Stacked_Interactions/Figure5"),pattern = paste0("FIXED_Predictedmu_VS_sqrtBAIj.plot.1.mean.J.I_for_3_level_of_logBAj.plot.1_rW.yO.sumO.mRY1.rds"),full.names = T)
assign("pLEGEND.",readRDS(A),envir = .GlobalEnv)
plegend <- pLEGEND.+theme(legend.direction ="horizontal",
                 legend.position = c(1, 1),
                 legend.justification = c("top"))+guides(colour = guide_legend(nrow = 2))

legend <- get_legend(plegend)
plot(legend)


A <- list.files(paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/paper2/Stacked_Interactions/"),pattern = paste0("Plotcat_rW.yO.sumO.mRY1.rds"),full.names = T)
# Load only the first half 

mapply(function(x,y){
  assign(paste0("p.",y),readRDS(x),envir = .GlobalEnv)
},x=A,y=1:8) ## be carreful with the order 





pall <- 
  p.7+theme(legend.position = "none")+
  labs(x="Relative drought index (SPEI)", y="Probability of recruitment (p)",title=paste0('a) Recruitment occurrence (p) VS SPEI along climatic marginality'))+
  
  p.3+theme(legend.position = "none")+
  labs(x="Relative drought index (SPEI)", y="Number of recruited trees (mu)",title=paste0('b) Recruitment count (mu) VS SPEI along climatic marginality'))+
  
  p.5+theme(legend.position = "none")+
  labs(x="INTRA (m².ha-1)", y="Probability of recruitment (p)",title=paste0('c) Recruitment occurrence (p) VS INTRA along climatic marginality'))+
  
  p.1+theme(legend.position = "none")+
  labs(x="INTRA (m².ha-1)", y="Number of recruited trees (mu)",title=paste0('d) Recruitment count (mu) VS INTRA along climatic marginality'))+
  
  p.8+theme(legend.position = "none")+
  labs(x="INTER (m².ha-1)", y="Probability of recruitment (p)",title=paste0('e) Recruitment occurrence (p) VS INTER along climatic marginality'))+
  
  p.4+theme(legend.position = "none")+
  labs(x="INTER (m².ha-1)", y="Number of recruited trees (mu)",title=paste0('f) Recruitment count (mu) VS INTER along climatic marginality'))+
  
  p.6+theme(legend.position = "none")+
  labs(x="M (trees.ha-1.yr-1)", y="Probability of recruitment (p)",title=paste0('g) Recruitment occurrence (p) VS M along climatic marginality'))+
  
  p.2+theme(legend.position = "none")+
  labs(x="M (trees.ha-1.yr-1)", y="Number of recruited trees (mu)",title=paste0('h) Recruitment count (mu) VS M along climatic marginality'))+
  plot_layout(ncol=2)


pall2 <- pall+legend+theme(plot.margin = unit(c(0,0,0,0), "cm"))+plot_layout(nrow=5,heights = c(1,1,1,1,0.2))

pall
save_plot(filename = paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/paper2/","Stacked.1.pdf"),
          plot = pall2,base_width = 12, base_height = 6, dpi = 300 ,units = "in",nrow=5,ncol=2)

## Figure 5


A <- list.files(paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/paper2/Stacked_Interactions/Figure5/"),pattern = paste0("rW.yO.sumO.mRY1.rds"),full.names = T)
# Load only the first half 
mapply(function(x,y){
  assign(paste0("p.",y),readRDS(x),envir = .GlobalEnv)
},x=A,y=1:8) ## be carreful with the order 


pall <- 
  p.6+theme(legend.position = "none")+
  labs(x="Relative drought index (SPEI)", y="Probability of recruitment (p)",title=paste0('a) Recruitment occurrence (p) VS SPEI along INTRA levels (m².ha-1)'))+
  
  p.2+theme(legend.position = "none")+
  labs(x="Relative drought index (SPEI)", y="Number of recruited trees (mu)",title=paste0('b) Recruitment count (mu) VS SPEI along INTRA levels (m².ha-1)'))+
  
  p.7+theme(legend.position = "none")+
  labs(x="Relative drought index (SPEI)", y="Probability of recruitment (p)",title=paste0('c) Recruitment occurrence (p) VS SPEI along INTER levels (m².ha-1)'))+
  
  p.3+theme(legend.position = "none")+
  labs(x="Relative drought index (SPEI)", y="Number of recruited trees (mu)",title=paste0('d) Recruitment count (mu) VS SPEI along INTER levels (m².ha-1)'))+
  
  p.8+theme(legend.position = "none")+
  labs(x="G (cm².ha-1.yr-1)", y="Probability of recruitment (p)",title=paste0('e) Recruitment occurrence (p) VS plot growth rate along INTRA levels (m².ha-1)'))+
  
  p.4+theme(legend.position = "none")+
  labs(x="G (cm².ha-1.yr-1)", y="Number of recruited trees (mu)",title=paste0('f) Recruitment count (mu) VS plot growth rate along INTRA levels (m².ha-1)'))+
  
  p.5+theme(legend.position = "none")+
  labs(x="DBH (cm)", y="Probability of recruitment (p)",title=paste0('g) Recruitment occurrence (p) VS plot DBH along INTER levels (m².ha-1)'))+
  
  p.1+theme(legend.position = "none")+
  labs(x="DBH (cm)", y="Number of recruited trees (mu)",title=paste0('h) Recruitment count (mu) VS plot DBH along INTER levels (m².ha-1)'))+
  plot_layout(ncol=2)


pall2 <- pall+legend+theme(plot.margin = unit(c(0,0,0,0), "cm"))+plot_layout(nrow=5,heights = c(1,1,1,1,0.2))

pall
save_plot(filename = paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/paper2/","Stacked.2.pdf"),
          plot = pall2,base_width = 12, base_height = 6, dpi = 300 ,units = "in",nrow=5,ncol=2)






## Figure 5 in two part
## First we take interactions with SPEI of Intra and Inter 
A <- list.files(paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/paper2/Stacked_Interactions/Figure5/"),pattern = paste0("rW.yO.sumO.mRY1.rds"),full.names = T)
# Load only the first half 
A <- A[c(2,3,6,7)]

mapply(function(x,y){
  assign(paste0("p.",y),readRDS(x),envir = .GlobalEnv)
},x=A,y=1:4) ## be carreful with the order 


pall <- 
  p.3+theme(legend.position = "none")+
  labs(x="Relative drought index (SPEI)", y="Probability of recruitment (p)",title=paste0('a) Recruitment occurrence (p) VS SPEI along INTRA levels (m².ha-1)'))+
  geom_point(size=0.1,stroke=0.1)+
  
  p.1+theme(legend.position = "none")+
  labs(x="Relative drought index (SPEI)", y="Number of recruited trees (mu)",title=paste0('b) Recruitment count (mu) VS SPEI along INTRA levels (m².ha-1)'))+
  geom_point(size=0.1,stroke=0.1)+
  
  p.4+theme(legend.position = "none")+
  labs(x="Relative drought index (SPEI)", y="Probability of recruitment (p)",title=paste0('c) Recruitment occurrence (p) VS SPEI along INTER levels (m².ha-1)'))+
  geom_point(size=0.1,stroke=0.1)+
  
  p.2+theme(legend.position = "none")+
  labs(x="Relative drought index (SPEI)", y="Number of recruited trees (mu)",title=paste0('d) Recruitment count (mu) VS SPEI along INTER levels (m².ha-1)'))+
  geom_point(size=0.1,stroke=0.1)+
  
  plot_layout(ncol=2)


pall2 <- pall+legend+theme(plot.margin = unit(c(0,0,0,0), "cm"))+plot_layout(nrow=3,heights = c(1,1,0.2))
pall
save_plot(filename = paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/paper2/","Figure5.Part1.2.pdf"),
          plot = pall2,base_width = 12, base_height = 5, dpi = 300 ,units = "in",nrow=3,ncol=2)



F## Part 2 Demography X competition
A <- list.files(paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/paper2/Stacked_Interactions/Figure5/"),pattern = paste0("rW.yO.sumO.mRY1.rds"),full.names = T)
# Load only the first half 
A <- A[c(1,4,5,8)]

mapply(function(x,y){
  assign(paste0("p.",y),readRDS(x),envir = .GlobalEnv)
},x=A,y=1:4) ## be carreful with the order 

pall <- 
  p.4+theme(legend.position = "none")+
  labs(x="G (cm².ha-1.yr-1)", y="Probability of recruitment (p)",title=paste0('a) Recruitment occurrence (p) VS plot growth rate along INTRA levels (m².ha-1)'))+
  geom_point(size=0.1,stroke=0.1)+
  p.2+theme(legend.position = "none")+
  labs(x="G (cm².ha-1.yr-1)", y="Number of recruited trees (mu)",title=paste0('b) Recruitment count (mu) VS plot growth rate along INTRA levels (m².ha-1)'))+
  geom_point(size=0.1,stroke=0.1)+
  
  p.3+theme(legend.position = "none")+
  labs(x="DBH (cm)", y="Probability of recruitment (p)",title=paste0('c) Recruitment occurrence (p) VS plot DBH along INTER levels (m².ha-1)'))+
  geom_point(size=0.1,stroke=0.1)+
  
  p.1+theme(legend.position = "none")+
  labs(x="DBH (cm)", y="Number of recruited trees (mu)",title=paste0('d) Recruitment count (mu) VS plot DBH along INTER levels (m².ha-1)'))+
  geom_point(size=0.1,stroke=0.1)+
  plot_layout(ncol=2)

pall2 <- pall+legend+theme(plot.margin = unit(c(0,0,0,0), "cm"))+plot_layout(nrow=3,heights = c(1,1,0.2))
save_plot(filename = paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/paper2/","Figure5.Part2.pdf"),
          plot = pall2,base_width = 12, base_height = 5, dpi = 300 ,units = "in",nrow=3,ncol=2)




