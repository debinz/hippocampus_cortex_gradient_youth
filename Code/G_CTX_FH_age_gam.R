library(mgcv)
library(lme4)
library(gamRR)
library(ggplot2)
library(itsadug)
library(tidyverse)
library(visreg)
library(R.matlab)
library(gratia)

#replace with absolute path of your work directory
setwd('E:\\Github\\hippocampus_cortex_gradient_youth')

dir <- "./Code"
result_dir <- "./Results/hipp_CTX_proj"
Fig_dir <- file.path(result_dir, "FH_dev_figure")
# Check if the figure directory exists
if (!dir.exists(Fig_dir)) {
  # If it doesn't exist, create the directory
  dir.create(Fig_dir)
  print(paste("Directory", Fig_dir, "created."))
}

Data <- paste(dir, "HCD_LS_2.0_subject_completeness_bh.csv", sep = "/") %>% read.csv()
head_motion <- paste(dir, "head_mv_par_all.csv", sep = "/") %>% read.csv()
Data <- rbind(Data,Data)
Data$mFD <- head_motion$mFD

L <- rep(0,652)
R <- rep(1,652)
Data$hemi <- c(L,R)

Measure_L <- paste(result_dir, "/indiv360/G_CTX_FH_all_L.csv", sep = "") %>% read.csv(header = FALSE)
Measure_R <- paste(result_dir, "/indiv360/G_CTX_FH_all_R.csv", sep = "") %>% read.csv(header = FALSE)
Measures <- rbind(Measure_L,Measure_L)
Measures <- data.frame(abs(Measures))

# we checked both the Pearson correlation ('*_p') and Spearman correlation ('*_s')
clist <- c('G1_p','G2_p','G3_p','G1_s','G2_s','G3_s') 
names(Measures) <- clist

Data <- cbind(Data, Measures)

Data$Age <- as.numeric(Data$Age)
Data$Sex <- as.factor(Data$Sex)
Data$site <- as.factor(Data$site)
Data$mFD<- as.numeric(Data$mFD)
Data$hemi <- as.factor(Data$hemi)

p_value<-matrix(NA,nrow=length(clist),ncol=1)
F_value<-matrix(NA,nrow=length(clist),ncol=1)
delta_R_sq<-matrix(NA,nrow=length(clist),ncol=1)
anova_pvalue<-matrix(NA,nrow=length(clist),ncol=1)

for (i in 1:length(clist)) {
  var_name<- clist[i]
  y1=Data[,var_name]
  # build GAM model
  fit_G <- gam(y1 ~ s(Age, k=4) + Sex + site +
                 s(hemi, bs="re") + mFD , data=Data, method="REML", na.action="na.omit")
  
  #GAM derivatives
  #Get derivatives of the smooth function using finite differences
  derv <- derivatives(fit_G, term = "s(Age)", interval = "simultaneous", unconditional = F) #derivative at 200 indices of smooth_var with a simultaneous CI
  #Identify derivative significance window(s)
  derv <- derv %>% #add "sig" column (TRUE/FALSE) to derv
    mutate(sig = !(0 > lower & 0 < upper)) #derivative is sig if the lower CI is not < 0 while the upper CI is > 0 (i.e., when the CI does not include 0)
  derv$sig_deriv = derv$derivative*derv$sig #add "sig_deriv derivatives column where non-significant derivatives are set to
  
  G_summary=summary(fit_G)
  p_value[i,] <- G_summary[["s.table"]]["s(Age)","p-value"]
  F_value[i,] <- G_summary[["s.table"]]["s(Age)","F"]
  
  #null model without s(Age)
  fit_G0 <- gam(y1 ~ Sex + site +
                  s(hemi, bs="re") + mFD , data=Data, method="REML", na.action="na.omit")
  G0_summary <- summary(fit_G0)
  
  delta_R_sq[i,] <- G_summary$r.sq-G0_summary$r.sq
  ### effect direction
  mean.derivative <- mean(derv$derivative)
  if(mean.derivative < 0){ #if the average derivative is less than 0, make the effect size estimate negative
    delta_R_sq[i,] <- delta_R_sq[i,]*-1}
  
  anova_pvalue[i,] <- anova.gam(fit_G0,fit_G,test='Chisq')$`Pr(>Chi)`[2]
  
  #visualization
  
  myplot <- visreg(fit_G, "Age", gg=T, type="conditional", scale = "response", 
                   overlay=T, partial=T, rug=F, ylab=paste(var_name, "_FH", sep = ""))
  myplot <- myplot + theme_bw() + 
    theme(axis.text = element_text(size=6, color = c("black")), 
          axis.title = element_text(size=6, color = c("black")), 
          axis.line = element_line(size=.22), axis.ticks = element_line(size=.22), 
          panel.grid=element_blank())  #remove background & grid
  myplot <- myplot + theme_classic() +
    scale_x_continuous(breaks=c(6, 8, 10, 12, 14, 16, 18, 20, 22), expand = c(0,.45))
  ggsave(filename=paste(var_name,"_FH_age.pdf", sep = ""), plot=myplot, 
         path=Fig_dir, 
         dpi = 600, width = 4 , height = 4)
}

GAM_Result<-cbind(F_value,p_value,delta_R_sq,anova_pvalue)
GAM_Result<-as.data.frame(GAM_Result)
colnames(GAM_Result) <- c("age_Fvalue", #GAM F-value for the age smooth term
                          "age_pvalue", #GAM p-value for the age smooth term
                          "age_deltaR2", #delta adjusted Rsq from age and age-null models
                          "Anova_age_pvalue") #Anova p-value comparing age and age-null models
                   
filepath<-paste(Fig_dir, "G_CTX_FH_age.csv", sep = "/")
write.csv(GAM_Result, filepath, row.names = F, quote = F)
