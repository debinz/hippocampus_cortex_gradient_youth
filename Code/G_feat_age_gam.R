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
result_dir <- "./Results/gradient_indiv/G_feat"
Fig_dir <- file.path(result_dir, "figure")

# Check if the figure directory exists
if (!dir.exists(Fig_dir)) {
  # If it doesn't exist, create the directory
  dir.create(Fig_dir)
  print(paste("Directory", Fig_dir, "created."))
}

# read the demographic information of HCP-D data
Data <- paste(dir, "HCD_LS_2.0_subject_completeness_bh.csv", sep = "/") %>% read.csv()
head_motion <- paste(dir, "head_mv_par_all.csv", sep = "/") %>% read.csv()

Data <- rbind(Data,Data)
Data$mFD <- head_motion$mFD

L <- rep(0,652)
R <- rep(1,652)
Data$hemi <- c(L,R)

Data$Age <- as.numeric(Data$Age)
Data$Sex <- as.factor(Data$Sex)
Data$site <- as.factor(Data$site)
Data$mFD<- as.numeric(Data$mFD)
Data$hemi <- as.factor(Data$hemi)

Age=Data$Age
Sex=Data$Sex
site=Data$site
mFD=Data$mFD
hemi=Data$hemi

for (g in 1:3) {
  
Measures <- paste(result_dir, "/feat_G", g, ".mat", sep = "") %>% readMat()
feat<-paste("feat.G",g,sep = "")
Feat_G_L<-Measures[feat][[1]][[1]][[1]]
Feat_G_R<-Measures[feat][[1]][[2]][[1]]

Feat_G <- rbind(Feat_G_L,Feat_G_R)
Feat_G<-data.frame(Feat_G)

clist <- c(paste0('Range',g), paste0('Std',g), paste0('Var',g), paste0('Expl',g), 
           paste0('Skew',g), paste0('Kurt',g))
names(Feat_G) <- clist
Data <- cbind(Data, Feat_G)

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

myplot <- visreg(fit_G, "Age", gg=T, type="conditional", scale = "response", col.point = "black",
                 overlay=T, partial=T, rug=F, ylab=paste("G", g, "_", strsplit(var_name,g), sep = ""))
myplot <- myplot + theme_bw() + theme(panel.grid=element_blank()) #remove background & grid
myplot <- myplot + theme_classic() +
  scale_x_continuous(breaks=c(6, 8, 10, 12, 14, 16, 18, 20, 22), expand = c(0,.45))
ggsave(filename=paste("G", g, "_", strsplit(var_name,g),"_age.pdf", sep = ""), plot=myplot,
       path=Fig_dir,
       width = 60, height = 60, units = "mm")
}

anova_pvalue_fdr=p.adjust(anova_pvalue,method = "fdr")

GAM_Result<-cbind(F_value,p_value,delta_R_sq,anova_pvalue,anova_pvalue_fdr)
GAM_Result<-as.data.frame(GAM_Result)
colnames(GAM_Result) <- c("age_Fvalue", #GAM F-value for the age smooth term
                          "age_pvalue", #GAM p-value for the age smooth term
                          "age_deltaR2", #delta adjusted Rsq from age and age-null models
                          "Anova_age_pvalue", #Anova p-value comparing age and age-null models
                          "Anova_age_adjpvalue") #fdr adjusted Anova p-value comparing age and age-null models
filepath<-paste(Fig_dir, "/G", g, "_feat_age.csv", sep = "")
write.csv(GAM_Result, filepath, row.names = F, quote = F)
}
