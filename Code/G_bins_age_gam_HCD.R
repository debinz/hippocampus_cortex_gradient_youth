library(mgcv)
library(lme4)
library(gamRR)
library(ggplot2)
library(itsadug)
library(tidyverse)
library(visreg)
library(openxlsx)
library(viridis)
library(gratia)

#replace with absolute path of your work directory
setwd('E:\\Github\\hippocampus_cortex_gradient_youth')

dir <- "./Code"
result_dir <- "./Results/gradient_indiv/G_bins"

Data <- paste(dir, "HCD_LS_2.0_subject_completeness_bh.csv", sep = "/") %>% read.csv()
head_motion <- paste(dir, "head_mv_par_all.csv", sep = "/") %>% read.csv()

Data <- rbind(Data,Data)
Data$mFD <- head_motion$mFD

L <- rep(0,652)
R <- rep(1,652)
Data$hemi <- c(L,R)
Data$Sex <- as.factor(Data$Sex)
Data$mFD<- as.numeric(Data$mFD)
Data$Age <- as.numeric(Data$Age)
Data$hemi <- as.factor(Data$hemi)
Data$site <- as.factor(Data$site)

Age=Data$Age
Sex=Data$Sex
mFD=Data$mFD
site=Data$site
hemi=Data$hemi

bin_nums=c(3,6,9,12,15)

label="mye"
Fig_dir <- file.path(result_dir, label, "figure")

# Check if the figure directory exists
if (!dir.exists(Fig_dir)) {
  # If it doesn't exist, create the directory
  dir.create(Fig_dir)
  print(paste("Directory", Fig_dir, "created."))
}


for (g_ind in 1:3) {

for (bin_num in bin_nums) {

MeasureL <- paste(result_dir, "/",label,"/G", g_ind, "_grp_", bin_num, "bins_",label,"_L.csv", sep = "") %>% read.csv(header = FALSE)
MeasureR <- paste(result_dir, "/",label,"/G", g_ind, "_grp_", bin_num, "bins_",label,"_R.csv", sep = "") %>% read.csv(header = FALSE)
Measures=rbind(MeasureL,MeasureR)

y_pred <- c()

p_value<-matrix(NA,nrow=bin_num,ncol=1)
F_value<-matrix(NA,nrow=bin_num,ncol=1)
delta_R_sq<-matrix(NA,nrow=bin_num,ncol=1)
anova_pvalue<-matrix(NA,nrow=bin_num,ncol=1)

for (bin in 1:bin_num) {
  Data$Measure <- Measures[,bin]
  Data$Measure <- as.numeric(Data$Measure)
  Measure=Data$Measure
  
  mod_gam <- gam(Measure ~ s(Age, bs="cs", k=4) + site + Sex #+ s(subj, bs="re")
                 +s(hemi, bs="re") + mFD , data=Data, method="REML") #, na.action="na.omit"
  
  #GAM derivatives
  #Get derivatives of the smooth function using finite differences
  derv <- derivatives(mod_gam, term = "s(Age)", interval = "simultaneous", unconditional = F) #derivative at 200 indices of smooth_var with a simultaneous CI
  #Identify derivative significance window(s)
  derv <- derv %>% #add "sig" column (TRUE/FALSE) to derv
    mutate(sig = !(0 > lower & 0 < upper)) #derivative is sig if the lower CI is not < 0 while the upper CI is > 0 (i.e., when the CI does not include 0)
  derv$sig_deriv = derv$derivative*derv$sig #add "sig_deriv derivatives column where non-significant derivatives are set to
  
  G_summary=summary(mod_gam)
  p_value[bin,] <- G_summary[["s.table"]]["s(Age)","p-value"]
  F_value[bin,] <- G_summary[["s.table"]]["s(Age)","F"]
  
  #null model without s(Age)
  mod_gam0 <- gam(Measure ~ Sex + site +
                  s(hemi, bs="re") + mFD , data=Data, method="REML", na.action="na.omit")
  G0_summary <- summary(mod_gam0)
  
  delta_R_sq[bin,] <- G_summary$r.sq-G0_summary$r.sq
  ### effect direction
  mean.derivative <- mean(derv$derivative)
  if(mean.derivative < 0){ #if the average derivative is less than 0, make the effect size estimate negative
    delta_R_sq[bin,] <- delta_R_sq[bin,]*-1}
  
  anova_pvalue[bin,] <- anova.gam(mod_gam0,mod_gam,test='Chisq')$`Pr(>Chi)`[2]
  # zero centered prediction  
  age_term = predict(mod_gam, type = 'terms')[, 's(Age)']
  #res_partial = residuals(mod_gam) + age_term
  Age = mod_gam$model$Age
  y_pred = cbind(y_pred,age_term[order(Age)])
  
  # orig prediction
  # Age = mod_gam$model$Age
  # Age_pred <- data.frame(Age=Age[order(Age)])
  # Age_pred <- fitted_values(mod_gam,
  #                         data=Age_pred,
  #                        #se.fit = TRUE,
  #                         #type = "link",
  #                          exclude = c("Sex", "s(hemi)","mFD","site"))
  # y_pred = cbind(y_pred,Age_pred$`fitted`)
}
p_value_fdr=p.adjust(p_value,method = "fdr")
anova_pvalue_fdr=p.adjust(anova_pvalue,method = "fdr")

GAM_Result<-cbind(F_value,p_value,delta_R_sq,anova_pvalue,p_value_fdr,anova_pvalue_fdr)
GAM_Result<-as.data.frame(GAM_Result)
colnames(GAM_Result) <- c("age_Fvalue", #GAM F-value for the age smooth term
                          "age_pvalue", #GAM p-value for the age smooth term
                          "age_deltaR2", #delta adjusted Rsq from age and age-null models
                          "Anova_age_pvalue", #Anova p-value comparing age and age-null models
                          "age_adjpvalue", #fdr adjusted p-value for the age smooth term
                          "Anova_age_adjpvalue") #fdr adjusted Anova p-value comparing age and age-null models
filepath<-paste(Fig_dir, "/G", g_ind, "_", bin_num, "bins_age.csv", sep = "")
write.csv(GAM_Result, filepath, row.names = F, quote = F)

cols <- plasma(bin_num)
names(cols) <- c(paste(label, "_bin", matrix(1:bin_num), sep=""))

# tiff(paste(dir, "/G_bins/",label,"/figure/G", g_ind, "_", bin_num, "bins_deltaRsq.tiff", sep = ""))
# barplot(height=GAM_Result$age_deltaR2,xlab="bin_num",ylab="Delta_adj_Rsq", col=cols)
# dev.off()

pdf(paste(Fig_dir, "/G", g_ind, "_", bin_num, "bins_deltaRsq.pdf", sep = ""))
barplot(height=GAM_Result$age_deltaR2,xlab="bin_num",ylab="Delta_adj_Rsq", col=cols)
dev.off()

myplot <- ggplot(data = mod_gam$model,aes(x=Age)) + 
  geom_line(aes(y=y_pred[,1], x = Age[order(Age)]), size=1, alpha = .8, colour=cols[1]) + 
  labs(y = paste("G", g_ind, sep=""), x = expression(Age))
theme_classic()

for (bin in 2:bin_num) {
  myplot <- myplot+
    geom_line(aes_(y=y_pred[,bin], x = Age[order(Age)]), size=1, alpha = .8, colour=cols[bin]) + 
    labs(y = paste("G", g_ind, sep=""), x = expression(Age))
  theme_classic()
}
myplot <- myplot+
  scale_colour_manual("",values=cols, breaks=names(cols))+
  theme_bw() + theme(axis.text = element_text(size=6, color = c("black")), 
                     axis.title = element_text(size=6, color = c("black")), 
                     axis.line = element_line(size=.22), axis.ticks = element_line(size=.22), 
                     panel.grid=element_blank())  #remove background & grid
myplot <- myplot + theme_classic() +
  scale_x_continuous(breaks=c(6, 8, 10, 12, 14, 16, 18, 20, 22), expand = c(0,.45))
myplot

filepath2<-paste(Fig_dir, "/G", g_ind, "_", bin_num, "bins_age.pdf", sep = "")
ggsave(filename = filepath2, device = "pdf", myplot, dpi = 500, width = 2.65 , height = 2.6)
filepath2<-paste(Fig_dir, "/G", g_ind, "_", bin_num, "bins_age.tiff", sep = "")
ggsave(filename = filepath2, device = "tiff", myplot, dpi = 500, width = 2.65 , height = 2.6)
}
}