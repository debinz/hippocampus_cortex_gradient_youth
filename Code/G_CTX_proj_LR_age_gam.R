library(mgcv)
library(lme4)
library(gamRR)
library(ggplot2)
library(itsadug)
library(dplyr)
library(visreg)
library(openxlsx)
library(viridis)
library(gratia)
library(Matrix)
extrafont::loadfonts()

#replace with absolute path of your work directory
setwd('E:\\Github\\hippocampus_cortex_gradient_youth')

dir <- "./Code"
Surf_temp <- "./Surf_temp"

result_dir <- "./Results/hipp_CTX_proj"
# result_dir <- paste(dir, "/G_CTX/figure360", sep = "")

Fig_dir <- file.path(result_dir, "CTX_proj_dev_figure")
# Check if the figure directory exists
if (!dir.exists(Fig_dir)) {
  # If it doesn't exist, create the directory
  dir.create(Fig_dir)
  print(paste("Directory", Fig_dir, "created."))
}

start_time <- Sys.time()

for (g_ind in 1:3) {

G_text=paste("G",g_ind,"_Projection", sep = "")

Data <- paste(dir, "HCD_LS_2.0_subject_completeness_bh.csv", sep = "/") %>% read.csv()
head_motion <- paste(dir, "head_mv_par_all.csv", sep = "/") %>% read.csv()
Measures_L <- paste(result_dir, "/indiv360/G", g_ind, "_CTX_all_L", ".csv", sep = "") %>% read.csv(header = FALSE)
Measures_R <- paste(result_dir, "/indiv360/G", g_ind, "_CTX_all_R", ".csv", sep = "") %>% read.csv(header = FALSE)
yeo7_label <- paste(Surf_temp, "/cortex_surf/glasser360_7networks.csv", sep = "") %>% read.csv(header = FALSE)

Measures <- rbind(Measures_L,Measures_R)
Data <- rbind(Data,Data)
Data$mFD <- head_motion$mFD

L <- rep(0,652)
R <- rep(1,652)
Data$hemi <- c(L,R)

Data$hemi <- as.factor(Data$hemi)
Data$Age <- as.numeric(Data$Age)
Data$Sex <- as.factor(Data$Sex)
Data$site <- as.factor(Data$site)
Data$mFD<- as.numeric(Data$mFD)

hemi=Data$hemi
Age=Data$Age
Sex=Data$Sex
mFD=Data$mFD
site=Data$site

y_pred <- c()
p_value<-matrix(NA,nrow=Region,ncol=1)
F_value<-matrix(NA,nrow=Region,ncol=1)
delta_R_sq<-matrix(NA,nrow=Region,ncol=1)
anova_pvalue<-matrix(NA,nrow=Region,ncol=1)

gam_estimated <- data.frame(matrix(data=NA, ncol=7))
gam_estimated <-setNames(gam_estimated,c("age","est","index","Yeo7_ind","Yeo7_name","age_deltaR2","derivative")) 

np=200
smooth_var="Age"

for (i in 1:Region) {
  Data$Measure <- Measures[,i]
  Data$Measure <- as.numeric(Data$Measure)
  Measure=Data$Measure
  
  mod_gam <- gam(Measure ~ s(Age, bs="cs", k=4) + Sex +site + #+ s(subj, bs="re")
                   s(hemi, bs="re") + mFD , data=Data, method="REML", na.action="na.omit") #, na.action="na.omit"
  
  #GAM derivatives
  #Get derivatives of the smooth function using finite differences
  derv <- derivatives(mod_gam, term = "s(Age)", n=np, interval = "simultaneous", unconditional = F) #derivative at 200 indices of smooth_var with a simultaneous CI
  #Identify derivative significance window(s)
  derv <- derv %>% #add "sig" column (TRUE/FALSE) to derv
    mutate(sig = !(0 > lower & 0 < upper)) #derivative is sig if the lower CI is not < 0 while the upper CI is > 0 (i.e., when the CI does not include 0)
  derv$sig_deriv = derv$derivative*derv$sig #add "sig_deriv derivatives column where non-significant derivatives are set to
  
  #get p_value and F_value of the s(Age)
  summary_model <- summary(mod_gam)
  p_value[i,] <- summary_model[["s.table"]]["s(Age)","p-value"]
  F_value[i,] <- summary_model[["s.table"]]["s(Age)","F"]
  R_sq <- summary_model$r.sq
  
  #null model without s(Age)
  mod_gam0 <- gam(Measure ~ Sex +site + #+ s(subj, bs="re")
                    s(hemi, bs="re") + mFD , data=Data, method="REML", na.action="na.omit") #, na.action="na.omit"
  summary_model0 <- summary(mod_gam0)
  R_sq0 <- summary_model0$r.sq
  
  delta_R_sq[i,] <- summary_model$r.sq-summary_model0$r.sq
  ### effect direction
  mean.derivative <- mean(derv$derivative)
  # if(mean.derivative < 0){ #if the average derivative is less than 0, make the effect size estimate negative
  #   delta_R_sq[i,] <- delta_R_sq[i,]*-1}

  anova_pvalue[i,] <- anova.gam(mod_gam0,mod_gam,test='Chisq')$`Pr(>Chi)`[2]
  
  #age_term = predict(mod_gam, type = 'terms')[, 's(Age)']
  #res_partial = residuals(mod_gam) + age_term
  #Age = mod_gam$model$Age
  #y_pred = cbind(y_pred,age_term[order(Age)])
  
  # Age = mod_gam$model$Age
  # Age_pred <- data.frame(Age=Age[order(Age)])
  # Age_pred <- fitted_values(mod_gam,
  #                           data=Age_pred,
  #                           #se.fit = TRUE,
  #                           #type = "link",
  #                           exclude = c("Sex", "s(hemi)","mFD","site"))
  # y_pred = cbind(y_pred,Age_pred$`fitted`)
  
  thisPred <- data.frame(init = rep(0,np))
  theseVars <- attr(mod_gam$terms,"term.labels") #gam model predictors (smooth_var + covariates)
  varClasses <- attr(mod_gam$terms,"dataClasses") #classes of the model predictors and y measure
  thisResp <- as.character(mod_gam$terms[[2]]) #the measure to predict
  for (v in c(1:length(theseVars))) { #fill the prediction df with data for predictions. These data will be used to predict the output measure (y) at np increments of the smooth_var, holding other model terms constant
    thisVar <- theseVars[[v]]
    thisClass <- varClasses[thisVar]
    if (thisVar == smooth_var) { 
      thisPred[,smooth_var] = seq(min(mod_gam$model[,smooth_var],na.rm = T),max(mod_gam$model[,smooth_var],na.rm = T), length.out = np) #generate a range of np data points, from minimum of smooth term to maximum of smooth term
    } else {
      switch (thisClass,
              "numeric" = {thisPred[,thisVar] = median(mod_gam$model[,thisVar])}, #make predictions based on median value
              "factor" = {thisPred[,thisVar] = levels(mod_gam$model[,thisVar])[[1]]}, #make predictions based on first level of factor 
              "ordered" = {thisPred[,thisVar] = levels(mod_gam$model[,thisVar])[[1]]} #make predictions based on first level of ordinal variable
      )
    }
  }
  pred <- thisPred %>% select(-init)
  #Estimate the smooth trajectory 
  prediction <- predict(object = mod_gam, newdata = pred, exclude = c("Sex", "s(hemi)","mFD","site"))
  
  GAM_ESTIMATES <- data.frame(matrix(data=NA, nrow=np, ncol=7))
  GAM_ESTIMATES <-setNames(GAM_ESTIMATES,c("age","est","index","Yeo7_ind","Yeo7_name","age_deltaR2","derivative"))
  GAM_ESTIMATES$age=pred$Age
  GAM_ESTIMATES$est=prediction
  GAM_ESTIMATES$index=i
  GAM_ESTIMATES$Yeo7_ind=yeo7_label$V3[i]
  GAM_ESTIMATES$Yeo7_name=yeo7_label$V4[i]
  GAM_ESTIMATES$age_deltaR2=delta_R_sq[i,]
  GAM_ESTIMATES$derivative=derv$derivative
  gam_estimated <- rbind(gam_estimated, GAM_ESTIMATES)
}

gam_estimated <- gam_estimated[-1,]

pvalue_fdr=p.adjust(p_value,method = "fdr")
anova_pvalue_fdr=p.adjust(anova_pvalue,method = "fdr")

GAM_Result<-cbind(F_value,p_value,delta_R_sq,anova_pvalue,anova_pvalue_fdr,pvalue_fdr)
GAM_Result<-as.data.frame(GAM_Result)
colnames(GAM_Result) <- c("age_Fvalue", #GAM F-value for the age smooth term
                          "age_pvalue", #GAM p-value for the age smooth term
                          "age_deltaR2", #delta adjusted Rsq from age and age-null models
                          "Anova_age_pvalue", #Anova p-value comparing age and age-null models
                          "Anova_age_adjpvalue", #fdr adjusted Anova p-value comparing age and age-null models
                          "age_adjpvalue") #fdr adjusted Gam age p-value
filepath<-paste(Fig_dir, "/G", g_ind, "_CTX_proj_dev.csv", sep = "")
write.csv(GAM_Result, filepath, row.names = F, quote = F)

filepath<-paste(Fig_dir, "/G", g_ind, "_gam_estimated.csv", sep = "")
write.csv(gam_estimated, filepath, row.names = F, quote = F)


## plot the cortical distribution of age_deltaR2
# ggseg(.data = GAM_Result, atlas = "glasser", mapping=aes(fill = age_deltaR2), position = c("stacked")) +
#   theme_void() +
#   labs(fill="") +
#   paletteer::scale_fill_paletteer_c("pals::ocean.matter", direction = -1, limits = c(-0.02, .03)) +
#   theme(legend.text = element_text(color = c("black")))

## plot region specific developmental trajectories (colored by age_deltaR2)
# gam_estimated_sort <- gam_estimated[order(gam_estimated$age_deltaR2),] #, decreasing = TRUE
myplot <- ggplot(gam_estimated,aes(age,est,group=age_deltaR2,color=age_deltaR2)) + 
  geom_line(size=.3, alpha = .8) + 
  scale_color_viridis(option = "D") +
  theme_classic() +
  labs(x = "\nAge", y = parse(text = G_text) ) +
  theme(legend.position = "none") +
  theme(axis.text = element_text(size=6, color = c("black")), axis.title = element_text(size=6, color = c("black")), axis.line = element_line(size=.22), axis.ticks = element_line(size=.22)) +
  scale_x_continuous(breaks=c(6, 8, 10, 12, 14, 16, 18, 20, 22), expand = c(0,.45)) 

filepath1<-paste(Fig_dir, "/G", g_ind, "_CTX_proj_dev.pdf", sep = "")
ggsave(filename = filepath1, device = "pdf", myplot, dpi = 500, width = 2.55 , height = 2.3)
filepath2<-paste(Fig_dir, "/G", g_ind, "_CTX_proj_dev.tiff", sep = "")
ggsave(filename = filepath2, device = "tiff", myplot, dpi = 500, width = 2.55 , height = 2.3)

# ## plot region specific developmental trajectories (colored by yeo7)
# myplot <- ggplot(gam_estimated,aes(age,est,group=age_deltaR2,color=as.factor(Yeo7_ind))) + 
#   geom_line(size=.75, alpha = .8, show.legend = TRUE) + 
#   scale_color_discrete(type = viridis(7),labels=c("Visual","Somatomotor","Dorsal Attention","Ventral Attention","Limbic","Frontoparietal","Default")) +
#   theme_classic() +
#   labs(x = "\nAge", y = parse(text = G_text) ) +
#   theme(legend.key.size = unit(0.5, 'cm'), legend.key.width = unit(0.4, "cm"),legend.title = element_text(size=6),legend.text = element_text(size=6)) +
#   theme(axis.text = element_text(size=6, color = c("black")), axis.title = element_text(size=6, color = c("black")), axis.line = element_line(size=.22), axis.ticks = element_line(size=.22)) +
#   scale_x_continuous(breaks=c(6, 8, 10, 12, 14, 16, 18, 20, 22), expand = c(0,.45))
# 
# filepath1<-paste(result_dir, "/ctx_dev_LR/G", g_ind, "_CTX_proj_dev_yeo7.pdf", sep = "")
# ggsave(filename = filepath1, device = "pdf", myplot, dpi = 500, width = 5.5 , height = 4)
# filepath2<-paste(result_dir, "/ctx_dev_LR/G", g_ind, "_CTX_proj_dev_yeo7.tiff", sep = "")
# ggsave(filename = filepath2, device = "tiff", myplot, dpi = 500, width = 5.5 , height = 4)
}

end_time <- Sys.time()
execution_time <- end_time - start_time
print(paste("Execution time: ", execution_time))