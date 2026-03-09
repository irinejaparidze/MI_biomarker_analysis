
library(ggplot2)
library(tidyverse)
library(tidyr)
library(dplyr)
library(rstatix)
library(ggpubr)
library(patchwork)
library(ROCit)
library(ggplotify)
library(xtable)


### STEMI - No Damage ROC curves###

DataA <- read.csv("/home/irina/Documents/Bioinformatics/Giorgi Cardiology Docs/AllGroupsForCardiology.csv")

DataA
################################################
#Crerated new version to use initial table later
################################################

DataB <- DataA
view(DataB)

#####################################################
#Getting all proteins and concentrations in one table
#####################################################

Data_Long <- DataB %>%
  pivot_longer(cols = c(Ubqng.ml, CXCR4ng.ml, SDF1ng.ml),
               names_to = "Protein",
               values_to = "Concentration"
  )

Data_Long

#####################
# Plotting histograms
#####################

ggplot(Data_Long, aes(x = Concentration)) +
  
  geom_histogram(aes(y = ..density..), bins = 5, fill = "lightblue", color = "black") +
  
  geom_density(alpha = .2, color = "black") +
 
  facet_wrap(~Protein, ncol = 1, scales = "free") +
  theme_classic() +
  labs(title = "Histogram and density distribution of proteins",
       y = "density") +
  theme(strip.background = element_blank(), 
        strip.text = element_text(hjust = 0, face = "bold"))

####################################
#SHAPIRO TEST TO CHECK THE NORMALITY
####################################

test1 <- shapiro.test(DataB$Ubqng.ml)
test1

test2 <- shapiro.test(DataB$CXCR4ng.ml)
test2

test3 <- shapiro.test(DataB$SDF1ng.ml)
test3

#####################################################################
#Changing NoDamage on Healthy (just to make clear version for artice)
#####################################################################

DataB$Diagnosis <- recode(DataB$Diagnosis,
                          "NoDamage"="Healthy")
DataB

##################################################
#### Subsetting Data to compare Stemi and  Healthy. Factorizing variables
##################################################

DataC<- subset(DataB, Diagnosis != "NSTEMI") #STEMI-Healthy

DataC$Diagnosis <- as.factor(DataC$Diagnosis)

DataE<- subset(DataB, Diagnosis != "STEMI") # NSTENI - Healthy

DataE$Diagnosis <- as.factor(DataE$Diagnosis)

DataC
DataE

str(DataE)
str(DataC)


#######################################################
#Creating tests
#######################################################

##############################################
#Wilcoxon test for STEMI preparation, pivoting
##############################################

#STEMI - Healthy

DC <- DataC[, 2:5]

df_longC <- DC %>%
  mutate(obs_id = row_number()) %>% 
  pivot_longer(
    cols = c(Ubqng.ml, CXCR4ng.ml, SDF1ng.ml),
    names_to = "Protein",
    values_to = "Concentration"
  )

df_longC$Protein <- recode(df_longC$Protein,
                           "Ubqng.ml"="Ubiquitin",
                           "CXCR4ng.ml"="CXCR4",
                           "SDF1ng.ml"="SDF1")


head(df_longC)
df_longC %>% group_by(Diagnosis, Protein) %>% summarize(Count = n())
tail(df_longC)


stat.test1 <-df_longC %>%
  group_by(Protein) %>%
  wilcox_test(Concentration ~ Diagnosis) %>%
  adjust_pvalue(method = "none")%>%
  add_significance()

stat.test1

stat.test2 <-df_longC %>%
  group_by(Protein) %>%
  wilcox_effsize(Concentration ~ Diagnosis) %>%
  adjust_pvalue(method = "none")%>%
  add_significance()

stat.test2

#######################
#Saving the results

write_csv(stat.test1, 
          "/home/irina/Documents/Bioinformatics/Giorgi Cardiology Docs/ROC files/Wilcoxon_STEMI.csv")



################################################
#Wilcoxon test for NSTEMI preparation, pivoting
################################################

DE <- DataE[, 2:5]

df_longE <- DE %>%
  mutate(obs_id = row_number()) %>% 
  pivot_longer(
    cols = c(Ubqng.ml, CXCR4ng.ml, SDF1ng.ml),
    names_to = "Protein",
    values_to = "Concentration"
  )

df_longE$Protein <- recode(df_longE$Protein,
                           "Ubqng.ml"="Ubiquitin",
                           "CXCR4ng.ml"="CXCR4",
                           "SDF1ng.ml"="SDF1")


head(df_longE)
df_longE %>% group_by(Diagnosis, Protein) %>% summarize(Count = n())
tail(df_longE)



stat.test1n <-df_longE %>%
  group_by(Protein) %>%
  wilcox_test(Concentration ~ Diagnosis) %>%
  adjust_pvalue(method = "none")%>%
  add_significance()

stat.test1n

stat.test2n <-df_longE %>%
  group_by(Protein) %>%
  wilcox_effsize(Concentration ~ Diagnosis) %>%
  adjust_pvalue(method = "none")%>%
  add_significance()

stat.test2n

##################
#Saving the results

write_csv(stat.test1n, 
          "/home/irina/Documents/Bioinformatics/Giorgi Cardiology Docs/ROC files/Wilcoxon_NSTEMI.csv")

#########################
#Wilcpoxon Visualizations
#########################

#############
#For STEMI
#############

stat.test1 <- df_longC %>%
  group_by(Protein) %>%
  wilcox_test(Concentration ~ Diagnosis) %>%
  adjust_pvalue(method = "none") %>%
  add_significance()

stat.test1$statistic <- round(stat.test1$statistic, 2)

colnames(stat.test1)[colnames(stat.test1) == "p"] <- "p"
colnames(stat.test1)[colnames(stat.test1) == "p.signif"] <- "p.signif"

stat.test1

stats <- stat.test1 %>%
  add_xy_position(
    x = "Diagnosis",
    fun = "max"
  )


my_bxp1 <- ggboxplot(
  df_longC,
  x = "Diagnosis",
  y = "Concentration",
  fill = "Diagnosis",
  add = "dotplot",
  palette = c("lightgreen", "#D55E00")
) +
  facet_wrap(~Protein, scales = "free_y") +
  stat_pvalue_manual(
    stats,
    label = "W = {statistic}, adj.p = {p.adj} {p.adj.signif}",
    tip.length = 0.02,
    size = 6
  ) +
  labs(
    x = "Diagnosis",
    y = "Concentration (ng/ml)",
    title = "a",
    subtitle = "Wilcoxon Test Comparison Between STEMI and Healthy Groups"
  ) +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 16, color = "darkblue"),
    axis.text = element_text(size = 14, color = "darkred", face = "bold"),
    plot.title = element_text(size = 18, face = "bold")
  )

my_bxp1

#############
#For NSTEMI
#############


stat.test1e <- df_longE %>%
  group_by(Protein) %>%
  wilcox_test(Concentration ~ Diagnosis) %>%
  adjust_pvalue(method = "none") %>%
  add_significance()

stat.test1e$statistic <- round(stat.test1e$statistic, 2)

colnames(stat.test1e)[colnames(stat.test1e) == "p"] <- "p"
colnames(stat.test1e)[colnames(stat.test1e) == "p.signif"] <- "p.signif"

stat.test1e

stats <- stat.test1e %>%
  add_xy_position(
    x = "Diagnosis",
    fun = "max"
  )


my_bxp1e <- ggboxplot(
  df_longE,
  x = "Diagnosis",
  y = "Concentration",
  fill = "Diagnosis",
  add = "dotplot",
  palette = c("lightgreen", "#D55E00")
) +
  facet_wrap(~Protein, scales = "free_y") +
  stat_pvalue_manual(
    stats,
    label = "W = {statistic}, adj.p = {p.adj} {p.adj.signif}",
    tip.length = 0.02,
    size = 6
  ) +
  labs(
    x = "Diagnosis",
    y = "Concentration (ng/ml)",
    title = "b",
    subtitle = "Wilcoxon Test Comparison Between NSTEMI and Healthy Groups"
  ) +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 16, color = "darkblue"),
    axis.text = element_text(size = 14, color = "darkred", face = "bold"),
    plot.title = element_text(size = 18, face = "bold")
  )

my_bxp1e


my_bxp1 <- my_bxp1 +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
  coord_cartesian(clip = "off")

my_bxp1e <- my_bxp1e +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
  coord_cartesian(clip = "off")

#############################
#Combining plots all together


combined_plot <- my_bxp1 / my_bxp1e &
  theme(plot.margin = margin(20, 60, 20, 60))

combined_plot

################################
#Saving plots in high resolution
ggsave(
  "/home/irina/Documents/Bioinformatics/Giorgi Cardiology Docs/ROC files/Protein_wilcoxon_comparison.tiff",
  combined_plot,
  width = 14,
  height = 10,
  dpi = 600
)

########################################
#Starting ROC calculations
########################################

summary(DataC)

######################
#Preping data

DataUbi <-  DataC[, 2:3]
DataCX <- DataC[, c(2, 4)]
DataSD <- DataC[, c(2, 5)]
DataTro <- DataC[, c(2, 6)]

DataUbi 
DataCX 
DataSD 
DataTro 


#################################################################
# For STEMI
#################################################################
#Saving Data in csv format
##########################

DataUbi$Health_Status <- ifelse(DataUbi$Diagnosis == 'STEMI', '+', '-')
DataUbi
write.csv(DataUbi, 'articleUbi_stemi.csv')

DataCX$Health_Status <- ifelse(DataCX$Diagnosis == 'STEMI', '+', '-')
DataCX
write.csv(DataCX, 'articleCX_stemi.csv')


DataSD$Health_Status <- ifelse(DataSD$Diagnosis == 'STEMI', '+', '-')
DataSD
write.csv(DataSD, 'articleSD_stemi.csv')

DataTro$Health_Status <- ifelse(DataTro$Diagnosis == 'STEMI', '+', '-')
DataTro
write.csv(DataTro, 'articleTRO_stemi.csv')


################################################################################
#STEMI-UBI
################################################################################

logistic.model_Ubi <- glm(as.factor(Diagnosis)~Ubqng.ml,
                      data = DataUbi,family = "binomial")


classUbi <- logistic.model_Ubi$y

scoreUbi <- logistic.model_Ubi$fitted.values

measureUBI <- measureit(score = scoreUbi, class = classUbi,
                     measure = c("ACC", "SENS", 'SPEC' ,"FSCR"), method = "non")
measureUBI

Cutoff_Ubi <-  cbind(Cutoff = measureUBI$Cutoff, Depth = measureUBI$Depth,
                     TP = measureUBI$TP, FP = measureUBI$FP, 
                     TN = measureUBI$TN, FN = measureUBI$FN,
                     ACC = measureUBI$ACC, SENS = measureUBI$SENS, 
                     SPEC = measureUBI$SPEC,
                     FSCR = measureUBI$FSCR)
Cutoff_Ubi

#Saving Cutoff data in csv

write.csv(Cutoff_Ubi, '/home/irina/Documents/Bioinformatics/Giorgi Cardiology Docs/ROC files/Cutpoints_ubi_STEMI.csv')





#ROC plots for ubi

summary(as.factor(DataUbi$Diagnosis))

roc_empirical <- rocit(score = DataUbi$Ubqng.ml, class = DataUbi$Health_Status,
                       negref = "-") 

roc_nonparametric_UBI <- rocit(score = DataUbi$Ubqng.ml, 
                               class = DataUbi$Health_Status,
                               negref = "-", 
                               method = "non") 
roc_nonparametric_UBI
summary(roc_nonparametric_UBI)



plot(roc_nonparametric_UBI, 
     legend = TRUE,
     legendpos = "bottomright",
     YIndex = TRUE,
     values = TRUE,
     main ="ROC Ubiquitin, AUC = 0.8891")+
  title("ROC Ubiquitin, AUC = 0.8891")



###################
######STEMI - CXCR4
###################

logistic.model_CX <- glm(as.factor(Diagnosis)~CXCR4ng.ml,
                          data = DataCX,family = "binomial", )
classCX <- logistic.model_CX$y
scoreCX <- logistic.model_CX$fitted.values



measureCX <- measureit(score = scoreCX, class = classCX,
                        measure = c("ACC", "SENS", "SENS", 'SPEC' ,"FSCR"), method = "non")
measureCX

Cutoff_CX <-  cbind(Cutoff = measureCX$Cutoff, Depth = measureCX$Depth,
                     TP = measureCX$TP, FP = measureCX$FP, 
                     TN = measureCX$TN, FN = measureCX$FN,
                     ACC = measureCX$ACC, SENS = measureCX$SENS, 
                     SPEC = measureCX$SPEC,
                     FSCR = measureCX$FSCR)

write.csv(Cutoff_CX, '/home/irina/Documents/Bioinformatics/Giorgi Cardiology Docs/ROC files/Cutpoints_CXCX4_STEMI.csv')


roc_nonparametric_CX <- rocit(score = DataCX$CXCR4ng.ml, 
                              class = DataCX$Health_Status,
                              negref = "-", 
                              method = "non") 
summary(roc_nonparametric_CX)


plot(roc_nonparametric_CX, 
     legend = TRUE,
     legendpos = "bottomright",
     YIndex = TRUE,
     values = TRUE) +
  title("ROC CXCR4, AUC = 0.8106")

###############
#####SDF1-STEMI
###############

logistic.model_SD <- glm(as.factor(Diagnosis)~SDF1ng.ml,
                         data = DataSD,family = "binomial")
classSD <- logistic.model_SD$y
scoreSD <- logistic.model_SD$fitted.values



measureSD <- measureit(score = scoreSD, class = classSD,
                       measure = c("ACC", "SENS", "SENS", 'SPEC' ,"FSCR"), method = "non")
measureSD

Cutoff_SD <-  cbind(Cutoff = measureSD$Cutoff, Depth = measureSD$Depth,
                    TP = measureSD$TP, FP = measureSD$FP, 
                    TN = measureSD$TN, FN = measureSD$FN,
                    ACC = measureSD$ACC, SENS = measureSD$SENS, 
                    SPEC = measureSD$SPEC,
                    FSCR = measureSD$FSCR)

write.csv(Cutoff_SD, '/home/irina/Documents/Bioinformatics/Giorgi Cardiology Docs/ROC files/Cutpoints_SDF1_STEMI.csv')


roc_nonparametric_SD <- rocit(score = DataSD$SDF1ng.ml, 
                              class = DataSD$Health_Status,
                              negref = "-", 
                              method = "non") 
summary(roc_nonparametric_SD)


plot(roc_nonparametric_SD, 
     legend = TRUE,
     legendpos = "bottomright",
     YIndex = TRUE,
     values = TRUE) +
  title("ROC SDF1, AUC = 0.8051")

##################
###Toponin_STEMI
#################

DataTro

logistic.model_TR <- glm(as.factor(Health_Status)~Troponing.ml,
                         data = DataTro, family = "binomial")

classTR <- logistic.model_TR$y
scoreTR <- logistic.model_TR$fitted.values



measureTR <- measureit(score = scoreTR, class = classTR,
                       measure = c("ACC", "SENS",'SPEC' ,"FSCR"))
measureTR

Cutoff_TR <-  cbind(Cutoff = measureTR$Cutoff, Depth = measureTR$Depth,
                    TP = measureTR$TP, FP = measureTR$FP, 
                    TN = measureTR$TN, FN = measureTR$FN,
                    ACC = measureTR$ACC, SENS = measureTR$SENS, 
                    SPEC = measureTR$SPEC,
                    FSCR = measureTR$FSCR)

write.csv(Cutoff_TR, '/home/irina/Documents/Bioinformatics/Giorgi Cardiology Docs/ROC files/Cutpoints_TR_STEMI.csv')


roc_nonparametric_TR <- rocit(score = DataTro$Troponing.ml, 
                              class = DataTro$Health_Status,
                              negref = "-") 
summary(roc_nonparametric_TR)


plot(roc_nonparametric_TR, 
     legend = TRUE,
     legendpos = "bottomright",
     YIndex = TRUE,
     values = TRUE) +
  title("ROC Troponin, AUC = 0.0.8619")

##############################################################################
#NSTEMI- NoDamage analysis
##############################################################################


DataAa <- DataA
DataAa
DataBa <- DataAa[, 2:6]
DataBa


#### Subsetting Data to compare Stemi and  no damage
#

DataCa<- subset(DataBa, Diagnosis != "STEMI")

DataCa$Diagnosis <- as.factor(DataCa$Diagnosis)


DataCa <- DataCa%>% 
  arrange(Diagnosis)
DataCa

write.csv(DataCa, 'NSTEMI.csv')

summary(DataCa)


DataUbia <-  DataCa[, 1:2]
DataCXa <- DataCa[, c(1, 3)]
DataSDa <- DataCa[, c(1, 4)]
DataTroa <- DataCa[, c(1, 5)]

DataCXa 
DataSDa 
DataTroa
DataUbia


#######################
#Subsetting Data
#######################

DataUbia$Health_Status <- ifelse(DataUbia$Diagnosis == 'NSTEMI', '+', '-')
DataUbia
write.csv(DataUbia, '/home/irina/Documents/Bioinformatics/Giorgi Cardiology Docs/ROC files/Ubi_Nstemi.csv')


DataCXa$Health_Status <- ifelse(DataCXa$Diagnosis == 'NSTEMI', '+', '-')
DataCXa
write.csv(DataCXa, '/home/irina/Documents/Bioinformatics/Giorgi Cardiology Docs/ROC files/CXCR4_Nstemi.csv')


DataSDa$Health_Status <- ifelse(DataSDa$Diagnosis == 'NSTEMI', '+', '-')
DataSDa
write.csv(DataSDa, '/home/irina/Documents/Bioinformatics/Giorgi Cardiology Docs/ROC files/SDf_Nstemi.csv')

DataTroa$Health_Status <- ifelse(DataTroa$Diagnosis == 'NSTEMI', '+', '-')
DataTroa
write.csv(DataTroa, '/home/irina/Documents/Bioinformatics/Giorgi Cardiology Docs/ROC files/articleTRO_Nstemi.csv')


####################
####Ubiquitin-NSTEMI
####################

logistic.model_Ubia <- glm(as.factor(Health_Status)~Ubqng.ml,
                          data = DataUbia,family = "binomial")



classUbia <- logistic.model_Ubia$y

scoreUbia <- logistic.model_Ubia$fitted.values

measureUBIa <- measureit(score = scoreUbia, class = classUbia,
                        measure = c("ACC", "SENS", 'SPEC' ,"FSCR"))
measureUBIa



Cutoff_Ubia <-  cbind(Cutoff = measureUBIa$Cutoff, Depth = measureUBIa$Depth,
                     TP = measureUBIa$TP, FP = measureUBIa$FP, 
                     TN = measureUBIa$TN, FN = measureUBIa$FN,
                     ACC = measureUBIa$ACC, SENS = measureUBIa$SENS, 
                     SPEC = measureUBIa$SPEC,
                     FSCR = measureUBIa$FSCR)

Cutoff_Ubia 

write.csv(Cutoff_Ubia, '/home/irina/Documents/Bioinformatics/Giorgi Cardiology Docs/ROC files/Cutpoints_ubi_nstemi_1.csv')




##################
#ROC plots for ubi
##################

summary(as.factor(DataUbia$Health_Status))

roc_empiricala <- rocit(score = DataUbia$Ubqng.ml, class = DataUbia$Health_Status,
                       negref = "-") 

roc_nonparametric_UBIa <- rocit(score = DataUbia$Ubqng.ml, 
                               class = DataUbia$Health_Status,
                               negref = "-", 
                               method = "non") 
summary(roc_nonparametric_UBIa)



a = cbind(Cutoff=roc_nonparametric_UBIa$Cutoff, 
      TPR=roc_nonparametric_UBIa$TPR, 
      FPR=roc_nonparametric_UBIa$FPR)
a


roc_nonparametric_UBIa

plot(roc_nonparametric_UBIa, 
     legend = TRUE,
     legendpos = "bottomright",
     YIndex = TRUE,
     values = TRUE) +
  title("ROC Ubiquitin, AUC = 0.9055")


##################
######CXCR4 NSTEMI
##################

logistic.model_CXa <- glm(as.factor(Health_Status)~CXCR4ng.ml,
                          data = DataCXa,family = "binomial", )
classCXa <- logistic.model_CXa$y
scoreCXa <- logistic.model_CXa$fitted.values


measureCXa <- measureit(score = scoreCXa, class = classCXa,
                       measure = c("ACC", "SENS", "SENS", 'SPEC' ,"FSCR"), method = "non")
measureCXa

Cutoff_CXa <-  cbind(Cutoff = measureCXa$Cutoff, Depth = measureCXa$Depth,
                    TP = measureCXa$TP, FP = measureCXa$FP, 
                    TN = measureCXa$TN, FN = measureCXa$FN,
                    ACC = measureCXa$ACC, SENS = measureCXa$SENS, 
                    SPEC = measureCXa$SPEC,
                    FSCR = measureCXa$FSCR)

write.csv(Cutoff_CXa, '/home/irina/Documents/Bioinformatics/Giorgi Cardiology Docs/ROC files/Cutpoints_CXCR4_nstemi.csv')

DataCXa

roc_nonparametric_CXa <- rocit(score = DataCXa$CXCR4ng.ml, 
                              class = DataCXa$Health_Status,
                              negref = "-", 
                              method = "non") 
summary(roc_nonparametric_CXa)


plot(roc_nonparametric_CXa, 
     legend = TRUE,
     legendpos = "bottomright",
     YIndex = TRUE,
     values = TRUE) + 
  title("ROC CXCR4, AUC = 0.6978 ")

################
#####SDF1-NSTEMI
################

logistic.model_SDa <- glm(as.factor(Health_Status)~SDF1ng.ml,
                         data = DataSDa,family = "binomial")
classSDa <- logistic.model_SDa$y
scoreSDa <- logistic.model_SDa$fitted.values

DataSDa

measureSDa <- measureit(score = scoreSDa, class = classSDa,
                       measure = c("ACC", "SENS", "SENS", 'SPEC' ,"FSCR"), method = "non")
measureSDa

Cutoff_SDa <-  cbind(Cutoff = measureSDa$Cutoff, Depth = measureSDa$Depth,
                    TP = measureSDa$TP, FP = measureSDa$FP, 
                    TN = measureSDa$TN, FN = measureSDa$FN,
                    ACC = measureSDa$ACC, SENS = measureSDa$SENS, 
                    SPEC = measureSDa$SPEC,
                    FSCR = measureSDa$FSCR)

write.csv(Cutoff_SDa, '/home/irina/Documents/Bioinformatics/Giorgi Cardiology Docs/ROC files/Cutpoints_SDF1_Nstemi.csv')

Cutoff_SDa

roc_nonparametric_SDa <- rocit(score = DataSDa$SDF1ng.ml, 
                              class = DataSDa$Health_Status,
                              negref = "-", 
                              method = "non") 
roc_nonparametric_SDa

summary(roc_nonparametric_SDa)


plot(roc_nonparametric_SDa, 
     legend = TRUE,
     legendpos = "bottomright",
     YIndex = TRUE,
     values = TRUE) +
title("ROC SDF1, AUC = 0.5618")


#################
#NSTEMI troponin
################

DataTroa

logistic.model_TRa <- glm(as.factor(Health_Status)~Troponing.ml,
                         data = DataTroa, family = "binomial")

classTRa <- logistic.model_TRa$y
scoreTRa <- logistic.model_TRa$fitted.values



measureTRa <- measureit(score = scoreTRa, class = classTRa,
                       measure = c("ACC", "SENS",'SPEC' ,"FSCR"))
measureTRa

Cutoff_TRa <-  cbind(Cutoff = measureTRa$Cutoff, Depth = measureTRa$Depth,
                    TP = measureTRa$TP, FP = measureTRa$FP, 
                    TN = measureTRa$TN, FN = measureTRa$FN,
                    ACC = measureTRa$ACC, SENS = measureTRa$SENS, 
                    SPEC = measureTRa$SPEC,
                    FSCR = measureTRa$FSCR)

write.csv(Cutoff_TRa, '/home/irina/Documents/Bioinformatics/Giorgi Cardiology Docs/ROC files/Cutpoints_TR_NSTEMI.csv')


roc_nonparametric_TRa <- rocit(score = DataTroa$Troponing.ml, 
                              class = DataTroa$Health_Status,
                              negref = "-"
                             ) 
summary(roc_nonparametric_TRa)


plot(roc_nonparametric_TRa, 
     legend = TRUE,
     legendpos = "bottomright",
     YIndex = TRUE,
     values = TRUE) +
  title("NSTEMI-Healthy,ROC Troponin, AUC = 1")


#######################################################
#Common final plots
# High resolution journal figure
#######################################################


#Ubiquitin
library(ggplotify)

ubi_plot <-  as.ggplot(function() {
  plot(roc_nonparametric_UBI, col = c(1), 
     main = "a",
     legend = FALSE, YIndex = FALSE)
  lines(roc_nonparametric_UBIa$TPR~roc_nonparametric_UBIa$FPR, 
      col = 2, lwd = 2)
  legend("bottomright", 
         legend = c("STEMI - Healthy ROC, Ubiquitin, AUC=0.8891", "NSTEMI - Healthy ROC, Ubiquitin, AUC=0.9055"),
         col = c(1,2), lwd = 1.6)
})


#CXCR4
CXCR4_plot <- as.ggplot(function() {
  plot(roc_nonparametric_CX, col = c(1), 
     main ="b",
     legend = FALSE, YIndex = FALSE)
  lines(roc_nonparametric_CXa$TPR~roc_nonparametric_CXa$FPR, 
      col = 2, lwd = 2)
  legend("bottomright", 
         legend = c("STEMI - Healthy ROC, CXCR4, AUC=0.8106", "NSTEMI - Healthy ROC, CXCR4, AUC=0.6978"),
         col = c(1,2), lwd = 1.6)
})


#SDF1

SDF1_plot <- as.ggplot(function() {
  plot(roc_nonparametric_SD, col = c(1),
     main = "c",
     legend = FALSE, YIndex = FALSE)
  lines(roc_nonparametric_SDa$TPR~roc_nonparametric_SDa$FPR, 
      col = 2, lwd = 2)
  legend("bottomright", 
         legend = c("STEMI - Healthy ROC, SDF1, AUC=0.8051", "NSTEMI - Healthy ROC, SDF1, AUC=0.5618"),
         col = c(1,2), lwd = 1.6)
})


#Troponin

Tro_plot <- as.ggplot(function() {
  plot(roc_nonparametric_TR, 
     col = c(1), 
     legend = FALSE, YIndex = FALSE,
     main = "d", )
  lines(roc_nonparametric_TRa$TPR~roc_nonparametric_TRa$FPR, 
      col = 2, lwd = 0.5)
  legend("bottomright", 
         legend = c("STEMI - Healthy ROC, Troponin, AUC=1", "NSTEMI - Healthy ROC, Troponin, AUC=1"),
         col = c(1,2), lwd = 1.6)
})




combined_plot <- (ubi_plot | CXCR4_plot) /
  (SDF1_plot | Tro_plot) +
  plot_annotation(tag_levels = "a")

combined_plot


ggsave(
  "/home/irina/Documents/Bioinformatics/Giorgi Cardiology Docs/ROC files/ROC_Combined.tiff",
  combined_plot,
  width = 14,
  height = 10,
  dpi = 600
)

