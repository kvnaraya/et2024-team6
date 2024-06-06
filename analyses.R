
############################ SETTING UP SHOP, READING IN DATA #################################

##Always start by clearing all existing variables.  
rm(list=ls(all=TRUE))

##Update packages
update.packages()

##Load relevant package_cs
library(reshape)
library(nlme)    #Non-linear mixed effects package
library(car)   #Regression package
library(lme4)    #Linear mixed effects package
library(ggplot2)   #Advanced plotting package
library(language_cR)  ##Package_cs containing useful functions for language research 
library(lattice)  ##Plotting package
library(lmerTest) 
library(dplyr)
library(tidyr)
library(readxl)
library(factoextra)
library(mice)
library(rstatix)
library(compareGroups)

##Read in data
setwd("/Users/keertananarayanan/Desktop/et2024")

#semantic data
SLdata1 = read.csv("SemanticData.csv", header = TRUE, sep = ",")
#cohort data
data1 = read.csv("CohortData.csv", header = TRUE, sep = ",")
#VW data 
VWdata = read.csv("VVWPdata.csv", header = TRUE, sep = ",")

SLdata <- na.omit(SLdata)
VWdata <- na.omit(VWdata)

colnames(SLdata1) <- paste("Semantic", colnames(SLdata1), sep = "_")
colnames(data1) <- paste("Cohort", colnames(data1), sep = "_")
colnames(VWdata) <- paste("VVWP", colnames(VWdata), sep = "_")
data<- cbind(SLdata1, data1)
data <- cbind(data, VWdata)

# transform target slope and crossover to get target timing per trial type

#Semantics 
data$Semantic_slope_z <-log10(data$Semantic_Slope_S)
data$Semantic_cross_z <-(data$Semantic_Crossover_S)    #don't need to log scale XO.
Semantic_meanslope  <- mean(data$Semantic_slope_z)
Semantic_meancross <- mean(data$Semantic_cross_z)
Semantic_sdslope  <- sd(data$Semantic_slope_z)
Semantic_sdcross <- sd(data$Semantic_cross_z)
data$Semantic_slope_z <- ((data$Semantic_slope_z)-(Semantic_meanslope))/(Semantic_sdslope)
data$Semantic_cross_z <- -1*(((data$Semantic_cross_z)-(Semantic_meancross))/(Semantic_sdcross))
data$Semantic_targettiming<- (data$Semantic_slope_z+data$Semantic_cross_z)/2

#Cohort
data$Cohort_slope_z <-log10(data$Cohort_Slope)
data$Cohort_cross_z <-(data$Cohort_Crossover)    #don't need to log scale XO.
Cohort_meanslope  <- mean(data$Cohort_slope_z)
Cohort_meancross <- mean(data$Cohort_cross_z)
Cohort_sdslope  <- sd(data$Cohort_slope_z)
Cohort_sdcross <- sd(data$Cohort_cross_z)
data$Cohort_slope_z <- ((data$Cohort_slope_z)-(Cohort_meanslope))/(Cohort_sdslope)
data$Cohort_cross_z <- -1*(((data$Cohort_cross_z)-(Cohort_meancross))/(Cohort_sdcross))
data$Cohort_targettiming<- (data$Cohort_slope_z+data$Cohort_cross_z)/2

#VVWP
data$VVWP_slope_z <-log10(data$VVWP_Slope)
data$VVWP_cross_z <-(data$VVWP_Crossover)    #don't need to log scale XO.
VVWP_meanslope  <- mean(na.omit(data$VVWP_slope_z))
VVWP_meancross <- mean(na.omit(data$VVWP_cross_z))
VVWP_sdslope  <- sd(na.omit(data$VVWP_slope_z))
VVWP_sdcross <- sd(na.omit(data$VVWP_cross_z))
data$VVWP_slope_z <- ((data$VVWP_slope_z)-(VVWP_meanslope))/(VVWP_sdslope)
data$VVWP_cross_z <- -1*(((data$VVWP_cross_z)-(VVWP_meancross))/(VVWP_sdcross))
data$VVWP_targettiming<- (data$VVWP_slope_z+data$VVWP_cross_z)/2

#target resolution
data$Semantic_resolution <- (data$Semantic_Max-((data$Semantic_CB2+data$Semantic_UB2)/2))
data$Cohort_resolution <- (data$Cohort_Max-((data$Cohort_CB2+data$Cohort_UB2)/2))
data$VVWP_resolution <- (data$VVWP_Max-((data$VVWP_CB2+data$VVWP_UB2)/2))

#impute missing data 
datax <- data[,c('Semantic_WJOC', 'Semantic_WJSR')]
dataxmissing <- mice(data=datax, meth=c("","pmm"))
dataxcomplete <- complete(dataxmissing,1)

#center continuous variables
data$WJOC_c <- data$Semantic_WJOC - mean(data$Semantic_WJOC)
data$WJSR_c <- data$Semantic_WJSR - mean(na.omit(data$Semantic_WJSR))
data$Cohort_age_c <- data$Cohort_age - mean(data$Cohort_age)
data$Semantic_age_c <- data$Semantic_age- mean(data$Semantic_age)
data$Cohort_acc_c <- data$Cohort_acc - mean(data$Cohort_acc)
data$Semantic_acc_c <- data$Semantic_acc- mean(data$Semantic_acc)
data$Cohort_rt_c <- data$Cohort_rt - mean(data$Cohort_rt)
data$Semantic_rt_c <- data$Semantic_rt- mean(data$Semantic_rt)
data$Semantic_gnrt_c <- data$Semantic_gnrt- mean(na.omit(data$Semantic_gnrt))
data$VVWP_eyemove_c <- data$VVWP_eyemove- mean(na.omit(data$VVWP_eyemove))
data$Semantic_ssrt_c <- data$Semantic_ssrt -mean(na.omit(data$Semantic_ssrt))
#Z scoring CC within age_c 
data <- group_by(data, Semantic_agegroup)
data$sscongruent_log<-log10(data$Semantic_sscongruent)
data$ssincongruent_log<- log10(data$Semantic_ssincongruent)
data$ssctotal_log<- log10(data$Semantic_sstotal)
data$ssscore <- (data$ssincongruent_log)-(data$sscongruent_log)
data$ssscore <- data$ssscore - mean(na.omit(data$ssscore))  #MAYBE DIVIDE BY SD within age_c?

# resdiualize CC within age 

#global stop
data4<- subset(data, !is.na(Cohort_sscongruent),)
resid <- lm(Cohort_globalstop ~ Cohort_age, data = data4)
   			resid.res <- (resid(resid))
   			resid.res <- as.data.frame(resid.res)
   			data4 <- cbind(data4, resid.res)
   			data4$Cohort_gs_res <- data4$resid.res
#responseinhib
resid2<- lm(Cohort_respinhib ~ Cohort_age, data = data4)
   			resid.res2<- (resid(resid2))
			resid.res2<- as.data.frame(resid.res2)
   			data4 <- cbind(data4, resid.res2)
   			data4$Cohort_ri_res <- data4$resid.res2
   
   	
data$nogoacc_z <-log10(data$Semantic_nogoacc)
meannogoacc  <- mean(na.omit(data$nogoacc_z))
sdnogoacc  <- sd(na.omit(data$nogoacc_z))
data$nogoacc_z <- ((data$nogoacc_z)-(meannogoacc))/(sdnogoacc)
data <- ungroup(data)

#exclusions
SLdata1 <- subset(SLdata1, Exclude!= "x")
data1 <- subset(data1, Exclude!= "x")
VWdata1 <- subset(VWdata, Exclude!= "x")

## descriptive data 

## language ability generally
WJOCall <- summary(data$Semantic_WJOC)
WJOCall
WJOCallsd <- sd(data$Semantic_WJOC)
WJOCallsd
WJSRall <- summary(data$Semantic_WJSR)
WJSRall
WJSRallsd <- sd(na.omit(data$Semantic_WJS))
WJSRallsd
below_85_count <- sum(data$Semantic_WJOC < 85)
below_85_count
below_85_count <- sum(data$Semantic_WJSR < 85)
below_85_count
# by age group 
WJOC <- tapply(data$Semantic_WJOC, data$Semantic_agegroup, summary)
WJOC
library(dplyr)
WJOC_summary <- data %>%
  group_by(Semantic_agegroup) %>%
  summarize(count_below_85 = sum(Semantic_WJOC < 85))
WJOC_summary
WJSR_summary <- data %>%
  group_by(Semantic_agegroup) %>%
  summarize(count_below_85 = sum(Semantic_WJSR < 85))
WJSR_summary
WJOCsd <- tapply(data$Semantic_WJOC, data$Semantic_agegroup, sd)
WJOCsd
WJSR <- tapply(data$Semantic_WJSR, data$Semantic_agegroup, summary)
WJSR
WJSRsd <- tapply(data$Semantic_WJSR, data$Semantic_agegroup, sd)
WJSRsd

#anova comparisons of language ability. 
language <- lm(Semantic_WJOC ~ Semantic_age_c, data = data)
summary(language)

language1 <- lm(Semantic_WJSR ~ Semantic_age_c, data = data)
summary(language1)

# average target fits 
R2statscohort <- summary(data$Cohort_R)
R2statscohort
R2cohortsd <- sd(data$Cohort_R)
R2cohortsd

R2statssem <- summary(data$Semantic_R)
R2statssem
R2semsd <- sd(data$Semantic_R)
R2semsd

#average competitor fits 
R2statscohortcomp <- summary(data$Cohort_CR)
R2statscohortcomp
R2cohortcompsd <- sd(data$Cohort_CR)
R2cohortcompsd

R2statssemcomp <- summary(data$Semantic_CR)
R2statssemcomp
R2semcompsd <- sd(data$Semantic_CR)
R2semcompsd

#average unrelated fits 
R2statscohortunr <- summary(data$Cohort_UR)
R2statscohortunr
R2cohortunrsd <- sd(data$Cohort_UR)
R2cohortunrsd

R2statssemunr <- summary(data$Semantic_UR)
R2statssemunr
R2semunrsd <- sd(data$Semantic_UR)
R2semunrsd


#mean and sd for accuracy
summary(data$Cohort_acc) #mean accuracy for cohort trials 
summary(data$Semantic_acc) # mean accuracy for semantic trials 
sd(data$Cohort_acc) # sd for cohort
sd(data$Semantic_acc) # sd for semantic

#accuracy and rt by trial type and age
datax <- read.table("data\\accrt.csv", header = TRUE, sep = ",")
datax <- read.table("accrt.csv", header = TRUE, sep = ",")
datax$rt_c <- datax$rt - mean(datax$rt)
datax$acc_c <- datax$acc- mean(datax$acc)

acc <-anova_test(data = datax, dv= acc, wid = subjectID, between = agegroup, within = trialtype) 
acc
rt <-anova_test(data = datax, dv= rt, wid = subjectID, between = agegroup, within = trialtype) 
rt
aggregate(acc ~ agegroup, datax, mean)
aggregate(rt ~ agegroup, datax, mean)
aggregate(rt ~ agegroup, datax, sd)
aggregate(rt ~ trialtype, datax, mean)
aggregate(rt ~ trialtype, datax, sd)

#mean and sd for RT 
summary(data$Cohort_rt) #mean rt for cohort trials 
summary(data$Semantic_rt) # mean rt for semantic trials 
sd(data$Cohort_rt) # sd for cohort
sd(data$Semantic_rt)# sd for semantic
data$Cohort_rt <- log10(data$Cohort_rt)
data$Semantic_rt <- log10(data$Semantic_rt)


############################### EVALUATING FIXED-EFFECTS QUESTIONS ##################################
## Age by cohort 

####TIMING#####

####language main effects 
#cohort
datac <- data[ !is.na(data$WJSR_c), ]   ###### NEED TO IMPUTE WOODCOCK SR
data2<- subset(data, !is.na(VVWP_targettiming),) # no NAs on VVWP
data3<- subset(data, !is.na(WJOC_c),) # no NAs on WJOC 
data3<- subset(data3, !is.na(WJSR_c),) # no NAS on WJSR
model1h <- lm(Cohort_targettiming ~ 1+ poly(Cohort_age_c, 2), data = data, na.action = na.omit)# age_c 
summary(model1h)
model1h <- lm(Cohort_targettiming ~ 1+ Cohort_age_c, data = data, na.action = na.omit)# age_c 
summary(model1h)
model1i <- lm(Cohort_targettiming ~ 1+ Cohort_age_c + WJOC_c + WJSR_c, data = data, na.action = na.omit) # main effect of language_c 
summary(model1i)
anova(model1h, model1i)
model1j <- lm(Cohort_targettiming ~ 1+ Cohort_age_c*WJOC_c +  Cohort_age_c*WJSR_c, data = data, na.action = na.omit)#interactions with language_c
summary(model1j)
anova(model1h, model1i, model1j)#interaction model is moderately significant over age_c
correlation <- cor(data$VVWP_eyemove, data$Cohort_targettiming, use = "complete.obs")
correlation
## CC main effects 
model1h <- lm(Cohort_targettiming ~ 1+ Cohort_age_c, data = data4,na.action = na.omit)# age_c 
summary(model1h)
model1i <- lm(Cohort_targettiming ~ 1+ Cohort_age_c + Cohort_gs_res + Cohort_ri_res, data = data4)
summary(model1i)
anova(model1h, model1i)
model1j <- lm(Cohort_targettiming ~ 1+ Cohort_age_c*Cohort_gs_res + Cohort_age_c*Cohort_ri_res, data = data4)#interactions with language_c
summary(model1j)
anova(model1h, model1i, model1j)#interaction model is not significant over age
ccttr2 <- (summary(model1i)$r.squared) - (summary(model1h)$r.squared)
ccttr2
ccttr3 <- (summary(model1j)$r.squared) - (summary(model1i)$r.squared)
ccttr3

datax <- subset(data, Cohort_UR>.8)
dataxx <- subset(data4, Cohort_UR>.8)

####RESOLUTION 
####### language main effects 
model2u <- lm(Cohort_resolution ~ 1+ poly(Cohort_age_c, 2), data = datax, na.action = na.omit) # age_c does not predict resolution 
summary(model2u) # no effect of age 
model2u <- lm(Cohort_resolution ~ 1+ Cohort_age_c, data = datax, na.action = na.omit) # age_c does not predict resolution 
summary(model2u) # no effect of age 
model2u1 <- lm(Cohort_resolution ~ 1+ Cohort_age_c + WJOC_c + WJSR_c, data = datax, na.action = na.omit)
summary(model2u1)
anova(model2u, model2u1)
model2v <- lm(Cohort_resolution ~ 1+ Cohort_age_c*WJOC_c +  Cohort_age_c*WJSR_c, data = datax, na.action = na.omit)
summary(model2v)
anova(model2u, model2u1, model2v)
#follow up on moderate interaction
#young 
datax <- subset(data, Cohort_agegroup == "young")
cor.test(datax$Cohort_resolution, datax$WJOC_c)
#middle 
datax <- subset(data, Cohort_agegroup == "middle")
cor.test(datax$Cohort_resolution, datax$WJOC_c)
#old 
datax <- subset(data, Cohort_agegroup == "teen")
cor.test(datax$Cohort_resolution, datax$WJOC_c)

#CC main effects 
# use jtools to compare effect sizes of CC on target timing and resolution. This will determine what the magnitude is of CC. 
model1h <- lm(Cohort_resolution ~ 1+ Cohort_age_c, data = dataxx)# age_c 
summary(model1h)
model1i <- lm(Cohort_resolution ~ 1+ Cohort_age_c + Cohort_gs_res + Cohort_ri_res, data = dataxx)
summary(model1i)
anova(model1h, model1i)
model1j <- lm(Cohort_resolution ~ 1+ Cohort_age_c*Cohort_gs_res + Cohort_age_c*Cohort_ri_res, data = dataxx)#interactions with language_c
summary(model1j)
anova(model1h, model1i, model1j)#interaction model is not significant over age
ccrr <- (summary(model1i)$r.squared) - (summary(model1h)$r.squared)
ccrr
ccrr2<- (summary(model1j)$r.squared) - (summary(model1i)$r.squared)
ccrr2

# compare effect sizes 
library(jtools)

# Create a linear regression model for each standardized dependent variable
lm_resolution <- lm(Cohort_resolution_std ~ Cohort_age_c + Cohort_gs_res + Cohort_ri_res, data = data4)
lm_targettiming <- lm(Cohort_targettiming_std ~ Cohort_age_c + Cohort_gs_res + Cohort_ri_res, data = data4)

# Use 'summ' to summarize the models with transformed response
summ(lm_resolution, scale = TRUE, transform.response = TRUE)
summ(lm_targettiming, scale = TRUE, transform.response = TRUE)


####PEAK HEIGHT 

###### language main effects
datax$Cohort_peak <- (datax$Cohort_CH)-(datax$Cohort_UH)
dataxx$Cohort_peak <- (dataxx$Cohort_CH) - (dataxx$Cohort_UH)
model1w <- lm((Cohort_CH-Cohort_UH) ~ 1+ poly(Cohort_age_c, 2), data = datax, na.action = na.omit)
summary(model1w) 
model1w <- lm((Cohort_CH-Cohort_UH) ~ 1+ Cohort_age_c, data = datax, na.action = na.omit)
summary(model1w)
model1w1 <- lm((Cohort_CH-Cohort_UH)~ 1+ Cohort_age_c + WJOC_c + WJSR_c, data = datax, na.action = na.omit)
summary(model1w1) 
anova(model1w, model1w1)
model1x <- lm((Cohort_CH-Cohort_UH)~ 1+ Cohort_age_c*WJOC_c + Cohort_age_c*WJSR_c, data = datax, na.action = na.omit)
summary(model1x) 
anova(model1w, model1w1, model1x)

## follow up on interactions
#young 
datax1 <- subset(datax, Cohort_agegroup == "young")
cor.test(datax1$Cohort_peak, datax1$WJOC_c)
#middle 
datax2 <- subset(datax, Cohort_agegroup == "middle")
cor.test(datax2$Cohort_peak, datax2$WJOC_c)
#old 
datax3 <- subset(datax, Cohort_agegroup == "teen")
cor.test(datax3$Cohort_peak, datax3$WJOC_c)

#young WJSR
cor.test(datax1$Cohort_peak, datax1$WJSR_c)
#middle 
cor.test(datax2$Cohort_peak, datax2$WJSR_c)
#old 
cor.test(datax3$Cohort_peak, datax3$WJSR_c)

#median split look at age effects 
summary(dataxx$Cohort_WJOC)
lowOC <- subset(dataxx, Cohort_WJOC<103.41)
highOC <- subset(dataxx, Cohort_WJOC>103.41)
ageOClow <- lm((Cohort_CH-Cohort_UH) ~ 1 + Cohort_age_c, data = lowOC)
summary(ageOClow)
ageOChigh <- lm((Cohort_CH-Cohort_UH) ~ 1 + Cohort_age_c, data = highOC)
summary(ageOChigh)

summary(dataxx$Cohort_WJSR)
lowSR <- subset(dataxx, Cohort_WJSR<98.49)
highSR <- subset(dataxx, Cohort_WJSR>98.49)
ageSRlow <- lm((Cohort_CH-Cohort_UH) ~ 1 + Cohort_age_c, data = lowSR)
summary(ageSRlow)
ageSRhigh <- lm((Cohort_CH-Cohort_UH) ~ 1 + Cohort_age_c, data = highSR)
summary(ageSRhigh)

#CC main effects 
model1h <- lm(Cohort_peak~Cohort_age_c, data = dataxx)
summary(model1h)
model1i <- lm(Cohort_peak ~ 1+ Cohort_age_c + Cohort_gs_res + Cohort_ri_res, data = dataxx)
summary(model1i)
anova(model1h, model1i)
model1j <- lm(Cohort_peak~ 1+ Cohort_age_c*Cohort_gs_res + Cohort_age_c*Cohort_ri_res, data = dataxx)
summary(model1j)
anova(model1h, model1i, model1j)#interaction model is not significant over age
ccp<- (summary(model1i)$r.squared) - (summary(model1h)$r.squared)
ccp
ccrr2<- (summary(model1j)$r.squared) - (summary(model1i)$r.squared)
ccrr2

##ONSET 
###### language_c main effects
model2r <- lm(Cohort_MinTime ~ 1+ poly(Cohort_age_c, 2), data = datax, na.action = na.omit) # age_c does not predict semantic extent 
summary(model2r) # moderate effect
model2r <- lm(Cohort_MinTime ~ 1+ Cohort_age_c, data = datax, na.action = na.omit) # age_c does not predict semantic extent 
summary(model2r) 
model2r1 <- lm(Cohort_MinTime ~ 1+ Cohort_age_c + WJOC_c + WJSR_c, data = datax, na.action = na.omit)
summary(model2r1)
anova(model2r, model2r1)
model2s <- lm(Cohort_MinTime ~ 1+ Cohort_age_c*WJOC_c + Cohort_age_c*WJSR_c, data = datax, na.action = na.omit)
summary(model2s)
anova(model2r, model2r1, model2s)

#CC main effects 
model1h <- lm(Cohort_MinTime~ 1 + Cohort_age_c, data = dataxx)
summary(model1h)
model1i <- lm(Cohort_MinTime ~ 1 + Cohort_age_c +Cohort_gs_res + Cohort_ri_res, data = dataxx)
summary(model1i)
anova(model1h, model1i)
model1j <- lm(Cohort_MinTime ~ 1 + Cohort_age_c*Cohort_gs_res + Cohort_age_c*Cohort_ri_res,data=dataxx)
summary(model1j)
anova(model1h, model1i, model1j)#interaction model is not significant over age
cco<-(summary(model1i)$r.squared) - (summary(model1h)$r.squared)
cco
ccrr2<- (summary(model1j)$r.squared) - (summary(model1i)$r.squared)
ccrr2

##OFFSET 
###### language_c main effects
model2r <- lm(Cohort_MaxTime ~ 1+ poly(Cohort_age_c, 2), data = datax, na.action = na.omit) # age_c does not predict semantic extent 
summary(model2r)
model2r <- lm(Cohort_MaxTime ~ 1+ Cohort_age_c, data = datax, na.action = na.omit) # age_c does not predict semantic extent 
summary(model2r)
model2r1 <- lm(Cohort_MaxTime ~ 1+ Cohort_age_c + WJOC_c + WJSR_c, data = datax, na.action = na.omit)
summary(model2r1)
anova(model2r, model2r1)
model2s <- lm(Cohort_MaxTime ~ 1+ Cohort_age_c*WJOC_c + Cohort_age_c*WJSR_c, data = datax, na.action = na.omit)
summary(model2s)
anova(model2r, model2r1, model2s)

#young
cor.test(datax1$Cohort_MaxTime, datax1$WJSR_c)
#middle 
cor.test(datax2$Cohort_MaxTime, datax2$WJSR_c)
#old 
cor.test(datax3$Cohort_MaxTime, datax3$WJSR_c)


#median split look at age effects 

ageOClow <- lm(Cohort_MaxTime ~ 1 + Cohort_age_c, data = lowOC)
summary(ageOClow)
ageOChigh <- lm(Cohort_Maxtime ~ 1 + Cohort_age_c, data = highOC)
summary(ageOChigh)

ageSRlow <- lm(Cohort_MaxTime ~ 1 + Cohort_age_c, data = lowSR)
summary(ageSRlow)
ageSRhigh <- lm(Cohort_MaxTime ~ 1 + Cohort_age_c, data = highSR)
summary(ageSRhigh)


#CC main effects 
model1h <- lm(Cohort_MaxTime~ 1 + Cohort_age_c, data = dataxx)
summary(model1h)
model1i <- lm(Cohort_MaxTime ~ 1 + Cohort_age_c +Cohort_gs_res + Cohort_ri_res, data = dataxx)
summary(model1i)
anova(model1h, model1i)
model1j <- lm(Cohort_MaxTime~1+ Cohort_age_c*Cohort_gs_res + Cohort_age_c*Cohort_ri_res, data=dataxx)
summary(model1j)
anova(model1h, model1i, model1j)#interaction model is not significant over age
cco<-(summary(model1i)$r.squared) - (summary(model1h)$r.squared)
cco
ccrr2<- (summary(model1j)$r.squared) - (summary(model1i)$r.squared)
ccrr2


datay <- subset(data, (!is.na(VVWP_CB1)))
datay <- subset(datay, (!is.na(Semantic_ssrt_c)))
datax <- subset(data, (!is.na(VVWP_CB1)))

## summary discussion about development of these parameters. Let's then account for VVWP for the ones that show age effects only.
#### VVWP 
#simple age effects 
model1h <- lm(VVWP_targettiming ~ 1+ Cohort_age_c, data = datax, na.action = na.omit)# age_c 
summary(model1h)
model1h <- lm(VVWP_resolution ~ 1+ Cohort_age_c, data = datax, na.action = na.omit)# age_c 
summary(model1h)
model1h <- lm((VVWP_CH-VVWP_UH) ~ 1+ Cohort_age_c, data = datax, na.action = na.omit)# age_c 
summary(model1h)
model1h <- lm(VVWP_MinTime ~ 1+ Cohort_age_c, data = datax, na.action = na.omit)# age_c 
summary(model1h)
model1h <- lm(VVWP_MaxTime ~ 1+ Cohort_age_c, data = datax, na.action = na.omit)# age_c 
summary(model1h)

#simple SWR effects 
model1h <- lm(VVWP_targettiming ~ 1+ Cohort_targettiming, data = datax, na.action = na.omit)# age_c 
summary(model1h)
model1h <- lm((VVWP_CH-VVWP_UH) ~ 1+ (Cohort_CH-Cohort_UH), data = datax, na.action = na.omit)# age_c 
summary(model1h)
model1h <- lm(VVWP_MinTime ~ 1+ Cohort_MinTime, data = datax, na.action = na.omit)# age_c 
summary(model1h)
model1h <- lm(VVWP_MaxTime ~ 1+ Cohort_MaxTime, data = datax, na.action = na.omit)# age_c 
summary(model1h)
#Cohort target timing 
  #C: Fixation Duration 
      modelzaz <- lm(Cohort_targettiming ~ 1 + VVWP_eyemove_c, data = datay)
      summary(modelzaz) 
  #D: Fixation Duration + Visual
      modelzz2 <- lm(Cohort_targettiming ~ 1+ VVWP_eyemove_c + (VVWP_targettiming+Semantic_ssrt_c),data=datay)
      summary(modelzz2) # significant effect, older children --> earlier target timing 
      mainvis <- (summary(modelzz2)$r.squared) - (summary(modelzaz)$r.squared)
        mainvis
  #E: Fixation Duration + Age  
      modelzbz <- lm(Cohort_targettiming ~ 1+ VVWP_eyemove_c + Cohort_age_c, data = datay) # both age_c and VVWP predicting TT 
      summary(modelzbz)
      mainage <- (summary(modelzbz)$r.squared) - (summary(modelzaz)$r.squared)
      mainage
  #F: Fixation Duration + Visual + Age 
      modelzbzz<- lm(Cohort_targettiming ~ 1+ VVWP_eyemove_c + Cohort_age_c + VVWP_targettiming+Semantic_ssrt_c,data=datay) 
      summary(modelzbzz)
    modelzbzz2<-lm(Cohort_targettiming ~ 1 + VVWP_targettiming+Semantic_ssrt_c+Cohort_age_c,data=datay) 
      summary(modelzbzz2)
        anova(modelzz2,modelzbzz)# does age predict
	  age<-(summary(modelzbzz)$r.squared) - (summary(modelzz2)$r.squared)
	  age
        anova(modelzbz,modelzbzz)# does VVWP predict
	  vis<-(summary(modelzbzz)$r.squared) - (summary(modelzbz)$r.squared)
	  vis
	  FD <- (summary(modelzaz)$r.squared)
	  shared<-(summary(modelzbzz)$r.squared)-(age+vis+FD)
	  shared
	  anova(modelzz2,modelzbzz)# does age predict

#correlations betweeen visual and cognitive 
	  data4$Vpeak <- (data4$VVWP_CH-data4$VVWP_UH)
	  data5 <- subset(data4,(!is.na(Vpeak)))
	  cor.test(data4$VVWP_MaxTime, data4$Cohort_gs)
	  cor.test(data4$VVWP_MinTime, data4$Cohort_gs)
	  cor.test(data4$VVWP_targettiming, data4$Cohort_gs)
	  cor.test(data4$Vpeak, data4$Cohort_gs)
	  cor.test(data4$VVWP_MaxTime, data4$Cohort_ri)
	  cor.test(data4$VVWP_MinTime, data4$Cohort_ri)
	  cor.test(data4$VVWP_targettiming, data4$Cohort_ri)
	  cor.test(data4$Vpeak, data4$Cohort_ri)
	  cor.test(data4$VVWP_MaxTime, data4$Semantic_ssrt)
	  cor.test(data4$VVWP_MinTime, data4$Semantic_ssrt)
	  cor.test(data4$VVWP_targettiming, data4$Semantic_ssrt)
	  cor.test(data4$Vpeak, data4$Semantic_ssrt)
	  cor.test(data4$VVWP_MaxTime, data4$VVWP_eyemove)
	  cor.test(data4$VVWP_MinTime, data4$VVWP_eyemove)
	  cor.test(data4$VVWP_targettiming, data4$VVWP_eyemove)
	  cor.test(data4$Vpeak, data4$VVWP_eyemove)

###### Cohort Peak
data3 <- subset(datay, Cohort_UR>.8)
data3$Cpeak <- (data3$Cohort_CH-data3$Cohort_UH)
data3$Vpeak <- (data3$VVWP_CH-data3$VVWP_UH)
  #C: Fixation Duration 
      modelzaz <- lm(Cpeak ~ 1 + VVWP_eyemove_c, data = datay)
      summary(modelzaz) 
  #D: Fixation Duration + Visual
      modelzz2 <- lm(Cpeak ~ 1+ VVWP_eyemove_c + Vpeak + Semantic_ssrt_c, data = datay)
      summary(modelzz2) # significant effect, older children --> earlier target timing 
  #E: Fixation Duration + Age  
      modelzbz <- lm(Cpeak ~ 1+ VVWP_eyemove_c + Cohort_age_c, data = datay) # both age_c and VVWP predicting TT 
      summary(modelzbz)
  #F: Fixation Duration + Visual + Age 
      modelzbzz<- lm(Cpeak ~ 1+ VVWP_eyemove_c + Vpeak + Cohort_age_c, data = datay) # both age_c and VVWP predicting TT 
      summary(modelzbzz)
        anova(modelzz2,modelzbzz)# does age predict
	  age<-(summary(modelzbzz)$r.squared) - (summary(modelzz2)$r.squared)
	  age
        anova(modelzbz,modelzbzz)# does VVWP predict
	  vis<-(summary(modelzbzz)$r.squared) - (summary(modelzbz)$r.squared)
	  vis
	  shared<-(summary(modelzbzz)$r.squared)-(age+vis)
	  shared

###### Cohort Onset Time

  #C: Fixation Duration 
      modelzaz <- lm(Cohort_MinTime~ 1 + VVWP_eyemove_c, data = datay)
      summary(modelzaz) 
  #D: Fixation Duration + Visual
      modelzz2 <- lm(Cohort_MinTime ~ 1+ VVWP_eyemove_c + (VVWP_MinTime+Semantic_ssrt_c),data=datay)
      summary(modelzz2) # significant effect, older children --> earlier target timing 
      mainvis <- (summary(modelzz2)$r.squared) - (summary(modelzaz)$r.squared)
      mainvis
  #E: Fixation Duration + Age  
      modelzbz <- lm(Cohort_MinTime ~ 1+ VVWP_eyemove_c + Cohort_age_c, data = datay) # both age_c and VVWP predicting TT 
      summary(modelzbz)
      mainage <- (summary(modelzbz)$r.squared) - (summary(modelzaz)$r.squared)
      mainage
  #F: Fixation Duration + Visual + Age 
      modelzbzz<- lm(Cohort_MinTime ~ 1+ VVWP_eyemove_c + VVWP_MinTime+Semantic_ssrt_c+Cohort_age_c,data=datay) 
      summary(modelzbzz)
        anova(modelzz2,modelzbzz)# does age predict
	  age<-(summary(modelzbzz)$r.squared) - (summary(modelzz2)$r.squared)
	  age
        anova(modelzbz,modelzbzz)# does VVWP predict
	  vis<-(summary(modelzbzz)$r.squared) - (summary(modelzbz)$r.squared)
	  vis
	  FD <- (summary(modelzaz)$r.squared)
	  shared<-(summary(modelzbzz)$r.squared)-(age+vis+FD)
	  shared
###### Cohort Offset Time

  #C: Fixation Duration 
      modelzaz <- lm(Cohort_MaxTime~ 1 + VVWP_eyemove_c, data = datay)
      summary(modelzaz) 
  #D: Fixation Duration + Visual
      modelzz2 <- lm(Cohort_MaxTime ~ 1+ VVWP_eyemove_c + (VVWP_MaxTime+Semantic_ssrt_c),data=datay)
      summary(modelzz2) # significant effect, older children --> earlier target timing 
      mainvis <- (summary(modelzz2)$r.squared) - (summary(modelzaz)$r.squared)
      mainvis
  #E: Fixation Duration + Age  
      modelzbz <- lm(Cohort_MaxTime ~ 1+ VVWP_eyemove_c + Cohort_age_c, data = datay) # both age_c and VVWP predicting TT 
      summary(modelzbz)
      mainage <- (summary(modelzbz)$r.squared) - (summary(modelzaz)$r.squared)
      mainage
  #F: Fixation Duration + Visual + Age 
      modelzbzz<- lm(Cohort_MaxTime ~ 1+ VVWP_eyemove_c + VVWP_MaxTime+Semantic_ssrt_c+Cohort_age_c,data=datay) 
      summary(modelzbzz)
        anova(modelzz2,modelzbzz)# does age predict
	  age<-(summary(modelzbzz)$r.squared) - (summary(modelzz2)$r.squared)
	  age
        anova(modelzbz,modelzbzz)# does VVWP predict
	  vis<-(summary(modelzbzz)$r.squared) - (summary(modelzbz)$r.squared)
	  vis
	  FD <- (summary(modelzaz)$r.squared)
	  shared<-(summary(modelzbzz)$r.squared)-(age+vis+FD)
	  shared

## 3.b. VVWP factor structure 

datay <- subset(dataxx, (!is.na(VVWP_CB1)))
datay <- subset(datay, (!is.na(Semantic_ssrt_c)))
#data2 <- subset(datay, Cohort_UR>.8)
#data2$Cpeak <- (data2$Cohort_CH-data2$Cohort_UH)
#data2$Vpeak <- (data2$VVWP_CH-data2$VVWP_UH)

## VVWP target timing 

## use the yhat for commonality 

library(yhat)

# Run regression and save it 
lm.out <- lm(VVWP_targettiming ~ 1 + Semantic_ssrt_c + VVWP_eyemove_c+ Cohort_gs_res + Cohort_ri_res + Cohort_age_c, data = datay)
# Perform commonality analysis
commonality_model <- regr(lm.out)
# Print the results
print(commonality_model)

# Load commonality analysis data
commonality_data <- commonality_model$Commonality_Data$CCTotalbyVar
commonality_data <- as.data.frame(commonality_data)
unique_effect_sizes <- data.frame(Variable = row.names(commonality_data),
                                  Unique = commonality_data$Unique)

# Get degrees of freedom from the LM_Output of commonality_model
df<- as.numeric(commonality_model$LM_Output$df[2]) # Adjust to the actual degrees of freedom value

# Calculate p-values using the t-distribution
p_values <- rep(NA, nrow(unique_effect_sizes))
for (i in seq_len(nrow(unique_effect_sizes))) {
  r_squared <- unique_effect_sizes$Unique[i]
  r <- sqrt(r_squared)  # Convert squared value back to r
  # Calculate the t-statistic
  t_statistic <- r * sqrt(df / (1 - r^2))
  p_positive <- pt(t_statistic, df = df, lower.tail = FALSE)
  p_negative <- pt(-t_statistic, df = df, lower.tail = TRUE)
  p_values[i] <- 2 * min(p_positive, p_negative)
}

# Add p-values to the dataframe
unique_effect_sizes$p_value <- p_values
# Print the result
print(unique_effect_sizes)

# Shared variance with FD 
# Load commonality analysis data
commonality_data <- commonality_model$Commonality_Data$CC
commonality_df <- as.data.frame(commonality_data)
# Add a column for row names
commonality_df$RowNames <- row.names(commonality_df)
# Words to search for in the column titles
target_words <- c("Common", "VVWP_eyemove_c")
# Filter rows with target words in the row names and exclude "Total"
filtered_rows <- commonality_df[
  grepl(paste(target_words, collapse = "|"), commonality_df$RowNames) & !grepl("Total", commonality_df$RowNames),
]
sum_of_coefficients <- sum(filtered_rows$Coefficient)
# Print the sum and the filtered rows
print(sum_of_coefficients)

# Shared variance Cognitive 
# Exclude rows containing "VVWP_eyemove_c" and "Total" in the row names
filtered_rows <- commonality_df[
  !grepl("VVWP_eyemove_c", commonality_df$RowNames) & !grepl("Total", commonality_df$RowNames),
]
sum_of_coefficients <- sum(filtered_rows$Coefficient)
# Print the sum and the filtered rows
print(sum_of_coefficients)



#ONSET 
# Run regression and save it 
lm.out <- lm(VVWP_MinTime ~ 1 + Semantic_ssrt_c + VVWP_eyemove_c+ Cohort_gs_res + Cohort_ri_res + Cohort_age_c, data = datay)
# Perform commonality analysis
commonality_model <- regr(lm.out)
# Print the results
print(commonality_model)

# Load commonality analysis data
commonality_data <- commonality_model$Commonality_Data$CCTotalbyVar
commonality_data <- as.data.frame(commonality_data)
unique_effect_sizes <- data.frame(Variable = row.names(commonality_data),
                                  Unique = commonality_data$Unique)

# Get degrees of freedom from the LM_Output of commonality_model
df<- as.numeric(commonality_model$LM_Output$df[2]) # Adjust to the actual degrees of freedom value

# Calculate p-values using the t-distribution
p_values <- rep(NA, nrow(unique_effect_sizes))
for (i in seq_len(nrow(unique_effect_sizes))) {
  r_squared <- unique_effect_sizes$Unique[i]
  r <- sqrt(r_squared)  # Convert squared value back to r
  # Calculate the t-statistic
  t_statistic <- r * sqrt(df / (1 - r^2))
  p_positive <- pt(t_statistic, df = df, lower.tail = FALSE)
  p_negative <- pt(-t_statistic, df = df, lower.tail = TRUE)
  p_values[i] <- 2 * min(p_positive, p_negative)
}

# Add p-values to the dataframe
unique_effect_sizes$p_value <- p_values
# Print the result
print(unique_effect_sizes)

# Shared variance with FD 
# Load commonality analysis data
commonality_data <- commonality_model$Commonality_Data$CC
commonality_df <- as.data.frame(commonality_data)
# Add a column for row names
commonality_df$RowNames <- row.names(commonality_df)
# Words to search for in the column titles
target_words <- c("Common", "VVWP_eyemove_c")
# Filter rows with target words in the row names and exclude "Total"
filtered_rows <- commonality_df[
  grepl(paste(target_words, collapse = "|"), commonality_df$RowNames) & !grepl("Total", commonality_df$RowNames),
]
sum_of_coefficients <- sum(filtered_rows$Coefficient)
# Print the sum and the filtered rows
print(sum_of_coefficients)

# Shared variance Cognitive 
# Exclude rows containing "VVWP_eyemove_c" and "Total" in the row names
filtered_rows <- commonality_df[
  !grepl("VVWP_eyemove_c", commonality_df$RowNames) & !grepl("Total", commonality_df$RowNames),
]
sum_of_coefficients <- sum(filtered_rows$Coefficient)
# Print the sum and the filtered rows
print(sum_of_coefficients)


#OFFSET 
# Run regression and save it 
lm.out <- lm(VVWP_MaxTime ~ 1 + Semantic_ssrt_c + VVWP_eyemove_c+ Cohort_gs_res + Cohort_ri_res + Cohort_age_c, data = datay)
# Perform commonality analysis
commonality_model <- regr(lm.out)
# Print the results
print(commonality_model)

# Load commonality analysis data
commonality_data <- commonality_model$Commonality_Data$CCTotalbyVar
commonality_data <- as.data.frame(commonality_data)
unique_effect_sizes <- data.frame(Variable = row.names(commonality_data),
                                  Unique = commonality_data$Unique)

# Get degrees of freedom from the LM_Output of commonality_model
df<- as.numeric(commonality_model$LM_Output$df[2]) # Adjust to the actual degrees of freedom value

# Calculate p-values using the t-distribution
p_values <- rep(NA, nrow(unique_effect_sizes))
for (i in seq_len(nrow(unique_effect_sizes))) {
  r_squared <- unique_effect_sizes$Unique[i]
  r <- sqrt(r_squared)  # Convert squared value back to r
  # Calculate the t-statistic
  t_statistic <- r * sqrt(df / (1 - r^2))
  p_positive <- pt(t_statistic, df = df, lower.tail = FALSE)
  p_negative <- pt(-t_statistic, df = df, lower.tail = TRUE)
  p_values[i] <- 2 * min(p_positive, p_negative)
}

# Add p-values to the dataframe
unique_effect_sizes$p_value <- p_values
# Print the result
print(unique_effect_sizes)

# Shared variance with FD 
# Load commonality analysis data
commonality_data <- commonality_model$Commonality_Data$CC
commonality_df <- as.data.frame(commonality_data)
# Add a column for row names
commonality_df$RowNames <- row.names(commonality_df)
# Words to search for in the column titles
target_words <- c("Common", "VVWP_eyemove_c")
# Filter rows with target words in the row names and exclude "Total"
filtered_rows <- commonality_df[
  grepl(paste(target_words, collapse = "|"), commonality_df$RowNames) & !grepl("Total", commonality_df$RowNames),
]
sum_of_coefficients <- sum(filtered_rows$Coefficient)
# Print the sum and the filtered rows
print(sum_of_coefficients)

# Shared variance Cognitive 
# Exclude rows containing "VVWP_eyemove_c" and "Total" in the row names
filtered_rows <- commonality_df[
  !grepl("VVWP_eyemove_c", commonality_df$RowNames) & !grepl("Total", commonality_df$RowNames),
]
sum_of_coefficients <- sum(filtered_rows$Coefficient)
# Print the sum and the filtered rows
print(sum_of_coefficients)

#A: all variables
modelzbzz <- lm(VVWP_targettiming ~ 1 + Semantic_ssrt_c + VVWP_eyemove_c+ Cohort_gs_res + Cohort_ri_res + Cohort_age_c, data = datay)
summary(modelzbzz) 
#B: Processing speed 
modelzz <- lm(VVWP_targettiming ~ 1 + Semantic_ssrt_c, data = datay)
summary(modelzz) 
modelzz <- lm(VVWP_targettiming ~ 1 + VVWP_eyemove_c+ Cohort_gs_res + Cohort_ri_res + Cohort_age_c, data = datay)
summary(modelzz)
anova(modelzbzz, modelzz)
cco1<-(summary(modelzbzz)$r.squared) - (summary(modelzz)$r.squared)
cco1

#C: Fixation duration
modelzz <- lm(VVWP_targettiming ~ 1 + VVWP_eyemove_c, data = datay)
summary(modelzz) 
modelzz <- lm(VVWP_targettiming ~ 1 + Semantic_ssrt_c+ Cohort_gs_res + Cohort_ri_res + Cohort_age_c, data = datay)
summary(modelzz)
anova(modelzbzz, modelzz)
cco2<- (summary(modelzz)$r.squared)
#D: Conflict Suppression
modelzz <- lm(VVWP_targettiming ~ 1 + Cohort_gs_res, data = datay)
summary(modelzz)
modelzz <- lm(VVWP_targettiming ~ 1 + Semantic_ssrt_c+ VVWP_eyemove_c + Cohort_ri_res + Cohort_age_c, data = datay)
summary(modelzz)
anova(modelzbzz, modelzz)
cco3<-(summary(modelzbzz)$r.squared) - (summary(modelzz)$r.squared)
cco3
#E: Conflict Inhibition
modelzz <- lm(VVWP_targettiming ~ 1 + Cohort_ri_res, data = datay)
summary(modelzz)
modelzz <- lm(VVWP_targettiming ~ 1 + Semantic_ssrt_c+ VVWP_eyemove_c + Cohort_gs_res + Cohort_age_c, data = datay)
summary(modelzz)
anova(modelzbzz, modelzz)
cco4<-(summary(modelzbzz)$r.squared) - (summary(modelzz)$r.squared)
cco4
#F: Age 
modelzz <- lm(VVWP_targettiming ~ 1 + Semantic_ssrt_c+ VVWP_eyemove_c + Cohort_gs_res + Cohort_ri_res, data = datay)
summary(modelzz)
anova(modelzbzz, modelzz)
cco5<-(summary(modelzbzz)$r.squared) - (summary(modelzz)$r.squared)
cco5
shared <- (summary(modelzbzz)$r.squared) - (cco1 + cco2 + cco3 + cco4 + cco5)
shared
## VVWP Onset Time 

#A: all variables
modelzbzz <- lm(VVWP_MinTime ~ 1 + Semantic_ssrt_c + VVWP_eyemove_c+ Cohort_gs_res + Cohort_ri_res + Cohort_age_c, data = datay)
summary(modelzbzz) # no main effect of fixation duration 
#B: Processing speed 
modelzz <- lm(VVWP_MinTime ~ 1 + Semantic_ssrt_c, data = datay)
summary(modelzz) 
modelzz <- lm(VVWP_MinTime ~ 1 + VVWP_eyemove_c+ Cohort_gs_res + Cohort_ri_res + Cohort_age_c, data = datay)
summary(modelzz)
anova(modelzbzz, modelzz)
cco1<-(summary(modelzbzz)$r.squared) - (summary(modelzz)$r.squared)
cco1
#C: Fixation duration
modelzz <- lm(VVWP_MinTime ~ 1 + VVWP_eyemove_c, data = datay)
summary(modelzz) 
modelzz <- lm(VVWP_MinTime ~ 1 + Semantic_ssrt_c+ Cohort_gs_res + Cohort_ri_res + Cohort_age_c, data = datay)
summary(modelzz)
anova(modelzbzz, modelzz)
cco2<-(summary(modelzz)$r.squared)
cco2
#D: Conflict Suppression 
modelzz <- lm(VVWP_MinTime ~ 1 + Cohort_gs_res, data = datay)
summary(modelzz)
modelzz <- lm(VVWP_MinTime ~ 1 + Semantic_ssrt_c+ VVWP_eyemove_c + Cohort_ri_res + Cohort_age_c, data = datay)
summary(modelzz)
anova(modelzbzz, modelzz)
cco3<-(summary(modelzbzz)$r.squared) - (summary(modelzz)$r.squared)
cco3
#E: Conflict Inhibition
modelzz <- lm(VVWP_MinTime ~ 1 + Cohort_ri_res, data = datay)
summary(modelzz)
modelzz <- lm(VVWP_MinTime ~ 1 + Semantic_ssrt_c+ VVWP_eyemove_c + Cohort_gs_res + Cohort_age_c, data = datay)
summary(modelzz)
anova(modelzbzz, modelzz)
cco4<-(summary(modelzbzz)$r.squared) - (summary(modelzz)$r.squared)
cco4
#F: Age 
modelzz <- lm(VVWP_MinTime ~ 1 + Semantic_ssrt_c+ VVWP_eyemove_c + Cohort_gs_res + Cohort_ri_res, data = datay)
summary(modelzz)
anova(modelzbzz, modelzz)
cco5<-(summary(modelzbzz)$r.squared) - (summary(modelzz)$r.squared)
cco5
shared <- (summary(modelzbzz)$r.squared) - (cco1 + cco2 + cco3 + cco4 + cco5)
shared


## VVWP Offset Time 

#A: all variables
modelzbzz <- lm(VVWP_MaxTime ~ 1 + Semantic_ssrt_c + VVWP_eyemove_c+ Cohort_gs_res + Cohort_ri_res + Cohort_age_c, data = datay)
summary(modelzbzz) 
#B: Processing speed 
modelzz <- lm(VVWP_MaxTime ~ 1 + VVWP_eyemove_c, data = datay)
summary(modelzz) 
modelzz <- lm(VVWP_MaxTime ~ 1 + Semantic_ssrt_c, data = datay)
summary(modelzz) 
modelzz <- lm(VVWP_MaxTime ~ 1 + VVWP_eyemove_c+ Cohort_gs_res + Cohort_ri_res + Cohort_age_c, data = datay)
summary(modelzz)
anova(modelzbzz, modelzz)
cco1<-(summary(modelzbzz)$r.squared) - (summary(modelzz)$r.squared)
cco1
#C: Fixation duration
modelzz <- lm(VVWP_MaxTime ~ 1 + Semantic_ssrt_c+ Cohort_gs_res + Cohort_ri_res + Cohort_age_c, data = datay)
summary(modelzz)
anova(modelzz, modelzbzz)
cco2<-(summary(modelzz)$r.squared)
cco2
#D: Conflict Suppression
modelzz <- lm(VVWP_MaxTime ~ 1 + Cohort_gs_res, data = datay)
summary(modelzz)
modelzz <- lm(VVWP_MaxTime ~ 1 + Semantic_ssrt_c+ VVWP_eyemove_c + Cohort_ri_res + Cohort_age_c, data = datay)
summary(modelzz)
anova(modelzbzz, modelzz)
cco3<-(summary(modelzbzz)$r.squared) - (summary(modelzz)$r.squared)
cco3
#E: Conflict Inhibition
modelzz <- lm(VVWP_MaxTime ~ 1 + Cohort_ri_res, data = datay)
summary(modelzz)
modelzz <- lm(VVWP_MaxTime ~ 1 + Semantic_ssrt_c+ VVWP_eyemove_c + Cohort_gs_res + Cohort_age_c, data = datay)
summary(modelzz)
anova(modelzbzz, modelzz)
cco4<-(summary(modelzbzz)$r.squared) - (summary(modelzz)$r.squared)
cco4
#F: Age 
modelzz <- lm(VVWP_MaxTime ~ 1 + Semantic_ssrt_c+ VVWP_eyemove_c + Cohort_gs_res + Cohort_ri_res, data = datay)
summary(modelzz)
anova(modelzbzz, modelzz)
cco5<-(summary(modelzbzz)$r.squared) - (summary(modelzz)$r.squared)
cco5
shared <- (summary(modelzbzz)$r.squared) - (cco1 + cco2 + cco3 + cco4 + cco5)
shared




############SEMANTIC COMP############## 

#####PEAK HEIGHT
####### language_c main effects
  dataz <- subset(data, Semantic_UR>.8)
  model1z <- lm((Semantic_CH-Semantic_UH) ~ 1+ Semantic_age_c, data = dataz, na.action = na.omit) # age_c does not significantly predict comp peak height
  summary(model1z) 
  model1z1 <- lm((Semantic_CH-Semantic_UH)~ 1+ Semantic_age_c+ WJOC_c + WJSR_c, data = dataz, na.action = na.omit)
  summary(model1z1)
  model1z1 <- lm((Semantic_CH-Semantic_UH)~ 1+ WJOC_c + WJSR_c, data = dataz, na.action = na.omit)
  summary(model1z1)
  anova(model1z, model1z1)
  model2a <- lm((Semantic_CH-Semantic_UH)~ 1+ Semantic_age_c*WJOC_c + Semantic_age_c*WJSR_c, data = dataz, na.action = na.omit)
  summary(model2a) 
  anova(model1z, model1z1, model2a)

#CC main effects 
model1h <- lm((Semantic_CH-Semantic_UH)~ 1 + Cohort_age_c, data = dataxx)
summary(model1h)
model1i <- lm((Semantic_CH-Semantic_UH) ~ 1 + Cohort_age_c +Cohort_gs_res + Cohort_ri_res, data = dataxx)
summary(model1i)
anova(model1h, model1i)
model1j <- lm((Semantic_CH-Semantic_UH)~1+ Cohort_age_c*Cohort_gs_res + Cohort_age_c*Cohort_ri_res, data=dataxx)
summary(model1j)
anova(model1h, model1i, model1j)#interaction model is not significant over age
cco<-(summary(model1i)$r.squared) - (summary(model1h)$r.squared)
cco
ccrr2<- (summary(model1j)$r.squared) - (summary(model1i)$r.squared)
ccrr2


##ONSET 
###### language_c main effects
   model2r <- lm(Semantic_MinTime ~ 1+ Semantic_age_c, data = dataz, na.action = na.omit) # age_c does not predict semantic extent 
   summary(model2r)
   model2r1 <- lm(Semantic_MinTime ~ 1+ Semantic_age_c + WJOC_c + WJSR_c, data = dataz, na.action = na.omit)
   summary(model2r1)
   anova(model2r, model2r1)
   model2s <- lm(Semantic_MinTime ~ 1+ Semantic_age_c*WJOC_c + Semantic_age_c*WJSR_c, data = dataz, na.action = na.omit)
   summary(model2s)
   anova(model2r, model2r1, model2s)

#CC main effects 
model1h <- lm(Semantic_MinTime~ 1 + Cohort_age_c, data = dataxx)
summary(model1h)
model1i <- lm(Semantic_MinTime ~ 1 + Cohort_age_c +Cohort_gs_res + Cohort_ri_res, data = dataxx)
summary(model1i)
anova(model1h, model1i)
model1j <- lm(Semantic_MinTime~1+ Cohort_age_c*Cohort_gs_res + Cohort_age_c*Cohort_ri_res, data=dataxx)
summary(model1j)
anova(model1h, model1i, model1j)#interaction model is not significant over age
cco<-(summary(model1i)$r.squared) - (summary(model1h)$r.squared)
cco
ccrr2<- (summary(model1j)$r.squared) - (summary(model1i)$r.squared)
ccrr2

##OFFSET 
###### language_c main effects
   model2r <- lm(Semantic_MaxTime ~ 1+ Semantic_age_c, data = dataz, na.action = na.omit) # age_c does not predict semantic extent 
   summary(model2r)
   model2r1 <- lm(Semantic_MaxTime ~ 1+ Semantic_age_c + WJOC_c + WJSR_c, data = dataz, na.action = na.omit)
   summary(model2r1)
   anova(model2r, model2r1)
   model2s <- lm(Semantic_MaxTime ~ 1+ Semantic_age_c*WJOC_c + Semantic_age_c*WJSR_c, data = dataz, na.action = na.omit)
   summary(model2s)
   anova(model2r, model2r1, model2s)

#CC main effects 
model1h <- lm(Semantic_MaxTime~ 1 + Cohort_age_c, data = dataxx)
summary(model1h)
model1i <- lm(Semantic_MaxTime ~ 1 + Cohort_age_c +Cohort_gs_res + Cohort_ri_res, data = dataxx)
summary(model1i)
anova(model1h, model1i)
model1j <- lm(Semantic_MaxTime~1+ Cohort_age_c*Cohort_gs_res + Cohort_age_c*Cohort_ri_res, data=dataxx)
summary(model1j)
anova(model1h, model1i, model1j)#interaction model is not significant over age
cco<-(summary(model1i)$r.squared) - (summary(model1h)$r.squared)
cco
ccrr2<- (summary(model1j)$r.squared) - (summary(model1i)$r.squared)
ccrr2

### CONTROLLING FOR COHORT TARGET TIMING
   ## Peak
   model4m <- lm((Semantic_CH-Semantic_UH) ~ 1+ Semantic_age_c + (Cohort_targettiming), data = dataz) # both age_c and VVWP do not predict peak
   summary(model4m)
   model4k <- lm((Semantic_CH-Semantic_UH) ~ 1+ Semantic_age_c, data = dataz) # age_c does not predict peak 
   summary(model4k)
   model4l <- lm((Semantic_CH-Semantic_UH) ~ 1 + Cohort_targettiming, data = dataz) # VVWP does not predict peak
   summary(model4l)
   anova(model4m, model4l) # age_c predicts over and above VVWP 
   anova(model4m, model4k) # VVWP predicts over and above age_c 
  
uniqueap <- (summary(model4m)$r.squared) - (summary(model4l)$r.squared) # R squared change for unique processing speed variance
uniqueap
uniquecp <- (summary(model4m)$r.squared) - (summary(model4k)$r.squared) # unique visual r squared over processing speed 
uniquecp
sharedacp <- (summary(model4m)$r.squared) - (uniqueap+uniquecp)
sharedacp

   ### ONSET Timing 
   model4m <- lm((Semantic_MinTime) ~ 1+ Semantic_age_c + (Cohort_targettiming), data = dataz) # both age_c and VVWP do not predict peak
   summary(model4m)
   model4k <- lm((Semantic_MinTime) ~ 1+ Semantic_age_c, data = dataz) # age_c does not predict peak 
   summary(model4k)
   model4l <- lm((Semantic_MinTime) ~ 1 + Cohort_targettiming, data = dataz) # VVWP does not predict peak
   summary(model4l)
   anova(model4l, model4m) # age_c predicts over and above VVWP 
   anova(model4k, model4m) # VVWP predicts over and above age_c 

uniqueap <- (summary(model4m)$r.squared) - (summary(model4l)$r.squared) # R squared change for unique processing speed variance
uniqueap
uniquecp <- (summary(model4m)$r.squared) - (summary(model4k)$r.squared) # unique visual r squared over processing speed 
uniquecp
sharedacp <- (summary(model4m)$r.squared) - (uniqueap+uniquecp)
sharedacp

  ### Offset Time
  model4m <- lm((Semantic_MaxTime) ~ 1+ Semantic_age_c + (Cohort_targettiming), data = dataz) # both age_c and VVWP do not predict peak
  summary(model4m)
  model4k <- lm((Semantic_MaxTime) ~ 1+ Semantic_age_c, data = dataz) # age_c does not predict peak 
  summary(model4k)
  model4l <- lm((Semantic_MaxTime) ~ 1 + Cohort_targettiming, data = dataz) # VVWP does not predict peak
  summary(model4l)
  anova(model4l, model4m) # age_c predicts over and above VVWP 
  anova(model4k, model4m) # VVWP predicts over and above age_c 


uniqueap <- (summary(model4m)$r.squared) - (summary(model4l)$r.squared) # R squared change for unique processing speed variance
uniqueap
uniquecp <- (summary(model4m)$r.squared) - (summary(model4k)$r.squared) # unique visual r squared over processing speed 
uniquecp
sharedacp <- (summary(model4m)$r.squared) - (uniqueap+uniquecp)
sharedacp

time <- seq(from=0, to=2000, by=4)
####cognitive control main effects 
#cohort target timing 
model1j <- lm(Cohort_targettiming ~ 1+ Cohort_globalstop + Cohort_respinhib, data = data4, na.action = na.omit)
summary(model1j) # global stop predicts target timing, but not response inhib
model1ja <- lm(Cohort_targettiming ~ 1 + Cohort_age_c, data = data4)
summary(model1ja) # age alone predicts tt 
model1jb <- lm(Cohort_targettiming ~ 1 + Cohort_age_c + Cohort_globalstop + Cohort_respinhib, data = data4)
summary(model1jb) # age and global stop both predicts tt 
anova(model1j, model1jb) # age predicts over and above CC 
anova(model1ja, model1jb) # CC predicts over and above age 
uniqueatt <- (summary(model1jb)$r.squared) - (summary(model1j)$r.squared) # R squared change for unique processing speed variance
uniqueatt
uniquecctt <- (summary(model1jb)$r.squared) - (summary(model1ja)$r.squared) # unique visual r squared over processing speed 
uniquecctt
sharedacctt <- (summary(model1jb)$r.squared) - (uniqueatt+uniquecctt)
sharedacctt


####Cohort resolution
model1jj <- lm(Cohort_resolution ~ 1+ Cohort_globalstop + Cohort_respinhib, data = data4, na.action = na.omit)
summary(model1jj) # global stop predicts target resolution, but not response inhib
model1jk <- lm(Cohort_resolution ~ 1 + Cohort_age_c, data = data4)
summary(model1jk) # age alone does not predict resolution
model1jl <- lm(Cohort_resolution ~ 1 + Cohort_age_c + Cohort_globalstop + Cohort_respinhib, data = data4)
summary(model1jl) # only global stop both predicts resolution 
anova(model1jl, model1jj) # age does not predict over and above CC 
anova(model1jl, model1jk) # CC moderately predicts over and above age 
uniqueatt <- (summary(model1jl)$r.squared) - (summary(model1jj)$r.squared) # R squared change for unique processing speed variance
uniqueatt
uniquecctt <- (summary(model1jl)$r.squared) - (summary(model1jk)$r.squared) # unique visual r squared over processing speed 
uniquecctt
sharedacctt <- (summary(model1jl)$r.squared) - (uniqueatt+uniquecctt)
sharedacctt

###### Peak 
data4 <- subset(data4, Cohort_UR>.8)
model1jo <- lm(Cohort_CH-Cohort_UH ~ 1+ Cohort_globalstop + Cohort_respinhib, data = data4, na.action = na.omit)
summary(model1jo) # neither global stop nor response inhib predict peak
model1jp <- lm(Cohort_CH-Cohort_UH ~ 1 + Cohort_age_c, data = data4)
summary(model1jp) # age alone does not predict peak
model1jl <- lm(Cohort_CH-Cohort_UH ~ 1 + Cohort_age_c + Cohort_globalstop + Cohort_respinhib, data = data4)
summary(model1jl)
anova(model1jo, model1jl) # age does not predict over and above CC 
anova(model1jp, model1jl) # CC

uniqueap <- (summary(model1jl)$r.squared) - (summary(model1jo)$r.squared) # R squared change for unique processing speed variance
uniqueap
uniqueccp <- (summary(model1jl)$r.squared) - (summary(model1jp)$r.squared) # unique visual r squared over processing speed 
uniqueccp
sharedaccp <- (summary(model1jl)$r.squared) - (uniqueap+uniqueccp)
sharedaccp

###### Min Time 
model1jy <- lm(Cohort_MinTime ~ 1+ Cohort_globalstop + Cohort_respinhib, data = data4, na.action = na.omit)
summary(model1jy) # neither global stop nor response inhib predict max
model1jz <- lm(Cohort_MinTime ~ 1 + Cohort_age_c, data = data4)
summary(model1jz) # age alone does not predict max
model1jl <- lm(Cohort_MinTime ~ 1 + Cohort_age_c + Cohort_globalstop + Cohort_respinhib, data = data4)
summary(model1jl)
anova(model1jy, model1jl) # age does not predict over and above CC 
anova(model1jz, model1jl) # CC

uniqueaon <- (summary(model1jl)$r.squared) - (summary(model1jy)$r.squared) # R squared change for unique processing speed variance
uniqueaon
uniqueccon <- (summary(model1jl)$r.squared) - (summary(model1jz)$r.squared) # unique visual r squared over processing speed 
uniqueccon
sharedaccon <- (summary(model1jl)$r.squared) - (uniqueaon+uniqueccon)
sharedaccon


###### Max Time 
model1jy <- lm(Cohort_MaxTime ~ 1+ Cohort_globalstop + Cohort_respinhib, data = data4, na.action = na.omit)
summary(model1jy) # neither global stop nor response inhib predict max
model1jz <- lm(Cohort_MaxTime ~ 1 + Cohort_age_c, data = data4)
summary(model1jz) # age alone does not predict max
model1jz1 <- lm(Cohort_MaxTime ~ 1 + Cohort_age_c + Cohort_globalstop + Cohort_respinhib, data = data4)
summary(model1jz1) # only global stop both predicts resolution 
anova(model1jz1, model1jy) # age does not predict over and above CC 
anova(model1jz1, model1jz) # CC does not predicts over and above age 

uniqueaoff <- (summary(model1jz1)$r.squared) - (summary(model1jy)$r.squared) # R squared change for unique processing speed variance
uniqueaoff
uniqueccoff <- (summary(model1jz1)$r.squared) - (summary(model1jz)$r.squared) # unique visual r squared over processing speed 
uniqueccoff
sharedaccoff <- (summary(model1jl)$r.squared) - (uniqueaoff+uniqueccoff)
sharedaccoff

### residuals visualization 
#baseline curve
	#Cohort 
		Cohort_meancross <- mean(data2$Cohort_Crossover) # mean cross
		Cohort_meancross <- as.data.frame(Cohort_meancross)
		data2 <- cbind(data2, Cohort_meancross)
		Cohort_meanslope <- mean(data2$Cohort_Slope) # mean slope
		Cohort_meanslope <- as.data.frame(Cohort_meanslope)
		data2 <- cbind(data2, Cohort_meanslope)
		Cohort_meanmax <- mean(data2$Cohort_Max) # mean max
		Cohort_meanmax <- as.data.frame(Cohort_meanmax)
		data2 <- cbind(data2, Cohort_meanmax)
		Cohort_meanmin <- mean(data2$Cohort_Min) # mean min 
		Cohort_meanmin <- as.data.frame(Cohort_meanmin)
		data2 <- cbind(data2, Cohort_meanmin)



#individual subjects add residuals
	#Cohort 
   		#Crossover
   			resid <- lm(Cohort_Crossover ~ VVWP_Crossover, data = data2)
   			resid.res <- (resid(resid))
			resid.coef <- as.data.frame(coef(resid))
   			resid.res <- as.data.frame(resid.res)
   			data2 <- cbind(data2, resid.res)
   			data2$Cohort_meancross <- as.numeric(data2$Cohort_meancross)
   			data2$Cohort_cross_res <- data2$Cohort_meancross + data2$resid.res
   
   		#Max 
   			resid1 <- lm(Cohort_Max~ VVWP_Max, data = data2)
  			resid.res1 <- (resid(resid1))
			resid.coef1 <- as.data.frame(coef(resid1))
   			resid.res1 <- as.data.frame(resid.res1)
   			data2 <- cbind(data2, resid.res1)
   			data2$Cohort_meanmax <- as.numeric(data2$Cohort_meanmax)
   			data2$Cohort_Max_res <- data2$Cohort_meanmax + data2$resid.res1
   
   		 #Min 
   			resid2 <- lm(Cohort_Min~ VVWP_Min, data = data2)
  			resid.res2 <- (resid(resid2))
			resid.coef2 <-  as.data.frame(coef(resid2))
   			resid.res2 <- as.data.frame(resid.res2)
   			data2 <- cbind(data2, resid.res2)
   			data2$Cohort_meanmin <- as.numeric(data2$Cohort_meanmin)
   			data2$Cohort_Min_res <- data2$Cohort_meanmin + data2$resid.res2
		#Slope 
   			resid3 <- lm(Cohort_Slope~ VVWP_Slope, data = data2)
   			resid.res3 <- (resid(resid3))
			resid.coef3 <- as.data.frame(coef(resid3))
   			resid.res3 <- as.data.frame(resid.res3)
   			data2 <- cbind(data2, resid.res3)
   			data2$Cohort_meanslope <- as.numeric(data2$Cohort_meanslope)
   			data2$Cohort_Slope_res <- data2$Cohort_meanslope + data2$resid.res3

   			
	#Cohort
   			cohlogistic <- function(x){
   			  (data2$Cohort_Max_res[x] - .0000000004)/ (1 + exp(4* 
   			                                                
		(data2$Cohort_Slope_res[x]/(data2$Cohort_Max_res[x]-.0000000004)) *(data2$Cohort_cross_res[x]- time))) + .0000000004}
   			
   			cohcurve <- data.frame(matrix(NA, 
   			                              nrow= 501,
   			                              ncol = 107))
   			rows <- as.numeric(nrow(data2))
   			for (x in 1:rows) {
   			  cohcurve[x] <- cohlogistic(x)}
   	
  

   			# Average logistics within age
   			agegroup <- data2$Semantic_agegroup
   			#Cohort
   			cohcurve <- t(cohcurve)
   			cohcurve <- cbind(cohcurve, agegroup)
   			cohcurve <- as.data.frame(cohcurve)
   			cohavg <- cohcurve %>%
   			  group_by(agegroup) %>%
   			  summarise_at(vars(V1:V501), list(mean))
   			cohavg <- t(cohavg)
   			cohavg <- cohavg[2:502,]
   			colnames(cohavg) <- c("(7-8)", "(11-12)", "(16-17)")
   			
   			
   			## effect of age without the residual visualization 	
   			# create logistic function 
   			# Cohort 
   			cohraw <- function(x){
   			  (data2$Cohort_Max[x] - .0000000004)/ (1 + exp(4* 
   			                                                  (data2$Cohort_Slope[x]/(data2$Cohort_Max[x]-.0000000004)) *(data2$Cohort_Crossover[x]- time))) + .0000000004}
   			
   			
   			
   			# logistic run 
   			# Cohort 
   			cohrawcurve <- data.frame(matrix(NA, 
   			                                 nrow= 501,
   			                                 ncol = 107))
   			rows <- as.numeric(nrow(data2))
   			for (x in 1:rows) {
   			  cohrawcurve[x] <- cohraw(x)}
   		
   			#Cohort by age
   			cohrawcurve <- t(cohrawcurve)
   			cohrawcurve <- cbind(cohrawcurve, agegroup)
   			cohrawcurve <- as.data.frame(cohrawcurve)
   			cohrawavg <- cohrawcurve %>%
   			  group_by(agegroup) %>%
   			  summarise_at(vars(V1:V501), list(mean))
   			cohrawavg <- t(cohrawavg)
   			cohrawavg <- cohrawavg[2:502,]
   			colnames(cohrawavg) <- c("(7-8)", "(11-12)", "(16-17)")
   			
   			#plot 
   	
   			#cohort without VVWP 
   			cohavg <- data.frame(cbind(cohavg, time))
   			ccohavgmelt <- melt(cohavg, id = ("time"))
   			ccohavgmelt$type <- c("B. Cohort Trials Controlled for Vis/Cog Skills")

			cohrawavg <- data.frame(cbind(cohrawavg, time))
			acohrawavgmelt <- melt(cohrawavg, id = ("time"))
   			acohrawavgmelt$type <- c("A. Cohort Trials")
   			
   	#just Cohort plot - USE THIS ONE
   			
   			plot <- data.frame(acohrawavgmelt, ccohavgmelt)
   			plot1 <- data.frame(value = c(plot[,"value"], plot[, "value.1"]))
   			plot1$age <-  c(plot[,"variable"], plot[, "variable.1"])
 	  		plot1$age <- recode(plot1$age, '1' = "7-8", '2' = "11-12", '3' = "16-17")
			  plot1$age <- factor(plot1$age, labels = c("7-8", "11-12", "16-17"))
   			plot1$type <- c(plot[,"type"], plot[, "type.1"])
   			plot1$time <- c(plot[,"time"], plot[,"time.1"])
			library(ggplot2)
   			ggplot(plot1, aes(x=time, y=value, colour=age)) + geom_line(size = 1.25) +
				scale_color_manual(name = "Age Group", values = c('7-8' = "#0F3BE4", '11-12'= "#E60202", '16-17' = "#ED7306")) + 
   			labs(x = "Time (msec)", y = "Proportion of Fixations") + 
				theme_minimal() +
			  facet_wrap(~type) +
  				theme(legend.justification=c(1, .5), legend.position=c(1, .5),
        			 	strip.text.x=element_text(size=14,family="Times New Roman"),
        				axis.title.x = element_text(size=14,family="Times New Roman"),
        				axis.text.x =element_text(size=12,family="Times New Roman"),
        				axis.title.y=element_text(size=14, family="Times New Roman"), 
					      axis.text.y=element_text(size=12, family="Times New Roman"), 
        				legend.title=element_text(size=14, family="Times New Roman"), 
        				legend.text=element_text(size=14, family="Times New Roman"), 
        				panel.grid.major = element_blank(),
        				panel.grid.minor = element_blank(), 
        				axis.line = element_line(),
        				axis.ticks = element_line(), 
					      panel.spacing.x = unit(6, "mm"), 
					      plot.margin = unit(c(.5,.5,.5,.5), "cm")) +
   			  geom_hline(yintercept=0) + 
   			  scale_y_continuous(expand = c(0, 0)) + 
   			  scale_x_continuous(expand = c(0, 0))
  			 
   			  
#just raw 
   			plot <- data.frame(acohrawavgmelt)
   			plot1 <- data.frame(value = c(plot[,"value"]))
   			plot1$age <-  c(plot[,"variable"])
   			plot1$age <- recode(plot1$age, '1' = "(7-8)", '2' = "(11-12)", '3' = "(16-17)")
   			plot1$type <- c(plot[,"type"])
   			plot1$time <- c(plot[,"time"])
   			ggplot(plot1, aes(x=time, y=value, color=age, fill=age)) + geom_line(size = 1.25) + 
   			  xlab("Time (msecs)") + ylab("Proportion of Fixations") + 
   			  scale_colour_discrete("Age Group", labels=c("(7-8)", "(11-12)", "(16-17)")) + 
   			  ggtitle("Proportion of Target Fixations on Cohort Trials")
 # just residualized 
   			
   			plot <- data.frame(ccohavgmelt)
   			plot1 <- data.frame(value = c(plot[, "value.1"]))
   			plot1$age <-  c(plot[, "variable.1"])
   			plot1$age <- recode(plot1$age, '1' = "(7-8)", '2' = "(11-12)", '3' = "(16-17)")
   			plot1$type <- c(plot[, "type.1"])
   			plot1$time <- c(plot[,"time.1"])
   			ggplot(plot1, aes(x=time, y=value, color=age, fill=age)) + geom_line(size = 1.25) +
   			  xlab("Time (msecs)") + ylab("Proportion of Fixations") + 
   			  scale_colour_discrete("Age Group", labels=c("(7-8)", "(11-12)", "(16-17)")) + 
   			  ggtitle("Proportion of Target Fixations on Cohort Trials Residualized by Vis")
   			
   			#Cohort and Semantic plot
   			
   			plot <- data.frame(acohrawavgmelt, bsemrawavgmelt, ccohavgmelt, dsemavgmelt)
   			plot1 <- data.frame(value = c(plot[,"value"], plot[,"value.1"], plot[, "value.2"], plot[,"value.3"]))
   			plot1$age <-  c(plot[,"variable"], plot[,"variable.1"], plot[, "variable.2"], plot[, "variable.3"])
   			plot1$age <- recode(plot1$age, '1' = "young", '2' = "middle", '3' = "teen")
   			plot1$type <- c(plot[,"type"], plot[, "type.1"], plot[, "type.2"], plot[,"type.3"])
   			plot1$time <- c(plot[,"time"], plot[, "time.1"], plot[,"time.2"], plot[,"time.3"])
   			ggplot(plot1, aes(x=time, y=value, color=age)) + geom_line(size = 1.25) +
   			  facet_wrap(~type) + 
   			  xlab("Time (msecs)") + ylab("Proportion of Fixations") + 
   			  scale_color_discrete("Age Group", limits=c("young", "middle", "teen"))
	
			#new plot
			library(ggplot2)
			ggplot(plot1, aes(x = time, y = value)) +
  				geom_line((aes(linetype= type)), size = 1) + 
  				scale_linetype_manual(name = "Age Group", values = c('(7-8)' = "solid", '(11-12)' = "dashed", '(16-17)' = "22" )) + 
  				labs(y = "Proportion of Fixations") +
  				theme_minimal() +
  				theme(legend.justification=c(1, .5), legend.position=c(1, .5),
        				strip.text.x=element_text(size=40,family="Times New Roman"),
        				axis.title.x = element_blank(),
        				axis.text.x =element_blank(),
        				axis.title.y=element_text(size=40,family="Times New Roman"),
        				axis.text.y=element_text(size=40, family="Times New Roman"), 
        				legend.title=element_text(size=40, family="Times New Roman"), 
        				legend.text=element_text(size=40, family="Times New Roman"), 
        				panel.grid.major = element_blank(),
        				panel.grid.minor = element_blank(), 
       				axis.line = element_line(),
        				axis.ticks = element_line(), 
        				axis.line.x = element_blank(), 
        				axis.ticks.x = element_blank()) +
  				geom_hline(yintercept=0)

#create dg function 
	#Semantic 
   			time <- seq(from=0, to=2000, by=4)
   semdg <- function(x){
	  if (time <= data$Semantic_CMu[x]) {
	      ((exp(((time-data$Semantic_CMu[x])^2)/(-2*(data$Semantic_CS1[x]^2))))* (data$Semantic_CH[x]-data$Semantic_CB1[x])) + data$Semantic_CB1[x]}
	  else {
	  ((exp(((time-data$Semantic_CMu[x])^2)/(-2*(data$Semantic_CS2[x]^2))))* (data$Semantic_CH[x]-data$Semantic_CB2[x])) + data$Semantic_CB2[x]}}                                                         
	
	semunrelateddg <- function(x){
	  if (time <= data$Semantic_UMu[x]) {
	    ((exp(((time-data$Semantic_UMu[x])^2)/(-2*(data$Semantic_US1[x]^2))))* (data$Semantic_UH[x]-data$Semantic_UB1[x])) + data$Semantic_UB1[x]}
	  else {
	    ((exp(((time-data$Semantic_UMu[x])^2)/(-2*(data$Semantic_US2[x]^2))))* (data$Semantic_UH[x]-data$Semantic_UB2[x])) + data$Semantic_UB2[x]}}                                                         
	
	#Cohort
	cohdg <- function(x){
	  if (time <= data$Cohort_CMu[x]) {
	    ((exp(((time-data$Cohort_CMu[x])^2)/(-2*(data$Cohort_CS1[x]^2))))* (data$Cohort_CH[x]-data$Cohort_CB1[x])) + data$Cohort_CB1[x]}
	  else {
	    ((exp(((time-data$Cohort_CMu[x])^2)/(-2*(data$Cohort_CS2[x]^2))))* (data$Cohort_CH[x]-data$Cohort_CB2[x])) + data$Cohort_CB2[x]}}                                                         
	
	cohunrelateddg <- function(x){
	  if (time <= data$Cohort_UMu[x]) {
	    ((exp(((time-data$Cohort_UMu[x])^2)/(-2*(data$Cohort_US1[x]^2))))* (data$Cohort_UH[x]-data$Cohort_UB1[x])) + data$Cohort_UB1[x]}
	  else {
	    ((exp(((time-data$Cohort_UMu[x])^2)/(-2*(data$Cohort_US2[x]^2))))* (data$Cohort_UH[x]-data$Cohort_UB2[x])) + data$Cohort_UB2[x]}}                                                         
	# VVWP 
	   vvdg <- function(x){
	  if (time <= data2$VVWP_CMu[x]) {
	    ((exp(((time-data2$VVWP_CMu[x])^2)/(-2*(data2$VVWP_CS1[x]^2))))* (data2$VVWP_CH[x]-data$VVWP_CB1[x])) + data2$VVWP_CB1[x]}
	  else {
	    ((exp(((time-data2$VVWP_CMu[x])^2)/(-2*(data2$VVWP_CS2[x]^2))))* (data2$VVWP_CH[x]-data$VVWP_CB2[x])) + data2$VVWP_CB2[x]}}                                                         
	
	vvunrelateddg <- function(x){
	   if (time <= data2$VVWP_UMu[x]) {
	    ((exp(((time-data2$VVWP_UMu[x])^2)/(-2*(data2$VVWP_US1[x]^2))))* (data2$VVWP_UH[x]-data2$VVWP_UB1[x])) + data2$VVWP_UB1[x]}
	  else {
	    ((exp(((time-data2$VVWP_UMu[x])^2)/(-2*(data2$VVWP_US2[x]^2))))* (data2$VVWP_UH[x]-data2$VVWP_UB2[x])) + data2$VVWP_UB2[x]}} 

# dg run 
	#Semantic 
	semdgdata <- data.frame(matrix(NA, 
						nrow= 501,
						ncol = 117))
	rows <- as.numeric(nrow(data))
	for (x in 1:rows) {
	  semdgdata[x] <- semdg(x)}
	data1 <- semdgdata
	colnames(data1) <- paste("sem", colnames(data1), sep = "_")
	#Semantic Unrelated 
  semunrelateddgdata <- data.frame(matrix(NA, 
                               nrow= 501,
                               ncol = 117))
  rows <- as.numeric(nrow(data))
  for (x in 1:rows) {
    semunrelateddgdata[x] <- semunrelateddg(x)}
  data2 <- semunrelateddgdata
  colnames(data2) <- paste("semunrelated", colnames(data2), sep = "_")
  semdiff <- (data1-data2)

#Cohort
  cohdgdata <- data.frame(matrix(NA, 
                               nrow= 501,
                               ncol = 117))
  rows <- as.numeric(nrow(data))
  for (x in 1:rows) {
    cohdgdata[x] <- cohdg(x)}
  data3 <- cohdgdata
  colnames(data1) <- paste("coh", colnames(data3), sep = "_")     
#Cohort Unrelated 
  cohunrelateddgdata <- data.frame(matrix(NA, 
                                        nrow= 501,
                                        ncol = 117))
  rows <- as.numeric(nrow(data))
  for (x in 1:rows) {
    cohunrelateddgdata[x] <- cohunrelateddg(x)}
  data4 <- cohunrelateddgdata
  colnames(data4) <- paste("cohunrelated", colnames(data4), sep = "_")
  cohdiff <- (data3-data4)

#VVWP
  vvdgdata <- data.frame(matrix(NA, 
                               nrow= 501,
                               ncol = 107))
  rows <- as.numeric(nrow(data2))
  for (x in 1:rows) {
    vvdgdata[x] <- vvdg(x)}
  data5 <- vvdgdata
    
#VVWP Unrelated 
  vvunrelateddgdata <- data.frame(matrix(NA, 
                                        nrow= 501,
                                        ncol = 107))
  rows <- as.numeric(nrow(data2))
  for (x in 1:rows) {
    vvunrelateddgdata[x] <- vvunrelateddg(x)}
  data6 <- vvunrelateddgdata
  colnames(data6) <- paste("vvunrelated", colnames(data6), sep = "_")
  vvdiff <- (data5-data6)
  
#normalize competitor curves 
  #Cohort 
  cohdgmax<- apply(cohdiff, 2, max)
  cohdgmax <- as.data.frame(cohdgmax)
  cohnorm1 <- rep(c(cohdgmax), each=501)
  cohnorm1 <- as.data.frame(cohnorm1)
  cohnorm1 <- t(cohnorm1)
  cohnorm1 <- as.data.frame(cohnorm1)
  cohnorm <- (cohdiff/cohnorm1)
#then export to Access to run Max, Min at .6 and threshold 
 write.csv(cohnorm, file= "cohnorm.csv")
                                                      
  #Semantic 
  semdgmax<- apply(semdiff, 2, max)
  semdgmax <- as.data.frame(semdgmax)
  semnorm1 <- rep(c(semdgmax), each=501)
  semnorm1 <- as.data.frame(semnorm1)
  semnorm1 <- t(semnorm1)
  semnorm1 <- as.data.frame(semnorm1)
  semnorm <- (semdiff/semnorm1)
#then export to Access to run Max, Min at .6 and threshold 
  write.csv(semnorm, file = "semnorm.csv")

  #VVWP
  vvdgmax<- apply(vvdiff, 2, max)
  vvdgmax <- as.data.frame(vvdgmax)
  vvnorm1 <- rep(c(vvdgmax), each=501)
  vvnorm1 <- as.data.frame(vvnorm1)
  vvnorm1 <- t(vvnorm1)
  vvnorm1 <- as.data.frame(vvnorm1)
  vvnorm <- (vvdiff/vvnorm1)
#then export to Access to run Max, Min at .6 and threshold 
  write.csv(vvnorm, file = "vvnorm.csv")






		
		ggplot(semrawavgmelt, aes(x=time, y=value, color=variable)) + geom_line(size = 1.25) + 
		  ggtitle("Semantic Trials Raw")
		
		ggplot(cohrawavgmelt, aes(x=time, y=value, color=variable)) + geom_line(size = 1.25) + 
		  ggtitle("Cohort Trials Raw")
		

