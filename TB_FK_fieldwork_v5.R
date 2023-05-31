###################################################################
######## Combining all FK and TB trophic level projects ########

#for Mac
setwd("~/Box Sync/R work/TB_FK_fieldwork")
search()

#loading packages
install.packages("vegan")
install.packages("tidyverse")
install.packages("codyn")
install.packages("ggplot2")
install.packages("reshape2")
install.packages("plotly")
install.packages("nmle")
install.packages("lme4")
install.packages("olsrr")
install.packages("car")
install.packages("patchwork")
install.packages("lmerTest")
install.packages("piecewiseSEM")
install.packages("multcomp")
install.packages("MuMIn")
install.packages("forcats")
install.packages("emmeans")

library(vegan)
library(tidyverse)
library(codyn)
library(ggplot2)
library(reshape2)
library(plotly)
library(nlme)
library(lme4)
library(olsrr)
library(car)
library(patchwork)
library(lmerTest)
library(piecewiseSEM)
library(multcomp)
library(MuMIn)
library(forcats)
library(emmeans)

#Set ggplot2 theme to black and white
theme_set(theme_bw())
#Update ggplot2 theme - make box around the x-axis title size 30, vertically justify x-axis title to 0.35, 
#Place a margin of 15 around the x-axis title.  
#Make the x-axis title size 30. For y-axis title, make the box size 30, put the writing at a 90 degree angle, and vertically justify the title to 0.5.  
#Add a margin of 15 and make the y-axis text size 25. Make the plot title size 30 and vertically justify it to 2.  Do not add any grid lines.  
#Do not add a legend title, and make the legend size 20
theme_update(axis.title.x=element_text(size=20, vjust=-0.35, margin=margin(t=12)),
             axis.text.x=element_text(size=20), axis.title.y=element_text(size=20, angle=90, vjust=0.5,
                                                                          margin=margin(r=15)), axis.text.y=element_text(size=20), plot.title =
               element_text(size=20, vjust=2), panel.grid.major=element_blank(),
             panel.grid.minor=element_blank(),
             legend.text=element_text(size=15))



#set colorblind friendly color palette
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#000000", "#CC79A7")

####################################################################



#########################################################################
########### Creating functions #########################

#repeated measures anova with 2 independent variables
anova_t3 <- function(IndVars=IndVars, DepVar=DepVar, RndForm=RndForm, Data=Data){
  anova_out <- {}
  IndVarMatrix <- matrix(nrow=length(IndVars),ncol=length(IndVars))
  IndVars2x <- c(IndVars,IndVars)
  
  for(REORDER in 1:length(IndVars)){
    IndVarMatrix[REORDER,] <- IndVars2x[REORDER:(length(IndVars)+(REORDER-1))]
  }
  rm(IndVars2x)
  
  for(RUN in 1:length(IndVars)){
    model_formula_temp <- paste0(DepVar,"~", paste0(IndVarMatrix[RUN,], collapse="*"))
    model_temp <- lme(as.formula(model_formula_temp)
                      , data=Data
                      , random = as.formula(RndForm)
                      , correlation=corCompSymm(form = as.formula(RndForm))
                      , control=lmeControl(returnObject=TRUE)
                      , na.action = na.omit)
    anova_out_temp <- anova(model_temp)
    
    if(length(IndVars)==2){ ## Pulls model output for variable that is last
      if(RUN==1){ ### This currently works for 
        anova_partial_temp <- rbind(anova_out_temp[(length(IndVars)+1),],
                                    anova_out_temp[4,]
        )
      }
      if(RUN %in% 2:length(IndVars)){
        anova_partial_temp <- anova_out_temp[(length(IndVars)+1),]
        
      }
    } # End if vars == 2 statement
    anova_out <- rbind(anova_out, anova_partial_temp)
  } # End reorder loop
  return(anova_out)
} # End function



#repeated measures anova with 3 independent variables
anova_t3_3ind <- function(IndVars=IndVars, DepVar=DepVar, RndForm=RndForm, Data=Data){
  anova_out <- {}
  IndVarMatrix <- matrix(nrow=length(IndVars),ncol=length(IndVars))
  IndVars2x <- c(IndVars,IndVars)
  
  for(REORDER in 1:length(IndVars)){
    IndVarMatrix[REORDER,] <- IndVars2x[REORDER:(length(IndVars)+(REORDER-1))]
  }
  rm(IndVars2x)
  
  for(RUN in 1:length(IndVars)){
    model_formula_temp <- paste0(DepVar,"~", paste0(IndVarMatrix[RUN,], collapse="*"))
    model_temp <- lme(as.formula(model_formula_temp)
                      , data=Data
                      , random = as.formula(RndForm)
                      , correlation=corCompSymm(form = as.formula(RndForm))
                      , control=lmeControl(returnObject=TRUE)
                      , na.action = na.omit)
    anova_out_temp <- anova(model_temp)
    
    if(length(IndVars)==3){ ## Pulls model output for variable that is last
      if(RUN==1){ ### This currently works for 
        anova_partial_temp <- rbind(anova_out_temp[(length(IndVars)+1),],
                                    anova_out_temp[7:8,]
        )
      }
      if(RUN %in% 2:length(IndVars)){
        anova_partial_temp <- rbind(anova_out_temp[(length(IndVars)+1),],
                                    anova_out_temp[7,]
        )
      }
    } # End if vars==3 statement
    
    if(length(IndVars)==2){ ## Pulls model output for variable that is last
      if(RUN==1){ ### This currently works for 
        anova_partial_temp <- rbind(anova_out_temp[(length(IndVars)+1),],
                                    anova_out_temp[4,]
        )
      }
      if(RUN %in% 2:length(IndVars)){
        anova_partial_temp <- anova_out_temp[(length(IndVars)+1),]
        
      }
    } # End if vars == 2 statement
    anova_out <- rbind(anova_out, anova_partial_temp)
  } # End reorder loop
  return(anova_out)
} # End function

##################################################




##################################################
####### Read in data ########

#metadata
FK_plotinfo <- read.csv("FK_plotinfo.csv", header = TRUE)
TB_plotinfo <- read.csv("TB_plotinfo.csv", header = TRUE)


#plant data
TBspcomp2019 <- read.csv("TB_2019_speciescomp.csv", header = TRUE)
FKspcomp2020 <- read.csv("FK_2020_speciescomp.csv", header = TRUE)
FKspcomp2021 <- read.csv("FK_2021_speciescomp.csv", header = TRUE)
FKspcomp2022 <- read.csv("FK_2022_speciescomp.csv", header = TRUE)

TBspecies <- read.csv("speciesinfo_TB.csv", header = TRUE)
FKspecies <- read.csv("speciesinfo_FK.csv", header = TRUE)


#insect data
TBinsID2019 <- read.csv("TB_2019_insectID.csv", header = TRUE, na.strings = "")
TBguild <- read.csv("TB_insect_feeding_guilds.csv", header = TRUE)
TBbio <- read.csv("TB_2019_insectbiomass.csv", header = TRUE)

FKinsID2020 <- read.csv("FK_2020_insectID.csv", header = TRUE, na.strings = "")
FKguild <- read.csv("FK_insect_feeding_guilds.csv", header = TRUE)
FK2020bio <- read.csv("FK_2020_insectbiomass.csv", header = TRUE)

FKinsID2021 <- read.csv("FK_2021_insectID.csv", header = TRUE, na.strings = "")
FK2021bio <- read.csv("FK_2021_insectbiomass.csv", header = TRUE)
FK2021herb <- read.csv("FK_2021_herbivory.csv", header = TRUE)

FK2022bio <- read.csv("FK_2022_insectbiomass.csv", header = TRUE)
FK2022herb <- read.csv("FK_2022_herbivory.csv", header = TRUE)

#microbe data
TB_mic_rich <- read.csv("TB_microbe_richness.csv", header = TRUE)  #observed features as richness
TBmicord <- read.delim("TBdistance-matrix.tsv", header = TRUE, row.names = NULL) %>%
  rename(id = X)
mapping <- read.csv("mapping.csv", header = TRUE)
mic_func <- read.delim("func_table3.tsv", header=TRUE, row.names=NULL) 
#FK plot 6 should always be removed due to low # of sequences
FK_mic_rich <- read.csv("FK_2020_microbe_richness.csv", header = TRUE)  #observed features as richness
FKmicord <- read.delim("FKdistance-matrix.tsv", header = TRUE, row.names = NULL) %>%
  rename(id = X)

###################################################################


###################################################################
######## Clean up sp comp ##########

#TB 2019
TBspcomp2019_2 <- TBspcomp2019 %>%
  filter(aerial_basal == "aerial") %>%
  dplyr::select(-c(date, est_tot, add_tot, add_tot_excel, cov_other_tot, cov_other_tot_excel, bareground, dung, lichen, moss, rock, litter))

TBspcomp2019_long <- TBspcomp2019_2 %>%
  group_by(location, year, grad_num, plot) %>%
  gather(symbol, cover, 6:ncol(TBspcomp2019_2)) %>%
  ungroup()

TB2019totcov <- TBspcomp2019_long %>%
  group_by(plot) %>%
  summarise(totcov = sum(cover)) %>%
  ungroup()

TB2019cover <- full_join(TBspcomp2019_long, TB2019totcov) %>%
  mutate(rel_cov = (cover/totcov) * 100) %>%
  full_join(TB_plotinfo) %>%
  dplyr::select(-c(aerial_basal, plot_type, plot_name, phone_lat, phone_long, garmin_lat, garmin_long))


#FK 2020
FKspcomp2020_2 <- FKspcomp2020 %>%
  filter(aerial_basal == "aerial") %>%
  dplyr::select(-c(date, est_tot, add_tot, add_tot_excel, cov_other_tot, cov_other_tot_excel, bareground, dung, lichen, litter, moss, rock))

FKspcomp2020_2[is.na(FKspcomp2020_2)] <- 0 #replace na with 0

FKspcomp2020_long <- FKspcomp2020_2 %>%
  group_by(location, year, grad_num, plot) %>%
  gather(symbol, cover, 6:ncol(FKspcomp2020_2)) %>%
  ungroup()

FK2020totcov <- FKspcomp2020_long %>%
  group_by(plot) %>%
  summarise(totcov = sum(cover)) %>%
  ungroup()

FK2020cover <- full_join(FKspcomp2020_long, FK2020totcov) %>%
  mutate(rel_cov = (cover/totcov) * 100) %>%
  full_join(FK_plotinfo) %>%
  dplyr::select(-c(aerial_basal, plot_type, plot_name, latitude, longitude))

#what plants are most dominant across plots
FK2020dom <- FK2020cover %>%
  group_by(symbol, invasion_percent) %>%
  summarise(avg = mean(cover)) %>%
  ungroup()  #HECO is most dominant across all plots/inv %

#FK 2021
FKspcomp2021_2 <- FKspcomp2021 %>%
  filter(aerial_basal == "aerial") %>%
  dplyr::select(-c(date, fried_litter, est_tot, add_tot, rock, dung, moss, lichen, mushroom, litter, bareground, final_tot))

FKspcomp2021_2[is.na(FKspcomp2021_2)] <- 0 #replace na with 0

FKspcomp2021_long <- FKspcomp2021_2 %>%
  group_by(location, year, grad_num, plot) %>%
  gather(symbol, cover, 6:ncol(FKspcomp2021_2)) %>%
  ungroup()

FK2021totcov <- FKspcomp2021_long %>%
  group_by(plot) %>%
  summarise(totcov = sum(cover)) %>%
  ungroup()

FK2021cover <- full_join(FKspcomp2021_long, FK2021totcov) %>%
  mutate(rel_cov = (cover/totcov) * 100) %>%
  full_join(FK_plotinfo) %>%
  filter(plot_type == "C") %>%
  dplyr::select(-c(aerial_basal, plot_type, plot_name, latitude, longitude))


#FK 2022
FKspcomp2022_2 <- FKspcomp2022 %>%
  filter(aerial_basal == "aerial") %>%
  dplyr::select(-c(date, est_tot, add_tot, rock, dung, moss, lichen, litter, bareground, final_tot))

FKspcomp2022_2[is.na(FKspcomp2022_2)] <- 0 #replace na with 0

FKspcomp2022_long <- FKspcomp2022_2 %>%
  group_by(location, year, grad_num, plot) %>%
  gather(symbol, cover, 6:ncol(FKspcomp2022_2)) %>%
  ungroup()

FK2022totcov <- FKspcomp2022_long %>%
  group_by(plot) %>%
  summarise(totcov = sum(cover)) %>%
  ungroup()

FK2022cover <- full_join(FKspcomp2022_long, FK2022totcov) %>%
  mutate(rel_cov = (cover/totcov) * 100) %>%
  full_join(FK_plotinfo) %>%
  filter(plot_type == "C") %>%
  dplyr::select(-c(aerial_basal, plot_type, plot_name, latitude, longitude))


#combine FK data
FKcover <- FK2020cover%>%
  full_join(FK2021cover) %>%
  full_join(FK2022cover)

###########################################################################


########################################################
######### Gradient effectiveness ########
#TB 2019
TB_BRAR_2019 <- TB2019cover %>%
  filter(invasive_type == "BRAR" & symbol == "BRAR") %>%
  group_by(year, invasion_percent) %>%
  summarise(avg_cov = mean(cover), se_cov = sd(cover)/sqrt(length(cover)), avg_tot = mean(totcov), se_tot = sd(totcov)/sqrt(length(totcov)), avg_rel = mean(rel_cov), se_rel = sd(rel_cov)/sqrt(length(rel_cov))) %>%
  ungroup()

TB_BRTE_2019 <- TB2019cover %>%
  filter(invasive_type == "BRTE" & symbol == "BRTE") %>%
  group_by(year, invasion_percent) %>%
  summarise(avg_cov = mean(cover), se_cov = sd(cover)/sqrt(length(cover)), avg_tot = mean(totcov), se_tot = sd(totcov)/sqrt(length(totcov)), avg_rel = mean(rel_cov), se_rel = sd(rel_cov)/sqrt(length(rel_cov))) %>%
  ungroup()

#FK
FK_BRAR <- FKcover %>%
  filter(symbol == "BRAR") %>%
  group_by(year, invasion_percent) %>%
  summarise(avg_cov = mean(cover), se_cov = sd(cover)/sqrt(length(cover)), avg_tot = mean(totcov), se_tot = sd(totcov)/sqrt(length(totcov)), avg_rel = mean(rel_cov), se_rel = sd(rel_cov)/sqrt(length(rel_cov))) %>%
  ungroup()



#boxplot
FK_BRAR_2020_box <- FKcover %>%
  filter(symbol == "BRAR") %>%
  filter(year == "2020")

TB_BRAR_2019_box <- TB2019cover %>%
  filter(invasive_type == "BRAR" & symbol == "BRAR") 

TB_BRTE_2019_box <- TB2019cover %>%
  filter(invasive_type == "BRTE" & symbol == "BRTE") 


FK_BRAR_box <- ggplot(data = FK_BRAR_2020_box, aes(x = factor(invasion_percent), y = rel_cov)) + 
  geom_boxplot(aes(group = invasion_percent)) + 
  #geom_point(size = 2, color = "blue", position = position_dodge(.05)) +
  ylim(0, 85) +
  labs(x = "Invasion Level (%)", y = "% Cover", title = "MT BRAR") +
  theme(plot.title = element_text(vjust = - 7, hjust = 0.07, size = 16))

TB_BRAR_box <- ggplot(data = TB_BRAR_2019_box, aes(x = factor(invasion_percent), y = rel_cov)) + 
  geom_boxplot(aes(group = invasion_percent)) + 
  #geom_point(size = 2, color = "blue", position = position_dodge(.05)) +
  ylim(0, 85) +
  labs(x = "Invasion Level (%)", y = "", title = "WY BRAR") +
  theme(plot.title = element_text(vjust = - 7, hjust = 0.07, size = 16))

TB_BRTE_box <- ggplot(data = TB_BRTE_2019_box, aes(x = factor(invasion_percent), y = rel_cov)) + 
  geom_boxplot(aes(group = invasion_percent)) + 
  ylim(0, 85) +
  #geom_point(size = 2, color = "blue", position = position_dodge(.05)) +
  labs(x = "Invasion Level (%)", y = "", title = "WY BRTE") +
  theme(plot.title = element_text(vjust = - 7, hjust = 0.07, size = 16))

FK_BRAR_box + TB_BRAR_box + TB_BRTE_box



#statistics

#TB 2019
#BRAR
TBBRAR <- TB2019cover %>%
  filter(invasive_type == "BRAR" & symbol == "BRAR") 

TBBRAR2 <- TBBRAR %>%
  mutate(invasion_percent = factor(invasion_percent))

TB_BRAR_stats <- lmerTest::lmer(data = TBBRAR2, rel_cov ~ invasion_percent +
                                  (1|grad_num))   #rel cov
anova(TB_BRAR_stats, type = 3) #p = 0.007062

summary(glht(TB_BRAR_stats, linfct = mcp(invasion_percent = "Tukey")), test = adjusted(type = "BH")) 
#0-50*, 0-75, 0-100, 25-75, 25-100

#check normality of residuals
resid1 <- lm(data = TBBRAR2, rel_cov ~ invasion_percent)
ols_plot_resid_hist(resid1) #right skew
ols_test_normality(resid1) #pass KS

TBBRAR3 <- TBBRAR2 %>%
  mutate(ln = log10(rel_cov + 1))

resid1 <- lm(data = TBBRAR3, ln ~ invasion_percent)
ols_plot_resid_hist(resid1) #more normal, go with transformed
ols_test_normality(resid1) #pass KS

TB_BRAR_stats2 <- lmerTest::lmer(data = TBBRAR3, ln ~ invasion_percent +
                                  (1|grad_num))   #rel cov
anova(TB_BRAR_stats2, type = 3) #p = 1.752e-06

summary(glht(TB_BRAR_stats2, linfct = mcp(invasion_percent = "Tukey")), test = adjusted(type = "BH")) 
#0-25, 0-50, 0-75, 0-100, 25-75, 25-100, 50-100


#linearity
plot(resid(TB_BRAR_stats2), TBBRAR3$ln) #looks linear
plot(TB_BRAR_stats2) #no pattern so indicates linearity

#homoscedascity
TBBRAR3$res <- residuals(TB_BRAR_stats2)
TBBRAR3$abs_res <- abs(TBBRAR3$res)
TBBRAR3$abs_res2 <- TBBRAR3$abs_res^2
levene_brar <- lm(data = TBBRAR3, abs_res2 ~ ln)
anova(levene_brar) #p = 0.8459 so > 0.05 so equal variance is met

#autocorrelation
acf(residuals(TB_BRAR_stats2, retype = "normalized"))  #lines dont go outside CI horizontal lines, so not autocorrelated
pacf(residuals(TB_BRAR_stats2, retype = "normalized"))




#BRTE
TBBRTE <- TB2019cover %>%
  filter(invasive_type == "BRTE" & symbol == "BRTE") 

TBBRTE2 <- TBBRTE %>%
  mutate(invasion_percent = factor(invasion_percent))

TB_BRTE_stats <- lmerTest::lmer(data = TBBRTE2, rel_cov ~ invasion_percent +
                                  (1|grad_num))   #rel cov
anova(TB_BRTE_stats, type = 3) #p = 4.862e-7

summary(glht(TB_BRTE_stats, linfct = mcp(invasion_percent = "Tukey")), test = adjusted(type = "BH")) 
#0-50, 0-75, 0-100, 25-75, 25-100, 50-75, 50-100, 75-100*


#check normality of residuals
resid1 <- lm(data = TBBRTE2, rel_cov ~ invasion_percent)
ols_plot_resid_hist(resid1) #normalish
ols_test_normality(resid1) #pass KS

TBBRTE3 <- TBBRTE2 %>%
  mutate(ln = log10(rel_cov + 1))

resid1 <- lm(data = TBBRTE3, ln ~ invasion_percent)
ols_plot_resid_hist(resid1) #more normal, go with transformed
ols_test_normality(resid1) #pass

TB_BRTE_stats2 <- lmerTest::lmer(data = TBBRTE3, ln ~ invasion_percent +
                                   (1|grad_num))   #rel cov
anova(TB_BRTE_stats2, type = 3) #p = 3.758e-11

summary(glht(TB_BRTE_stats2, linfct = mcp(invasion_percent = "Tukey")), test = adjusted(type = "BH")) 
#0-25, 0-50, 0-75, 0-100, 25-50, 25-75, 25-100, 50-75, 50-100


#linearity
plot(resid(TB_BRTE_stats2), TBBRTE3$ln) #looks linear
plot(TB_BRTE_stats2) #no pattern so indicates linearity

#homoscedascity
TBBRTE3$res <- residuals(TB_BRTE_stats2)
TBBRTE3$abs_res <- abs(TBBRTE3$res)
TBBRTE3$abs_res2 <- TBBRTE3$abs_res^2
levene_brte <- lm(data = TBBRTE3, abs_res2 ~ ln)
anova(levene_brte) #p = 0.4019 so > 0.05 so equal variance is met

#autocorrelation
acf(residuals(TB_BRTE_stats2, retype = "normalized"))  #lines dont go outside CI horizontal lines, so not autocorrelated
pacf(residuals(TB_BRTE_stats2, retype = "normalized"))



#FK  
FKBRAR <- FKcover %>%
  filter(symbol == "BRAR")

FKBRAR2 <- FKBRAR %>%
  mutate(invasion_percent = factor(invasion_percent))

FK_BRAR_stats <- anova_t3(IndVars = c("year", "invasion_percent"), 
                          DepVar = "rel_cov",       #rel cov
                          RndForm = "~1 | grad_num", 
                          Data = FKBRAR) #inv% p < 0.0001, year p < 0.0001, year*inv% p < 0.0001

FK_BRAR_stats2 <- anova_t3(IndVars = c("year", "invasion_percent"), 
                           DepVar = "cover",       #abs cov
                           RndForm = "~1 | grad_num", 
                           Data = FKBRAR) #inv% p < 0.0001, year p < 0.0001, year*inv% p < 0.0001


FK_BRAR2020_stats <- lmerTest::lmer(data = subset(FKBRAR2, year == "2020"), rel_cov ~ invasion_percent +
                                      (1|grad_num))   #rel cov
anova(FK_BRAR2020_stats, type = 3) #p = 1.017e-13

summary(glht(FK_BRAR2020_stats, linfct = mcp(invasion_percent = "Tukey")), test = adjusted(type = "BH")) 
#all sig different

#check normality of residuals
resid1 <- lm(data = subset(FKBRAR2, year == "2020"), rel_cov ~ invasion_percent)
ols_plot_resid_hist(resid1) #normalish
ols_test_normality(resid1) #pass KS


FKBRAR3 <- FKBRAR2 %>%
  filter(year == "2020") %>%
  mutate(ln = log10(rel_cov + 1))

resid1 <- lm(data = FKBRAR3, ln ~ invasion_percent)
ols_plot_resid_hist(resid1) #more normal, go with transformed
ols_test_normality(resid1) #pass KS

FK_BRAR_stats2 <- lmerTest::lmer(data = FKBRAR3, ln ~ invasion_percent +
                                   (1|grad_num))   #rel cov
anova(FK_BRAR_stats2, type = 3) #p < 2.2e-16

summary(glht(FK_BRAR_stats2, linfct = mcp(invasion_percent = "Tukey")), test = adjusted(type = "BH")) 
#0-50, 0-100, 50-100


#linearity
plot(resid(FK_BRAR_stats2), FKBRAR3$ln) #looks linear
plot(FK_BRAR_stats2) #no pattern so indicates linearity

#homoscedascity
FKBRAR3$res <- residuals(FK_BRAR_stats2)
FKBRAR3$abs_res <- abs(FKBRAR3$res)
FKBRAR3$abs_res2 <- FKBRAR3$abs_res^2
levene_brar2 <- lm(data = FKBRAR3, abs_res2 ~ ln)
anova(levene_brar2) #p = 0.08751 so > 0.05 so equal variance is met

#autocorrelation
acf(residuals(FK_BRAR_stats2, retype = "normalized"))  #lines dont go outside CI horizontal lines, so not autocorrelated
pacf(residuals(FK_BRAR_stats2, retype = "normalized"))



FK_BRAR2021_stats <- lmerTest::lmer(data = subset(FKBRAR, year == "2021"), rel_cov ~ invasion_percent +
                                      (1|grad_num))   #rel cov
anova(FK_BRAR2021_stats, type = 3) #p = 1.227e-5

#check normality of residuals
resid1 <- lm(data = subset(FKBRAR, year == "2021"), rel_cov ~ invasion_percent)
ols_plot_resid_hist(resid1) #normal
ols_test_normality(resid1) #pass


FK_BRAR2022_stats <- lmerTest::lmer(data = subset(FKBRAR, year == "2022"), rel_cov ~ invasion_percent +
                                      (1|grad_num))   #rel cov
anova(FK_BRAR2022_stats, type = 3) #p = 0.0002784

#check normality of residuals
resid1 <- lm(data = subset(FKBRAR, year == "2022"), rel_cov ~ invasion_percent)
ols_plot_resid_hist(resid1) #normal
ols_test_normality(resid1) #pass



## yes, gradients do hold; use untransformed

###################################################################



####################################################################################
############ Plant, insect, soil microbes species richness #############
#Thunder Basin
#Plant species richness
TB_comm <- community_structure(df = TB2019cover,
                               abundance.var = "cover", 
                               replicate.var = "plot", 
                               metric = "Evar") %>%
  full_join(TB_plotinfo)

#add in BRAR and BRTE cover columns
TB_commBRAR <- TB2019cover %>%
  filter(invasive_type == "BRAR" & symbol == "BRAR") %>%
  mutate(BRARcov = cover) %>%
  mutate(BRARrel = rel_cov) %>%
  dplyr::select(c(plot, invasion_percent, BRARcov, BRARrel)) %>%
  full_join(TB_comm) %>%
  drop_na(BRARcov) %>%
  rename(plant_rich = richness, plant_even = Evar)

TB_commBRTE <- TB2019cover %>%
  filter(invasive_type == "BRTE" & symbol == "BRTE") %>%
  mutate(BRTEcov = cover) %>%
  mutate(BRTErel = rel_cov) %>%
  dplyr::select(c(plot, invasion_percent, BRTEcov, BRTErel)) %>%
  full_join(TB_comm) %>%
  drop_na(BRTEcov) %>%
  rename(plant_rich = richness, plant_even = Evar)

#statistics
#TB richness
TBBRARrichstats <- lmerTest::lmer(data = TB_commBRAR, plant_rich ~ BRARrel +
                                    (1|grad_num)) 

anova(TBBRARrichstats, type = 3) #p = 0.04586
rsquared(TBBRARrichstats) #marg = 0.1661731, cond = 0.1806797

TBBRTErichstats <- lmerTest::lmer(data = TB_commBRTE, plant_rich ~ BRTErel +
                                    (1|grad_num)) 

anova(TBBRTErichstats, type = 3) #p = 0.02282
rsquared(TBBRTErichstats)  #r2 = 0.1987397

#normality check
res_a <- lm(data = TB_commBRAR, plant_rich ~ BRARrel)
ols_plot_resid_hist(res_a) #normalish
ols_test_normality(res_a) #pass

res_b <- lm(data = TB_commBRTE, plant_rich ~ BRTErel)
ols_plot_resid_hist(res_b) #normalish
ols_test_normality(res_b) #pass

#linearity
plot(resid(TBBRARrichstats), TB_commBRAR$BRARrel) #looks linear
plot(TBBRARrichstats) #no pattern so indicates linearity

#homoscedascity
TB_commBRAR$res <- residuals(TBBRARrichstats)
TB_commBRAR$abs_res <- abs(TB_commBRAR$res)
TB_commBRAR$abs_res2 <- TB_commBRAR$abs_res^2
levene_brar_rich <- lm(data = TB_commBRAR, abs_res2 ~ BRARrel)
anova(levene_brar_rich) #p = 0.8898 so > 0.05 so equal variance is met

#autocorrelation
acf(residuals(TBBRARrichstats, retype = "normalized"))  #lines dont go outside CI horizontal lines, so not autocorrelated
pacf(residuals(TBBRARrichstats, retype = "normalized"))


#linearity
plot(resid(TBBRTErichstats), TB_commBRTE$BRTErel) #looks linear
plot(TBBRTErichstats) #no pattern so indicates linearity

#homoscedascity
TB_commBRTE$res <- residuals(TBBRTErichstats)
TB_commBRTE$abs_res <- abs(TB_commBRTE$res)
TB_commBRTE$abs_res2 <- TB_commBRTE$abs_res^2
levene_brte_rich <- lm(data = TB_commBRTE, abs_res2 ~ BRTErel)
anova(levene_brte_rich) #p = 0.4564 so > 0.05 so equal variance is met

#autocorrelation
acf(residuals(TBBRTErichstats, retype = "normalized"))  #lines dont go outside CI horizontal lines, so not autocorrelated
pacf(residuals(TBBRTErichstats, retype = "normalized"))


#Insect family richness
TBinsID <- TBinsID2019 %>%
  group_by(location, year, plot, order, family) %>%
  summarise(count = length(sample)) %>%
  ungroup()
TBinsID2 <- TBinsID %>%
  drop_na(family)

TB_inscomm <- community_structure(df = TBinsID2,
                                  abundance.var = "count", 
                                  replicate.var = "plot", 
                                  metric = "Evar") %>%
  full_join(TB_plotinfo) %>%
  drop_na(richness) %>%
  rename(ins_rich = richness, ins_even = Evar)

TB_inscomm_BRAR <- TB_inscomm %>%
  filter(invasive_type == "BRAR") %>%
  full_join(TB_commBRAR)

TB_inscomm_BRTE <- TB_inscomm %>%
  filter(invasive_type == "BRTE") %>%
  full_join(TB_commBRTE)

#statistics
TBBRARinsrichstats <- lmerTest::lmer(data = TB_inscomm_BRAR, ins_rich ~ BRARrel +
                                       (1|grad_num)) 

anova(TBBRARinsrichstats, type = 3) #p = 0.5815

TBBRTEinsrichstats <- lmerTest::lmer(data = TB_inscomm_BRTE, ins_rich ~ BRTErel +
                                       (1|grad_num)) 

anova(TBBRTEinsrichstats, type = 3) #p = 0.05783
rsquared(TBBRTEinsrichstats)  #marg = 0.1168212, cond = 0.3361565

#normality check - use untransformed data
res_a <- lm(data = TB_inscomm_BRAR, ins_rich ~ BRARrel)
ols_plot_resid_hist(res_a) #normalish
ols_test_normality(res_a) #pass

res_b <- lm(data = TB_inscomm_BRTE, ins_rich ~ BRTErel)
ols_plot_resid_hist(res_b) #normalish
ols_test_normality(res_b) #pass


#linearity
TB_inscomm_BRAR <- TB_inscomm_BRAR %>%
  drop_na(ins_rich)
plot(resid(TBBRARinsrichstats), TB_inscomm_BRAR$ins_rich) #looks linear
plot(TBBRARinsrichstats) #no pattern so indicates linearity

#homoscedascity
TB_inscomm_BRAR$res <- residuals(TBBRARinsrichstats)
TB_inscomm_BRAR$abs_res <- abs(TB_inscomm_BRAR$res)
TB_inscomm_BRAR$abs_res2 <- TB_inscomm_BRAR$abs_res^2
levene_brar_insrich <- lm(data = TB_inscomm_BRAR, abs_res2 ~ BRARrel)
anova(levene_brar_insrich) #p = 0.256 so > 0.05 so equal variance is met

#autocorrelation
acf(residuals(TBBRARinsrichstats, retype = "normalized"))  #lines dont go outside CI horizontal lines, so not autocorrelated
pacf(residuals(TBBRARinsrichstats, retype = "normalized"))

#linearity
plot(resid(TBBRTEinsrichstats), TB_inscomm_BRTE$ins_rich) #looks linear
plot(TBBRTEinsrichstats) #no pattern so indicates linearity

#homoscedascity
TB_inscomm_BRTE$res <- residuals(TBBRTEinsrichstats)
TB_inscomm_BRTE$abs_res <- abs(TB_inscomm_BRTE$res)
TB_inscomm_BRTE$abs_res2 <- TB_inscomm_BRTE$abs_res^2
levene_brte_insrich <- lm(data = TB_inscomm_BRTE, abs_res2 ~ BRTErel)
anova(levene_brte_insrich) #p = 0.232 so > 0.05 so equal variance is met

#autocorrelation
acf(residuals(TBBRARinsrichstats, retype = "normalized"))  #lines dont go outside CI horizontal lines, so not autocorrelated
pacf(residuals(TBBRARinsrichstats, retype = "normalized"))



###Microbe ASV richness 
TB_micrich_BRAR <- TB_mic_rich %>%
  filter(invasive_type == "BRAR") %>%
  full_join(TB_inscomm_BRAR)

TB_micrich_BRTE <- TB_mic_rich %>%
  filter(invasive_type == "BRTE") %>%
  full_join(TB_inscomm_BRTE)

#statistics
TB_micrich_BRARstats <- lmerTest::lmer(data = TB_micrich_BRAR, observed_features ~ BRARrel +
                                         (1|grad_num)) 

anova(TB_micrich_BRARstats, type = 3) #p = 0.5933

TB_micrich_BRTEstats <- lmerTest::lmer(data = TB_micrich_BRTE, observed_features ~ BRTErel +
                                         (1|grad_num)) 

anova(TB_micrich_BRTEstats, type = 3) #p = 0.06406
rsquared(TB_micrich_BRTEstats)  #r2 = 0.136205

#normality check - use untransformed data
res_a <- lm(data = TB_micrich_BRAR, observed_features ~ BRARrel)
ols_plot_resid_hist(res_a) #normalish
ols_test_normality(res_a) #pass

res_b <- lm(data = TB_micrich_BRTE, observed_features ~ BRTErel)
ols_plot_resid_hist(res_b) #normalish
ols_test_normality(res_b) #pass


#linearity
TB_micrich_BRAR <- TB_micrich_BRAR %>%
  drop_na(observed_features) %>%
  drop_na(BRARrel)
plot(resid(TB_micrich_BRARstats), TB_micrich_BRAR$observed_features) #looks linear
plot(TB_micrich_BRARstats) #no pattern so indicates linearity

#homoscedascity
TB_micrich_BRAR$res <- residuals(TB_micrich_BRARstats)
TB_micrich_BRAR$abs_res <- abs(TB_micrich_BRAR$res)
TB_micrich_BRAR$abs_res2 <- TB_micrich_BRAR$abs_res^2
levene_brar_micrich <- lm(data = TB_micrich_BRAR, abs_res2 ~ BRARrel)
anova(levene_brar_micrich) #p = 0.5796 so > 0.05 so equal variance is met

#autocorrelation
acf(residuals(TB_micrich_BRARstats, retype = "normalized"))  #lines dont go outside CI horizontal lines, so not autocorrelated
pacf(residuals(TB_micrich_BRARstats, retype = "normalized"))

#linearity
TB_micrich_BRTE <- TB_micrich_BRTE %>%
  drop_na(observed_features) %>%
  drop_na(BRTErel)
plot(resid(TB_micrich_BRTEstats), TB_micrich_BRTE$observed_features) #looks linear
plot(TB_micrich_BRTEstats) #no pattern so indicates linearity

#homoscedascity
TB_micrich_BRTE$res <- residuals(TB_micrich_BRTEstats)
TB_micrich_BRTE$abs_res <- abs(TB_micrich_BRTE$res)
TB_micrich_BRTE$abs_res2 <- TB_micrich_BRTE$abs_res^2
levene_brte_micrich <- lm(data = TB_micrich_BRTE, abs_res2 ~ BRTErel)
anova(levene_brte_micrich) #p = 0.9663 so > 0.05 so equal variance is met

#autocorrelation
acf(residuals(TB_micrich_BRTEstats, retype = "normalized"))  #lines dont go outside CI horizontal lines, so not autocorrelated
pacf(residuals(TB_micrich_BRTEstats, retype = "normalized"))



#Fort Keogh
#Plant species richness
FK_comm <- community_structure(df = FK2020cover,
                               abundance.var = "cover", 
                               replicate.var = "plot", 
                               metric = "Evar") %>%
  full_join(FK_plotinfo)

#add in BRAR and BRTE cover columns
FK_commBRAR <- FK2020cover %>%
  filter(symbol == "BRAR") %>%
  mutate(BRARcov = cover) %>%
  mutate(BRARrel = rel_cov) %>%
  dplyr::select(c(plot, invasion_percent, BRARcov, BRARrel)) %>%
  full_join(FK_comm) %>%
  drop_na(BRARcov) %>%
  rename(plant_rich = richness, plant_even = Evar)

#statistics
#FK richness
FKBRARrichstats <- lmerTest::lmer(data = FK_commBRAR, plant_rich ~ BRARrel +
                                    (1|grad_num)) 

anova(FKBRARrichstats, type = 3) #p = 0.0001627
rsquared(FKBRARrichstats) #marg = 0.2436837, cond = 0.3859954

#normality check
res_a <- lm(data = FK_commBRAR, plant_rich ~ BRARrel)
ols_plot_resid_hist(res_a) #normalish
ols_test_normality(res_a) #pass


#linearity
plot(resid(FKBRARrichstats), FK_commBRAR$plant_rich) #looks linear
plot(FKBRARrichstats) #no pattern so indicates linearity

#homoscedascity
FK_commBRAR$res <- residuals(FKBRARrichstats)
FK_commBRAR$abs_res <- abs(FK_commBRAR$res)
FK_commBRAR$abs_res2 <- FK_commBRAR$abs_res^2
levene_brar_plrich2 <- lm(data = FK_commBRAR, abs_res2 ~ BRARrel)
anova(levene_brar_plrich2) #p = 0.8788 so > 0.05 so equal variance is met

#autocorrelation
acf(residuals(FKBRARrichstats, retype = "normalized"))  #lines dont go outside CI horizontal lines, so not autocorrelated
pacf(residuals(FKBRARrichstats, retype = "normalized"))


#Insect family richness
FKinsID <- FKinsID2020 %>%
  group_by(location, year, plot, order, family) %>%
  summarise(count = length(sample)) %>%
  ungroup()
FKinsID2 <- FKinsID %>%
  drop_na(family)

FK_inscomm <- community_structure(df = FKinsID2,
                                  abundance.var = "count", 
                                  replicate.var = "plot", 
                                  metric = "Evar") %>%
  full_join(FK_plotinfo) %>%
  drop_na(richness) %>%
  rename(ins_rich = richness, ins_even = Evar)

FK_inscomm_BRAR <- FK_inscomm %>%
  full_join(FK_commBRAR)

#statistics
FKBRARinsrichstats <- lmerTest::lmer(data = FK_inscomm_BRAR, ins_rich ~ BRARrel +
                                       (1|grad_num)) 

anova(FKBRARinsrichstats, type = 3) #p = 0.5961

#normality check - use untransformed data
res_a <- lm(data = FK_inscomm_BRAR, ins_rich ~ BRARrel)
ols_plot_resid_hist(res_a) #normalish
ols_test_normality(res_a) #pass

#linearity
plot(resid(FKBRARinsrichstats), FK_inscomm_BRAR$ins_rich) #looks linear
plot(FKBRARinsrichstats) #no pattern so indicates linearity

#homoscedascity
FK_inscomm_BRAR$res <- residuals(FKBRARinsrichstats)
FK_inscomm_BRAR$abs_res <- abs(FK_inscomm_BRAR$res)
FK_inscomm_BRAR$abs_res2 <- FK_inscomm_BRAR$abs_res^2
levene_brar_insrich2 <- lm(data = FK_inscomm_BRAR, abs_res2 ~ BRARrel)
anova(levene_brar_insrich2) #p = 0.9355 so > 0.05 so equal variance is met

#autocorrelation
acf(residuals(FKBRARinsrichstats, retype = "normalized"))  #lines dont go outside CI horizontal lines, so not autocorrelated
pacf(residuals(FKBRARinsrichstats, retype = "normalized"))



###Microbe ASV richness 
FK_micrich_BRAR <- FK_mic_rich %>%
  full_join(FK_inscomm_BRAR)

#statistics
FK_micrich_BRARstats <- lmerTest::lmer(data = FK_micrich_BRAR, observed_features ~ BRARrel +
                                         (1|grad_num)) 

anova(FK_micrich_BRARstats, type = 3) #p = 0.7286

#normality check - use untransformed data
res_a <- lm(data = FK_micrich_BRAR, observed_features ~ BRARrel)
ols_plot_resid_hist(res_a) #right skew
ols_test_normality(res_a) #fail

FK_micrich_BRAR_trans <- FK_micrich_BRAR %>%
  mutate(ln = log10(observed_features))

res_a <- lm(data = FK_micrich_BRAR_trans, ln ~ BRARrel)
ols_plot_resid_hist(res_a) #right skew
ols_test_normality(res_a) #pass some, better

#redo stats on trans data
FK_micrich_BRARstats2 <- lmerTest::lmer(data = FK_micrich_BRAR_trans, ln ~ BRARrel +
                                          (1|grad_num)) 

anova(FK_micrich_BRARstats2, type = 3) #p = 0.8558

#take outlier out of plot 41
FK_micrich_BRAR_out <- FK_mic_rich %>%
  full_join(FK_inscomm_BRAR) %>%
  filter(plot != "41")

#statistics
FK_micrich_BRARstats2 <- lmerTest::lmer(data = FK_micrich_BRAR_out, observed_features ~ BRARrel +
                                          (1|grad_num)) 

anova(FK_micrich_BRARstats2, type = 3) #p = 0.6659
AICc(FK_micrich_BRARstats2)  #590.9003

#normality check - use untransformed data
res_a <- lm(data = FK_micrich_BRAR_out, observed_features ~ BRARrel)
ols_plot_resid_hist(res_a) #normal
ols_test_normality(res_a) #pass


#linearity
FK_micrich_BRAR_out <- FK_micrich_BRAR_out %>%
  drop_na(observed_features)
plot(resid(FK_micrich_BRARstats2), FK_micrich_BRAR_out$observed_features) #looks linear
plot(FK_micrich_BRARstats2) #no pattern so indicates linearity

#homoscedascity
FK_micrich_BRAR_out$res <- residuals(FK_micrich_BRARstats2)
FK_micrich_BRAR_out$abs_res <- abs(FK_micrich_BRAR_out$res)
FK_micrich_BRAR_out$abs_res2 <- FK_micrich_BRAR_out$abs_res^2
levene_brar_micrich2 <- lm(data = FK_micrich_BRAR_out, abs_res2 ~ BRARrel)
anova(levene_brar_micrich2) #p = 0.7025 so > 0.05 so equal variance is met

#autocorrelation
acf(residuals(FK_micrich_BRARstats2, retype = "normalized"))  #lines dont go outside CI horizontal lines, so not autocorrelated
pacf(residuals(FK_micrich_BRARstats2, retype = "normalized"))



#############################################################################################################


#########################################################################
############ Combined FK and TB richness figure #############
#Plant
title1 <- expression(paste("MT ", italic("B. arvensis")))

FK_BRAR_rich <- ggplot(data = FK_commBRAR, aes(x = BRARrel, y = plant_rich)) + 
  geom_point(size = 3) + 
  xlim(0, 80) +
  ylim(0, 20) +
  geom_smooth(method = "lm", se = FALSE) + 
  labs(x = "", y = "Plant species richness\n(# species/plot)", title = title1) +
  theme(plot.title = element_text(vjust = - 7, hjust = 0.07, size = 16))

#Insect
FK_BRAR_insrich <- ggplot(data = FK_inscomm_BRAR, aes(x = BRARrel, y = ins_rich)) + 
  geom_point(size = 3) +
  xlim(0, 80) +
  ylim(0, 20) +
  labs(x = "", y = "Insect family richness\n(# families/plot)", title = title1) +
  theme(plot.title = element_text(vjust = - 7, hjust = 0.07, size = 16))

#Microbe - no outlier
FK_micrich_BRAR3 <- ggplot(data = FK_micrich_BRAR_out, aes(x = BRARrel, y = observed_features)) + 
  geom_point(size = 3) + 
  xlim(0, 80) +
  ylim(500, 2000) +
  labs(x = "Cover (%)", y = "Microbe OTU richness\n(# OTUs/plot)", title = title1) +
  theme(plot.title = element_text(vjust = - 7, hjust = 0.07, size = 16))

#Plant
title2 <- expression(paste("WY ", italic("B. arvensis")))

TB_BRAR_rich <- ggplot(data = TB_commBRAR, aes(x = BRARrel, y = plant_rich)) + 
  geom_point(size = 3) + 
  xlim(0, 80) +
  ylim(0, 20) +
  geom_smooth(method = "lm", se = FALSE) + 
  labs(x = "", y = "", title = title2) +
  theme(plot.title = element_text(vjust = - 7, hjust = 0.07, size = 16))

title3 <- expression(paste("WY ", italic("B. tectorum")))

TB_BRTE_rich <- ggplot(data = TB_commBRTE, aes(x = BRTErel, y = plant_rich)) + 
  geom_point(size = 3) + 
  xlim(0, 80) +
  ylim(0, 20) +
  geom_smooth(method = "lm", se = FALSE) + 
  labs(x = "", y = "", title = title3) +
  theme(plot.title = element_text(vjust = - 7, hjust = 0.07, size = 16))

#Insect
TB_BRAR_insrich <- ggplot(data = TB_inscomm_BRAR, aes(x = BRARrel, y = ins_rich)) + 
  geom_point(size = 3) +
  xlim(0, 80) +
  ylim(0, 20) +
  labs(x = "", y = "", title = title2) +
  theme(plot.title = element_text(vjust = - 7, hjust = 0.07, size = 16))

TB_BRTE_insrich <- ggplot(data = TB_inscomm_BRTE, aes(x = BRTErel, y = ins_rich)) + 
  geom_point(size = 3) + 
  xlim(0, 80) +
  ylim(0, 20) +
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed") + 
  labs(x = "", y = "", title = title3) +
  theme(plot.title = element_text(vjust = - 7, hjust = 0.07, size = 16))

#Microbe
TB_micrich_BRAR2 <- ggplot(data = TB_micrich_BRAR, aes(x = BRARrel, y = observed_features)) + 
  geom_point(size = 3) + 
  xlim(0, 80) +
  ylim(500, 2000) +
  labs(x = "Cover (%)", y = "", title = title2) +
  theme(plot.title = element_text(vjust = - 7, hjust = 0.07, size = 16))

TB_micrich_BRTE2 <- ggplot(data = TB_micrich_BRTE, aes(x = BRTErel, y = observed_features)) + 
  geom_point(size = 3) + 
  xlim(0, 80) +
  ylim(500, 2000) +
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed") + 
  labs(x = "Cover (%)", y = "", title = title3) +
  theme(plot.title = element_text(vjust = - 7, hjust = 0.07, size = 16))


FK_BRAR_rich + TB_BRAR_rich + TB_BRTE_rich + FK_BRAR_insrich + TB_BRAR_insrich + TB_BRTE_insrich + FK_micrich_BRAR3 + TB_micrich_BRAR2 + TB_micrich_BRTE2 + plot_layout(ncol = 3)

#rich 1400 x 1200




###############################################################################



#####################################################################################
############ Plant, insect, soil microbes species ordination ############
#Thunder Basin
#Plant community NMDS - non brome
#BRAR
BRARTBspcomp2019_2_join <- TB_plotinfo %>%
  full_join(TB_BRAR_2019) %>%
  full_join(TBspcomp2019_2) %>%
  dplyr::select(-c(plot_type, plot_name, phone_lat, phone_long, garmin_lat, garmin_long, aerial_basal, BRAR)) #no BRAR

BRARspcomp <- BRARTBspcomp2019_2_join %>%
  filter(invasive_type == "BRAR")

BRAR_spcompdist <- BRARspcomp %>%
  dplyr::select(-c(1:12))
BRAR_spcomp_meta <- BRARspcomp %>%
  dplyr::select(c(1:12))

BRARspcompNMDS <- metaMDS(BRAR_spcompdist, autotransform = FALSE, shrink = FALSE, distance = "bray")
BRARspcompNMDS #stress = 0.081 fair

BRARspcompscores <- data.frame(scores(BRARspcompNMDS, display = "sites"))
BRARspcompscores2 <- cbind(BRAR_spcomp_meta, BRARspcompscores)
TB_BRARspcomptoplot <- BRARspcompscores2

#BRTE
BRTETBspcomp2019_2_join <- TB_plotinfo %>%
  full_join(TB_BRTE_2019) %>%
  full_join(TBspcomp2019_2) %>%
  dplyr::select(-c(plot_type, plot_name, phone_lat, phone_long, garmin_lat, garmin_long, aerial_basal, BRTE)) #no BRTE

BRTEspcomp <- BRTETBspcomp2019_2_join %>%
  filter(invasive_type == "BRTE")

BRTE_spcompdist <- BRTEspcomp %>%
  dplyr::select(-c(1:12))
BRTE_spcomp_meta <- BRTEspcomp %>%
  dplyr::select(c(1:12))

BRTEspcompNMDS <- metaMDS(BRTE_spcompdist, autotransform = FALSE, shrink = FALSE, distance = "bray")
BRTEspcompNMDS #stress = 0.14 fair

BRTEspcompscores <- data.frame(scores(BRTEspcompNMDS, display = "sites"))
BRTEspcompscores2 <- cbind(BRTE_spcomp_meta, BRTEspcompscores)
TB_BRTEspcomptoplot <- BRTEspcompscores2

#stats
TBBRARcompstats <- adonis2(formula = BRAR_spcompdist ~ avg_rel + (1|grad_num), data = BRAR_spcomp_meta, method = "bray")
TBBRARcompstats #inv% p = 0.031

RVAideMemoire::pairwise.perm.manova(dist(BRAR_spcompdist, "euclidean"), BRAR_spcomp_meta$invasion_percent, nperm = 999, p.method = "BH")
#none

TBBRTEcompstats <- adonis2(formula = BRTE_spcompdist ~ avg_rel + (1|grad_num), data = BRTE_spcomp_meta, method = "bray")
TBBRTEcompstats #inv% p = 0.001

RVAideMemoire::pairwise.perm.manova(dist(BRTE_spcompdist, "euclidean"), BRTE_spcomp_meta$invasion_percent, nperm = 999, p.method = "BH")



#Insect community NMDS
TBinsID3 <- TBinsID2 %>%
  mutate(order_fam = paste(order, family, sep = "_")) %>%
  dplyr::select(-c(order, family)) %>%
  spread(key = order_fam, value = count)
TBinsID3[is.na(TBinsID3)] <- 0

#BRAR
BRARTBins <- TB_BRAR_2019 %>%
  full_join(TB_plotinfo) %>%
  full_join(TBinsID3) %>%
  filter(invasive_type == "BRAR") %>%
  filter(plot != 8)

BRAR_insdist <- BRARTBins %>%
  dplyr::select(-c(1:18))
BRAR_ins_meta <- BRARTBins %>%
  dplyr::select(c(1:18))

BRARinsNMDS <- metaMDS(BRAR_insdist, autotransform = FALSE, shrink = FALSE, distance = "bray")
BRARinsNMDS #stress = 0.098 fair

BRARinsscores <- data.frame(scores(BRARinsNMDS, display = "sites"))
BRARinsscores2 <- cbind(BRAR_ins_meta, BRARinsscores)
TB_BRARinstoplot <- BRARinsscores2

#BRTE
BRTETBins <- TB_BRTE_2019 %>%
  full_join(TB_plotinfo) %>%
  full_join(TBinsID3) %>%
  filter(invasive_type == "BRTE") %>%
  filter(plot != 8)

BRTE_insdist <- BRTETBins %>%
  dplyr::select(-c(1:18))
BRTE_ins_meta <- BRTETBins %>%
  dplyr::select(c(1:18))

BRTEinsNMDS <- metaMDS(BRTE_insdist, autotransform = FALSE, shrink = FALSE, distance = "bray")
BRTEinsNMDS #stress = 0.097 fair

BRTEinsscores <- data.frame(scores(BRTEinsNMDS, display = "sites"))
BRTEinsscores2 <- cbind(BRTE_ins_meta, BRTEinsscores)
TB_BRTEinstoplot <- BRTEinsscores2

#stats
TBBRARinsstats <- adonis2(formula = BRAR_insdist ~ avg_rel + (1|grad_num), data = BRAR_ins_meta, method = "bray")
TBBRARinsstats #inv% p = 0.026

RVAideMemoire::pairwise.perm.manova(dist(BRAR_insdist, "euclidean"), BRAR_ins_meta$invasion_percent, nperm = 999, p.method = "BH")
#none

TBBRTEinsstats <- adonis2(formula = BRTE_insdist ~ avg_rel + (1|grad_num), data = BRTE_ins_meta, method = "bray")
TBBRTEinsstats #inv% p = 0.009

RVAideMemoire::pairwise.perm.manova(dist(BRTE_insdist, "euclidean"), BRTE_ins_meta$invasion_percent, nperm = 999, p.method = "BH")
#0-50, 0-75, 0-100

#Microbe OTU NMDS
TBmicord2 <- TB_mic_rich %>%
  full_join(TBmicord)
#BRAR
TBmicord2_BRAR <- TB_BRAR_2019 %>%
  full_join(TBmicord2) %>%
  filter(invasive_type == "BRAR")

TBmicord2_BRAR_dist <- TBmicord2_BRAR %>%
  dplyr::select(-c(1:15))
TBmicord2_BRAR_meta <- TBmicord2_BRAR %>%
  dplyr::select(c(1:15))

TBmicBRARNMDS <- metaMDS(TBmicord2_BRAR_dist, autotransform = FALSE, shrink = FALSE, distance = "bray")
TBmicBRARNMDS #stress = 0.0483951 good

TBBRARmicscores <- data.frame(scores(TBmicBRARNMDS, display = "sites"))
TBBRARmicscores2 <- cbind(TBmicord2_BRAR_meta, TBBRARmicscores)
BRARtoplot <- TBBRARmicscores2

#BRTE
TBmicord2_BRTE <- TB_BRTE_2019 %>%
  full_join(TBmicord2) %>%
  filter(invasive_type == "BRTE")

TBmicord2_BRTE_dist <- TBmicord2_BRTE %>%
  dplyr::select(-c(1:15))
TBmicord2_BRTE_meta <- TBmicord2_BRTE %>%
  dplyr::select(c(1:15))

TBmicBRTENMDS <- metaMDS(TBmicord2_BRTE_dist, autotransform = FALSE, shrink = FALSE, distance = "bray")
TBmicBRTENMDS #stress = 0.1063307 fair

TBBRTEmicscores <- data.frame(scores(TBmicBRTENMDS, display = "sites"))
TBBRTEmicscores2 <- cbind(TBmicord2_BRTE_meta, TBBRTEmicscores)
BRTEtoplot <- TBBRTEmicscores2

#stats
TBBRARmicstats <- adonis2(formula = TBmicord2_BRAR_dist ~ avg_rel + (1|grad_num), data = TBmicord2_BRAR_meta, method = "bray")
TBBRARmicstats #inv% p = 0.425

TBBRTEmicstats <- adonis2(formula = TBmicord2_BRTE_dist ~ avg_rel + (1|grad_num), data = TBmicord2_BRTE_meta, method = "bray")
TBBRTEmicstats #inv% p = 0.001

RVAideMemoire::pairwise.perm.manova(dist(TBmicord2_BRTE_dist, "euclidean"), TBmicord2_BRTE_meta$invasion_percent, nperm = 999, p.method = "BH")
#0-25, 0-50, 0-75, 0-100, 25-75, 25-100






#Fort Keogh
#Plant community NMDS - non brome
#BRAR
BRARFKspcomp2020_2_join <- FK_plotinfo %>%
  full_join(FK_BRAR) %>%
  full_join(FKspcomp2020_2) %>%
  filter(year == "2020") %>%
  dplyr::select(-c(plot_type, latitude, longitude, plot_name, aerial_basal, BRAR)) #no BRAR

BRARFKspcomp2020_2_join[is.na(BRARFKspcomp2020_2_join)] <- 0

BRAR_spcompdist2 <- BRARFKspcomp2020_2_join %>%
  dplyr::select(-c(1:12))
BRAR_spcomp_meta2 <- BRARFKspcomp2020_2_join %>%
  dplyr::select(c(1:12))

BRARspcompNMDS2 <- metaMDS(BRAR_spcompdist2, autotransform = FALSE, shrink = FALSE, distance = "bray")
BRARspcompNMDS2 #stress = 0.216 poor

BRARspcompscoresFK <- data.frame(scores(BRARspcompNMDS2, display = "sites"))
BRARspcompscores2FK <- cbind(BRAR_spcomp_meta2, BRARspcompscoresFK)
FK_BRARspcomptoplot <- BRARspcompscores2FK

#stats
FKBRARcompstats <- adonis2(formula = BRAR_spcompdist2 ~ avg_rel + (1|grad_num), data = BRAR_spcomp_meta2, method = "bray")
FKBRARcompstats #inv% p = 0.001

RVAideMemoire::pairwise.perm.manova(dist(BRAR_spcompdist2, "euclidean"), BRAR_spcomp_meta2$invasion_percent, nperm = 999, p.method = "BH")
#0-50, 0-100


#Insect community NMDS
FKinsID3 <- FKinsID2 %>%
  mutate(order_fam = paste(order, family, sep = "_")) %>%
  dplyr::select(-c(order, family)) %>%
  spread(key = order_fam, value = count)
FKinsID3[is.na(FKinsID3)] <- 0

#BRAR
BRARFKins <- FK_BRAR %>%
  full_join(FK_plotinfo) %>%
  full_join(FKinsID3) %>%
  filter(year == "2020") 

BRAR_insdist2 <- BRARFKins %>%
  dplyr::select(-c(1:16))
BRAR_ins_meta2 <- BRARFKins %>%
  dplyr::select(c(1:16))

BRARinsNMDS2 <- metaMDS(BRAR_insdist2, autotransform = FALSE, shrink = FALSE, distance = "bray")
BRARinsNMDS2 #stress = 0.15 poor

BRARinsscoresFK <- data.frame(scores(BRARinsNMDS2, display = "sites"))
BRARinsscores2FK <- cbind(BRAR_ins_meta2, BRARinsscoresFK)
FK_BRARinstoplot <- BRARinsscores2FK

#stats
FKBRARinsstats <- adonis2(formula = BRAR_insdist2 ~ avg_rel + (1|grad_num), data = BRAR_ins_meta2, method = "bray")
FKBRARinsstats #inv% p = 0.005

RVAideMemoire::pairwise.perm.manova(dist(BRAR_insdist2, "euclidean"), BRAR_ins_meta2$invasion_percent, nperm = 999, p.method = "BH")
#none

#Microbe OTU NMDS
FKmicord2 <- FK_mic_rich %>%
  full_join(FKmicord)
#BRAR
FKmicord2_BRAR <- FK_BRAR %>%
  filter(year == "2020") %>%
  full_join(FKmicord2) 

FKmicord2_BRAR_dist <- FKmicord2_BRAR %>%
  dplyr::select(-c(1:15))
FKmicord2_BRAR_meta <- FKmicord2_BRAR %>%
  dplyr::select(c(1:15))

FKmicBRARNMDS <- metaMDS(FKmicord2_BRAR_dist, autotransform = FALSE, shrink = FALSE, distance = "bray")
FKmicBRARNMDS #stress = 0.113 fair

FKBRARmicscores <- data.frame(scores(FKmicBRARNMDS, display = "sites"))
FKBRARmicscores2 <- cbind(FKmicord2_BRAR_meta, FKBRARmicscores)
FKBRARtoplot <- FKBRARmicscores2

#stats
FKBRARmicstats <- adonis2(formula = FKmicord2_BRAR_dist ~ avg_rel + (1|grad_num), data = FKmicord2_BRAR_meta, method = "bray")
FKBRARmicstats #inv% p = 0.151



#######################################################################


#######################################################################
############ Combined FK and TB species NMDS #######
cbPalette_nmds <- c("0" = "#999999", "25" = "#E69F00", "50" = "#56B4E9", "75" = "#009E73", "100" = "#F0E442")
cbPalette_nmds2 <- c("0" = "navy", "25" = "royalblue1", "50" = "cadetblue2", "75" = "aquamarine3", "100" = "seagreen4")

#Plant
FKBRARcompNMDS1 <- ggplot(FK_BRARspcomptoplot, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(NMDS1, NMDS2, colour = factor(invasion_percent), size = invasion_percent)) +
  xlim(-1, 1) + 
  ylim(-0.8, 0.8) + 
  labs(x = "", colour = "Invasion Percent", size = "Invasion Percent", title = title1) +
  theme(plot.title = element_text(vjust = - 7, hjust = 0.07, size = 16), legend.position = "none") + 
  scale_colour_manual(values = cbPalette_nmds2) + 
  scale_size_continuous(range = c(4, 8))

#Insect
FKBRARinsNMDS1 <- ggplot(FK_BRARinstoplot, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(NMDS1, NMDS2, colour = factor(invasion_percent), size = invasion_percent)) +
  xlim(-1, 1.3) + 
  ylim(-0.5, 1) +
  labs(x = "", colour = "Invasion Percent", size = "Invasion Percent", title = title1) +
  theme(plot.title = element_text(vjust = - 7, hjust = 0.07, size = 16), legend.position = "none") + 
  scale_colour_manual(values = cbPalette_nmds2) + 
  scale_size_continuous(range = c(4, 8))

#Microbe
FKBRARmicNMDS1 <- ggplot(FKBRARtoplot, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(NMDS1, NMDS2, colour = factor(invasion_percent), size = invasion_percent)) +
  xlim(-0.05, 0.1) + 
  ylim(-0.05, 0.1) +
  labs(colour = "Invasion Percent", size = "Invasion Percent", title = title1) +
  theme(plot.title = element_text(vjust = - 7, hjust = 0.07, size = 16), legend.position = "none") + 
  scale_colour_manual(values = cbPalette_nmds2) + 
  scale_size_continuous(range = c(4, 8))


#Plant
BRARcompNMDS1 <- ggplot(TB_BRARspcomptoplot, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(NMDS1, NMDS2, colour = factor(invasion_percent), size = invasion_percent)) +
  xlim(-0.5, 1.8) + 
  ylim(-0.8, 0.8) + 
  labs(y = "", x = "", colour = "Invasion Percent", size = "Invasion Percent", title = title2) +
  theme(plot.title = element_text(vjust = - 7, hjust = 0.07, size = 16), legend.position = "none") + 
  scale_colour_manual(values = cbPalette_nmds2) + 
  scale_size_continuous(range = c(4, 8))

BRTEcompNMDS1 <- ggplot(TB_BRTEspcomptoplot, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(NMDS1, NMDS2, colour = factor(invasion_percent), size = invasion_percent)) +
  xlim(-0.6, 1.8) + 
  ylim(-1.5, 1.5) +
  labs(x = "", y = "", colour = "Invasion Percent", size = "Invasion Percent", title = title3) +
  theme(plot.title = element_text(vjust = - 7, hjust = 0.07, size = 16), legend.position = "none") + 
  scale_colour_manual(values = cbPalette_nmds2) + 
  scale_size_continuous(range = c(4, 8))

#Insect
BRARinsNMDS1 <- ggplot(TB_BRARinstoplot, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(NMDS1, NMDS2, colour = factor(invasion_percent), size = invasion_percent)) +
  xlim(-2.5, 2) + 
  ylim(-1.5, 2.5) +
  labs(y = "", x = "", colour = "Invasion Percent", size = "Invasion Percent", title = title2) +
  theme(plot.title = element_text(vjust = - 7, hjust = 0.07, size = 16), legend.position = "none") + 
  scale_colour_manual(values = cbPalette_nmds2) + 
  scale_size_continuous(range = c(4, 8))

BRTEinsNMDS1 <- ggplot(TB_BRTEinstoplot, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(NMDS1, NMDS2, colour = factor(invasion_percent), size = invasion_percent)) +
  xlim(-1.5, 0.9) + 
  ylim(-1, 2) +
  labs(x = "", y = "", colour = "Invasion Percent", size = "Invasion Percent", title = title3) +
  theme(plot.title = element_text(vjust = - 7, hjust = 0.07, size = 16), legend.position = "none") + 
  scale_colour_manual(values = cbPalette_nmds2) + 
  scale_size_continuous(range = c(4, 8))

#Microbe
BRARmicNMDS1 <- ggplot(BRARtoplot, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(NMDS1, NMDS2, colour = factor(invasion_percent), size = invasion_percent)) +
  xlim(-0.2, 0.1) + 
  ylim(-0.05, 0.08) +
  labs(y = "", colour = "Invasion Percent", size = "Invasion Percent", title = title2) +
  theme(plot.title = element_text(vjust = - 7, hjust = 0.07, size = 16), legend.position = "none") + 
  scale_colour_manual(values = cbPalette_nmds2) + 
  scale_size_continuous(range = c(4, 8))

BRTEmicNMDS1 <- ggplot(BRTEtoplot, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(NMDS1, NMDS2, colour = factor(invasion_percent), size = invasion_percent)) +
  xlim(-0.1, 0.05) + 
  ylim(-0.06, 0.08) + 
  labs(y = "", colour = "Invasion Percent", size = "Invasion Percent", title = title3) +
  theme(plot.title = element_text(vjust = - 7, hjust = 0.07, size = 16)) + 
  scale_colour_manual(values = cbPalette_nmds2) + 
  scale_size_continuous(range = c(4, 8))


FKBRARcompNMDS1 + BRARcompNMDS1 + BRTEcompNMDS1 + FKBRARinsNMDS1 + BRARinsNMDS1 + BRTEinsNMDS1 + FKBRARmicNMDS1 +BRARmicNMDS1 + BRTEmicNMDS1 + plot_layout(ncol = 3)

#NMDScom 1900x1400


###############################################################################




######################################################################################################
############ Plant, insect, soil microbes functional ordination ############
#Thunbder Basin
#Plants

TB2019cover_sp <- TB2019cover %>%
  full_join(TBspecies)
TB2019cover_sp <- TB2019cover_sp[c(1:5, 7:14, 15, 16, 17, 6)]

#BRAR
TB2019cover_sp_BRAR <- TB2019cover_sp %>%
  filter(invasive_type == "BRAR") %>%
  filter(symbol != "BRAR") %>%
  dplyr::select(-c(location, symbol, totcov, rel_cov, scientific_name, functional_group, perennial_annual, native_invasive, family, spec_funct_group)) %>%
  group_by(plot, invasion_percent, grad_num, funct2) %>%
  summarise(avg_cov = mean(cover)) %>%
  spread(key = funct2, value = avg_cov) %>%
  ungroup() %>%
  full_join(TB_BRAR_2019)

TB2019cover_sp_BRAR <- TB2019cover_sp_BRAR[c(1:3, 10:16, 4:9)]

TB_BRAR_NMDS_pl_func_dist <- TB2019cover_sp_BRAR %>%
  dplyr::select(c(11:16))
TB_BRAR_NMDS_pl_func_meta <- TB2019cover_sp_BRAR %>%
  dplyr::select(-c(11:16))

TB_BRAR_NMDS_pl_func <- metaMDS(TB_BRAR_NMDS_pl_func_dist, autotransform = FALSE, shrink = FALSE, distance = "bray")
TB_BRAR_NMDS_pl_func #stress = 0.06 fair

TB_BRAR_NMDS_pl_func_scores <- data.frame(scores(TB_BRAR_NMDS_pl_func, display = "sites"))
TB_BRAR_NMDS_pl_func_scores2 <- cbind(TB_BRAR_NMDS_pl_func_meta, TB_BRAR_NMDS_pl_func_scores)
TB_BRAR_NMDS_pl_func_plot <- TB_BRAR_NMDS_pl_func_scores2

#BRTE
TB2019cover_sp_BRTE <- TB2019cover_sp %>%
  filter(invasive_type == "BRTE") %>%
  filter(symbol != "BRTE") %>%
  dplyr::select(-c(location, symbol, totcov, rel_cov, scientific_name, functional_group, perennial_annual, native_invasive, family, spec_funct_group)) %>%
  group_by(plot, invasion_percent, grad_num, funct2) %>%
  summarise(avg_cov = mean(cover)) %>%
  spread(key = funct2, value = avg_cov) %>%
  ungroup() %>%
  full_join(TB_BRTE_2019)

TB2019cover_sp_BRTE <- TB2019cover_sp_BRTE[c(1:3, 10:16, 4:9)]

TB_BRTE_NMDS_pl_func_dist <- TB2019cover_sp_BRTE %>%
  dplyr::select(c(11:16))
TB_BRTE_NMDS_pl_func_meta <- TB2019cover_sp_BRTE %>%
  dplyr::select(-c(11:16))

TB_BRTE_NMDS_pl_func <- metaMDS(TB_BRTE_NMDS_pl_func_dist, autotransform = FALSE, shrink = FALSE, distance = "bray")
TB_BRTE_NMDS_pl_func #stress = 0.11 fair

TB_BRTE_NMDS_pl_func_scores <- data.frame(scores(TB_BRTE_NMDS_pl_func, display = "sites"))
TB_BRTE_NMDS_pl_func_scores2 <- cbind(TB_BRTE_NMDS_pl_func_meta, TB_BRTE_NMDS_pl_func_scores)
TB_BRTE_NMDS_pl_func_plot <- TB_BRTE_NMDS_pl_func_scores2

#stats
TBBRARplstats2 <- adonis2(formula = TB_BRAR_NMDS_pl_func_dist ~ avg_rel + (1|grad_num), data = TB_BRAR_NMDS_pl_func_meta, method = "bray")
TBBRARplstats2 #inv% p = 0.051

RVAideMemoire::pairwise.perm.manova(dist(TB_BRAR_NMDS_pl_func_dist, "euclidean"), TB_BRAR_NMDS_pl_func_meta$invasion_percent, nperm = 999, p.method = "BH")
#none


TBBRTEplstats2 <- adonis2(formula = TB_BRTE_NMDS_pl_func_dist ~ avg_rel + (1|grad_num), data = TB_BRTE_NMDS_pl_func_meta, method = "bray")
TBBRTEplstats2 #inv% p = 0.001

RVAideMemoire::pairwise.perm.manova(dist(TB_BRTE_NMDS_pl_func_dist, "euclidean"), TB_BRTE_NMDS_pl_func_meta$invasion_percent, nperm = 999, p.method = "BH")
#0-75, 0-100, 25-75, 25-100, 50-75, 50-100



#Insects
TBinsID_guild <- TBguild %>%
  rename(order = order2, family = family2) %>%
  full_join(TBinsID2) %>%
  group_by(plot, guild) %>%
  summarise(avg_n = mean(count)) %>%
  ungroup() %>%
  full_join(TB_plotinfo)

#BRAR
TBinsID_guild_BRAR <- TBinsID_guild %>%
  filter(invasive_type == "BRAR") %>%
  full_join(TB_BRAR_2019) %>%
  dplyr::select(-c(5, 7, 10:13)) %>%
  filter(plot != 8)

TBinsID_guild_BRAR <- TBinsID_guild_BRAR[c(1, 4:14, 2, 3)]

TBinsID_guild_BRAR2 <- TBinsID_guild_BRAR %>%
  spread(key = guild, value = avg_n)

TBinsID_guild_BRAR2[is.na(TBinsID_guild_BRAR2)] <- 0

TB_BRAR_NMDS_ins_func_dist <- TBinsID_guild_BRAR2 %>%
  dplyr::select(c(13:19))
TB_BRAR_NMDS_ins_func_meta <- TBinsID_guild_BRAR2 %>%
  dplyr::select(-c(13:19))

TB_BRAR_NMDS_ins_func <- metaMDS(TB_BRAR_NMDS_ins_func_dist, autotransform = FALSE, shrink = FALSE, distance = "bray")
TB_BRAR_NMDS_ins_func #stress = 0.06 good

TB_BRAR_NMDS_ins_func_scores <- data.frame(scores(TB_BRAR_NMDS_ins_func, display = "sites"))
TB_BRAR_NMDS_ins_func_scores2 <- cbind(TB_BRAR_NMDS_ins_func_meta, TB_BRAR_NMDS_ins_func_scores)
TB_BRAR_NMDS_ins_func_plot <- TB_BRAR_NMDS_ins_func_scores2

#BRTE
TBinsID_guild_BRTE <- TBinsID_guild %>%
  filter(invasive_type == "BRTE") %>%
  full_join(TB_BRTE_2019) %>%
  dplyr::select(-c(5, 7, 10:13)) 

TBinsID_guild_BRTE <- TBinsID_guild_BRTE[c(1, 4:14, 2, 3)]

TBinsID_guild_BRTE2 <- TBinsID_guild_BRTE %>%
  spread(key = guild, value = avg_n)

TBinsID_guild_BRTE2[is.na(TBinsID_guild_BRTE2)] <- 0

TB_BRTE_NMDS_ins_func_dist <- TBinsID_guild_BRTE2 %>%
  dplyr::select(c(13:18))
TB_BRTE_NMDS_ins_func_meta <- TBinsID_guild_BRTE2 %>%
  dplyr::select(-c(13:18))

TB_BRTE_NMDS_ins_func <- metaMDS(TB_BRTE_NMDS_ins_func_dist, autotransform = FALSE, shrink = FALSE, distance = "bray")
TB_BRTE_NMDS_ins_func #stress = 0.08 good

TB_BRTE_NMDS_ins_func_scores <- data.frame(scores(TB_BRTE_NMDS_ins_func, display = "sites"))
TB_BRTE_NMDS_ins_func_scores2 <- cbind(TB_BRTE_NMDS_ins_func_meta, TB_BRTE_NMDS_ins_func_scores)
TB_BRTE_NMDS_ins_func_plot <- TB_BRTE_NMDS_ins_func_scores2

#stats
TBBRARinsstats2 <- adonis2(formula = TB_BRAR_NMDS_ins_func_dist ~ avg_rel + (1|grad_num), data = TB_BRAR_NMDS_ins_func_meta, method = "bray")
TBBRARinsstats2 #inv% p = 0.177

TBBRTEinsstats2 <- adonis2(formula = TB_BRTE_NMDS_ins_func_dist ~ avg_rel + (1|grad_num), data = TB_BRTE_NMDS_ins_func_meta, method = "bray")
TBBRTEinsstats2 #inv% p = 0.033

RVAideMemoire::pairwise.perm.manova(dist(TB_BRTE_NMDS_ins_func_dist, "euclidean"), TB_BRTE_NMDS_ins_func_meta$invasion_percent, nperm = 999, p.method = "BH")
#0-50, 0-75, 0-100

#Microbes
mic_func_long <- mic_func %>%
  rename(func = OUT_ID) %>%
  group_by(func) %>%
  gather(index, count, 2:ncol(mic_func)) %>%
  ungroup() %>%
  full_join(mapping) 

mic_func_long <- mic_func_long[c(2, 4:10, 1, 3)] #reorder columns

TBmic_func_long <- mic_func_long %>%
  filter(location == "TB") %>%
  filter(count != 0)

TBmic_func_wide <- TBmic_func_long %>% 
  spread(key = func, value = count)
TBmic_func_wide[is.na(TBmic_func_wide)] <- 0

#BRAR
TBmic_func_wide_BRAR <- TBmic_func_wide %>%
  filter(invasive_type == "BRAR")

TB_BRAR_NMDS_mic_func_dist <- TBmic_func_wide_BRAR %>%
  dplyr::select(c(9:ncol(TBmic_func_wide_BRAR)))
TB_BRAR_NMDS_mic_func_meta <- TBmic_func_wide_BRAR %>%
  dplyr::select(-c(9:ncol(TBmic_func_wide_BRAR))) %>%
  full_join(TB_BRTE_2019)

TB_BRAR_NMDS_mic_func <- metaMDS(TB_BRAR_NMDS_mic_func_dist, autotransform = FALSE, shrink = FALSE, distance = "bray")
TB_BRAR_NMDS_mic_func #stress = 0.06 good

TB_BRAR_NMDS_mic_func_scores <- data.frame(scores(TB_BRAR_NMDS_mic_func, display = "sites"))
TB_BRAR_NMDS_mic_func_scores2 <- cbind(TB_BRAR_NMDS_mic_func_meta, TB_BRAR_NMDS_mic_func_scores)
TB_BRAR_NMDS_mic_func_plot <- TB_BRAR_NMDS_mic_func_scores2 %>%
  full_join(TB_BRAR_2019)

#BRTE
TBmic_func_wide_BRTE <- TBmic_func_wide %>%
  filter(invasive_type == "BRTE")

TB_BRTE_NMDS_mic_func_dist <- TBmic_func_wide_BRTE %>%
  dplyr::select(c(9:ncol(TBmic_func_wide_BRTE)))
TB_BRTE_NMDS_mic_func_meta <- TBmic_func_wide_BRTE %>%
  dplyr::select(-c(9:ncol(TBmic_func_wide_BRTE))) %>%
  full_join(TB_BRTE_2019)

TB_BRTE_NMDS_mic_func <- metaMDS(TB_BRTE_NMDS_mic_func_dist, autotransform = FALSE, shrink = FALSE, distance = "bray")
TB_BRTE_NMDS_mic_func #stress = 0.06 good

TB_BRTE_NMDS_mic_func_scores <- data.frame(scores(TB_BRTE_NMDS_mic_func, display = "sites"))
TB_BRTE_NMDS_mic_func_scores2 <- cbind(TB_BRTE_NMDS_mic_func_meta, TB_BRTE_NMDS_mic_func_scores)
TB_BRTE_NMDS_mic_func_plot <- TB_BRTE_NMDS_mic_func_scores2 %>%
  full_join(TB_BRTE_2019)

#stats
TBBRARmicstats2 <- adonis2(formula = TB_BRAR_NMDS_mic_func_dist ~ avg_rel + (1|grad_num), data = TB_BRAR_NMDS_mic_func_meta, method = "bray")
TBBRARmicstats2 #inv% p = 0.125

TBBRTEmicstats2 <- adonis2(formula = TB_BRTE_NMDS_mic_func_dist ~ avg_rel + (1|grad_num), data = TB_BRTE_NMDS_mic_func_meta, method = "bray")
TBBRTEmicstats2 #inv% p = 0.016

RVAideMemoire::pairwise.perm.manova(dist(TB_BRTE_NMDS_mic_func_dist, "euclidean"), TB_BRTE_NMDS_mic_func_meta$invasion_percent, nperm = 999, p.method = "BH")
#none





#Fort Keogh
#Plants

FK2020cover_sp <- FK2020cover %>%
  full_join(FKspecies)
FK2020cover_sp <- FK2020cover_sp[c(1:5, 7:15, 16, 6)]

FK_BRAR_2020 <- FK_BRAR %>%
  filter(year == "2020")

#BRAR
FK2020cover_sp_BRAR <- FK2020cover_sp %>%
  filter(symbol != "BRAR") %>%
  dplyr::select(-c(location, symbol, totcov, rel_cov, scientific_name, functional_group, perennial_annual, native_invasive, spec_funct_group)) %>%
  group_by(plot, invasion_percent, grad_num, funct2) %>%
  summarise(avg_cov = mean(cover)) %>%
  spread(key = funct2, value = avg_cov) %>%
  ungroup() %>%
  full_join(FK_BRAR_2020) %>%
  drop_na(plot)

FK2020cover_sp_BRAR <- FK2020cover_sp_BRAR[c(1:3, 10:16, 4:9)]

FK_BRAR_NMDS_pl_func_dist <- FK2020cover_sp_BRAR %>%
  dplyr::select(c(11:16))
FK_BRAR_NMDS_pl_func_meta <- FK2020cover_sp_BRAR %>%
  dplyr::select(-c(11:16))

FK_BRAR_NMDS_pl_func <- metaMDS(FK_BRAR_NMDS_pl_func_dist, autotransform = FALSE, shrink = FALSE, distance = "bray")
FK_BRAR_NMDS_pl_func #stress = 0.15 ok

FK_BRAR_NMDS_pl_func_scores <- data.frame(scores(FK_BRAR_NMDS_pl_func, display = "sites"))
FK_BRAR_NMDS_pl_func_scores2 <- cbind(FK_BRAR_NMDS_pl_func_meta, FK_BRAR_NMDS_pl_func_scores)
FK_BRAR_NMDS_pl_func_plot <- FK_BRAR_NMDS_pl_func_scores2

#stats
FKBRARplstats2 <- adonis2(formula = FK_BRAR_NMDS_pl_func_dist ~ avg_rel + (1|grad_num), data = FK_BRAR_NMDS_pl_func_meta, method = "bray")
FKBRARplstats2 #inv% p = 0.001

RVAideMemoire::pairwise.perm.manova(dist(FK_BRAR_NMDS_pl_func_dist, "euclidean"), FK_BRAR_NMDS_pl_func_meta$invasion_percent, nperm = 999, p.method = "BH")
#0-50, 0-100


#Insects
FKinsID_guild <- FKguild %>%
  rename(order = order2, family = family2) %>%
  full_join(FKinsID2) %>%
  group_by(plot, guild) %>%
  summarise(avg_n = mean(count)) %>%
  ungroup() %>%
  full_join(FK_plotinfo)

#BRAR
FKinsID_guild_BRAR <- FKinsID_guild %>%
  full_join(FK_BRAR_2020) %>%
  dplyr::select(-c(6, 8, 10, 11)) 

FKinsID_guild_BRAR <- FKinsID_guild_BRAR[c(1, 4:14, 2, 3)]

FKinsID_guild_BRAR2 <- FKinsID_guild_BRAR %>%
  spread(key = guild, value = avg_n)

FKinsID_guild_BRAR2[is.na(FKinsID_guild_BRAR2)] <- 0

FK_BRAR_NMDS_ins_func_dist <- FKinsID_guild_BRAR2 %>%
  dplyr::select(c(13:19))
FK_BRAR_NMDS_ins_func_meta <- FKinsID_guild_BRAR2 %>%
  dplyr::select(-c(13:19))

FK_BRAR_NMDS_ins_func <- metaMDS(FK_BRAR_NMDS_ins_func_dist, autotransform = FALSE, shrink = FALSE, distance = "bray")
FK_BRAR_NMDS_ins_func #stress = 0.09 good

FK_BRAR_NMDS_ins_func_scores <- data.frame(scores(FK_BRAR_NMDS_ins_func, display = "sites"))
FK_BRAR_NMDS_ins_func_scores2 <- cbind(FK_BRAR_NMDS_ins_func_meta, FK_BRAR_NMDS_ins_func_scores)
FK_BRAR_NMDS_ins_func_plot <- FK_BRAR_NMDS_ins_func_scores2

#stats
FKBRARinsstats2 <- adonis2(formula = FK_BRAR_NMDS_ins_func_dist ~ avg_rel + (1|grad_num), data = FK_BRAR_NMDS_ins_func_meta, method = "bray")
FKBRARinsstats2 #inv% p = 0.043

RVAideMemoire::pairwise.perm.manova(dist(FK_BRAR_NMDS_ins_func_dist, "euclidean"), FK_BRAR_NMDS_ins_func_meta$invasion_percent, nperm = 999, p.method = "BH")
#none

#Microbes
mic_func_long <- mic_func %>%
  rename(func = OUT_ID) %>%
  group_by(func) %>%
  gather(index, count, 2:ncol(mic_func)) %>%
  ungroup() %>%
  full_join(mapping) 

mic_func_long <- mic_func_long[c(2, 4:10, 1, 3)] #reorder columns

FKmic_func_long <- mic_func_long %>%
  filter(location == "FK") %>%
  filter(count != 0)

FKmic_func_wide <- FKmic_func_long %>% 
  spread(key = func, value = count)
FKmic_func_wide[is.na(FKmic_func_wide)] <- 0

#BRAR
FK_BRAR_NMDS_mic_func_dist <- FKmic_func_wide %>%
  dplyr::select(c(9:ncol(FKmic_func_wide)))
FK_BRAR_NMDS_mic_func_meta <- FKmic_func_wide %>%
  dplyr::select(-c(9:ncol(FKmic_func_wide))) %>%
  full_join(FK_BRAR_2020) 

FK_BRAR_NMDS_mic_func <- metaMDS(FK_BRAR_NMDS_mic_func_dist, distance = "bray")
FK_BRAR_NMDS_mic_func #stress = 0.113 fair

FK_BRAR_NMDS_mic_func_scores <- data.frame(scores(FK_BRAR_NMDS_mic_func, display = "sites"))
FK_BRAR_NMDS_mic_func_scores2 <- cbind(FK_BRAR_NMDS_mic_func_meta, FK_BRAR_NMDS_mic_func_scores)
FK_BRAR_NMDS_mic_func_plot <- FK_BRAR_NMDS_mic_func_scores2 %>%
  full_join(FK_BRAR_2020)


#stats
FKBRARmicstats2 <- adonis2(formula = FK_BRAR_NMDS_mic_func_dist ~ avg_rel + (1|grad_num), data = FK_BRAR_NMDS_mic_func_meta, method = "bray")
FKBRARmicstats2 #inv% p = 0.315


###########################################################################################



###############################################################################
############ Combined FK and TB functional ordination #########
#using less plant functional groupings


#Plant
FKBRARcompNMDS2 <- ggplot(FK_BRAR_NMDS_pl_func_plot, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(NMDS1, NMDS2, colour = factor(invasion_percent), size = invasion_percent)) +
  xlim(-1, 1.5) + 
  ylim(-0.8, 0.8) + 
  labs(x = "", colour = "Invasion Percent", size = "Invasion Percent", title = title1) +
  theme(plot.title = element_text(vjust = - 7, hjust = 0.07, size = 16), legend.position = "none") + 
  scale_colour_manual(values = cbPalette_nmds2) + 
  scale_size_continuous(range = c(4, 8))

#Insect
FKBRARinsNMDS2 <- ggplot(FK_BRAR_NMDS_ins_func_plot, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(NMDS1, NMDS2, colour = factor(invasion_percent), size = invasion_percent)) +
  xlim(-1.2, 1.5) + 
  ylim(-1, 1) +
  labs(x = "", colour = "Invasion Percent", size = "Invasion Percent", title = title1) +
  theme(plot.title = element_text(vjust = - 7, hjust = 0.07, size = 16), legend.position = "none") + 
  scale_colour_manual(values = cbPalette_nmds2) + 
  scale_size_continuous(range = c(4, 8))

#Microbe
FKBRARmicNMDS2 <- ggplot(FK_BRAR_NMDS_mic_func_plot, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(NMDS1, NMDS2, colour = factor(invasion_percent), size = invasion_percent)) +
  xlim(-0.5, 0.5) + 
  ylim(-0.3, 0.3) +
  labs(colour = "Invasion Percent", size = "Invasion Percent", title = title1) +
  theme(plot.title = element_text(vjust = - 7, hjust = 0.07, size = 16), legend.position = "none") + 
  scale_colour_manual(values = cbPalette_nmds2) + 
  scale_size_continuous(range = c(4, 8))

#Plant
BRARcompNMDS2 <- ggplot(TB_BRAR_NMDS_pl_func_plot, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(NMDS1, NMDS2, colour = factor(invasion_percent), size = invasion_percent)) +
  xlim(-1, 2.7) + 
  ylim(-0.8, 0.8) + 
  labs(x = "", y = "", colour = "Invasion Percent", size = "Invasion Percent", title = title2) +
  theme(plot.title = element_text(vjust = - 7, hjust = 0.07, size = 16), legend.position = "none") + 
  scale_colour_manual(values = cbPalette_nmds2) + 
  scale_size_continuous(range = c(4, 8))

BRTEcompNMDS2 <- ggplot(TB_BRTE_NMDS_pl_func_plot, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(NMDS1, NMDS2, colour = factor(invasion_percent), size = invasion_percent)) +
  xlim(-1, 2.1) + 
  ylim(-1.2, 1.5) +
  labs(x = "", y = "", colour = "Invasion Percent", size = "Invasion Percent", title = title3) +
  theme(plot.title = element_text(vjust = - 7, hjust = 0.07, size = 16), legend.position = "none") + 
  scale_colour_manual(values = cbPalette_nmds2) + 
  scale_size_continuous(range = c(4, 8))

#Insect
BRARinsNMDS2 <- ggplot(TB_BRAR_NMDS_ins_func_plot, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(NMDS1, NMDS2, colour = factor(invasion_percent), size = invasion_percent)) +
  xlim(-2, 1.5) + 
  ylim(-2, 2) +
  labs(x = "", y = "", colour = "Invasion Percent", size = "Invasion Percent", title = title2) +
  theme(plot.title = element_text(vjust = - 7, hjust = 0.07, size = 16), legend.position = "none") + 
  scale_colour_manual(values = cbPalette_nmds2) + 
  scale_size_continuous(range = c(4, 8))

BRTEinsNMDS2 <- ggplot(TB_BRTE_NMDS_ins_func_plot, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(NMDS1, NMDS2, colour = factor(invasion_percent), size = invasion_percent)) +
  xlim(-1, 2.3) + 
  ylim(-1, 1.1) +
  labs(x = "", y = "", colour = "Invasion Percent", size = "Invasion Percent", title = title3) +
  theme(plot.title = element_text(vjust = - 7, hjust = 0.07, size = 16), legend.position = "none") + 
  scale_colour_manual(values = cbPalette_nmds2) + 
  scale_size_continuous(range = c(4, 8))

#Microbe
BRARmicNMDS2 <- ggplot(TB_BRAR_NMDS_mic_func_plot, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(NMDS1, NMDS2, colour = factor(invasion_percent), size = invasion_percent)) +
  xlim(-0.9, 0.9) + 
  ylim(-0.6, 0.6) +
  labs(y = "", colour = "Invasion Percent", size = "Invasion Percent", title = title2) +
  theme(plot.title = element_text(vjust = - 7, hjust = 0.07, size = 16), legend.position = "none") + 
  scale_colour_manual(values = cbPalette_nmds2) + 
  scale_size_continuous(range = c(4, 8))

BRTEmicNMDS2 <- ggplot(TB_BRTE_NMDS_mic_func_plot, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(NMDS1, NMDS2, colour = factor(invasion_percent), size = invasion_percent)) +
  xlim(-0.6, 0.9) + 
  ylim(-0.3, 0.3) + 
  labs(y = "", colour = "Invasion Percent", size = "Invasion Percent", title = title3) +
  theme(plot.title = element_text(vjust = - 7, hjust = 0.07, size = 16)) + 
  scale_colour_manual(values = cbPalette_nmds2) + 
  scale_size_continuous(range = c(4, 8))

FKBRARcompNMDS2 + BRARcompNMDS2 + BRTEcompNMDS2 + FKBRARinsNMDS2 + BRARinsNMDS2 + BRTEinsNMDS2 + FKBRARmicNMDS2 + BRARmicNMDS2 + BRTEmicNMDS2 + plot_layout(ncol = 3)


#NMDScom2 1900x1400


############################################################################



#####################################################################################
################ Multiple regression lines of functional groups #######
cbPalette_rac <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "red3", "#CC79A7", "#0072B2", "grey24", "purple", "burlywood4" )

#Fort Keogh
#Plants
FK2020cover_sp_BRAR3 <- FK2020cover_sp  %>%
  filter(symbol != "BRAR") %>%
  full_join(FK_BRAR_2020) %>%
  full_join(FKspecies) %>%
  drop_na(location) %>%
  group_by(plot, funct2) %>%
  summarise(avg_cov = sum(cover)) %>%
  ungroup() 

FKcov_func <- FK2020cover_sp_BRAR3 %>%
  full_join(FK_commBRAR)

FK_BRAR_func_reg <- ggplot(data = FKcov_func, aes(x = BRARrel, y = avg_cov, color = funct2, shape = funct2)) +
  #geom_point(size = 3) +
  geom_smooth(data = subset(FKcov_func, funct2 == "C3 Annual Grass"), aes(group = funct2), method = "lm", se = TRUE, linetype = "dashed") +
  geom_smooth(data = subset(FKcov_func, funct2 == "C3 Perennial Grass"), aes(group = funct2), method = "lm", se = TRUE) +
  geom_smooth(data = subset(FKcov_func, funct2 == "C4 Perennial Grass"), aes(group = funct2), method = "lm", se = TRUE) +
  geom_smooth(data = subset(FKcov_func, funct2 == "Cactus"), aes(group = funct2), method = "lm", se = FALSE, linetype = "dotted") +
  geom_smooth(data = subset(FKcov_func, funct2 == "Forb"), aes(group = funct2), method = "lm", se = TRUE, linetype = "dashed") +
  geom_smooth(data = subset(FKcov_func, funct2 == "Sub-Shrub/Shrub"), aes(group = funct2), method = "lm", se = TRUE) +
  scale_colour_manual(values = cbPalette_rac) +
  #scale_shape_manual(values = c(19, 18, 17, 15, 8, 7)) +
  labs(x = "BRAR Cover (%)", y = "Cover (%)", title = "MT BRAR", color = "Functional Group", shape = "Functional Group") + 
  theme(plot.title = element_text(vjust = - 7, hjust = 0.07, size = 16)) 

#C3 AG
c3ag_fk <- lmerTest::lmer(avg_cov ~ BRARrel + (1|grad_num), data = subset(FKcov_func, funct2 == "C3 Annual Grass"))
summary(c3ag_fk) #p = 0.0626; slope = 0.012081
anova(c3ag_fk)

#check normality
resid1.0 <- lm(data = subset(FKcov_func, funct2 == "C3 Annual Grass"), avg_cov ~ BRARrel)
ols_plot_resid_hist(resid1.0) #right skew
ols_test_normality(resid1.0) #fail
#log trans
FKcov_func_c3ag_trans <- FKcov_func %>%
  filter(funct2 == "C3 Annual Grass") %>%
  mutate(ln = log10(avg_cov + 0.001)) 
resid1.0 <- lm(data = FKcov_func_c3ag_trans, ln ~ BRARrel)
ols_plot_resid_hist(resid1.0) #right skew, but better; go with log
ols_test_normality(resid1.0) #fail

c3ag_fk2 <- lmerTest::lmer(ln ~ BRARrel + (1|grad_num), data = FKcov_func_c3ag_trans)
summary(c3ag_fk2) #p = 0.3957; slope = 0.007547
anova(c3ag_fk2)

#linearity
plot(resid(c3ag_fk2), FKcov_func_c3ag_trans$ln) #looks more linear
plot(c3ag_fk2) #not much pattern so indicates linearity

#homoscedascity
FKcov_func_c3ag_trans$res <- residuals(c3ag_fk2)
FKcov_func_c3ag_trans$abs_res <- abs(FKcov_func_c3ag_trans$res)
FKcov_func_c3ag_trans$abs_res2 <- FKcov_func_c3ag_trans$abs_res^2
levene_brar <- lm(data = FKcov_func_c3ag_trans, abs_res2 ~ BRARrel)
anova(levene_brar) #p = 0.1024 so > 0.05 so equal variance is met

#autocorrelation
acf(residuals(c3ag_fk2, retype = "normalized"))  #lines dont go outside CI horizontal lines, so not autocorrelated
pacf(residuals(c3ag_fk2, retype = "normalized"))

#C3 PG
c3pg_fk <- lmerTest::lmer(avg_cov ~ BRARrel + (1|grad_num), data = subset(FKcov_func, funct2 == "C3 Perennial Grass"))
summary(c3pg_fk) #p = 0.00039; slope = -0.2265
anova(c3pg_fk)

#check normality
resid1.0 <- lm(data = subset(FKcov_func, funct2 == "C3 Perennial Grass"), avg_cov ~ BRARrel)
ols_plot_resid_hist(resid1.0) #normal
ols_test_normality(resid1.0) #pass

#linearity
FKcov_func_c3pg <- FKcov_func %>%
  filter(funct2 == "C3 Perennial Grass")
plot(resid(c3pg_fk), FKcov_func_c3pg$avg_cov) #looks linear
plot(c3pg_fk) #not much pattern so indicates linearity

#homoscedascity
FKcov_func_c3pg$res <- residuals(c3pg_fk)
FKcov_func_c3pg$abs_res <- abs(FKcov_func_c3pg$res)
FKcov_func_c3pg$abs_res2 <- FKcov_func_c3pg$abs_res^2
levene_brar <- lm(data = FKcov_func_c3pg, abs_res2 ~ BRARrel)
anova(levene_brar) #p = 0.5444 so > 0.05 so equal variance is met

#autocorrelation
acf(residuals(c3pg_fk, retype = "normalized"))  #lines dont go outside CI horizontal lines, so not autocorrelated
pacf(residuals(c3pg_fk, retype = "normalized"))


#C4 PG
c4pg_fk <- lmerTest::lmer(avg_cov ~ BRARrel + (1|grad_num), data = subset(FKcov_func, funct2 == "C4 Perennial Grass"))
summary(c4pg_fk) #p = 0.000726; slope = -0.12763
anova(c4pg_fk)

#check normality
resid1.0 <- lm(data = subset(FKcov_func, funct2 == "C4 Perennial Grass"), avg_cov ~ BRARrel)
ols_plot_resid_hist(resid1.0) #normal
ols_test_normality(resid1.0) #pass

#linearity
FKcov_func_c4pg <- FKcov_func %>%
  filter(funct2 == "C4 Perennial Grass")
plot(resid(c4pg_fk), FKcov_func_c4pg$avg_cov) #looks linear
plot(c4pg_fk) #no pattern so indicates linearity

#homoscedascity
FKcov_func_c4pg$res <- residuals(c4pg_fk)
FKcov_func_c4pg$abs_res <- abs(FKcov_func_c4pg$res)
FKcov_func_c4pg$abs_res2 <- FKcov_func_c4pg$abs_res^2
levene_brar <- lm(data = FKcov_func_c4pg, abs_res2 ~ BRARrel)
anova(levene_brar) #p = 0.1255 so > 0.05 so equal variance is met

#autocorrelation
acf(residuals(c4pg_fk, retype = "normalized"))  #lines dont go outside CI horizontal lines, so not autocorrelated
pacf(residuals(c4pg_fk, retype = "normalized"))

#Cactus
cc_fk <- lmerTest::lmer(avg_cov ~ BRARrel + (1|grad_num), data = subset(FKcov_func, funct2 == "Cactus"))
summary(cc_fk) #not sig
anova(cc_fk)

#check normality
resid1.0 <- lm(data = subset(FKcov_func, funct2 == "Cactus"), avg_cov ~ BRARrel)
ols_plot_resid_hist(resid1.0) #right skew 
ols_test_normality(resid1.0) #fail
#trans
FKcov_func_cc_trans <- FKcov_func %>%
  filter(funct2 == "Cactus") %>%
  mutate(ln = log10(avg_cov + 0.001))
resid1.0 <- lm(data = FKcov_func_cc_trans, ln ~ BRARrel)
ols_plot_resid_hist(resid1.0) #better
ols_test_normality(resid1.0) #fail
#redo stats
cc_fk2 <- lmerTest::lmer(ln ~ BRARrel + (1|grad_num), data = FKcov_func_cc_trans)
summary(cc_fk2) #not sig
anova(cc_fk2)

#linearity
plot(resid(cc_fk2), FKcov_func_cc_trans$ln) #looks linear
plot(cc_fk2) #no pattern so indicates linearity

#homoscedascity
FKcov_func_cc_trans$res <- residuals(cc_fk2)
FKcov_func_cc_trans$abs_res <- abs(FKcov_func_cc_trans$res)
FKcov_func_cc_trans$abs_res2 <- FKcov_func_cc_trans$abs_res^2
levene_brar <- lm(data = FKcov_func_cc_trans, abs_res2 ~ BRARrel)
anova(levene_brar) #p = 0.2365 so > 0.05 so equal variance is met

#autocorrelation
acf(residuals(cc_fk2, retype = "normalized"))  #lines dont go outside CI horizontal lines, so not autocorrelated
pacf(residuals(cc_fk2, retype = "normalized"))


#Forb
cf_fk <- lmerTest::lmer(avg_cov ~ BRARrel + (1|grad_num), data = subset(FKcov_func, funct2 == "Forb"))
summary(cf_fk) #p = 0.0876; slope = -0.06572
anova(cf_fk)

#check normality
resid1.0 <- lm(data = subset(FKcov_func, funct2 == "Forb"), avg_cov ~ BRARrel)
ols_plot_resid_hist(resid1.0) #normal
ols_test_normality(resid1.0) #pass

#linearity
FKcov_func_cf <- FKcov_func %>%
  filter(funct2 == "Forb")
plot(resid(cf_fk), FKcov_func_cf$avg_cov) #looks linear
plot(cf_fk) #no pattern so indicates linearity

#homoscedascity
FKcov_func_cf$res <- residuals(cf_fk)
FKcov_func_cf$abs_res <- abs(FKcov_func_cf$res)
FKcov_func_cf$abs_res2 <- FKcov_func_cf$abs_res^2
levene_brar <- lm(data = FKcov_func_cf, abs_res2 ~ BRARrel)
anova(levene_brar) #p = 0.9831 so > 0.05 so equal variance is met

#autocorrelation
acf(residuals(cf_fk, retype = "normalized"))  #lines dont go outside CI horizontal lines, so not autocorrelated
pacf(residuals(cf_fk, retype = "normalized"))

#Shrub
cs_fk <- lmerTest::lmer(avg_cov ~ BRARrel + (1|grad_num), data = subset(FKcov_func, funct2 == "Sub-Shrub/Shrub"))
summary(cs_fk) #p = 0.000365; slope = -0.16169
anova(cs_fk)

#check normality
resid1.0 <- lm(data = subset(FKcov_func, funct2 == "Sub-Shrub/Shrub"), avg_cov ~ BRARrel)
ols_plot_resid_hist(resid1.0) #right skew
ols_test_normality(resid1.0) #fail
#trans
FKcov_func_cs_trans <- FKcov_func %>%
  filter(funct2 == "Sub-Shrub/Shrub") %>%
  mutate(ln = log10(avg_cov + 0.001))
resid1.0 <- lm(data = FKcov_func_cs_trans, ln ~ BRARrel)
ols_plot_resid_hist(resid1.0) #better
ols_test_normality(resid1.0) #pass some
#redo stats
cs_fk2 <- lmerTest::lmer(ln ~ BRARrel + (1|grad_num), data = FKcov_func_cs_trans)
summary(cs_fk2) #p = 0.001917; slope = -0.03603
anova(cs_fk2)

#linearity
plot(resid(cs_fk2), FKcov_func_cs_trans$ln) #looks linear
plot(cs_fk2) #no pattern so indicates linearity

#homoscedascity
FKcov_func_cs_trans$res <- residuals(cs_fk2)
FKcov_func_cs_trans$abs_res <- abs(FKcov_func_cs_trans$res)
FKcov_func_cs_trans$abs_res2 <- FKcov_func_cs_trans$abs_res^2
levene_brar <- lm(data = FKcov_func_cs_trans, abs_res2 ~ BRARrel)
anova(levene_brar) #p = 0.8295 so > 0.05 so equal variance is met

#autocorrelation
acf(residuals(cs_fk2, retype = "normalized"))  #lines dont go outside CI horizontal lines, so not autocorrelated
pacf(residuals(cs_fk2, retype = "normalized"))



#Insects
FKinsID_guild_BRAR3 <- FKinsID_guild_BRAR %>%
  group_by(plot, guild) %>%
  summarise(avg = sum(avg_n)) %>%
  ungroup()

FKinsfunc <- FKinsID_guild_BRAR3 %>%
  full_join(FK_commBRAR)

FK_BRAR_func_reg_ins <- ggplot(data = FKinsfunc, aes(x = BRARrel, y = avg, color = guild, shape = guild)) +
  #geom_point(size = 3) +
  geom_smooth(data = subset(FKinsfunc, guild == "Leaf Chewing Herbivore"), aes(group = guild), method = "lm", se = TRUE) +
  geom_smooth(data = subset(FKinsfunc, guild == "Predator"), aes(group = guild), method = "lm", se = FALSE, linetype = "dotted") +
  geom_smooth(data = subset(FKinsfunc, guild == "Sap Sucking Herbivore"), aes(group = guild), method = "lm", se = TRUE, linetype = "dashed") +
  geom_smooth(data = subset(FKinsfunc, guild == "Parasitoid"), aes(group = guild), method = "lm", se = FALSE, linetype = "dotted") +
  scale_colour_manual(values = cbPalette_rac) +
  #scale_shape_manual(values = c(19, 18, 17, 15, 8, 7)) +
  labs(x = "BRAR Cover (%)", y = "Abundance (count)", title = "MT BRAR", color = "Functional Group", shape = "Functional Group") + 
  theme(plot.title = element_text(vjust = - 7, hjust = 0.07, size = 16)) 


det_fk <- lmerTest::lmer(avg ~ BRARrel + (1|grad_num), data = subset(FKinsfunc, guild == "Detritivore"))
summary(det_fk) #can't include; only 1 found
anova(det_fk)

#Leaf chewers
lch_fk <- lmerTest::lmer(avg ~ BRARrel + (1|grad_num), data = subset(FKinsfunc, guild == "Leaf Chewing Herbivore"))
summary(lch_fk) #p = 0.00439; slope = -0.07776
anova(lch_fk)

#check normality
resid1.0 <- lm(data = subset(FKinsfunc, guild == "Leaf Chewing Herbivore"), avg ~ BRARrel)
ols_plot_resid_hist(resid1.0) #right skew
ols_test_normality(resid1.0) #fail

#trans
FKinsfunctrans <- FKinsfunc %>%
  filter(guild == "Leaf Chewing Herbivore") %>%
  mutate(ln = log10(avg + 0.001))
resid1.0 <- lm(data = FKinsfunctrans, ln ~ BRARrel)
ols_plot_resid_hist(resid1.0) #better
ols_test_normality(resid1.0) #pass
#redo stats
lch_fk2 <- lmerTest::lmer(ln ~ BRARrel + (1|grad_num), data = FKinsfunctrans)
summary(lch_fk2) #p = 0.000512; slope = -0.006031
anova(lch_fk2)

#linearity
plot(resid(lch_fk2), FKinsfunctrans$ln) #looks linear
plot(lch_fk2) #no pattern so indicates linearity

#homoscedascity
FKinsfunctrans$res <- residuals(lch_fk2)
FKinsfunctrans$abs_res <- abs(FKinsfunctrans$res)
FKinsfunctrans$abs_res2 <- FKinsfunctrans$abs_res^2
levene_brar <- lm(data = FKinsfunctrans, abs_res2 ~ BRARrel)
anova(levene_brar) #p = 0.5358 so > 0.05 so equal variance is met

#autocorrelation
acf(residuals(lch_fk2, retype = "normalized"))  #lines dont go outside CI horizontal lines, so not autocorrelated
pacf(residuals(lch_fk2, retype = "normalized"))


#Predator
pre_fk <- lmerTest::lmer(avg ~ BRARrel + (1|grad_num), data = subset(FKinsfunc, guild == "Predator"))
summary(pre_fk) #no sig
anova(pre_fk)

#check normality
resid1.0 <- lm(data = subset(FKinsfunc, guild == "Predator"), avg ~ BRARrel)
ols_plot_resid_hist(resid1.0) #right skew
ols_test_normality(resid1.0) #fail

#trans
FKinsfunc_predtrans <- FKinsfunc %>%
  filter(guild == "Predator") %>%
  mutate(ln = log10(avg + 0.001))
resid1.0 <- lm(data = FKinsfunc_predtrans, ln ~ BRARrel)
ols_plot_resid_hist(resid1.0) #better
ols_test_normality(resid1.0) #pass
#redo stats
pre_fk2 <- lmerTest::lmer(ln ~ BRARrel + (1|grad_num), data = FKinsfunc_predtrans)
summary(pre_fk2) #p = 0.3462; slope = -0.001289
anova(pre_fk2)

#linearity
plot(resid(pre_fk2), FKinsfunc_predtrans$ln) #looks linear
plot(pre_fk2) #no pattern so indicates linearity

#homoscedascity
FKinsfunc_predtrans$res <- residuals(pre_fk2)
FKinsfunc_predtrans$abs_res <- abs(FKinsfunc_predtrans$res)
FKinsfunc_predtrans$abs_res2 <- FKinsfunc_predtrans$abs_res^2
levene_brar <- lm(data = FKinsfunc_predtrans, abs_res2 ~ BRARrel)
anova(levene_brar) #p = 0.9283 so > 0.05 so equal variance is met

#autocorrelation
acf(residuals(pre_fk2, retype = "normalized"))  #lines dont go outside CI horizontal lines, so not autocorrelated
pacf(residuals(pre_fk2, retype = "normalized"))


#Sap sucking
ssh_fk <- lmerTest::lmer(avg ~ BRARrel + (1|grad_num), data = subset(FKinsfunc, guild == "Sap Sucking Herbivore"))
summary(ssh_fk) #p = 0.0599; slope = -0.13499
anova(ssh_fk)

#check normality
resid1.0 <- lm(data = subset(FKinsfunc, guild == "Sap Sucking Herbivore"), avg ~ BRARrel)
ols_plot_resid_hist(resid1.0) #right skew
ols_test_normality(resid1.0) #fail

#trans
FKinsfunc_sshtrans <- FKinsfunc %>%
  filter(guild == "Sap Sucking Herbivore") %>%
  mutate(ln = log10(avg + 0.001))
resid1.0 <- lm(data = FKinsfunc_sshtrans, ln ~ BRARrel)
ols_plot_resid_hist(resid1.0) #better
ols_test_normality(resid1.0) #pass
#redo stats
ssh_fk2 <- lmerTest::lmer(ln ~ BRARrel + (1|grad_num), data = FKinsfunc_sshtrans)
summary(ssh_fk2) #p = 0.02978; slope = -0.006083
anova(ssh_fk2)

#linearity
plot(resid(ssh_fk2), FKinsfunc_sshtrans$ln) #looks linear
plot(ssh_fk2) #no pattern so indicates linearity

#homoscedascity
FKinsfunc_sshtrans$res <- residuals(ssh_fk2)
FKinsfunc_sshtrans$abs_res <- abs(FKinsfunc_sshtrans$res)
FKinsfunc_sshtrans$abs_res2 <- FKinsfunc_sshtrans$abs_res^2
levene_brar <- lm(data = FKinsfunc_sshtrans, abs_res2 ~ BRARrel)
anova(levene_brar) #p = 0.7443 so > 0.05 so equal variance is met

#autocorrelation
acf(residuals(ssh_fk2, retype = "normalized"))  #lines dont go outside CI horizontal lines, so not autocorrelated
pacf(residuals(ssh_fk2, retype = "normalized"))


#Parasitoid
par_fk <- lmerTest::lmer(avg ~ BRARrel + (1|grad_num), data = subset(FKinsfunc, guild == "Parasitoid"))
summary(par_fk) #no sig
anova(par_fk)

#check normality
resid1.0 <- lm(data = subset(FKinsfunc, guild == "Parasitoid"), avg ~ BRARrel)
ols_plot_resid_hist(resid1.0) #right skew
ols_test_normality(resid1.0) #fail - stick with untransformed; that was no better

#linearity
FKinsfunc_par <- FKinsfunc %>%
  filter(guild == "Parasitoid")
plot(resid(par_fk), FKinsfunc_partrans$ln) #looks linear
plot(par_fk2) #no pattern so indicates linearity

#homoscedascity
FKinsfunc_par$res <- residuals(par_fk)
FKinsfunc_par$abs_res <- abs(FKinsfunc_par$res)
FKinsfunc_par$abs_res2 <- FKinsfunc_par$abs_res^2
levene_brar <- lm(data = FKinsfunc_par, abs_res2 ~ BRARrel)
anova(levene_brar) #p = 0.3209 so > 0.05 so equal variance is met

#autocorrelation
acf(residuals(par_fk, retype = "normalized"))  #lines dont go outside CI horizontal lines, so not autocorrelated
pacf(residuals(par_fk, retype = "normalized"))


pneh_fk <- lmerTest::lmer(avg ~ BRARrel + (1|grad_num), data = subset(FKinsfunc, guild == "Pollen/Nectar Eating Herbivore"))
summary(pneh_fk) #not enough data to include
anova(pneh_fk)

oh_fk <- lmerTest::lmer(avg ~ BRARrel + (1|grad_num), data = subset(FKinsfunc, guild == "Other Herbivore"))
summary(oh_fk) #not enough data to include



#Microbes
#FKmicfunc <- BRARFKmic_func_long2 %>%
  dplyr::select(func) %>%
  unique()
#write.csv(FKmicfunc, "fkmic.csv")

BRARFKmic_func_long2 <- FKmic_func_long %>%
  full_join(FK_BRAR_2020) %>%
  group_by(plot, func) %>%
  summarise(avg = sum(count)) %>%
  ungroup() %>%
  full_join(FK_commBRAR) %>%
  drop_na(func)

BRARFKmic_func_long3 <- BRARFKmic_func_long2 %>%
  filter(func != "other")

BRARFKmic_func_long4 <- BRARFKmic_func_long2 %>%
  group_by(func) %>%
  summarise(avg2 = mean(avg)) %>%
  ungroup()
#top 7 most abundant functions (not other): chemoheterotrophy, aerobic_chemoheterotrophy, nitrate_reduction, nitrification, aerobic_ammonia_oxidation, nitrogen_fixation, manganese_oxidation 


FK_BRAR_func_reg_mic <- ggplot(data = BRARFKmic_func_long3, aes(x = BRARrel, y = avg, color = func)) +
  #geom_point(size = 3) +
  geom_smooth(data = subset(BRARFKmic_func_long3, func == "oxygenic_photoautotrophy"), aes(group = func), method = "lm", se = TRUE, linetype = "dashed", position = position_dodge(0.5)) +
  geom_smooth(data = subset(BRARFKmic_func_long3, func == "photosynthetic_cyanobacteria"), aes(group = func), method = "lm", se = TRUE, linetype = "dashed", position = position_dodge(0.5)) +
  #geom_smooth(data = subset(BRARFKmic_func_long3, func != "oxygenic_photoautotrophy"), aes(group = func), method = "lm", se = FALSE, linetype = "dotted") +
  #scale_colour_manual(values = cbPalette_rac) +
  #scale_shape_manual(values = c(19, 18, 17, 15, 8, 7)) +
  labs(x = "BRAR Cover (%)", y = "Abundance (count)", title = "MT BRAR", color = "Functional Group", shape = "Functional Group") + 
  theme(plot.title = element_text(vjust = - 7, hjust = 0.07, size = 16)) 

FK_BRAR_func_reg_mic2 <- ggplot(data = BRARFKmic_func_long3, aes(x = BRARrel, y = avg, color = func)) +
  #geom_point(size = 3) +
  geom_smooth(data = subset(BRARFKmic_func_long3, func == "chemoheterotrophy"), aes(group = func), method = "lm", se = FALSE, linetype = "dotted", position = position_dodge(0.5)) +
  geom_smooth(data = subset(BRARFKmic_func_long3, func == "aerobic_chemoheterotrophy"), aes(group = func), method = "lm", se = FALSE, linetype = "dotted", position = position_dodge(0.5)) +
  geom_smooth(data = subset(BRARFKmic_func_long3, func == "nitrate_reduction"), aes(group = func), method = "lm", se = FALSE, linetype = "dotted") +
  geom_smooth(data = subset(BRARFKmic_func_long3, func == "nitrification"), aes(group = func), method = "lm", se = FALSE, linetype = "dotted") +
  geom_smooth(data = subset(BRARFKmic_func_long3, func == "aerobic_ammonia_oxidation"), aes(group = func), method = "lm", se = FALSE, linetype = "dotted") +
  geom_smooth(data = subset(BRARFKmic_func_long3, func == "nitrogen_fixation"), aes(group = func), method = "lm", se = FALSE, linetype = "dotted") +
  geom_smooth(data = subset(BRARFKmic_func_long3, func == "manganese_oxidation"), aes(group = func), method = "lm", se = FALSE, linetype = "dotted") +
  scale_colour_manual(values = cbPalette_rac) +
  #scale_shape_manual(values = c(19, 18, 17, 15, 8, 7)) +
  labs(x = "BRAR Cover (%)", y = "Abundance (count)", title = "MT BRAR", color = "Functional Group", shape = "Functional Group") + 
  theme(plot.title = element_text(vjust = - 7, hjust = 0.07, size = 16)) 


#Aerobic ammonia oxidation
fk_mic1 <- lmerTest::lmer(avg ~ BRARrel + (1|grad_num), data = subset(BRARFKmic_func_long2, func == "aerobic_ammonia_oxidation"))
summary(fk_mic1) #not sig
anova(fk_mic1)

#check normality
resid1.0 <- lm(data = subset(BRARFKmic_func_long2, func == "aerobic_ammonia_oxidation"), avg ~ BRARrel)
ols_plot_resid_hist(resid1.0) #right skew
ols_test_normality(resid1.0) #fail

#trans
BRARFKmic_func_mic1 <- BRARFKmic_func_long2 %>%
  filter(func == "aerobic_ammonia_oxidation") %>%
  mutate(ln = log10(avg + 0.001))
resid1.0 <- lm(data = BRARFKmic_func_mic1, ln ~ BRARrel)
ols_plot_resid_hist(resid1.0) #better
ols_test_normality(resid1.0) #pass
#redo stats
fk_mic1a <- lmerTest::lmer(ln ~ BRARrel + (1|grad_num), data = BRARFKmic_func_mic1)
summary(fk_mic1a) #p = 0.8392; slope = -0.0003897
anova(fk_mic1a)

#linearity
plot(resid(fk_mic1a), BRARFKmic_func_mic1$ln) #looks linear
plot(fk_mic1a) #no pattern so indicates linearity

#homoscedascity
BRARFKmic_func_mic1$res <- residuals(fk_mic1a)
BRARFKmic_func_mic1$abs_res <- abs(BRARFKmic_func_mic1$res)
BRARFKmic_func_mic1$abs_res2 <- BRARFKmic_func_mic1$abs_res^2
levene_brar <- lm(data = BRARFKmic_func_mic1, abs_res2 ~ BRARrel)
anova(levene_brar) #p = 0.6529 so > 0.05 so equal variance is met

#autocorrelation
acf(residuals(fk_mic1a, retype = "normalized"))  #lines dont go outside CI horizontal lines, so not autocorrelated
pacf(residuals(fk_mic1a, retype = "normalized"))


#Aerobic chemo
fk_mic2 <- lmerTest::lmer(avg ~ BRARrel + (1|grad_num), data = subset(BRARFKmic_func_long2, func == "aerobic_chemoheterotrophy"))
summary(fk_mic2) #not sig
anova(fk_mic2)

#check normality
resid1.0 <- lm(data = subset(BRARFKmic_func_long2, func == "aerobic_chemoheterotrophy"), avg ~ BRARrel)
ols_plot_resid_hist(resid1.0) #right skew
ols_test_normality(resid1.0) #fail

#trans
BRARFKmic_func_mic2 <- BRARFKmic_func_long2 %>%
  filter(func == "aerobic_chemoheterotrophy") %>%
  mutate(ln = log10(avg + 0.001))
resid1.0 <- lm(data = BRARFKmic_func_mic2, ln ~ BRARrel)
ols_plot_resid_hist(resid1.0) #better
ols_test_normality(resid1.0) #pass
#redo stats
fk_mic2a <- lmerTest::lmer(ln ~ BRARrel + (1|grad_num), data = BRARFKmic_func_mic2)
summary(fk_mic2a) #p = 0.848; slope = 3.364e-4
anova(fk_mic2a)

#linearity
plot(resid(fk_mic2a), BRARFKmic_func_mic2$ln) #looks linear
plot(fk_mic2a) #no pattern so indicates linearity

#homoscedascity
BRARFKmic_func_mic2$res <- residuals(fk_mic2a)
BRARFKmic_func_mic2$abs_res <- abs(BRARFKmic_func_mic2$res)
BRARFKmic_func_mic2$abs_res2 <- BRARFKmic_func_mic2$abs_res^2
levene_brar <- lm(data = BRARFKmic_func_mic2, abs_res2 ~ BRARrel)
anova(levene_brar) #p = 0.6373 so > 0.05 so equal variance is met

#autocorrelation
acf(residuals(fk_mic2a, retype = "normalized"))  #lines dont go outside CI horizontal lines, so not autocorrelated
pacf(residuals(fk_mic2a, retype = "normalized"))


#Aerobic nitrite oxid
fk_mic3 <- lmerTest::lmer(avg ~ BRARrel + (1|grad_num), data = subset(BRARFKmic_func_long2, func == "aerobic_nitrite_oxidation"))
summary(fk_mic3) #not sig
anova(fk_mic3)

#check normality
resid1.0 <- lm(data = subset(BRARFKmic_func_long2, func == "aerobic_nitrite_oxidation"), avg ~ BRARrel)
ols_plot_resid_hist(resid1.0) #right skew
ols_test_normality(resid1.0) #fail

#trans
BRARFKmic_func_mic3 <- BRARFKmic_func_long2 %>%
  filter(func == "aerobic_nitrite_oxidation") %>%
  mutate(ln = log10(avg + 0.001))
resid1.0 <- lm(data = BRARFKmic_func_mic3, ln ~ BRARrel)
ols_plot_resid_hist(resid1.0) #better
ols_test_normality(resid1.0) #pass
#redo stats
fk_mic3a <- lmerTest::lmer(ln ~ BRARrel + (1|grad_num), data = BRARFKmic_func_mic3)
summary(fk_mic3a) #p = 0.4751; slope = 0.002111
anova(fk_mic3a)

#linearity
plot(resid(fk_mic3a), BRARFKmic_func_mic3$ln) #looks linear
plot(fk_mic3a) #no pattern so indicates linearity

#homoscedascity
BRARFKmic_func_mic3$res <- residuals(fk_mic3a)
BRARFKmic_func_mic3$abs_res <- abs(BRARFKmic_func_mic3$res)
BRARFKmic_func_mic3$abs_res2 <- BRARFKmic_func_mic3$abs_res^2
levene_brar <- lm(data = BRARFKmic_func_mic3, abs_res2 ~ BRARrel)
anova(levene_brar) #p = 0.4791 so > 0.05 so equal variance is met

#autocorrelation
acf(residuals(fk_mic3a, retype = "normalized"))  #lines dont go outside CI horizontal lines, so not autocorrelated
pacf(residuals(fk_mic3a, retype = "normalized"))


#Animal parasites...
fk_mic4 <- lmerTest::lmer(avg ~ BRARrel + (1|grad_num), data = subset(BRARFKmic_func_long2, func == "animal_parasites_or_symbionts"))
summary(fk_mic4) #not sig
anova(fk_mic4)

#check normality
resid1.0 <- lm(data = subset(BRARFKmic_func_long2, func == "animal_parasites_or_symbionts"), avg ~ BRARrel)
ols_plot_resid_hist(resid1.0) #right skew
ols_test_normality(resid1.0) #fail

#trans
BRARFKmic_func_mic4 <- BRARFKmic_func_long2 %>%
  filter(func == "animal_parasites_or_symbionts") %>%
  mutate(ln = log10(avg + 0.001))
resid1.0 <- lm(data = BRARFKmic_func_mic4, ln ~ BRARrel)
ols_plot_resid_hist(resid1.0) #better
ols_test_normality(resid1.0) #pass
#redo stats
fk_mic4a <- lmerTest::lmer(ln ~ BRARrel + (1|grad_num), data = BRARFKmic_func_mic4)
summary(fk_mic4a) #p = 0.3501; slope = 0.002794
anova(fk_mic4a)

#linearity
plot(resid(fk_mic4a), BRARFKmic_func_mic4$ln) #looks linear
plot(fk_mic4a) #no pattern so indicates linearity

#homoscedascity
BRARFKmic_func_mic4$res <- residuals(fk_mic4a)
BRARFKmic_func_mic4$abs_res <- abs(BRARFKmic_func_mic4$res)
BRARFKmic_func_mic4$abs_res2 <- BRARFKmic_func_mic4$abs_res^2
levene_brar <- lm(data = BRARFKmic_func_mic4, abs_res2 ~ BRARrel)
anova(levene_brar) #p = 0.9953 so > 0.05 so equal variance is met

#autocorrelation
acf(residuals(fk_mic4a, retype = "normalized"))  #lines dont go outside CI horizontal lines, so not autocorrelated
pacf(residuals(fk_mic4a, retype = "normalized"))


#Anoxygen photo
fk_mic5 <- lmerTest::lmer(avg ~ BRARrel + (1|grad_num), data = subset(BRARFKmic_func_long2, func == "anoxygenic_photoautotrophy"))
summary(fk_mic5) #not sig 
anova(fk_mic5)

#check normality
resid1.0 <- lm(data = subset(BRARFKmic_func_long2, func == "anoxygenic_photoautotrophy"), avg ~ BRARrel)
ols_plot_resid_hist(resid1.0) #right skew
ols_test_normality(resid1.0) #fail

#trans
BRARFKmic_func_mic5 <- BRARFKmic_func_long2 %>%
  filter(func == "anoxygenic_photoautotrophy") %>%
  mutate(ln = log10(avg + 0.001))
resid1.0 <- lm(data = BRARFKmic_func_mic5, ln ~ BRARrel)
ols_plot_resid_hist(resid1.0) #better
ols_test_normality(resid1.0) #pass
#redo stats
fk_mic5a <- lmerTest::lmer(ln ~ BRARrel + (1|grad_num), data = BRARFKmic_func_mic5)
summary(fk_mic5a) #p = 0.7558; slope = 7.918e-4
anova(fk_mic5a)

#linearity
plot(resid(fk_mic5a), BRARFKmic_func_mic5$ln) #looks linear
plot(fk_mic5a) #no pattern so indicates linearity

#homoscedascity
BRARFKmic_func_mic5$res <- residuals(fk_mic5a)
BRARFKmic_func_mic5$abs_res <- abs(BRARFKmic_func_mic5$res)
BRARFKmic_func_mic5$abs_res2 <- BRARFKmic_func_mic5$abs_res^2
levene_brar <- lm(data = BRARFKmic_func_mic5, abs_res2 ~ BRARrel)
anova(levene_brar) #p = 0.9726 so > 0.05 so equal variance is met

#autocorrelation
acf(residuals(fk_mic5a, retype = "normalized"))  #lines dont go outside CI horizontal lines, so not autocorrelated
pacf(residuals(fk_mic5a, retype = "normalized"))


#Anox photo S
fk_mic6 <- lmerTest::lmer(avg ~ BRARrel + (1|grad_num), data = subset(BRARFKmic_func_long2, func == "anoxygenic_photoautotrophy_S_oxidizing"))
summary(fk_mic6) #not sig 
anova(fk_mic6)

#check normality
resid1.0 <- lm(data = subset(BRARFKmic_func_long2, func == "anoxygenic_photoautotrophy_S_oxidizing"), avg ~ BRARrel)
ols_plot_resid_hist(resid1.0) #right skew
ols_test_normality(resid1.0) #fail

#trans
BRARFKmic_func_mic6 <- BRARFKmic_func_long2 %>%
  filter(func == "anoxygenic_photoautotrophy_S_oxidizing") %>%
  mutate(ln = log10(avg + 0.001))
resid1.0 <- lm(data = BRARFKmic_func_mic6, ln ~ BRARrel)
ols_plot_resid_hist(resid1.0) #better
ols_test_normality(resid1.0) #pass
#redo stats
fk_mic6a <- lmerTest::lmer(ln ~ BRARrel + (1|grad_num), data = BRARFKmic_func_mic6)
summary(fk_mic6a) #p = 0.7558; slope = 7.918e-4
anova(fk_mic6a)

#linearity
plot(resid(fk_mic6a), BRARFKmic_func_mic6$ln) #looks linear
plot(fk_mic6a) #no pattern so indicates linearity

#homoscedascity
BRARFKmic_func_mic6$res <- residuals(fk_mic6a)
BRARFKmic_func_mic6$abs_res <- abs(BRARFKmic_func_mic6$res)
BRARFKmic_func_mic6$abs_res2 <- BRARFKmic_func_mic6$abs_res^2
levene_brar <- lm(data = BRARFKmic_func_mic6, abs_res2 ~ BRARrel)
anova(levene_brar) #p = 0.9726 so > 0.05 so equal variance is met

#autocorrelation
acf(residuals(fk_mic6a, retype = "normalized"))  #lines dont go outside CI horizontal lines, so not autocorrelated
pacf(residuals(fk_mic6a, retype = "normalized"))


#Aromatic comp deg
fk_mic7 <- lmerTest::lmer(avg ~ BRARrel + (1|grad_num), data = subset(BRARFKmic_func_long2, func == "aromatic_compound_degradation"))
summary(fk_mic7) #not sig 
anova(fk_mic7)

#check normality
resid1.0 <- lm(data = subset(BRARFKmic_func_long2, func == "aromatic_compound_degradation"), avg ~ BRARrel)
ols_plot_resid_hist(resid1.0) #right skew
ols_test_normality(resid1.0) #fail

#trans
BRARFKmic_func_mic7 <- BRARFKmic_func_long2 %>%
  filter(func == "aromatic_compound_degradation") %>%
  mutate(ln = log10(avg + 0.001))
resid1.0 <- lm(data = BRARFKmic_func_mic7, ln ~ BRARrel)
ols_plot_resid_hist(resid1.0) #better
ols_test_normality(resid1.0) #pass
#redo stats
fk_mic7a <- lmerTest::lmer(ln ~ BRARrel + (1|grad_num), data = BRARFKmic_func_mic7)
summary(fk_mic7a) #p = 0.1986; slope = -0.002333
anova(fk_mic7a)

#linearity
plot(resid(fk_mic7a), BRARFKmic_func_mic7$ln) #looks linear
plot(fk_mic7a) #no pattern so indicates linearity

#homoscedascity
BRARFKmic_func_mic7$res <- residuals(fk_mic7a)
BRARFKmic_func_mic7$abs_res <- abs(BRARFKmic_func_mic7$res)
BRARFKmic_func_mic7$abs_res2 <- BRARFKmic_func_mic7$abs_res^2
levene_brar <- lm(data = BRARFKmic_func_mic7, abs_res2 ~ BRARrel)
anova(levene_brar) #p = 0.5987 so > 0.05 so equal variance is met

#autocorrelation
acf(residuals(fk_mic7a, retype = "normalized"))  #lines dont go outside CI horizontal lines, so not autocorrelated
pacf(residuals(fk_mic7a, retype = "normalized"))


#Chemohetero
fk_mic8 <- lmerTest::lmer(avg ~ BRARrel + (1|grad_num), data = subset(BRARFKmic_func_long2, func == "chemoheterotrophy"))
summary(fk_mic8) #not sig 
anova(fk_mic8)

#check normality
resid1.0 <- lm(data = subset(BRARFKmic_func_long2, func == "chemoheterotrophy"), avg ~ BRARrel)
ols_plot_resid_hist(resid1.0) #right skew
ols_test_normality(resid1.0) #fail

#trans
BRARFKmic_func_mic8 <- BRARFKmic_func_long2 %>%
  filter(func == "chemoheterotrophy") %>%
  mutate(ln = log10(avg + 0.001))
resid1.0 <- lm(data = BRARFKmic_func_mic8, ln ~ BRARrel)
ols_plot_resid_hist(resid1.0) #better
ols_test_normality(resid1.0) #pass
#redo stats
fk_mic8a <- lmerTest::lmer(ln ~ BRARrel + (1|grad_num), data = BRARFKmic_func_mic8)
summary(fk_mic8a) #p = 0.8477; slope = 0.000337
anova(fk_mic8a)

#linearity
plot(resid(fk_mic8a), BRARFKmic_func_mic8$ln) #looks linear
plot(fk_mic8a) #no pattern so indicates linearity

#homoscedascity
BRARFKmic_func_mic8$res <- residuals(fk_mic8a)
BRARFKmic_func_mic8$abs_res <- abs(BRARFKmic_func_mic8$res)
BRARFKmic_func_mic8$abs_res2 <- BRARFKmic_func_mic8$abs_res^2
levene_brar <- lm(data = BRARFKmic_func_mic8, abs_res2 ~ BRARrel)
anova(levene_brar) #p = 0.638 so > 0.05 so equal variance is met

#autocorrelation
acf(residuals(fk_mic8a, retype = "normalized"))  #lines dont go outside CI horizontal lines, so not autocorrelated
pacf(residuals(fk_mic8a, retype = "normalized"))


#Chitinolysis
fk_mic9 <- lmerTest::lmer(avg ~ BRARrel + (1|grad_num), data = subset(BRARFKmic_func_long2, func == "chitinolysis"))
summary(fk_mic9) #not sig
anova(fk_mic9)

#check normality
resid1.0 <- lm(data = subset(BRARFKmic_func_long2, func == "chitinolysis"), avg ~ BRARrel)
ols_plot_resid_hist(resid1.0) #right skew
ols_test_normality(resid1.0) #fail

#trans
BRARFKmic_func_mic9 <- BRARFKmic_func_long2 %>%
  filter(func == "chitinolysis") %>%
  mutate(ln = log10(avg + 0.001))
resid1.0 <- lm(data = BRARFKmic_func_mic9, ln ~ BRARrel)
ols_plot_resid_hist(resid1.0) #better
ols_test_normality(resid1.0) #pass
#redo stats
fk_mic9a <- lmerTest::lmer(ln ~ BRARrel + (1|grad_num), data = BRARFKmic_func_mic9)
summary(fk_mic9a) #p = 0.8936; slope = -0.0003734
anova(fk_mic9a)

#linearity
plot(resid(fk_mic9a), BRARFKmic_func_mic9$ln) #looks linear
plot(fk_mic9a) #no pattern so indicates linearity

#homoscedascity
BRARFKmic_func_mic9$res <- residuals(fk_mic9a)
BRARFKmic_func_mic9$abs_res <- abs(BRARFKmic_func_mic9$res)
BRARFKmic_func_mic9$abs_res2 <- BRARFKmic_func_mic9$abs_res^2
levene_brar <- lm(data = BRARFKmic_func_mic9, abs_res2 ~ BRARrel)
anova(levene_brar) #p = 0.3698 so > 0.05 so equal variance is met

#autocorrelation
acf(residuals(fk_mic9a, retype = "normalized"))  #lines dont go outside CI horizontal lines, so not autocorrelated
pacf(residuals(fk_mic9a, retype = "normalized"))


#Dark oxid sulf
fk_mic10 <- lmerTest::lmer(avg ~ BRARrel + (1|grad_num), data = subset(BRARFKmic_func_long2, func == "dark_oxidation_of_sulfur_compounds"))
summary(fk_mic10) #not sig
anova(fk_mic10)

#check normality
resid1.0 <- lm(data = subset(BRARFKmic_func_long2, func == "dark_oxidation_of_sulfur_compounds"), avg ~ BRARrel)
ols_plot_resid_hist(resid1.0) #right skew
ols_test_normality(resid1.0) #fail

#trans
BRARFKmic_func_mic10 <- BRARFKmic_func_long2 %>%
  filter(func == "dark_oxidation_of_sulfur_compounds") %>%
  mutate(ln = log10(avg + 0.001))
resid1.0 <- lm(data = BRARFKmic_func_mic10, ln ~ BRARrel)
ols_plot_resid_hist(resid1.0) #better
ols_test_normality(resid1.0) #pass
#redo stats
fk_mic10a <- lmerTest::lmer(ln ~ BRARrel + (1|grad_num), data = BRARFKmic_func_mic10)
summary(fk_mic10a) #p = 0.5772; slope = 0.003558
anova(fk_mic10a)

#linearity
plot(resid(fk_mic10a), BRARFKmic_func_mic10$ln) #looks linear
plot(fk_mic10a) #no pattern so indicates linearity

#homoscedascity
BRARFKmic_func_mic10$res <- residuals(fk_mic10a)
BRARFKmic_func_mic10$abs_res <- abs(BRARFKmic_func_mic10$res)
BRARFKmic_func_mic10$abs_res2 <- BRARFKmic_func_mic10$abs_res^2
levene_brar <- lm(data = BRARFKmic_func_mic10, abs_res2 ~ BRARrel)
anova(levene_brar) #p = 0.597 so > 0.05 so equal variance is met

#autocorrelation
acf(residuals(fk_mic10a, retype = "normalized"))  #lines dont go outside CI horizontal lines, so not autocorrelated
pacf(residuals(fk_mic10a, retype = "normalized"))


#Denitification
fk_mic11 <- lmerTest::lmer(avg ~ BRARrel + (1|grad_num), data = subset(BRARFKmic_func_long2, func == "denitrification"))
summary(fk_mic11) #not sig
anova(fk_mic11)

#check normality
resid1.0 <- lm(data = subset(BRARFKmic_func_long2, func == "denitrification"), avg ~ BRARrel)
ols_plot_resid_hist(resid1.0) #right skew
ols_test_normality(resid1.0) #fail

#trans
BRARFKmic_func_mic11 <- BRARFKmic_func_long2 %>%
  filter(func == "denitrification") %>%
  mutate(ln = log10(avg + 0.001))
resid1.0 <- lm(data = BRARFKmic_func_mic11, ln ~ BRARrel)
ols_plot_resid_hist(resid1.0) #better
ols_test_normality(resid1.0) #pass
#redo stats
fk_mic11a <- lmerTest::lmer(ln ~ BRARrel + (1|grad_num), data = BRARFKmic_func_mic11)
summary(fk_mic11a) #p = 0.7558; slope = 7.918e-04
anova(fk_mic11a)

#linearity
plot(resid(fk_mic11a), BRARFKmic_func_mic11$ln) #looks linear
plot(fk_mic11a) #no pattern so indicates linearity

#homoscedascity
BRARFKmic_func_mic11$res <- residuals(fk_mic11a)
BRARFKmic_func_mic11$abs_res <- abs(BRARFKmic_func_mic11$res)
BRARFKmic_func_mic11$abs_res2 <- BRARFKmic_func_mic11$abs_res^2
levene_brar <- lm(data = BRARFKmic_func_mic11, abs_res2 ~ BRARrel)
anova(levene_brar) #p = 0.9726 so > 0.05 so equal variance is met

#autocorrelation
acf(residuals(fk_mic11a, retype = "normalized"))  #lines dont go outside CI horizontal lines, so not autocorrelated
pacf(residuals(fk_mic11a, retype = "normalized"))


#Fermentation
fk_mic12 <- lmerTest::lmer(avg ~ BRARrel + (1|grad_num), data = subset(BRARFKmic_func_long2, func == "fermentation"))
summary(fk_mic12) #not sig
anova(fk_mic12)

#check normality
resid1.0 <- lm(data = subset(BRARFKmic_func_long2, func == "fermentation"), avg ~ BRARrel)
ols_plot_resid_hist(resid1.0) #right skew
ols_test_normality(resid1.0) #fail

#trans
BRARFKmic_func_mic12 <- BRARFKmic_func_long2 %>%
  filter(func == "fermentation") %>%
  mutate(ln = log10(avg + 0.001))
resid1.0 <- lm(data = BRARFKmic_func_mic12, ln ~ BRARrel)
ols_plot_resid_hist(resid1.0) #better
ols_test_normality(resid1.0) #pass
#redo stats
fk_mic12a <- lmerTest::lmer(ln ~ BRARrel + (1|grad_num), data = BRARFKmic_func_mic12)
summary(fk_mic12a) #p = 0.809; slope = -0.0006924
anova(fk_mic12a)

#linearity
plot(resid(fk_mic12a), BRARFKmic_func_mic12$ln) #looks linear
plot(fk_mic12a) #no pattern so indicates linearity

#homoscedascity
BRARFKmic_func_mic12$res <- residuals(fk_mic12a)
BRARFKmic_func_mic12$abs_res <- abs(BRARFKmic_func_mic12$res)
BRARFKmic_func_mic12$abs_res2 <- BRARFKmic_func_mic12$abs_res^2
levene_brar <- lm(data = BRARFKmic_func_mic12, abs_res2 ~ BRARrel)
anova(levene_brar) #p = 0.8681 so > 0.05 so equal variance is met

#autocorrelation
acf(residuals(fk_mic12a, retype = "normalized"))  #lines dont go outside CI horizontal lines, so not autocorrelated
pacf(residuals(fk_mic12a, retype = "normalized"))


#Manganese oxid
fk_mic13 <- lmerTest::lmer(avg ~ BRARrel + (1|grad_num), data = subset(BRARFKmic_func_long2, func == "manganese_oxidation"))
summary(fk_mic13) #not sig
anova(fk_mic13)

#check normality
resid1.0 <- lm(data = subset(BRARFKmic_func_long2, func == "manganese_oxidation"), avg ~ BRARrel)
ols_plot_resid_hist(resid1.0) #right skew
ols_test_normality(resid1.0) #fail

#trans
BRARFKmic_func_mic13 <- BRARFKmic_func_long2 %>%
  filter(func == "manganese_oxidation") %>%
  mutate(ln = log10(avg + 0.001))
resid1.0 <- lm(data = BRARFKmic_func_mic13, ln ~ BRARrel)
ols_plot_resid_hist(resid1.0) #better
ols_test_normality(resid1.0) #pass
#redo stats
fk_mic13a <- lmerTest::lmer(ln ~ BRARrel + (1|grad_num), data = BRARFKmic_func_mic13)
summary(fk_mic13a) #p = 0.9885; slope = -4.935e-5
anova(fk_mic13a)

#linearity
plot(resid(fk_mic13a), BRARFKmic_func_mic13$ln) #looks linear
plot(fk_mic13a) #no pattern so indicates linearity

#homoscedascity
BRARFKmic_func_mic13$res <- residuals(fk_mic13a)
BRARFKmic_func_mic13$abs_res <- abs(BRARFKmic_func_mic13$res)
BRARFKmic_func_mic13$abs_res2 <- BRARFKmic_func_mic13$abs_res^2
levene_brar <- lm(data = BRARFKmic_func_mic13, abs_res2 ~ BRARrel)
anova(levene_brar) #p = 0.7545 so > 0.05 so equal variance is met

#autocorrelation
acf(residuals(fk_mic13a, retype = "normalized"))  #lines dont go outside CI horizontal lines, so not autocorrelated
pacf(residuals(fk_mic13a, retype = "normalized"))


#Methanol oxidation
fk_mic14 <- lmerTest::lmer(avg ~ BRARrel + (1|grad_num), data = subset(BRARFKmic_func_long2, func == "methanol_oxidation"))
summary(fk_mic14) #not sig 
anova(fk_mic14)

#check normality
resid1.0 <- lm(data = subset(BRARFKmic_func_long2, func == "methanol_oxidation"), avg ~ BRARrel)
ols_plot_resid_hist(resid1.0) #right skew
ols_test_normality(resid1.0) #fail

#trans
BRARFKmic_func_mic14 <- BRARFKmic_func_long2 %>%
  filter(func == "methanol_oxidation") %>%
  mutate(ln = log10(avg + 0.001))
resid1.0 <- lm(data = BRARFKmic_func_mic14, ln ~ BRARrel)
ols_plot_resid_hist(resid1.0) #better
ols_test_normality(resid1.0) #pass
#redo stats
fk_mic14a <- lmerTest::lmer(ln ~ BRARrel + (1|grad_num), data = BRARFKmic_func_mic14)
summary(fk_mic14a) #p = 0.2491; slope = 0.002414
anova(fk_mic14a)

#linearity
plot(resid(fk_mic14a), BRARFKmic_func_mic14$ln) #looks linear
plot(fk_mic14a) #no pattern so indicates linearity

#homoscedascity
BRARFKmic_func_mic14$res <- residuals(fk_mic14a)
BRARFKmic_func_mic14$abs_res <- abs(BRARFKmic_func_mic14$res)
BRARFKmic_func_mic14$abs_res2 <- BRARFKmic_func_mic14$abs_res^2
levene_brar <- lm(data = BRARFKmic_func_mic14, abs_res2 ~ BRARrel)
anova(levene_brar) #p = 0.6732 so > 0.05 so equal variance is met

#autocorrelation
acf(residuals(fk_mic14a, retype = "normalized"))  #lines dont go outside CI horizontal lines, so not autocorrelated
pacf(residuals(fk_mic14a, retype = "normalized"))


#Methylotrophy
fk_mic15 <- lmerTest::lmer(avg ~ BRARrel + (1|grad_num), data = subset(BRARFKmic_func_long2, func == "methylotrophy"))
summary(fk_mic15) #not sig
anova(fk_mic15)

#check normality
resid1.0 <- lm(data = subset(BRARFKmic_func_long2, func == "methylotrophy"), avg ~ BRARrel)
ols_plot_resid_hist(resid1.0) #right skew
ols_test_normality(resid1.0) #fail

#trans
BRARFKmic_func_mic15 <- BRARFKmic_func_long2 %>%
  filter(func == "methylotrophy") %>%
  mutate(ln = log10(avg + 0.001))
resid1.0 <- lm(data = BRARFKmic_func_mic15, ln ~ BRARrel)
ols_plot_resid_hist(resid1.0) #better
ols_test_normality(resid1.0) #pass
#redo stats
fk_mic15a <- lmerTest::lmer(ln ~ BRARrel + (1|grad_num), data = BRARFKmic_func_mic15)
summary(fk_mic15a) #p = 0.2491; slope = 0.002414
anova(fk_mic15a)

#linearity
plot(resid(fk_mic15a), BRARFKmic_func_mic15$ln) #looks linear
plot(fk_mic15a) #no pattern so indicates linearity

#homoscedascity
BRARFKmic_func_mic15$res <- residuals(fk_mic15a)
BRARFKmic_func_mic15$abs_res <- abs(BRARFKmic_func_mic15$res)
BRARFKmic_func_mic15$abs_res2 <- BRARFKmic_func_mic15$abs_res^2
levene_brar <- lm(data = BRARFKmic_func_mic15, abs_res2 ~ BRARrel)
anova(levene_brar) #p = 0.6732 so > 0.05 so equal variance is met

#autocorrelation
acf(residuals(fk_mic15a, retype = "normalized"))  #lines dont go outside CI horizontal lines, so not autocorrelated
pacf(residuals(fk_mic15a, retype = "normalized"))


#Nitrate denitrif
fk_mic16 <- lmerTest::lmer(avg ~ BRARrel + (1|grad_num), data = subset(BRARFKmic_func_long2, func == "nitrate_denitrification"))
summary(fk_mic16) #not sig
anova(fk_mic16)

#check normality
resid1.0 <- lm(data = subset(BRARFKmic_func_long2, func == "nitrate_denitrification"), avg ~ BRARrel)
ols_plot_resid_hist(resid1.0) #right skew
ols_test_normality(resid1.0) #fail

#trans
BRARFKmic_func_mic16 <- BRARFKmic_func_long2 %>%
  filter(func == "nitrate_denitrification") %>%
  mutate(ln = log10(avg + 0.001))
resid1.0 <- lm(data = BRARFKmic_func_mic16, ln ~ BRARrel)
ols_plot_resid_hist(resid1.0) #better
ols_test_normality(resid1.0) #pass
#redo stats
fk_mic16a <- lmerTest::lmer(ln ~ BRARrel + (1|grad_num), data = BRARFKmic_func_mic16)
summary(fk_mic16a) #p = 0.756; slope = 7.198e-4
anova(fk_mic16a)

#linearity
plot(resid(fk_mic16a), BRARFKmic_func_mic16$ln) #looks linear
plot(fk_mic16a) #no pattern so indicates linearity

#homoscedascity
BRARFKmic_func_mic16$res <- residuals(fk_mic16a)
BRARFKmic_func_mic16$abs_res <- abs(BRARFKmic_func_mic16$res)
BRARFKmic_func_mic16$abs_res2 <- BRARFKmic_func_mic16$abs_res^2
levene_brar <- lm(data = BRARFKmic_func_mic16, abs_res2 ~ BRARrel)
anova(levene_brar) #p = 0.6732 so > 0.05 so equal variance is met

#autocorrelation
acf(residuals(fk_mic16a, retype = "normalized"))  #lines dont go outside CI horizontal lines, so not autocorrelated
pacf(residuals(fk_mic16a, retype = "normalized"))


#Nitrate reduction
fk_mic17 <- lmerTest::lmer(avg ~ BRARrel + (1|grad_num), data = subset(BRARFKmic_func_long2, func == "nitrate_reduction"))
summary(fk_mic17) #not sig
anova(fk_mic17)

#check normality
resid1.0 <- lm(data = subset(BRARFKmic_func_long2, func == "nitrate_reduction"), avg ~ BRARrel)
ols_plot_resid_hist(resid1.0) #right skew
ols_test_normality(resid1.0) #fail

#trans
BRARFKmic_func_mic17 <- BRARFKmic_func_long2 %>%
  filter(func == "nitrate_reduction") %>%
  mutate(ln = log10(avg + 0.001))
resid1.0 <- lm(data = BRARFKmic_func_mic17, ln ~ BRARrel)
ols_plot_resid_hist(resid1.0) #better
ols_test_normality(resid1.0) #pass
#redo stats
fk_mic17a <- lmerTest::lmer(ln ~ BRARrel + (1|grad_num), data = BRARFKmic_func_mic17)
summary(fk_mic17a) #p = 0.779; slope = -0.0005477
anova(fk_mic17a)

#linearity
plot(resid(fk_mic17a), BRARFKmic_func_mic17$ln) #looks linear
plot(fk_mic17a) #no pattern so indicates linearity

#homoscedascity
BRARFKmic_func_mic17$res <- residuals(fk_mic17a)
BRARFKmic_func_mic17$abs_res <- abs(BRARFKmic_func_mic17$res)
BRARFKmic_func_mic17$abs_res2 <- BRARFKmic_func_mic17$abs_res^2
levene_brar <- lm(data = BRARFKmic_func_mic17, abs_res2 ~ BRARrel)
anova(levene_brar) #p = 0.6732 so > 0.05 so equal variance is met

#autocorrelation
acf(residuals(fk_mic17a, retype = "normalized"))  #lines dont go outside CI horizontal lines, so not autocorrelated
pacf(residuals(fk_mic17a, retype = "normalized"))


#Nitrate respiration
fk_mic18 <- lmerTest::lmer(avg ~ BRARrel + (1|grad_num), data = subset(BRARFKmic_func_long2, func == "nitrate_respiration"))
summary(fk_mic18) #not sig 
anova(fk_mic18)

#check normality
resid1.0 <- lm(data = subset(BRARFKmic_func_long2, func == "nitrate_respiration"), avg ~ BRARrel)
ols_plot_resid_hist(resid1.0) #right skew
ols_test_normality(resid1.0) #fail

#trans
BRARFKmic_func_mic18 <- BRARFKmic_func_long2 %>%
  filter(func == "nitrate_respiration") %>%
  mutate(ln = log10(avg + 0.001))
resid1.0 <- lm(data = BRARFKmic_func_mic18, ln ~ BRARrel)
ols_plot_resid_hist(resid1.0) #better
ols_test_normality(resid1.0) #pass
#redo stats
fk_mic18a <- lmerTest::lmer(ln ~ BRARrel + (1|grad_num), data = BRARFKmic_func_mic18)
summary(fk_mic18a) #p = 0.7627; slope = 7.717e-4
anova(fk_mic18a)

#linearity
plot(resid(fk_mic18a), BRARFKmic_func_mic18$ln) #looks linear
plot(fk_mic18a) #no pattern so indicates linearity

#homoscedascity
BRARFKmic_func_mic18$res <- residuals(fk_mic18a)
BRARFKmic_func_mic18$abs_res <- abs(BRARFKmic_func_mic18$res)
BRARFKmic_func_mic18$abs_res2 <- BRARFKmic_func_mic18$abs_res^2
levene_brar <- lm(data = BRARFKmic_func_mic18, abs_res2 ~ BRARrel)
anova(levene_brar) #p = 0.9956 so > 0.05 so equal variance is met

#autocorrelation
acf(residuals(fk_mic18a, retype = "normalized"))  #lines dont go outside CI horizontal lines, so not autocorrelated
pacf(residuals(fk_mic18a, retype = "normalized"))


#Nitrificaiton
fk_mic19 <- lmerTest::lmer(avg ~ BRARrel + (1|grad_num), data = subset(BRARFKmic_func_long2, func == "nitrification"))
summary(fk_mic19) #not sig
anova(fk_mic19)

#check normality
resid1.0 <- lm(data = subset(BRARFKmic_func_long2, func == "nitrification"), avg ~ BRARrel)
ols_plot_resid_hist(resid1.0) #right skew
ols_test_normality(resid1.0) #fail

#trans
BRARFKmic_func_mic19 <- BRARFKmic_func_long2 %>%
  filter(func == "nitrification") %>%
  mutate(ln = log10(avg + 0.001))
resid1.0 <- lm(data = BRARFKmic_func_mic19, ln ~ BRARrel)
ols_plot_resid_hist(resid1.0) #better
ols_test_normality(resid1.0) #pass
#redo stats
fk_mic19a <- lmerTest::lmer(ln ~ BRARrel + (1|grad_num), data = BRARFKmic_func_mic19)
summary(fk_mic19a) #p = 0.865; slope = -0.0003264
anova(fk_mic19a)

#linearity
plot(resid(fk_mic19a), BRARFKmic_func_mic19$ln) #looks linear
plot(fk_mic19a) #no pattern so indicates linearity

#homoscedascity
BRARFKmic_func_mic19$res <- residuals(fk_mic19a)
BRARFKmic_func_mic19$abs_res <- abs(BRARFKmic_func_mic19$res)
BRARFKmic_func_mic19$abs_res2 <- BRARFKmic_func_mic19$abs_res^2
levene_brar <- lm(data = BRARFKmic_func_mic19, abs_res2 ~ BRARrel)
anova(levene_brar) #p = 0.6656 so > 0.05 so equal variance is met

#autocorrelation
acf(residuals(fk_mic19a, retype = "normalized"))  #lines dont go outside CI horizontal lines, so not autocorrelated
pacf(residuals(fk_mic19a, retype = "normalized"))


#Nitrite denitrif
fk_mic20 <- lmerTest::lmer(avg ~ BRARrel + (1|grad_num), data = subset(BRARFKmic_func_long2, func == "nitrite_denitrification"))
summary(fk_mic20) #not sig
anova(fk_mic20)

#check normality
resid1.0 <- lm(data = subset(BRARFKmic_func_long2, func == "nitrite_denitrification"), avg ~ BRARrel)
ols_plot_resid_hist(resid1.0) #right skew
ols_test_normality(resid1.0) #fail

#trans
BRARFKmic_func_mic20 <- BRARFKmic_func_long2 %>%
  filter(func == "nitrite_denitrification") %>%
  mutate(ln = log10(avg + 0.001))
resid1.0 <- lm(data = BRARFKmic_func_mic20, ln ~ BRARrel)
ols_plot_resid_hist(resid1.0) #better
ols_test_normality(resid1.0) #pass
#redo stats
fk_mic20a <- lmerTest::lmer(ln ~ BRARrel + (1|grad_num), data = BRARFKmic_func_mic20)
summary(fk_mic20a) #p = 0.756; slope = 7.918e-4
anova(fk_mic20a)


#Nitrite resp
fk_mic21 <- lmerTest::lmer(avg ~ BRARrel + (1|grad_num), data = subset(BRARFKmic_func_long2, func == "nitrite_respiration"))
summary(fk_mic21) #not sig
anova(fk_mic21)

#check normality
resid1.0 <- lm(data = subset(BRARFKmic_func_long2, func == "nitrite_respiration"), avg ~ BRARrel)
ols_plot_resid_hist(resid1.0) #right skew
ols_test_normality(resid1.0) #fail

#trans
BRARFKmic_func_mic21 <- BRARFKmic_func_long2 %>%
  filter(func == "nitrite_respiration") %>%
  mutate(ln = log10(avg + 0.001))
resid1.0 <- lm(data = BRARFKmic_func_mic21, ln ~ BRARrel)
ols_plot_resid_hist(resid1.0) #better
ols_test_normality(resid1.0) #pass
#redo stats
fk_mic21a <- lmerTest::lmer(ln ~ BRARrel + (1|grad_num), data = BRARFKmic_func_mic21)
summary(fk_mic21a) #p = 0.7558; slope = 7.918e-4
anova(fk_mic21a)


#Nitrogen fixation
fk_mic22 <- lmerTest::lmer(avg ~ BRARrel + (1|grad_num), data = subset(BRARFKmic_func_long2, func == "nitrogen_fixation"))
summary(fk_mic22) #not sig
anova(fk_mic22)

#check normality
resid1.0 <- lm(data = subset(BRARFKmic_func_long2, func == "nitrogen_fixation"), avg ~ BRARrel)
ols_plot_resid_hist(resid1.0) #right skew
ols_test_normality(resid1.0) #fail

#trans
BRARFKmic_func_mic22 <- BRARFKmic_func_long2 %>%
  filter(func == "nitrogen_fixation") %>%
  mutate(ln = log10(avg + 0.001))
resid1.0 <- lm(data = BRARFKmic_func_mic22, ln ~ BRARrel)
ols_plot_resid_hist(resid1.0) #better
ols_test_normality(resid1.0) #pass
#redo stats
fk_mic22a <- lmerTest::lmer(ln ~ BRARrel + (1|grad_num), data = BRARFKmic_func_mic22)
summary(fk_mic22a) #p = 0.7344; slope = 5.818e-4
anova(fk_mic22a)

#linearity
plot(resid(fk_mic22a), BRARFKmic_func_mic22$ln) #looks linear
plot(fk_mic22a) #no pattern so indicates linearity

#homoscedascity
BRARFKmic_func_mic22$res <- residuals(fk_mic22a)
BRARFKmic_func_mic22$abs_res <- abs(BRARFKmic_func_mic22$res)
BRARFKmic_func_mic22$abs_res2 <- BRARFKmic_func_mic22$abs_res^2
levene_brar <- lm(data = BRARFKmic_func_mic22, abs_res2 ~ BRARrel)
anova(levene_brar) #p = 0.9417 so > 0.05 so equal variance is met

#autocorrelation
acf(residuals(fk_mic22a, retype = "normalized"))  #lines dont go outside CI horizontal lines, so not autocorrelated
pacf(residuals(fk_mic22a, retype = "normalized"))


#Nitrogen respiration
fk_mic23 <- lmerTest::lmer(avg ~ BRARrel + (1|grad_num), data = subset(BRARFKmic_func_long2, func == "nitrogen_respiration"))
summary(fk_mic23) #not sig
anova(fk_mic23)


#Nitrous oxide denit
fk_mic24 <- lmerTest::lmer(avg ~ BRARrel + (1|grad_num), data = subset(BRARFKmic_func_long2, func == "nitrous_oxide_denitrification"))
summary(fk_mic24) #not sig
anova(fk_mic24)


#Nonphoto cyano
fk_mic25 <- lmerTest::lmer(avg ~ BRARrel + (1|grad_num), data = subset(BRARFKmic_func_long2, func == "nonphotosynthetic_cyanobacteria"))
summary(fk_mic25) #not sig
anova(fk_mic25)

#check normality
resid1.0 <- lm(data = subset(BRARFKmic_func_long2, func == "nonphotosynthetic_cyanobacteria"), avg ~ BRARrel)
ols_plot_resid_hist(resid1.0) #right skew
ols_test_normality(resid1.0) #fail

#trans
BRARFKmic_func_mic25 <- BRARFKmic_func_long2 %>%
  filter(func == "nonphotosynthetic_cyanobacteria") %>%
  mutate(ln = log10(avg + 0.001))
resid1.0 <- lm(data = BRARFKmic_func_mic25, ln ~ BRARrel)
ols_plot_resid_hist(resid1.0) #better
ols_test_normality(resid1.0) #pass
#redo stats
fk_mic25a <- lmerTest::lmer(ln ~ BRARrel + (1|grad_num), data = BRARFKmic_func_mic25)
summary(fk_mic25a) #p = 0.8576; slope = -0.0005561
anova(fk_mic25a)

#linearity
plot(resid(fk_mic25a), BRARFKmic_func_mic25$ln) #looks linear
plot(fk_mic25a) #no pattern so indicates linearity

#homoscedascity
BRARFKmic_func_mic25$res <- residuals(fk_mic25a)
BRARFKmic_func_mic25$abs_res <- abs(BRARFKmic_func_mic25$res)
BRARFKmic_func_mic25$abs_res2 <- BRARFKmic_func_mic25$abs_res^2
levene_brar <- lm(data = BRARFKmic_func_mic25, abs_res2 ~ BRARrel)
anova(levene_brar) #p = 0.7691 so > 0.05 so equal variance is met

#autocorrelation
acf(residuals(fk_mic25a, retype = "normalized"))  #lines dont go outside CI horizontal lines, so not autocorrelated
pacf(residuals(fk_mic25a, retype = "normalized"))


#Other
fk_mic26 <- lmerTest::lmer(avg ~ BRARrel + (1|grad_num), data = subset(BRARFKmic_func_long2, func == "other"))
summary(fk_mic26) #not sig
anova(fk_mic26)

#check normality
resid1.0 <- lm(data = subset(BRARFKmic_func_long2, func == "other"), avg ~ BRARrel)
ols_plot_resid_hist(resid1.0) #right skew
ols_test_normality(resid1.0) #fail

#trans
BRARFKmic_func_mic26 <- BRARFKmic_func_long2 %>%
  filter(func == "other") %>%
  mutate(ln = log10(avg + 0.001))
resid1.0 <- lm(data = BRARFKmic_func_mic26, ln ~ BRARrel)
ols_plot_resid_hist(resid1.0) #better
ols_test_normality(resid1.0) #pass
#redo stats
fk_mic26a <- lmerTest::lmer(ln ~ BRARrel + (1|grad_num), data = BRARFKmic_func_mic26)
summary(fk_mic26a) #p = 0.8753; slope = -0.000292
anova(fk_mic26a)

#linearity
plot(resid(fk_mic26a), BRARFKmic_func_mic26$ln) #looks linear
plot(fk_mic26a) #no pattern so indicates linearity

#homoscedascity
BRARFKmic_func_mic26$res <- residuals(fk_mic26a)
BRARFKmic_func_mic26$abs_res <- abs(BRARFKmic_func_mic26$res)
BRARFKmic_func_mic26$abs_res2 <- BRARFKmic_func_mic26$abs_res^2
levene_brar <- lm(data = BRARFKmic_func_mic26, abs_res2 ~ BRARrel)
anova(levene_brar) #p = 0.6732 so > 0.05 so equal variance is met

#autocorrelation
acf(residuals(fk_mic26a, retype = "normalized"))  #lines dont go outside CI horizontal lines, so not autocorrelated
pacf(residuals(fk_mic26a, retype = "normalized"))


#Oxygenic photo
fk_mic27 <- lmerTest::lmer(avg ~ BRARrel + (1|grad_num), data = subset(BRARFKmic_func_long2, func == "oxygenic_photoautotrophy"))
summary(fk_mic27) #p = 0.064801; slope = -1.1677
anova(fk_mic27)

#check normality
resid1.0 <- lm(data = subset(BRARFKmic_func_long2, func == "oxygenic_photoautotrophy"), avg ~ BRARrel)
ols_plot_resid_hist(resid1.0) #right skew
ols_test_normality(resid1.0) #fail

#trans
BRARFKmic_func_mic27 <- BRARFKmic_func_long2 %>%
  filter(func == "oxygenic_photoautotrophy") %>%
  mutate(ln = log10(avg + 0.001))
resid1.0 <- lm(data = BRARFKmic_func_mic27, ln ~ BRARrel)
ols_plot_resid_hist(resid1.0) #better
ols_test_normality(resid1.0) #pass
#redo stats
fk_mic27a <- lmerTest::lmer(ln ~ BRARrel + (1|grad_num), data = BRARFKmic_func_mic27)
summary(fk_mic27a) #p = 0.07983; slope = -0.00754
anova(fk_mic27a)

#linearity
plot(resid(fk_mic27a), BRARFKmic_func_mic27$ln) #looks linear
plot(fk_mic27a) #no pattern so indicates linearity

#homoscedascity
BRARFKmic_func_mic27$res <- residuals(fk_mic27a)
BRARFKmic_func_mic27$abs_res <- abs(BRARFKmic_func_mic27$res)
BRARFKmic_func_mic27$abs_res2 <- BRARFKmic_func_mic27$abs_res^2
levene_brar <- lm(data = BRARFKmic_func_mic27, abs_res2 ~ BRARrel)
anova(levene_brar) #p = 0.4031 so > 0.05 so equal variance is met

#autocorrelation
acf(residuals(fk_mic27a, retype = "normalized"))  #lines dont go outside CI horizontal lines, so not autocorrelated
pacf(residuals(fk_mic27a, retype = "normalized"))


#Photoauto
fk_mic28 <- lmerTest::lmer(avg ~ BRARrel + (1|grad_num), data = subset(BRARFKmic_func_long2, func == "photoautotrophy"))
summary(fk_mic28) #not sig
anova(fk_mic28)

#check normality
resid1.0 <- lm(data = subset(BRARFKmic_func_long2, func == "photoautotrophy"), avg ~ BRARrel)
ols_plot_resid_hist(resid1.0) #right skew
ols_test_normality(resid1.0) #fail

#trans
BRARFKmic_func_mic28 <- BRARFKmic_func_long2 %>%
  filter(func == "photoautotrophy") %>%
  mutate(ln = log10(avg + 0.001))
resid1.0 <- lm(data = BRARFKmic_func_mic28, ln ~ BRARrel)
ols_plot_resid_hist(resid1.0) #better
ols_test_normality(resid1.0) #pass
#redo stats
fk_mic28a <- lmerTest::lmer(ln ~ BRARrel + (1|grad_num), data = BRARFKmic_func_mic28)
summary(fk_mic28a) #p = 0.8033; slope = -0.0005558
anova(fk_mic28a)

#linearity
plot(resid(fk_mic28a), BRARFKmic_func_mic28$ln) #looks linear
plot(fk_mic28a) #no pattern so indicates linearity

#homoscedascity
BRARFKmic_func_mic28$res <- residuals(fk_mic28a)
BRARFKmic_func_mic28$abs_res <- abs(BRARFKmic_func_mic28$res)
BRARFKmic_func_mic28$abs_res2 <- BRARFKmic_func_mic28$abs_res^2
levene_brar <- lm(data = BRARFKmic_func_mic28, abs_res2 ~ BRARrel)
anova(levene_brar) #p = 0.7277 so > 0.05 so equal variance is met

#autocorrelation
acf(residuals(fk_mic28a, retype = "normalized"))  #lines dont go outside CI horizontal lines, so not autocorrelated
pacf(residuals(fk_mic28a, retype = "normalized"))


#Photoheterotrophy
fk_mic29 <- lmerTest::lmer(avg ~ BRARrel + (1|grad_num), data = subset(BRARFKmic_func_long2, func == "photoheterotrophy"))
summary(fk_mic29) #not sig
anova(fk_mic29)

#check normality
resid1.0 <- lm(data = subset(BRARFKmic_func_long2, func == "photoheterotrophy"), avg ~ BRARrel)
ols_plot_resid_hist(resid1.0) #right skew
ols_test_normality(resid1.0) #fail

#trans
BRARFKmic_func_mic29 <- BRARFKmic_func_long2 %>%
  filter(func == "photoheterotrophy") %>%
  mutate(ln = log10(avg + 0.001))
resid1.0 <- lm(data = BRARFKmic_func_mic29, ln ~ BRARrel)
ols_plot_resid_hist(resid1.0) #better
ols_test_normality(resid1.0) #pass
#redo stats
fk_mic29a <- lmerTest::lmer(ln ~ BRARrel + (1|grad_num), data = BRARFKmic_func_mic29)
summary(fk_mic29a) #p = 0.8194; slope = 5.183e-4
anova(fk_mic29a)

#linearity
plot(resid(fk_mic29a), BRARFKmic_func_mic29$ln) #looks linear
plot(fk_mic29a) #no pattern so indicates linearity

#homoscedascity
BRARFKmic_func_mic29$res <- residuals(fk_mic29a)
BRARFKmic_func_mic29$abs_res <- abs(BRARFKmic_func_mic29$res)
BRARFKmic_func_mic29$abs_res2 <- BRARFKmic_func_mic29$abs_res^2
levene_brar <- lm(data = BRARFKmic_func_mic29, abs_res2 ~ BRARrel)
anova(levene_brar) #p = 0.861 so > 0.05 so equal variance is met

#autocorrelation
acf(residuals(fk_mic29a, retype = "normalized"))  #lines dont go outside CI horizontal lines, so not autocorrelated
pacf(residuals(fk_mic29a, retype = "normalized"))


#Photo cyano
fk_mic30 <- lmerTest::lmer(avg ~ BRARrel + (1|grad_num), data = subset(BRARFKmic_func_long2, func == "photosynthetic_cyanobacteria"))
summary(fk_mic30) #p = 00.064801; slope = -1.1677
anova(fk_mic30)



#Phototrophy
fk_mic31 <- lmerTest::lmer(avg ~ BRARrel + (1|grad_num), data = subset(BRARFKmic_func_long2, func == "phototrophy"))
summary(fk_mic31) #not sig
anova(fk_mic31)

#check normality
resid1.0 <- lm(data = subset(BRARFKmic_func_long2, func == "phototrophy"), avg ~ BRARrel)
ols_plot_resid_hist(resid1.0) #right skew
ols_test_normality(resid1.0) #fail

#trans
BRARFKmic_func_mic31 <- BRARFKmic_func_long2 %>%
  filter(func == "phototrophy") %>%
  mutate(ln = log10(avg + 0.001))
resid1.0 <- lm(data = BRARFKmic_func_mic31, ln ~ BRARrel)
ols_plot_resid_hist(resid1.0) #better
ols_test_normality(resid1.0) #pass
#redo stats
fk_mic31a <- lmerTest::lmer(ln ~ BRARrel + (1|grad_num), data = BRARFKmic_func_mic31)
summary(fk_mic31a) #p = 0.8567; slope = -0.000375
anova(fk_mic31a)

#linearity
plot(resid(fk_mic31a), BRARFKmic_func_mic31$ln) #looks linear
plot(fk_mic31a) #no pattern so indicates linearity

#homoscedascity
BRARFKmic_func_mic31$res <- residuals(fk_mic31a)
BRARFKmic_func_mic31$abs_res <- abs(BRARFKmic_func_mic31$res)
BRARFKmic_func_mic31$abs_res2 <- BRARFKmic_func_mic31$abs_res^2
levene_brar <- lm(data = BRARFKmic_func_mic31, abs_res2 ~ BRARrel)
anova(levene_brar) #p = 0.7515 so > 0.05 so equal variance is met

#autocorrelation
acf(residuals(fk_mic31a, retype = "normalized"))  #lines dont go outside CI horizontal lines, so not autocorrelated
pacf(residuals(fk_mic31a, retype = "normalized"))


#Pred or exo
fk_mic32 <- lmerTest::lmer(avg ~ BRARrel + (1|grad_num), data = subset(BRARFKmic_func_long2, func == "predatory_or_exoparasitic"))
summary(fk_mic32) #not sig 
anova(fk_mic32)

#check normality
resid1.0 <- lm(data = subset(BRARFKmic_func_long2, func == "predatory_or_exoparasitic"), avg ~ BRARrel)
ols_plot_resid_hist(resid1.0) #right skew
ols_test_normality(resid1.0) #fail

#trans
BRARFKmic_func_mic32 <- BRARFKmic_func_long2 %>%
  filter(func == "predatory_or_exoparasitic") %>%
  mutate(ln = log10(avg + 0.001))
resid1.0 <- lm(data = BRARFKmic_func_mic32, ln ~ BRARrel)
ols_plot_resid_hist(resid1.0) #better
ols_test_normality(resid1.0) #pass
#redo stats
fk_mic32a <- lmerTest::lmer(ln ~ BRARrel + (1|grad_num), data = BRARFKmic_func_mic32)
summary(fk_mic32a) #p = 0.6024; slope = -0.001127
anova(fk_mic32a)

#linearity
plot(resid(fk_mic32a), BRARFKmic_func_mic32$ln) #looks linear
plot(fk_mic32a) #no pattern so indicates linearity

#homoscedascity
BRARFKmic_func_mic32$res <- residuals(fk_mic32a)
BRARFKmic_func_mic32$abs_res <- abs(BRARFKmic_func_mic32$res)
BRARFKmic_func_mic32$abs_res2 <- BRARFKmic_func_mic32$abs_res^2
levene_brar <- lm(data = BRARFKmic_func_mic32, abs_res2 ~ BRARrel)
anova(levene_brar) #p = 0.7515 so > 0.05 so equal variance is met

#autocorrelation
acf(residuals(fk_mic32a, retype = "normalized"))  #lines dont go outside CI horizontal lines, so not autocorrelated
pacf(residuals(fk_mic32a, retype = "normalized"))


#Resp sulf
fk_mic33 <- lmerTest::lmer(avg ~ BRARrel + (1|grad_num), data = subset(BRARFKmic_func_long2, func == "respiration_of_sulfur_compounds"))
summary(fk_mic33) #not sig
anova(fk_mic33)

#check normality
resid1.0 <- lm(data = subset(BRARFKmic_func_long2, func == "respiration_of_sulfur_compounds"), avg ~ BRARrel)
ols_plot_resid_hist(resid1.0) #right skew
ols_test_normality(resid1.0) #fail

#trans
BRARFKmic_func_mic33 <- BRARFKmic_func_long2 %>%
  filter(func == "respiration_of_sulfur_compounds") %>%
  mutate(ln = log10(avg + 0.001))
resid1.0 <- lm(data = BRARFKmic_func_mic33, ln ~ BRARrel)
ols_plot_resid_hist(resid1.0) #better
ols_test_normality(resid1.0) #pass
#redo stats
fk_mic33a <- lmerTest::lmer(ln ~ BRARrel + (1|grad_num), data = BRARFKmic_func_mic33)
summary(fk_mic33a) #p = 0.6024; slope = -0.001127
anova(fk_mic33a)

#linearity
plot(resid(fk_mic33a), BRARFKmic_func_mic33$ln) #looks linear
plot(fk_mic33a) #no pattern so indicates linearity

#homoscedascity
BRARFKmic_func_mic33$res <- residuals(fk_mic33a)
BRARFKmic_func_mic33$abs_res <- abs(BRARFKmic_func_mic33$res)
BRARFKmic_func_mic33$abs_res2 <- BRARFKmic_func_mic33$abs_res^2
levene_brar <- lm(data = BRARFKmic_func_mic33, abs_res2 ~ BRARrel)
anova(levene_brar) #p = 0.6111 so > 0.05 so equal variance is met

#autocorrelation
acf(residuals(fk_mic33a, retype = "normalized"))  #lines dont go outside CI horizontal lines, so not autocorrelated
pacf(residuals(fk_mic33a, retype = "normalized"))


#Sulfate resp
fk_mic34 <- lmerTest::lmer(avg ~ BRARrel + (1|grad_num), data = subset(BRARFKmic_func_long2, func == "sulfate_respiration"))
summary(fk_mic34) #not sig
anova(fk_mic34)


#Ureolysis
fk_mic35 <- lmerTest::lmer(avg ~ BRARrel + (1|grad_num), data = subset(BRARFKmic_func_long2, func == "ureolysis"))
summary(fk_mic35) #not sig
anova(fk_mic35)

#check normality
resid1.0 <- lm(data = subset(BRARFKmic_func_long2, func == "ureolysis"), avg ~ BRARrel)
ols_plot_resid_hist(resid1.0) #right skew
ols_test_normality(resid1.0) #fail

#trans
BRARFKmic_func_mic35 <- BRARFKmic_func_long2 %>%
  filter(func == "ureolysis") %>%
  mutate(ln = log10(avg + 0.001))
resid1.0 <- lm(data = BRARFKmic_func_mic35, ln ~ BRARrel)
ols_plot_resid_hist(resid1.0) #better
ols_test_normality(resid1.0) #pass
#redo stats
fk_mic35a <- lmerTest::lmer(ln ~ BRARrel + (1|grad_num), data = BRARFKmic_func_mic35)
summary(fk_mic35a) #p = 0.3727; slope = 0.001641
anova(fk_mic35a)

#linearity
plot(resid(fk_mic35a), BRARFKmic_func_mic35$ln) #looks linear
plot(fk_mic35a) #no pattern so indicates linearity

#homoscedascity
BRARFKmic_func_mic35$res <- residuals(fk_mic35a)
BRARFKmic_func_mic35$abs_res <- abs(BRARFKmic_func_mic35$res)
BRARFKmic_func_mic35$abs_res2 <- BRARFKmic_func_mic35$abs_res^2
levene_brar <- lm(data = BRARFKmic_func_mic35, abs_res2 ~ BRARrel)
anova(levene_brar) #p = 0.823 so > 0.05 so equal variance is met

#autocorrelation
acf(residuals(fk_mic35a, retype = "normalized"))  #lines dont go outside CI horizontal lines, so not autocorrelated
pacf(residuals(fk_mic35a, retype = "normalized"))


#Cellulo
fk_mic36 <- lmerTest::lmer(avg ~ BRARrel + (1|grad_num), data = subset(BRARFKmic_func_long2, func == "cellulolysis"))
summary(fk_mic36) #not sig
anova(fk_mic36)

#check normality
resid1.0 <- lm(data = subset(BRARFKmic_func_long2, func == "cellulolysis"), avg ~ BRARrel)
ols_plot_resid_hist(resid1.0) #right skew
ols_test_normality(resid1.0) #fail

#trans
BRARFKmic_func_mic36 <- BRARFKmic_func_long2 %>%
  filter(func == "cellulolysis") %>%
  mutate(ln = log10(avg + 0.001))
resid1.0 <- lm(data = BRARFKmic_func_mic36, ln ~ BRARrel)
ols_plot_resid_hist(resid1.0) #better
ols_test_normality(resid1.0) #pass
#redo stats
fk_mic36a <- lmerTest::lmer(ln ~ BRARrel + (1|grad_num), data = BRARFKmic_func_mic36)
summary(fk_mic36a) #p = 0.6941; slope = 0.001489
anova(fk_mic36a)

#linearity
plot(resid(fk_mic36a), BRARFKmic_func_mic36$ln) #looks linear
plot(fk_mic36a) #no pattern so indicates linearity

#homoscedascity
BRARFKmic_func_mic36$res <- residuals(fk_mic36a)
BRARFKmic_func_mic36$abs_res <- abs(BRARFKmic_func_mic36$res)
BRARFKmic_func_mic36$abs_res2 <- BRARFKmic_func_mic36$abs_res^2
levene_brar <- lm(data = BRARFKmic_func_mic36, abs_res2 ~ BRARrel)
anova(levene_brar) #p = 0.823 so > 0.05 so equal variance is met

#autocorrelation
acf(residuals(fk_mic36a, retype = "normalized"))  #lines dont go outside CI horizontal lines, so not autocorrelated
pacf(residuals(fk_mic36a, retype = "normalized"))


#Human associated
fk_mic37 <- lmerTest::lmer(avg ~ BRARrel + (1|grad_num), data = subset(BRARFKmic_func_long2, func == "human_associated"))
summary(fk_mic37) #not sig
anova(fk_mic37)

#check normality
resid1.0 <- lm(data = subset(BRARFKmic_func_long2, func == "human_associated"), avg ~ BRARrel)
ols_plot_resid_hist(resid1.0) #right skew
ols_test_normality(resid1.0) #fail

#trans
BRARFKmic_func_mic37 <- BRARFKmic_func_long2 %>%
  filter(func == "human_associated") %>%
  mutate(ln = log10(avg + 0.001))
resid1.0 <- lm(data = BRARFKmic_func_mic37, ln ~ BRARrel)
ols_plot_resid_hist(resid1.0) #better
ols_test_normality(resid1.0) #pass
#redo stats
fk_mic37a <- lmerTest::lmer(ln ~ BRARrel + (1|grad_num), data = BRARFKmic_func_mic37)
summary(fk_mic37a) #p = 0.7019; slope = 0.001937
anova(fk_mic37a)

#linearity
plot(resid(fk_mic37a), BRARFKmic_func_mic37$ln) #looks linear
plot(fk_mic37a) #no pattern so indicates linearity

#homoscedascity
BRARFKmic_func_mic37$res <- residuals(fk_mic37a)
BRARFKmic_func_mic37$abs_res <- abs(BRARFKmic_func_mic37$res)
BRARFKmic_func_mic37$abs_res2 <- BRARFKmic_func_mic37$abs_res^2
levene_brar <- lm(data = BRARFKmic_func_mic37, abs_res2 ~ BRARrel)
anova(levene_brar) #p = 0.8756 so > 0.05 so equal variance is met

#autocorrelation
acf(residuals(fk_mic37a, retype = "normalized"))  #lines dont go outside CI horizontal lines, so not autocorrelated
pacf(residuals(fk_mic37a, retype = "normalized"))


#Human pathogens all
fk_mic38 <- lmerTest::lmer(avg ~ BRARrel + (1|grad_num), data = subset(BRARFKmic_func_long2, func == "human_pathogens_all"))
summary(fk_mic38) #not sig
anova(fk_mic38)

#check normality
resid1.0 <- lm(data = subset(BRARFKmic_func_long2, func == "human_pathogens_all"), avg ~ BRARrel)
ols_plot_resid_hist(resid1.0) #right skew
ols_test_normality(resid1.0) #fail

#trans
BRARFKmic_func_mic38 <- BRARFKmic_func_long2 %>%
  filter(func == "human_pathogens_all") %>%
  mutate(ln = log10(avg + 0.001))
resid1.0 <- lm(data = BRARFKmic_func_mic38, ln ~ BRARrel)
ols_plot_resid_hist(resid1.0) #better
ols_test_normality(resid1.0) #pass
#redo stats
fk_mic38a <- lmerTest::lmer(ln ~ BRARrel + (1|grad_num), data = BRARFKmic_func_mic38)
summary(fk_mic38a) #p = 0.2993; slope = 0.005339
anova(fk_mic38a)

#linearity
plot(resid(fk_mic38a), BRARFKmic_func_mic38$ln) #looks linear
plot(fk_mic38a) #no pattern so indicates linearity

#homoscedascity
BRARFKmic_func_mic38$res <- residuals(fk_mic38a)
BRARFKmic_func_mic38$abs_res <- abs(BRARFKmic_func_mic38$res)
BRARFKmic_func_mic38$abs_res2 <- BRARFKmic_func_mic38$abs_res^2
levene_brar <- lm(data = BRARFKmic_func_mic38, abs_res2 ~ BRARrel)
anova(levene_brar) #p = 0.8756 so > 0.05 so equal variance is met

#autocorrelation
acf(residuals(fk_mic38a, retype = "normalized"))  #lines dont go outside CI horizontal lines, so not autocorrelated
pacf(residuals(fk_mic38a, retype = "normalized"))


#Iron resp
fk_mic39 <- lmerTest::lmer(avg ~ BRARrel + (1|grad_num), data = subset(BRARFKmic_func_long2, func == "iron_respiration"))
summary(fk_mic39) #not sig
anova(fk_mic39)

#check normality
resid1.0 <- lm(data = subset(BRARFKmic_func_long2, func == "iron_respiration"), avg ~ BRARrel)
ols_plot_resid_hist(resid1.0) #right skew
ols_test_normality(resid1.0) #fail

#trans
BRARFKmic_func_mic39 <- BRARFKmic_func_long2 %>%
  filter(func == "iron_respiration") %>%
  mutate(ln = log10(avg + 0.001))
resid1.0 <- lm(data = BRARFKmic_func_mic39, ln ~ BRARrel)
ols_plot_resid_hist(resid1.0) #better
ols_test_normality(resid1.0) #pass
#redo stats
fk_mic39a <- lmerTest::lmer(ln ~ BRARrel + (1|grad_num), data = BRARFKmic_func_mic39)
summary(fk_mic39a) #p = 0.5517; slope = 0.0051253
anova(fk_mic39a)

#linearity
plot(resid(fk_mic39a), BRARFKmic_func_mic39$ln) #looks linear
plot(fk_mic39a) #no pattern so indicates linearity

#homoscedascity
BRARFKmic_func_mic39$res <- residuals(fk_mic39a)
BRARFKmic_func_mic39$abs_res <- abs(BRARFKmic_func_mic39$res)
BRARFKmic_func_mic39$abs_res2 <- BRARFKmic_func_mic39$abs_res^2
levene_brar <- lm(data = BRARFKmic_func_mic39, abs_res2 ~ BRARrel)
anova(levene_brar) #p = 0.8756 so > 0.05 so equal variance is met

#autocorrelation
acf(residuals(fk_mic39a, retype = "normalized"))  #lines dont go outside CI horizontal lines, so not autocorrelated
pacf(residuals(fk_mic39a, retype = "normalized"))


#Intra para
fk_mic40 <- lmerTest::lmer(avg ~ BRARrel + (1|grad_num), data = subset(BRARFKmic_func_long2, func == "intracellular_parasites"))
summary(fk_mic40) #not sig
anova(fk_mic40)

#check normality
resid1.0 <- lm(data = subset(BRARFKmic_func_long2, func == "intracellular_parasites"), avg ~ BRARrel)
ols_plot_resid_hist(resid1.0) #right skew
ols_test_normality(resid1.0) #fail

#trans
BRARFKmic_func_mic40 <- BRARFKmic_func_long2 %>%
  filter(func == "intracellular_parasites") %>%
  mutate(ln = log10(avg + 0.001))
resid1.0 <- lm(data = BRARFKmic_func_mic40, ln ~ BRARrel)
ols_plot_resid_hist(resid1.0) #better
ols_test_normality(resid1.0) #pass
#redo stats
fk_mic40a <- lmerTest::lmer(ln ~ BRARrel + (1|grad_num), data = BRARFKmic_func_mic40)
summary(fk_mic40a) #p = 0.7076; slope = -0.002465
anova(fk_mic40a)

#linearity
plot(resid(fk_mic40a), BRARFKmic_func_mic40$ln) #looks linear
plot(fk_mic40a) #no pattern so indicates linearity

#homoscedascity
BRARFKmic_func_mic40$res <- residuals(fk_mic40a)
BRARFKmic_func_mic40$abs_res <- abs(BRARFKmic_func_mic40$res)
BRARFKmic_func_mic40$abs_res2 <- BRARFKmic_func_mic40$abs_res^2
levene_brar <- lm(data = BRARFKmic_func_mic40, abs_res2 ~ BRARrel)
anova(levene_brar) #p = 0.1777 so > 0.05 so equal variance is met

#autocorrelation
acf(residuals(fk_mic40a, retype = "normalized"))  #lines dont go outside CI horizontal lines, so not autocorrelated
pacf(residuals(fk_mic40a, retype = "normalized"))


#Alipath
fk_mic41 <- lmerTest::lmer(avg ~ BRARrel + (1|grad_num), data = subset(BRARFKmic_func_long2, func == "aliphatic_non_methane_hydrocarbon_degradation"))
summary(fk_mic41) #not sig
anova(fk_mic41)

#check normality
resid1.0 <- lm(data = subset(BRARFKmic_func_long2, func == "aliphatic_non_methane_hydrocarbon_degradation"), avg ~ BRARrel)
ols_plot_resid_hist(resid1.0) #right skew
ols_test_normality(resid1.0) #fail

#trans
BRARFKmic_func_mic41 <- BRARFKmic_func_long2 %>%
  filter(func == "aliphatic_non_methane_hydrocarbon_degradation") %>%
  mutate(ln = log10(avg + 0.001))
resid1.0 <- lm(data = BRARFKmic_func_mic41, ln ~ BRARrel)
ols_plot_resid_hist(resid1.0) #better
ols_test_normality(resid1.0) #pass
#redo stats
fk_mic41a <- lmerTest::lmer(ln ~ BRARrel + (1|grad_num), data = BRARFKmic_func_mic41)
summary(fk_mic41a) #p = 0.9283; slope = 0.0002906
anova(fk_mic41a)

#linearity
plot(resid(fk_mic41a), BRARFKmic_func_mic41$ln) #looks linear
plot(fk_mic41a) #no pattern so indicates linearity

#homoscedascity
BRARFKmic_func_mic41$res <- residuals(fk_mic41a)
BRARFKmic_func_mic41$abs_res <- abs(BRARFKmic_func_mic41$res)
BRARFKmic_func_mic41$abs_res2 <- BRARFKmic_func_mic41$abs_res^2
levene_brar <- lm(data = BRARFKmic_func_mic41, abs_res2 ~ BRARrel)
anova(levene_brar) #p = 0.9322 so > 0.05 so equal variance is met

#autocorrelation
acf(residuals(fk_mic41a, retype = "normalized"))  #lines dont go outside CI horizontal lines, so not autocorrelated
pacf(residuals(fk_mic41a, retype = "normalized"))


#Aromatic
fk_mic42 <- lmerTest::lmer(avg ~ BRARrel + (1|grad_num), data = subset(BRARFKmic_func_long2, func == "aromatic_hydrocarbon_degradation"))
summary(fk_mic42) #not sig 
anova(fk_mic42)

#Hydro deg
fk_mic43 <- lmerTest::lmer(avg ~ BRARrel + (1|grad_num), data = subset(BRARFKmic_func_long2, func == "hydrocarbon_degradation"))
summary(fk_mic43) #not sig
anova(fk_mic43)

#Human gut
fk_mic44 <- lmerTest::lmer(avg ~ BRARrel + (1|grad_num), data = subset(BRARFKmic_func_long2, func == "human_gut"))
summary(fk_mic44) #not enough data

#Mammal gut
fk_mic45 <- lmerTest::lmer(avg ~ BRARrel + (1|grad_num), data = subset(BRARFKmic_func_long2, func == "mammal_gut"))
summary(fk_mic45) #not enough data

#Plant path
fk_mic46 <- lmerTest::lmer(avg ~ BRARrel + (1|grad_num), data = subset(BRARFKmic_func_long2, func == "plant_pathogen"))
summary(fk_mic46) #not enough data

#Xylano
fk_mic47 <- lmerTest::lmer(avg ~ BRARrel + (1|grad_num), data = subset(BRARFKmic_func_long2, func == "xylanolysis"))
summary(fk_mic47) #not enough data 











#Thunder Basin
#BRAR Plants

TB2019cover_sp_BRAR3 <- TB2019cover_sp %>%
  filter(invasive_type == "BRAR") %>%
  filter(symbol != "BRAR") %>%
  full_join(TB_BRAR_2019) %>%
  group_by(plot, funct2) %>%
  summarise(avg_cov = sum(cover)) %>%
  ungroup() 


TBBRARcov_func <- TB2019cover_sp_BRAR3 %>%
  full_join(TB_commBRAR)

TB_BRAR_func_reg <- ggplot(data = TBBRARcov_func, aes(x = BRARrel, y = avg_cov, color = funct2, shape = funct2)) +
  #geom_point(size = 3) +
  geom_smooth(data = subset(TBBRARcov_func, funct2 == "C3 Annual Grass"), aes(group = funct2), method = "lm", se = TRUE, linetype = "dashed") +
  geom_smooth(data = subset(TBBRARcov_func, funct2 == "C3 Perennial Grass"), aes(group = funct2), method = "lm", se = FALSE, linetype = "dotted") +
  geom_smooth(data = subset(TBBRARcov_func, funct2 == "C4 Perennial Grass"), aes(group = funct2), method = "lm", se = TRUE) +
  geom_smooth(data = subset(TBBRARcov_func, funct2 == "Cactus"), aes(group = funct2), method = "lm", se = FALSE, linetype = "dotted") +
  geom_smooth(data = subset(TBBRARcov_func, funct2 == "Forb"), aes(group = funct2), method = "lm", se = FALSE, linetype = "dotted") +
  geom_smooth(data = subset(TBBRARcov_func, funct2 == "Sub-Shrub/Shrub"), aes(group = funct2), method = "lm", se = FALSE, linetype = "dotted") +
  scale_colour_manual(values = cbPalette_rac) +
  #scale_shape_manual(values = c(19, 18, 17, 15, 8, 7)) +
  labs(x = "BRAR Cover (%)", y = "Cover (%)", title = "WY BRAR", color = "Functional Group", shape = "Functional Group") + 
  theme(plot.title = element_text(vjust = - 7, hjust = 0.07, size = 16)) 

#C3 AG
c3ag_tb <- lmerTest::lmer(avg_cov ~ BRARrel + (1|grad_num), data = subset(TBBRARcov_func, funct2 == "C3 Annual Grass"))
summary(c3ag_tb) #p = 0.0513; slope = -0.00373
anova(c3ag_tb)

#check normality
resid1.0 <- lm(data = subset(TBBRARcov_func, funct2 == "C3 Annual Grass"), avg_cov ~ BRARrel)
ols_plot_resid_hist(resid1.0) #normalish
ols_test_normality(resid1.0) #fail
#log trans
TBBRARcov_func_c3ag_trans <- TBBRARcov_func %>%
  filter(funct2 == "C3 Annual Grass") %>%
  mutate(ln = log10(avg_cov + 0.001)) 
resid1.0 <- lm(data = TBBRARcov_func_c3ag_trans, ln ~ BRARrel)
ols_plot_resid_hist(resid1.0) #still skew, stick with untrans
ols_test_normality(resid1.0) #fail

#linearity
plot(resid(c3ag_tb), TBBRARcov_func_c3ag_trans$avg_cov) #looks more linear
plot(c3ag_tb) #not much pattern so indicates linearity

#homoscedascity
TBBRARcov_func_c3ag_trans$res <- residuals(c3ag_tb)
TBBRARcov_func_c3ag_trans$abs_res <- abs(TBBRARcov_func_c3ag_trans$res)
TBBRARcov_func_c3ag_trans$abs_res2 <- TBBRARcov_func_c3ag_trans$abs_res^2
levene_brar <- lm(data = TBBRARcov_func_c3ag_trans, abs_res2 ~ BRARrel)
anova(levene_brar) #p = 0.628 so > 0.05 so equal variance is met

#autocorrelation
acf(residuals(c3ag_tb, retype = "normalized"))  #lines dont go outside CI horizontal lines, so not autocorrelated
pacf(residuals(c3ag_tb, retype = "normalized"))


#C3 perennial
c3pg_tb <- lmerTest::lmer(avg_cov ~ BRARrel + (1|grad_num), data = subset(TBBRARcov_func, funct2 == "C3 Perennial Grass"))
summary(c3pg_tb) #not sig
anova(c3pg_tb)

#check normality
resid1.0 <- lm(data = subset(TBBRARcov_func, funct2 == "C3 Perennial Grass"), avg_cov ~ BRARrel)
ols_plot_resid_hist(resid1.0) #normalish
ols_test_normality(resid1.0) #pass
#log trans
TBBRARcov_func_c3pg_trans <- TBBRARcov_func %>%
  filter(funct2 == "C3 Perennial Grass") 

#linearity
plot(resid(c3pg_tb), TBBRARcov_func_c3pg_trans$avg_cov) #looks more linear
plot(c3pg_tb) #not much pattern so indicates linearity

#homoscedascity
TBBRARcov_func_c3pg_trans$res <- residuals(c3pg_tb)
TBBRARcov_func_c3pg_trans$abs_res <- abs(TBBRARcov_func_c3pg_trans$res)
TBBRARcov_func_c3pg_trans$abs_res2 <- TBBRARcov_func_c3pg_trans$abs_res^2
levene_brar <- lm(data = TBBRARcov_func_c3pg_trans, abs_res2 ~ BRARrel)
anova(levene_brar) #p = 0.8435 so > 0.05 so equal variance is met

#autocorrelation
acf(residuals(c3pg_tb, retype = "normalized"))  #lines dont go outside CI horizontal lines, so not autocorrelated
pacf(residuals(c3pg_tb, retype = "normalized"))


#C4 PG
c4pg_tb <- lmerTest::lmer(avg_cov ~ BRARrel + (1|grad_num), data = subset(TBBRARcov_func, funct2 == "C4 Perennial Grass"))
summary(c4pg_tb) #p = 2.86e-5; slope = -0.53154
anova(c4pg_tb)

#check normality
resid1.0 <- lm(data = subset(TBBRARcov_func, funct2 == "C4 Perennial Grass"), avg_cov ~ BRARrel)
ols_plot_resid_hist(resid1.0) #normalish
ols_test_normality(resid1.0) #pass
#log trans
TBBRARcov_func_c4pg_trans <- TBBRARcov_func %>%
  filter(funct2 == "C4 Perennial Grass") 

#linearity
plot(resid(c4pg_tb), TBBRARcov_func_c4pg_trans$avg_cov) #looks more linear
plot(c4pg_tb) #not much pattern so indicates linearity

#homoscedascity
TBBRARcov_func_c4pg_trans$res <- residuals(c4pg_tb)
TBBRARcov_func_c4pg_trans$abs_res <- abs(TBBRARcov_func_c4pg_trans$res)
TBBRARcov_func_c4pg_trans$abs_res2 <- TBBRARcov_func_c4pg_trans$abs_res^2
levene_brar <- lm(data = TBBRARcov_func_c4pg_trans, abs_res2 ~ BRARrel)
anova(levene_brar) #p = 0.6299 so > 0.05 so equal variance is met

#autocorrelation
acf(residuals(c4pg_tb, retype = "normalized"))  #lines dont go outside CI horizontal lines, so not autocorrelated
pacf(residuals(c4pg_tb, retype = "normalized"))


#Cactus
cc_tb <- lmerTest::lmer(avg_cov ~ BRARrel + (1|grad_num), data = subset(TBBRARcov_func, funct2 == "Cactus"))
summary(cc_tb) #not sig
anova(cc_tb)

#check normality
resid1.0 <- lm(data = subset(TBBRARcov_func, funct2 == "Cactus"), avg_cov ~ BRARrel)
ols_plot_resid_hist(resid1.0) #right skew
ols_test_normality(resid1.0) #fail
#log trans
TBBRARcov_func_cc_trans <- TBBRARcov_func %>%
  filter(funct2 == "Cactus") %>%
  mutate(ln = log10(avg_cov + 0.001))

resid1.0 <- lm(data = TBBRARcov_func_cc_trans, ln ~ BRARrel)
ols_plot_resid_hist(resid1.0) #right skew
ols_test_normality(resid1.0) #fail some but better

cc_tb2 <- lmerTest::lmer(ln ~ BRARrel + (1|grad_num), data = TBBRARcov_func_cc_trans)
summary(cc_tb2) #not sig
anova(cc_tb2)

#linearity
plot(resid(cc_tb2), TBBRARcov_func_cc_trans$ln) #looks more linear
plot(cc_tb2) #not much pattern so indicates linearity

#homoscedascity
TBBRARcov_func_cc_trans$res <- residuals(cc_tb2)
TBBRARcov_func_cc_trans$abs_res <- abs(TBBRARcov_func_cc_trans$res)
TBBRARcov_func_cc_trans$abs_res2 <- TBBRARcov_func_cc_trans$abs_res^2
levene_brar <- lm(data = TBBRARcov_func_cc_trans, abs_res2 ~ BRARrel)
anova(levene_brar) #p = 0.4839 so > 0.05 so equal variance is met

#autocorrelation
acf(residuals(cc_tb, retype = "normalized"))  #lines dont go outside CI horizontal lines, so not autocorrelated
pacf(residuals(cc_tb, retype = "normalized"))


#Forb
cf_tb <- lmerTest::lmer(avg_cov ~ BRARrel + (1|grad_num), data = subset(TBBRARcov_func, funct2 == "Forb"))
summary(cf_tb) #not sig
anova(cf_tb)

#check normality
resid1.0 <- lm(data = subset(TBBRARcov_func, funct2 == "Forb"), avg_cov ~ BRARrel)
ols_plot_resid_hist(resid1.0) #right skew
ols_test_normality(resid1.0) #fail
#log trans
TBBRARcov_func_cf_trans <- TBBRARcov_func %>%
  filter(funct2 == "Forb") %>%
  mutate(ln = log10(avg_cov + 0.001))

resid1.0 <- lm(data = TBBRARcov_func_cf_trans, ln ~ BRARrel)
ols_plot_resid_hist(resid1.0) #right skew
ols_test_normality(resid1.0) #fail some but better

cf_tb2 <- lmerTest::lmer(ln ~ BRARrel + (1|grad_num), data = TBBRARcov_func_cf_trans)
summary(cf_tb2) #not sig
anova(cf_tb2)

#linearity
plot(resid(cf_tb2), TBBRARcov_func_cf_trans$ln) #looks more linear
plot(cf_tb2) #not much pattern so indicates linearity

#homoscedascity
TBBRARcov_func_cf_trans$res <- residuals(cf_tb2)
TBBRARcov_func_cf_trans$abs_res <- abs(TBBRARcov_func_cf_trans$res)
TBBRARcov_func_cf_trans$abs_res2 <- TBBRARcov_func_cf_trans$abs_res^2
levene_brar <- lm(data = TBBRARcov_func_cf_trans, abs_res2 ~ BRARrel)
anova(levene_brar) #p = 0.4839 so > 0.05 so equal variance is met

#autocorrelation
acf(residuals(cf_tb, retype = "normalized"))  #lines dont go outside CI horizontal lines, so not autocorrelated
pacf(residuals(cf_tb, retype = "normalized"))

#Shrub
cs_tb <- lmerTest::lmer(avg_cov ~ BRARrel + (1|grad_num), data = subset(TBBRARcov_func, funct2 == "Sub-Shrub/Shrub"))
summary(cs_tb) #not sig
anova(cs_tb)

#check normality
resid1.0 <- lm(data = subset(TBBRARcov_func, funct2 == "Sub-Shrub/Shrub"), avg_cov ~ BRARrel)
ols_plot_resid_hist(resid1.0) #right skew
ols_test_normality(resid1.0) #fail
#log trans
TBBRARcov_func_cs_trans <- TBBRARcov_func %>%
  filter(funct2 == "Sub-Shrub/Shrub") %>%
  mutate(ln = log10(avg_cov + 0.001))

resid1.0 <- lm(data = TBBRARcov_func_cs_trans, ln ~ BRARrel)
ols_plot_resid_hist(resid1.0) #right skew
ols_test_normality(resid1.0) #fail and worse

#linearity
plot(resid(cs_tb), TBBRARcov_func_cs_trans$avg_cov) #looks more linear
plot(cs_tb2) #not much pattern so indicates linearity

#homoscedascity
TBBRARcov_func_cs_trans$res <- residuals(cs_tb)
TBBRARcov_func_cs_trans$abs_res <- abs(TBBRARcov_func_cs_trans$res)
TBBRARcov_func_cs_trans$abs_res2 <- TBBRARcov_func_cs_trans$abs_res^2
levene_brar <- lm(data = TBBRARcov_func_cs_trans, abs_res2 ~ BRARrel)
anova(levene_brar) #p = 0.6103 so > 0.05 so equal variance is met

#autocorrelation
acf(residuals(cs_tb, retype = "normalized"))  #lines dont go outside CI horizontal lines, so not autocorrelated
pacf(residuals(cs_tb, retype = "normalized"))






#BRAR Insects
TBinsID_guild_BRAR3 <- TBinsID_guild_BRAR %>%
  group_by(plot, guild) %>%
  summarise(avg = sum(avg_n)) %>%
  ungroup()

TBBRARinsfunc <- TBinsID_guild_BRAR3 %>%
  full_join(TB_commBRAR) %>%
  drop_na(guild)

TB_BRAR_func_reg_ins <- ggplot(data = TBBRARinsfunc, aes(x = BRARrel, y = avg, color = guild, shape = guild)) +
  #geom_point(size = 3) +
  geom_smooth(data = subset(TBBRARinsfunc, guild == "Leaf Chewing Herbivore"), aes(group = guild), method = "lm", se = FALSE, linetype = "dotted") +
  geom_smooth(data = subset(TBBRARinsfunc, guild == "Predator"), aes(group = guild), method = "lm", se = FALSE, linetype = "dotted") +
  geom_smooth(data = subset(TBBRARinsfunc, guild == "Sap Sucking Herbivore"), aes(group = guild), method = "lm", se = FALSE, linetype = "dotted") +
  geom_smooth(data = subset(TBBRARinsfunc, guild == "Parasitoid"), aes(group = guild), method = "lm", se = TRUE) +
  scale_colour_manual(values = cbPalette_rac) +
  #scale_shape_manual(values = c(19, 18, 17, 15, 8, 7)) +
  labs(x = "BRAR Cover (%)", y = "Abundance (count)", title = "WY BRAR", color = "Functional Group", shape = "Functional Group") + 
  theme(plot.title = element_text(vjust = - 7, hjust = 0.07, size = 16)) 

#Fungivore
fun_tb <- lmerTest::lmer(avg ~ BRARrel + (1|grad_num), data = subset(TBBRARinsfunc, guild == "Fungivore"))
summary(fun_tb) #can't include; only 1 found
anova(fun_tb)

#Leaf chewing
lch_tb <- lmerTest::lmer(avg ~ BRARrel + (1|grad_num), data = subset(TBBRARinsfunc, guild == "Leaf Chewing Herbivore"))
summary(lch_tb) #no sig
anova(lch_tb)

#check normality
resid1.0 <- lm(data = subset(TBBRARinsfunc, guild == "Leaf Chewing Herbivore"), avg ~ BRARrel)
ols_plot_resid_hist(resid1.0) #right skew
ols_test_normality(resid1.0) #fail
#log trans
TBBRARins_lchfunc <- TBBRARinsfunc %>%
  filter(guild == "Leaf Chewing Herbivore") %>%
  mutate(ln = log10(avg + 0.001))

resid1.0 <- lm(data = TBBRARins_lchfunc, ln ~ BRARrel)
ols_plot_resid_hist(resid1.0) #right skew
ols_test_normality(resid1.0) #fail some but better

lch_tb2 <- lmerTest::lmer(ln ~ BRARrel + (1|grad_num), data = TBBRARins_lchfunc)
summary(lch_tb2) #no sig
anova(lch_tb2)

#linearity
plot(resid(lch_tb2), TBBRARins_lchfunc$ln) #looks more linear
plot(lch_tb2) #not much pattern so indicates linearity

#homoscedascity
TBBRARins_lchfunc$res <- residuals(lch_tb2)
TBBRARins_lchfunc$abs_res <- abs(TBBRARins_lchfunc$res)
TBBRARins_lchfunc$abs_res2 <- TBBRARins_lchfunc$abs_res^2
levene_brar <- lm(data = TBBRARins_lchfunc, abs_res2 ~ BRARrel)
anova(levene_brar) #p = 0.8394 so > 0.05 so equal variance is met

#autocorrelation
acf(residuals(lch_tb2, retype = "normalized"))  #lines dont go outside CI horizontal lines, so not autocorrelated
pacf(residuals(lch_tb2, retype = "normalized"))


#Predator
pre_tb <- lmerTest::lmer(avg ~ BRARrel + (1|grad_num), data = subset(TBBRARinsfunc, guild == "Predator"))
summary(pre_tb) #no sig
anova(pre_tb)

#check normality
resid1.0 <- lm(data = subset(TBBRARinsfunc, guild == "Predator"), avg ~ BRARrel)
ols_plot_resid_hist(resid1.0) #right skew
ols_test_normality(resid1.0) #fail
#log trans
TBBRARins_prefunc <- TBBRARinsfunc %>%
  filter(guild == "Predator") %>%
  mutate(ln = log10(avg + 0.001))

resid1.0 <- lm(data = TBBRARins_prefunc, ln ~ BRARrel)
ols_plot_resid_hist(resid1.0) #right skew
ols_test_normality(resid1.0) #fail some but better

pre_tb2 <- lmerTest::lmer(ln ~ BRARrel + (1|grad_num), data = TBBRARins_prefunc)
summary(pre_tb2) #no sig
anova(pre_tb2)

#linearity
plot(resid(pre_tb2), TBBRARins_prefunc$ln) #looks more linear
plot(pre_tb2) #not much pattern so indicates linearity

#homoscedascity
TBBRARins_prefunc$res <- residuals(pre_tb2)
TBBRARins_prefunc$abs_res <- abs(TBBRARins_prefunc$res)
TBBRARins_prefunc$abs_res2 <- TBBRARins_prefunc$abs_res^2
levene_brar <- lm(data = TBBRARins_prefunc, abs_res2 ~ BRARrel)
anova(levene_brar) #p = 0.383 so > 0.05 so equal variance is met

#autocorrelation
acf(residuals(pre_tb2, retype = "normalized"))  #lines dont go outside CI horizontal lines, so not autocorrelated
pacf(residuals(pre_tb2, retype = "normalized"))


#Sap sucking
ssh_tb <- lmerTest::lmer(avg ~ BRARrel + (1|grad_num), data = subset(TBBRARinsfunc, guild == "Sap Sucking Herbivore"))
summary(ssh_tb) #no sig
anova(ssh_tb)

#check normality
resid1.0 <- lm(data = subset(TBBRARinsfunc, guild == "Sap Sucking Herbivore"), avg ~ BRARrel)
ols_plot_resid_hist(resid1.0) #right skew
ols_test_normality(resid1.0) #fail
#log trans
TBBRARins_sshfunc <- TBBRARinsfunc %>%
  filter(guild == "Sap Sucking Herbivore") %>%
  mutate(ln = log10(avg + 0.001))

resid1.0 <- lm(data = TBBRARins_sshfunc, ln ~ BRARrel)
ols_plot_resid_hist(resid1.0) #right skew
ols_test_normality(resid1.0) #pass

ssh_tb2 <- lmerTest::lmer(ln ~ BRARrel + (1|grad_num), data = TBBRARins_sshfunc)
summary(ssh_tb2) #p = 0.02093, slope = 0.008859
anova(ssh_tb2)

#linearity
plot(resid(ssh_tb2), TBBRARins_sshfunc$ln) #looks more linear
plot(ssh_tb2) #not much pattern so indicates linearity

#homoscedascity
TBBRARins_sshfunc$res <- residuals(ssh_tb2)
TBBRARins_sshfunc$abs_res <- abs(TBBRARins_sshfunc$res)
TBBRARins_sshfunc$abs_res2 <- TBBRARins_sshfunc$abs_res^2
levene_brar <- lm(data = TBBRARins_sshfunc, abs_res2 ~ BRARrel)
anova(levene_brar) #p = 0.2847 so > 0.05 so equal variance is met

#autocorrelation
acf(residuals(ssh_tb2, retype = "normalized"))  #lines dont go outside CI horizontal lines, so not autocorrelated
pacf(residuals(ssh_tb2, retype = "normalized"))


#Parasitoid
par_tb <- lmerTest::lmer(avg ~ BRARrel + (1|grad_num), data = subset(TBBRARinsfunc, guild == "Parasitoid"))
summary(par_tb) #p = 0.00505; slope = 0.03809
anova(par_tb)

#check normality
resid1.0 <- lm(data = subset(TBBRARinsfunc, guild == "Parasitoid"), avg ~ BRARrel)
ols_plot_resid_hist(resid1.0) #normal
ols_test_normality(resid1.0) #pass
#log trans
TBBRARins_parfunc <- TBBRARinsfunc %>%
  filter(guild == "Parasitoid") %>%
  mutate(ln = log10(avg + 0.001))

#linearity
plot(resid(par_tb), TBBRARins_parfunc$avg) #looks more linear
plot(par_tb) #not much pattern so indicates linearity

#homoscedascity
TBBRARins_parfunc$res <- residuals(par_tb)
TBBRARins_parfunc$abs_res <- abs(TBBRARins_parfunc$res)
TBBRARins_parfunc$abs_res2 <- TBBRARins_parfunc$abs_res^2
levene_brar <- lm(data = TBBRARins_parfunc, abs_res2 ~ BRARrel)
anova(levene_brar) #p = 0.5451 so > 0.05 so equal variance is met

#autocorrelation
acf(residuals(par_tb, retype = "normalized"))  #lines dont go outside CI horizontal lines, so not autocorrelated
pacf(residuals(par_tb, retype = "normalized"))



pneh_tb <- lmerTest::lmer(avg ~ BRARrel + (1|grad_num), data = subset(TBBRARinsfunc, guild == "Pollen/Nectar Eating Herbivore"))
summary(pneh_tb) #not enough data to include

oh_tb <- lmerTest::lmer(avg ~ BRARrel + (1|grad_num), data = subset(TBBRARinsfunc, guild == "Other Herbivore"))
summary(oh_tb) #not enough data to include




#BRAR Microbes
#TBmicfunc <- BRARTBmic_func_long2 %>%
dplyr::select(func) %>%
  unique()
#write.csv(TBmicfunc, "tbmic.csv")

BRARTBmic_func_long2 <- TBmic_func_long %>%
  filter(invasive_type == "BRAR") %>%
  full_join(TB_BRAR_2019) %>%
  group_by(plot, func) %>%
  summarise(avg = sum(count)) %>%
  ungroup() %>%
  drop_na(func) %>%
  full_join(TB_commBRAR)

BRARTBmic_func_long3 <- BRARTBmic_func_long2 %>%
  filter(func != "other")

#
BRARFKmic_func_long4 <- BRARFKmic_func_long2 %>%
  group_by(func) %>%
  summarise(avg2 = mean(avg)) %>%
  ungroup()
#top 7 most abundant functions (not other): chemoheterotrophy, aerobic_chemoheterotrophy, nitrate_reduction, nitrification, aerobic_ammonia_oxidation, nitrogen_fixation, manganese_oxidation 

#only whats significant
TB_BRAR_func_reg_mic <- ggplot(data = BRARTBmic_func_long3, aes(x = BRARrel, y = avg, color = func)) +
  #geom_point(size = 3) +
  geom_smooth(data = subset(BRARTBmic_func_long3, func == "aerobic_ammonia_oxidation"), aes(group = func), method = "lm", se = TRUE, linetype = "dashed") +
  geom_smooth(data = subset(BRARTBmic_func_long3, func == "chitinolysis"), aes(group = func), method = "lm", se = TRUE) +
  geom_smooth(data = subset(BRARTBmic_func_long3, func == "dark_oxidation_of_sulfur_compounds"), aes(group = func), method = "lm", se = TRUE) +
  geom_smooth(data = subset(BRARTBmic_func_long3, func == "nitrate_reduction"), aes(group = func), method = "lm", se = TRUE) +
  geom_smooth(data = subset(BRARTBmic_func_long3, func == "oxygenic_photoautotrophy"), aes(group = func), method = "lm", se = TRUE, linetype = "dashed") +
  geom_smooth(data = subset(BRARTBmic_func_long3, func == "photoautotrophy"), aes(group = func), method = "lm", se = TRUE, linetype = "dashed") +
  geom_smooth(data = subset(BRARTBmic_func_long3, func == "photosynthetic_cyanobacteria"), aes(group = func), method = "lm", se = TRUE, linetype = "dashed") +
  geom_smooth(data = subset(BRARTBmic_func_long3, func == "phototrophy"), aes(group = func), method = "lm", se = TRUE, linetype = "dashed") +
  geom_smooth(data = subset(BRARTBmic_func_long3, func == "human_associated"), aes(group = func), method = "lm", se = TRUE) +
  geom_smooth(data = subset(BRARTBmic_func_long3, func == "human_pathogens_all"), aes(group = func), method = "lm", se = TRUE, linetype = "dashed") +
  geom_smooth(data = subset(BRARTBmic_func_long3, func == "iron_respiration"), aes(group = func), method = "lm", se = TRUE, linetype = "dashed") +
  geom_smooth(data = subset(BRARTBmic_func_long3, func == "aliphatic_non_methane_hydrocarbon_degradation"), aes(group = func), method = "lm", se = TRUE) +
  geom_smooth(data = subset(BRARTBmic_func_long3, func == "aromatic_hydrocarbon_degradation"), aes(group = func), method = "lm", se = TRUE) +
  geom_smooth(data = subset(BRARTBmic_func_long3, func == "hydrocarbon_degradation"), aes(group = func), method = "lm", se = TRUE) +
  #geom_smooth(data = subset(BRARFKmic_func_long3, func != "oxygenic_photoautotrophy"), aes(group = func), method = "lm", se = FALSE, linetype = "dotted") +
  #scale_colour_manual(values = cbPalette_rac) +
  #scale_shape_manual(values = c(19, 18, 17, 15, 8, 7)) +
  ylim(0, 2500) +
  labs(x = "BRAR Cover (%)", y = "Abundance (count)", title = "WY BRAR", color = "Functional Group", shape = "Functional Group") + 
  theme(plot.title = element_text(vjust = - 7, hjust = 0.07, size = 16)) 

#
FK_BRAR_func_reg_mic2 <- ggplot(data = BRARFKmic_func_long3, aes(x = BRARrel, y = avg, color = func)) +
  #geom_point(size = 3) +
  geom_smooth(data = subset(BRARFKmic_func_long3, func == "chemoheterotrophy"), aes(group = func), method = "lm", se = FALSE, linetype = "dotted", position = position_dodge(0.5)) +
  geom_smooth(data = subset(BRARFKmic_func_long3, func == "aerobic_chemoheterotrophy"), aes(group = func), method = "lm", se = FALSE, linetype = "dotted", position = position_dodge(0.5)) +
  geom_smooth(data = subset(BRARFKmic_func_long3, func == "nitrate_reduction"), aes(group = func), method = "lm", se = FALSE, linetype = "dotted") +
  geom_smooth(data = subset(BRARFKmic_func_long3, func == "nitrification"), aes(group = func), method = "lm", se = FALSE, linetype = "dotted") +
  geom_smooth(data = subset(BRARFKmic_func_long3, func == "aerobic_ammonia_oxidation"), aes(group = func), method = "lm", se = FALSE, linetype = "dotted") +
  geom_smooth(data = subset(BRARFKmic_func_long3, func == "nitrogen_fixation"), aes(group = func), method = "lm", se = FALSE, linetype = "dotted") +
  geom_smooth(data = subset(BRARFKmic_func_long3, func == "manganese_oxidation"), aes(group = func), method = "lm", se = FALSE, linetype = "dotted") +
  scale_colour_manual(values = cbPalette_rac) +
  #scale_shape_manual(values = c(19, 18, 17, 15, 8, 7)) +
  labs(x = "BRAR Cover (%)", y = "Abundance (count)", title = "MT BRAR", color = "Functional Group", shape = "Functional Group") + 
  theme(plot.title = element_text(vjust = - 7, hjust = 0.07, size = 16)) 
#

#Aerobic ammon oxid
tb_mic1 <- lmerTest::lmer(avg ~ BRARrel + (1|grad_num), data = subset(BRARTBmic_func_long2, func == "aerobic_ammonia_oxidation"))
summary(tb_mic1) #p = 0.0979; slope = -10.819
anova(tb_mic1)

#check normality
resid1.0 <- lm(data = subset(BRARTBmic_func_long2, func == "aerobic_ammonia_oxidation"), avg ~ BRARrel)
ols_plot_resid_hist(resid1.0) 
ols_test_normality(resid1.0) 

#Aerobic chemo
tb_mic2 <- lmerTest::lmer(avg ~ BRARrel + (1|grad_num), data = subset(BRARTBmic_func_long2, func == "aerobic_chemoheterotrophy"))
summary(tb_mic2) #not sig
anova(tb_mic2)

#check normality
resid1.0 <- lm(data = subset(BRARTBmic_func_long2, func == "aerobic_chemoheterotrophy"), avg ~ BRARrel)
ols_plot_resid_hist(resid1.0) 
ols_test_normality(resid1.0) 

#Aerobic nitrite oxid
tb_mic3 <- lmerTest::lmer(avg ~ BRARrel + (1|grad_num), data = subset(BRARTBmic_func_long2, func == "aerobic_nitrite_oxidation"))
summary(tb_mic3) #not sig
anova(tb_mic3)

#check normality
resid1.0 <- lm(data = subset(BRARTBmic_func_long2, func == "aerobic_nitrite_oxidation"), avg ~ BRARrel)
ols_plot_resid_hist(resid1.0) 
ols_test_normality(resid1.0) 

#Animal para or sym
tb_mic4 <- lmerTest::lmer(avg ~ BRARrel + (1|grad_num), data = subset(BRARTBmic_func_long2, func == "animal_parasites_or_symbionts"))
summary(tb_mic4) #not sig
anova(tb_mic4)

#check normality
resid1.0 <- lm(data = subset(BRARTBmic_func_long2, func == "animal_parasites_or_symbionts"), avg ~ BRARrel)
ols_plot_resid_hist(resid1.0) 
ols_test_normality(resid1.0) 
#log trans
BRARTBmic_func_4 <- BRARTBmic_func_long2 %>%
  filter(func == "animal_parasites_or_symbionts") %>%
  mutate(ln = log10(avg + 0.001)) 

resid1.0 <- lm(data = BRARTBmic_func_4, ln ~ BRARrel)
ols_plot_resid_hist(resid1.0) 
ols_test_normality(resid1.0) 

tb_mic4a <- lmerTest::lmer(ln ~ BRARrel + (1|grad_num), data = BRARTBmic_func_4)
summary(tb_mic4a) #not sig
anova(tb_mic4a)

#ANOXY PHOTO
tb_mic5 <- lmerTest::lmer(avg ~ BRARrel + (1|grad_num), data = subset(BRARTBmic_func_long2, func == "anoxygenic_photoautotrophy"))
summary(tb_mic5) #not sig 
anova(tb_mic5)

#check normality
resid1.0 <- lm(data = subset(BRARTBmic_func_long2, func == "anoxygenic_photoautotrophy"), avg ~ BRARrel)
ols_plot_resid_hist(resid1.0) 
ols_test_normality(resid1.0) 

#Anoxy S
tb_mic6 <- lmerTest::lmer(avg ~ BRARrel + (1|grad_num), data = subset(BRARTBmic_func_long2, func == "anoxygenic_photoautotrophy_S_oxidizing"))
summary(tb_mic6) #not sig 
anova(tb_mic6)

#Aromatic comp deg
tb_mic7 <- lmerTest::lmer(avg ~ BRARrel + (1|grad_num), data = subset(BRARTBmic_func_long2, func == "aromatic_compound_degradation"))
summary(tb_mic7) #not sig
anova(tb_mic7)

#check normality
resid1.0 <- lm(data = subset(BRARTBmic_func_long2, func == "aromatic_compound_degradation"), avg ~ BRARrel)
ols_plot_resid_hist(resid1.0) 
ols_test_normality(resid1.0) 

#Chemohetero
tb_mic8 <- lmerTest::lmer(avg ~ BRARrel + (1|grad_num), data = subset(BRARTBmic_func_long2, func == "chemoheterotrophy"))
summary(tb_mic8) #not sig 
anova(tb_mic8)

#check normality
resid1.0 <- lm(data = subset(BRARTBmic_func_long2, func == "chemoheterotrophy"), avg ~ BRARrel)
ols_plot_resid_hist(resid1.0) 
ols_test_normality(resid1.0) 
#log trans
BRARTBmic_func_8 <- BRARTBmic_func_long2 %>%
  filter(func == "chemoheterotrophy") %>%
  mutate(ln = log10(avg + 0.001)) 

resid1.0 <- lm(data = BRARTBmic_func_8, ln ~ BRARrel)
ols_plot_resid_hist(resid1.0) 
ols_test_normality(resid1.0) 

tb_mic8a <- lmerTest::lmer(ln ~ BRARrel + (1|grad_num), data = BRARTBmic_func_8)
summary(tb_mic8a) #not sig
anova(tb_mic8a)

#Chitino
tb_mic9 <- lmerTest::lmer(avg ~ BRARrel + (1|grad_num), data = subset(BRARTBmic_func_long2, func == "chitinolysis"))
summary(tb_mic9) #p = 0.0473; slope = 0.6791
anova(tb_mic9)

#check normality
resid1.0 <- lm(data = subset(BRARTBmic_func_long2, func == "chitinolysis"), avg ~ BRARrel)
ols_plot_resid_hist(resid1.0) 
ols_test_normality(resid1.0) 
#log trans
BRARTBmic_func_9 <- BRARTBmic_func_long2 %>%
  filter(func == "chitinolysis") %>%
  mutate(ln = log10(avg + 0.001)) 

resid1.0 <- lm(data = BRARTBmic_func_9, ln ~ BRARrel)
ols_plot_resid_hist(resid1.0) 
ols_test_normality(resid1.0) 

tb_mic9a <- lmerTest::lmer(ln ~ BRARrel + (1|grad_num), data = BRARTBmic_func_9)
summary(tb_mic9a) #not sig
anova(tb_mic9a)

#Dark oxid S
tb_mic10 <- lmerTest::lmer(avg ~ BRARrel + (1|grad_num), data = subset(BRARTBmic_func_long2, func == "dark_oxidation_of_sulfur_compounds"))
summary(tb_mic10) #p = 0.03625; slope = 0.25304
anova(tb_mic10)

#check normality
resid1.0 <- lm(data = subset(BRARTBmic_func_long2, func == "dark_oxidation_of_sulfur_compounds"), avg ~ BRARrel)
ols_plot_resid_hist(resid1.0) 
ols_test_normality(resid1.0) 
#log trans
BRARTBmic_func_10 <- BRARTBmic_func_long2 %>%
  filter(func == "dark_oxidation_of_sulfur_compounds") %>%
  mutate(ln = log10(avg + 0.001)) 

resid1.0 <- lm(data = BRARTBmic_func_10, ln ~ BRARrel)
ols_plot_resid_hist(resid1.0) 
ols_test_normality(resid1.0) 

tb_mic10a <- lmerTest::lmer(ln ~ BRARrel + (1|grad_num), data = BRARTBmic_func_10)
summary(tb_mic10a) #not sig
anova(tb_mic10a)

#Denitrif
tb_mic11 <- lmerTest::lmer(avg ~ BRARrel + (1|grad_num), data = subset(BRARTBmic_func_long2, func == "denitrification"))
summary(tb_mic11) #not sig
anova(tb_mic11)

#check normality
resid1.0 <- lm(data = subset(BRARTBmic_func_long2, func == "denitrification"), avg ~ BRARrel)
ols_plot_resid_hist(resid1.0) 
ols_test_normality(resid1.0) 

#Fermentation
tb_mic12 <- lmerTest::lmer(avg ~ BRARrel + (1|grad_num), data = subset(BRARTBmic_func_long2, func == "fermentation"))
summary(tb_mic12) #not sig
anova(tb_mic12)

#check normality
resid1.0 <- lm(data = subset(BRARTBmic_func_long2, func == "fermentation"), avg ~ BRARrel)
ols_plot_resid_hist(resid1.0) 
ols_test_normality(resid1.0) 
#log trans
BRARTBmic_func_12 <- BRARTBmic_func_long2 %>%
  filter(func == "fermentation") %>%
  mutate(ln = log10(avg + 0.001)) 

resid1.0 <- lm(data = BRARTBmic_func_12, ln ~ BRARrel)
ols_plot_resid_hist(resid1.0) 
ols_test_normality(resid1.0) 

tb_mic12a <- lmerTest::lmer(ln ~ BRARrel + (1|grad_num), data = BRARTBmic_func_12)
summary(tb_mic12a) #not sig
anova(tb_mic12a)

#Mang oxid
tb_mic13 <- lmerTest::lmer(avg ~ BRARrel + (1|grad_num), data = subset(BRARTBmic_func_long2, func == "manganese_oxidation"))
summary(tb_mic13) #not sig
anova(tb_mic13)

#check normality
resid1.0 <- lm(data = subset(BRARTBmic_func_long2, func == "manganese_oxidation"), avg ~ BRARrel)
ols_plot_resid_hist(resid1.0) 
ols_test_normality(resid1.0) 
#log trans
BRARTBmic_func_13 <- BRARTBmic_func_long2 %>%
  filter(func == "manganese_oxidation") %>%
  mutate(ln = log10(avg + 0.001)) 

resid1.0 <- lm(data = BRARTBmic_func_13, ln ~ BRARrel)
ols_plot_resid_hist(resid1.0) 
ols_test_normality(resid1.0) 

tb_mic13a <- lmerTest::lmer(ln ~ BRARrel + (1|grad_num), data = BRARTBmic_func_13)
summary(tb_mic13a) #not sig
anova(tb_mic13a)

#Meth oxid
tb_mic14 <- lmerTest::lmer(avg ~ BRARrel + (1|grad_num), data = subset(BRARTBmic_func_long2, func == "methanol_oxidation"))
summary(tb_mic14) #not sig 
anova(tb_mic14)

#check normality
resid1.0 <- lm(data = subset(BRARTBmic_func_long2, func == "methanol_oxidation"), avg ~ BRARrel)
ols_plot_resid_hist(resid1.0) 
ols_test_normality(resid1.0) 
#log trans
BRARTBmic_func_14 <- BRARTBmic_func_long2 %>%
  filter(func == "methanol_oxidation") %>%
  mutate(ln = log10(avg + 0.001)) 

resid1.0 <- lm(data = BRARTBmic_func_14, ln ~ BRARrel)
ols_plot_resid_hist(resid1.0) 
ols_test_normality(resid1.0) 

tb_mic14a <- lmerTest::lmer(ln ~ BRARrel + (1|grad_num), data = BRARTBmic_func_14)
summary(tb_mic14a) #not sig
anova(tb_mic14a)

#Methlyo
tb_mic15 <- lmerTest::lmer(avg ~ BRARrel + (1|grad_num), data = subset(BRARTBmic_func_long2, func == "methylotrophy"))
summary(tb_mic15) #not sig
anova(tb_mic15)

#Nitrate denit
tb_mic16 <- lmerTest::lmer(avg ~ BRARrel + (1|grad_num), data = subset(BRARTBmic_func_long2, func == "nitrate_denitrification"))
summary(tb_mic16) #not sig
anova(tb_mic16)

#Nitrate reduc
tb_mic17 <- lmerTest::lmer(avg ~ BRARrel + (1|grad_num), data = subset(BRARTBmic_func_long2, func == "nitrate_reduction"))
summary(tb_mic17) #p = 0.027309; slope = -38.293
anova(tb_mic17)

#check normality
resid1.0 <- lm(data = subset(BRARTBmic_func_long2, func == "nitrate_reduction"), avg ~ BRARrel)
ols_plot_resid_hist(resid1.0) 
ols_test_normality(resid1.0) 
#log trans
BRARTBmic_func_17 <- BRARTBmic_func_long2 %>%
  filter(func == "nitrate_reduction") %>%
  mutate(ln = log10(avg + 0.001)) 

resid1.0 <- lm(data = BRARTBmic_func_17, ln ~ BRARrel)
ols_plot_resid_hist(resid1.0) 
ols_test_normality(resid1.0) 

tb_mic17a <- lmerTest::lmer(ln ~ BRARrel + (1|grad_num), data = BRARTBmic_func_17)
summary(tb_mic17a) #not sig
anova(tb_mic17a)

#Nitrate resp
tb_mic18 <- lmerTest::lmer(avg ~ BRARrel + (1|grad_num), data = subset(BRARTBmic_func_long2, func == "nitrate_respiration"))
summary(tb_mic18) #not sig 
anova(tb_mic18)

#check normality
resid1.0 <- lm(data = subset(BRARTBmic_func_long2, func == "nitrate_respiration"), avg ~ BRARrel)
ols_plot_resid_hist(resid1.0) 
ols_test_normality(resid1.0) 

#Nitrification
tb_mic19 <- lmerTest::lmer(avg ~ BRARrel + (1|grad_num), data = subset(BRARTBmic_func_long2, func == "nitrification"))
summary(tb_mic19) #not sig
anova(tb_mic19)

#check normality
resid1.0 <- lm(data = subset(BRARTBmic_func_long2, func == "nitrification"), avg ~ BRARrel)
ols_plot_resid_hist(resid1.0) 
ols_test_normality(resid1.0) 

#Nitrite denit
tb_mic20 <- lmerTest::lmer(avg ~ BRARrel + (1|grad_num), data = subset(BRARTBmic_func_long2, func == "nitrite_denitrification"))
summary(tb_mic20) #not sig
anova(tb_mic20)

#Nitrite resp
tb_mic21 <- lmerTest::lmer(avg ~ BRARrel + (1|grad_num), data = subset(BRARTBmic_func_long2, func == "nitrite_respiration"))
summary(tb_mic21) #not sig
anova(tb_mic21)

#Nitrogen fix
tb_mic22 <- lmerTest::lmer(avg ~ BRARrel + (1|grad_num), data = subset(BRARTBmic_func_long2, func == "nitrogen_fixation"))
summary(tb_mic22) #not sig
anova(tb_mic22)

#check normality
resid1.0 <- lm(data = subset(BRARTBmic_func_long2, func == "nitrogen_fixation"), avg ~ BRARrel)
ols_plot_resid_hist(resid1.0) 
ols_test_normality(resid1.0) 

#Nitrogen resp
tb_mic23 <- lmerTest::lmer(avg ~ BRARrel + (1|grad_num), data = subset(BRARTBmic_func_long2, func == "nitrogen_respiration"))
summary(tb_mic23) #not sig
anova(tb_mic23)

#Nitrous oxide den
tb_mic24 <- lmerTest::lmer(avg ~ BRARrel + (1|grad_num), data = subset(BRARTBmic_func_long2, func == "nitrous_oxide_denitrification"))
summary(tb_mic24) #not sig
anova(tb_mic24)

#Nonphoto cyano
tb_mic25 <- lmerTest::lmer(avg ~ BRARrel + (1|grad_num), data = subset(BRARTBmic_func_long2, func == "nonphotosynthetic_cyanobacteria"))
summary(tb_mic25) #not sig
anova(tb_mic25)

#check normality
resid1.0 <- lm(data = subset(BRARTBmic_func_long2, func == "nonphotosynthetic_cyanobacteria"), avg ~ BRARrel)
ols_plot_resid_hist(resid1.0) 
ols_test_normality(resid1.0) 
#log trans
BRARTBmic_func_25 <- BRARTBmic_func_long2 %>%
  filter(func == "nonphotosynthetic_cyanobacteria") %>%
  mutate(ln = log10(avg + 0.001)) 

resid1.0 <- lm(data = BRARTBmic_func_25, ln ~ BRARrel)
ols_plot_resid_hist(resid1.0) 
ols_test_normality(resid1.0) 

tb_mic25a <- lmerTest::lmer(ln ~ BRARrel + (1|grad_num), data = BRARTBmic_func_25)
summary(tb_mic25a) #not sig
anova(tb_mic25a)

#Other
tb_mic26 <- lmerTest::lmer(avg ~ BRARrel + (1|grad_num), data = subset(BRARTBmic_func_long2, func == "other"))
summary(tb_mic26) #not sig
anova(tb_mic26)

#check normality
resid1.0 <- lm(data = subset(BRARTBmic_func_long2, func == "other"), avg ~ BRARrel)
ols_plot_resid_hist(resid1.0) 
ols_test_normality(resid1.0) 

#Oxy photo
tb_mic27 <- lmerTest::lmer(avg ~ BRARrel + (1|grad_num), data = subset(BRARTBmic_func_long2, func == "oxygenic_photoautotrophy"))
summary(tb_mic27) #p = 0.064801; slope = -1.1677
anova(tb_mic27)

#check normality
resid1.0 <- lm(data = subset(BRARTBmic_func_long2, func == "oxygenic_photoautotrophy"), avg ~ BRARrel)
ols_plot_resid_hist(resid1.0) 
ols_test_normality(resid1.0) 
#log trans
BRARTBmic_func_27 <- BRARTBmic_func_long2 %>%
  filter(func == "oxygenic_photoautotrophy") %>%
  mutate(ln = log10(avg + 0.001)) 

resid1.0 <- lm(data = BRARTBmic_func_27, ln ~ BRARrel)
ols_plot_resid_hist(resid1.0) 
ols_test_normality(resid1.0) 

tb_mic27a <- lmerTest::lmer(ln ~ BRARrel + (1|grad_num), data = BRARTBmic_func_27)
summary(tb_mic27a)
anova(tb_mic27a)

#Photoauto
tb_mic28 <- lmerTest::lmer(avg ~ BRARrel + (1|grad_num), data = subset(BRARTBmic_func_long2, func == "photoautotrophy"))
summary(tb_mic28) #p = 0.08493; slope = -46.052
anova(tb_mic28)

#check normality
resid1.0 <- lm(data = subset(BRARTBmic_func_long2, func == "photoautotrophy"), avg ~ BRARrel)
ols_plot_resid_hist(resid1.0) 
ols_test_normality(resid1.0) 
#log trans
BRARTBmic_func_28 <- BRARTBmic_func_long2 %>%
  filter(func == "photoautotrophy") %>%
  mutate(ln = log10(avg + 0.001)) 

resid1.0 <- lm(data = BRARTBmic_func_28, ln ~ BRARrel)
ols_plot_resid_hist(resid1.0) 
ols_test_normality(resid1.0) 

tb_mic28a <- lmerTest::lmer(ln ~ BRARrel + (1|grad_num), data = BRARTBmic_func_28)
summary(tb_mic28a)
anova(tb_mic28a)

#Photohetero
tb_mic29 <- lmerTest::lmer(avg ~ BRARrel + (1|grad_num), data = subset(BRARTBmic_func_long2, func == "photoheterotrophy"))
summary(tb_mic29) #not sig
anova(tb_mic29)

#check normality
resid1.0 <- lm(data = subset(BRARTBmic_func_long2, func == "photoheterotrophy"), avg ~ BRARrel)
ols_plot_resid_hist(resid1.0) 
ols_test_normality(resid1.0) 

#Photo cyano
tb_mic30 <- lmerTest::lmer(avg ~ BRARrel + (1|grad_num), data = subset(BRARTBmic_func_long2, func == "photosynthetic_cyanobacteria"))
summary(tb_mic30) #p = 0.064801; slope = -1.1677
anova(tb_mic30)

#Phototrophy
tb_mic31 <- lmerTest::lmer(avg ~ BRARrel + (1|grad_num), data = subset(BRARTBmic_func_long2, func == "phototrophy"))
summary(tb_mic31) #p = 0.08599; slope = -46.970
anova(tb_mic31)

#check normality
resid1.0 <- lm(data = subset(BRARTBmic_func_long2, func == "phototrophy"), avg ~ BRARrel)
ols_plot_resid_hist(resid1.0) 
ols_test_normality(resid1.0) 
#log trans
BRARTBmic_func_31 <- BRARTBmic_func_long2 %>%
  filter(func == "phototrophy") %>%
  mutate(ln = log10(avg + 0.001)) 

resid1.0 <- lm(data = BRARTBmic_func_31, ln ~ BRARrel)
ols_plot_resid_hist(resid1.0) 
ols_test_normality(resid1.0) 

tb_mic31a <- lmerTest::lmer(ln ~ BRARrel + (1|grad_num), data = BRARTBmic_func_31)
summary(tb_mic31a)
anova(tb_mic31a)

#Pred or exo
tb_mic32 <- lmerTest::lmer(avg ~ BRARrel + (1|grad_num), data = subset(BRARTBmic_func_long2, func == "predatory_or_exoparasitic"))
summary(tb_mic32) #not sig 
anova(tb_mic32)

#check normality
resid1.0 <- lm(data = subset(BRARTBmic_func_long2, func == "predatory_or_exoparasitic"), avg ~ BRARrel)
ols_plot_resid_hist(resid1.0) 
ols_test_normality(resid1.0) 
#log trans
BRARTBmic_func_32 <- BRARTBmic_func_long2 %>%
  filter(func == "predatory_or_exoparasitic") %>%
  mutate(ln = log10(avg + 0.001)) 

resid1.0 <- lm(data = BRARTBmic_func_32, ln ~ BRARrel)
ols_plot_resid_hist(resid1.0) 
ols_test_normality(resid1.0) 

tb_mic32a <- lmerTest::lmer(ln ~ BRARrel + (1|grad_num), data = BRARTBmic_func_32)
summary(tb_mic32a)
anova(tb_mic32a)

#Resp of sulfur
tb_mic33 <- lmerTest::lmer(avg ~ BRARrel + (1|grad_num), data = subset(BRARTBmic_func_long2, func == "respiration_of_sulfur_compounds"))
summary(tb_mic33) #not sig
anova(tb_mic33)

#check normality
resid1.0 <- lm(data = subset(BRARTBmic_func_long2, func == "respiration_of_sulfur_compounds"), avg ~ BRARrel)
ols_plot_resid_hist(resid1.0) 
ols_test_normality(resid1.0) 
#log trans
BRARTBmic_func_33 <- BRARTBmic_func_long2 %>%
  filter(func == "respiration_of_sulfur_compounds") %>%
  mutate(ln = log10(avg + 0.001)) 

resid1.0 <- lm(data = BRARTBmic_func_33, ln ~ BRARrel)
ols_plot_resid_hist(resid1.0) 
ols_test_normality(resid1.0) 

tb_mic33a <- lmerTest::lmer(ln ~ BRARrel + (1|grad_num), data = BRARTBmic_func_33)
summary(tb_mic33a)
anova(tb_mic33a)

#Sulf resp
tb_mic34 <- lmerTest::lmer(avg ~ BRARrel + (1|grad_num), data = subset(BRARTBmic_func_long2, func == "sulfate_respiration"))
summary(tb_mic34) #not sig
anova(tb_mic34)

#Ureolysis
tb_mic35 <- lmerTest::lmer(avg ~ BRARrel + (1|grad_num), data = subset(BRARTBmic_func_long2, func == "ureolysis"))
summary(tb_mic35) #not sig
anova(tb_mic35)

#check normality
resid1.0 <- lm(data = subset(BRARTBmic_func_long2, func == "ureolysis"), avg ~ BRARrel)
ols_plot_resid_hist(resid1.0) 
ols_test_normality(resid1.0) 

#Cellulo
tb_mic36 <- lmerTest::lmer(avg ~ BRARrel + (1|grad_num), data = subset(BRARTBmic_func_long2, func == "cellulolysis"))
summary(tb_mic36) #not sig
anova(tb_mic36)

#check normality
resid1.0 <- lm(data = subset(BRARTBmic_func_long2, func == "cellulolysis"), avg ~ BRARrel)
ols_plot_resid_hist(resid1.0) 
ols_test_normality(resid1.0) 
#log trans
BRARTBmic_func_36 <- BRARTBmic_func_long2 %>%
  filter(func == "cellulolysis") %>%
  mutate(ln = log10(avg + 0.001)) 

resid1.0 <- lm(data = BRARTBmic_func_36, ln ~ BRARrel)
ols_plot_resid_hist(resid1.0) 
ols_test_normality(resid1.0) 

tb_mic36a <- lmerTest::lmer(ln ~ BRARrel + (1|grad_num), data = BRARTBmic_func_36)
summary(tb_mic36a)
anova(tb_mic36a)

#Human assoc
tb_mic37 <- lmerTest::lmer(avg ~ BRARrel + (1|grad_num), data = subset(BRARTBmic_func_long2, func == "human_associated"))
summary(tb_mic37) #p = 0.0155; slope = -3.362
anova(tb_mic37)

#check normality
resid1.0 <- lm(data = subset(BRARTBmic_func_long2, func == "human_associated"), avg ~ BRARrel)
ols_plot_resid_hist(resid1.0) 
ols_test_normality(resid1.0) 
#log trans
BRARTBmic_func_37 <- BRARTBmic_func_long2 %>%
  filter(func == "human_associated") %>%
  mutate(ln = log10(avg + 0.001)) 

resid1.0 <- lm(data = BRARTBmic_func_37, ln ~ BRARrel)
ols_plot_resid_hist(resid1.0) 
ols_test_normality(resid1.0) 

tb_mic37a <- lmerTest::lmer(ln ~ BRARrel + (1|grad_num), data = BRARTBmic_func_37)
summary(tb_mic37a)
anova(tb_mic37a)

#Human path all
tb_mic38 <- lmerTest::lmer(avg ~ BRARrel + (1|grad_num), data = subset(BRARTBmic_func_long2, func == "human_pathogens_all"))
summary(tb_mic38) #p = 0.0764; slope = -2.773
anova(tb_mic38)

#check normality
resid1.0 <- lm(data = subset(BRARTBmic_func_long2, func == "human_pathogens_all"), avg ~ BRARrel)
ols_plot_resid_hist(resid1.0) 
ols_test_normality(resid1.0) 
#log trans
BRARTBmic_func_38 <- BRARTBmic_func_long2 %>%
  filter(func == "human_pathogens_all") %>%
  mutate(ln = log10(avg + 0.001)) 

resid1.0 <- lm(data = BRARTBmic_func_38, ln ~ BRARrel)
ols_plot_resid_hist(resid1.0) 
ols_test_normality(resid1.0) 

tb_mic38a <- lmerTest::lmer(ln ~ BRARrel + (1|grad_num), data = BRARTBmic_func_38)
summary(tb_mic38a)
anova(tb_mic38a)

#Iron resp
tb_mic39 <- lmerTest::lmer(avg ~ BRARrel + (1|grad_num), data = subset(BRARTBmic_func_long2, func == "iron_respiration"))
summary(tb_mic39) #p = 0.0658; slope = -1.078
anova(tb_mic39)

#check normality
resid1.0 <- lm(data = subset(BRARTBmic_func_long2, func == "iron_respiration"), avg ~ BRARrel)
ols_plot_resid_hist(resid1.0) 
ols_test_normality(resid1.0) 

#Intracell para
tb_mic40 <- lmerTest::lmer(avg ~ BRARrel + (1|grad_num), data = subset(BRARTBmic_func_long2, func == "intracellular_parasites"))
summary(tb_mic40) #not sig
anova(tb_mic40)

#check normality
resid1.0 <- lm(data = subset(BRARTBmic_func_long2, func == "intracellular_parasites"), avg ~ BRARrel)
ols_plot_resid_hist(resid1.0) 
ols_test_normality(resid1.0) 

#Alipath nonmeth
tb_mic41 <- lmerTest::lmer(avg ~ BRARrel + (1|grad_num), data = subset(BRARTBmic_func_long2, func == "aliphatic_non_methane_hydrocarbon_degradation"))
summary(tb_mic41) #p = 0.0159; slope = -2.27736
anova(tb_mic41)

#check normality
resid1.0 <- lm(data = subset(BRARTBmic_func_long2, func == "aliphatic_non_methane_hydrocarbon_degradation"), avg ~ BRARrel)
ols_plot_resid_hist(resid1.0) 
ols_test_normality(resid1.0) 
#log trans
BRARTBmic_func_41 <- BRARTBmic_func_long2 %>%
  filter(func == "aliphatic_non_methane_hydrocarbon_degradation") %>%
  mutate(ln = log10(avg + 0.001)) 

resid1.0 <- lm(data = BRARTBmic_func_41, ln ~ BRARrel)
ols_plot_resid_hist(resid1.0) 
ols_test_normality(resid1.0) 

tb_mic41a <- lmerTest::lmer(ln ~ BRARrel + (1|grad_num), data = BRARTBmic_func_41)
summary(tb_mic41a)
anova(tb_mic41a)

#Aromatic hydrocar
tb_mic42 <- lmerTest::lmer(avg ~ BRARrel + (1|grad_num), data = subset(BRARTBmic_func_long2, func == "aromatic_hydrocarbon_degradation"))
summary(tb_mic42) #p = 0.0159; slope = -2.27736
anova(tb_mic42)

#Hydrocarbon deg
tb_mic43 <- lmerTest::lmer(avg ~ BRARrel + (1|grad_num), data = subset(BRARTBmic_func_long2, func == "hydrocarbon_degradation"))
summary(tb_mic43) #p = 0.0159; slope = -2.27736
anova(tb_mic43)

#Human gut
tb_mic44 <- lmerTest::lmer(avg ~ BRARrel + (1|grad_num), data = subset(BRARTBmic_func_long2, func == "human_gut"))
summary(tb_mic44) #not sig
anova(tb_mic44)

#check normality
resid1.0 <- lm(data = subset(BRARTBmic_func_long2, func == "human_gut"), avg ~ BRARrel)
ols_plot_resid_hist(resid1.0) 
ols_test_normality(resid1.0) 
#log trans
BRARTBmic_func_44 <- BRARTBmic_func_long2 %>%
  filter(func == "human_gut") %>%
  mutate(ln = log10(avg + 0.001)) 

resid1.0 <- lm(data = BRARTBmic_func_44, ln ~ BRARrel)
ols_plot_resid_hist(resid1.0) 
ols_test_normality(resid1.0) 

tb_mic44a <- lmerTest::lmer(ln ~ BRARrel + (1|grad_num), data = BRARTBmic_func_44)
summary(tb_mic44a)
anova(tb_mic44a)

#Mammal gut
tb_mic45 <- lmerTest::lmer(avg ~ BRARrel + (1|grad_num), data = subset(BRARTBmic_func_long2, func == "mammal_gut"))
summary(tb_mic45) #not sig
anova(tb_mic45)

#Plant path
tb_mic46 <- lmerTest::lmer(avg ~ BRARrel + (1|grad_num), data = subset(BRARTBmic_func_long2, func == "plant_pathogen"))
summary(tb_mic46) #not enough data
anova(tb_mic46)

#Xyano
tb_mic47 <- lmerTest::lmer(avg ~ BRARrel + (1|grad_num), data = subset(BRARTBmic_func_long2, func == "xylanolysis"))
summary(tb_mic47) #not enough data
anova(tb_mic47)











#BRTE Plants
TB2019cover_sp_BRTE3 <- TB2019cover_sp %>%
  filter(invasive_type == "BRTE") %>%
  filter(symbol != "BRTE") %>%
  full_join(TB_BRTE_2019) %>%
  group_by(plot, funct2) %>%
  summarise(avg_cov = sum(cover)) %>%
  ungroup() 


TBBRTEcov_func <- TB2019cover_sp_BRTE3 %>%
  full_join(TB_commBRTE)

TB_BRTE_func_reg <- ggplot(data = TBBRTEcov_func, aes(x = BRTErel, y = avg_cov, color = funct2, shape = funct2)) +
  #geom_point(size = 3) +
  geom_smooth(data = subset(TBBRTEcov_func, funct2 == "C3 Annual Grass"), aes(group = funct2), method = "lm", se = TRUE) +
  geom_smooth(data = subset(TBBRTEcov_func, funct2 == "C3 Perennial Grass"), aes(group = funct2), method = "lm", se = TRUE) +
  geom_smooth(data = subset(TBBRTEcov_func, funct2 == "C4 Perennial Grass"), aes(group = funct2), method = "lm", se = TRUE) +
  geom_smooth(data = subset(TBBRTEcov_func, funct2 == "Cactus"), aes(group = funct2), method = "lm", se = FALSE, linetype = "dotted") +
  geom_smooth(data = subset(TBBRTEcov_func, funct2 == "Forb"), aes(group = funct2), method = "lm", se = FALSE, linetype = "dotted") +
  geom_smooth(data = subset(TBBRTEcov_func, funct2 == "Sub-Shrub/Shrub"), aes(group = funct2), method = "lm", se = FALSE, linetype = "dotted") +
  scale_colour_manual(values = cbPalette_rac) +
  #scale_shape_manual(values = c(19, 18, 17, 15, 8, 7)) +
  labs(x = "BRTE Cover (%)", y = "Cover (%)", title = "WY BRTE", color = "Functional Group", shape = "Functional Group") + 
  theme(plot.title = element_text(vjust = - 7, hjust = 0.07, size = 16)) 


#C3 AG
c3ag_tb2 <- lmerTest::lmer(avg_cov ~ BRTErel + (1|grad_num), data = subset(TBBRTEcov_func, funct2 == "C3 Annual Grass"))
summary(c3ag_tb2) #p = 0.04873; slope = -0.004728
anova(c3ag_tb2, type = 3)

#check normality
resid1.0 <- lm(data = subset(TBBRTEcov_func, funct2 == "C3 Annual Grass"), avg_cov ~ BRTErel)
ols_plot_resid_hist(resid1.0) #normalish
ols_test_normality(resid1.0) #pass some
#log trans
TBBRTEcov_func_c3ag_trans <- TBBRTEcov_func %>%
  filter(funct2 == "C3 Annual Grass") %>%
  mutate(ln = log10(avg_cov + 0.001)) 
resid1.0 <- lm(data = TBBRTEcov_func_c3ag_trans, ln ~ BRTErel)
ols_plot_resid_hist(resid1.0) #still skew, stick with untrans
ols_test_normality(resid1.0) #fail

#linearity
plot(resid(c3ag_tb2), TBBRTEcov_func_c3ag_trans$avg_cov) #looks more linear
plot(c3ag_tb2) #not much pattern so indicates linearity

#homoscedascity
TBBRTEcov_func_c3ag_trans$res <- residuals(c3ag_tb2)
TBBRTEcov_func_c3ag_trans$abs_res <- abs(TBBRTEcov_func_c3ag_trans$res)
TBBRTEcov_func_c3ag_trans$abs_res2 <- TBBRTEcov_func_c3ag_trans$abs_res^2
levene_brte <- lm(data = TBBRTEcov_func_c3ag_trans, abs_res2 ~ BRTErel)
anova(levene_brte) #p = 0.7086 so > 0.05 so equal variance is met

#autocorrelation
acf(residuals(c3ag_tb2, retype = "normalized"))  #lines dont go outside CI horizontal lines, so not autocorrelated
pacf(residuals(c3ag_tb2, retype = "normalized"))


#C3 PG
c3pg_tb2 <- lmerTest::lmer(avg_cov ~ BRTErel + (1|grad_num), data = subset(TBBRTEcov_func, funct2 == "C3 Perennial Grass"))
summary(c3pg_tb2) #p = 0.0468; slope = -0.05450
anova(c3pg_tb2)

#check normality
resid1.0 <- lm(data = subset(TBBRTEcov_func, funct2 == "C3 Perennial Grass"), avg_cov ~ BRTErel)
ols_plot_resid_hist(resid1.0) #normalish
ols_test_normality(resid1.0) #pass
#log trans
TBBRTEcov_func_c3pg_trans <- TBBRTEcov_func %>%
  filter(funct2 == "C3 Perennial Grass") %>%
  mutate(ln = log10(avg_cov + 0.001)) 

#linearity
plot(resid(c3pg_tb2), TBBRTEcov_func_c3pg_trans$avg_cov) #looks more linear
plot(c3pg_tb2) #not much pattern so indicates linearity

#homoscedascity
TBBRTEcov_func_c3pg_trans$res <- residuals(c3pg_tb2)
TBBRTEcov_func_c3pg_trans$abs_res <- abs(TBBRTEcov_func_c3pg_trans$res)
TBBRTEcov_func_c3pg_trans$abs_res2 <- TBBRTEcov_func_c3pg_trans$abs_res^2
levene_brte <- lm(data = TBBRTEcov_func_c3pg_trans, abs_res2 ~ BRTErel)
anova(levene_brte) #p = 0.3288 so > 0.05 so equal variance is met

#autocorrelation
acf(residuals(c3pg_tb2, retype = "normalized"))  #lines dont go outside CI horizontal lines, so not autocorrelated
pacf(residuals(c3pg_tb2, retype = "normalized"))


#C4 PG
c4pg_tb2 <- lmerTest::lmer(avg_cov ~ BRTErel + (1|grad_num), data = subset(TBBRTEcov_func, funct2 == "C4 Perennial Grass"))
summary(c4pg_tb2) #p = 5.07e-8; slope = -0.4311
anova(c4pg_tb2)

#check normality
resid1.0 <- lm(data = subset(TBBRTEcov_func, funct2 == "C4 Perennial Grass"), avg_cov ~ BRTErel)
ols_plot_resid_hist(resid1.0) #normalish
ols_test_normality(resid1.0) #pass
#log trans
TBBRTEcov_func_c4pg_trans <- TBBRTEcov_func %>%
  filter(funct2 == "C4 Perennial Grass") %>%
  mutate(ln = log10(avg_cov + 0.001)) 

#linearity
plot(resid(c4pg_tb2), TBBRTEcov_func_c4pg_trans$avg_cov) #looks more linear
plot(c4pg_tb2) #not much pattern so indicates linearity

#homoscedascity
TBBRTEcov_func_c4pg_trans$res <- residuals(c4pg_tb2)
TBBRTEcov_func_c4pg_trans$abs_res <- abs(TBBRTEcov_func_c4pg_trans$res)
TBBRTEcov_func_c4pg_trans$abs_res2 <- TBBRTEcov_func_c4pg_trans$abs_res^2
levene_brte <- lm(data = TBBRTEcov_func_c4pg_trans, abs_res2 ~ BRTErel)
anova(levene_brte) #p = 0.7208 so > 0.05 so equal variance is met

#autocorrelation
acf(residuals(c4pg_tb2, retype = "normalized"))  #lines dont go outside CI horizontal lines, so not autocorrelated
pacf(residuals(c4pg_tb2, retype = "normalized"))


#Cactus
cc_tb2 <- lmerTest::lmer(avg_cov ~ BRTErel + (1|grad_num), data = subset(TBBRTEcov_func, funct2 == "Cactus"))
summary(cc_tb2) #not sig
anova(cc_tb2)

#check normality
resid1.0 <- lm(data = subset(TBBRTEcov_func, funct2 == "Cactus"), avg_cov ~ BRTErel)
ols_plot_resid_hist(resid1.0) #normalish
ols_test_normality(resid1.0) #pass some
#log trans
TBBRTEcov_func_cc_trans <- TBBRTEcov_func %>%
  filter(funct2 == "Cactus") %>%
  mutate(ln = log10(avg_cov + 0.001)) 

resid1.0 <- lm(data = TBBRTEcov_func_cc_trans, ln ~ BRTErel)
ols_plot_resid_hist(resid1.0) #normalish
ols_test_normality(resid1.0) #worse, stick with untrans

#linearity
plot(resid(cc_tb2), TBBRTEcov_func_cc_trans$avg_cov) #looks more linear
plot(cc_tb2) #not much pattern so indicates linearity

#homoscedascity
TBBRTEcov_func_cc_trans$res <- residuals(cc_tb2)
TBBRTEcov_func_cc_trans$abs_res <- abs(TBBRTEcov_func_cc_trans$res)
TBBRTEcov_func_cc_trans$abs_res2 <- TBBRTEcov_func_cc_trans$abs_res^2
levene_brte <- lm(data = TBBRTEcov_func_cc_trans, abs_res2 ~ BRTErel)
anova(levene_brte) #p = 0.7071 so > 0.05 so equal variance is met

#autocorrelation
acf(residuals(cc_tb2, retype = "normalized"))  #lines dont go outside CI horizontal lines, so not autocorrelated
pacf(residuals(cc_tb2, retype = "normalized"))


#Forb
cf_tb2 <- lmerTest::lmer(avg_cov ~ BRTErel + (1|grad_num), data = subset(TBBRTEcov_func, funct2 == "Forb"))
summary(cf_tb2) #not sig
anova(cf_tb2)

#check normality
resid1.0 <- lm(data = subset(TBBRTEcov_func, funct2 == "Forb"), avg_cov ~ BRTErel)
ols_plot_resid_hist(resid1.0) #normalish
ols_test_normality(resid1.0) #pass some
#log trans
TBBRTEcov_func_cf_trans <- TBBRTEcov_func %>%
  filter(funct2 == "Forb") %>%
  mutate(ln = log10(avg_cov + 0.001)) 

resid1.0 <- lm(data = TBBRTEcov_func_cf_trans, ln ~ BRTErel)
ols_plot_resid_hist(resid1.0) #normalish
ols_test_normality(resid1.0) #better

cf_tb2a <- lmerTest::lmer(ln ~ BRTErel + (1|grad_num), data = TBBRTEcov_func_cf_trans)
summary(cf_tb2a) #not sig
anova(cf_tb2a)

#linearity
plot(resid(cf_tb2a), TBBRTEcov_func_cf_trans$ln) #looks more linear
plot(cf_tb2a) #not much pattern so indicates linearity

#homoscedascity
TBBRTEcov_func_cf_trans$res <- residuals(cf_tb2a)
TBBRTEcov_func_cf_trans$abs_res <- abs(TBBRTEcov_func_cf_trans$res)
TBBRTEcov_func_cf_trans$abs_res2 <- TBBRTEcov_func_cf_trans$abs_res^2
levene_brte <- lm(data = TBBRTEcov_func_cf_trans, abs_res2 ~ BRTErel)
anova(levene_brte) #p = 0.4564 so > 0.05 so equal variance is met

#autocorrelation
acf(residuals(cf_tb2a, retype = "normalized"))  #lines dont go outside CI horizontal lines, so not autocorrelated
pacf(residuals(cf_tb2a, retype = "normalized"))


#Shrub
cs_tb2 <- lmerTest::lmer(avg_cov ~ BRTErel + (1|grad_num), data = subset(TBBRTEcov_func, funct2 == "Sub-Shrub/Shrub"))
summary(cs_tb2) #not sig
anova(cs_tb2)

#check normality
resid1.0 <- lm(data = subset(TBBRTEcov_func, funct2 == "Sub-Shrub/Shrub"), avg_cov ~ BRTErel)
ols_plot_resid_hist(resid1.0) #right skew
ols_test_normality(resid1.0) #fail
#log trans
TBBRTEcov_func_cs_trans <- TBBRTEcov_func %>%
  filter(funct2 == "Sub-Shrub/Shrub") %>%
  mutate(ln = log10(avg_cov + 0.001)) 

resid1.0 <- lm(data = TBBRTEcov_func_cs_trans, ln ~ BRTErel)
ols_plot_resid_hist(resid1.0) #normalish
ols_test_normality(resid1.0) #better

cs_tb2a <- lmerTest::lmer(ln ~ BRTErel + (1|grad_num), data = TBBRTEcov_func_cs_trans)
summary(cs_tb2a) #not sig
anova(cs_tb2a)

#linearity
plot(resid(cs_tb2a), TBBRTEcov_func_cs_trans$ln) #looks more linear
plot(cs_tb2a) #not much pattern so indicates linearity

#homoscedascity
TBBRTEcov_func_cs_trans$res <- residuals(cs_tb2a)
TBBRTEcov_func_cs_trans$abs_res <- abs(TBBRTEcov_func_cs_trans$res)
TBBRTEcov_func_cs_trans$abs_res2 <- TBBRTEcov_func_cs_trans$abs_res^2
levene_brte <- lm(data = TBBRTEcov_func_cs_trans, abs_res2 ~ BRTErel)
anova(levene_brte) #p = 0.3336 so > 0.05 so equal variance is met

#autocorrelation
acf(residuals(cs_tb2a, retype = "normalized"))  #lines dont go outside CI horizontal lines, so not autocorrelated
pacf(residuals(cs_tb2a, retype = "normalized"))





#BRTE Insects
TBinsID_guild_BRTE3 <- TBinsID_guild_BRTE %>%
  group_by(plot, guild) %>%
  summarise(avg = sum(avg_n)) %>%
  ungroup()

TBBRTEinsfunc <- TBinsID_guild_BRTE3 %>%
  full_join(TB_commBRTE) %>%
  drop_na(guild)

TB_BRTE_func_reg_ins <- ggplot(data = TBBRTEinsfunc, aes(x = BRTErel, y = avg, color = guild, shape = guild)) +
  #geom_point(size = 3) +
  geom_smooth(data = subset(TBBRTEinsfunc, guild == "Leaf Chewing Herbivore"), aes(group = guild), method = "lm", se = FALSE, linetype = "dotted") +
  geom_smooth(data = subset(TBBRTEinsfunc, guild == "Predator"), aes(group = guild), method = "lm", se = FALSE, linetype = "dotted") +
  geom_smooth(data = subset(TBBRTEinsfunc, guild == "Sap Sucking Herbivore"), aes(group = guild), method = "lm", se = FALSE, linetype = "dotted") +
  geom_smooth(data = subset(TBBRTEinsfunc, guild == "Parasitoid"), aes(group = guild), method = "lm", se = FALSE, linetype = "dotted") +
  scale_colour_manual(values = cbPalette_rac) +
  #scale_shape_manual(values = c(19, 18, 17, 15, 8, 7)) +
  labs(x = "BRTE Cover (%)", y = "Abundance (count)", title = "WY BRTE", color = "Functional Group", shape = "Functional Group") + 
  theme(plot.title = element_text(vjust = - 7, hjust = 0.07, size = 16)) 

#Leaf chewing
lch_tb2 <- lmerTest::lmer(avg ~ BRTErel + (1|grad_num), data = subset(TBBRTEinsfunc, guild == "Leaf Chewing Herbivore"))
summary(lch_tb2) #no sig
anova(lch_tb2)

#check normality
resid1.0 <- lm(data = subset(TBBRTEinsfunc, guild == "Leaf Chewing Herbivore"), avg ~ BRTErel)
ols_plot_resid_hist(resid1.0) #right skew
ols_test_normality(resid1.0) #fail
#log trans
TBBRTEins_lchfunc <- TBBRTEinsfunc %>%
  filter(guild == "Leaf Chewing Herbivore") %>%
  mutate(ln = log10(avg + 0.001)) 

resid1.0 <- lm(data = TBBRTEins_lchfunc, ln ~ BRTErel)
ols_plot_resid_hist(resid1.0) #normalish
ols_test_normality(resid1.0) #better

lch_tb2a <- lmerTest::lmer(ln ~ BRTErel + (1|grad_num), data = TBBRTEins_lchfunc)
summary(lch_tb2a) #not sig
anova(lch_tb2a)

#linearity
plot(resid(lch_tb2a), TBBRTEins_lchfunc$ln) #looks more linear
plot(lch_tb2a) #not much pattern so indicates linearity

#homoscedascity
TBBRTEins_lchfunc$res <- residuals(lch_tb2a)
TBBRTEins_lchfunc$abs_res <- abs(TBBRTEins_lchfunc$res)
TBBRTEins_lchfunc$abs_res2 <- TBBRTEins_lchfunc$abs_res^2
levene_brte <- lm(data = TBBRTEins_lchfunc, abs_res2 ~ BRTErel)
anova(levene_brte) #p = 0.9628 so > 0.05 so equal variance is met

#autocorrelation
acf(residuals(lch_tb2a, retype = "normalized"))  #lines dont go outside CI horizontal lines, so not autocorrelated
pacf(residuals(lch_tb2a, retype = "normalized"))


#Predator
pre_tb2 <- lmerTest::lmer(avg ~ BRTErel + (1|grad_num), data = subset(TBBRTEinsfunc, guild == "Predator"))
summary(pre_tb2) #no sig
anova(pre_tb2)

#check normality
resid1.0 <- lm(data = subset(TBBRTEinsfunc, guild == "Predator"), avg ~ BRTErel)
ols_plot_resid_hist(resid1.0) #right skew
ols_test_normality(resid1.0) #fail
#log trans
TBBRTEins_prefunc <- TBBRTEinsfunc %>%
  filter(guild == "Predator") %>%
  mutate(ln = log10(avg + 0.001)) 

resid1.0 <- lm(data = TBBRTEins_prefunc, ln ~ BRTErel)
ols_plot_resid_hist(resid1.0) #normalish
ols_test_normality(resid1.0) #better

pre_tb2a <- lmerTest::lmer(ln ~ BRTErel + (1|grad_num), data = TBBRTEins_prefunc)
summary(pre_tb2a) #not sig
anova(pre_tb2a)

#linearity
plot(resid(pre_tb2a), TBBRTEins_prefunc$ln) #looks more linear
plot(pre_tb2a) #not much pattern so indicates linearity

#homoscedascity
TBBRTEins_prefunc$res <- residuals(pre_tb2a)
TBBRTEins_prefunc$abs_res <- abs(TBBRTEins_prefunc$res)
TBBRTEins_prefunc$abs_res2 <- TBBRTEins_prefunc$abs_res^2
levene_brte <- lm(data = TBBRTEins_prefunc, abs_res2 ~ BRTErel)
anova(levene_brte) #p = 0.2331 so > 0.05 so equal variance is met

#autocorrelation
acf(residuals(pre_tb2a, retype = "normalized"))  #lines dont go outside CI horizontal lines, so not autocorrelated
pacf(residuals(pre_tb2a, retype = "normalized"))


#Sap sucking
ssh_tb2 <- lmerTest::lmer(avg ~ BRTErel + (1|grad_num), data = subset(TBBRTEinsfunc, guild == "Sap Sucking Herbivore"))
summary(ssh_tb2) #no sig
anova(ssh_tb2)

#check normality
resid1.0 <- lm(data = subset(TBBRTEinsfunc, guild == "Sap Sucking Herbivore"), avg ~ BRTErel)
ols_plot_resid_hist(resid1.0) #right skew
ols_test_normality(resid1.0) #fail
#log trans
TBBRTEins_sshfunc <- TBBRTEinsfunc %>%
  filter(guild == "Sap Sucking Herbivore") %>%
  mutate(ln = log10(avg + 0.001)) 

resid1.0 <- lm(data = TBBRTEins_sshfunc, ln ~ BRTErel)
ols_plot_resid_hist(resid1.0) #normalish
ols_test_normality(resid1.0) #pass

ssh_tb2a <- lmerTest::lmer(ln ~ BRTErel + (1|grad_num), data = TBBRTEins_sshfunc)
summary(ssh_tb2a) #not sig
anova(ssh_tb2a)

#linearity
plot(resid(ssh_tb2a), TBBRTEins_sshfunc$ln) #looks more linear
plot(ssh_tb2a) #not much pattern so indicates linearity

#homoscedascity
TBBRTEins_sshfunc$res <- residuals(ssh_tb2a)
TBBRTEins_sshfunc$abs_res <- abs(TBBRTEins_sshfunc$res)
TBBRTEins_sshfunc$abs_res2 <- TBBRTEins_sshfunc$abs_res^2
levene_brte <- lm(data = TBBRTEins_sshfunc, abs_res2 ~ BRTErel)
anova(levene_brte) #p = 0.9781 so > 0.05 so equal variance is met

#autocorrelation
acf(residuals(ssh_tb2a, retype = "normalized"))  #lines dont go outside CI horizontal lines, so not autocorrelated
pacf(residuals(ssh_tb2a, retype = "normalized"))



par_tb2 <- lm(avg ~ BRTErel, data = subset(TBBRTEinsfunc, guild == "Parasitoid"))
summary(par_tb2) #no sig
anova(par_tb2)
par_tb2 <- lmerTest::lmer(avg ~ BRTErel + (1|grad_num), data = subset(TBBRTEinsfunc, guild == "Parasitoid"))
summary(par_tb2) #cant be run

#check normality
resid1.0 <- lm(data = subset(TBBRTEinsfunc, guild == "Parasitoid"), avg ~ BRTErel)
ols_plot_resid_hist(resid1.0) #right skew
ols_test_normality(resid1.0) #fail some, use untrans
#log trans
TBBRTEins_parfunc <- TBBRTEinsfunc %>%
  filter(guild == "Parasitoid") %>%
  mutate(ln = log10(avg + 0.001)) 



pneh_tb2 <- lm(avg ~ BRTErel, data = subset(TBBRTEinsfunc, guild == "Pollen/Nectar Eating Herbivore"))
summary(pneh_tb2) #not enough data to include
pneh_tb2 <- lmerTest::lmer(avg ~ BRTErel + (1|grad_num), data = subset(TBBRTEinsfunc, guild == "Pollen/Nectar Eating Herbivore"))
summary(pneh_tb2) #not enough data to include


oh_tb2 <- lm(avg ~ BRTErel, data = subset(TBBRTEinsfunc, guild == "Other Herbivore"))
summary(oh_tb2) #not enough data to include
oh_tb2 <- lmerTest::lmer(avg ~ BRTErel + (1|grad_num), data = subset(TBBRTEinsfunc, guild == "Other Herbivore"))
summary(oh_tb2) #not enough data to include





#BRTE Microbes
#TBmicfunc2 <- BRTETBmic_func_long2 %>%
dplyr::select(func) %>%
  unique()
#write.csv(TBmicfunc2, "tbmic2.csv")

BRTETBmic_func_long2 <- TBmic_func_long %>%
  filter(invasive_type == "BRTE") %>%
  full_join(TB_BRTE_2019) %>%
  group_by(plot, func) %>%
  summarise(avg = sum(count)) %>%
  ungroup() %>%
  drop_na(func) %>%
  full_join(TB_commBRTE)

BRTETBmic_func_long3 <- BRTETBmic_func_long2 %>%
  filter(func != "other")

#
BRTETBmic_func_long4 <- BRTETBmic_func_long2 %>%
  group_by(func) %>%
  summarise(avg2 = mean(avg)) %>%
  ungroup()
#top 7 most abundant functions (not other): chemoheterotrophy, aerobic_chemoheterotrophy, nitrate_reduction, nitrification, aerobic_ammonia_oxidation, nitrogen_fixation, manganese_oxidation 

#only whats significant
TB_BRTE_func_reg_mic <- ggplot(data = BRTETBmic_func_long3, aes(x = BRTErel, y = avg, color = func)) +
  #geom_point(size = 3) +
  geom_smooth(data = subset(BRTETBmic_func_long3, func == "anoxygenic_photoautotrophy"), aes(group = func), method = "lm", se = TRUE) +
  geom_smooth(data = subset(BRTETBmic_func_long3, func == "anoxygenic_photoautotrophy_S_oxidizing"), aes(group = func), method = "lm", se = TRUE) +
  geom_smooth(data = subset(BRTETBmic_func_long3, func == "aromatic_compound_degradation"), aes(group = func), method = "lm", se = TRUE) +
  geom_smooth(data = subset(BRTETBmic_func_long3, func == "denitrification"), aes(group = func), method = "lm", se = TRUE) +
  geom_smooth(data = subset(BRTETBmic_func_long3, func == "methanol_oxidation"), aes(group = func), method = "lm", se = TRUE) +
  geom_smooth(data = subset(BRTETBmic_func_long3, func == "methylotrophy"), aes(group = func), method = "lm", se = TRUE) +
  geom_smooth(data = subset(BRTETBmic_func_long3, func == "nitrate_denitrification"), aes(group = func), method = "lm", se = TRUE) +
  geom_smooth(data = subset(BRTETBmic_func_long3, func == "nitrate_respiration"), aes(group = func), method = "lm", se = TRUE) +
  geom_smooth(data = subset(BRTETBmic_func_long3, func == "nitrite_denitrification"), aes(group = func), method = "lm", se = TRUE) +
  geom_smooth(data = subset(BRTETBmic_func_long3, func == "nitrite_respiration"), aes(group = func), method = "lm", se = TRUE) +
  geom_smooth(data = subset(BRTETBmic_func_long3, func == "nitrogen_respiration"), aes(group = func), method = "lm", se = TRUE) +
  geom_smooth(data = subset(BRTETBmic_func_long3, func == "nitrous_oxide_denitrification"), aes(group = func), method = "lm", se = TRUE) +
  geom_smooth(data = subset(BRTETBmic_func_long3, func == "oxygenic_photoautotrophy"), aes(group = func), method = "lm", se = TRUE) +
  geom_smooth(data = subset(BRTETBmic_func_long3, func == "photoautotrophy"), aes(group = func), method = "lm", se = TRUE) +
  geom_smooth(data = subset(BRTETBmic_func_long3, func == "photoheterotrophy"), aes(group = func), method = "lm", se = TRUE) +
  geom_smooth(data = subset(BRTETBmic_func_long3, func == "photosynthetic_cyanobacteria"), aes(group = func), method = "lm", se = TRUE) +
  geom_smooth(data = subset(BRTETBmic_func_long3, func == "phototrophy"), aes(group = func), method = "lm", se = TRUE) +
  geom_smooth(data = subset(BRTETBmic_func_long3, func == "predatory_or_exoparasitic"), aes(group = func), method = "lm", se = TRUE, linetype = "dashed") +
  geom_smooth(data = subset(BRTETBmic_func_long3, func == "ureolysis"), aes(group = func), method = "lm", se = TRUE) +
  #geom_smooth(data = subset(BRARFKmic_func_long3, func != "oxygenic_photoautotrophy"), aes(group = func), method = "lm", se = FALSE, linetype = "dotted") +
  #scale_colour_manual(values = cbPalette_rac) +
  #scale_shape_manual(values = c(19, 18, 17, 15, 8, 7)) +
  ylim(0, 2500) +
  labs(x = "BRTE Cover (%)", y = "Abundance (count)", title = "WY BRTE", color = "Functional Group", shape = "Functional Group") + 
  theme(plot.title = element_text(vjust = - 7, hjust = 0.07, size = 16)) 

#
FK_BRAR_func_reg_mic2 <- ggplot(data = BRARFKmic_func_long3, aes(x = BRARrel, y = avg, color = func)) +
  #geom_point(size = 3) +
  geom_smooth(data = subset(BRARFKmic_func_long3, func == "chemoheterotrophy"), aes(group = func), method = "lm", se = FALSE, linetype = "dotted", position = position_dodge(0.5)) +
  geom_smooth(data = subset(BRARFKmic_func_long3, func == "aerobic_chemoheterotrophy"), aes(group = func), method = "lm", se = FALSE, linetype = "dotted", position = position_dodge(0.5)) +
  geom_smooth(data = subset(BRARFKmic_func_long3, func == "nitrate_reduction"), aes(group = func), method = "lm", se = FALSE, linetype = "dotted") +
  geom_smooth(data = subset(BRARFKmic_func_long3, func == "nitrification"), aes(group = func), method = "lm", se = FALSE, linetype = "dotted") +
  geom_smooth(data = subset(BRARFKmic_func_long3, func == "aerobic_ammonia_oxidation"), aes(group = func), method = "lm", se = FALSE, linetype = "dotted") +
  geom_smooth(data = subset(BRARFKmic_func_long3, func == "nitrogen_fixation"), aes(group = func), method = "lm", se = FALSE, linetype = "dotted") +
  geom_smooth(data = subset(BRARFKmic_func_long3, func == "manganese_oxidation"), aes(group = func), method = "lm", se = FALSE, linetype = "dotted") +
  scale_colour_manual(values = cbPalette_rac) +
  #scale_shape_manual(values = c(19, 18, 17, 15, 8, 7)) +
  labs(x = "BRAR Cover (%)", y = "Abundance (count)", title = "MT BRAR", color = "Functional Group", shape = "Functional Group") + 
  theme(plot.title = element_text(vjust = - 7, hjust = 0.07, size = 16)) 
#

#Aerobic ammon oxid
tb_mic1a <- lmerTest::lmer(avg ~ BRTErel + (1|grad_num), data = subset(BRTETBmic_func_long2, func == "aerobic_ammonia_oxidation"))
summary(tb_mic1a) #not sig
anova(tb_mic1a)

#check normality
resid1.0 <- lm(data = subset(BRTETBmic_func_long2, func == "aerobic_ammonia_oxidation"), avg ~ BRTErel)
ols_plot_resid_hist(resid1.0) #right skew
ols_test_normality(resid1.0) #fail some, use untrans
#log trans
BRTETBmic_func_1 <- BRTETBmic_func_long2 %>%
  filter(func == "aerobic_ammonia_oxidation") %>%
  mutate(ln = log10(avg + 0.001)) 

resid1.0 <- lm(data = BRTETBmic_func_1, ln ~ BRTErel)
ols_plot_resid_hist(resid1.0) 
ols_test_normality(resid1.0) 

tb_mic1ab <- lmerTest::lmer(ln ~ BRTErel + (1|grad_num), data = BRTETBmic_func_1)
summary(tb_mic1ab) #not sig
anova(tb_mic1ab)


#Aerobic chemo
tb_mic2a <- lmerTest::lmer(avg ~ BRTErel + (1|grad_num), data = subset(BRTETBmic_func_long2, func == "aerobic_chemoheterotrophy"))
summary(tb_mic2a) #not sig
anova(tb_mic2a)

#check normality
resid1.0 <- lm(data = subset(BRTETBmic_func_long2, func == "aerobic_chemoheterotrophy"), avg ~ BRTErel)
ols_plot_resid_hist(resid1.0) 
ols_test_normality(resid1.0) 
#log trans
BRTETBmic_func_2 <- BRTETBmic_func_long2 %>%
  filter(func == "aerobic_chemoheterotrophy") %>%
  mutate(ln = log10(avg + 0.001)) 

resid1.0 <- lm(data = BRTETBmic_func_2, ln ~ BRTErel)
ols_plot_resid_hist(resid1.0) #right skew
ols_test_normality(resid1.0) #fail some, use untrans

tb_mic2ab <- lmerTest::lmer(ln ~ BRTErel + (1|grad_num), data = BRTETBmic_func_2)
summary(tb_mic2ab) #p = 0.07886; slope = -0.001961
anova(tb_mic2ab)


#Aerobic nitrite oxid
tb_mic3a <- lmerTest::lmer(avg ~ BRTErel + (1|grad_num), data = subset(BRTETBmic_func_long2, func == "aerobic_nitrite_oxidation"))
summary(tb_mic3a) #not sig
anova(tb_mic3a)

#check normality
resid1.0 <- lm(data = subset(BRTETBmic_func_long2, func == "aerobic_nitrite_oxidation"), avg ~ BRTErel)
ols_plot_resid_hist(resid1.0) 
ols_test_normality(resid1.0) 

#Animal par 
tb_mic4a <- lmerTest::lmer(avg ~ BRTErel + (1|grad_num), data = subset(BRTETBmic_func_long2, func == "animal_parasites_or_symbionts"))
summary(tb_mic4a) #not sig
anova(tb_mic4a)

#check normality
resid1.0 <- lm(data = subset(BRTETBmic_func_long2, func == "animal_parasites_or_symbionts"), avg ~ BRTErel)
ols_plot_resid_hist(resid1.0) 
ols_test_normality(resid1.0) 
#log trans
BRTETBmic_func_4 <- BRTETBmic_func_long2 %>%
  filter(func == "animal_parasites_or_symbionts") %>%
  mutate(ln = log10(avg + 0.001)) 

resid1.0 <- lm(data = BRTETBmic_func_4, ln ~ BRTErel)
ols_plot_resid_hist(resid1.0) 
ols_test_normality(resid1.0) 

tb_mic4ab <- lmerTest::lmer(ln ~ BRTErel + (1|grad_num), data = BRTETBmic_func_4)
summary(tb_mic4ab) #p = 0.07886; slope = -0.001961
anova(tb_mic4ab)

#Anoxy photo
tb_mic5a <- lmerTest::lmer(avg ~ BRTErel + (1|grad_num), data = subset(BRTETBmic_func_long2, func == "anoxygenic_photoautotrophy"))
summary(tb_mic5a) #p = 0.0488; slope = -4.724
anova(tb_mic5a)

#check normality
resid1.0 <- lm(data = subset(BRTETBmic_func_long2, func == "anoxygenic_photoautotrophy"), avg ~ BRTErel)
ols_plot_resid_hist(resid1.0) 
ols_test_normality(resid1.0) 
#log trans
BRTETBmic_func_5 <- BRTETBmic_func_long2 %>%
  filter(func == "anoxygenic_photoautotrophy") %>%
  mutate(ln = log10(avg + 0.001)) 

resid1.0 <- lm(data = BRTETBmic_func_5, ln ~ BRTErel)
ols_plot_resid_hist(resid1.0) 
ols_test_normality(resid1.0) 

tb_mic5ab <- lmerTest::lmer(ln ~ BRTErel + (1|grad_num), data = BRTETBmic_func_5)
summary(tb_mic5ab) #p = 0.009843; slope = -0.003979
anova(tb_mic5ab)


#Anoxy photo S
tb_mic6a <- lmerTest::lmer(avg ~ BRTErel + (1|grad_num), data = subset(BRTETBmic_func_long2, func == "anoxygenic_photoautotrophy_S_oxidizing"))
summary(tb_mic6a) #p = 0.0488; slope = -4.724
anova(tb_mic6a)

#Aromatic comp deg
tb_mic7a <- lmerTest::lmer(avg ~ BRTErel + (1|grad_num), data = subset(BRTETBmic_func_long2, func == "aromatic_compound_degradation"))
summary(tb_mic7a) #p = 0.00523; slope = -8.759 
anova(tb_mic7a)

#check normality
resid1.0 <- lm(data = subset(BRTETBmic_func_long2, func == "aromatic_compound_degradation"), avg ~ BRTErel)
ols_plot_resid_hist(resid1.0) 
ols_test_normality(resid1.0) 


#Chemohetro
tb_mic8a <- lmerTest::lmer(avg ~ BRTErel + (1|grad_num), data = subset(BRTETBmic_func_long2, func == "chemoheterotrophy"))
summary(tb_mic8a) #not sig 
anova(tb_mic8a)

#check normality
resid1.0 <- lm(data = subset(BRTETBmic_func_long2, func == "chemoheterotrophy"), avg ~ BRTErel)
ols_plot_resid_hist(resid1.0) 
ols_test_normality(resid1.0) 
#log trans
BRTETBmic_func_8 <- BRTETBmic_func_long2 %>%
  filter(func == "chemoheterotrophy") %>%
  mutate(ln = log10(avg + 0.001)) 

resid1.0 <- lm(data = BRTETBmic_func_8, ln ~ BRTErel)
ols_plot_resid_hist(resid1.0) 
ols_test_normality(resid1.0) 

tb_mic8ab <- lmerTest::lmer(ln ~ BRTErel + (1|grad_num), data = BRTETBmic_func_8)
summary(tb_mic8ab) #p = 0.06063; slope = -0.002043
anova(tb_mic8ab)

#Chitin
tb_mic9a <- lmerTest::lmer(avg ~ BRTErel + (1|grad_num), data = subset(BRTETBmic_func_long2, func == "chitinolysis"))
summary(tb_mic9a) #not sig
anova(tb_mic9a)

#check normality
resid1.0 <- lm(data = subset(BRTETBmic_func_long2, func == "chitinolysis"), avg ~ BRTErel)
ols_plot_resid_hist(resid1.0) 
ols_test_normality(resid1.0) 
#log trans
BRTETBmic_func_9 <- BRTETBmic_func_long2 %>%
  filter(func == "chitinolysis") %>%
  mutate(ln = log10(avg + 0.001)) 

resid1.0 <- lm(data = BRTETBmic_func_9, ln ~ BRTErel)
ols_plot_resid_hist(resid1.0) 
ols_test_normality(resid1.0) 

tb_mic9ab <- lmerTest::lmer(ln ~ BRTErel + (1|grad_num), data = BRTETBmic_func_9)
summary(tb_mic9ab) #p = 0.06063; slope = -0.002043
anova(tb_mic9ab)

#Dark oxid S
tb_mic10a <- lmerTest::lmer(avg ~ BRTErel + (1|grad_num), data = subset(BRTETBmic_func_long2, func == "dark_oxidation_of_sulfur_compounds"))
summary(tb_mic10a) #not sig
anova(tb_mic10a)

#check normality
resid1.0 <- lm(data = subset(BRTETBmic_func_long2, func == "dark_oxidation_of_sulfur_compounds"), avg ~ BRTErel)
ols_plot_resid_hist(resid1.0) 
ols_test_normality(resid1.0) 


#Denitrificaiton
tb_mic11a <- lmerTest::lmer(avg ~ BRTErel + (1|grad_num), data = subset(BRTETBmic_func_long2, func == "denitrification"))
summary(tb_mic11a) #p = 0.0499; slope = -4.706
anova(tb_mic11a)

#check normality
resid1.0 <- lm(data = subset(BRTETBmic_func_long2, func == "denitrification"), avg ~ BRTErel)
ols_plot_resid_hist(resid1.0) 
ols_test_normality(resid1.0) 
#log trans
BRTETBmic_func_11 <- BRTETBmic_func_long2 %>%
  filter(func == "denitrification") %>%
  mutate(ln = log10(avg + 0.001)) 

resid1.0 <- lm(data = BRTETBmic_func_11, ln ~ BRTErel)
ols_plot_resid_hist(resid1.0) 
ols_test_normality(resid1.0) 

tb_mic11ab <- lmerTest::lmer(ln ~ BRTErel + (1|grad_num), data = BRTETBmic_func_11)
summary(tb_mic11ab) #p = 0.0105; slope = -0.003954
anova(tb_mic11ab)


#Fermentation
tb_mic12a <- lmerTest::lmer(avg ~ BRTErel + (1|grad_num), data = subset(BRTETBmic_func_long2, func == "fermentation"))
summary(tb_mic12a) #not sig
anova(tb_mic12a)

#check normality
resid1.0 <- lm(data = subset(BRTETBmic_func_long2, func == "fermentation"), avg ~ BRTErel)
ols_plot_resid_hist(resid1.0) 
ols_test_normality(resid1.0) 
#log trans
BRTETBmic_func_12 <- BRTETBmic_func_long2 %>%
  filter(func == "fermentation") %>%
  mutate(ln = log10(avg + 0.001)) 

resid1.0 <- lm(data = BRTETBmic_func_12, ln ~ BRTErel)
ols_plot_resid_hist(resid1.0) 
ols_test_normality(resid1.0) 

tb_mic12ab <- lmerTest::lmer(ln ~ BRTErel + (1|grad_num), data = BRTETBmic_func_12)
summary(tb_mic12ab) #p = 0.3361; slope = -0.002060
anova(tb_mic12ab)

#Mang oxid
tb_mic13a <- lmerTest::lmer(avg ~ BRTErel + (1|grad_num), data = subset(BRTETBmic_func_long2, func == "manganese_oxidation"))
summary(tb_mic13a) #not sig
anova(tb_mic13a)

#check normality
resid1.0 <- lm(data = subset(BRTETBmic_func_long2, func == "manganese_oxidation"), avg ~ BRTErel)
ols_plot_resid_hist(resid1.0) 
ols_test_normality(resid1.0) 
#log trans
BRTETBmic_func_13 <- BRTETBmic_func_long2 %>%
  filter(func == "manganese_oxidation") %>%
  mutate(ln = log10(avg + 0.001)) 

resid1.0 <- lm(data = BRTETBmic_func_13, ln ~ BRTErel)
ols_plot_resid_hist(resid1.0) 
ols_test_normality(resid1.0) 

tb_mic13ab <- lmerTest::lmer(ln ~ BRTErel + (1|grad_num), data = BRTETBmic_func_13)
summary(tb_mic13ab) #p = 0.3324; slope = -0.004080
anova(tb_mic13ab)

#Meth oxid
tb_mic14a <- lmerTest::lmer(avg ~ BRTErel + (1|grad_num), data = subset(BRTETBmic_func_long2, func == "methanol_oxidation"))
summary(tb_mic14a) #p = 0.0102; slope = -4.117 
anova(tb_mic14a)

#check normality
resid1.0 <- lm(data = subset(BRTETBmic_func_long2, func == "methanol_oxidation"), avg ~ BRTErel)
ols_plot_resid_hist(resid1.0) 
ols_test_normality(resid1.0) 
#log trans
BRTETBmic_func_14 <- BRTETBmic_func_long2 %>%
  filter(func == "methanol_oxidation") %>%
  mutate(ln = log10(avg + 0.001)) 

resid1.0 <- lm(data = BRTETBmic_func_14, ln ~ BRTErel)
ols_plot_resid_hist(resid1.0) 
ols_test_normality(resid1.0) 

tb_mic14ab <- lmerTest::lmer(ln ~ BRTErel + (1|grad_num), data = BRTETBmic_func_14)
summary(tb_mic14ab) #p = 0.00451; slope = -0.005561
anova(tb_mic14ab)

#Methylo
tb_mic15a <- lmerTest::lmer(avg ~ BRTErel + (1|grad_num), data = subset(BRTETBmic_func_long2, func == "methylotrophy"))
summary(tb_mic15a) #p = 0.0102; slope = -4.117
anova(tb_mic15a)

#Nitrate denit
tb_mic16a <- lmerTest::lmer(avg ~ BRTErel + (1|grad_num), data = subset(BRTETBmic_func_long2, func == "nitrate_denitrification"))
summary(tb_mic16a) #p = 0.0499; slope = -4.706
anova(tb_mic16a)

#Nitrate reduc
tb_mic17a <- lmerTest::lmer(avg ~ BRTErel + (1|grad_num), data = subset(BRTETBmic_func_long2, func == "nitrate_reduction"))
summary(tb_mic17a) #not sig
anova(tb_mic17a)

#check normality
resid1.0 <- lm(data = subset(BRTETBmic_func_long2, func == "nitrate_reduction"), avg ~ BRTErel)
ols_plot_resid_hist(resid1.0) 
ols_test_normality(resid1.0) 
#log trans
BRTETBmic_func_17 <- BRTETBmic_func_long2 %>%
  filter(func == "nitrate_reduction") %>%
  mutate(ln = log10(avg + 0.001)) 

resid1.0 <- lm(data = BRTETBmic_func_17, ln ~ BRTErel)
ols_plot_resid_hist(resid1.0) 
ols_test_normality(resid1.0) 

tb_mic17ab <- lmerTest::lmer(ln ~ BRTErel + (1|grad_num), data = BRTETBmic_func_17)
summary(tb_mic17ab) #p = 0.1333; slope = -0.002524
anova(tb_mic17ab)

#Nitrate resp
tb_mic18a <- lmerTest::lmer(avg ~ BRTErel + (1|grad_num), data = subset(BRTETBmic_func_long2, func == "nitrate_respiration"))
summary(tb_mic18a) #p = 0.0488; slope = -4.722
anova(tb_mic18a)

#check normality
resid1.0 <- lm(data = subset(BRTETBmic_func_long2, func == "nitrate_respiration"), avg ~ BRTErel)
ols_plot_resid_hist(resid1.0) 
ols_test_normality(resid1.0) 
#log trans
BRTETBmic_func_18 <- BRTETBmic_func_long2 %>%
  filter(func == "nitrate_respiration") %>%
  mutate(ln = log10(avg + 0.001)) 

resid1.0 <- lm(data = BRTETBmic_func_18, ln ~ BRTErel)
ols_plot_resid_hist(resid1.0) 
ols_test_normality(resid1.0) 

tb_mic18ab <- lmerTest::lmer(ln ~ BRTErel + (1|grad_num), data = BRTETBmic_func_18)
summary(tb_mic18ab) #p = 0.0104; slope = -0.003955
anova(tb_mic18ab)

#Nitrification
tb_mic19a <- lmerTest::lmer(avg ~ BRTErel + (1|grad_num), data = subset(BRTETBmic_func_long2, func == "nitrification"))
summary(tb_mic19a) #not sig
anova(tb_mic19a)

#check normality
resid1.0 <- lm(data = subset(BRTETBmic_func_long2, func == "nitrification"), avg ~ BRTErel)
ols_plot_resid_hist(resid1.0) 
ols_test_normality(resid1.0) 
#log trans
BRTETBmic_func_19 <- BRTETBmic_func_long2 %>%
  filter(func == "nitrification") %>%
  mutate(ln = log10(avg + 0.001)) 

resid1.0 <- lm(data = BRTETBmic_func_19, ln ~ BRTErel)
ols_plot_resid_hist(resid1.0) 
ols_test_normality(resid1.0) 

tb_mic19ab <- lmerTest::lmer(ln ~ BRTErel + (1|grad_num), data = BRTETBmic_func_19)
summary(tb_mic19ab) #p = 0.145; slope = -0.002207
anova(tb_mic19ab)

#Nitrate denitif
tb_mic20a <- lmerTest::lmer(avg ~ BRTErel + (1|grad_num), data = subset(BRTETBmic_func_long2, func == "nitrite_denitrification"))
summary(tb_mic20a) #p = 0.0499; slope = -4.706
anova(tb_mic20a)

#Nitrate resp
tb_mic21a <- lmerTest::lmer(avg ~ BRTErel + (1|grad_num), data = subset(BRTETBmic_func_long2, func == "nitrite_respiration"))
summary(tb_mic21a) #p = 0.0499; slope = -4.706
anova(tb_mic21a)

#Nitrogen fix
tb_mic22a <- lmerTest::lmer(avg ~ BRTErel + (1|grad_num), data = subset(BRTETBmic_func_long2, func == "nitrogen_fixation"))
summary(tb_mic22a) #not sig
anova(tb_mic22a)

#check normality
resid1.0 <- lm(data = subset(BRTETBmic_func_long2, func == "nitrogen_fixation"), avg ~ BRTErel)
ols_plot_resid_hist(resid1.0) 
ols_test_normality(resid1.0) 
#log trans
BRTETBmic_func_22 <- BRTETBmic_func_long2 %>%
  filter(func == "nitrogen_fixation") %>%
  mutate(ln = log10(avg + 0.001)) 

resid1.0 <- lm(data = BRTETBmic_func_22, ln ~ BRTErel)
ols_plot_resid_hist(resid1.0) 
ols_test_normality(resid1.0) 

tb_mic22ab <- lmerTest::lmer(ln ~ BRTErel + (1|grad_num), data = BRTETBmic_func_22)
summary(tb_mic22ab) #p = 0.473; slope = -0.0006978
anova(tb_mic22ab)

#Nitrogen resp
tb_mic23a <- lmerTest::lmer(avg ~ BRTErel + (1|grad_num), data = subset(BRTETBmic_func_long2, func == "nitrogen_respiration"))
summary(tb_mic23a) #p = 0.0488; slope = -4.722
anova(tb_mic23a)

#Nitrous oxide denit
tb_mic24a <- lmerTest::lmer(avg ~ BRTErel + (1|grad_num), data = subset(BRTETBmic_func_long2, func == "nitrous_oxide_denitrification"))
summary(tb_mic24a) #p = 0.0499; slope = -4.706
anova(tb_mic24a)

#Nonphoto cyano
tb_mic25a <- lmerTest::lmer(avg ~ BRTErel + (1|grad_num), data = subset(BRTETBmic_func_long2, func == "nonphotosynthetic_cyanobacteria"))
summary(tb_mic25a) #not sig
anova(tb_mic25a)

#check normality
resid1.0 <- lm(data = subset(BRTETBmic_func_long2, func == "nonphotosynthetic_cyanobacteria"), avg ~ BRTErel)
ols_plot_resid_hist(resid1.0) 
ols_test_normality(resid1.0) 
#log trans
BRTETBmic_func_25 <- BRTETBmic_func_long2 %>%
  filter(func == "nonphotosynthetic_cyanobacteria") %>%
  mutate(ln = log10(avg + 0.001)) 

resid1.0 <- lm(data = BRTETBmic_func_25, ln ~ BRTErel)
ols_plot_resid_hist(resid1.0) 
ols_test_normality(resid1.0) 

tb_mic25ab <- lmerTest::lmer(ln ~ BRTErel + (1|grad_num), data = BRTETBmic_func_25)
summary(tb_mic25ab) #p = 0.1821; slope = -0.004106
anova(tb_mic25ab)

#Other
tb_mic26a <- lmerTest::lmer(avg ~ BRTErel + (1|grad_num), data = subset(BRTETBmic_func_long2, func == "other"))
summary(tb_mic26a) #p = 0.0471; slope = -239.0
anova(tb_mic26a)

#check normality
resid1.0 <- lm(data = subset(BRTETBmic_func_long2, func == "other"), avg ~ BRTErel)
ols_plot_resid_hist(resid1.0) 
ols_test_normality(resid1.0) 


#Oxygen photo
tb_mic27a <- lmerTest::lmer(avg ~ BRTErel + (1|grad_num), data = subset(BRTETBmic_func_long2, func == "oxygenic_photoautotrophy"))
summary(tb_mic27a) #p = 0.00190; slope = -22.596
anova(tb_mic27a)

#check normality
resid1.0 <- lm(data = subset(BRTETBmic_func_long2, func == "oxygenic_photoautotrophy"), avg ~ BRTErel)
ols_plot_resid_hist(resid1.0) 
ols_test_normality(resid1.0) 
#log trans
BRTETBmic_func_27 <- BRTETBmic_func_long2 %>%
  filter(func == "oxygenic_photoautotrophy") %>%
  mutate(ln = log10(avg + 0.001)) 

resid1.0 <- lm(data = BRTETBmic_func_27, ln ~ BRTErel)
ols_plot_resid_hist(resid1.0) 
ols_test_normality(resid1.0) 

tb_mic27ab <- lmerTest::lmer(ln ~ BRTErel + (1|grad_num), data = BRTETBmic_func_27)
summary(tb_mic27ab) #p = 0.000411; slope = -0.019326
anova(tb_mic27ab)

#Photoauto
tb_mic28a <- lmerTest::lmer(avg ~ BRTErel + (1|grad_num), data = subset(BRTETBmic_func_long2, func == "photoautotrophy"))
summary(tb_mic28a) #p = 0.000285; slope = -27.287
anova(tb_mic28a)

#check normality
resid1.0 <- lm(data = subset(BRTETBmic_func_long2, func == "photoautotrophy"), avg ~ BRTErel)
ols_plot_resid_hist(resid1.0) 
ols_test_normality(resid1.0) 

#Photohetero
tb_mic29a <- lmerTest::lmer(avg ~ BRTErel + (1|grad_num), data = subset(BRTETBmic_func_long2, func == "photoheterotrophy"))
summary(tb_mic29a) #p = 0.0462; slope = -5.201
anova(tb_mic29a)

#check normality
resid1.0 <- lm(data = subset(BRTETBmic_func_long2, func == "photoheterotrophy"), avg ~ BRTErel)
ols_plot_resid_hist(resid1.0) 
ols_test_normality(resid1.0) 
#log trans
BRTETBmic_func_29 <- BRTETBmic_func_long2 %>%
  filter(func == "photoheterotrophy") %>%
  mutate(ln = log10(avg + 0.001)) 

resid1.0 <- lm(data = BRTETBmic_func_29, ln ~ BRTErel)
ols_plot_resid_hist(resid1.0) 
ols_test_normality(resid1.0) 

tb_mic29ab <- lmerTest::lmer(ln ~ BRTErel + (1|grad_num), data = BRTETBmic_func_29)
summary(tb_mic29ab) #p = 0.01113; slope = -0.003919
anova(tb_mic29ab)

#Photo cyano
tb_mic30a <- lmerTest::lmer(avg ~ BRTErel + (1|grad_num), data = subset(BRTETBmic_func_long2, func == "photosynthetic_cyanobacteria"))
summary(tb_mic30a) #p = 0.00190; slope = -22.596
anova(tb_mic30a)

#Phototrophy
tb_mic31a <- lmerTest::lmer(avg ~ BRTErel + (1|grad_num), data = subset(BRTETBmic_func_long2, func == "phototrophy"))
summary(tb_mic31a) #p = 0.00022; slope = -27.876
anova(tb_mic31a)

#check normality
resid1.0 <- lm(data = subset(BRTETBmic_func_long2, func == "phototrophy"), avg ~ BRTErel)
ols_plot_resid_hist(resid1.0) 
ols_test_normality(resid1.0) 

#Pred or exo
tb_mic32a <- lmerTest::lmer(avg ~ BRTErel + (1|grad_num), data = subset(BRTETBmic_func_long2, func == "predatory_or_exoparasitic"))
summary(tb_mic32a) #p = 0.0903; slope = -6.323
anova(tb_mic32a)

#check normality
resid1.0 <- lm(data = subset(BRTETBmic_func_long2, func == "predatory_or_exoparasitic"), avg ~ BRTErel)
ols_plot_resid_hist(resid1.0) 
ols_test_normality(resid1.0) 

#Resp sulfur
tb_mic33a <- lmerTest::lmer(avg ~ BRTErel + (1|grad_num), data = subset(BRTETBmic_func_long2, func == "respiration_of_sulfur_compounds"))
summary(tb_mic33a) #not sig
anova(tb_mic33a)

#check normality
resid1.0 <- lm(data = subset(BRTETBmic_func_long2, func == "respiration_of_sulfur_compounds"), avg ~ BRTErel)
ols_plot_resid_hist(resid1.0) 
ols_test_normality(resid1.0) 

#Sulfur resp
tb_mic34a <- lmerTest::lmer(avg ~ BRTErel + (1|grad_num), data = subset(BRTETBmic_func_long2, func == "sulfate_respiration"))
summary(tb_mic34a) #not sig
anova(tb_mic34a)

#Ureolysis
tb_mic35a <- lmerTest::lmer(avg ~ BRTErel + (1|grad_num), data = subset(BRTETBmic_func_long2, func == "ureolysis"))
summary(tb_mic35a) #p = 0.0171; slope = -6.604
anova(tb_mic35a)

#check normality
resid1.0 <- lm(data = subset(BRTETBmic_func_long2, func == "ureolysis"), avg ~ BRTErel)
ols_plot_resid_hist(resid1.0) 
ols_test_normality(resid1.0) 

#Cellulo
tb_mic36a <- lmerTest::lmer(avg ~ BRTErel + (1|grad_num), data = subset(BRTETBmic_func_long2, func == "cellulolysis"))
summary(tb_mic36a) #not sig
anova(tb_mic36a)

#check normality
resid1.0 <- lm(data = subset(BRTETBmic_func_long2, func == "cellulolysis"), avg ~ BRTErel)
ols_plot_resid_hist(resid1.0) 
ols_test_normality(resid1.0) 
#log trans
BRTETBmic_func_36 <- BRTETBmic_func_long2 %>%
  filter(func == "cellulolysis") %>%
  mutate(ln = log10(avg + 0.001)) 

resid1.0 <- lm(data = BRTETBmic_func_36, ln ~ BRTErel)
ols_plot_resid_hist(resid1.0) 
ols_test_normality(resid1.0) 

tb_mic36ab <- lmerTest::lmer(ln ~ BRTErel + (1|grad_num), data = BRTETBmic_func_36)
summary(tb_mic36ab) #p = 0.544; slope = 0.001506
anova(tb_mic36ab)

#Human associ
tb_mic37a <- lmerTest::lmer(avg ~ BRTErel + (1|grad_num), data = subset(BRTETBmic_func_long2, func == "human_associated"))
summary(tb_mic37a) #not sig
anova(tb_mic37a)

#check normality
resid1.0 <- lm(data = subset(BRTETBmic_func_long2, func == "human_associated"), avg ~ BRTErel)
ols_plot_resid_hist(resid1.0) 
ols_test_normality(resid1.0) 
#log trans
BRTETBmic_func_37 <- BRTETBmic_func_long2 %>%
  filter(func == "human_associated") %>%
  mutate(ln = log10(avg + 0.001)) 

resid1.0 <- lm(data = BRTETBmic_func_37, ln ~ BRTErel)
ols_plot_resid_hist(resid1.0) 
ols_test_normality(resid1.0) 

tb_mic37ab <- lmerTest::lmer(ln ~ BRTErel + (1|grad_num), data = BRTETBmic_func_37)
summary(tb_mic37ab) #p = 0.3577; slope = -0.003378
anova(tb_mic37ab)

#Human path
tb_mic38a <- lmerTest::lmer(avg ~ BRTErel + (1|grad_num), data = subset(BRTETBmic_func_long2, func == "human_pathogens_all"))
summary(tb_mic38a) #not sig
anova(tb_mic38a)

#check normality
resid1.0 <- lm(data = subset(BRTETBmic_func_long2, func == "human_pathogens_all"), avg ~ BRTErel)
ols_plot_resid_hist(resid1.0) 
ols_test_normality(resid1.0) 
#log trans
BRTETBmic_func_38 <- BRTETBmic_func_long2 %>%
  filter(func == "human_pathogens_all") %>%
  mutate(ln = log10(avg + 0.001)) 

resid1.0 <- lm(data = BRTETBmic_func_38, ln ~ BRTErel)
ols_plot_resid_hist(resid1.0) 
ols_test_normality(resid1.0) 

tb_mic38ab <- lmerTest::lmer(ln ~ BRTErel + (1|grad_num), data = BRTETBmic_func_38)
summary(tb_mic38ab) #p = 0.5395; slope = -0.002751
anova(tb_mic38ab)

#Iron resp
tb_mic39a <- lmerTest::lmer(avg ~ BRTErel + (1|grad_num), data = subset(BRTETBmic_func_long2, func == "iron_respiration"))
summary(tb_mic39a) #not sig
anova(tb_mic39a)

#check normality
resid1.0 <- lm(data = subset(BRTETBmic_func_long2, func == "iron_respiration"), avg ~ BRTErel)
ols_plot_resid_hist(resid1.0) 
ols_test_normality(resid1.0) 
#log trans
BRTETBmic_func_39 <- BRTETBmic_func_long2 %>%
  filter(func == "iron_respiration") %>%
  mutate(ln = log10(avg + 0.001)) 

resid1.0 <- lm(data = BRTETBmic_func_39, ln ~ BRTErel)
ols_plot_resid_hist(resid1.0) 
ols_test_normality(resid1.0) 

tb_mic39ab <- lmerTest::lmer(ln ~ BRTErel + (1|grad_num), data = BRTETBmic_func_39)
summary(tb_mic39ab) #p = 0.458; slope = -0.00329
anova(tb_mic39ab)

#Intracel para
tb_mic40a <- lmerTest::lmer(avg ~ BRTErel + (1|grad_num), data = subset(BRTETBmic_func_long2, func == "intracellular_parasites"))
summary(tb_mic40a) #not sig
anova(tb_mic40a)

#check normality
resid1.0 <- lm(data = subset(BRTETBmic_func_long2, func == "intracellular_parasites"), avg ~ BRTErel)
ols_plot_resid_hist(resid1.0) 
ols_test_normality(resid1.0) 
#log trans
BRTETBmic_func_40 <- BRTETBmic_func_long2 %>%
  filter(func == "intracellular_parasites") %>%
  mutate(ln = log10(avg + 0.001)) 

resid1.0 <- lm(data = BRTETBmic_func_40, ln ~ BRTErel)
ols_plot_resid_hist(resid1.0) 
ols_test_normality(resid1.0) 

tb_mic40ab <- lmerTest::lmer(ln ~ BRTErel + (1|grad_num), data = BRTETBmic_func_40)
summary(tb_mic40ab) #p = 0.1243; slope = 0.004074
anova(tb_mic40ab)

#Alipath non meth
tb_mic41a <- lmerTest::lmer(avg ~ BRTErel + (1|grad_num), data = subset(BRTETBmic_func_long2, func == "aliphatic_non_methane_hydrocarbon_degradation"))
summary(tb_mic41a) #not sig
anova(tb_mic41a)

#check normality
resid1.0 <- lm(data = subset(BRTETBmic_func_long2, func == "aliphatic_non_methane_hydrocarbon_degradation"), avg ~ BRTErel)
ols_plot_resid_hist(resid1.0) 
ols_test_normality(resid1.0) 
#log trans
BRTETBmic_func_41 <- BRTETBmic_func_long2 %>%
  filter(func == "aliphatic_non_methane_hydrocarbon_degradation") %>%
  mutate(ln = log10(avg + 0.001)) 

resid1.0 <- lm(data = BRTETBmic_func_41, ln ~ BRTErel)
ols_plot_resid_hist(resid1.0) 
ols_test_normality(resid1.0) 

tb_mic41ab <- lmerTest::lmer(ln ~ BRTErel + (1|grad_num), data = BRTETBmic_func_41)
summary(tb_mic41ab) #p = 0.4405; slope = -0.006543
anova(tb_mic41ab)

#Aromat hydro
tb_mic42a <- lmerTest::lmer(avg ~ BRTErel + (1|grad_num), data = subset(BRTETBmic_func_long2, func == "aromatic_hydrocarbon_degradation"))
summary(tb_mic42a) #not sig
anova(tb_mic42a)

#check normality
resid1.0 <- lm(data = subset(BRTETBmic_func_long2, func == "aromatic_hydrocarbon_degradation"), avg ~ BRTErel)
ols_plot_resid_hist(resid1.0) 
ols_test_normality(resid1.0) 
#log trans
BRTETBmic_func_42 <- BRTETBmic_func_long2 %>%
  filter(func == "aromatic_hydrocarbon_degradation") %>%
  mutate(ln = log10(avg + 0.001)) 

resid1.0 <- lm(data = BRTETBmic_func_42, ln ~ BRTErel)
ols_plot_resid_hist(resid1.0) 
ols_test_normality(resid1.0) 

tb_mic42ab <- lmerTest::lmer(ln ~ BRTErel + (1|grad_num), data = BRTETBmic_func_42)
summary(tb_mic42ab) #p = 0.6412; slope = -0.003151
anova(tb_mic42ab)

#Hydro deg
tb_mic43a <- lmerTest::lmer(avg ~ BRTErel + (1|grad_num), data = subset(BRTETBmic_func_long2, func == "hydrocarbon_degradation"))
summary(tb_mic43a) #not sig
anova(tb_mic43a)

#Human gut
tb_mic44a <- lmerTest::lmer(avg ~ BRTErel + (1|grad_num), data = subset(BRTETBmic_func_long2, func == "human_gut"))
summary(tb_mic44a) #not sig
anova(tb_mic44a)

#check normality
resid1.0 <- lm(data = subset(BRTETBmic_func_long2, func == "human_gut"), avg ~ BRTErel)
ols_plot_resid_hist(resid1.0) 
ols_test_normality(resid1.0) 
#log trans
BRTETBmic_func_44 <- BRTETBmic_func_long2 %>%
  filter(func == "human_gut") %>%
  mutate(ln = log10(avg + 0.001)) 

resid1.0 <- lm(data = BRTETBmic_func_44, ln ~ BRTErel)
ols_plot_resid_hist(resid1.0) 
ols_test_normality(resid1.0) 

tb_mic44ab <- lmerTest::lmer(ln ~ BRTErel + (1|grad_num), data = BRTETBmic_func_44)
summary(tb_mic44ab) #p = 0.9496; slope = -0.0003283
anova(tb_mic44ab)

#Mammal gut
tb_mic45a <- lmerTest::lmer(avg ~ BRTErel + (1|grad_num), data = subset(BRTETBmic_func_long2, func == "mammal_gut"))
summary(tb_mic45a) #not sig
anova(tb_mic45a)

#Plant pathogen
tb_mic46a <- lmerTest::lmer(avg ~ BRTErel + (1|grad_num), data = subset(BRTETBmic_func_long2, func == "plant_pathogen"))
summary(tb_mic46a) #not sig
anova(tb_mic46a)

#check normality
resid1.0 <- lm(data = subset(BRTETBmic_func_long2, func == "plant_pathogen"), avg ~ BRTErel)
ols_plot_resid_hist(resid1.0) 
ols_test_normality(resid1.0) 


#Xyanolysis
tb_mic47a <- lmerTest::lmer(avg ~ BRTErel + (1|grad_num), data = subset(BRTETBmic_func_long2, func == "xylanolysis"))
summary(tb_mic47a) #not enough data 
anova(tb_mic47a)

#Dark thiosul
tb_mic48a <- lmerTest::lmer(avg ~ BRTErel + (1|grad_num), data = subset(BRTETBmic_func_long2, func == "dark_thiosulfate_oxidation"))
summary(tb_mic48a) #not sig 
anova(tb_mic48a)

#check normality
resid1.0 <- lm(data = subset(BRTETBmic_func_long2, func == "dark_thiosulfate_oxidation"), avg ~ BRTErel)
ols_plot_resid_hist(resid1.0) 
ols_test_normality(resid1.0) 

#Human path
tb_mic49a <- lmerTest::lmer(avg ~ BRTErel + (1|grad_num), data = subset(BRTETBmic_func_long2, func == "human_pathogens_pneumonia"))
summary(tb_mic49a) #not sig 
anova(tb_mic49a)
#check normality
resid1.0 <- lm(data = subset(BRTETBmic_func_long2, func == "human_pathogens_pneumonia"), avg ~ BRTErel)
ols_plot_resid_hist(resid1.0) 
ols_test_normality(resid1.0) 




tb_mic50a <- lmerTest::lmer(avg ~ BRTErel + (1|grad_num), data = subset(BRTETBmic_func_long2, func == "human_pathogens_septicemia"))
summary(tb_mic50a) #not sig 
anova(tb_mic50a)

tb_mic51a <- lmerTest::lmer(avg ~ BRTErel + (1|grad_num), data = subset(BRTETBmic_func_long2, func == "plastic_degradation"))
summary(tb_mic51a) #not sig
anova(tb_mic51a)






##############################################################################################################




#############################################################################################################
############# Multiple regression figures ##########


#### put all these graphs together

#FK
#plant
FK_BRAR_func_reg <- ggplot(data = FKcov_func, aes(x = BRARrel, y = avg_cov, color = funct2)) +
  geom_smooth(data = subset(FKcov_func, funct2 == "C3 Annual Grass"), aes(group = funct2), method = "lm", se = FALSE, linetype = "dotted") +
  geom_smooth(data = subset(FKcov_func, funct2 == "C3 Perennial Grass"), aes(group = funct2), method = "lm", se = TRUE) +
  geom_smooth(data = subset(FKcov_func, funct2 == "C4 Perennial Grass"), aes(group = funct2), method = "lm", se = TRUE) +
  geom_smooth(data = subset(FKcov_func, funct2 == "Cactus"), aes(group = funct2), method = "lm", se = FALSE, linetype = "dotted") +
  geom_smooth(data = subset(FKcov_func, funct2 == "Forb"), aes(group = funct2), method = "lm", se = TRUE, linetype = "dashed") +
  geom_smooth(data = subset(FKcov_func, funct2 == "Sub-Shrub/Shrub"), aes(group = funct2), method = "lm", se = TRUE) +
  scale_colour_manual(values = c("black", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "red3")) +
  ylim(0, 45) + 
  labs(x = "", y = "Cover (%)", title = title1, color = "Functional Group", shape = "Functional Group") + 
  theme(plot.title = element_text(vjust = - 7, hjust = 0.07, size = 16)) 

#insect
FK_BRAR_func_reg_ins <- ggplot(data = FKinsfunc, aes(x = BRARrel, y = avg, color = guild)) +
  geom_smooth(data = subset(FKinsfunc, guild == "Leaf Chewing Herbivore"), aes(group = guild), method = "lm", se = TRUE) +
  geom_smooth(data = subset(FKinsfunc, guild == "Predator"), aes(group = guild), method = "lm", se = FALSE, linetype = "dotted") +
  geom_smooth(data = subset(FKinsfunc, guild == "Sap Sucking Herbivore"), aes(group = guild), method = "lm", se = TRUE) +
  geom_smooth(data = subset(FKinsfunc, guild == "Parasitoid"), aes(group = guild), method = "lm", se = FALSE, linetype = "dotted") +
  scale_colour_manual(values = c("black", "#E69F00", "#56B4E9", "#009E73")) +
  ylim(0, 15) +
  labs(x = "", y = "Abundance (count)", title = title1, color = "Functional Group", shape = "Functional Group") + 
  theme(plot.title = element_text(vjust = - 7, hjust = 0.07, size = 16)) 

#microbe
FK_BRAR_func_reg_mic <- ggplot(data = BRARFKmic_func_long3, aes(x = BRARrel, y = avg, color = func)) +
  geom_smooth(data = subset(BRARFKmic_func_long3, func == "oxygenic_photoautotrophy"), aes(group = func), method = "lm", se = TRUE, linetype = "dashed", position = position_dodge(0.5)) +
  geom_smooth(data = subset(BRARFKmic_func_long3, func == "photosynthetic_cyanobacteria"), aes(group = func), method = "lm", se = TRUE, linetype = "dashed", position = position_dodge(0.5)) +
  scale_colour_manual(values = c("oxygenic_photoautotrophy" = "black", "photosynthetic_cyanobacteria" = "#009E73")) +
  ylim(0, 65) +
  labs(x = "Cover (%)", y = "Abundance (count)", title = title1, color = "Functional Group", shape = "Functional Group") + 
  theme(plot.title = element_text(vjust = - 7, hjust = 0.07, size = 16)) 


#TB
#BRAR
#plant
TB_BRAR_func_reg <- ggplot(data = TBBRARcov_func, aes(x = BRARrel, y = avg_cov, color = funct2)) +
  geom_smooth(data = subset(TBBRARcov_func, funct2 == "C3 Annual Grass"), aes(group = funct2), method = "lm", se = TRUE, linetype = "dashed") +
  geom_smooth(data = subset(TBBRARcov_func, funct2 == "C3 Perennial Grass"), aes(group = funct2), method = "lm", se = FALSE, linetype = "dotted") +
  geom_smooth(data = subset(TBBRARcov_func, funct2 == "C4 Perennial Grass"), aes(group = funct2), method = "lm", se = TRUE) +
  geom_smooth(data = subset(TBBRARcov_func, funct2 == "Cactus"), aes(group = funct2), method = "lm", se = FALSE, linetype = "dotted") +
  geom_smooth(data = subset(TBBRARcov_func, funct2 == "Forb"), aes(group = funct2), method = "lm", se = FALSE, linetype = "dotted") +
  geom_smooth(data = subset(TBBRARcov_func, funct2 == "Sub-Shrub/Shrub"), aes(group = funct2), method = "lm", se = FALSE, linetype = "dotted") +
  scale_colour_manual(values = c("black", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "red3")) +
  ylim(0, 45) +
  labs(x = "", y = "", title = title2, color = "Functional Group", shape = "Functional Group") + 
  theme(plot.title = element_text(vjust = - 7, hjust = 0.07, size = 16)) 

#insect
TB_BRAR_func_reg_ins <- ggplot(data = TBBRARinsfunc, aes(x = BRARrel, y = avg, color = guild)) +
  geom_smooth(data = subset(TBBRARinsfunc, guild == "Leaf Chewing Herbivore"), aes(group = guild), method = "lm", se = FALSE, linetype = "dotted") +
  geom_smooth(data = subset(TBBRARinsfunc, guild == "Predator"), aes(group = guild), method = "lm", se = FALSE, linetype = "dotted") +
  geom_smooth(data = subset(TBBRARinsfunc, guild == "Sap Sucking Herbivore"), aes(group = guild), method = "lm", se = TRUE) +
  geom_smooth(data = subset(TBBRARinsfunc, guild == "Parasitoid"), aes(group = guild), method = "lm", se = TRUE) +
  scale_colour_manual(values = c("black", "#E69F00", "#56B4E9", "#009E73")) +
  ylim(0, 20) +
  labs(x = "", y = "", title = title2, color = "Functional Group", shape = "Functional Group") + 
  theme(plot.title = element_text(vjust = - 7, hjust = 0.07, size = 16)) 

#microbe
TB_BRAR_func_reg_mic <- ggplot(data = BRARTBmic_func_long3, aes(x = BRARrel, y = avg, color = func)) +
  geom_smooth(data = subset(BRARTBmic_func_long3, func == "aerobic_ammonia_oxidation"), aes(group = func), method = "lm", se = TRUE, linetype = "dashed") +
  geom_smooth(data = subset(BRARTBmic_func_long3, func == "chitinolysis"), aes(group = func), method = "lm", se = TRUE, linetype = "dashed") +
  geom_smooth(data = subset(BRARTBmic_func_long3, func == "dark_oxidation_of_sulfur_compounds"), aes(group = func), method = "lm", se = TRUE, linetype = "dashed") +
  geom_smooth(data = subset(BRARTBmic_func_long3, func == "nitrate_reduction"), aes(group = func), method = "lm", se = TRUE, linetype = "dashed") +
  geom_smooth(data = subset(BRARTBmic_func_long3, func == "oxygenic_photoautotrophy"), aes(group = func), method = "lm", se = TRUE) +
  geom_smooth(data = subset(BRARTBmic_func_long3, func == "photoautotrophy"), aes(group = func), method = "lm", se = TRUE) +
  geom_smooth(data = subset(BRARTBmic_func_long3, func == "photosynthetic_cyanobacteria"), aes(group = func), method = "lm", se = TRUE) +
  geom_smooth(data = subset(BRARTBmic_func_long3, func == "phototrophy"), aes(group = func), method = "lm", se = TRUE) +
  geom_smooth(data = subset(BRARTBmic_func_long3, func == "iron_respiration"), aes(group = func), method = "lm", se = TRUE, linetype = "dashed") +
  scale_colour_manual(values = c("oxygenic_photoautotrophy" = "black", "photoautotrophy" = "#56B4E9", "photosynthetic_cyanobacteria" = "#009E73", 
                                 "phototrophy" = "#F0E442", "aerobic_ammonia_oxidation" = "red3", "chitinolysis" = "purple", "dark_oxidation_of_sulfur_compounds" = "green", "nitrate_reduction" = "navy", 
                                 "human_associated" = "lightpink", "human_pathogens_all" = "burlywood", "iron_respiration" = "royalblue", "aliphatic_non_methane_hydrocarbon_degradation" = "cyan2", 
                                 "aromatic_hydrocarbon_degradation" = "gold4", "hydrocarbon_degradation" = "springgreen2")) +
  ylim(0, 2500) +
  labs(x = "Cover (%)", y = "", title = title2, color = "Functional Group", shape = "Functional Group") + 
  theme(plot.title = element_text(vjust = - 7, hjust = 0.07, size = 16)) 


#BRTE
#plant
TB_BRTE_func_reg <- ggplot(data = TBBRTEcov_func, aes(x = BRTErel, y = avg_cov, color = funct2)) +
  geom_smooth(data = subset(TBBRTEcov_func, funct2 == "C3 Annual Grass"), aes(group = funct2), method = "lm", se = TRUE) +
  geom_smooth(data = subset(TBBRTEcov_func, funct2 == "C3 Perennial Grass"), aes(group = funct2), method = "lm", se = TRUE) +
  geom_smooth(data = subset(TBBRTEcov_func, funct2 == "C4 Perennial Grass"), aes(group = funct2), method = "lm", se = TRUE) +
  geom_smooth(data = subset(TBBRTEcov_func, funct2 == "Cactus"), aes(group = funct2), method = "lm", se = FALSE, linetype = "dotted") +
  geom_smooth(data = subset(TBBRTEcov_func, funct2 == "Forb"), aes(group = funct2), method = "lm", se = FALSE, linetype = "dotted") +
  geom_smooth(data = subset(TBBRTEcov_func, funct2 == "Sub-Shrub/Shrub"), aes(group = funct2), method = "lm", se = FALSE, linetype = "dotted") +
  scale_colour_manual(values = c("black", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "red3")) +
  ylim(0, 35) +
  labs(x = "", y = "", title = title3, color = "Functional Group", shape = "Functional Group") + 
  theme(plot.title = element_text(vjust = - 7, hjust = 0.07, size = 16)) 

#insect
TB_BRTE_func_reg_ins <- ggplot(data = TBBRTEinsfunc, aes(x = BRTErel, y = avg, color = guild)) +
  geom_smooth(data = subset(TBBRTEinsfunc, guild == "Leaf Chewing Herbivore"), aes(group = guild), method = "lm", se = FALSE, linetype = "dotted") +
  geom_smooth(data = subset(TBBRTEinsfunc, guild == "Predator"), aes(group = guild), method = "lm", se = FALSE, linetype = "dotted") +
  geom_smooth(data = subset(TBBRTEinsfunc, guild == "Sap Sucking Herbivore"), aes(group = guild), method = "lm", se = FALSE, linetype = "dotted") +
  geom_smooth(data = subset(TBBRTEinsfunc, guild == "Parasitoid"), aes(group = guild), method = "lm", se = FALSE, linetype = "dotted") +
  scale_colour_manual(values = c("black", "#E69F00", "#56B4E9", "#009E73")) +
  ylim(0, 10) +
  labs(x = "", y = "", title = title3, color = "Functional Group", shape = "Functional Group") + 
  theme(plot.title = element_text(vjust = - 7, hjust = 0.07, size = 16)) 

#microbe
TB_BRTE_func_reg_mic <- ggplot(data = BRTETBmic_func_long3, aes(x = BRTErel, y = avg, color = func)) +
  geom_smooth(data = subset(BRTETBmic_func_long3, func == "aerobic_chemoheterotrophy"), aes(group = func), method = "lm", se = TRUE, linetype = "dashed") +
  geom_smooth(data = subset(BRTETBmic_func_long3, func == "anoxygenic_photoautotrophy"), aes(group = func), method = "lm", se = TRUE) +
  geom_smooth(data = subset(BRTETBmic_func_long3, func == "anoxygenic_photoautotrophy_S_oxidizing"), aes(group = func), method = "lm", se = TRUE) +
  geom_smooth(data = subset(BRTETBmic_func_long3, func == "aromatic_compound_degradation"), aes(group = func), method = "lm", se = TRUE) +
  geom_smooth(data = subset(BRTETBmic_func_long3, func == "chemoheterotrophy"), aes(group = func), method = "lm", se = TRUE, linetype = "dashed") +
  geom_smooth(data = subset(BRTETBmic_func_long3, func == "denitrification"), aes(group = func), method = "lm", se = TRUE) +
  geom_smooth(data = subset(BRTETBmic_func_long3, func == "methanol_oxidation"), aes(group = func), method = "lm", se = TRUE) +
  geom_smooth(data = subset(BRTETBmic_func_long3, func == "methylotrophy"), aes(group = func), method = "lm", se = TRUE) +
  geom_smooth(data = subset(BRTETBmic_func_long3, func == "nitrate_denitrification"), aes(group = func), method = "lm", se = TRUE) +
  geom_smooth(data = subset(BRTETBmic_func_long3, func == "nitrate_respiration"), aes(group = func), method = "lm", se = TRUE) +
  geom_smooth(data = subset(BRTETBmic_func_long3, func == "nitrite_denitrification"), aes(group = func), method = "lm", se = TRUE) +
  geom_smooth(data = subset(BRTETBmic_func_long3, func == "nitrite_respiration"), aes(group = func), method = "lm", se = TRUE) +
  geom_smooth(data = subset(BRTETBmic_func_long3, func == "nitrogen_respiration"), aes(group = func), method = "lm", se = TRUE) +
  geom_smooth(data = subset(BRTETBmic_func_long3, func == "nitrous_oxide_denitrification"), aes(group = func), method = "lm", se = TRUE) +
  geom_smooth(data = subset(BRTETBmic_func_long3, func == "oxygenic_photoautotrophy"), aes(group = func), method = "lm", se = TRUE) +
  geom_smooth(data = subset(BRTETBmic_func_long3, func == "photoautotrophy"), aes(group = func), method = "lm", se = TRUE) +
  geom_smooth(data = subset(BRTETBmic_func_long3, func == "photoheterotrophy"), aes(group = func), method = "lm", se = TRUE) +
  geom_smooth(data = subset(BRTETBmic_func_long3, func == "photosynthetic_cyanobacteria"), aes(group = func), method = "lm", se = TRUE) +
  geom_smooth(data = subset(BRTETBmic_func_long3, func == "phototrophy"), aes(group = func), method = "lm", se = TRUE) +
  geom_smooth(data = subset(BRTETBmic_func_long3, func == "predatory_or_exoparasitic"), aes(group = func), method = "lm", se = TRUE, linetype = "dashed") +
  geom_smooth(data = subset(BRTETBmic_func_long3, func == "ureolysis"), aes(group = func), method = "lm", se = TRUE) +
  scale_colour_manual(values = c("aerobic_chemoheterotrophy" = "darkgreen","oxygenic_photoautotrophy" = "black", "photoautotrophy" = "#56B4E9", "photosynthetic_cyanobacteria" = "#009E73", 
                                 "phototrophy" = "#F0E442", "anoxygenic_photoautotrophy" = "red3", "anoxygenic_photoautotrophy_S_oxidizing" = "purple", "aromatic_compound_degradation" = "green", 
                                "chemoheterotrophy" = "purple4", "denitrification" = "navy", "methanol_oxidation" = "lightpink", "methylotrophy" = "burlywood", "nitrate_denitrification" = "royalblue", "nitrate_respiration" = "cyan2", 
                                 "nitrite_denitrification" = "gold4", "nitrite_respiration" = "springgreen2", "nitrogen_respiration" = "mediumorchid3", "nitrous_oxide_denitrification" = "lightsalmon1", 
                                 "photoheterotrophy" = "grey36", "predatory_or_exoparasitic" = "tomato", "ureolysis" = "#E69F00")) +
  ylim(0, 2500) +
  labs(x = "Cover (%)", y = "", title = title3, color = "Functional Group", shape = "Functional Group") + 
  theme(plot.title = element_text(vjust = - 7, hjust = 0.07, size = 16)) 


FK_BRAR_func_reg + TB_BRAR_func_reg + TB_BRTE_func_reg + 
FK_BRAR_func_reg_ins + TB_BRAR_func_reg_ins + TB_BRTE_func_reg_ins +
FK_BRAR_func_reg_mic + TB_BRAR_func_reg_mic + TB_BRTE_func_reg_mic + plot_layout(ncol = 3)


#reg 3300 x 1700



###########################################################################







##########################################################################
############# Insect biomass ##########
#FK 2020
FK2020bio_join <- FK2020bio %>%
  full_join(FK_plotinfo) %>%
  full_join(FK_commBRAR) %>%
  drop_na(total_biomass)

#FK2021
FK2021bio_join <- FKcover %>%
  filter(year == "2021") %>%
  filter(symbol == "BRAR") %>%
  mutate(BRARrel = rel_cov) %>%
  full_join(FK2021bio) %>%
  drop_na(total_biomass) %>%
  drop_na(BRARrel)

#FK2022
FK2022bio_join <- FKcover %>%
  filter(year == "2022") %>%
  filter(symbol == "BRAR") %>%
  mutate(BRARrel = rel_cov) %>%
  full_join(FK2022bio) %>%
  drop_na(total_biomass) %>%
  drop_na(BRARrel)

#TB
TBbio_join <- TBbio %>%
  full_join(TB_plotinfo) 

TBbio_join_BRAR <- TBbio_join %>%
  filter(invasive_type == "BRAR") %>%
  full_join(TB_commBRAR) %>%
  filter(plot != "8")

TBbio_join_BRTE <- TBbio_join %>%
  filter(invasive_type == "BRTE") %>%
  full_join(TB_commBRTE) 




#Figures
#FK
FK2020bio_fig <- ggplot(data = FK2020bio_join, aes(x = BRARrel, y = total_biomass)) + 
  geom_point(size = 3) + 
  #xlim(0, 80) +
  ylim(0, 1.5) +
  labs(x = "BRAR Cov (%)", y = "Insect biomass (g)", title = "MT BRAR 2020") +
  theme(plot.title = element_text(vjust = - 7, hjust = 0.07, size = 16))

FK2021bio_fig <- ggplot(data = FK2021bio_join, aes(x = BRARrel, y = total_biomass)) + 
  geom_point(size = 3) + 
  #xlim(0, 80) +
  ylim(0, 0.5) +
  labs(x = "BRAR Cov (%)", y = "Insect biomass (g)", title = "MT BRAR 2021") +
  theme(plot.title = element_text(vjust = - 7, hjust = 0.07, size = 16))

FK2022bio_fig <- ggplot(data = FK2022bio_join, aes(x = BRARrel, y = total_biomass)) + 
  geom_point(size = 3) + 
  #xlim(0, 80) +
  ylim(0, 0.8) +
  labs(x = "BRAR Cov (%)", y = "Insect biomass (g)", title = "MT BRAR 2022") +
  theme(plot.title = element_text(vjust = - 7, hjust = 0.07, size = 16))

#TB
TBBRARbio_fig <- ggplot(data = TBbio_join_BRAR, aes(x = BRARrel, y = total_biomass)) + 
  geom_point(size = 3) + 
  #xlim(0, 80) +
  ylim(0, 0.2) +
  labs(x = "BRAR Cov (%)", y = "Insect biomass (g)", title = "WY BRAR 2019") +
  theme(plot.title = element_text(vjust = - 7, hjust = 0.07, size = 16))

TBBRTEbio_fig <- ggplot(data = TBbio_join_BRTE, aes(x = BRTErel, y = total_biomass)) + 
  geom_point(size = 3) + 
  #xlim(0, 80) +
  ylim(0, 0.1) +
  labs(x = "BRTE Cov (%)", y = "Insect biomass (g)", title = "WY BRTE 2019") +
  theme(plot.title = element_text(vjust = - 7, hjust = 0.07, size = 16))

FK2020bio_fig + FK2021bio_fig + FK2022bio_fig + TBBRARbio_fig + TBBRTEbio_fig + plot_layout(nrow = 1)

#bio 2300 x 400

#stats
#FK 2020
FKbio2020stats <- lmerTest::lmer(data = FK2020bio_join, total_biomass ~ BRARrel +
                                   (1|grad_num)) 

anova(FKbio2020stats, type = 3) #p = 0.313

#normality check
res_a <- lm(data = FK2020bio_join, total_biomass ~ BRARrel)
ols_plot_resid_hist(res_a) #right skew
ols_test_normality(res_a) #pass KS

#try trans
FK2020bio_join_trans <- FK2020bio_join %>%
  mutate(log = log10(total_biomass))

res_a <- lm(data = FK2020bio_join_trans, log ~ BRARrel)
ols_plot_resid_hist(res_a) #right skew
ols_test_normality(res_a) #pass KS and AD

#redo stats with log
FKbio2020stats2 <- lmerTest::lmer(data = FK2020bio_join_trans, log ~ BRARrel +
                                    (1|grad_num)) 

anova(FKbio2020stats2, type = 3) #p = 0.1459

#linearity
plot(resid(FKbio2020stats2), FK2020bio_join_trans$log) #looks linear
plot(FKbio2020stats2) #no pattern so indicates linearity

#homoscedascity
FK2020bio_join_trans$res <- residuals(FKbio2020stats2)
FK2020bio_join_trans$abs_res <- abs(FK2020bio_join_trans$res)
FK2020bio_join_trans$abs_res2 <- FK2020bio_join_trans$abs_res^2
levene_brar_bio20 <- lm(data = FK2020bio_join_trans, abs_res2 ~ BRARrel)
anova(levene_brar_bio20) #p = 0.2878 so > 0.05 so equal variance is met

#autocorrelation
acf(residuals(FKbio2020stats2, retype = "normalized"))  #lines dont go outside CI horizontal lines, so not autocorrelated
pacf(residuals(FKbio2020stats2, retype = "normalized"))


#FK 2021
FKbio2021stats <- lmerTest::lmer(data = FK2021bio_join, total_biomass ~ BRARrel +
                                   (1|grad_num)) 

anova(FKbio2021stats, type = 3) #p = 0.8061

#normality check
res_a <- lm(data = FK2021bio_join, total_biomass ~ BRARrel)
ols_plot_resid_hist(res_a) #close to normal
ols_test_normality(res_a) #pass KS and AD

#try trans
FK2021bio_join_trans <- FK2021bio_join %>%
  mutate(log = log10(total_biomass))

res_a <- lm(data = FK2021bio_join_trans, log ~ BRARrel)
ols_plot_resid_hist(res_a) #normal
ols_test_normality(res_a) #pass

#redo stats with log
FKbio2021stats2 <- lmerTest::lmer(data = FK2021bio_join_trans, log ~ BRARrel +
                                    (1|grad_num)) 

anova(FKbio2021stats2, type = 3) #p = 0.3043

#linearity
plot(resid(FKbio2021stats2), FK2021bio_join_trans$log) #looks linear
plot(FKbio2021stats2) #no pattern so indicates linearity

#homoscedascity
FK2021bio_join_trans$res <- residuals(FKbio2021stats2)
FK2021bio_join_trans$abs_res <- abs(FK2021bio_join_trans$res)
FK2021bio_join_trans$abs_res2 <- FK2021bio_join_trans$abs_res^2
levene_brar_bio21 <- lm(data = FK2021bio_join_trans, abs_res2 ~ BRARrel)
anova(levene_brar_bio21) #p = 0.1215 so > 0.05 so equal variance is met

#autocorrelation
acf(residuals(FKbio2021stats2, retype = "normalized"))  #lines dont go outside CI horizontal lines, so not autocorrelated
pacf(residuals(FKbio2021stats2, retype = "normalized"))

#FK 2022
FKbio2022stats <- lmerTest::lmer(data = FK2022bio_join, total_biomass ~ BRARrel +
                                   (1|grad_num)) 

anova(FKbio2022stats, type = 3) #p = 0.3225

#normality check
res_a <- lm(data = FK2022bio_join, total_biomass ~ BRARrel)
ols_plot_resid_hist(res_a) #skew
ols_test_normality(res_a) #fail

#try trans
FK2022bio_join_trans <- FK2022bio_join %>%
  mutate(log = log10(total_biomass))

res_a <- lm(data = FK2022bio_join_trans, log ~ BRARrel)
ols_plot_resid_hist(res_a) #normaler
ols_test_normality(res_a) #pass 

#redo stats with log
FKbio2022stats2 <- lmerTest::lmer(data = FK2022bio_join_trans, log ~ BRARrel +
                                    (1|grad_num)) 

anova(FKbio2022stats2, type = 3) #p = 0.8967

#linearity
plot(resid(FKbio2022stats2), FK2022bio_join_trans$log) #looks linear
plot(FKbio2022stats2) #no pattern so indicates linearity

#homoscedascity
FK2022bio_join_trans$res <- residuals(FKbio2022stats2)
FK2022bio_join_trans$abs_res <- abs(FK2022bio_join_trans$res)
FK2022bio_join_trans$abs_res2 <- FK2022bio_join_trans$abs_res^2
levene_brar_bio22 <- lm(data = FK2022bio_join_trans, abs_res2 ~ BRARrel)
anova(levene_brar_bio22) #p = 0.07 so > 0.05 so equal variance is met

#autocorrelation
acf(residuals(FKbio2022stats2, retype = "normalized"))  #lines dont go outside CI horizontal lines, so not autocorrelated
pacf(residuals(FKbio2022stats2, retype = "normalized"))


#TB BRAR
TBbio_join_BRARstats <- lmerTest::lmer(data = TBbio_join_BRAR, total_biomass ~ BRARrel +
                                         (1|grad_num)) 

anova(TBbio_join_BRARstats, type = 3) #p = 0.4155

#normality check
res_a <- lm(data = TBbio_join_BRAR, total_biomass ~ BRARrel)
ols_plot_resid_hist(res_a) # normal
ols_test_normality(res_a) #pass 

#linearity
plot(resid(TBbio_join_BRARstats), TBbio_join_BRAR$total_biomass) #looks linear
plot(TBbio_join_BRARstats) #no pattern so indicates linearity

#homoscedascity
TBbio_join_BRAR$res <- residuals(TBbio_join_BRARstats)
TBbio_join_BRAR$abs_res <- abs(TBbio_join_BRAR$res)
TBbio_join_BRAR$abs_res2 <- TBbio_join_BRAR$abs_res^2
levene_brar_biotb <- lm(data = TBbio_join_BRAR, abs_res2 ~ BRARrel)
anova(levene_brar_biotb) #p = 0.2654 so > 0.05 so equal variance is met

#autocorrelation
acf(residuals(TBbio_join_BRARstats, retype = "normalized"))  #lines dont go outside CI horizontal lines, so not autocorrelated
pacf(residuals(TBbio_join_BRARstats, retype = "normalized"))


#TB BRTE
TBbio_join_BRTEstats <- lmerTest::lmer(data = TBbio_join_BRTE, total_biomass ~ BRTErel +
                                         (1|grad_num)) 

anova(TBbio_join_BRTEstats, type = 3) #p = 0.3385

#normality check
res_a <- lm(data = TBbio_join_BRTE, total_biomass ~ BRTErel)
ols_plot_resid_hist(res_a) # normal
ols_test_normality(res_a) #pass 

#try trans
TBbio_join_BRTE_trans <- TBbio_join_BRTE %>%
  mutate(log = log10(total_biomass))

res_a <- lm(data = TBbio_join_BRTE_trans, log ~ BRTErel)
ols_plot_resid_hist(res_a) #normaler
ols_test_normality(res_a) #pass 

#redo stats with log
TBbio_join_BRTEstats2 <- lmerTest::lmer(data = TBbio_join_BRTE_trans, log ~ BRTErel +
                                          (1|grad_num)) 

anova(TBbio_join_BRTEstats2, type = 3) #p = 0.8939

#linearity
plot(resid(TBbio_join_BRTEstats2), TBbio_join_BRTE_trans$log) #looks linear
plot(TBbio_join_BRTEstats2) #no pattern so indicates linearity

#homoscedascity
TBbio_join_BRTE_trans$res <- residuals(TBbio_join_BRTEstats2)
TBbio_join_BRTE_trans$abs_res <- abs(TBbio_join_BRTE_trans$res)
TBbio_join_BRTE_trans$abs_res2 <- TBbio_join_BRTE_trans$abs_res^2
levene_brte_biotb <- lm(data = TBbio_join_BRTE_trans, abs_res2 ~ BRTErel)
anova(levene_brte_biotb) #p = 0.2654 so > 0.05 so equal variance is met

#autocorrelation
acf(residuals(TBbio_join_BRTEstats2, retype = "normalized"))  #lines dont go outside CI horizontal lines, so not autocorrelated
pacf(residuals(TBbio_join_BRTEstats2, retype = "normalized"))






#######################################################################




####################################################################
################### Insect herbivory ##########
#FK herbivory 2021 and 2022
#2021
FK2021herb_join <- FKcover %>%
  filter(year == "2021") %>%
  filter(symbol == "BRAR") %>%
  mutate(BRARrel = rel_cov) %>%
  full_join(FK2021herb) %>%
  filter(native_invasive != "introduced")

#2022
FK2022herb_join <- FKcover %>%
  filter(year == "2022") %>%
  filter(symbol == "BRAR") %>%
  mutate(BRARrel = rel_cov) %>%
  full_join(FK2022herb) %>%
  filter(native_invasive != "introduced")
FK2022herb_join[is.na(FK2022herb_join)] <- 0


FK2021_sp_box <- ggplot(data = FK2021herb_join, aes(x = plant_sp, y = per_herb_plant, color = native_invasive)) +
  geom_boxplot() +
  labs(title = "MT 2021")

FK2022_sp_box <- ggplot(data = FK2022herb_join, aes(x = plant_sp, y = per_herb_plant, color = native_invasive)) +
  geom_boxplot() +
  labs(title = "MT 2022")

FK2021_sp_box / FK2022_sp_box

#sp1  2500 x 500   


FK2021herb_join2 <- FK2021herb_join %>%
  rowwise() %>%
  mutate(avg_leaf_herb = mean(c(per_herb_leaf1, per_herb_leaf2, per_herb_leaf3, per_herb_leaf4, per_herb_leaf5, per_herb_leaf6, per_herb_leaf7, per_herb_leaf8, per_herb_leaf9, per_herb_leaf10), na.rm = TRUE)) %>%
  group_by(native_invasive, plot, BRARrel, invasion_percent, grad_num) %>%
  summarise(avg_leaf = mean(avg_leaf_herb), avg_plant = mean(per_herb_plant)) %>%
  ungroup()

FK2022herb_join2 <- FK2022herb_join %>%
  rowwise() %>%
  mutate(avg_leaf_herb = mean(c(per_herb_leaf1, per_herb_leaf2, per_herb_leaf3, per_herb_leaf4, per_herb_leaf5, per_herb_leaf6, per_herb_leaf7, per_herb_leaf8, per_herb_leaf9, per_herb_leaf10), na.rm = TRUE)) %>%
  group_by(native_invasive, plot, BRARrel, invasion_percent, grad_num) %>%
  summarise(avg_leaf = mean(avg_leaf_herb), avg_plant = mean(per_herb_plant)) %>%
  ungroup()



#stats
#2021
#whole plant herbivory
#normality
res_a <- lm(data = FK2021herb_join2, avg_plant ~ BRARrel * native_invasive)
ols_plot_resid_hist(res_a) # normalish
ols_test_normality(res_a) #pass KS
#same results w/inv %

#leaf herbivory
#normality
res_a <- lm(data = FK2021herb_join2, avg_leaf ~ invasion_percent * native_invasive)
ols_plot_resid_hist(res_a) # normalish
ols_test_normality(res_a) #pass KS
#same results w/inv %

#transformed
FK2021herb_join2_trans <- FK2021herb_join2 %>%
  mutate(ln = log10(avg_plant + 0.1)) %>%
  mutate(ln_lf = log10(avg_leaf + 0.1))

res_a <- lm(data = FK2021herb_join2_trans, ln ~ BRARrel * native_invasive)
ols_plot_resid_hist(res_a) # normal
ols_test_normality(res_a) #pass
#same results w/inv %

res_a <- lm(data = FK2021herb_join2_trans, ln_lf ~ BRARrel * native_invasive)
ols_plot_resid_hist(res_a) # normal
ols_test_normality(res_a) #pass
#same results w/inv %


#use log transformed

#factored - plant
FK2021herb_join2_trans2 <- FK2021herb_join2_trans %>%
  mutate(invasion_percent = as.factor(invasion_percent))

FK2021_pl_stats2 <- lmerTest::lmer(data = FK2021herb_join2_trans2, ln ~ invasion_percent * native_invasive + (1|grad_num))
anova(FK2021_pl_stats2, type = 3) #BRAR % p = 0.001917, nat_inv p = 2.172e-6

tuk21 <- emmeans(FK2021_pl_stats2, specs = pairwise ~ invasion_percent : native_invasive)
tuk21$contrasts

#linearity
plot(resid(FK2021_pl_stats2), FK2021herb_join2_trans2$ln) #looks linear
plot(FK2021_pl_stats2) #no pattern so indicates linearity

#homoscedascity
FK2021herb_join2_trans2$res <- residuals(FK2021_pl_stats2)
FK2021herb_join2_trans2$abs_res <- abs(FK2021herb_join2_trans2$res)
FK2021herb_join2_trans2$abs_res2 <- FK2021herb_join2_trans2$abs_res^2
levene_brar_herb21 <- lm(data = FK2021herb_join2_trans2, abs_res2 ~ invasion_percent)
anova(levene_brar_herb21) #p = 0.2568 so > 0.05 so equal variance is met

#autocorrelation
acf(residuals(FK2021_pl_stats2, retype = "normalized"))  #lines dont go outside CI horizontal lines, so not autocorrelated
pacf(residuals(FK2021_pl_stats2, retype = "normalized"))

#factored - leaf
FK2021_lf_stats2 <- lmerTest::lmer(data = FK2021herb_join2_trans2, ln_lf ~ invasion_percent * native_invasive + (1|grad_num))
anova(FK2021_lf_stats2, type = 3) #BRAR % p = 0.001009, nat_inv p = 0.004826

tuk21_2 <- emmeans(FK2021_lf_stats2, specs = pairwise ~ invasion_percent : native_invasive)
tuk21_2$contrasts




#2022
#whole plant herbivory
#normality
res_a <- lm(data = FK2022herb_join2, avg_plant ~ invasion_percent * native_invasive)
ols_plot_resid_hist(res_a) # normalish
ols_test_normality(res_a) #pass KS
#same results w/inv %

#leaf herbivory
#normality
res_a <- lm(data = FK2022herb_join2, avg_leaf ~ BRARrel * native_invasive)
ols_plot_resid_hist(res_a) # normalish
ols_test_normality(res_a) #pass KS
#same results w/inv %

#transformed
FK2022herb_join2_trans <- FK2022herb_join2 %>%
  mutate(ln = log10(avg_plant + 0.1)) %>%
  mutate(ln_lf = log10(avg_leaf + 0.1))

res_a <- lm(data = FK2022herb_join2_trans, ln ~ invasion_percent * native_invasive)
ols_plot_resid_hist(res_a) # normal
ols_test_normality(res_a) #pass
#same results w/inv %

res_a <- lm(data = FK2022herb_join2_trans, ln_lf ~ invasion_percent * native_invasive)
ols_plot_resid_hist(res_a) # normal
ols_test_normality(res_a) #pass
#same results w/inv %

#use log transformed
#factored - plant
FK2022herb_join2_trans2 <- FK2022herb_join2_trans %>%
  mutate(invasion_percent = as.factor(invasion_percent))

FK2022_pl_stats2 <- lmerTest::lmer(data = FK2022herb_join2_trans2, ln ~ invasion_percent * native_invasive + (1|grad_num))
anova(FK2022_pl_stats2, type = 3) #BRAR % p = 0.1313, nat_inv p = 3.674e-7

tuk22 <- emmeans(FK2022_pl_stats2, specs = pairwise ~ invasion_percent : native_invasive)
tuk22$contrasts

#linearity
plot(resid(FK2022_pl_stats2), FK2022herb_join2_trans2$ln) #looks linear
plot(FK2022_pl_stats2) #no pattern so indicates linearity

#homoscedascity
FK2022herb_join2_trans2$res <- residuals(FK2022_pl_stats2)
FK2022herb_join2_trans2$abs_res <- abs(FK2022herb_join2_trans2$res)
FK2022herb_join2_trans2$abs_res2 <- FK2022herb_join2_trans2$abs_res^2
levene_brar_herb22 <- lm(data = FK2022herb_join2_trans2, abs_res2 ~ invasion_percent)
anova(levene_brar_herb22) #p = 0.7633 so > 0.05 so equal variance is met

#autocorrelation
acf(residuals(FK2022_pl_stats2, retype = "normalized"))  #lines dont go outside CI horizontal lines, so not autocorrelated
pacf(residuals(FK2022_pl_stats2, retype = "normalized"))


#factored - leaf
FK2022_lf_stats2 <- lmerTest::lmer(data = FK2022herb_join2_trans2, ln_lf ~ invasion_percent * native_invasive + (1|grad_num))
anova(FK2022_lf_stats2, type = 3) #BRAR % p = 0.06684, nat_inv p = 1.125e-8

tuk22_2 <- emmeans(FK2022_lf_stats2, specs = pairwise ~ invasion_percent : native_invasive)
tuk22_2$contrasts


#on average, how much more herbivory is on natives than invasives regardless of inv level?
FK2021herb_join3 <- FK2021herb_join2 %>%
  group_by(native_invasive) %>%
  summarise(avg = mean(avg_plant), se = sd(avg_plant)/sqrt(length(avg_plant))) %>%
  ungroup()   #inv: 0.8402778 +/- 0.1405994; nat: 4.0616319 +/- 0.7811700

FK2022herb_join3 <- FK2022herb_join2 %>%
  group_by(native_invasive) %>%
  summarise(avg = mean(avg_plant), se = sd(avg_plant)/sqrt(length(avg_plant))) %>%
  ungroup()   #inv: 1.019097 +/- 0.2855309; nat: 4.359028 +/- 0.6510514



#figures
#whole plant
FK2021herb_fig3 <- ggplot(data = FK2021herb_join2, aes(x = factor(invasion_percent), y = avg_plant, color = native_invasive)) + 
  geom_boxplot() + 
  #xlim(0, 60) +
  ylim(0, 15) +
  scale_colour_manual(values = cbPalette) +
  labs(x = "Invasion Level (%)", y = "Total plant herbivory (%)", title = "MT BRAR 2021", color = "Plant Status") +
  theme(plot.title = element_text(vjust = - 7, hjust = 0.07, size = 16), legend.position = "none")

FK2022herb_fig3 <- ggplot(data = FK2022herb_join2, aes(x = factor(invasion_percent), y = avg_plant, color = native_invasive)) + 
  geom_boxplot() + 
  #xlim(0, 60) +
  ylim(0, 15) +
  scale_colour_manual(values = cbPalette) +
  labs(x = "Invasion Level (%)", y = "", title = "MT BRAR 2022", color = "Plant Status") +
  theme(plot.title = element_text(vjust = - 7, hjust = 0.07, size = 16))

FK2021herb_fig3 + FK2022herb_fig3

#herb_pl 1000 x 400


#leaf
FK2021herb_fig4 <- ggplot(data = FK2021herb_join2, aes(x = factor(invasion_percent), y = avg_leaf, color = native_invasive)) + 
  geom_boxplot() + 
  #xlim(0, 60) +
  ylim(0, 21) +
  scale_colour_manual(values = cbPalette) +
  labs(x = "Invasion Level (%)", y = "Leaf Herbivory (%)", title = "MT BRAR 2021", color = "Plant Status") +
  theme(plot.title = element_text(vjust = - 7, hjust = 0.07, size = 16), legend.position = "none")

FK2022herb_fig4 <- ggplot(data = FK2022herb_join2, aes(x = factor(invasion_percent), y = avg_leaf, color = native_invasive)) + 
  geom_boxplot() + 
  #xlim(0, 60) +
  ylim(0, 21) +
  scale_colour_manual(values = cbPalette) +
  labs(x = "Invasion Level (%)", y = "", title = "MT BRAR 2022", color = "Plant Status") +
  theme(plot.title = element_text(vjust = - 7, hjust = 0.07, size = 16))

FK2021herb_fig4 + FK2022herb_fig4

#herb_lf 1000 x 400





########################################################################################




######################################################################################
############## FK 2021/2022 Insect Biomass and Herbivory Combined Figure ########
#biomass
title4 <- expression(paste("MT 2021 ", italic("B. arvensis")))

FK2021bio_fig <- ggplot(data = FK2021bio_join, aes(x = BRARrel, y = total_biomass)) + 
  geom_point(size = 3) + 
  #xlim(0, 80) +
  ylim(0, 0.5) +
  labs(x = "Cover (%)", y = "Insect biomass (g)", title = title4) +
  theme(plot.title = element_text(vjust = - 7, hjust = 0.07, size = 16))

title5 <- expression(paste("MT 2022 ", italic("B. arvensis")))

FK2022bio_fig <- ggplot(data = FK2022bio_join, aes(x = BRARrel, y = total_biomass)) + 
  geom_point(size = 3) + 
  #xlim(0, 80) +
  ylim(0, 0.8) +
  labs(x = "Cover (%)", y = "", title = title5) +
  theme(plot.title = element_text(vjust = - 7, hjust = 0.07, size = 16))

#whole plant herb
FK2021herb_fig3 <- ggplot(data = FK2021herb_join2, aes(x = factor(invasion_percent), y = avg_plant, color = native_invasive)) + 
  geom_boxplot() + 
  #xlim(0, 60) +
  ylim(0, 15) +
  scale_colour_manual(values = cbPalette) +
  labs(x = "Invasion Level (%)", y = "Total plant herbivory (%)", title = title4, color = "Plant Status") +
  theme(plot.title = element_text(vjust = - 7, hjust = 0.07, size = 16), legend.position = "none")

FK2022herb_fig3 <- ggplot(data = FK2022herb_join2, aes(x = factor(invasion_percent), y = avg_plant, color = native_invasive)) + 
  geom_boxplot() + 
  #xlim(0, 60) +
  ylim(0, 15) +
  scale_colour_manual(values = cbPalette) +
  labs(x = "Invasion Level (%)", y = "", title = title5, color = "Plant Status") +
  theme(plot.title = element_text(vjust = - 7, hjust = 0.07, size = 16))

FK2021bio_fig + FK2022bio_fig + FK2021herb_fig3 + FK2022herb_fig3 + plot_layout(ncol = 2)

#bioherb 1000 x 800










##########################################################################################


