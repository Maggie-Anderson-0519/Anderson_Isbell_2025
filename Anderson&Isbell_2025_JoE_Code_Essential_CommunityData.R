##############
# Title: Complete data and code associated with "Experimental warming drives local grassland plant species loss"
# Author: Maggie Anderson
# Date: 3/21/2025
##############
#
# Copyright (c) [2025] [Maggie Anderson]
#
# This code is the intellectual property of the author. Unauthorized use, distribution, or reproduction of this code, in whole or in part, is prohibited without the express written permission of the author. 
#
##############
# INDEX
##############
# (1) Setup
# (2) Load data
# (3) Demonstration: calculate diversity & biomass responses
# (4) Calculate SPEI 
#   - Figure 1: Climate and SPEI trends during study 
# (5) Prepare data for summary figures and analyses/models
# (6) Calculate Coefficient of Variation (CV)
# (7) Model selection
# (8) Plot model coefficients
#   - Figure 2: Warming effects on community responses
#   - Figure 3: Dominance & A. gerardii responses to warming
#   - Figure 4: Functional group responses to warming (no year)
# (9) Supplementary figures
##############

################################################################################
# (1) Setup ####################################################################
################################################################################
# Clear workspace
rm(list = ls())

# Load pacman package if needed
# install.packages("pacman"); library(pacman)

# Load necessary packages
pacman::p_load(
  DescTools,    # %like% operator
  emmeans,      # model evaluation
  fossil,       # ecol.dist for distance matrices
  ggpubr,       # ggarrange() for multiple plots
  ggtext,       # subscripts for plot text
  ggthemes,     # "theme_few()" in ggplot
  grid,         # commo axis labels for multiple plots
  insight,      # mmrm model assessment
  lme4,         # modeling
  lmerTest,     # modeling
  mmrm,         # mixed models for repeated measures. 
  # More on covariance structures in mmrm found here: https://cran.r-project.org/web/packages/mmrm/vignettes/covariance.html
  ncddf,        # netcf data manipulation (SPEI)
  nlme,         # modeling
  raster,       # SPEI
  RColorBrewer, # nice plot colors
  tidyverse,    # data wrangling
  vegan)        # diversity indices

# Set working directory
setwd(" YOUR DIRECTORY NAME HERE") # store all data files cloned from the "Anderson_Isbell_2025" GitHub repo here

################################################################################
# (2) Load data ################################################################
################################################################################

# Community-level data
df <- read.csv("CommunityLevelDataClean_2017_22.csv", row.names =  1) # first column should be "Plot"
# Growing season temperature and precipitation data from Cedar Creek weather station
tp1 <- read.csv("e080_DailyClimateSummary.csv")
# SPEI data for Cedar Creek
spei_df <- read.csv("SPEI_data.csv")
# Soil moisture and temperature data for Cedar Creek (for supplement)
smt1722.mod <- read.csv("SoilMoisture_2017-22_Summarized.csv", row.names = 1)

################################################################################
# (3) Demonstration: calculate diversity & biomass responses ###################
################################################################################
# Because I provide summarized community data *not raw data* for this project, I will demonstrate how I calculated diversity using the "dune" dataset from the vegan package in R
# Source URL: https://cran.r-project.org/web/packages/BiodiversityR/readme/README.html
data(dune) # community data
    # "Data set dune is a community data set, where variables (columns) typically correspond to different species and data represents abundance of each species. Species names were abbreviated to eight characters, with for example Agrostol representing Agrostis stolonifera."
data(dune.env) # environmental data
    # "Data set dune.env is an environmental data set, where variables (columns) correspond to different descriptors (typically continuous and categorical variables) of the sample sites. One of the variables is Management, a categorical variable that describes different management categories, coded as BF (an abbreviation for biological farming), HF (hobby farming), NM (nature conservation management) and SF (standard farming)."

dune.demo <- cbind(dune.env, dune) # combine 
str(dune.demo) # view data
# species are in cols 6:35

# Create a copy of the dataframe 'df' and assign it to 'pa'
pa <- dune.demo
# Loop through columns 6 to 35 (all species) and convert each column to logical (TRUE/FALSE)
for(i in 6:35) {
  pa[,i] <- as.logical(pa[,i])
}

# Calculate the number of species present and store it in a new column 'S' 
dune.demo$S <- rowSums(pa[,6:35], na.rm = TRUE)

# Reshape the dataframe 'dune.demo' from wide to long format, pivoting the species columns
dune.demo.l <- dune.demo %>%
  pivot_longer(cols = Achimill:Callcusp)
dune.demo.l <- dune.demo.l[!is.na(dune.demo.l$value) & dune.demo.l$value > 0,] ## Remove rows where 'value' is NA or less than or equal to 0

# New dataframe for diversity indices
dune.demo.p <- dune.demo.l %>%
  group_by(Management, Use, Manure) %>% # dune.env grouping variables. In my data, I used Year, Plot, and Subplot here. 
  dplyr::summarise(Dq0=sum(value^0), # Richness
                   Dq1=exp(-sum(value*log(value))), # Exponential of Shannon Entropy: measures the uncertainty of a random process
                   Dq2=sum(value^2)^(1/(1-2)), # Inverse Simpson Diversity
                   BP=max(value), # Berger-Parker Dominance Index (from a proportion)
                   p.Ag=sum(value[name == 'Agrostol']), # Calculate proportional biomass of "Agrostol" from dune data. In my data, I replace this with 'Andropogon.gerardii' to get total Andropogon gerardii biomass
                   biomass=sum(value)) # proportional biomass of Andropogon gerardii

dune.demo <- dune.demo %>%
  left_join(dune.demo.p) # join diversity information to original dataframe

dune.demo$invSimpEven <- dune.demo$Dq2/dune.demo$Dq0 # evenness based on inverse simpson

# reorder dataframe 
dune.demo <- dune.demo[order(dune.demo$Management,dune.demo$Use,dune.demo$Manure),] # dune.env grouping variables. In my data, I used Year, Plot, and Subplot here.

# Calculate proportions of each functional group
# For dune, we will do this with legumes, forbs, C3 grasses, and shrubs. Ther eare no C4 grasses in the dune dataset. 
# Note that in my data, I calculate biomass proportion of C4 grasses instead of shrubs, as shown below

# Calculate proportions of each functional group
dune.demo <- dune.demo %>%
  mutate(
    Yo = rowSums(across(6:35), na.rm=TRUE), # total biomass
    LegumeBiomass = ifelse(Yo > 0, rowSums(across(c("Trifprat", "Trifrepe", "Vicilath")), na.rm=TRUE) / Yo, 0),
    ForbBiomass = ifelse(Yo > 0, rowSums(across(c("Achimill", "Anthodor", "Bellpere", "Chenalbu", "Cirsarve", 
                                                  "Empenigr", "Hyporadi", "Planlanc", "Ranuflam", "Rumeacet", 
                                                  "Sagiproc")), na.rm=TRUE) / Yo, 0),
    C3Biomass = ifelse(Yo > 0, rowSums(across(c("Agrostol", "Airaprae", "Alopgeni", "Bromhord", "Elymrepe", 
                                                "Lolipere", "Poaprat", "Poatriv")), na.rm=TRUE) / Yo, 0),
    ShrubBiomass = ifelse(Yo > 0, rowSums(across(c("Salirepe", "Callcusp")), na.rm=TRUE) / Yo, 0) # Substitute for C4 biomass in my data
  ) %>%
  dplyr::select(-(Achimill:Callcusp)) # drop species columns we no longer need

# End demonstration. 
# note that our dune.demo dataset now has the following response variables:
str(dune.demo)

# note the response variable columns present in the community dataset available for this project:
str(df)

################################################################################
# (4) Summarize temperature data & calculate SPEI ##############################
################################################################################

# reformat date column
tp1$date <- as.Date(tp1$Date, "%m/%d/%Y")

# make year/month/day columns
tp1$year <- as.numeric(str_sub(tp1$date, end=4))
tp1$month <- as.numeric(str_sub(tp1$date, start=6, end=7))
tp1$day <- as.numeric(str_sub(tp1$date, start=9, end=10))

# alter variables
tp1$MinTemp.degC. <- (tp1$MinTemp.degF. - 32) * 5/9
tp1$MaxTemp.degC. <- (tp1$MaxTemp.degF. - 32) * 5/9
tp1$Precip.mm. <- tp1$Precip.inches. * 25.4

# filter to data from 2017 - 2020 and May - September
# convert temperature to deg. C and precipitation to mm
# summarize to get average daily temperature and 
tp2 <- tp1 %>%
  filter(year %in% c(2017,2018,2019,2020,2021,2022)) %>%
  filter(month %in% c(5,6,7,8)) %>%
  dplyr::select(-c(Date,date,day,MinTemp.degF.,MaxTemp.degF.,Precip.inches.)) # cols we no longer need

# get summaries of min and max temperature
tp2.temp <- tp2 %>%
  dplyr::select(-Precip.mm.) %>%
  group_by(year, month, TempSource, PrecipSource) %>%
  summarise_all(mean) %>%
  pivot_longer(cols = MinTemp.degC.:MaxTemp.degC., 
               names_to = "Measure",
               values_to = "Temp.degC.") %>%
  mutate(Measure = ifelse(Measure == "MinTemp.degC.", "MinTemp",
                          ifelse(Measure == "MaxTemp.degC.", "MaxTemp",NA)))

# max precip 
tp2.precip <- tp2 %>%
  dplyr::select(-c(MinTemp.degC.,MaxTemp.degC.,month)) %>%
  group_by(year, TempSource, PrecipSource) %>%
  summarise_all(sum)

# summaries of growing season historical temps
# summer high temps
tp3.temp.gs <- tp1 %>%
  filter(year >= 1966 & year <= 2016) %>% # use data from 50 years prior to study
  filter(month %in% c(5,6,7,8)) %>% # May - August (instead of September)
  dplyr::select(MaxTemp.degC., MinTemp.degC.) %>%
  dplyr::summarize(mean.max.summer.temp = mean(MaxTemp.degC.,na.rm=T),
                   n.max=n(), 
                   sd.max = sd(MaxTemp.degC.,na.rm=T), 
                   se.max = sd.max/sqrt(n.max),
                   mean.min.summer.temp = mean(MinTemp.degC.,na.rm=T),
                   n.min=n(), 
                   sd.min = sd(MinTemp.degC.,na.rm=T), 
                   se.min = sd.min/sqrt(n.min));tp3.temp.gs

# winter low temps
tp3.temp.low <- tp1 %>%
  filter(year >= 1966 & year <= 2016) %>%
  filter(month %in% c(12,1,2)) %>% # December - February
  dplyr::select(MinTemp.degC.) %>%
  dplyr::summarize(mean.min.winter.temp = mean(MinTemp.degC.,na.rm=T),
                   n=n(), 
                   sd = sd(MinTemp.degC.,na.rm=T), 
                   se = sd/sqrt(n));tp3.temp.low

# yearly precip
tp3.precip <- tp1 %>%
  filter(year >= 1966 & year <= 2016) %>%
  dplyr::select(year,Precip.mm.) %>%
  group_by(year) %>%
  summarise_all(sum) %>%
  dplyr::select(-year) %>%
  dplyr::summarize(mean.precip = mean(Precip.mm.,na.rm=T),
                   n=n(), 
                   sd = sd(Precip.mm.,na.rm=T), 
                   se = sd/sqrt(n));tp3.precip

# growing season precip
tp3.precip.gs <- tp1 %>%
  filter(year >= 1966 & year <= 2016) %>%
  filter(month %in% c(5,6,7,8)) %>% # May - August (instead of September)
  dplyr::select(year,Precip.mm.) %>%
  dplyr::group_by(year) %>%
  summarise_all(sum) %>%
  dplyr::select(-year) %>%
  dplyr::summarize(mean.precip = mean(Precip.mm.,na.rm=T),
                   n=n(), 
                   sd = sd(Precip.mm.,na.rm=T), 
                   se = sd/sqrt(n));tp3.precip.gs

# SPEI Calculations
spei1.1 <- spei_df %>%
  data.frame() %>%
  dplyr::filter(Month == 8) # Extract data from August (when biomass clipping occurred)

spei1.1$SPEI_4cat <- "Normal"
qnorm(c(0.1,0.25,0.75,0.9)) # norm function provides the quantile of the normal distribution at a specified cumulative density (Isbell et al. 2015, Nature)

spei1.1$SPEI_4cat[spei1.1$SPEI >= 1.2815516] <- "Extremely wet" 
spei1.1$SPEI_4cat[spei1.1$SPEI >= 0.6744898 & spei1.1$SPEI < 1.2815516] <- "Wet"
spei1.1$SPEI_4cat[spei1.1$SPEI <= -1.2815516] <- "Extremely dry"
spei1.1$SPEI_4cat[spei1.1$SPEI < -0.6744898 & spei1.1$SPEI > -1.2815516] <- "Dry"

unique(spei1.1$SPEI_4cat)

spei2 <- spei1.1 %>%
  rename(SPEI04 = SPEI) %>%
  filter(2017 <= Year & Year <= 2023)   # filter data to years of interest

# Figure 1: Climate and SPEI trends during study ###############################

palette1<- c("red4", "coral1")

# Temp and precip over time (Cedar Creek)
tp2.temp$year <- as.factor(tp2.temp$year)
t<-ggplot(tp2.temp, aes(x = year, y = Temp.degC., col = Measure)) +
  stat_summary(fun=mean,geom="point", position = position_dodge(width = 0.5), size=2) +
  stat_summary(fun.data = mean_se,  geom = "errorbar", width = 0.25, size = 0.75, position = position_dodge(width = 0.5)) +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  ylab(bquote('Growing season temperature,'~degree*C)) +
  xlab("") +
  geom_hline(yintercept = tp3.temp.gs$mean.max.summer.temp,col="red4",linetype="dashed") +
  geom_hline(yintercept = tp3.temp.gs$mean.min.summer.temp,col="coral1",linetype="dashed") +
  scale_color_manual(values=palette1, labels=c("MaxTemp" = "Maximum", "MinTemp" = "Minimum")) +
  theme(strip.background=element_rect(fill="white"),
        legend.position = c(0.75,0.5),
        legend.title = element_blank(),
        legend.background = element_rect(fill="white",
                                         size=0.5, linetype="solid", 
                                         colour ="grey95"));t

tp2.precip$year <- as.factor(tp2.precip$year)
p<-ggplot(tp2.precip, aes(x = year, y = Precip.mm.)) +
  stat_summary(fun=mean,geom="bar", col = "black", fill = "steelblue", position = position_dodge(width = 0.5), size=0.5, width=0.25) +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  geom_hline(yintercept = tp3.precip.gs$mean.precip,col="steelblue4",linetype="dashed") +
  ylab(bquote('Total growing season precipitation,'~mm)) +
  xlab("") +
  theme(strip.background=element_rect(fill="white"),legend.position = "top")

spei.plot <- ggplot(spei2, aes(x = Year, y = SPEI04)) +
  stat_summary(fun=mean,geom="bar", col = "black", fill = "grey40", position = position_dodge(width = 0.5), size=0.5, width=0.25) +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  #geom_hline(yintercept = spei_mean$SPEI04, col="grey20",linetype="dashed") +
  geom_hline(yintercept = 0, col="black",linetype="solid", size = 0.5) +
  ylab(bquote('SPEI (May - August)')) +
  xlab("Year") +
  theme(strip.background=element_rect(fill="white"),legend.position = "top");spei.plot
ggarrange(t+rremove("xlab"),p+rremove("xlab"),spei.plot,ncol=1)

################################################################################
# (5) Prepare data for summary figures and analyses/models #####################
################################################################################

# Function for calculating standard error
calcSE<-function(x){
  sd(x)/sqrt(length(x))}

# change df to long format
dfl <- df %>%
  dplyr::select(Year:Sorghastrum.nutans) %>%
  pivot_longer(cols = Achillea.millefolium:Sorghastrum.nutans, 
               names_to = "Species", values_to = "TotalBiomass.g.m2") %>%
  na.omit() %>%
  mutate(Trtmt = ifelse(HeatTrt=="Heated" & WaterTrt=="Control", "+Heat",
                        ifelse(HeatTrt=="Control" & WaterTrt=="Drought","-Water",
                               ifelse(HeatTrt=="Heated" & WaterTrt=="Drought","+Heat-Water","Control")))) # make special combined treatment column (will make more sense later)

# make functional group column for new df (dfl)
dfl$FunctionalGroup <- ifelse(dfl$Species %in% c("Lespedeza.capitata", "Lupinus.perennis","Amorpha.Dalea"), "Legume",
                              ifelse(dfl$Species %in% c("Achillea.millefolium","Monarda.Oligoneuron","Asclepias.tuberosa","Liatris.aspera"), "Forb",
                                     ifelse(dfl$Species %in% c("Pascopyrum.smithii","Elymus.canadensis","Koeleria.macrantha","Poa.pratensis"), "C3",
                                            ifelse(dfl$Species %in% c("Andropogon.gerardii","Schizachyrium.scoparium","Sorghastrum.nutans","Panicum.virgatum"), "C4",NA))))

# Log-transform relative & total biomass
dfl <- dfl %>%
  mutate(TotalBiomassLog = log1p(TotalBiomass.g.m2))

# get mean species biomass per treatment
Means <- dfl %>%
  dplyr::filter(RichnessTrt == 16) %>%
  dplyr::group_by(Trtmt, HeatTrt, WaterTrt, Species, FunctionalGroup, Year) %>%
  dplyr::select(-c(Subplot,RichnessTrt)) %>%
  summarise_each(funs(mean,sd,se=sd(.)/sqrt(n())))

Means.all <- dfl %>%
  filter(RichnessTrt == 16) %>%
  group_by(Trtmt, HeatTrt, WaterTrt, Species, FunctionalGroup) %>%
  dplyr::select(-c(Subplot,RichnessTrt,Year)) %>%
  summarise_each(funs(mean,sd,se=sd(.)/sqrt(n())))

Means.trt.yr <- df %>%
  filter(RichnessTrt == 16) %>%
  dplyr::select(Year,HeatTrt,WaterTrt,Dq0,Dq1,Dq2,BP,p.Ag,invSimpEven,LegumeBiomass,ForbBiomass,C3Biomass,C4Biomass) %>%
  group_by(HeatTrt, WaterTrt, Year) %>%
  summarise_each(funs(mean,sd,se=sd(.)/sqrt(n())))

# Re-name and recategorize treatment responses for plotting
df$Trtmt<- as.factor(ifelse(df$WaterTrt == "Control" & df$HeatTrt == "Control", "Control",
                            ifelse(df$WaterTrt == "Drought" & df$HeatTrt == "Control","Drought",
                                   ifelse(df$WaterTrt == "Control" & df$HeatTrt == "Heated", "Warming",
                                          ifelse(df$WaterTrt == "Drought" & df$HeatTrt == "Heated","Warming*Drought",NA)))))
df$WaterTrt <- ifelse(df$WaterTrt == "Drought","-Water","Control")
df$HeatTrt <- ifelse(df$HeatTrt == "Heated","+Heat","Control")
df$WaterTrt <- factor(df$WaterTrt, levels = c("Control","-Water"))
df$HeatTrt <- factor(df$HeatTrt, levels = c("Control","+Heat"))

# # make summary dataset for responses
# df.sum <- df %>%
#   dplyr::filter(RichnessTrt == 16) %>%
#   dplyr::select(Year,HeatTrt,WaterTrt,Dq0,Dq1,Dq2,BP,p.Ag,invSimpEven,LegumeBiomass,ForbBiomass,C3Biomass,C4Biomass) %>%
#   dplyr::group_by(HeatTrt, WaterTrt, Year) %>%
#   dplyr::summarise_each(funs(mean,sd,se=sd(.)/sqrt(n())))
# 
# # make summary dataset for WARMING ONLY across years
# df.sum.heat <- df %>%
#   dplyr::filter(RichnessTrt == 16) %>%
#   dplyr::select(Year,HeatTrt,Dq0,Dq1,Dq2,BP,p.Ag,invSimpEven,LegumeBiomass,ForbBiomass,C3Biomass,C4Biomass) %>%
#   dplyr::group_by(HeatTrt, Year) %>%
#   dplyr::summarise_each(funs(mean,sd,se=sd(.)/sqrt(n())))
# 
# # make log response ratio dataset 
# df.lrr <- df %>%
#   dplyr::filter(RichnessTrt == 16) %>%
#   dplyr::select(Year,Plot,Trtmt,Dq0,Dq1,Dq2,BP,p.Ag,invSimpEven,LegumeBiomass,ForbBiomass,C3Biomass,C4Biomass) %>%
#   pivot_longer(cols = Dq0:C4Biomass, names_to = "var", values_to = "response") %>%
#   pivot_wider(names_from = Trtmt, values_from = response) %>%
#   dplyr::mutate(DroughtLRR = log(Drought / Control),
#                 WarmingLRR = log(Warming / Control),
#                 `Warming*DroughtLRR` = log(`Warming*Drought` / Control)) %>%
#   dplyr::select(Year, Plot, var, DroughtLRR,WarmingLRR,`Warming*DroughtLRR`) %>%
#   pivot_longer(cols = DroughtLRR:`Warming*DroughtLRR`, names_to = "TrtResp", values_to = "LRR") %>%
#   pivot_wider(names_from = var, values_from = "LRR")

################################################################################
# (6) Calculate Coefficient of Variation (CV) ##################################
################################################################################
# Reformat df for CV calculations
df.cv <- df %>%
  dplyr::filter(RichnessTrt == 16) %>%
  dplyr::select(HeatTrt,WaterTrt,Plot,Subplot,Year,Dq0,Dq1,Dq2,BP,p.Ag,invSimpEven,LegumeBiomass,ForbBiomass,C3Biomass,C4Biomass) %>%
  pivot_longer(cols=Dq0:C4Biomass,
               names_to = "Measure",
               values_to = "Value")

# Calculate mean and standard deviation for each Treatment-Year combination
cv_results <- df.cv %>%
  dplyr::group_by(HeatTrt, WaterTrt, Plot, Subplot, Measure) %>%
  dplyr::summarize(mean_response = mean(Value),
                   sd_response = sd(Value)) %>%
  dplyr::mutate(cv = (sd_response / mean_response))

# Print the results
head(cv_results)

################################################################################
# (7) Model selection ##########################################################
################################################################################

# Filter dataframe to focus on 16-species mixtures, since this allows us to evaluate community-level diversity responses
df <- df[df$RichnessTrt==16,] 

# Treat year as a factor in fixed effects, this helps get DenDF correct and the repeated variable must be a factor according to Phil Dixon: https://pdixon.stat.iastate.edu/stat571/labs/6%20Nov%20Rep%20meas/asparagus.pdf
df$YrFac <- as.factor(as.character(df$Year)) # factor

# Create a new subplot variable that does not repeat values
df$sp <- as.factor(paste(df$Plot,df$Subplot,sep=".")) 

# Initialize an empty dataframe to store the combined results
all_coefs <- data.frame() # coefficients for main figures
all_coefs_sup <- data.frame() # coefficietns for supplemental figures
aic_all <- data.frame() # AICc values for all tested models

# Our results depend on temporal autocorrelation structure, therefore, need to evaluate which one is most parsimonious in addition to evaluating which fixed effects relationship is most parsimonious
# We will do this for each response variable in our study:

#### Models for community responses
{
# Begin model selection
  # Note that in some cases the models are commented out because they throw the optimizer error "L-BFGS-B needs finite values of 'fn'" indicating that model parameters are converging to values that return non-finite output. In these cases, we assume a poor model fit and move on to select a different model which converges.
  ################################################################################
  # Richness (Dq0)
  ################################################################################
  #first, find most parsimonious temporal autocorrelation structure
  fit1.Dq0.factorYr.mmrm <- mmrm(Dq0~HeatTrt*WaterTrt*YrFac  + ar1(YrFac|sp), data=df )   # ar1 = homogeneous autoregressive covariance strs assume that "the correlation between any two elements is equal to r for adjacent elements". In other words, covariances between time points decline exponentially.
  #fit2.Dq0.factorYr.mmrm <- mmrm(Dq0~HeatTrt*WaterTrt*YrFac  + cs(YrFac|sp), data=df )    # Optimizer error# cs = homogeneous compound symmetry assumes a constant correlation between time points
  fit2.5.Dq0.factorYr.mmrm <- mmrm(Dq0~HeatTrt*WaterTrt*YrFac + csh(YrFac|sp), data=df ) # csh = heterogeneous compound symmetry 
  fit3.Dq0.factorYr.mmrm <- mmrm(Dq0~HeatTrt*WaterTrt*YrFac + us(YrFac|sp), data=df )    # us = unstructured covariance str
  fit4.Dq0.factorYr.mmrm <- mmrm(Dq0~HeatTrt*WaterTrt*YrFac + ad(YrFac|sp), data=df )    # ad = homogeneous covariance str, observations are approx. equally spaced in time 
  fit5.Dq0.factorYr.mmrm <- mmrm(Dq0~HeatTrt*WaterTrt*YrFac  + ar1h(YrFac|sp), data=df )  # ar1h = same as hr 1, but assumes heterogeneous variance
  AIC(fit1.Dq0.factorYr.mmrm) 
  #AIC(fit2.Dq0.factorYr.mmrm)
  AIC(fit2.5.Dq0.factorYr.mmrm) # best fit
  AIC(fit3.Dq0.factorYr.mmrm)
  AIC(fit4.Dq0.factorYr.mmrm)
  AIC(fit5.Dq0.factorYr.mmrm)
  # get delta AIC
  aic_covariance <- data.frame(
    Model = character(),
    AIC = numeric(),
    Delta_AIC = numeric(),
    stringsAsFactors = FALSE
  )
  
  # List of models
  models_covariance <- list(
    fit1.Dq0.factorYr.mmrm,
    #fit2.Dq0.factorYr.mmrm,
    fit2.5.Dq0.factorYr.mmrm,
    fit3.Dq0.factorYr.mmrm,
    fit4.Dq0.factorYr.mmrm,
    fit5.Dq0.factorYr.mmrm
  )
  
  # Model names
  model_names_covariance <- c(
    "AR1",
    #"CS",
    "CSH",
    "US",
    "AD",
    "AR1H"
  )
  
  # Calculate AIC and delta-AIC for each model
  aic_values <- sapply(models_covariance, AIC)
  min_aic <- min(aic_values)
  delta_aic <- aic_values - min_aic
  
  # Populate the dataframe
  aic_covariance <- data.frame(
    Model = model_names_covariance,
    AIC = aic_values,
    Delta_AIC = delta_aic,
    response = "Richness",  
    selection_stage = "v-cov"
  )
  
  aic_covariance <- aic_covariance %>% arrange(AIC)  # Sort by AIC
  aic_all <- bind_rows(aic_all, aic_covariance)  # bind into main table
  
  # do this for the csh model (new lowest AIC)
  # look at the anova table and double check NumDF and DenDF
  Dq0.emm2.5 <- emmeans(fit2.5.Dq0.factorYr.mmrm, c('HeatTrt','WaterTrt','YrFac'))
  joint_tests(Dq0.emm2.5) 
  
  Dq0.emm4 <- emmeans(fit2.5.Dq0.factorYr.mmrm,c('HeatTrt'))
  Dq0.postHoc <- multcomp::cld(object = Dq0.emm4,Letters = letters,alpha = 0.05);Dq0.postHoc 
  
  # get treatment effect:
  Dq0.contrasts <- data.frame(summary(contrast(Dq0.emm2.5, method = "trt.vs.ctrl", by = c("WaterTrt", "YrFac"))))
  Dq0.contrasts.all <- data.frame(summary(contrast(Dq0.emm2.5, method = "trt.vs.ctrl", by = c("YrFac"))))
  
  # Get most parsimonious model form
  fit2.5.Dq0.factorYr.mmrm.a <- mmrm(Dq0~HeatTrt*WaterTrt*YrFac +  + csh(YrFac|sp), data=df )
  fit2.5.Dq0.factorYr.mmrm.b <- mmrm(Dq0~HeatTrt*WaterTrt+YrFac + csh(YrFac|sp), data=df )
  fit2.5.Dq0.factorYr.mmrm.c <- mmrm(Dq0~HeatTrt+WaterTrt+YrFac + csh(YrFac|sp), data=df )
  fit2.5.Dq0.factorYr.mmrm.d <- mmrm(Dq0~HeatTrt*WaterTrt+YrFac - YrFac:HeatTrt:WaterTrt  + csh(YrFac|sp), data=df )
  fit2.5.Dq0.factorYr.mmrm.e <- mmrm(Dq0~HeatTrt+WaterTrt+YrFac - YrFac:HeatTrt:WaterTrt  + csh(YrFac|sp), data=df )
  fit2.5.Dq0.factorYr.mmrm.f <- mmrm(Dq0~HeatTrt*WaterTrt*YrFac - YrFac:HeatTrt:WaterTrt + csh(YrFac|sp), data=df )
  AIC(fit2.5.Dq0.factorYr.mmrm.a) # fully-interactive model fits best
  AIC(fit2.5.Dq0.factorYr.mmrm.b)
  AIC(fit2.5.Dq0.factorYr.mmrm.c)
  AIC(fit2.5.Dq0.factorYr.mmrm.d)
  AIC(fit2.5.Dq0.factorYr.mmrm.e)
  AIC(fit2.5.Dq0.factorYr.mmrm.f)
  
  ########################################################################
  # Delta-AIC
  aic_parsimony <- data.frame(
    Model = character(),
    AIC = numeric(),
    Delta_AIC = numeric(),
    stringsAsFactors = FALSE
  )
  
  # List of models
  models_parsimony <- list(
    fit2.5.Dq0.factorYr.mmrm.a,
    fit2.5.Dq0.factorYr.mmrm.b,
    fit2.5.Dq0.factorYr.mmrm.c,
    fit2.5.Dq0.factorYr.mmrm.d,
    fit2.5.Dq0.factorYr.mmrm.e,
    fit2.5.Dq0.factorYr.mmrm.f
    
  )
  
  # Model names
  model_names_parsimony <- c(
    "Full Interactive",
    "Year additive",
    "Fully additive",
    "Full Interactive - no 3-way",
    "Year additive - no 3-way",
    "Fully additive - no 3-way"
  )
  
  # Calculate AIC and delta-AIC for each model
  aic_values_parsimony <- sapply(models_parsimony, AIC)
  min_aic_parsimony <- min(aic_values_parsimony)
  delta_aic_parsimony <- aic_values_parsimony - min_aic_parsimony
  
  # Populate the dataframe
  aic_parsimony <- data.frame(
    Model = model_names_parsimony,
    AIC = aic_values_parsimony,
    Delta_AIC = delta_aic_parsimony,
    response = "Richness",  # Repeat "Dq0" for each row
    selection_stage = "parsimony"
  )
  
  # Sort by AIC
  aic_parsimony <- aic_parsimony %>% arrange(AIC)
  
  # bind into main table
  aic_all <- bind_rows(aic_all, aic_parsimony)
  # Best fit for Dq0 is a CSH fully-interactive model (tested above)
  # extract coefficients
  fit_summary.Dq0 <- data.frame(summary(fit2.5.Dq0.factorYr.mmrm.a)$coefficients) %>% mutate(response = rep("Dq0"))
  
  ################################################################################
  # Exponential of Shannon's entropy (Dq1)
  ################################################################################
  fit1.Dq1.factorYr.mmrm <- mmrm(Dq1~HeatTrt*WaterTrt*YrFac  + ar1(YrFac|sp), data=df )
  fit2.Dq1.factorYr.mmrm <- mmrm(Dq1~HeatTrt*WaterTrt*YrFac + cs(YrFac|sp), data=df ) #best fit = lowest AIC
  fit2.5.Dq1.factorYr.mmrm <- mmrm(Dq1~HeatTrt*WaterTrt*YrFac + csh(YrFac|sp), data=df ) 
  fit3.Dq1.factorYr.mmrm <- mmrm(Dq1~HeatTrt*WaterTrt*YrFac + us(YrFac|sp), data=df )
  fit4.Dq1.factorYr.mmrm <- mmrm(Dq1~HeatTrt*WaterTrt*YrFac + ad(YrFac|sp), data=df )
  fit5.Dq1.factorYr.mmrm <- mmrm(Dq1~HeatTrt*WaterTrt*YrFac + ar1h(YrFac|sp), data=df )
  AIC(fit1.Dq1.factorYr.mmrm)
  AIC(fit2.Dq1.factorYr.mmrm) # best fit
  AIC(fit2.5.Dq1.factorYr.mmrm) 
  AIC(fit3.Dq1.factorYr.mmrm)
  AIC(fit4.Dq1.factorYr.mmrm)
  AIC(fit5.Dq1.factorYr.mmrm)
  # get delta AIC
  aic_covariance <- data.frame(
    Model = character(),
    AIC = numeric(),
    Delta_AIC = numeric(),
    stringsAsFactors = FALSE
  )
  
  # List of models
  models_covariance <- list(
    fit1.Dq1.factorYr.mmrm,
    fit2.Dq1.factorYr.mmrm,
    fit2.5.Dq1.factorYr.mmrm,
    fit3.Dq1.factorYr.mmrm,
    fit4.Dq1.factorYr.mmrm,
    fit5.Dq1.factorYr.mmrm
  )
  
  # Model names
  model_names_covariance <- c(
    "AR1",
    "CS",
    "CSH",
    "US",
    "AD",
    "AR1H"
  )
  
  # Calculate AIC and delta-AIC for each model
  aic_values <- sapply(models_covariance, AIC)
  min_aic <- min(aic_values)
  delta_aic <- aic_values - min_aic
  
  # Populate the dataframe
  aic_covariance <- data.frame(
    Model = model_names_covariance,
    AIC = aic_values,
    Delta_AIC = delta_aic,
    response = "Shannon entropy",  
    selection_stage = "v-cov"
  )
  
  aic_covariance <- aic_covariance %>% arrange(AIC)  # Sort by AIC
  aic_all <- bind_rows(aic_all, aic_covariance)  # bind into main table
  
  Dq1.emm1 <- emmeans(fit2.Dq1.factorYr.mmrm,c('HeatTrt','WaterTrt','YrFac'))
  joint_tests(Dq1.emm1) 
  Dq1.emm2 <- emmeans(fit2.Dq1.factorYr.mmrm,c('YrFac','HeatTrt')) #warming significantly decreased diversity in some years (2017, 2019)
  Dq1.postHoc <- multcomp::cld(object = Dq1.emm2,Letters = letters,alpha = 0.05);Dq1.postHoc # add letters to each mean
  
  # get treatment effect:
  Dq1.contrasts <- data.frame(summary(contrast(Dq1.emm1, method = "trt.vs.ctrl", by = c("WaterTrt", "YrFac"))))
  Dq1.contrasts.all <- data.frame(summary(contrast(Dq1.emm1, method = "trt.vs.ctrl", by = c("YrFac"))))
  
  # Get most parsimonious model form
  fit2.Dq1.factorYr.mmrm.a <- mmrm(Dq1~HeatTrt*WaterTrt*YrFac + cs(YrFac|sp), data=df )
  fit2.Dq1.factorYr.mmrm.b <- mmrm(Dq1~HeatTrt*WaterTrt+YrFac + cs(YrFac|sp), data=df )
  fit2.Dq1.factorYr.mmrm.c <- mmrm(Dq1~HeatTrt+WaterTrt+YrFac + cs(YrFac|sp), data=df )
  fit2.Dq1.factorYr.mmrm.d <- mmrm(Dq1~HeatTrt*WaterTrt+YrFac - YrFac:HeatTrt:WaterTrt  + cs(YrFac|sp), data=df )
  fit2.Dq1.factorYr.mmrm.e <- mmrm(Dq1~HeatTrt+WaterTrt+YrFac - YrFac:HeatTrt:WaterTrt  + cs(YrFac|sp), data=df )
  fit2.Dq1.factorYr.mmrm.f <- mmrm(Dq1~HeatTrt*WaterTrt*YrFac - YrFac:HeatTrt:WaterTrt + cs(YrFac|sp), data=df )
  AIC(fit2.Dq1.factorYr.mmrm.a) # fully-interactive model fits best
  AIC(fit2.Dq1.factorYr.mmrm.b)
  AIC(fit2.Dq1.factorYr.mmrm.c)
  AIC(fit2.Dq1.factorYr.mmrm.d)
  AIC(fit2.Dq1.factorYr.mmrm.e)
  AIC(fit2.Dq1.factorYr.mmrm.f)
  
  # Delta-AIC
  aic_parsimony <- data.frame(
    Model = character(),
    AIC = numeric(),
    Delta_AIC = numeric(),
    stringsAsFactors = FALSE
  )
  
  # List of models
  models_parsimony <- list(
    fit2.Dq1.factorYr.mmrm.a,
    fit2.Dq1.factorYr.mmrm.b,
    fit2.Dq1.factorYr.mmrm.c,
    fit2.Dq1.factorYr.mmrm.d,
    fit2.Dq1.factorYr.mmrm.e,
    fit2.Dq1.factorYr.mmrm.f
  )
  
  # Model names
  model_names_parsimony <- c(
    "Full Interactive",
    "Year additive",
    "Fully additive",
    "Full Interactive - no 3-way",
    "Year additive - no 3-way",
    "Fully additive - no 3-way"
  )
  
  # Calculate AIC and delta-AIC for each model
  aic_values_parsimony <- sapply(models_parsimony, AIC)
  min_aic_parsimony <- min(aic_values_parsimony)
  delta_aic_parsimony <- aic_values_parsimony - min_aic_parsimony
  
  # Populate the dataframe
  aic_parsimony <- data.frame(
    Model = model_names_parsimony,
    AIC = aic_values_parsimony,
    Delta_AIC = delta_aic_parsimony,
    response = "Shannon entropy",  # Repeat "Dq0" for each row
    selection_stage = "parsimony"
  )
  
  # Sort by AIC
  aic_parsimony <- aic_parsimony %>% arrange(AIC)
  
  # bind into main table
  aic_all <- bind_rows(aic_all, aic_parsimony)
  # Best fit for Dq1 is a CS fully-interactive model (tested above)
  # extract coefficients
  fit_summary.Dq1 <- data.frame(summary(fit2.Dq1.factorYr.mmrm.a)$coefficients) %>% mutate(response = rep("Dq1")) 
  
  Dq1.emm1 <- emmeans(fit2.Dq1.factorYr.mmrm.a,c('HeatTrt','WaterTrt','YrFac'))
  joint_tests(Dq1.emm1) 
  Dq1.emm2 <- emmeans(fit2.Dq1.factorYr.mmrm.a,c('HeatTrt','YrFac')) 
  Dq1.postHoc <- multcomp::cld(object = Dq1.emm2,Letters = letters,alpha = 0.05); Dq1.postHoc 
  
  ################################################################################
  # Inverse Simpson's diversity (Dq2)
  ################################################################################
  fit1.Dq2.factorYr.mmrm <- mmrm(Dq2~HeatTrt*WaterTrt*YrFac  + ar1(YrFac|sp), data=df)
  fit2.Dq2.factorYr.mmrm <- mmrm(Dq2~HeatTrt*WaterTrt*YrFac  + cs(YrFac|sp), data=df) #best fit = lowest AIC
  fit2.5.Dq2.factorYr.mmrm <- mmrm(Dq2~HeatTrt*WaterTrt*YrFac  + csh(YrFac|sp), data=df)
  fit3.Dq2.factorYr.mmrm <- mmrm(Dq2~HeatTrt*WaterTrt*YrFac + us(YrFac|sp), data=df)
  fit4.Dq2.factorYr.mmrm <- mmrm(Dq2~HeatTrt*WaterTrt*YrFac  + ad(YrFac|sp), data=df)
  fit5.Dq2.factorYr.mmrm <- mmrm(Dq2~HeatTrt*WaterTrt*YrFac  + ar1h(YrFac|sp), data=df)
  AIC(fit1.Dq2.factorYr.mmrm)
  AIC(fit2.Dq2.factorYr.mmrm) # best fit
  AIC(fit2.5.Dq2.factorYr.mmrm) 
  AIC(fit3.Dq2.factorYr.mmrm)
  AIC(fit4.Dq2.factorYr.mmrm)
  AIC(fit5.Dq2.factorYr.mmrm)
  # get delta AIC
  aic_covariance <- data.frame(
    Model = character(),
    AIC = numeric(),
    Delta_AIC = numeric(),
    stringsAsFactors = FALSE
  )
  
  # List of models
  models_covariance <- list(
    fit1.Dq2.factorYr.mmrm,
    fit2.Dq2.factorYr.mmrm,
    fit2.5.Dq2.factorYr.mmrm,
    fit3.Dq2.factorYr.mmrm,
    fit4.Dq2.factorYr.mmrm,
    fit5.Dq2.factorYr.mmrm
  )
  
  # Model names
  model_names_covariance <- c(
    "AR1",
    "CS",
    "CSH",
    "US",
    "AD",
    "AR1H"
  )
  
  # Calculate AIC and delta-AIC for each model
  aic_values <- sapply(models_covariance, AIC)
  min_aic <- min(aic_values)
  delta_aic <- aic_values - min_aic
  
  # Populate the dataframe
  aic_covariance <- data.frame(
    Model = model_names_covariance,
    AIC = aic_values,
    Delta_AIC = delta_aic,
    response = "Inverse Simpson diversity",
    selection_stage = "v-cov"
  )
  
  aic_covariance <- aic_covariance %>% arrange(AIC)  # Sort by AIC
  aic_all <- bind_rows(aic_all, aic_covariance)  # bind into main table
  
  
  # Get most parsimonious model form
  fit2.Dq2.factorYr.mmrm.a <- mmrm(Dq2~HeatTrt*WaterTrt*YrFac + cs(YrFac|sp), data=df )
  fit2.Dq2.factorYr.mmrm.b <- mmrm(Dq2~HeatTrt*WaterTrt+YrFac + cs(YrFac|sp), data=df )
  fit2.Dq2.factorYr.mmrm.c <- mmrm(Dq2~HeatTrt+WaterTrt+YrFac + cs(YrFac|sp), data=df )
  fit2.Dq2.factorYr.mmrm.d <- mmrm(Dq2~HeatTrt*WaterTrt+YrFac - YrFac:HeatTrt:WaterTrt  + cs(YrFac|sp), data=df )
  fit2.Dq2.factorYr.mmrm.e <- mmrm(Dq2~HeatTrt+WaterTrt+YrFac - YrFac:HeatTrt:WaterTrt  + cs(YrFac|sp), data=df )
  fit2.Dq2.factorYr.mmrm.f <- mmrm(Dq2~HeatTrt*WaterTrt*YrFac - YrFac:HeatTrt:WaterTrt + cs(YrFac|sp), data=df )
  AIC(fit2.Dq2.factorYr.mmrm.a) # fully-interactive model fits best
  AIC(fit2.Dq2.factorYr.mmrm.b)
  AIC(fit2.Dq2.factorYr.mmrm.c)
  AIC(fit2.Dq2.factorYr.mmrm.d)
  AIC(fit2.Dq2.factorYr.mmrm.e)
  AIC(fit2.Dq2.factorYr.mmrm.f)
  
  Dq2.emm1 <- emmeans(fit2.Dq2.factorYr.mmrm,c('HeatTrt','WaterTrt','YrFac'))
  joint_tests(Dq2.emm1) 
  Dq2.emm2 <- emmeans(fit2.Dq2.factorYr.mmrm,'HeatTrt')
  Dq2.postHoc <- multcomp::cld(object = Dq2.emm2,Letters = letters,alpha = 0.05); Dq2.postHoc 
  
  # get treatment effect:
  Dq2.contrasts <- data.frame(summary(contrast(Dq2.emm1, method = "trt.vs.ctrl", by = c("WaterTrt", "YrFac"))))
  Dq2.contrasts.all <- data.frame(summary(contrast(Dq2.emm1, method = "trt.vs.ctrl", by = c("YrFac"))))
  
  
  # Delta-AIC
  aic_parsimony <- data.frame(
    Model = character(),
    AIC = numeric(),
    Delta_AIC = numeric(),
    stringsAsFactors = FALSE
  )
  
  # List of models
  models_parsimony <- list(
    fit2.Dq2.factorYr.mmrm.a,
    fit2.Dq2.factorYr.mmrm.b,
    fit2.Dq2.factorYr.mmrm.c,
    fit2.Dq2.factorYr.mmrm.d,
    fit2.Dq2.factorYr.mmrm.e,
    fit2.Dq2.factorYr.mmrm.f
  )
  
  # Model names
  model_names_parsimony <- c(
    "Full Interactive",
    "Year additive",
    "Fully additive",
    "Full Interactive - no 3-way",
    "Year additive - no 3-way",
    "Fully additive - no 3-way"
  )
  
  # Calculate AIC and delta-AIC for each model
  aic_values_parsimony <- sapply(models_parsimony, AIC)
  min_aic_parsimony <- min(aic_values_parsimony)
  delta_aic_parsimony <- aic_values_parsimony - min_aic_parsimony
  
  # Populate the dataframe
  aic_parsimony <- data.frame(
    Model = model_names_parsimony,
    AIC = aic_values_parsimony,
    Delta_AIC = delta_aic_parsimony,
    response = "Inverse Simpson diversity",  # Repeat "Dq0" for each row
    selection_stage = "parsimony"
  )
  
  # Sort by AIC
  aic_parsimony <- aic_parsimony %>% arrange(AIC)
  
  # bind into main table
  aic_all <- bind_rows(aic_all, aic_parsimony)
  # Best fit for Dq2 is a CS fully-interactive model (tested above)
  # extract coefficients
  fit_summary.Dq2 <- data.frame(summary(fit2.Dq2.factorYr.mmrm.a)$coefficients) %>% mutate(response = rep("Dq2"))
  
  ################################################################################
  # Berger-Parker dominance (BP)
  ################################################################################
  fit1.BP.factorYr.mmrm <- mmrm(BP~HeatTrt*WaterTrt*YrFac + ar1(YrFac|sp), data=df )
  fit2.BP.factorYr.mmrm <- mmrm(BP~HeatTrt*WaterTrt*YrFac  + cs(YrFac|sp), data=df ) #best fit = lowest AIC
  # fit2.5.BP.factorYr.mmrm <- mmrm(BP~HeatTrt*WaterTrt*YrFac  + csh(YrFac|sp), data=df ) 
  fit3.BP.factorYr.mmrm <- mmrm(BP~HeatTrt*WaterTrt*YrFac  + us(YrFac|sp), data=df )
  fit4.BP.factorYr.mmrm <- mmrm(BP~HeatTrt*WaterTrt*YrFac + ad(YrFac|sp), data=df )
  fit5.BP.factorYr.mmrm <- mmrm(BP~HeatTrt*WaterTrt*YrFac + ar1h(YrFac|sp), data=df )
  AIC(fit1.BP.factorYr.mmrm)
  AIC(fit2.BP.factorYr.mmrm) # best fit
  #AIC(fit2.5.BP.factorYr.mmrm)
  AIC(fit3.BP.factorYr.mmrm)
  AIC(fit4.BP.factorYr.mmrm)
  AIC(fit5.BP.factorYr.mmrm)
  
  # get delta AIC
  aic_covariance <- data.frame(
    Model = character(),
    AIC = numeric(),
    Delta_AIC = numeric(),
    stringsAsFactors = FALSE
  )
  
  # List of models
  models_covariance <- list(
    fit1.BP.factorYr.mmrm,
    fit2.BP.factorYr.mmrm,
    #fit2.5.BP.factorYr.mmrm,
    fit3.BP.factorYr.mmrm,
    fit4.BP.factorYr.mmrm,
    fit5.BP.factorYr.mmrm
  )
  
  # Model names
  model_names_covariance <- c(
    "AR1",
    "CS",
    #"CSH",
    "US",
    "AD",
    "AR1H"
  )
  
  # Calculate AIC and delta-AIC for each model
  aic_values <- sapply(models_covariance, AIC)
  min_aic <- min(aic_values)
  delta_aic <- aic_values - min_aic
  
  # Populate the dataframe
  aic_covariance <- data.frame(
    Model = model_names_covariance,
    AIC = aic_values,
    Delta_AIC = delta_aic,
    response = "Berger-Parker Index",
    selection_stage = "v-cov"
  )
  
  aic_covariance <- aic_covariance %>% arrange(AIC)  # Sort by AIC
  aic_all <- bind_rows(aic_all, aic_covariance)  # bind into main table
  
  BP.emm1 <- emmeans(fit3.BP.factorYr.mmrm,c('HeatTrt','WaterTrt','YrFac'))
  joint_tests(BP.emm1)
  
  BP.emm2 <- emmeans(fit3.BP.factorYr.mmrm,'HeatTrt')
  BP.postHoc <- multcomp::cld(object = BP.emm2,Letters = letters,alpha = 0.05); BP.postHoc
  
  # get treatment effect:
  BP.contrasts <- data.frame(summary(contrast(BP.emm1, method = "trt.vs.ctrl", by = c("WaterTrt", "YrFac"))))
  BP.contrasts.all <- data.frame(summary(contrast(BP.emm1, method = "trt.vs.ctrl", by = c("YrFac"))))
  
  # Get most parsimonious model form
  fit2.BP.factorYr.mmrm.a <- mmrm(BP~HeatTrt*WaterTrt*YrFac + cs(YrFac|sp), data=df )
  fit2.BP.factorYr.mmrm.b <- mmrm(BP~HeatTrt*WaterTrt+YrFac + cs(YrFac|sp), data=df )
  fit2.BP.factorYr.mmrm.c <- mmrm(BP~HeatTrt+WaterTrt+YrFac + cs(YrFac|sp), data=df )
  fit2.BP.factorYr.mmrm.d <- mmrm(BP~HeatTrt*WaterTrt+YrFac - YrFac:HeatTrt:WaterTrt  + cs(YrFac|sp), data=df )
  fit2.BP.factorYr.mmrm.e <- mmrm(BP~HeatTrt+WaterTrt+YrFac - YrFac:HeatTrt:WaterTrt  + cs(YrFac|sp), data=df )
  fit2.BP.factorYr.mmrm.f <- mmrm(BP~HeatTrt*WaterTrt*YrFac - YrFac:HeatTrt:WaterTrt + cs(YrFac|sp), data=df )
  AIC(fit2.BP.factorYr.mmrm.a) # fully-interactive model fits best
  AIC(fit2.BP.factorYr.mmrm.b)
  AIC(fit2.BP.factorYr.mmrm.c)
  AIC(fit2.BP.factorYr.mmrm.d)
  AIC(fit2.BP.factorYr.mmrm.e)
  AIC(fit2.BP.factorYr.mmrm.f)
  
  # Delta-AIC
  aic_parsimony <- data.frame(
    Model = character(),
    AIC = numeric(),
    Delta_AIC = numeric(),
    stringsAsFactors = FALSE
  )
  
  # List of models
  models_parsimony <- list(
    fit2.BP.factorYr.mmrm.a,
    fit2.BP.factorYr.mmrm.b,
    fit2.BP.factorYr.mmrm.c,
    fit2.BP.factorYr.mmrm.d,
    fit2.BP.factorYr.mmrm.e,
    fit2.BP.factorYr.mmrm.f
  )
  
  # Model names
  model_names_parsimony <- c(
    "Full Interactive",
    "Year additive",
    "Fully additive",
    "Full Interactive - no 3-way",
    "Year additive - no 3-way",
    "Fully additive - no 3-way"
  )
  
  # Calculate AIC and delta-AIC for each model
  aic_values_parsimony <- sapply(models_parsimony, AIC)
  min_aic_parsimony <- min(aic_values_parsimony)
  delta_aic_parsimony <- aic_values_parsimony - min_aic_parsimony
  
  # Populate the dataframe
  aic_parsimony <- data.frame(
    Model = model_names_parsimony,
    AIC = aic_values_parsimony,
    Delta_AIC = delta_aic_parsimony,
    response = "Berger-Parker Index",  # Repeat "Dq0" for each row
    selection_stage = "parsimony"
  )
  
  # Sort by AIC
  aic_parsimony <- aic_parsimony %>% arrange(AIC)
  
  # bind into main table
  aic_all <- bind_rows(aic_all, aic_parsimony)
  # Best fit for Berger-Parker dominance is a CS fully-interactive model (tested above)
  # extract coefficients
  fit_summary.BP <- data.frame(summary(fit2.BP.factorYr.mmrm.a)$coefficients) %>% mutate(response = rep("BP"))
  
  ################################################################################
  # Inverse Simpson Evenness (invSimpEven)
  ################################################################################ 
  fit1.invSimpEven.factorYr.mmrm <- mmrm(invSimpEven~HeatTrt*WaterTrt*YrFac + ar1(YrFac|sp), data=df )
  fit2.invSimpEven.factorYr.mmrm <- mmrm(invSimpEven~HeatTrt*WaterTrt*YrFac + cs(YrFac|sp), data=df ) #best fit = lowest AIC
  fit2.5.invSimpEven.factorYr.mmrm <- mmrm(invSimpEven~HeatTrt*WaterTrt*YrFac + csh(YrFac|sp), data=df ) 
  fit3.invSimpEven.factorYr.mmrm <- mmrm(invSimpEven~HeatTrt*WaterTrt*YrFac + us(YrFac|sp), data=df )
  fit4.invSimpEven.factorYr.mmrm <- mmrm(invSimpEven~HeatTrt*WaterTrt*YrFac + ad(YrFac|sp), data=df )
  fit5.invSimpEven.factorYr.mmrm <- mmrm(invSimpEven~HeatTrt*WaterTrt*YrFac + ar1h(YrFac|sp), data=df )
  AIC(fit1.invSimpEven.factorYr.mmrm)
  AIC(fit2.invSimpEven.factorYr.mmrm)
  AIC(fit2.5.invSimpEven.factorYr.mmrm) #  best fit
  AIC(fit3.invSimpEven.factorYr.mmrm)
  AIC(fit4.invSimpEven.factorYr.mmrm)
  AIC(fit5.invSimpEven.factorYr.mmrm)
  # get delta AIC
  aic_covariance <- data.frame(
    Model = character(),
    AIC = numeric(),
    Delta_AIC = numeric(),
    stringsAsFactors = FALSE
  )
  
  # List of models
  models_covariance <- list(
    fit1.invSimpEven.factorYr.mmrm,
    fit2.invSimpEven.factorYr.mmrm,
    fit2.5.invSimpEven.factorYr.mmrm,
    fit3.invSimpEven.factorYr.mmrm,
    fit4.invSimpEven.factorYr.mmrm,
    fit5.invSimpEven.factorYr.mmrm
  )
  
  # Model names
  model_names_covariance <- c(
    "AR1",
    "CS",
    "CSH",
    "US",
    "AD",
    "AR1H"
  )
  
  # Calculate AIC and delta-AIC for each model
  aic_values <- sapply(models_covariance, AIC)
  min_aic <- min(aic_values)
  delta_aic <- aic_values - min_aic
  
  # Populate the dataframe
  aic_covariance <- data.frame(
    Model = model_names_covariance,
    AIC = aic_values,
    Delta_AIC = delta_aic,
    response = "Evenness", 
    selection_stage = "v-cov"
  )
  
  aic_covariance <- aic_covariance %>% arrange(AIC)  # Sort by AIC
  aic_all <- bind_rows(aic_all, aic_covariance)  # bind into main table
  
  invSimpEven.emm1 <- emmeans(fit3.invSimpEven.factorYr.mmrm,c('HeatTrt','WaterTrt','YrFac'))
  joint_tests(invSimpEven.emm1) #HeatTrt:YrFac p=0.0888
  
  invSimpEven.emm2 <- emmeans(fit3.invSimpEven.factorYr.mmrm,'HeatTrt') 
  invSimpEven.postHoc <- multcomp::cld(object = invSimpEven.emm2,Letters = letters,alpha = 0.05);invSimpEven.postHoc # warming alters evenness (reduces it?) in dry years
  # get treatment effect:
  invSimpEven.contrasts <- data.frame(summary(contrast(invSimpEven.emm1, method = "trt.vs.ctrl", by = c("WaterTrt", "YrFac"))))
  invSimpEven.contrasts.all <- data.frame(summary(contrast(invSimpEven.emm1, method = "trt.vs.ctrl", by = c("YrFac"))))
  
  # Get most parsimonious model form
  # most of the following models do not converge 
  fit2.5.invSimpEven.factorYr.mmrm.a <- mmrm(invSimpEven~HeatTrt*WaterTrt*YrFac + csh(YrFac|sp), data=df )
  # fit2.5.invSimpEven.factorYr.mmrm.b <- mmrm(invSimpEven~HeatTrt*WaterTrt+YrFac + csh(YrFac|sp), data=df )
  # fit2.5.invSimpEven.factorYr.mmrm.c <- mmrm(invSimpEven~HeatTrt+WaterTrt+YrFac + csh(YrFac|sp), data=df )
  # fit2.5.invSimpEven.factorYr.mmrm.d <- mmrm(invSimpEven~HeatTrt*WaterTrt+YrFac - YrFac:HeatTrt:WaterTrt  + csh(YrFac|sp), data=df )
  # fit2.5.invSimpEven.factorYr.mmrm.e <- mmrm(invSimpEven~HeatTrt+WaterTrt+YrFac - YrFac:HeatTrt:WaterTrt  + csh(YrFac|sp), data=df )
  # fit2.5.invSimpEven.factorYr.mmrm.f <- mmrm(invSimpEven~HeatTrt*WaterTrt*YrFac - YrFac:HeatTrt:WaterTrt + csh(YrFac|sp), data=df )
  AIC(fit2.5.invSimpEven.factorYr.mmrm.a) # fully-interactive model fits "best" since no other models converge
  # AIC(fit2.5.invSimpEven.factorYr.mmrm.b)
  # AIC(fit2.5.invSimpEven.factorYr.mmrm.c)
  # AIC(fit2.5.invSimpEven.factorYr.mmrm.d)
  # AIC(fit2.5.invSimpEven.factorYr.mmrm.e)
  # AIC(fit2.5.invSimpEven.factorYr.mmrm.f)
  
  # Delta-AIC
  aic_parsimony <- data.frame(
    Model = character(),
    AIC = numeric(),
    Delta_AIC = numeric(),
    stringsAsFactors = FALSE
  )
  
  # List of models
  models_parsimony <- list(
    fit2.5.invSimpEven.factorYr.mmrm.a #,
    #fit2.5.invSimpEven.factorYr.mmrm.b#,
    # fit2.5.invSimpEven.factorYr.mmrm.c,
    #fit2.5.invSimpEven.factorYr.mmrm.d
  )
  
  # Model names
  model_names_parsimony <- c(
    "Full Interactive",
    "Year condition additive" #,
    # "Year & year condition additive",
    # "Fully additive"
  )
  
  # Calculate AIC and delta-AIC for each model
  aic_values_parsimony <- sapply(models_parsimony, AIC)
  min_aic_parsimony <- min(aic_values_parsimony)
  delta_aic_parsimony <- aic_values_parsimony - min_aic_parsimony
  
  # Populate the dataframe
  aic_parsimony <- data.frame(
    Model = model_names_parsimony,
    AIC = aic_values_parsimony,
    Delta_AIC = delta_aic_parsimony,
    response = "Evenness",  # Repeat "Dq0" for each row
    selection_stage = "parsimony"
  )
  
  # Sort by AIC
  aic_parsimony <- aic_parsimony %>% arrange(AIC)
  
  # bind into main table
  aic_all <- bind_rows(aic_all, aic_parsimony)
  # Best fit for Inverse Simpson Evenness is a CSH fully-interactive model (tested above), because other models fail to converge
  # extract coefficients
  fit_summary.invSimpEven <- data.frame(summary(fit2.5.invSimpEven.factorYr.mmrm.a)$coefficients) %>% mutate(response = rep("invSimpEven"))
  
  ################################################################################
  # Andropogon gerardii proportional biomass (p.Ag)
  ################################################################################
  #Andro relative biomass
  fit1.p.Ag.factorYr.mmrm <- mmrm(p.Ag~HeatTrt*WaterTrt*YrFac + ar1(YrFac|sp), data=df )
  fit2.p.Ag.factorYr.mmrm <- mmrm(p.Ag~HeatTrt*WaterTrt*YrFac + cs(YrFac|sp), data=df ) #best fit = lowest AIC
  fit2.5.p.Ag.factorYr.mmrm <- mmrm(p.Ag~HeatTrt*WaterTrt*YrFac + csh(YrFac|sp), data=df ) 
  fit3.p.Ag.factorYr.mmrm <- mmrm(p.Ag~HeatTrt*WaterTrt*YrFac + us(YrFac|sp), data=df )
  fit4.p.Ag.factorYr.mmrm <- mmrm(p.Ag~HeatTrt*WaterTrt*YrFac + ad(YrFac|sp), data=df )
  fit5.p.Ag.factorYr.mmrm <- mmrm(p.Ag~HeatTrt*WaterTrt*YrFac + ar1h(YrFac|sp), data=df )
  AIC(fit1.p.Ag.factorYr.mmrm)
  AIC(fit2.p.Ag.factorYr.mmrm) #best fit = lowest AIC
  AIC(fit2.5.p.Ag.factorYr.mmrm) 
  AIC(fit3.p.Ag.factorYr.mmrm)
  AIC(fit4.p.Ag.factorYr.mmrm)
  AIC(fit5.p.Ag.factorYr.mmrm)
  
  # get delta AIC
  aic_covariance <- data.frame(
    Model = character(),
    AIC = numeric(),
    Delta_AIC = numeric(),
    stringsAsFactors = FALSE
  )
  
  # List of models
  models_covariance <- list(
    fit1.p.Ag.factorYr.mmrm,
    fit2.p.Ag.factorYr.mmrm,
    fit2.5.p.Ag.factorYr.mmrm,
    fit3.p.Ag.factorYr.mmrm,
    fit4.p.Ag.factorYr.mmrm,
    fit5.p.Ag.factorYr.mmrm
  )
  
  # Model names
  model_names_covariance <- c(
    "AR1",
    "CS",
    "CSH",
    "US",
    "AD",
    "AR1H"
  )
  
  # Calculate AIC and delta-AIC for each model
  aic_values <- sapply(models_covariance, AIC)
  min_aic <- min(aic_values)
  delta_aic <- aic_values - min_aic
  
  # Populate the dataframe
  aic_covariance <- data.frame(
    Model = model_names_covariance,
    AIC = aic_values,
    Delta_AIC = delta_aic,
    response = "A. gerardii proportion", 
    selection_stage = "v-cov"
  )
  
  aic_covariance <- aic_covariance %>% arrange(AIC)  # Sort by AIC
  aic_all <- bind_rows(aic_all, aic_covariance)  # bind into main table
  
  p.Ag.emm1 <- emmeans(fit3.p.Ag.factorYr.mmrm,c('HeatTrt','WaterTrt','YrFac'))
  joint_tests(p.Ag.emm1) 
  
  p.Ag.emm2 <- emmeans(fit3.p.Ag.factorYr.mmrm,'HeatTrt') 
  multcomp::cld(object = p.Ag.emm2,Letters = letters,alpha = 0.05)
  
 # get treatment effect:
  p.Ag.contrasts <- data.frame(summary(contrast(p.Ag.emm1, method = "trt.vs.ctrl", by = c("WaterTrt", "YrFac"))))
  p.Ag.contrasts.all <- data.frame(summary(contrast(p.Ag.emm1, method = "trt.vs.ctrl", by = c("YrFac"))))
  
  
  # MAGGIE ADDED: get most parsimonious model form
  fit2.p.Ag.factorYr.mmrm.a <- mmrm(p.Ag~HeatTrt*WaterTrt*YrFac + cs(YrFac|sp), data=df )
  fit2.p.Ag.factorYr.mmrm.b <- mmrm(p.Ag~HeatTrt*WaterTrt+YrFac + cs(YrFac|sp), data=df )
  fit2.p.Ag.factorYr.mmrm.c <- mmrm(p.Ag~HeatTrt+WaterTrt+YrFac + cs(YrFac|sp), data=df )
  fit2.p.Ag.factorYr.mmrm.d <- mmrm(p.Ag~HeatTrt*WaterTrt+YrFac - YrFac:HeatTrt:WaterTrt  + cs(YrFac|sp), data=df )
  fit2.p.Ag.factorYr.mmrm.e <- mmrm(p.Ag~HeatTrt+WaterTrt+YrFac - YrFac:HeatTrt:WaterTrt  + cs(YrFac|sp), data=df )
  fit2.p.Ag.factorYr.mmrm.f <- mmrm(p.Ag~HeatTrt*WaterTrt*YrFac - YrFac:HeatTrt:WaterTrt + cs(YrFac|sp), data = df )
  AIC(fit2.p.Ag.factorYr.mmrm.a) 
  AIC(fit2.p.Ag.factorYr.mmrm.b)
  AIC(fit2.p.Ag.factorYr.mmrm.c) # fully-ADDITIVE model fits best? WTF?
  AIC(fit2.p.Ag.factorYr.mmrm.d)
  AIC(fit2.p.Ag.factorYr.mmrm.e)
  AIC(fit2.p.Ag.factorYr.mmrm.f)
  
  # Delta-AIC
  aic_parsimony <- data.frame(
    Model = character(),
    AIC = numeric(),
    Delta_AIC = numeric(),
    stringsAsFactors = FALSE
  )
  
  # List of models
  models_parsimony <- list(
    fit2.p.Ag.factorYr.mmrm.a,
    fit2.p.Ag.factorYr.mmrm.b,
    fit2.p.Ag.factorYr.mmrm.c,
    fit2.p.Ag.factorYr.mmrm.d,
    fit2.p.Ag.factorYr.mmrm.e,
    fit2.p.Ag.factorYr.mmrm.f
  )
  
  # Model names
  model_names_parsimony <- c(
    "Full Interactive",
    "Year additive",
    "Fully additive",
    "Full Interactive - no 3-way",
    "Year additive - no 3-way",
    "Fully additive - no 3-way"
  )
  
  # Calculate AIC and delta-AIC for each model
  aic_values_parsimony <- sapply(models_parsimony, AIC)
  min_aic_parsimony <- min(aic_values_parsimony)
  delta_aic_parsimony <- aic_values_parsimony - min_aic_parsimony
  
  # Populate the dataframe
  aic_parsimony <- data.frame(
    Model = model_names_parsimony,
    AIC = aic_values_parsimony,
    Delta_AIC = delta_aic_parsimony,
    response = "A. gerardii proportion",  # Repeat "Dq0" for each row
    selection_stage = "parsimony"
  )
  
  # Sort by AIC
  aic_parsimony <- aic_parsimony %>% arrange(AIC)
  
  # bind into main table
  aic_all <- bind_rows(aic_all, aic_parsimony)
  # Best fit for proportion of Andropogon gerardii is a CS ADDITIVE model
  
  # now test new best-fit CS additive model
  p.Ag.emm1c <- emmeans(fit2.p.Ag.factorYr.mmrm.c,c('HeatTrt','WaterTrt','YrFac')) # for models
  p.Ag.emm1 <- emmeans(fit2.p.Ag.factorYr.mmrm,c('HeatTrt','WaterTrt','YrFac')) # for statistical output
  joint_tests(p.Ag.emm1c) #HeatTrt: p=0.0862, YrFac p=0.0032
  #p.Ag.emm2c <- emmeans(fit2.p.Ag.factorYr.mmrm.c) 
  #p.Ag.postHoc <- multcomp::cld(object = p.Ag.emm2c,Letters = letters,alpha = 0.05); p.Ag.postHoc # A. gerardii biomass is significantly different in dry years (decreases?)
  # extract coefficients
  fit_summary.p.Ag <- data.frame(summary(fit2.p.Ag.factorYr.mmrm.c)$coefficients) %>% mutate(response = rep("Ag"))
  
}
 
#### Models for functional group proportion responses      
{        
  ################################################################################
  # Legume biomass proportion
  ################################################################################
  fit1.LegumeBiomass.factorYr.mmrm <- mmrm(LegumeBiomass~HeatTrt*WaterTrt*YrFac + ar1(YrFac|sp), data=df )
  fit2.LegumeBiomass.factorYr.mmrm <- mmrm(LegumeBiomass~HeatTrt*WaterTrt*YrFac + cs(YrFac|sp), data=df ) #best fit = lowest AIC
  fit2.5.LegumeBiomass.factorYr.mmrm <- mmrm(LegumeBiomass~HeatTrt*WaterTrt*YrFac + csh(YrFac|sp), data=df ) 
  fit3.LegumeBiomass.factorYr.mmrm <- mmrm(LegumeBiomass~HeatTrt*WaterTrt*YrFac + us(YrFac|sp), data=df )
  fit4.LegumeBiomass.factorYr.mmrm <- mmrm(LegumeBiomass~HeatTrt*WaterTrt*YrFac + ad(YrFac|sp), data=df )
  fit5.LegumeBiomass.factorYr.mmrm <- mmrm(LegumeBiomass~HeatTrt*WaterTrt*YrFac + ar1h(YrFac|sp), data=df )
  AIC(fit1.LegumeBiomass.factorYr.mmrm)
  AIC(fit2.LegumeBiomass.factorYr.mmrm) 
  AIC(fit2.5.LegumeBiomass.factorYr.mmrm) 
  AIC(fit3.LegumeBiomass.factorYr.mmrm)
  AIC(fit4.LegumeBiomass.factorYr.mmrm)
  AIC(fit5.LegumeBiomass.factorYr.mmrm)
  # get delta AIC
  aic_covariance <- data.frame(
    Model = character(),
    AIC = numeric(),
    Delta_AIC = numeric(),
    stringsAsFactors = FALSE
  )
  
  # List of models
  models_covariance <- list(
    fit1.LegumeBiomass.factorYr.mmrm,
    fit2.LegumeBiomass.factorYr.mmrm,
    fit2.5.LegumeBiomass.factorYr.mmrm,
    fit3.LegumeBiomass.factorYr.mmrm,
    fit4.LegumeBiomass.factorYr.mmrm,
    fit5.LegumeBiomass.factorYr.mmrm
  )
  
  # Model names
  model_names_covariance <- c(
    "AR1",
    "CS",
    "CSH",
    "US",
    "AD",
    "AR1H"
  )
  
  # Calculate AIC and delta-AIC for each model
  aic_values <- sapply(models_covariance, AIC)
  min_aic <- min(aic_values)
  delta_aic <- aic_values - min_aic
  
  # Populate the dataframe
  aic_covariance <- data.frame(
    Model = model_names_covariance,
    AIC = aic_values,
    Delta_AIC = delta_aic,
    response = "Legume biomass prop.", 
    selection_stage = "v-cov"
  )
  
  aic_covariance <- aic_covariance %>% arrange(AIC)  # Sort by AIC
  aic_all <- bind_rows(aic_all, aic_covariance)  # bind into main table
  
  LegumeBiomass.emm1 <- emmeans(fit4.LegumeBiomass.factorYr.mmrm,c('HeatTrt','WaterTrt','YrFac'))
  joint_tests(LegumeBiomass.emm1) 
  LegumeBiomass.emm2 <- emmeans(fit4.LegumeBiomass.factorYr.mmrm,'HeatTrt') 
  LegumeBiomass.postHoc <- multcomp::cld(object = LegumeBiomass.emm2,Letters = letters,alpha = 0.05); LegumeBiomass.postHoc
  # get treatment effect:
  LegumeBiomass.contrasts <- data.frame(summary(contrast(LegumeBiomass.emm1, method = "trt.vs.ctrl", by = c("WaterTrt", "YrFac"))))
  LegumeBiomass.contrasts.all <- data.frame(summary(contrast(LegumeBiomass.emm1, method = "trt.vs.ctrl", by = c("YrFac"))))
  
  # get most parsimonious model form
  fit4.LegumeBiomass.factorYr.mmrm.a <- mmrm(LegumeBiomass~HeatTrt*WaterTrt*YrFac + ad(YrFac|sp), data=df )
  fit4.LegumeBiomass.factorYr.mmrm.b <- mmrm(LegumeBiomass~HeatTrt*WaterTrt+YrFac + ad(YrFac|sp), data=df )
  fit4.LegumeBiomass.factorYr.mmrm.c <- mmrm(LegumeBiomass~HeatTrt+WaterTrt+YrFac + ad(YrFac|sp), data=df )
  fit4.LegumeBiomass.factorYr.mmrm.d <- mmrm(LegumeBiomass~HeatTrt*WaterTrt+YrFac - YrFac:HeatTrt:WaterTrt  + ad(YrFac|sp), data=df )
  fit4.LegumeBiomass.factorYr.mmrm.e <- mmrm(LegumeBiomass~HeatTrt+WaterTrt+YrFac - YrFac:HeatTrt:WaterTrt  + ad(YrFac|sp), data=df )
  fit4.LegumeBiomass.factorYr.mmrm.f <- mmrm(LegumeBiomass~HeatTrt*WaterTrt*YrFac - YrFac:HeatTrt:WaterTrt + ad(YrFac|sp), data=df )
  AIC(fit4.LegumeBiomass.factorYr.mmrm.a) 
  AIC(fit4.LegumeBiomass.factorYr.mmrm.b)
  AIC(fit4.LegumeBiomass.factorYr.mmrm.c) # fully ADDITIVE fits best
  AIC(fit4.LegumeBiomass.factorYr.mmrm.d)
  AIC(fit4.LegumeBiomass.factorYr.mmrm.e)
  AIC(fit4.LegumeBiomass.factorYr.mmrm.f)
  
  # Best fit for proportion of Legume proportional biomass is a CS additive model
  # Delta-AIC
  aic_parsimony <- data.frame(
    Model = character(),
    AIC = numeric(),
    Delta_AIC = numeric(),
    stringsAsFactors = FALSE
  )
  
  # List of models
  models_parsimony <- list(
    fit4.LegumeBiomass.factorYr.mmrm.a,
    fit4.LegumeBiomass.factorYr.mmrm.b,
    fit4.LegumeBiomass.factorYr.mmrm.c,
    fit4.LegumeBiomass.factorYr.mmrm.d,
    fit4.LegumeBiomass.factorYr.mmrm.e,
    fit4.LegumeBiomass.factorYr.mmrm.f
  )
  
  # Model names
  model_names_parsimony <- c(
    "Full Interactive",
    "Year additive",
    "Fully additive",
    "Full Interactive - no 3-way",
    "Year additive - no 3-way",
    "Fully additive - no 3-way"
  )
  
  # Calculate AIC and delta-AIC for each model
  aic_values_parsimony <- sapply(models_parsimony, AIC)
  min_aic_parsimony <- min(aic_values_parsimony)
  delta_aic_parsimony <- aic_values_parsimony - min_aic_parsimony
  
  # Populate the dataframe
  aic_parsimony <- data.frame(
    Model = model_names_parsimony,
    AIC = aic_values_parsimony,
    Delta_AIC = delta_aic_parsimony,
    response = "Legume biomass prop.",  # Repeat for each row
    selection_stage = "parsimony"
  )
  
  # Sort by AIC
  aic_parsimony <- aic_parsimony %>% arrange(AIC)
  
  # bind into main table
  aic_all <- bind_rows(aic_all, aic_parsimony)
  
  # now test new best-fit CS additive model
  LegumeBiomass.emm1c <- emmeans(fit4.LegumeBiomass.factorYr.mmrm.c,c('HeatTrt','WaterTrt','YrFac')) 
  LegumeBiomass.emm1 <- emmeans(fit4.LegumeBiomass.factorYr.mmrm.c,c('HeatTrt','WaterTrt','YrFac')) 
  joint_tests(LegumeBiomass.emm1c)  
  LegumeBiomass.emm2c <- emmeans(fit4.LegumeBiomass.factorYr.mmrm.a,'HeatTrt') 
  LegumeBiomass.postHoc <- multcomp::cld(object = LegumeBiomass.emm2c,Letters = letters,alpha = 0.05); LegumeBiomass.postHoc 
  # extract coefficients
  fit_summary.LegumeBiomass <- data.frame(summary(fit2.LegumeBiomass.factorYr.mmrm)$coefficients) %>% mutate(response = rep("LegumeBiomass"))
  
  
  ################################################################################
  # Forb proportional biomass
  ################################################################################
  fit1.ForbBiomass.factorYr.mmrm <- mmrm(ForbBiomass~HeatTrt*WaterTrt*YrFac + ar1(YrFac|sp), data=df )
  fit2.ForbBiomass.factorYr.mmrm <- mmrm(ForbBiomass~HeatTrt*WaterTrt*YrFac + cs(YrFac|sp), data=df ) #best fit = lowest AIC
  #fit2.5.ForbBiomass.factorYr.mmrm <- mmrm(ForbBiomass~HeatTrt*WaterTrt*YrFac + csh(YrFac|sp), data=df ) 
  fit3.ForbBiomass.factorYr.mmrm <- mmrm(ForbBiomass~HeatTrt*WaterTrt*YrFac + us(YrFac|sp), data=df )
  fit4.ForbBiomass.factorYr.mmrm <- mmrm(ForbBiomass~HeatTrt*WaterTrt*YrFac + ad(YrFac|sp), data=df )
  fit5.ForbBiomass.factorYr.mmrm <- mmrm(ForbBiomass~HeatTrt*WaterTrt*YrFac + ar1h(YrFac|sp), data=df )
  AIC(fit1.ForbBiomass.factorYr.mmrm)
  AIC(fit2.ForbBiomass.factorYr.mmrm)
  #AIC(fit2.5.ForbBiomass.factorYr.mmrm) 
  AIC(fit3.ForbBiomass.factorYr.mmrm)  # best fit
  AIC(fit4.ForbBiomass.factorYr.mmrm)
  AIC(fit5.ForbBiomass.factorYr.mmrm)
  
  # get delta AIC
  aic_covariance <- data.frame(
    Model = character(),
    AIC = numeric(),
    Delta_AIC = numeric(),
    stringsAsFactors = FALSE
  )
  
  # List of models
  models_covariance <- list(
    fit1.ForbBiomass.factorYr.mmrm,
    fit2.ForbBiomass.factorYr.mmrm,
    #fit2.5.ForbBiomass.factorYr.mmrm,
    fit3.ForbBiomass.factorYr.mmrm,
    fit4.ForbBiomass.factorYr.mmrm,
    fit5.ForbBiomass.factorYr.mmrm
  )
  
  # Model names
  model_names_covariance <- c(
    "AR1",
    "CS",
    #"CSH",
    "US",
    "AD",
    "AR1H"
  )
  
  # Calculate AIC and delta-AIC for each model
  aic_values <- sapply(models_covariance, AIC)
  min_aic <- min(aic_values)
  delta_aic <- aic_values - min_aic
  
  # Populate the dataframe
  aic_covariance <- data.frame(
    Model = model_names_covariance,
    AIC = aic_values,
    Delta_AIC = delta_aic,
    response = "Forb biomass prop.", 
    selection_stage = "v-cov"
  )
  
  aic_covariance <- aic_covariance %>% arrange(AIC)  # Sort by AIC
  aic_all <- bind_rows(aic_all, aic_covariance)  # bind into main table
  
  ForbBiomass.emm1 <- emmeans(fit4.ForbBiomass.factorYr.mmrm,c('HeatTrt','WaterTrt','YrFac'))
  joint_tests(ForbBiomass.emm1) 
  # get treatment effect:
  ForbBiomass.contrasts <- data.frame(summary(contrast(ForbBiomass.emm1, method = "trt.vs.ctrl", by = c("WaterTrt", "YrFac"))))
  ForbBiomass.contrasts.all <- data.frame(summary(contrast(ForbBiomass.emm1, method = "trt.vs.ctrl", by = c("YrFac"))))
  
  # get most parsimonious model form
  fit4.ForbBiomass.factorYr.mmrm.a <- mmrm(ForbBiomass~HeatTrt*WaterTrt*YrFac + ad(YrFac|sp), data=df )
  fit4.ForbBiomass.factorYr.mmrm.b <- mmrm(ForbBiomass~HeatTrt*WaterTrt+YrFac + ad(YrFac|sp), data=df )
  fit4.ForbBiomass.factorYr.mmrm.c <- mmrm(ForbBiomass~HeatTrt+WaterTrt+YrFac + ad(YrFac|sp), data=df )
  fit4.ForbBiomass.factorYr.mmrm.d <- mmrm(ForbBiomass~HeatTrt*WaterTrt+YrFac - YrFac:HeatTrt:WaterTrt  + ad(YrFac|sp), data=df )
  fit4.ForbBiomass.factorYr.mmrm.e <- mmrm(ForbBiomass~HeatTrt+WaterTrt+YrFac - YrFac:HeatTrt:WaterTrt  + ad(YrFac|sp), data=df )
  fit4.ForbBiomass.factorYr.mmrm.f <- mmrm(ForbBiomass~HeatTrt*WaterTrt*YrFac - YrFac:HeatTrt:WaterTrt + ad(YrFac|sp), data=df )
  AIC(fit4.ForbBiomass.factorYr.mmrm.a) 
  AIC(fit4.ForbBiomass.factorYr.mmrm.b)
  AIC(fit4.ForbBiomass.factorYr.mmrm.c) # fully-additive model best fit again
  AIC(fit4.ForbBiomass.factorYr.mmrm.d)
  AIC(fit4.ForbBiomass.factorYr.mmrm.e)
  AIC(fit4.ForbBiomass.factorYr.mmrm.f)
  
  # Best fit for proportion of Forb proportional biomass is a US ADDITIVE model
  # Delta-AIC
  aic_parsimony <- data.frame(
    Model = character(),
    AIC = numeric(),
    Delta_AIC = numeric(),
    stringsAsFactors = FALSE
  )
  
  # List of models
  models_parsimony <- list(
    fit4.ForbBiomass.factorYr.mmrm.a,
    fit4.ForbBiomass.factorYr.mmrm.b,
    fit4.ForbBiomass.factorYr.mmrm.c,
    fit4.ForbBiomass.factorYr.mmrm.d,
    fit4.ForbBiomass.factorYr.mmrm.e,
    fit4.ForbBiomass.factorYr.mmrm.f
  )
  
  # Model names
  model_names_parsimony <- c(
    "Full Interactive",
    "Year additive",
    "Fully additive",
    "Full Interactive - no 3-way",
    "Year additive - no 3-way",
    "Fully additive - no 3-way"
    
  )
  
  # Calculate AIC and delta-AIC for each model
  aic_values_parsimony <- sapply(models_parsimony, AIC)
  min_aic_parsimony <- min(aic_values_parsimony)
  delta_aic_parsimony <- aic_values_parsimony - min_aic_parsimony
  
  # Populate the dataframe
  aic_parsimony <- data.frame(
    Model = model_names_parsimony,
    AIC = aic_values_parsimony,
    Delta_AIC = delta_aic_parsimony,
    response = "Forb biomass prop.",  # Repeat "Dq0" for each row
    selection_stage = "parsimony"
  )
  
  # Sort by AIC
  aic_parsimony <- aic_parsimony %>% arrange(AIC)
  
  # bind into main table
  aic_all <- bind_rows(aic_all, aic_parsimony)
  
  # now test new best-fit US additive model
  ForbBiomass.emm1c <- emmeans(fit4.ForbBiomass.factorYr.mmrm.c,c('HeatTrt','WaterTrt','YrFac'))
  joint_tests(ForbBiomass.emm1c) 
  # extract coefficients
  fit_summary.ForbBiomass <- data.frame(summary(fit4.ForbBiomass.factorYr.mmrm.c)$coefficients) %>% mutate(response = rep("ForbBiomass"))
  ForbBiomass.emm1 <- emmeans(fit4.ForbBiomass.factorYr.mmrm.c,c('HeatTrt','WaterTrt','YrFac'))
  
  ################################################################################
  # C3 grass biomass proportion
  ################################################################################
  fit1.C3Biomass.factorYr.mmrm <- mmrm(C3Biomass~HeatTrt*WaterTrt*YrFac + ar1(YrFac|sp), data=df )
  fit2.C3Biomass.factorYr.mmrm <- mmrm(C3Biomass~HeatTrt*WaterTrt*YrFac + cs(YrFac|sp), data=df ) #best fit = lowest AIC
  fit2.5.C3Biomass.factorYr.mmrm <- mmrm(C3Biomass~HeatTrt*WaterTrt*YrFac + csh(YrFac|sp), data=df )
  fit3.C3Biomass.factorYr.mmrm <- mmrm(C3Biomass~HeatTrt*WaterTrt*YrFac + us(YrFac|sp), data=df )
  fit4.C3Biomass.factorYr.mmrm <- mmrm(C3Biomass~HeatTrt*WaterTrt*YrFac + ad(YrFac|sp), data=df )
  fit5.C3Biomass.factorYr.mmrm <- mmrm(C3Biomass~HeatTrt*WaterTrt*YrFac + ar1h(YrFac|sp), data=df )
  AIC(fit1.C3Biomass.factorYr.mmrm)
  AIC(fit2.C3Biomass.factorYr.mmrm)
  AIC(fit2.5.C3Biomass.factorYr.mmrm) 
  AIC(fit3.C3Biomass.factorYr.mmrm)
  AIC(fit4.C3Biomass.factorYr.mmrm)
  AIC(fit5.C3Biomass.factorYr.mmrm) #best fit = lowest AIC
  # get delta AIC
  aic_covariance <- data.frame(
    Model = character(),
    AIC = numeric(),
    Delta_AIC = numeric(),
    stringsAsFactors = FALSE
  )
  
  # List of models
  models_covariance <- list(
    fit1.C3Biomass.factorYr.mmrm,
    fit2.C3Biomass.factorYr.mmrm,
    fit2.5.C3Biomass.factorYr.mmrm,
    fit3.C3Biomass.factorYr.mmrm,
    fit4.C3Biomass.factorYr.mmrm,
    fit5.C3Biomass.factorYr.mmrm
  )
  
  # Model names
  model_names_covariance <- c(
    "AR1",
    "CS",
    "CSH",
    "US",
    "AD",
    "AR1H"
  )
  
  # Calculate AIC and delta-AIC for each model
  aic_values <- sapply(models_covariance, AIC)
  min_aic <- min(aic_values)
  delta_aic <- aic_values - min_aic
  
  # Populate the dataframe
  aic_covariance <- data.frame(
    Model = model_names_covariance,
    AIC = aic_values,
    Delta_AIC = delta_aic,
    response = "C3 grass biomass prop.", 
    selection_stage = "v-cov"
  )
  
  aic_covariance <- aic_covariance %>% arrange(AIC)  # Sort by AIC
  aic_all <- bind_rows(aic_all, aic_covariance)  # bind into main table
  
  C3Biomass.emm1 <- emmeans(fit5.C3Biomass.factorYr.mmrm,c('HeatTrt','WaterTrt','YrFac'))
  joint_tests(C3Biomass.emm1) 
  C3Biomass.emm2 <- emmeans(fit5.C3Biomass.factorYr.mmrm,c('YrFac')) 
  multcomp::cld(object = C3Biomass.emm2,Letters = letters,alpha = 0.05) 
  # get treatment effect:
  C3Biomass.contrasts <- data.frame(summary(contrast(C3Biomass.emm1, method = "trt.vs.ctrl", by = c("WaterTrt", "YrFac"))))
  C3Biomass.contrasts.all <- data.frame(summary(contrast(C3Biomass.emm1, method = "trt.vs.ctrl", by = c("YrFac"))))
  
  # Get most parsimonious model form
  fit5.C3Biomass.factorYr.mmrm.a <- mmrm(C3Biomass~HeatTrt*WaterTrt*YrFac + ar1h(YrFac|sp), data=df )
  fit5.C3Biomass.factorYr.mmrm.b <- mmrm(C3Biomass~HeatTrt*WaterTrt+YrFac + ar1h(YrFac|sp), data=df )
  fit5.C3Biomass.factorYr.mmrm.c <- mmrm(C3Biomass~HeatTrt+WaterTrt+YrFac + ar1h(YrFac|sp), data=df )
  fit5.C3Biomass.factorYr.mmrm.d <- mmrm(C3Biomass~HeatTrt*WaterTrt+YrFac - YrFac:HeatTrt:WaterTrt  + ar1h(YrFac|sp), data=df )
  fit5.C3Biomass.factorYr.mmrm.e <- mmrm(C3Biomass~HeatTrt+WaterTrt+YrFac - YrFac:HeatTrt:WaterTrt  + ar1h(YrFac|sp), data=df )
  fit5.C3Biomass.factorYr.mmrm.f <- mmrm(C3Biomass~HeatTrt*WaterTrt*YrFac - YrFac:HeatTrt:WaterTrt + ar1h(YrFac|sp), data=df )
  AIC(fit5.C3Biomass.factorYr.mmrm.a) 
  AIC(fit5.C3Biomass.factorYr.mmrm.b)
  AIC(fit5.C3Biomass.factorYr.mmrm.c) # additive model fits best 
  AIC(fit5.C3Biomass.factorYr.mmrm.d)
  AIC(fit5.C3Biomass.factorYr.mmrm.e)
  AIC(fit5.C3Biomass.factorYr.mmrm.f)
  
  # MAGGIE ADDED: now test new best-fit AR1H additive model
  C3Biomass.emm1c <- emmeans(fit5.C3Biomass.factorYr.mmrm.c,c('HeatTrt','WaterTrt','YrFac'))
  joint_tests(C3Biomass.emm1c) 
  C3Biomass.emm2d <- emmeans(fit5.C3Biomass.factorYr.mmrm.d,c('WaterTrt')) 
  multcomp::cld(object = C3Biomass.emm2d,Letters = letters,alpha = 0.05) 
  # extract coefficients
  fit_summary.C3Biomass <- data.frame(summary(fit5.C3Biomass.factorYr.mmrm.d)$coefficients) %>% mutate(response = rep("C3Biomass"))
  
  # Delta-AIC
  aic_parsimony <- data.frame(
    Model = character(),
    AIC = numeric(),
    Delta_AIC = numeric(),
    stringsAsFactors = FALSE
  )
  
  # List of models
  models_parsimony <- list(
    fit5.C3Biomass.factorYr.mmrm.a,
    fit5.C3Biomass.factorYr.mmrm.b,
    fit5.C3Biomass.factorYr.mmrm.c,
    fit5.C3Biomass.factorYr.mmrm.d,
    fit5.C3Biomass.factorYr.mmrm.e,
    fit5.C3Biomass.factorYr.mmrm.f
  )
  
  # Model names
  model_names_parsimony <- c(
    "Full Interactive",
    "Year additive",
    "Fully additive",
    "Full Interactive - no 3-way",
    "Year additive - no 3-way",
    "Fully additive - no 3-way"
  )
  
  # Calculate AIC and delta-AIC for each model
  aic_values_parsimony <- sapply(models_parsimony, AIC)
  min_aic_parsimony <- min(aic_values_parsimony)
  delta_aic_parsimony <- aic_values_parsimony - min_aic_parsimony
  
  # Populate the dataframe
  aic_parsimony <- data.frame(
    Model = model_names_parsimony,
    AIC = aic_values_parsimony,
    Delta_AIC = delta_aic_parsimony,
    response = "C3 grass biomass prop.",  # Repeat "Dq0" for each row
    selection_stage = "parsimony"
  )
  
  # Sort by AIC
  aic_parsimony <- aic_parsimony %>% arrange(AIC)
  
  # bind into main table
  aic_all <- bind_rows(aic_all, aic_parsimony)
  
  ################################################################################
  # C4 grass biomass
  ################################################################################
  fit1.C4Biomass.factorYr.mmrm <- mmrm(C4Biomass~HeatTrt*WaterTrt*YrFac + ar1(YrFac|sp), data=df )
  fit2.C4Biomass.factorYr.mmrm <- mmrm(C4Biomass~HeatTrt*WaterTrt*YrFac + cs(YrFac|sp), data=df ) #best fit = lowest AIC
  fit2.5.C4Biomass.factorYr.mmrm <- mmrm(C4Biomass~HeatTrt*WaterTrt*YrFac + csh(YrFac|sp), data=df ) 
  fit3.C4Biomass.factorYr.mmrm <- mmrm(C4Biomass~HeatTrt*WaterTrt*YrFac + us(YrFac|sp), data=df )
  fit4.C4Biomass.factorYr.mmrm <- mmrm(C4Biomass~HeatTrt*WaterTrt*YrFac + ad(YrFac|sp), data=df )
  fit5.C4Biomass.factorYr.mmrm <- mmrm(C4Biomass~HeatTrt*WaterTrt*YrFac + ar1h(YrFac|sp), data=df )
  AIC(fit1.C4Biomass.factorYr.mmrm)
  AIC(fit2.C4Biomass.factorYr.mmrm) 
  AIC(fit2.5.C4Biomass.factorYr.mmrm) 
  AIC(fit3.C4Biomass.factorYr.mmrm)#best fit = lowest AIC
  AIC(fit4.C4Biomass.factorYr.mmrm)
  AIC(fit5.C4Biomass.factorYr.mmrm)
  
  # get delta AIC
  aic_covariance <- data.frame(
    Model = character(),
    AIC = numeric(),
    Delta_AIC = numeric(),
    stringsAsFactors = FALSE
  )
  
  # List of models
  models_covariance <- list(
    fit1.C4Biomass.factorYr.mmrm,
    fit2.C4Biomass.factorYr.mmrm,
    fit2.5.C4Biomass.factorYr.mmrm,
    fit3.C4Biomass.factorYr.mmrm,
    fit4.C4Biomass.factorYr.mmrm,
    fit5.C4Biomass.factorYr.mmrm
  )
  
  # Model names
  model_names_covariance <- c(
    "AR1",
    "CS",
    "CSH",
    "US",
    "AD",
    "AR1H"
  )
  
  # Calculate AIC and delta-AIC for each model
  aic_values <- sapply(models_covariance, AIC)
  min_aic <- min(aic_values)
  delta_aic <- aic_values - min_aic
  
  # Populate the dataframe
  aic_covariance <- data.frame(
    Model = model_names_covariance,
    AIC = aic_values,
    Delta_AIC = delta_aic,
    response = "C4 grass biomass prop.", 
    selection_stage = "v-cov"
  )
  
  aic_covariance <- aic_covariance %>% arrange(AIC)  # Sort by AIC
  aic_all <- bind_rows(aic_all, aic_covariance)  # bind into main table
  
  
  C4Biomass.emm1 <- emmeans(fit3.C4Biomass.factorYr.mmrm,c('HeatTrt','WaterTrt','YrFac'))
  joint_tests(C4Biomass.emm1) 
  # get treatment effect:
  C4Biomass.contrasts <- data.frame(summary(contrast(C4Biomass.emm1, method = "trt.vs.ctrl", by = c("WaterTrt", "YrFac"))))
  C4Biomass.contrasts.all <- data.frame(summary(contrast(C4Biomass.emm1, method = "trt.vs.ctrl", by = c("YrFac"))))
  
  # get most parsimonious model form
  fit3.C4Biomass.factorYr.mmrm.a <- mmrm(C4Biomass~HeatTrt*WaterTrt*YrFac + us(YrFac|sp), data=df )
  fit3.C4Biomass.factorYr.mmrm.b <- mmrm(C4Biomass~HeatTrt*WaterTrt+YrFac + us(YrFac|sp), data=df )
  fit3.C4Biomass.factorYr.mmrm.c <- mmrm(C4Biomass~HeatTrt+WaterTrt+YrFac + us(YrFac|sp), data=df )
  fit3.C4Biomass.factorYr.mmrm.d <- mmrm(C4Biomass~HeatTrt*WaterTrt+YrFac - YrFac:HeatTrt:WaterTrt + us(YrFac|sp), data=df )
  fit3.C4Biomass.factorYr.mmrm.e <- mmrm(C4Biomass~HeatTrt+WaterTrt+YrFac - YrFac:HeatTrt:WaterTrt + us(YrFac|sp), data=df )
  fit3.C4Biomass.factorYr.mmrm.f <- mmrm(C4Biomass~HeatTrt*WaterTrt*YrFac - YrFac:HeatTrt:WaterTrt + us(YrFac|sp), data=df )
  AIC(fit3.C4Biomass.factorYr.mmrm.a) 
  AIC(fit3.C4Biomass.factorYr.mmrm.b)
  AIC(fit3.C4Biomass.factorYr.mmrm.c) # additive model fits best
  AIC(fit3.C4Biomass.factorYr.mmrm.d)
  AIC(fit3.C4Biomass.factorYr.mmrm.e)
  AIC(fit3.C4Biomass.factorYr.mmrm.f)
  
  # Best fit for proportion of C4 proportional biomass is a CS ADDITIVE model
  
  # now test new best-fit CS additive model
  C4Biomass.emm1c <- emmeans(fit3.C4Biomass.factorYr.mmrm.c,c('HeatTrt','WaterTrt','YrFac')) 
  joint_tests(C4Biomass.emm1c) 

  # Delta-AIC
  aic_parsimony <- data.frame(
    Model = character(),
    AIC = numeric(),
    Delta_AIC = numeric(),
    stringsAsFactors = FALSE
  )
  
  # List of models
  models_parsimony <- list(
    fit3.C4Biomass.factorYr.mmrm.a,
    fit3.C4Biomass.factorYr.mmrm.b,
    fit3.C4Biomass.factorYr.mmrm.c,
    fit3.C4Biomass.factorYr.mmrm.d,
    fit3.C4Biomass.factorYr.mmrm.e,
    fit3.C4Biomass.factorYr.mmrm.f
  )
  
  # Model names
  model_names_parsimony <- c(
    "Full Interactive",
    "Year additive",
    "Fully additive",
    "Full Interactive - no 3-way",
    "Year additive - no 3-way",
    "Fully additive - no 3-way"
  )
  
  # Calculate AIC and delta-AIC for each model
  aic_values_parsimony <- sapply(models_parsimony, AIC)
  min_aic_parsimony <- min(aic_values_parsimony)
  delta_aic_parsimony <- aic_values_parsimony - min_aic_parsimony
  
  # Populate the dataframe
  aic_parsimony <- data.frame(
    Model = model_names_parsimony,
    AIC = aic_values_parsimony,
    Delta_AIC = delta_aic_parsimony,
    response = "C4 grass biomass prop.",  # Repeat "Dq0" for each row
    selection_stage = "parsimony"
  )
  
  # Sort by AIC
  aic_parsimony <- aic_parsimony %>% arrange(AIC)
  
  # bind into main table
  aic_all <- bind_rows(aic_all, aic_parsimony)

################################################################################
  
  # Assemble dataframe for p-values and other stats
  mmrm.stats <- bind_rows(as.data.frame(joint_tests(Dq0.emm2.5)) %>% mutate(response = rep("Richness")),
                          as.data.frame(joint_tests(Dq1.emm1)) %>% mutate(response = rep("Shannon diversity")),
                          as.data.frame(joint_tests(Dq2.emm1)) %>% mutate(response = rep("Inv. Simpson diversity")),
                          as.data.frame(joint_tests(BP.emm1)) %>% mutate(response = rep("Dominance")),
                          as.data.frame(joint_tests(invSimpEven.emm1)) %>% mutate(response = rep("Evenness")),
                          as.data.frame(joint_tests(p.Ag.emm1)) %>% mutate(response = rep("proportion A. gerardii")),
                          as.data.frame(joint_tests(LegumeBiomass.emm1c)) %>% mutate(response = rep("Legume biomass")),
                          as.data.frame(joint_tests(ForbBiomass.emm1c)) %>% mutate(response = rep("Forb biomass")),
                          as.data.frame(joint_tests(C3Biomass.emm1c)) %>% mutate(response = rep("C3 biomass")),
                          as.data.frame(joint_tests(C4Biomass.emm1c)) %>% mutate(response = rep("C4 biomass"))) %>%
    # make column for significance
    mutate(sig = ifelse(p.value < 0.001, "***",
                        ifelse(p.value < 0.01 & p.value >= 0.001, "**",
                               ifelse(p.value < 0.05 & p.value >= 0.01, "*","")))) %>%
    dplyr::select(response, `model term`, df1, df2, F.ratio, p.value, sig)
  
  # Assemble dataframe for coefficients
  mmrm.coefs <- bind_rows(as.data.frame(Dq0.emm2.5) %>% mutate(response = rep("Dq0")),
                          as.data.frame(Dq1.emm1) %>% mutate(response = rep("Dq1")),
                          as.data.frame(Dq2.emm1) %>% mutate(response = rep("Dq2")),
                          as.data.frame(BP.emm1) %>% mutate(response = rep("BP")),
                          as.data.frame(invSimpEven.emm1) %>% mutate(response = rep("invSimpEven")),
                          as.data.frame(p.Ag.emm1) %>% mutate(response = rep("p.Ag")),
                          as.data.frame(LegumeBiomass.emm1) %>% mutate(response = rep("LegumeBiomass")),
                          as.data.frame(ForbBiomass.emm1) %>% mutate(response = rep("ForbBiomass")),
                          as.data.frame(C3Biomass.emm1) %>% mutate(response = rep("C3Biomass")),
                          as.data.frame(C4Biomass.emm1) %>% mutate(response = rep("C4Biomass"))) %>%
    filter(!is.na(df))
  
  # view dataframe of model coefficients
  print(head(mmrm.coefs))
  
  # Assemble dataframes for differences in coefficients
  # Warming effects only
  mmrm.diff.coefs <- bind_rows(as.data.frame(Dq0.contrasts) %>% mutate(response = rep("Dq0")),
                               as.data.frame(Dq1.contrasts) %>% mutate(response = rep("Dq1")),
                               as.data.frame(Dq2.contrasts) %>% mutate(response = rep("Dq2")),
                               as.data.frame(BP.contrasts) %>% mutate(response = rep("BP")),
                               as.data.frame(invSimpEven.contrasts) %>% mutate(response = rep("invSimpEven")),
                               as.data.frame(p.Ag.contrasts) %>% mutate(response = rep("p.Ag")),
                               as.data.frame(LegumeBiomass.contrasts) %>% mutate(response = rep("LegumeBiomass")),
                               as.data.frame(ForbBiomass.contrasts) %>% mutate(response = rep("ForbBiomass")),
                               as.data.frame(C3Biomass.contrasts) %>% mutate(response = rep("C3Biomass")),
                               as.data.frame(C4Biomass.contrasts) %>% mutate(response = rep("C4Biomass"))) %>%
    filter(!is.na(df))  %>%
    mutate(significance = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01  ~ "**",
      p.value < 0.05  ~ "*",
      p.value < 0.1 ~ "˙",
      TRUE            ~ ""
    ))
  
  # Warming and drought effects 
  mmrm.diff.coefs.all <- bind_rows(as.data.frame(Dq0.contrasts.all) %>% mutate(response = rep("Dq0")),
                                   as.data.frame(Dq1.contrasts.all) %>% mutate(response = rep("Dq1")),
                                   as.data.frame(Dq2.contrasts.all) %>% mutate(response = rep("Dq2")),
                                   as.data.frame(BP.contrasts.all) %>% mutate(response = rep("BP")),
                                   as.data.frame(invSimpEven.contrasts.all) %>% mutate(response = rep("invSimpEven")),
                                   as.data.frame(p.Ag.contrasts.all) %>% mutate(response = rep("p.Ag")),
                                   as.data.frame(LegumeBiomass.contrasts.all) %>% mutate(response = rep("LegumeBiomass")),
                                   as.data.frame(ForbBiomass.contrasts.all) %>% mutate(response = rep("ForbBiomass")),
                                   as.data.frame(C3Biomass.contrasts.all) %>% mutate(response = rep("C3Biomass")),
                                   as.data.frame(C4Biomass.contrasts.all) %>% mutate(response = rep("C4Biomass"))) %>%
    filter(!is.na(df)) %>%
    mutate(contrast = ifelse(contrast == "(+Heat Control) - Control Control","+Heat",
                             ifelse(contrast == "(Control -Water) - Control Control", "-Water",
                                    ifelse(contrast == "(+Heat -Water) - Control Control","+Heat:-Water", contrast))))  %>%
    mutate(significance = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01  ~ "**",
      p.value < 0.05  ~ "*",
      p.value < 0.1 ~ "˙",
      TRUE            ~ ""
    ))
  
  
}

################################################################################
# (8) Plot model coefficients ##################################################
################################################################################

# Ensure factor levels are ordered properly
df$HeatTrt <- factor(df$HeatTrt, levels = c("+Heat", "Control"))
df$WaterTrt <- factor(df$WaterTrt, levels = c("-Water", "Control"))

# Run all best-fit models for plotting
fit2.5.Dq0.factorYr.mmrm <- mmrm(Dq0~HeatTrt*WaterTrt*YrFac + csh(YrFac|sp), data=df )
fit2.Dq1.factorYr.mmrm <- mmrm(Dq1~HeatTrt*WaterTrt*YrFac + cs(YrFac|sp), data=df )
fit2.Dq2.factorYr.mmrm <- mmrm(Dq2~HeatTrt*WaterTrt*YrFac  + cs(YrFac|sp), data=df)
fit2.BP.factorYr.mmrm <- mmrm(BP~HeatTrt*WaterTrt*YrFac  + cs(YrFac|sp), data=df )
fit2.5.invSimpEven.factorYr.mmrm <- mmrm(invSimpEven~HeatTrt*WaterTrt*YrFac + csh(YrFac|sp), data=df )
fit2.p.Ag.factorYr.mmrm <- mmrm(p.Ag~HeatTrt*WaterTrt*YrFac + cs(YrFac|sp), data=df )
fit4.LegumeBiomass.factorYr.mmrm <- mmrm(LegumeBiomass~HeatTrt*WaterTrt*YrFac + ad(YrFac|sp), data=df)
fit3.ForbBiomass.factorYr.mmrm <- mmrm(ForbBiomass~HeatTrt*WaterTrt*YrFac + us(YrFac|sp), data=df)
fit5.C3Biomass.factorYr.mmrm <- mmrm(C3Biomass~HeatTrt*WaterTrt*YrFac + ar1h(YrFac|sp), data=df)
fit3.C4Biomass.factorYr.mmrm <- mmrm(C4Biomass~HeatTrt*WaterTrt*YrFac + us(YrFac|sp), data=df )

# plots for warming effects only
models <- list(
  Dq0 = fit2.5.Dq0.factorYr.mmrm,
  Dq1 = fit2.Dq1.factorYr.mmrm,
  Dq2 = fit2.Dq2.factorYr.mmrm,
  BP = fit2.BP.factorYr.mmrm,
  invSimpEven = fit2.5.invSimpEven.factorYr.mmrm,
  pAg = fit2.p.Ag.factorYr.mmrm
)

# Function to extract heat treatment effects for each model
extract_treatment_effects <- function(model, response_name) {
  # Get marginal means
  em_results <- emmeans(model, ~ HeatTrt | YrFac)
  
  # Compute Heat - Control effect
  heat_response <- as.data.frame(contrast(em_results, method = "pairwise", by = "YrFac"))
  
  # Ensure correct year ordering and add response name
  heat_response <- heat_response %>%
    mutate(YrFac = factor(YrFac, levels = c("2017", "2018", "2019", "2020", "2021", "2022")),
           Response = response_name,
           LCL = estimate - 1.96 * SE,  # 95% CI Lower
           UCL = estimate + 1.96 * SE)  # 95% CI Upper
  
  return(heat_response)
}

# Loop through all models and extract results
results_list <- map2(models, names(models), extract_treatment_effects)

# Combine all responses into a single dataframe
all_results <- bind_rows(results_list)
all_results$Response <- factor(all_results$Response, levels = c("Dq0", "Dq1", "Dq2", "invSimpEven", "BP", "pAg"))

# Define custom labels for facets
response_labels <- c(
  Dq0 = "Richness (q=0)",
  Dq1 = "Shannon diversity (q=1)",
  Dq2 = "Inverse Simpson diversity (q=2)",
  invSimpEven = "Evenness",
  BP = "Dominance",
  pAg = "<i>Andropogon gerardii</i> biomass prop."  # Uses HTML tags for italics
)

# # Extract warming responses for plotting
# # List of models and responses
# models <- list(
#   Dq0 = fit2.5.Dq0.factorYr.mmrm,
#   Dq1 = fit2.Dq1.factorYr.mmrm,
#   Dq2 = fit2.Dq2.factorYr.mmrm,
#   BP = fit2.BP.factorYr.mmrm,
#   invSimpEven = fit2.5.invSimpEven.factorYr.mmrm,
#   pAg = fit2.p.Ag.factorYr.mmrm
# )
# 
# # Function to extract treatment effects for each model
# extract_treatment_effects <- function(model, response_name) {
#   # Get marginal means
#   em_results <- emmeans(model, ~ HeatTrt | YrFac)
#   
#   # Compute Heat - Control effect
#   heat_response <- as.data.frame(contrast(em_results, method = "pairwise", by = "YrFac"))
#   
#   # Ensure correct year ordering and add response name
#   heat_response <- heat_response %>%
#     mutate(YrFac = factor(YrFac, levels = c("2017", "2018", "2019", "2020", "2021", "2022")),
#            Response = response_name,
#            LCL = estimate - 1.96 * SE,  # 95% CI Lower
#            UCL = estimate + 1.96 * SE)  # 95% CI Upper
#   
#   return(heat_response)
# }
# 
# # Loop through all models and extract results
# results_list <- map2(models, names(models), extract_treatment_effects)
# 
# # Combine all responses into a single dataframe
# all_results <- bind_rows(results_list)
# all_results$Response <- factor(all_results$Response, levels = c("Dq0", "Dq1", "Dq2", "invSimpEven", "BP", "pAg"))

# # Define custom labels for facets
# response_labels <- c(
#   Dq0 = "Richness (q=0)",
#   Dq1 = "Shannon diversity (q=1)",
#   Dq2 = "Inverse Simpson diversity (q=2)",
#   invSimpEven = "Evenness",
#   BP = "Dominance",
#   pAg = "<i>Andropogon gerardii</i> biomass prop."  # Uses HTML tags for italics
# )

# Filter for diversity responses
fig2_data <- all_results %>% filter(Response %in% c("Dq0", "Dq1", "Dq2", "invSimpEven"))

# Figure 2: Warming effects on community responses #############################
ggplot(fig2_data, aes(x = YrFac, y = estimate)) +
  geom_errorbar(aes(ymin = estimate - SE, ymax = estimate + SE),
                position = position_dodge(0.5), width = 0, size = 0.75) +
  geom_point(position = position_dodge(0.5), size = 2) +
  facet_wrap(~Response, scales = "free", ncol = 2, labeller = labeller(Response = response_labels)) +
  theme_few() +
  theme(strip.text = element_markdown(size = 10)) +  # Enables HTML rendering
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Year", y = "Community response to warming")

# Filter for BP & pAg
fig3_data <- all_results %>% filter(Response %in% c("BP", "pAg"))

# Figure 3: Dominance & A. gerardii responses to warming #######################
ggplot(fig3_data, aes(x = YrFac, y = estimate)) +
  geom_errorbar(aes(ymin = estimate - SE, ymax = estimate + SE),
                position = position_dodge(0.5), width = 0, size = 0.75) +
  geom_point(position = position_dodge(0.5), size = 2) +
  facet_wrap(~Response, scales = "free", ncol = 2, labeller = labeller(Response = response_labels)) +
  theme_few() +
  theme(strip.text = element_markdown(size = 10)) + 
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Year", y = "Community response to warming")

# Figure 4: Functional group responses to warming (no year) ###################

# list functional group models which do not include year effects
# we're doing this because we found no Year x Treatment interactive effects on legume responses
fit3.C4Biomass.factorYr.mmrm.c.noYr <- mmrm(C4Biomass~HeatTrt + us(YrFac|sp), data=df )
fit5.C3Biomass.factorYr.mmrm.c.noYr <- mmrm(C3Biomass~HeatTrt + ar1h(YrFac|sp), data=df )
fit4.ForbBiomass.factorYr.mmrm.c.noYr <- mmrm(ForbBiomass~HeatTrt + ad(YrFac|sp), data=df )
fit4.LegumeBiomass.factorYr.mmrm.c.noYr <- mmrm(LegumeBiomass~HeatTrt + ad(YrFac|sp), data=df )

# replicate across all FGs
C4Biomass.emm1 <- emmeans(fit3.C4Biomass.factorYr.mmrm.c.noYr,c('HeatTrt'))
joint_tests(C4Biomass.emm1) 
multcomp::cld(object = C4Biomass.emm1,Letters = letters,alpha = 0.05) 
C4Biomass.contrasts.noYr <- data.frame(summary(contrast(C4Biomass.emm1, method = "trt.vs.ctrl")))

C3Biomass.emm1 <- emmeans(fit5.C3Biomass.factorYr.mmrm.c.noYr,c('HeatTrt'))
joint_tests(C3Biomass.emm1) 
multcomp::cld(object = C3Biomass.emm1,Letters = letters,alpha = 0.05) 
C3Biomass.contrasts.noYr <- data.frame(summary(contrast(C3Biomass.emm1, method = "trt.vs.ctrl")))

ForbBiomass.emm1 <- emmeans(fit4.ForbBiomass.factorYr.mmrm.c.noYr,c('HeatTrt'))
joint_tests(ForbBiomass.emm1) 
multcomp::cld(object = ForbBiomass.emm1,Letters = letters,alpha = 0.05) 
ForbBiomass.contrasts.noYr <- data.frame(summary(contrast(ForbBiomass.emm1, method = "trt.vs.ctrl")))

LegumeBiomass.emm1 <- emmeans(fit4.LegumeBiomass.factorYr.mmrm.c.noYr,c('HeatTrt'))
joint_tests(LegumeBiomass.emm1) 
multcomp::cld(object = LegumeBiomass.emm1,Letters = letters,alpha = 0.05)
LegumeBiomass.contrasts.noYr <- data.frame(summary(contrast(LegumeBiomass.emm1, method = "trt.vs.ctrl")))

# Bind into dataframe
mmrm.diff.coefs.fg <- bind_rows(as.data.frame(LegumeBiomass.contrasts.noYr) %>% mutate(response = rep("LegumeBiomass")),
                                as.data.frame(ForbBiomass.contrasts.noYr) %>% mutate(response = rep("ForbBiomass")),
                                as.data.frame(C3Biomass.contrasts.noYr) %>% mutate(response = rep("C3Biomass")),
                                as.data.frame(C4Biomass.contrasts.noYr) %>% mutate(response = rep("C4Biomass"))) %>%
  filter(!is.na(df))  %>%
  mutate(significance = case_when(
    p.value < 0.001 ~ "***",
    p.value < 0.01  ~ "**",
    p.value < 0.05  ~ "*",
    p.value < 0.1 ~ "˙",
    TRUE            ~ ""
  ))

fig3b.supp.df <- mmrm.diff.coefs.fg %>%
  mutate(response = fct_relevel(response,
                                "LegumeBiomass","ForbBiomass","C3Biomass","C4Biomass")) 

ggplot(fig3b.supp.df, aes(x = response, y = estimate)) +
  geom_errorbar(aes(ymin = estimate - SE, ymax = estimate + SE),
                position = position_dodge(0.5), width = 0, size = 0.75) +
  geom_point(position = position_dodge(0.5), size = 2) +
  scale_shape_manual(breaks = c("Dry", "Normal", "Wet"),
                     values = c(15, 16, 17), name = "Precipitation") +
  theme_few() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Functional group", y = expression("Biomass proportion response to warming (difference from Control,  g/m"^phantom(0)^2*", %)")) + 
  scale_x_discrete(labels = c("LegumeBiomass" = expression("Legumes"),
                              "ForbBiomass" = expression("Forbs"),
                              "C3Biomass" = expression("C"["3"]~"grasses"),
                              "C4Biomass" = expression("C"["4"]~"grasses")))

################################################################################
# (9) Supplementary figures ####################################################
################################################################################
# Figure A1: Soil Moisture Response to Warming and Drought #####################

# view data
str(smt1722.mod) # note that there is no data from 2018 due to measurements not being taken before 11am that year.

# summarize
smt1722.mod.sum <- smt1722.mod %>%
  group_by(LTERYEAR,Exp,Plot,response,measure,YrFac,Trtmt,Subplot,sp) %>%
  summarise_all(mean) %>%
  as.data.frame() %>%
  mutate_if(is.double, as.factor)

# change characters to factors for plotting
smt1722.mod.sum <- smt1722.mod.sum %>%
  as.data.frame() %>%
  mutate(YrFac = as.factor(YrFac)) %>%
  mutate_if(is.character, as.factor) 

# find most parsimonious temporal autocorrelation structure for soil moisture responses
fit1.soilM_m3m3.factorYr.mmrm <- mmrm(value~ 0+response*YrFac + ar1(YrFac|sp), data=smt1722.mod.sum)   # ar1 = homogeneous autoregressive covariance strs assume that "the correlation between any two elements is equal to r for adjacent elements". In other words, covariances between time points decline exponentially.
fit2.soilM_m3m3.factorYr.mmrm <- mmrm(value~0+response*YrFac + cs(YrFac|sp), data=smt1722.mod.sum) # cs = homogeneous compound symmetry assumes a constant correlation between time points
#fit2.5.soilM_m3m3.factorYr.mmrm <- mmrm(value~0+response*YrFac + csh(YrFac|sp), data=smt1722.mod.sum) # MAGGIE ADDED # csh = heterogeneous compound symmetry # BEST fit otherwise
fit3.soilM_m3m3.factorYr.mmrm <- mmrm(value~0+response*YrFac + us(YrFac|sp), data=smt1722.mod.sum)    # us = unstructured covariance str
fit4.soilM_m3m3.factorYr.mmrm <- mmrm(value~0+response*YrFac + ad(YrFac|sp), data=smt1722.mod.sum)    # ad = homogeneous covariance str, observations are approx. equally spaced in time 
fit5.soilM_m3m3.factorYr.mmrm <- mmrm(value~0+response*YrFac + ar1h(YrFac|sp), data=smt1722.mod.sum)  # ar1h = same as hr 1, but assumes heterogeneous variance
AIC(fit1.soilM_m3m3.factorYr.mmrm)
AIC(fit2.soilM_m3m3.factorYr.mmrm)
#AIC(fit2.5.soilM_m3m3.factorYr.mmrm) 
AIC(fit3.soilM_m3m3.factorYr.mmrm) # best fit
AIC(fit4.soilM_m3m3.factorYr.mmrm)
AIC(fit5.soilM_m3m3.factorYr.mmrm)

# get most parsimonious model form
fit3.soilM_m3m3.factorYr.mmrm.a <- mmrm(value~response*YrFac + us(YrFac|sp), data=smt1722.mod.sum)
fit3.soilM_m3m3.factorYr.mmrm.b <- mmrm(value~response+YrFac + us(YrFac|sp), data=smt1722.mod.sum)
AIC(fit3.soilM_m3m3.factorYr.mmrm.a) # fully-interactive model fits best
AIC(fit3.soilM_m3m3.factorYr.mmrm.b)
# Best fit for soilM_m3m3 is a US fully-interactive model (tested above)
summary(fit3.soilM_m3m3.factorYr.mmrm.a)

# Prepare model coefficients for plotting
smt1722.mod.lmer <- smt1722.mod %>%
  dplyr::mutate(YrFac = ifelse(YrFac == 1, 2017,
                               ifelse(YrFac == 2, 2018,
                                      ifelse(YrFac == 3, 2019, 
                                             ifelse(YrFac == 4, 2020, 
                                                    ifelse(YrFac == 5, 2021, 
                                                           ifelse(YrFac == 6, 2022, NA)))))),
                YrFac = as.factor(YrFac),
                Plot = as.factor(Plot)) 

fit3.soilM_m3m3.factorYr.lmer <- lmer(value ~ response * YrFac + (1 |Plot/Subplot), data = smt1722.mod.lmer)

# Extract the estimated marginal means (EMMs) and standard errors for the interaction
soilm_coefs <- as.data.frame(emmeans(fit3.soilM_m3m3.factorYr.lmer, ~ response * YrFac, method = "pairwise", reverse = TRUE))

# Plot coefficients
soilm_coefs %>%
  ggplot(., aes(x = as.factor(YrFac), y = emmean, col = response)) +
  geom_errorbar(aes(ymin=emmean-SE, ymax=emmean+SE),
                position=position_dodge(0.5), width = 0.25, size = 0.75) +
  geom_point(position=position_dodge(0.5), size = 2.5) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_few() +
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank()) +
  scale_color_manual(values = c("#92C5DE","#F4A582","#CA0020"), name = "Treatment", labels = c("-Water", "+Heat", "+Heat -Water")) +
  geom_hline(yintercept = 0, linetype="dashed") +
  ylab(bquote('Soil moisture response')) +
  scale_shape_manual(values=c(15,16,17), name = "Yearly precip.") +
  labs(x = "Year", y = "Soil moisture change (treatment response)")

# Figure A2: Biomass responses by species ######################################
fills <- c("#117733","#117733","#117733","#117733", 
           "#D0CC77", "#D0CC77", "#D0CC77", "#D0CC77", 
           "#000088", "#000088", "#000088", "#000088",
           "#882255", "#882255","#882255") 


{# Control - -Water, all years    
  Means1 <- Means %>%
    filter(Trtmt %in% c("Control", "-Water")) 
  cw<-ggplot() +
    geom_errorbar(data = Means1, aes(x=factor(Trtmt, level = c("Control","-Water", "+Heat", "+Heat-Water")),
                                     ymin=TotalBiomassLog_mean-TotalBiomassLog_se,
                                     ymax=TotalBiomassLog_mean+TotalBiomassLog_se,
                                     col = interaction(Species, FunctionalGroup), group = interaction(Species)), width=0,
                  size = 1, position = position_dodge(width = 0.5), alpha=0.2)+
    geom_point(data = Means1, aes(x=factor(Trtmt, level = c("Control","-Water", "+Heat", "+Heat-Water")),
                                  y=TotalBiomassLog_mean,
                                  col = interaction(Species, FunctionalGroup), group = interaction(Species),
                                  shape = interaction(Species, FunctionalGroup)), 
               position = position_dodge(width = 0.5), size=1)+
    geom_line(data = Means1, aes(x = factor (Trtmt, level = c("Control","-Water", "+Heat", "+Heat-Water")), 
                                 y = TotalBiomassLog_mean,
                                 col = interaction(Species, FunctionalGroup), group = interaction(Species)), alpha = 0.5,show.legend = F,
              position = position_dodge(width = 0.5)) +
    theme_bw() +
    labs(x = "", y= "",
         color = "Species & functional group",
         shape = "Species & functional group") +
    # ylab(bquote('log(1 + biomass '~gm^-2~')'))+
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid.minor.x = element_blank()) +
    facet_grid(.~Year) + 
    scale_fill_manual(values = fills) +
    scale_shape_manual(values=c(1:4,1:4,1:4,1:3)) +
    scale_color_manual(values = fills) +
    theme(strip.background=element_rect(fill="white"))
  
  # Control - +Heat, all years    
  Means2 <- Means %>%
    dplyr::filter(Trtmt %in% c("Control", "+Heat")) 
  ch<-ggplot() +
    geom_errorbar(data = Means2, aes(x=factor(Trtmt, level = c("Control","-Water", "+Heat", "+Heat-Water")),
                                     ymin=TotalBiomassLog_mean-TotalBiomassLog_se,
                                     ymax=TotalBiomassLog_mean+TotalBiomassLog_se,
                                     col = interaction(Species, FunctionalGroup), group = interaction(Species)), width=0,
                  size = 1, position = position_dodge(width = 0.5), alpha=0.2)+
    geom_point(data = Means2, aes(x=factor(Trtmt, level = c("Control","-Water", "+Heat", "+Heat-Water")),
                                  y=TotalBiomassLog_mean,
                                  col = interaction(Species, FunctionalGroup), group = interaction(Species),
                                  shape = interaction(Species, FunctionalGroup)), 
               position = position_dodge(width = 0.5), size=1)+
    geom_line(data = Means2, aes(x = factor (Trtmt, level = c("Control","-Water", "+Heat", "+Heat-Water")), 
                                 y = TotalBiomassLog_mean,
                                 col = interaction(Species, FunctionalGroup), group = interaction(Species)), alpha = 0.5,show.legend = F,
              position = position_dodge(width = 0.5)) +
    theme_bw() +
    labs(x = "", fill = "Species & functional group") +
    ylab(bquote('log(1 + biomass '~gm^-2~')'))+
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid.minor.x = element_blank()) +
    facet_grid(.~Year) +
    scale_fill_manual(values = fills) +
    scale_shape_manual(values=c(1:4,1:4,1:4,1:3)) +
    scale_color_manual(values = fills) +
    theme(strip.background=element_rect(fill="white"))
  
  # Control - +Heat-Water, all years    
  Means3 <- Means %>%
    filter(Trtmt %in% c("Control", "+Heat-Water")) 
  cwh<-ggplot() +
    geom_errorbar(data = Means3, aes(x=factor(Trtmt, level = c("Control","-Water", "+Heat", "+Heat-Water")),
                                     ymin=TotalBiomassLog_mean-TotalBiomassLog_se,
                                     ymax=TotalBiomassLog_mean+TotalBiomassLog_se,
                                     col = interaction(Species, FunctionalGroup), group = interaction(Species)), width=0,
                  size = 1, position = position_dodge(width = 0.5), alpha=0.2)+
    geom_point(data = Means3, aes(x=factor(Trtmt, level = c("Control","-Water", "+Heat", "+Heat-Water")),
                                  y=TotalBiomassLog_mean,
                                  col = interaction(Species, FunctionalGroup), group = interaction(Species),
                                  shape = interaction(Species, FunctionalGroup)), 
               position = position_dodge(width = 0.5), size=1)+
    geom_line(data = Means3, aes(x = factor (Trtmt, level = c("Control","-Water", "+Heat", "+Heat-Water")), 
                                 y = TotalBiomassLog_mean,
                                 col = interaction(Species, FunctionalGroup), group = interaction(Species)), alpha = 0.5,show.legend = F,
              position = position_dodge(width = 0.5)) +
    theme_bw() +
    labs(x = "", y = " ") +
    theme(axis.text.x = element_text(angle = 30, hjust = 1),
          panel.grid.minor.x = element_blank()) +
    facet_grid(.~Year) + 
    scale_fill_manual(values = fills) +
    scale_shape_manual(values=c(1:4,1:4,1:4,1:3)) +
    scale_color_manual(values = fills) +
    theme(strip.background=element_rect(fill="white"));cwh
}
ggarrange(cw+rremove("xlab"),ch+rremove("xlab"),cwh+rremove("xlab"), ncol= 1, common.legend = TRUE, legend = "right") # LogBiomassSps_allYears_20231009
ggsave('LogBiomassSps_allYears_20231009.png',width = 8, height = 8, dpi=300)
# Function to extract ALL treatment effects for each model
extract_treatment_effects <- function(model, response_name) {
  # Get marginal means
  em_results <- emmeans(model, ~ HeatTrt * WaterTrt | YrFac) 
  
  # Compute Heat - Control effect
  trt_response <- as.data.frame(contrast(em_results, method = "pairwise", by = "YrFac"))
  
  # Ensure correct year ordering and add response name
  trt_response <- trt_response %>%
    mutate(YrFac = factor(YrFac, levels = c("2017", "2018", "2019", "2020", "2021", "2022")),
           Response = response_name,
           LCL = estimate - 1.96 * SE,  # 95% CI Lower
           UCL = estimate + 1.96 * SE)  # 95% CI Upper
  
  return(trt_response)
}

# Loop through all models and extract results
results_list <- map2(models, names(models), extract_treatment_effects)

# Combine all responses into a single dataframe
all_results <- bind_rows(results_list)

# rename factor levels of all_results$contrast
all_results<- all_results %>%
  dplyr::filter(grepl("Control Control", contrast)) %>%
  mutate(contrast = factor(contrast, levels=c("Control Control - (Control -Water)","Control Control - (+Heat Control)", "Control Control - (+Heat -Water)")))

all_results$Response <- factor(all_results$Response, levels = c("Dq0", "Dq1", "Dq2", "invSimpEven", "BP", "pAg", "Legume","Forb","C3","C4"))

# Define custom labels for facets
response_labels <- c(
  Dq0 = "Richness (q=0)",
  Dq1 = "Shannon diversity (q=0)",
  Dq2 = "Inverse Simpson diversity (q=0)",
  invSimpEven = "Evenness",
  BP = "Dominance",
  pAg = "<i>A. gerardii</i> biomass",  # Italicized species name
  Legume = "Legumes",
  Forb = "Forbs",
  C3 = "C<sub>3</sub> grasses",
  C4 = "C<sub>4</sub> grasses"
)

# Filter for diversity responses
plot1_data <- all_results %>% filter(Response %in% c("Dq0", "Dq1", "Dq2", "invSimpEven"))

# reset factor levels
df$HeatTrt <- factor(df$HeatTrt, levels = c("+Heat", "Control"))
df$WaterTrt <- factor(df$WaterTrt, levels = c("-Water", "Control"))

# List of models and responses
models <- list(
  Dq0 = fit2.5.Dq0.factorYr.mmrm,
  Dq1 = fit2.Dq1.factorYr.mmrm,
  Dq2 = fit2.Dq2.factorYr.mmrm,
  BP = fit2.BP.factorYr.mmrm,
  invSimpEven = fit2.5.invSimpEven.factorYr.mmrm,
  pAg = fit2.p.Ag.factorYr.mmrm,
  Legume = fit4.LegumeBiomass.factorYr.mmrm,
  Forb = fit3.ForbBiomass.factorYr.mmrm,
  C3 = fit5.C3Biomass.factorYr.mmrm,
  C4 = fit3.C4Biomass.factorYr.mmrm
)

# Function to extract treatment effects for each model
extract_treatment_effects <- function(model, response_name) {
  # Get marginal means
  em_results <- emmeans(model, ~ HeatTrt * WaterTrt | YrFac) 
  
  # Compute Heat - Control effect
  trt_response <- as.data.frame(contrast(em_results, method = "pairwise", by = "YrFac"))
  
  # Ensure correct year ordering and add response name
  trt_response <- trt_response %>%
    mutate(YrFac = factor(YrFac, levels = c("2017", "2018", "2019", "2020", "2021", "2022")),
           Response = response_name,
           LCL = estimate - 1.96 * SE,  # 95% CI Lower
           UCL = estimate + 1.96 * SE)  # 95% CI Upper
  
  return(trt_response)
}

# Loop through all models and extract results
results_list <- map2(models, names(models), extract_treatment_effects)

# Combine all responses into a single dataframe
all_results <- bind_rows(results_list)

# rename factor levels of all_results$contrast
all_results<- all_results %>%
  filter(grepl("Control Control", contrast)) %>%
  #dplyr::mutate(contrast = sub(".*\\(([^)]+)\\).*", "\\1", contrast),
  mutate(contrast = factor(contrast, levels=c("(Control -Water) - Control Control","(+Heat Control) - Control Control", "(+Heat -Water) - Control Control")))

all_results$Response <- factor(all_results$Response, levels = c("Dq0", "Dq1", "Dq2", "invSimpEven", "BP", "pAg", "Legume","Forb","C3","C4"))

# Define custom labels for facets
response_labels <- c(
  Dq0 = "Richness (q=0)",
  Dq1 = "Shannon diversity (q=0)",
  Dq2 = "Inverse Simpson diversity (q=0)",
  invSimpEven = "Evenness",
  BP = "Dominance",
  pAg = "<i>A. gerardii</i> biomass",  # Italicized species name
  Legume = "Legumes",
  Forb = "Forbs",
  C3 = "C<sub>3</sub> grasses",
  C4 = "C<sub>4</sub> grasses"
)

# Figure A3: All treatment effects on community responses ######################
# Filter for diversity responses
figA3_data <- all_results %>% filter(Response %in% c("Dq0", "Dq1", "Dq2", "invSimpEven"))

# Generate Plot 1
ggplot(figA3_data, aes(x = YrFac, y = estimate, col = contrast)) +
  geom_errorbar(aes(ymin = estimate - SE, ymax = estimate + SE),
                position = position_dodge(0.5), width = 0, size = 0.75) +
  geom_point(position = position_dodge(0.5), size = 2) +
  facet_wrap(~Response, scales = "free", ncol = 2, labeller = labeller(Response = response_labels)) +
  theme_few() +
  theme(strip.text = element_markdown(size = 10)) +  # Enables HTML rendering
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Year", y = "Community response to warming") +
  scale_color_manual(values = c("#92C5DE", "#F39E78", "#CA0020"),
                     name = "Treatment effect",labels=c("(Control -Water) - Control Control"="-Water",
                                                        "(+Heat Control) - Control Control"="+Heat",
                                                        "(+Heat -Water) - Control Control" ="+Heat-Water")) 

# Figure A4: All treatment effects on dominance/A. gerardii biomass ############
# Filter for BP & pAg
figA4_data <- all_results %>% filter(Response %in% c("BP", "pAg"))
ggplot(figA4_data, aes(x = YrFac, y = estimate, col = contrast)) +
  geom_errorbar(aes(ymin = estimate - SE, ymax = estimate + SE),
                position = position_dodge(0.5), width = 0, size = 0.75) +
  geom_point(position = position_dodge(0.5), size = 2) +
  facet_wrap(~Response, scales = "free", ncol = 2, labeller = labeller(Response = response_labels)) +
  theme_few() +
  theme(strip.text = element_markdown(size = 10)) +  # Enables HTML rendering
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Year", y = "Community response to warming") +
  scale_color_manual(values = c("#92C5DE", "#F39E78", "#CA0020"),
                     name = "Treatment effect",labels=c("(Control -Water) - Control Control"="-Water",
                                                        "(+Heat Control) - Control Control"="+Heat",
                                                        "(+Heat -Water) - Control Control"="-Water+Heat"))

# Figure A5: All treatment effects on functional group biomass proportions #####
# Filter for functional groups
FigA5_data <- all_results %>% filter(Response %in% c("Legume","Forb","C3","C4"))

ggplot(FigA5_data, aes(x = YrFac, y = estimate, col = contrast)) +
  geom_errorbar(aes(ymin = estimate - SE, ymax = estimate + SE),
                position = position_dodge(0.5), width = 0, size = 0.75) +
  geom_point(position = position_dodge(0.5), size = 2) +
  facet_wrap(~Response, scales = "free", ncol = 2, labeller = labeller(Response = response_labels)) +
  theme_few() +
  theme(strip.text = element_markdown(size = 10)) +  # Enables HTML rendering
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Year", y = "Community response to warming") +
  scale_color_manual(values = c("#92C5DE", "#F39E78", "#CA0020"),
                     name = "Treatment effect",labels=c("(Control -Water) - Control Control"="-Water",
                                                        "(+Heat Control) - Control Control"="+Heat",
                                                        "(+Heat -Water) - Control Control"="-Water+Heat")) +
  labs(x="Year", y = expression("Biomass proportion response to warming (difference from Control,  g/m"^phantom(0)^2*", %)"))

# Figure A7: Model assessment ##################################################
par(mfrow=c(3,2))
# Dq0
get_varcov(fit2.5.Dq0.factorYr.mmrm) # variance-covariance matrix
residuals <- get_residuals(fit2.5.Dq0.factorYr.mmrm) # get residuals
fitted.values <- get_predicted(fit2.5.Dq0.factorYr.mmrm) # get fitted values (i.e., predictions for the response)
j <- data.frame(cbind(fitted.values,residuals))
scatter.smooth(j$residuals ~ j$fitted.values, span = 2/3, degree = 1,col="blue", main="Dq0, mmrm")

# Inverse Simpson Evenness
get_varcov(fit2.5.invSimpEven.factorYr.mmrm) # variance-covariance matrix
residuals <- get_residuals(fit2.5.invSimpEven.factorYr.mmrm) # residuals
fitted.values <- get_predicted(fit2.5.invSimpEven.factorYr.mmrm) # fitted values 
j <- data.frame(cbind(fitted.values,residuals))
scatter.smooth(j$residuals ~ j$fitted.values, span = 2/3, degree = 1,col="blue", main="Evenness, mmrm")

# Dq1
get_varcov(fit2.Dq1.factorYr.mmrm) # variance-covariance matrix
residuals <- get_residuals(fit2.Dq1.factorYr.mmrm) # residuals
fitted.values <- get_predicted(fit2.Dq1.factorYr.mmrm) # fitted values 
j <- data.frame(cbind(fitted.values,residuals))
scatter.smooth(j$residuals ~ j$fitted.values, span = 2/3, degree = 1,col="blue", main="Dq1, mmrm")

# BP
get_varcov(fit2.BP.factorYr.mmrm) # variance-covariance matrix
residuals <- get_residuals(fit2.BP.factorYr.mmrm) # residuals
fitted.values <- get_predicted(fit2.BP.factorYr.mmrm) # fitted values 
j <- data.frame(cbind(fitted.values,residuals))
scatter.smooth(j$residuals ~ j$fitted.values, span = 2/3, degree = 1,col="blue", main="BP, mmrm")

# Dq2
get_varcov(fit2.Dq2.factorYr.mmrm) # variance-covariance matrix
residuals <- get_residuals(fit2.Dq2.factorYr.mmrm) # residuals
fitted.values <- get_predicted(fit2.Dq2.factorYr.mmrm) # fitted values 
j <- data.frame(cbind(fitted.values,residuals))
scatter.smooth(j$residuals ~ j$fitted.values, span = 2/3, degree = 1,col="blue", main="Dq2, mmrm")

# p.Ag
get_varcov(fit2.p.Ag.factorYr.mmrm) # variance-covariance matrix
residuals <- get_residuals(fit2.5.p.Ag.factorYr.mmrm) # residuals
fitted.values <- get_predicted(fit2.5.p.Ag.factorYr.mmrm) # fitted values 
j <- data.frame(cbind(fitted.values,residuals))
scatter.smooth(j$residuals ~ j$fitted.values, span = 2/3, degree = 1,col="blue", main="Prop. Andge, mmrm")

par(mfrow=c(2,2))
# Legume
get_varcov(fit2.LegumeBiomass.factorYr.mmrm) # variance-covariance matrix
residuals <- get_residuals(fit2.LegumeBiomass.factorYr.mmrm) # residuals
fitted.values <- get_predicted(fit2.LegumeBiomass.factorYr.mmrm) # fitted values 
j <- data.frame(cbind(fitted.values,residuals))
scatter.smooth(j$residuals ~ j$fitted.values, span = 2/3, degree = 1,col="blue", main="Legume, mmrm")

# Forb
get_varcov(fit3.ForbBiomass.factorYr.mmrm) # variance-covariance matrix
residuals <- get_residuals(fit3.ForbBiomass.factorYr.mmrm) # residuals
fitted.values <- get_predicted(fit3.ForbBiomass.factorYr.mmrm) # fitted values 
j <- data.frame(cbind(fitted.values,residuals))
scatter.smooth(j$residuals ~ j$fitted.values, span = 2/3, degree = 1,col="blue", main="Forb, mmrm")

# C3 grass
get_varcov(fit5.C3Biomass.factorYr.mmrm) # variance-covariance matrix
residuals <- get_residuals(fit5.C3Biomass.factorYr.mmrm) # residuals
fitted.values <- get_predicted(fit5.C3Biomass.factorYr.mmrm) # fitted values 
j <- data.frame(cbind(fitted.values,residuals))
scatter.smooth(j$residuals ~ j$fitted.values, span = 2/3, degree = 1,col="blue", main="C3 grass, mmrm")

# C4 grass
get_varcov(fit2.C4Biomass.factorYr.mmrm) # variance-covariance matrix
residuals <- get_residuals(fit2.C4Biomass.factorYr.mmrm) # residuals
fitted.values <- get_predicted(fit2.C4Biomass.factorYr.mmrm) # fitted values 
j <- data.frame(cbind(fitted.values,residuals))
scatter.smooth(j$residuals ~ j$fitted.values, span = 2/3, degree = 1,col="blue", main="C4 grass, mmrm")

par(mfrow=c(1,1))
# Figure A8: Coefficient of Variation ##########################################
# Dataframe for CV plot
comm.cv <- cv_results %>% 
  dplyr::filter(Measure %in% c("BP","Dq0","Dq1","Dq2","invSimpEven","p.Ag")) %>%
  mutate(Measure = factor(Measure, 
                          levels = c("Dq0", "Dq1", "invSimpEven", "Dq2", "BP", "p.Ag")))

# Create a mapping of labels
measure_labels <- c("Dq0" = "Richness~(q==0)",
                    "Dq1" = "Shannon~diversity~(q==1)",
                    "invSimpEven" = "Evenness",
                    "Dq2" = "Inv.~Simpson~diversity~(q==2)",
                    "BP" = "Dominance",
                    "p.Ag" = "italic(A.~gerardii)~biomass~proportion")

comm_cv_plot <- ggplot(comm.cv, aes(x = interaction(WaterTrt, HeatTrt), y = cv, fill = interaction(WaterTrt, HeatTrt), col = interaction(WaterTrt, HeatTrt))) + 
  stat_summary(fun=mean, geom="point", position = position_dodge(width = 0.5), size=2) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.25, size = 0.75, position = position_dodge(width = 0.5)) + 
  theme_bw() + 
  ggtitle("A") +
  ylab('Coefficient of Variation') +
  scale_fill_manual(values = c("#0571B0", "#92C5DE", "#F39E78", "#CA0020"), 
                    name = "Treatment", 
                    labels = c("Control.Control" = "Control", 
                               "Drought.Control" = "-Water", 
                               "Control.Heated" = "+Heat", 
                               "Drought.Heated" = "+Heat-Water")) +
  scale_color_manual(values = c("#0571B0", "#92C5DE", "#F39E78", "#CA0020"), 
                     name = "Treatment", 
                     labels = c("Control.Control" = "Control", 
                                "Drought.Control" = "-Water", 
                                "Control.Heated" = "+Heat", 
                                "Drought.Heated" = "+Heat-Water")) + 
  theme(strip.background = element_rect(fill="white"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right") +
  facet_wrap(~ Measure, labeller = as_labeller(measure_labels, label_parsed), nrow = 3); comm_cv_plot


# Create factor without expressions
fg.cv <- cv_results %>% 
  filter(Measure %in% c("LegumeBiomass","ForbBiomass","C3Biomass","C4Biomass")) %>%
  mutate(Measure = factor(Measure, 
                          levels = c("LegumeBiomass", "ForbBiomass", "C3Biomass", "C4Biomass")))

# Mapping of labels with expressions
measure_labels_fg <- c("LegumeBiomass" = "Legumes",
                       "ForbBiomass" = "Forbs",
                       "C3Biomass" = "C[3]~grasses",
                       "C4Biomass" = "C[4]~grasses")

# Adjust ggplot code to include parsing
fg_cv_plot<-ggplot(fg.cv, aes(x = interaction(WaterTrt, HeatTrt), y = cv, fill = interaction(WaterTrt, HeatTrt), col = interaction(WaterTrt, HeatTrt))) + 
  stat_summary(fun=mean, geom="point", position = position_dodge(width = 0.5), size=2) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.25, size = 0.75, position = position_dodge(width = 0.5)) + 
  theme_bw() + 
  ylab('Coefficient of Variation') +
  scale_fill_manual(values = c("#0571B0", "#92C5DE", "#F39E78", "#CA0020"), 
                    name = "Treatment", 
                    labels = c("Control.Control" = "Control", 
                               "Drought.Control" = "-Water", 
                               "Control.Heated" = "+Heat", 
                               "Drought.Heated" = "+Heat-Water")) +
  scale_color_manual(values = c("#0571B0", "#92C5DE", "#F39E78", "#CA0020"), 
                     name = "Treatment", 
                     labels = c("Control.Control" = "Control", 
                                "Drought.Control" = "-Water", 
                                "Control.Heated" = "+Heat", 
                                "Drought.Heated" = "+Heat-Water")) + 
  theme(strip.background = element_rect(fill="white"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right") +
  facet_wrap(~ Measure, labeller = as_labeller(measure_labels_fg, label_parsed)) +
  ggtitle("B")

# Arrange these plots with the previous ones using ggarrange
ggarrange(comm_cv_plot, fg_cv_plot, ncol = 1, common.legend = TRUE, legend = "right", heights = c(3,2))

# Celebratory frog if you made it to the end of my code
frog <- "\U0001F438"; cat(frog,  "\n") 
