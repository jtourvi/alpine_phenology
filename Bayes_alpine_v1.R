#Code for Tourville et al. 2026 Northeastern North America (NNA) 
#alpine phenology analysis and figures

#set local working directory here

##### Load required and useful packages ######

require(ggplot2)
require(sf)
require(sp)
require(rgdal)
require(grid)
require(gridExtra)
require(egg)
require(rnaturalearth)
require(stars)
require(viridis)
require(R2jags) 
require(plotly)
require(dplyr)
require(tidyr)
require(tidyverse)
require(AICcmodavg)
require(cowplot)
require(car)
require(patchwork)
require(reshape)
require(sjPlot)
require(sjstats)
require(lme4)
require(chillR)
require(lubridate)
require(ggpmisc)
require(loo)
require(brms)
require(ggcorrplot)
require(rnaturalearth)
require(rnaturalearthdata)
require(terra)
require(tidyterra)
require(ggspatial)
require(figpatch)

##### Read in data ######

#Need DOY.df.csv
#Need fig2_data.csv
#Need year_change.csv

# See metadata for variable descriptions.

DOY.df <- read.csv('DOY.df.csv',
                    header = TRUE, sep = ",")

#DOY.df$diff = DOY.df$doy_first - DOY.df$SMOD #make DOY diff between SMOD and observation date
#DOY.df = merge(DOY.df, lc3, by = c("ids", "year_first"), all.x = TRUE)#merge in local climate data for subset of close observations if not done yet

#subseting of data for useful analyses
DOY.df2 = subset(DOY.df, select = c("ids", "scientific_name", "group", "type", "code", "source", "peak", "year_first", "doy_first",
                                   "TPI", "Tmean_SMOD", "Tmean_SMOD_l", "SMOD_class", "diff",
                                   "Tmean_winter"))
DOY.df2 = na.omit(DOY.df2)
DOY.df = subset(DOY.df, select = c("ids", "scientific_name", "group", "type", "code", "source", "peak", "year_first", "doy_first",
                                   "TPI", "Tmean_SMOD", "SMOD_class", "diff",
                                   "Tmean_winter"))
DOY.df = na.omit(DOY.df) # need to calculate GDD just between SMOD and DOY for every observation

#Clean data and add climate, topographic, and snow variables

#Scale predictor variables (Temp, SMOD, Elevation, TSI, TPI, SR)
#Turn on and off annoatation based on which variables are used

#DOY.df$Tmean_spring = scale(DOY.df$Tmean_spring)
#DOY.df$SMOD = scale(DOY.df$SMOD)
DOY.df$Tmean_winter = scale(DOY.df$Tmean_winter)
#DOY.df2$Tmean_winter = scale(DOY.df2$Tmean_winter)
#DOY.df$gdd_0_spring_accum = scale(DOY.df$gdd_0_spring_accum)
#DOY.df$WindShelter = scale(DOY.df$WindShelter)
DOY.df$TPI = scale(DOY.df$TPI)
#DOY.df2$TPI = scale(DOY.df2$TPI)
DOY.df$Prcp_spring = scale(DOY.df$Prcp_spring)
DOY.df$Swe_spring = scale(DOY.df$Swe_spring)
#DOY.df$Tmax_spring = scale(DOY.df$Tmax_spring)
#DOY.df$Tmin_spring = scale(DOY.df$Tmin_spring)
DOY.df$Tmean_SMOD = scale(DOY.df$Tmean_SMOD)
#DOY.df2$Tmean_SMOD = scale(DOY.df2$Tmean_SMOD)
DOY.df$gdd_0_SMOD = scale(DOY.df$gdd_0_SMOD)
#DOY.df2$Tmean_SMOD_l = scale(DOY.df2$Tmean_SMOD_l)
#DOY.df$gdd_0_SMOD_l = scale(DOY.df$gdd_0_SMOD_l)
DOY.df$SMOD_class = scale(DOY.df$SMOD_class)
DOY.df$Vp_spring = scale(DOY.df$Vp_spring)
DOY.df$elevation = scale(DOY.df$elevation)
DOY.df$aspect = scale(DOY.df$aspect)
DOY.df$slope = scale(DOY.df$slope)
#DOY.df2$SMOD_class = scale(DOY.df2$SMOD_class)
#DOY.df$diff = scale(DOY.df$diff)

#Species only df for later
#Alpine
df.AL.F1 = subset(DOY.df, scientific_name == "Carex bigelowii")
df.AL.F2 = subset(DOY.df, scientific_name == "Clintonia borealis")
df.AL.F3 = subset(DOY.df, scientific_name == "Cornus canadensis")
df.AL.F4 = subset(DOY.df, scientific_name == "Diapensia lapponica")
df.AL.F5 = subset(DOY.df, scientific_name == "Geum peckii")
df.AL.F6 = subset(DOY.df, scientific_name == "Kalmia procumbens")
df.AL.F7 = subset(DOY.df, scientific_name == "Maianthemum canadense")
df.AL.F8 = subset(DOY.df, scientific_name == "Rhododendron groenlandicum")
df.AL.F9 = subset(DOY.df, scientific_name == "Rhododendron lapponicum")
df.AL.F10 = subset(DOY.df, scientific_name == "Sibbaldiopsis tridentata")
df.AL.F11 = subset(DOY.df, scientific_name == "Vaccinium uliginosum")
df.AL.F12 = subset(DOY.df, scientific_name == "Vaccinium vitis-idaea")

#look at quick relationships between important variables and vif

DOY.df1 = subset(DOY.df, select = c("doy_first", "Tmean_SMOD", "Tmean_spring", "Tmax_spring", "Tmin_spring",
                                    "Tmean_winter", "gdd_0_spring_accum", "gdd_0_SMOD", "Prcp_spring", "Swe_spring",
                                    "SMOD", "SMOD_class", "TPI", "WindShelter", "VRML", "elevation", "slope", "aspect"))
DOY.df1 = na.omit(DOY.df1)

DOY.df1$Tmean_spring = scale(DOY.df1$Tmean_spring)
DOY.df1$SMOD = scale(DOY.df1$SMOD)
DOY.df1$Tmean_winter = scale(DOY.df1$Tmean_winter)
DOY.df1$gdd_0_spring_accum = scale(DOY.df1$gdd_0_spring_accum)
DOY.df1$WindShelter = scale(DOY.df1$WindShelter)
DOY.df1$TPI = scale(DOY.df1$TPI)
DOY.df1$Prcp_spring = scale(DOY.df1$Prcp_spring)
DOY.df1$Swe_spring = scale(DOY.df1$Swe_spring)
DOY.df1$Tmax_spring = scale(DOY.df1$Tmax_spring)
DOY.df1$Tmin_spring = scale(DOY.df1$Tmin_spring)
DOY.df1$Tmean_SMOD = scale(DOY.df1$Tmean_SMOD)
DOY.df1$gdd_0_SMOD = scale(DOY.df1$gdd_0_SMOD)
DOY.df1$SMOD_class = scale(DOY.df1$SMOD_class)
DOY.df1$elevation = scale(DOY.df1$elevation)
DOY.df1$VRML = scale(DOY.df1$VRML)
DOY.df1$aspect = scale(DOY.df1$aspect)
DOY.df1$slope = scale(DOY.df1$slope)

#standardize colum names
colnames(DOY.df1)[colnames(DOY.df1) == 'slope'] <- 'Slope'
colnames(DOY.df1)[colnames(DOY.df1) == 'aspect'] <- 'Aspect'
colnames(DOY.df1)[colnames(DOY.df1) == 'elevation'] <- 'Elevation'
colnames(DOY.df1)[colnames(DOY.df1) == 'doy_first'] <- 'DOY'
colnames(DOY.df1)[colnames(DOY.df1) == 'gdd_0_SMOD'] <- 'GDD_SMOD'
colnames(DOY.df1)[colnames(DOY.df1) == 'gdd_0_spring_accum'] <- 'GDD_S_Accum'

#variable correlation plot
corr <- round(cor(DOY.df1, method = "spearman"), 2)

corr = as.data.frame(corr)

cor_p = ggcorrplot(corr, hc.order = FALSE,
           type = "lower",
           outline.color = "white",
           lab = TRUE,
           legend.title = "Spearman's \u03C1", lab_size = 3) +
           theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=1))
cor_p #export

tiff(file = "cor_p.tiff", width = 8, height = 7, units = 'in', res = 600, pointsize = 11)
cor_p
dev.off()

#quick frequentist look at full model with important variables

lma = lmer(doy_first ~ Tmean_winter + TPI + Tmean_SMOD + SMOD_class + (1|scientific_name) + (1|peak) + (1|source), data = DOY.df)
summary(lma)
plot_model(lma, sort.est = TRUE, title = "DOY Flowering", show.values = TRUE,
           show.p = TRUE, axis.title = "Standardized Regression Coefficient")
vif(lma)
tab_model(lma)

lma101 = lmer(doy_first ~ Tmean_winter + TPI + Tmean_SMOD + SMOD_class + (1|scientific_name) + (1|peak) + (1|source), data = DOY.df)
summary(lma)
vif(lma)
tab_model(lma)

aictab(lma, lma101)

#multiple model exploratory comparisons

lma11 = lmer(doy_first ~ Tmean_SMOD + (1|scientific_name) + (1|peak) + (1|source), data = DOY.df)
summary(lma11)
#plot_model(lma11, sort.est = TRUE, title = "Day of Year Flowering") #export
#vif(lma11)
tab_model(lma11)

lma12 = lmer(doy_first ~ Tmean_winter + (1|scientific_name) + (1|peak) + (1|source), data = DOY.df)
summary(lma12)
#plot_model(lma12, sort.est = TRUE, title = "Day of Year Flowering") #export
#vif(lma12)
tab_model(lma12)

lma13 = lmer(doy_first ~ TPI + (1|scientific_name) + (1|peak) + (1|source), data = DOY.df)
summary(lma13)
#plot_model(lma13, sort.est = TRUE, title = "Day of Year Flowering") #export
#vif(lma13)
tab_model(lma13)

lma14 = lmer(doy_first ~ SMOD_class + (1|scientific_name) + (1|peak) + (1|source), data = DOY.df)
summary(lma14)
#plot_model(lma14, sort.est = TRUE, title = "Day of Year Flowering") #export
#vif(lma14)
tab_model(lma14)

lma2 = lmer(doy_first ~ Tmean_winter + TPI + Tmean_SMOD + SMOD_class + (1|scientific_name) + (1|peak) + (1|source), data = DOY.df)

summary(lma2)

plot_model(lma2, sort.est = TRUE)

vif(lma2)
tab_model(lma2)
###

#Topography

l8 = lmer(doy_first ~ TPI + elevation + aspect + slope + (1|scientific_name) + (1|peak) + (1|source), data = DOY.df)
l9 = lmer(doy_first ~ TPI + elevation + aspect + (1|scientific_name) + (1|peak) + (1|source), data = DOY.df)
l10 = lmer(doy_first ~ TPI + elevation + (1|scientific_name) + (1|peak) + (1|source), data = DOY.df)
l11 = lmer(doy_first ~ TPI + (1|scientific_name) + (1|peak) + (1|source), data = DOY.df)
l12 = lmer(doy_first ~ TPI + elevation + slope + (1|scientific_name) + (1|peak) + (1|source), data = DOY.df)
l13 = lmer(doy_first ~ TPI + slope + (1|scientific_name) + (1|peak) + (1|source), data = DOY.df)
l14 = lmer(doy_first ~ slope + (1|scientific_name) + (1|peak) + (1|source), data = DOY.df)
l15 = lmer(doy_first ~ TPI + aspect + slope + (1|scientific_name) + (1|peak) + (1|source), data = DOY.df)
l16 = lmer(doy_first ~ aspect + slope + (1|scientific_name) + (1|peak) + (1|source), data = DOY.df)
l17 = lmer(doy_first ~ aspect + (1|scientific_name) + (1|peak) + (1|source), data = DOY.df)
l18 = lmer(doy_first ~ elevation + (1|scientific_name) + (1|peak) + (1|source), data = DOY.df)
l19 = lmer(doy_first ~ elevation + aspect + (1|scientific_name) + (1|peak) + (1|source), data = DOY.df)
l20 = lmer(doy_first ~ elevation + slope + (1|scientific_name) + (1|peak) + (1|source), data = DOY.df)
l21 = lmer(doy_first ~ TPI + aspect + (1|scientific_name) + (1|peak) + (1|source), data = DOY.df)
l22 = lmer(doy_first ~ elevation + aspect + slope + (1|scientific_name) + (1|peak) + (1|source), data = DOY.df)

models <- list(l8,l9,l10,l11,l12,l13,l14,l15,l16,l17,l18,l19,l20,l21,l22)

#specify model names
mod.names <- c('l8', 'l9', 'l10', 'l11', 'l12', 'l13', 'l14', 'l15', 'l16', 'l17', 'l18', 'l19', 'l20', 'l21', 'l22')

#calculate AIC of each model
aictab(cand.set = models, modnames = mod.names)
#Climate 

l1 = lmer(doy_first ~ Tmean_SMOD + Swe_spring + Tmean_winter + (1|scientific_name) + (1|peak) + (1|source), data = DOY.df)
l2 = lmer(doy_first ~ Tmean_SMOD + Swe_spring + (1|scientific_name) + (1|peak) + (1|source), data = DOY.df)
l3 = lmer(doy_first ~ Tmean_SMOD + (1|scientific_name) + (1|peak) + (1|source), data = DOY.df)
l4 = lmer(doy_first ~ Tmean_SMOD + Tmean_winter + (1|scientific_name) + (1|peak) + (1|source), data = DOY.df)
l5 = lmer(doy_first ~ Swe_spring + Tmean_winter + (1|scientific_name) + (1|peak) + (1|source), data = DOY.df)
l6 = lmer(doy_first ~ Swe_spring + (1|scientific_name) + (1|peak) + (1|source), data = DOY.df)
l7 = lmer(doy_first ~ Tmean_winter + (1|scientific_name) + (1|peak) + (1|source), data = DOY.df)

#

models <- list(l1,l2,l3,l4,l5,l6,l7)

#specify model names
mod.names <- c('l1', 'l2', 'l3', 'l4', 'l5', 'l6', 'l7')

#calculate AIC of each model
aictab(cand.set = models, modnames = mod.names)

#All

l23 = lmer(doy_first ~ Tmean_SMOD + Tmean_winter + TPI + elevation + SMOD_class + (1|scientific_name) + (1|peak) + (1|source), data = DOY.df)
l24 = lmer(doy_first ~ Tmean_SMOD + Tmean_winter + TPI + SMOD_class + (1|scientific_name) + (1|peak) + (1|source), data = DOY.df)
l25 = lmer(doy_first ~ Tmean_SMOD + Tmean_winter + elevation + SMOD_class + (1|scientific_name) + (1|peak) + (1|source), data = DOY.df)
l26 = lmer(doy_first ~ Tmean_SMOD + TPI + elevation + SMOD_class + (1|scientific_name) + (1|peak) + (1|source), data = DOY.df)
l27 = lmer(doy_first ~ Tmean_SMOD + elevation + SMOD_class + (1|scientific_name) + (1|peak) + (1|source), data = DOY.df)
l28 = lmer(doy_first ~ Tmean_SMOD + TPI + SMOD_class + (1|scientific_name) + (1|peak) + (1|source), data = DOY.df)
l29 = lmer(doy_first ~ Tmean_winter + TPI + elevation + SMOD_class + (1|scientific_name) + (1|peak) + (1|source), data = DOY.df)
l30 = lmer(doy_first ~ Tmean_winter + elevation + SMOD_class + (1|scientific_name) + (1|peak) + (1|source), data = DOY.df)
l31 = lmer(doy_first ~ Tmean_winter + TPI + SMOD_class + (1|scientific_name) + (1|peak) + (1|source), data = DOY.df)

models <- list(l23,l24,l25,l26,l27,l28,l29,l30,l31)

#specify model names
mod.names <- c('l23', 'l24', 'l25', 'l26', 'l27', 'l28', 'l29', 'l30', 'l31')

#calculate AIC of each model
aictab(cand.set = models, modnames = mod.names)

vif(l23)
tab_model(l23)
summary(l23)

vif(l24)
tab_model(l24)
###

###brms - play around with brms fitting of models - WAIC and LOO comparisons

fit1 <- brm(doy_first ~ Tmean_SMOD + TPI + SMOD_class + Tmean_winter
            + (1|scientific_name) + (1|source) + (1|peak),
            data = DOY.df, family = gaussian())
summary(fit1)
bayes_R2(fit1)
plot(conditional_effects(fit1, effects = "TPI:SMOD_class"))

fit111 <- brm(doy_first ~ Tmean_SMOD + TPI * SMOD_class + Tmean_winter
            + (1|scientific_name) + (1|source) + (1|peak),
            data = DOY.df, family = gaussian())
summary(fit111)
bayes_R2(fit111)

fit2 <- brm(doy_first ~ Tmean_SMOD
            + (1|scientific_name) + (1|source) + (1|peak),
            data = DOY.df, family = gaussian())

fit3 <- brm(doy_first ~ Tmean_SMOD + Tmean_winter
            + (1|scientific_name) + (1|source) + (1|peak),
            data = DOY.df, family = gaussian())

fit4 <- brm(doy_first ~ TPI + SMOD_class
            + (1|scientific_name) + (1|source) + (1|peak),
            data = DOY.df, family = gaussian())

fit5 <- brm(doy_first ~ TPI
            + (1|scientific_name) + (1|source) + (1|peak),
            data = DOY.df, family = gaussian())

fit6 <- brm(doy_first ~ SMOD_class
            + (1|scientific_name) + (1|source) + (1|peak),
            data = DOY.df, family = gaussian())

fit7 <- brm(doy_first ~ TPI * SMOD_class
            + (1|scientific_name) + (1|source) + (1|peak),
            data = DOY.df, family = gaussian())

fit8 <- brm(doy_first ~ (1|scientific_name) + (1|source) + (1|peak),
            data = DOY.df, family = gaussian())

fit9 <- brm(doy_first ~ Tmean_winter
            + (1|scientific_name) + (1|source) + (1|peak),
            data = DOY.df, family = gaussian())

fit10 <- brm(doy_first ~ Tmean_SMOD_l + TPI + SMOD_class + Tmean_winter
            + (1|scientific_name) + (1|source) + (1|peak),
            data = DOY.df2, family = gaussian())
summary(fit10)

fit11 <- brm(doy_first ~ Tmean_SMOD + TPI + SMOD_class + Tmean_winter
             + (1|scientific_name) + (1|source) + (1|peak),
             data = DOY.df2, family = gaussian())
summary(fit11)

waic(fit1, fit2, fit3, fit4, fit5, fit6, fit7, fit8, fit9, fit111)
waic(fit10, fit11)
loo(fit10, fit11)

fitS1 <- brm(doy_first ~ Tmean_SMOD + TPI + SMOD_class + Tmean_winter
            + (1|source) + (1|peak),
            data = df.AL.F1, family = gaussian())

fitS2 <- brm(doy_first ~ Tmean_SMOD + TPI + SMOD_class + Tmean_winter
             + (1|source) + (1|peak),
             data = df.AL.F2, family = gaussian())

fitS3 <- brm(doy_first ~ Tmean_SMOD + TPI + SMOD_class + Tmean_winter
             + (1|source) + (1|peak),
             data = df.AL.F3, family = gaussian())

fitS4 <- brm(doy_first ~ Tmean_SMOD + TPI + SMOD_class + Tmean_winter
             + (1|source) + (1|peak),
             data = df.AL.F4, family = gaussian())

fitS5 <- brm(doy_first ~ Tmean_SMOD + TPI + SMOD_class + Tmean_winter
             + (1|source) + (1|peak),
             data = df.AL.F5, family = gaussian())

fitS6 <- brm(doy_first ~ Tmean_SMOD + TPI + SMOD_class + Tmean_winter
             + (1|source) + (1|peak),
             data = df.AL.F6, family = gaussian())

fitS7 <- brm(doy_first ~ Tmean_SMOD + TPI + SMOD_class + Tmean_winter
             + (1|source) + (1|peak),
             data = df.AL.F7, family = gaussian())

fitS8 <- brm(doy_first ~ Tmean_SMOD + TPI + SMOD_class + Tmean_winter
             + (1|source) + (1|peak),
             data = df.AL.F8, family = gaussian())

fitS9 <- brm(doy_first ~ Tmean_SMOD + TPI + SMOD_class + Tmean_winter
             + (1|source) + (1|peak),
             data = df.AL.F9, family = gaussian())

fitS10 <- brm(doy_first ~ Tmean_SMOD + TPI + SMOD_class + Tmean_winter
             + (1|source) + (1|peak),
             data = df.AL.F10, family = gaussian())

fitS11 <- brm(doy_first ~ Tmean_SMOD + TPI + SMOD_class + Tmean_winter
             + (1|source) + (1|peak),
             data = df.AL.F11, family = gaussian())

fitS12 <- brm(doy_first ~ Tmean_SMOD + TPI + SMOD_class + Tmean_winter
             + (1|source) + (1|peak),
             data = df.AL.F12, family = gaussian())

bayes_R2(fitS1)
bayes_R2(fitS2)
bayes_R2(fitS3)
bayes_R2(fitS4)
bayes_R2(fitS5)
bayes_R2(fitS6)
bayes_R2(fitS7)
bayes_R2(fitS8)
bayes_R2(fitS9)
bayes_R2(fitS10)
bayes_R2(fitS11)
bayes_R2(fitS12)

###

fitY0 <- brm(doy_first ~ year_first
            + (1|scientific_name) + (1|source) + (1|peak),
            data = DOY.df, family = gaussian())

fitY1 <- brm(doy_first ~ year_first
             + (1|source) + (1|peak),
             data = df.AL.F1, family = gaussian())

fitY2 <- brm(doy_first ~ year_first
             + (1|source) + (1|peak),
             data = df.AL.F2, family = gaussian())

fitY3 <- brm(doy_first ~ year_first
             + (1|source) + (1|peak),
             data = df.AL.F3, family = gaussian())

fitY4 <- brm(doy_first ~ year_first
             + (1|source) + (1|peak),
             data = df.AL.F4, family = gaussian())

fitY5 <- brm(doy_first ~ year_first
             + (1|source) + (1|peak),
             data = df.AL.F5, family = gaussian())

fitY6 <- brm(doy_first ~ year_first
             + (1|source) + (1|peak),
             data = df.AL.F6, family = gaussian())

fitY7 <- brm(doy_first ~ year_first
             + (1|source) + (1|peak),
             data = df.AL.F7, family = gaussian())

fitY8 <- brm(doy_first ~ year_first
             + (1|source) + (1|peak),
             data = df.AL.F8, family = gaussian())

fitY9 <- brm(doy_first ~ year_first
             + (1|source) + (1|peak),
             data = df.AL.F9, family = gaussian())

fitY10 <- brm(doy_first ~ year_first
             + (1|source) + (1|peak),
             data = df.AL.F10, family = gaussian())

fitY11 <- brm(doy_first ~ year_first
             + (1|source) + (1|peak),
             data = df.AL.F11, family = gaussian())

fitY12 <- brm(doy_first ~ year_first
             + (1|source) + (1|peak),
             data = df.AL.F12, family = gaussian())

summary(fitY0)
bayes_R2(fitY0)
summary(fitY1)
bayes_R2(fitY1)
summary(fitY2)
bayes_R2(fitY2)
summary(fitY3)
bayes_R2(fitY3)
summary(fitY4)
bayes_R2(fitY4)
summary(fitY5)
bayes_R2(fitY5)
summary(fitY6)
bayes_R2(fitY6)
summary(fitY7)
bayes_R2(fitY7)
summary(fitY8)
bayes_R2(fitY8)
summary(fitY9)
bayes_R2(fitY9)
summary(fitY10)
bayes_R2(fitY10)
summary(fitY11)
bayes_R2(fitY11)
summary(fitY12)
bayes_R2(fitY12)

##### JAGS models ##### Main fits

DOY.df = subset(DOY.df, select = c("ids", "scientific_name", "group", "type", "code", "source", "peak", "year_first", "doy_first",
                                   "TPI", "Tmean_SMOD", "SMOD_class", "diff", "Tmean_winter", "Tmean_SMOD_l"))
DOY.df = na.omit(DOY.df)

# Model with common intercept, fixed effects for temperature(s), landscape(s), SMOD, and species random 
# effects.

#Broad models

# Model:
mod <- "model
{
  for(k in 1:Ntot) {
  FLD[k] ~ dnorm(mu[k],tau)
  mu[k] <- b0 + b2*TPI[k] + b3*clim_smod[k] + b4*SMOD_c[k] + b5*clim_w[k] + a1[spnumber[k]] + a2[datan[k]] + a3[peakn[k]]
  } #End k loop
  
  for(l in 1:N.spp){
    a1[l] ~ dnorm(0, 0.001)
  } #End l loop
  
  for(i in 1:dataset){
    a2[i] ~ dnorm(0, 0.001)
  } #End i loop
  
  for(j in 1:peak){
    a3[j] ~ dnorm(0, 0.001)
  } #End j loop
  
  b0 ~ dnorm(0, 0.001)
  b2 ~ dnorm(0, 0.001)
  b3 ~ dnorm(0, 0.001)
  b4 ~ dnorm(0, 0.001)
  b5 ~ dnorm(0, 0.001)
  tau <- sigma^-2
  sigma ~ dunif(0,100)
}"
  write(mod, "model.txt")
  
  params = c("a1", "a2", "a3", "b0", "b2", "b3", "b4", "b5", "tau", "sigma", "mu") #parameters to monitor
  
  # List of observations by country and lifeform
  model.vect <- list(DOY.df)
  model.vect2 <- list("DOY.df")
  
  # Model loop
  
  SPT.model.latents <- data.frame(Species_Type = rep(c('Alpine'), each = 3),
                                  n = NA, DIC = NA, r2 = NA, 
                                  a1.mean = NA, a1.sd = NA, a1.CI2.5 = NA, a1.CI97.5 = NA,
                                  a2.mean = NA, a2.sd = NA, a2.CI2.5 = NA, a2.CI97.5 = NA,
                                  a3.mean = NA, a3.sd = NA, a3.CI2.5 = NA, a3.CI97.5 = NA,
                                  b2.mean = NA, b2.sd = NA, b2.CI2.5 = NA, b2.CI97.5 = NA,
                                  b3.mean = NA, b3.sd = NA, b3.CI2.5 = NA, b3.CI97.5 = NA,
                                  b4.mean = NA, b4.sd = NA, b4.CI2.5 = NA, b4.CI97.5 = NA,
                                  b5.mean = NA, b5.sd = NA, b5.CI2.5 = NA, b5.CI97.5 = NA,
                                  b0.mean = NA, b0.sd = NA, b0.CI2.5 = NA, b0.CI97.5 = NA)
  
  for(i in 1:1)
    print(Sys.time()) # Track run time
    print(i) # Track iterations
    
    df <- model.vect[[i]]
    
    spnumber = as.numeric(factor(DOY.df$scientific_name))
    N.spp = max(spnumber)
    
    Ntot = dim(DOY.df)[1]
    
    SPT.model.latents$n <- Ntot
    
    datan = as.numeric(factor(DOY.df$source))
    dataset = max(datan)
    
    peakn = as.numeric(factor(DOY.df$peak))
    peak = max(peakn)
    
    #input for JAGS
    input = list(Ntot = Ntot, N.spp = N.spp, spnumber = as.numeric(factor(DOY.df$scientific_name)), datan = as.numeric(factor(DOY.df$source)),
                 dataset = max(datan), peakn = as.numeric(factor(DOY.df$peak)), peak = max(peakn), FLD = DOY.df$doy_first, 
                 TPI = as.vector(DOY.df$TPI), clim_smod = as.vector(DOY.df$Tmean_SMOD), SMOD_c = as.vector(DOY.df$SMOD_class),
                 clim_w = as.vector(DOY.df$Tmean_winter)) #input data
    
    mod.FLD.Spring <- jags(model = "model.txt", data = input, param = params,
                           n.chains = 3, #number of separate MCMC chains
                           n.iter = 5000, #number of iterations per chain
                           n.burnin = 1000) #number of initial iterations to discard
    SPT.model.latents$DIC <- DIC(mod.FLD.Spring)
    
    attach.jags(mod.FLD.Spring)
    
    muFLD = colMeans(mu)
    SPT.model.latents$r2 <- summary(lm(muFLD ~ DOY.df$doy_first)$r.squared) 
    
    # Latent posteriors
    SPT.model.latents$b2.mean = mean(b2)
    SPT.model.latents$b2.sd = sd(b2)
    SPT.model.latents$b2.CI2.5 = quantile(b2, probs = 0.025, names = F)
    SPT.model.latents$b2.CI97.5 = quantile(b2, probs = 0.975, names = F)
    
    SPT.model.latents$b3.mean = mean(b3)
    SPT.model.latents$b3.sd = sd(b3)
    SPT.model.latents$b3.CI2.5 = quantile(b3, probs = 0.025, names = F)
    SPT.model.latents$b3.CI97.5 = quantile(b3, probs = 0.975, names = F)
    
    SPT.model.latents$b4.mean = mean(b4)
    SPT.model.latents$b4.sd = sd(b4)
    SPT.model.latents$b4.CI2.5 = quantile(b4, probs = 0.025, names = F)
    SPT.model.latents$b4.CI97.5 = quantile(b4, probs = 0.975, names = F)
    
    SPT.model.latents$b5.mean = mean(b5)
    SPT.model.latents$b5.sd = sd(b5)
    SPT.model.latents$b5.CI2.5 = quantile(b5, probs = 0.025, names = F)
    SPT.model.latents$b5.CI97.5 = quantile(b5, probs = 0.975, names = F)
    
    SPT.model.latents$b0.mean = mean(b0)
    SPT.model.latents$b0.sd = sd(b0)
    SPT.model.latents$b0.CI2.5 = quantile(b0, probs = 0.025, names = F)
    SPT.model.latents$b0.CI97.5 = quantile(b0, probs = 0.975, names = F)
    
    SPT.model.latents$a1.mean = mean(a1)
    SPT.model.latents$a1.sd = sd(a1)
    SPT.model.latents$a1.CI2.5 = quantile(a1, probs = 0.025, names = F)
    SPT.model.latents$a1.CI97.5 = quantile(a1, probs = 0.975, names = F)
    
    SPT.model.latents$a2.mean = mean(a2)
    SPT.model.latents$a2.sd = sd(a2)
    SPT.model.latents$a2.CI2.5 = quantile(a2, probs = 0.025, names = F)
    SPT.model.latents$a2.CI97.5 = quantile(a2, probs = 0.975, names = F)
    
    SPT.model.latents$a3.mean = mean(a3)
    SPT.model.latents$a3.sd = sd(a3)
    SPT.model.latents$a3.CI2.5 = quantile(a3, probs = 0.025, names = F)
    SPT.model.latents$a3.CI97.5 = quantile(a3, probs = 0.975, names = F)
    
    # Species-level hierarchical posteriors
    df.sp <- data.frame(Species_Type = SPT.model.latents$Species_Type, 
                        Sp = 1:N.spp,
                        a1.mean = NA, a1.sd = NA, a1.CI2.5 = NA, a1.CI97.5 = NA) 
    
    for(j in 1:N.spp){
      
      df.sp$a1.mean[j] <- mean(a1[,j])
      df.sp$a1.sd[j] = sd(a1[,j])
      df.sp$a1.CI2.5[j] = quantile(a1[,j], probs = 0.025, names = F)
      df.sp$a1.CI97.5[j] = quantile(a1[,j], probs = 0.975, names = F)
    }
    
    # Dataset-level hierarchical posteriors
    df.set <- data.frame(Species_Type = SPT.model.latents$Species_Type, 
                         Set = 1:dataset,
                         a2.mean = NA, a2.sd = NA, a2.CI2.5 = NA, a2.CI97.5 = NA) 
    
    for(j in 1:dataset){
      
      df.set$a2.mean[j] <- mean(a2[,j])
      df.set$a2.sd[j] = sd(a2[,j])
      df.set$a2.CI2.5[j] = quantile(a2[,j], probs = 0.025, names = F)
      df.set$a2.CI97.5[j] = quantile(a2[,j], probs = 0.975, names = F)
    }
    
    df.peak <- data.frame(Species_Type = SPT.model.latents$Species_Type,
                         P = 1:peak,
                         a3.mean = NA, a3.sd = NA, a3.CI2.5 = NA, a3.CI97.5 = NA) 
    
    for(l in 1:peak){
      
      df.peak$a3.mean[l] <- mean(a3[,l])
      df.peak$a3.sd[l] = sd(a3[,l])
      df.peak$a3.CI2.5[l] = quantile(a3[,l], probs = 0.025, names = F)
      df.peak$a3.CI97.5[l] = quantile(a3[,l], probs = 0.975, names = F)
    }
    
    write.csv(df.sp, paste0(model.vect2[i],".Sp.Post.csv"))
    
    write.csv(df.set, paste0(model.vect2[i],".Set.Post.csv"))
    
    write.csv(df.peak, paste0(model.vect2[i],".Peak.Post.csv"))
    detach.jags()
  
  write.csv(SPT.model.latents, paste0("SpT.model.Latent.Post.csv"))
  
  ###end run here
  
  ##### JAGS models ##### per species
  
  #Models per species if needed # do not run
  
  # Model:
  mod1 <- "model
{
  for(k in 1:Ntot) {
  FLD[k] ~ dnorm(mu[k],tau)
  mu[k] <- b0 + b2*TPI[k] + b3*clim_smod[k] + b4*SMOD_c[k] + b5*clim_w[k] + a2[datan[k]] + a3[peakn[k]]
  } #End k loop
  
  for(i in 1:dataset){
    a2[i] ~ dnorm(0, 0.001)
  } #End i loop
  
  for(j in 1:peak){
    a3[j] ~ dnorm(0, 0.001)
  } #End j loop
  
  b0 ~ dnorm(0, 0.001)
  b2 ~ dnorm(0, 0.001)
  b3 ~ dnorm(0, 0.001)
  b4 ~ dnorm(0, 0.001)
  b5 ~ dnorm(0, 0.001)
  tau <- sigma^-2
  sigma ~ dunif(0,100)
}"
  write(mod1, "model1.txt")
  
  params1 = c("a2", "a3", "b0", "b2", "b3", "b4", "b5", "tau", "sigma", "mu") #parameters to monitor
  
  # List of observations by country and lifeform ##########start here
  model.vect <- list(df.AL.F1, df.AL.F2, df.AL.F3, df.AL.F4, df.AL.F5, df.AL.F6,
                     df.AL.F7, df.AL.F8, df.AL.F9, df.AL.F10, df.AL.F11, df.AL.F12)
  model.vect2 <- list("df.AL.F1", "df.AL.F2", "df.AL.F3", "df.AL.F4", "df.AL.F5", "df.AL.F6",
                      "df.AL.F7", "df.AL.F8", "df.AL.F9", "df.AL.F10", "df.AL.F11", "df.AL.F12")
  
  # Model loop
  
  SPT.model.latents1 <- data.frame(Species_Type = c('CABI', 'CLBO', 'COCA', 'DILA', 'GEPE', 'KAPR',
                                                   'MACA', 'RHGR', 'RHLA', 'SITR', 'VAUL', 'VAVI'),
                                  Growth_Habit = c('alpine', 'alpine', 'alpine', 'alpine', 'alpine', 'alpine',
                                                   'alpine', 'alpine', 'alpine', 'alpine', 'alpine', 'alpine'),
                                  Type = c("herb", "herb", "herb", "shrub", "herb", "shrub", "herb", "shrub",
                                           "shrub", "herb", "shrub", "shrub"),
                                  Form = c("arctic", "trans", "trans", "arctic", "arctic", "arctic", "trans",
                                           "trans", "arctic", "arctic", "arctic", "arctic"),
                                  n = NA, DIC = NA, r2 = NA,
                                  a2.mean = NA, a2.sd = NA, a2.CI2.5 = NA, a2.CI97.5 = NA,
                                  a3.mean = NA, a3.sd = NA, a3.CI2.5 = NA, a3.CI97.5 = NA,
                                  b2.mean = NA, b2.sd = NA, b2.CI2.5 = NA, b2.CI97.5 = NA,
                                  b3.mean = NA, b3.sd = NA, b3.CI2.5 = NA, b3.CI97.5 = NA,
                                  b4.mean = NA, b4.sd = NA, b4.CI2.5 = NA, b4.CI97.5 = NA,
                                  b5.mean = NA, b5.sd = NA, b5.CI2.5 = NA, b5.CI97.5 = NA,
                                  b0.mean = NA, b0.sd = NA, b0.CI2.5 = NA, b0.CI97.5 = NA)
  
  
  for(i in 1:12){
    print(Sys.time()) # Track run time
    print(i) # Track iterations
    
    df <- model.vect[[i]]
    
    spnumber = as.numeric(factor(df$code))
    N.spp = max(spnumber)
    
    Ntot = dim(df)[1]
    
    SPT.model.latents1$n[i] <- Ntot
    
    datan = as.numeric(factor(df$source))
    dataset = max(datan)
    
    peakn = as.numeric(factor(df$peak))
    peak = max(peakn)
    
    #input for JAGS
    input = list(Ntot = Ntot, N.spp = N.spp, spnumber = as.numeric(factor(df$scientific_name)), datan = as.numeric(factor(df$source)),
                 dataset = max(datan), peakn = as.numeric(factor(df$peak)), peak = max(peakn), FLD = df$doy_first, 
                 TPI = as.vector(df$TPI), clim_smod = as.vector(df$Tmean_SMOD), SMOD_c = as.vector(df$SMOD_class),
                 clim_w = as.vector(df$Tmean_winter)) #input data
    
    mod.FLD.Spring <- jags(model = "model1.txt", data = input, param = params1,
                           n.chains = 3, #number of separate MCMC chains
                           n.iter = 5000, #number of iterations per chain
                           n.burnin = 1000) #number of initial iterations to discard
    SPT.model.latents1$DIC[i] <- DIC(mod.FLD.Spring)
    
    attach.jags(mod.FLD.Spring)
    
    muFLD = colMeans(mu)
    SPT.model.latents1$r2[i] <- summary(lm(muFLD ~ df$doy_first)$r.squared) 
    
    # Latent posteriors
    SPT.model.latents1$b2.mean[i] = mean(b2)
    SPT.model.latents1$b2.sd[i] = sd(b2)
    SPT.model.latents1$b2.CI2.5[i] = quantile(b2, probs = 0.025, names = F)
    SPT.model.latents1$b2.CI97.5[i] = quantile(b2, probs = 0.975, names = F)
    
    SPT.model.latents1$b3.mean[i] = mean(b3)
    SPT.model.latents1$b3.sd[i] = sd(b3)
    SPT.model.latents1$b3.CI2.5[i] = quantile(b3, probs = 0.025, names = F)
    SPT.model.latents1$b3.CI97.5[i] = quantile(b3, probs = 0.975, names = F)
    
    SPT.model.latents1$b4.mean[i] = mean(b4)
    SPT.model.latents1$b4.sd[i] = sd(b4)
    SPT.model.latents1$b4.CI2.5[i] = quantile(b4, probs = 0.025, names = F)
    SPT.model.latents1$b4.CI97.5[i] = quantile(b4, probs = 0.975, names = F)
    
    SPT.model.latents1$b5.mean[i] = mean(b5)
    SPT.model.latents1$b5.sd[i] = sd(b5)
    SPT.model.latents1$b5.CI2.5[i] = quantile(b5, probs = 0.025, names = F)
    SPT.model.latents1$b5.CI97.5[i] = quantile(b5, probs = 0.975, names = F)
    
    SPT.model.latents1$b0.mean[i] = mean(b0)
    SPT.model.latents1$b0.sd[i] = sd(b0)
    SPT.model.latents1$b0.CI2.5[i] = quantile(b0, probs = 0.025, names = F)
    SPT.model.latents1$b0.CI97.5[i] = quantile(b0, probs = 0.975, names = F)
    
    SPT.model.latents1$a2.mean[i] = mean(a2)
    SPT.model.latents1$a2.sd[i] = sd(a2)
    SPT.model.latents1$a2.CI2.5[i] = quantile(a2, probs = 0.025, names = F)
    SPT.model.latents1$a2.CI97.5[i] = quantile(a2, probs = 0.975, names = F)
    
    SPT.model.latents1$a3.mean[i] = mean(a3)
    SPT.model.latents1$a3.sd[i] = sd(a3)
    SPT.model.latents1$a3.CI2.5[i] = quantile(a3, probs = 0.025, names = F)
    SPT.model.latents1$a3.CI97.5[i] = quantile(a3, probs = 0.975, names = F)
    
    # Species-level hierarchical posteriors

    # Dataset-level hierarchical posteriors
    df.set <- data.frame(Species_Type = SPT.model.latents1$Species_Type[i], Growth_Habit = SPT.model.latents1$Growth_Habit[i], 
                         Type = SPT.model.latents1$Type[i], Form = SPT.model.latents1$Form[i],
                         Set = 1:dataset,
                         a2.mean = NA, a2.sd = NA, a2.CI2.5 = NA, a2.CI97.5 = NA) 
    
    for(j in 1:dataset){
      
      df.set$a2.mean[j] <- mean(a2[,j])
      df.set$a2.sd[j] = sd(a2[,j])
      df.set$a2.CI2.5[j] = quantile(a2[,j], probs = 0.025, names = F)
      df.set$a2.CI97.5[j] = quantile(a2[,j], probs = 0.975, names = F)
    }
    
    df.peak <- data.frame(Species_Type = SPT.model.latents1$Species_Type[i], Growth_Habit = SPT.model.latents1$Growth_Habit[i], 
                          Type = SPT.model.latents1$Type[i], Form = SPT.model.latents1$Form[i],
                          P = 1:peak,
                          a3.mean = NA, a3.sd = NA, a3.CI2.5 = NA, a3.CI97.5 = NA) 
    
    for(j in 1:peak){
      
      df.peak$a3.mean[j] <- mean(a3[,j])
      df.peak$a3.sd[j] = sd(a3[,j])
      df.peak$a3.CI2.5[j] = quantile(a3[,j], probs = 0.025, names = F)
      df.peak$a3.CI97.5[j] = quantile(a3[,j], probs = 0.975, names = F)
    }
    
    write.csv(df.set, paste0(model.vect2[i],".Sp.Post_sp1.csv"))
    
    write.csv(df.peak, paste0(model.vect2[i],".Sp.Post_sp1.csv"))
    detach.jags()
  }
  
  write.csv(SPT.model.latents1, paste0("SpT.model.Latent.Post_spclass.csv"))
  
####Figures

#Figure 3

fig2_d = read.csv("fig2_data.csv", header = TRUE)
fig2_d$Type = as.factor(fig2_d$Type)
fig2_d$order1 = with(fig2_d, reorder(Species, Spring_Temperature, max))
  
  spt_plot1 = ggplot(data = fig2_d, aes(x = Spring_Temperature, y = order1, color = Type)) +
    geom_errorbar(aes(xmin = Spring_Temperature - ST_95, xmax = Spring_Temperature + ST_95), width = 0, position = position_nudge(x = 0, y = 0.3)) +
    geom_point(aes(colour = Type), size = 4, position = position_nudge(x = 0, y = 0.3)) +
    geom_errorbar(aes(xmin = Winter_Temperature - WT_95, xmax = Winter_Temperature + WT_95), width = 0, position = position_nudge(x = 0, y = -0.3)) +
    geom_point(aes(x = Winter_Temperature, y = order1, colour = Type), pch = 17, size = 4, position = position_nudge(x = 0, y = -0.3)) +
    #geom_point(aes(x = Mean_Net_Temperature_Sensitivity, y = order1), color = "red", pch = 18, size = 5) +
    geom_vline(xintercept = 0, linetype = "dotted") +
    scale_color_manual(values = c("#AFB83B", "#800080", "#F0C807", "#1656AD", "#FF0000"), name = "Functional Type") +
    ylab("") +
    xlab("Sensitivity (days/\u00b0C)") +
    ggtitle("Mean Temperature") +
    annotate('rect', xmin = -Inf, xmax = Inf, 
             ymin=c(0.5, 2.5, 4.5, 6.5, 8.5, 10.5, 12.5), 
             ymax=c(1.5, 3.5, 5.5, 7.5, 9.5, 11.5, 13.5), alpha=0.2, fill="grey50") +
    theme_classic() +
    #facet_grid(rows = vars(Type), cols = vars(parameter_type), scales = "free") +
    theme(axis.text.y = element_text(size = 14),
          plot.title = element_text(size = 16),
          axis.title = element_text(size = 16),
          legend.title = element_text(size = 16),
          legend.text = element_text(size = 14),
          axis.text.x = element_text(size = 14),
          legend.position="none",
          strip.placement = "outside")
  #strip.background = element_blank(),
  #strip.text.x = element_blank())
  spt_plot1 
  
  spt_plot5 = ggplot(data = fig2_d, aes(x = Mean_Net_Temperature_Sensitivity, y = order1)) +
    #geom_errorbar(aes(xmin = Spring_Temperature - ST_95, xmax = Spring_Temperature + ST_95), width = 0, position = position_nudge(x = 0, y = 0.3)) +
    geom_point(color = "black", fill = "red", pch = 23, size = 4) +
    #geom_errorbar(aes(xmin = Winter_Temperature - WT_95, xmax = Winter_Temperature + WT_95), width = 0, position = position_nudge(x = 0, y = -0.3)) +
    #geom_point(aes(x = Winter_Temperature, y = order1, colour = Type), pch = 17, size = 4, position = position_nudge(x = 0, y = -0.3)) +
    #geom_point(aes(x = Mean_Net_Temperature_Sensitivity, y = order1), color = "red", pch = 18, size = 5) +
    geom_vline(xintercept = 0, linetype = "dotted") +
    #scale_color_viridis(name = "Type", discrete = TRUE) +
    ylab("") +
    xlab("Net Sensitivity (days/\u00b0C)") +
    ggtitle("Spring - Winter Mean Temperature") +
    annotate('rect', xmin = -Inf, xmax = Inf, 
             ymin=c(0.5, 2.5, 4.5, 6.5, 8.5, 10.5, 12.5), 
             ymax=c(1.5, 3.5, 5.5, 7.5, 9.5, 11.5, 13.5), alpha=0.2, fill="grey50") +
    theme_classic() +
    #facet_grid(rows = vars(Type), cols = vars(parameter_type), scales = "free") +
    theme(axis.text.y = element_text(size = 14),
          plot.title = element_text(size = 16),
          axis.title = element_text(size = 16),
          legend.title = element_text(size = 16),
          legend.text = element_text(size = 14),
          axis.text.x = element_text(size = 14),
          legend.position="none",
          strip.placement = "outside")
  #strip.background = element_blank(),
  #strip.text.x = element_blank())
  spt_plot5 
  
  fig2_d$order2 = with(fig2_d, reorder(Species, TPI, max))
  
  spt_plot2 = ggplot(data = fig2_d, aes(x = TPI, y = order2, color = Type)) +
    geom_errorbar(aes(xmin = TPI - TPI_95, xmax = TPI + TPI_95), width = 0, position = "dodge") +
    geom_point(aes(colour = Type), size = 4, position = position_dodge(width = 0.3)) +
    geom_vline(xintercept = 0, linetype = "dotted") +
    scale_color_manual(values = c("#AFB83B", "#800080", "#F0C807", "#1656AD", "#FF0000"), name = "Functional Type") +
    xlab("Sensitivity (days/TPI Unit)") +
    ylab("") +
    ggtitle("Topographic Position Index") +
    annotate('rect', xmin = -Inf, xmax = Inf, 
             ymin=c(0.5, 2.5, 4.5, 6.5, 8.5, 10.5, 12.5), 
             ymax=c(1.5, 3.5, 5.5, 7.5, 9.5, 11.5, 13.5), alpha=0.2, fill="grey50") +
    theme_classic() +
    #facet_grid(cols = vars(Type), scales = "free") +
    theme(axis.text.x = element_text(size = 14),
          plot.title = element_text(size = 16),
          axis.title = element_text(size = 16),
          legend.title = element_text(size = 16),
          legend.text = element_text(size = 14),
          axis.text.y = element_text(size = 14),
          legend.position="bottom",
          strip.placement = "outside",
  legend.justification = "left")
  #strip.background = element_blank(),
  #strip.text.x = element_blank())
  spt_plot2 
  
  fig2_d$order3 = with(fig2_d, reorder(Species, SMOD_Class, max))
  
  spt_plot3 = ggplot(data = fig2_d, aes(x = SMOD_Class, y = order3, color = Type)) +
    geom_errorbar(aes(xmin = SMOD_Class - SMOD_95, xmax = SMOD_Class + SMOD_95), width = 0, position = "dodge") +
    geom_point(aes(colour = Type), size = 4, position = position_dodge(width = 0.3)) +
    geom_vline(xintercept = 0, linetype = "dotted") +
    scale_color_manual(values = c("#AFB83B", "#800080", "#F0C807", "#1656AD", "#FF0000"), name = "Functional Type") +
    xlab("Sensitivity (days/SMOD)") +
    ylab("") +
    ggtitle("Snowmelt Date") +
    annotate('rect', xmin = -Inf, xmax = Inf, 
             ymin=c(0.5, 2.5, 4.5, 6.5, 8.5, 10.5, 12.5), 
             ymax=c(1.5, 3.5, 5.5, 7.5, 9.5, 11.5, 13.5), alpha=0.2, fill="grey50") +
    theme_classic() +
    #facet_grid(cols = vars(Type), scales = "free") +
    theme(axis.text.x = element_text(size = 14),
          plot.title = element_text(size = 16),
          axis.title = element_text(size = 16),
          legend.title = element_text(size = 16),
          legend.text = element_text(size = 14),
          axis.text.y = element_text(size = 14),
          legend.position="none")
  #strip.background = element_blank(),
  #strip.text.x = element_blank())
  spt_plot3 

  spt_p = (spt_plot1 + spt_plot5) / (spt_plot2 + spt_plot3) + plot_annotation(tag_levels = "A")
  spt_p

  tiff(file = "coeffs.tiff", width = 12, height = 11, units = 'in', res = 600, pointsize = 11)
  spt_p
  dev.off()
  
###More exploratory data visalization - main and SI figures
  
data.all3 = read.csv("data_all.csv", header = TRUE)  
id = ddply(data.all3, ~scientific_name, summarise, sum = length(scientific_name), med = median(doy_first))
  unique(data.all3$scientific_name)
ids = ddply(data.all3, ~scientific_name, summarise, sum = length(scientific_name), sum1 = length(unique(peak)), med = median(doy_first))
  unique(data.all3$scientific_name)
ids1 = ddply(data.all3, ~peak, summarise, sum = length(scientific_name), med = median(doy_first),
             sum1 = length(unique(scientific_name)), smod1 = min(SMOD, na.rm = TRUE), smod2 = max(SMOD, na.rm = TRUE))

hist(data.all3$doy_first, xlab = "DOY Flowering", main = "(A) Histogram of All Observation DOY")
group_medians <- data.all3 %>%
  group_by(scientific_name) %>%
  summarize(median_value = median(doy_first), .groups = 'drop')

p2 <- ggplot(data=data.all3, aes(x=doy_first, group=scientific_name, fill=scientific_name)) +
  geom_density(adjust=1.5, alpha=0.5) +
  xlab("Day of Year Flower (DOY)") +
  ylab("Density") +
  scale_fill_manual(values = c("#F0C807","#FF0000", "#FF0000",
                              "#800080", "#F0C807", "#800080", "#FF0000",
                              "#FF0000", "#1656AD", "#F0C807",
                              "#1656AD", "#1656AD"), name = "Functional Type") +
  facet_wrap(vars(scientific_name), scales = "fixed", ncol = 4) +
  theme_bw() +
  geom_vline(data = group_medians, aes(xintercept = median_value),
             color = "red",
             linetype = "dashed",
             size = 1) +
  theme(legend.position = "none")
p2

tiff(file = "med1.tiff", width = 9, height = 7, units = 'in', res = 600, pointsize = 11)
p2
dev.off()

hist(data.all3$year_first,
     xlab = "Year", main = "(B) Histogram of Annual Flowering Observations")#export figure

hist(DOY.df$elevation)
a_s_plot3 = ggplot(DOY.df, aes(x = elevation, y = doy_first, color = function1)) +
  geom_point(size = 2, alpha = 0.3) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  scale_color_manual(values = c("#800080", "#F0C807", "#1656AD", "#FF0000"), name = "Functional Type") +
  xlab("Elevation (m)") +
  ylab("DOY Flower") +
  theme_bw() +
  facet_wrap(vars(scientific_name), scales = "fixed", ncol = 4) +
  theme(legend.position = "bottom")
a_s_plot3

tiff(file = "spp_elev.tiff", width = 8, height = 5, units = 'in', res = 600, pointsize = 11)
a_s_plot3
dev.off()

#Figure 6

time_plot = ggplot(DOY.df, aes(x = year_first, y = doy_first, color = function1)) +
  geom_point(size = 2.5, alpha = 0.3, color = "grey20") +
  #scale_x_continuous(labels = label_number(accuracy = 1)) +
  geom_smooth(method = "lm", se = FALSE, lwd = 1.75) +
  scale_color_manual(values = c("#800080", "#F0C807", "#1656AD", "#FF0000"), name = "Functional Type") +
  ylab("Day of Year Flower (DOY)") +
  xlab("Year") +
  theme_bw() +
  #facet_wrap(vars(scientific_name), scales = "free", ncol = 4) +
  theme(strip.text.x = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 11),
        axis.text.x = element_text(size = 14),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        legend.position = "none")
time_plot

#checking normality of species observation data
hist(df.AL.F1$doy_first1)#CABI
df.AL.F1$doy_first1 = log10(df.AL.F1$doy_first)
shapiro.test(df.AL.F1$doy_first1)
ggdensity(df.AL.F1, x = "doy_first", fill = "lightgray", title = "doy_first1") +
  #scale_x_continuous(limits = c(3, 12)) +
  stat_overlay_normal_density(color = "red", linetype = "dashed")

hist(df.AL.F2$doy_first)#CLBO
shapiro.test(df.AL.F2$doy_first)
ggdensity(df.AL.F2, x = "doy_first", fill = "lightgray", title = "doy_first") +
  #scale_x_continuous(limits = c(3, 12)) +
  stat_overlay_normal_density(color = "red", linetype = "dashed")
hist(df.AL.F3$doy_first)#COCA
ggdensity(df.AL.F3, x = "doy_first", fill = "lightgray", title = "doy_first1") +
  #scale_x_continuous(limits = c(3, 12)) +
  stat_overlay_normal_density(color = "red", linetype = "dashed")
shapiro.test(df.AL.F3$doy_first)
hist(df.AL.F4$doy_first)#DILA
shapiro.test(df.AL.F4$doy_first)
hist(df.AL.F5$doy_first)#GEPE
shapiro.test(df.AL.F5$doy_first)
hist(df.AL.F6$doy_first)#KAPR
shapiro.test(df.AL.F6$doy_first)
hist(df.AL.F7$doy_first)#MACA
shapiro.test(df.AL.F7$doy_first)
hist(df.AL.F8$doy_first)#RHGR
shapiro.test(df.AL.F8$doy_first)
hist(df.AL.F9$doy_first)#RHLA
shapiro.test(df.AL.F9$doy_first)
hist(df.AL.F10$doy_first)#SITR
shapiro.test(df.AL.F10$doy_first)
hist(df.AL.F11$doy_first)#VAUL
shapiro.test(df.AL.F11$doy_first)
hist(df.AL.F12$doy_first)#VAVI
shapiro.test(df.AL.F12$doy_first)

supp_den = ggdensity(DOY.df, x = "doy_first", fill = "lightgray", title = "") +
  #scale_x_continuous(limits = c(3, 12)) +
  stat_overlay_normal_density(color = "red", linetype = "dashed") +
  xlab("DOY Flower") +
  ylab("Density") +
  facet_wrap(vars(scientific_name), scales = "free", ncol = 4)
supp_den

tiff(file = "supp_den.tiff", width = 10, height = 6, units = 'in', res = 600, pointsize = 11)
supp_den
dev.off()

time_plot3 = ggplot(DOY.df, aes(x = year_first, y = doy_first, color = function1)) +
  geom_point(size = 2, alpha = 0.3, color = "grey20") +
  #scale_x_continuous(labels = label_number(accuracy = 1)) +
  geom_smooth(method = "lm", se = FALSE, lwd = 1.75) +
  scale_color_manual(values = c("#800080", "#F0C807", "#1656AD", "#FF0000"), name = "Functional Type") +
  ylab("Day of Year Flower (DOY)") +
  xlab("Year") +
  facet_wrap(vars(name), scales = "fixed", ncol = 4) +
  theme_bw() +
  theme(strip.text.x = element_text(size = 11),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.title = element_text(size = 22),
        legend.text = element_text(size = 22),
        axis.text.y = element_text(size = 14),
        axis.text.x = element_text(size = 14, angle = 45, vjust = 1, hjust = 1),
        legend.position = "bottom",
        legend.justification = "left",
        legend.spacing.x = unit(3.0, 'in'))
time_plot3

lm2 = lm(doy_m ~ year_first, data = data.all4)
summary(lm2)

lm3 = lm(doy_first ~ year_first, data = DOY.df)
summary(lm3)

time_plot_t = ggplot(DOY.df, aes(x = year_first, y = Tmean_SMOD)) +
  geom_point() +
  #geom_line() +
  geom_smooth(method = "lm", se = TRUE) +
  ylab("Tmean Spring (\u00b0C)")
  #xlab("Year")
  #facet_wrap(vars(scientific_name), scales = "free", ncol = 4)
time_plot_t

time_plot_w = ggplot() +
  geom_point(data = DOY.df, aes(x = year_first, y = Tmean_winter), color = "navy", size = 3, alpha = 0.2, pch = 1) +
  geom_point(data = DOY.df, aes(x = year_first, y = Tmean_SMOD), color = "red4", size = 3, alpha = 0.2, pch = 1) +
  #geom_line() +
  geom_smooth(data = DOY.df, aes(x = year_first, y = Tmean_winter), method = "lm", se = TRUE, color = "navy", lwd = 1.75) +
  geom_smooth(data = DOY.df, aes(x = year_first, y = Tmean_SMOD), method = "lm", se = TRUE, color = "red4", lwd = 1.75) +
  ylab("Tmean (\u00b0C)") +
  xlab("Year") +
  theme_bw() +
  theme(axis.text = element_text(size = 14),
        axis.title.y = element_text(size = 14)) +
  geom_label(aes(x = 2010, y = 16, label = 'Spring Slope = 0.053, R\u00B2 = 0.54'), parse = FALSE,
           color = "black", size = 5, fontface = "bold") +
  geom_label(aes(x = 2010, y = -3, label = 'Winter Slope = 0.082, R\u00B2 = 0.56'), parse = FALSE,
           color = "black", size = 5, fontface = "bold") +
  theme(axis.title.x = element_text(size = 14),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1))
#facet_wrap(vars(scientific_name), scales = "free", ncol = 4)
time_plot_w

yc = read.csv("year_change.csv", header = TRUE)
yc$order2 = with(yc, reorder(Species, Slope, min))

time_plot2 = ggplot(data = yc, aes(x = Slope, y = order2, color = Functional_Type)) +
  geom_errorbar(aes(xmin = CI2.5, xmax = CI97.5), width = 0, position = "dodge") +
  geom_point(aes(colour = Functional_Type), size = 4, position = position_dodge(width = 0.3)) +
  geom_vline(xintercept = 0, linetype = "dotted") +
  scale_color_manual(values = c("#800080", "#F0C807", "#1656AD", "#FF0000"), name = "Functional Type") +
  xlab("\u03B2*Year (95% CI)") +
  ylab("") +
  ggtitle("") +
  annotate('rect', xmin = -Inf, xmax = Inf, 
           ymin=c(0.5, 2.5, 4.5, 6.5, 8.5, 10.5), 
           ymax=c(1.5, 3.5, 5.5, 7.5, 9.5, 11.5), alpha=0.2, fill="grey50") +
  theme_classic() +
  #facet_grid(cols = vars(Type), scales = "free") +
  theme(axis.text.x = element_text(size = 14),
        plot.title = element_text(size = 16),
        axis.title = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        legend.position="none",
        strip.placement = "outside") +
  guides(color=guide_legend(nrow=2, byrow=TRUE)) 
#strip.background = element_blank(),
#strip.text.x = element_blank())
time_plot2 

t_p = time_plot_w + free(time_plot) + time_plot3 + time_plot2 +
         plot_annotation(tag_levels = 'A') + plot_layout(ncol = 2) +
  plot_layout(widths = c(1.75, 1))
t_p

tiff(file = "time_m.tiff", width = 12.5, height = 10, units = 'in', res = 600, pointsize = 11)
t_p
dev.off()

#Visualizing SMOD data
DOY.df1 = subset(DOY.df, SMOD_class > 3)
DOY.df1$SMOD_class = as.factor(DOY.df1$SMOD_class)
#DOY.df1 = na.omit(DOY.df)

#Figure 5

smod_plot4 = ggplot(DOY.df1, aes(x = SMOD_class, y = diff, fill = function1)) +
  geom_violin(drop = FALSE, alpha = 0.7) +
  stat_summary(fun.y=median, geom="point", pch = 19, size=2, color = "black") +
  #geom_smooth(method = "lm", se = FALSE) +
  ylab("Flowering After SMOD (days)") +
  xlab("SMOD Class") +
  scale_fill_manual(values = c("#800080", "#F0C807", "#1656AD", "#FF0000"), name = "Functional Type") +
  facet_wrap(vars(scientific_name), scales = "fixed", ncol = 4) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black")) +
  theme(strip.text.x = element_text(size = 11),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 11),
        legend.position = "bottom")
smod_plot4

tiff(file = "SMOD_diff.tiff", width = 9.5, height = 6, units = 'in', res = 600, pointsize = 11)
smod_plot4
dev.off()

###Maps

mtn_points <- st_read("mtns.shp")
mtn_points$value = c(9, 8, 7, 6, 4, 5, 3, 2, 1, 13, 12, 11, 10)
alp_area <- st_read("alp_bound.shp")
NE_alp <- st_read("NE_alpine.shp")
na_dem <- rast("C:/Users/jtourville/OneDrive - Appalachian Mountain Club/Desktop/Alp_Phe_Man/WMNF_comm_sci/na_clip.tif")#change to directory - spatial data for maps available upon request
na_dem = as.numeric(na_dem)

NNA <- ne_countries(scale = 50, continent = "North America")
pal <- wes_palette("Zissou1", 100, type = "continuous")

#Figure 1

medium_scale_map <- ggplot() +
  geom_spatraster(data = na_dem, maxcell = Inf) +
  scale_fill_gradientn(colours = pal, na.value = "white") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(fill = "Elevation (m)") +
  geom_sf(data = NNA, alpha = 0.1, lwd = 1.1, color = "black") +
  geom_sf_label(data = NNA,
    aes(label = brk_name),
    color = "black",
    size = 3) +
  geom_sf(data = mtn_points, color = "red", fill = NA, pch = 24, size = 4, stroke = 2) +
  geom_text_repel(data = mtn_points,
                    aes(label = value, geometry = geometry),
                    stat = "sf_coordinates", min.segment.length = 0.2, box.padding = 0.5, max.overlaps = 13, size = 5) +
  annotation_scale(location = "br") + # br = bottom right
  annotation_north_arrow(style = north_arrow_nautical, location = "tr") +
  annotate("text", x = -68.8, y = 46.7, label = "United States", 
           color = "black", size = 3, fontface = "bold") +
  annotate("text", x = -68, y = 47.7, label = "Canada", 
           color = "black", size = 3, fontface = "bold") +
  coord_sf(xlim = c(-66, -77), ylim = c(42.5, 48)) +
  ggtitle("") +
  xlab("Longitude") +
  ylab("Latitude") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 12),
axis.title = element_text(size = 14),
legend.title = element_text(size = 14),
legend.text = element_text(size = 12),
legend.position = "bottom",
legend.spacing.x = unit(3, 'cm'),
axis.text.y = element_text(size = 12),
panel.border = element_rect(colour = "black", fill=NA, linewidth=1)) +
  guides(fill = guide_colorbar(frame.colour = "black",
                               barwidth = 15,
                               barheight = 1.5))
medium_scale_map

#img1 <- fig("C:/Users/jtourville/Downloads/Picture1.jpg")#to read in images for map figure
#img2 <- fig("C:/Users/jtourville/Downloads/Picture2.jpg")
#img3 <- fig("C:/Users/jtourville/Downloads/Picture3.jpg")
#img4 <- fig("C:/Users/jtourville/Downloads/Picture4.jpg")

map1 = medium_scale_map + (img1 + img2 + img3 + img4) + plot_annotation(tag_levels = "A") +
  plot_layout(ncol = 2, nrow = 1)
map1

tiff(file = "map1.tiff", width = 11, height = 5.5, units = 'in', res = 600, pointsize = 11)
map1
dev.off()

#SI Map

NE_scale_map <- ggplot() +
  geom_sf(data = NNA, alpha = 0.9, lwd = 1.5) +
  geom_sf_label(data = NNA,
                aes(label = brk_name),
                color = "black",
                size = 3) +
  geom_sf(data = NE_alp, alpha = 0.9, lwd = 2, color = "red") +
  annotation_scale(location = "br") + # br = bottom right
  annotation_north_arrow(style = north_arrow_nautical, location = "tr") +
  annotate("text", x = -68.8, y = 46.7, label = "United States", 
           color = "black", size = 3.5, fontface = "bold") +
  annotate("text", x = -68, y = 47.7, label = "Canada", 
           color = "black", size = 3.5, fontface = "bold") +
  coord_sf(xlim = c(-64, -76), ylim = c(43, 49.5)) +
  ggtitle("") +
  xlab("Longitude") +
  ylab("Latitude") +
  theme(axis.text.x = element_text(size = 12),
        axis.title = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1))
NE_scale_map

tiff(file = "map_supp.tiff", width = 8.5, height = 7, units = 'in', res = 600, pointsize = 11)
NE_scale_map
dev.off()

#Figure 4 SMOD Map

pres_smod <- rast("C:/Users/jtourville/OneDrive - Appalachian Mountain Club/Desktop/Alp_Phe_Man/WMNF_comm_sci/planet_smod_2017_2024_mean.tif")
melt_area <- st_read("SMOD_early.shp")
pres_p <- st_read("pres_points.shp")
pres_tpi <- rast("C:/Users/jtourville/OneDrive - Appalachian Mountain Club/Desktop/Alp_Phe_Man/WMNF_comm_sci/pres_tpi.tif")
pres_p$SMOD = as.numeric(pres_p$SMOD)
pres_p$diff = as.numeric(pres_p$diff)
pres_p$doy = pres_p$SMOD + pres_p$diff

pres_map <- ggplot() +
  geom_spatraster(data = pres_smod) +
  scale_fill_viridis(option = "C", direction = -1, na.value = "grey90", limits = c(100,200)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  geom_sf(data = melt_area, lwd = 0.5, fill = NA, color = "black") +
  labs(fill = "SMOD") +
  ggtitle("") +
  xlab("Longitude") +
  ylab("Latitude") +
  theme(axis.text.x = element_text(size = 12),
        axis.title = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.position = c(0.1,0.85),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1))
pres_map

tiff(file = "map_smod.tiff", width = 6, height = 6.5, units = 'in', res = 600, pointsize = 11)
pres_map
dev.off()

###END