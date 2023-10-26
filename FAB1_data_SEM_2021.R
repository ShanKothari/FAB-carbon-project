setwd("C:/Users/querc/Dropbox/FABCarbonProject/")
library(lavaan)
library(tidySEM)
library(piecewiseSEM)
library(ggplot2)
library(lme4)
library(lmerTest)

C_agg<-read.csv("ProcessedData/Cseq.csv")

## drop rows with missing data
C_agg_sub<-C_agg[-which(is.na(C_agg$soilC) | is.na(C_agg$macro250)),]

## z-standardize important variables
C_agg_standard<-C_agg_sub
standard_cols<-c("species_richness","woodyC","soilC",
                 "macro250","percentAM","percentCon")
C_agg_standard[,standard_cols]<-scale(C_agg_standard[,standard_cols])

##################################
## global estimation of SEMs

full_model_lavaan<-'
woodyC~1+species_richness+percentAM+percentCon
macro250~1+percentAM+percentCon
soilC~1+species_richness+percentAM+woodyC+macro250
'

fit_full_model <- sem(full_model_lavaan,
                      data=C_agg_standard)
# fitMeasures(fit_full_model)

## factors that influence leaf litter quantity (woody C)
## or quality (% needleleaf) don't seem to influence soil C
sub_model_lavaan<-'
woodyC~1+species_richness+percentAM+percentCon
macro250~1+percentAM+percentCon
soilC~1+species_richness+percentAM
'

fit_sub_model <- sem(sub_model_lavaan,
                     data=C_agg_standard)
# fitMeasures(fit_sub_model)

min_model_lavaan<-'
woodyC~1+species_richness+percentAM+percentCon
macro250~1+percentCon
soilC~1+species_richness+percentAM
'

fit_min_model <- sem(min_model_lavaan,
                     data=C_agg_standard)
# fitMeasures(fit_min_model)

## note: maybe we should cite my spectra paper to show that conifer
## litter is more recalcitrant...

###################################
## local estimation?

localfit_full_model<-psem(
  lmer(woodyC~species_richness+percentCon+percentAM+(1|block),data=C_agg_standard),
  lmer(macro250~percentAM+percentCon+(1|block),data=C_agg_standard),
  lm(soilC~species_richness+percentAM+woodyC+macro250,
     data=C_agg_standard)
)

## checking individual components
# check_model(lmer(woodyC~species_richness+percentCon+percentAM+(1|block),data=C_agg_standard))
# check_model(lmer(macro250~percentAM+percentCon+(1|block),data=C_agg_standard))
# check_model(lm(soilC~species_richness+percentAM+woodyC+macro250,
#                data=C_agg_standard))

localfit_sub_model<-psem(
  lmer(woodyC~species_richness+percentCon+percentAM+(1|block),data=C_agg_standard),
  lmer(macro250~percentAM+percentCon+(1|block),data=C_agg_standard),
  lm(soilC~species_richness+percentAM+macro250,
     data=C_agg_standard)
)

localfit_min_model<-psem(
  lmer(woodyC~species_richness+percentCon+percentAM+(1|block),data=C_agg_standard),
  lmer(macro250~percentCon+(1|block),data=C_agg_standard),
  lm(soilC~species_richness+percentAM,
     data=C_agg_standard)
)
