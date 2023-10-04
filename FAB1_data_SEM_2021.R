library(lavaan)
library(tidySEM)
library(piecewiseSEM)
library(ggplot2)
library(lme4)
library(lmerTest)

C_agg<-read.csv("ProcessedData/Cseq.csv")

## drop rows with missing data
C_agg<-C_agg[-which(is.na(C_agg$soilC) | is.na(C_agg$macro250)),]

## z-standardize important variables
C_agg_standard<-C_agg
standard_cols<-c("species_richness","woodyC","soilC",
                 "macro250","percentAM","percentCon")
C_agg_standard[,standard_cols]<-scale(C_agg_standard[,standard_cols])

##################################
## local estimation of SEMs

full_model_lavaan<-'
woodyC~1+species_richness+percentCon
macro250~1+woodyC+percentAM+percentCon
soilC~1+species_richness+percentAM+woodyC+percentCon+macro250
'

fit_full_model <- sem(full_model_lavaan,
                      data=C_agg_standard)
fitMeasures(fit_full_model)

sub_model_lavaan<-'
woodyC~1+species_richness+percentCon
macro250~1+woodyC+percentAM+percentCon
soilC~1+species_richness+percentAM
'

fit_sub_model <- sem(sub_model_lavaan,
                     data=C_agg_standard)
fitMeasures(fit_sub_model)
