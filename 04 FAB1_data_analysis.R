setwd("C:/Users/querc/Dropbox/FABCarbonProject/")

library(lme4)
library(lmerTest)
library(MuMIn)
library(car)
library(performance)
library(piecewiseSEM)
library(MASS)

carbon_seq<-read.csv("ProcessedData/Cseq.csv")

#############################
## AIC-based model selection for predictors of C stocks and OY

####
## total woody C
mwt<-lm(woodyC~species_richness+FDis+PSV,data=carbon_seq)
mwt_block<-lmer(woodyC~species_richness+FDis+PSV+(1|block),
                data=carbon_seq, REML=F)
mwt_robust<-rlm(woodyC~species_richness+FDis+PSV,data=carbon_seq)

options(na.action = "na.fail") # required for dredge to run
mwt_dredge <- dredge(mwt,beta = "none",evaluate = T, rank = AICc)
mwt_block_dredge <- dredge(mwt_block,beta = "none",evaluate = T,rank = AICc)
mwt_robust_dredge <- dredge(mwt_robust,beta = "none",evaluate = T,rank = AICc)
#summary(model.avg(mwt_dredge, subset = delta <= 2))
options(na.action = "na.omit")

mwt_model<-get.models(mwt_dredge, subset = 2)[[1]]
## high heteroskedasticity, non-normality of variance
## but robust regression yields similar results
## check_model(mwt_model)

####
## woody C overyielding
carbon_seq_woodyOY<-carbon_seq[which(!is.na(carbon_seq$woodyOY)),]
mw<-lm(woodyOY~species_richness+FDis+PSV,
       data=carbon_seq_woodyOY)
mw_block<-lmer(woodyOY~species_richness+FDis+PSV+(1|block),
               data=carbon_seq_woodyOY,REML=F)
mw_robust<-rlm(woodyOY~species_richness+FDis+PSV,data=carbon_seq_woodyOY)

options(na.action = "na.fail")
mw_dredge <- dredge(mw, beta = "none", evaluate = T, rank = AICc)
mw_block_dredge <- dredge(mw_block, beta = "none", evaluate = T, rank = AICc)
mw_robust_dredge <- dredge(mw_robust, beta = "none", evaluate = T, rank = AICc)
options(na.action = "na.omit")

mw_model<-get.models(mw_dredge, subset = 1)[[1]]
check_model(mw_model)

####
## complementarity effects
mce<-lm(woodyCE~species_richness+FDis+PSV,data=carbon_seq_woodyOY)
mce_block<-lmer(woodyCE~species_richness+FDis+PSV+(1|block),
                data=carbon_seq_woodyOY, REML=F)
mce_robust<-rlm(woodyCE~species_richness+FDis+PSV,data=carbon_seq_woodyOY)

options(na.action = "na.fail")
mce_dredge <- dredge(mce, beta = "none", evaluate = T, rank = AICc)
mce_block_dredge <- dredge(mce_block, beta = "none", evaluate = T, rank = AICc)
mce_robust_dredge <- dredge(mce_robust, beta = "none", evaluate = T, rank = AICc)
options(na.action = "na.omit")

mce_model<-get.models(mce_dredge, subset = 2)[[1]]
check_model(mce_model)

## because of persistent disputes about the meaningfulness
## of the magnitudes of CE and SE, we might just want
## to see whether they differ from 0
## the distribution is skewed so we use a non-parametric test

wilcox.test(carbon_seq_woodyOY$woodyCE)

####
## selection effects
mse<-lm(woodySE~species_richness+FDis+PSV,data=carbon_seq_woodyOY)
mse_block<-lmer(woodySE~species_richness+FDis+PSV+(1|block),
                data=carbon_seq_woodyOY,REML=F)
mse_robust<-rlm(woodySE~species_richness+FDis+PSV,data=carbon_seq_woodyOY)

options(na.action = "na.fail")
mse_dredge <- dredge(mse, beta = "none", evaluate = T, rank = AICc)
mse_block_dredge <- dredge(mse_block, beta = "none", evaluate = T, rank = AICc)
mse_robust_dredge <- dredge(mse_robust, beta = "none", evaluate = T, rank = AICc)
options(na.action = "na.omit")

mse_model<-get.models(mse_dredge, subset = 1)[[1]]
check_model(mse_model)

## because of persistent disputes about the meaningfulness
## of the magnitudes of CE and SE, we might just want
## to see whether they differ from 0
## the distribution is skewed so we use a non-parametric test

wilcox.test(carbon_seq_woodyOY$woodySE)

####
## total change in soil C
carbon_seq_soilC<-carbon_seq[which(!is.na(carbon_seq$soilC)),]
mst<-lm(soilC~species_richness+FDis+PSV,data=carbon_seq_soilC)
mst_block<-lmer(soilC~species_richness+FDis+PSV+(1|block),
                data=carbon_seq_soilC,REML=F)

options(na.action = "na.fail")
mst_dredge <- dredge(mst, beta = "none", evaluate = T, rank = AICc)
mst_block_dredge <- dredge(mst_block, beta = "none", evaluate = T, rank = AICc)
options(na.action = "na.omit")

mst_model<-get.models(mst_dredge, subset = 1)[[1]]
check_model(mst_model)

####
## overyielding in (change in) soil C
carbon_seq_soilOY<-carbon_seq[which(!is.na(carbon_seq$soilOY)),]
ms<-lm(soilOY~species_richness+FDis+PSV,data=carbon_seq_soilOY)
ms_block<-lmer(soilOY~species_richness+FDis+PSV+(1|block),
               data=carbon_seq_soilOY,REML=F)

options(na.action = "na.fail")
ms_dredge <- dredge(ms, beta = "none", evaluate = T, rank = AICc)
ms_block_dredge <- dredge(ms_block, beta = "none", evaluate = T, rank = AICc)
options(na.action = "na.omit")

## the 'selected' model is actually the second
ms_model<-get.models(ms_dredge, subset = 1)[[1]]
check_model(ms_model)

####
## total fine root C
carbon_seq_rootC<-carbon_seq[which(!is.na(carbon_seq$rootC)),]
mrt<-lm(rootC~species_richness+FDis+PSV,data=carbon_seq_rootC)
mrt_block<-lmer(rootC~species_richness+FDis+PSV+(1|block),
                data=carbon_seq_rootC,REML=F)

options(na.action = "na.fail")
mrt_dredge <- dredge(mrt, beta = "none", evaluate = T, rank = AICc)
mrt_block_dredge <- dredge(mrt_block, beta = "none", evaluate = T, rank = AICc)
options(na.action = "na.omit")

## the 'selected' model is actually the second
mrt_model<-get.models(mrt_dredge, subset = 1)[[1]]
check_model(mrt_model)

####
## overyielding in fine root C
carbon_seq_rootOY<-carbon_seq[which(!is.na(carbon_seq$rootOY)),]
mr<-lm(rootOY~species_richness+FDis+PSV,data=carbon_seq_rootOY)
mr_block<-lmer(rootOY~species_richness+FDis+PSV+(1|block),
               data=carbon_seq_rootOY,REML=F)

options(na.action = "na.fail")
mr_dredge <- dredge(mr, beta = "none", evaluate = T, rank = AICc)
mr_block_dredge <- dredge(mr_block, beta = "none", evaluate = T, rank = AICc)
options(na.action = "na.omit")

mr_model<-get.models(mr_dredge, subset = 1)[[1]]
check_model(mr_model)

####
## macroaggregates
## we add block here because it seems to consistently
## improve models + change parameter estimates
carbon_seq_macro<-carbon_seq[which(!is.na(carbon_seq$macro250)),]
mma<-lm(macro250~species_richness+FDis+PSV,data=carbon_seq_macro)
mma_block<-lmer(macro250~species_richness+FDis+PSV+(1|block),
                data=carbon_seq_macro,REML=F)

options(na.action = "na.fail")
mma_dredge <- dredge(mma, beta = "none", evaluate = T, rank = AICc)
mma_block_dredge <- dredge(mma_block, beta = "none", evaluate = T, rank = AICc)
options(na.action = "na.omit")

mma_model<-get.models(mma_block_dredge, subset = 1)[[1]]
check_model(mma_model)

#####################################
## structural equation models
## here we use local estimation in piecewiseSEM

## drop rows with missing data
carbon_seq_sub<-carbon_seq[-which(is.na(carbon_seq$soilC) | is.na(carbon_seq$macro250)),]

## z-standardize important variables
carbon_seq_standard<-carbon_seq_sub
standard_cols<-c("species_richness","woodyC","soilC",
                 "macro250","percentAM","percentCon")
carbon_seq_standard[,standard_cols]<-scale(carbon_seq_standard[,standard_cols])

## this was the original proposed model
## before examining any SEM output
localfit_orig_model<-psem(
  lmer(woodyC~species_richness+percentCon+(1|block),data=carbon_seq_standard),
  lmer(macro250~percentAM+percentCon+woodyC+(1|block),data=carbon_seq_standard),
  lmer(soilC~species_richness+percentAM+percentCon+woodyC+macro250+(1|block),
     data=carbon_seq_standard)
)

## a full model
localfit_full_model<-psem(
  lmer(woodyC~species_richness+percentCon+percentAM+(1|block),data=carbon_seq_standard),
  lmer(macro250~percentAM+percentCon+(1|block),data=carbon_seq_standard),
  lm(soilC~species_richness+percentAM+woodyC+macro250,
     data=carbon_seq_standard)
)

## checking individual components
# check_model(lmer(woodyC~species_richness+percentCon+percentAM+(1|block),data=carbon_seq_standard))
# check_model(lmer(macro250~percentAM+percentCon+(1|block),data=carbon_seq_standard))
# check_model(lm(soilC~species_richness+percentAM+woodyC+macro250,
#                data=carbon_seq_standard))

## factors that influence leaf litter quantity (woody C)
## or quality (% needleleaf) don't seem to influence soil C
localfit_sub_model<-psem(
  lmer(woodyC~species_richness+percentCon+percentAM+(1|block),data=carbon_seq_standard),
  lmer(macro250~percentAM+percentCon+(1|block),data=carbon_seq_standard),
  lm(soilC~species_richness+percentAM+macro250,
     data=carbon_seq_standard)
)

localfit_min_model<-psem(
  lmer(woodyC~species_richness+percentCon+percentAM+(1|block),data=carbon_seq_standard),
  lmer(macro250~percentCon+(1|block),data=carbon_seq_standard),
  lm(soilC~species_richness+percentAM,
     data=carbon_seq_standard)
)

