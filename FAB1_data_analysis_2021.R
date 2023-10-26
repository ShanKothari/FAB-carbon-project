setwd("C:/Users/querc/Dropbox/FABCarbonProject/")

library(reshape2)
library(ggplot2)
library(ggpubr)
library(pdiv) ## Pascal Niklaus's pdiv package
library(lme4)
library(lmerTest)
library(MuMIn)
library(car)
library(performance)

FABdata<-read.csv("ProcessedData/FAB_Cestimate.csv")

## aggregate aboveground woody carbon by plot
## *10000/(16*1000) changes units from kg per plot to Mg per hectare
C_agg<-aggregate(C_estimate~plot,data=FABdata,
                 FUN=function(x) sum(x,na.rm=T)*10/16)
colnames(C_agg)<-c("plot","woodyC")

plot_guide<-read.csv("OriginalData/plotkey_biomass.csv")
C_agg$species_richness<-plot_guide$Treatment[match(C_agg$plot,plot_guide$Plot)]
C_agg$block<-plot_guide$Block[match(C_agg$plot,plot_guide$Plot)]

###################################################
## overyielding calculations

## this is to make sure that dead trees are included as 0s for means,
## but only up to one per position
FABdata<-FABdata[FABdata$surveyed_2019=="Yes",]
FABdata$C_estimate[which(is.na(FABdata$C_estimate))]<-0

## calculate overyielding matching to the monoculture of the same block
mono.means<-aggregate(C_estimate~species_code+block,
                      data=FABdata[FABdata$species_richness==1,],
                      FUN=mean)
FABdata$mono.means<-apply(FABdata,1,
                          function(x) {
                            mono.means$C_estimate[mono.means$species_code==x["species_code"] & mono.means$block==x["block"]]
                          })
FABdata$ind.OY<-FABdata$C_estimate-FABdata$mono.means

OY.agg<-aggregate(ind.OY~plot,data=FABdata,
                  FUN=function(x) sum(x,na.rm=T)*10/16)
OY.agg$species_richness<-plot_guide$Treatment[match(OY.agg$plot,plot_guide$Plot)]
OY.agg$block<-plot_guide$Block[match(OY.agg$plot,plot_guide$Plot)]
OY.agg<-OY.agg[-which(OY.agg$species_richness==1),]

############################
## aboveground woody CE/SE calculations

## for calculations to work, need to fix misplant in ACRU monoculture 147
FABdata_mod<-FABdata
## these two lines don't seem to be doing anything?
## due to misplants CE + SE don't quite add up to OY as calculated above
replacement<-mean(FABdata_mod$C_estimate[FABdata_mod$species_code=="ACRU" & FABdata_mod$plot==147],na.rm=T)
FABdata_mod$C_estimate[FABdata_mod$species_code=="TIAM" & FABdata_mod$plot==147]<-replacement
FABdata_mod$species_code[FABdata_mod$species_code=="TIAM" & FABdata_mod$plot==147]<-"ACRU"

## generate indicators of species composition
FABplot_list<-split(FABdata_mod,f = FABdata_mod$plot)
FABplot_comp<-unlist(lapply(FABplot_list,function(plot) paste(unique(plot$species_code),collapse="|")))

## get estimates of total woody carbon per species per plot
C_sp_plot<-aggregate(C_estimate~species_code+plot+block,
                     data=FABdata_mod, FUN=function(x) sum(x)*10/16)
## add composition indicators
C_sp_plot$sp_comp<-FABplot_comp[match(C_sp_plot$plot,names(FABplot_comp))]

## estimate fractions of each species in each plot
## fulcrum_id is just a haphazardly chosen (and irrelevant) variable
counts<-aggregate(fulcrum_id~species_code+plot,
                  data=FABdata_mod,
                  FUN=length)
counts$fractions<-counts$fulcrum_id/64

## check that they have the same order of species and plots
## this should return the number of rows in both counts and C_sp_plot
sum(counts$species_code==C_sp_plot$species_code & counts$plot==C_sp_plot$plot)
C_sp_plot$fractions<-counts$fractions

C_partition<-addpart(C_estimate~sp_comp/species_code+plot,
                     fractions= ~fractions,
                     groups= ~block,
                     data=C_sp_plot)

# ## trying calculations manually
# C_sp_plot$mono.means<-apply(C_sp_plot,1,
#                             function(x) {
#                               C_sp_plot$C_estimate[C_sp_plot$sp_comp==x["species_code"] & C_sp_plot$block==x["block"]]
#                             })
# C_sp_plot$mono.exp<-C_sp_plot$mono.means*C_sp_plot$fractions
# C_sp_plot$RY_Oi<-C_sp_plot$C_estimate/C_sp_plot$mono.means
# C_sp_plot$deltaRY<-C_sp_plot$RY_Oi-C_sp_plot$fractions
# 
# ## population covariance, borrowed from Pascal Niklaus
# covp <- function(x, y) {
#   sum((x - mean(x)) * (y - mean(y))) / length(x);
# }
# 
# ## for a plot, sp richness*covp(mono.means,deltaRY) is selection
# ## and sp richness*mean(mono.means)*mean(deltaRY) is complementarity

#############################
## proportions of each species planted

counts$fulcrum_id<-NULL
FAB_planted<-reshape(counts,
                         idvar = "plot",
                         timevar = "species_code",
                         direction = "wide")

colnames(FAB_planted)<-gsub(pattern="fractions.",
                                replacement="",
                                x=colnames(FAB_planted))

FAB_planted[which(is.na(FAB_planted),arr.ind=T)]<-0
FAB_planted$percentAM<-rowSums(FAB_planted[,c("ACNE","ACRU","JUVI")],na.rm=T)
FAB_planted$percentCon<-rowSums(FAB_planted[,c("JUVI","PIBA","PIRE","PIST")],na.rm=T)

#############################
## light data from 2018 (just for context)

ground_light<-read.csv("OriginalData/ground_light_processed_2018.csv",
                       fileEncoding="latin1")

ground_light$Annotation<-gsub(pattern="PLOT",
                              replacement="",
                              x=ground_light$Annotation)
ground_light$Annotation<-gsub(pattern="SK",
                              replacement="",
                              x=ground_light$Annotation)

C_agg$ground_light_2018<-ground_light$Tau....[match(C_agg$plot,
                                               ground_light$Annotation)]

############################
## belowground C

belowground<-read.csv("OriginalData/FAB5 'FINAL' Data - Data.csv")
belowground$above_C<-C_agg$woodyC[match(belowground$Plot,C_agg$plot)]
belowground$aboveOY<-OY.agg$ind.OY[match(belowground$Plot,OY.agg$plot)]

belowground$block<-plot_guide$Block[match(belowground$Plot,plot_guide$Plot)]

## estimate bulk density per block
BD<-read.csv("OriginalData/FAB_BD.csv")
BD_wide<-dcast(Plot.ID~Depth,
               data=BD[,c("Plot.ID","Depth","Bulk.density.g.cm3")])
BD_wide$BD<-3/4*BD_wide$`0-15`+1/4*BD_wide$`15-30`
BD_wide$block<-plot_guide$Block[match(BD_wide$Plot.ID,plot_guide$Plot)]

belowground$BD<-NA
belowground$BD[belowground$block==1]<-mean(BD_wide$BD[BD_wide$block==1])
belowground$BD[belowground$block==2]<-mean(BD_wide$BD[BD_wide$block==2])
belowground$BD[belowground$block==3]<-mean(BD_wide$BD[BD_wide$block==3])

## estimate accumulation of soil carbon
## first take the difference in percent
## (which assumes bulk density remains unchanged)
belowground$soilC_diff<-belowground$X..C_2019-belowground$X.C_2013
belowground$soilN_diff<-belowground$X..N_2019-belowground$X.N_2013
## to g C accumulated / cm^3 soil
belowground$soilC_diff_vol<-belowground$soilC_diff/100*belowground$BD
belowground$soilN_diff_vol<-belowground$soilN_diff/100*belowground$BD
## 200000 cm^3 per m^2 sampled, 10000 m^2 per ha, 1000000 g per Mg
belowground$soilC_diff_plot<-belowground$soilC_diff_vol*200000*10000/1000000
belowground$soilN_diff_plot<-belowground$soilN_diff_vol*200000*10000/1000000

## soil C pool at the end (not change)
belowground$soilC_pool_vol<-belowground$X..C_2019/100*belowground$BD
belowground$soilC_pool_plot<-belowground$soilC_pool_vol*200000*10000/1000000

## root carbon, assuming roots are 50% carbon
## sampled from 2 in diameter root cores, 5 cores per plot
## 10000 cm2 to m2, 10000 m2 per ha, 1/1000000 Mg per g
belowground$rootC<-belowground$Roots.Dry.Mass/(5*2.54^2*pi)*10000*10000*0.5/1000000
## these plots weren't actually sampled
belowground$rootC[belowground$Plot %in% c(60,76,96,114)]<-NA

sp_comp<-belowground[,c("Plot","Sp.Richness","block","soilC_diff_plot","rootC",
                        "ACNE","ACRU","BEPA","JUVI",
                        "PIBA","PIRE","PIST","QUAL",
                        "QUEL","QUMA","QURU","TIAM")]
sp_comp_long<-melt(sp_comp,id.vars = c("Plot","Sp.Richness","block","soilC_diff_plot","rootC"))

## attach monoculture
sp_comp_long$mono_soil<-unlist(apply(sp_comp_long,1,function(x) {
  mono_block<-which(sp_comp_long$value>0.99 & sp_comp_long$variable==x["variable"] & sp_comp_long$block==x["block"])
  return(sp_comp_long$soilC_diff_plot[mono_block])
}))

sp_comp_long$mono_root<-unlist(apply(sp_comp_long,1,function(x) {
  mono_block<-which(sp_comp_long$value>0.99 & sp_comp_long$variable==x["variable"] & sp_comp_long$block==x["block"])
  return(sp_comp_long$rootC[mono_block])
}))

sp_comp_long$mono_soil_exp<-sp_comp_long$mono_soil*sp_comp_long$value
sp_comp_long$mono_root_exp<-sp_comp_long$mono_root*sp_comp_long$value

tot_mono_soil_exp<-aggregate(mono_soil_exp~Plot,data=sp_comp_long,FUN=sum)
tot_mono_root_exp<-aggregate(mono_root_exp~Plot,data=sp_comp_long,FUN=sum)

belowground$mono_soil_exp<-tot_mono_soil_exp$mono_soil_exp[match(belowground$Plot,tot_mono_soil_exp$Plot)]
belowground$soil_OY<-belowground$soilC_diff_plot-belowground$mono_soil_exp
belowground$mono_root_exp<-tot_mono_root_exp$mono_root_exp[match(belowground$Plot,tot_mono_root_exp$Plot)]
belowground$root_OY<-belowground$rootC-belowground$mono_root_exp

belowground_sub<-belowground[-which(belowground$Sp.Richness==1),]
# plot(belowground_sub$soil_OY~belowground_sub$Sp.Richness)
# plot(belowground_sub$root_OY~belowground_sub$Sp.Richness)

###############################
## combine

C_agg$block<-plot_guide$Block[match(C_agg$plot,plot_guide$Plot)]

## AG woody OY, CE, and SE
C_agg$woodyOY<-OY.agg$ind.OY[match(C_agg$plot,OY.agg$plot)]
C_agg$woodyCE<-C_partition$CE.C_estimate[match(C_agg$plot,C_partition$plot)]
C_agg$woodySE<-C_partition$SE.C_estimate[match(C_agg$plot,C_partition$plot)]

## percentages of soil C and N
belowground_match<-match(C_agg$plot,belowground$Plot)
C_agg$perC_init<-belowground$X.C_2013[belowground_match]
C_agg$perC<-belowground$X..C_2019[belowground_match]
C_agg$perN<-belowground$X..N_2019[belowground_match]

## change in soil C at plot scale
C_agg$soilC<-belowground$soilC_diff_plot[belowground_match]

## total soil C at scale
C_agg$soilC_pool<-belowground$soilC_pool_plot[belowground_match]

## soil C overyielding
C_agg$soilOY<-belowground_sub$soil_OY[match(C_agg$plot,belowground_sub$Plot)]

## root C and overyielding and plot scale
C_agg$rootC<-belowground$rootC[belowground_match]
C_agg$rootOY<-belowground_sub$root_OY[match(C_agg$plot,belowground_sub$Plot)]

## adding together AG wood, soil, and roots
C_agg$totalC<-C_agg$woodyC+C_agg$soilC+C_agg$rootC
C_agg$totalOY<-C_agg$woodyOY+C_agg$soilOY+C_agg$rootOY

## macroaggregates
C_agg$macro250<-belowground$X..Mass.of.250[belowground_match]

## soil moisture
C_agg$soil_moisture<-belowground$Soil.moisture[belowground_match]

## planted proportions of AM and coniferous trees
C_agg$percentAM<-FAB_planted$percentAM[match(C_agg$plot,FAB_planted$plot)]
C_agg$percentCon<-FAB_planted$percentCon[match(C_agg$plot,FAB_planted$plot)]

## change in soil N
C_agg$soilN<-belowground$soilN_diff_plot[belowground_match]

## log (root C/AG woody C)
C_agg$logBA<-log(C_agg$rootC/C_agg$woodyC)

## functional and phylogenetic diversity
FD<-read.csv("OriginalData/ecy1958-sup-0003-tables1.csv")
C_agg$PSV<-FD$PSV[match(C_agg$plot,FD$Plot)]
C_agg$FDis<-FD$FDis[match(C_agg$plot,FD$Plot)]

## leaf type and mycotype
## the 0.97 here is for plot 147, which has one Tilia in it
## (by mistake)
C_agg$mycotype<-ifelse(C_agg$percentAM==0,"E",
                       ifelse(C_agg$percentAM>0.97,"A","B"))
C_agg$leaf_type<-ifelse(C_agg$percentCon==0,"D",
                        ifelse(C_agg$percentCon==1,"C","B"))

## optional step to set PSV to 0 in monocultures
C_agg$PSV[C_agg$species_richness==1]<-0

write.csv(C_agg,"ProcessedData/Cseq.csv",row.names = F)

##############################
## simple light analyses

C_agg$BEPA<-FAB_planted$BEPA[match(C_agg$plot,FAB_planted$plot)]

quantile(C_agg$ground_light_2018[C_agg$leaf_type=="C"],
         probs=c(0.1,0.5,0.9),na.rm=T)
quantile(C_agg$ground_light_2018[C_agg$leaf_type=="B"],
         probs=c(0.1,0.5,0.9),na.rm=T)
quantile(C_agg$ground_light_2018[C_agg$leaf_type=="D" & C_agg$BEPA>0],
         probs=c(0.1,0.5,0.9),na.rm=T)
quantile(C_agg$ground_light_2018[C_agg$leaf_type=="D" & is.na(C_agg$BEPA)],
         probs=c(0.1,0.5,0.9),na.rm=T)

#############################
## plot monoculture biomass

C_agg_mono<-C_agg[which(C_agg$species_richness==1),]
C_agg_mono$species<-C_sp_plot$sp_comp[match(C_agg_mono$plot,C_sp_plot$plot)]

wood_mono<-ggplot(C_agg_mono,aes(x=species,y=woodyC))+
  geom_boxplot()+theme_bw()+
  theme(text=element_text(size=20),
        axis.text.x = element_blank(),
        axis.title.x = element_blank())+
  labs(y="Woody C (kg per plot)",x="Species")

soil_mono<-ggplot(C_agg_mono,aes(x=species,y=soilC))+
  geom_boxplot()+theme_bw()+
  theme(text=element_text(size=20),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  labs(y="Change in soil C (kg per plot)",x="Species")

png("Images/wood_soil_mono.png",width=8,height=8,units = "in",res=150)
ggarrange(plotlist = list(wood_mono,soil_mono),
          ncol=1,heights=c(1,1.3))
dev.off()

#############################
## various basic plots

png("Images/soilC_div.png",width=5,height=5,units = "in",res=150)
ggplot(C_agg,aes(x=species_richness,y=soilC))+
  geom_point(size=2)+theme_bw()+
  geom_smooth(method="lm",se=F)+
  theme(text = element_text(size=20),
        plot.margin = margin(0.1,0.2,0,0,"in"))+
  labs(x="Species richness",y=expression("Soil C accrual (kg plot"^-1*")"))
dev.off()

png("Images/soilOY_div.png",width=5,height=5,units = "in",res=150)
ggplot(C_agg,aes(x=species_richness,y=soilOY))+
  geom_point(size=2)+theme_bw()+
  geom_smooth(method="lm",se=F)+
  theme(text = element_text(size=20),
        plot.margin = margin(0.1,0.2,0,0,"in"))+
  labs(x="Species richness",y=expression("Soil C overyielding (kg plot"^-1*")"))
dev.off()

png("Images/woodyC_div.png",width=5,height=5,units = "in",res=150)
ggplot(C_agg,aes(x=species_richness,y=woodyC))+
  geom_point(size=2)+theme_bw()+
  theme(text = element_text(size=20),
        plot.margin = margin(0.1,0.2,0,0,"in"))+
  labs(x="Species richness",y=expression("Wood C (kg plot"^-1*")"))
dev.off()

png("Images/woodyOY_div.png",width=5,height=5,units = "in",res=150)
ggplot(C_agg,aes(x=species_richness,y=woodyOY))+
  geom_point(size=2)+
  geom_smooth(method="lm",se=F)+
  theme_bw()+
  theme(text = element_text(size=20),
        plot.margin = margin(0.1,0.2,0,0,"in"))+
  labs(x="Species richness",y=expression("Wood C overyielding (kg plot"^-1*")"))
dev.off()

png("Images/soil_wood_C.png",width=5,height=5,units = "in",res=150)
ggplot(C_agg,aes(x=soilC,y=woodyC))+
  geom_point(size=2)+
  theme_bw()+
  theme(text = element_text(size=20),
        plot.margin = margin(0.1,0.2,0,0,"in"))+
  labs(x=expression("Soil C accrual (kg plot"^-1*")"),
       y=expression("Wood C (kg plot"^-1*")"))
dev.off()

#############################
## AIC-based model selection for predictors of C stocks and OY

library(lme4)
library(lmerTest)

####
## total woody C
mwt<-lm(woodyC~species_richness+FDis+PSV,data=C_agg)
mwt_block<-lmer(woodyC~species_richness+FDis+PSV+(1|block),
                data=C_agg, REML=F)
mwt_robust<-rlm(woodyC~species_richness+FDis+PSV,data=C_agg)

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
C_agg_woodyOY<-C_agg[which(!is.na(C_agg$woodyOY)),]
mw<-lm(woodyOY~species_richness+FDis+PSV,
       data=C_agg_woodyOY)
mw_block<-lmer(woodyOY~species_richness+FDis+PSV+(1|block),
               data=C_agg_woodyOY,REML=F)
mw_robust<-rlm(woodyOY~species_richness+FDis+PSV,data=C_agg_woodyOY)

options(na.action = "na.fail")
mw_dredge <- dredge(mw, beta = "none", evaluate = T, rank = AICc)
mw_block_dredge <- dredge(mw_block, beta = "none", evaluate = T, rank = AICc)
mw_robust_dredge <- dredge(mw_robust, beta = "none", evaluate = T, rank = AICc)
options(na.action = "na.omit")

mw_model<-get.models(mw_dredge, subset = 1)[[1]]
check_model(mw_model)

####
## complementarity effects
mce<-lm(woodyCE~species_richness+FDis+PSV,data=C_agg_woodyOY)
mce_block<-lmer(woodyCE~species_richness+FDis+PSV+(1|block),
                data=C_agg_woodyOY, REML=F)
mce_robust<-rlm(woodyCE~species_richness+FDis+PSV,data=C_agg_woodyOY)

options(na.action = "na.fail")
mce_dredge <- dredge(mce, beta = "none", evaluate = T, rank = AICc)
mce_block_dredge <- dredge(mce_block, beta = "none", evaluate = T, rank = AICc)
mce_robust_dredge <- dredge(mce_robust, beta = "none", evaluate = T, rank = AICc)
options(na.action = "na.omit")

mce_model<-get.models(mce_dredge, subset = 2)[[1]]
check_model(mce_model)

####
## selection effects
mse<-lm(woodySE~species_richness+FDis+PSV,data=C_agg_woodyOY)
mse_block<-lmer(woodySE~species_richness+FDis+PSV+(1|block),
                data=C_agg_woodyOY,REML=F)
mse_robust<-rlm(woodySE~species_richness+FDis+PSV,data=C_agg_woodyOY)

options(na.action = "na.fail")
mse_dredge <- dredge(mse, beta = "none", evaluate = T, rank = AICc)
mse_block_dredge <- dredge(mse_block, beta = "none", evaluate = T, rank = AICc)
mse_robust_dredge <- dredge(mse_robust, beta = "none", evaluate = T, rank = AICc)
options(na.action = "na.omit")

mse_model<-get.models(mse_dredge, subset = 1)[[1]]
check_model(mse_model)

####
## total change in soil C
C_agg_soilC<-C_agg[which(!is.na(C_agg$soilC)),]
mst<-lm(soilC~species_richness+FDis+PSV,data=C_agg_soilC)
mst_block<-lmer(soilC~species_richness+FDis+PSV+(1|block),
                data=C_agg_soilC,REML=F)

options(na.action = "na.fail")
mst_dredge <- dredge(mst, beta = "none", evaluate = T, rank = AICc)
mst_block_dredge <- dredge(mst_block, beta = "none", evaluate = T, rank = AICc)
options(na.action = "na.omit")

mst_model<-get.models(mst_dredge, subset = 1)[[1]]
check_model(mst_model)

####
## overyielding in (change in) soil C
C_agg_soilOY<-C_agg[which(!is.na(C_agg$soilOY)),]
ms<-lm(soilOY~species_richness+FDis+PSV,data=C_agg_soilOY)
ms_block<-lmer(soilOY~species_richness+FDis+PSV+(1|block),
               data=C_agg_soilOY,REML=F)

options(na.action = "na.fail")
ms_dredge <- dredge(ms, beta = "none", evaluate = T, rank = AICc)
ms_block_dredge <- dredge(ms_block, beta = "none", evaluate = T, rank = AICc)
options(na.action = "na.omit")

## the 'selected' model is actually the second
ms_model<-get.models(ms_dredge, subset = 1)[[1]]
check_model(ms_model)

####
## total fine root C
C_agg_rootC<-C_agg[which(!is.na(C_agg$rootC)),]
mrt<-lm(rootC~species_richness+FDis+PSV,data=C_agg_rootC)
mrt_block<-lmer(rootC~species_richness+FDis+PSV+(1|block),
                data=C_agg_rootC,REML=F)

options(na.action = "na.fail")
mrt_dredge <- dredge(mrt, beta = "none", evaluate = T, rank = AICc)
mrt_block_dredge <- dredge(mrt_block, beta = "none", evaluate = T, rank = AICc)
options(na.action = "na.omit")

## the 'selected' model is actually the second
mrt_model<-get.models(mrt_dredge, subset = 1)[[1]]
check_model(mrt_model)

####
## overyielding in fine root C
C_agg_rootOY<-C_agg[which(!is.na(C_agg$rootOY)),]
mr<-lm(rootOY~species_richness+FDis+PSV,data=C_agg_rootOY)
mr_block<-lmer(rootC~species_richness+FDis+PSV+(1|block),
                data=C_agg_rootOY,REML=F)

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
C_agg_macro<-C_agg[which(!is.na(C_agg$macro250)),]
mma<-lm(macro250~species_richness+FDis+PSV,data=C_agg_macro)
mma_block<-lmer(macro250~species_richness+FDis+PSV+(1|block),
                data=C_agg_macro,REML=F)

options(na.action = "na.fail")
mma_dredge <- dredge(mma, beta = "none", evaluate = T, rank = AICc)
mma_block_dredge <- dredge(mma_block, beta = "none", evaluate = T, rank = AICc)
options(na.action = "na.omit")

mma_model<-get.models(mma_block_dredge, subset = 1)[[1]]
check_model(mma_model)
