setwd("C:/Users/querc/Dropbox/FABCarbonProject/")

library(reshape2)
library(ggplot2)
library(pdiv) ## Pascal Niklaus's pdiv package

FABdata<-read.csv("ProcessedData/FAB_Cestimate.csv")

## aggregate carbon by plot
C_agg<-aggregate(FABdata$C_estimate,by=list(FABdata$plot),FUN=sum,na.rm=T)
colnames(C_agg)<-c("plot","woodyC")
## these plots aren't actually surveyed
C_agg<-C_agg[-which(C_agg$plot %in% c(112,119)),]

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

OY.agg<-aggregate(ind.OY~plot,data=FABdata,FUN=sum)
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
                     data=FABdata_mod,
                     FUN=sum)
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

## estimate accumulation of soil carbon -- first take the difference in percent
## assumes bulk density remains unchanged
belowground$soilC_diff<-belowground$X..C_2019-belowground$X.C_2013
belowground$soilN_diff<-belowground$X..N_2019-belowground$X.N_2013
## to g C accumulated / cm^3 soil
belowground$soilC_diff_vol<-belowground$soilC_diff/100*belowground$BD
belowground$soilN_diff_vol<-belowground$soilN_diff/100*belowground$BD
## 200000 cm^3 per m^2 sampled, 16 m^2 per plot, 1000 g per kg
belowground$soilC_diff_plot<-belowground$soilC_diff_vol*200000*16/1000
belowground$soilN_diff_plot<-belowground$soilN_diff_vol*200000*16/1000

## soil C pool (not change)
belowground$soilC_pool_vol<-belowground$X..C_2019/100*belowground$BD
belowground$soilC_pool_plot<-belowground$soilC_pool_vol*200000*16/1000

## root carbon, assuming roots are 50% carbon
## sampled from 2 in diameter root cores, 5 cores per plot
belowground$rootC<-belowground$Roots.Dry.Mass/(5*2.54^2*pi)*10000*16*0.5/1000
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
plot(belowground_sub$soil_OY~belowground_sub$Sp.Richness)
plot(belowground_sub$root_OY~belowground_sub$Sp.Richness)

###############################
## combine

C_agg$block<-plot_guide$Block[match(C_agg$plot,plot_guide$Plot)]
C_agg$woodyOY<-OY.agg$ind.OY[match(C_agg$plot,OY.agg$plot)]
C_agg$woodyCE<-C_partition$CE.C_estimate[match(C_agg$plot,C_partition$plot)]
C_agg$woodySE<-C_partition$SE.C_estimate[match(C_agg$plot,C_partition$plot)]
C_agg$soilC<-belowground$soilC_diff_plot[match(C_agg$plot,belowground$Plot)]
C_agg$soilC_pool<-belowground$soilC_pool_plot[match(C_agg$plot,belowground$Plot)]
C_agg$soilOY<-belowground_sub$soil_OY[match(C_agg$plot,belowground_sub$Plot)]
C_agg$rootC<-belowground$rootC[match(C_agg$plot,belowground$Plot)]
C_agg$rootOY<-belowground_sub$root_OY[match(C_agg$plot,belowground_sub$Plot)]
C_agg$totalC<-C_agg$woodyC+C_agg$soilC+C_agg$rootC
C_agg$totalOY<-C_agg$woodyOY+C_agg$soilOY+C_agg$rootOY

C_agg$soilN<-belowground$soilN_diff_plot[match(C_agg$plot,belowground$Plot)]
C_agg$logBA<-log(C_agg$rootC/C_agg$woodyC)

#############################
## plot monoculture biomass

C_agg_mono<-C_agg[which(C_agg$species_richness==1),]
C_agg_mono$species<-C_sp_plot$sp_comp[match(C_agg_mono$plot,C_sp_plot$plot)]

ggplot(C_agg_mono,aes(x=species,y=log10(woodyC)))+
  geom_boxplot()+theme_bw()+
  labs(y="Woody C (kg per plot)")

ggplot(C_agg_mono,aes(x=species,y=soilC))+
  geom_boxplot()+theme_bw()+
  labs(y="Change in soil C (kg per plot)")

#############################
## various basic plots

png("Images/soilC_div.png",width=5,height=5,units = "in",res=150)
ggplot(C_agg,aes(x=species_richness,y=soilC))+
  geom_point(size=2)+theme_bw()+
  theme(text = element_text(size=20),
        plot.margin = margin(0.1,0.2,0,0,"in"))+
  labs(x="Species richness",y=expression("Soil C accrual (kg plot"^-1*")"))
dev.off()

ggplot(C_agg,aes(x=species_richness,y=soilOY))+
  geom_point(size=2)+theme_bw()+
  theme(text = element_text(size=20),
        plot.margin = margin(0.1,0.2,0,0,"in"))+
  labs(x="Species richness",y=expression("Soil C overyielding (kg plot"^-1*")"))

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
## add functional and phylogenetic diversity

FD<-read.csv("OriginalData/ecy1958-sup-0003-tables1.csv")
C_agg$PSV<-FD$PSV[match(C_agg$plot,FD$Plot)]
C_agg$FDis<-FD$FDis[match(C_agg$plot,FD$Plot)]

write.csv(C_agg,"ProcessedData/Cseq.csv",row.names = F)

#############################
## AIC-based model selection for predictors of OY

library(lme4)
library(lmerTest)

mw_null<-lm(woodyOY~1,data=C_agg)
mw1<-lm(woodyOY~species_richness,data=C_agg)
mw2<-lm(woodyOY~FDis,data=C_agg)
mw3<-lm(woodyOY~PSV,data=C_agg)
mw4<-lm(woodyOY~species_richness+FDis,data=C_agg)
mw5<-lm(woodyOY~species_richness+PSV,data=C_agg)
mw6<-lm(woodyOY~species_richness+FDis+PSV,data=C_agg)
AIC(mw_null)
AIC(mw1)
AIC(mw2)
AIC(mw3)
AIC(mw4)
AIC(mw5)
AIC(mw6)

mce_null<-lm(woodyCE~1,data=C_agg)
mce1<-lm(woodyCE~species_richness,data=C_agg)
mce2<-lm(woodyCE~FDis,data=C_agg)
mce3<-lm(woodyCE~PSV,data=C_agg)
mce4<-lm(woodyCE~species_richness+FDis,data=C_agg)
mce5<-lm(woodyCE~species_richness+PSV,data=C_agg)
mce6<-lm(woodyCE~species_richness+FDis+PSV,data=C_agg)
AIC(mce_null)
AIC(mce1)
AIC(mce2)
AIC(mce3)
AIC(mce4)
AIC(mce5)
AIC(mce6)

mse_null<-lm(woodySE~1,data=C_agg)
mse1<-lm(woodySE~species_richness,data=C_agg)
mse2<-lm(woodySE~FDis,data=C_agg)
mse3<-lm(woodySE~PSV,data=C_agg)
mse4<-lm(woodySE~species_richness+FDis,data=C_agg)
mse5<-lm(woodySE~species_richness+PSV,data=C_agg)
mse6<-lm(woodySE~species_richness+FDis+PSV,data=C_agg)
AIC(mse_null)
AIC(mse1)
AIC(mse2)
AIC(mse3)
AIC(mse4)
AIC(mse5)
AIC(mse6)

ms_null<-lm(soilOY~1,data=C_agg)
ms1<-lm(soilOY~species_richness,data=C_agg)
ms2<-lm(soilOY~FDis,data=C_agg)
ms3<-lm(soilOY~PSV,data=C_agg)
ms4<-lm(soilOY~species_richness+FDis,data=C_agg)
ms5<-lm(soilOY~species_richness+PSV,data=C_agg)
ms6<-lm(soilOY~species_richness+FDis+PSV,data=C_agg)
AIC(ms_null)
AIC(ms1)
AIC(ms2)
AIC(ms3)
AIC(ms4)
AIC(ms5)
AIC(ms6)

mr_null<-lm(rootOY~1,data=C_agg)
mr1<-lm(rootOY~species_richness,data=C_agg)
mr2<-lm(rootOY~FDis,data=C_agg)
mr3<-lm(rootOY~PSV,data=C_agg)
mr4<-lm(rootOY~species_richness+FDis,data=C_agg)
mr5<-lm(rootOY~species_richness+PSV,data=C_agg)
mr6<-lm(rootOY~species_richness+FDis+PSV,data=C_agg)
AIC(mr_null)
AIC(mr1)
AIC(mr2)
AIC(mr3)
AIC(mr4)
AIC(mr5)
AIC(mr6)