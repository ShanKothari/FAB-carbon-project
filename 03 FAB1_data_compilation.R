setwd("C:/Users/querc/Dropbox/FABCarbonProject/")

library(reshape2)
library(pdiv) ## Pascal Niklaus's pdiv package
## installed via (with devtools):
# install_github("pascal-niklaus/pdiv/pdiv")

## as a reminder, this data set has the FAB inventory data
## where each row is one individual. in the second step (02)
## we estimated the carbon content of each tree in 2019
## based on measured height and measured/estimated basal diameter
FABdata<-read.csv("ProcessedData/FAB_Cestimate.csv")

## examine mortality rates across plots
# FABdata_alive<-FABdata[which(FABdata$deadmissing_2019=="No"),]
# max((64-table(FABdata_alive$plot))/64)
# min((64-table(FABdata_alive$plot))/64)
# median((64-table(FABdata_alive$plot))/64)

####################################
## aggregate aboveground woody carbon by plot

## *10000/(16*1000) changes units from kg per plot to Mg per hectare
C_agg<-aggregate(C_estimate~plot,data=FABdata,
                 FUN=function(x) sum(x,na.rm=T)*10/16)
colnames(C_agg)<-c("plot","woodyC")

## get the species richness and block for each plot
plot_guide<-read.csv("OriginalData/plotkey_biomass.csv")
C_agg$species_richness<-plot_guide$Treatment[match(C_agg$plot,plot_guide$Plot)]
C_agg$block<-plot_guide$Block[match(C_agg$plot,plot_guide$Plot)]

###################################################
## overyielding calculations

## this is to make sure that dead trees are included as 0s for means;
## first we have to get rid of trees that were replaced during replanting
FABdata<-FABdata[FABdata$surveyed_2019=="Yes",]
FABdata$C_estimate[which(is.na(FABdata$C_estimate))]<-0

## calculate individual overyielding by matching each
## individual to the monoculture mean of the same species
## in the same block
mono.means<-aggregate(C_estimate~species_code+block,
                      data=FABdata[FABdata$species_richness==1,],
                      FUN=mean)
FABdata$mono.means<-apply(FABdata,1,
                          function(x) {
                            mono.means$C_estimate[mono.means$species_code==x["species_code"] & mono.means$block==x["block"]]
                          })
## the amount by which the individual exceeds monoculture-based expectations
FABdata$ind.OY<-FABdata$C_estimate-FABdata$mono.means

## sum individual overyielding to the plot level
OY.agg<-aggregate(ind.OY~plot,data=FABdata,
                  FUN=function(x) sum(x,na.rm=T)*10/16)
OY.agg$species_richness<-plot_guide$Treatment[match(OY.agg$plot,plot_guide$Plot)]
OY.agg$block<-plot_guide$Block[match(OY.agg$plot,plot_guide$Plot)]
## overyielding should be 0 for monocultures (can be verified)
## we then proceed to delete rows corresponding to monocultures
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

## get estimates of total woody carbon per species per plot
C_sp_plot<-aggregate(C_estimate~species_code+plot+block,
                     data=FABdata_mod, FUN=function(x) sum(x)*10/16)

## generate indicators of species composition
FABplot_list<-split(FABdata_mod,f = FABdata_mod$plot)
FABplot_comp<-unlist(lapply(FABplot_list,function(plot) paste(unique(plot$species_code),collapse="|")))
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

## now we use pdiv to carry out the additive partitioning
C_partition<-addpart(C_estimate~sp_comp/species_code+plot,
                     fractions= ~fractions,
                     groups= ~block,
                     data=C_sp_plot)

## pdiv's calculations can be verified manually
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
## calculate proportions of conifers and AM species planted

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
## add light data from 2018 (just for context)

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
belowground$block<-plot_guide$Block[match(belowground$Plot,plot_guide$Plot)]

## estimate bulk density based on block means
## we calculate 0-20 cm as a weighted avg of 0-15 and 15-30
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
belowground$perC_change<-belowground$X..C_2019-belowground$X.C_2013
## convert to g C accumulated / cm^3 soil using BD
## 20 cm * 10000 cm^2 per m^2, 10000 m^2 per ha, 1000000 g per Mg
belowground$soilC<-(belowground$perC_change/100*belowground$BD)*200000*10000/1000000

## do the same for soil nitrogen, even though we don't use it
belowground$perN_change<-belowground$X..N_2019-belowground$X.N_2013
belowground$soilN<-(belowground$perN_change/100*belowground$BD)*200000*10000/1000000

## soil C pool at the beginning and the end end (not change)
belowground$soilC_2013<-(belowground$X.C_2013/100*belowground$BD)*200000*10000/1000000
belowground$soilC_2019<-(belowground$X..C_2019/100*belowground$BD)*200000*10000/1000000

## fine root carbon, assuming roots are 50% carbon
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

## attach monoculture values from the same species and block
sp_comp_long$mono_soil<-unlist(apply(sp_comp_long,1,function(x) {
  mono_block<-which(sp_comp_long$value>0.99 & sp_comp_long$variable==x["variable"] & sp_comp_long$block==x["block"])
  return(sp_comp_long$soilC_diff_plot[mono_block])
}))

sp_comp_long$mono_root<-unlist(apply(sp_comp_long,1,function(x) {
  mono_block<-which(sp_comp_long$value>0.99 & sp_comp_long$variable==x["variable"] & sp_comp_long$block==x["block"])
  return(sp_comp_long$rootC[mono_block])
}))

## multiply monoculture values by planted fractions of the species
sp_comp_long$mono_soil_exp<-sp_comp_long$mono_soil*sp_comp_long$value
sp_comp_long$mono_root_exp<-sp_comp_long$mono_root*sp_comp_long$value

## and add together within plots to get plot-level expectations
tot_mono_soil_exp<-aggregate(mono_soil_exp~Plot,data=sp_comp_long,FUN=sum)
tot_mono_root_exp<-aggregate(mono_root_exp~Plot,data=sp_comp_long,FUN=sum)

## subtract expectations from observed to get overyielding
belowground$mono_soil_exp<-tot_mono_soil_exp$mono_soil_exp[match(belowground$Plot,tot_mono_soil_exp$Plot)]
belowground$soil_OY<-belowground$soilC_diff_plot-belowground$mono_soil_exp
belowground$mono_root_exp<-tot_mono_root_exp$mono_root_exp[match(belowground$Plot,tot_mono_root_exp$Plot)]
belowground$root_OY<-belowground$rootC-belowground$mono_root_exp

belowground_sub<-belowground[-which(belowground$Sp.Richness==1),]
# plot(belowground_sub$soil_OY~belowground_sub$Sp.Richness)
# plot(belowground_sub$root_OY~belowground_sub$Sp.Richness)

###############################
## combine all data

C_agg$block<-plot_guide$Block[match(C_agg$plot,plot_guide$Plot)]

## AG woody OY, CE, and SE
C_agg$woodyOY<-OY.agg$ind.OY[match(C_agg$plot,OY.agg$plot)]
C_agg$woodyCE<-C_partition$CE.C_estimate[match(C_agg$plot,C_partition$plot)]
C_agg$woodySE<-C_partition$SE.C_estimate[match(C_agg$plot,C_partition$plot)]

## percentages of soil C and N at beginning and end
belowground_match<-match(C_agg$plot,belowground$Plot)
C_agg$perC_init<-belowground$X.C_2013[belowground_match]
C_agg$perC<-belowground$X..C_2019[belowground_match]
C_agg$perN_init<-belowground$X.N_2013[belowground_match]
C_agg$perN<-belowground$X..N_2019[belowground_match]

## bulk density at the block level
C_agg$BD<-belowground$BD[belowground_match]

## change in soil C at plot scale
C_agg$soilC<-belowground$soilC[belowground_match]

## total soil C at plot scale at beginning and end
C_agg$soilC_2013<-belowground$soilC_2013[belowground_match]
C_agg$soilC_2019<-belowground$soilC_2019[belowground_match]

## soil C overyielding
C_agg$soilOY<-belowground_sub$soil_OY[match(C_agg$plot,belowground_sub$Plot)]

## root C and overyielding and plot scale
C_agg$rootC<-belowground$rootC[belowground_match]
C_agg$rootOY<-belowground_sub$root_OY[match(C_agg$plot,belowground_sub$Plot)]

## adding together aboveground wood, soil, and fine roots
C_agg$totalC<-C_agg$woodyC+C_agg$soilC+C_agg$rootC
C_agg$totalOY<-C_agg$woodyOY+C_agg$soilOY+C_agg$rootOY

## macroaggregates
C_agg$macro250<-belowground$X..Mass.of.250[belowground_match]

## soil moisture
C_agg$soil_moisture<-belowground$Soil.moisture[belowground_match]

## planted proportions of AM and coniferous trees
C_agg$percentAM<-FAB_planted$percentAM[match(C_agg$plot,FAB_planted$plot)]
C_agg$percentCon<-FAB_planted$percentCon[match(C_agg$plot,FAB_planted$plot)]

## log (root C/AG woody C)
C_agg$logBA<-log(C_agg$rootC/C_agg$woodyC)

## functional and phylogenetic diversity
FD<-read.csv("OriginalData/ecy1958-sup-0003-tables1.csv")
C_agg$PSV<-FD$PSV[match(C_agg$plot,FD$Plot)]
C_agg$FDis<-FD$FDis[match(C_agg$plot,FD$Plot)]

## categorical leaf type and mycotype
## the 0.97 here is for plot 147
## which has one Tilia in it (by mistake)
C_agg$mycotype<-ifelse(C_agg$percentAM==0,"E",
                       ifelse(C_agg$percentAM>0.97,"A","B"))
C_agg$leaf_type<-ifelse(C_agg$percentCon==0,"D",
                        ifelse(C_agg$percentCon==1,"C","B"))

## optional step to set PSV to 0 in monocultures
C_agg$PSV[C_agg$species_richness==1]<-0

write.csv(C_agg,"ProcessedData/Cseq.csv",row.names = F)