setwd("C:/Users/querc/Dropbox/FABCarbonProject/")
FABdata<-read.csv("ProcessedData/FAB_edited.csv")

library(reshape2)
library(ggplot2)
library(pdiv) ## Pascal Niklaus's pdiv package

###########################
## estimation models

FABsub<-FABdata[,c("individual_id","position","species_code",
                   "dbh_2016","dbh_2017","dbh_2018","dbh_2019","dbh_2020",
                   "diameter_2016","diameter_2017","diameter_2018","diameter_2019","diameter_2020",
                   "height_2016","height_2017","height_2018","height_2019","height_2020")]
FABsub_long<-melt(FABsub,id.vars = c("individual_id","position","species_code"))

measurement_split<-strsplit(as.character(FABsub_long$variable),split="_")
FABsub_long$variable<-unlist(lapply(measurement_split,function(x) x[[1]]))
FABsub_long$year<-unlist(lapply(measurement_split,function(x) x[[2]]))

FABmodel<-dcast(individual_id+position+species_code+year~variable,
                data=FABsub_long,
                fun.aggregate = mean)

## infer basal diameter based on linear models that
## predict it from dbh and height in the same year,
## using all individual x years when all three were measured
## NOTE: if basal diameter was measured in the year to be predicted,
## this function will just return the measured value rather than
## estimating it from the model
infer.diam.dbh<-function(diameter_ind,dbh_ind,height_ind,species_code,df){
  if(is.na(height_ind)){
    return(NA)
  }
  if(is.na(diameter_ind) & is.na(dbh_ind)){
    return(NA)
  }
  if(!is.na(diameter_ind)){
    return(diameter_ind)
  }
  model<-lm(diameter~dbh+height,data=df[df$species_code==species_code,])
  return(predict(model,newdata=data.frame(dbh=dbh_ind,height=height_ind)))
}

## infer basal diameter based on linear models that
## predict it from just height in the same year,
## using all individual x years when all three were measured
## the models produced using this approach do kind of look nicer
## although that 'niceness' might be illusory, due solely
## to the presence of smaller individuals that anchor the
## lower end of the regressions

## NOTE: if basal diameter was measured in the year to be predicted,
## this function will just return the measured value rather than
## estimating it from the model
infer.diam<-function(diameter_ind,height_ind,species_code,df){
  if(is.na(height_ind)){
    return(NA)
  }
  if(!is.na(diameter_ind)){
    return(diameter_ind)
  }
  model<-lm(diameter~height,data=df[df$species_code==species_code,])
  return(predict(model,newdata=data.frame(height=height_ind)))
}

## infer basal diameter using models that incorporate height only
FABdata$diameter_2019_inf<-NA
for(i in 1:nrow(FABdata)){
  FABdata$diameter_2019_inf[i]<-with(FABdata,infer.diam(diameter_2019[i],
                                                        height_2019[i],
                                                        species_code[i],
                                                        df=FABmodel))
}

## infer basal diameter using models that incorporate height and DBH
FABdata$diameter_2019_inf_dbh<-NA
for(i in 1:nrow(FABdata)){
  FABdata$diameter_2019_inf_dbh[i]<-with(FABdata,infer.diam.dbh(diameter_2019[i],
                                                                dbh_2019[i],
                                                                height_2019[i],
                                                                species_code[i],
                                                                df=FABmodel))
}

## calculate percent RMSD for models
# max<-range(FABmodel$diameter[FABmodel$species_code=="TIAM" & !is.na(FABmodel$dbh)],na.rm=T)[2]
# min<-range(FABmodel$diameter[FABmodel$species_code=="TIAM" & !is.na(FABmodel$dbh)],na.rm=T)[1]
# perRMSD<-sqrt(mean(lm(diameter~height,data=FABmodel[FABmodel$species_code=="TIAM",])$residuals^2))/(max-min)

##########################
## wood density from Jenkins
## carbon content from Lamlom and Savidge

species_code<-c("ACNE","ACRU","BEPA","JUVI","PIBA","PIRE",
            "PIST","QUAL","QUEL","QUMA","QURU","TIAM")
density<-c(0.44,0.49,0.48,0.44,0.40,0.41,
             0.34,0.60,0.58,0.58,0.56,0.32)
C_content<-c(0.4934,0.4864,0.4837,0.5214,0.5040,0.5328,
             0.4974,0.4957,0.4963,0.4957,0.4963,0.4643)
wood_df<-data.frame(species_code=species_code,
                    density=density,
                    C_content=C_content)

##########################
## allometric equations for BEPA, JUVI, PIBA, PIRE, PIST

#######
## BEPA

## Fatemi PV
# FABdata$BEPAestimate1<-NA
# for(i in 1:nrow(FABdata)){
# #  if(!is.numeric(FABdata$dbh_2019[i])){FABdata$BEPAestimate1[i]<-NA}
#   PV<-0.5*pi*(FABdata$dbh_2019[i]/20)^2*FABdata$height_2019[i]
#   stem_wood<- 10^(-0.753+1.096*log10(PV))
#   stem_bark<- 10^(-1.627+1.010*log10(PV))
#   FABdata$BEPAestimate1[i]<-(stem_wood+stem_bark)/1000
#   if((is.na(FABdata$dbh_2019[i]) || FABdata$dbh_2019[i] < 20) & !is.na(FABdata$diameter_2019[i])){
#     FABdata$BEPAestimate1[i]<-(PV*BEPAdensity)/1000
#   }
#   if(FABdata$species_code[i]!="BEPA"){FABdata$BEPAestimate1[i]<-NA}
# }

## Fatemi DBH
# FABdata$BEPAestimate1X<-NA
# for(i in 1:nrow(FABdata)){
#   #  if(!is.numeric(FABdata$dbh_2019[i])){FABdata$BEPAestimate1X[i]<-NA}
#   stem_wood<- 10^(1.739+2.638*log10(FABdata$dbh_2019[i]/10))
#   stem_bark<- 10^(0.823+2.711*log10(FABdata$dbh_2019[i]/10))
#   FABdata$BEPAestimate1X[i]<-(stem_wood+stem_bark)/1000
#   if((is.na(FABdata$dbh_2019[i]) || FABdata$dbh_2019[i] < 20) & !is.na(FABdata$diameter_2019[i])){
#     FABdata$BEPAestimate1X[i]<-(0.5*pi*(FABdata$diameter_2019[i]/20)^2*FABdata$height_2019[i]*BEPAdensity)/1000
#   }
#   if(FABdata$species_code[i]!="BEPA"){FABdata$BEPAestimate1X[i]<-NA}
# }

## Lambert
FABdata$BEPAestimate2<-NA
for(i in 1:nrow(FABdata)){
  #  if(!is.numeric(FABdata$dbh_2019[i])){FABdata$BEPAestimate2[i]<-NA}
  stem_wood<-0.0338*(FABdata$dbh_2019[i]/10)^2.0702*(FABdata$height_2019[i]/100)^0.6876
  stem_bark<-0.0080*(FABdata$dbh_2019[i]/10)^1.9754*(FABdata$height_2019[i]/100)^0.6659
  branch<-0.0257*(FABdata$dbh_2019[i]/10)^3.1754*(FABdata$height_2019[i]/100)^-0.9417
  FABdata$BEPAestimate2[i]<-stem_wood+stem_bark+branch
  BEPAdensity<-wood_df$density[wood_df$species_code=="BEPA"]
  if(is.na(FABdata$dbh_2019[i]) || FABdata$dbh_2019[i] < 15 || FABdata$height_2019[i] < 260){
    FABdata$BEPAestimate2[i]<-(1/3*FABdata$height_2019[i]*(FABdata$diameter_2019_inf_dbh[i]/20)^2*pi*BEPAdensity)/1000
  }
  if(FABdata$species_code[i]!="BEPA"){FABdata$BEPAestimate2[i]<-NA}
}

## Ter-Mikaelian 1997 -- Ker 1984
# FABdata$BEPAestimate3<-NA
# for(i in 1:nrow(FABdata)){
#   #  if(!is.numeric(FABdata$dbh_2019[i])){FABdata$BEPAestimate3[i]<-NA}
#   FABdata$BEPAestimate3[i]<-0.0847*(FABdata$dbh_2019[i]/10)^2.4029
#   if(is.na(FABdata$dbh_2019[i]) & !is.na(FABdata$diameter_2019[i])){
#     FABdata$BEPAestimate3[i]<-(1/3*FABdata$height_2019[i]*(FABdata$diameter_2019[i]/20)^2*pi*BEPAdensity)/1000
#   }
#   if(FABdata$species_code[i]!="BEPA"){FABdata$BEPAestimate3[i]<-NA}
# }
# 
# ## Ter-Mikaelian 1997 -- Schmitt and Grigal
# FABdata$BEPAestimate4<-NA
# for(i in 1:nrow(FABdata)){
#   #  if(!is.numeric(FABdata$dbh_2019[i])){FABdata$BEPAestimate4[i]<-NA}
#   FABdata$BEPAestimate4[i]<-0.0923*(FABdata$dbh_2019[i]/10)^2.4800
#   if(is.na(FABdata$dbh_2019[i]) & !is.na(FABdata$diameter_2019[i])){
#     FABdata$BEPAestimate4[i]<-(1/3*FABdata$height_2019[i]*(FABdata$diameter_2019[i]/20)^2*pi*BEPAdensity)/1000
#   }
#   if(FABdata$species_code[i]!="BEPA"){FABdata$BEPAestimate4[i]<-NA}
# }

#######
## JUVI

## Lambert
FABdata$JUVIestimate1<-NA
for(i in 1:nrow(FABdata)){
  #  if(!is.numeric(FABdata$dbh_2019[i])){FABdata$JUVIestimate1[i]<-NA}
  stem_wood<-0.1277*(FABdata$dbh_2019[i]/10)^1.9778
  stem_bark<-0.0377*(FABdata$dbh_2019[i]/10)^1.6064
  branch<-0.0254*(FABdata$dbh_2019[i]/10)^2.2884
  FABdata$JUVIestimate1[i]<-stem_wood+stem_bark+branch
  JUVIdensity<-wood_df$density[wood_df$species_code=="JUVI"]
  if(is.na(FABdata$dbh_2019[i]) || FABdata$dbh_2019[i] < 20 || FABdata$height_2019[i] < 270){
    FABdata$JUVIestimate1[i]<-(1/3*FABdata$height_2019[i]*(FABdata$diameter_2019_inf_dbh[i]/20)^2*pi*JUVIdensity)/1000
  }
  if(FABdata$species_code[i]!="JUVI"){FABdata$JUVIestimate1[i]<-NA}
}

#######
## PIBA

# ## Ter-Mikaelian 1997 -- Ker 1984
# FABdata$PIBAestimate1<-NA
# for(i in 1:nrow(FABdata)){
#   #  if(!is.numeric(FABdata$dbh_2019[i])){FABdata$PIBAestimate1[i]<-NA}
#   FABdata$PIBAestimate1[i]<-0.1470*(FABdata$dbh_2019[i]/10)^2.1673
#   if(is.na(FABdata$dbh_2019[i]) & !is.na(FABdata$diameter_2019[i])){
#     FABdata$PIBAestimate1[i]<-(1/3*FABdata$height_2019[i]*(FABdata$diameter_2019[i]/20)^2*pi*PIBAdensity)/1000
#   }
#   if(FABdata$species_code[i]!="PIBA"){FABdata$PIBAestimate1[i]<-NA}
# }

## Lambert
FABdata$PIBAestimate2<-NA
for(i in 1:nrow(FABdata)){
  #  if(!is.numeric(FABdata$dbh_2019[i])){FABdata$PIBAestimate2[i]<-NA}
  stem_wood<-0.0199*(FABdata$dbh_2019[i]/10)^1.6883*(FABdata$height_2019[i]/100)^1.2456
  stem_bark<-0.0141*(FABdata$dbh_2019[i]/10)^1.5994*(FABdata$height_2019[i]/100)^0.5957
  branch<-0.0185*(FABdata$dbh_2019[i]/10)^3.0584*(FABdata$height_2019[i]/100)^-0.9816
  FABdata$PIBAestimate2[i]<-stem_wood+stem_bark+branch
  PIBAdensity<-wood_df$density[wood_df$species_code=="PIBA"]
  if(is.na(FABdata$dbh_2019[i]) || FABdata$dbh_2019[i] < 17 || FABdata$height_2019[i] < 290){
    FABdata$PIBAestimate2[i]<-(1/3*FABdata$height_2019[i]*(FABdata$diameter_2019_inf_dbh[i]/20)^2*pi*PIBAdensity)/1000
  }
  if(FABdata$species_code[i]!="PIBA"){FABdata$PIBAestimate2[i]<-NA}
}

#######
## PIRE

# ## Ter-Mikaelian 1997 -- Ker 1980
# ## seems like it overestimates
# FABdata$PIREestimate1<-NA
# for(i in 1:nrow(FABdata)){
#   #  if(!is.numeric(FABdata$dbh_2019[i])){FABdata$PIREestimate1[i]<-NA}
#   FABdata$PIREestimate1[i]<-0.0586*(FABdata$dbh_2019[i]/10)^2.3892
#   if(is.na(FABdata$dbh_2019[i]) & !is.na(FABdata$diameter_2019[i])){
#     FABdata$PIREestimate1[i]<-(1/3*FABdata$height_2019[i]*(FABdata$diameter_2019[i]/20)^2*pi*PIREdensity)/1000
#   }
#   if(FABdata$species_code[i]!="PIRE"){FABdata$PIREestimate1[i]<-NA}
# }

## Lambert
FABdata$PIREestimate2<-NA
for(i in 1:nrow(FABdata)){
  #  if(!is.numeric(FABdata$dbh_2019[i])){FABdata$PIREestimate2[i]<-NA}
  # stem_wood<-0.0106*(FABdata$dbh_2019[i]/10)^1.7725*(FABdata$height_2019[i]/100)^1.3285
  # stem_bark<-0.0277*(FABdata$dbh_2019[i]/10)^1.5192*(FABdata$height_2019[i]/100)^0.4645
  # branch<-0.0125*(FABdata$dbh_2019[i]/10)^3.3865*(FABdata$height_2019[i]/100)^-1.1939
  stem_wood<-0.0564*(FABdata$dbh_2019[i]/10)^2.4465
  stem_bark<-0.0188*(FABdata$dbh_2019[i]/10)^2.0527
  branch<-0.0033*(FABdata$dbh_2019[i]/10)^2.7515
  FABdata$PIREestimate2[i]<-stem_wood+stem_bark+branch
  PIREdensity<-wood_df$density[wood_df$species_code=="PIRE"]
  if(is.na(FABdata$dbh_2019[i]) || FABdata$dbh_2019[i] < 13 || FABdata$height_2019[i] < 180){
    FABdata$PIREestimate2[i]<-(1/3*FABdata$height_2019[i]*(FABdata$diameter_2019_inf_dbh[i]/20)^2*pi*PIREdensity)/1000
  }
  if(FABdata$species_code[i]!="PIRE"){FABdata$PIREestimate2[i]<-NA}
}

#######
## PIST

## Ter-Mikaelian 1997 -- Ker 1980
# FABdata$PISTestimate1<-NA
# for(i in 1:nrow(FABdata)){
#   #  if(!is.numeric(FABdata$dbh_2019[i])){FABdata$PISTestimate1[i]<-NA}
#   FABdata$PISTestimate1[i]<-0.0414*(FABdata$dbh_2019[i]/10)^2.5360
#   if(is.na(FABdata$dbh_2019[i]) & !is.na(FABdata$diameter_2019[i])){
#     FABdata$PISTestimate1[i]<-(1/3*FABdata$height_2019[i]*(FABdata$diameter_2019[i]/20)^2*pi*PISTdensity)/1000
#   }
#   if(FABdata$species_code[i]!="PIST"){FABdata$PISTestimate1[i]<-NA}
# }

## Lambert
FABdata$PISTestimate2<-NA
for(i in 1:nrow(FABdata)){
  #  if(!is.numeric(FABdata$dbh_2019[i])){FABdata$PISTestimate2[i]<-NA}
  # stem_wood<-0.0170*(FABdata$dbh_2019[i]/10)^1.7779*(FABdata$height_2019[i]/100)^1.1370
  # stem_bark<-0.0069*(FABdata$dbh_2019[i]/10)^1.6589*(FABdata$height_2019[i]/100)^0.9582
  # branch<-0.0184*(FABdata$dbh_2019[i]/10)^3.1968*(FABdata$height_2019[i]/100)^-1.0876
  stem_wood<-0.0997*(FABdata$dbh_2019[i]/10)^2.2709
  stem_bark<-0.0192*(FABdata$dbh_2019[i]/10)^2.2038
  branch<-0.0056*(FABdata$dbh_2019[i]/10)^2.6011
  FABdata$PISTestimate2[i]<-stem_wood+stem_bark+branch
  PISTdensity<-wood_df$density[wood_df$species_code=="PIST"]
  if(is.na(FABdata$dbh_2019[i]) || FABdata$dbh_2019[i] < 15 || FABdata$height_2019[i] < 230){
    FABdata$PISTestimate2[i]<-(1/3*FABdata$height_2019[i]*(FABdata$diameter_2019_inf_dbh[i]/20)^2*pi*PISTdensity)/1000
  }
  if(FABdata$species_code[i]!="PIST"){FABdata$PISTestimate2[i]<-NA}
}

####################
## all together now
## can change diameter_2019_inf to diameter_2019_inf_dbh
## to switch to models that use both height and dbh

FABdata$biomass_estimate<-NA
for(i in 1:nrow(FABdata)){
  ## Lambert equations
  if(FABdata$species_code[i]=="BEPA"){
    stem_wood<-0.0338*(FABdata$dbh_2019[i]/10)^2.0702*(FABdata$height_2019[i]/100)^0.6876
    stem_bark<-0.0080*(FABdata$dbh_2019[i]/10)^1.9754*(FABdata$height_2019[i]/100)^0.6659
    branch<-0.0257*(FABdata$dbh_2019[i]/10)^3.1754*(FABdata$height_2019[i]/100)^-0.9417
    FABdata$biomass_estimate[i]<-stem_wood+stem_bark+branch
    BEPAdensity<-wood_df$density[wood_df$species_code=="BEPA"]
    if(is.na(FABdata$dbh_2019[i]) || FABdata$dbh_2019[i] < 15 || FABdata$height_2019[i] < 260){
      FABdata$biomass_estimate[i]<-(1/3*FABdata$height_2019[i]*(FABdata$diameter_2019_inf_dbh[i]/20)^2*pi*BEPAdensity)/1000
    }
  }
  if(FABdata$species_code[i]=="JUVI"){
    stem_wood<-0.1277*(FABdata$dbh_2019[i]/10)^1.9778
    stem_bark<-0.0377*(FABdata$dbh_2019[i]/10)^1.6064
    branch<-0.0254*(FABdata$dbh_2019[i]/10)^2.2884
    FABdata$biomass_estimate[i]<-stem_wood+stem_bark+branch
    JUVIdensity<-wood_df$density[wood_df$species_code=="JUVI"]
    if(is.na(FABdata$dbh_2019[i]) || FABdata$dbh_2019[i] < 20 || FABdata$height_2019[i] < 270){
      FABdata$biomass_estimate[i]<-(1/3*FABdata$height_2019[i]*(FABdata$diameter_2019_inf_dbh[i]/20)^2*pi*JUVIdensity)/1000
    }
  }
  if(FABdata$species_code[i]=="PIBA"){
    stem_wood<-0.0199*(FABdata$dbh_2019[i]/10)^1.6883*(FABdata$height_2019[i]/100)^1.2456
    stem_bark<-0.0141*(FABdata$dbh_2019[i]/10)^1.5994*(FABdata$height_2019[i]/100)^0.5957
    branch<-0.0185*(FABdata$dbh_2019[i]/10)^3.0584*(FABdata$height_2019[i]/100)^-0.9816
    FABdata$biomass_estimate[i]<-stem_wood+stem_bark+branch
    PIBAdensity<-wood_df$density[wood_df$species_code=="PIBA"]
    if(is.na(FABdata$dbh_2019[i]) || FABdata$dbh_2019[i] < 17 || FABdata$height_2019[i] < 290){
      FABdata$biomass_estimate[i]<-(1/3*FABdata$height_2019[i]*(FABdata$diameter_2019_inf_dbh[i]/20)^2*pi*PIBAdensity)/1000
    }
  }
  if(FABdata$species_code[i]=="PIRE"){
    stem_wood<-0.0564*(FABdata$dbh_2019[i]/10)^2.4465
    stem_bark<-0.0188*(FABdata$dbh_2019[i]/10)^2.0527
    branch<-0.0033*(FABdata$dbh_2019[i]/10)^2.7515
    FABdata$biomass_estimate[i]<-stem_wood+stem_bark+branch
    PIREdensity<-wood_df$density[wood_df$species_code=="PIRE"]
    if(is.na(FABdata$dbh_2019[i]) || FABdata$dbh_2019[i] < 13 || FABdata$height_2019[i] < 180){
      FABdata$biomass_estimate[i]<-(1/3*FABdata$height_2019[i]*(FABdata$diameter_2019_inf_dbh[i]/20)^2*pi*PIREdensity)/1000
    }
  }
  if(FABdata$species_code[i]=="PIST"){
    stem_wood<-0.0997*(FABdata$dbh_2019[i]/10)^2.2709
    stem_bark<-0.0192*(FABdata$dbh_2019[i]/10)^2.2038
    branch<-0.0056*(FABdata$dbh_2019[i]/10)^2.6011
    FABdata$biomass_estimate[i]<-stem_wood+stem_bark+branch
    PISTdensity<-wood_df$density[wood_df$species_code=="PIST"]
    if(is.na(FABdata$dbh_2019[i]) || FABdata$dbh_2019[i] < 15 || FABdata$height_2019[i] < 230){
      FABdata$biomass_estimate[i]<-(1/3*FABdata$height_2019[i]*(FABdata$diameter_2019_inf_dbh[i]/20)^2*pi*PISTdensity)/1000
    }
  }
  if(FABdata$species_code[i] %in% c("ACNE","ACRU","QUAL","QUEL","QUMA","QURU","TIAM")){
    sp_density<-wood_df$density[which(wood_df$species_code==FABdata$species_code[i])]
    FABdata$biomass_estimate[i]<-(1/3*FABdata$height_2019[i]*(FABdata$diameter_2019_inf_dbh[i]/20)^2*pi*sp_density)/1000
  }
}

## purely hardwood vs. softwood
# FABdata$biomass_estimate<-NA
# for(i in 1:nrow(FABdata)){
#   ## Lambert equations
#   if(FABdata$species_code[i] %in% c("JUVI","PIBA","PIRE","PIST")){
#     # stem_wood<-0.0284*(FABdata$dbh_2019[i]/10)^1.6894*(FABdata$height_2019[i]/100)^1.0857
#     # stem_bark<-0.0100*(FABdata$dbh_2019[i]/10)^1.8463*(FABdata$height_2019[i]/100)^0.5616
#     # branch<-0.0301*(FABdata$dbh_2019[i]/10)^3.0038*(FABdata$height_2019[i]/100)^-1.0520
#     stem_wood<-0.0648*(FABdata$dbh_2019[i]/10)^2.3927
#     stem_bark<-0.0162*(FABdata$dbh_2019[i]/10)^2.1959
#     branch<-0.0156*(FABdata$dbh_2019[i]/10)^2.2916
#     FABdata$biomass_estimate[i]<-stem_wood+stem_bark+branch
#     density<-wood_df$density[wood_df$species_code==FABdata$species_code[i]]
#     if(is.na(FABdata$dbh_2019[i]) || FABdata$dbh_2019[i] < 13 || FABdata$height_2019[i] < 180){
#       FABdata$biomass_estimate[i]<-(1/3*FABdata$height_2019[i]*(FABdata$diameter_2019_inf[i]/20)^2*pi*density)/1000
#     }
#   }
#   if(FABdata$species_code[i] %in% c("ACNE","ACRU","BEPA","QUAL","QUEL","QUMA","QURU","TIAM")){
#     stem_wood<-0.0359*(FABdata$dbh_2019[i]/10)^2.0263*(FABdata$height_2019[i]/100)^0.6987
#     stem_bark<-0.0094*(FABdata$dbh_2019[i]/10)^1.8677*(FABdata$height_2019[i]/100)^0.6985
#     branch<-0.0433*(FABdata$dbh_2019[i]/10)^2.6817*(FABdata$height_2019[i]/100)^-0.5731
#     FABdata$biomass_estimate[i]<-stem_wood+stem_bark+branch
#     density<-wood_df$density[which(wood_df$species_code==FABdata$species_code[i])]
#     if(is.na(FABdata$dbh_2019[i]) || FABdata$dbh_2019[i] < 11 || FABdata$height_2019[i] < 250){
#       FABdata$biomass_estimate[i]<-(1/3*FABdata$height_2019[i]*(FABdata$diameter_2019_inf[i]/20)^2*pi*density)/1000
#     }
#   }
# }

## fill in estimated carbon using species specific wood carbon concentration
FABdata$C_content<-wood_df$C_content[match(FABdata$species_code,wood_df$species_code)]
FABdata$C_estimate<-FABdata$biomass_estimate*FABdata$C_content

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

## for calculations to work, fix issues about misplanted species
## in most of these plots, there's only one, but plot 31 had more severe issues
FABdata_mod<-FABdata
FABdata_mod$species_code[FABdata_mod$species_code=="ACRU" & FABdata_mod$plot==18]<-"ACNE"
# FABdata_mod<-FABdata_mod[-which(FABdata_mod$plot==31),]
FABdata_mod$species_code[FABdata_mod$species_code=="PIST" & FABdata_mod$plot==34]<-"PIRE"
FABdata_mod$species_code[FABdata_mod$species_code=="PIST" & FABdata_mod$plot==37]<-"PIRE"
FABdata_mod$species_code[FABdata_mod$species_code=="TIAM" & FABdata_mod$plot==147]<-"ACRU"

## generate indicators of species composition
FABplot_list<-split(FABdata_mod,f = FABdata_mod$plot)
FABplot_comp<-unlist(lapply(FABplot_list,function(plot) paste(unique(plot$species_code),collapse="|")))

## get estimates of total woody carbon per species per plot
C_sp_plot<-aggregate(C_estimate~species_code+plot+block,
                     data=FABdata_mod,
                     FUN=sum)
C_sp_plot$sp_comp<-FABplot_comp[match(C_sp_plot$plot,names(FABplot_comp))]

## estimate fractions
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