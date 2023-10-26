setwd("C:/Users/querc/Dropbox/FABCarbonProject/")

library(reshape2)

FABdata<-read.csv("ProcessedData/FAB_edited.csv")

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

write.csv(FABdata,"ProcessedData/FAB_Cestimate.csv")