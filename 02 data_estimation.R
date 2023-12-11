setwd("C:/Users/querc/Dropbox/FABCarbonProject/")

library(reshape2)

## read in the cleaned FAB data
FABdata<-read.csv("ProcessedData/FAB_cleaned.csv")

###########################
## estimation models

## this next series of steps transforms the FAB data
## such that each individual x year is a row with DBH,
## basal diameter, and height among the columns
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

## here, we fit species-specific linear models that predict
## basal diameter based on dbh and height in the same year
## among these individual x year observations

## if basal diameter was measured for an individual x year,
## this function will just return the measured value
## otherwise it will return the model-estimated value
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

## this is the same as the above, but only using height
## rather than height and DBH to predict basal diameter

## the models produced using this approach do kind of look nicer
## although that 'niceness' might be illusory, due solely
## to the presence of smaller individuals that anchor the
## lower end of the regressions
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

## infer basal diameter using models that incorporate height and DBH
FABdata$diameter_2019_inf_dbh<-NA
for(i in 1:nrow(FABdata)){
  FABdata$diameter_2019_inf_dbh[i]<-with(FABdata,infer.diam.dbh(diameter_2019[i],
                                                                dbh_2019[i],
                                                                height_2019[i],
                                                                species_code[i],
                                                                df=FABmodel))
}

## infer basal diameter using models that incorporate height only
FABdata$diameter_2019_inf<-NA
for(i in 1:nrow(FABdata)){
  FABdata$diameter_2019_inf[i]<-with(FABdata,infer.diam(diameter_2019[i],
                                                        height_2019[i],
                                                        species_code[i],
                                                        df=FABmodel))
}

## calculate percent RMSD for models
# max<-range(FABmodel$diameter[FABmodel$species_code=="TIAM" & !is.na(FABmodel$dbh)],na.rm=T)[2]
# min<-range(FABmodel$diameter[FABmodel$species_code=="TIAM" & !is.na(FABmodel$dbh)],na.rm=T)[1]
# perRMSD<-sqrt(mean(lm(diameter~height,data=FABmodel[FABmodel$species_code=="TIAM",])$residuals^2))/(max-min)

##########################
## create a table of wood density (from Jenkins)
## and carbon content (from Lamlom and Savidge)

species_code<-c("ACNE","ACRU","BEPA","JUVI","PIBA","PIRE",
            "PIST","QUAL","QUEL","QUMA","QURU","TIAM")
density<-c(0.44,0.49,0.48,0.44,0.40,0.41,
             0.34,0.60,0.58,0.58,0.56,0.32)
C_content<-c(0.4934,0.4864,0.4837,0.5214,0.5040,0.5328,
             0.4974,0.4957,0.4963,0.4957,0.4963,0.4643)
wood_df<-data.frame(species_code=species_code,
                    density=density,
                    C_content=C_content)

####################
## estimating biomass for all species

## here we use the height-based estimates of unmeasured
## basal diameter when needed, but we can change the variable
## 'diameter_2019_inf' to 'diameter_2019_inf_dbh'
## to switch to models that use both height and DBH

FABdata$biomass_estimate<-NA
for(i in 1:nrow(FABdata)){
  
  ## Lambert's allometric equations for BEPA, JUVI, PIBA, PIRE, PIST
  ## previous versions of this script on GitHub include various intercomparisons
  ## of allometric equations to test which are the most consistent
  ## I've removed those here for readability
  
  if(FABdata$species_code[i]=="BEPA"){
    stem_wood<-0.0338*(FABdata$dbh_2019[i]/10)^2.0702*(FABdata$height_2019[i]/100)^0.6876
    stem_bark<-0.0080*(FABdata$dbh_2019[i]/10)^1.9754*(FABdata$height_2019[i]/100)^0.6659
    branch<-0.0257*(FABdata$dbh_2019[i]/10)^3.1754*(FABdata$height_2019[i]/100)^-0.9417
    FABdata$biomass_estimate[i]<-stem_wood+stem_bark+branch
    
    BEPAdensity<-wood_df$density[wood_df$species_code=="BEPA"]
    
    ## for trees that fall short of Lambert et al's size thresholds
    ## we calculate the conical volume and multiply it by wood density
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
  
  ## these species generally fall short of Lambert et al's size thresholds
  if(FABdata$species_code[i] %in% c("ACNE","ACRU","QUAL","QUEL","QUMA","QURU","TIAM")){
    sp_density<-wood_df$density[which(wood_df$species_code==FABdata$species_code[i])]
    FABdata$biomass_estimate[i]<-(1/3*FABdata$height_2019[i]*(FABdata$diameter_2019_inf_dbh[i]/20)^2*pi*sp_density)/1000
  }
}

## fill in estimated carbon using species specific wood carbon concentration
FABdata$C_content<-wood_df$C_content[match(FABdata$species_code,wood_df$species_code)]
FABdata$C_estimate<-FABdata$biomass_estimate*FABdata$C_content

write.csv(FABdata,"ProcessedData/FAB_Cestimate.csv")