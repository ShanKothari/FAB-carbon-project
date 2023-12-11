## this script represents my best attempt at
## cleaning the FAB survey data

#######################################
## biomass estimates from 2019 data

setwd("C:/Users/querc/Dropbox/FABCarbonProject/")
FAB_CL<-read.csv("OriginalData/fab1_growth_CL.csv")

########################
## fixing apparent errors

#######
## ACNE

FAB_CL$diameter_2019[FAB_CL$position==4226]<-16 ## off by a cm?

## position 167 grew too much between 2019 and 2020?

#######
## ACRU

## something is up with these individuals, but hard to tell what exactly
## all have basal diameter < dbh in 2016 and dbh 2016 > dbh 2017
## I'm going to assume basal diameter and dbh were switched
FAB_CL$diameter_2016[FAB_CL$position==75]<-22
FAB_CL$diameter_2016[FAB_CL$position==2033]<-14
FAB_CL$diameter_2016[FAB_CL$position==7220]<-15
FAB_CL$dbh_2016[FAB_CL$position==75]<-15.5
FAB_CL$dbh_2016[FAB_CL$position==2033]<-8.5 ## 8.3 seems impossible?
FAB_CL$dbh_2016[FAB_CL$position==7220]<-13.5

FAB_CL$diameter_2020[FAB_CL$position==118]<-20 ## cm rather than mm?

FAB_CL$diameter_2020[FAB_CL$position==5876]<-7.5 ## read as 1.5 rather than 7.5?

## switched D and DBH?
FAB_CL$diameter_2018[FAB_CL$position==2561]<-8
FAB_CL$dbh_2018[FAB_CL$position==2561]<-NA

## position 3023? DBH in  2018 -- too short

#######
## BEPA

FAB_CL$dbh_2018[FAB_CL$position==6828]<-16 ## off by a cm?

FAB_CL$diameter_2018[FAB_CL$position==560]<-NA ## basal diameter doesn't make sense

FAB_CL$dbh_2020[FAB_CL$position==6805]<-21 ## off by a cm?

## not sure what's going on here, just interpolated between 2018 and 2020
FAB_CL$dbh_2019[FAB_CL$position==4805]<-44

FAB_CL$dbh_2019[FAB_CL$position==1132]<-19 ## missing, interpolated

#######
## JUVI

FAB_CL$dbh_2018[FAB_CL$position==2029]<-17.5 ## off by a cm?

FAB_CL$dbh_2019[FAB_CL$position==1029]<-14.5 ## missing, assumed constant

#######
## PIBA

FAB_CL$height_2019[FAB_CL$position==3213]<-392 ## off by 100 cm?
FAB_CL$dbh_2018[FAB_CL$position==5921]<-29 ## off by 2 cm?
FAB_CL$dbh_2018[FAB_CL$position==7435]<-22.5 ## off by a cm?

## missing, assumed constant
FAB_CL$diameter_2019[FAB_CL$position==5901]<-14

## marked as alive but without measurements--so interpolating
FAB_CL$height_2019[FAB_CL$position==3608]<-241
FAB_CL$dbh_2019[FAB_CL$position==3608]<-19

#######
## PIRE

FAB_CL$dbh_2018[FAB_CL$position==1142]<-24.5 ## off by two cm?

## PIRE 7153 seems wrong in 2020 but not sure how to fix

## switched
FAB_CL$diameter_2019[FAB_CL$position==7317]<-NA
FAB_CL$dbh_2019[FAB_CL$position==7317]<-2.5

## 8858 seems off between 2019 and 2020

#######
## PIST

FAB_CL$diameter_2018[FAB_CL$position==244]<-15.5 ## off by a cm?

## switched?
FAB_CL$diameter_2017[FAB_CL$position==630]<-NA
FAB_CL$dbh_2017[FAB_CL$position==630]<-12
FAB_CL$diameter_2018[FAB_CL$position==630]<-NA
FAB_CL$dbh_2018[FAB_CL$position==630]<-20.5

## switched?
FAB_CL$diameter_2019[FAB_CL$position==895]<-NA
FAB_CL$dbh_2019[FAB_CL$position==895]<-7.5

FAB_CL$dbh_2018[FAB_CL$position==3552]<-23.5 ## off by 2 cm?

FAB_CL$dbh_2020[FAB_CL$position==373]<-52 ## off by 4 cm?

FAB_CL$height_2018[FAB_CL$position==1846]<-236 ## off by a meter?

FAB_CL$height_2018[FAB_CL$position==1216]<-277 ## bad handwriting?

FAB_CL$height_2018[FAB_CL$position==3160]<-313 ## off by a meter?

FAB_CL$height_2019[FAB_CL$position==198]<-331 ## off by a meter?

## missing and interpolated
FAB_CL$dbh_2019[FAB_CL$position==2574]<-32
FAB_CL$diameter_2019[FAB_CL$individual_id=="PIST-2015-8105"]<-32.5
FAB_CL$dbh_2019[FAB_CL$individual_id=="PIST-2015-8105"]<-14.5

#######
## QUAL

FAB_CL$deadmissing_2019[FAB_CL$position==1904]<-"Yes"

## switched
FAB_CL$diameter_2018[FAB_CL$position==1663]<-NA
FAB_CL$dbh_2018[FAB_CL$position==1663]<-18.5

FAB_CL$diameter_2019[FAB_CL$position==8250]<-15 ## off by a cm?

FAB_CL$diameter_2018[FAB_CL$position==6264]<-17 ## off by a cm?

## this individual died before
FAB_CL$surveyed_2019[FAB_CL$individual_id=="QUAL-2015-1904"]<-"No"

#######
## QUEL

## looks like duplicated from basal diameter; interpolated
FAB_CL$dbh_2017[FAB_CL$position==1909]<-3.5

FAB_CL$diameter_2019[FAB_CL$position==6423]<-10 ## cm rather than mm?

FAB_CL$diameter_2019[FAB_CL$position==5095]<-14 ## off by a cm?

FAB_CL$diameter_2019[FAB_CL$position==6304]<-15 ## off by a cm?

FAB_CL$diameter_2020[FAB_CL$position==1311]<-15 ## off by a cm?

FAB_CL$diameter_2018[FAB_CL$position==1313]<-17 ## off by a cm?

FAB_CL$dbh_2018[FAB_CL$position==4980]<-NA ## can't have DBH
FAB_CL$dbh_2018[FAB_CL$position==5101]<-NA ## can't have DBH

## plot 105?

#######
## QUMA

FAB_CL$height_2019[FAB_CL$position==6724]<-101.5 ## missing digit?

FAB_CL$diameter_2020[FAB_CL$position==6366]<-24 ## off by 2 cm?

FAB_CL$dbh_2020[FAB_CL$position==8660]<-NA ## can't have DBH

FAB_CL$diameter_2019[FAB_CL$position==6042]<-16.5 ## missing; interpolated

## unsure about this
FAB_CL$height_2019[FAB_CL$position==3107]<-102 ## can't tell what happened here, but perhaps misread?

#######
## TIAM

## switched?
FAB_CL$diameter_2018[FAB_CL$position==4591]<-14
FAB_CL$dbh_2018[FAB_CL$position==4591]<-3

FAB_CL$dbh_2018[FAB_CL$position==1143]<-16.5 ## off by a cm?

FAB_CL$dbh_2020[FAB_CL$position==2825]<-16 ## off by a cm?

FAB_CL$height_2019[FAB_CL$position==6545]<-211 ## off by a meter?

#############################
## write data

write.csv(FAB_CL,"ProcessedData/FAB_cleaned.csv",row.names = F)
