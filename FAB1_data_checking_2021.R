setwd("C:/Users/querc/Dropbox/FABCarbonProject/")
FAB_CL<-read.csv("OriginalData/fab1_growth_CL.csv")

with(FAB_CL,which(height_2018>250 & diameter_2018<22 & species_code=="BEPA"))

with(FAB_CL[FAB_CL$species_code=="BEPA",],plot(dbh_2018~diameter_2018))
with(FAB_CL[FAB_CL$species_code=="BEPA",],plot(dbh_2019~diameter_2019))
with(FAB_CL[FAB_CL$species_code=="BEPA",],plot(dbh_2020~diameter_2020))
with(FAB_CL[FAB_CL$species_code=="BEPA",],plot(diameter_2018~height_2018))
with(FAB_CL[FAB_CL$species_code=="BEPA",],plot(diameter_2019~height_2019))
with(FAB_CL[FAB_CL$species_code=="BEPA",],plot(diameter_2020~height_2020))
with(FAB_CL[FAB_CL$species_code=="BEPA",],plot(dbh_2018~height_2018))
with(FAB_CL[FAB_CL$species_code=="BEPA",],plot(dbh_2019~height_2019))
with(FAB_CL[FAB_CL$species_code=="BEPA",],plot(dbh_2020~height_2020))

with(FAB_CL,which(dbh_2019>22 & dbh_2018<13 & species_code=="BEPA"))

with(FAB_CL[FAB_CL$species_code=="BEPA",],plot(height_2019~height_2018))
with(FAB_CL[FAB_CL$species_code=="BEPA",],plot(height_2020~height_2019))
with(FAB_CL[FAB_CL$species_code=="BEPA",],plot(diameter_2019~diameter_2018))
with(FAB_CL[FAB_CL$species_code=="BEPA",],plot(diameter_2020~diameter_2019))
with(FAB_CL[FAB_CL$species_code=="BEPA",],plot(dbh_2019~dbh_2018))
with(FAB_CL[FAB_CL$species_code=="BEPA",],plot(dbh_2020~dbh_2019))
