# clean up HZ19 data
library(dplyr)

HZ <- read.csv("c:\\Documents and Settings/exemp/Downloads/EmanuelData.csv")
names(HZ)[names(HZ) == "X_Longit"] <- "Longitude"
names(HZ)[names(HZ) == "Y_Latit"] <- "Latitude"
names(HZ)[names(HZ) == "PIN"] <- "EH_ID"

names(HZ)[1:34]
HZ <- HZ[,1:34]

write.csv(x = HZ, "~/GitHub/Mouse_Eimeria_Field/data/Field_data/HZ10-19_Genotypes.csv")
