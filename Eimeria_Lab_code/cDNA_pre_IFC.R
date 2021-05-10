library(dplyr)
library(tidyverse)

DNA <- read.csv("~/Documents/cDNA_done.csv")

OV14_17p <- read.csv("~/Mouse_Eimeria_Databasing/data/Eimeria_detection/HZ14-17_negatives.csv")
OV14_17n <- read.csv("~/Mouse_Eimeria_Databasing/data/Eimeria_detection/HZ14-17_positives.csv")
OV14_17 <- rbind(OV14_17n, OV14_17p)
OV14_17 <- select(OV14_17, Mouse_ID, Year, delta_ct_cewe_MminusE, eimeriaSpecies)
rm(OV14_17n, OV14_17p)
OV14_17$MC <- NA

OV18 <- read.csv("~/Mouse_Eimeria_Databasing/data/Eimeria_detection/Svenja/table_ct_and_more.csv")

names(OV18)[names(OV18) == "deltaCtMmE_tissue"] <- "delta_ct_cewe_MminusE"
names(OV18)[names(OV18) == "Eimeria.subspecies"] <- "eimeriaSpecies"
names(OV18)[names(OV18) == "Eimeria.presence.in.Caecum"] <- "MC"
OV18$Year <- 2018
OV18 <- select(OV18, Mouse_ID, Year, delta_ct_cewe_MminusE, eimeriaSpecies, MC)

OV14_18 <- rbind(OV14_17, OV18)

OV19 <- read.csv("~/Mouse_Eimeria_Databasing/data/Eimeria_detection/HZ19_CEWE_qPCR.csv")
names(OV19)[names(OV19) == "delta"] <- "delta_ct_cewe_MminusE"
OV19$Year <- 2019
OV19$eimeriaSpecies <- NA

OV19 <- select(OV19, Mouse_ID, Year, delta_ct_cewe_MminusE, eimeriaSpecies, MC)
t19 <- subset(OV19, OV19$MC == TRUE)

OV14_19 <- rbind(OV19, OV14_18)
rm(OV14_17, OV14_18, OV18, OV19)
OV14_19$MC[OV14_19$eimeriaSpecies == "E_ferrisi"] <- TRUE
OV14_19$MC[OV14_19$eimeriaSpecies == "E_falciformis"] <- TRUE
OV14_19$MC[OV14_19$eimeriaSpecies == "Double_tbd"] <- TRUE
OV14_19$MC[OV14_19$eimeriaSpecies == "Double_ferrisi_vermiformis"] <- TRUE
OV14_19$MC[OV14_19$eimeriaSpecies == "Other"] <- TRUE
OV14_19$MC[OV14_19$eimeriaSpecies == "Negative"] <- FALSE

OV14_19$eimeriaSpecies[OV14_19$eimeriaSpecies == "Double_ferrisi_vermiformis"] <- "Double"
OV14_19$eimeriaSpecies[OV14_19$eimeriaSpecies == "Double_tbd"] <- "Double"
OV14_19$eimeriaSpecies[OV14_19$eimeriaSpecies == "E_falciformis"] <- "E.falciformis"
OV14_19$eimeriaSpecies[OV14_19$eimeriaSpecies == "E_ferrisi"] <- "E.ferrisi"
OV14_19$eimeriaSpecies[OV14_19$eimeriaSpecies == "E. falciformis"] <- "E.falciformis"
OV14_19$eimeriaSpecies[OV14_19$eimeriaSpecies == "non infected"] <- "Negative"
OV14_19$eimeriaSpecies[OV14_19$eimeriaSpecies == "E. ferrisi"] <- "E.ferrisi"
OV14_19$eimeriaSpecies[OV14_19$eimeriaSpecies == "Eimeria sp."] <- "Other"

levels(factor(OV14_19$eimeriaSpecies))



processed <- merge(DNA, OV14_19, all.x =  T)

