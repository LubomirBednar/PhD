set.seed(123)
# get the samples we already have
DNA <- read.csv("~/Documents/cDNA_done.csv")
################## make a big table of all the samples and remove the ones we already have
OV16_17 <- read.csv("~/Mouse_Eimeria_Databasing/data/Eimeria_detection/FINALqpcrData_2016_2017_threshold5.csv")
levels(factor(OV16_17$qPCRstatus)) # just checking for pos and neg only
# assign Eimeria species
OVsp <- read.csv("~/Mouse_Eimeria_Databasing/data/Eimeria_detection/Eimeria_species_assignment_14_17.csv")
OVsp <- select(OVsp, Mouse_ID, Species)

a <- sapply(OVsp, is.factor)
OVsp[a] <- lapply(OVsp[a], as.character)
b <- sapply(OV16_17, is.factor)
OV16_17[b] <- lapply(OV16_17[b], as.character)
OVsp <- data.frame(OVsp)
OV16_17 <- data.frame(OV16_17)

OV16_17 <- merge(OV16_17, OVsp, by = "Mouse_ID", all.x = T)
OV16_17n <- subset(OV16_17, OV16_17$qPCRstatus == "negative")
OV16_17p <- subset(OV16_17, OV16_17$qPCRstatus == "positive")
# split into years and equalize positive to negative samples 50-50 with random sampling
n16 <- subset(OV16_17n, OV16_17n$year == 2016)
t16 <- subset(OV16_17p, OV16_17p$year == 2016)
DNA16t <- merge(DNA, t16)
DNA16n <- merge(DNA, n16)
n16 <- n16[ !(n16$Mouse_ID %in% DNA16n$Mouse_ID), ]
t16 <- t16[ !(t16$Mouse_ID %in% DNA16t$Mouse_ID), ]
n16 <- sample_n(n16, 20)
OV16 <- rbind(n16, t16)
rm(n16, t16, DNA16n, DNA16t)
OV16 <- select(OV16, Mouse_ID, year, delta_ct_cewe_MminusE, eimeriaSpecies)
OV16$MC <- NA

n17 <- subset(OV16_17n, OV16_17n$Year == 2017)
t17 <- subset(OV16_17p, OV16_17p$Year == 2017)
n17 <- sample_n(n17, 59)
OV17 <- rbind(n17, t17)
rm(n17, t17)
OV17 <- select(OV17, Mouse_ID, Year, delta_ct_cewe_MminusE, eimeriaSpecies)
OV17$MC <- NA

OV18 <- read.csv("~/Mouse_Eimeria_Databasing/data/Eimeria_detection/Svenja/table_ct_and_more.csv")

names(OV18)[names(OV18) == "deltaCtMmE_tissue"] <- "delta_ct_cewe_MminusE"
names(OV18)[names(OV18) == "Eimeria.subspecies"] <- "eimeriaSpecies"
names(OV18)[names(OV18) == "Eimeria.presence.in.Caecum"] <- "MC"
OV18$Year <- 2018
n18 <- subset(OV18, OV18$MC == FALSE)
t18 <- subset(OV18, OV18$MC == TRUE)
n18 <- sample_n(n18, 24)
OV18 <- rbind(t18, n18)
rm(t18, n18)
OV18 <- select(OV18, Mouse_ID, Year, delta_ct_cewe_MminusE, eimeriaSpecies, MC)

rm(OV14_17n, OV14_17p)

OV19 <- read.csv("~/Mouse_Eimeria_Databasing/data/Eimeria_detection/HZ19_CEWE_qPCR.csv")
names(OV19)[names(OV19) == "delta"] <- "delta_ct_cewe_MminusE"
OV19$Year <- 2019
OV19$eimeriaSpecies <- NA
OV19 <- select(OV19, Mouse_ID, Year, delta_ct_cewe_MminusE, eimeriaSpecies, MC)
t19 <- subset(OV19, OV19$MC == TRUE)
n19 <- subset(OV19, OV19$MC == FALSE)
n19 <- sample_n(n19, 68)
OV19 <- rbind(n19, t19)
rm(t19, n19)

# blob everything together
OV <- rbind(OV14, OV15)
OV <- rbind(OV, OV16)
OV <- rbind(OV, OV17)
OV <- rbind(OV, OV18)
OV <- rbind(OV, OV19)
rm(OV14, OV15, OV16, OV17, OV18, OV19)



# standardise eimeria species columns and confitionally fill in missing MCs
levels(factor(OV$eimeriaSpecies))
i <- sapply(OV, is.factor)
OV[i] <- lapply(OV[i], as.character)
OV$eimeriaSpecies[OV$eimeriaSpecies == "Double_ferrisi_vermiformis"] <- "Double"
OV$eimeriaSpecies[OV$eimeriaSpecies == "Double_tbd"] <- "Double"
OV$eimeriaSpecies[OV$eimeriaSpecies == "E_falciformis"] <- "E.falciformis"
OV$eimeriaSpecies[OV$eimeriaSpecies == "E_ferrisi"] <- "E.ferrisi"
OV$eimeriaSpecies[OV$eimeriaSpecies == "E. falciformis"] <- "E.falciformis"
OV$eimeriaSpecies[OV$eimeriaSpecies == "non infected"] <- "Negative"
OV$eimeriaSpecies[OV$eimeriaSpecies == "E. ferrisi"] <- "E.ferrisi"
OV$eimeriaSpecies[OV$eimeriaSpecies == "Eimeria sp."] <- "Other"
levels(factor(OV$eimeriaSpecies)) # success, no NA warnings and only 5 levels of factor

OV$MC[OV$eimeriaSpecies == "E.ferrisi"] <- TRUE
OV$MC[OV$eimeriaSpecies == "E.falciformis"] <- TRUE
OV$MC[OV$eimeriaSpecies == "Double"] <- TRUE
OV$MC[OV$eimeriaSpecies == "Other"] <- TRUE
OV$MC[OV$eimeriaSpecies == "Negative"] <- FALSE
levels(factor(OV$MC)) # success (no NAs)

processed <- merge(DNA, OV, all.x =  T)
missing <- merge(DNA, OV, all.y = T)
