set.seed(123)
# get the samples we already have
DNA <- read.csv("https://raw.githubusercontent.com/LubomirBednar/PhD/master/cDNA_done.csv")
##################
################## make a big table of all the samples and remove the ones we already have
##################
OV16_17 <- read.csv("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Field/master/data/Eimeria_detection/FINALqpcrData_2016_2017_threshold5.csv")
levels(factor(OV16_17$qPCRstatus)) # just checking for pos and neg only
# assign Eimeria species
OVsp <- read.csv("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Field/master/data/Eimeria_detection/Eimeria_species_assignment_14_17.csv")
OVsp <- dplyr::select(OVsp, Mouse_ID, Species)
# remove invisible spaces
OVsp$Mouse_ID <- gsub(" ", "", x = OVsp$Mouse_ID)
# create overall 2016 - 17 table
OV16_17 <- merge(OV16_17, OVsp, by = "Mouse_ID")
names(OV16_17)[names(OV16_17) == "qPCRstatus"] <- "MC"
OV16_17n <- subset(OV16_17, OV16_17$MC == "negative")
OV16_17n$MC[OV16_17n$MC == "negative"] <- FALSE
OV16_17p <- subset(OV16_17, OV16_17$MC == "positive")
OV16_17p$MC[OV16_17p$MC == "positive"] <- TRUE
########## split into years and equalize positive to negative samples 50-50 with random sampling
# create negative group for 2016
n16 <- subset(OV16_17n, OV16_17n$year == 2016)
# create positive group for 2016
t16 <- subset(OV16_17p, OV16_17p$year == 2016)
# check what we already have from positive 2016
DNA16t <- merge(DNA, t16)
# check what we already have from negative 2016
DNA16n <- merge(DNA, n16)
# remove already processed samples 
n16 <- n16[ !(n16$Mouse_ID %in% DNA16n$Mouse_ID), ]
t16 <- t16[ !(t16$Mouse_ID %in% DNA16t$Mouse_ID), ]
# because t16 is 20 samples, randomly sample 20 n16 to have 50-50 spread
n16 <- dplyr::sample_n(n16, 20)
OV16 <- rbind(n16, t16)
rm(n16, t16, DNA16n, DNA16t)
OV16 <- dplyr::select(OV16, Mouse_ID, year, delta_ct_cewe_MminusE, Species, MC)

########### repeat for 2017
n17 <- subset(OV16_17n, OV16_17n$year == 2017)
t17 <- subset(OV16_17p, OV16_17p$year == 2017)
DNA17t <- merge(DNA, t17)
DNA17n <- merge(DNA, n17)

n17 <- n17[ !(n17$Mouse_ID %in% DNA17n$Mouse_ID), ]
t17 <- t17[ !(t17$Mouse_ID %in% DNA17t$Mouse_ID), ]
# t17 is 47 samples so "sample_n" 47 n17
n17 <- dplyr::sample_n(n17, 47)
OV17 <- rbind(n17, t17)
rm(n17, t17, DNA17n, DNA17t)
OV17 <- dplyr::select(OV17, Mouse_ID, year, delta_ct_cewe_MminusE, Species, MC)

# same for 2018
OV18 <- read.csv("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Field/master/data/Eimeria_detection/Svenja/table_ct_and_more.csv")

names(OV18)[names(OV18) == "deltaCtMmE_tissue"] <- "delta_ct_cewe_MminusE"
names(OV18)[names(OV18) == "Eimeria.subspecies"] <- "Species"
names(OV18)[names(OV18) == "Eimeria.presence.in.Caecum"] <- "MC"
names(OV18)[names(OV18) == "Name"] <- "Mouse_ID"
OV18$year <- 2018
OV18 <- dplyr::select(OV18, Mouse_ID, year, delta_ct_cewe_MminusE, Species, MC)

n18 <- subset(OV18, OV18$MC == FALSE)
t18 <- subset(OV18, OV18$MC == TRUE)
DNA18t <- merge(DNA, t18)
DNA18n <- merge(DNA, n18)

n18 <- n18[ !(n18$Mouse_ID %in% DNA18n$Mouse_ID), ]
t18 <- t18[ !(t18$Mouse_ID %in% DNA18t$Mouse_ID), ]
# looks like all samples are processed so no need for sampling more negatives
rm(OV18, n18, t18, DNA18n, DNA18t)
# removed everything because all 2018 samples have already been extracted

OV19 <- read.csv("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Field/master/data/Eimeria_detection/HZ19_CEWE_qPCR.csv")
names(OV19)[names(OV19) == "delta"] <- "delta_ct_cewe_MminusE"
OV19$year <- 2019
OV19$Species <- NA # because we don't know those yet
OV19 <- dplyr::select(OV19, Mouse_ID, year, delta_ct_cewe_MminusE, Species, MC)
t19 <- subset(OV19, OV19$MC == TRUE)
n19 <- subset(OV19, OV19$MC == FALSE)

DNA19t <- merge(DNA, t19)
DNA19n <- merge(DNA, n19)

n19 <- n19[ !(n19$Mouse_ID %in% DNA19n$Mouse_ID), ]
t19 <- t19[ !(t19$Mouse_ID %in% DNA19t$Mouse_ID), ]

# remaining t19 after subtracting samples we already have is 31 so match that in negatives
n19 <- dplyr::sample_n(n19, 31)
OV19 <- rbind(n19, t19)
rm(t19, n19, DNA19n, DNA19t)
# match colnames


# blob everything together
OV <- rbind(OV16, OV17)
OV <- subset(OV, OV$Species != "E_vermiformis") # remove vermiformis
OV <- rbind(OV, OV19)
rm(OV14, OV15, OV16, OV17, OV18, OV19, OV16_17, OV16_17n, OV16_17p)
# check the integrity
levels(factor(OV$MC))
levels(factor(OV$Species))

# write out
write.csv(OV, "~/GitHub/PhD/missing_IFC_samples_to_process.csv")


