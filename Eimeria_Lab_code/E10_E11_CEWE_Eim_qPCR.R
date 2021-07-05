# E10 and E11 qPCRs
library(dplyr)
library(reshape2)
library(tidyr)

# load in files
a <- read.csv("C:/Users/exemp/OneDrive/Desktop/E10_E11_qPCRs/E10_E11_CEWE_qPCR_01_repeat.csv")
b <- read.csv("C:/Users/exemp/OneDrive/Desktop/E10_E11_qPCRs/E10_E11_CEWE_qPCR_02_Franzi.csv")
b$Sample <- paste0("LM0", b$Sample)
c <- read.csv("C:/Users/exemp/OneDrive/Desktop/E10_E11_qPCRs/E10_E11_CEWE_qPCR.csv")
d <- read.csv("C:/Users/exemp/OneDrive/Desktop/E10_E11_qPCRs/E10_E11_CEWE_qPCR_01_01.csv")
e <- read.csv("C:/Users/exemp/OneDrive/Desktop/E10_E11_qPCRs/E10_E11_qPCR_CEWE_05.csv")
f <- read.csv("C:/Users/exemp/OneDrive/Desktop/E10_E11_qPCRs/E10_E11_qPCR_CEWE_06.csv")
# add MCs
g <- read.csv("C:/Users/exemp/OneDrive/Desktop/E10_E11_qPCRs/")

# combine
DF <- bind_rows(a, b, c, d, e, f)
DF <- select(DF, Sample, Cq, Target)
DF$Target[DF$Target == "Mus"] <- "mouse"
DF$Target[DF$Target == "CDC42"] <- "mouse"
DF$Target[DF$Target == "Eim"] <- "eimeria"
DF$Target[DF$Target == "COI"] <- "eimeria"
DF <- subset(DF, !is.na(DF$Cq))

DF <- DF %>% group_by(Sample, Target) %>% summarise(across(everything(), list(mean)))
DF <- pivot_wider(data = DF, names_from = Target, values_from = Cq_1)
DF$delta <- DF$mouse - DF$eimeria

design1 <- read.csv("~/GitHub/Eimeria_Lab/data/Experimental_design/E10_112020_Eim_DESIGN.csv")
design1$birthday <- NA
design1$EH_ID <- gsub("_", "", design1$EH_ID)
design2 <- read.csv("~/GitHub/Eimeria_Lab/data/Experimental_design/E11_012021_Eim_DESIGN.csv")

design <- bind_rows(design1, design2)
names(DF)[names(DF) == "Sample"] <- "EH_ID"

DF1 <- merge(design, DF, all.x = T)

missing <- subset(DF1, is.na(DF1$delta))
