# E10 and E11 qPCRs
library(dplyr)
library(reshape2)
library(tidyr)

# load in files
a <- read.csv("C:/Users/exemp/Downloads/0101.csv")
b <- read.csv("C:/Users/exemp/Downloads/1.csv")
c <- read.csv("C:/Users/exemp/Downloads/2.csv")
c$Sample <- paste0("LM0", c$Sample)
d <- read.csv("C:/Users/exemp/Downloads/5.csv")
e <- read.csv("C:/Users/exemp/Downloads/6.csv")
f <- read.csv("C:/Users/exemp/Downloads/Luk2.csv")
g <- read.csv("C:/Users/exemp/Downloads/Luke1.csv")


# combine
DF <- bind_rows(a, b, c, d, e, f, g)
DF <- subset(DF, DF$Omit == "false")

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
# add melting curves
MC <- read.csv("C:/Users/exemp/OneDrive/Desktop/E10_E11_qPCRs/E10_E11_CEWE_qPCR_MC.csv")

DF1 <- merge(MC, DF1)
DF1$Eim_MC <- gsub("_", "", DF1$Eim_MC)


library(ggplot2)

ggplot(DF1, aes(x = infection_history, y = delta, color = Eim_MC)) +
  geom_boxplot()

library(data.table)
setDT(DF1)[infection_history=="UNI:UNI", delta:=paste0(delta,NA)]


DF1$delta <- ifelse(DF1$infection_history == "UNI:UNI", NA,no = T)
DF2 <- subset(DF1, DF1$experiment == "E11")
DF3 <- subset(DF1, DF1$experiment == "E10")

write.csv(DF3, "C:/Users/exemp/OneDrive/Documents/GitHub/Eimeria_Lab/data/Experiment_results/E10_112020_Eim_CEWE_qPCR.csv")
write.csv(DF2, "C:/Users/exemp/OneDrive/Documents/GitHub/Eimeria_Lab/data/Experiment_results/E11_012021_Eim_CEWE_qPCR.csv")
