# P4 FACS
library(httr)
library(RCurl)
library(Rmisc)
library(tidyverse)
library(readxl)
library(dplyr)
# load in xlsx by sheet + remove NA rows 
FACSraw1 <- read_xlsx("~/GitHub/Eimeria_Lab/data/Experiment_results/Lubomir FACS data_02.10.2020.xlsx", sheet = 1)
FACSraw1 <- FACSraw1[-c(14,15),]
FACSraw2 <- read_xlsx("~/GitHub/Eimeria_Lab/data/Experiment_results/Lubomir FACS data_02.10.2020.xlsx", sheet = 2)
FACSraw2 <- FACSraw2[-c(14,15),]
# separate Sample into EH_ID and organ (aka. position) + combine
FACSraw1 <- separate(FACSraw1, Sample, into = c("EH_ID", "Position"))
FACSraw2 <- separate(FACSraw2, Sample, into = c("EH_ID", "Position"))
FACS <- full_join(FACSraw1, FACSraw2)
FACS$EH_ID <- sub("LM", "LM_", FACS$EH_ID)
# create R and .csv friendly column names
colnames(FACS)[3]<- "CD4"
colnames(FACS)[4]<- "Treg"
colnames(FACS)[5]<- "Div_Treg"
colnames(FACS)[6]<- "Treg17"
FACS$`FSC-A, SSC-A subset/single/live/CD4+/Foxp3-,Freq. of Parent` <- NULL
colnames(FACS)[7]<- "Th1"
colnames(FACS)[8]<- "Div_Th1"
colnames(FACS)[9]<- "Th17"
colnames(FACS)[10]<- "Div_Th17"
colnames(FACS)[11]<- "CD8"
colnames(FACS)[12]<- "Act_CD8"
colnames(FACS)[13]<- "Div_Act_CD8"
colnames(FACS)[14]<- "IFNy_CD4"
colnames(FACS)[15]<- "IFNy_CD8"
# remove "%" signs and convert to numeric
FACS[] <- lapply(FACS, gsub, pattern='%', replacement='')
FACS[, 3:15] <- sapply(FACS[, 3:15], as.numeric)
FACS$Position <- as.factor(FACS$Position)
# add labels and dpi reference columns
design <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experimental_design/P4b_082020_Eim_design.csv"))
design <- subset(design, design$dpi == 8)
FACS <- merge(FACS, design)
colnames(FACS)[17] <- "label"
FACS$label <- sub("^", "P4b", FACS$label)


# write out 
write.csv(FACS, "~/GitHub/Eimeria_Lab/data/Experiment_results/P4_082020_Eim_FACS.csv")
