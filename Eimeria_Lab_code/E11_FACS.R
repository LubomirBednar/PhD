# E11 FACS
library(tidyverse)
library(readxl)
library(dplyr)
library(ggplot2)
library(ggpubr)



# load in xlsx by sheet + remove NA rows 
FACSraw1 <- read_xlsx("~/Downloads/Lubomir FACS data_03.03.2021.xlsx", sheet = 1)
FACSraw1 <- na.omit(FACSraw1)
FACSraw2 <- read_xlsx("~/Downloads/Lubomir FACS data_03.03.2021.xlsx", sheet = 2)
FACSraw2 <- na.omit(FACSraw2)

# separate Sample into EH_ID and organ (aka. position) + combine
FACSraw1 <- separate(FACSraw1, Sample, into = c("EH_ID", "Position"))
FACSraw2 <- separate(FACSraw2, Sample, into = c("EH_ID", "Position"))
FACS <- full_join(FACSraw1, FACSraw2)
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
FACS$dpi <- 8
FACS$batch <- "b"
FACS.long <- pivot_longer(FACS, 
             !c("EH_ID", "dpi", "batch", "Position"), names_to = "Pop", values_to = "count")

label <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/E11_012021_Eim_record.csv")
FACS <- merge(FACS, label)
FACS.long <- merge(FACS.long, label)


# write out 
write.csv(FACS, "~/Eimeria_Lab/data/Experiment_results/E11_012021_Eim_FACS.csv")

ggplot(FACS.long, aes(x = challenge_infection, y = count, color = challenge_infection)) +
  geom_violin() + 
  stat_comp
  facet_wrap(~Pop, scales = "free_y")
