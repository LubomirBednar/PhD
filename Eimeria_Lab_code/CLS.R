library(ggplot2)
library(Rmisc)
library(httr)
library(RCurl)
library(tidyverse)
library(ggpubr)
library(dplyr)

# load in the COMPLETE tables (should contain everything at LABEL level)
P3 <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/P3_112019_Eim_COMPLETE.csv")
P3$X <- NULL
P4 <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/P4_082020_Eim_COMPLETE.csv")
P4$X <- NULL
P4$experiment <- "P4"
P4$labels.1 <- NULL
P4$label <- NULL
# unify columns by renaming and selecting
### remake with underscores for clarity
colnames(P3)[18] <- "infection_history"
colnames(P4)[4] <- "infection_history"

### add strain to P3 and P4 + other necessary columns for later
P3$Strain <- "SWISS"
P4$Strain <- "SWISS"
P3$EXP_type <- "CLS"
P4$EXP_type <- "CLS"
P4$Eim_MC[P4$Eim_MC == "TRUE"] <- "pos"
P4$Eim_MC[P4$Eim_MC == "FALSE"] <- "neg"


P3$AVG <- NULL
P3$faeces_weight <- NULL
P3$Strain <- "SWISS"
P4$EXP_type <- "CLS"
P4$faeces_weight <- NULL
P4$Strain <- "SWISS"

complete <- rbind.fill(P3, P4)

write.csv(complete, "../GitHub/Eimeria_Lab/data/Experiment_results/CLS_complete.csv")

# CLS <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/CLS_complete.csv")

CLS <- complete
CLS$Eimeria1[CLS$primary == "E64"] <- "E.ferrisi"
CLS$Eimeria1[CLS$primary == "E139"] <- "E.ferrisi"
CLS$Eimeria1[CLS$primary == "E88"] <- "E.falciformis"
CLS$Eimeria1[CLS$primary == "Eflab"] <- "E.falciformis"
CLS$Eimeria1[CLS$primary == "UNI"] <- "Uninfected"


CLS$Eimeria2[CLS$challenge == "E64"] <- "E.ferrisi"
CLS$Eimeria2[CLS$challenge == "E139"] <- "E.ferrisi"
CLS$Eimeria2[CLS$challenge == "E88"] <- "E.falciformis"
CLS$Eimeria2[CLS$challenge == "UNI"] <- "Uninfected"
# primary
ggplot(subset(CLS, CLS$OPG > 0 & !is.na(CLS$primary)), aes(Eimeria1, OPG, color = Eimeria1)) +
  geom_boxplot() +
  geom_jitter()

ggplot(subset(CLS, !is.na(CLS$primary)), aes( x = dpi, y = relative_weight, color = Eimeria1)) + 
  geom_smooth() + 
  geom_jitter() + 
  facet_wrap(~Eimeria1) +
  ylim(82, 110)

ggplot(subset(CLS, !is.na(CLS$primary)), aes( x = dpi, y = relative_weight, color = primary)) + 
  geom_smooth() + 
  geom_jitter() + 
  facet_wrap(~primary) + 
  ylim(70, 110)

# challenge
ggplot(subset(CLS, CLS$OPG > 0 & !is.na(CLS$challenge)), aes(Eimeria2, OPG, color = Eimeria2)) +
  geom_boxplot() +
  geom_jitter()

ggplot(subset(CLS, !is.na(CLS$challenge)), aes( x = dpi, y = relative_weight, color = Eimeria2)) + 
  geom_smooth() + 
  geom_jitter() + 
  facet_wrap(~Eimeria2)

ggplot(subset(CLS, !is.na(CLS$challenge)), aes( x = dpi, y = relative_weight, color = challenge)) + 
  geom_smooth() + 
  geom_jitter() + 
  facet_wrap(~challenge) + 
  ylim(82, 110)

# batches
ggplot(subset(CLS, CLS$OPG > 0), aes(x = dpi, y = OPG, color = batch)) +
  geom_jitter() +
  facet_wrap(~infection_history)

ggplot(subset(CLS, !is.na(CLS$challenge)), aes(x = dpi, y = relative_weight, color = challenge)) +
  geom_jitter() +
  geom_smooth() +
  facet_wrap(~infection_history)

ggplot(subset(CLS, !is.na(CLS$primary)), aes(x = dpi, y = relative_weight, color = primary)) +
  geom_jitter() +
  geom_smooth() +
  facet_wrap(~primary)

ggplot(subset(CLS, !is.na(CLS$primary)), aes(x = OPG, y = relative_weight, color = primary)) +
  geom_jitter() +
  geom_smooth() +
  facet_wrap(~primary)

weight <- select(CLS, relative_weight, dpi)

ggplot(subset(CLS, !is.na(CLS$primary)), aes(dpi, log(OPG))) + 
  geom_line() +
  geom_line(data = CLS, color = "red")
