library(ggplot2)
library(Rmisc)
library(httr)
library(RCurl)
library(tidyverse)
library(ggpubr)

# load in the COMPLETE tables (should contain everything at LABEL level)
P3 <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/P3_112019_Eim_COMPLETE.csv"))
P3$X <- NULL
P4 <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/P4_082020_Eim_COMPLETE.csv"))
P4$X <- NULL
P4$experiment <- "P4"
P4$label.1 <- NULL
# unify columns by renaming and selecting
### remake with underscores for clarity
colnames(P3)[10] <- "total_oocysts"
colnames(P4)[29] <- "total_oocysts"

colnames(P3)[19] <- "infection_history"
colnames(P4)[4] <- "infection_history"

colnames(P3)[15] <- "weight_change"
colnames(P4)[8] <- "weight_change"
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
