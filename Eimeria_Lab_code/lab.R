# P3, P4, E7 and E10 (CLS (classic laboratory strains), WDS(wild derived strains))
library(ggplot2)
library(Rmisc)
library(httr)
library(RCurl)
library(tidyverse)
# load in the COMPLETE tables (should contain everything at LABEL level)
P3 <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/P3_112019_Eim_COMPLETE.csv"))
P4 <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/P4_082020_Eim_COMPLETE.csv"))
E7 <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/E7_112018_Eim_COMPLETE.csv"))
# E10 coming soon
# unify columns by renaming and selecting
### remake with underscores for clarity
colnames(P3)[11] <- "total_oocysts"
colnames(P4)[29] <- "total_oocysts"
colnames(E7)[22] <- "total_oocysts"
colnames(P3)[20] <- "infection_history"
colnames(P4)[5] <- "infection_history"
colnames(E7)[4] <- "infection_history"
colnames(P3)[16] <- "weight_change"
colnames(P4)[9] <- "weight_change"
colnames(E7)[12] <- "weight_change"
colnames(E7)[2] <- "label"
### add strain to P3 and P4 + other necessary columns for later
P3$Strain <- "SWISS"
P4$Strain <- "SWISS"
P4$IFNy_CEWE <- NA
P3$EXP_type <- "CLS"
P4$EXP_type <- "CLS"
E7$EXP_type <- "WDS"
P4$Eim_MC[P4$Eim_MC == "TRUE"] <- "pos"
P4$Eim_MC[P4$Eim_MC == "FALSE"] <- "neg"

### pick down to bare minimum of 13 columns
P3 <- select(P3, EH_ID, dpi, label, delta, IFNy_CEWE, Eim_MC, total_oocysts, infection_history, weight_change, 
             primary, challenge, OPG, IFNy_CEWE, Strain, EXP_type)
P4 <- select(P4, EH_ID, dpi, label, delta, IFNy_CEWE, Eim_MC, total_oocysts, infection_history, weight_change, 
             primary, challenge, OPG, IFNy_CEWE, Strain, EXP_type)
E7 <- select(E7, EH_ID, dpi, label, delta, IFNy_CEWE, Eim_MC, total_oocysts, infection_history, weight_change, 
             primary, challenge, OPG, IFNy_CEWE, Strain, EXP_type)
# add FACS data
# P3_FACS - need this from Hongwei
# E10_FACS - coming soon
P4_FACS <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/P4_082020_Eim_FACS.csv"))
E7_FACS <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/E7_112018_Eim_FACS.csv"))
colnames(E7_FACS)[3] <- "label"
E7_FACS$Position <- "mLN"

# select only FACS columns
P4_FACS <- select(P4_FACS, EH_ID, dpi, label, Position, CD4, Treg, Div_Treg, Treg17, Th1, Div_Th1, Th17, Div_Th17, CD8,
                  Act_CD8, Div_Act_CD8, IFNy_CD4, IFNy_CD8)
E7_FACS <- select(E7_FACS, EH_ID, dpi, label, Position, CD4, Treg, Div_Treg, Treg17, Th1, Div_Th1, Th17, Div_Th17, CD8,
                  Act_CD8, Div_Act_CD8, IFNy_CD4, IFNy_CD8)
# use FACS and completes to create end point immuno datasets
P3_immuno <- subset(P3, !is.na(P3$Eim_MC))
# P3_immuno <- merge(P3_immuno, P3_FACS)
write.csv(P3_immuno, "~/GitHub/Eimeria_Lab/data/Experiment_results/P3_112019_Eim_immuno.csv")
P4_immuno <- subset(P4, !is.na(P4$Eim_MC))
P4_immuno <- merge(P4_immuno, P4_FACS)
write.csv(P4_immuno, "~/GitHub/Eimeria_Lab/data/Experiment_results/P4_082020_Eim_immuno.csv")
E7_immuno <- subset(E7, !is.na(E7$Eim_MC))
E7_immuno <- merge(E7_immuno, E7_FACS)
write.csv(E7_immuno, "~/GitHub/Eimeria_Lab/data/Experiment_results/E7_112018_Eim_immuno.csv")
# P3_immuno <- merge(P3_immuno, P3_FACS)

# purely explorative now

lab <- full_join(P3_immuno, P4_immuno)
lab <- full_join(lab, E7_immuno)
