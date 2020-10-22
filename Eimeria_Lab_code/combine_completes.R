library(tidyverse)
library(ggplot2)
library(Rmisc)
library(httr)
library(RCurl)

E7q <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/E7_112018_Eim_CEWE_qPCR.csv"))
E7q$X <- NULL
E7e <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/E7_112018_Eim_CEWE_ELISA.csv"))
E7e$X <- NULL
E7f <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/E7_112018_Eim_FACS.csv"))
E7f$X <- NULL
E7p <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/E7_112018_Eim_weight_oocyst.csv"))
E7p$X <- NULL

E7 <- merge(E7q, E7e)
E7 <- merge(E7, E7f)
E7 <- merge(E7, E7p)
E7$EXP_type <- "WDS"


P3q <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/P3_112019_Eim_CEWE_qPCR.csv"))
P3q$X <- NULL
MC <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/E7%26P3_Eim_MCs.csv"), sep = ";")
MC$X <- NULL
MC$X.1 <- NULL
P3q <- merge(P3q, MC)
P3e <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/P3_112019_Eim_CEWE_ELISA.csv"))
P3e$X <- NULL
P3p <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/P3_112019_Eim_Weight%26Oocyst_complete.csv"))
P3p$X <- NULL

P3 <- merge(P3q, P3e)
colnames(P3p)[1] <- "label"
P3 <- merge(P3, P3p)
P3$EXP_type <- "CLS"
colnames(E7)[2] <- "label"

# add P4
P4 <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/P4_082020_Eim_COMPLETE.csv"))
P4$X <- NULL
# order columns for easier navigation
P3 <- P3[,order(colnames(P3))]
E7 <- E7[,order(colnames(E7))]
P4 <- P4[,order(colnames(P4))]

P3$AVG <- NULL
P3$batch <- NULL
E7$average <- NULL
colnames(E7)[2] <- "Eim_MC"
E7$comment <- NULL
E7$CXCR3 <- NULL
E7$IRG6 <- NULL
E7$IFNy_FEC <- NULL
E7$IL.12 <- NULL
E7$label.1 <- NULL
E7$oocyst_1 <- NULL
E7$oocyst_2 <- NULL
E7$oocyst_3 <- NULL
E7$oocyst_4 <- NULL
E7$volume_PBS_mL <- NULL
P3$faeces_weight <- NULL
E7$fecweight <- NULL
P3$Strain <- "SWISS"
colnames(P3)[7] <- "IFNy_CEWE"
colnames(P3)[15] <- "Wchange"
E7$IFNy <- NULL
P4$batch <- NULL
P4$EXP_type <- "CLS"
P4$IFNy_CEWE <- NA
P4$faeces_weight <- NULL
P4$Strain <- "SWISS"
colnames(P4)[27] <- "Wchange"
P4 <- subset(P4, !is.na(P4$Eim_MC))
P4$Eim_MC[P4$Eim_MC == "TRUE"] <- "pos"
P4$Eim_MC[P4$Eim_MC == "FALSE"] <- "neg"

complete <- rbind.fill(P3, E7)
complete <- rbind.fill(complete, P4)


write.csv(complete, "../GitHub/Eimeria_Lab/data/Experiment_results/E7_P3_P4_Eim_complete.csv")

write.csv(complete, "../Eimeria_Lab/data/Experiment_results/E7_P3_P$_Eim_complete.csv")
