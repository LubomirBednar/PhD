# P3 combining script for weight, oocysts, qPCR, RT-qPCR, ELISA and hopefully FACS
library(httr)
library(RCurl)
library(dplyr)
library(Rmisc)

# load in weight and oocysts

P3_oocyst1 <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/P3_112019_Eim_oocyst.csv")
P3_oocyst1$X <- NULL


P3_oocyst2 <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/P3_112019_Eim_oocysts.csv")
P3_oocyst2$X <- NULL
names(P3_oocyst2)[names(P3_oocyst2) == "oocyst_1"] <- "oocyst_sq1"
names(P3_oocyst2)[names(P3_oocyst2) == "oocyst_2"] <- "oocyst_sq2"
names(P3_oocyst2)[names(P3_oocyst2) == "oocyst_3"] <- "oocyst_sq3"
names(P3_oocyst2)[names(P3_oocyst2) == "oocyst_4"] <- "oocyst_sq4"
names(P3_oocyst2)[names(P3_oocyst2) == "AVG"] <- "oocyst_mean"

P3_oocyst <- merge(P3_oocyst1, P3_oocyst2, all = T)


P3a_record <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/P3a_112019_Eim_Record.csv")
P3b_record <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/P3b_112019_Eim_Record.csv")
P3b_record$X <- NULL
P3a_record$labels <- sub("^", "P3a", P3a_record$labels)
P3a_record$batch <- "a"
P3b_record$labels <- sub("^", "P3b", P3b_record$labels)
P3b_record$batch <- "b"
P3_record <- rbind(P3a_record, P3b_record)
P3_para <- merge(P3_record, P3_oocyst)
write.csv(P3_para, "C:/Users/exemp/Documents/P3_para.csv")
P3_para <- read.csv("C:/Users/exemp/Documents/P3_para.csv")

P3_design <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experimental_design/P3_112019_Eim_design.csv")

P3_para <- merge(P3_para, P3_design, all.x = T) 
P3_para$day_change <- NULL
# load in qPCRs
P3_qPCR <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/P3_112019_Eim_CEWE_qPCR.csv")
P3_qPCR$X <- NULL
P3_qPCR$dpi <- 8
P3_qPCR$batch <- "b"
P3 <- merge(P3_para, P3_qPCR, all.x = T)

# load in RT-qPCRs
# P3_RT <- "https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/P3_112019_Eim_CEWE_RTqPCR.csv"
# P3_RT <- read.csv(text = getURL(P3_RT))
# P3_RT$X <- NULL

# load in CEWE ELISA (important to merge CEWE ELISAs with qPCR and RTqPCR to give them labels)
P3_CEWE_ELISA <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/P3_112019_Eim_CEWE_ELISA.csv")
P3_CEWE_ELISA$X <- NULL
colnames(P3_CEWE_ELISA)[2] <- "IFNy_CEWE"
P3 <- merge(P3, P3_CEWE_ELISA, all.x = T)


# # load in FEC ELISA
# P3_FEC_ELISA <- "https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/P3_112019_Eim_FEC_ELISAs/P3_112019_Eim_FEC_ELISA1_complete.csv"
# P3_FEC_ELISA <- read.csv(text = getURL(P3_FEC_ELISA))
# P3_FEC_ELISA$X <- NULL
# colnames(P3_FEC_ELISA)[2] <- "IFNy_FEC"

# load in qPCR MCs for Eimeria
P3_MC <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/E7%26P3_Eim_MCs.csv")
P3_MC$X <- NULL
P3_MC$X.1 <- NULL
P3_MC$batch <- "b"
P3_MC$dpi <- 8

P3 <- merge(P3, P3_MC, all.x = T)

# up to date most complete P3 dataset
write.csv(P3, "./Eimeria_Lab/data/Experiment_results/P3_112019_Eim_COMPLETE.csv")
write.csv(P3, "../GitHub/Eimeria_Lab/data/Experiment_results/P3_112019_Eim_COMPLETE.csv")
