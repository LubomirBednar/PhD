# P3 RT-qPCR and infection intensity qPCR analysis
library(tidyverse)
library(RCurl)
library(reshape2)
library(ggpubr)
library(ggplot2)
library(naniar)
library(data.table)

# load in raw

RT1 <- "https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/P3_112019_Eim_RT-qPCRs/P3_112019_Eim_RT-qPCR1.CSV"
RT1 <- read.csv(text = getURL(RT1))

RT2 <- "https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/P3_112019_Eim_RT-qPCRs/P3_112019_Eim_RT-qPCR3.CSV"
RT2 <- read.csv(text = getURL(RT2))

RT3 <- "https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/P3_112019_Eim_RT-qPCRs/P3_112019_Eim_RT-qPCR5.CSV"
RT3 <- read.csv(text = getURL(RT3))

RT4 <- "https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/P3_112019_Eim_RT-qPCRs/P3_112019_Eim_RT-qPCR7.CSV"
RT4 <- read.csv(text = getURL(RT4))

RT <- rbind(RT1, RT2)
RT <- rbind(RT, RT3)
RT <- rbind(RT, RT4)
RT$Pos <- NULL
RT$Amount.SYBR <- NULL
RT$Ct.Dev..SYBR <- NULL

names(RT)[names(RT) == "Target.SYBR"] <- "Target"
names(RT)[names(RT) == "Ct.SYBR"] <- "RT.Ct"
names(RT)[names(RT) == "Name"] <- "EH_ID"
RT <- RT %>% dplyr::group_by(EH_ID, Target) %>% dplyr::summarise(RT.Ct = mean(RT.Ct, na.rm = T))
RT <- ungroup(RT)
RT.wide <- reshape2::dcast(RT, EH_ID ~ Target, value.var = "RT.Ct")

refGenes <- c("B-actin", "GAPDH")
targetGenes <- c("CXCR3", "IL-12", "IRG6")

RT.wide <- RT.wide %>% mutate(refMean = rowMeans(na.rm = TRUE, dplyr::select(RT.wide, refGenes)))
RT.wide <- data.frame(RT.wide)
# Ref - Target = lower value (bigger minus means less expression)
RT.wide$CXCR3 <- (RT.wide$refMean - RT.wide$CXCR3)
RT.wide$IRG6 <- (RT.wide$refMean - RT.wide$IRG6)
RT.wide$IL.12 <- (RT.wide$refMean - RT.wide$IL.12)
RT.wide$refMean <- NULL
RT.wide$B.actin <- NULL
RT.wide$GAPDH <- NULL

#### write out
write.csv(RT.wide, "./Eimeria_Lab/data/EXperiment_results/P3_112019_Eim_CEWE_RTqPCR.csv")
write.csv(RT.wide, "D:/Eimeria_Lab/data/Experiment_results/P3_112019_Eim_CEWE_RTqPCR.csv")
