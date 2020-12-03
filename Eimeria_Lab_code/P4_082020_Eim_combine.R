# P4 processing script into clean tables

library(tidyverse)
library(ggplot2)
library(Rmisc)
library(RCurl)
library(httr)

# add design P4a (add batch, amend labels)
P4a_design <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experimental_design/P4a_082020_Eim_design.csv"))
# add and clean oocyst P4a
P4a_oocyst <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/P4a_082020_Eim_oocysts1.csv"))
# add and clean record P4a
P4a_record <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/P4a_082020_Eim_Record.csv"))
P4a_record$relative_weight <- as.numeric(as.character(P4a_record$relative_weight))
P4a_record$relative_weight[which(P4a_record$relative_weight == "#VALUE!")] <- NA
# add design P4b
P4b_design <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experimental_design/P4b_082020_Eim_design.csv"))
# add and clean oocyst P4b
P4b_oocyst <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/P4b_082020_Eim_oocysts.csv"))
# add and clean record P4b
P4b_record <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/P4b_082020_Eim_Record.csv"))
P4b_record$relative_weight <- as.numeric(as.character(P4b_record$relative_weight))
# add cleaned qPCR P4b
P4b_qPCR <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/P4_082020_Eim_CEWE_qPCR.csv"))
P4b_qPCR$dpi <- 8
# add P4b FACS
P4b_FACS <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/P4_082020_Eim_FACS.csv"))
P4b_FACS$dpi <- 8
P4b_FACS$X <- NULL
# add P4b ELISA
P4b_ELISA <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/P4_082020_Eim_ELISA.csv"))
P4b_ELISA$X <- NULL
P4b_ELISA$dpi <- 8

# make infection history
P4ainf <- dplyr::select(P4a_design, EH_ID, primary)
P4ainf <- distinct(P4ainf)
P4binf <- dplyr::select(P4b_design, EH_ID, challenge)
P4binf <- distinct(P4binf)
P4_inf <- merge(P4ainf, P4binf)
P4_inf$infHistory <- paste(P4_inf$primary, P4_inf$challenge, sep =  ":")
# add inf history to P4a and P4b designs
P4_inf$primary <- NULL
P4_inf$challenge <- NULL
P4a_design <- merge(P4a_design, P4_inf, all = T)
P4b_design <- merge(P4b_design, P4_inf, all = T)
# create P4a
P4a <- merge(P4a_design, P4a_oocyst, all = T)
P4a <- merge(P4a, P4a_record, all = T)
# add batch and label uniques
P4a$batch <- "a"
P4a$labels <- sub("^", "P4a", P4a$labels)

# merge all into P4b and then P4
P4b <- merge(P4b_design, P4b_oocyst, all = T)
P4b <- merge(P4b, P4b_record, all = T)
P4b <- merge(P4b, P4b_qPCR, all = T)
P4b$batch <- "b"
P4b$labels <- sub("^", "P4b", P4b$labels)
P4b <- merge(P4b, P4b_FACS, all = T)
P4b <- merge(P4b, P4b_ELISA, all = T)
# merge P4a and P4b
P4 <- merge(P4a, P4b, all = T)

# calculate oocyst AVG, total and OPG
P4$total_oocysts <- ((P4$oocyst_sq1 
                          + P4$oocyst_sq2 
                          + P4$oocyst_sq3 
                          + P4$oocyst_sq4) / 4) * 
  10000 * # because volume chamber
  2
P4$OPG <- P4$total_oocysts/ P4$feces_weight

# replace !VALUE# with NAs
P4$day_change[P4$day_change == "#VALUE!"] <- NA
# final clean to have same columns as P3:
# label	EH_ID	IFNy_CEWE	delta	CXCR3	IL.12	IRG6	Eim_MC	AVG	totalOocysts	dpi	weight	
# weight_dpi0	faeces_weight	wloss	batch	primary	challenge	infHistory	OPG	IFNy_FEC

P4$day_change <- NULL
P4$experiment <- "P4"
P4$dilution <- 1


write.csv(P4, "~/GitHub/Eimeria_Lab/data/Experiment_results/P4_082020_Eim_COMPLETE.csv")
P4o <- select(P4, labels, experiment, oocyst_sq1, oocyst_sq2, oocyst_sq3, oocyst_sq4, oocyst_mean, OPG, dilution)
write.csv(P4o, "~/GitHub/Eimeria_Lab/data/Experiment_results/P4_082020_Eim_oocyst.csv")
P4w <- select(P4, labels, experiment, weight, weight_dpi0, relative_weight, feces_weight, dpi)
write.csv(P4w, "~/GitHub/Eimeria_Lab/data/Experiment_results/P4_082020_Eim_record.csv")

ggplot(subset(P4, !is.na(P4$Eim_MC)), aes(x = wloss, y = delta, color = Eim_MC)) +
  geom_jitter() +
  facet_wrap(~infHistory)

ggplot(subset(P4, !is.na(P4$Eim_MC)), aes(x = dpi, y = OPG, color = challenge)) + 
  geom_jitter() +
  geom_smooth()

ggplot(subset(P4, !is.na(P4$primary)), aes(x = dpi, y = OPG, color = primary))  +
  geom_jitter(width = 0.2, size = 2) + 
  geom_smooth(se = F) +
  facet_wrap(~primary)

ggplot(P4, aes(x = dpi, y = OPG, color = batch))  +
  geom_jitter(width = 0.2, size = 2) + 
  geom_smooth(se = F) +
  facet_wrap(~infHistory)

