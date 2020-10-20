# P4 processing script into clean tables

library(tidyverse)
library(ggplot2)
library(Rmisc)
library(RCurl)
library(httr)

# add design P4a (add batch, amend labels)
P4a_design <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experimental_design/P4a_082020_Eim_design.csv"))
# add and clean oocyst P4a
P4a_oocyst <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/P4a_082020_Eim_oocysts.csv"))
# add and clean record P4a
P4a_record <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/P4a_082020_Eim_Record.csv"))
# add design P4b
P4b_design <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experimental_design/P4b_082020_Eim_design.csv"))
# add and clean oocyst P4b
P4b_oocyst <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/P4b_082020_Eim_oocysts.csv"))
# add and clean record P4b
P4b_record <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/P4b_082020_Eim_Record.csv"))
# add cleaned qPCR P4b
P4b_qPCR <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/P4_082020_Eim_CEWE_qPCR.csv"))
P4b_qPCR$dpi <- 8

# make infection history
P4ainf <- select(P4a_design, EH_ID, primary)
P4ainf <- distinct(P4ainf)
P4binf <- select(P4b_design, EH_ID, challenge)
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
P4b <- merge(P4b_design, P4b_oocyst, all = T) # MLZ twice but both 0
P4b <- merge(P4b, P4b_record, all = T)
P4b <- merge(P4b, P4b_qPCR, all = T)
P4b$batch <- "b"
P4b$labels <- sub("^", "P4b", P4b$labels)
# merge P4a and P4b
P4 <- merge(P4a, P4b, all = T)

# calculate oocyst AVG, total and OPG
P4$totalOocysts <- ((P4$oocyst_1 
                          + P4$oocyst_2 
                          + P4$oocyst_3 
                          + P4$oocyst_4) / 4) * 
  10000 * # because volume chamber
  2
P4$OPG <- P4$totalOocysts/ P4$faeces_weight


