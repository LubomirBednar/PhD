# P4b processing into separate clean tables for weightloss, oocysts, qPCR, etc.

library(httr)
library(RCurl)
library(dplyr)
library(ggplot2)
library(Rmisc)

# Design first
P4b_design <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experimental_design/P4b_082020_Eim_design.csv"))
# Record next
P4b_record <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/P4b_082020_Eim_Record.csv"))
# add oocysts
P4b_oocyst <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/P4b_082020_Eim_oocysts.csv"), sep = "")
# add qPCR results
P4b_qPCR <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/P4_082020_Eim_CEWE_qPCR.csv"))
P4b_qPCR$dpi <- 8
# start combining
P4b <- merge(P4b_design, P4b_record, all.y = T)
P4b <- merge(P4b, P4b_oocyst, all = T)
P4b <- merge(P4b, P4b_qPCR, all = T)
# add P4b to lables

# calculate OPG and total oocysts

# remove useless columns

# add batch

# make sure columns are: label, EH_ID, IFNy_CEWE,	delta, Eim_MC,	AVG,	totalOocysts,	dpi,	weight,	weight_dpi0,	
# faeces_weight,	wloss, batch, primary,	challenge,	infHistory,	OPG.
