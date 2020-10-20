# P4a processing into separate clean tables for weightloss, oocysts, qPCR, etc.

library(httr)
library(RCurl)
library(dplyr)
library(ggplot2)
library(Rmisc)

# Design first
P4a_design <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experimental_design/P4a_082020_Eim_design.csv"))

# Record next
P4a_record <- read.csv(text = getURL(""))
P4a_record$X <- NULL
# start combining
P4a <- merge(P4a_design, P4a_record)
