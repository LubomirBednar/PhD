# get libraries
library(RCurl)
library(ggplot2)
library(dplyr)
library(tidyverse)

# get raw files of records
E11aR <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/E11a_012021_Eim_record.csv")
E11aR$X <- NULL
E11aR$batch <- "a"
E11bR <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/E11b_012021_Eim_record.csv")
E11bR$X <- NULL 
E11bR$batch <- "b"
# get raw design file
E11D <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experimental_design/E11_012021_Eim_DESIGN.csv")

# merge together into E10R + remove redundants
E11aR$experiment <- "E11"
E11bR$experiment <- "E11"
E11aR$dpi_dissect <- "NA"
E11aR <- dplyr::select(E11aR, EH_ID, dpi, labels, weight, relative_weight, weight_dpi0, feces_weight, experiment, dpi_dissect, batch)
E11aRD <- merge(E11aR, E11D)
E11aRD$challenge_infection <- NA

E11bRD <- merge(E11bR, E11D)
E11bRD$primary_infection <- NA
rm(E11aR, E11bR)
# merge together into E10D and remove redundants
E11 <- rbind(E11aRD, E11bRD)
rm(E11bRD, E11aRD)
# look at primary infection
ggplot(subset(E11, !is.na(E11$primary_infection)), aes(x = dpi, y = relative_weight, color = primary_infection)) +
  geom_jitter() +
  geom_smooth() +
  geom_hline(yintercept = 82, color = "red")# + facet_wrap(~primary)
# look at challenge infection
ggplot(subset(E11, !is.na(E11$challenge_infection)), aes(x = dpi, y = relative_weight, color = challenge_infection)) +
  geom_jitter() +
  geom_smooth() +
  geom_hline(yintercept = 82, color = "red") #+ facet_wrap(~challenge)
# look at challenge with infection history
ggplot(subset(E11, E11$batch == "b"), aes(x = dpi, y = relative_weight, color = challenge_infection)) +
  geom_jitter() +
  geom_smooth(se = F) +
  facet_wrap(~infection_history) +
  geom_hline(yintercept = 90) +
  geom_hline(yintercept = 82, color = "red")

write.csv(E11, "~/GitHub/Eimeria_Lab/data/Experiment_results/E11_012021_Eim_record.csv")
