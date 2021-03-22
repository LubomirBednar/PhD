# get libraries
library(RCurl)
library(ggplot2)
library(dplyr)
library(tidyverse)

# get raw files of records
E10aR <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/E10a_112020_Eim_record.csv")
E10aR$X <- NULL
E10aR$batch <- "a"
E10bR <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/E10b_112020_Eim_record.csv")
E10bR$X <- NULL 
E10bR$batch <- "b"
# get raw design file
E10D <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experimental_design/E10_112020_Eim_DESIGN.csv")

# merge together into E10R + remove redundants
E10aR <- dplyr::select(E10aR, EH_ID, dpi, labels, weight, batch, relative_weight, weight_dpi0, feces_weight, mouse_strain, experiment)
E10aRD <- merge(E10aR, E10D)

E10bRD <- merge(E10bR, E10D)
rm(E10aR, E10bR)
# merge together into E10D and remove redundants
E10 <- rbind(E10aRD, E10bRD)
rm(E10bRD, E10aRD)
# look at primary infection
ggplot(subset(E10, !is.na(E10$primary_infection)), aes(x = dpi, y = relative_weight, color = primary_infection)) +
  geom_jitter() +
  geom_smooth() +
  geom_hline(yintercept = 82, color = "red")# + facet_wrap(~primary)
# look at challenge infection
ggplot(subset(E10, !is.na(E10$challenge_infection)), aes(x = dpi, y = relative_weight, color = challenge_infection)) +
  geom_jitter() +
  geom_smooth() +
  geom_hline(yintercept = 82, color = "red") #+ facet_wrap(~challenge)
# look at challenge with infection history
ggplot(subset(E10, E10$batch == "b"), aes(x = dpi, y = relative_weight, color = challenge_infection)) +
  geom_jitter() +
  geom_smooth(se = F) +
  facet_wrap(~infection_history) +
  geom_hline(yintercept = 90) +
  geom_hline(yintercept = 82, color = "red")

write.csv(E10, "~/GitHub/Eimeria_Lab/data/Experiment_results/E10_112020_Eim_record.csv")
