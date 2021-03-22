# WDS

library(ggplot2)
library(dplyr)
library(tidyverse)

# E7 parasitology
E7 <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/E7_112018_Eim_weight_oocyst.csv")
E7$X <- NULL
E7$experiment <- "E7"
E7 <- select(E7, "EH_ID", "labels", "experiment", "OPG", "dpi", "relative_weight", "mouse_strain", "challenge_infection", "infection_history")
E7$primary_infection <- "NA"
E7$EH_ID <- as.character(sub("_", "", E7$EH_ID))

# E6 parasitology
# take E7 infection_history and add to E6
inf <- select(E7, infection_history, EH_ID)
inf <- distinct(inf)

E6 <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/E5_062018_Eim_weight_oocyst.csv")
E6 <- select(E6, "EH_ID", "labels", "experiment", "OPG", "dpi", "relative_weight", "Strain", "primary_infection")
E6 <- distinct(E6)
names(E6)[names(E6) == "Strain"] <- "mouse_strain"
names(E6)[names(E6) == "Batch"] <- "batch"
E6a <- subset(E6, E6$experiment == "Expe005_2a")
E6b <- subset(E6, E6$experiment == "Expe005_2b")
E6 <- rbind(E6a, E6b)
E6$experiment <- "E6"
E6$challenge_infection <- "NA"
E6 <- merge(E6, inf, all.x = T)

E6 <- E6 %>% separate(infection_history, c("primary_infection", "challenge_infection"), sep = ":")
E6 <- within(E6, infection_history <- paste(primary_infection,challenge_infection, sep=":"))
E6$challenge_infection <- NA

# E10 parasitology

E10 <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/E10_112020_Eim_record.csv")
E10$X <- NULL
E10 <- E10 %>% mutate(challenge_infection = ifelse(batch == "a" , "NA", challenge_infection))
E10 <- E10 %>% mutate(primary_infection = ifelse(batch == "b" , "NA", primary_infection))
E10 <- select(E10, "EH_ID", "labels", "experiment", "dpi", "relative_weight", "mouse_strain", 
              "challenge_infection", "primary_infection", "infection_history")
E10$EH_ID <- as.character(sub("_", "", E10$EH_ID))


# E11 parasitology
E11 <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/E11_012021_Eim_record.csv")
E11$X <- NULL
E11 <- E11 %>% mutate(challenge_infection = ifelse(batch == "a" , "NA", challenge_infection))
E11 <- E11 %>% mutate(primary_infection = ifelse(batch == "b" , "NA", primary_infection))
E11 <- select(E11, "EH_ID", "labels", "experiment", "dpi", "relative_weight", "mouse_strain", 
              "challenge_infection", "primary_infection", "infection_history")
# join them together
WDS <- rbind(E6, E7)
# add mock OPG to bind
E10$OPG <- "NA"
WDS <- rbind(WDS, E10)
E11$OPG <- "NA"
WDS <- rbind(WDS, E11)

WDS[WDS == "NA"] <- NA 
WDS$relative_weight <- as.numeric(WDS$relative_weight)

WDS_primary <- subset(WDS, !is.na(WDS$primary_infection))
WDS_challenge <- subset(WDS, !is.na(WDS$challenge_infection))
# shuld have a dataset with all weight data and some OPG data

ggplot(WDS_primary, aes(x = dpi, y = relative_weight, color = primary_infection)) +
  geom_jitter(width = 0.2) +
  geom_smooth() +
  facet_wrap(~experiment)

ggplot(WDS_challenge, aes(x = dpi, y = relative_weight, color = challenge_infection)) +
  geom_jitter(width = 0.2) +
  geom_smooth() +
  facet_wrap(~infection_history)

