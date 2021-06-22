# E10 complete
E10_record <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/E10_112020_Eim_record.csv")
E10_record$X <- NULL
E10_oocyst <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/E10_112020_Eim_oocyst.csv")
E10_oocyst$X <- NULL

E10 <- merge(E10_record, E10_oocyst)
# for now this is complete
write.csv(E10, "C:/Users/exemp/OneDrive/Documents/GitHub/Eimeria_Lab/data/Experiment_results/E10_112020_Eim_COMPLETE.csv")

# primary shedding

ggplot(subset(E10, !is.na(E10$primary_infection)), aes(x = dpi, y = OPG, color = primary_infection, group = EH_ID)) +
  geom_point() +
  geom_line() +
  facet_wrap(~EH_ID)

ggplot(subset(E10, !is.na(E10$primary_infection)), aes(x = dpi, y = relative_weight, color = primary_infection)) +
  geom_smooth() +
  geom_jitter() +
  facet_wrap(~primary_infection)     
