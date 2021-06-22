# E10 oocysts
E10a <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/E10a_112020_Eim_oocyst.csv")
E10a$X <- NULL
E10a$batch <- "a"
E10b <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/E10b_112020_Eim_oocyst.csv")
E10b$X <- NULL
E10b$batch <- "b"

E10 <- rbind(E10a, E10b)
# dilution is always 1mL in processing so set that flat
E10$dilution <- 1
# calculcate OPG
E10$total_oocysts <- ((E10$oocyst_sq1 
                      + E10$oocyst_sq2 
                      + E10$oocyst_sq3 
                      + E10$oocyst_sq4) / 4) * 
  10000 * # because volume chamber
  E10$dilution
# need record for feces weight
E10_record <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/E10_112020_Eim_record.csv")
E10_record$X <- NULL

E10_oocyst <- merge(E10, E10_record, all.y = T)


E10_oocyst$OPG <- E10_oocyst$total_oocysts / E10_oocyst$feces_weight 
# comply with documentation labels (unique timepoint identifier), experiment (unique experiment identifier), 
# oocyst_sq1, oocys_sq2, oocyst_sq3, oocyst_sq4, oocyst_mean, OPG, dilution.
E10_oocyst <- select(E10_oocyst, labels, experiment, oocyst_sq1, oocyst_sq2, oocyst_sq3, oocyst_sq4, oocyst_mean, OPG, dilution)

write.csv(E10_oocyst, "C:/Users/exemp/OneDrive/Documents/GitHub/Eimeria_Lab/data/Experiment_results/E10_112020_Eim_oocyst.csv")
