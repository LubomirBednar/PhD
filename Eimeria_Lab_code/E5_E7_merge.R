library(dplyr)
# weight and oocysts
E5O <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/E5_062018_Eim_oocyst.csv")
E5O$EH_ID <- NULL
E5O$OPG <- NULL
E5O$X <- NULL
E5W <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/E5_062018_Eim_record.csv")
E5W$weightloss <- NULL
E5W$X <- NULL

E7O <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/E7_112018_Eim_oocyst.csv")
E7O$OPG <- NULL
E7W <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/E7_112018_Eim_Record.csv")
E7W$EH_ID <- gsub("_", "", E7W$EH_ID)

EO <- rbind(E5O, E7O)
EO$oocyst_mean <- rowMeans(EO[,c("oocyst_sq1", "oocyst_sq2", "oocyst_sq3", "oocyst_sq4")], na.rm=TRUE)


EW <- rbind(E5W, E7W)

write.csv(EO, "C:/Users/exemp/OneDrive/Documents/GitHub/Eimeria_Lab/data/Experiment_results/E5_062018_Eim_oocyst.csv")
write.csv(EW, "C:/Users/exemp/OneDrive/Documents/GitHub/Eimeria_Lab/data/Experiment_results/E5_062018_Eim_record.csv")
# design
E5D <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experimental_design/E5_062018_Eim_DESIGN.csv")
E5D$experiment <- 
E7D <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experimental_design/E7_112018_Eim_DESIGN.csv")


ED <- merge(E5D, E7D, all.x = T)


# qPCRs
q <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/E7_112018_Eim_CEWE_qPCR.csv")
