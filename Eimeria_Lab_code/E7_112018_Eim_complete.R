# E7 complete
# E7 parasitology
E7 <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/E7_112018_Eim_weight_oocyst.csv")
E7$X <- NULL
E7$experiment <- "E7b"
E7 <- select(E7, "EH_ID", "labels", "experiment", "OPG", "dpi", "relative_weight", "mouse_strain", "challenge_infection", "infection_history")
E7$primary_infection <- "NA"
E7$EH_ID <- as.character(sub("_", "", E7$EH_ID))
E7$batch <- "b"

# E6 parasitology
# take E7 infection_history and add to E6
inf <- select(E7, infection_history, EH_ID)
inf <- distinct(inf)

E6 <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/E6_062018_Eim_weight_oocyst.csv")
E6$challenge_infection <- NA
names(E6)[names(E6) == "Mouse_strain"] <- "mouse_strain"
E6 <- merge(E6, inf, all = T)
E6$X <- NULL
E6$experiment <- "E7a"
E6$batch <- "a"

E7 <- rbind(E6, E7)

# add immune factors
E7immuno <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/E7_112018_Eim_immuno.csv")
names(E7immuno)[names(E7immuno) == "label"] <- "labels"
E7immuno <- select(E7immuno, EH_ID, dpi, labels, delta, IFNy_CEWE, Eim_MC, Position, CD4, Treg, Div_Treg, Treg17, Th1,
                   Div_Th1, Th17, Div_Th17, CD8, Act_CD8, Div_Act_CD8, IFNy_CD4, IFNy_CD8)
E7immuno$EH_ID <- gsub( "_", "", E7immuno$EH_ID)


E7 <- merge(E7, E7immuno, all.x = T)
write.csv(E7, "C:/Users/exemp/OneDrive/Documents/GitHub/Eimeria_Lab/data/Experiment_results/E7_112018_Eim_COMPLETE.csv")

#E7gene <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/E7_112018_Eim_CEWE_RT-qPCR.csv")

# E7gene.long <- pivot_longer(E7gene, cols = c("CXCR3", "IRG6", "IL.12"))
# E7gene.long$dpi <- 8
# E7gene.long$EH_ID <- as.character(sub("_", "", E7gene.long$EH_ID))
# E7gene.long <- merge(E7gene.long, E7)
# 
# ggplot(E7gene.long, aes(x = name, y = value, color = challenge_infection)) +
#   geom_boxplot() +
#   geom_smooth() +
#   ggtitle("E7 gene expression")

E7immuno <- select(E7immuno, EH_ID, dpi, labels, delta, IFNy_CEWE, Eim_MC, EXP_type, CD4, Treg, Div_Treg, Position,
                   Treg17, Th1, Div_Th1, Th17, Div_Th17, CD8, Act_CD8, Div_Act_CD8, IFNy_CD4, IFNy_CD8)
E7immuno$EH_ID <- as.character(sub("_", "", E7immuno$EH_ID))
E7immuno$experiment <- "E7"