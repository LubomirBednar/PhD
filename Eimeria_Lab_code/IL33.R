library(dplyr)

q <- read.csv("C:/Users/exemp/Desktop/E10_E11_qPCRs/E10_E11_CEWE_Eim_qPCR.csv")
IL <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/E11_012021_Eim_IL33.csv")
IL <- select(IL, EH_ID, Gene.expr.IL.33, Gene.expr.MCP.1)

ILq <- merge(q, IL, by = "EH_ID")
ILq88 <- subset(ILq, ILq$challenge_infection == "E88")

model1 <- lm(Gene.expr.IL.33 ~ delta, ILq88)

summary(model1)
