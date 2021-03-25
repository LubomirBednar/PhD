library(Hmisc)
library(tidyverse)
library(dplyr)
library(corrplot)
library(PerformanceAnalytics)

# combined gene expression
complete <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/E7_P3_P4_Eim_complete.csv")
complete$X <- NULL
complete$label.1 <- NULL
complete <- distinct(complete)
complete <- subset(complete, !is.na(IFNy_CEWE))
complete$Position <- replace_na(complete$Position, "mLN")
complete <- subset(complete, Position!="Spleen")
# rewrite MCs from previous experiment as per observations (amp + MC) 334, 335, 340
complete$Eim_MC[3] <- "neg"
complete$Eim_MC[4] <- "neg"
complete$Eim_MC[9] <- "neg"
complete$Eim_MC[6] <- "neg"
complete$Eim_MC[15] <- "neg"
# subset for bad elisa
complete <- data.frame(complete)
complete2 <- complete[-c(70:94),]
complete1 <- complete
# if there are oocysts, it is positive
complete$Eim_MC[complete$OPG > 0] <- "infected"
complete$Eim_MC[complete$Eim_MC == "pos"] <- "infected"
complete$Eim_MC[complete$Eim_MC == "neg"] <- "uninfected"
complete1$Eim_MC[complete1$Eim_MC == "pos"] <- "infected"
complete1$Eim_MC[complete1$Eim_MC == "neg"] <- "uninfected"
# lab IFNy
complete1$Eimeria[complete1$challenge == "E64"] <- "E.ferrisi"
complete1$Eimeria[complete1$challenge == "E88"] <- "E.falciformis"
complete1$Eimeria[complete1$challenge == "UNI"] <- "Uninfected"

complete1$Eimeria[complete1$Eim_MC == "uninfected"] <- "Uninfected"

E7 <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/E7_112018_Eim_CEWE_RT-qPCR.csv")
E7$X <- NULL
E7$EXP_type <- "WDS"
P3 <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/P3_112019_Eim_CEWE_RTqPCR.csv")
P3$X <- NULL
P3$EXP_type <- "CLS"
names(E7)[names(E7) == "Mouse_ID"] <- "EH_ID"

lab_RT <- rbind(E7, P3)
RT_merge <- select(complete1, EH_ID, delta, Eim_MC, Eimeria)
lab_RT <- merge(lab_RT, RT_merge)

HZ18 <- read.csv("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Field/master/data/Gene_expression/HZ18_complete.csv")
HZ18 <- select(HZ18, Mouse_ID, Target, deltaCtMmE_tissue, Eimeria.subspecies, NE, inf)
HZ18$inf[HZ18$inf == "TRUE"] <- "infected"
HZ18$inf[HZ18$inf == "FALSE"] <- "uninfected"
# rename columns
names(HZ18)[names(HZ18) == "deltaCtMmE_tissue"] <- "delta"

HZ18 <- pivot_wider(HZ18, names_from = Target, values_from = NE)

# name columns to combine, drop genes not in common
names(HZ18)[names(HZ18) == "Mouse_ID"] <- "EH_ID"
names(HZ18)[names(HZ18) == "IL-12b"] <- "IL.12"
HZ18$GBP2 <- NULL
HZ18$`IL-6` <- NULL
names(HZ18)[names(HZ18) == "inf"] <- "Eim_MC"
names(HZ18)[names(HZ18) == "Eimeria.subspecies"] <- "Eimeria"
HZ18$Eimeria[HZ18$Eimeria == "non infected"] <- "Uninfected"
HZ18$EXP_type <- "wild"

foo <- rbind(HZ18, lab_RT)
# finally, start working on correlations
foo <- column_to_rownames(.data = foo, var = "EH_ID")

foo_num <- foo[,c(1,4,5,6)]
res <- cor(foo_num, use = "complete.obs")
round(res, 2)
res2 <- rcorr(as.matrix(foo_num))
# borrowing this lovely function from STHDA
# ++++++++++++++++++++++++++++
# flattenCorrMatrix
# ++++++++++++++++++++++++++++
# cormat : matrix of the correlation coefficients
# pmat : matrix of the correlation p-values
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}

flattenCorrMatrix(res2$r, res2$P)
symnum(res, abbr.colnames = F)

corrplot(res, type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45)

corrplot(res2$r, type="upper", order="hclust", 
         p.mat = res2$P, sig.level = 0.01, insig = "blank")

chart.Correlation(foo_num, histogram=TRUE, pch=19)

col<- colorRampPalette(c("blue", "white", "red"))(20)
heatmap(x = res, col = col, symm = TRUE)
