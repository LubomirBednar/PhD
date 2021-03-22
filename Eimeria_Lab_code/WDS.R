# WDS
library(ggplot2)
library(dplyr)
library(tidyverse)
library(magrittr)
library(ggpubr)

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

E6 <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/E6_062018_Eim_weight_oocyst.csv")
E6$challenge_infection <- NA
names(E6)[names(E6) == "Mouse_strain"] <- "mouse_strain"
E6 <- merge(E6, inf, all = T)
E6$X <- NULL
# E10 parasitology

E10 <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/E10_112020_Eim_record.csv")
E10$X <- NULL
E10 <- E10 %>% mutate(challenge_infection = ifelse(batch == "a" , "NA", challenge_infection))
E10 <- E10 %>% mutate(primary_infection = ifelse(batch == "b" , "NA", primary_infection))
E10 <- select(E10, "EH_ID", "labels", "experiment", "dpi", "relative_weight", "mouse_strain", 
              "challenge_infection", "primary_infection", "infection_history")
E10$EH_ID <- as.character(sub("_", "", E10$EH_ID))
E10$mouse_strain <- as.character(sub("PWD", "PWD_PWD", E10$mouse_strain))
E10$mouse_strain <- as.character(sub("SCHUNT", "SCHUNT_SCHUNT", E10$mouse_strain))

# E11 parasitology
E11 <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/E11_012021_Eim_record.csv")
E11$X <- NULL
E11 <- E11 %>% mutate(challenge_infection = ifelse(batch == "a" , "NA", challenge_infection))
E11 <- E11 %>% mutate(primary_infection = ifelse(batch == "b" , "NA", primary_infection))
E11 <- select(E11, "EH_ID", "labels", "experiment", "dpi", "relative_weight", "mouse_strain", 
              "challenge_infection", "primary_infection", "infection_history")
E11$mouse_strain <- as.character(sub("PWD", "PWD_PWD", E11$mouse_strain))
E11$mouse_strain <- as.character(sub("SCHUNT", "SCHUNT_SCHUNT", E11$mouse_strain))

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
  geom_point() +
  geom_smooth() + 
  facet_wrap(~mouse_strain)

ggplot(WDS_challenge, aes(x = dpi, y = relative_weight, color = challenge_infection)) +
  geom_jitter(width = 0.2) +
  geom_smooth() +
  facet_wrap(~mouse_strain)

####
# add immune factors
E7immuno <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/E7_112018_Eim_immuno.csv")
names(E7immuno)[names(E7immuno) == "label"] <- "labels"
E7immuno <- select(E7immuno, EH_ID, dpi, labels, delta, IFNy_CEWE, Eim_MC, EXP_type, CD4, Treg, Div_Treg, Position,
                   Treg17, Th1, Div_Th1, Th17, Div_Th17, CD8, Act_CD8, Div_Act_CD8, IFNy_CD4, IFNy_CD8)
E7immuno$EH_ID <- as.character(sub("_", "", E7immuno$EH_ID))


#needs fixing after delta, IFNy cEWE and such exist (Position is causing havoc so remove it for now)
E11immuno <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/E11_012021_Eim_FACS.csv")
E11immuno <- select(E11immuno, EH_ID, dpi, labels, CD4, Treg, Div_Treg, Position,
                    Treg17, Th1, Div_Th1, Th17, Div_Th17, CD8, Act_CD8, Div_Act_CD8, IFNy_CD4, IFNy_CD8)
E11immuno$delta <- "NA"
E11immuno$IFNy_CEWE <- "NA"
E11immuno$Eim_MC <- "NA"
E11immuno$EXP_type <- "WDS"

E_immuno <- rbind(E7immuno, E11immuno)


WDS <- merge(WDS, E_immuno, all.x = T)

############### try MDS
# select only weight and immune aspects for MDS in the mesentery
WDS_MDS <- subset(WDS, WDS$dpi == 8)
WDS_MDS <- subset(WDS_MDS, WDS_MDS$Position == "mLN")
WDS_MDS <- select(WDS_MDS, EH_ID, relative_weight, CD4, Treg, Div_Treg, Treg17,
                  Th1, Div_Th1, Th17, Div_Th17, CD8, Act_CD8, Div_Act_CD8, IFNy_CD4, IFNy_CD8)
WDS_MDS <- na.omit(WDS_MDS)
WDS_MDS <- WDS_MDS %>% remove_rownames %>% column_to_rownames(var="EH_ID")

distance.matrix <- dist(scale(WDS_MDS,center = T, scale = T),
                              method = "euclidean")
mds.stuff <- cmdscale(distance.matrix, eig = T, x.ret = T)
mds.var.per <- round(mds.stuff$eig/sum(mds.stuff$eig)*100, 1)
mds.values <- mds.stuff$points
mds.data <- data.frame(Sample=rownames(mds.values),
                       X=mds.values[,1],
                       Y=mds.values[,2])
ggplot(data = mds.data, aes(x = X, y = Y, label=Sample))+
  geom_text()+
  theme_bw()+
  xlab(paste("MDS1 - ", mds.var.per[1], "%", sep = ""))+
  ylab(paste("MDS2 - ", mds.var.per[2], "%", sep = ""))+
ggtitle("MDS plot using Euclidean distance")
# add variables back to MDS for graphing
names(mds.data)[names(mds.data) == "Sample"] <- "EH_ID"

WDS_MDS <- subset(WDS, WDS$dpi == 8)
WDS_MDS <- subset(WDS_MDS, WDS_MDS$Position == "mLN")
WDS_MDS <- select(WDS_MDS, EH_ID, challenge_infection, mouse_strain, infection_history, relative_weight, CD4, Treg, Div_Treg, Treg17,
                  Th1, Div_Th1, Th17, Div_Th17, CD8, Act_CD8, Div_Act_CD8, IFNy_CD4, IFNy_CD8)
WDS_MDS <- na.omit(WDS_MDS)
mds.data <- merge(mds.data, WDS_MDS)
# graph with looking at infection isolate, mouse strain, etc
ggplot(data = mds.data, aes(x = X, y = Y, label=EH_ID, color = challenge_infection))+
  geom_point()+
  theme_bw()+
  xlab(paste("MDS1 - ", mds.var.per[1], "%", sep = ""))+
  ylab(paste("MDS2 - ", mds.var.per[2], "%", sep = ""))+
  ggtitle("MDS plot coloured by infection")

ggplot(data = mds.data, aes(x = X, y = Y, label=EH_ID, color = mouse_strain))+
  geom_point()+
  theme_bw()+
  xlab(paste("MDS1 - ", mds.var.per[1], "%", sep = ""))+
  ylab(paste("MDS2 - ", mds.var.per[2], "%", sep = ""))+
  ggtitle("MDS plot coloured by mouse strain")

ggplot(data = mds.data, aes(x = X, y = Y, label=EH_ID, color = infection_history))+
  geom_point()+
  theme_bw()+
  xlab(paste("MDS1 - ", mds.var.per[1], "%", sep = ""))+
  ylab(paste("MDS2 - ", mds.var.per[2], "%", sep = ""))+
  ggtitle("MDS plot coloured by infection history")

################ explore PCA
WDS_PCA <- subset(WDS, WDS$dpi == 8)
WDS_PCA <- subset(WDS_PCA, WDS_PCA$Position == "mLN")
WDS_PCA <- select(WDS_PCA, EH_ID, #challenge_infection, mouse_strain, infection_history, 
                  relative_weight, CD4, Treg, Div_Treg, Treg17,
                  Th1, Div_Th1, Th17, Div_Th17, CD8, Act_CD8, Div_Act_CD8, IFNy_CD4, IFNy_CD8)
WDS_PCA <- na.omit(WDS_PCA)
WDS_PCA <- WDS_PCA %>% remove_rownames %>% column_to_rownames(var="EH_ID")
pca <- prcomp(WDS_PCA, scale = T, center = T)
pca.var <- pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100, 1)
barplot(pca.var.per, main = "Scree Plot", xlab = "Principal component",
        ylab = "Percent Variation")
pca.data <- data.frame(Sample=rownames(pca$x),
                       X = pca$x[,1],
                       Y = pca$x[,2])
names(pca.data)[names(pca.data) == "Sample"] <- "EH_ID"
pca.data <- remove_rownames(pca.data)
WDS_PCA <- subset(WDS, WDS$dpi == 8)
WDS_PCA <- subset(WDS_PCA, WDS_PCA$Position == "mLN")
WDS_PCA <- select(WDS_PCA, EH_ID, infection_history, mouse_strain, challenge_infection, relative_weight, CD4, Treg, Div_Treg, Treg17,
                  Th1, Div_Th1, Th17, Div_Th17, CD8, Act_CD8, Div_Act_CD8, IFNy_CD4, IFNy_CD8)
pca.data <- merge(pca.data, WDS_PCA, by = "EH_ID")
loading_scores <- pca$rotation[,1]
scores <- abs(loading_scores)
score_ranked <- sort(scores, decreasing = T)
top10 <- names(score_ranked[1:10])
pca$rotation[top10,1]

ggplot(data = pca.data, aes(x = X, y = Y, label = EH_ID, color = challenge_infection)) +
  geom_point(size = 3)+
  xlab(paste("PC1 - ", pca.var.per[1], "%", sep = "")) +
  ylab(paste("PC2 - ", pca.var.per[2], "%", sep = "")) +
  theme_bw() +
  ggtitle("PCA by challenge")

ggplot(data = pca.data, aes(x = X, y = Y, label = EH_ID, color = mouse_strain)) +
  geom_point(size = 3)+
  xlab(paste("PC1 - ", pca.var.per[1], "%", sep = "")) +
  ylab(paste("PC2 - ", pca.var.per[2], "%", sep = "")) +
  theme_bw() +
  ggtitle("PCA by mouse strain")

ggplot(data = pca.data, aes(x = X, y = Y, label = EH_ID, color = infection_history)) +
  geom_point(size = 3)+
  xlab(paste("PC1 - ", pca.var.per[1], "%", sep = "")) +
  ylab(paste("PC2 - ", pca.var.per[2], "%", sep = "")) +
  theme_bw() +
  ggtitle("PCA by history")

######################### effect of cell populations on oocysts and weight
WDS_8 <- subset(WDS, WDS$dpi == 8)
WDS_8 <- subset(WDS_8, !is.na(WDS_8$CD4))


                