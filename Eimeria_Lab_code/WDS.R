# WDS
library(ggplot2)
library(dplyr)
library(tidyverse)
library(magrittr)
library(ggpubr)
library(ggeffects)
library(effects)
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

# CLS
CLS <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/CLS_complete.csv")
CLS$X <- NULL 
CLS$EH_ID <- as.character(sub("_", "", CLS$EH_ID))
names(CLS)[names(CLS) == "Strain"] <- "mouse_strain"
names(CLS)[names(CLS) == "primary"] <- "primary_infection"
names(CLS)[names(CLS) == "challenge"] <- "challenge_infection"
CLSimmuno <- select(CLS, EH_ID, dpi, labels, delta, IFNy_CEWE, Eim_MC, EXP_type, CD4, Treg, Div_Treg, Position,
                    Treg17, Th1, Div_Th1, Th17, Div_Th17, CD8, Act_CD8, Div_Act_CD8, IFNy_CD4, IFNy_CD8)
CLSimmuno <- subset(CLSimmuno, CLSimmuno$dpi == 8)
CLSimmuno <- subset(CLSimmuno, CLSimmuno$Position == "mLN")

CLS <- select(CLS, "EH_ID", "labels", "experiment", "dpi", "relative_weight", "mouse_strain", 
              "challenge_infection", "primary_infection", "infection_history", "OPG")
CLS <- distinct(CLS)

# join them together
WDS <- rbind(E6, E7)
# add mock OPG to bind
E10$OPG <- "NA"
WDS <- rbind(WDS, E10)
E11$OPG <- "NA"
WDS <- rbind(WDS, E11)
WDS <- rbind(WDS, CLS)

WDS[WDS == "NA"] <- NA 
WDS$relative_weight <- as.numeric(WDS$relative_weight)
# add columns with parasite names instead of strains + infection history
WDS$eimeria_species_primary[WDS$primary_infection == "E64"] <- "FER"
WDS$eimeria_species_challenge[WDS$challenge_infection == "E64"] <- "FER"
WDS$eimeria_species_primary[WDS$primary_infection == "E88"] <- "FAL"
WDS$eimeria_species_challenge[WDS$challenge_infection == "E88"] <- "FAL"
WDS$eimeria_species_primary[WDS$primary_infection == "Eflab"] <- "FAL"
WDS$eimeria_species_challenge[WDS$challenge_infection == "Eflab"] <- "FAL"
WDS$eimeria_species_primary[WDS$primary_infection == "E139"] <- "FER"
WDS$eimeria_species_challenge[WDS$challenge_infection == "E139"] <- "FER"
WDS$eimeria_species_primary[WDS$primary_infection == "UNI"] <- "UNI"
WDS$eimeria_species_challenge[WDS$challenge_infection == "UNI"] <- "UNI"

WDS$eimeria_species_history[WDS$infection_history == "UNI:UNI"] <- "UNI:UNI"
WDS$eimeria_species_history[WDS$infection_history == "E88:E88"] <- "FAL:FAL"
WDS$eimeria_species_history[WDS$infection_history == "E88:E64"] <- "FAL:FER"
WDS$eimeria_species_history[WDS$infection_history == "E88:UNI"] <- "FAL:UNI"
WDS$eimeria_species_history[WDS$infection_history == "Eflab:UNI"] <- "FAL:UNI"
WDS$eimeria_species_history[WDS$infection_history == "Eflab:E88"] <- "FAL:FAL"
WDS$eimeria_species_history[WDS$infection_history == "E64:E64"] <- "FER:FER"
WDS$eimeria_species_history[WDS$infection_history == "E64:E88"] <- "FER:FAL"
WDS$eimeria_species_history[WDS$infection_history == "E64:UNI"] <- "FER:UNI"
WDS$eimeria_species_history[WDS$infection_history == "UNI:E88"] <- "UNI:FAL"
WDS$eimeria_species_history[WDS$infection_history == "UNI:E64"] <- "UNI:FER"










WDS_primary <- subset(WDS, !is.na(WDS$primary_infection))
WDS_challenge <- subset(WDS, !is.na(WDS$challenge_infection))
# shuld have a dataset with all weight data and some OPG data

ggplot(WDS_primary, aes(x = dpi, y = relative_weight, color = primary_infection)) +
  geom_jitter(width = 0.2) +
  geom_smooth() + 
  facet_wrap(~mouse_strain)

ggplot(WDS_primary, aes(x = dpi, y = relative_weight, color = primary_infection)) +
  geom_jitter(width = 0.2) +
  geom_smooth() +
  facet_wrap(~primary_infection) +
  ggtitle("WDS weightloss in primary infection")

ggplot(WDS_challenge, aes(x = dpi, y = relative_weight, color = challenge_infection)) +
  geom_jitter(width = 0.2) +
  geom_smooth() +
  facet_wrap(~mouse_strain)

ggplot(WDS_challenge, aes(x = dpi, y = relative_weight, color = challenge_infection)) +
  geom_jitter(width = 0.2) +
  geom_smooth() +
  ggtitle("WDS weightloss in challenge infection")

####
# add immune factors
E7immuno <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/E7_112018_Eim_immuno.csv")
names(E7immuno)[names(E7immuno) == "label"] <- "labels"
E7gene <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/E7_112018_Eim_CEWE_RT-qPCR.csv")

E7gene.long <- pivot_longer(E7gene, cols = c("CXCR3", "IRG6", "IL.12"))
E7gene.long$dpi <- 8
E7gene.long$EH_ID <- as.character(sub("_", "", E7gene.long$EH_ID))
E7gene.long <- merge(E7gene.long, E7)

ggplot(E7gene.long, aes(x = name, y = value, color = challenge_infection)) +
  geom_boxplot() +
  geom_smooth() +
  ggtitle("E7 gene expression")

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
E_immuno <- rbind(E_immuno, CLSimmuno)

WDS_challenge <- merge(WDS_challenge, E_immuno, all.x = T)

############### try MDS
# select only weight and immune aspects for MDS in the mesentery
WDS_MDS <- subset(WDS_challenge, WDS_challenge$dpi == 8)
WDS_MDS <- subset(WDS_MDS, WDS_MDS$Position == "mLN")
WDS_MDS <- select(WDS_MDS, EH_ID, CD4, Treg, Div_Treg, Treg17,
                  Th1, Div_Th1, Th17, Div_Th17, CD8, Act_CD8, Div_Act_CD8, IFNy_CD4, IFNy_CD8)
WDS_MDS <- na.omit(WDS_MDS)
WDS_MDS <- distinct(WDS_MDS)
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

WDS_MDS <- subset(WDS_challenge, WDS_challenge$dpi == 8)
WDS_MDS <- subset(WDS_MDS, WDS_MDS$Position == "mLN")
WDS_MDS <- select(WDS_MDS, EH_ID, experiment, challenge_infection, mouse_strain, infection_history,
                  eimeria_species_challenge, eimeria_species_history,
                  relative_weight, CD4, Treg, Div_Treg, Treg17, Th1, Div_Th1, Th17, Div_Th17, CD8,
                  Act_CD8, Div_Act_CD8, IFNy_CD4, IFNy_CD8)
WDS_MDS <- na.omit(WDS_MDS)
WDS_MDS <- distinct(WDS_MDS)
mds.data <- merge(mds.data, WDS_MDS)
# graph with looking at infection isolate, mouse strain, etc
ggplot(data = mds.data, aes(x = X, y = Y, label=EH_ID, color = eimeria_species_challenge))+
  geom_point(size = 3)+
  xlab(label = "infection history") +
  ylab(label = bquote(atop("maximum weight loss (%)",
                           "lower = higher impact"))) +
  labs(fill = "PCA1 independent predictor") +
  theme_bw()+
  xlab(paste("MDS1 - ", mds.var.per[1], "%", sep = ""))+
  ylab(paste("MDS2 - ", mds.var.per[2], "%", sep = ""))+
  ggtitle("immune cell MDS by infection")

ggplot(data = mds.data, aes(x = X, y = Y, label=EH_ID, color = mouse_strain))+
  geom_point()+
  theme_bw()+
  xlab(paste("MDS1 - ", mds.var.per[1], "%", sep = ""))+
  ylab(paste("MDS2 - ", mds.var.per[2], "%", sep = ""))+
  ggtitle("MDS plot coloured by mouse strain")

ggplot(data = mds.data, aes(x = X, y = Y, label=EH_ID, color = eimeria_species_history))+
  geom_point()+
  theme_bw()+
  xlab(paste("MDS1 - ", mds.var.per[1], "%", sep = ""))+
  ylab(paste("MDS2 - ", mds.var.per[2], "%", sep = ""))+
  ggtitle("MDS plot coloured by infection history")
# and check for batch effect
ggplot(data = mds.data, aes(x = X, y = Y, label=EH_ID, color = experiment))+
  geom_point(size = 3)+
  theme_bw()+
  xlab(paste("MDS1 - ", mds.var.per[1], "%", sep = ""))+
  ylab(paste("MDS2 - ", mds.var.per[2], "%", sep = ""))+
  ggtitle("immune cell MDS by infection")

################ explore PCA
WDS_PCA <- subset(WDS_challenge, WDS_challenge$dpi == 8)
WDS_PCA <- subset(WDS_PCA, WDS_PCA$Position == "mLN")
WDS_PCA <- select(WDS_PCA, EH_ID, #challenge_infection, mouse_strain, infection_history, 
                  CD4, Treg, Div_Treg, Treg17,
                  Th1, Div_Th1, Th17, Div_Th17, CD8, Act_CD8, Div_Act_CD8, IFNy_CD4, IFNy_CD8)
WDS_PCA <- na.omit(WDS_PCA)
WDS_PCA <- distinct(WDS_PCA)
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
WDS_PCA <- subset(WDS_challenge, WDS_challenge$dpi == 8)
WDS_PCA <- subset(WDS_PCA, WDS_PCA$Position == "mLN")
WDS_PCA <- select(WDS_PCA, EH_ID, experiment, challenge_infection, mouse_strain, infection_history,
                  eimeria_species_challenge, eimeria_species_history,
                  relative_weight, CD4, Treg, Div_Treg, Treg17, Th1, Div_Th1, Th17, Div_Th17, CD8,
                  Act_CD8, Div_Act_CD8, IFNy_CD4, IFNy_CD8)
pca.data <- merge(pca.data, WDS_PCA, by = "EH_ID")
loading_scores <- pca$rotation[,1]
scores <- abs(loading_scores)
score_ranked <- sort(scores, decreasing = T)
top10 <- names(score_ranked[1:13])
pca$rotation[top10,1]
pca.var.per.df <- cbind(pca.var.per, top10)
pca.var.per.df <- data.frame(pca.var.per.df)
pca.var.per.df$pca.var.per <- as.numeric(as.character(pca.var.per.df$pca.var.per))
pca.var.per.df$top10 <- factor(pca.var.per.df$top10,levels = top10)

ggplot(pca.var.per.df, aes(top10, pca.var.per, fill = top10)) +
  geom_col() +
  theme(strip.text.x = element_text(size = 14, face = "bold"),
        strip.text.y = element_text(size = 14, face = "bold")) +
  scale_x_discrete(name ="Cell populations") +
  scale_y_discrete(name ="% of variation explained") +
  labs(fill = "") +
  
  geom_text(aes(label = pca.var.per)) +
  theme_bw()

ggplot(data = pca.data, aes(x = X, y = Y, label = EH_ID, color = eimeria_species_challenge)) +
  geom_point(size = 3)+
  xlab(paste("PC1 - ", pca.var.per[1], "%", sep = "")) +
  ylab(paste("PC2 - ", pca.var.per[2], "%", sep = "")) +
  theme_bw() +
  theme(axis.text=element_text(size=12, face = "bold"),
        title = element_text(size = 16, face = "bold"),
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        strip.text.y = element_text(size = 14, face = "bold"),
        legend.text = element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"))+
  ggtitle("T-cell populations PCA by challenge infection") +
  labs(color='infection species') 
# looks good, test for batch independence, add experiment column
ggplot(data = pca.data, aes(x = X, y = Y, label = EH_ID, color = experiment)) +
  geom_point(size = 3)+
  xlab(paste("PC1 - ", pca.var.per[1], "%", sep = "")) +
  ylab(paste("PC2 - ", pca.var.per[2], "%", sep = "")) +
  theme_bw() +
  theme(axis.text=element_text(size=12, face = "bold"),
        title = element_text(size = 16, face = "bold"),
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        strip.text.y = element_text(size = 14, face = "bold"),
        legend.text = element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"))+
  ggtitle("Checking for batch effect") +
  labs(color='Experiment ID') 


# no clustering, moving on
ggplot(data = pca.data, aes(x = X, y = Y, label = EH_ID, color = eimeria_species_history)) +
  geom_point(size = 3)+
  xlab(paste("PC1 - ", pca.var.per[1], "%", sep = "")) +
  ylab(paste("PC2 - ", pca.var.per[2], "%", sep = "")) +
  theme_bw() +
  ggtitle("PCA by history")

######################### effect of cell populations on oocysts and weight
pca.data$primary <- unlist(lapply(str_split(string = pca.data$infection_history, pattern = ":"), "[[", 1))
pca.data$secondary <- unlist(lapply(str_split(string = pca.data$infection_history, pattern = ":"), "[[", 2))
pca.data$primary_species <- unlist(lapply(str_split(string = pca.data$eimeria_species_history, pattern = ":"), "[[", 1))
pca.data$secondary_species <- unlist(lapply(str_split(string = pca.data$eimeria_species_history, pattern = ":"), "[[", 2))

pca.data$infection_history <- factor(pca.data$infection_history,
                                      levels = c("UNI:UNI", "E88:E88", "E88:E64", "E88:UNI", "Eflab:UNI",
                                                 "Eflab:E88", "E64:E64", "E64:E88", "E64:UNI", "UNI:E88", "UNI:E64"))

pca.data$eimeria_species_history <- factor(pca.data$eimeria_species_history,
                                      levels = c("UNI:UNI", "FAL:FAL", "FAL:FER", "FAL:UNI", "FER:FER", 
                                                 "FER:FAL", "FER:UNI", "UNI:FAL", "UNI:FER"))


mod1 <- lm(relative_weight~primary_species*secondary_species+X+Y, pca.data)
mod2 <- lm(relative_weight~secondary_species+X+Y, pca.data)
mod3 <- lm(relative_weight~X+Y, pca.data)
                
summary(mod1) # not enough power, groups go NA = rank deficient
summary(mod2) # strong X effect
summary(mod3) # strong X effect

# this is very cool but we're using relative weight from day 8, which won't tell us much
# make maximum OPG (OPG when available) and maximum weightloss columns ad test those
# confirmed, no primary challenges
max_weight_primary <- WDS_primary %>% select(EH_ID, relative_weight, primary_infection, 
                                               challenge_infection, infection_history) 
max_weight_challenge <- WDS_challenge %>% select(EH_ID, relative_weight, primary_infection, 
                                               challenge_infection, infection_history)
# except that this uses only data from challenge for primary:uni
max_weight_primary <- max_weight_primary[!is.na(max_weight_primary$relative_weight), ]
max_weight_challenge <- max_weight_challenge[!is.na(max_weight_challenge$relative_weight),]
# creates weird 1 sample with 0 so remove it
# max_weight_challenge <- max_weight_challenge[ !(max_weight_challenge$EH_ID %in% "LM0350"),]
# max_OPG<- WDS %>% select(EH_ID, OPG)

# max_OPG$OPG <- as.numeric(as.character(max_OPG$OPG))
max_weight_primary <- max_weight_primary %>% dplyr::group_by(EH_ID) %>% dplyr::summarise(relative_weight = min(relative_weight, na.rm = T))
max_weight_challenge <- max_weight_challenge %>% dplyr::group_by(EH_ID) %>% dplyr::summarise(relative_weight = min(relative_weight, na.rm = T))
# max_OPG <- max_OPG %>% dplyr::group_by(EH_ID) %>% dplyr::summarise(OPG = max(OPG, na.rm = T))
names(max_weight_primary)[names(max_weight_primary) == "relative_weight"] <- "maximum_weight_loss_primary"
names(max_weight_challenge)[names(max_weight_challenge) == "relative_weight"] <- "maximum_weight_loss_challenge"
# now this gets complicated. infection_history needs maximum weightloss, but tor eflect primary infections which
# are in challenge. and current infections need chalenge weight loss. Otherwise infected:uninfected mess up models

pca.data1 <- merge(pca.data, max_weight_primary, by = "EH_ID")
pca.data1 <- merge(pca.data1, max_weight_challenge, by = "EH_ID")

# anything that is infected:uninfected needs to be the in the final maximum weightloss column from primary_max_WLoss
# can do min across both columns but that didnt work because its non disciminating to our criteria so scrapping
# pca.challenge <- transform(pca.challenge, min = pmin(maximum_weight_loss_primary, maximum_weight_loss_challenge))
############## THIS WILL ONLY WORK FOR infection_history MODELS! (pca.data.swap1)
pca.data.swap1 <- subset(pca.data1, infection_history == "E64:UNI" | 
                          infection_history == "E88:UNI"|
                          infection_history == "Eflab:UNI")
# drop challenge weightloss column to avoic UNI weight in infection history
pca.data.swap1$maximum_weight_loss_challenge <- NULL
names(pca.data.swap1)[names(pca.data.swap1) == "maximum_weight_loss_primary"] <- "maximum_weight_loss"
# name the columns the same and merge
pca.data.swap2 <- subset(pca.data1, infection_history == "UNI:UNI" | 
                      infection_history == "E88:E88"|
                      infection_history == "E88:E64"|
                      infection_history == "Eflab:E88"|
                      infection_history == "E64:E64"|
                      infection_history == "E64:E88"|
                      infection_history == "UNI:E88"|
                      infection_history == "UNI:E64")
pca.data.swap2$maximum_weight_loss_primary <- NULL
names(pca.data.swap2)[names(pca.data.swap2) == "maximum_weight_loss_challenge"] <- "maximum_weight_loss"

pca.data.swap <- rbind(pca.data.swap1, pca.data.swap2)

# for this the FAL:UNI and UNI:FAL are the same, join under new column and model
pca.data.swap$relevant_history[pca.data.swap$eimeria_species_history == "UNI:UNI"] <- "UNI"
pca.data.swap$relevant_history[pca.data.swap$eimeria_species_history == "UNI:FAL"] <- "FAL"
pca.data.swap$relevant_history[pca.data.swap$eimeria_species_history == "FAL:UNI"] <- "FAL"
pca.data.swap$relevant_history[pca.data.swap$eimeria_species_history == "FER:UNI"] <- "FER"
pca.data.swap$relevant_history[pca.data.swap$eimeria_species_history == "UNI:FER"] <- "FER"
pca.data.swap$relevant_history[pca.data.swap$eimeria_species_history == "FAL:FAL"] <- "FAL:FAL"
pca.data.swap$relevant_history[pca.data.swap$eimeria_species_history == "FER:FER"] <- "FER:FER"
pca.data.swap$relevant_history[pca.data.swap$eimeria_species_history == "FER:FAL"] <- "FER:FAL"
pca.data.swap$relevant_history[pca.data.swap$eimeria_species_history == "FAL:FER"] <- "FAL:FER"

pca.data.swap$relevant_history <- factor(pca.data.swap$relevant_history,
                                         levels = c("UNI", "FAL", "FER","FER:FER", 
                                                    "FER:FAL", "FAL:FAL", "FAL:FER"))

mod4 <- lm(maximum_weight_loss~relevant_history+X+Y, pca.data.swap) 
summary(mod4)
# percentage error
sigma(mod4)*100/mean(pca.data.swap$maximum_weight_loss)

mod4gg <- ggpredict(mod4, c("relevant_history", "X"))
plot(mod4gg) + 
  geom_point() +
  xlab(label = "infection history") +
  ylab(label = bquote(atop("maximum weight loss (%)",
                   "lower = higher impact"))) +
  theme_bw() +
  ggtitle("Predicted weight loss by infection history")


mod4.1 <- lm(maximum_weight_loss~relevant_history*X, pca.data.swap) 
summary(mod4.1)

mod4.1gg <- ggpredict(mod4.1, c("relevant_history", "X"))
plot(mod4.1gg) + 
  geom_point() +
  xlab(label = "infection history") +
  ylab(label = bquote(atop("maximum weight loss (%)",
                           "lower = higher impact"))) +
  labs(fill = "PCA1 independent predictor") +
  theme_bw() +
  theme(axis.text = element_text(size = 15), 
        axis.title = element_text(size = 18), 
        plot.title = element_text(size = 18),
        axis.text.x = element_text()) +
  ggtitle("Predicted weight loss by infection history")



ggplot(mod4.1, aes(y = predict(mod4.1,terms = c("relevant_history", "X")), 
                   x = X, colour = relevant_history, group = relevant_history)) +
  # geom_jitter(size = 2,width =1, aes(y = predict(mod4.1,terms = c("relevant_history", "X"), 
  # colour = group ,pca.data.swap))) +
  geom_jitter(width = 1)+
  geom_smooth() +
  theme_bw() + 
  ggtitle("") + 
  xlab("PCA1 coordinates") + 
  ylab(bquote(atop("maximum weight loss %",
                   "lower = worse"))) +
  theme(axis.text = element_text(size = 15), 
        axis.title = element_text(size = 18), 
        plot.title = element_text(size = 18),
        legend.text = element_text(size = 12),
        legend.title = element_blank(),
        axis.text.x = element_text())


mod4.3 <- lm(maximum_weight_loss~mouse_strain + X, pca.data.swap)
summary(mod4.3)
mod4.3gg <- ggpredict(mod4.3, c("mouse_strain", "X"))
plot(mod4.3gg) + 
  geom_point() +
  xlab(label = "mouse_strain") +
  ylab(label = bquote(atop("maximum weight loss (%)",
                           "lower = higher impact"))) +
  labs(fill = "PCA1 independent predictor") +
  theme_bw() +
  theme(axis.text = element_text(size = 15), 
        axis.title = element_text(size = 18), 
        plot.title = element_text(size = 18),
        axis.text.x = element_text()) +
  ggtitle("Predicted weight loss by infection history")





# interesting, check strains that might be pulling it


# continue challenge analysis with secondary only
# model 5 exploration

mod5 <- lm(maximum_weight_loss_challenge~secondary+X, pca.data1) 
summary(mod5)
mod5gg <- ggpredict(mod5, c("secondary", "X"))
plot(mod5gg)


ggplot(pca.data1, aes(y = maximum_weight_loss_challenge, x = secondary, colour = X)) +
  geom_boxplot(aes(y = predict(mod5,terms = c("secondary", "X"), pca.data1))) +
  geom_jitter(width = 0.2, size = 3) +
  theme_bw() + 
  ggtitle("") + 
  xlab("challenge infection") + 
  ylab("maximum weightloss") +
  theme(axis.text = element_text(size = 15), 
        axis.title = element_text(size = 18), 
        plot.title = element_text(size = 18))
# mod6 <- lm(maximum_weight_loss_cxhallenge~X+Y, pca.data1)
# summary(mod6)
# mod6gg <- ggpredict(mod6, c("X"))
# plot(mod6gg)
# add comparisons





mod5gg1 <- ggeffect(mod5, c("secondary", "X")) 
# ggplot(mod5gg1, aes(x, predicted, color = group)) +
#   geom_boxplot() +
#   geom_point() +
#   stat_compare_means(comparisons = my_comparisons,
#                      method = "wilcox.test", 
#                      aes(label = ..p.signif..), 
#                      size = 3, label.y.npc =0.95) +
#   geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1) +
#   xlab("P. berghei dataset node degree")
  
# make mod 4 intercept on uninfected
mod7 <- lm(maximum_weight_loss~infection_history+X+Y, pca.challenge) 
summary(mod7)
mod7gg <- ggpredict(mod7, c("infection_history", "X"))
plot(mod7gg)

# add raw cell counts to the  pca.challenge and test that
immuno <- select(WDS_8, Eim_MC, Position, EH_ID, CD4, Treg, Div_Treg, Treg17,
                 Th1, Div_Th1, Th17, Div_Th17, CD8, Act_CD8, Div_Act_CD8, IFNy_CD4, IFNy_CD8)
immuno <- subset(immuno, immuno$Position == "mLN")

pca.challenge <- merge(pca.challenge, immuno)

mod8 <- lm(maximum_weight_loss~infection_history+Th1, pca.challenge)
summary(mod8)
mod8gg <- ggpredict(mod8, c("infection_history", "Th1"))
plot(mod8gg)

ggimmuno <- pivot_longer(immuno, cols = c("Div_Treg", "Treg", "CD4", "Treg17", "Th1", "Div_Th1", "Th17", "Div_Th17",
                                          "CD8", "Act_CD8", "Div_Act_CD8", "IFNy_CD4", "IFNy_CD8"))
ggpara <- select(pca.challenge, EH_ID, infection_history, mouse_strain, primary, secondary, maximum_weight_loss, X, Y)
ggimmuno <- merge(ggimmuno, ggpara)

ggplot(ggimmuno, aes(secondary, value, color = secondary)) +
  geom_boxplot() +
  geom_jitter(width = 0.2) +
  facet_wrap(~name, scales = "free")

ggplot(ggimmuno, aes(maximum_weight_loss, value, group = infection_history, color = infection_history)) +
  geom_point(size = 2) +
  geom_smooth() +
  facet_wrap(~name, scales = "free")

mod10 <- lm(maximum_weight_loss~Eim_MC+X+Y, subset(pca.challenge, !is.na(pca.challenge$Eim_MC)))
summary(mod10)

ggplot(subset(pca.challenge, !is.na(pca.challenge$Eim_MC)), aes(y = maximum_weight_loss, 
                          x = X, fill = Eim_MC, colour = Eim_MC)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_smooth(aes(y = predict(mod10, pca.challenge)), se = F) +
  theme_bw() + 
  ggtitle("") + 
  xlab("X") + 
  ylab("maximum weightloss") +
  theme(axis.text = element_text(size = 15), 
        axis.title = element_text(size = 18), 
        plot.title = element_text(size = 18))

ggplot(subset(ggimmuno, !ggimmuno$Eim_MC == "NA"), aes(x = maximum_weight_loss, y = value, colour = Eim_MC)) +
  geom_point(size = 2, alpha = 0.7) +
  theme_bw() + 
  ggtitle("") + 
  xlab("maximum weightloss") + 
  ylab("cell populations %") +
  facet_wrap(~name, scales = "free") +
  geom_smooth() +
  theme(axis.text = element_text(size = 15), 
        axis.title = element_text(size = 18), 
        plot.title = element_text(size = 18))
# look st mosue strains

mod4.1 <- lm(maximum_weight_loss~mouse_strain+X+Y, pca.challenge)
summary(mod4.1)
gg1.1 <- ggpredict(mod4.1, c("mouse_strain", "X"))
plot(gg1.1)
# # check health of the data, redundant now after fixing
# mod4_matrix <- model.matrix(maximum_weight_loss~primary*secondary+X+Y, pca.data1)
# length(mod4$coefficients) < mod4$rank
# library(Matrix)
# cat(rankMatrix(mod4_matrix),"\n")
# cat(rankMatrix(mod4_matrix), "\n")

# make model for homologou, heterologous, unprotected and naive
pca.challenge$infection_type = case_when(
   pca.challenge$infection_history == "E88:E88" ~ "homologous",
   pca.challenge$infection_history == "E88:E64" ~ "heterologous",
   pca.challenge$infection_history == "E88:UNI" ~ "unprotected",
   pca.challenge$infection_history == "E64:E64" ~ "homologous",
   pca.challenge$infection_history == "E64:E88" ~ "heterologous",  
   pca.challenge$infection_history == "E64:UNI" ~ "unprotected",
   pca.challenge$infection_history == "UNI:E64" ~ "unprotected",
   pca.challenge$infection_history == "UNI:E88" ~ "unprotected",
   pca.challenge$infection_history == "UNI:UNI" ~ "anaive"
   )

mod7 <- lm(maximum_weight_loss~infection_type+X+Y, pca.challenge)
summary(mod7)

gg1 <- ggpredict(mod7, c("infection_type", "X"))
plot(gg1)

ggplot(gg1, aes(x, predicted, colour = group)) + 
  geom_line() + 
  geom_jitter(width = 0.2)


gg2 <- ggpredict(mod7, terms = c("infection_type"))
gg3 <- plot(gg2) 
my_comparisons <-  list(c("homologous", "anaive"), 
                        c("heterologous", "anaive"), 
                        c("homologous" , "heterologous"),
                        c("homologous", "unprotected"),
                        c("heterologous", "unprotected"),
                        c("unprotected", "anaive"))

