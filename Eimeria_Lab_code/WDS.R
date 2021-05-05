# WDS
library(ggplot2)
library(dplyr)
library(tidyverse)
library(magrittr)
library(ggpubr)
library(ggeffects)
library(effects)
library(car)
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
CLSimmuno <- select(CLS, EH_ID, experiment, dpi, labels, delta, IFNy_CEWE, Eim_MC, EXP_type, CD4, Treg, Div_Treg, Position,
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
E7immuno$experiment <- "E7"
#needs fixing after delta, IFNy cEWE and such exist (Position is causing havoc so remove it for now)
E11immuno <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/E11_012021_Eim_FACS.csv")


E11immuno <- select(E11immuno, EH_ID, dpi, labels, experiment, CD4, Treg, Div_Treg, Position,
                  Treg17, Th1, Div_Th1, Th17, Div_Th17, CD8, Act_CD8, Div_Act_CD8, IFNy_CD4, IFNy_CD8)

E11immuno$delta <- "NA"
E11immuno$IFNy_CEWE <- "NA"
E11immuno$Eim_MC <- "NA"
E11immuno$EXP_type <- "WDS"

E_immuno <- rbind(E7immuno, E11immuno)
E_immuno <- rbind(E_immuno, CLSimmuno)

WDS_challenge <- merge(WDS_challenge, E_immuno, all.x = T)

ggimmuno <- pivot_longer(E_immuno, cols = c("Div_Treg", "Treg", "CD4", "Treg17", "Th1", "Div_Th1", "Th17", "Div_Th17",
                                          "CD8", "Act_CD8", "Div_Act_CD8", "IFNy_CD4", "IFNy_CD8"))

ggplot(ggimmuno, aes(name, value, color = experiment)) +
  geom_boxplot() + 
  facet_wrap(~Eim_MC)

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


pca.data$mouse_strain <- factor(pca.data$mouse_strain,
                                     levels = c("SCHUNT_SCHUNT", "SWISS", "PWD_PWD", "BUSNA_STRA", "STRA_BUSNA",
                                                "PWD_SCHUNT", "STRA_STRA", "STRA_SCHUNT", "PWD_BUSNA", "SCHUNT_PWD",
                                                "SCHUNT_STRA", "BUSNA_BUSNA", "BUSNA_PWD"))

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
# organize to compare to uninfected
pca.data.swap$relevant_history <- factor(pca.data.swap$relevant_history,
                                       levels = c("UNI", "FAL", "FER","FER:FER", 
                                                  "FER:FAL", "FAL:FAL", "FAL:FER"))
# organize to compare to SWISS
strains <- unique(pca.data.swap$mouse_strain)
pca.data.swap$mouse_strain <- factor(pca.data.swap$mouse_strain,
                                       levels = strains)


# check distribution for lm or glm

############################### make basic models and check for residuals with anova + AIC
library(AICcmodavg)
############## most basic ones
# without X and only history
wloss_history_mod <- lm(maximum_weight_loss~relevant_history, pca.data.swap)
summary(wloss_history_mod) # gives strong effect on primary infections and X (immune aspects)
Anova(wloss_history_mod) # confirmed, without significant residuals
sigma(wloss_history_mod)*100/mean(pca.data.swap$maximum_weight_loss)# percentage error =  6.315369

# without X and only mouse strain 
wloss_history_mouse_mod <- lm(maximum_weight_loss~relevant_history + mouse_strain , pca.data.swap)
summary(wloss_history_mouse_mod) # strain effects remain more or less + strong infection effect
Anova(wloss_history_mouse_mod) # no residual significance, history has strongest effect
sigma(wloss_history_mouse_mod)*100/mean(pca.data.swap$maximum_weight_loss)# percentage error = 5.703082


# weightloss and relevant history with pca components
wloss_history_X_Y_mod <- lm(maximum_weight_loss~relevant_history+X+Y, pca.data.swap) 
summary(wloss_history_X_Y_mod) # gives strong effect on primary infections and X (immune aspects)
Anova(wloss_history_X_Y_mod) # confirmed, without significant residuals
sigma(wloss_history_X_Y_mod)*100/mean(pca.data.swap$maximum_weight_loss)# percentage error = 6.084542

# weightloss and relevant history with X as interaction
wloss_history_xX_mod <- lm(maximum_weight_loss~relevant_history*X, pca.data.swap) 
summary(wloss_history_xX_mod) # no significance (too many groups?)
Anova(wloss_history_xX_mod) # confirmed, without significant residuals
sigma(wloss_history_xX_mod)*100/mean(pca.data.swap$maximum_weight_loss)# percentage error = 6.060462

# weightloss and relevant history with X as interaction (not necessary)
wloss_history_xX_xY_mod <- lm(maximum_weight_loss~relevant_history * X * Y, pca.data.swap) 
summary(wloss_history_xX_xY_mod) # no significance (too many groups?)
Anova(wloss_history_xX_xY_mod) # confirmed, without significant residuals
sigma(wloss_history_xX_xY_mod)*100/mean(pca.data.swap$maximum_weight_loss)# percentage error = 5.420395


# weightloss and relevant history with pca X
wloss_history_X_mod <- lm(maximum_weight_loss~relevant_history + X, pca.data.swap)
summary(wloss_history_X_mod) # some strains have significant effect
Anova(wloss_history_X_mod) # X is only borderline significant with house strains
sigma(wloss_history_X_mod)*100/mean(pca.data.swap$maximum_weight_loss)# percentage error = 6.046203

# weightloss and relevant history with mouse strain and X
wloss_history_mouse_X_mod <- lm(maximum_weight_loss~relevant_history + mouse_strain + X, pca.data.swap)
summary(wloss_history_mouse_X_mod) # strain effects remain more or less + strong infection effect
Anova(wloss_history_mouse_X_mod) # no residual significance, history has strongest effect
sigma(wloss_history_mouse_X_mod)*100/mean(pca.data.swap$maximum_weight_loss)# percentage error = 5.451805


# compare weightloss~history (anova tests in order specified so separate test for now)
# to w~h to w~h + X
anova(wloss_history_mod, wloss_history_X_mod) # X better of course
# w~h to w~h + mouse
anova(wloss_history_mod, wloss_history_mouse_mod) # mouse better of course
# to w~h + mouse + X
anova(wloss_history_mod, wloss_history_mouse_X_mod) # mouse + X better of course
#  w~h+mouse to w~h+mouse+X
anova(wloss_history_mouse_mod, wloss_history_mouse_X_mod)
#  w~h+X to w~h+mouse+X + load mouse only behind combined model to compare (because only mouse and only X are incomparable)
anova(wloss_history_X_mod, wloss_history_mouse_X_mod, wloss_history_mouse_mod)
anova(wloss_history_mouse_mod, wloss_history_mouse_X_mod, wloss_history_X_mod)

# mouse better than X but combined best


###################### make X response variable to strains (between mouse strain variance 
pca.data.swap$mouse_strain <- relevel(pca.data.swap$mouse_strain, "SCHUNT_SCHUNT")
library(RColorBrewer)
nb.cols <- 13
mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(nb.cols)
# correlate weightloss and X (immune)
shapiro.test(pca.data.swap$X)
shapiro.test(pca.data.swap$maximum_weight_loss)
cor(pca.data.swap$X, pca.data.swap$maximum_weight_loss)

#this one only on infected
X_mouse_mod <- lm(X~mouse_strain, pca.data.swap)
summary(X_mouse_mod) # strain effects remain more or less + strong infection effect
Anova(X_mouse_mod) # no residual significance, history has strongest effect
sigma(X_mouse_mod)*100/mean(pca.data.swap$maximum_weight_loss)# percentage error = 2.275198
X_mouse_modgg <- ggpredict(X_mouse_mod, c("mouse_strain"))
plot(X_mouse_modgg) + 
  geom_point() +
  xlab(label = "mouse_strain") +
  ylab(label = bquote(atop("immune response intensity"))) +
  theme_bw() +
  ggtitle("Predicted immune response intensity by mouse strain")

# add weightloss
X_mouse_weight_mod <- lm(X~mouse_strain + maximum_weight_loss, pca.data.swap)
summary(X_mouse_weight_mod) # strain effects remain more or less + strong infection effect
Anova(X_mouse_weight_mod) # no residual significance, history has strongest effect
sigma(X_mouse_weight_mod)*100/mean(pca.data.swap$maximum_weight_loss)# percentage error = 2.287074
X_mouse_weight_modgg <- ggpredict(X_mouse_weight_mod, c("maximum_weight_loss", "mouse_strain"))
plot(X_mouse_weight_modgg) + 
  scale_fill_manual(values = mycolors) +
  geom_point() +
  xlab(label = "mouse_strain") +
  ylab(label = bquote(atop("immune response intensity"))) +
  theme_bw() +
  ggtitle("Predicted immune response intensity by mouse strain")

X_history_weight_mod <- lm(X~relevant_history + maximum_weight_loss, pca.data.swap)
summary(X_history_weight_mod) # strain effects remain more or less + strong infection effect
Anova(X_history_weight_mod) # no residual significance, history has strongest effect
sigma(X_history_weight_mod)*100/mean(pca.data.swap$maximum_weight_loss)# percentage error = 2.287074
X_history_weight_modgg <- ggpredict(X_history_weight_mod, c("relevant_history", "maximum_weight_loss"))
plot(X_history_weight_modgg) + 
  geom_point() +
  xlab(label = "mouse_strain") +
  ylab(label = bquote(atop("immune response intensity"))) +
  theme_bw() +
  ggtitle("Predicted immune response intensity by mouse strain")

# winner?
X_mouse_history_mod <- lm(X~mouse_strain + relevant_history, pca.data.swap)
summary(X_mouse_history_mod)
Anova(X_mouse_history_mod) # no residual significance, history has strongest effect
sigma(X_mouse_history_mod)*100/mean(pca.data.swap$maximum_weight_loss)# percentage error = 1.85589

X_mouse_history_modgg <- ggpredict(X_mouse_history_mod, c("mouse_strain","relevant_history"))
plot(X_mouse_history_modgg) + 
  geom_point() +
  xlab(label = "mouse strain") +
  ylab(label = bquote(atop("change in immune cell populations (%)"))) +
  theme_bw() +
  ggtitle("")



# most complex for immune
X_weight_mouse_history_mod <- lm(X~maximum_weight_loss + mouse_strain + relevant_history, pca.data.swap)
summary(X_weight_mouse_history_mod) # strain effects remain more or less + strong infection effect
Anova(X_weight_mouse_history_mod) # no residual significance, history has strongest effect
sigma(X_weight_mouse_history_mod)*100/mean(pca.data.swap$maximum_weight_loss)# percentage error = 1.774119

library(RColorBrewer)
# Define the number of colors you want
nb.cols <- 18
mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(nb.cols)
# Create a ggplot with 18 colors 
# Use scale_fill_manual
ggplot(df) + 
  geom_col(aes(name, Sepal.Length, fill = factor(Sepal.Length))) +
  scale_fill_manual(values = mycolors) +
  theme_minimal() +
  theme(legend.position = "top")


X_weight_mouse_history_modgg <- ggpredict(X_weight_mouse_history_mod, c("maximum_weight_loss","mouse_strain", "relevant_history"))
plot(X_weight_mouse_history_modgg) + 
  geom_point() +
  scale_fill_manual(values = mycolors) +
  xlab(label = "weight loss") +
  ylab(label = bquote(atop("change in immune cell populations (%)"))) +
  theme_bw() +
  ggtitle("")


# anova this
anova(X_mouse_history_mod, X_mouse_weight_history_mod,X_mouse_weight_mod)

testDispersion(X_mouse_mod)
testDispersion(X_mouse_weight_history_mod)
testDispersion(X_mouse_weight_mod)

X_mouse_mod_simulated <- simulateResiduals(fittedModel = X_mouse_mod, plot = T)
X_mouse_weight_history_mod_simulated <- simulateResiduals(fittedModel = X_mouse_weight_history_mod, plot = T)
X_mouse_weight_mod_simulated <- simulateResiduals(fittedModel = X_mouse_weight_mod, plot = T)

plot(wloss_history_X_mod_simulated) # comes back with quantile deviations
plotResiduals(wloss_history_X_mod_simulated, form = pca.data.swap$X) # significant on weightloss
plot(wloss_history_mouse_X_mod_simulated)
plotResiduals(wloss_history_mouse_X_mod_simulated, form = pca.data.swap$X) # also significant but not in plot()
plotResiduals(wloss_history_mouse_X_mod_simulated, form = pca.data.swap$mouse_strain)
plot(wloss_history_mouse_mod_simulated) # comes back with quantile deviations
plotResiduals(wloss_history_mouse_mod_simulated, form = pca.data.swap$mouse_strain)
# or within mouse strain variance in X ( what drives wloss more))

# check for worst weightloss mouse strain, rearrange intercept and look again
#######################################




# wloss_history_mouse_X_mod has all the factors we are interested in + has the lowest percentage of error on weight
# compare with wloss_history_X_mod, wloss_history_mouse_mod

# graph wloss_history_X_mod
wloss_history_X_modgg <- ggpredict(wloss_history_X_mod, c("relevant_history", "X"))
plot(wloss_history_X_modgg) + 
  geom_point() +
  xlab(label = "infection history") +
  ylab(label = bquote(atop("maximum weight loss (%)",
                   "lower = higher impact"))) +
  theme_bw() +
  ggtitle("Predicted weightloss by infection history and principal immune parameters")

wloss_history_modgg <- ggpredict(wloss_history_mod, c("relevant_history"))
plot(wloss_history_modgg) + 
  geom_point() +
  xlab(label = "infection history") +
  ylab(label = bquote(atop("maximum weight loss (%)",
                           "lower = higher impact"))) +
  theme_bw() +
  ggtitle("Predicted weightloss by infection history")




wloss_X <- lm(maximum_weight_loss~X, pca.data.swap)
summary(wloss_X)
wloss_Xgg <- ggpredict(wloss_X, c("X"))
plot(wloss_Xgg) + 
  geom_jitter() +
  xlab(label = "mouse strain") +
  ylab(label = bquote(atop("change in immune cell populations (%)"))) +
  theme_bw() +
  ggtitle("")

ggplot(pca.data.swap, aes(maximum_weight_loss, X)) +
  geom_point()



summary(X_mouse_mod)
X_mouse_modgg <- ggpredict(X_mouse_mod, c("mouse_strain"))
plot(X_mouse_modgg) + 
  geom_point() +
  xlab(label = "mouse strain") +
  ylab(label = bquote(atop("change in immune cell populations (%)"))) +
  theme_bw() +
  ggtitle("differences in immune cell populations between strains")

X_history_mod <- lm(X~relevant_history, pca.data.swap)
summary(X_history_mod)
X_history_modgg <- ggpredict(X_history_mod, c("relevant_history"))
plot(X_history_modgg) + 
  geom_point() +
  xlab(label = "mouse strain") +
  ylab(label = bquote(atop("change in immune cell populations (%)"))) +
  theme_bw() +
  ggtitle("differences in immune cell populations between strains")



wloss_X_modgg <- ggpredict(wloss_X_mod, c("X"))
plot(wloss_X_modgg) + 
  geom_point() +
  xlab(label = "immune cells") +
  ylab(label = bquote(atop("weight loss (%)"))) +
  theme_bw() +
  ggtitle("")








# graph wloss_history_mouse_mod



wloss_history_mouse_modgg <- ggpredict(wloss_history_mouse_mod, c("mouse_strain", "relevant_history"))
plot(wloss_history_mouse_modgg) + 
  geom_point() +
  scale_fill_manual(values = mycolors) +
  xlab(label = "mouse strain") +
  ylab(label = bquote(atop("maximum weight loss (%)",
                           "lower = higher impact"))) +
  theme_bw() +
  ggtitle("Predicted weightloss by infection history and mouse strain")

# graph  wloss_history_mouse_X_mod
wloss_history_mouse_X_modgg <- ggpredict(wloss_history_mouse_X_mod, c("relevant_history", "mouse_strain", "X"))
plot(wloss_history_mouse_X_modgg) + 
  geom_point() +
  scale_fill_manual(values = mycolors) +
  xlab(label = "infection history") +
  ylab(label = bquote(atop("maximum weight loss (%)",
                           "lower = higher impact"))) +
  theme_bw() +
  ggtitle("Predicted weightloss by infection history and principal immune parameters")


# compare more than visualy, AIC tab for the relevant models
models <- list(wloss_history_X_mod, wloss_history_mouse_X_mod, wloss_history_mouse_mod)
model.names <- c("wloss_history_X_mod", "wloss_history_mouse_X_mod", "wloss_history_mouse_mod")
aictab(cand.set = models, modnames = model.names, sort = T)
# add likelyhood tests

# check with DHARMa for wloss_history_X_mod, wloss_history_mouse_X_mod, wloss_history_mouse_mod
library(DHARMa)
testDispersion(wloss_history_X_mod)
testDispersion(wloss_history_mouse_X_mod)
testDispersion(wloss_history_mouse_mod)

wloss_history_X_mod_simulated <- simulateResiduals(fittedModel = wloss_history_X_mod, plot = T)
wloss_history_mouse_X_mod_simulated <- simulateResiduals(fittedModel = wloss_history_mouse_X_mod, plot = T)
wloss_history_mouse_mod_simulated <- simulateResiduals(fittedModel = wloss_history_mouse_mod, plot = T)

plot(wloss_history_X_mod_simulated) # comes back with quantile deviations
plotResiduals(wloss_history_X_mod_simulated, form = pca.data.swap$X) # significant on weightloss
plot(wloss_history_mouse_X_mod_simulated)
plotResiduals(wloss_history_mouse_X_mod_simulated, form = pca.data.swap$X) # also significant but not in plot()
plotResiduals(wloss_history_mouse_X_mod_simulated, form = pca.data.swap$mouse_strain)
plot(wloss_history_mouse_mod_simulated) # comes back with quantile deviations
plotResiduals(wloss_history_mouse_mod_simulated, form = pca.data.swap$mouse_strain) # significant on weightloss

#######################################################################################################################
################### Can all weightloss be explained by the PCA loading components?
# Scatter plot: Visualize the linear relationship between the predictor and response
scatter.smooth(x=pca.data.swap$maximum_weight_loss, y=pca.data.swap$X, main="Weight loss ~ PCA1")
# linearly increasing but messy
# boxplot for outliers
par(mfrow=c(1, 2))  # divide graph area in 2 columns
boxplot(pca.data.swap$maximum_weight_loss, main="weight loss", sub=paste("Outlier rows: ", 
        boxplot.stats(pca.data.swap$maximum_weight_loss)$out))  # box plot for 'speed'
boxplot(pca.data.swap$X, main="PCA1", sub=paste("Outlier rows: ", 
        boxplot.stats(pca.data.swap$X)$out))  # box plot for 'distance'

library(e1071)
par(mfrow=c(1, 2))  # divide graph area in 2 columns
plot(density(pca.data.swap$maximum_weight_loss), main="Density Plot: weight loss",
     ylab="Frequency", sub=paste("Skewness:", round(e1071::skewness(pca.data.swap$maximum_weight_loss), 2)))
polygon(density(pca.data.swap$maximum_weight_loss), col="red")

plot(density(pca.data.swap$X), main="Density Plot: PCA1", 
     ylab="Frequency", sub=paste("Skewness:", round(e1071::skewness(pca.data.swap$X), 2)))  # density plot for 'dist'
polygon(density(pca.data.swap$X), col="red")
cor(pca.data.swap$maximum_weight_loss, pca.data.swap$X)
# quite weak correlation at 0.2380415
##############
wloss_X_mod <- lm(maximum_weight_loss ~ X, pca.data.swap)
summary(wloss_X_mod) # we have a linear relationship
# test the assumptions for linear model:

#### 1) check residuals
mean(wloss_X_mod$residuals) #  1.488036e-16

### 2) Homoscedasticity of residuals or equal variance
par(mfrow=c(2,2))  # set 2 rows and 2 column plot layout
plot(wloss_X_mod) # equal variance, check on subset
wloss_X_mod20 <- lm(maximum_weight_loss ~ X, pca.data.swap[1:20, ])  #  linear model
plot(wloss_X_mod20) # still holds

#### 3) autocorrelation
# Durbin-Watson test
lmtest::dwtest(wloss_X_mod)
# data:  wloss_X_mod
# DW = 1.2149, p-value = 6.761e-05
# alternative hypothesis: true autocorrelation is greater than 0
# How to rectify?
library(DataCombine)
pca.data.swap1 <- data.frame(pca.data.swap, resid_mod1=wloss_X_mod$residuals)
pca.data.swap1_1<- slide(pca.data.swap1, Var="resid_mod1", NewVar = "lag1", slideBy = -1)
pca.data.swap1_2 <- na.omit(pca.data.swap1_1)
wloss_X_mod2 <- lm(maximum_weight_loss ~ X + lag1, data=pca.data.swap1_2)
acf(wloss_X_mod2$residuals)
library(snpar)
runs.test(wloss_X_mod2$residuals)  # runs test
# Approximate runs rest
# 
# data:  wloss_X_mod2$residuals
# Runs = 45, p-value = 0.5798
# alternative hypothesis: two.sided
lmtest::dwtest(wloss_X_mod2)
# Durbin-Watson test
# 
# data:  wloss_X_mod2
# DW = 2.1992, p-value = 0.8009
# alternative hypothesis: true autocorrelation is greater than 0
# passed above 2 so the residuals are not autocorrelated

### 4) The X variables and residuals are uncorrelated
cor.test(pca.data.swap$X, wloss_X_mod$residuals)
# Pearson's product-moment correlation
# 
# data:  pca.data.swap$X and wloss_X_mod$residuals
# t = 2.0589e-16, df = 83, p-value = 1
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.2131241  0.2131241
# sample estimates:
#          cor 
# 2.259943e-17
# X variables and residuals are uncorrelated

### 5) The number of observations must be greater than number of Xs = YES

### 6) The variability in X values is positive
var(pca.data.swap$X) 
# 5.434293 satisfied

### 7) No perfect multicollinearity
vif(wloss_X_mod) # VIF = 1, strict cutoff is below 2 so satisifed.

### 8) Normality of residuals
library(gvlma)
gvlma(wloss_X_mod)
# Call:
#   lm(formula = maximum_weight_loss ~ X, data = pca.data.swap)
# 
# Coefficients:
#   (Intercept)            X  
# 92.2108       0.6729  
# 
# 
# ASSESSMENT OF THE LINEAR MODEL ASSUMPTIONS
# USING THE GLOBAL TEST ON 4 DEGREES-OF-FREEDOM:
#   Level of Significance =  0.05 
# 
# Call:
#   gvlma(x = wloss_X_mod) 
# 
# Value  p-value                   Decision
# Global Stat        10.51951 0.032529 Assumptions NOT satisfied!
#   Skewness            7.99686 0.004686 Assumptions NOT satisfied!
#   Kurtosis            0.02621 0.871377    Assumptions acceptable.
# Link Function       0.02597 0.871977    Assumptions acceptable.
# Heteroscedasticity  2.47047 0.116003    Assumptions acceptable.

wloss_X_mod_influence <- influence.measures(wloss_X_mod)

gvlma(wloss_history_mouse_X_mod)
gvlma(wloss_history_X_mod)
gvlma(wloss_history_mod)
gvlma(wloss_X_mod)
gvlma(X_mouse_mod)

bla <- lm(maximum_weight_loss ~ X, pca.data.swap)
summary(bla)
gvlma(bla)

################ training the model
# # Residuals:
# Min       1Q   Median       3Q      Max 
# -17.5444  -3.8127   0.7998   5.0119   9.6131 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  92.2108     0.6984 132.032   <2e-16 ***
#   X             0.6729     0.3014   2.233   0.0283 *  
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 6.439 on 83 degrees of freedom
# Multiple R-squared:  0.05666,	Adjusted R-squared:  0.0453 
# F-statistic: 4.986 on 1 and 83 DF,  p-value: 0.02825
AIC(wloss_X_mod)  # AIC => 561.7967
BIC(wloss_X_mod)  # BIC => 569.1247
# Step 1: Create the training (development) and test (validation) data samples from original data.
set.seed(100)  # setting seed to reproduce results of random sampling
trainingRowIndex <- sample(1:nrow(pca.data.swap), 0.8*nrow(pca.data.swap))  # row indices for training data
trainingData <- pca.data.swap[trainingRowIndex, ]  # model training data
testData  <- pca.data.swap[-trainingRowIndex, ]   # test data
# Step 2: Develop the model on the training data and use it to predict the distance on test data
# Build the model on training data -
lmMod <- lm(maximum_weight_loss ~ X, data=trainingData)  # build the model
distPred <- predict(lmMod, testData)  # predict distance
# Step 3: Review diagnostic measures.
summary(lmMod) # compare to full data model
# Call:
#   lm(formula = maximum_weight_loss ~ X, data = trainingData)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -17.2541  -2.9820   0.5073   4.7311   9.5918 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  92.1011     0.7784 118.314   <2e-16 ***
#   X             0.7725     0.3616   2.136   0.0364 *  
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 6.316 on 66 degrees of freedom
# Multiple R-squared:  0.06468,	Adjusted R-squared:  0.0505 
# F-statistic: 4.564 on 1 and 66 DF,  p-value: 0.03637
summary(wloss_X_mod)
# Call:
# lm(formula = maximum_weight_loss ~ X, data = pca.data.swap)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -17.5444  -3.8127   0.7998   5.0119   9.6131 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  92.2108     0.6984 132.032   <2e-16 ***
#   X             0.6729     0.3014   2.233   0.0283 *  
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 6.439 on 83 degrees of freedom
# Multiple R-squared:  0.05666,	Adjusted R-squared:  0.0453 
# F-statistic: 4.986 on 1 and 83 DF,  p-value: 0.02825

# Step 4: Calculate prediction accuracy and error rates
## higher correlation accuracy implies that the actuals and predicted values have similar directional movement
actuals_preds <- data.frame(cbind(actuals=testData$maximum_weight_loss, predicteds=distPred))  # make actuals_predicteds dataframe.
correlation_accuracy <- cor(actuals_preds)  
head(actuals_preds)
# MeanAbsolutePercentageError (MAPE)
min_max_accuracy <- mean(apply(actuals_preds, 1, min) / apply(actuals_preds, 1, max))  
# => 93.96%, min_max accuracy
mape <- mean(abs((actuals_preds$predicteds - actuals_preds$actuals))/actuals_preds$actuals)  
# => 6.39%, mean absolute percentage deviation
library(DAAG)
cvResults <- suppressWarnings(CVlm(data = pca.data.swap, form.lm=maximum_weight_loss ~ X, m=5,
                                   dots=FALSE, seed=29, legend.pos="topleft",  printit=FALSE, 
                                   main="Small symbols are predicted values while bigger ones are actuals."));  # performs the CV
attr(cvResults, 'ms')  # => 251.2783 mean squared error


# code infection history differently (only present infection) 
# immunized or not + acute strain

# continue challenge analysis with secondary only
# model 5 exploration

# make pca.challenge (only current infection weightloss as compared to maximum weightloss)


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

