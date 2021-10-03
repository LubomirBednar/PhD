# WDS
library(ggplot2)
library(dplyr)
library(tidyverse)
library(magrittr)
library(ggpubr)
library(ggeffects)
library(effects)
library(car)
library(AICcmodavg)
library(RColorBrewer)
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

E11q <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/E11_012021_Eim_CEWE_qPCR.csv")
E11q <- select(E11q, EH_ID, delta, Eim_MC)
E11q$dpi <- 8
E11q$primary_infection <- NA

E11 <- merge(E11, E11q, all.x = T)


# CLS
CLS <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/CLS_complete.csv")
CLS$X <- NULL 
CLS$EH_ID <- as.character(sub("_", "", CLS$EH_ID))
names(CLS)[names(CLS) == "Strain"] <- "mouse_strain"
names(CLS)[names(CLS) == "primary"] <- "primary_infection"
names(CLS)[names(CLS) == "challenge"] <- "challenge_infection"
CLSimmuno <- select(CLS, EH_ID, experiment, dpi, labels, delta, IFNy_CEWE, Eim_MC, EXP_type, CD4, Treg, Div_Treg, Position,
                    Treg17, Th1, Div_Th1, Th17, Div_Th17, CD8, Act_CD8, Div_Act_CD8, IFNy_CD4, IFNy_CD8, batch)
CLSimmuno <- subset(CLSimmuno, CLSimmuno$dpi == 8)
CLSimmuno <- subset(CLSimmuno, CLSimmuno$batch == "b")

CLS_primary <- subset(CLS, CLS$batch == "a")
CLS_primary$challenge_infection <- NA

CLS_challenge <- subset(CLS, CLS$batch == "b")
CLS_challenge$primary_infection <- NA

CLS <- rbind(CLS_primary, CLS_challenge)

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
WDS_primary$OPG <- as.numeric(WDS_primary$OPG)

# shuld have a dataset with all weight data and some OPG data

ggplot(WDS_primary, aes(x = dpi, y = relative_weight, color = eimeria_species_primary)) +
  geom_jitter(width = 0.2) +
  geom_smooth() +
  labs(y = "weight relative to dpi 0", x = "days post infection", color = "infecting species") +
  theme(axis.text=element_text(size=12, face = "bold"),
        title = element_text(size = 16, face = "bold"),
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        strip.text.y = element_text(size = 14, face = "bold"),
        legend.text = element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold")) +
   ggtitle("Weight loss during primary infections")
  
ggplot(WDS_primary, aes(x = dpi, y = OPG, color = eimeria_species_primary)) +
  geom_jitter(width = 0.2) +
  geom_smooth() +
  labs(y = "oocysts per gram of feces", x = "days post infection", color = "infecting species") +
  theme(axis.text=element_text(size=12, face = "bold"),
        title = element_text(size = 16, face = "bold"),
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        strip.text.y = element_text(size = 14, face = "bold"),
        legend.text = element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold")) +
  ggtitle("Oocyst shedding during primary infection")




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

############ not finished IFNy CEWE
# E11inf <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/E10_E11_Eim_CEWE_ELISA.csv")
# E10inf <- 

E11immuno$delta <- "NA"
E11immuno$IFNy_CEWE <- "NA"
E11immuno$Eim_MC <- "NA"
E11immuno$EXP_type <- "WDS"

E_immuno <- rbind(E7immuno, E11immuno)
E_immuno <- rbind(E_immuno, CLSimmuno)


#################### add and process IFC runs
IFC1 <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/IFC1.csv")
IFC2 <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/IFC2.csv")
IFC3 <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/IFC3.csv")
IFC4 <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/IFC4.csv")
IFC5 <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/IFC5.csv")

IFC <- bind_rows(IFC1, IFC2, IFC3, IFC4, IFC5)
IFC$EH_ID <- as.character(sub("A0729", "AA0729", IFC$EH_ID))
# remove unsuccessful amplifications
IFC <- subset(IFC, IFC$Value != 999)

IFC <- IFC %>% dplyr::group_by(EH_ID, Target) %>% dplyr::summarise(Ct = mean(Value, na.rm = T))

IFC <- distinct(IFC)

IFC.wide <- pivot_wider(IFC, names_from = Target, values_from = Ct)
IFC.wide <- data.frame(IFC.wide)
# write.csv(IFC.wide, "~/GitHub/PhD/IFC.wide.csv")

# IFC.wide <- read.csv("C:/Users/exemp/OneDrive/Documents/GitHub/PhD/IFC.wide.csv")
# IFC.wide$X <- NULL

IFC.wide[2:ncol(IFC.wide)] <- IFC.wide[2:ncol(IFC.wide)]-IFC.wide[,18]

# refGenes <- c("PPIB", "GAPDH")
# targetGenes <- c("CASP1","CXCL9","CXCR3","IDO1", "IFNG", "IL10", "IL12A","IL13",
#                  "IL1RN","IRGM1","MPO","MUC2","MUC5AC","MYD88","NCR1", "PRF1", "RETNLB", 
#                  "SOCS1","TICAM1","TNF","IL6","IL17A")
# 
# IFC.wide <- IFC.wide %>% mutate(refMean = rowMeans(select(., refGenes)))
# IFC.wide <- data.frame(IFC.wide)
# refMean <- as.numeric(IFC.wide$refMean)

IFC.wide$PPIB <- NULL
IFC.wide$GAPDH <- NULL
IFC.wide$dpi <- 8

# merge into WDS
WDS_challenge <- merge(WDS_challenge, E_immuno, all.x = T)
WDS_challenge <- merge(WDS_challenge, IFC.wide, all.x = T)



############################################################################################
IFC.long <- pivot_longer(IFC.wide, cols = c("CASP1","CXCL9","CXCR3","IDO1", "IFNG", "IL10", "IL12A","IL13",
                                            "IL1RN","IRGM1","MPO","MUC2","MUC5AC","MYD88","NCR1", "PRF1", "RETNLB", 
                                            "SOCS1","TICAM1","TNF","IL6","IL17A"))

names(IFC.long)[names(IFC.long) == "name"] <- "Target"
names(IFC.long)[names(IFC.long) == "value"] <- "Ct"

ggimmuno <- pivot_longer(E_immuno, cols = c("Div_Treg", "Treg", "CD4", "Treg17", "Th1", "Div_Th1", "Th17", "Div_Th17",
                                          "CD8", "Act_CD8", "Div_Act_CD8", "IFNy_CD4", "IFNy_CD8"))


ggimmuno <- merge(ggimmuno, IFC.long)

ggplot(ggimmuno, aes(x = Target, y = Ct, color = Target)) +
  geom_boxplot() +
  facet_wrap(~Eim_MC)

################ explore PCA
WDS_PCA <- subset(WDS_challenge, WDS_challenge$dpi == 8)
WDS_PCA <- subset(WDS_PCA, WDS_PCA$Position == "mLN")
WDS_PCA <- select(WDS_PCA, EH_ID, #challenge_infection, mouse_strain, infection_history, 
                CD4, Treg, Div_Treg, Treg17, Th1, Div_Th1, Th17, Div_Th17, CD8, Act_CD8, Div_Act_CD8, IFNy_CD4, IFNy_CD8, 
                # CASP1, CXCL9, CXCR3, IDO1,  IFNG,  IL10,  IL12A, IL13, IL1RN, IRGM1, MPO, MUC2, MUC5AC, MYD88, NCR1, PRF1,  RETNLB,  
                # SOCS1, TICAM1, TNF, IL6, IL17A)
)

# WDS_PCA <- na.omit(WDS_PCA)
WDS_PCA <- distinct(WDS_PCA)
WDS_PCA <- WDS_PCA %>% remove_rownames %>% column_to_rownames(var="EH_ID")
pca <- prcomp(na.omit(WDS_PCA), center = TRUE, scale = TRUE)
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
                eimeria_species_challenge, eimeria_species_history, Eim_MC,
                relative_weight, CD4, Treg, Div_Treg, Treg17, Th1, Div_Th1, Th17, Div_Th17, CD8, Act_CD8, Div_Act_CD8, 
                IFNy_CD4, IFNy_CD8#, CASP1, CXCL9, CXCR3, IDO1,  IFNG,  IL10, IL12A, IL13, IL1RN, IRGM1, MPO, MUC2, 
                #MUC5AC, MYD88, NCR1, PRF1,  RETNLB, SOCS1, TICAM1, TNF, IL6, IL17A)
)
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

ggplot(data = pca.data, aes(x = X, y = Y, label = EH_ID, color = Eim_MC)) +
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
labs(color='infection') 




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
pca.data$secondary_species <- factor(pca.data$secondary_species,
                                     levels =  c("UNI", "FER", "FAL"))

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

#################### same for oocysts
max_oocyst_primary <- WDS_primary %>% select(EH_ID, OPG, primary_infection, 
                                             challenge_infection, infection_history) 
max_oocyst_challenge <- WDS_challenge %>% select(EH_ID, OPG, primary_infection, 
                                                 challenge_infection, infection_history)
# except that this uses only data from challenge for primary:uni
max_oocyst_primary <- max_oocyst_primary[!is.na(max_oocyst_primary$OPG), ]
max_oocyst_challenge <- max_oocyst_challenge[!is.na(max_oocyst_challenge$OPG),]
# creates weird 1 sample with 0 so remove it
# max_weight_challenge <- max_weight_challenge[ !(max_weight_challenge$EH_ID %in% "LM0350"),]
# max_OPG<- WDS %>% select(EH_ID, OPG)

# max_OPG$OPG <- as.numeric(as.character(max_OPG$OPG))
max_oocyst_primary <- max_oocyst_primary %>% dplyr::group_by(EH_ID) %>% dplyr::summarise(OPG = max(OPG, na.rm = T))
max_oocyst_challenge <- max_oocyst_challenge %>% dplyr::group_by(EH_ID) %>% dplyr::summarise(OPG = max(OPG, na.rm = T))
# max_OPG <- max_OPG %>% dplyr::group_by(EH_ID) %>% dplyr::summarise(OPG = max(OPG, na.rm = T))
names(max_oocyst_primary)[names(max_oocyst_primary) == "OPG"] <- "maximum_OPG_primary"
names(max_oocyst_challenge)[names(max_oocyst_challenge) == "OPG"] <- "maximum_OPG_challenge"


pca.data1 <- merge(pca.data, max_weight_primary, by = "EH_ID")
pca.data1 <- merge(pca.data1, max_weight_challenge, by = "EH_ID")
pca.data1 <- merge(pca.data1, max_oocyst_primary, by = "EH_ID")
pca.data1 <- merge(pca.data1, max_oocyst_challenge, by = "EH_ID")

# anything that is infected:uninfected needs to be the in the final maximum weightloss column from primary_max_WLoss
# can do min across both columns but that didnt work because its non disciminating to our criteria so scrapping
# pca.challenge <- transform(pca.challenge, min = pmin(maximum_weight_loss_primary, maximum_weight_loss_challenge))
############## THIS WILL ONLY WORK FOR infection_history MODELS! (pca.data.swap1)
pca.data.swap1 <- subset(pca.data1, infection_history == "E64:UNI" | 
                        infection_history == "E88:UNI"|
                        infection_history == "Eflab:UNI")
# drop challenge weightloss column to avoic UNI weight in infection history

# name the columns the same and merge
pca.data.swap2 <- subset(pca.data1, infection_history == "UNI:UNI" | 
                    infection_history == "E88:E88"|
                    infection_history == "E88:E64"|
                    infection_history == "Eflab:E88"|
                    infection_history == "E64:E64"|
                    infection_history == "E64:E88"|
                    infection_history == "UNI:E88"|
                    infection_history == "UNI:E64")


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

# Define the number of colors you want
nb.cols <- 18
mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(nb.cols)
#####################################################################################################
#####################compare more than visualy, AIC tab for the relevant models######################
# EXPLORATION so full dataset







############### weight loss explanation
# use dataset without SWISS mice (validation and training)
# pca.data.swapval <- subset(pca.data.swap, pca.data.swap$mouse_strain == "SWISS")
# pca.data.swap <- subset(pca.data.swap, !pca.data.swap$mouse_strain == "SWISS")

# now build these to reflext wild better (only current infections and hybrid, mmd and mmm mouse strains)
# already got secondary_species so let's code the mouse strains to simpler groups
pca.data.swap$simple_strain[pca.data.swap$mouse_strain == "SWISS"] <- "SWISS"
pca.data.swap$simple_strain[pca.data.swap$mouse_strain == "SCHUNT_SCHUNT"] <- "Mmd"
pca.data.swap$simple_strain[pca.data.swap$mouse_strain == "PWD_PWD"] <- "Mmm"
pca.data.swap$simple_strain[pca.data.swap$mouse_strain == "BUSNA_STRA"] <- "Hybrid"
pca.data.swap$simple_strain[pca.data.swap$mouse_strain == "STRA_BUSNA"] <- "Hybrid"
pca.data.swap$simple_strain[pca.data.swap$mouse_strain == "PWD_SCHUNT"] <- "Hybrid"
pca.data.swap$simple_strain[pca.data.swap$mouse_strain == "STRA_STRA"] <- "Mmd"
pca.data.swap$simple_strain[pca.data.swap$mouse_strain == "STRA_SCHUNT"] <- "Mmd"
pca.data.swap$simple_strain[pca.data.swap$mouse_strain == "PWD_BUSNA"] <- "Mmm"
pca.data.swap$simple_strain[pca.data.swap$mouse_strain == "SCHUNT_PWD"] <- "Hybrid"
pca.data.swap$simple_strain[pca.data.swap$mouse_strain == "SCHUNT_STRA"] <- "Mmd"
pca.data.swap$simple_strain[pca.data.swap$mouse_strain == "BUSNA_BUSNA"] <- "Mmm"
pca.data.swap$simple_strain[pca.data.swap$mouse_strain == "BUSNA_PWD"] <- "Mmm"

write.csv(pca.data.swap, "C:/Users/exemp/OneDrive/Documents/pca.data.swap.csv")



































# ggplot(pca.data.swap, aes(x = dpi, y = relative_weight, color = eimeria_species_primary)) +
#   geom_jitter(width = 0.2) +
#   geom_smooth() +
#   labs(y = "weight relative to dpi 0", x = "days post infection", color = "infecting species") +
#   theme(axis.text=element_text(size=12, face = "bold"),
#         title = element_text(size = 16, face = "bold"),
#         axis.title=element_text(size=14,face="bold"),
#         strip.text.x = element_text(size = 14, face = "bold"),
#         strip.text.y = element_text(size = 14, face = "bold"),
#         legend.text = element_text(size=12, face = "bold"),
#         legend.title = element_text(size = 12, face = "bold")) +
#   ggtitle("infection intensities in wild and wild-derived mice")
# order factors first
pca.data.swap$simple_strain <- factor(pca.data.swap$simple_strain,
                                      levels = c("SWISS", "Mmm", "Mmd", "Hybrid"))
pca.data.swap$secondary_species <- factor(pca.data.swap$secondary_species,
                                          levels = c("UNI", "FER", "FAL"))

# make function for testing all the potential effects
#interests <- names(select(pca.data.swap, X, Y, challenge_infection, mouse_strain, infection_history, 
                          # eimeria_species_challenge, eimeria_species_history, relative_weight, CD4, Treg, 
                          # Div_Treg, Treg17, Th1, Div_Th1, Th17, Div_Th17,CD8, Act_CD8, Div_Act_CD8, IFNy_CD4,
                          # IFNy_CD8, secondary_species, relevant_history, simple_strain))

pca.data.swap$maximum_OPG_challenge <- as.numeric(pca.data.swap$maximum_OPG_challenge)
wloss_X_current_simpstrain_OPG_mod <- lm(maximum_weight_loss_challenge ~ X + secondary_species + simple_strain + 
                                         maximum_OPG_challenge, pca.data.swap)
summary(wloss_X_current_simpstrain_OPG_mod)

OPG_X_current_mouse_mod <- lm(maximum_OPG_challenge ~ X + secondary_species + mouse_strain, pca.data.swap)
summary(OPG_X_current_mouse_mod)

wloss_X_mod <- lm(maximum_weight_loss_challenge ~ X, pca.data.swap)
summary(wloss_X_mod)
wloss_X_history_mod <- lm(maximum_weight_loss_challenge ~ X + eimeria_species_history, pca.data.swap)
summary(wloss_X_history_mod)
wloss_X_history_mouse_mod <- lm(maximum_weight_loss_challenge ~ X +eimeria_species_history + mouse_strain, pca.data.swap)
summary(wloss_X_history_mouse_mod)
anova(wloss_X_mod, wloss_X_history_mod, wloss_X_history_mouse_mod)

wloss_X_current_mod <- lm(maximum_weight_loss_challenge ~ X + secondary_species, pca.data.swap)
summary(wloss_X_current_mod)
wloss_X_current_simpstrain_mod <- lm(maximum_weight_loss_challenge ~ X + secondary_species + simple_strain, pca.data.swap)
summary(wloss_X_current_simpstrain_mod)

# include validation dataset
# TEMPORARY remove UNI:FAL as it is unique for this dataset
pca.data.swapval <- subset(pca.data.swap, !pca.data.swap$eimeria_species_history == "UNI:FAL")
predict(wloss_X_history_mod, newdata = pca.data.swapval, interval = "prediction")
predict(wloss_X_current_mod, newdata = pca.data.swapval, interval = "prediction")

# compare with full info model CANT USE ANOVA BECAUSE NOT NESTED
models <- list(wloss_X_history_mouse_mod, wloss_X_current_simpstrain_mod, 
               wloss_X_mod, wloss_X_history_mod)
model.names <- c("wloss_X_history_mouse_mod", "wloss_X_current_simpstrain_mod", "wloss_X_mod", 
                  "wloss_X_history_mod")
aictab(cand.set = models, modnames = model.names, sort = T)



# 1. Add predictions 
pred.int <- predict(wloss_X_history_mod, newdata = pca.data.swapval, interval = "prediction")
mydata <- cbind(pca.data.swapval, pred.int)
# 2. Regression line + confidence intervals
library("ggplot2")
p <- ggplot(mydata, aes(maximum_weight_loss_challenge, fit)) +
  geom_point(aes(color = "predicted")) +
  stat_smooth(method = lm)
# 3. Add prediction intervals

p + geom_point(data = mydata, 
               mapping = aes(x = fit, y = maximum_weight_loss_challenge, color = "original")) +
  stat_smooth(method = lm) +
  geom_line(aes(y = lwr), color = "red", linetype = "dashed")+
  geom_line(aes(y = upr), color = "red", linetype = "dashed")


########################## add random forest











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

