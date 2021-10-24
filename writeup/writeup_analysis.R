# import libraries
library(ggplot2)
library(missMDA)
library(FactoMineR)
library(dplyr)
library(tidyr)
library(tidyverse)
library(magrittr)
library(ggpubr)
library(ggeffects)
library(effects)
library(car)
library(AICcmodavg)
library(RColorBrewer)
# library(textshape)
library(tibble)

# read in the babies
IT <- read.csv("./Documents/GitHub/PhD/writeup/IT.csv")
IT$X <- NULL
IMCH <- read.csv("./Documents/GitHub/PhD/writeup/IMCH.csv")
IMCH$X <- NULL
names(IMCH)[names(IMCH) == "OOC"] <- "OPG"

# make more general columns
IT$eimeria_species_primary[IT$primary_infection == "E64"] <- "FER"
IT$eimeria_species_challenge[IT$challenge_infection == "E64"] <- "FER"
IT$eimeria_species_primary[IT$primary_infection == "E88"] <- "FAL"
IT$eimeria_species_challenge[IT$challenge_infection == "E88"] <- "FAL"
IT$eimeria_species_primary[IT$primary_infection == "Eflab"] <- "FAL"
IT$eimeria_species_challenge[IT$challenge_infection == "Eflab"] <- "FAL"
IT$eimeria_species_primary[IT$primary_infection == "E139"] <- "FER"
IT$eimeria_species_challenge[IT$challenge_infection == "E139"] <- "FER"
IT$eimeria_species_primary[IT$primary_infection == "UNI"] <- "UNI"
IT$eimeria_species_challenge[IT$challenge_infection == "UNI"] <- "UNI"

IT$eimeria_species_history[IT$infection_history == "UNI:UNI"] <- "UNI:UNI"
IT$eimeria_species_history[IT$infection_history == "E88:E88"] <- "FAL:FAL"
IT$eimeria_species_history[IT$infection_history == "E88:E64"] <- "FAL:FER"
IT$eimeria_species_history[IT$infection_history == "E88:UNI"] <- "FAL:UNI"
IT$eimeria_species_history[IT$infection_history == "Eflab:UNI"] <- "FAL:UNI"
IT$eimeria_species_history[IT$infection_history == "Eflab:E88"] <- "FAL:FAL"
IT$eimeria_species_history[IT$infection_history == "E64:E64"] <- "FER:FER"
IT$eimeria_species_history[IT$infection_history == "E64:E88"] <- "FER:FAL"
IT$eimeria_species_history[IT$infection_history == "E64:UNI"] <- "FER:UNI"
IT$eimeria_species_history[IT$infection_history == "UNI:E88"] <- "UNI:FAL"
IT$eimeria_species_history[IT$infection_history == "UNI:E64"] <- "UNI:FER"

##################################################      start PCA on FACS     ##########################################################

IT_PCA <- subset(IT, IT$dpi == 8)
IT_PCA <- subset(IT_PCA, IT_PCA$infection == "challenge")
IT_PCA <- select(IT_PCA, EH_ID, #challenge_infection, mouse_strain, infection_history, 
                  CD4, Treg, Div_Treg, Treg17, Th1, Div_Th1, Th17, Div_Th17, CD8, Act_CD8, Div_Act_CD8, IFNy_CD4, IFNy_CD8, 
                  # CASP1, CXCL9, CXCR3, IDO1,  IFNG,  IL10,  IL12A, IL13, IL1RN, IRGM1, MPO, MUC2, MUC5AC, MYD88, NCR1, PRF1,  RETNLB,  
                  # SOCS1, TICAM1, TNF, IL6, IL17A)
)

IT_PCA <- distinct(IT_PCA)
unique(IT_PCA$EH_ID)
IT_PCA <- IT_PCA %>% remove_rownames %>% column_to_rownames(var="EH_ID")
IT_PCA <- na.omit(IT_PCA)

PCA(IT_PCA, scale.unit=TRUE, ncp=5, graph=T)

pca <- prcomp(IT_PCA, center = TRUE, scale = TRUE)
pca.var <- pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100, 1)
barplot(pca.var.per, main = "Scree Plot", xlab = "Principal component",
        ylab = "Percent Variation")

pca.data <- data.frame(Sample=rownames(pca$x),
                       X = pca$x[,1],
                       Y = pca$x[,2])
names(pca.data)[names(pca.data) == "Sample"] <- "EH_ID"
pca.data <- remove_rownames(pca.data)
IT_PCA <- select(IT, EH_ID, experiment, challenge_infection, mouse_strain, infection_history,
                  eimeria_species_challenge, eimeria_species_history, Eim_MC,
                  relative_weight, CD4, Treg, Div_Treg, Treg17, Th1, Div_Th1, Th17, Div_Th17, CD8, Act_CD8, Div_Act_CD8, 
                  IFNy_CD4, IFNy_CD8#, CASP1, CXCL9, CXCR3, IDO1,  IFNG,  IL10, IL12A, IL13, IL1RN, IRGM1, MPO, MUC2, 
                  #MUC5AC, MYD88, NCR1, PRF1,  RETNLB, SOCS1, TICAM1, TNF, IL6, IL17A)
)
IT_PCA <- subset(IT, IT$dpi == 8)
IT_PCA <- subset(IT_PCA, IT_PCA$infection == "challenge")
pca.data <- merge(pca.data, IT_PCA, by = "EH_ID")

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

##################################################      start PCA on IFCs     ##########################################################
IT_PCA1 <- subset(IT, IT$dpi == 8)
IT_PCA1 <- subset(IT_PCA1, IT_PCA1$infection == "challenge")
IT_PCA1 <- select(IT_PCA1, EH_ID, #challenge_infection, mouse_strain, infection_history, 
                   #CD4, Treg, Div_Treg, Treg17, Th1, Div_Th1, Th17, Div_Th17, CD8, Act_CD8, Div_Act_CD8, IFNy_CD4, IFNy_CD8, 
                   CASP1, CXCL9, CXCR3, IDO1,  IFNG,  IL10,  IL12A, IL13, IL1RN, IRGM1, MPO, MUC2, MUC5AC, MYD88, NCR1, PRF1,  RETNLB,  
                   SOCS1, TICAM1, TNF, IL6, IL17A)

IT_PCA1 <- distinct(IT_PCA1)
EH_ID_list <- unique(IT_PCA1$EH_ID)
IT_PCA1 <- IT_PCA1 %>% remove_rownames %>% column_to_rownames(var="EH_ID")
pca1 <- imputePCA(IT_PCA1[,c(1:22)], ncp = 1, method = "Regularized")
pca1 <- data.frame(matrix(unlist(pca1['completeObs']), nrow=154),stringsAsFactors=FALSE)
gene_names <- c("CASP1", "CXCL9", "CXCR3", "IDO1",  "IFNG",  "IL10",  "IL12A", "IL13", "IL1RN", "IRGM1", "MPO", "MUC2", "MUC5AC", "MYD88",
                "NCR1", "PRF1",  "RETNLB", "SOCS1", "TICAM1", "TNF", "IL6", "IL17A")
row.names(pca1) <- EH_ID_list
colnames(pca1) <- gene_names
gene_values <- pca1

PCA(pca1, scale.unit=TRUE, ncp=5, graph=T)

pca1 <- prcomp(pca1, center = TRUE, scale = TRUE)
pca.var1 <- pca1$sdev^2
pca.var.per1 <- round(pca.var1/sum(pca.var1)*100, 1)
barplot(pca.var.per1, main = "Scree Plot", xlab = "Principal component",
        ylab = "Percent Variation")

pca.data1 <- data.frame(Sample=rownames(pca1$x),
                       X = pca1$x[,1],
                       Y = pca1$x[,2])
names(pca.data1)[names(pca.data1) == "Sample"] <- "EH_ID"
pca.data1 <- remove_rownames(pca.data1)


IT_PCA1 <- select(IT, EH_ID, experiment, challenge_infection, mouse_strain, infection_history,
                   eimeria_species_challenge, eimeria_species_history, Eim_MC,
                   relative_weight, dpi, infection #CD4, Treg, Div_Treg, Treg17, Th1, Div_Th1, Th17, Div_Th17, CD8, Act_CD8, Div_Act_CD8, 
                   # IFNy_CD4, IFNy_CD8, 
                   # CASP1, CXCL9, CXCR3, IDO1,  IFNG,  IL10, IL12A, IL13, IL1RN, IRGM1, MPO, MUC2, 
                   # MUC5AC, MYD88, NCR1, PRF1,  RETNLB, SOCS1, TICAM1, TNF, IL6, IL17A
)

IT_PCA1 <- subset(IT_PCA1, IT_PCA1$dpi == 8)
IT_PCA1 <- subset(IT_PCA1, IT_PCA1$infection == "challenge")
IT_PCA1 <- distinct(IT_PCA1)
IT_PCA1 <- cbind(IT_PCA1, gene_values)

pca.data1 <- merge(pca.data1, IT_PCA1, by = "EH_ID")

loading_scores1 <- pca1$rotation[,1]
scores1 <- abs(loading_scores1)
score_ranked1 <- sort(scores1, decreasing = T)
top101 <- names(score_ranked1[1:22])
pca1$rotation[top101,1]

pca.var.per.df1 <- cbind(pca.var.per1, top101)

pca.var.per.df1 <- data.frame(pca.var.per.df1)
pca.var.per.df1$pca.var.per1 <- as.numeric(as.character(pca.var.per.df1$pca.var.per1))
pca.var.per.df1$top101 <- factor(pca.var.per.df1$top101,levels = top101)

ggplot(pca.var.per.df1, aes(top101, pca.var.per1, fill = top101)) +
  geom_col() +
  theme(strip.text.x = element_text(size = 14, face = "bold"),
        strip.text.y = element_text(size = 14, face = "bold")) +
  scale_x_discrete(name ="Genes") +
  scale_y_discrete(name ="% of variation explained") +
  labs(fill = "") +
  
  geom_text(aes(label = pca.var.per1)) +
  theme_bw()

ggplot(data = subset(pca.data1, !is.na(pca.data1$relative_weight)), aes(x = X, y = Y, label = EH_ID, color = challenge_infection)) +
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
  ggtitle("Gene expression PCA by challenge infection") +
  labs(color='infection') 
# keep only the important bits for main table
PCA_keep <- distinct(pca.data[,c("EH_ID", "X", "Y", "dpi")])
PCA1_keep <- distinct(pca.data1[,c("EH_ID", "X", "Y", "dpi")])

################################################# transform the data
# this is very cool but we're using relative weight from day 8, which won't tell us much
# make maximum OPG (OPG when available) and maximum weightloss columns ad test those
# confirmed, no primary challenges
max_weight_primary <- subset(IT, IT$infection == "primary")
max_weight_primary <- max_weight_primary %>% select(EH_ID, relative_weight, primary_infection, 
                                             challenge_infection, infection_history) 
max_weight_challenge <- subset(IT, IT$infection == "challenge")
max_weight_challenge <- max_weight_challenge %>% select(EH_ID, relative_weight, primary_infection, 
                                                 challenge_infection, infection_history)
# except that this uses only data from challenge for primary:uni
max_weight_primary <- max_weight_primary[!is.na(max_weight_primary$relative_weight), ]
max_weight_challenge <- max_weight_challenge[!is.na(max_weight_challenge$relative_weight),]
# creates weird 1 sample with 0 so remove it
# max_weight_challenge <- max_weight_challenge[ !(max_weight_challenge$EH_ID %in% "LM0350"),]
# max_OPG<- WDS %>% select(EH_ID, OPG)

# check we don't lose anything here
# unique(max_weight_challenge$EH_ID) 215
# unique(max_weight_primary$EH_ID) 127

# max_OPG$OPG <- as.numeric(as.character(max_OPG$OPG))
max_weight_primary <- max_weight_primary %>% dplyr::group_by(EH_ID) %>% dplyr::summarise(relative_weight = min(relative_weight, na.rm = T))
max_weight_challenge <- max_weight_challenge %>% dplyr::group_by(EH_ID) %>% dplyr::summarise(relative_weight = min(relative_weight, na.rm = T))
# max_OPG <- max_OPG %>% dplyr::group_by(EH_ID) %>% dplyr::summarise(OPG = max(OPG, na.rm = T))
names(max_weight_primary)[names(max_weight_primary) == "relative_weight"] <- "maximum_weight_loss_primary"
names(max_weight_challenge)[names(max_weight_challenge) == "relative_weight"] <- "maximum_weight_loss_challenge"
# now this gets complicated. infection_history needs maximum weightloss, but tor eflect primary infections which
# are in challenge. and current infections need chalenge weight loss. Otherwise infected:uninfected mess up models

#################### same for oocysts
max_oocyst_primary <- subset(IT, IT$infection == "primary")
max_oocyst_primary <-  max_oocyst_primary %>% select(EH_ID, OPG, primary_infection, 
                                             challenge_infection, infection_history) 

max_oocyst_challenge <- subset(IT, IT$infection == "challenge")
max_oocyst_challenge <- max_oocyst_challenge %>% select(EH_ID, OPG, primary_infection, 
                                                 challenge_infection, infection_history)
# except that this uses only data from challenge for primary:uni
max_oocyst_primary <- max_oocyst_primary[!is.na(max_oocyst_primary$OPG), ]
max_oocyst_challenge <- max_oocyst_challenge[!is.na(max_oocyst_challenge$OPG),]

# max_OPG$OPG <- as.numeric(as.character(max_OPG$OPG))
max_oocyst_primary <- max_oocyst_primary %>% dplyr::group_by(EH_ID) %>% dplyr::summarise(OPG = max(OPG, na.rm = T))
max_oocyst_challenge <- max_oocyst_challenge %>% dplyr::group_by(EH_ID) %>% dplyr::summarise(OPG = max(OPG, na.rm = T))
names(max_oocyst_primary)[names(max_oocyst_primary) == "OPG"] <- "maximum_OPG_primary"
names(max_oocyst_challenge)[names(max_oocyst_challenge) == "OPG"] <- "maximum_OPG_challenge"

ITa <- merge(IT, max_weight_primary, by = "EH_ID")
ITa <- merge(ITa, max_weight_challenge, by = "EH_ID")
ITa <- merge(ITa, max_oocyst_primary, by = "EH_ID")
ITa <- merge(ITa, max_oocyst_challenge, by = "EH_ID")

# anything that is infected:uninfected needs to be the in the final maximum weightloss column from primary_max_WLoss
# can do min across both columns but that didnt work because its non disciminating to our criteria so scrapping
############## THIS WILL ONLY WORK FOR infection_history MODELS! (pca.data.swap1)
ITa$eimeria_species_history <- paste(ITa$eimeria_species_primary, ITa$eimeria_species_challenge, sep = ":")

# IT1 <- subset(ITa, infection_history == "E64:UNI" | 
#                            infection_history == "E88:UNI"|
#                            infection_history == "Eflab:UNI")
# drop challenge weightloss column to avoic UNI weight in infection history

# name the columns the same and merge
# IT2 <- subset(ITa, infection_history == "UNI:UNI" | 
#                            infection_history == "E88:E88"|
#                            infection_history == "E88:E64"|
#                            infection_history == "Eflab:E88"|
#                            infection_history == "E64:E64"|
#                            infection_history == "E64:E88"|
#                            infection_history == "UNI:E88"|
#                            infection_history == "UNI:E64")
# 
# 
# ITb <- rbind(IT1, IT2)
ITchallenge <- subset(ITa, ITa$infection == "challenge")

# for this the FAL:UNI and UNI:FAL are the same, join under new column and model
ITchallenge$relevant_history[ITchallenge$eimeria_species_history == "UNI:UNI"] <- "UNI"
ITchallenge$relevant_history[ITchallenge$eimeria_species_history == "UNI:FAL"] <- "FAL"
ITchallenge$relevant_history[ITchallenge$eimeria_species_history == "FAL:UNI"] <- "FAL"
ITchallenge$relevant_history[ITchallenge$eimeria_species_history == "FER:UNI"] <- "FER"
ITchallenge$relevant_history[ITchallenge$eimeria_species_history == "UNI:FER"] <- "FER"
ITchallenge$relevant_history[ITchallenge$eimeria_species_history == "FAL:FAL"] <- "FAL:FAL"
ITchallenge$relevant_history[ITchallenge$eimeria_species_history == "FER:FER"] <- "FER:FER"
ITchallenge$relevant_history[ITchallenge$eimeria_species_history == "FER:FAL"] <- "FER:FAL"
ITchallenge$relevant_history[ITchallenge$eimeria_species_history == "FAL:FER"] <- "FAL:FER"

ITchallenge$relevant_history <- factor(ITchallenge$relevant_history,
                                         levels = c("UNI", "FAL", "FER","FER:FER", 
                                                    "FER:FAL", "FAL:FAL", "FAL:FER"))
# organize to compare to SWISS
strains <- unique(ITchallenge$mouse_strain)
ITchallenge$mouse_strain <- factor(ITchallenge$mouse_strain,
                                     levels = strains)

# Define the number of colors you want
nb.cols <- 18
mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(nb.cols)

############### weight loss explanation
# use dataset without SWISS mice (validation and training)
# ITchallengeval <- subset(ITchallenge, ITchallenge$mouse_strain == "SWISS")
# ITchallenge <- subset(ITchallenge, !ITchallenge$mouse_strain == "SWISS")

# now build these to reflext wild better (only current infections and hybrid, mmd and mmm mouse strains)
# already got secondary_species so let's code the mouse strains to simpler groups
ITchallenge$simple_strain[ITchallenge$mouse_strain == "NMRI"] <- "SWISS"
ITchallenge$simple_strain[ITchallenge$mouse_strain == "SCHUNT_SCHUNT"] <- "Mmd"
ITchallenge$simple_strain[ITchallenge$mouse_strain == "PWD_PWD"] <- "Mmm"
ITchallenge$simple_strain[ITchallenge$mouse_strain == "BUSNA_STRA"] <- "Hybrid"
ITchallenge$simple_strain[ITchallenge$mouse_strain == "STRA_BUSNA"] <- "Hybrid"
ITchallenge$simple_strain[ITchallenge$mouse_strain == "PWD_SCHUNT"] <- "Hybrid"
ITchallenge$simple_strain[ITchallenge$mouse_strain == "STRA_STRA"] <- "Mmd"
ITchallenge$simple_strain[ITchallenge$mouse_strain == "STRA_SCHUNT"] <- "Mmd"
ITchallenge$simple_strain[ITchallenge$mouse_strain == "PWD_BUSNA"] <- "Mmm"
ITchallenge$simple_strain[ITchallenge$mouse_strain == "SCHUNT_PWD"] <- "Hybrid"
ITchallenge$simple_strain[ITchallenge$mouse_strain == "SCHUNT_STRA"] <- "Mmd"
ITchallenge$simple_strain[ITchallenge$mouse_strain == "BUSNA_BUSNA"] <- "Mmm"
ITchallenge$simple_strain[ITchallenge$mouse_strain == "BUSNA_PWD"] <- "Mmm"

challenge <- select(ITchallenge, EH_ID, experiment, relative_weight, dpi, primary_infection, challenge_infection, mouse_strain, 
                    infection_history, OPG, IFNy_CEWE, CD4, IFNy_CEWE, CD4, Treg , Div_Treg, Treg17 ,Th1 ,Div_Th1 ,Th17, Div_Th17, 
                    CD8, Act_CD8, Div_Act_CD8, IFNy_CD4, IFNy_CD8, CASP1, CXCL9, CXCR3 ,IDO1 , IFNG , IL10 ,IL12A ,IL13 , IL1RN ,
                    IRGM1 , MPO ,  MUC2 , MUC5AC ,MYD88 , NCR1 , PRF1 , RETNLB ,SOCS1 , TICAM1 ,TNF ,  IL6 ,  IL17A , Eim_MC ,delta ,
                    eimeria_species_primary, eimeria_species_challenge, eimeria_species_history, maximum_weight_loss_primary,
                    maximum_weight_loss_challenge, maximum_OPG_primary, maximum_OPG_challenge, relevant_history, simple_strain )
# add PCA components (renamed)
colnames(PCA_keep) <- c("EH_ID", "FACS_PCA1", "FACS_PCA2", "dpi")
colnames(PCA1_keep) <- c("EH_ID", "GE_PCA1", "GE_PCA2", "dpi")

challenge <- full_join(challenge, PCA_keep)
challenge <- full_join(challenge, PCA1_keep)

# save
write.csv(challenge, "C:/Users/exemp/OneDrive/Documents/challenge.csv")
write.csv(challenge, "./Documents/GitHub/PhD/writeup/challenge.csv")
# cleanup
rm(list = c("df", "IMCH", "IT", "ITa", "ITchallenge", "max_oocyst_challenge", "loading_scores", "loading_scores1", "score_ranked", 
            "score_ranked1", "max_oocyst_primary", "max_weight_challenge", "max_weight_primary", "scores", "scores1", "top10", "top101"))
rm(list=ls(pattern="PCA"))
rm(list=ls(pattern="pca"))

#####################################################################################################
#####################compare more than visualy, AIC tab for the relevant models######################
# EXPLORATION so full dataset









# order factors first
challenge$simple_strain <- factor(challenge$simple_strain,
                                      levels = c("SWISS", "Mmm", "Mmd", "Hybrid"))
challenge$secondary_species <- factor(challenge$eimeria_species_challenge,
                                          levels = c("UNI", "FER", "FAL"))
challenge8 <- subset(challenge, challenge$dpi == 8)
# build basic models
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


