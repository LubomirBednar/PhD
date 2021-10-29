# import libraries
library(ggplot2)
library(dplyr)
# library(tidyr)
library(tidyverse)
library(AICcmodavg)
# library(RColorBrewer)
# library(textshape)
# library(tibble)
library(caret)


set.seed(8)

# read in the babies
IT <- read.csv("./Documents/GitHub/PhD/writeup/IT.csv")
IT$X <- NULL
# not sure where the inf came from but lets remove it
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
library(missMDA)
IT_PCA <- subset(IT, IT$dpi == 8)
IT_PCA <- subset(IT_PCA, IT_PCA$infection == "challenge")
IT_PCA <- select(IT_PCA, EH_ID, #challenge_infection, mouse_strain, infection_history, 
                  CD4, Treg, Div_Treg, Treg17, Th1, Div_Th1, Th17, Div_Th17, CD8, Act_CD8, Div_Act_CD8, IFNy_CD4, IFNy_CD8 
                  # CASP1, CXCL9, CXCR3, IDO1,  IFNG,  IL10,  IL12A, IL13, IL1RN, IRGM1, MPO, MUC2, MUC5AC, MYD88, NCR1, PRF1,  RETNLB,  
                  # SOCS1, TICAM1, TNF, IL6, IL17A)
)

IT_PCA <- distinct(IT_PCA)
unique(IT_PCA$EH_ID)
IT_PCA <- IT_PCA %>% remove_rownames %>% column_to_rownames(var="EH_ID")
IT_PCA <- na.omit(IT_PCA)
library(FactoMineR)
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
facs_names <- c("CD4", "Treg", "Div_Treg", "Treg17", "Th1", "Div_Th1", "Th17", "Div_Th17", "CD8", "Act_CD8", "Div_Act_CD8", 
                "IFNy_CD4", "IFNy_CD8")
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

ggplot(data = CHALLENGE, aes(x = FACS_PCA1, y = FACS_PCA2, label = EH_ID, color = maximum_weight_loss_challenge)) +
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
  labs(color='weight retained') 

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

ggplot(data = subset(CHALLENGE, !is.na(CHALLENGE$relative_weight)), aes(x = GE_PCA1, y = GE_PCA2, label = EH_ID, color = maximum_weight_loss_challenge)) +
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
  ggtitle("Gene expression PCA by weight loss") +
  labs(color='weight retained') 
# keep only the important bits for main table
PCA_keep <- distinct(pca.data[,c("EH_ID", "X", "Y", "dpi", "CD4", "Treg", "Div_Treg", "Treg17", "Th1", "Div_Th1", "Th17", "Div_Th17",
                                 "CD8", "Act_CD8", "Div_Act_CD8", "IFNy_CD4", "IFNy_CD8")])
PCA1_keep <- distinct(pca.data1[,c("EH_ID", "X", "Y", "dpi", "CASP1", "CXCL9", "CXCR3", "IDO1",  "IFNG",  "IL10",  "IL12A", "IL13", 
                                   "IL1RN", "IRGM1", "MPO", "MUC2", "MUC5AC", "MYD88","NCR1", "PRF1",  "RETNLB", "SOCS1", "TICAM1", 
                                   "TNF", "IL6", "IL17A")])

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

# # Define the number of colors you want
# nb.cols <- 18
# mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(nb.cols)

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

CHALLENGE <- select(ITchallenge, EH_ID, experiment, relative_weight, dpi, primary_infection, challenge_infection, mouse_strain, 
                    infection_history, OPG, IFNy_CEWE, Eim_MC, delta, eimeria_species_primary, eimeria_species_challenge, 
                    eimeria_species_history, maximum_weight_loss_primary, maximum_weight_loss_challenge, maximum_OPG_primary, 
                    maximum_OPG_challenge, relevant_history, simple_strain)

ggplot(subset(CHALLENGE, !is.na(CHALLENGE$eimeria_species_challenge)), aes(x = dpi, y = OPG, color = eimeria_species_challenge)) +
  geom_jitter(width = 0.2) +
#  geom_smooth() +
  labs(y = "weight relative to dpi 0", x = "days post infection", color = "infecting species") +
  theme(axis.text=element_text(size=12, face = "bold"),
        title = element_text(size = 16, face = "bold"),
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        strip.text.y = element_text(size = 14, face = "bold"),
        legend.text = element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold")) +
  ggtitle("Oocyst shedding during challenge infections")



# add PCA components (renamed)
names(PCA_keep)[names(PCA_keep) == "X"] <- "FACS_PCA1"
names(PCA_keep)[names(PCA_keep) == "Y"] <- "FACS_PCA2"

names(PCA1_keep)[names(PCA1_keep) == "X"] <- "GE_PCA1"
names(PCA1_keep)[names(PCA1_keep) == "Y"] <- "GE_PCA2"

CHALLENGE <- full_join(CHALLENGE, PCA_keep)
CHALLENGE <- full_join(CHALLENGE, PCA1_keep)

# save
write.csv(CHALLENGE, "C:/Users/exemp/OneDrive/Documents/challenge.csv")
write.csv(CHALLENGE, "./Documents/GitHub/PhD/writeup/challenge.csv")
# cleanup
rm(list = c("df", "IMCH", "IT", "ITa", "max_oocyst_challenge", "loading_scores", "loading_scores1", "score_ranked", 
            "score_ranked1", "max_oocyst_primary", "max_weight_challenge", "max_weight_primary", "scores", "scores1", "top10", "top101"))
rm(list=ls(pattern="PCA"))
rm(list=ls(pattern="pca"))

#####################################################################################################
##################### compare more than visually, AIC tab for the relevant models ######################
# order factors first
CHALLENGE$simple_strain <- factor(CHALLENGE$simple_strain,
                                      levels = c("SWISS", "Mmm", "Mmd", "Hybrid"))
CHALLENGE$OPG <- as.numeric(CHALLENGE$OPG)

# HINT: 
# In statistics, linear regression is used to model a relationship between a continuous dependent variable and one or more 
# independent variables. The independent variable can be either categorical or numerical. The case when we have only one 
# independent variable then it is called as simple linear regression. If we have more than one independent variable, 
# then it is called as multivariate regression.

# might need a naming convention:
# P1 = FACS_PCA1, P2 = FACS_PCA2, P3 = GE_PCA1, P4 = GE_PCA2,
# M1 = simple_strain, M2 = mouse_strain,
# W1 = maximum_weight_loss_challenge,  W2 = relative_weight,
# O1 = maximum_OPG_challenge, O2 = OPG,
# E1 = eimeria_species_challenge, E2 = challenge_infection, 
# Build order: dependent variable ~ main independent variable + additional regressors

######### build basic models looking at what affects weight loss and maximum weight loss
# We think the weight loss or inversely weight retention will be most influenced by 
# 1) immunological parameters
# 2) parasite species
# 3) mouse strain
# 4) parasite success (OPG)
# 5) infection history

# since weight loss is of interest:
# but can't keep relative weight because of correlation problems
CHALLENGE$relative_weight <- NULL
CHALLENGE$EH_ID <- NULL
CHALLENGE$maximum_OPG_primary <- NULL
CHALLENGE$maximum_weight_loss_primary <- NULL
CHALLENGE$OPG <- NULL

CHALLENGE_immuno <- subset(CHALLENGE, CHALLENGE$dpi == 8)

CHALLENGE_immuno <- mutate_if(CHALLENGE_immuno, is.character, as.factor)
CHALLENGE_immuno <- distinct(CHALLENGE_immuno)
CHALLENGE_immuno <- data.frame(CHALLENGE_immuno)

# remove all NAs, this cuts out huge amount of samples because of missing FACS
CHALLENGE_immuno <- na.omit(CHALLENGE_immuno)
# so let's make other data sets for training and testing to split
CHALLENGE_facs <- subset(CHALLENGE, CHALLENGE$dpi == 8 & !is.na(CHALLENGE$FACS_PCA1))
CHALLENGE_facs <- distinct(CHALLENGE_facs)
CHALLENGE_ge <- subset(CHALLENGE, CHALLENGE$dpi == 8 & !is.na(CHALLENGE$GE_PCA1))
CHALLENGE_ge <- distinct(CHALLENGE_ge)
# Split data into train and test
num_names <- c("maximum_weight_loss_challenge", "IFNy_CEWE", "delta", "maximum_OPG_challenge")
other_names <- c("experiment", "challenge_infection", "mouse_strain", "infection_history", "Eim_MC", "simple_strain", "eimeria_species_challenge",
                 "eimeria_species_history", "relevant_history")
pca_names <- c("FACS_PCA1", "FACS_PCA2", "GE_PCA1", "GE_PCA2")
selected_names <- c("eimeria_species_challenge", "relevant_history", "simple_strain")

immuno_index <- sample(2, nrow(CHALLENGE_immuno), replace = TRUE, prob = c(0.7, 0.3))
train_immuno <- CHALLENGE_immuno[immuno_index==1,]
test_immuno <- CHALLENGE_immuno[immuno_index==2,]

train_immunorf <- dplyr::select(train_immuno, maximum_weight_loss_challenge, all_of(selected_names), all_of(pca_names), all_of(gene_names), 
                              all_of(facs_names))
test_immunorf <- dplyr::select(test_immuno, maximum_weight_loss_challenge, all_of(selected_names), all_of(pca_names), all_of(gene_names), 
                               all_of(facs_names))
############################################### add random forest
library(randomForest)
library(e1071)
library(caret)
# Define the control for F-fold cross-validation
trControl <- trainControl(method = "cv",
                          number = 10, # 10++
                          search = "grid")
# values of mtry that are close to the total number of variables in the model
# may weaken the forest by making the individual decision trees more correlated
# this influences time consumption A LOT!
tuneGrid <- expand.grid(.mtry = c(1 : 100))
# Run the model
rf_default <- caret::train(maximum_weight_loss_challenge ~ .,
                    data = train_immunorf,
                    method = "rf",
                    trControl = trControl)
print(rf_default)
def_predict <- predict(rf_default, test_immunorf)
def_predict <- as.vector(def_predict)
# mtry  RMSE      Rsquared   MAE     
# 2     5.135472  0.4411085  4.406684
# 53    5.036686  0.4586127  4.400776
# 105   5.067056  0.4356753  4.396825
# 
# RMSE was used to select the optimal model using the smallest value.
# The final value used for the model was mtry = 53.
error_def_predict <- def_predict - test_immunorf$maximum_weight_loss_challenge
error_def_predict <- sqrt(mean((error_def_predict)^2))
print(error_def_predict) # 4.137446
print(varImp(rf_default))

# top 10
# TNF                    100.00
# TICAM1                  96.12
# RETNLB                  88.66
# CD4                     32.91
# CD8                     31.65
# Treg                    25.87
# maximum_OPG_challenge   24.70
# CASP1                   22.42
# CXCR3                   21.69
# Treg17                  18.81

varImp(rf_default)$importance %>% 
  as.data.frame() %>%
  rownames_to_column() %>%
  arrange(Overall) %>%
  mutate(rowname = forcats::fct_inorder(rowname )) %>%
  ggplot()+
  geom_col(aes(x = rowname, y = Overall, fill = rowname))+
  coord_flip()+
  theme_bw() + 
  theme(legend.position = "none")

RFD_importance <- varImp(rf_default)$importance %>% 
  as.data.frame() %>%
  rownames_to_column() %>%
  arrange(Overall) %>%
  mutate(rowname = forcats::fct_inorder(rowname )) 
RFD_20 <- RFD_importance[c(50:30),]

ggplot(RFD_20) +
  geom_col(aes(x = rowname, y = Overall, fill = rowname))+
  labs(x= "weight loss influencing factors", y = "degree of influence") +
  coord_flip()+
  theme_bw() + 
  theme(legend.position = "none")


# pick best mtry value by changing tuneGrid
rf_mtry <- train(maximum_weight_loss_challenge ~ .,
                 data = train_immunorf,
                 method = "rf",
                 tuneGrid = tuneGrid,
                 trControl = trControl,
                 importance = TRUE,
                 nodesize = 14,
                 ntree = 300)
print(rf_mtry)
rf_mtry$bestTune$mtry # mtry = 5, 92???
best_mtry <- rf_mtry$bestTune$mtry
store_maxnode <- list()
tuneGrid <- expand.grid(.mtry = best_mtry)
for (maxnodes in c(5: 25)) {
  set.seed(8)
  rf_maxnode <- train(maximum_weight_loss_challenge ~ .,
                      data = train_immunorf,
                      method = "rf",
                      tuneGrid = tuneGrid,
                      trControl = trControl,
                      importance = TRUE,
                      nodesize = 14,
                      maxnodes = maxnodes,
                      ntree = 300)
  current_iteration <- toString(maxnodes)
  store_maxnode[[current_iteration]] <- rf_maxnode
}
results_mtry <- resamples(store_maxnode)
summary(results_mtry)
# maxnode 7 has lowest RSME and 9 for MAE

store_maxtrees <- list()
for (ntree in c(250, 300, 350, 400, 450, 500, 550, 600, 800, 1000, 2000, 5000)) {
  set.seed(8)
  rf_maxtrees <- train(maximum_weight_loss_challenge ~ .,
                       data = train_immunorf,
                       method = "rf",
                       tuneGrid = tuneGrid,
                       trControl = trControl,
                       importance = TRUE,
                       nodesize = 14,
                       maxnodes = 7,
                       ntree = ntree)
  key <- toString(ntree)
  store_maxtrees[[key]] <- rf_maxtrees
}
results_tree <- resamples(store_maxtrees)
summary(results_tree)
# ntree smallest for RMSE at 300, and 350 for MAE, stick with 350

# R^2 used for explaining how well the independent variables in the linear regression model 
# explains the variability in the dependent variable
# The lower value of MAE and RMSE implies higher accuracy of a regression model 
# based on low MAE, low RMSE and high R^2 => mtry = 5, node = 9, ntree = 350

RF1 <- train(maximum_weight_loss_challenge ~ .,
                data = train_immunorf,
                method = "rf",
                tuneGrid = tuneGrid,
                trControl = trControl,
                importance = TRUE,
                nodesize = 14,
                ntree = 5000,
                maxnodes = 7)
print(RF1)
RF1predict <- predict(RF1, test_immunorf)
# RMSE      Rsquared   MAE     
# 4.637729  0.3814265  4.069071
error_RF1_predict <- RF1predict - test_immunorf$maximum_weight_loss_challenge
error_RF1_predict <- sqrt(mean((error_RF1_predict)^2))
print(error_RF1_predict) # 2.643995
print(varImp(RF1))
# top 10
# TNF                             100.00
# RETNLB                           98.71
# TICAM1                           97.66
# Treg                             96.65
# MYD88                            86.57
# IL13                             78.29
# GE_PCA1                          77.48
# maximum_OPG_challenge            74.30
# IL12A                            71.54
# IL6                              69.51

RF1_importance <- varImp(RF1)$importance %>% 
  as.data.frame() %>%
  rownames_to_column() %>%
  arrange(Overall) %>%
  mutate(rowname = forcats::fct_inorder(rowname )) 
RF1_20 <- RF1_importance[c(50:30),]

ggplot(RF1_20) +
  geom_col(aes(x = rowname, y = Overall, fill = rowname))+
  labs(x= "weight loss influencing factors", y = "degree of influence") +
  coord_flip()+
  theme_bw() + 
  theme(legend.position = "none")
  

# new dataset after LMs
RF2_train <- dplyr::select(train_immuno, maximum_weight_loss_challenge, all_of(selected_names), all_of(pca_names), 
                           all_of(gene_names), all_of(facs_names))
RF2_test <- dplyr::select(test_immuno, maximum_weight_loss_challenge, all_of(selected_names), all_of(pca_names), 
                          all_of(gene_names), all_of(facs_names))

RF2 <- train(maximum_weight_loss_challenge ~ .,
                    data = RF2_train,
                    method = "rf",
                    tuneGrid = tuneGrid,
                    trControl = trControl,
                    importance = TRUE,
                    nodesize = 10,
                    ntree = 5000,
                    maxnodes = 10)
print(RF2)
RF2predict <- predict(RF2, RF2_test)
error_RF2_predict <- RF2predict - RF2_test$maximum_weight_loss_challenge
error_RF2_predict <- sqrt(mean((error_RF2_predict)^2))
print(error_RF2_predict) # 4.077054
print(varImp(RF2))

RF2_importance <- varImp(RF2)$importance %>% 
  as.data.frame() %>%
  rownames_to_column() %>%
  arrange(Overall) %>%
  mutate(rowname = forcats::fct_inorder(rowname )) 
RF2_20 <- RF2_importance[c(50:30),]

ggplot(RF2_20) +
  geom_col(aes(x = rowname, y = Overall, fill = rowname))+
  labs(x= "weight loss influencing factors", y = "degree of influence") +
  coord_flip()+
  theme_bw() + 
  theme(legend.position = "none") 

# get mean square error 
# take test_immuno$maximum_weight_loss_challenge - prediction
# square the =
# compare the errors with LM against RF1 and rf_default
# RMSE      Rsquared   MAE     
# 2.347099  0.8498383  1.565854
# test direction?


# further random forest insights on interactions
library(randomForestExplainer)
forest <- randomForest::randomForest(maximum_weight_loss_challenge ~ ., data = train_immuno,
                                     ntree = 600, localImp = TRUE, mtry =9, maxnodes = 16)
min_depth_frame <- min_depth_distribution(forest)
head(min_depth_frame, n = 10)
plot_min_depth_distribution(forest) # gives the same result as below but takes longer
plot_min_depth_distribution(min_depth_frame)
plot_min_depth_distribution(min_depth_frame, mean_sample = "relevant_trees", k = 20)



forest2 <- randomForest(maximum_weight_loss_challenge ~ ., data = train_immuno, localImp = TRUE)
min_depth_frame2 <- min_depth_distribution(forest2)
head(min_depth_frame2, n = 10)
plot_min_depth_distribution(forest2) # gives the same result as below but takes longer
plot_min_depth_distribution(min_depth_frame2)
plot_min_depth_distribution(min_depth_frame2, mean_sample = "relevant_trees", k = 20)


######################################## LMs for weight loss
library(MASS)
library(leaps)

# Training model
CHALLENGE_facs$challenge_infection <- factor(CHALLENGE_facs$challenge_infection,
                                             levels = c("UNI", "E64", "E88"))
CHALLENGE_facs$mouse_strain <- factor(CHALLENGE_facs$mouse_strain,
                                             levels = c("NMRI", "PWD_PWD", "BUSNA_PWD", "BUSNA_BUSNA", "SCHUNT_STRA", "SCHUNT_PWD", "PWD_BUSNA",
                                                        "STRA_SCHUNT", "STRA_STRA", "PWD_SCHUNT", "SCHUNT_SCHUNT", "STRA_BUSNA", "BUSNA_STRA"))

CHALLENGE_facs <- mutate_if(CHALLENGE_facs, is.character, as.factor)
CHALLENGE_ge <- mutate_if(CHALLENGE_ge, is.character, as.factor)
CHALLENGE <- mutate_if(CHALLENGE, is.character, as.factor)

facs_index <- sample(2, nrow(CHALLENGE_facs), replace = TRUE, prob = c(0.7, 0.3))
train_facs <- CHALLENGE_facs[facs_index==1,]
test_facs <- CHALLENGE_facs[facs_index==2,]

ge_index <- sample(2, nrow(CHALLENGE_ge), replace = TRUE, prob = c(0.7, 0.3))
train_ge <- CHALLENGE_ge[ge_index==1,]
test_ge <- CHALLENGE_ge[ge_index==2,]


rankifremoved <- sapply(1:ncol(your.matrix), function (x) qr(your.matrix[,-x])$rank)
which(rankifremoved == max(rankifremoved))
########### start with the most basic categories
# linearize <- function(IV, regressors, dataset) {
#   LM <- lm(IV ~ regressors, dataset)
#   print(summary(LM))
# }
# linearize("maximum_weight_loss_challenge", "mouse_strain", train_facs)
#  figure out the passing of arguments
CHALLENGE_para <- subset(ITchallenge, !is.na(ITchallenge$maximum_weight_loss_challenge))
CHALLENGE_para <- subset(CHALLENGE_para, CHALLENGE_para$dpi == 8)
CHALLENGE_para <- distinct(CHALLENGE_para)
para_index <- sample(2, nrow(CHALLENGE_para), replace = TRUE, prob = c(0.7, 0.3))
train_para <- CHALLENGE_para[para_index==1,]
test_para <- CHALLENGE_para[para_index==2,]

LM1 <- lm(maximum_weight_loss_challenge ~ challenge_infection, train_para)
summary(LM1) # ** challenge_infectionUNI    4.268
LM2 <- lm(maximum_weight_loss_challenge ~ infection_history, train_para)
summary(LM2) # * infection_historyE64:UNI   12.93014, infection_historyE88:UNI   14.47939, * infection_historyUNI:UNI
LM3 <- lm(maximum_weight_loss_challenge ~ mouse_strain, train_para)
summary(LM3) # nope
LM4 <- lm(maximum_weight_loss_challenge ~ eimeria_species_challenge, train_para)
summary(LM4) # *** eimeria_species_challengeUNI    6.240
LM5 <- lm(maximum_weight_loss_challenge ~ relevant_history, train_para)
summary(LM5) # . relevant_historyFER:FAL   -5.399 of interest
LM6 <- lm(maximum_weight_loss_challenge ~ simple_strain, train_para)
summary(LM6) # * simple_strainSWISS

selected_names <- c("eimeria_species_challenge", "relevant_history", "simple_strain")

basic_models <- list("LM1" = LM1, "LM2" = LM2, "LM3" = LM3, "LM4" = LM4, "LM5" = LM5, "LM6" = LM6)
aictab(basic_models) # LM1/4 for infection, LM2 for history, LM6 for strain
LM7_data <- dplyr::select(train_facs, maximum_weight_loss_challenge, all_of(selected_names), FACS_PCA1, FACS_PCA2)
LM7_test <- dplyr::select(train_facs, maximum_weight_loss_challenge, all_of(selected_names), FACS_PCA1, FACS_PCA2)
LM7 <- lm(maximum_weight_loss_challenge ~ ., LM7_data)
LM7_step <- stepAIC(LM7, direction = "both")
summary(LM7)


LM8_data <- dplyr::select(train_facs, maximum_weight_loss_challenge, all_of(selected_names), all_of(facs_names))
LM8_test <- dplyr::select(train_facs, maximum_weight_loss_challenge, all_of(selected_names), all_of(facs_names))
LM8 <- lm(maximum_weight_loss_challenge ~ ., LM8_data)
LM8_step <- stepAIC(LM8, direction =  "backward")  
summary(LM8_step)

facs_models <- list("PCA" = LM7_step, "AIC all" = LM8_step)
aictab(facs_models)
# > facs_models <- list("PCA" = LM7, "AIC all" = LM8_step)
# > aictab(facs_models)
# Model selection based on AICc:
#   
#   K   AICc Delta_AICc AICcWt Cum.Wt      LL
# AIC all 14 347.88        0.0   0.93   0.93 -154.69
# PCA      8 353.08        5.2   0.07   1.00 -166.97



LM9_data <- dplyr::select(train_ge, GE_PCA1, GE_PCA2, all_of(selected_names), maximum_weight_loss_challenge)
LM9_data <- na.omit(LM9_data)
LM9_test <- dplyr::select(test_ge, GE_PCA1, GE_PCA2, all_of(selected_names), maximum_weight_loss_challenge)
LM9_test <- na.omit(LM9_test)
LM9 <- lm(maximum_weight_loss_challenge ~ ., LM9_data)
LM9_step <- stepAIC(LM9, direction =  "backward")
summary(LM9_step)

LM10_data <- dplyr::select(train_ge, all_of(gene_names), all_of(selected_names), maximum_weight_loss_challenge)
LM10_data <- na.omit(LM10_data)
LM10_test <- dplyr::select(test_ge, all_of(gene_names), all_of(selected_names), maximum_weight_loss_challenge)
LM10_test <- na.omit(LM10_test)
LM10 <- lm(maximum_weight_loss_challenge ~ ., LM10_data)
LM10_step <- stepAIC(LM10, direction = "backward")
summary(LM10_step)
control <- trainControl(method = "cv", number = 10)
LM10_steptr <- train(maximum_weight_loss_challenge ~ CXCL9 + CXCR3 + 
                     IDO1 + IL1RN + IRGM1 + MPO + MYD88 + RETNLB + IL17A + eimeria_species_challenge, 
                   data = LM10_data, trControl = control, method = "lm")

ge_models <- list("PCA" = LM9_step, "AIC all" = LM10_step)
aictab(facs_models)
# > ge_models <- list("PCA" = LM9_step, "AIC all" = LM10_step)
# > aictab(facs_models)
# Model selection based on AICc:
#   
#   K   AICc Delta_AICc AICcWt Cum.Wt      LL
# AIC all 14 347.88        0.0   0.93   0.93 -154.69
# PCA      8 353.08        5.2   0.07   1.00 -166.97



LM11_data <- dplyr::select(train_immuno, maximum_weight_loss_challenge, all_of(selected_names), CD4 , Treg , Div_Treg , 
                             Treg17 , Th1 , Div_Th1 , Th17 , Act_CD8 , Div_Act_CD8 , IFNy_CD8, CXCL9 , CXCR3 , 
                             IDO1 , IFNG , IL13 , IL1RN , IRGM1 , MYD88 , RETNLB , TICAM1 , 
                             TNF , IL17A)
LM11_test <- dplyr::select(test_immuno, maximum_weight_loss_challenge, all_of(selected_names), CD4 , Treg , Div_Treg , 
                          Treg17 , Th1 , Div_Th1 , Th17 , Act_CD8 , Div_Act_CD8 , IFNy_CD8, CXCL9 , CXCR3 , 
                          IDO1 , IFNG , IL13 , IL1RN , IRGM1 , MYD88 , RETNLB , TICAM1 , 
                          TNF , IL17A)
LM11 <- lm(maximum_weight_loss_challenge ~ ., LM11_data)
summary(LM11)
# combined LARGE
LM11_prediction <- predict(LM11, newdata = LM11_test)
LM11_test$maximum_weight_loss_challenge  
LM11_prediction  
LM11_error <- LM11_test$maximum_weight_loss_challenge - LM11_prediction  
LM11_RMSE <- sqrt(mean(LM11_error^2))  
LM11_RMSE # 6.2525


# FACS only
LM8_prediction <- predict(LM8_step, newdata = LM8_test)  
LM8_test$maximum_weight_loss_challenge  
LM8_prediction  
LM8_error <- LM8_test$maximum_weight_loss_challenge - LM8_prediction  
LM8_RMSE <- sqrt(mean(LM8_error^2))  
LM8_RMSE # 3.911566

# GE only
LM10_prediction <- predict(LM10_step, newdata = LM10_test)  
LM10_test$maximum_weight_loss_challenge  
LM10_prediction  
LM10_error <- LM10_test$maximum_weight_loss_challenge - LM10_prediction  
LM10_RMSE <- sqrt(mean(LM10_error^2)) 
LM10_RMSE # 6.056478

LM7_prediction <- predict(LM7_step, newdata = LM7_test)  
LM7_test$maximum_weight_loss_challenge  
LM7_prediction  
LM7_error <- LM7_test$maximum_weight_loss_challenge - LM7_prediction  
LM7_RMSE <- sqrt(mean(LM7_error^2)) 
LM7_RMSE # 6.056478

LM9_prediction <- predict(LM9_step, newdata = LM9_test)  
LM9_test$maximum_weight_loss_challenge  
LM9_prediction  
LM9_error <- LM9_test$maximum_weight_loss_challenge - LM9_prediction  
LM9_RMSE <- sqrt(mean(LM9_error^2)) 
LM9_RMSE # 6.056478


# RF1 with all was 4.029967

best_LMs <- list("LM11" = LM11, "LM8" = LM8_step, "LM9" = LM9_step, "LM10" = LM10_step, "LM7" = LM7_step)
aictab(best_LMs)
# Model selection based on AICc:
#   
#   K   AICc Delta_AICc AICcWt Cum.Wt      LL
# LM8  14 347.88       0.00   0.93   0.93 -154.69
# LM7   8 353.08       5.20   0.07   1.00 -166.97
# LM10 14 511.62     163.74   0.00   1.00 -238.77
# LM9  14 544.26     196.38   0.00   1.00 -255.09

# LM8_RMSE   3.911566
# LM10_RMSE  6.056478
# LM7_RMSE   5.037743
# LM9_RMSE   5.8303

# RF1 RMSE   4.052919
# RFD_RMSE   4.137446

LM7_prediction <- as.data.frame(LM7_prediction)
LM8_prediction <- as.data.frame(LM8_prediction)
FACS_real <- as.data.frame(train_facs$maximum_weight_loss_challenge)
names(FACS_real)[names(FACS_real) == "train_facs$maximum_weight_loss_challenge"] <- "original"

LM9_prediction <- as.data.frame(LM9_prediction)
LM10_prediction <- as.data.frame(LM10_prediction)
GE_real <- as.data.frame(LM10_test$maximum_weight_loss_challenge)
names(GE_real)[names(GE_real) == "LM10_test$maximum_weight_loss_challenge"] <- "original"

RFD_prediction <- as.data.frame(def_predict)
RF1_prediction <- as.data.frame(RF1predict)
RF2_prediction <- as.data.frame(RF2predict)
RF_real <- as.data.frame(test_immunorf$maximum_weight_loss_challenge)
names(RF_real)[names(RF_real) == "test_immunorf$maximum_weight_loss_challenge"] <- "original"

PF <- bind_cols(LM7_prediction, LM8_prediction, FACS_real)
PG <- bind_cols(LM9_prediction, LM10_prediction, GE_real)
PR <- bind_cols(RFD_prediction, RF1_prediction, RF2_prediction, RF_real)


ggplot(PF, aes(x = original)) + 
  geom_jitter(width = 0.5,aes(y = original, colour = "original")) + 
  geom_line(aes(y = original, colour = "original")) +
  geom_smooth(aes(y = LM7_prediction, colour = "PCA"), se = F, method = "lm") +
  geom_smooth(aes(y = LM8_prediction, colour = "stepwise"), se = F, method = "lm") +
  geom_jitter(width = 0.5,aes(y = LM7_prediction, colour = "PCA")) +
  geom_jitter(width = 0.5,aes(y = LM8_prediction, colour = "stepwise")) +
  scale_x_continuous(breaks=c(10,20,30,40,50,60,70,80,90,100)) +
  ggtitle("Linear model predictions based on cell counts") +
  theme_bw()

ggplot(PG, aes(x = original)) + 
  geom_jitter(width = 0.5,aes(y = original, colour = "original")) + 
  geom_line(aes(y = original, colour = "original")) +
  geom_smooth(aes(y = LM9_prediction, colour = "PCA"), se = F, method = "lm") +
  geom_smooth(aes(y = LM10_prediction, colour = "stepwise"), se = F, method = "lm") +
  geom_jitter(width = 0.5,aes(y = LM9_prediction, colour = "PCA")) +
  geom_jitter(width = 0.5,aes(y = LM10_prediction, colour = "stepwise")) +
  scale_x_continuous(breaks=c(10,20,30,40,50,60,70,80,90,100)) +
  ggtitle("Linear model predictions based on Gene expression") +
  theme_bw()

ggplot(PR, aes(x = original)) + 
  geom_jitter(width = 0.5, aes(y = original, colour = "original")) + 
  geom_line(aes(y = original, colour = "original")) +
  geom_smooth(aes(y = def_predict, colour = "default"), se = F, method = "lm") +
  geom_smooth(aes(y = RF1predict, colour = "RF1"), se = F, method = "lm") +
  geom_jitter(width = 0.5, aes(y = def_predict, colour = "default")) +
  geom_jitter(width = 0.5, aes(y = RF1predict, colour = "RF1")) +
  scale_x_continuous(breaks=c(10,20,30,40,50,60,70,80,90,100)) +
  ggtitle("Linear model predictions based on random forest") +
  theme_bw()


