library(dplyr)
library(ggplot2)
library(lme4)
library(lmerTest)
library(ggpubr)
library(tidyverse)

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

# Read in, set names as characters for later changes and average replicates 
efficiencies <- read.csv("./E1_RT-qPCRs/efficiencies.csv")
# efficiency is calculated by Q (E^(-Ct)) == efficiency factor ^ - Cq


############################################################################################
CDC_1 <- read.csv("./E1_RT-qPCRs/runs_select/CDC24_1.csv")
CDC_1$Sample <- as.character(CDC_1$Sample)

CDC_2 <- read.csv("./E1_RT-qPCRs/runs_select/CDC24_2.csv")
CDC_2$Sample <- as.character(CDC_2$Sample)

CDC <- rbind(CDC_1, CDC_2)

CDC <- CDC %>% dplyr::group_by(Sample, Target) %>% dplyr::summarise(Cq = gm_mean(Cq, na.rm = T))

rm(CDC_1, CDC_2)
###################################################################################################
CXCL_1 <- read.csv("./E1_RT-qPCRs/runs_select/CXCL9_1.csv")
CXCL_1$Sample <- as.character(CXCL_1$Sample)

CXCL_2 <- read.csv("./E1_RT-qPCRs/runs_select/CXCL9_2.csv")
CXCL_2$Sample <- as.character(CXCL_2$Sample)

CXCL_3 <- read.csv("./E1_RT-qPCRs/runs_select/CXCL9_3.csv")
CXCL_3$Sample <- as.character(CXCL_3$Sample)

CXCL <- rbind(CXCL_1, CXCL_2)
CXCL <- rbind(CXCL, CXCL_3)

CXCL <- CXCL %>% dplyr::group_by(Sample, Target) %>% dplyr::summarise(Cq = gm_mean(Cq, na.rm = T))
rm(CXCL_1, CXCL_2, CXCL_3)
###################################################################################################
IFN_1 <- read.csv("./E1_RT-qPCRs/runs_select/IFN3_1.csv")
IFN_1$Sample <- as.character(IFN_1$Sample)

IFN_2 <- read.csv("./E1_RT-qPCRs/runs_select/IFN3_2.csv")
IFN_2$Sample <- as.character(IFN_2$Sample)

IFN_3 <- read.csv("./E1_RT-qPCRs/runs_select/IFN3_3.csv")
IFN_3$Sample <- as.character(IFN_3$Sample)

IFN <- rbind(IFN_1, IFN_2)
IFN <- rbind(IFN, IFN_3)
IFN <- IFN %>% dplyr::group_by(Sample, Target) %>% dplyr::summarise(Cq = gm_mean(Cq, na.rm = T))
rm(IFN_1, IFN_2, IFN_3)
###################################################################################################
IL10_1 <- read.csv("./E1_RT-qPCRs/runs_select/IL10_1.csv")
IL10_1$Sample <- as.character(IL10_1$Sample)

IL10_2 <- read.csv("./E1_RT-qPCRs/runs_select/IL10_2.csv")
IL10_2$Sample <- as.character(IL10_2$Sample)

IL10_3 <- read.csv("./E1_RT-qPCRs/runs_select/IL10_3.csv")
IL10_3$Sample <- as.character(IL10_3$Sample)

IL10_4 <- read.csv("./E1_RT-qPCRs/runs_select/IL10_4.csv")
IL10_4$Sample <- as.character(IL10_4$Sample)

IL10 <- rbind(IL10_1, IL10_2)
IL10 <- rbind(IL10, IL10_3)
IL10 <- rbind(IL10, IL10_4)
IL10 <- IL10 %>% dplyr::group_by(Sample, Target) %>% dplyr::summarise(Cq = gm_mean(Cq, na.rm = T))
rm(IL10_1, IL10_2, IL10_3, IL10_4)
###################################################################################################
IL12_1 <- read.csv("./E1_RT-qPCRs/runs_select/IL12_1.csv")
IL12_1$Sample <- as.character(IL12_1$Sample)

IL12_2 <- read.csv("./E1_RT-qPCRs/runs_select/IL12_2.csv")
IL12_2$Sample <- as.character(IL12_2$Sample)

IL12_3 <- read.csv("./E1_RT-qPCRs/runs_select/IL12_3.csv")
IL12_3$Sample <- as.character(IL12_3$Sample)

IL12_4 <- read.csv("./E1_RT-qPCRs/runs_select/IL12_4.csv")
IL12_4$Sample <- as.character(IL12_4$Sample)

IL12_5 <- read.csv("./E1_RT-qPCRs/runs_select/IL12_5.csv")
IL12_5$Sample <- as.character(IL12_5$Sample)

IL12_6 <- read.csv("./E1_RT-qPCRs/runs_select/IL12_6.csv")
IL12_6$Sample <- as.character(IL12_6$Sample)

IL12 <- rbind(IL12_1, IL12_2)
IL12 <- rbind(IL12, IL12_3)
IL12 <- rbind(IL12, IL12_4)
IL12 <- rbind(IL12, IL12_5)
IL12 <- rbind(IL12, IL12_6)

IL12 <- IL12 %>% dplyr::group_by(Sample, Target) %>% dplyr::summarise(Cq = gm_mean(Cq, na.rm = T))
rm(IL12_1, IL12_2, IL12_3, IL12_4, IL12_5, IL12_6)
###################################################################################################
Ppia_1 <- read.csv("./E1_RT-qPCRs/runs_select/Ppia_1.csv")
Ppia_1$Sample <- as.character(Ppia_1$Sample)

Ppia_2 <- read.csv("./E1_RT-qPCRs/runs_select/Ppia_2.csv")
Ppia_2$Sample <- as.character(Ppia_2$Sample)
Ppia <- rbind(Ppia_1, Ppia_2)
Ppia <- Ppia %>% dplyr::group_by(Sample, Target) %>% dplyr::summarise(Cq = gm_mean(Cq, na.rm = T))
rm(Ppia_1, Ppia_2)
###################################################################################################
IL6_1 <- read.csv("./E1_RT-qPCRs/runs_select/IL6_1.csv")
IL6_1$Sample <- as.character(IL6_1$Sample)

IL6_2 <- read.csv("./E1_RT-qPCRs/runs_select/IL6_2.csv")
IL6_2$Sample <- as.character(IL6_2$Sample)

IL6_3 <- read.csv("./E1_RT-qPCRs/runs_select/IL6_3.csv")
IL6_3$Sample <- as.character(IL6_3$Sample)

IL6 <- rbind(IL6_1, IL6_2)
IL6 <- rbind(IL6, IL6_3)

IL6 <- IL6 %>% dplyr::group_by(Sample, Target) %>% dplyr::summarise(Cq = gm_mean(Cq, na.rm = T))
rm(IL6_1, IL6_2, IL6_3)
###################################################################################################
Ppib_1 <- read.csv("./E1_RT-qPCRs/runs_select/Ppib_1.csv")
Ppib_1$Sample <- as.character(Ppib_1$Sample)

Ppib_2 <- read.csv("./E1_RT-qPCRs/runs_select/Ppib_2.csv")
Ppib_2$Sample <- as.character(Ppib_2$Sample)

Ppib_3 <- read.csv("./E1_RT-qPCRs/runs_select/Ppib_3.csv")
Ppib_3$Sample <- as.character(Ppib_3$Sample)


Ppib <- rbind(Ppib_1, Ppib_2)
Ppib <- rbind(Ppib, Ppib_3)

Ppib <- Ppib %>% dplyr::group_by(Sample, Target) %>% dplyr::summarise(Cq = gm_mean(Cq, na.rm = T))
rm(Ppib_1, Ppib_2, Ppib_3)
###################################################################################################
STAT6_1 <- read.csv("./E1_RT-qPCRs/runs_select/STAT6_1.csv")
STAT6_1$Sample <- as.character(STAT6_1$Sample)

STAT6_2 <- read.csv("./E1_RT-qPCRs/runs_select/STAT6_2.csv")
STAT6_2$Sample <- as.character(STAT6_2$Sample)

STAT6_3 <- read.csv("./E1_RT-qPCRs/runs_select/STAT6_3.csv")
STAT6_3$Sample <- as.character(STAT6_3$Sample)

STAT6 <- rbind(STAT6_1, STAT6_2)
STAT6 <- rbind(STAT6, STAT6_3)

STAT6 <- STAT6 %>% dplyr::group_by(Sample, Target) %>% dplyr::summarise(Cq = gm_mean(Cq, na.rm = T))
rm(STAT6_1, STAT6_2, STAT6_3)
###################################################################################################
TGF_1 <- read.csv("./E1_RT-qPCRs/runs_select/TGF_1.csv")
TGF_1$Sample <- as.character(TGF_1$Sample)

TGF_2 <- read.csv("./E1_RT-qPCRs/runs_select/TGF_2.csv")
TGF_2$Sample <- as.character(TGF_2$Sample)

TGF_3 <- read.csv("./E1_RT-qPCRs/runs_select/TGF_3.csv")
TGF_3$Sample <- as.character(TGF_3$Sample)

TGF <- rbind(TGF_1, TGF_2)
TGF <- rbind(TGF, TGF_3)

TGF <- TGF %>% dplyr::group_by(Sample, Target) %>% dplyr::summarise(Cq = gm_mean(Cq, na.rm = T))
rm(TGF_1, TGF_2, TGF_3)
###################################################################################################

RT <- rbind(CDC, CXCL)
RT <- rbind(RT, IFN)
RT <- rbind(RT, IL10)
RT <- rbind(RT, IL12)
RT <- rbind(RT, IL6)
RT <- rbind(RT, Ppia)
RT <- rbind(RT, Ppib)
RT <- rbind(RT, STAT6)
RT <- rbind(RT, TGF)

RT.long <- RT %>% dplyr::group_by(Sample, Target) %>% dplyr::summarise(Cq = mean(Cq, na.rm = T))
RT.long <- data.frame(RT.long)
RT.wide <- reshape(RT.long[, c("Target", "Sample","Cq")],
                   timevar = "Target", idvar = "Sample", direction = "wide")

# compare graphically becaus I'm just disabled like that
HKG1 <- dplyr::filter(RT.long, Target == "Ppia")
HKG2 <- dplyr::filter(RT.long, Target ==  "Ppib")
HKG3 <- dplyr::filter(RT.long, Target ==  "CDC24")
HKG <- rbind(HKG1, HKG2)
HKG <- rbind(HKG, HKG3)
ggplot(HKG, aes(x = Target, y = Cq, color = Target)) +
  geom_boxplot() +
  geom_jitter() +
  labs(y="raw expression", x = "target", colour = "Target") +
  theme(title = element_text(size = 16, face = "bold"),
        axis.text=element_text(size=12, face = "bold"),
        axis.title=element_text(size=14,face="bold"),
        legend.position = "none",
        strip.text.x = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_blank()) +
  ggtitle("HKG differences E1")
HKG$EXP <- "E1"
# add infection to test for HKG suitability (basic visual exploration)
E1 <- read.csv("https://raw.githubusercontent.com/derele/Jan2017Exp/master/Experiment_Table_raw_NMRI_Jan2017.csv")
names(E1)[names(E1) == "mouseID"] <- "EH_ID"
E1$dpi.diss <- gsub('7dip', '7dpi', E1$dpi.diss)
E1$dpi.diss <- factor(E1$dpi.diss, levels = c("3dpi", "5dpi", "7dpi", "9dpi", "11dpi"))
E1$inf.strain <- as.factor(E1$inf.strain)
HKG$Sample <- sub("^", "LM00", HKG$Sample)
names(HKG)[names(HKG) == "Sample"] <- "EH_ID"
HKG <- merge(HKG, E1, by = "EH_ID")

ggplot(HKG, aes(x = Target, y = Cq, color = Target)) +
  geom_boxplot() +
  geom_jitter() +
  labs(y="raw expression", x = "target", colour = "Target") +
  theme(title = element_text(size = 16, face = "bold"),
        axis.text=element_text(size=12, face = "bold"),
        axis.title=element_text(size=14,face="bold"),
        legend.position = "none",
        strip.text.x = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_blank()) +
  facet_wrap(~inf.strain) +
  ggtitle("HKG differences E1 CEWE")
# test the HKGs 
summary(HKG)
HKG_explore <- subset(HKG, select = c("EH_ID", "Target", "Cq", "inf.strain"))
# check distributions
ggplot(HKG_explore, aes(x=Cq)) +
  geom_density() +
  facet_wrap(~Target)

# check wilcox
HKG_explore$Target <- ordered(HKG_explore$Target,
                         levels = c("CDC24", "Ppia", "Ppib"))

group_by(HKG_explore, Target) %>%
  dplyr::summarise(
    count = n(),
    mean = mean(Cq, na.rm = TRUE),
    sd = sd(Cq, na.rm = TRUE),
    median = median(Cq, na.rm = TRUE),
    IQR = IQR(Cq, na.rm = TRUE)
  )

ggviolin(HKG_explore, x = "Target", y = "Cq", 
          color = "Target", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
          order = c("CDC24", "Ppia", "Ppib"),
          ylab = "Cq", xlab = "Target")

ggline(HKG_explore, x = "Target", y = "Cq", 
       add = c("mean_se", "jitter"), 
       order = c("CDC24", "Ppia", "Ppib"),
       ylab = "Cq", xlab = "Target")
# kruskal test for differences between the HKG targets
kruskal.test(Cq ~ Target, data = HKG_explore) # less that 0.05 =  significantly different
# pairwise wilcox to find out where the differences are
pairwise.wilcox.test(HKG_explore$Cq, HKG_explore$Target,
                     p.adjust.method = "BH")

compare_means(Cq ~ Target,  data = HKG_explore, method = "wilcox")
#specify comparisons to use
my_comparisons <- list( c("CDC24", "Ppia"), c("CDC24", "Ppib"), c("Ppia", "Ppib") )
ggboxplot(HKG_explore, x = "Target", y = "Cq",
          color = "Target", palette = "jco") + 
  # Add pairwise comparisons p-value as stars
  stat_compare_means(comparisons = my_comparisons, label = "p.signif") + 
  stat_compare_means(label.y = 50)     # Add global p-value


# try with all genes against Cq mean
ggboxplot(RT, x = "Target", y = "Cq", color = "Target", 
          add = "jitter", legend = "none") +
  rotate_x_text(angle = 45)+
  geom_hline(yintercept = mean(RT$Cq), linetype = 2)+ # Add horizontal line at base mean
  stat_compare_means(method = "kruskal.test", label.y = 40)+        # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "wilcox",
                     ref.group = ".all.")    

# test variability of HKGs under infection (uninf vs strains)
HKG$inf <- ifelse(HKG$inf.strain == "Uninf", "UNI", "INF")

ggboxplot(HKG, x = "Target", y = "Cq", color = "inf", 
          add = "jitter", legend = "none") +
  rotate_x_text(angle = 45)+
  geom_hline(yintercept = mean(RT$Cq), linetype = 2)+ # Add horizontal line at base mean
  stat_compare_means(method = "kruskal.test", label.y = 40)+        # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "wilcox",
                     ref.group = ".all.")  

# test cmpare means between infected and uninfected in HKGs
CDC42 <- subset(HKG, HKG$Target == "CDC24")
compare_means(Cq ~ inf,  data = CDC42)

Ppia <- subset(HKG, HKG$Target == "Ppia")
compare_means(Cq ~ inf,  data = Ppia)

Ppib <- subset(HKG, HKG$Target == "Ppib")
compare_means(Cq ~ inf,  data = Ppib)


# make wide for working with ranks across 3 variables
HKG_wide <- pivot_wider(data = HKG, names_from = Target, values_from = Cq)

HKG_wide_names <- paste("rank", names(HKG_wide)[6:8], sep="_")
HKG_wide[HKG_wide_names] <-  mutate_each(HKG_wide[6:8],funs(rank(., ties.method="first")))

HKG_rank <- pivot_longer(data = HKG_wide,
                         names_to = "rank",
                         cols = c("rank_Ppib", "rank_CDC24", "rank_Ppia"))
HKG_rank$Ppip <- NULL
HKG_rank$Ppia <- NULL
HKG_rank$CDC42 <- NULL


HKG_rank_Ppia <- subset(HKG_rank, HKG_rank$rank == "rank_Ppia")
compare_means(value ~ inf,  data = HKG_rank_Ppia)

HKG_rank_Ppip <- subset(HKG_rank, HKG_rank$rank == "rank_Ppib")
compare_means(value ~ inf,  data = HKG_rank_Ppip)

HKG_rank_CDC42 <- subset(HKG_rank, HKG_rank$rank == "rank_CDC24")
compare_means(value ~ inf,  data = HKG_rank_CDC42)

ggplot(HKG_rank,aes(x = rank, y = value, color = EH_ID, group = EH_ID)) +
  geom_point() + 
  geom_line() +
  facet_wrap(~EH_ID) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5,))



# name columns to match other data sets (Mouse_ID, Target) + make RT.CT numeric
names(RT.wide)[names(RT.wide) == "Cq.CDC24"] <- "CDC24"
names(RT.wide)[names(RT.wide) == "Cq.CXCL9"] <- "CXCL9"
names(RT.wide)[names(RT.wide) == "Cq.IFNy"] <- "IFNy"
names(RT.wide)[names(RT.wide) == "Cq.IL-10"] <- "IL10"
names(RT.wide)[names(RT.wide) == "Cq.IL-12"] <- "IL12"
names(RT.wide)[names(RT.wide) == "Cq.IL-6"] <- "IL6"
names(RT.wide)[names(RT.wide) == "Cq.Ppia"] <- "Ppia"
names(RT.wide)[names(RT.wide) == "Cq.Ppib"] <- "Ppib"
names(RT.wide)[names(RT.wide) == "Cq.STAT6"] <- "STAT6"
names(RT.wide)[names(RT.wide) == "Cq.TGF-b"] <- "TGFb"

# set ref and target genes
refGenes <- c("Ppia", "Ppib", "CDC24")
targetGenes <- c("CXCL9", "IFNy", "IL10", "IL12", "IL6", "STAT6", "TGFb")
# ### add primer efficiency
# CDC_eff <- 100.2^subset(E1.wide, E1.wide)
# 
# CE.eff <-  eff.factor^(CE.wide[, c(refGenes, targetGenes)] * -1)
# 
# normIDX <- apply(CE.eff[, refGenes], 1, prod)^
#   (1/length(refGenes))
# 
# CE.norm <- CE.eff[, targetGenes] / normIDX
# 
# names(CE.norm) <- gsub("CqM", "NE", names(CE.norm))
# 
# ## fix some very odd outlier numbers
# CE.norm[CE.norm > 1] <- NA
# 
# ## dropping everything but IDs and normalized values... look into SDs,
# ## non-normalized etc... if needed!!
# CE.norm <- cbind(EH_ID=CE.wide[, "EH_ID"], CE.norm)
# 
# ## too lazy to write this more concisely...
# CE.long <- reshape(CE.norm,
#                    direction = "long",
#                    idvar = "EH_ID", ids = EH_ID,
#                    varying = list(grep("^NE\\.", colnames(CE.norm))),
#                    times = grep("^NE\\.", colnames(CE.norm), value=TRUE))
# 
# ## too confused to write this concisely
# rownames(CE.long) <- NULL
# CE.long$time <-  gsub("NE\\.", "", CE.long$time)
# names(CE.long) <- c("EH_ID", "Gene", "NE")
# 
# CE.final <- merge(CE.long, stab, all.y=TRUE)


# calculate ref genes in new column and subtract targets from HKG average, create new columns
require(dplyr)
RT.wide <- RT.wide %>% mutate(refMean = rowMeans(na.rm = TRUE, dplyr::select(RT.wide, refGenes)))
RT.wide <- data.frame(RT.wide)
refMean <- as.numeric(RT.wide$refMean)

# continue with averaging refgenes and subtracting targets from them
RT.wide$CXCL9 <- (RT.wide$refMean - RT.wide$CXCL9)
RT.wide$IFNy <- (RT.wide$refMean - RT.wide$IFNy)
RT.wide$IL10 <- (RT.wide$refMean - RT.wide$IL10)
RT.wide$IL12 <- (RT.wide$refMean - RT.wide$IL12)
RT.wide$IL6 <- (RT.wide$refMean - RT.wide$IL6)
RT.wide$STAT6 <- (RT.wide$refMean - RT.wide$STAT6)
RT.wide$TGFb <- (RT.wide$refMean - RT.wide$TGFb)
RT.wide[RT.wide=="NaN"]<-NA
# remove HKGs after normalization
RT.wide$CDC24 <- NULL
RT.wide$Ppia <- NULL
RT.wide$Ppib <- NULL
RT.wide$refMean <- NULL


RT.long <- reshape(RT.wide, 
                   direction = "long",
                   varying = list(names(RT.wide)[2:8]),
                   v.names = "NE",
                   times = c("CXCL9", "IFNy", "IL10", "IL12", "IL6", "STAT6", "TGFb"),
                   timevar = "Target",
                   idvar = "Sample")
rownames(RT.long) <- NULL
# add mouse names properly and rename Sample to EH_ID
RT.wide$Sample <- sub("^", "LM00", RT.wide$Sample)
names(RT.wide)[names(RT.wide) == "Sample"] <- "EH_ID"
RT.long$Sample <- sub("^", "LM00", RT.long$Sample)
names(RT.long)[names(RT.long) == "Sample"] <- "EH_ID"
# load in mouse history to check how the infections explain expressions
E1 <- read.csv("https://raw.githubusercontent.com/derele/Jan2017Exp/master/Experiment_Table_raw_NMRI_Jan2017.csv")
names(E1)[names(E1) == "mouseID"] <- "EH_ID"
E1$dpi.diss <- gsub('7dip', '7dpi', E1$dpi.diss)
E1.wide <- merge(E1, RT.wide, by = "EH_ID")
E1.long <- merge(E1, RT.long, by = "EH_ID")
# order dpi.diss
E1.long$dpi.diss <- factor(E1.long$dpi.diss, levels = c("3dpi", "5dpi", "7dpi", "9dpi", "11dpi"))
E1.long$inf.strain <- as.factor(E1.long$inf.strain)

ggplot(subset(E1.long, !is.na(E1.long$NE)), aes(x = as.numeric(dpi.diss), y = NE, color = inf.strain)) +
  geom_point() +
  geom_smooth(se = F) +
  facet_wrap(~Target, scales = "free") +
scale_colour_brewer("infection\nisolate", palette = "Dark2")

# add genes as factors to order
E1.long$Target = factor(E1.long$Target, levels=c('CXCL9','IL6','IL10','IL12', "IFNy", "TGFb", "STAT6"))

# plot Spleen
# CytokinesSP <- 
ggplot(subset(E1.long, !is.na(E1.long$NE)), aes(x = dpi.diss, y = NE, color=inf.strain, group = inf.strain)) +
  geom_point() +
  geom_smooth(se=FALSE) +
  facet_wrap(~Target, scales="free", nrow=2)+
  scale_colour_brewer("infection\nisolate", palette = "Dark2") +
  scale_y_continuous("normalized mRNA expression")+
  theme_bw()

write.csv(RT.wide, "./E1_RT-qPCRs/repeats.csv")
