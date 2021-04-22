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
CXCL9_1 <- read.csv("https://raw.githubusercontent.com/derele/Jan2017Exp/master/RT-qPCRs_2020/results/MES/E1_CXCL9_MES1.csv")
CXCL9_1$Sample <- as.character(CXCL9_1$Sample)

CXCL9_2 <- read.csv("https://raw.githubusercontent.com/derele/Jan2017Exp/master/RT-qPCRs_2020/results/MES/E1_CXCL9_MES2.csv")
CXCL9_2$Sample <- as.character(CXCL9_2$Sample)
CXCL9 <- rbind(CXCL9_1, CXCL9_2)

CXCL9 <- CXCL9 %>% dplyr::group_by(Sample, Target) %>% dplyr::summarise(Cq = gm_mean(Cq, na.rm = T))

rm(CXCL9_1, CXCL9_2)
###################################################################################################
IFNy_1 <- read.csv("https://raw.githubusercontent.com/derele/Jan2017Exp/master/RT-qPCRs_2020/results/MES/E1_IFNy_MES1.csv")
IFNy_1$Sample <- as.character(IFNy_1$Sample)

IFNy_2 <- read.csv("https://raw.githubusercontent.com/derele/Jan2017Exp/master/RT-qPCRs_2020/results/MES/E1_IFNy_MES2.csv")
IFNy_2$Sample <- as.character(IFNy_2$Sample)
IFNy <- rbind(IFNy_1, IFNy_2)

IFNy <- IFNy %>% dplyr::group_by(Sample, Target) %>% dplyr::summarise(Cq = gm_mean(Cq, na.rm = T))
rm(IFNy_1, IFNy_2)
###################################################################################################
IL10_1 <- read.csv("https://raw.githubusercontent.com/derele/Jan2017Exp/master/RT-qPCRs_2020/results/MES/E1_IL-10_MES1.csv")
IL10_1$Sample <- as.character(IL10_1$Sample)

IL10 <- IL10_1 %>% dplyr::group_by(Sample, Target) %>% dplyr::summarise(Cq = gm_mean(Cq, na.rm = T))
rm(IL10_1)
###################################################################################################
IL12_1 <- read.csv("https://raw.githubusercontent.com/derele/Jan2017Exp/master/RT-qPCRs_2020/results/MES/E1_IL-12_MES1.csv")
IL12_1$Sample <- as.character(IL12_1$Sample)

IL12_2 <- read.csv("https://raw.githubusercontent.com/derele/Jan2017Exp/master/RT-qPCRs_2020/results/MES/E1_IL-12_MES2.csv")
IL12_2$Sample <- as.character(IL12_2$Sample)

IL12 <- rbind(IL12_1, IL12_2)

IL12 <- IL12 %>% dplyr::group_by(Sample, Target) %>% dplyr::summarise(Cq = gm_mean(Cq, na.rm = T))
rm(IL12_1, IL12_2)
###################################################################################################
IL6_1 <- read.csv("https://raw.githubusercontent.com/derele/Jan2017Exp/master/RT-qPCRs_2020/results/MES/E1_IL-6_MES1.csv")
IL6_1$Sample <- as.character(IL6_1$Sample)

IL6_2 <- read.csv("https://raw.githubusercontent.com/derele/Jan2017Exp/master/RT-qPCRs_2020/results/MES/E1_IL-6_MES2.csv")
IL6_2$Sample <- as.character(IL6_2$Sample)

IL6 <- rbind(IL6_1, IL6_2)

IL6 <- IL6 %>% dplyr::group_by(Sample, Target) %>% dplyr::summarise(Cq = gm_mean(Cq, na.rm = T))
rm(IL6_1, IL6_2)
###################################################################################################
Ppia_1 <- read.csv("https://raw.githubusercontent.com/derele/Jan2017Exp/master/RT-qPCRs_2020/results/MES/E1_Ppia_MES1.csv")
Ppia_1$Sample <- as.character(Ppia_1$Sample)


Ppia <- Ppia_1 %>% dplyr::group_by(Sample, Target) %>% dplyr::summarise(Cq = gm_mean(Cq, na.rm = T))
rm(Ppia_1)
###################################################################################################
Ppip_1 <- read.csv("https://raw.githubusercontent.com/derele/Jan2017Exp/master/RT-qPCRs_2020/results/MES/E1_Ppip_MES1.csv")
Ppip_1$Sample <- as.character(Ppip_1$Sample)

Ppip_2 <- read.csv("https://raw.githubusercontent.com/derele/Jan2017Exp/master/RT-qPCRs_2020/results/MES/E1_Ppip_MES2.csv")
Ppip_2$Sample <- as.character(Ppip_2$Sample)

Ppip <- rbind(Ppip_1, Ppip_2)

Ppip <- Ppip %>% dplyr::group_by(Sample, Target) %>% dplyr::summarise(Cq = gm_mean(Cq, na.rm = T))
rm(Ppip_1, Ppip_2)
###################################################################################################
CDC42_1 <- read.csv("https://raw.githubusercontent.com/derele/Jan2017Exp/master/RT-qPCRs_2020/results/MES/E1_CDC42_MES1.csv")
CDC42_1$Sample <- as.character(CDC42_1$Sample)

CDC42 <- CDC42_1 %>% dplyr::group_by(Sample, Target) %>% dplyr::summarise(Cq = gm_mean(Cq, na.rm = T))
rm(CDC42_1)
###################################################################################################
STAT6_1 <- read.csv("https://raw.githubusercontent.com/derele/Jan2017Exp/master/RT-qPCRs_2020/results/MES/E1_STAT6_MES1.csv")
STAT6_1$Sample <- as.character(STAT6_1$Sample)

STAT6_2 <- read.csv("https://raw.githubusercontent.com/derele/Jan2017Exp/master/RT-qPCRs_2020/results/MES/E1_STAT6_MES2.csv")
STAT6_2$Sample <- as.character(STAT6_2$Sample)

STAT6 <- rbind(STAT6_1, STAT6_2)

STAT6 <- STAT6%>% dplyr::group_by(Sample, Target) %>% dplyr::summarise(Cq = gm_mean(Cq, na.rm = T))
rm(STAT6_1, STAT6_2)
###################################################################################################
TGF_1 <- read.csv("https://raw.githubusercontent.com/derele/Jan2017Exp/master/RT-qPCRs_2020/results/MES/E1_TGFb_MES1.csv")
TGF_1$Sample <- as.character(TGF_1$Sample)

TGF_2 <- read.csv("https://raw.githubusercontent.com/derele/Jan2017Exp/master/RT-qPCRs_2020/results/MES/E1_TGFb_MES2.csv")
TGF_2$Sample <- as.character(TGF_2$Sample)

TGF <- rbind(TGF_1, TGF_2)

TGF <- TGF %>% dplyr::group_by(Sample, Target) %>% dplyr::summarise(Cq = gm_mean(Cq, na.rm = T))
rm(TGF_1, TGF_2)
###################################################################################################

RT <- rbind(CXCL9, IFNy)
RT <- rbind(RT, IL10)
RT <- rbind(RT, IL12)
RT <- rbind(RT, IL6)
RT <- rbind(RT, Ppip)
RT <- rbind(RT, STAT6)
RT <- rbind(RT, TGF)
RT <- rbind(RT, Ppia)
RT <- rbind(RT, CDC42)

RT.long <- RT %>% dplyr::group_by(Sample, Target) %>% dplyr::summarise(Cq = mean(Cq, na.rm = T))
RT.long <- data.frame(RT.long)
RT.wide <- reshape(RT.long[, c("Target", "Sample","Cq")],
                   timevar = "Target", idvar = "Sample", direction = "wide")

# compare graphically becaus I'm just disabled like that
# HKG1 <- dplyr::filter(RT.long, Target == "Ppia")
# HKG2 <- dplyr::filter(RT.long, Target ==  "Ppip")
# HKG3 <- dplyr::filter(RT.long, Target ==  "CDC42")
# HKG <- rbind(HKG1, HKG2)
# HKG <- rbind(HKG, HKG3)
# ggplot(HKG, aes(x = Target, y = Cq, color = Target)) +
#   geom_boxplot() +
#   geom_jitter() +
#   labs(y="raw expression", x = "target", colour = "Target") +
#   theme(title = element_text(size = 16, face = "bold"),
#         axis.text=element_text(size=12, face = "bold"),
#         axis.title=element_text(size=14,face="bold"),
#         legend.position = "none",
#         strip.text.x = element_text(size = 14, face = "bold"),
#         legend.text=element_text(size=12, face = "bold"),
#         legend.title = element_blank()) +
#   ggtitle("HKG differences E1")
# HKG$EXP <- "E1"
# # add infection to test for HKG suitability (basic visual exploration)
# E1 <- read.csv("https://raw.githubusercontent.com/derele/Jan2017Exp/master/Experiment_Table_raw_NMRI_Jan2017.csv")
# names(E1)[names(E1) == "mouseID"] <- "EH_ID"
# E1$dpi.diss <- gsub('7dip', '7dpi', E1$dpi.diss)
# E1$dpi.diss <- factor(E1$dpi.diss, levels = c("3dpi", "5dpi", "7dpi", "9dpi", "11dpi"))
# E1$inf.strain <- as.factor(E1$inf.strain)
# HKG$Sample <- sub("^", "LM00", HKG$Sample)
# names(HKG)[names(HKG) == "Sample"] <- "EH_ID"
# HKG <- merge(HKG, E1, by = "EH_ID")
# 
# ggplot(HKG, aes(x = Target, y = Cq, color = Target)) +
#   geom_boxplot() +
#   geom_jitter() +
#   labs(y="raw expression", x = "target", colour = "Target") +
#   theme(title = element_text(size = 16, face = "bold"),
#         axis.text=element_text(size=12, face = "bold"),
#         axis.title=element_text(size=14,face="bold"),
#         legend.position = "none",
#         strip.text.x = element_text(size = 14, face = "bold"),
#         legend.text=element_text(size=12, face = "bold"),
#         legend.title = element_blank()) +
#   facet_wrap(~inf.strain) +
#   ggtitle("HKG differences E1 CEWE")
# # test the HKGs 
# summary(HKG)
# HKG_explore <- subset(HKG, select = c("EH_ID", "Target", "Cq", "inf.strain"))
# # check distributions
# ggplot(HKG_explore, aes(x=Cq)) +
#   geom_density() +
#   facet_wrap(~Target)
# 
# # check wilcox
# HKG_explore$Target <- ordered(HKG_explore$Target,
#                               levels = c("CDC42", "Ppia", "Ppip"))
# 
# group_by(HKG_explore, Target) %>%
#   dplyr::summarise(
#     count = n(),
#     mean = mean(Cq, na.rm = TRUE),
#     sd = sd(Cq, na.rm = TRUE),
#     median = median(Cq, na.rm = TRUE),
#     IQR = IQR(Cq, na.rm = TRUE)
#   )
# 
# ggviolin(HKG_explore, x = "Target", y = "Cq", 
#          color = "Target", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
#          order = c("CDC42", "Ppia", "Ppip"),
#          ylab = "Cq", xlab = "Target")
# 
# ggline(HKG_explore, x = "Target", y = "Cq", 
#        add = c("mean_se", "jitter"), 
#        order = c("CDC42", "Ppia", "Ppip"),
#        ylab = "Cq", xlab = "Target")
# # kruskal test for differences between the HKG targets
# kruskal.test(Cq ~ Target, data = HKG_explore) # less that 0.05 =  significantly different
# # pairwise wilcox to find out where the differences are
# pairwise.wilcox.test(HKG_explore$Cq, HKG_explore$Target,
#                      p.adjust.method = "BH")
# 
# compare_means(Cq ~ Target,  data = HKG_explore, method = "wilcox")
# #specify comparisons to use
# my_comparisons <- list( c("CDC42", "Ppia"), c("CDC42", "Ppip"), c("Ppia", "Ppip") )
# ggboxplot(HKG_explore, x = "Target", y = "Cq",
#           color = "Target", palette = "jco") + 
#   facet_wrap(~inf.strain) +
#   # Add pairwise comparisons p-value as stars
#   stat_compare_means(comparisons = my_comparisons, label = "p.signif") + 
#   stat_compare_means(label.y = 50)     # Add global p-value
# 
# 
# # try with all genes against Cq mean
# ggboxplot(RT, x = "Target", y = "Cq", color = "Target", 
#           add = "jitter", legend = "none") +
#   rotate_x_text(angle = 45)+
#   geom_hline(yintercept = mean(RT$Cq), linetype = 2)+ # Add horizontal line at base mean
#   stat_compare_means(method = "kruskal.test", label.y = 40)+        # Add global annova p-value
#   stat_compare_means(label = "p.signif", method = "wilcox",
#                      ref.group = ".all.")    
# 
# # test variability of HKGs under infection (uninf vs strains)
# HKG$inf <- ifelse(HKG$inf.strain == "Uninf", "UNI", "INF")
# 
# ggboxplot(HKG, x = "Target", y = "Cq", color = "inf", 
#           add = "jitter", legend = "none") +
#   rotate_x_text(angle = 45)+
#   geom_hline(yintercept = mean(RT$Cq), linetype = 2)+ # Add horizontal line at base mean
#   stat_compare_means(method = "kruskal.test", label.y = 40)+        # Add global annova p-value
#   stat_compare_means(label = "p.signif", method = "wilcox",
#                      ref.group = ".all.")  
# 
# # test cmpare means between infected and uninfected in HKGs
# CDC42 <- subset(HKG, HKG$Target == "CDC42")
# compare_means(Cq ~ inf,  data = CDC42)
# 
# Ppia <- subset(HKG, HKG$Target == "Ppia")
# compare_means(Cq ~ inf,  data = Ppia)
# 
# Ppip <- subset(HKG, HKG$Target == "Ppip")
# compare_means(Cq ~ inf,  data = Ppip)
# 
# 
# # make wide for working with ranks across 3 variables
# HKG_wide <- pivot_wider(data = HKG, names_from = Target, values_from = Cq)
# 
# HKG_wide_names <- paste("rank", names(HKG_wide)[6:8], sep="_")
# HKG_wide[HKG_wide_names] <-  mutate_each(HKG_wide[6:8],funs(rank(., ties.method="first")))
# 
# HKG_rank <- pivot_longer(data = HKG_wide,
#                          names_to = "rank",
#                          cols = c("rank_Ppip", "rank_CDC42", "rank_Ppia"))
# HKG_rank$Ppip <- NULL
# HKG_rank$Ppia <- NULL
# HKG_rank$CDC42 <- NULL
# 
# 
# HKG_rank_Ppia <- subset(HKG_rank, HKG_rank$rank == "rank_Ppia")
# compare_means(value ~ inf,  data = HKG_rank_Ppia)
# 
# HKG_rank_Ppip <- subset(HKG_rank, HKG_rank$rank == "rank_Ppip")
# compare_means(value ~ inf,  data = HKG_rank_Ppip)
# 
# HKG_rank_CDC42 <- subset(HKG_rank, HKG_rank$rank == "rank_CDC42")
# compare_means(value ~ inf,  data = HKG_rank_CDC42)
# 
# ggplot(HKG_rank,aes(x = rank, y = value, color = EH_ID, group = EH_ID)) +
#   geom_point() + 
#   geom_line() +
#   facet_wrap(~EH_ID) +
#   theme(axis.text.x = element_text(angle = 45, vjust = 0.5,))
# 


# name columns to match other data sets (Mouse_ID, Target) + make RT.CT numeric
names(RT.wide)[names(RT.wide) == "Cq.CDC42"] <- "CDC42"
names(RT.wide)[names(RT.wide) == "Cq.CXCL9"] <- "CXCL9"
names(RT.wide)[names(RT.wide) == "Cq.IFN-y"] <- "IFNy"
names(RT.wide)[names(RT.wide) == "Cq.IL-10"] <- "IL10"
names(RT.wide)[names(RT.wide) == "Cq.IL-12"] <- "IL12"
names(RT.wide)[names(RT.wide) == "Cq.IL-6"] <- "IL6"
names(RT.wide)[names(RT.wide) == "Cq.Ppia"] <- "Ppia"
names(RT.wide)[names(RT.wide) == "Cq.Ppip"] <- "Ppip"
names(RT.wide)[names(RT.wide) == "Cq.STAT6"] <- "STAT6"
names(RT.wide)[names(RT.wide) == "Cq.TGF-b"] <- "TGFb"

# set ref and target genes
refGenes <- c("Ppia", "Ppip", "CDC42")
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
RT.wide$CDC42 <- NULL
RT.wide$Ppia <- NULL
RT.wide$Ppip <- NULL
RT.wide$refMean <- NULL

RT.long <- pivot_longer(RT.wide, cols = c("CXCL9", "IFNy", "IL10", "IL12", "IL6", "STAT6", "TGFb"))
names(RT.long)[names(RT.long) == "name"] <- "Target"
names(RT.long)[names(RT.long) == "value"] <- "NE"

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
  facet_wrap(~Target, scales = "free_y") +
  scale_colour_brewer("infection\nisolate", palette = "Dark2") +
  scale_x_continuous(labels=c("dpi3", "dpi5", "dpi7", "dpi9", "dpi11"))

# add genes as factors to order
E1.long$Target = factor(E1.long$Target, levels=c('CXCL9','IL6','IL10','IL12', "IFNy", "TGFb", "STAT6"))

## Now contrasting against negative control
# l3v3l setting introduces only NAs
E1.long$inf.strain = factor(E1.long$inf.strain, levels(E1.long$inf.strain)[c(4,1:3)])

modCXCL9.cu <- lmer(NE~inf.strain + (1|dpi.diss), data=subset(E1.long, E1.long$Target%in%"CXCL9"))
summary(modCXCL9.cu)

modIL10.cu <- lmer(NE~inf.strain + (1|dpi.diss), data=subset(E1.long, E1.long$Target%in%"IL10"))
summary(modIL10.cu)

modIL12.cu <- lmer(NE~inf.strain + (1|dpi.diss), data=subset(E1.long, E1.long$Target%in%"IL12"))
summary(modIL12.cu)

modIL6.cu <- lmer(NE~inf.strain + (1|dpi.diss), data=subset(E1.long, E1.long$Target%in%"IL6"))
summary(modIL6.cu)

modIFNG.cu <- lmer(NE~inf.strain + (1|dpi.diss), data=subset(E1.long, E1.long$Target%in%"IFNy"))
summary(modIFNG.cu)

modSTAT6.cu <- lmer(NE~inf.strain  +(1|dpi.diss), data=subset(E1.long, E1.long$Target%in%"STAT6"))
summary(modSTAT6.cu)

modTGFB.cu <- lmer(NE~inf.strain + (1|dpi.diss), data=subset(E1.long, E1.long$Target%in%"TGFb"))
summary(modTGFB.cu)

ggplot(fortify(modTGFB.cu), aes(inf.strain , NE, color=dpi.diss)) +
  stat_summary(fun.data=mean_se, geom="pointrange") +
  stat_summary(aes(y=.fitted), fun.y=mean, geom="line")


tab_model(modCXCL9.cu, modIL10.cu, modIL12.cu, modIL6.cu,
          modIFNG.cu, modSTAT6.cu, modTGFB.cu,
          file="CEtable_VS_non_infected(itercept).html",
          dv.labels=c("CXCL9", "IL10", "IL12", "IL6",
                      "IFNG", "STAT6", "TGFB"))
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