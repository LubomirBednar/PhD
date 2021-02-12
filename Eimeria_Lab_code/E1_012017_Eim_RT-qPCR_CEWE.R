library(dplyr)
library(ggplot2)


CDC_1 <- read.csv("./E1_RT-qPCRs/runs_select/CDC24_1.csv")
CDC_1$Sample <- as.character(CDC_1$Sample)
CDC_1 <- CDC_1 %>% dplyr::group_by(Sample, Target) %>% dplyr::summarise(Cq = mean(Cq, na.rm = T))
###################################################################################################
CXCL_1 <- read.csv("./E1_RT-qPCRs/runs_select/CXCL9_1.csv")
CXCL_1$Sample <- as.character(CXCL_1$Sample)

CXCL_2 <- read.csv("./E1_RT-qPCRs/runs_select/CXCL9_2.csv")
CXCL_2$Sample <- as.character(CXCL_2$Sample)
CXCL <- rbind(CXCL_1, CXCL_2)
CXCL <- CXCL %>% dplyr::group_by(Sample, Target) %>% dplyr::summarise(Cq = mean(Cq, na.rm = T))
rm(CXCL_1, CXCL_2)
###################################################################################################
IFN <- read.csv("./E1_RT-qPCRs/runs_select/IFN3_1.csv")
IFN$Sample <- as.character(IFN$Sample)
IFN <- IFN %>% dplyr::group_by(Sample, Target) %>% dplyr::summarise(Cq = mean(Cq, na.rm = T))
###################################################################################################
IL10_1 <- read.csv("./E1_RT-qPCRs/runs_select/IL10_1.csv")
IL10_1$Sample <- as.character(IL10_1$Sample)

IL10_2 <- read.csv("./E1_RT-qPCRs/runs_select/IL10_2.csv")
IL10_2$Sample <- as.character(IL10_2$Sample)
IL10 <- rbind(IL10_1, IL10_2)
IL10 <- IL10 %>% dplyr::group_by(Sample, Target) %>% dplyr::summarise(Cq = mean(Cq, na.rm = T))
rm(IL10_1, IL10_2)
###################################################################################################
IL12_1 <- read.csv("./E1_RT-qPCRs/runs_select/IL12_1.csv")
IL12_1$Sample <- as.character(IL12_1$Sample)

IL12_2 <- read.csv("./E1_RT-qPCRs/runs_select/IL12_2.csv")
IL12_2$Sample <- as.character(IL12_2$Sample)

IL12_3 <- read.csv("./E1_RT-qPCRs/runs_select/IL12_3.csv")
IL12_3$Sample <- as.character(IL12_3$Sample)
IL12 <- rbind(IL12_1, IL12_2)
IL12 <- rbind(IL12, IL12_3)
IL12 <- IL12 %>% dplyr::group_by(Sample, Target) %>% dplyr::summarise(Cq = mean(Cq, na.rm = T))
rm(IL12_1, IL12_2, IL12_3)
###################################################################################################
Ppia_1 <- read.csv("./E1_RT-qPCRs/runs_select/Ppia_1.csv")
Ppia_1$Sample <- as.character(Ppia_1$Sample)

Ppia_2 <- read.csv("./E1_RT-qPCRs/runs_select/Ppia_2.csv")
Ppia_2$Sample <- as.character(Ppia_2$Sample)
Ppia <- rbind(Ppia_1, Ppia_2)
Ppia <- Ppia %>% dplyr::group_by(Sample, Target) %>% dplyr::summarise(Cq = mean(Cq, na.rm = T))
rm(Ppia_1, Ppia_2)
###################################################################################################
IL6_1 <- read.csv("./E1_RT-qPCRs/runs_select/IL6_1.csv")
IL6_1$Sample <- as.character(IL6_1$Sample)

IL6_2 <- read.csv("./E1_RT-qPCRs/runs_select/IL6_2.csv")
IL6_2$Sample <- as.character(IL6_2$Sample)
IL6 <- rbind(IL6_1, IL6_2)
IL6 <- IL6 %>% dplyr::group_by(Sample, Target) %>% dplyr::summarise(Cq = mean(Cq, na.rm = T))
rm(IL6_1, IL6_2)
###################################################################################################
Ppib_1 <- read.csv("./E1_RT-qPCRs/runs_select/Ppib_1.csv")
Ppib_1$Sample <- as.character(Ppib_1$Sample)
Ppib <- Ppib_1
Ppib <- Ppib %>% dplyr::group_by(Sample, Target) %>% dplyr::summarise(Cq = mean(Cq, na.rm = T))
rm(Ppib_1)
###################################################################################################
STAT6_1 <- read.csv("./E1_RT-qPCRs/runs_select/STAT6_1.csv")
STAT6_1$Sample <- as.character(STAT6_1$Sample)

STAT6_2 <- read.csv("./E1_RT-qPCRs/runs_select/STAT6_2.csv")
STAT6_2$Sample <- as.character(STAT6_2$Sample)
STAT6 <- rbind(STAT6_1, STAT6_2)
STAT6 <- STAT6 %>% dplyr::group_by(Sample, Target) %>% dplyr::summarise(Cq = mean(Cq, na.rm = T))
rm(STAT6_1, STAT6_2)
###################################################################################################
TGF_1 <- read.csv("./E1_RT-qPCRs/runs_select/TGF_1.csv")
TGF_1$Sample <- as.character(TGF_1$Sample)

TGF_2 <- read.csv("./E1_RT-qPCRs/runs_select/TGF_2.csv")
TGF_2$Sample <- as.character(TGF_2$Sample)
TGF <- rbind(TGF_1, TGF_2)
TGF <- TGF %>% dplyr::group_by(Sample, Target) %>% dplyr::summarise(Cq = mean(Cq, na.rm = T))
rm(TGF_1, TGF_2)
###################################################################################################

RT <- rbind(CDC_1, CXCL)
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
  labs(y="deltaCT = Target - HKG", x = "deltaCT = Mouse - Eimeria", colour = "infection") +
  theme(title = element_text(size = 16, face = "bold"),
        axis.text=element_text(size=12, face = "bold"),
        axis.title=element_text(size=14,face="bold"),
        legend.position = "none",
        strip.text.x = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_blank()) +
  ggtitle("HKG differences E7")
HKG$EXP <- "E7"

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
E1.wide <- merge(E1, RT.wide, by = "EH_ID")
E1.long <- merge(E1, RT.long, by = "EH_ID")
# order dpi.diss
E1.long$dpi.diss <- factor(E1.long$dpi.diss, levels = c("3dpi", "5dpi", "7dip", "9dpi", "11dpi"))
E1.long$inf.strain <- as.factor(E1.long$inf.strain)

ggplot(subset(E1.long, !is.na(E1.long$NE)), aes(x = as.numeric(dpi.diss), y = NE, color = inf.strain)) +
  geom_point() +
  geom_smooth(se = F) +
  facet_wrap(~Target, scales = "free")

