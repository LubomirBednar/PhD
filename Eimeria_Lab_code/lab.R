# P3, P4, E7 and E10 (CLS (classic laboratory strains), WDS(wild derived strains))
library(ggplot2)
library(Rmisc)
library(httr)
library(RCurl)
library(tidyverse)
library(ggpubr)

# load in the COMPLETE tables (should contain everything at LABEL level)
P3 <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/P3_112019_Eim_COMPLETE.csv"))
P3$X <- NULL
P4 <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/P4_082020_Eim_COMPLETE.csv"))
P4$X <- NULL
P4$experiment <- "P4"
E7 <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/E7_112018_Eim_COMPLETE.csv"))
E7$X <- NULL
E7$experiment <- "E7"
# E10 coming soon
# unify columns by renaming and selecting
### remake with underscores for clarity
colnames(P3)[10] <- "total_oocysts"
colnames(P4)[30] <- "total_oocysts"
colnames(E7)[21] <- "total_oocysts"

colnames(P3)[19] <- "infection_history"
colnames(P4)[4] <- "infection_history"
colnames(E7)[3] <- "infection_history"

colnames(P3)[15] <- "weight_change"
colnames(P4)[8] <- "weight_change"
colnames(E7)[11] <- "weight_change"

colnames(E7)[1] <- "label"
### add strain to P3 and P4 + other necessary columns for later
P3$Strain <- "SWISS"
P4$Strain <- "SWISS"
P3$EXP_type <- "CLS"
P4$EXP_type <- "CLS"
E7$EXP_type <- "WDS"
P4$Eim_MC[P4$Eim_MC == "TRUE"] <- "pos"
P4$Eim_MC[P4$Eim_MC == "FALSE"] <- "neg"

### pick down to bare minimum of 13 columns
P3 <- dplyr::select(P3, EH_ID, dpi, label, delta, IFNy_CEWE, Eim_MC, total_oocysts, infection_history, weight_change, 
             primary, challenge, OPG, IFNy_CEWE, Strain, EXP_type)
P4 <- dplyr::select(P4, EH_ID, dpi, label, delta, IFNy_CEWE, Eim_MC, total_oocysts, infection_history, weight_change, 
             primary, challenge, OPG, IFNy_CEWE, Strain, EXP_type)
E7 <- dplyr::select(E7, EH_ID, dpi, label, delta, IFNy_CEWE, Eim_MC, total_oocysts, infection_history, weight_change, 
             primary, challenge, OPG, IFNy_CEWE, Strain, EXP_type)

lab <- rbind(P3, P4)
lab <- rbind(lab, E7)

write.csv(lab, "~/GitHub/Eimeria_Lab/data/Experiment_results/Lab_COMPLETE.csv")
# add FACS data
# P3_FACS - need this from Hongwei
# E10_FACS - coming soon
P4_FACS <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/P4_082020_Eim_FACS.csv"))
E7_FACS <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/E7_112018_Eim_FACS.csv"))
colnames(E7_FACS)[3] <- "label"
E7_FACS$Position <- "mLN"

# select only FACS columns
P4_FACS <- dplyr::select(P4_FACS, EH_ID, dpi, label, Position, CD4, Treg, Div_Treg, Treg17, Th1, Div_Th1, Th17, Div_Th17, CD8,
                  Act_CD8, Div_Act_CD8, IFNy_CD4, IFNy_CD8)
E7_FACS <- dplyr::select(E7_FACS, EH_ID, dpi, label, Position, CD4, Treg, Div_Treg, Treg17, Th1, Div_Th1, Th17, Div_Th17, CD8,
                  Act_CD8, Div_Act_CD8, IFNy_CD4, IFNy_CD8)
# use FACS and completes to create end point immuno datasets
P3_immuno <- subset(P3, !is.na(P3$Eim_MC))
# P3_immuno <- merge(P3_immuno, P3_FACS)
write.csv(P3_immuno, "~/GitHub/Eimeria_Lab/data/Experiment_results/P3_112019_Eim_immuno.csv")
P4_immuno <- subset(P4, !is.na(P4$Eim_MC))
P4_immuno <- merge(P4_immuno, P4_FACS)
write.csv(P4_immuno, "~/GitHub/Eimeria_Lab/data/Experiment_results/P4_082020_Eim_immuno.csv")
E7_immuno <- subset(E7, !is.na(E7$Eim_MC))
E7_immuno <- merge(E7_immuno, E7_FACS)
write.csv(E7_immuno, "~/GitHub/Eimeria_Lab/data/Experiment_results/E7_112018_Eim_immuno.csv")
# P3_immuno <- merge(P3_immuno, P3_FACS)

# purely explorative now

lab_immuno <- full_join(P3_immuno, P4_immuno)
lab_immuno <- full_join(lab_immuno, E7_immuno)
lab_immuno <- distinct(lab_immuno)
write.csv(lab_immuno, "~/GitHub/Eimeria_Lab/data/Experiment_results/lab_immuno.csv")

# IFny vs delta exploration
ggscatter(lab_immuno, x = "IFNy_CEWE", y = "delta", color = "Eim_MC", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "kendall",
          xlab = "IFN", ylab = "delta") +
  facet_wrap(~Eim_MC)

ggscatter(subset(lab_immuno, lab_immuno$Eim_MC == "neg"), x = "IFNy_CEWE", y = "delta", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "kendall",
          xlab = "IFN", ylab = "delta")

ggscatter(lab_immuno, x = "delta", y = "IFNy_CEWE", add = "reg.line", color = "Eim_MC") +
  facet_grid(EXP_type~Eim_MC) +
  labs(y = "IFN-y (pg/mL)", x = "deltaCT = Mouse - Eimeria", color = "infection status", fill = "infection status") +
  theme(axis.text=element_text(size=12, face = "bold"),
        title = element_text(size = 16, face = "bold"),
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        strip.text.y = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold")) +
  ylim(0, 510) +
  ggtitle("infection intensity effect on IFN-y abundance")

# transform into long to use cell populations
# remove dpi and just across mLN
lab_immuno$dpi <- NULL
lab_immuno <- subset(lab_immuno, lab_immuno$Position == "mLN")
lab_FACS <- dplyr::select(lab_immuno, EH_ID, delta, CD4, Treg, Div_Treg, Treg17, Th1, Div_Th1,
                          Th17, Div_Th17, CD8, Act_CD8, Div_Act_CD8, IFNy_CD4, IFNy_CD8)
lab_FACS <- distinct(lab_FACS)
lab_FACS_long <- reshape2::melt(lab_FACS,
                                direction = "long",
                                varying = list(names(lab_immuno)[2:14]),
                                v.names = "cell.pop",
                                na.rm = T, value.name = "counts", 
                                id.vars = c("EH_ID", "delta"))
names(lab_FACS_long)[names(lab_FACS_long) == "variable"] <- "pop"
# add needed columns back to FACS_long
necessary <- dplyr::select(lab, EH_ID, Eim_MC, IFNy_CEWE, challenge, EXP_type, OPG)
necessary <- subset(necessary, !is.na(necessary$Eim_MC))
necessary <- distinct(necessary)
lab_FACS_long <- distinct(lab_FACS_long)
lab_FACS_long <- merge(lab_FACS_long, necessary)
lab_immuno_long <- lab_FACS_long 
write.csv(lab_immuno_long, "~/GitHub/Eimeria_Lab/data/Experiment_results/lab_immuno_long.csv")

# cells

ggplot(subset(lab_FACS_long, lab_FACS_long$EXP_type == "WDS"), aes(x = Eim_MC, y = counts, color = Eim_MC)) +
  geom_boxplot(outlier.shape=NA) + 
  geom_point(position=position_jitterdodge()) +
  facet_wrap(~pop, scales = "free_y") +
  stat_compare_means(method = "wilcox.test", aes(label = ..p.signif..), size = 6, label.y.npc =0.95) +
  ggtitle("WDS cell counts")

ggplot(subset(lab_FACS_long, lab_FACS_long$EXP_type == "CLS"), aes(x = Eim_MC, y = counts, color = Eim_MC)) +
  geom_boxplot(outlier.shape=NA) + 
  geom_point(position=position_jitterdodge()) +
  facet_wrap(~pop, scales = "free_y") +
  stat_compare_means(method = "wilcox.test", aes(label = ..p.signif..), size = 6, label.y.npc =0.95) +
  ggtitle("CLS cell counts")

ggplot(lab_FACS_long, aes(x = Eim_MC, y = counts, color = Eim_MC)) +
  geom_boxplot(outlier.shape=NA) + 
  geom_point(position=position_jitterdodge()) +
  facet_wrap(~pop, scales = "free_y") +
  stat_compare_means(method = "wilcox.test", aes(label = ..p.signif..), size = 6, label.y.npc =0.95) +
  ggtitle("Lab mice cell counts")

# combine completes
# order columns for easier navigation
P3 <- P3[,order(colnames(P3))]
E7 <- E7[,order(colnames(E7))]
P4 <- P4[,order(colnames(P4))]

P3$AVG <- NULL
P3$batch <- NULL
E7$average <- NULL
colnames(E7)[1] <- "Eim_MC"
E7$comment <- NULL
E7$CXCR3 <- NULL
E7$IRG6 <- NULL
E7$IFNy_FEC <- NULL
E7$IL.12 <- NULL
E7$label.1 <- NULL
E7$oocyst_1 <- NULL
E7$oocyst_2 <- NULL
E7$oocyst_3 <- NULL
E7$oocyst_4 <- NULL
E7$volume_PBS_mL <- NULL
P3$faeces_weight <- NULL
E7$fecweight <- NULL
P3$Strain <- "SWISS"
E7$IFNy <- NULL
P4$batch <- NULL
P4$EXP_type <- "CLS"
P4$faeces_weight <- NULL
P4$Strain <- "SWISS"
P4 <- subset(P4, !is.na(P4$Eim_MC))

complete <- rbind.fill(P3, P4)


ggplot(subset(lab, !is.na(lab$challenge) & lab$OPG > 0), aes(x = infection_history, y = OPG, color = infection_history)) + 
  geom_boxplot() + 
  geom_jitter()


  summarize(mean_size = mean(OPG, na.rm = TRUE))

ggplot(OPG, aes(x = infection_history, y = mean_size, color = infection_history)) + 
  geom_boxplot()
