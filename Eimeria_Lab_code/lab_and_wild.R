library(ggplot2)
library(Rmisc)
library(httr)
library(RCurl)
library(plyr)
library(dplyr)
library(reshape2)
library(tidyverse)
library(ggpubr)
library(stats)

lab_long <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/lab_immuno_long.csv"))
wild_long <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/data/HZ19_immuno_long.csv"))

# rename and delet columns for merging
lab_long$X <- NULL
wild_long$X <- NULL
wild_long$Ct.Eimeria <- NULL
wild_long$Ct.Mus <- NULL
names(wild_long)[names(wild_long) == "Mouse_ID"] <- "EH_ID"
names(wild_long)[names(wild_long) == "MC"] <- "Eim_MC"
lab_long$Eim_MC[lab_long$Eim_MC == "pos"] <- "infected"
lab_long$Eim_MC[lab_long$Eim_MC == "neg"] <- "uninfected"
wild_long$Eim_MC[wild_long$Eim_MC == TRUE] <- "infected"
wild_long$Eim_MC[wild_long$Eim_MC == FALSE] <- "uninfected"
names(wild_long)[names(wild_long) == "IFNy"] <- "IFNy_CEWE"
wild_long$Position <- NULL
lab_long$OPG <- NULL
lab_long$challenge <- NULL
# join
long <- rbind(lab_long, wild_long)
names(long)[names(long) == "EXP_type"] <- "strain_type"
long$EXP_type <- NA
long$EXP_type[long$strain_type == "CLS"] <- "lab"
long$EXP_type[long$strain_type == "WDS"] <- "lab"
long$EXP_type[long$strain_type == "wild"] <- "wild"

# start graphing
# infection(delta) dependent IFNy increase
ggscatter(subset(long, !is.na(long$delta)), x = "delta", y = "IFNy_CEWE", add = "reg.line", color = "Eim_MC") +
  facet_grid(EXP_type~Eim_MC) +
  labs(y = "IFN-y (pg/mL)", x = "deltaCT = Mouse - Eimeria", color = "infection status", fill = "infection status") +
  theme(axis.text=element_text(size=12, face = "bold"),
        title = element_text(size = 16, face = "bold"),
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        strip.text.y = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold")) +
  ylim(0, 1000) +
  ggtitle("infection intensity effect on IFN-y abundance")

# cell populations of wild
ggplot(subset(long, long$EXP_type == "wild" & !is.na(long$pop)), aes(x = Eim_MC, y = counts, color = Eim_MC)) +
  geom_boxplot(outlier.shape=NA) + 
  geom_point(position=position_jitterdodge()) +
  facet_wrap(~pop, scales = "free_y") +
  stat_compare_means(method = "wilcox.test", aes(label = ..p.signif..), size = 6, label.y.npc =0.95) +
  ggtitle("wild mice cell counts")
# cell populations of lab
ggplot(subset(long, long$EXP_type == "lab"), aes(x = Eim_MC, y = counts, color = Eim_MC)) +
  geom_boxplot(outlier.shape=NA) + 
  geom_point(position=position_jitterdodge()) +
  facet_wrap(~pop, scales = "free_y") +
  stat_compare_means(method = "wilcox.test", aes(label = ..p.signif..), size = 6, label.y.npc =0.95) +
  ggtitle("lab mice cell counts")

# after checking the above, only IFNy_CD4, IFNy_CD8, 
# IFNy_CD4 cells between lab and wild
ggplot(distinct(subset(long, long$pop == "IFNy_CD4" & !is.na(long$delta))), aes(x = EXP_type, y = counts, color = Eim_MC)) +
  geom_boxplot(outlier.shape=NA) + 
  geom_point(position=position_jitterdodge()) +
  facet_wrap(~pop, scales = "free_y") +
  stat_compare_means(method = "wilcox.test", aes(label = ..p.signif..), size = 8, label.y.npc =0.95) +
  labs(y = "cell counts %", x = "origin") +
  theme(axis.text=element_text(size=12, face = "bold"),
        title = element_text(size = 16, face = "bold"),
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        strip.text.y = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold")) +
  ggtitle("IFN-y producing CD4+ cells")
# IFNy_CD8 cells between lab and wild
ggplot(distinct(subset(long, long$pop == "IFNy_CD8" & !is.na(long$delta))), aes(x = EXP_type, y = counts, color = Eim_MC)) +
  geom_boxplot(outlier.shape=NA) + 
  geom_point(position=position_jitterdodge()) +
  facet_wrap(~pop, scales = "free_y") +
  stat_compare_means(method = "wilcox.test", aes(label = ..p.signif..), size = 8, label.y.npc =0.95) +
  labs(y = "cell counts %", x = "origin") +
  theme(axis.text=element_text(size=12, face = "bold"),
        title = element_text(size = 16, face = "bold"),
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        strip.text.y = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold")) +
  ggtitle("IFN-y producing CD8+ cells")

# remake lab cells flow chart explanation of Th1 response
ggplot(distinct(subset(long, long$pop == "CD4" & long$EXP_type == "lab")), aes(x = Eim_MC, y = counts, color = Eim_MC)) +
  geom_boxplot(outlier.shape=NA) + 
  geom_point(position=position_jitterdodge()) +
  stat_compare_means(method = "wilcox.test", aes(label = ..p.signif..), size = 8, label.y.npc =0.95) +
  labs(y = "cell counts %", x = "origin") +
  theme(axis.text=element_text(size=12, face = "bold"),
        title = element_text(size = 16, face = "bold"),
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        strip.text.y = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold")) +
  ggtitle("CD4+ in laboratory mice")

ggplot(distinct(subset(long, long$pop == "Th1" & long$EXP_type == "lab")), aes(x = Eim_MC, y = counts, color = Eim_MC)) +
  geom_boxplot(outlier.shape=NA) + 
  geom_point(position=position_jitterdodge()) +
  stat_compare_means(method = "wilcox.test", aes(label = ..p.signif..), size = 8, label.y.npc =0.95) +
  labs(y = "cell counts %", x = "origin") +
  theme(axis.text=element_text(size=12, face = "bold"),
        title = element_text(size = 16, face = "bold"),
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        strip.text.y = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold")) +
  ggtitle("Th1 in laboratory mice")

ggplot(distinct(subset(long, long$pop == "Div_Th1" & long$EXP_type == "lab")), aes(x = Eim_MC, y = counts, color = Eim_MC)) +
  geom_boxplot(outlier.shape=NA) + 
  geom_point(position=position_jitterdodge()) +
  stat_compare_means(method = "wilcox.test", aes(label = ..p.signif..), size = 8, label.y.npc =0.95) +
  labs(y = "cell counts %", x = "origin") +
  theme(axis.text=element_text(size=12, face = "bold"),
        title = element_text(size = 16, face = "bold"),
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        strip.text.y = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold")) +
  ggtitle("Actively dividing Th1 in laboratory mice")
# remake lab cells flow chart explanation of CD8 response
ggplot(distinct(subset(long, long$pop == "CD8" & long$EXP_type == "lab")), aes(x = Eim_MC, y = counts, color = Eim_MC)) +
  geom_boxplot(outlier.shape=NA) + 
  geom_point(position=position_jitterdodge()) +
  stat_compare_means(method = "wilcox.test", aes(label = ..p.signif..), size = 8, label.y.npc =0.95) +
  labs(y = "cell counts %", x = "origin") +
  theme(axis.text=element_text(size=12, face = "bold"),
        title = element_text(size = 16, face = "bold"),
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        strip.text.y = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold")) +
  ggtitle("CD8+ in laboratory mice")

ggplot(distinct(subset(long, long$pop == "Act_CD8" & long$EXP_type == "lab")), aes(x = Eim_MC, y = counts, color = Eim_MC)) +
  geom_boxplot(outlier.shape=NA) + 
  geom_point(position=position_jitterdodge()) +
  stat_compare_means(method = "wilcox.test", aes(label = ..p.signif..), size = 8, label.y.npc =0.95) +
  labs(y = "cell counts %", x = "origin") +
  theme(axis.text=element_text(size=12, face = "bold"),
        title = element_text(size = 16, face = "bold"),
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        strip.text.y = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold")) +
  ggtitle("T-bet+ (activated) CD8+ in laboratory mice")

ggplot(distinct(subset(long, long$pop == "Div_Act_CD8" & long$EXP_type == "lab")), aes(x = Eim_MC, y = counts, color = Eim_MC)) +
  geom_boxplot(outlier.shape=NA) + 
  geom_point(position=position_jitterdodge()) +
  stat_compare_means(method = "wilcox.test", aes(label = ..p.signif..), size = 8, label.y.npc =0.95) +
  labs(y = "cell counts %", x = "origin") +
  theme(axis.text=element_text(size=12, face = "bold"),
        title = element_text(size = 16, face = "bold"),
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        strip.text.y = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold")) +
  ggtitle("Dividing T-bet+ (activated) CD8+ in laboratory mice")

