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
library(lmerTest)
library(modelr)
library(sjPlot)

lab_long <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/lab_immuno_long.csv"))
wild_long <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/data/HZ19_immuno_long.csv"))
lab <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/Lab_COMPLETE.csv"))
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
wild_long$challenge <- NA
# join
long <- rbind(lab_long, wild_long)
names(long)[names(long) == "EXP_type"] <- "strain_type"
long$EXP_type <- NA
long$EXP_type[long$strain_type == "CLS"] <- "lab"
long$EXP_type[long$strain_type == "WDS"] <- "lab"
long$EXP_type[long$strain_type == "wild"] <- "wild"

# start graphing
# delta by itself
long_delta <- dplyr::select(long, EH_ID, delta, Eim_MC, EXP_type, challenge)
long_delta <- distinct(long_delta)
ggplot((subset(long_delta, !is.na(long_delta$delta))), aes(x = Eim_MC, y = delta, color = Eim_MC)) +
  geom_violin() +
  facet_grid(~EXP_type, drop = T) +
  geom_jitter(stat = "identity") +
  stat_compare_means(method = "wilcox.test", aes(label = ..p.signif..), label.x = 1.5, label.y = 5, comparisons = list(c("infected", "uninfected"))) +
  labs(y = "deltaCT = Mouse - Eimeria", x = "infection status", color = "infection status") +
  theme(axis.text=element_text(size=12, face = "bold"),
        title = element_text(size = 16, face = "bold"),
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        strip.text.y = element_text(size = 14, face = "bold"),
        legend.text = element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold")) +
  ggtitle("infection intensities in wild and wild-derived mice")

# infection(delta) dependent IFNy increase (rewrite for lab where needed)
ggscatter(subset(long, long$EXP_type == "lab" & long$challenge == "E64"), x = "delta", y = "IFNy_CEWE", 
          add = "reg.line", color = "Eim_MC") +
  stat_cor(aes(label= paste(..rr.label.., ..p.label.., sep = "~`,`~")),label.x = -10, label.y = 800) +
  stat_regline_equation(label.x = -10, label.y = 750) +
  facet_grid(EXP_type~Eim_MC) +
  labs(y = "IFN-y (pg/mL)", x = "deltaCT = Mouse - Eimeria", color = "infection status", fill = "infection status") +
  theme(axis.text=element_text(size=12, face = "bold"),
        title = element_text(size = 16, face = "bold"),
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        strip.text.y = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold")) +
  ggtitle("infection intensity effect on IFN-y abundance in the lab")

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

# check ferrisi

ggplot(subset(long, !is.na(long$challenge) & long$Eim_MC == "infected"),aes(x = challenge, y = counts, color = challenge)) +
  geom_boxplot(outlier.shape=NA) + 
  geom_jitter() +
  stat_compare_means(method = "wilcox.test", aes(label = ..p.signif..), size = 6, label.y.npc =0.95, label.x = 1.5) +
  labs(y = "cell counts %", x = "origin") +
  facet_wrap(~pop, scales = "free_y") +
  theme(axis.text=element_text(size=12, face = "bold"),
        title = element_text(size = 16, face = "bold"),
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        strip.text.y = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold")) +
  ggtitle("")
  
# test ferrisi in lab

modIFNyCEWElab <- lmerTest::lmer(IFNy_CEWE~challenge + (1|EXP_type), data = lab)
summary(modIFNyCEWElab)

tab_model(modIFNyCEWElab, 
          file="IFNtable_VS_E64(itercept)lab.html",
          dv.labels=c("CEWE"))
# in wild
wild <- reshape(data = wild_long, direction = "wide", timevar = "pop")
modIFNyCEWEwild <- lmerTest::lmer(IFNy_CEWE~challenge + (1|EXP_type), data = wild)
summary(modIFNyCEWEwild)





# try this
facs.measure.cols <- c("ThCD4p", "TcCD8p", "Th1IFNgp_in_CD4p", "Th17IL17Ap_in_CD4p", 
                       "Tc1IFNgp_in_CD8p", "Treg_Foxp3_in_CD4p", "Dividing_Ki67p_in_Foxp3p", 
                       "RORgtp_in_Foxp3p", "ThCD4p_Foxp3n", "Th1Tbetp_in_CD4pFoxp3n", "Dividing_Ki67p_in_Tbetp", 
                       "Th17RORgp_in_CD4pFoxp3n", "Dividing_Ki67p_in_RORgtp")

mods.l <- lapply(facs.measure.cols, function (x) {
  lm(get(x) ~ (primary * challenge) + Position,
     data=E7)
})
names(mods.l) <- facs.measure.cols
lapply(mods.l, summary)
