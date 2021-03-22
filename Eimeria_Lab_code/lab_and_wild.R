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

lab_long <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/lab_immuno_long.csv")
wild_long <- read.csv("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/data/HZ19_immuno_long.csv")
lab <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/Lab_COMPLETE.csv")
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

my_comparisons <-  list(c("E.falciformis", "E.ferrisi"), 
                        c("E.falciformis", "Uninfected"), 
                        c("E.ferrisi" , "Uninfected"))

my_comparisons1 <- list(c("infected", "uninfected"))

my_comparisons3 <-  list(c("E. falciformis", "non infected"), 
                        c("E. ferrisi" , "non infected"),
                        c("Eimeria sp.", "non infected"))

give.n <- function(x){
  return(c(y = median(x)*0.5, label = length(x)))
}
# start graphing
# delta by itself
long_delta <- dplyr::select(long, EH_ID, delta, Eim_MC, EXP_type, challenge)
long_delta <- distinct(long_delta)
ggplot((subset(long_delta, !is.na(long_delta$delta))), aes(x = Eim_MC, y = delta, color = Eim_MC)) +
  geom_violin() +
  facet_grid(~EXP_type, drop = T) +
  geom_jitter(stat = "identity") +
  stat_compare_means(method = "wilcox.test", aes(label = ..p.signif..), label.x = 1.5, 
                     label.y = 5, comparisons = list(c("infected", "uninfected"))) +
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
complete <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/E7_P3_P4_Eim_complete.csv")
complete$X <- NULL
complete$label.1 <- NULL
complete <- distinct(complete)
complete <- subset(complete, !is.na(IFNy_CEWE))
complete$Position <- replace_na(complete$Position, "mLN")
complete <- subset(complete, Position!="Spleen")
# rewrite MCs from previous experiment as per observations (amp + MC) 334, 335, 340
complete$Eim_MC[3] <- "neg"
complete$Eim_MC[4] <- "neg"
complete$Eim_MC[9] <- "neg"
complete$Eim_MC[6] <- "neg"
complete$Eim_MC[15] <- "neg"
# subset for bad elisa
complete <- data.frame(complete)
complete2 <- complete[-c(70:94),]
complete1 <- complete
# if there are oocysts, it is positive
complete$Eim_MC[complete$OPG > 0] <- "infected"
complete$Eim_MC[complete$Eim_MC == "pos"] <- "infected"
complete$Eim_MC[complete$Eim_MC == "neg"] <- "uninfected"
complete1$Eim_MC[complete1$Eim_MC == "pos"] <- "infected"
complete1$Eim_MC[complete1$Eim_MC == "neg"] <- "uninfected"
# lab IFNy
complete1$Eimeria[complete1$challenge == "E64"] <- "E.ferrisi"
complete1$Eimeria[complete1$challenge == "E88"] <- "E.falciformis"
complete1$Eimeria[complete1$challenge == "UNI"] <- "Uninfected"

complete1$Eimeria[complete1$Eim_MC == "uninfected"] <- "Uninfected"
# run MDS on complete for now




###### add E10 and E11 to complete

E10W <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/E10_112020_Eim_record.csv")
E10W$X <- NULL
E10W <- select(E10W, "EH_ID", "labels", "dpi", "batch", "relative_weight", "challenge_infection", "infection_history")

E11W <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/E11_012021_Eim_record.csv")
E11W$X <- NULL
E11W <- select(E11W, "EH_ID", "labels", "dpi", "batch", "relative_weight", "challenge_infection", "infection_history")


names(long)[names(long) == "Mouse_ID"] <- "EH_ID"





# test and add to ggplot (remove stat cor and regline)
IFN <- select(complete1, EH_ID, delta, IFNy_CEWE, Eimeria)
# graph before reordering for models
ggplot(IFN, aes(x = delta, y = IFNy_CEWE, color = Eimeria)) +
  geom_point(size = 2, show.legend = F) + 
  geom_smooth(method = "lm", show.legend = F) + 
  facet_wrap(~Eimeria) +
  labs(y = "IFN-y (pg/mL)", x = "infection intensity") +
  theme(axis.text=element_text(size=12, face = "bold"),
        title = element_text(size = 16, face = "bold"),
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        strip.text.y = element_text(size = 14, face = "bold")) +
  scale_color_manual(breaks = c("Uninfected", "E.falciformis", "E.ferrisi"),
                     values=c("#009999", "#FF6666", "#339933")) +
  ylim(0, 880) +
  ggtitle("infection intensity effect on IFN-y abundance in lab")

#continue reordering for models
IFN$Eimeria <- factor(IFN$Eimeria, levels = c("Uninfected", "E.falciformis", "E.ferrisi"))

IFNU <- subset(IFN, IFN$Eimeria == "Uninfected")
IFNfer <- subset(IFN, IFN$Eimeria == "E.ferrisi")
IFNfal <- subset(IFN, IFN$Eimeria == "E.falciformis")

# unifected model lab IFN
U <- lm(IFNy_CEWE~delta, IFNU)
summary(U)
tab_model(U, 
          file="IFN_delta_lab_U.html",
          dv.labels=c("IFN-y"))
# ferrisi model
fer <- lm(IFNy_CEWE~delta, IFNfer)
summary(fer)
tab_model(fer, 
          file="IFN_delta_lab_fer.html",
          dv.labels=c("IFN-y"))
# falciformis model
fal <- lm(IFNy_CEWE~ + delta, IFNfal)
summary(fal)
tab_model(fal, 
          file="IFN_delta_lab_fal.html",
          dv.labels=c("IFN-y"))

############################# wild IFNy
HZ19 <- read.csv("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Field/master/data/HZ19_immuno.csv")
HZ19$MC.Eimeria[HZ19$MC.Eimeria == "TRUE"] <- "infected"
HZ19$MC.Eimeria[HZ19$MC.Eimeria == "FALSE"] <- "uninfected"
HZ19 <- select(HZ19, Mouse_ID, delta, IFNy, MC.Eimeria)
HZ19 <- distinct(HZ19)


# models and graph IFN wild
ggplot(subset(HZ19, !is.na(HZ19$IFNy) & !is.na(HZ19$delta)), aes(x = delta, y = IFNy, color = MC.Eimeria)) +
  geom_point(size = 2, show.legend = F) + 
  geom_smooth(method = "lm", show.legend = F) + 
  facet_wrap(~MC.Eimeria) +
  labs(y = "IFN-y (pg/mL)", x = "infection intensity") +
  theme(axis.text=element_text(size=12, face = "bold"),
        title = element_text(size = 16, face = "bold"),
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        strip.text.y = element_text(size = 14, face = "bold")) +
  scale_color_manual(breaks = c("uninfected", "infected"),
                     values=c("#009999", "#FF6666")) +
  ggtitle("infection intensity effect on IFN-y abundance")






# wild IFN
HZ19$MC.Eimeria <- factor(HZ19$MC.Eimeria, levels = c("uninfected", "infected"))

IFNUw <- subset(HZ19, HZ19$MC.Eimeria == "uninfected")
IFNIw <- subset(HZ19, HZ19$MC.Eimeria == "infected")

# unifected model wild
Uw<- lm(IFNy~delta, IFNUw)
summary(Uw)
tab_model(Uw, 
          file="IFN_delta_wild_U.html",
          dv.labels=c("IFN-y"))
# infected model wild
Iw <- lm(IFNy~delta, IFNIw)
summary(Iw)
tab_model(Iw, 
          file="IFN_delta_wild_I.html",
          dv.labels=c("IFN-y"))

########################################## long fix uninfected challenges
long$Eimeria[long$challenge == "E64"] <- "E.ferrisi"
long$Eimeria[long$challenge == "E88"] <- "E.falciformis"
long$Eimeria[long$challenge == "UNI"] <- "Uninfected"

long$Eimeria[long$Eim_MC == "uninfected"] <- "Uninfected"

########################################## cell populations of wild

ggplot(subset(long, long$EXP_type == "wild"),
       aes(x = Eim_MC, y = counts, color = Eim_MC)) +
  geom_boxplot(outlier.shape=NA, show.legend = F) + 
  geom_jitter(size = 3, width = 0.3, show.legend = F) +
  stat_compare_means(comparisons = my_comparisons1,
                     method = "wilcox.test", 
                     aes(label = ..p.signif..), 
                     size = 3, label.y.npc =0.95) +
  facet_wrap(~pop, scales = "free") +
  labs(y = "cell counts %", x = "") + 
  theme(axis.text=element_text(size=12, face = "bold"),
        title = element_text(size = 16, face = "bold"),
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        strip.text.y = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold")) +
  ggtitle("")





############################################ cell populations of lab

ggplot(subset(long, long$EXP_type == "lab"),
       aes(x = Eimeria, y = counts, color = Eimeria)) +
  geom_boxplot(outlier.shape=NA, show.legend = F) + 
  geom_jitter(size = 3, width = 0.3, show.legend = F) +
  stat_compare_means(comparisons = my_comparisons,
                     method = "wilcox.test", 
                     aes(label = ..p.signif..), 
                     size = 3, label.y.npc =0.95) +
  facet_wrap(~pop, scales = "free") +
  labs(y = "cell counts %", x = "") + 
  theme(axis.text=element_text(size=12, face = "bold"),
        title = element_text(size = 16, face = "bold"),
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        strip.text.y = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold")) +
  ggtitle("")



####### compare lab and wild

ggplot(subset(long, !is.na(long$Eim_MC)),
       aes(x = Eim_MC, y = counts, color = Eim_MC)) +
  geom_boxplot(outlier.shape=NA, show.legend = F) + 
  geom_jitter(size = 3, width = 0.3, show.legend = F) +
  stat_compare_means(comparisons = my_comparisons,
                     method = "wilcox.test", 
                     aes(label = ..p.signif..), 
                     size = 3, label.y.npc =0.95) +
  facet_grid(EXP_type~pop, scales = "free") +
  labs(y = "cell counts %", x = "") + 
  theme(axis.text=element_text(size=12, face = "bold"),
        title = element_text(size = 16, face = "bold"),
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        strip.text.y = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold")) +
  ggtitle("")

# remake lab cells flow chart explanation of Th1 response

ggplot(complete1, 
       aes(x = Eimeria, y = CD4, color = Eimeria)) +
  geom_boxplot(outlier.shape=NA, show.legend = F) + 
  geom_jitter(size = 3, width = 0.3, show.legend = F) +
  stat_compare_means(comparisons = my_comparisons,
                                        method = "wilcox.test", 
                                        aes(label = ..p.signif..), 
                                        size = 6, label.y.npc =0.95) +
  labs(y = "cell counts %", x = "") + 
  theme(axis.text=element_text(size=12, face = "bold"),
        title = element_text(size = 16, face = "bold"),
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        strip.text.y = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold")) +
  scale_color_manual(breaks = c("Uninfected", "E.falciformis", "E.ferrisi"),
                     values=c("#009999", "#FF6666", "#339933")) +
  ggtitle("CD4+ in laboratory mice")

ggplot(complete1, 
       aes(x = Eimeria, y = Th1, color = Eimeria)) +
  geom_boxplot(outlier.shape=NA, show.legend = F) + 
  geom_jitter(size = 3, width = 0.3, show.legend = F) +
  stat_compare_means(comparisons = my_comparisons,
                     method = "wilcox.test", 
                     aes(label = ..p.signif..), 
                     size = 6, label.y.npc =0.95) +
  labs(y = "cell counts %", x = "") + 
  theme(axis.text=element_text(size=12, face = "bold"),
        title = element_text(size = 16, face = "bold"),
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        strip.text.y = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold")) +
  scale_color_manual(breaks = c("Uninfected", "E.falciformis", "E.ferrisi"),
                     values=c("#009999", "#FF6666", "#339933")) +
  ggtitle("Th1 in laboratory mice")

ggplot(complete1, 
       aes(x = Eimeria, y = Div_Th1, color = Eimeria)) +
  geom_boxplot(outlier.shape=NA, show.legend = F) + 
  geom_jitter(size = 3, width = 0.3, show.legend = F) +
  stat_compare_means(comparisons = my_comparisons,
                     method = "wilcox.test", 
                     aes(label = ..p.signif..), 
                     size = 6, label.y.npc =0.95) +
  labs(y = "cell counts %", x = "") + 
  theme(axis.text=element_text(size=12, face = "bold"),
        title = element_text(size = 16, face = "bold"),
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        strip.text.y = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold")) +
  scale_color_manual(breaks = c("Uninfected", "E.falciformis", "E.ferrisi"),
                     values=c("#009999", "#FF6666", "#339933")) +
  ggtitle("Actively dividing Th1 in laboratory mice")


# remake lab cells flow chart explanation of CD8 response
ggplot(complete1, 
       aes(x = Eimeria, y = CD8, color = Eimeria)) +
  geom_boxplot(outlier.shape=NA, show.legend = F) + 
  geom_jitter(size = 3, width = 0.3, show.legend = F) +
  stat_compare_means(comparisons = my_comparisons,
                     method = "wilcox.test", 
                     aes(label = ..p.signif..), 
                     size = 6, label.y.npc =0.95) +
  labs(y = "cell counts %", x = "") + 
  theme(axis.text=element_text(size=12, face = "bold"),
        title = element_text(size = 16, face = "bold"),
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        strip.text.y = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold")) +
  scale_color_manual(breaks = c("Uninfected", "E.falciformis", "E.ferrisi"),
                     values=c("#009999", "#FF6666", "#339933")) +
  ggtitle("CD8+ in laboratory mice")

ggplot(complete1, 
       aes(x = Eimeria, y = Act_CD8, color = Eimeria)) +
  geom_boxplot(outlier.shape=NA, show.legend = F) + 
  geom_jitter(size = 3, width = 0.3, show.legend = F) +
  stat_compare_means(comparisons = my_comparisons,
                     method = "wilcox.test", 
                     aes(label = ..p.signif..), 
                     size = 6, label.y.npc =0.95) +
  labs(y = "cell counts %", x = "") + 
  theme(axis.text=element_text(size=12, face = "bold"),
        title = element_text(size = 16, face = "bold"),
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        strip.text.y = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold")) +
  scale_color_manual(breaks = c("Uninfected", "E.falciformis", "E.ferrisi"),
                     values=c("#009999", "#FF6666", "#339933")) +
  ggtitle("T-bet+ (activated) CD8+ in laboratory mice")

ggplot(complete1, 
       aes(x = Eimeria, y =Div_Act_CD8, color = Eimeria)) +
  geom_boxplot(outlier.shape=NA, show.legend = F) + 
  geom_jitter(size = 3, width = 0.3, show.legend = F) +
  stat_compare_means(comparisons = my_comparisons,
                     method = "wilcox.test", 
                     aes(label = ..p.signif..), 
                     size = 6, label.y.npc =0.95) +
  labs(y = "cell counts %", x = "") + 
  theme(axis.text=element_text(size=12, face = "bold"),
        title = element_text(size = 16, face = "bold"),
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        strip.text.y = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold")) +
  scale_color_manual(breaks = c("Uninfected", "E.falciformis", "E.ferrisi"),
                     values=c("#009999", "#FF6666", "#339933")) +
  ggtitle("Dividing T-bet+ (activated) CD8+ in laboratory mice")

# IFNy_CD4 
ggplot(complete1, 
       aes(x = Eimeria, y =IFNy_CD4, color = Eimeria)) +
  geom_boxplot(outlier.shape=NA, show.legend = F) + 
  geom_jitter(size = 3, width = 0.3, show.legend = F) +
  stat_compare_means(comparisons = my_comparisons,
                     method = "wilcox.test", 
                     aes(label = ..p.signif..), 
                     size = 6, label.y.npc =0.95) +
  labs(y = "cell counts %", x = "") + 
  theme(axis.text=element_text(size=12, face = "bold"),
        title = element_text(size = 16, face = "bold"),
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        strip.text.y = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold")) +
  scale_color_manual(breaks = c("Uninfected", "E.falciformis", "E.ferrisi"),
                     values=c("#009999", "#FF6666", "#339933")) +
  ggtitle("IFN-y producing CD4+ cells in laboratory mice")

# IFNy_CD8
ggplot(complete1, 
       aes(x = Eimeria, y =IFNy_CD8, color = Eimeria)) +
  geom_boxplot(outlier.shape=NA, show.legend = F) + 
  geom_jitter(size = 3, width = 0.3, show.legend = F) +
  stat_compare_means(comparisons = my_comparisons,
                     method = "wilcox.test", 
                     aes(label = ..p.signif..), 
                     size = 6, label.y.npc =0.95) +
  labs(y = "cell counts %", x = "") + 
  theme(axis.text=element_text(size=12, face = "bold"),
        title = element_text(size = 16, face = "bold"),
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        strip.text.y = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold")) +
  scale_color_manual(breaks = c("Uninfected", "E.falciformis", "E.ferrisi"),
                     values=c("#009999", "#FF6666", "#339933")) +
  ggtitle("IFN-y producing CD8+ cells in laboratory mice")




# check wild
ggplot(subset(long, long$EXP_type == "wild"), 
       aes(x = factor(Eim_MC, levels = c("infected", "uninfected")), y = counts, color = Eim_MC)) +
  geom_boxplot(outlier.shape=NA) + 
  facet_wrap(~pop) +
  geom_point(position=position_jitterdodge()) +
  stat_compare_means(comparisons = list(c("infected", "uninfected")),
                     method = "wilcox.test", 
                     aes(label = ..p.signif..), 
                     size = 8, label.y.npc =0.95, label.y = 75) +
  labs(y = "cell counts %", x = "") +
  theme(axis.text=element_text(size=12, face = "bold"),
        title = element_text(size = 16, face = "bold"),
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        strip.text.y = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold")) +
  ggtitle("")
# IFNy_CD4, IFNy_CD8, TH17 and Treg17

ggplot(subset(long, long$pop == "IFNy_CD4" & !is.na(long$Eim_MC)), 
       aes(x = factor(Eim_MC, levels = c("infected", "uninfected")), y = counts, color = Eim_MC)) +
  geom_boxplot(outlier.shape=NA, show.legend = F) + 
  facet_wrap(~EXP_type) +
  geom_point(position=position_jitterdodge(), show.legend = F) +
  stat_compare_means(comparisons = list(c("infected", "uninfected")),
                     method = "wilcox.test", 
                     aes(label = ..p.signif..), 
                     size = 6, label.y.npc =0.95, label.y = 15) +
  labs(y = "cell counts %", x = "") +
  theme(axis.text=element_text(size=12, face = "bold"),
        title = element_text(size = 16, face = "bold"),
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        strip.text.y = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold")) +
  ggtitle("IFNy producing CD4+ populations comparison")

ggplot(subset(long, long$pop == "IFNy_CD8" & !is.na(long$Eim_MC)), 
       aes(x = factor(Eim_MC, levels = c("infected", "uninfected")), y = counts, color = Eim_MC)) +
  geom_boxplot(outlier.shape=NA, show.legend = F) + 
  facet_wrap(~EXP_type) +
  geom_point(position=position_jitterdodge(), show.legend = F) +
  stat_compare_means(comparisons = list(c("infected", "uninfected")),
                     method = "wilcox.test", 
                     aes(label = ..p.signif..), 
                     size = 6, label.y.npc =0.95, label.y = 53) +
  labs(y = "cell counts %", x = "") +
  theme(axis.text=element_text(size=12, face = "bold"),
        title = element_text(size = 16, face = "bold"),
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        strip.text.y = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold")) +
  ylim(0, 60) +
  ggtitle("IFNy producing CD8+ populations comparison")

ggplot(subset(long, long$pop == "Th17" & !is.na(long$Eim_MC)), 
       aes(x = factor(Eim_MC, levels = c("infected", "uninfected")), y = counts, color = Eim_MC)) +
  geom_boxplot(outlier.shape=NA, show.legend = F) + 
  facet_wrap(~EXP_type) +
  geom_point(position=position_jitterdodge(), show.legend = F) +
  stat_compare_means(comparisons = list(c("infected", "uninfected")),
                     method = "wilcox.test", 
                     aes(label = ..p.signif..), 
                     size = 6, label.y.npc =0.95, label.y = 10) +
  labs(y = "cell counts %", x = "") +
  theme(axis.text=element_text(size=12, face = "bold"),
        title = element_text(size = 16, face = "bold"),
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        strip.text.y = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold")) +
  ylim(0, 12) +
  ggtitle("Th17 populations comparison")

ggplot(subset(long, long$pop == "Treg17" & !is.na(long$Eim_MC)), 
       aes(x = factor(Eim_MC, levels = c("infected", "uninfected")), y = counts, color = Eim_MC)) +
  geom_boxplot(outlier.shape=NA, show.legend = F) + 
  facet_wrap(~EXP_type) +
  geom_point(position=position_jitterdodge(), show.legend = F) +
  stat_compare_means(comparisons = list(c("infected", "uninfected")),
                     method = "wilcox.test", 
                     aes(label = ..p.signif..), 
                     size = 6, label.y.npc =0.95, label.y = 27) +
  labs(y = "cell counts %", x = "") +
  theme(axis.text=element_text(size=12, face = "bold"),
        title = element_text(size = 16, face = "bold"),
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        strip.text.y = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold")) +
  ylim(0,30) +
  ggtitle("Treg17 populations comparison")

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
# Th17 cells between lab and wild
ggplot(distinct(subset(long, long$pop == "Th17" & !is.na(long$delta))), aes(x = EXP_type, y = counts, color = Eim_MC)) +
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

############################ wild species and genes
HZ18 <- read.csv("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Field/master/data/Gene_expression/HZ18_complete.csv")
HZ18 <- select(HZ18, Mouse_ID, Target, deltaCtMmE_tissue, Eimeria.subspecies, NE, inf)
HZ18$inf[HZ18$inf == "TRUE"] <- "infected"
HZ18$inf[HZ18$inf == "FALSE"] <- "uninfected"
# rename columns
names(HZ18)[names(HZ18) == "deltaCtMmE_tissue"] <- "delta"
# graph
ggplot(HZ18,
       aes(x = Eimeria.subspecies, y = NE, color = Eimeria.subspecies)) +
  geom_boxplot(outlier.shape=NA) + 
  geom_jitter(size = 3, width = 0.3) +
  stat_compare_means(comparisons = my_comparisons3,
                     method = "wilcox.test", 
                     aes(label = ..p.signif..), 
                     size = 3, label.y.npc =0.95) +
  facet_wrap(~Target, scales = "free") +
  labs(y = "Normalised expression", x = "") + 
  theme(axis.text=element_text(size=12, face = "bold"),
        title = element_text(size = 16, face = "bold"),
        axis.title=element_text(size=14,face="bold"),
        axis.title.x = element_blank(),
        axis.text.x=element_blank(),
        legend.position = c(0.8, 0.3),
        strip.text.x = element_text(size = 14, face = "bold"),
        strip.text.y = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold")) +
  scale_color_manual(breaks = c("non infected", "E. falciformis", "E. ferrisi", "Eimeria sp."),
                     values=c("#009999", "#FF6666", "#339933", "#993399")) +
  ggtitle("Gene expression in the wild")
# add intensity for comparison
ggplot(HZ18, aes(x = Eimeria.subspecies, y = delta, color = Eimeria.subspecies)) +
  geom_jitter(size = 2) +
  facet_wrap(~Target) + 
  geom_hline(yintercept=-4, linetype="dashed", color = "red", size = 1) +
  labs(y = "infection intensity", x = "infecting species") +
  theme(axis.text=element_text(size=12, face = "bold"),
        title = element_text(size = 16, face = "bold"),
        axis.title=element_text(size=14,face="bold"),
        axis.title.x = element_blank(),
        axis.text.x=element_blank(),
        legend.position = c(0.8, 0.3),
        strip.text.x = element_text(size = 14, face = "bold"),
        strip.text.y = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold")) +
  scale_color_manual(breaks = c("non infected", "E. falciformis", "E. ferrisi", "Eimeria sp."),
                     values=c("#009999", "#FF6666", "#339933", "#993399")) +
  ggtitle("Gene expression in the wild")

############# lab species and genes
E7 <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/E7_112018_Eim_CEWE_RT-qPCR.csv")
E7$X <- NULL
P3 <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/P3_112019_Eim_CEWE_RTqPCR.csv")
P3$X <- NULL
names(E7)[names(E7) == "Mouse_ID"] <- "EH_ID"

lab_RT <- rbind(E7, P3)
RT_merge <- select(complete1, EH_ID, delta, Eim_MC, Eimeria)
lab_RT <- merge(lab_RT, RT_merge)
RT_long <- gather(lab_RT, Target, NE, CXCR3:IL.12, factor_key=TRUE)

ggplot(RT_long,
         aes(x = Eimeria, y = NE, color = Eimeria)) +
  geom_boxplot(outlier.shape=NA, show.legend = F) + 
  geom_jitter(size = 3, width = 0.3, show.legend = F) +
  stat_compare_means(comparisons = my_comparisons,
                     method = "wilcox.test", 
                     aes(label = ..p.signif..), 
                     size = 3, label.y.npc =0.95) +
  facet_wrap(~Target, scales = "free") +
  labs(y = "Normalised expression", x = "") + 
  theme(axis.text=element_text(size=12, face = "bold"),
        title = element_text(size = 16, face = "bold"),
        axis.title=element_text(size=14,face="bold"),
        axis.text.x=element_text(size = 10, face = "bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        strip.text.y = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold")) +
  ggtitle("Gene expression in the lab")

# well shit...
# wild species and genes
HZRT <- read.csv(text = getURL())



################### toy with lab


CLS <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/CLS_complete.csv")
CLS$Eimeria1[CLS$primary == "E64"] <- "E.ferrisi"
CLS$Eimeria1[CLS$primary == "E139"] <- "E.ferrisi"
CLS$Eimeria1[CLS$primary == "E88"] <- "E.falciformis"
CLS$Eimeria1[CLS$primary == "Eflab"] <- "E.falciformis"
CLS$Eimeria1[CLS$primary == "UNI"] <- "Uninfected"


CLS$Eimeria2[CLS$challenge == "E64"] <- "E.ferrisi"
CLS$Eimeria2[CLS$challenge == "E139"] <- "E.ferrisi"
CLS$Eimeria2[CLS$challenge == "E88"] <- "E.falciformis"
CLS$Eimeria2[CLS$challenge == "UNI"] <- "Uninfected"
# primary
ggplot(subset(CLS, CLS$OPG > 0 & !is.na(CLS$primary)), aes(Eimeria1, OPG, color = Eimeria1)) +
  geom_boxplot() +
  geom_jitter()

ggplot(subset(CLS, !is.na(CLS$primary)), aes( x = dpi, y = relative_weight, color = Eimeria1)) + 
  geom_smooth() + 
  geom_jitter() + 
  facet_wrap(~Eimeria1)

ggplot(subset(CLS, !is.na(CLS$primary)), aes( x = dpi, y = relative_weight, color = primary)) + 
  geom_smooth() + 
  geom_jitter() + 
  facet_wrap(~primary) + 
  ylim(70, 110)

# challenge
ggplot(subset(CLS, CLS$OPG > 0 & !is.na(CLS$challenge)), aes(Eimeria2, OPG, color = Eimeria2)) +
  geom_boxplot() +
  geom_jitter()

ggplot(subset(CLS, !is.na(CLS$challenge)), aes( x = dpi, y = relative_weight, color = Eimeria2)) + 
  geom_smooth() + 
  geom_jitter() + 
  facet_wrap(~Eimeria2)

ggplot(subset(CLS, !is.na(CLS$challenge)), aes( x = dpi, y = relative_weight, color = challenge)) + 
  geom_smooth() + 
  geom_jitter() + 
  facet_wrap(~challenge) + 
  ylim(70, 110)

# batches
ggplot(subset(CLS, CLS$OPG > 0), aes(x = dpi, y = OPG, color = batch)) +
  geom_jitter() +
  facet_wrap(~infection_history)

ggplot(CLS, aes(x = dpi, y = relative_weight, color = batch)) +
  geom_jitter() +
  geom_smooth() +
  facet_wrap(~infection_history) +
  ylim(70,120)
