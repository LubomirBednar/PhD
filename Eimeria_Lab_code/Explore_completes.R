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
library(psych)
library(PerformanceAnalytics)

# good graphing
# ggplot(HZgraph, 
#        aes(x = MC, y = NE, color = MC)) +
#   geom_boxplot(outlier.shape = NA) +
#   geom_jitter() +
#   facet_wrap("Target", scales = "free") +
#   labs(y="deltaCT = Target - HKG", x = "infected", colour = "infected") +
#   theme(axis.text=element_text(size=12, face = "bold"),
#         title = element_text(size = 16, face = "bold"),
#         axis.title=element_text(size=14,face="bold"),
#         strip.text.x = element_text(size = 14, face = "bold"),
#         legend.text=element_text(size=12, face = "bold"),
#         legend.title = element_text(size = 12, face = "bold"))+
#   ggtitle("Gene expression in wild samples")
#
# good modeling
# summary(lm(formula = IFNy~counts, data = drop_na(subset(wild_immuno_compare, wild_immuno_compare$pop == "IFNy_CD4"&wild_immuno_compare$Eim_MC == "infected"))))

###### 
complete <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/E7_P3_E6_complete.csv"))
# make negative MCs into NAs in a new column
complete$delta_clean <- complete$delta
complete <- mutate(complete, delta_clean = ifelse(Eim_MC == "neg", -30, delta_clean))
complete$dpi <- as.factor(complete$dpi)
complete$X <- NULL
########## make column with E. ferrisi and E. falciformis only
complete$Eimeria_p <- gsub("E64|E139", replacement = "E.ferrisi", complete$primary)
complete$Eimeria_p <- paste(gsub("E88|Eflab", replacement = "E.falciformis", complete$Eimeria_p))

complete$Eimeria_c <- gsub("E64|E139", replacement = "E.ferrisi", complete$challenge)
complete$Eimeria_c <- paste(gsub("E88|Eflab", replacement = "E.falciformis", complete$Eimeria_c))
complete$EXP_type <- "lab"
########## start creating lab_immuno from complete to avoid NAs in wrong places
##### first delta then genes then FACS then para
### 1) delta
lab_delta <- dplyr::select(complete, EH_ID, delta, delta_clean, Eim_MC, EXP_type, challenge)
lab_delta <- dplyr::distinct(lab_delta)
### 2) genes
lab_genes <- dplyr::select(complete, EH_ID, CXCR3, IL.12, IRG6)
lab_genes <- dplyr::distinct(lab_genes)
#lab_genes <- subset(lab_genes, !is.na(lab_genes$CXCR3))
# tranform into long
lab_genes_long <- reshape2::melt(lab_genes,

                        direction = "long",
                        varying = list(names(lab_genes)[2:4]),
                        v.names = "NE",
                        na.rm = T, value.name = "NE", 
                        id.vars = c("EH_ID"))
names(lab_genes_long)[names(lab_genes_long) == "variable"] <- "Target"
### 3) FACS
lab_FACS <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/E7_112018_Eim_FACS.csv"))
lab_FACS$X <- NULL
lab_FACS <- dplyr::select(lab_FACS, EH_ID, infHistory, CD4, Treg, Div_Treg, Treg17, Treg_prop, Th1, Div_Th1,
                           Th17, Div_Th17, CD8, Act_CD8, Div_Act_CD8, IFNy_CD4, IL17A_CD4, IFNy_CD8)
lab_FACS <- dplyr::distinct(lab_FACS)
# transform into long
lab_FACS_long <- reshape2::melt(lab_FACS,
             direction = "long",
             varying = list(names(lab_FACS)[19:34]),
             v.names = "cell.pop",
             na.rm = T, value.name = "counts", 
             id.vars = c("EH_ID", "infHistory"))
names(lab_FACS_long)[names(lab_FACS_long) == "variable"] <- "pop"
# make a summary for comparing with wild (no Pos or Ant difference)
lab_FACS_long <- lab_FACS_long %>% dplyr::group_by(EH_ID, pop, infHistory) %>% dplyr::summarise(counts = mean(counts, na.rm = T))
# write.csv(lab_FACS_long, "~/Eimeria_Lab/data/3_recordingTables/E7_112018_Eim_lab_FACS_long.csv")

# IFNy data ELISA
lab_ELISA_CEWE <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/E7_112018_Eim_CEWE_ELISA.csv"))
lab_ELISA_CEWE$X <- NULL
lab_ELISA_CEWE$label <- NULL

### finally join into lab_immuno

lab_immuno <- merge(lab_genes_long, lab_delta)
lab_immuno <- subset(lab_immuno, !is.na(lab_immuno$delta))
lab_immuno <- merge(lab_immuno, lab_FACS_long)
lab_immuno <- merge(lab_immuno, lab_ELISA_CEWE)

########## create wild_immuno in it's own script
##### add and process so that both tables can be rbound
wild_immuno <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/data/HZ19_immuno_long.csv"))
lab_immuno_compare <- dplyr::select(lab_immuno, EH_ID, delta, Eim_MC, EXP_type, pop, counts, IFNy)
wild_immuno_compare <- dplyr::select(wild_immuno, Mouse_ID, delta, MC, EXP_type, IFNy, pop, counts)
names(lab_immuno_compare)[names(lab_immuno_compare) == "EH_ID"] <- "Mouse_ID"
names(wild_immuno_compare)[names(wild_immuno_compare) == "MC"] <- "Eim_MC"
wild_immuno_compare$Eim_MC[wild_immuno_compare$Eim_MC == "TRUE"] <- "infected"
wild_immuno_compare$Eim_MC[wild_immuno_compare$Eim_MC == "FALSE"] <- "uninfected"
lab_immuno_compare$Eim_MC <- as.character(lab_immuno_compare$Eim_MC)
lab_immuno_compare$Eim_MC[lab_immuno_compare$Eim_MC == "pos"] <- "infected"
lab_immuno_compare$Eim_MC[lab_immuno_compare$Eim_MC == "neg"] <- "uninfected"
wild_immuno_compare <- subset(wild_immuno_compare, !wild_immuno_compare$Mouse_ID == "AA_0791")
lab_immuno_compare <- subset(lab_immuno_compare, !lab_immuno_compare$Mouse_ID == "LM_0294")
lab_immuno_compare$EXP_type[lab_immuno_compare$EXP_type == "lab"] <- "wild-derived"

immuno <- rbind(lab_immuno_compare, wild_immuno_compare)
immuno <- subset(immuno, !immuno$pop == "Treg_prop")
immuno <- subset(immuno, !immuno$Mouse_ID == "LM_0294")
immuno <- subset(immuno, !immuno$Mouse_ID == "AA_0791")

# just infection intensities to start
lab_immuno_delta <- lab_immuno_compare
wild_immuno_delta <- wild_immuno_compare
lab_immuno_delta <- dplyr::select(lab_immuno_delta, Mouse_ID, delta, Eim_MC, EXP_type)
wild_immuno_delta <- dplyr::select(wild_immuno_delta, Mouse_ID, delta, Eim_MC, EXP_type)
lab_immuno_delta <- distinct(lab_immuno_delta)
wild_immuno_delta <- distinct(wild_immuno_delta)

immuno_delta <- rbind(lab_immuno_delta, wild_immuno_delta)
immuno_delta <- distinct(immuno_delta)
immuno_delta <- subset(immuno_delta, !immuno_delta$Mouse_ID == "LM_0294")
immuno_delta <- subset(immuno_delta, !immuno_delta$Mouse_ID == "AA_0791")

# convert Eim_MC to factor otherwise ggplot complains
immuno_delta$Eim_MC <- as.factor(immuno_delta$Eim_MC)
ggplot(subset(immuno_delta, !is.na(immuno_delta$delta)), aes(x = Eim_MC, y = delta, color = Eim_MC)) +
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

# graph out infection intensity effect on IFN-y abundance 
summary(lm(formula = IFNy~counts, data = drop_na(subset(wild_immuno_compare, wild_immuno_compare$pop == "IFNy_CD4"&wild_immuno_compare$Eim_MC == "infected"))))
# ggqqplot(immuno$delta, ylab = "delta")
# ggqqplot(immuno$IFNy, ylab = "IFN")
# qPCRxIFN <- cor.test(immuno$delta, 
#                      immuno$IFNy, 
#                      method = "spearman")
# qPCRxIFN1 <- cor.test(immuno$delta, 
#                      immuno$IFNy, 
#                      method = "pearson")
# qPCRxIFN2 <- cor.test(immuno$delta, 
#                      immuno$IFNy, 
#                      method = "kendall")
# plot(lm(formula = IFNy~counts, data = drop_na(subset(wild_immuno_compare, wild_immuno_compare$pop == "IFNy_CD4" & wild_immuno_compare$Eim_MC == "infected"))))
# summary(lm(formula = IFNy~delta, data = drop_na(subset(wild_immuno_compare, wild_immuno_compare$Eim_MC == "infected"))))

summary(lm(formula = IFNy~delta, data = subset(immuno, immuno$Eim_MC == "infected" & immuno$EXP_type == "lab")))


cor(immuno[,unlist(lapply(immuno, is.numeric))])
pairs.panels(immuno[c(2,6,7)])

infected_immuno <- subset(immuno, immuno$Eim_MC == "infected")
infected_immuno <- distinct(infected_immuno)
chart.Correlation(subset(distinct(infected_immuno[c(2,6,7)])))

uninfected_immuno <- subset(immuno, immuno$Eim_MC == "uninfected")
uninfected_immuno <- distinct(uninfected_immuno)
chart.Correlation(subset(distinct(uninfected_immuno[c(2,6,7)])))
# finihs this at some point   
ggplot(subset(immuno, !is.na(immuno$delta)), aes(x = delta, y = IFNy, colour = Eim_MC)) +
  facet_grid(EXP_type~Eim_MC) +
  stat_summary(fun=mean) + 
  geom_smooth(method='lm')

  



ggscatter(subset(immuno, !is.na(immuno$delta)), x = "delta", y = "IFNy", add = "reg.line", color = "Eim_MC") +
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




# see what wild mice cell populations are relevant during infection 
# (Treg up, Treg 17 down, Th17 down, IFNy CD4 down, IFNy CD8 down)
ggplot(subset(immuno, immuno$EXP_type == "wild"), aes(x = Eim_MC, y = counts, color = Eim_MC)) +
  geom_boxplot(outlier.shape=NA) + 
  geom_point(position=position_jitterdodge()) +
  facet_wrap(~pop, scales = "free_y") +
  stat_compare_means(method = "wilcox.test", aes(label = ..p.signif..), size = 8, label.y.npc =0.95) +
  ggtitle("wild mice cell counts")
# same for lab mice
# (CD4 down, Treg up, Th1 up, Div_Th1 up, Div_Th17 up, Act_CD8 up, Div_Act_CD8 up, IFNy_CD4 up, IFNy_CD8 up)
ggplot(distinct(subset(immuno, immuno$EXP_type == "wild-derived")), aes(x = Eim_MC, y = counts, color = Eim_MC)) +
  geom_boxplot(outlier.shape=NA) + 
  geom_point(position=position_jitterdodge()) +
  facet_wrap(~pop, scales = "free_y") +
  stat_compare_means(method = "wilcox.test", aes(label = ..p.signif..), size = 8, label.y.npc =0.95) +
  ggtitle("lab mice cell counts")
# that one in common population
ggplot(distinct(subset(immuno, immuno$pop == "IFNy_CD4" & !is.na(immuno$delta))), aes(x = EXP_type, y = counts, color = Eim_MC)) +
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



# now let's compare groups of interest between wild and lab
# wild
# IFN  and delta dependent behaviour of these cells wild (IFNy CD4, IFNy CD8, Th17 and Treg17)
ggscatter(drop_na(subset(wild_immuno_compare, wild_immuno_compare$pop == "IFNy_CD4")), x = "counts", y = "IFNy", add = "reg.line", color = "Eim_MC") +
  facet_wrap(~Eim_MC, scales = "free")+
  stat_cor(method = "spearman",label.x =-5, label.y = 600) +
  stat_regline_equation(label.x = -5, label.y = 500) + 
  labs(y = "IFN-y (pg/mL)", x = "cell counts %") +
  theme(axis.text=element_text(size=12, face = "bold"),
        title = element_text(size = 16, face = "bold"),
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        strip.text.y = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold")) +
  ggtitle("wild IFNy CD4 cell count effect on IFN-y abundance")
ggscatter(drop_na(subset(wild_immuno_compare, wild_immuno_compare$pop == "IFNy_CD4")), y = "counts", x = "delta", add = "reg.line", color = "Eim_MC") +
  facet_grid(~Eim_MC, scales = "free")+
  stat_cor(method = "spearman",label.x = -10, label.y = 7) +
  stat_regline_equation(label.x = -10, label.y = 5) + 
  labs(x = "delta", y = "cell counts %") +
  theme(axis.text=element_text(size=12, face = "bold"),
        title = element_text(size = 16, face = "bold"),
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        strip.text.y = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold")) +
  ggtitle("wild delta effect on IFNy_CD4 cell counts")
# IFNy producing CD8s
ggscatter(drop_na(subset(wild_immuno_compare, wild_immuno_compare$pop == "IFNy_CD8")), x = "counts", y = "IFNy", add = "reg.line", color = "Eim_MC") +
  facet_wrap(~Eim_MC, scales = "free")+
  stat_cor(method = "spearman",label.x =-5, label.y = 600) +
  stat_regline_equation(label.x = -5, label.y = 500) + 
  labs(y = "IFN-y (pg/mL)", x = "cell counts %") +
  theme(axis.text=element_text(size=12, face = "bold"),
        title = element_text(size = 16, face = "bold"),
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        strip.text.y = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold")) +
  ggtitle("wild IFNy CD8 cell count effect on IFN-y abundance")
ggscatter(drop_na(subset(wild_immuno_compare, wild_immuno_compare$pop == "IFNy_CD8")), y = "counts", x = "delta", add = "reg.line", color = "Eim_MC") +
  facet_grid(~Eim_MC, scales = "free")+
  stat_cor(method = "spearman",label.x = -10, label.y = 7) +
  stat_regline_equation(label.x = -10, label.y = 5) + 
  labs(x = "delta", y = "cell counts %") +
  theme(axis.text=element_text(size=12, face = "bold"),
        title = element_text(size = 16, face = "bold"),
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        strip.text.y = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold")) +
  ggtitle("wild delta effect on IFNy_CD8 cell counts")
# Th17s 
ggscatter(drop_na(subset(wild_immuno_compare, wild_immuno_compare$pop == "Th17")), y = "counts", x = "IFNy", add = "reg.line", color = "Eim_MC") +
  facet_grid(~Eim_MC, scales = "free")+
  stat_cor(method = "spearman",label.x =0, label.y = 6) +
  stat_regline_equation(label.x = 0, label.y = 5) + 
  labs(x = "IFN-y (pg/mL)", y = "cell counts %") +
  ylim(0,10) +
  theme(axis.text=element_text(size=12, face = "bold"),
        title = element_text(size = 16, face = "bold"),
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        strip.text.y = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold")) +
  ggtitle("wild IFN-y effect on Th17 cell counts")
ggscatter(drop_na(subset(wild_immuno_compare, wild_immuno_compare$pop == "Th17")), y = "counts", x = "delta", add = "reg.line", color = "Eim_MC") +
  facet_grid(~Eim_MC, scales = "free")+
  stat_cor(method = "spearman",label.x = -10, label.y = 7) +
  stat_regline_equation(label.x = -10, label.y = 5) + 
  labs(x = "delta", y = "cell counts %") +
  ylim(0,10) +
  theme(axis.text=element_text(size=12, face = "bold"),
        title = element_text(size = 16, face = "bold"),
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        strip.text.y = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold")) +
  ggtitle("wild delta effect on Th17 cell counts")
# Treg 17s
ggscatter(drop_na(subset(wild_immuno_compare, wild_immuno_compare$pop == "Treg17")), y = "counts", x = "IFNy", add = "reg.line", color = "Eim_MC") +
  facet_grid(~Eim_MC, scales = "free")+
  stat_cor(method = "spearman",label.x =0, label.y = 30) +
  stat_regline_equation(label.x = 0, label.y = 25) + 
  labs(x = "IFN-y (pg/mL)", y = "cell counts %") +
  theme(axis.text=element_text(size=12, face = "bold"),
        title = element_text(size = 16, face = "bold"),
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        strip.text.y = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold")) +
  ggtitle("wild IFN-y effect on Treg17 cell counts")
ggscatter(drop_na(subset(wild_immuno_compare, wild_immuno_compare$pop == "Treg17")), y = "counts", x = "delta", add = "reg.line", color = "Eim_MC") +
  facet_grid(~Eim_MC, scales = "free")+
  stat_cor(method = "spearman",label.x = -10, label.y = 7) +
  stat_regline_equation(label.x = -10, label.y = 5) + 
  labs(x = "delta", y = "cell counts %") +
  theme(axis.text=element_text(size=12, face = "bold"),
        title = element_text(size = 16, face = "bold"),
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        strip.text.y = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold")) +
  ggtitle("wild delta effect on Treg17 cell counts")
# lab
# IFN and delta dependent behaviour of these cells lab (Th1, Div_Th1, Act_CD8, Div_Act_CD8, IFNy_CD4)
ggscatter(drop_na(subset(lab_immuno_compare, lab_immuno_compare$pop == "IFNy_CD4")), y = "counts", x = "IFNy", add = "reg.line", color = "Eim_MC") +
  facet_grid(~Eim_MC, scales = "free")+
  stat_cor(method = "spearman",label.x =0, label.y = 30) +
  stat_regline_equation(label.x = 0, label.y = 25) + 
  labs(x = "IFN-y (pg/mL)", y = "cell counts %") +
  theme(axis.text=element_text(size=12, face = "bold"),
        title = element_text(size = 16, face = "bold"),
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        strip.text.y = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold")) +
  ggtitle("lab IFN-y CD4 cell counts effect on IFN-y abundance")
ggscatter(drop_na(subset(lab_immuno_compare, lab_immuno_compare$pop == "IFNy_CD4")), y = "counts", x = "delta", add = "reg.line", color = "Eim_MC") +
  facet_grid(~Eim_MC, scales = "free")+
  stat_cor(method = "spearman",label.x = -10, label.y = 7) +
  stat_regline_equation(label.x = -10, label.y = 5) + 
  labs(x = "delta", y = "cell counts %") +
  theme(axis.text=element_text(size=12, face = "bold"),
        title = element_text(size = 16, face = "bold"),
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        strip.text.y = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold")) +
  ggtitle("lab delta effect on IFNy CD4 cell counts")
# Th1s
ggscatter(drop_na(subset(lab_immuno_compare, lab_immuno_compare$pop == "Th1")), x = "counts", y = "IFNy", add = "reg.line", color = "Eim_MC") +
  facet_grid(~Eim_MC, scales = "free")+
  stat_cor(method = "spearman",label.x =0, label.y = 30) +
  stat_regline_equation(label.x = 0, label.y = 25) + 
  labs(y = "IFN-y (pg/mL)", x = "cell counts %") +
  theme(axis.text=element_text(size=12, face = "bold"),
        title = element_text(size = 16, face = "bold"),
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        strip.text.y = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold")) +
  ggtitle("lab Th1 cell counts effect on IFN-y abundance")
ggscatter(drop_na(subset(lab_immuno_compare, lab_immuno_compare$pop == "Th1")), y = "counts", x = "delta", add = "reg.line", color = "Eim_MC") +
  facet_grid(~Eim_MC, scales = "free")+
  stat_cor(method = "spearman",label.x = -10, label.y = 7) +
  stat_regline_equation(label.x = -10, label.y = 5) + 
  labs(x = "delta", y = "cell counts %") +
  theme(axis.text=element_text(size=12, face = "bold"),
        title = element_text(size = 16, face = "bold"),
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        strip.text.y = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold")) +
  ggtitle("lab delta effect on Th1 cell counts")
# Div_Th1
ggscatter(drop_na(subset(lab_immuno_compare, lab_immuno_compare$pop == "Div_Th1")), y = "counts", x = "IFNy", add = "reg.line", color = "Eim_MC") +
  facet_grid(~Eim_MC, scales = "free")+
  stat_cor(method = "spearman",label.x =0, label.y = 30) +
  stat_regline_equation(label.x = 0, label.y = 25) + 
  labs(x = "IFN-y (pg/mL)", y = "cell counts %") +
  theme(axis.text=element_text(size=12, face = "bold"),
        title = element_text(size = 16, face = "bold"),
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        strip.text.y = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold")) +
  ggtitle("lab Div_Th1 cell counts effect on IFN-y abundance")
ggscatter(drop_na(subset(lab_immuno_compare, lab_immuno_compare$pop == "Div_Th1")), y = "counts", x = "delta", add = "reg.line", color = "Eim_MC") +
  facet_grid(~Eim_MC, scales = "free")+
  stat_cor(method = "spearman",label.x = -10, label.y = 7) +
  stat_regline_equation(label.x = -10, label.y = 5) + 
  labs(x = "delta", y = "cell counts %") +
  theme(axis.text=element_text(size=12, face = "bold"),
        title = element_text(size = 16, face = "bold"),
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        strip.text.y = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold")) +
  ggtitle("lab delta effect on Div_Th1 cell counts")
# Act_CD8s
ggscatter(drop_na(subset(lab_immuno_compare, lab_immuno_compare$pop == "Act_CD8")), y = "counts", x = "IFNy", add = "reg.line", color = "Eim_MC") +
  facet_grid(~Eim_MC, scales = "free")+
  stat_cor(method = "spearman",label.x =0, label.y = 30) +
  stat_regline_equation(label.x = 0, label.y = 25) + 
  labs(x = "IFN-y (pg/mL)", y = "cell counts %") +
  theme(axis.text=element_text(size=12, face = "bold"),
        title = element_text(size = 16, face = "bold"),
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        strip.text.y = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold")) +
  ggtitle("lab Act_CD8 cell counts effect on IFN-y abundance")
ggscatter(drop_na(subset(lab_immuno_compare, lab_immuno_compare$pop == "Act_CD8")), y = "counts", x = "delta", add = "reg.line", color = "Eim_MC") +
  facet_grid(~Eim_MC, scales = "free")+
  stat_cor(method = "spearman",label.x = -10, label.y = 7) +
  stat_regline_equation(label.x = -10, label.y = 5) + 
  labs(x = "delta", y = "cell counts %") +
  theme(axis.text=element_text(size=12, face = "bold"),
        title = element_text(size = 16, face = "bold"),
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        strip.text.y = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold")) +
  ggtitle("lab delta effect on Act_CD8 cell counts")
# Div_Act_CD8
ggscatter(drop_na(subset(lab_immuno_compare, lab_immuno_compare$pop == "Div_Act_CD8")), y = "counts", x = "IFNy", add = "reg.line", color = "Eim_MC") +
  facet_grid(~Eim_MC, scales = "free")+
  stat_cor(method = "spearman",label.x =0, label.y = 30) +
  stat_regline_equation(label.x = 0, label.y = 25) + 
  labs(x = "IFN-y (pg/mL)", y = "cell counts %") +
  theme(axis.text=element_text(size=12, face = "bold"),
        title = element_text(size = 16, face = "bold"),
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        strip.text.y = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold")) +
  ggtitle("lab Div_Act_CD8 cell counts effect on IFN-y abundance")
ggscatter(drop_na(subset(lab_immuno_compare, lab_immuno_compare$pop == "Div_Act_CD8")), y = "counts", x = "delta", add = "reg.line", color = "Eim_MC") +
  facet_grid(~Eim_MC, scales = "free")+
  stat_cor(method = "spearman",label.x = -10, label.y = 7) +
  stat_regline_equation(label.x = -10, label.y = 5) + 
  labs(x = "delta", y = "cell counts %") +
  theme(axis.text=element_text(size=12, face = "bold"),
        title = element_text(size = 16, face = "bold"),
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        strip.text.y = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold")) +
  ggtitle("lab delta effect on Div_Act_CD8 cell counts")







# looking at dofferences to confirm with IFNy
ggscatter(subset(immuno, immuno$Eim_MC == "infected"), x = "counts", y = "IFNy", add = "reg.line", color = "Eim_MC") +
  facet_grid(pop~EXP_type, scales = "free")+
  stat_cor(method = "spearman",label.x =0, label.y = 450) +
  stat_regline_equation(label.x = 0, label.y = 350) + 
  labs(y = "IFN-y (pg/mL)", x = "cell counts %") +
  theme(axis.text=element_text(size=12, face = "bold"),
        title = element_text(size = 16, face = "bold"),
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        strip.text.y = element_text(size = 8, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold")) +
  ggtitle("")

CD4 <- subset(immuno, immuno$pop == "CD4")
Th1 <- subset(immuno, immuno$pop == "Th1")
Div_Th1 <- subset(immuno, immuno$pop == "Div_Th1")
Div_Th17 <- subset(immuno, immuno$pop == "Div_Th17")
Act_CD8 <- subset(immuno, immuno$pop == "Act_CD8")
Div_Act_CD8 <- subset(immuno, immuno$pop == "Div_Act_CD8")
IFNy_CD4 <- subset(immuno, immuno$pop == "IFNy_CD4")
IFNy_CD8 <- subset(immuno, immuno$pop == "IFNy_CD8")
IL17A_CD4 <- subset(immuno, immuno$pop == "IL17A_CD4")
Th17 <- subset(immuno, immuno$pop == "Th17")

# cytokine producing cells
ggscatter(drop_na(IFNy_CD4), x = "counts", y = "IFNy", add = "reg.line", color = "Eim_MC") +
  facet_grid(Eim_MC~EXP_type, scales = "free")+
  stat_cor(method = "spearman",label.x =0, label.y = 450) +
  stat_regline_equation(label.x = 0, label.y = 350) + 
  labs(y = "IFN-y (pg/mL)", x = "cell counts %") +
  theme(axis.text=element_text(size=12, face = "bold"),
        title = element_text(size = 16, face = "bold"),
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        strip.text.y = element_text(size = 12, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold")) +
  ggtitle("")
ggscatter(drop_na(IFNy_CD8), x = "counts", y = "IFNy", add = "reg.line", color = "Eim_MC") +
  facet_grid(Eim_MC~EXP_type, scales = "free")+
  stat_cor(method = "spearman",label.x =0, label.y = 450) +
  stat_regline_equation(label.x = 0, label.y = 350) + 
  labs(y = "IFN-y (pg/mL)", x = "cell counts %") +
  theme(axis.text=element_text(size=12, face = "bold"),
        title = element_text(size = 16, face = "bold"),
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        strip.text.y = element_text(size = 12, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold")) +
  ggtitle("")
ggscatter(drop_na(IL17A_CD4), y = "counts", x = "IFNy", add = "reg.line", color = "Eim_MC") +
  facet_grid(Eim_MC~EXP_type, scales = "free")+
  stat_cor(method = "spearman",label.x =0, label.y = 10) +
  stat_regline_equation(label.x = 0, label.y = 5) + 
  labs(x = "IFN-y (pg/mL)", y = "cell counts %") +
  theme(axis.text=element_text(size=12, face = "bold"),
        title = element_text(size = 16, face = "bold"),
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        strip.text.y = element_text(size = 12, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold")) +
  ggtitle("")

# Th1s
ggscatter(drop_na(Th1), x = "counts", y = "IFNy", add = "reg.line", color = "Eim_MC") +
  facet_grid(Eim_MC~EXP_type, scales = "free")+
  stat_cor(method = "spearman",label.x =0, label.y = 450) +
  stat_regline_equation(label.x = 0, label.y = 350) + 
  labs(y = "IFN-y (pg/mL)", x = "cell counts %") +
  theme(axis.text=element_text(size=12, face = "bold"),
        title = element_text(size = 16, face = "bold"),
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        strip.text.y = element_text(size = 12, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold")) +
  ggtitle("")
ggscatter(drop_na(Th1), x = "delta", y = "counts", add = "reg.line", color = "Eim_MC") +
  facet_grid(Eim_MC~EXP_type, scales = "free")+
  stat_cor(method = "spearman",label.x =0, label.y = 50) +
  stat_regline_equation(label.x = 0, label.y = 40) + 
  labs(y = "cell counts %", x = "delta") +
  theme(axis.text=element_text(size=12, face = "bold"),
        title = element_text(size = 16, face = "bold"),
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        strip.text.y = element_text(size = 12, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold")) +
  ggtitle("")

ggscatter(drop_na(Div_Th1), x = "counts", y = "IFNy", add = "reg.line", color = "Eim_MC") +
  facet_grid(Eim_MC~EXP_type, scales = "free")+
  stat_cor(method = "spearman",label.x =0, label.y = 450) +
  stat_regline_equation(label.x = 0, label.y = 350) + 
  labs(y = "IFN-y (pg/mL)", x = "cell counts %") +
  theme(axis.text=element_text(size=12, face = "bold"),
        title = element_text(size = 16, face = "bold"),
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        strip.text.y = element_text(size = 12, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold")) +
  ggtitle("")
ggscatter(drop_na(Div_Th1), x = "delta", y = "counts", add = "reg.line", color = "Eim_MC") +
  facet_grid(Eim_MC~EXP_type, scales = "free")+
  stat_cor(method = "spearman",label.x =0, label.y = 50) +
  stat_regline_equation(label.x = 0, label.y = 40) + 
  labs(y = "cell counts %", x = "delta") +
  theme(axis.text=element_text(size=12, face = "bold"),
        title = element_text(size = 16, face = "bold"),
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        strip.text.y = element_text(size = 12, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold")) +
  ggtitle("")

# Th17s
ggscatter(drop_na(Th17), y = "counts", x = "IFNy", add = "reg.line", color = "Eim_MC") +
  facet_grid(Eim_MC~EXP_type, scales = "free")+
  stat_cor(method = "spearman",label.x =0, label.y = 30) +
  stat_regline_equation(label.x = 0, label.y = 20) + 
  labs(x = "IFN-y (pg/mL)", y = "cell counts %") +
  theme(axis.text=element_text(size=12, face = "bold"),
        title = element_text(size = 16, face = "bold"),
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        strip.text.y = element_text(size = 12, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold")) +
  ggtitle("")
ggscatter(drop_na(Th17), x = "delta", y = "counts", add = "reg.line", color = "Eim_MC") +
  facet_grid(Eim_MC~EXP_type, scales = "free")+
  stat_cor(method = "spearman",label.x =0, label.y = 50) +
  stat_regline_equation(label.x = 0, label.y = 40) + 
  labs(y = "cell counts %", x = "delta") +
  theme(axis.text=element_text(size=12, face = "bold"),
        title = element_text(size = 16, face = "bold"),
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        strip.text.y = element_text(size = 12, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold")) +
  ggtitle("")


ggscatter(drop_na(Div_Th17), y = "counts", x = "IFNy", add = "reg.line", color = "Eim_MC") +
  facet_grid(Eim_MC~EXP_type, scales = "free")+
  stat_cor(method = "spearman",label.x =0, label.y = 450) +
  stat_regline_equation(label.x = 0, label.y = 350) + 
  labs(x = "IFN-y (pg/mL)", y = "cell counts %") +
  theme(axis.text=element_text(size=12, face = "bold"),
        title = element_text(size = 16, face = "bold"),
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        strip.text.y = element_text(size = 12, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold")) +
  ggtitle("")






ggscatter(drop_na(Act_CD8), x = "counts", y = "IFNy", add = "reg.line", color = "Eim_MC") +
  facet_grid(Eim_MC~EXP_type, scales = "free")+
  stat_cor(method = "spearman",label.x =0, label.y = 450) +
  stat_regline_equation(label.x = 0, label.y = 350) + 
  labs(y = "IFN-y (pg/mL)", x = "cell counts %") +
  theme(axis.text=element_text(size=12, face = "bold"),
        title = element_text(size = 16, face = "bold"),
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        strip.text.y = element_text(size = 12, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold")) +
  ggtitle("")

ggscatter(drop_na(Div_Act_CD8), x = "counts", y = "IFNy", add = "reg.line", color = "Eim_MC") +
  facet_grid(Eim_MC~EXP_type, scales = "free")+
  stat_cor(method = "spearman",label.x =0, label.y = 450) +
  stat_regline_equation(label.x = 0, label.y = 350) + 
  labs(y = "IFN-y (pg/mL)", x = "cell counts %") +
  theme(axis.text=element_text(size=12, face = "bold"),
        title = element_text(size = 16, face = "bold"),
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        strip.text.y = element_text(size = 12, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold")) +
  ggtitle("")





