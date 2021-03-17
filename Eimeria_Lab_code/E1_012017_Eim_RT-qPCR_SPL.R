library(dplyr)
library(ggplot2)
library(lme4)
library(lmerTest)
library(ggpubr)
library(tidyverse)

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

HKG <- read.csv("https://raw.githubusercontent.com/derele/Jan2017Exp/master/raw_xls_data/raw_HKGs.csv")
# got the Cq mean so remove Cq
HKG$Cq <- NULL
names(HKG)[names(HKG) == "Cq.Mean"] <- "Cq"
names(HKG)[names(HKG) == "Sample"] <- "EH_ID"
names(HKG)[names(HKG) == "Gene"] <- "Target"
HKG <- HKG %>% dplyr::group_by(EH_ID, Target) %>% dplyr::summarise(Cq = gm_mean(Cq, na.rm = T))

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
  ggtitle("HKG differences E1 SPL")

# add infection info
E1 <- read.csv("https://raw.githubusercontent.com/derele/Jan2017Exp/master/Experiment_Table_raw_NMRI_Jan2017.csv")
names(E1)[names(E1) == "mouseID"] <- "EH_ID"
E1$dpi.diss <- gsub('7dip', '7dpi', E1$dpi.diss)
E1$dpi.diss <- factor(E1$dpi.diss, levels = c("3dpi", "5dpi", "7dpi", "9dpi", "11dpi"))
E1$inf.strain <- as.factor(E1$inf.strain)
HKG$EH_ID <- toupper(HKG$EH_ID)
HKG <- merge(HKG, E1, by = "EH_ID")

HKG_explore <- subset(HKG, select = c("EH_ID", "Target", "Cq", "inf.strain"))
# check distributions
ggplot(HKG_explore, aes(x=Cq)) +
  geom_density() +
  facet_wrap(~Target)

# check wilcox
HKG_explore$Target <- ordered(HKG_explore$Target,
                              levels = c("CDC42", "Ppia", "Ppip"))

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
         order = c("CDC42", "Ppia", "Ppip"),
         ylab = "Cq", xlab = "Target")

ggline(HKG_explore, x = "Target", y = "Cq", 
       add = c("mean_se", "jitter"), 
       order = c("CDC42", "Ppia", "Ppip"),
       ylab = "Cq", xlab = "Target")
# kruskal test for differences between the HKG targets
kruskal.test(Cq ~ Target, data = HKG_explore) # less that 0.05 =  significantly different
# pairwise wilcox to find out where the differences are
pairwise.wilcox.test(HKG_explore$Cq, HKG_explore$Target,
                     p.adjust.method = "BH")

compare_means(Cq ~ Target,  data = HKG_explore)
#specify comparisons to use
my_comparisons <- list( c("CDC42", "Ppia"), c("CDC42", "Ppip"), c("Ppia", "Ppip") )
ggboxplot(HKG_explore, x = "Target", y = "Cq",
          color = "Target", palette = "jco") + 
  # Add pairwise comparisons p-value as stars
  stat_compare_means(comparisons = my_comparisons, label = "p.signif") + 
  stat_compare_means(label.y = 30)     # Add global p-value


# try with all genes against Cq mean
ggboxplot(HKG, x = "Target", y = "Cq", color = "Target", 
          add = "jitter", legend = "none") +
  rotate_x_text(angle = 45)+
  geom_hline(yintercept = mean(HKG$Cq), linetype = 2)+ # Add horizontal line at base mean
  stat_compare_means(method = "kruskal.test", label.y = 30)+        # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "wilcox",
                     ref.group = ".all.")    

# test variability of HKGs under infection (uninf vs strains)
HKG$inf <- ifelse(HKG$inf.strain == "Uninf", "UNI", "INF")

ggboxplot(HKG, x = "Target", y = "Cq", color = "inf", 
          add = "jitter", legend = "none") +
  rotate_x_text(angle = 45)+
  geom_hline(yintercept = mean(HKG$Cq), linetype = 2)+ # Add horizontal line at base mean
  stat_compare_means(method = "kruskal.test", label.y = 30)+        # Add global kruskal p-value
  stat_compare_means(label = "p.signif", method = "wilcox",
                     ref.group = ".all.")  


# CDC42 <- subset(HKG, HKG$Target == "CDC42")
# compare_means(Cq ~ inf,  data = CDC42)
# 
# Ppia <- subset(HKG, HKG$Target == "Ppia")
# compare_means(Cq ~ inf,  data = Ppia)
# 
# Ppib <- subset(HKG, HKG$Target == "Ppip")
# compare_means(Cq ~ inf,  data = Ppib)
# 
# ggplot(HKG, aes(x = Cq, y = EH_ID, color = inf)) +
#   geom_point() +
#   geom_smooth()


# make wide for working with ranks across 3 variables
HKG_wide <- pivot_wider(data = HKG, names_from = Target, values_from = Cq)

HKG_wide_names <- paste("rank", names(HKG_wide)[5:7], sep="_")
HKG_wide[HKG_wide_names] <-  mutate_each(HKG_wide[5:7],funs(rank(., ties.method="first")))

HKG_rank <- pivot_longer(data = HKG_wide,
                         names_to = "rank",
                         cols = c("rank_Ppip", "rank_CDC42", "rank_Ppia"))
HKG_rank$Ppip <- NULL
HKG_rank$Ppia <- NULL
HKG_rank$CDC42 <- NULL


HKG_rank_Ppia <- subset(HKG_rank, HKG_rank$rank == "rank_Ppia")
compare_means(value ~ inf,  data = HKG_rank_Ppia)

HKG_rank_Ppip <- subset(HKG_rank, HKG_rank$rank == "rank_Ppip")
compare_means(value ~ inf,  data = HKG_rank_Ppip)

HKG_rank_CDC42 <- subset(HKG_rank, HKG_rank$rank == "rank_CDC42")
compare_means(value ~ inf,  data = HKG_rank_CDC42)

my_comparisons2 <- list( c("rank_CDC42", "rank_Ppia"), c("rank_CDC42", "rank_Ppip"), c("rank_Ppia", "rank_Ppip") )

ggplot(HKG_rank,aes(x = rank, y = value, color = EH_ID, group = EH_ID)) +
  geom_point() + 
  geom_line() +
  facet_wrap(~EH_ID) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5,))

kruskal.test(HKG_rank$value, HKG$EH_ID)
pairwise.wilcox.test(HKG_rank$value, HKG_rank$rank,
                     p.adjust.method = "BH")

