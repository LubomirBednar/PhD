library(Rmisc)
library(httr)
library(RCurl)
library(plyr)
library(dplyr)
library(reshape2)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(stats)
library(ggsignif)
library(ggpmisc)

qPCR <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/data/Eimeria_detection/HZ19_CEWE_qPCR.csv"))
qPCR$X <- NULL
# doesn"t exist yet
#RT <- read.csv(text = getURL(""))

ELISA_CEWE <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/data/HZ19_CEWE_ELISA.csv"))
ELISA_CEWE$X <- NULL


FACS <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/data/HZ19_MES_FACS.csv"))
FACS$X <- NULL


immuno <- merge(qPCR, ELISA_CEWE, all = T)
immuno <- merge(immuno, FACS, all = T)
#make into long
SpleenDF <- filter(FACS, Position == "spleen")
mLNDF <- filter(FACS, Position == "mLN")
FACS.long.mln <- melt(mLNDF,
                direction = "long",
                varying = list(names(mLNDF)[1:14]),
                v.names = "cell.pop",
                na.rm = T, value.name = "counts", 
                id.vars = c("Mouse_ID", "Position"))
FACS.long.mln <- na.omit(FACS.long.mln)
names(FACS.long.mln)[names(FACS.long.mln) == "variable"] <- "pop"

immuno.long <- merge(qPCR, ELISA_CEWE)
immuno.long <- merge(immuno.long, FACS.long.mln)

write.csv(immuno.long, "/Users/Luke Bednar/Mouse_Eimeria_Databasing/data/HZ19_immuno_long.csv")

# IFNy vs delta
ggscatter(immuno.long, x = "delta", y = "IFNy", add = "reg.line", color = "MC") +
  facet_grid(~MC, scales = "free")+
  stat_cor(method = "spearman", label.x =-5, label.y = 600) +
  stat_regline_equation(label.x = -5, label.y = 500) + 
  labs(y = "IFN-y (pg/mL)", x = "deltaCT = Mouse - Eimeria") +
  theme(axis.text=element_text(size=12, face = "bold"),
        title = element_text(size = 16, face = "bold"),
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        strip.text.y = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold")) +
  ggtitle("infection intensity effect on IFN-y abundance")

# let's have a look
ggscatter(immuno.long, x = "IFNy", y = "delta", add = "reg.line", color = "MC") +
  facet_wrap(~MC)+
  stat_cor(label.x = 50, label.y = 0) +
  stat_regline_equation(label.x = 50, label.y = 2) + 
  ggtitle("HZ19 infections vs IFNy")

# now FACS look at the delta vs populations
ggscatter(immuno.long, x = "delta", y = "counts", add = "reg.line", color = "MC") +
  facet_grid(~pop)+
  stat_cor(label.x = -20, label.y = 0) +
  stat_regline_equation(label.x = -20, label.y = 5) + 
  ggtitle("HZ19 infections vs IFNy")

ggscatter(subset(immuno.long, immuno.long$MC == "TRUE"), x = "counts", y = "delta", add = "reg.line", color = "MC") +
  facet_wrap(~pop)+
  stat_cor(label.x = -20, label.y = 0) +
  stat_regline_equation(label.x = -20, label.y = 5) + 
  ggtitle("HZ19 infections vs IFNy")

ggscatter(subset(immuno.long, immuno.long$MC == "FALSE"), x = "counts", y = "delta", add = "reg.line", color = "MC") +
  facet_wrap(~pop)+
  stat_cor(label.x = -20, label.y = 0) +
  stat_regline_equation(label.x = -20, label.y = 5) + 
  ggtitle("HZ19 infections vs IFNy")
