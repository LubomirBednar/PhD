library(httr)
library(RCurl)
library(Rmisc)
library(tidyverse)
library(purrr)
library(ggplot2)
library(reshape2)
library(naniar)

# load in runs
RT1 <- "https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/data/Gene_expression/HZ16-17_RT-qPCRs/HZ16-17_RT-qPCR1.CSV"
RT1 <- read.csv(text = getURL(RT1))

RT2 <- "https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/data/Gene_expression/HZ16-17_RT-qPCRs/HZ16-17_RT-qPCR2.CSV"
RT2 <- read.csv(text = getURL(RT2))

RT3 <- "https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/data/Gene_expression/HZ16-17_RT-qPCRs/HZ16-17_RT-qPCR3.CSV"
RT3 <- read.csv(text = getURL(RT3))

RT4 <- "https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/data/Gene_expression/HZ16-17_RT-qPCRs/HZ16-17_RT-qPCR4.CSV"
RT4 <- read.csv(text = getURL(RT4))

RT5 <- "https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/data/Gene_expression/HZ16-17_RT-qPCRs/HZ16-17_RT-qPCR5.CSV"
RT5 <- read.csv(text = getURL(RT5))

RT6 <- "https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/data/Gene_expression/HZ16-17_RT-qPCRs/HZ16-17_RT-qPCR6.CSV"
RT6 <- read.csv(text = getURL(RT6))

RT7 <- "https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/data/Gene_expression/HZ16-17_RT-qPCRs/HZ16-17_RT-qPCR7.CSV"
RT7 <- read.csv(text = getURL(RT7))

RT8 <- "https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/data/Gene_expression/HZ16-17_RT-qPCRs/HZ16-17_RT-qPCR8.CSV"
RT8 <- read.csv(text = getURL(RT8))

RT9 <- "https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/data/Gene_expression/HZ16-17_RT-qPCRs/HZ16-17_RT-qPCR9.CSV"
RT9 <- read.csv(text = getURL(RT9))
  
RT10 <- "https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/data/Gene_expression/HZ16-17_RT-qPCRs/HZ16-17_RT-qPCR10.CSV" 
RT10 <- read.csv(text = getURL(RT10))  
  
RT11 <- "https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/data/Gene_expression/HZ16-17_RT-qPCRs/HZ16-17_RT-qPCR11.CSV" 
RT11 <- read.csv(text = getURL(RT11)) 

# clean and merge it all
RT1$Ct.Mean.SYBR <- NULL
RT1$Ct.Dev..SYBR <- NULL
RT4$Ct.Mean.SYBR <- NULL
RT4$Ct.Dev..SYBR <- NULL
RT1$Amount.SYBR <- NULL
RT2$Amount.SYBR <- NULL
RT3$Amount.SYBR <- NULL
RT4$Amount.SYBR <- NULL
RT5$Amount.SYBR <- NULL
RT6$Amount.SYBR <- NULL
RT7$Amount.SYBR <- NULL
RT8$Amount.SYBR <- NULL
RT9$Amount.SYBR <- NULL
RT10$Amount.SYBR <- NULL
RT11$Amount.SYBR <- NULL

RT1$Pos <- NULL
RT2$Pos <- NULL
RT3$Pos <- NULL
RT4$Pos <- NULL
RT5$Pos <- NULL
RT6$Pos <- NULL
RT7$Pos <- NULL
RT8$Pos <- NULL
RT9$Pos <- NULL
RT10$Pos <- NULL
RT11$Pos <- NULL

RT <- rbind(RT1, RT2)
RT <- rbind(RT, RT3)
RT <- rbind(RT, RT4)
RT <- rbind(RT, RT5)
RT <- rbind(RT, RT6)
RT <- rbind(RT, RT7)
RT <- rbind(RT, RT8)
RT <- rbind(RT, RT9)
RT <- rbind(RT, RT10)
RT <- rbind(RT, RT11)

#remove negative controls
RT <- RT[!grepl("IRG6A", RT$Name),]
RT <- RT[!grepl("CXCR3", RT$Name),]
RT <- RT[!grepl("IL-12rb1", RT$Name),]
RT <- RT[!grepl("B-actin", RT$Name),]
RT <- RT[!grepl("GAPDH", RT$Name),]

# name columns to match other data sets (Mouse_ID, Target) + make RT.CT numeric
names(RT)[names(RT) == "Target.SYBR"] <- "Target"
names(RT)[names(RT) == "Ct.SYBR"] <- "RT.Ct"
names(RT)[names(RT) == "Name"] <- "Mouse_ID"
RT$RT.Ct <- as.numeric(as.character(RT$RT.Ct))

# remove NAs, calculate averages + save long
RT <- na.omit(RT)
RT.long <- RT %>% dplyr::group_by(Mouse_ID, Target) %>% dplyr::summarise(RT.Ct = mean(RT.Ct), .drop = FALSE)
RT.long <- data.frame(RT.long)
RT.wide <- reshape(RT.long[, c("Target", "Mouse_ID","RT.Ct")],
                   timevar = "Target", idvar = "Mouse_ID", direction = "wide")
# set ref and target genes
refGenes <- c("RT.Ct.B-actin", "RT.Ct.GAPDH")
targetGenes <- c("RT.Ct.CXCR3", "RT.Ct.IL-12", "RT.Ct.IRG6")
# compare graphically becaus I'm just disabled like that
HKG1 <- dplyr::filter(RT.long, Target == "B-actin")
HKG2 <- dplyr::filter(RT.long, Target ==  "GAPDH")
HKG <- rbind(HKG1, HKG2)
ggplot(HKG, aes(x = Target, y = RT.Ct, color = Target)) +
  geom_point() +
  geom_boxplot() +
  stat_compare_means(aes(label = ..p.signif..), size = 8, label.y.npc =1) +
  labs(y="deltaCT = Target - HKG", x = "deltaCT = Mouse - Eimeria", colour = "infection") +
  theme(title = element_text(size = 16, face = "bold"),
        axis.text=element_text(size=12, face = "bold"),
        axis.title=element_text(size=14,face="bold"),
        legend.position = "none",
        strip.text.x = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_blank()) +
  ggtitle("HKG differences HZ16-17")
HKG$EXP <- "HZ16-17"
write.csv(HKG, "C:/Users/Luke Bednar/Eimeria_Lab/data/3_recordingTables/HKG_HZ16-17.csv")


# calculate ref genes in new column and subtract targets from HKG average, create new columns
require(dplyr)
RT.wide <- RT.wide %>% mutate(refMean = rowMeans(na.rm = TRUE, dplyr::select(RT.wide, refGenes)))
RT.wide <- data.frame(RT.wide)
refMean <- as.numeric(RT.wide$refMean)

# continue with averaging refgenes and subtracting targets from them
RT.wide$CXCR3 <- (RT.wide$refMean - RT.wide$RT.Ct.CXCR3)
RT.wide$IRG6 <- (RT.wide$refMean - RT.wide$RT.Ct.IRG6)
RT.wide$IL.12 <- (RT.wide$refMean - RT.wide$RT.Ct.IL.12)
# remove non normalized expressions
RT.wide$RT.Ct.CXCR3 <- NULL
RT.wide$RT.Ct.IRG6 <- NULL
RT.wide$RT.Ct.IL.12 <- NULL
RT.wide$RT.Ct.GAPDH <- NULL
RT.wide$RT.Ct.B.actin <- NULL
RT.wide$refMean <- NULL

RT.long <- melt(RT.wide, id.vars = "Mouse_ID")
names(RT.long)[names(RT.long) == "variable"] <- "Target"
names(RT.long)[names(RT.long) == "value"] <- "NE"
# correct names + add sample column
# RT <- separate(RT, c("Mouse_ID"), into = c("AA", "Number"))
# RT$Mouse_ID <- sub("^", "AA_", RT$Number )
# RT$AA <- NULL
# RT$Number <- NULL
# calculate averages
# RT <- RT %>% dplyr::group_by(Mouse_ID, Target) %>% dplyr::summarise(RT.Ct = mean(RT.Ct))
#basic graph
ggplot(RT.long, aes(x = NE, y = Mouse_ID)) +
  geom_point() +
  facet_wrap("Target")

# ggplot(POS, aes(x = HI, y = NE)) +
#   geom_point() +
#   geom_smooth() +
#   facet_wrap("Target", scales = "free_y") +
#   theme(axis.text=element_text(size=12, face = "bold"), 
#         axis.title=element_text(size=14,face="bold"),
#         strip.text.x = element_text(size = 14, face = "bold"),
#         legend.text=element_text(size=12, face = "bold"),
#         legend.title = element_text(size = 12, face = "bold"))+
#   ggtitle("Eimeria positive samples gene expression")

HZ17genotype <- "https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/data/Field_data/HZ10_HZ17_Genotypes_47-510.csv"
HZ17genotype <- read.csv(text = getURL(HZ17genotype))

HZ18genotype <- "https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/data/Field_data/HZ18_Genotypes.csv"
HZ18genotype <- read.csv(text = getURL(HZ18genotype))

HZgenotype <- rbind.fill(HZ17genotype, HZ18genotype)
# subest by HI
HImus <- select(HZgenotype, HI, Mouse_ID, Year)
HImus <- distinct(HImus)
#load in dissections
# HZ18dissection <- "https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/data/Field_data/HZ18_Dissections.csv"
# HZ18dissection <- read.csv(text = getURL(HZ18dissection))
# 
# HZ16dissection <- "https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/data/Field_data/HZ16_Dissection_1-211.csv"
# HZ16dissection <- read.csv(text = getURL(HZ16dissection))
# 
# HZ17dissection <- "https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/data/Field_data/HZ17_Dissections_212-236.csv"
# HZ17dissection <- read.csv(text = getURL(HZ17dissection))
# 
# HZ17dissection2 <- "https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/data/Field_data/HZ17_Dissections_237-523.csv"
# HZ17dissection2 <- read.csv(text = getURL(HZ17dissection2))
# 
# HZdissections <- rbind.fill(HZ18dissection, HZ16dissection)
# HZdissections <- rbind.fill(HZdissections, HZ17dissection)
# HZdissections <- rbind.fill(HZdissections, HZ17dissection2)
# 
# #subset by columns relevant for mapping and gene expression 
# diss <- select(HZdissections, Mouse_ID, Latitude, Longitude, Year, Sex, Status, Body_weight, Spleen, ASP, SYP, HET, MART, CP, HD, HM, MM, TM)
# diss <- unique(diss)
# # merge HImus and diss (MAKES MESS WITH YEARS)
# HImus <- merge(HImus, diss, by = "Mouse_ID")

# merge HImus and RT (rename from HZ18 because it is NOT!°!!!!!!!!!!!!!)
HZ16_17 <- merge(RT.long, HImus)
HZ16_17$Year <- as.factor(HZ16_17$Year)
HZ18G <- select(HZ18genotype, Mouse_ID, Year, HI)
HZ18 <- "https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/data/Gene_expression/HZ18_RT-qPCR_RTlong.csv"
HZ18 <- read.csv(text = getURL(HZ18))
HZ18$inf <- NULL
HZ18$deltaCtMmE_tissue <- NULL
HZ18 <- merge(HZ18G, HZ18)
HZ18$Target <- as.character(HZ18$Target)
HZ18 <- HZ18[!(HZ18$Target == "IL-6"),]
HZ18 <- HZ18[!(HZ18$Target == "GBP2"),]
HZ18$Target[HZ18$Target == "IL-12b"] <- "IL.12"
HZ18$Target <- as.factor(HZ18$Target)
HZ <- rbind(HZ18, HZ16_17)
HZ <- distinct(HZ,.keep_all = TRUE)

ggplot(HZ, aes(x = NE, y = HI)) +
  geom_point() +
  coord_flip() +
  facet_wrap("Target", scales = "free_y") +
  theme(axis.text=element_text(size=12, face = "bold"), 
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"))+
  ggtitle("Overall wild gene expression")

# ### Load in known MC positives and merge together
HZ16_17MC <- "https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/data/Eimeria_detection/HZ16-17_InfInt_MC_Lorenzo%26Mert.csv"
HZ16_17MC <- read.csv(text = getURL(HZ16_17MC), na.strings = c("", "NA"))
TruePositives1 <- subset(HZ16_17MC, Caecum == "pos")
#just to assign years, fix later by rewriting the table
years_add <- "https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/data/Eimeria_detection/FINALqpcrData_2016_2017_threshold3.75.csv"
years_add <- read.csv(text = getURL(years_add))
years_add <- select(years_add, year, Mouse_ID)
HZ16_17MC <- merge(HZ16_17MC, years_add)
colnames(HZ16_17MC)[7] <- "Year"
TruePositives1 <- select(HZ16_17MC, Mouse_ID, Year, Caecum)
# write.csv(TruePositives, file = "~/Mouse_Eimeria_Databasing/Mouse_Eimeria_Databasing/data/Eimeria_detection/MC_verified_positives.csv")
HZ18MC <- "https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/data/Eimeria_detection/Svenja/table_ct_and_more.csv"
HZ18MC <- read.csv(text = getURL(HZ18MC))
colnames(HZ18MC)[4] <- "Caecum"
HZ18MC$Year <- "2018"
HZ18MC$Caecum <- as.character(HZ18MC$Caecum)
HZ18MC$Caecum[HZ18MC$Caecum == "TRUE"] <- "pos"
HZ18MC$Caecum[HZ18MC$Caecum == "FALSE"] <- "neg"
TruePositives2 <- subset(HZ18MC, Caecum == "pos")
TruePositives1 <- select(TruePositives1, Mouse_ID, Caecum, Year)
colnames(TruePositives2)[1] <- "Mouse_ID"
TruePositives2 <- select(TruePositives2, Mouse_ID, Caecum, Year)
TruePositives2 <- separate(TruePositives2, c("Mouse_ID"), into = c("Tissue", "AA", "Mouse_ID"))
TruePositives2$Mouse_ID <- sub("^", "AA_0", TruePositives2$Mouse_ID)
TruePositives2$Tissue <- NULL
TruePositives2$AA <- NULL
TruePositives <- rbind(TruePositives1, TruePositives2)
TruePositives <- subset(TruePositives, Caecum == "pos")

# load in all known negatives and merge together
TrueNegatives1 <- subset(HZ16_17MC, Caecum == "neg")
TrueNegatives2 <- subset(HZ18MC, Caecum == "neg")
TrueNegatives1 <- select(TrueNegatives1, Mouse_ID, Caecum, Year)
TrueNegatives2 <- separate(TrueNegatives2, c("Name"), into = c("Tissue", "AA", "Mouse_ID"))
TrueNegatives2$Mouse_ID <- sub("^", "AA_0", TrueNegatives2$Mouse_ID)
TrueNegatives2$Tissue <- NULL
TrueNegatives2$AA <- NULL
TrueNegatives2 <- select(TrueNegatives2, Mouse_ID, Caecum, Year)
TrueNegatives <- rbind(TrueNegatives1, TrueNegatives2)
Trues <- rbind(TrueNegatives, TruePositives)

HZ1 <- merge(HZ, Trues)

# HZ1 <- subset(HZ1, Caecum == "pos" & Caecum == "neg")

#HZ1 <- subset(HZ1, Caecum == "pos" & Caecum == "neg")

# HZ1$Caecum <- replace_na(HZ1$Caecum, "neg")
colnames(HZ1)[6] <- "MC"
# compare in one big DF and ggplot to see POS vs NEG
ggplot(data=subset(HZ1, !is.na(x = HZ1$Target)), aes(x = HI, y = NE, color = MC)) +
  geom_point() +
  geom_smooth() +
  facet_wrap("Target", scales = "free_y") +
  theme(axis.text=element_text(size=12, face = "bold"), 
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"))+
  ggtitle("Overall wild gene expression vs HI")

##################################################################################
# now load in intensity data and add it to HZ1
int1 <- "https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/data/Eimeria_detection/FINALqpcrData_2016_2017_threshold3.75.csv"
int1 <- read.csv(text = getURL(int1))

int1 <- select(int1, Mouse_ID, delta_ct_cewe_MminusE, year)

int1 <- dplyr::select(int1, Mouse_ID, delta_ct_cewe_MminusE, year)

colnames(int1)[2] <- "delta"
colnames(int1)[3] <- "Year"

int2 <- select(HZ18MC, Name, deltaCtMmE_tissue, Year)
colnames(int2)[2] <- "delta"
int2 <- separate(int2, c("Name"), into = c("Tissue", "AA", "Mouse_ID"))
int2$Mouse_ID <- sub("^", "AA_0", int2$Mouse_ID)
int2$Tissue <- NULL
int2$AA <- NULL


int <- rbind(int1, int2)

HZ1 <- merge(int, HZ1)

ggplot(data=subset(HZ1, !is.na(x = HZ1$Target) & !is.na(x = HZ1$MC)), aes(x = delta, y = NE, color = MC)) +
  geom_point() +
  geom_smooth() +
  facet_wrap("Target", scales = "free_y") +
  theme(axis.text=element_text(size=12, face = "bold"), 
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"))+
  ggtitle("Overall wild gene expression vs delta")
# write out the RT-qPCRs that need to be done to fill in missing positives

missing <- HZ1[is.na(HZ1$NE),]
write.csv(missing, "~/Mouse_Eimeria_Databasing/data/Gene_expression/MC_identified_extra_samples_to_process.csv")

# pick out high delta, MC negative samples
High_delta_negs <- subset(HZ1, MC == "neg" & delta > -5)
write.csv(High_delta_negs, "~/Mouse_Eimeria_Databasing/data/Eimeria_detection/HZ16-18_high_delta_negatives.csv")

## Add NA group for negative MCs as Emanuel suggested (converts numbers to chr so X axis looks wild)
HZ1$Eimstatus <- ifelse(HZ1$MC == "neg", "NA", HZ1$delta)
HZ1$Eimstatus <- as.numeric(HZ1$Eimstatus)


ggplot(HZ1, aes(x = Eimstatus, y = NE, color = MC)) +
  geom_point() +
  geom_smooth() +
  geom_miss_point() +
  facet_wrap("Target", scales = "free") +
  theme(axis.text=element_text(size=12, face = "bold"), 
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"))+
  ggtitle("Overall wild gene expression vs delta")

############################## Add oocyst data
oocysts <- "https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/data/Eimeria_detection/Eimeria_oocysts_2015%262017_Lorenzo.csv"
oocysts <- read.csv(text = getURL(oocysts))
oocysts <- select(oocysts, Mouse_ID, mean_neubauer, OPG)
HZ16and17 <- merge(oocysts, HZ16and17, by = "Mouse_ID", all.y = TRUE)
# check for oocyst positive and qPCR negative flukes (keeps reducing the overall number... maybe stick to oocyst and HZ merge)
DoublePositive <- dplyr::filter(HZ16and17, qPCRstatus == "positive", OPG > 0)
DoublePositive$positive <- "double"
OocystPositive <- dplyr::filter(HZ16and17, qPCRstatus == "negative", OPG > 0)
OocystPositive$positive <- "oocyst"
qPCRPositive <- dplyr::filter(HZ16and17, qPCRstatus == "positive", OPG == 0)
qPCRPositive$positive <- "qPCR"
PositiveInvestigate <- rbind(DoublePositive, OocystPositive, qPCRPositive)
PositiveInvestigate <- distinct(PositiveInvestigate)
HZ16and17 <- merge(HZ16and17, PositiveInvestigate, by = "Mouse_ID")

################### graph like Svenja's results #######################
HZ2 <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/data/Gene_expression/HZ18_RT-qPCR_RTlong.csv"))
HZ1 <- select(HZ1, Mouse_ID, Target, delta, NE, HI, MC)
HZ2 <- select(HZ2, Mouse_ID, Target, deltaCtMmE_tissue, inf, NE, HI)
#### name columns to match
colnames(HZ2)[3] <- "delta"
colnames(HZ2)[4] <- "MC"
HZ1$MC <- as.character(HZ1$MC)
HZ1$MC[HZ1$MC == "pos"] <- "TRUE"
HZ1$MC[HZ1$MC == "neg"] <- "FALSE"
HZ2$Target <- as.character(HZ2$Target)
HZ2$Target[HZ2$Target == "IL-12b"] <- "IL.12"
HZgraph <- rbind(HZ1, HZ2)


ggplot(HZgraph, 
       aes(x = MC, y = NE, color = MC)) +
  geom_jitter() +
  geom_boxplot() +
  facet_wrap("Target", scales = "free") +
  theme(axis.text=element_text(size=12, face = "bold"), 
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"))+
  ggtitle("HZ16-18_gene_downregulation")

ggplot(HZgraph, 
       aes(x = MC, y = NE, color = MC)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter() +
  stat_compare_means(aes(label = ..p.signif..)) +
  facet_wrap("Target", scales = "free") +
  labs(y="deltaCT = Target - HKG", x = "infected", colour = "infected") +
  theme(axis.text=element_text(size=12, face = "bold"),
        title = element_text(size = 16, face = "bold"),
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"))+
  ggtitle("Gene expression in wild samples")

write.csv(HZgraph, "C:/Users/Luke Bednar/Mouse_Eimeria_Databasing/data/Gene_expression/HZ16-18_gene_expression.csv")
