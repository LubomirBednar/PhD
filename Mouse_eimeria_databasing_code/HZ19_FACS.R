library(httr)
library(RCurl)
library(Rmisc)
library(tidyverse)
library(readxl)
library(dplyr)
# guess this one will have to be hard coded (win)
FACSraw1 <- read_xlsx("/Users/Luke Bednar/Mouse_Eimeria_Databasing/data/Field_data/HZ19_FACS_raw.xlsx", sheet = 1)
FACSraw2 <- read_xlsx("/Users/Luke Bednar/Mouse_Eimeria_Databasing/data/Field_data/HZ19_FACS_raw.xlsx", sheet = 2)
FACSraw3 <- read_xlsx("/Users/Luke Bednar/Mouse_Eimeria_Databasing/data/Field_data/HZ19_FACS_raw.xlsx", sheet = 3)
FACSraw4 <- read_xlsx("/Users/Luke Bednar/Mouse_Eimeria_Databasing/data/Field_data/HZ19_FACS_raw.xlsx", sheet = 4)
# guess this one will have to be hard coded (linux)
FACSraw1 <- read_xlsx("./Mouse_Eimeria_Databasing/data/Field_data/HZ19_FACS_raw.xlsx", sheet = 1)
FACSraw2 <- read_xlsx("./Mouse_Eimeria_Databasing/data/Field_data/HZ19_FACS_raw.xlsx", sheet = 2)
FACSraw3 <- read_xlsx("./Mouse_Eimeria_Databasing/data/Field_data/HZ19_FACS_raw.xlsx", sheet = 3)
FACSraw4 <- read_xlsx("./Mouse_Eimeria_Databasing/data/Field_data/HZ19_FACS_raw.xlsx", sheet = 4)




# extract sample names and position 
FACSraw1$Mouse_ID <-gsub("\\d+: (mLN|spleen)_(\\d{3})_\\d{3}.fcs", "AA_0\\2", FACSraw1$sample)
FACSraw1$Position <- gsub("\\d+: (mLN|spleen)_(\\d{3})_\\d{3}.fcs", "\\1", FACSraw1$sample)
FACSraw2$Mouse_ID <-gsub("\\d+: (mLN|spleen)_(\\d{3})_\\d{3}.fcs", "AA_0\\2", FACSraw2$sample)
FACSraw2$Position <- gsub("\\d+: (mLN|spleen)_(\\d{3})_\\d{3}.fcs", "\\1", FACSraw2$sample)
FACSraw3$Mouse_ID <-gsub("\\d+: (mLN|spleen)_(\\d{3})_\\d{3}.fcs", "AA_0\\2", FACSraw3$sample)
FACSraw3$Position <- gsub("\\d+: (mLN|spleen)_(\\d{3})_\\d{3}.fcs", "\\1", FACSraw3$sample)

FACSraw3.1 <- subset(FACSraw3[1:45,])
FACSraw3.2 <- subset(FACSraw3[46:90,])
FACSraw3.2$Mouse_ID <-paste(gsub("(mLN|spleen)_(\\d{3})_\\d{3}.fcs", "AA_0\\2", FACSraw3.2$sample))
FACSraw3.2$Position <-paste(gsub("(mLN|spleen)_(\\d{3})_\\d{3}.fcs", "\\1", FACSraw3.2$sample))

FACSraw4$Mouse_ID <-gsub("(mLN|spleen)_(\\d{3})_\\d{3}.fcs", "AA_0\\2", FACSraw4$Sample)
FACSraw4$Position <- gsub("(mLN|spleen)_(\\d{3})_\\d{3}.fcs", "\\1", FACSraw4$Sample)
# remove that strage Sample column
FACSraw1 <- FACSraw1[,-c(1)]
FACSraw2 <- FACSraw2[,-c(1)]
FACSraw3.1 <- FACSraw3.1[,-c(1)]
FACSraw3.2 <- FACSraw3.2[,-c(1)]
FACSraw4 <- FACSraw4[,-c(1)]
# combine into one and remove wrong sample (Hongwei said) makes some NAs
FACS <- full_join(FACSraw1, FACSraw2)
FACS <- full_join(FACS, FACSraw3.1)
FACS <- full_join(FACS, FACSraw3.2)
FACS <- full_join(FACS, FACSraw4)

#####################################################################################################################################

#create R and .csv friendly column names
colnames(FACS)[1]<- "CD4"
colnames(FACS)[2]<- "Treg"
colnames(FACS)[3]<- "Div_Treg"
colnames(FACS)[4]<- "Treg17"
#remove this until known
FACS$`FSC-A, SSC-A subset/single/live/CD4+/Foxp3-,Freq. of Parent` <- NULL
colnames(FACS)[5]<- "Th1"
colnames(FACS)[6]<- "Div_Th1"
colnames(FACS)[7]<- "Th17"
colnames(FACS)[8]<- "Div_Th17"
colnames(FACS)[9]<- "CD8"
colnames(FACS)[10]<- "Act_CD8"
colnames(FACS)[11]<- "Div_Act_CD8"
colnames(FACS)[12]<- "IFNy_CD4"
colnames(FACS)[13]<- "IL17A_CD4"
colnames(FACS)[14]<- "IFNy_CD8"

# remove "%" signs and convert to numeric
FACS[] <- lapply(FACS, gsub, pattern='%', replacement='')
FACS[, 1:14] <- sapply(FACS[, 1:14], as.numeric)
FACS$Position <- as.factor(FACS$Position)
# write out
write.csv(FACS, "/Users/Luke Bednar/Mouse_Eimeria_Databasing/data/HZ19_MES_FACS.csv")

###########################################################################
FACS <- read.csv("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Field/master/data_input/HZ19_immuno.csv")
FACSmln <- subset(FACS, FACS$Position == "mLN")

FACSmln <- select(FACSmln, Mouse_ID, CD4, Treg, Div_Treg, Treg17, Th1, Div_Th1, Th17, Div_Th17, CD8, Act_CD8, Div_Act_CD8, 
                  IFNy_CD4, IFNy_CD8)

FACSmln <- FACSmln %>% remove_rownames %>% column_to_rownames(var="Mouse_ID")

pca <- prcomp(na.omit(FACSmln), center = TRUE, scale = TRUE)
pca.var <- pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100, 1)
barplot(pca.var.per, main = "Scree Plot", xlab = "Principal component",
        ylab = "Percent Variation")

pca.data <- data.frame(Sample=rownames(pca$x),
                       X = pca$x[,1],
                       Y = pca$x[,2])
names(pca.data)[names(pca.data) == "Sample"] <- "Mouse_ID"
pca.data <- remove_rownames(pca.data)
FACSmln <- select(FACS, Mouse_ID, EXP_type, Position, MC.Eimeria, delta, IFNy, CD4, Treg, Div_Treg, Treg17, Th1, 
                  Div_Th1, Th17, Div_Th17, CD8, Act_CD8, Div_Act_CD8, 
                  IFNy_CD4, IFNy_CD8#, CASP1, CXCL9, CXCR3, IDO1,  IFNG,  IL10, IL12A, IL13, IL1RN, IRGM1, MPO, MUC2, 
                  #MUC5AC, MYD88, NCR1, PRF1,  RETNLB, SOCS1, TICAM1, TNF, IL6, IL17A)
)
FACSmln <- subset(FACSmln, FACSmln$Position == "mLN")
FACSmln <- distinct(FACSmln)

pca.data <- merge(pca.data, FACSmln, by = "Mouse_ID")
loading_scores <- pca$rotation[,1]
scores <- abs(loading_scores)
score_ranked <- sort(scores, decreasing = T)
top10 <- names(score_ranked[1:13])
pca$rotation[top10,1]
pca.var.per.df <- cbind(pca.var.per, top10)
pca.var.per.df <- data.frame(pca.var.per.df)
pca.var.per.df$pca.var.per <- as.numeric(as.character(pca.var.per.df$pca.var.per))
pca.var.per.df$top10 <- factor(pca.var.per.df$top10,levels = top10)

ggplot(pca.var.per.df, aes(top10, pca.var.per, fill = top10)) +
  geom_col() +
  theme(strip.text.x = element_text(size = 14, face = "bold"),
        strip.text.y = element_text(size = 14, face = "bold")) +
  scale_x_discrete(name ="Cell populations") +
  scale_y_discrete(name ="% of variation explained") +
  labs(fill = "") +
  
  geom_text(aes(label = pca.var.per)) +
  theme_bw()

# remove AA0791 and 0699 because it's very strange (probably bad measurement)
pca.data <- pca.data[pca.data$Mouse_ID != "AA_0791", ]
pca.data <- pca.data[pca.data$Mouse_ID != "AA_0699", ]

ggplot(subset(pca.data, !is.na(pca.data$MC.Eimeria)), aes(x = X, y = Y, label = Mouse_ID, color = MC.Eimeria)) +
  geom_point(size = 3)+
  xlab(paste("PC1 - ", pca.var.per[1], "%", sep = "")) +
  ylab(paste("PC2 - ", pca.var.per[2], "%", sep = "")) +
  theme_bw() +
  theme(axis.text=element_text(size=12, face = "bold"),
        title = element_text(size = 16, face = "bold"),
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        strip.text.y = element_text(size = 14, face = "bold"),
        legend.text = element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"))+
  ggtitle("T-cell populations PCA by infection status") +
  labs(color='infected') 
# looks good, test for batch independence, add experiment column
ggplot(data = pca.data, aes(x = X, y = Y, label = EH_ID, color = experiment)) +
  geom_point(size = 3)+
  xlab(paste("PC1 - ", pca.var.per[1], "%", sep = "")) +
  ylab(paste("PC2 - ", pca.var.per[2], "%", sep = "")) +
  theme_bw() +
  theme(axis.text=element_text(size=12, face = "bold"),
        title = element_text(size = 16, face = "bold"),
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        strip.text.y = element_text(size = 14, face = "bold"),
        legend.text = element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"))+
  ggtitle("Checking for batch effect") +
  labs(color='Experiment ID') 

# # write out long for HZ19 immuno
# SpleenDF <- filter(FACS, Position == "spleen")
# mLNDF <- filter(FACS, Position == "mLN")
# FACS.long.mln <- melt(mLNDF,
#              direction = "long",
#              varying = list(names(mLNDF)[1:14]),
#              v.names = "cell.pop",
#              na.rm = T, value.name = "counts", 
#              id.vars = c("Mouse_ID", "Position"))
# FACS.long.mln <- na.omit(FACS.long.mln)
# names(FACS.long.mln)[names(FACS.long.mln) == "variable"] <- "pop"
# write.csv(FACS.long.mln, "/Users/Luke Bednar/Mouse_Eimeria_Databasing/data/Field_data/HZ19_FACS_long_mln.csv")
# #################### process like E7 
# FACS <- dplyr::select(HZ19, Mouse_ID, CD4p, CD8p, Th1IFNgp_in_CD4p, Th17IL17Ap_in_CD4p, Tc1IFNgp_in_CD8p, Treg_Foxp3_in_CD4p,
#                       Dividing_Ki67p_in_Foxp3p, RORgtp_in_Foxp3p, Th1Tbetp_in_CD4pFoxp3n, Dividing_Ki67p_in_Tbetp,
#                       Th17RORgp_in_CD4pFoxp3n, Dividing_Ki67p_in_RORgtp, Position, infHistory)
# 
# FACS <- dplyr::distinct(FACS)
# 
# 
# 
# 
# # tranform into long
# 
# FACS <- melt(FACS,
#              direction = "long",
#              varying = list(names(FACS)[2:13]),
#              v.names = "cell.pop",
#              na.rm = T, value.name = "counts", 
#              id.vars = c("EH_ID", "Position", "infHistory"))
# FACS <- na.omit(FACS)
# names(FACS)[names(FACS) == "variable"] <- "pop"
# 
# #################################
# ggplot(HZ19, aes(y =  , x = , color = )) +
#   geom_point() +
#   # ylim(2, -17) +
#   facet_grid(Target~infHistory, scales = "free") +
#   theme(axis.text=element_text(size=12, face = "bold"), 
#         axis.title=element_text(size=14,face="bold"),
#         strip.text.x = element_text(size = 14, face = "bold"),
#         legend.text=element_text(size=12, face = "bold"),
#         legend.title = element_text(size = 12, face = "bold"))+
#   ggtitle("")
# 
# 
# #test for normality
# 
# ## tabulate  medians for 
# # cell.medians <- lapply(facs.measure.cols, function (x){
# #   tapply(HZ19[, x], list(HZ19$Sex, as.factor(HZ19$tissue)), median)
# # })
# # names(cell.medians) <- facs.measure.cols
# # cell.medians
# 
# # cell.means <- lapply(facs.measure.cols, function (x){
# #   tapply(HZ19[, x], list(HZ19$infHistory, as.factor(HZ19$Position)), mean)
# # })
# # names(cell.means) <- facs.measure.cols
# # cell.means
# 
# #check distribution with histogram
# # histogram(~Sex | facs.measure.cols, data = HZ19)
# # histogram(~tissue | facs.measure.cols, data = HZ19)
# 
# ## #check overall populations between tissues
# # plotCells.tissue <- function (col){
# #   ggplot(HZ19, aes(tissue, get(col))) +
# #     geom_boxplot() +
# #     geom_jitter(width=0.2) +
# #     ggtitle(col)
# # }
# # 
# # facs_boxplots.tissue <- lapply(facs.measure.cols, plotCells.tissue)
# # names(facs_boxplots.tissue) <-  facs.measure.cols
# # 
# # for(i in seq_along(facs_boxplots.tissue)){
# #   pdf(paste0(names(facs_boxplots.tissue)[[i]], ".tissue.pdf"))
# #   plot(facs_boxplots.tissue[[i]])
# #   dev.off()
# # }
# 
# ### raw counts are modeled either as poisson or negative binomial in
# ### either case one could use the overall count (cell_counts) as
# ### "offset" to specify the "duration of observation" (normally
# ### offsets are used as a ratio, counto over time). I tried that, but
# ### then figured out that I know too little about how to interprete
# ### counts... expecially because the overall cell numbers are varying
# ### SO MUCH that this changes the results completely!!!
# 
# # distribution testing before modeling
# hist(HZ19$CD4p)
# descdist(HZ19$CD4p)
# descdist(HZ19$`Foxp3p_in_CD4p(Treg)`)
# 
# # subest by spleen
# SpleenDF <- filter(HZ19, tissue == "spleen")
# mLNDF <- filter(HZ19, tissue == "mLN")
# 
# # model interaction of cell populations with primary and secondary infection + constant position direction (PRIMARY : SECONDARY + POSITION)
# mods.l <- lapply(facs.measure.cols, function (x) {
#   lm(get(x) ~ (Body_weight * Spleen),
#      data=SpleenDF)
# })
# names(mods.l) <- facs.measure.cols
# lapply(mods.l, summary)
# 
# for(i in seq_along(facs.measure.cols)){
#   eff <- ggpredict(mods.l[[i]], terms=c("Body_weight", "Spleen"))
#   plot <-  plot(eff, rawdata=TRUE) +
#     scale_y_continuous(paste("percent", facs.measure.cols[[i]])) +
#     ggtitle(paste("predicted values of", facs.measure.cols[[i]]))
#   pdf(paste0(facs.measure.cols[[i]], ".spleen.pdf"))
#   print(plot)
#   dev.off()
# }
# 
# # model interaction of cell populations with primary, secondary infection and position (PRIMARY : SECONDARY : POSITION)
# mods.i <- lapply(facs.measure.cols, function (x) {
#   lm(get(x) ~ primary * challenge * Position,
#      data=HZ19)
# })
# names(mods.i) <- facs.measure.cols
# lapply(mods.i, summary)
# 
# for(i in seq_along(facs.measure.cols)){
#   eff <- ggpredict(mods.i[[i]], terms=c("primary", "challenge", "Position"))
#   plot <-  plot(eff, rawdata=TRUE) +
#     scale_y_continuous(paste("percent", facs.measure.cols[[i]])) +
#     ggtitle(paste("predicted values of", facs.measure.cols[[i]]))
#   pdf(paste0(facs.measure.cols[[i]], ".priXchaXpos.pdf"))
#   print(plot)
#   dev.off()
# }
# 
# # comparison of models
# lapply(seq_along(mods.i), function(i) anova(mods.i[[i]], mods.l[[i]]))
# 
# #check the model when using HI categories as well
# 
# modsHY.l <- lapply(facs.measure.cols, function (x) {
#   lm(get(x) ~ (primary * challenge) + Position + HybridStatus,
#      data=HZ19)
# })
# 
# names(modsHY.l) <- facs.measure.cols
# 
# lapply(modsHY.l, summary)
# 
# # comparison of all 3 models
# lapply(seq_along(mods.i), function(i) anova(mods.i[[i]], mods.l[[i]], modsHY.l[[i]]))
# ## And WOW (I reall wrote the above A PRIORY, otherwise... mayor
# ## fishing excursion ;-)...), but Tc1IFNgp_in_CD8p are lower in
# ## HYBRIDS look at THIS!!
# summary(modsHY.l[["Tc1IFNgp_in_CD8p"]])
# 
# WOW <- ggpredict(modsHY.l[["Tc1IFNgp_in_CD8p"]],
#                  terms=c("primary", "challenge", "HybridStatus"))
# 
# summary(modsHY.l[["Tc1IFNgp_in_CD8p"]])
# 
# pdf("WINNER_Tc1IFNgp_in_CD8p.effects.pdf")
# plot(WOW)
# dev.off()
# 
# 
# ## Now... I fooled myself a bit to that enthusiasm, as I expected
# ## "HybridStatusoutbred hybrids" to be ... well ... hybrids. Turns out
# ## this are the within subspecies outbreds. Let's do some PostHoc
# ## comparison. 
# 
# summary(glht(modsHY.l[["Tc1IFNgp_in_CD8p"]], mcp(HybridStatus="Tukey")))
# 
# ## nothing too shocking here, just that "outbred hybrids" have a trend
# ## towards lower cell proportions compared to "inter subsp. hybrids"
# 
# # ---------------------------------------------------------- Make connections between facets of models--------
# # transform data for graphing
# # HZ19.long <- reshape(data = HZ19, timevar = "infHistory", idvar = "EH_ID", direction = "long", varying = facs.measure.cols)
# 
# # HZ19.melt <- melt(setDT(HZ19), measure=patterns(facs.measure.cols), 
# #     value.name = facs.measure.cols, variable.name='EH_ID')
