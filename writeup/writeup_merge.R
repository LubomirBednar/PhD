# processing all available data for PhD thesis

# import libraries
library(tidyr) # data warngling
library(dplyr) # data warngling
library(ggplot2) # graphs
library(RColorBrewer) # custom colours
library(AICcmodavg) # Akaines information criterion

# Start by getting all tables into the environment.
# LAB: design tables, weight loss, oocysts, ELISA, intensity, gene expression, FACS
# WILD: Hybrid Index (HI), oocysts, ELISA, intensity, gene expression, FACS


# 
# LAB DESIGN
design_E5 <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experimental_design/E57_xxxxx_Eim_DESIGN.csv")
design_E10 <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experimental_design/E10_112020_Eim_DESIGN.csv")
design_E11 <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experimental_design/E11_012021_Eim_DESIGN.csv")
design_P3 <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experimental_design/P3_112019_Eim_design.csv")
design_P4 <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experimental_design/P4_082020_Eim_DESIGN.csv")

# #standardise to in common columns
# design_E5 <- select(design_E5, 
#                     EH_ID, unique mouse identifier
#                     experiment, identifier of experiment
#                     mouse_strain, genetic background, equivalent to HI in wild
#                     primary_infection, 
#                     challenge_infection,
#                     infection_history combination of primary and challenge infection
#                     )

design_E10 <- select(design_E10, EH_ID, experiment, mouse_strain, primary_infection, challenge_infection, infection_history)
design_E11 <- select(design_E11, EH_ID, experiment, mouse_strain, primary_infection, challenge_infection, infection_history)

##add infection_history to P3 and P4
design_P3$infection_history <- paste(design_P3$primary_infection, design_P3$challenge_infection, sep = ":")
design_P4$infection_history <- paste(design_P4$primary_infection, design_P4$challenge_infection, sep = ":")

#bind designs into 1 table and remove clutter (rm function here, caution)
DESIGN <- bind_rows(design_P4, design_P3, design_E5, design_E11, design_E10) #32+32+117+16+20=217
DESIGN$EH_ID <- gsub(" ", "", DESIGN$EH_ID)
DESIGN$EH_ID <- gsub("_", "", DESIGN$EH_ID)
rm(list = ls("pattern" = "design_"))
write.csv(DESIGN, "../GitHub/PhD/writeup/DESIGN.csv")
write.csv(DESIGN, "~/Documents/GitHub/PhD/writeup/DESIGN.csv")



# LAB WEIGHT
weight_E10 <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/E10_112020_Eim_record.csv")
weight_E11 <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/E11_012021_Eim_record.csv")
weight_E5 <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/E57_xxxxx_Eim_record.csv")
weight_P3 <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/P3_112019_Eim_record.csv")
weight_P4 <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/P4_082020_Eim_record.csv")

#standardise to in common columns
WEIGHT <- bind_rows(weight_P4, weight_P3, weight_E5, weight_E11, weight_E10)
WEIGHT$X <- NULL
WEIGHT$EH_ID <- gsub(" ", "", WEIGHT$EH_ID)
WEIGHT$EH_ID <- gsub("_", "", WEIGHT$EH_ID)
rm(list = ls("pattern" = "weight_"))
write.csv(WEIGHT, "../GitHub/PhD/writeup/WEIGHT.csv")
write.csv(WEIGHT, "~/Documents/GitHub/PhD/writeup/WEIGHT.csv")


# 
# LAB OOCYST
oocyst_E10 <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/E10_112020_Eim_oocyst.csv")
oocyst_E11 <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/E11_Oocyst_updateCSV.csv")
oocyst_E5 <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/E57_xxxxx_Eim_oocyst.csv")
oocyst_P3 <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/P3_112019_Eim_oocyst.csv")
oocyst_P4 <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/P4_082020_Eim_oocyst.csv")
## standardise to in common columns
oocyst_P3 <- select(oocyst_P3, labels, oocyst_sq1, oocyst_sq2, oocyst_sq3, oocyst_sq4, dilution)
oocyst_P4 <- select(oocyst_P4, labels, oocyst_sq1, oocyst_sq2, oocyst_sq3, oocyst_sq4, dilution)
OOCYST <- bind_rows(oocyst_E10, oocyst_E11, oocyst_E5, oocyst_P3, oocyst_P4)

# join OOCYST, WEIGHT and DEISGN together into PARA
PARA_1 <- full_join(OOCYST, WEIGHT)
PARA_2 <- full_join(WEIGHT, DESIGN)
PARA <- full_join(PARA_1, PARA_2)

# calculate oocyst AVG, total and OPG
PARA$total_oocysts <- ((PARA$oocyst_sq1 
                        + PARA$oocyst_sq2 
                        + PARA$oocyst_sq3 
                        + PARA$oocyst_sq4) / 4) * 
  10000 * # because volume chamber
  2

PARA <- subset(PARA, PARA$labels != "E57aQSW")

PARA$OPG <- PARA$total_oocysts/ PARA$feces_weight

rm(list = ls("pattern" = "oocyst_"))
rm(list = ls("pattern" = "PARA_"))
write.csv(OOCYST, "../GitHub/PhD/writeup/OOCYST.csv")
write.csv(OOCYST, "~/Documents/GitHub/PhD/writeup/OOCYST.csv")

write.csv(PARA, "../GitHub/PhD/writeup/PARA.csv")
write.csv(PARA, "~/Documents/GitHub/PhD/writeup/PARA.csv")

# LAB ELISA
elisa_E10_E11 <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/E10_E11_Eim_CEWE_ELISA.csv")
elisa_E10_E11$X <- NULL
names(elisa_E10_E11)[names(elisa_E10_E11) == "IFNy"] <- "IFNy_CEWE"

elisa_P4 <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/P4_082020_Eim_CEWE_ELISA.csv")
elisa_P4$X <- NULL
names(elisa_P4)[names(elisa_P4) == "IFNy"] <- "IFNy_CEWE"
elisa_P3 <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/P3_112019_Eim_CEWE_ELISA.csv")
elisa_P3$X <- NULL
elisa_E7 <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/E7_112018_Eim_CEWE_ELISA.csv")
elisa_E7$X <- NULL
names(elisa_E7)[names(elisa_E7) == "IFNy"] <- "IFNy_CEWE"

elisa_P3 <- select(elisa_P3, EH_ID, IFNy_CEWE)
elisa_E7 <- select(elisa_E7, EH_ID, IFNy_CEWE)

ELISA <- bind_rows(elisa_E10_E11, elisa_E7, elisa_P3, elisa_P4)
ELISA$EH_ID <- gsub(" ", "", ELISA$EH_ID)
ELISA$EH_ID <- gsub("_", "", ELISA$EH_ID)
ELISA$dpi <- 8
rm(list = ls("pattern" = "elisa_"))
write.csv(ELISA, "../GitHub/PhD/writeup/ELISA.csv")
write.csv(ELISA, "~/Documents/GitHub/PhD/writeup/ELISA.csv")


# LAB INTENSITY
int_E10 <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/E10_112020_Eim_CEWE_qPCR.csv")
int_E10 <- select(int_E10, EH_ID, Eim_MC, delta)
int_E11 <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/E11_012021_Eim_CEWE_qPCR.csv")
int_E11 <- select(int_E11, EH_ID, Eim_MC, delta)
int_E7 <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/E7_112018_Eim_CEWE_qPCR.csv")
int_E7 <- select(int_E7, EH_ID, Eim_MC, delta) 
int_P3 <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/P3_112019_Eim_CEWE_qPCR.csv")
int_P3a <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/E7%26P3_Eim_MCs.csv")
int_P3a <- select(int_P3a, EH_ID, Eim_MC)
int_P3 <- left_join(int_P3, int_P3a, by.x = T)
int_P3$X <- NULL
int_P4 <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/P4_082020_Eim_CEWE_qPCR.csv")
int_P4$Eim_MC[int_P4$Eim_MC == "TRUE"] <- "pos"
int_P4$Eim_MC[int_P4$Eim_MC == "FALSE"] <- "neg"

INTENSITY <- bind_rows(int_E10, int_E11, int_E7, int_P3, int_P4)
INTENSITY$EH_ID <- gsub(" ", "", INTENSITY$EH_ID)
INTENSITY$EH_ID <- gsub("_", "", INTENSITY$EH_ID)
INTENSITY$dpi <- 8
rm(list = ls("pattern" = "int_"))
write.csv(INTENSITY, "../GitHub/PhD/writeup/INTENSITY.csv")
write.csv(INTENSITY, "~/Documents/GitHub/PhD/writeup/INTENSITY.csv")

# LAB GENE EXPRESSION
IFC1 <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/IFC1.csv")
IFC2 <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/IFC2.csv")
IFC3 <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/IFC3.csv")
IFC4 <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/IFC4.csv")
IFC5 <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/IFC5.csv")

IFC <- bind_rows(IFC1, IFC2, IFC3, IFC4, IFC5)
IFC$EH_ID <- as.character(sub("A0729", "AA0729", IFC$EH_ID))
# remove unsuccessful amplifications
IFC <- subset(IFC, IFC$Value != 999)

IFC <- IFC %>% dplyr::group_by(EH_ID, Target) %>% dplyr::summarise(Ct = mean(Value, na.rm = T))

IFC <- distinct(IFC)

IFC.wide <- pivot_wider(IFC, names_from = Target, values_from = Ct)
# IFC.wide <- data.frame(IFC.wide)
# write.csv(IFC.wide, "~/GitHub/PhD/IFC.wide.csv")
write.csv(IFC.wide, "./Documents/GitHub/PhD/writeup/IFC_wide.csv")
# IFC.wide <- read.csv("C:/Users/exemp/OneDrive/Documents/GitHub/PhD/IFC.wide.csv")
IFC.wide <- read.csv("./Documents/GitHub/PhD/writeup/IFC_wide.csv")
IFC.wide$X <- NULL
# IFC.wide[2:ncol(IFC.wide)] <- IFC.wide[2:ncol(IFC.wide)]-IFC.wide[,18]
GAPDH_GM <- na.omit(IFC.wide$GAPDH)
GAPDH_GM <- exp(mean(log(GAPDH_GM))) # geometric mean of GAPDH

PPIB_GM <- na.omit(IFC.wide$PPIB)
PPIB_GM <- exp(mean(log(PPIB_GM))) # geometric mean of PPIB
# > var(IFC.wide$PPIB, na.rm = T)
# [1] 15.62563
# > var(IFC.wide$GAPDH, na.rm = T)
# [1] 11.49809
IFC.wide$GAPDH <- 11.49809
IFC.wide$PPIB <- NULL
refGenes <- c("GAPDH")
targetGenes <- c("CASP1","CXCL9","CXCR3","IDO1", "IFNG", "IL10", "IL12A","IL13",
                  "IL1RN","IRGM1","MPO","MUC2","MUC5AC","MYD88","NCR1", "PRF1", "RETNLB", 
                  "SOCS1","TICAM1","TNF","IL6","IL17A")
 
IFC.wide <- IFC.wide %>% mutate(refMean = rowMeans(select(., refGenes)))
IFC.wide <- data.frame(IFC.wide)
refMean <- as.numeric(IFC.wide$refMean)

IFC.wide$PPIB <- NULL
IFC.wide$GAPDH <- NULL
IFC.wide[2:ncol(IFC.wide)] <- IFC.wide[2:ncol(IFC.wide)]-IFC.wide[,24]
IFC.wide$refMean <- NULL
IFC.wide$dpi <- 8
IFC_LAB <- filter(IFC.wide, !grepl("AA", EH_ID))
IFC_LAB <- filter(IFC_LAB, !grepl("NTC", EH_ID))
IFC_WILD <- filter(IFC.wide, !grepl("LM", EH_ID))

############################################################################################
IFC.long <- pivot_longer(IFC.wide, cols = c("CASP1","CXCL9","CXCR3","IDO1", "IFNG", "IL10", "IL12A","IL13",
                                            "IL1RN","IRGM1","MPO","MUC2","MUC5AC","MYD88","NCR1", "PRF1", "RETNLB", 
                                            "SOCS1","TICAM1","TNF","IL6","IL17A"))

names(IFC.long)[names(IFC.long) == "name"] <- "Target"
names(IFC.long)[names(IFC.long) == "value"] <- "Ct"
rm(list = "IFC1", "IFC2", "IFC3", "IFC4", "IFC5", "IFC", "PPIN_GM", "refGenes", "refMean")
# LAB FACS
E7_FACS <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/E7_112018_Eim_immuno.csv")
E7_FACS <- select(E7_FACS, EH_ID, CD4, Treg, Div_Treg, Treg17, Th1, Div_Th1, Th17, Div_Th17, CD8, Act_CD8, Div_Act_CD8, IFNy_CD4, IFNy_CD8)

CLS <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/CLS_complete.csv")
CLS$X <- NULL 
CLS$EH_ID <- as.character(sub("_", "", CLS$EH_ID))
names(CLS)[names(CLS) == "Strain"] <- "mouse_strain"
names(CLS)[names(CLS) == "primary"] <- "primary_infection"
names(CLS)[names(CLS) == "challenge"] <- "challenge_infection"
CLSimmuno <- select(CLS, EH_ID, experiment, dpi, labels, delta, IFNy_CEWE, Eim_MC, EXP_type, CD4, Treg, Div_Treg, Position,
                    Treg17, Th1, Div_Th1, Th17, Div_Th17, CD8, Act_CD8, Div_Act_CD8, IFNy_CD4, IFNy_CD8, batch)
CLSimmuno <- subset(CLSimmuno, CLSimmuno$dpi == 8)
CLSimmuno <- subset(CLSimmuno, CLSimmuno$batch == "b")
CLSimmuno <- subset(CLSimmuno, CLSimmuno$Position == "mLN")
CLS_FACS <- select(CLSimmuno, EH_ID, CD4, Treg, Div_Treg, Treg17, Th1, Div_Th1, Th17, Div_Th17, CD8, Act_CD8, Div_Act_CD8, IFNy_CD4, IFNy_CD8)
rm(list = "CLS", "CLSimmuno")
CLS_FACS <- distinct(CLS_FACS)
E11_FACS <- read.csv("./Documents/GitHub/Eimeria_Lab/data/Experiment_results/E11_012021_Eim_FACS.csv")
E11_FACS <- subset(E11_FACS, E11_FACS$Position == "mLN")
E11_FACS <- select(E11_FACS, EH_ID, CD4, Treg, Div_Treg, Treg17, Th1, Div_Th1, Th17, Div_Th17, CD8, Act_CD8, Div_Act_CD8, IFNy_CD4, IFNy_CD8)
FACS <- bind_rows(E7_FACS, E11_FACS, CLS_FACS)
rm(list = "CLS_FACS", "E11_FACS", "E7_FACS", "DESIGN", "OOCYST", "WEIGHT")
FACS$dpi <- 8
FACS$infection <- "challenge"
FACS$EH_ID <- gsub(" ", "", FACS$EH_ID)
FACS$EH_ID <- gsub("_", "", FACS$EH_ID)

# dpi 8 and infection challenge to all terminal measurements
IFC_LAB$infection <- "challenge"
INTENSITY$infection <- "challenge"
ELISA$infection <- "challenge"

IMMUNO <- full_join(ELISA, FACS)
IMMUNO <- full_join(IMMUNO, IFC_LAB)
IMMUNO <- full_join(IMMUNO, INTENSITY)

CHALLENGE <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data_products/Challenge_infections.csv")
IMCH <- full_join(CHALLENGE, IMMUNO)
IT <- full_join(PARA, IMMUNO)
write.csv(IMCH, "./Documents/GitHub/PhD/writeup/IMCH.csv")
write.csv(IT, "./Documents/GitHub/PhD/writeup/IT.csv")
