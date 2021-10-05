# processing all available data for PhD thesis

# import libraries
library(tidyverse) # data warngling
library(dplyr) # data warngling
library(ggplot2) # graphs
library(RColorBrewer) # custom colours
library(AICcmodavg) # Akaines information criterion

# Start by getting all tables into the environment.
# LAB: design tables, weight loss, oocysts, ELISA, intensity, gene expression, FACS
# WILD: Hybrid Index (HI), oocysts, ELISA, intensity, gene expression, FACS



# LAB DESIGN
design_E5 <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experimental_design/E57_xxxxx_Eim_DESIGN.csv")
design_E10 <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experimental_design/E10_112020_Eim_DESIGN.csv")
design_E11 <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experimental_design/E11_012021_Eim_DESIGN.csv")
design_P3 <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experimental_design/P3_112019_Eim_design.csv")
design_P4 <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experimental_design/P4_082020_Eim_DESIGN.csv")

## standardise to in common columns
design_E5 <- select(design_E5, 
                    EH_ID, # unique mouse identifier
                    experiment, # identifier of experiment
                    mouse_strain, # genetic background, equivalent to HI in wild
                    primary_infection, 
                    challenge_infection,
                    infection_history # combination of primary and challenge infection
                    )

design_E10 <- select(design_E10, EH_ID, experiment, mouse_strain, primary_infection, challenge_infection, infection_history)
design_E11 <- select(design_E11, EH_ID, experiment, mouse_strain, primary_infection, challenge_infection, infection_history)

### add infection_history to P3 and P4
design_P3$infection_history <- paste(design_P3$primary_infection, design_P3$challenge_infection, sep = ":")
design_P4$infection_history <- paste(design_P4$primary_infection, design_P4$challenge_infection, sep = ":")

## bind designs into 1 table and remove clutter (rm function here, caution)
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

## standardise to in common columns
WEIGHT <- bind_rows(weight_P4, weight_P3, weight_E5, weight_E11, weight_E10)
WEIGHT$X <- NULL
WEIGHT$EH_ID <- gsub(" ", "", WEIGHT$EH_ID)
WEIGHT$EH_ID <- gsub("_", "", WEIGHT$EH_ID)
rm(list = ls("pattern" = "weight_"))
write.csv(WEIGHT, "../GitHub/PhD/writeup/WEIGHT.csv")
write.csv(WEIGHT, "~/Documents/GitHub/PhD/writeup/WEIGHT.csv")



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

# LAB FACS


