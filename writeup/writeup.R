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
rm(list = ls("pattern" = "design_"))
write.csv(DESIGN, "../GitHub/PhD/writeup/DESIGN.csv")



# LAB WEIGHT
weight_E10 <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/E10_112020_Eim_record.csv")
weight_E11 <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/E11_012021_Eim_record.csv")
weight_E5 <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/E57_xxxxx_Eim_record.csv")
weight_P3 <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/P3_112019_Eim_record.csv")
weight_P4 <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/P4_082020_Eim_record.csv")

## standardise to in common columns
WEIGHT <- bind_rows(weight_P4, weight_P3, weight_E5, weight_E11, weight_E10)
rm(list = ls("pattern" = "weight_"))
write.csv(DESIGN, "../GitHub/PhD/writeup/WEIGHT.csv")



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
rm(list = ls("pattern" = "oocyst_"))
write.csv(OOCYST, "../GitHub/PhD/writeup/OOCYST.csv")

# LAB ELISA


# LAB INTENSITY

# LAB GENE EXPRESSION

# LAB FACS


