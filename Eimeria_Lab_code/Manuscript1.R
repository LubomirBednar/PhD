# Manuscript 1:

# load libraries
library(ggplot2)
library(dplyr)


##### create full WDS 
# E6 
E6 <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/E5_062018_Eim_complete.csv")
E6 <- select(E6, "EH_ID", "labels", "experiment", "dpi", "OPG", "weightloss", "primary_infection", "Strain",
             "Batch")
E6$OPG <- as.numeric(as.character(E6$OPG))


E7 <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experimental_design/E7_112018_Eim_DESIGN.csv")
E10 <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experimental_design/E10_112020_Eim_DESIGN.csv")
E11 <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experimental_design/E11_012021_Eim_DESIGN.csv")








P3 <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experimental_design/P3_112019_Eim_design.csv")
P4 <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experimental_design/P4_082020_Eim_DESIGN.csv")
# check colnames and select to combine
names(E6)[names(E6) == "Batch"] <- "batch"
E6 <- select(E6, "EH_ID", "mouse_strain", "Sex", "primary_infection", "experiment", "batch")



##### join all existing records

##### join all existing oocysts

##### join all existing FACS

##### join all existing ELISAs

##### join all existing genes































# should contain P3 and P4
CLS <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/CLS_complete.csv")
CLS$X <- NULL
names(CLS)[names(CLS) == "Strain"] <- "mouse_strain"
names(CLS)[names(CLS) == "labels"] <- "label"
names(CLS)[names(CLS) == "primary"] <- "primary_infection"
names(CLS)[names(CLS) == "challenge"] <- "challenge_infection"
names(CLS)[names(CLS) == "Position"] <- "FACS_tissue"
# select for constant column names
CLS <- select(CLS, "label", "EH_ID", "dpi", "relative_weight", "batch", "primary_infection", 
              "challenge_infection", "infection_history","OPG","experiment", "mouse_strain",
              "EXP_type", "Eim_MC", "delta", "FACS_tissue", "CD4", "Treg", "Div_Treg", "Treg17",
              "Th1", "Div_Th1", "Th17", "Div_Th17", "CD8", "Act_CD8", "Div_Act_CD8", "IFNy_CD4",
              "IFNy_CD8", "IFNy_MES", "IFNy_CEWE")

##### WDS table  E6, E7, E10, E11
E7immuno <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/E7_112018_Eim_immuno.csv")
E7immuno$X <- NULL
E7immuno$batch <- "b"
E7immuno$experiment <- "E7"
names(E7immuno)[names(E7immuno) == "primary"] <- "primary_infection"
names(E7immuno)[names(E7immuno) == "challenge"] <- "challenge_infection"
names(E7immuno)[names(E7immuno) == "Position"] <- "FACS_tissue"
names(E7immuno)[names(E7immuno) == "weight_change"] <- "relative_weight"
names(E7immuno)[names(E7immuno) == "Strain"] <- "mouse_strain"



E7immuno<- select(E7immuno, "label", "EH_ID", "dpi", "relative_weight", "batch", "primary_infection", 
                  "challenge_infection", "infection_history","OPG","experiment", "mouse_strain",
                  "EXP_type", "Eim_MC", "delta", "FACS_tissue", "CD4", "Treg", "Div_Treg", "Treg17",
                  "Th1", "Div_Th1", "Th17", "Div_Th17", "CD8", "Act_CD8", "Div_Act_CD8", "IFNy_CD4",
                  "IFNy_CD8", "IFNy_CEWE")

E7para <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/E7_112018_Eim_COMPLETE.csv")
E7para$X <- NULL
E7para$batch <- "b"
E7para$experiment <- "E7"
E7para$EXP_type <- "WDS"

names(E7para)[names(E7para) == "labels"] <- "label"
names(E7para)[names(E7para) == "Strain"] <- "mouse_strain"
names(E7para)[names(E7para) == "weight_change"] <- "relative_weight"
names(E7para)[names(E7para) == "primary"] <- "primary_infection"
names(E7para)[names(E7para) == "challenge"] <- "challenge_infection"
names(E7para)[names(E7para) == "IL.12"] <- "IL12"

E7para <- select(E7para, "label", "EH_ID", "dpi", "relative_weight", "batch", "primary_infection", 
                  "challenge_infection", "infection_history","OPG","experiment", "mouse_strain",
                  "EXP_type", "Eim_MC", "delta", "IFNy_CEWE", "IL12", "CXCR3", "IRG6")




E10 <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/E10_112020_Eim_record.csv")
E10$X <- NULL
names(E10)[names(E10) == "labels"] <- "label"
E10$EXP_type <- "WDS"
E10<- select(E10, "label", "EH_ID", "dpi", "relative_weight", "batch", "primary_infection", 
                 "challenge_infection", "infection_history","experiment", "mouse_strain",
                 "EXP_type")

E11 <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/E11_012021_Eim_record.csv")
E11$X <- NULL
names(E11)[names(E11) == "labels"] <- "label"
E11$EXP_type <- "WDS"
E11<- select(E11, "label", "EH_ID", "dpi", "relative_weight", "batch", "primary_infection", 
             "challenge_infection", "infection_history","experiment", "mouse_strain",
             "EXP_type")


WDS <- E7immuno %>%
  full_join(E11) %>%
  full_join(E10) %>%
  full_join(E7para)

