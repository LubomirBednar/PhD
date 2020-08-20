library(Rmisc)
library(httr)
library(RCurl)
library(arsenal)

# read in record tables
P4a <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/P4a_082020_Eim_Record.csv"))
P4b <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/P4b_082020_Eim_Record.csv"))

#make random labels
random_labels <- do.call(paste0, replicate(3, sample(LETTERS, 480, TRUE), FALSE))
RLa <- random_labels[1:240]
RLb <- random_labels[241:480]

#assign to label columns
P4a$labels <- RLa
P4b$labels <- RLb

#check if labels are unique
comparedf(P4a, P4b)
summary(comparedf(P4a, P4b))

#write out
write.csv(P4a, "../Eimeria_Lab/data/Experiment_results/P4a_082020_Eim_Record.csv")
write.csv(P4b, "../Eimeria_Lab/data/Experiment_results/P4b_082020_Eim_Record.csv")
