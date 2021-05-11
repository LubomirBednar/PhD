library(Rmisc)
library(httr)
library(RCurl)
library(arsenal)

# read in record tables
E10a <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experimental_design/E10a_112020_Eim_design.csv"))

#make random labels use make.unique
random_labels <- do.call(paste0, replicate(3, sample(LETTERS, 384, TRUE), FALSE))
RLa <- random_labels[1:384]

#assign to label columns
E10a$labels <- RLa

#check if labels are unique
E10a$labels <- make.unique(E10a$labels, sep = "1")

#write out
write.csv(E10a, "../luke/Eimeria_Lab/data/Experimental_design/E10a_112020_Eim_design.csv")
