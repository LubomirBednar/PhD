# create design table
# Load information table
library(RCurl)
library(httr)

# load in initial dataset from GitHub (must be raw.)
infoTable <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experimental_design/E11_012021_Eim_info.csv"))
# check last experiment and get highest EH_ID
lastEH_ID <- "last used EH_ID"
# divide dataset into groups as desired
## e.g.: 100 mice would divide into c(rep("isolate1", 25), rep("isolate2", 25), 
##                                  rep("isolate3", 25), rep("uninfected", 25))
# for challenge infection plans this has to be spread further
infection_history <- c(rep("UNI:UNI", 3), rep("E88:UNI", 3), rep("E64:UNI", 3), rep("UNI:E64", 3), rep("UNI:E88", 3), 
                       rep("E88:E88", 4), rep("E64:E64", 4), rep("E88:E64", 4), rep("E64:E88", 4))
# number of mice
Nmice = nrow(infoTable)
#Give EH_IDs
num = as.numeric(sub("LM", "", lastEH_ID))
num = num + (1:(Nmice))
EH_ID = paste0("LM", sprintf("%04d", num))
#Assign infection isolate
designTable <- data.frame(primary_infection = primary_infection,
                          challenge_infection = challenge_infection,
                          EH_ID= EH_ID)
# Spread names and infections randomly among mice (restrict by total amount of mice)
infoTable$EH_ID <- sample(EH_ID, size = 54)
infoTable$primary_infection <- sample(primary_infection, size = 54)
infoTable$challenge_infection <- sample(challenge_infection, size = 54)
# merge
finaldesignTable <- merge(infoTable, designTable, all.x = T)
# write out
write.csv(finaldesignTable,
          "location on local machine",
          row.names = F, quote = F )