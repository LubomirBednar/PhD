# Load information table
library(RCurl)
# load in initial dataset from GitHub
infoTable = read.csv(text = getURL("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experimental_design/E7_112018_Eim_DESIGN.csv"))
# check last experiment and get highest EH_ID
lastEH_ID <- "LM0226"
# divide dataset into groups as desired
## e.g.: 100 mice would divide into c(rep("isolate1", 25), rep("isolate2", 25), 
##                                  rep("isolate3", 25), rep("uninfected", 25))
infection_isolate <- c(rep("E64", 27), rep("E88", 27))
# number of mice
Nmice = nrow(infoTable)
#Give EH_IDs
num = as.numeric(sub("LM", "", lastEH_ID))
num = num + (1:(Nmice))
EH_ID = paste0("LM", sprintf("%04d", num))
#Assign infection isolate
designTable <- data.frame(infection_isolate = infection_isolate,
                           EH_ID= EH_ID)
#Spread names randomly among mice
infoTable$EH_ID <- sample(EH_ID)
#merge
finaldesignTable <- merge(infoTable, designTable)
# write out
write.csv(finaldesignTable,
         "location on local machine",
         row.names = F, quote = F )
# 
# library(experiment)
# 
# expe <- randomize(data = finaldesignTable, group = c("E64", "E88"),
#                    indx = finaldesignTable$EH_ID, block = finaldesignTable$Strain)
# 
#  trt <- data.frame(infection_isolate = expe$treatment)
#  trt$EH_ID <- rownames(trt)
#  rownames(trt) <- NULL





