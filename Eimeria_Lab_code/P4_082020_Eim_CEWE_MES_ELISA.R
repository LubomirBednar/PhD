library(httr)
library(RCurl)
library(Rmisc)
library(dplyr)
library(drc)
library(data.table)
# load in data
P4_std <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/P4_082020_Eim_CEWE_MES_ELISA_std.csv"))
P4_MES <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/P4_082020_Eim_MES_ELISA_samples.csv"))
P4_CEWE <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/P4_082020_Eim_CEWE_ELISA_samples.csv"))
# use drc to construct standard curve and pinpoint protein content
model1<-drm(OD~Conc,
            fct=LL.4(names=c("Slope", "Lower", "Upper", "ED50")),
            data=P4_std)
plot(model1)

P4C <- ED(model1, P4_CEWE$OD, type="absolute", display=F)
row.names(P4C) <- P4_CEWE$EH_ID

points(y=P4_CEWE$OD,x=P4C[,1],col="lightblue",pch=19,cex=2)
text(y =P4_CEWE$OD, x = P4C[,1], labels=P4_CEWE$EH_ID, data=P4C, cex=0.9, font=2)

P4C <- data.frame(P4C)
colnames(P4C)[1] <- "IFNy"
P4C <- dplyr::select(P4C, IFNy)
setDT(P4C, keep.rownames = TRUE)[]
colnames(P4C)[1] <- "EH_ID"
P4C <- merge(P4C, P4_CEWE)
P4C$OD <- NULL
# CEWE ELISA complete, write out
write.csv(P4C, "~/GitHub/Eimeria_Lab/data/Experiment_results/P4_082020_Eim_CEWE_ELISA.csv")

# moving on to MES, using same standard curve

P4_MES <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/P4_082020_Eim_MES_ELISA_samples.csv"))
P4M <- ED(model1, P4_MES$OD, type = "absolute", display = F)
row.names(P4M) <- P4_MES$EH_ID

points(y=P4_MES$OD,x=P4M[,1],col="lightblue",pch=19,cex=2)
text(y =P4_MES$OD, x = P4M[,1], labels=P4_MES$EH_ID, data=P4M, cex=0.9, font=2)

P4M <- data.frame(P4M)
colnames(P4M)[1] <- "IFNy"
P4M <- dplyr::select(P4M, IFNy)
setDT(P4M, keep.rownames = TRUE)[]
colnames(P4M)[1] <- "EH_ID"
P4M <- merge(P4M, P4_CEWE)
P4M$OD <- NULL

write.csv(P4M, "~/GitHub/Eimeria_Lab/data/Experiment_results/P4_082020_Eim_MES_ELISA.csv")
# merge and write out
colnames(P4C)[2] <- "IFNy_CEWE"
colnames(P4M)[2] <- "IFNy_MES"
P4E <- merge(P4M, P4C)

write.csv(P4E, "~/GitHub/Eimeria_Lab/data/Experiment_results/P4_082020_Eim_ELISA.csv")

