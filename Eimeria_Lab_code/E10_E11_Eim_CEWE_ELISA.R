library(dplyr)
library(tidyr)
library(reshape2)
library(drc)
library(data.table)

E_std <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/E10_E11_Eim_CEWE_ELISA_std.csv")

E_samples <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/E10_E11_Eim_CEWE_ELISA_samples.csv")

model1<-drm(OD~Conc,
            fct=LL.4(names=c("Slope", "Lower", "Upper", "ED50")),
            data=E_std)
plot(model1)

E<-ED(model1, E_samples$OD, type="absolute", display=F)
row.names(E) <- E_samples$label

points(y=E_samples$OD,x=E[,1],col="lightblue",pch=19,cex=2)
text(y =E_samples$OD, x = E[,1], labels=E_samples$EH_ID, data=E, cex=0.9, font=2)

E <- data.frame(E)
colnames(E)[1] <- "IFNy"
E <- dplyr::select(E, IFNy)
setDT(E, keep.rownames = TRUE)[]
colnames(E)[1] <- "label"
E <- cbind(E_samples, E)
E$OD <- NULL
E$label <- NULL
E$EH_ID <- paste0("LM0", E$EH_ID)

write.csv(E, "./GitHub/Eimeria_Lab/data/Experiment_results/E10_E11_Eim_CEWE_ELISA.csv")
