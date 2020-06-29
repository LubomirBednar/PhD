# fecal ELISas HZ19 CEWE

library(httr)
library(RCurl)
library(Rmisc)
library(dplyr)
library(ggplot2)
library(tidyr)
library(reshape2)
library(drc)
library(data.table)

##### add clean tables ELISA 1
HZ_std1 <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/data/ELISAs/HZ19_CEWE_ELISA1_std.csv"))

HZ19_samples1 <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/data/ELISAs/HZ19_CEWE_ELISA1_samples.csv"))

###### use drc to construct standard curve and pinpointprotein content

model1<-drm(OD~Conc,
            fct=LL.4(names=c("Slope", "Lower", "Upper", "ED50")),
            data=HZ_std1)
plot(model1)

HZ1<-ED(model1, HZ19_samples1$OD, type="absolute", display=F)
row.names(HZ1) <- HZ19_samples1$Mouse_ID

points(y=HZ19_samples1$OD,x=HZ19_samples1[,1],col="lightblue",pch=19,cex=2)
text(y =HZ19_samples1$OD, x = HZ1[,1], labels=HZ19_samples1$Mouse_ID, data=HZ1, cex=0.9, font=2)

HZ1 <- data.frame(HZ1)
colnames(HZ1)[1] <- "IFNy"
HZ1 <- dplyr::select(HZ1, IFNy)
setDT(HZ1, keep.rownames = TRUE)[]
colnames(HZ1)[1] <- "Mouse_ID"
HZ1<- merge(HZ1, HZ19_samples1)
HZ1$OD <- NULL

# write.csv(HZ1, "./Documents/Mouse_Eimeria_Databasing/data/ELISAs/HZ19_CEWE_ELISA1_complete.csv")
write.csv(HZ1, "C:/Users/Luke Bednar/Mouse_Eimeria_Databasing/data/ELISAs/HZ19_CEWE_ELISA1_complete.csv")

############### load in ELISA2

E2_std <- "https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/data/ELISAs/HZ19_CEWE_ELISA2_std.csv"
E2_std <- read.csv(text = getURL(E2_std))

E2_samples <- "https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/data/ELISAs/HZ19_CEWE_ELISA2_samples.csv"
E2_samples <- read.csv(text = getURL(E2_samples))

###### use drc to construct standard curve and pinpointprotein content

model2<-drm(OD~Conc,
            fct=LL.4(names=c("Slope", "Lower", "Upper", "ED50")),
            data=E2_std)
plot(model2)

E2<-ED(model2, E2_samples$OD, type="absolute", display=F)
row.names(E2) <- E2_samples$Mouse_ID

points(y=E2_samples$OD,x=E2[,1],col="lightblue",pch=19,cex=2)
text(y =E2_samples$OD, x = E2[,1], labels=E2_samples$Mouse_ID, data=E2, cex=0.9, font=2)

E2 <- data.frame(E2)
colnames(E2)[1] <- "IFNy"
E2 <- dplyr::select(E2, IFNy)
setDT(E2, keep.rownames = TRUE)[]
colnames(E2)[1] <- "Mouse_ID"


write.csv(E2, "C:/Users/Luke Bednar/Mouse_Eimeria_Databasing/data/ELISAs/HZ19_CEWE_ELISA2_complete.csv")
# write.csv(E2, "C:/Users/Luke Bednar/Documents/Eimeria_Lab/data/3_recordingTables/P3_112019_Eim_CEWE_ELISAs/P3_112019_Eim_CEWE_ELISA2_complete.csv")


############### load in ELISA3

E3_std <- "https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/data/ELISAs/HZ19_CEWE_ELISA3_std.csv"
E3_std <- read.csv(text = getURL(E3_std))

E3_samples <- "https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/data/ELISAs/HZ19_CEWE_ELISA3_samples.csv"
E3_samples <- read.csv(text = getURL(E3_samples))

###### use drc to construct standard curve and pinpointprotein content

model3<-drm(OD~Conc,
            fct=LL.4(names=c("Slope", "Lower", "Upper", "ED50")),
            data=E3_std)
plot(model3)

E3<-ED(model3, E3_samples$OD, type="absolute", display=F)
row.names(E3) <- E3_samples$Mouse_ID

points(y=E3_samples$OD,x=E3[,1],col="lightblue",pch=19,cex=2)
text(y =E3_samples$OD, x = E3[,1], labels=E3_samples$Mouse_ID, data=E3, cex=0.9, font=2)

E3 <- data.frame(E3)
colnames(E3)[1] <- "IFNy"
E3 <- dplyr::select(E3, IFNy)
setDT(E3, keep.rownames = TRUE)[]
colnames(E3)[1] <- "Mouse_ID"

# write out ELISA 3
write.csv(E3, "C:/Users/Luke Bednar/Mouse_Eimeria_Databasing/data/ELISAs/HZ19_CEWE_ELISA3_complete.csv")





############### load in ELISA4

E4_std <- "https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/data/ELISAs/HZ19_CEWE_ELISA4_std.csv"
E4_std <- read.csv(text = getURL(E4_std))

E4_samples <- "https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/data/ELISAs/HZ19_CEWE_ELISA4_samples.csv"
E4_samples <- read.csv(text = getURL(E4_samples))

###### use drc to construct standard curve and pinpointprotein content

model4<-drm(OD~Conc,
            fct=LL.4(names=c("Slope", "Lower", "Upper", "ED50")),
            data=E4_std)
plot(model4)

E4<-ED(model4, E4_samples$OD, type="absolute", display=F)
row.names(E4) <- E4_samples$Mouse_ID

points(y=E4_samples$OD,x=E4[,1],col="lightblue",pch=19,cex=2)
text(y =E4_samples$OD, x = E4[,1], labels=E4_samples$Mouse_ID, data=E4, cex=0.9, font=2)

E4 <- data.frame(E4)
colnames(E4)[1] <- "IFNy"
E4 <- dplyr::select(E4, IFNy)
setDT(E4, keep.rownames = TRUE)[]
colnames(E4)[1] <- "Mouse_ID"

#write out ELISA 4
write.csv(E4, "C:/Users/Luke Bednar/Mouse_Eimeria_Databasing/data/ELISAs/HZ19_CEWE_ELISA4_complete.csv")


############## Load in ELISA 5
E5_std <- "https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/data/ELISAs/HZ19_CEWE_ELISA5_std.csv"
E5_std <- read.csv(text = getURL(E5_std))

E5_samples <- "https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/data/ELISAs/HZ19_CEWE_ELISA5_samples.csv"
E5_samples <- read.csv(text = getURL(E5_samples))

###### use drc to construct standard curve and pinpointprotein content

model5<-drm(OD~Conc,
            fct=LL.4(names=c("Slope", "Lower", "Upper", "ED50")),
            data=E5_std)
plot(model5)

E5<-ED(model5, E5_samples$OD, type="absolute", display=F)
row.names(E5) <- E5_samples$Mouse_ID

points(y=E5_samples$OD,x=E5[,1],col="lightblue",pch=19,cex=2)
text(y =E5_samples$OD, x = E5[,1], labels=E5_samples$Mouse_ID, data=E5, cex=0.9, font=2)

E5 <- data.frame(E5)
colnames(E5)[1] <- "IFNy"
E5 <- dplyr::select(E5, IFNy)
setDT(E5, keep.rownames = TRUE)[]
colnames(E5)[1] <- "Mouse_ID"

#write out ELISA 4
write.csv(E5, "C:/Users/Luke Bednar/Mouse_Eimeria_Databasing/data/ELISAs/HZ19_CEWE_ELISA5_complete.csv")

######################### merge all and write out
E <- rbind(HZ1, E2)
E <- rbind(E, E3)
E <- rbind(E, E4)
E <- rbind(E, E5)


write.csv(E, "C:/Users/Luke Bednar/Mouse_Eimeria_Databasing/data/ELISAs/HZ19_CEWE_ELISA.csv")

#################
