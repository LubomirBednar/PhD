# clean up HZ19 data
library(dplyr)

HZ <- read.csv("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Field/master/data/Field_data/HZ10-19_Genotypes.csv")


## Calculate the hybrid index
get.HI <- function (x){
  dom <- nchar(x) - nchar(gsub("d", "", x))
  mus <- nchar(x) - nchar(gsub("m", "", x))
  mus/(mus + dom)
}

HIHZ <- HZ %>%
  select(c(2,13:20, 30:35))
lapply(HIHZ[2:15], get.HI)
