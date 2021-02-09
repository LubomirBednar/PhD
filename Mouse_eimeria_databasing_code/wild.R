# HMHZ wild mice 
library(ggplot2)
library(Rmisc)
library(httr)
library(RCurl)
library(tidyverse)

# load in past 4 years
HZ19 <- read.csv(text = getURL(""))
HZ18 <- read.csv(text = getURL(""))
HZ17 <- read.csv(text = getURL(""))
HZ16 <- read.csv(text = getURL(""))
# load in wild_immuno_long
wild_immuno_long <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/data/HZ19_immuno_long.csv"))
