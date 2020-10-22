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