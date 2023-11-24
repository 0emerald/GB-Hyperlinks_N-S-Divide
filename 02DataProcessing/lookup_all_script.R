library(tidyverse)
library(stringr)
library(data.table)
library(plyr)
library(dplyr)
library(stats) 
library(urltools)
library(utils)
library(feather)

################################################
### SCRIPT TO OUTPUT A lookup_all FILE
### LOOKUP FILE THAT IS NOT FILTERED BY POSTCODE
################################################

# This is the file of the columns pc and host extracted from all.1.csv (1.9GB) (20s read-in)
hosts_pc <- read.csv("hosts-pc.csv")

head(hosts_pc, 100)

# Unique values in hosts_BS
# I hope this looks for unique pairs -- it does
hosts_pc_unique <- unique(hosts_pc)

# 553 kB files
lookup_all <- write.csv(hosts_pc_unique, "lookup_all.csv", row.names = FALSE, col.names = TRUE)
# outputs filel lookup_all (1GB)