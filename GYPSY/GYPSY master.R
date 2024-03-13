rm(list=ls())
setwd("./GYPSY")
packages <- c("parallel","R.utils","dplyr","dtplyr","data.table")
lapply(packages,require,character.only=T)
input.dir <- "H:/Shared drives/Growth & Yield Lab/Data Sets/PGYI/2023-09-08 PGYI Export" # Replace with directory for PGYI tables
cluster.num <- 20 # Replace with a number that's less than the number of cores your CPU has
#scripts <- list.files("GYPSY subscripts",full.names=T)
#for(x in scripts){
#  source(x, echo=verbose)
#}