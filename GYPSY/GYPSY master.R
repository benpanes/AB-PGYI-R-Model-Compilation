rm(list=ls())
setwd("./GYPSY")
packages <- c("parallel","R.utils","dplyr","data.table")
lapply(packages,require,character.only=T)
input.dir <- "C:/Users/Ben/Desktop/Projects/R Compilation Code/2. Input Data" # Replace with directory for PGYI tables
cluster.num <- 20 # Replace with a number that's less than the number of cores your CPU has
scripts <- list.files("1. Scripts",full.names=T)
for(x in scripts){
  source(x, echo=verbose)
}

test