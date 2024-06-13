rm(list = ls())
suppressPackageStartupMessages({
  library(dplyr)
})

source("/home/qizhan/others/PhD/projects/intervention/natComRevision/analysis/files/deriveImmRate/utils.R")
readDir <- "/project2/pascualmm/QZ/PhD/projects/intervention/natComRevision/files/deriveImmRate/mat/"
prefix <- "survey"
nums <- 1:2

for (i in 1:length(nums)) {
  num <- nums[i]
  file <- paste0(readDir, prefix, "_", num, ".RData")
  load(file)
  mut <- 5.64971751*10^(-7)*60*59/2/60 # primarily ectopic recombination rate, measured only in vitro. Alternatively, mut = 5.64971751*10^(-7)*60*59/2, but then mut/D/60 (migration brings in 60 genes at once).  
  D <- estD(matAll)
  mig <- mut/D 
  print(mig)
}
