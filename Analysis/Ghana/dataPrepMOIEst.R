# Survey 1 (Oct 2012) >> Survey 1
# Survey 4 (Oct 2014) >> Survey 2
# Survey 5 (Oct 2015) >> Survey 3

rm(list=ls())
suppressPackageStartupMessages({
  library(dplyr)
})

saveDir <- "/project2/pascualmm/QZ/PhD/projects/intervention/natComRevision/files/figures/main/Fig5-290923_NYU/MOI/MOIEstInputs/"
if (!dir.exists(saveDir)) {
  dir.create(saveDir)
}

readDir <- "/project2/pascualmm/QZ/PhD/projects/intervention/natComRevision/files/figures/main/Fig5-290923_NYU/GhanaIndividualSurvey/UpsBC/"
prefix <- "survey"
nums  <- c(1,4,5)
files <- paste0(prefix, "_", nums, ".csv")

for (i in 1:length(files)) {
  file <- files[i]
  num <- nums[i]
  dataT <- read.csv(paste0(readDir, file), header = T, row.names = 1)
  print(dim(dataT))
  table(dataT > 1)
  isolateSizes <- rowSums(dataT > 0)
  MOIInput <- data.frame("HostID" = names(isolateSizes), "DBLa_upsBC_rep_size" = unname(isolateSizes))
  write.csv(MOIInput, file = paste0(saveDir, prefix, "_", num, ".csv"), row.names = F)
}

