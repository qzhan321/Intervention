rm(list = ls())
suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
})
readDir <- "/project2/pascualmm/QZ/PhD/projects/intervention/natComRevision/files/figures/main/Fig5-290923_NYU/GhanaIndividualSurvey/Alltypes/"
prefix <- "survey"
nums <- 1:2

readDirEpi <- "/project2/pascualmm/QZ/PhD/projects/intervention/natComRevision/files/figures/main/Fig5-290923_NYU/GhanaSurvey/"
epi <- read.csv(paste0(readDirEpi, "Ghana_Survey_Merged_Epi_MOI_S1_S7_290923_NYU_130324.csv"), header = T, row.names = NULL)
epi$SeqID <- str_replace(epi$SeqID, "-", ".")
areas <- c("VeaGowrie", "Soe") 

saveDir <- "/project2/pascualmm/QZ/PhD/projects/intervention/natComRevision/files/deriveImmRate/mat/"
if (!dir.exists(saveDir)) {
  dir.create(saveDir)
}
  
for (i in 1:length(nums)) {
  num <- nums[i]
  file <- paste0(readDir, prefix, "_", num, ".csv")
  filecsv <- read.csv(file)
  locations <- data.frame("SeqID" = rownames(filecsv)) %>% left_join(epi %>% select(SeqID, CatchmentArea), by = "SeqID")
  stopifnot(nrow(filecsv) == nrow(locations))
  stopifnot(!any(is.na(locations$CatchmentArea)))
  
  matAll <- NULL
  indexSum <- 0
  for (j in 1:length(areas)) {
    area <- areas[j]
    location <- locations %>% filter(CatchmentArea %in% area)
    index <- which(rownames(filecsv) %in% location$SeqID)
    indexSum <- indexSum + length(index)
    mat <- filecsv[index, , drop=F]
    mat2 <- as.data.frame(colSums(mat))
    colnames(mat2) <- area
    matAll <- bind_cols(matAll, mat2)
  }
  stopifnot(indexSum == nrow(filecsv))
  save(matAll, file = paste0(saveDir, prefix, "_", num, ".RData"))
}
