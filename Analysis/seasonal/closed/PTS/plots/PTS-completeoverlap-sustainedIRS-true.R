rm(list = ls())
suppressPackageStartupMessages({
  library(RSQLite)
  library(reshape2)
  library(vegan)
  library(tidyverse)
  library(cowplot)
  library(gridExtra)
  library(zoo)
})
run <- 4
readDir <- paste0("/project2/pascualmm/QZ/PhD/projects/intervention/files", run, "/actualRuns/PTS/")
state <- "true"
types <- c("dir", "nondir")
seasonality <- "seasonal"
openness <- "closed"
nums <- c(101:109)
prefix <- "sim"
IRSDur <- 10
preIRS <- 200
postIRS <- 0
T_YEAR <- 360
timesSelected <- c(preIRS - 1, preIRS + 1, preIRS + 2)*T_YEAR + 300

saveDir0 <- paste0("/project2/pascualmm/QZ/PhD/projects/intervention/files", run, "/actualRuns/PTSCompleteOverlap/")
dir.create(saveDir0)
saveDir1 <- paste0(saveDir0, seasonality, "/")
dir.create(saveDir1)
saveDir2 <- paste0(saveDir1, openness, "/")
dir.create(saveDir2)

reps <- 0
for (i in 1:length(nums)) {
  num <- nums[i]
  print(num)
  
  for (r in reps) {
    for (j in 1:length(types)) {
      type <- types[j]
      
      saveDir3 <- paste0(saveDir2, type, "/")
      dir.create(saveDir3)
      saveDir4 <- paste0(saveDir3, state, "/")
      dir.create(saveDir4)
      
      load(paste0(readDir, seasonality, "/", openness, "/", type, "/", state, "/", prefix, "_", num, "_r", r, ".RData"))
      
      if (type == "dir") {
        PTSMatrixList <- PTSDirMatrixList
      } else {
        PTSMatrixList <- PTSNonDirMatrixList
      }
      
      dfAll <- NULL
      for (k in 1:length(timesSelected)) {
        t <- timesSelected[k]
        print(t)
        index <- which(names(PTSMatrixList) == paste0(t, "-", t))
        if (length(index) == 1) {
          PTS <- PTSMatrixList[[index]]
          if (nrow(PTS) > 1) {
            df <- data.frame("PTS" = PTS[row(PTS)!=col(PTS)], "time" = t, "num" = num)
            dfAll <- rbind(dfAll, df)
          }
        }
      }
      save(dfAll, file = paste0(saveDir4, "PTS_selectedTimes_", num, "_r", r, ".RData"))
    }
  }
}
