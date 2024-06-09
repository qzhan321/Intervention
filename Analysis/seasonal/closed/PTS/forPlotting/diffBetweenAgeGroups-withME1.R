rm(list = ls())
suppressPackageStartupMessages({
  library(dplyr)
})
run <- 4
readDirPTS <- paste0("/project2/pascualmm/QZ/PhD/projects/intervention/files", run, "/actualRuns/PTS/")
seasonality <- "seasonal"
openness <- "closed"
nums <- c(101:109, 115:121)
prefix <- "sim"
IRSDur <- 10
preIRS <- 200
T_YEAR <- 360
layers <- c(preIRS - 1, preIRS + 1)*T_YEAR + 300
state <- "withME"

readDirMOI <- paste0("/project2/pascualmm/QZ/PhD/projects/intervention/files", run, "/actualRuns/MOI/")

saveDir0 <- paste0("/project2/pascualmm/QZ/PhD/projects/intervention/files", run, "/actualRuns/PTSAgeDiff1/")
dir.create(saveDir0)
saveDir1 <- paste0(saveDir0, seasonality, "/")
dir.create(saveDir1)
saveDir2 <- paste0(saveDir1, openness, "/")
dir.create(saveDir2)

MOI_considered <- 1:20
ageLow <- 10
ageHigh <- 20

types <- c("dir", "nondir")
nums_w_reps <- NULL
for (i in 1:length(nums)) {
  num <- nums[i]
  print(num)
  if (num %in% nums_w_reps) {
    reps <- 0:2
  } else {
    reps <- 0
  }
  
  load(paste0(readDirMOI, seasonality, "/", openness, "/", state, "/MOI_", num, ".RData"))
  for (r in reps) {
    for (m in 1:length(types)) {
      type <- types[m]
      saveDir3 <- paste0(saveDir2, type, "/")
      dir.create(saveDir3)
      saveDir4 <- paste0(saveDir3, state, "/")
      dir.create(saveDir4)
      
      PTSAll <- NULL
      load(paste0(readDirPTS, seasonality, "/", openness, "/", type, "/", state, "/", prefix, "_", num, "_r", r, ".RData"))
      
      if (type == "dir") {
        PTSMatrixList <- PTSDirMatrixList
      } else {
        PTSMatrixList <- PTSNonDirMatrixList
      }
      for (j in 1:length(layers)) {
        layer <- layers[j]
        print(layer)
        MOIsub <- MOIAll %>% filter(time == layer, rep == r) %>% mutate(host_id = paste0(time, "-", host_id))
        timeLabel <- paste0(layer, "-", layer)
        index <- which(names(PTSMatrixList) == timeLabel)
        PTS <- PTSMatrixList[[index]]
        index <- which(names(ageList) == timeLabel)
        Age <- ageList[[index]]
        
        for (k in 1:length(MOI_considered)) {
          MOI_considered_sub <- MOI_considered[k]
          print(MOI_considered_sub)
          MOIAge <- MOIsub %>% left_join(Age, by = "host_id")
          hosts1 <- MOIAge %>% filter(age <= ageLow, MOI %in% MOI_considered_sub)
          if (nrow(hosts1) > 1) {
            hosts1 <- unique(hosts1$host_id)
            print(length(hosts1))
            PTS1 <- PTS[hosts1, hosts1]
            pts1 <- data.frame("PTS" = PTS1[row(PTS1)!=col(PTS1)], "time" = layer, "num" = num, "ageGroup" = "1-10", "MOI" = MOI_considered_sub)
            PTSAll <- rbind(PTSAll, pts1)
          }
          hosts2 <- MOIAge %>% filter(age >= ageHigh, MOI %in% MOI_considered_sub)
          if (nrow(hosts2) > 1) {
            hosts2 <- unique(hosts2$host_id)
            print(length(hosts2))
            PTS2 <- PTS[hosts2, hosts2]
            pts2 <- data.frame("PTS" = PTS2[row(PTS2)!=col(PTS2)], "time" = layer, "num" = num, "ageGroup" = ">=20", "MOI" = MOI_considered_sub)
            PTSAll <- rbind(PTSAll, pts2)
          }
        }
      } 
      save(PTSAll, file = paste0(saveDir4, "PTSAcrossAgeGroups_", num, "_r", r, ".RData"))
    }
  }
}

