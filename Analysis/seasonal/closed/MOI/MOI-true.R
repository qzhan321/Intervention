rm(list=ls())
suppressPackageStartupMessages({
  library(dplyr)
  library(RSQLite)
})
run <- 4
wd <- paste0("/scratch/midway2/qizhan/PhD/projects/intervention/simulation", run, "/actualRuns/")
seasonality <- "seasonal"
openness <- "closed"
saveDir0 <- paste0("/project2/pascualmm/QZ/PhD/projects/intervention/files", run, "/actualRuns/MOI/")
dir.create(saveDir0)
saveDir1 <- paste0(saveDir0, seasonality, "/")
dir.create(saveDir1)
saveDir2 <- paste0(saveDir1, openness, "/")
dir.create(saveDir2)
fName <- "true"
saveDir3 <- paste0(saveDir2, fName, "/")
dir.create(saveDir3)

readDir <- paste0("/project2/pascualmm/QZ/PhD/projects/intervention/files", run, "/actualRuns/infTablePresence/")
prefix <- "sim"
nums <- 101:110

preIRS <- 200
IRSDur <- 10
T_YEAR <- 360
layers <- c(preIRS - 1, preIRS + 1)*T_YEAR  + 300

for (i in 1:length(nums)) {
  num <- nums[i]
  load(paste0(readDir, seasonality, "/", openness, "/", prefix, "_", num, "_r0.RData"))
  infStrain_pre_sub <- infStrain_pre %>% filter(time %in% layers)
  MOIAll <- infStrain_pre_sub %>% group_by(time, host_id) %>% summarise(MOI = length(unique(uniqStrain))) %>% mutate(num = nums[i])
  save(MOIAll, file = paste0(saveDir3, "MOI_", num, ".RData"))
}
