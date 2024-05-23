rm(list = ls())
suppressPackageStartupMessages({
  library(dplyr)
})

run <- 4
seasonality <- "seasonal"
openness <- "closed"
saveDir0 <- paste0("/project2/pascualmm/QZ/PhD/projects/intervention/files", run, "/actualRuns/geneFreqTime/")
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

preIRS <- 200
IRSDur <- 10
T_YEAR <- 360

IRSTypes <- "10yIRS"
numsList <- list(101:110)
nums_w_reps <- 107:110
for (i in 1:length(seasonality)) {
  s <- seasonality[i]
  for (j in 1:length(openness)) {
    o <- openness[j]
    for (k in 1:length(IRSTypes)) {
      IRSType <- IRSTypes[k]
      nums <- numsList[[k]]
      
      for (p in 1:length(nums)) {
        num <- nums[p]
        if (num %in% nums_w_reps) {
          reps <- 0:2
        } else {
          reps <- 0
        }
        for (r in reps) {
          load(paste0(readDir, s, "/", o, "/", prefix, "_", num, "_r", r, ".RData"))
          df1 <- infStrain_pre %>% filter(is_functional == 1) %>% dplyr::group_by(time, gene_id) %>% summarise(n=n())
          df2 <- infStrain_pre %>% filter(is_functional == 1) %>% dplyr::group_by(time, gene_id) %>% summarise(n=n()) %>% group_by(time) %>% summarise(nTotal = sum(n))
          df12 <- df1 %>% left_join(df2, by = "time")
          df12 <- df12 %>% mutate(p = n/nTotal)
          
          df12 <- df12 %>% mutate(num = num, r = r)  
          
          save(df12, file = paste0(saveDir3, prefix, "_", num, "_rep_", r, "_", IRSType))
        }
      }
    }
  }
}
