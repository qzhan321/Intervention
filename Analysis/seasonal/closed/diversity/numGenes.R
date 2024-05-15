rm(list = ls())
suppressPackageStartupMessages({
  library(dplyr)
})
run <- 4
wd <- paste0("/project2/pascualmm/QZ/PhD/projects/intervention/files", run, "/actualRuns/infTablePresence/")
seasonality <- "seasonal"
openness <- "closed"
prefix<-"sim"
nums <- 101:110
saveDir0 <- paste0("/project2/pascualmm/QZ/PhD/projects/intervention/files", run, "/actualRuns/diversity/numGenes/")
dir.create(saveDir0)
saveDir1 <- paste0(saveDir0, seasonality, "/")
dir.create(saveDir1)
saveDir2 <- paste0(saveDir1, openness, "/")
dir.create(saveDir2)
fName <- "true"
saveDir3 <- paste0(saveDir2, fName, "/")
dir.create(saveDir3)

samplingPeriod<-30
preIRS <- 200
IRSDur <- 10
T_YEAR <- 360
postIRS <- 0
times <- c(seq(preIRS-3, preIRS+IRSDur+postIRS-1, 1)*T_YEAR+300)

nums_w_reps <- 107:110
for (a in 1:length(nums)) {
  num <- nums[a]
  print(num)
  if (num %in% nums_w_reps) {
    reps <- 0:2
    files <- paste0(prefix, "_", num, "_r", reps, ".RData")
  } else {
    reps <- 0
    files <- paste0(prefix, "_", num, "_r", reps, ".RData")
  }
  for (b in 1:length(files)) {
    file <- files[b]
    r <- reps[b]
    load(paste0(wd, seasonality, "/", openness, "/", file))
    infStrain_pre_sub <- infStrain_pre %>% filter(is_functional > 0, time %in% times)
    
    timesInf <- sort(unique(infStrain_pre_sub$time))
    nGenes <- infStrain_pre_sub %>% group_by(time) %>% summarise(numGenes=length(unique(gene_id))) 
    df <- nGenes %>% mutate("num" = num, "rep" = r)  
    save(df, file = paste0(saveDir3, "numGenes-", num, "-r", r, ".RData"))
  }
}
