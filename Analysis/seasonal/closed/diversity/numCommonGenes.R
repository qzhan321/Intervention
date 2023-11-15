rm(list = ls())
suppressPackageStartupMessages({
  library(RSQLite)
  library(reshape2)
  library(vegan)
  library(tidyverse)
  library(cowplot)
  library(gridExtra)
  library(data.table)
})
run <- 4
wd <- paste0("/project2/pascualmm/QZ/PhD/projects/intervention/files", run, "/actualRuns/infTablePresence/")
seasonality <- "seasonal"
openness <- "closed"
prefix<-"sim"
nums <- c(101:110, 115:123)
saveDir0 <- paste0("/project2/pascualmm/QZ/PhD/projects/intervention/files", run, "/actualRuns/diversity/numCommonGenes/")
dir.create(saveDir0)
saveDir1 <- paste0(saveDir0, seasonality, "/")
dir.create(saveDir1)
saveDir2 <- paste0(saveDir1, openness, "/")
dir.create(saveDir2)
fName <- "true"
saveDir3 <- paste0(saveDir2, fName, "/")
dir.create(saveDir3)

fetchdb<-function(dbname,query,numQuery = 20000000) {
  r<-dbSendQuery(conn=dbname, query)
  er<-dbFetch(r,numQuery)
  while(!dbHasCompleted(r)){
    er <- rbind(er, dbFetch(r, numQuery))
    print(nrow(er))
  }
  dbClearResult(r)
  return(er)
}
samplingPeriod<-30
preIRS <- 200
IRSDur <- 10
T_YEAR <- 360
postIRS <- 0
times <- c(seq(preIRS-3, preIRS+IRSDur+postIRS-1, 1)*T_YEAR+300)

nums_w_reps <- c(107:110,120:123)
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
    nCommonGenes <- sapply(timesInf, function(x){infStrain_pre_sub %>%
        filter(time==x) %>% group_by(gene_id) %>% summarise(n=n()) %>% ungroup() %>% select(n) %>% prop.table(.) %>% summarise(nCommonGenes = sum(.^2))})
    nCommonGenes <- 1/unlist(nCommonGenes)
    df <- data.frame("num" = num, "numCommon" = nCommonGenes, "time" = timesInf, "rep" = r)  
    
    save(df, file = paste0(saveDir3, "numCommonGenes-", num, "-r", r, ".RData"))
  }
}
