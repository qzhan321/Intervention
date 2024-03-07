rm(list = ls())
suppressPackageStartupMessages({
  library(RSQLite)
  library(reshape2)
  library(vegan)
  library(tidyverse)
  library(cowplot)
  library(scales)
})
run <- 4
wd <- paste0("/scratch/midway2/qizhan/PhD/projects/intervention/simulation", run, "/actualRuns/")
seasonality <- "seasonal"
openness <- "closed"
saveDir0 <- paste0("/project2/pascualmm/QZ/PhD/projects/intervention/files", run, "/actualRuns/PTS/")
dir.create(saveDir0)
saveDir1 <- paste0(saveDir0, seasonality, "/")
dir.create(saveDir1)
saveDir2 <- paste0(saveDir1, openness, "/")
dir.create(saveDir2)
saveDir3 <- paste0(saveDir2, "dir/")
dir.create(saveDir3)
saveDir4 <- paste0(saveDir2, "nondir/")
dir.create(saveDir4)
fName <- "true"
saveDir5 <- paste0(saveDir3, fName, "/")
dir.create(saveDir5)
saveDir6 <- paste0(saveDir4, fName, "/")
dir.create(saveDir6)

readDir <- paste0("/project2/pascualmm/QZ/PhD/projects/intervention/files", run, "/actualRuns/infTablePresence/")
prefix<-"sim"
nums <- c(101:110)

PTSDir<-function(mat){
  newmat<-tcrossprod(mat>0)
  newmat<-newmat/rowSums(mat>0)
  return(newmat)
}

PTSNonDir<-function(mat){
  newmat1<-2*tcrossprod(mat>0)
  newmat2<-matrix(rep(rowSums(mat), nrow(mat)), ncol = nrow(mat), byrow = F)
  newmat<-newmat1/(newmat2 + t(newmat2))
  return(newmat)
}

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

tempfunction <- function(x,infStrain) {
  temp1 <- infStrain %>% filter(time %in% x) 
  paraCb <- acast(temp1, host_id~gene_id, length, value.var = "host_id")
  outMat1<-PTSDir(as.matrix(paraCb))
  outMat2<-PTSNonDir(as.matrix(paraCb))
  return(list("dir" = outMat1, "nondir" = outMat2))
}

IRSDur <- 10
preIRS <- 200
T_YEAR <- 360
postIRS <- 0
times <- seq(preIRS - 2, preIRS + 3, 1)*T_YEAR + 300
nums_w_reps <- NULL
for (a in 1:length(nums)) {
  num <- nums[a]
  print(num)
  if (num %in% nums_w_reps) {
    reps <- 0:2
  } else {
    reps <- 0
  }
  for (r in reps) {
    sampleSqlFile<-paste(wd, seasonality, "/", openness, "/", prefix, "_", num, "/sqlitesDir/", prefix, "_", num, "_r0_sd.sqlite",sep="")
    print(sampleSqlFile)
    db<-dbConnect(SQLite(),dbname = sampleSqlFile)
    
    sc <- "select id, birth_time, death_time from hosts"
    hosts <- fetchdb(db, sc)
    colnames(hosts)[1] <- "host_id"
    
    load(paste0(readDir, seasonality, "/", openness, "/", prefix, "_", num, "_r", r, ".RData"))
    infStrain <- infStrain_pre %>% filter(is_functional == 1)
    
    infStrain <- infStrain %>% left_join(hosts, by = "host_id")
    infStrain <- infStrain %>% mutate(age = (time-birth_time)/T_YEAR,
                                      host_id = paste0(time, "-", host_id))
    dbDisconnect(db)
    
    times_len <- length(times)
    PTSDirMatrixList <- vector("list", (1 + times_len)/2*times_len)
    PTSNonDirMatrixList <- vector("list", (1 + times_len)/2*times_len)
    ageList <- vector("list", (1 + times_len)/2*times_len)
    counter <- 1
    for (m in 1:length(times)) {
      t1 <- times[m]
      infStrainSub <- infStrain %>% filter(time %in% c(t1))
      hostIds1 <- unique(infStrainSub$host_id)
      for (n in m:length(times)) {
        t2 <- times[n]
        infStrainSub <- infStrain %>% filter(time %in% c(t2)) 
        hostIds2 <- unique(infStrainSub$host_id)
        if (length(hostIds1) == 0 | length(hostIds2) == 0) {
          names(PTSDirMatrixList)[counter] <- paste0(t1, "-", t2)
          names(PTSNonDirMatrixList)[counter] <- paste0(t1, "-", t2)
          counter <- counter + 1
          break
        } else {
          PTS <- tempfunction(c(t1, t2), infStrain)
          PTS_dir <- PTS[["dir"]]
          PTS_nondir <- PTS[["nondir"]]
          stopifnot(hostIds1 %in% rownames(PTS_dir))
          stopifnot(hostIds2 %in% colnames(PTS_dir))
          
          stopifnot(hostIds1 %in% rownames(PTS_nondir))
          stopifnot(hostIds2 %in% colnames(PTS_nondir))
          
          PTSDirSub <- PTS_dir[hostIds1, hostIds2]
          PTSNonDirSub <- PTS_nondir[hostIds1, hostIds2]
          
          PTSDirMatrixList[[counter]] <- PTSDirSub
          PTSNonDirMatrixList[[counter]] <- PTSNonDirSub
          
          names(PTSDirMatrixList)[counter] <- paste0(t1, "-", t2)
          names(PTSNonDirMatrixList)[counter] <- paste0(t1, "-", t2)
          
          if (t1 == t2) {
            Age <- infStrain %>% filter(host_id %in% hostIds1) %>% select(host_id, age) %>% distinct(host_id, age)
          } else {
            Age <- infStrain %>% filter(host_id %in% c(hostIds1, hostIds2)) %>% select(host_id, age) %>% distinct(host_id, age)
          }
          ageList[[counter]] <- Age
          names(ageList)[counter] <- paste0(t1, "-", t2)
          counter <- counter + 1
        }
      }
    }
    save(PTSDirMatrixList, ageList, file = paste0(saveDir5, prefix, "_", num, "_r", r, ".RData"))
    save(PTSNonDirMatrixList, ageList, file = paste0(saveDir6, prefix, "_", num, "_r", r, ".RData"))
  }
}
