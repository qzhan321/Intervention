rm(list = ls())
suppressPackageStartupMessages({
  library(dplyr)
  library(RSQLite)
})

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
run <- 4
wd <- paste0("/scratch/midway2/qizhan/PhD/projects/intervention/simulation", run, "/actualRuns/")
seasonality <- "seasonal"
openness <- "closed"
saveDir0 <- paste0("/project2/pascualmm/QZ/PhD/projects/intervention/files", run, "/actualRuns/infectionDur/")
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

IRSTypes <- c("10yIRS", "2yIRS")
numsList <- list(c(101:110, 115:123),
                 c(1:14,15:23))
nums_w_reps <- c(107:110,120:123,
                 10:14,20:23)
for (i in 1:length(seasonality)) {
  s <- seasonality[i]
  for (j in 1:length(openness)) {
    o <- openness[j]
    for (k in c(1)) {
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
          sampleSqlFile <- paste(wd, s, "/", o, "/", prefix, "_", num,"/sqlitesDir/", prefix, "_", num, "_r", r, "_sd.sqlite",sep="")
          db<-dbConnect(SQLite(),dbname = sampleSqlFile)
          sc<-"select * from sampled_duration"   
          durInfo <- fetchdb(db, sc)
          durInfo <- durInfo %>% filter((time >= (preIRS - 3) * T_YEAR & time < preIRS * T_YEAR)| 
                                          (time >= (preIRS + 0) * T_YEAR & time - duration > preIRS * T_YEAR + 14 & time < (preIRS + 3) * T_YEAR)|
                                          (time >= (preIRS + 7) * T_YEAR))
          durTable <- durInfo 
          durTable$rep <- r
          durTable$num <- num
          durTable <- durTable %>% mutate(IRS = case_when(
            time < preIRS * T_YEAR ~ "Pre-IRS",
            time >= preIRS*T_YEAR & time < (preIRS + 3) * T_YEAR ~ paste0("Early Stage of IRS"),
            time >= (preIRS + 7)*T_YEAR ~ paste0("Late Stage of IRS")
          ))
          
          dbDisconnect(db)
          
          save(durTable, file = paste0(saveDir3, prefix, "_", num, "_rep_", r, "_", IRSType, ".RData"))
        }
      }
    }
  }
}