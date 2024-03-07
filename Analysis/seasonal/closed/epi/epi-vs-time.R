rm(list = ls())
suppressPackageStartupMessages({
  library(RSQLite)
  library(reshape2)
  library(vegan)
  library(tidyverse)
  library(cowplot)
  library(gridExtra)
})
run <- 4
wd <- paste0("/scratch/midway2/qizhan/PhD/projects/intervention/simulation", run, "/actualRuns/")
seasonality <- "seasonal"
openness <- "closed"

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
IRSType <- "2yIRS"

if (IRSType == "10yIRS") {
  numsList <- list(101:110)
  nums_w_reps <- c(107:110)
} else {
  numsList <- list(1:14)
  nums_w_reps <- c(10:14)
}

prefix <- "sim"
IRSDurs <- rep(10, length(numsList))
preIRSs <- rep(200, length(numsList))
T_YEAR <- 360
samplingPeriod<-30

SIsummaryTableCombined <- NULL

tStarts <- preIRSs-3
for (i in 1:length(numsList)) {
  nums <- numsList[[i]]
  preIRS <- preIRSs[i]
  IRSDur <- IRSDurs[i]
  tStart <- tStarts[i]
  for (j in 1:length(nums)) {
    num <- nums[j]
    if (num %in% nums_w_reps) {
      reps <- 0:2
    } else {
      reps <- 0
    }
    for (r in reps) {
      sampleSqlFile<-paste(wd, seasonality, "/", openness, "/", prefix, "_", num,"/sqlitesDir/", prefix, "_", num, "_r", r, "_sd.sqlite",sep="")
      print(sampleSqlFile)
      db<-dbConnect(SQLite(),dbname = sampleSqlFile)
      
      sc<-"select * from summary"   
      summaryInfo<-fetchdb(db, sc)
      summaryInfo<-summaryInfo %>% mutate(num = num) %>% filter(time > tStart*T_YEAR) %>% mutate(year = ceiling(time/T_YEAR))
      summaryTable<-summaryInfo%>%mutate(Prevalence = n_infected/10000,
                                         MOI = n_infections/n_infected)
      eirs <- summaryInfo %>% 
        select(time, num, year, n_infected_bites) %>% group_by(num, year) %>% 
        summarise(EIR = sum(n_infected_bites)/10000)
      crs <- summaryInfo %>% 
        select(time, num, year, n_total_bites) %>% group_by(num, year) %>% 
        summarise(CR = sum(n_total_bites)/10000)
      
      summaryTable <- summaryTable %>% left_join(eirs, by = c("num", "year"))
      summaryTable <- summaryTable %>% left_join(crs, by = c("num", "year"))
      summaryTable$rep <- r
      
      summaryTable <- summaryTable %>% mutate(IRS = case_when(
        time<preIRS*T_YEAR ~ "preIRS",
        time>=preIRS*T_YEAR ~ paste0("I-", num - min(nums) + 1)
      ))
      
      summaryTableApp <- summaryTable[summaryTable$time==preIRS * T_YEAR,,drop=F]
      summaryTableApp$IRS = "preIRS"
      summaryTable <- rbind(summaryTable, summaryTableApp)
      if (IRSType == "2yIRS" & all(summaryTable$Prevalence > 0)) {
        SIsummaryTableCombined <- rbind(SIsummaryTableCombined, summaryTable)
      }
      dbDisconnect(db)
    }
  }
}

summaryTableCombined <- SIsummaryTableCombined
saveDir0 <- paste0("/project2/pascualmm/QZ/PhD/projects/intervention/files", run, "/actualRuns/epi/")
dir.create(saveDir0)
saveDir1 <- paste0(saveDir0, seasonality, "/")
dir.create(saveDir1)
saveDir2 <- paste0(saveDir1, openness, "/")
dir.create(saveDir2)

sfileName <- "summaryInfoTable"
save(summaryTableCombined, file = paste0(saveDir2, sfileName, "-", IRSType, ".RData"))
