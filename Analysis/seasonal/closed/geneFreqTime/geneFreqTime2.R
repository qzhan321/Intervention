rm(list = ls())
suppressPackageStartupMessages({
  library(RSQLite)
  library(reshape2)
  library(vegan)
  library(tidyverse)
  library(cowplot)
  library(gridExtra)
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

IRSTypes <- c("10yIRS", "2yIRS")
numsList <- list(c(101:110),
                 c(1:14))
nums_w_reps <- c(107:110,
                 10:14)
for (i in 1:length(seasonality)) {
  s <- seasonality[i]
  for (j in 1:length(openness)) {
    o <- openness[j]
    for (k in c(2)) {
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
          
          sc<-paste0("select time, host_id, strain_id from sampled_infections")
          sampledInf<-fetchdb(db,sc)
          infStrain<-sampledInf %>% mutate(uniqStrain = 1:nrow(sampledInf))
          
          strain_ids <- unique(infStrain$strain_id)
          sc<-paste0("select id, ind, gene_id from sampled_strains where id IN (", paste(noquote(strain_ids), collapse = ","), ")")
          strainInfo<-fetchdb(db, sc)
          colnames(strainInfo)[1]<-"strain_id"
          
          sc<-"select * from genes"
          genes<-fetchdb(db, sc)
          colnames(genes)[1]<-"gene_id"
          
          sc<-"select * from alleles"
          alleleInfo<-fetchdb(db, sc)
          geneAllele<-dcast(alleleInfo, gene_id~locus, value.var="allele")
          colnames(geneAllele)[2:ncol(geneAllele)]<-paste("l",1:(ncol(geneAllele)-1),sep="")
          
          infStrain<-left_join(infStrain, strainInfo, by="strain_id")
          infStrain<-left_join(infStrain, genes, by="gene_id")
          infStrain<-left_join(infStrain, geneAllele, by="gene_id")
          infStrain <- infStrain %>% mutate(gene_id = paste0(l1, "-", l2))
          
          df1 <- infStrain %>% filter(is_functional == 1) %>% dplyr::group_by(time, gene_id) %>% summarise(n=n())
          df2 <- infStrain %>% filter(is_functional == 1) %>% dplyr::group_by(time, gene_id) %>% summarise(n=n()) %>% group_by(time) %>% summarise(nTotal = sum(n))
          df12 <- df1 %>% left_join(df2, by = "time")
          df12 <- df12 %>% mutate(p = n/nTotal)
          
          df12 <- df12 %>% mutate(num = num, r = r)   
          
          save(df12, file = paste0(saveDir3, prefix, "_", num, "_rep_", r, "_", IRSType))
        }
      }
    }
  }
}
