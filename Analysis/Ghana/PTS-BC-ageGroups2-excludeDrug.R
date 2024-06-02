rm(list = ls())
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyverse)
})
pairwiseTypeSharing<-function(mat){
  newmat<-tcrossprod(mat>0)
  newmat<-newmat/rowSums(mat>0)
  return(newmat)
}

readDir <- "/project2/pascualmm/QZ/PhD/projects/FOI/files/Ghana/surveys/BC/"
surveys <- c("survey_1", "survey_4", "survey_5")
files <- paste0(surveys, ".csv")
years <- c(2012, 2014, 2015)

readDirMOI <- "/project2/pascualmm/QZ/PhD/projects/intervention/natComRevision/files/figures/main/Fig5/MOI/MOIEst/"
readDirEpi <- "/home/qizhan/others/PhD/projects/intervention/natComRevision/utils/"
epi <- read.csv(paste0(readDirEpi, "Ghana_Survey_Merged_Epi_MOI_S1_S7_070721_UChicago_080822.csv"), header = T, row.names = 1)
epi$SeqID <- str_replace(epi$SeqID, "-", ".")
epi <- epi %>% mutate(AgeGroups3 = case_when(
  (AgeGroups2 == "Children (1-5 years)") ~ "1-10",
  (AgeGroups2 == "Children (6-10 years)") ~ "1-10",
  (AgeGroups2 == "Adolescents (11-20 years)") ~ "11-20",
  (AgeGroups2 == "Adults (>20 years)") ~ ">20"
))

PTSdfAll <- NULL
MOIInclude <- 1:20
for (i in 1:length(files)) {
  survey <- surveys[i]
  load(paste0(readDirMOI, survey, ".RData"))
  MOI <- outputList$indLevelMOI
  MOIAge <- MOI %>% left_join(epi %>% select("SeqID", "AgeGroups3", "MalTreat2"), by = c("HostID" = "SeqID"))
  file <- files[i]
  isolateType <- read.csv(paste0(readDir, file), header = T, row.names = 1)
  PTS <- pairwiseTypeSharing(isolateType)
  
  stopifnot(identical(as.character(MOI$HostID), rownames(PTS)))
  
  ageGroupV <- MOIAge$AgeGroups3
  ageGroup <- unique(ageGroupV)
  for (j in 1:length(ageGroup)) {
    ageGroup_single <- ageGroup[j]
    print(ageGroup_single)
    for (k in 1:length(MOIInclude)) {
      MOIInclude_single <- MOIInclude[k]
      hostsSub <- MOIAge %>% filter(AgeGroups3 == ageGroup_single, MalTreat2 == "No", MOI %in% MOIInclude_single) %>% .$HostID
      if (length(hostsSub) > 1) {
        PTSsub <- PTS[c(hostsSub), c(hostsSub)]
        PTSSubdf <- data.frame("PTS" = PTSsub[row(PTSsub)!=col(PTSsub)], "age" = ageGroup_single, "survey" = survey, "MOI" = MOIInclude_single)
        PTSdfAll <- rbind(PTSdfAll, PTSSubdf)
      }
    }
  }
}
PTSdfAll <- PTSdfAll %>% mutate(year = case_when((survey == "survey_1")~2012,
                          (survey == "survey_4")~2014,
                          (survey == "survey_5")~2015))
PTSdfAll <- PTSdfAll %>% mutate(state = case_when((survey == "survey_1")~"pre IRS",
                                           (survey == "survey_4")~ "2y into IRS",
                                           (survey == "survey_5")~"right post IRS (2015)"))
PTSdfAll$age <- factor(PTSdfAll$age, levels = c("1-10", "11-20", ">20"))
PTSdfAll$state <- factor(PTSdfAll$state, levels = c("pre IRS", "2y into IRS", "right post IRS (2015)"))
saveDir <- "/project2/pascualmm/QZ/PhD/projects/intervention/natComRevision/files/figures/main/Fig5/PTSAgeGroups1-excludeDrug/"
if (!dir.exists(saveDir)) {
  dir.create(saveDir)
}
save(PTSdfAll, file = paste0(saveDir, "PTS-BC-byAge-separateMOI-1-to-20-combineYoungAgeGroups-add2015-2"))
