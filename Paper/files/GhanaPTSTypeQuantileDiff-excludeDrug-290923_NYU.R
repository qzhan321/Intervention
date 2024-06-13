rm(list = ls())
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyverse)
})
MOIdfAllGhana <- NULL
prefix <- "survey"
nums <- c(1,4,5)
date <- "290923_NYU"
readDir <- paste0("/project2/pascualmm/QZ/PhD/projects/intervention/natComRevision/files/figures/main/Fig5-", date, "/MOI/MOIEst/")
readDirEpi <- paste0("/project2/pascualmm/QZ/PhD/projects/intervention/natComRevision/files/figures/main/Fig5-", date, "/GhanaSurvey/")
epi <- read.csv(paste0(readDirEpi, "Ghana_Survey_Merged_Epi_MOI_S1_S7_290923_NYU_130324.csv"), header = T, row.names = NULL)
epi$SeqID <- str_replace(epi$SeqID, "-", ".")
epi <- epi %>% mutate(AgeGroups3 = case_when(
  (AgeGroups2 == "Children (1-5 years)") ~ "1-10",
  (AgeGroups2 == "Children (6-10 years)") ~ "1-10",
  (AgeGroups2 == "Adolescents (11-20 years)") ~ "11-20",
  (AgeGroups2 == "Adults (>20 years)") ~ ">20"
))
for (i in 1:length(nums)) {
  num <- nums[i]
  file <- paste0(readDir, prefix, "_", num, ".RData")
  load(file)
  MOIEstTemp <- outputList$indLevelMOI %>% left_join(epi %>% select("SeqID", "MalTreat2"), by = c("HostID" = "SeqID")) %>%
    filter(MalTreat2 == "No")
  MOIEst <- MOIEstTemp %>% mutate(survey = paste0(prefix, "_", num))
  MOIdf <- MOIEst %>% group_by(survey) %>% 
    summarise(PropMonoclonal = sum(MOI %in% c(1,2))/n())
  MOIdfAllGhana <- bind_rows(MOIdfAllGhana, MOIdf)
}

PTSdfAllGhana <- NULL
states <- c("pre IRS", "2y into IRS", "right post IRS (2015)")
file = paste0("/project2/pascualmm/QZ/PhD/projects/intervention/natComRevision/files/figures/main/Fig5-", date, "/PTSAgeGroups1-excludeDrug/PTS-BC-byAge-separateMOI-1-to-20-combineYoungAgeGroups-add2015")
load(file)
for (i in 1:length(states)) {
  s <- states[i]
  PTSdfAllSub <- PTSdfAll %>% filter(state == s) 
  PTSdfTemp <- PTSdfAllSub %>% group_by(state, age) %>%
    summarise(quantiles = seq(0,1,0.001), PTSQuantiles = quantile(PTS, probs = seq(0,1,0.001)))
  PTSdfTemp1 <- PTSdfTemp %>% filter(age == "1-10")
  PTSdfTemp2 <- PTSdfTemp %>% filter(age == ">20") 
  PTSdf <- PTSdfTemp1 %>% ungroup() %>% select(-age) %>% 
    left_join(PTSdfTemp2 %>% ungroup() %>% select(-age), by = c("state", "quantiles")) %>%
    mutate(PTSQuantilesDiff = PTSQuantiles.x - PTSQuantiles.y) %>% select(-c(PTSQuantiles.x, PTSQuantiles.y)) %>%
    mutate(state = case_when(state == "2y into IRS" ~ "IRS",
                             state == "pre IRS" ~ "Pre-IRS",
                             state == "right post IRS (2015)" ~ "Post-IRS"))
  PTSdfAllGhana <- bind_rows(PTSdfAllGhana, PTSdf)
}
PTSType <- "quantileDiff"
saveDir <- "/project2/pascualmm/QZ/PhD/projects/intervention/natComRevision/files/phaseDiagram/"
save(MOIdfAllGhana, PTSdfAllGhana, file = paste0(saveDir, "Ghana", "_PTSType_", PTSType, "-excludeDrug-", date, ".RData"))

