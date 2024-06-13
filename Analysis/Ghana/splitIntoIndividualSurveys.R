# 2012: S1, wet
# 2013: S2, dry
# 2014: S3, dry
# 2014: S4, wet
# 2015: S5, wet
# 2016: S6, dry
# 2017: S7, wet
# SMC starts from 2016 wet season, so S7 is under its impact.
# pre-IRS: S1, S2; IRS: S3, S4; right post-IRS (IRS sort of): S5, S6; post-IRS: 2017
rm(list=ls())
suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
})

readDir <- "/project2/pascualmm/QZ/PhD/projects/intervention/natComRevision/files/figures/main/Fig5-290923_NYU/GhanaSurvey/"
epi <- read.csv(paste0(readDir, "Ghana_Survey_Merged_Epi_MOI_S1_S7_290923_NYU_130324.csv"), header = T, row.names = NULL)
epi$SeqID <- str_replace(epi$SeqID, "-", ".")
nums <- 1:7
types <- c("Alltypes", "UpsBC", "UpsA")
prefix <- "survey"
saveDir0 <- "/project2/pascualmm/QZ/PhD/projects/intervention/natComRevision/files/figures/main/Fig5-290923_NYU/GhanaIndividualSurvey/"
if (!dir.exists(saveDir0)) {
  dir.create(saveDir0)
}
for (i in 1:length(types)) {
  type <- types[i]
  saveDir1 <- paste0(saveDir0, type, "/")
  if (!dir.exists(saveDir1)) {
    dir.create(saveDir1)
  }
  surveysAll <- read.csv(paste0(readDir, "Ghana_Survey_Final_S1_S7_MatrixPTS_", type, "_290923.csv"), header = T, row.names = 1)
  print(dim(surveysAll))
  for (j in 1:length(nums)) {
    num <- nums[j]
    isolatesSurvey <- epi %>% filter(Survey == num) %>% .$SeqID
    isolatesSurveySub <- intersect(colnames(surveysAll), isolatesSurvey)
    survey <- t(surveysAll[, isolatesSurveySub])
    print(dim(survey))
    write.table(survey, file = paste0(saveDir1, prefix, "_", num, ".csv"), sep = ",", row.names = T, col.names = T)
  }
}
