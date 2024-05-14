rm(list=ls())
suppressPackageStartupMessages({
  library(dplyr)
  library(RSQLite)
})
run <- 4
wd <- paste0("/scratch/midway2/qizhan/PhD/projects/intervention/simulation", run, "/actualRuns/")
seasonality <- "seasonal"
openness <- "closed"
saveDir0 <- paste0("/project2/pascualmm/QZ/PhD/projects/intervention/files", run, "/actualRuns/MOI/")
dir.create(saveDir0)
saveDir1 <- paste0(saveDir0, seasonality, "/")
dir.create(saveDir1)
saveDir2 <- paste0(saveDir1, openness, "/")
dir.create(saveDir2)
fName <- "withME"
saveDir3 <- paste0(saveDir2, fName, "/")
dir.create(saveDir3)

sampledHostsSaveDir0 <- paste0("/project2/pascualmm/QZ/PhD/projects/intervention/files", run, "/actualRuns/MicroscopyPosHosts/")
dir.create(sampledHostsSaveDir0)
sampledHostsSaveDir1 <- paste0(sampledHostsSaveDir0, seasonality, "/")
dir.create(sampledHostsSaveDir1)
sampledHostsSaveDir2 <- paste0(sampledHostsSaveDir1, openness, "/")
dir.create(sampledHostsSaveDir2)
sampledHostsSaveDir3 <- paste0(sampledHostsSaveDir2, fName, "/")
dir.create(sampledHostsSaveDir3)

readDir <- paste0("/project2/pascualmm/QZ/PhD/projects/intervention/files", run, "/actualRuns/infTablePresence/")
prefix <- "sim"
nums <- 101:109

load(paste0("~/others/papersOfficial/FOI/analysis/utils/s_givenMOI_list"))
p_s_givenMOI <- s_givenMOI_list
MOI_max <- 20
names(p_s_givenMOI) <- as.character(1:MOI_max)

MOI_priors <- rep(1/MOI_max, MOI_max)
names(MOI_priors) <- as.character(1:MOI_max)

f <- read.table(paste0("~/others/papersOfficial/FOI/analysis/utils/NbVarGenes_MOI1_upsBC_AllSurveys_Weight.txt"), header = T)

a <- min(f$DBLa_upsBC_rep_size)
b <- max(f$DBLa_upsBC_rep_size)

s_single <- max(45, b)

preIRS <- 200
IRSDur <- 10
T_YEAR <- 360
layers <- c(preIRS - 1, preIRS + 1)*T_YEAR  + 300

p_microscopy <- 0.57

for (i in 1:length(nums)) {
  num <- nums[i]
  sampledHostsSaveDir4 <- paste0(sampledHostsSaveDir3, num, "/")
  dir.create(sampledHostsSaveDir4)
  MOIAll <- NULL
  load(paste0(readDir, seasonality, "/", openness, "/", prefix, "_", num, "_r0.RData"))
  for (ii in 1:length(layers)) {
    layer <- layers[ii]
    infStrain_pre_sub <- infStrain_pre %>% filter(presence == 1, time == layer, is_functional == 1)
    hostsSub <- sample(unique(infStrain_pre_sub$host_id), round(p_microscopy*length(unique(infStrain_pre_sub$host_id))))
    microscopyPosHosts <- hostsSub
    # save(microscopyPosHosts, file = paste0(sampledHostsSaveDir4, "rep_0", "time_", layer, ".RData"))
    infStrain_pre_sub <- infStrain_pre_sub %>% filter(host_id %in% hostsSub)
    s_temp_dat <- infStrain_pre_sub %>% group_by(host_id) %>% summarise(n = n_distinct(gene_id))
      
    s_temp <- s_temp_dat$n
    names <- s_temp_dat$host_id
    for (j in 1:length(s_temp)) {
      s <- s_temp[j]
      if (s >= a & s <= s_single*MOI_max) {
        numerator <- 0
        for (b in 1:MOI_max) {
          if (!is.na(p_s_givenMOI[[as.character(b)]][as.character(s)])) {
            numerator <- numerator + p_s_givenMOI[[as.character(b)]][as.character(s)]*MOI_priors[as.character(b)]
          }
        }
          
        if (numerator == 0) {
          p_c_givens_all <- rep(0, length(1:MOI_max))
        } else {
          p_c_givens_all <- rep(NA, length(1:MOI_max))
          for (c in 1:MOI_max) {
            if (!is.na(p_s_givenMOI[[as.character(c)]][as.character(s)])) {
              temp <- p_s_givenMOI[[as.character(c)]][as.character(s)]*MOI_priors[as.character(c)]/numerator
              names(temp) <- NULL
            } else {
              temp <- 0
            }
            p_c_givens_all[c] <- temp
          }
        }
          
        print(sum(p_c_givens_all))  
        df <- data.frame("p" = max(p_c_givens_all), "MOI" = c(1:MOI_max)[which.max(p_c_givens_all)], "numDBLaTypes" = s, 
                         "host_id" = names[j], "num" = num, "time" = layer)
        MOIAll <- rbind(MOIAll, df)
      }
    }
  }
  save(MOIAll, file = paste0(saveDir3, "MOI_", num, ".RData"))
}
