rm(list = ls())
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})
# simulation
readDir0 <- "/project2/pascualmm/QZ/PhD/projects/intervention/files4/actualRuns/"
folders <- c("MOI", "PTSAgeDiff1")
seasonality <- c("seasonal", "non-seasonal")
openness <- c("closed", "semi-open", "regionally-open")
states <- c("withME", "true")
numsList <- list(
  "seasonal closed" = list("Seasonal \nClosed \n \n" = 101:108),
  "seasonal semi-open" = list("Seasonal \nSemi-open \nBaseline Migration \n" = 201:208, 
                              "Seasonal \nSemi-open \nHigh Migration \n" = 301:307),
  "seasonal regionally-open" = list("Seasonal \nRegionally-open \nBaseline Migration \nMedium Pool" = 201:207, 
                                    "Seasonal \nRegionally-open \nBaseline Migration \nLarge Pool" = 219:227,
                                    "Seasonal \nRegionally-open \nHigh Migration \nMedium Pool" = 237:243, 
                                    "Seasonal \nRegionally-open \nHigh Migration \nLarge Pool" = 255:263),
  "non-seasonal closed" = list("Non-seasonal \nClosed \n \n" = 101:108),
  "non-seasonal semi-open" = list("Non-seasonal \nSemi-open \nBaseline Migration \n" = 201:208, 
                                  "Non-seasonal \nSemi-open \nHigh Migration \n" = 301:308),
  "non-seasonal regionally-open" = list("Non-seasonal \nRegionally-open \nBaseline Migration \nMedium Pool" = 201:207, 
                                        "Non-seasonal \nRegionally-open \nBaseline Migration \nLarge Pool" = 219:227,
                                        "Non-seasonal \nRegionally-open \nHigh Migration \nMedium Pool" = 237:243, 
                                        "Non-seasonal \nRegionally-open \nHigh Migration \nLarge Pool" = 255:263)
)
regimesList <- list(
  "seasonal closed" = list(c("Fast-rebound", "Fast-rebound", "Fast-rebound", "Fast-rebound", 
                             "Approaching Transition", 
                             "Transition", "Transition", 
                             "Slow-rebound")),
  "seasonal semi-open" = list(c("Fast-rebound", "Fast-rebound", "Fast-rebound", "Fast-rebound", 
                                "Approaching Transition",
                                "Transition", "Transition", "Transition"), 
                              c("Fast-rebound", "Fast-rebound", "Fast-rebound", "Fast-rebound", 
                                "Approaching Transition",
                                "Transition", "Transition")),
  "seasonal regionally-open" = list(c("Fast-rebound", "Fast-rebound", "Fast-rebound", 
                                      "Approaching Transition", 
                                      "Transition", "Transition", 
                                      "Slow-rebound"), 
                                    c("Fast-rebound", "Fast-rebound", "Fast-rebound",
                                      "Approaching Transition", 
                                      "Transition", "Transition", "Transition", 
                                      "Slow-rebound", "Slow-rebound"),
                                    c("Fast-rebound", "Fast-rebound", "Fast-rebound",  
                                      "Approaching Transition", 
                                      "Transition", "Transition", 
                                      "Slow-rebound"), 
                                    c("Fast-rebound", "Fast-rebound", "Fast-rebound", "Fast-rebound",
                                      "Approaching Transition", 
                                      "Transition", "Transition", 
                                      "Slow-rebound", "Slow-rebound")),
  "non-seasonal closed" = list(c("Fast-rebound", "Fast-rebound", "Fast-rebound", 
                                 "Approaching Transition",
                                 "Transition", "Transition", "Transition", 
                                 "Slow-rebound")),
  "non-seasonal semi-open" = list(c("Fast-rebound", "Fast-rebound", "Fast-rebound",
                                    "Approaching Transition", 
                                    "Transition", "Transition", "Transition", 
                                    "Slow-rebound"), 
                                  c("Fast-rebound", "Fast-rebound", "Fast-rebound",
                                    "Approaching Transition", 
                                    "Transition", "Transition", "Transition", 
                                    "Slow-rebound")),
  "non-seasonal regionally-open" = list(c("Fast-rebound", "Fast-rebound", "Fast-rebound",
                                          "Approaching Transition", 
                                          "Transition", "Transition", 
                                          "Slow-rebound"), 
                                        c("Fast-rebound", "Fast-rebound", "Fast-rebound",
                                          "Approaching Transition", 
                                          "Transition", "Transition", 
                                          "Slow-rebound", "Slow-rebound", "Slow-rebound"),
                                        c("Fast-rebound", "Fast-rebound", "Fast-rebound",
                                          "Approaching Transition", 
                                          "Transition", "Transition", 
                                          "Slow-rebound"), 
                                        c("Fast-rebound", "Fast-rebound",
                                          "Approaching Transition", 
                                          "Transition", "Transition", "Transition",
                                          "Slow-rebound", "Slow-rebound", "Slow-rebound"))
)
T_YEAR <- 360
type <- "dir"
PTSType <- "quantileDiff"
saveDir <- "/project2/pascualmm/QZ/PhD/projects/intervention/natComRevision/files/phaseDiagram/"
MOIdfAllSim <- NULL
PTSdfAllSim <- NULL
for (i in 1:length(folders)) {
  folder <- folders[i]
  readDir1 <- paste0(readDir0, folder, "/")
  for (j in 1:length(seasonality)) {
    s <- seasonality[j]
    readDir2 <- paste0(readDir1, s, "/")
    for (k in 1:length(openness)) {
      o <- openness[k]
      if (o == "closed") {
        preIRS <- 200
      } else {
        preIRS <- 150
      }
      if (s == "seasonal") {
        times <- c((preIRS - 1)*T_YEAR + 300, (preIRS + 1)*T_YEAR + 300)
      } else {
        times <- c(preIRS*T_YEAR, (preIRS + 2)*T_YEAR)
      }
      readDir3 <- paste0(readDir2, o, "/")
      nums <- numsList[[paste0(s, " ", o)]]
      regimes <- regimesList[[paste0(s, " ", o)]]
      for (m in 1:length(nums)) {
        numsSingle <- nums[[m]]
        label <- names(nums)[m]
        regimesSingle <- regimes[[m]]
        for (n in 1:length(states)) {
          state <- states[n]
          if (folder == "MOI") {
            readDir4 <- paste0(readDir3, state, "/")
          } else {
            readDir4 <- paste0(readDir3, type, "/", state, "/")
          }
          for (p in 1:length(numsSingle)) {
            numSingle <- numsSingle[p]
            regimeSingle <- regimesSingle[p]
            if (folder == "MOI") {
              file <- paste0(readDir4, "MOI_", numSingle, ".RData")
              load(file)
              MOIAllSub <- MOIAll %>% filter(time %in% times) %>% 
                mutate(IRSPhase = ifelse(time <= preIRS * T_YEAR, "Pre-IRS", "IRS")) %>%
                mutate(Indicator = "MOI", Regime = ifelse(IRSPhase == "Pre-IRS", "Pre-IRS", regimeSingle),
                       State = state)
              if (numSingle != min(numsSingle)) {
                MOIAllSub <- MOIAllSub %>% filter(IRSPhase == "IRS")
              }
              MOIdf <- MOIAllSub %>% group_by(IRSPhase, Indicator, Regime, State) %>% 
                summarise(PropMonoclonal = sum(MOI %in% c(1,2))/n())
              MOIdfAllSim <- bind_rows(MOIdfAllSim, MOIdf %>% mutate(
                Seasonality = ifelse(s == "seasonal", "Seasonal", "Non-seasonal"),
                Openness = case_when((o == "closed")~"Closed",
                                     (o == "semi-open")~"Semi-open",
                                     (o == "regionally-open")~"Regionally-open"),
                Label = label,
                IRSLabel = ifelse(IRSPhase == "Pre-IRS", "Pre-IRS", paste0("I-", numSingle - min(numsSingle) + 1))
              ))
            } else {
              file <- paste0(readDir4, "PTSAcrossAgeGroups_", numSingle, "_r0.RData")
              load(file)
              PTSAllSub <- PTSAll %>% filter(time %in% times) %>% 
                mutate(IRSPhase = ifelse(time <= preIRS * T_YEAR, "Pre-IRS", "IRS")) %>%
                mutate(Indicator = "PTS", Regime = ifelse(IRSPhase == "Pre-IRS", "Pre-IRS", regimeSingle), State = state)
              PTSdfTemp <- PTSAllSub %>% group_by(IRSPhase, Indicator, Regime, State, ageGroup) %>%
                summarise(quantiles = seq(0,1,0.001), PTSQuantiles = quantile(PTS, probs = seq(0,1,0.001)))
              PTSdfTemp1 <- PTSdfTemp %>% filter(ageGroup == "1-10")
              PTSdfTemp2 <- PTSdfTemp %>% filter(ageGroup == ">=20")
              PTSdfTemp12 <- PTSdfTemp1 %>% ungroup() %>% select(-ageGroup) %>% 
                left_join(PTSdfTemp2 %>% ungroup() %>% select(-ageGroup), by = c("IRSPhase", "Indicator", "Regime", "State", "quantiles")) %>%
                mutate(PTSQuantilesDiff = PTSQuantiles.x - PTSQuantiles.y) %>% select(-c(PTSQuantiles.x, PTSQuantiles.y))
              PTSdf <- PTSdfTemp12 
              if (numSingle != min(numsSingle)) {
                PTSdf <- PTSdf %>% filter(IRSPhase == "IRS")
              }
              PTSdfAllSim <- bind_rows(PTSdfAllSim, PTSdf %>% mutate(
                Seasonality = ifelse(s == "seasonal", "Seasonal", "Non-seasonal"),
                Openness = case_when((o == "closed")~"Closed",
                                     (o == "semi-open")~"Semi-open",
                                     (o == "regionally-open")~"Regionally-open"),
                Label = label,
                IRSLabel = ifelse(IRSPhase == "Pre-IRS", "Pre-IRS", paste0("I-", numSingle - min(numsSingle) + 1))
              ))
            }
          }
        }
      }
    }
  }
}
save(MOIdfAllSim, PTSdfAllSim, file = paste0(saveDir, "simulation", "_PTSType_", PTSType, ".RData"))

