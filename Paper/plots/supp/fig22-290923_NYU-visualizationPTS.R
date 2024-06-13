rm(list = ls())
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggpubr)
})
sizeV <- 20
sizeVFactor <- 1.15
readDir0 <- "/project2/pascualmm/QZ/PhD/projects/intervention/natComRevision/files/phaseDiagram/"
PTSTypes <- "quantileDiff"
PTSQuantileThresholds <- c(0.40, 0.30)
quantilesChosen <- seq(0,1,0.001)
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}
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
date <- "290923_NYU"
saveDir0 <- paste0("/project2/pascualmm/QZ/PhD/projects/intervention/natComRevision/plots/figures/supp/fig22-", date, "/PTSQuantileDiff/")
if (!dir.exists(saveDir0)) {
  dir.create(saveDir0)
}
for (i in 1:length(PTSTypes)) {
  PTSType <- PTSTypes[i]
  saveDir1 <- paste0(saveDir0, PTSType, "/")
  if (!dir.exists(saveDir1)) {
    dir.create(saveDir1)
  }
  saveDir2 <- paste0(saveDir1, "sim/")
  if (!dir.exists(saveDir2)) {
    dir.create(saveDir2)
  }
  saveDir3 <- paste0(saveDir1, "Ghana/")
  if (!dir.exists(saveDir3)) {
    dir.create(saveDir3)
  }
  file1 <- paste0(readDir0, "simulation", "_PTSType_", PTSType, ".RData")
  load(file1)
  file2 <- paste0(readDir0, "Ghana", "_PTSType_", PTSType, "-", date, ".RData")
  load(file2)
  for (j in 1:length(states)) {
    state <- states[j]
    saveDir4 <- paste0(saveDir2, state, "/")
    if (!dir.exists(saveDir4)) {
      dir.create(saveDir4)
    }
    PTSQuantileThreshold <- PTSQuantileThresholds[j]
    PTSdfAllGhana <- PTSdfAllGhana %>% mutate(Regime = paste0("Ghana ", state))
    
    for (k in 1:length(seasonality)) {
      s <- seasonality[k]
      saveDir5 <- paste0(saveDir4, s, "/")
      if (!dir.exists(saveDir5)) {
        dir.create(saveDir5)
      }
      for (m in 1:length(openness)) {
        o <- openness[m]
        nums <- numsList[[paste0(s, " ", o)]]
        saveDir6 <- paste0(saveDir5, o, "/")
        if (!dir.exists(saveDir6)) {
          dir.create(saveDir6)
        }
        for (n in 1:length(nums)) {
          numsSingle <- nums[[n]]
          label <- names(numsList[[paste0(s, " ", o)]])[n]
          xtext <- seq(0, PTSQuantileThreshold, length.out = 3)
          dfPlotSim <- PTSdfAllSim %>% filter(
            quantiles <= PTSQuantileThreshold, 
            quantiles %in% quantilesChosen,
            Seasonality == firstup(s),
            Openness == firstup(o),
            Label == label,
            State == state) %>% 
            select(quantiles, PTSQuantilesDiff, Regime)
          
          for (p in 1:length(unique(dfPlotSim$Regime))) {
            RegimeSingle <- unique(dfPlotSim$Regime)[p]
            if (RegimeSingle == "Pre-IRS") {
              ytitleC <- "black"
            } else {
              ytitleC <- "white"
            }
            p1 <- ggplot(dfPlotSim %>% filter(Regime == RegimeSingle), 
                         aes(x = quantiles, y = PTSQuantilesDiff)) + 
              geom_point() +
              ggtitle(RegimeSingle) + 
              xlab("Quantiles") + ylab("Difference in PTS (Children - Adults)") +
              theme_bw() + theme(
                plot.title = element_text(color="black", size=sizeV, hjust = 0.5),
                axis.title.x = element_text(color="black", size=sizeV),
                axis.title.y = element_text(color=ytitleC, size=sizeV),
                axis.text.x = element_text(color="black", size=sizeV, angle=0),
                axis.text.y = element_text(color="black", size=sizeV, angle=0),
                legend.text = element_text(color="black", size=sizeV, angle=0),
                legend.title = element_text(color="black", size=sizeV, angle=0),
                strip.text = element_text(color="black", size=sizeV, angle=0)) +
              scale_x_continuous(breaks=xtext, labels=as.character(xtext))
            p1_rmL <- p1 + guides(fill = "none", col = "none") 
            ggsave(p1_rmL, filename = paste0(saveDir6, label, "_", RegimeSingle, ".pdf"), height = 5, width = 5)
          }
          
          if (s == "seasonal" & o == "closed" & state == "withME") {
            # p1_legend <- get_legend(p1)
            # ggsave(paste0(saveDir2, "lg.pdf"), p1_legend, width = 5, height = 5)
            
            dfPlotGhana <- PTSdfAllGhana %>% filter(
              quantiles <= PTSQuantileThreshold, 
              quantiles %in% quantilesChosen) %>% 
              select(quantiles, PTSQuantilesDiff, Regime)
            
            for (q in 1:length(unique(dfPlotGhana$Regime))) {
              RegimeSingle <- unique(dfPlotGhana$Regime)[q]
              if (RegimeSingle == "Ghana Pre-IRS") {
                ytitleC <- "black"
              } else {
                ytitleC <- "white"
              }
              p2 <- ggplot(dfPlotGhana %>% filter(Regime == RegimeSingle), 
                           aes(x = quantiles, y = PTSQuantilesDiff)) + 
                geom_point() +
                ggtitle(RegimeSingle) + 
                xlab("Quantiles") + ylab("Difference in PTS (Children - Adults)") +
                theme_bw() + theme(
                  plot.title = element_text(color="black", size=sizeV, hjust = 0.5),
                  axis.title.x = element_text(color="black", size=sizeV),
                  axis.title.y = element_text(color=ytitleC, size=sizeV),
                  axis.text.x = element_text(color="black", size=sizeV, angle=0),
                  axis.text.y = element_text(color="black", size=sizeV, angle=0),
                  legend.text = element_text(color="black", size=sizeV, angle=0),
                  legend.title = element_text(color="black", size=sizeV, angle=0),
                  strip.text = element_text(color="black", size=sizeV, angle=0)) +
                scale_x_continuous(breaks=xtext, labels=as.character(xtext))
              p2_rmL <- p2 + guides(fill = "none", col = "none") 
              ggsave(p2_rmL, filename = paste0(saveDir3, RegimeSingle, ".pdf"), height = 5, width = 5)
            }
            # p2_legend <- get_legend(p2)
            # ggsave(paste0(saveDir3, "lg.pdf"), p2_legend, width = 5, height = 5)
          }
        }
      }
    }
  }
}
