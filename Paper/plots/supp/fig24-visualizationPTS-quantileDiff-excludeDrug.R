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
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}
PTSQuantileThreshold <- 0.40
quantilesChosen <- seq(0,1,0.001)
xtext <- seq(0, PTSQuantileThreshold, length.out = 3)
saveDir0 <- "/project2/pascualmm/QZ/PhD/projects/intervention/natComRevision/plots/figures/supp/fig24/PTSQuantileDiffExcludeDrug/"
if (!dir.exists(saveDir0)) {
  dir.create(saveDir0)
}
for (i in 1:length(PTSTypes)) {
  PTSType <- PTSTypes[i]
  saveDir1 <- paste0(saveDir0, PTSType, "/")
  if (!dir.exists(saveDir1)) {
    dir.create(saveDir1)
  }
  saveDir2 <- paste0(saveDir1, "Ghana/")
  if (!dir.exists(saveDir2)) {
    dir.create(saveDir2)
  }
  file1 <- paste0(readDir0, "Ghana", "_PTSType_", PTSType, "-excludeDrug.RData")
  load(file1)
  PTSdfAllGhana <- PTSdfAllGhana %>% mutate(Regime = paste0("Ghana ", state))
  dfPlotGhana <- PTSdfAllGhana %>% filter(
    quantiles <= PTSQuantileThreshold, 
    quantiles %in% quantilesChosen) %>% 
    select(quantiles, PTSQuantilesDiff, Regime)
  for (j in 1:length(unique(dfPlotGhana$Regime))) {
    RegimeSingle <- unique(dfPlotGhana$Regime)[j]
    if (RegimeSingle == "Ghana Pre-IRS") {
      ytitleC <- "black"
    } else {
      ytitleC <- "white"
    }
    p1 <- ggplot(dfPlotGhana %>% filter(Regime == RegimeSingle), aes(x = quantiles, y = PTSQuantilesDiff)) + 
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
    ggsave(p1_rmL, filename = paste0(saveDir2, RegimeSingle, "-excludeDrug.pdf"), height = 5, width = 5)
  }
}
