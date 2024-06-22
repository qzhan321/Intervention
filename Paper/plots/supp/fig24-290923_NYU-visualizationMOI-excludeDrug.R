rm(list = ls())
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})
statsSummaryf <- function(x) {
  r <- quantile(x, probs = c(0.00, 0.05, 0.5, 0.95, 1.00))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}
date <- "290923_NYU"
sizeV <- 12.5
sizeVFactor <- 1.15
readDir0 <- "/project2/pascualmm/QZ/PhD/projects/intervention/natComRevision/files/phaseDiagram/"
PTSTypes <- "quantileDiff"
states <- c("withME")
saveDir0 <- paste0("/project2/pascualmm/QZ/PhD/projects/intervention/natComRevision/plots/figures/supp/fig24-", date, "/MOIExcludeDrug/")
if (!dir.exists(saveDir0)) {
  dir.create(saveDir0)
}
for (i in 1:length(PTSTypes)) {
  PTSType <- PTSTypes[i]
  file1 <- paste0(readDir0, "simulation", "_PTSType_", PTSType, ".RData")
  load(file1)
  file2 <- paste0(readDir0, "Ghana", "_PTSType_", PTSType, "-excludeDrug-", date, ".RData")
  load(file2)
  for (j in 1:length(states)) {
    state <- states[j]
    saveDir1 <- paste0(saveDir0, state, "/")
    if (!dir.exists(saveDir1)) {
      dir.create(saveDir1)
    }
    dfPlotSim <- MOIdfAllSim %>% filter(State == state, Seasonality == "Seasonal") %>% 
      ungroup() %>%
      select(Regime, PropMonoclonal) 
    dfPlotGhana <- MOIdfAllGhana %>% mutate(Regime = case_when(
      (survey == "survey_1") ~ "Ghana Pre-IRS",
      (survey == "survey_4") ~ "Ghana IRS",
      (survey == "survey_5") ~ "Ghana Post-IRS"
    )) %>% select(-survey)
    dfPlot <- bind_rows(dfPlotSim, dfPlotGhana)
    dfPlot$Regime <- factor(dfPlot$Regime, levels = c("Pre-IRS", "Fast-rebound", "Approaching Transition",
                                                      "Transition", "Slow-rebound",
                                                      "Ghana Pre-IRS", "Ghana IRS", "Ghana Post-IRS"))
    p1 <- ggplot(dfPlot, aes(x = Regime, y = PropMonoclonal)) + 
      stat_summary(fun.data = statsSummaryf, geom = "boxplot", position = position_dodge(1)) +
      xlab("Regime") + ylab("Proportion of MOI = 1 or 2") +
      theme_bw() + theme(
        plot.title = element_text(color="black", size=sizeV, hjust = 0.5),
        axis.title.x = element_text(color="black", size=sizeV),
        axis.title.y = element_text(color="black", size=sizeV),
        axis.text.x = element_text(color="black", size=sizeV, angle=-45),
        axis.text.y = element_text(color="black", size=sizeV, angle=0),
        legend.text = element_text(color="black", size=sizeV, angle=0),
        legend.title = element_text(color="black", size=sizeV, angle=0),
        strip.text = element_text(color="black", size=sizeV, angle=0)) 
    p1
    ggsave(p1, filename = paste0(saveDir1, "MOI.pdf"), height = 6, width = 6)
  }
}
