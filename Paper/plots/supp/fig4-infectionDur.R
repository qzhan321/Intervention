rm(list = ls())
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(ggpubr)
})
run <- 4
readDir <- paste0("/project2/pascualmm/QZ/PhD/projects/intervention/files", run, "/actualRuns/infectionDur/")
openness <- c("closed", "semi-open", "regionally-open")
seasonality <- c("seasonal", "non-seasonal")
preIRSList <- list(
  "seasonal closed" = 200,
  "seasonal semi-open" = 150, 
  "seasonal regionally-open" = 150,
  "non-seasonal closed" = 200,
  "non-seasonal semi-open" = 150,
  "non-seasonal regionally-open" = 150
)
IRSDur <- 10
IRSType <- paste0(IRSDur, "yIRS")
postIRS <- 10 - IRSDur

numsPlotList <- list(
  "seasonal closed" = list(104),
  "seasonal semi-open" = list(204, 304),
  "seasonal regionally-open" = list(204, 222, 240, 258),
  "non-seasonal closed" = list(104),
  "non-seasonal semi-open" = list(204, 304),
  "non-seasonal regionally-open" = list(204, 222, 240, 258)
)
labelsList <- list(
  "seasonal closed" = list("Seasonal \nClosed \n \n", "Seasonal \nClosed \n \n"),
  "seasonal semi-open" = list("Seasonal \nSemi-open \nBaseline Migration \n", 
                              "Seasonal \nSemi-open \nHigh Migration \n"),
  "seasonal regionally-open" = list("Seasonal \nRegionally-open \nBaseline Migration \nMedium Pool",
                                    "Seasonal \nRegionally-open \nBaseline Migration \nLarge Pool",
                                    "Seasonal \nRegionally-open \nHigh Migration \nMedium Pool",
                                    "Seasonal \nRegionally-open \nHigh Migration \nLarge Pool"),
  "non-seasonal closed" = list("Non-seasonal \nClosed \n \n"),
  "non-seasonal semi-open" = list("Non-seasonal \nSemi-open \nBaseline Migration \n",   
                                  "Non-seasonal \nSemi-open \nHigh Migration \n"),
  "non-seasonal regionally-open" = list("Non-seasonal \nRegionally-open \nBaseline Migration \nMedium Pool",
                                        "Non-seasonal \nRegionally-open \nBaseline Migration \nLarge Pool",
                                        "Non-seasonal \nRegionally-open \nHigh Migration \nMedium Pool",
                                        "Non-seasonal \nRegionally-open \nHigh Migration \nLarge Pool")
)
immsList <- list(
  "seasonal closed" = list("Three-year"),
  "seasonal semi-open" = list("Three-year", "Three-year"),
  "seasonal regionally-open" = list("Three-year", "Three-year", "Three-year", "Three-year"),
  "non-seasonal closed" = list("Three-year"),
  "non-seasonal semi-open" = list("Three-year", "Three-year"),
  "non-seasonal regionally-open" = list("Three-year", "Three-year", "Three-year", "Three-year")
)

T_YEAR <- 360
N <- 10000
sizeV <- 30

saveDir0 <- "/project2/pascualmm/QZ/PhD/projects/intervention/natComRevision/plots/figures/supp/fig4/infectionDur/"
if (!dir.exists(saveDir0)) {
  dir.create(saveDir0)
}
type <- "true"
prefix <- "sim"
rep <- 0
for (i in 1:length(seasonality)) {
  s <- seasonality[i]
  saveDir1 <- paste0(saveDir0, s, "/")
  if (!dir.exists(saveDir1)) {
    dir.create(saveDir1)
  }
  for (j in 1:length(openness)) {
    o <- openness[j]
    saveDir2 <- paste0(saveDir1, o, "/")
    if (!dir.exists(saveDir2)) {
      dir.create(saveDir2)
    }
    preIRS <- preIRSList[[paste0(s, " ", o)]]
    numsPlot <- numsPlotList[[paste0(s, " ", o)]]
    labels <- labelsList[[paste0(s, " ", o)]]
    imms <- immsList[[paste0(s, " ", o)]]
    for (k in 1:length(numsPlot)) {
      nums <- numsPlot[[k]]
      label <- labels[[k]]
      imm <- imms[[k]]
      saveDir3 <- paste0(saveDir2, imm, "/")
      if (!dir.exists(saveDir3)) {
        dir.create(saveDir3)
      }
      
      for (m in 1:length(nums)) {
        num <- nums[m]
        file <- paste0(readDir, s, "/", o, "/", type, "/", prefix, "_", num, "_rep_", rep, "_", IRSType, ".RData")
        load(file)
        durTable$IRS <- factor(durTable$IRS, levels = c("Pre-IRS", "Early Stage of IRS", "Late Stage of IRS"))
        if (s == "seasonal" & o == "closed") {
          ytitleC <- "black"
        } else if (s == "seasonal" & o == "regionally-open" & num %in% 237:245) {
          ytitleC <- "black"
        } else if (s == "non-seasonal" & o == "regionally-open" & num %in% 201:210) {
          ytitleC <- "black"
        } else {
          ytitleC <- "white"
        }
        if (s == "seasonal" & o == "semi-open" & num %in% 301:311) {
          xtitleC <- "black"
        } else if (s == "non-seasonal" & o == "closed" & num %in% 101:111) {
          xtitleC <- "black"
        } else if (s == "non-seasonal" & o == "regionally-open" & num %in% 237:246) {
          xtitleC <- "black"
        } else {
          xtitleC <- "white"
        }
        p <- ggplot(durTable, aes(duration, col = IRS, fill = IRS, alpha = 0.5)) + 
          geom_density()+
          ggtitle(label) +
          xlab("Infection Duration (Days)") + ylab("Density") +
          theme_bw() + theme(
            plot.title = element_text(color="black", size=sizeV/1.25, hjust = 0.5),
            axis.title.x = element_text(color=xtitleC, size=sizeV),
            axis.title.y = element_text(color=ytitleC, size=sizeV),
            axis.text.x = element_text(color="black", size=sizeV, angle=0),
            axis.text.y = element_text(color="black", size=sizeV, angle=0),
            legend.text = element_text(color="black", size=sizeV, angle=0),
            legend.title = element_text(color="black", size=sizeV, angle=0),
            strip.text = element_text(color="black", size=sizeV, angle=0),
            plot.margin = unit(c(5,10,5,5), "points")) + 
          scale_color_manual(values = c(alpha("#0000a7", 0.5), alpha("#008176", 0.5), alpha("#8d1f17", 0.5)), name = NULL) +
          scale_fill_manual(values = c(alpha("#0000a7", 0.5), alpha("#008176", 0.5), alpha("#8d1f17", 0.5)), name = NULL) +
          guides(fill = "none", col = "none", alpha = "none") 
        ggsave(paste0(saveDir3, prefix, "_", num, "_rep", rep, "_durInfo-", IRSType, ".pdf"), p, width = 9.25, height = 9.25)
      }
    }
  }
}

p <- ggplot(durTable, aes(duration, col = IRS, fill = IRS, alpha = 0.5)) + 
  geom_density()+
  ggtitle(label) +
  xlab("Infection Duration (Days)") + ylab("Density") +
  theme_bw() + theme(
    plot.title = element_text(color="black", size=sizeV/1.25, hjust = 0.5),
    axis.title.x = element_text(color=xtitleC, size=sizeV),
    axis.title.y = element_text(color=ytitleC, size=sizeV),
    axis.text.x = element_text(color="black", size=sizeV, angle=0),
    axis.text.y = element_text(color="black", size=sizeV, angle=0),
    legend.text = element_text(color="black", size=sizeV, angle=0),
    legend.title = element_text(color="black", size=sizeV, angle=0),
    strip.text = element_text(color="black", size=sizeV, angle=0),
    plot.margin = unit(c(5,10,5,5), "points")) + 
  scale_color_manual(values = c(alpha("#0000a7", 0.5), alpha("#008176", 0.5), alpha("#8d1f17", 0.5)), name = NULL) +
  scale_fill_manual(values = c(alpha("#0000a7", 0.5), alpha("#008176", 0.5), alpha("#8d1f17", 0.5)), name = NULL) +
  guides(alpha = "none") 
p_legend <- get_legend(p)
ggsave(paste0(saveDir0, "lg.pdf"), p_legend, width = 9.25, height = 9.25)
