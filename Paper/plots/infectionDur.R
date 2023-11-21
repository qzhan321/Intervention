rm(list = ls())
suppressPackageStartupMessages({
  library(ggplot2)
  library(gridExtra)
  library(scales)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(ggpubr)
  library(grid)
  library(cowplot)
})
run <- 4
source(paste0("/home/qizhan/others/PhD/projects/intervention/writings/simulation", run, "/plots/colorFunc.R"))
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
postIRS <- 10-IRSDur

numsPlotList <- list(
  "seasonal closed" = list(101:110,115:123),
  "seasonal semi-open" = list(201:211, 216:225,  
                              301:311, 316:325),
  "seasonal regionally-open" = list(201:209,219:229,
                                    237:245,255:265,
                                    301:308,319:328,
                                    337:344,355:364),
  "non-seasonal closed" = list(101:111, 115:124),
  "non-seasonal semi-open" = list(201:211,
                                  216:225,301:311,316:326),
  "non-seasonal regionally-open" = list(201:210,219:228,
                                        237:246,255:264,
                                        301:309,319:327,
                                        337:345,355:363)
)
labelsList <- list(
  "seasonal closed" = list("Seasonal \nClosed \n \n", "Seasonal \nClosed \n \n"),
  "seasonal semi-open" = list("Seasonal \nSemi-open \nBaseline Migration \n", "Seasonal \nSemi-open \nBaseline Migration \n", 
                              "Seasonal \nSemi-open \nHigh Migration \n", "Seasonal \nSemi-open \nHigh Migration \n"),
  "seasonal regionally-open" = list("Seasonal \nRegionally-open \nBaseline Migration \nMedium Pool",
                                    "Seasonal \nRegionally-open \nBaseline Migration \nLarge Pool",
                                    "Seasonal \nRegionally-open \nHigh Migration \nMedium Pool",
                                    "Seasonal \nRegionally-open \nHigh Migration \nLarge Pool",
                                    "Seasonal \nRegionally-open \nBaseline Migration \nMedium Pool",
                                    "Seasonal \nRegionally-open \nBaseline Migration \nLarge Pool",
                                    "Seasonal \nRegionally-open \nHigh Migration \nMedium Pool",
                                    "Seasonal \nRegionally-open \nHigh Migration \nLarge Pool"),
  "non-seasonal closed" = list("Non-seasonal \nClosed \n \n", "Non-seasonal \nClosed \n \n"),
  "non-seasonal semi-open" = list("Non-seasonal \nSemi-open \nBaseline Migration \n", "Non-seasonal \nSemi-open \nBaseline Migration \n",  
                                  "Non-seasonal \nSemi-open \nHigh Migration \n", "Non-seasonal \nSemi-open \nHigh Migration \n"),
  "non-seasonal regionally-open" = list("Non-seasonal \nRegionally-open \nBaseline Migration \nMedium Pool",
                                        "Non-seasonal \nRegionally-open \nBaseline Migration \nLarge Pool",
                                        "Non-seasonal \nRegionally-open \nHigh Migration \nMedium Pool",
                                        "Non-seasonal \nRegionally-open \nHigh Migration \nLarge Pool",
                                        "Non-seasonal \nRegionally-open \nBaseline Migration \nMedium Pool",
                                        "Non-seasonal \nRegionally-open \nBaseline Migration \nLarge Pool",
                                        "Non-seasonal \nRegionally-open \nHigh Migration \nMedium Pool",
                                        "Non-seasonal \nRegionally-open \nHigh Migration \nLarge Pool")
  
)
immsList <- list(
  "seasonal closed" = list("Three-year", "Four-year"),
  "seasonal semi-open" = list("Three-year", "Four-year",
                              "Three-year", "Four-year"),
  "seasonal regionally-open" = list("Three-year", "Three-year", "Three-year", "Three-year",
                                    "Four-year", "Four-year", "Four-year", "Four-year"),
  "non-seasonal closed" = list("Three-year", "Four-year"),
  "non-seasonal semi-open" = list("Three-year", "Four-year",
                                  "Three-year", "Four-year"),
  "non-seasonal regionally-open" = list("Three-year", "Three-year", "Three-year", "Three-year",
                                        "Four-year", "Four-year", "Four-year", "Four-year")
)

T_YEAR <- 360

N <- 10000
pAll1 <- list()
pAll2 <- list()
sizeV <- 31

saveDir0 <- paste0("/project2/pascualmm/QZ/PhD/projects/intervention/writings/simulation", run, "/plots/")
dir.create(saveDir0)
saveDir1 <- paste0(saveDir0, "infectionDur/")
dir.create(saveDir1)
type <- "true"
prefix <- "sim"
for (i in 1:length(seasonality)) {
  s <- seasonality[i]
  saveDir2 <- paste0(saveDir1, s, "/")
  dir.create(saveDir2)
  
  for (j in 1:length(openness)) {
    o <- openness[j]
    saveDir3 <- paste0(saveDir2, o, "/")
    dir.create(saveDir3)
  
    preIRS <- preIRSList[[paste0(s, " ", o)]]
    numsPlot <- numsPlotList[[paste0(s, " ", o)]]
    labels <- labelsList[[paste0(s, " ", o)]]
    imms <- immsList[[paste0(s, " ", o)]]
    for (k in 1:length(numsPlot)) {
      nums <- numsPlot[[k]]
      label <- labels[[k]]
      imm <- imms[[k]]
      
      saveDir4 <- paste0(saveDir3, imm, "/")
      dir.create(saveDir4)
      
      for (m in 1:length(nums)) {
        num <- nums[m]
        rep <- 0
        file <- paste0(readDir, s, "/", o, "/", type, "/", prefix, "_", num, "_rep_", rep, "_", IRSType, ".RData")
        load(file)
        durTable$IRS <- factor(durTable$IRS, levels = c("Pre-IRS", paste0("Early Stage of IRS"), paste0("Late Stage of IRS")))
        p <- ggplot(durTable, aes(duration, col = IRS, fill = IRS, alpha = 0.5)) + 
          geom_density()+
          xlab("Infection Duration (Days)") + ylab("Density") +
          theme_bw() + theme(
            plot.title = element_text(color="black", size=sizeV, hjust = 0.5),
            axis.title.x = element_text(color="black", size=sizeV),
            axis.title.y = element_text(color="black", size=sizeV),
            axis.text.x = element_text(color="black", 
                                       size=sizeV, angle=0),
            axis.text.y = element_text(color="black", 
                                       size=sizeV, angle=0),
            legend.text = element_text(color="black", 
                                       size=sizeV, angle=0),
            legend.title = element_text(color="black", 
                                        size=sizeV, angle=0),
            strip.text = element_text(color="black", 
                                      size=sizeV, angle=0),
            plot.margin = unit(c(5,20,5,5), "points")) + 
          scale_color_manual(values = c(alpha("#0000a7", 0.5), alpha("#008176", 0.5), alpha("#8d1f17", 0.5)), name = NULL) +
          scale_fill_manual(values = c(alpha("#0000a7", 0.5), alpha("#008176", 0.5), alpha("#8d1f17", 0.5)), name = NULL) +
          guides(fill = "none", col = "none", alpha = "none") 
        ggsave(paste0(saveDir4, prefix, "_", num, "_rep", rep, "_durInfo-", IRSType, ".pdf"), p, width = 6.5, height = 6.5)
      }
    }
  }
}

# p <- ggplot(durTable, aes(duration, col = IRS, fill = IRS, alpha = 0.5)) + 
#   geom_density()+
#   xlab("Duration of Infection (Days)") + ylab("Density") +
#   theme_bw() + theme(
#     plot.title = element_text(color="black", size=sizeV, hjust = 0.5),
#     axis.title.x = element_text(color="white", size=sizeV),
#     axis.title.y = element_text(color="white", size=sizeV),
#     axis.text.x = element_text(color="black", 
#                                size=sizeV, angle=0),
#     axis.text.y = element_text(color="black", 
#                                size=sizeV, angle=0),
#     legend.text = element_text(color="black", 
#                                size=sizeV, angle=0),
#     legend.title = element_text(color="black", 
#                                 size=sizeV, angle=0),
#     strip.text = element_text(color="black", 
#                               size=sizeV, angle=0)) + 
#   scale_color_manual(values = c(alpha("#0000a7", 0.5), alpha("#008176", 0.5), alpha("#8d1f17", 0.5)), name = NULL) +
#   scale_fill_manual(values = c(alpha("#0000a7", 0.5), alpha("#008176", 0.5), alpha("#8d1f17", 0.5)), name = NULL) +
#   guides(alpha = "none") 
# p_legend <- get_legend(p)
# ggsave(paste0(saveDir1, "lg.pdf"), p_legend, width = 6.5, height = 6.5)
