rm(list = ls())
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(cowplot)
})
run <- 4
source(paste0("/home/qizhan/others/PhD/projects/intervention/writings/simulation", run, "/plots/colorFunc.R"))
readDir <- paste0("/project2/pascualmm/QZ/PhD/projects/intervention/files", run, "/actualRuns/epi/")
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
IRSDur <- 2
IRSType <- paste0(IRSDur, "yIRS")
postIRS <- 10 - IRSDur

numsPlotList <- list(
  "seasonal closed" = list(1:14),
  "seasonal semi-open" = list(1:14,101:115),
  "seasonal regionally-open" = list(1:12,19:34,37:48,55:70),
  "non-seasonal closed" = list(1:13),
  "non-seasonal semi-open" = list(1:15,101:115),
  "non-seasonal regionally-open" = list(1:13,19:34,37:49,55:70)
)

labelsList <- list(
  "seasonal closed" = list("Seasonal \nClosed \n \n"),
  "seasonal semi-open" = list("Seasonal \nSemi-open \nBaseline Migration \n", "Seasonal \nSemi-open \nHigh Migration \n"),
  "seasonal regionally-open" = list("Seasonal \nRegionally-open \nBaseline Migration \nMedium Pool",
                                    "Seasonal \nRegionally-open \nBaseline Migration \nLarge Pool",
                                    "Seasonal \nRegionally-open \nHigh Migration \nMedium Pool",
                                    "Seasonal \nRegionally-open \nHigh Migration \nLarge Pool"),
  "non-seasonal closed" = list("Non-seasonal \nClosed \n \n"),
  "non-seasonal semi-open" = list("Non-seasonal \nSemi-open \nBaseline Migration \n", "Non-seasonal \nSemi-open \nHigh Migration \n"),
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
sizeV <- 25

saveDir0 <- "/project2/pascualmm/QZ/PhD/projects/intervention/natComRevision/plots/figures/supp/fig25/"
if (!dir.exists(saveDir0)) {
  dir.create(saveDir0)
}
saveDir1 <- paste0(saveDir0, "prevalenceDynamicsTransientIRS/")
if (!dir.exists(saveDir1)) {
  dir.create(saveDir1)
}
pointSize <- 1.25
for (i in 1:length(seasonality)) {
  s <- seasonality[i]
  saveDir2 <- paste0(saveDir1, s, "/")
  if (!dir.exists(saveDir2)) {
    dir.create(saveDir2)
  }
  
  for (j in 1:length(openness)) {
    o <- openness[j]
    saveDir3 <- paste0(saveDir2, o, "/")
    if (!dir.exists(saveDir3)) {
      dir.create(saveDir3)
    }
    
    if (s == "seasonal" & o == "closed") next
    
    file <- paste0(readDir, s, "/", o, "/", "summaryInfoTable-", IRSType, ".RData")
    load(file)
    
    preIRS <- preIRSList[[paste0(s, " ", o)]]
    numsPlot <- numsPlotList[[paste0(s, " ", o)]]
    labels <- labelsList[[paste0(s, " ", o)]]
    imms <- immsList[[paste0(s, " ", o)]]
    for (k in 1:length(numsPlot)) {
      nums <- numsPlot[[k]]
      if (s == "seasonal" & o == "semi-open" & identical(nums, 1:14)) next
      if (s == "seasonal" & o == "regionally-open" & identical(nums, 1:12)) next
      label <- labels[[k]]
      imm <- imms[[k]]
      
      saveDir4 <- paste0(saveDir3, imm, "/")
      if (!dir.exists(saveDir4)) {
        dir.create(saveDir4)
      }
      
      summaryTable1 <- summaryTableCombined %>% filter(num %in% nums)
      summaryTable2 <- summaryTable1 %>% filter(time > (preIRS-3)*T_YEAR, time <= (preIRS + IRSDur + postIRS)*T_YEAR)
      summaryTable3 <- summaryTable2 %>% group_by(time, IRS) %>%
        summarise(mP = mean(Prevalence), vP = sd(Prevalence),
                  maxP = max(Prevalence), minP = min(Prevalence),
                  mD = mean(n_circulating_genes), vD = sd(n_circulating_genes),
                  maxD = max(n_circulating_genes), minD = min(n_circulating_genes),
                  mMOI = mean(MOI, na.rm = T), vMOI = sd(MOI, na.rm = T))
      summaryTable4 <- summaryTable3 %>% mutate(vP = ifelse(is.na(vP), 0, vP),
                                                vD = ifelse(is.na(vD), 0, vD),
                                                vMOI = ifelse(is.na(vMOI), 0, vMOI))
      summaryTable5 <- summaryTable4 %>% mutate(timePlot = time/T_YEAR - preIRS)
      summaryTable5 <- summaryTable5 %>% mutate(IRS = ifelse(IRS == "preIRS", "Pre-IRS", IRS))
      summaryTable5$IRS <- factor(summaryTable5$IRS, levels = c("Pre-IRS", paste0("I-", 1:16)))
      if (s == "seasonal" & o == "semi-open" & identical(nums, 101:115)) {
        ytitleC <- "black"
        ytextC <- "black"
      } else if (s == "non-seasonal" & o == "closed") {
        ytitleC <- "black"
        ytextC <- "black"
      } else if (s == "non-seasonal" & o == "regionally-open" & identical(nums, 19:34)) {
        ytitleC <- "black"
        ytextC <- "black"
      } else {
        ytitleC <- "white"
        ytextC <- "white"
      }
      p1 <- ggplot(summaryTable5, aes(timePlot, mP, col = IRS, fill = IRS))+
        geom_line(size = pointSize)+
        geom_ribbon(aes(y = mP, ymin = minP, ymax = maxP), alpha = .4, linetype = 0) +
        xlab("Years Since IRS Starts \n(Transient IRS)") + ylab("Prevalence") +
        ggtitle(label) +  
        theme_bw() + theme(
          plot.title = element_text(color="black", size=sizeV, hjust = 0.5),
          axis.title.x = element_text(color="black", size=sizeV),
          axis.title.y = element_text(color=ytitleC, size=sizeV),
          axis.text.x = element_text(color="black", size=sizeV, angle=0),
          axis.text.y = element_text(color=ytextC, size=sizeV, angle=0),
          legend.text = element_text(color="black", size=sizeV, angle=0),
          legend.title = element_text(color="black", size=sizeV, angle=0),
          strip.text = element_text(color="black", size=sizeV, angle=0)) +
        coord_cartesian(ylim = c(0, 1)) +
        scale_color_manual(values = Turbo(out.colors = 17)[1:(length(nums)+1)], name = "") +
        scale_fill_manual(values = Turbo(out.colors = 17)[1:(length(nums)+1)], name = "") +
        scale_x_continuous(breaks = seq(-2, 10, 2)) 
      print(p1)
      
      p1_rmL <- p1 + guides(fill = "none", col = "none") 
      ggsave(paste0(saveDir4, "num-", min(nums), "-", max(nums), "-time.pdf"), p1_rmL, width = 6, height = 6)
      
      if (s == "non-seasonal" & o == "regionally-open" & identical(nums, 19:34)) {
        p1_legend <- get_legend(p1)
        ggsave(paste0(saveDir1, "lg.pdf"), p1_legend, width = 6, height = 6)
      }
    }
  }
}

