rm(list = ls())
suppressPackageStartupMessages({
  library(RSQLite)
  library(reshape2)
  library(vegan)
  library(tidyverse)
  library(igraph)
  library(cowplot)
  library(scales)
  # library(adegenet)
  library(factoextra)
  library(ggplot2)
  library(ggpubr)
  library(dplyr)
  library(RSQLite)
  library(reshape2)
  library(vegan)
  library(tidyverse)
  library(cowplot)
  library(gridExtra)
  library(grid)
})
run <- 4
source("/home/qizhan/others/PhD/projects/intervention/writings/simulation3/plots/colorFunc.R")
saveDir0 <- paste0("/project2/pascualmm/QZ/PhD/projects/intervention/writings/simulation", run, "/plots/prevalenceDynamicsShortTerm/")
dir.create(saveDir0)
saveDir1 <- paste0("/project2/pascualmm/QZ/PhD/projects/intervention/writings/simulation", run, "/plots/prevalenceDynamicsShortTermVsLongTerm/")
dir.create(saveDir1)
readDir <- paste0("/project2/pascualmm/QZ/PhD/projects/intervention/files", run, "/actualRuns/epi/")
numsPlotList <- list(
  "seasonal closed" = list(101:110),
  "seasonal semi-open" = list(201:211,  
                              301:311),
  "seasonal regionally-open" = list(201:209,219:229,
                                    237:245,255:265),
  "non-seasonal closed" = list(101:111),
  "non-seasonal semi-open" = list(201:211,
                                  301:311),
  "non-seasonal regionally-open" = list(201:210,219:228,
                                        237:246,255:264)
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
  "seasonal semi-open" = list("Three-year",
                              "Three-year"),
  "seasonal regionally-open" = list("Three-year", "Three-year", "Three-year", "Three-year"),
  "non-seasonal closed" = list("Three-year"),
  "non-seasonal semi-open" = list("Three-year",
                                  "Three-year"),
  "non-seasonal regionally-open" = list("Three-year", "Three-year", "Three-year", "Three-year")
)
T_YEAR <- 360
IRSDur <- 10
postIRS <- 10 - IRSDur
IRSType <- paste0(IRSDur, "yIRS")
seasonality <- c("seasonal", "non-seasonal")
openness <- c("closed", "semi-open", "regionally-open")
preIRSList <- list(
  "seasonal closed" = 200,
  "seasonal semi-open" = 150, 
  "seasonal regionally-open" = 150,
  "non-seasonal closed" = 200,
  "non-seasonal semi-open" = 150,
  "non-seasonal regionally-open" = 150
)
sizeV <- 30
pAll1 <- list()
pAll2 <- list()
for (i in 1:length(seasonality)) {
  s <- seasonality[i]
  saveDir2 <- paste0(saveDir0, s, "/")
  dir.create(saveDir2)
  saveDir3 <- paste0(saveDir1, s, "/")
  dir.create(saveDir3)
  
  for (j in 1:length(openness)) {
    o <- openness[j]
    saveDir4 <- paste0(saveDir2, o, "/")
    dir.create(saveDir4)
    saveDir5 <- paste0(saveDir3, o, "/")
    dir.create(saveDir5)
    
    file <- paste0(readDir, s, "/", o, "/", "summaryInfoTable-", IRSType, ".RData")
    load(file)
    preIRS <- preIRSList[[paste0(s, " ", o)]]
    numsPlot <- numsPlotList[[paste0(s, " ", o)]]
    labels <- labelsList[[paste0(s, " ", o)]]
    imms <- immsList[[paste0(s, " ", o)]]
    for (k in 1:length(numsPlot)) {
      nums <- numsPlot[[k]]
      label <- labels[[k]]
      imm <- imms[[k]]
      
      saveDir6 <- paste0(saveDir4, imm, "/")
      dir.create(saveDir6)
      saveDir7 <- paste0(saveDir5, imm, "/")
      dir.create(saveDir7)
      
      summaryTable1 <- summaryTableCombined %>% filter(num %in% nums)
      summaryTable2 <- summaryTable1 %>% filter(time > (preIRS-3)*T_YEAR, time <= (preIRS + IRSDur + postIRS)*T_YEAR)
      summaryTable3 <- summaryTable2 %>% dplyr::group_by(time, IRS) %>%
        dplyr::summarise(mP = mean(Prevalence), vP = sd(Prevalence), 
                         maxP = max(Prevalence), minP = min(Prevalence),
                         mD = mean(n_circulating_genes), vD = sd(n_circulating_genes),
                         maxD = max(n_circulating_genes), minD = min(n_circulating_genes),
                         mMOI = mean(MOI, na.rm = T), vMOI = sd(MOI, na.rm = T))
      summaryTable4 <- summaryTable3 %>% mutate(vP = ifelse(is.na(vP), 0, vP), 
                                                vD = ifelse(is.na(vD), 0, vD),
                                                vMOI = ifelse(is.na(vMOI), 0, vMOI))
      summaryTable5 <- summaryTable4 %>% mutate(timePlot = time/T_YEAR - preIRS) 
      summaryTable5$IRS <- factor(summaryTable5$IRS, levels = c("preIRS", paste0("I-", 1:18)))
      
      if (s == "seasonal") {
        t1 <- (preIRS + 1)*T_YEAR + 300
        t2 <- (preIRS + 9)*T_YEAR + 300
      } else {
        t1 <- (preIRS + 2)*T_YEAR 
        t2 <- (preIRS + 10)*T_YEAR 
      }
      df1 <- summaryTable5 %>% dplyr::filter(time == t1) %>% dplyr::ungroup() %>%
        dplyr::select(IRS, mP, minP, maxP)
      df1 <- df1 %>% group_by(IRS) %>% summarise(mP = unique(mP),
                                                 minP = unique(minP),
                                                 maxP = unique(maxP))
      colnames(df1)[2:4] <- paste0("short_term_rebound_", colnames(df1)[2:4])
      
      df2 <- summaryTable5 %>% dplyr::filter(time == t2) %>% dplyr::ungroup() %>%
        dplyr::select(IRS, mP, minP, maxP)
      df2 <- df2 %>% group_by(IRS) %>% summarise(mP = unique(mP),
                                                 minP = unique(minP),
                                                 maxP = unique(maxP))
      colnames(df2)[2:4] <- paste0("longer_term_rebound_", colnames(df2)[2:4])
      
      df12 <- left_join(df1, df2, by = c("IRS"))
      
      p1 <- ggplot(df12, aes(x=short_term_rebound_mP, y=longer_term_rebound_mP, 
                             col = IRS, fill = IRS))+ 
        geom_point(size = 4)+
        geom_errorbarh(aes(xmin=short_term_rebound_minP, xmax=short_term_rebound_maxP)) +
        geom_pointrange(aes(ymin=longer_term_rebound_minP, ymax=longer_term_rebound_maxP)) +
        xlab("Prevalence \n(Short-term Rebound)") + ylab("Prevalence \n(Long-term Rebound)") +
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
                                    size=sizeV, angle=0)) +
        scale_color_manual(values = Turbo(out.colors = 17)[2:(length(nums) + 1)]) +
        scale_fill_manual(values = Turbo(out.colors = 17)[2:(length(nums) + 1)]) +
        coord_cartesian(xlim = c(0,1), ylim = c(0,1)) + 
        guides(fill = "none", col = "none") +
        scale_x_continuous(expand = expansion(mult = 0.1)) + 
        scale_y_continuous(expand = expansion(mult = 0.1))
      print(p1)
      index <- length(pAll1)
      pAll1[[index + 1]] <- p1
      
      forwardYs <- 1
      if (s == "seasonal") {
        df1 <- summaryTable5 %>% filter(time == (preIRS + forwardYs)*T_YEAR + 300)
      } else {
        df1 <- summaryTable5 %>% filter(time == (preIRS + forwardYs + 1)*T_YEAR)
      }
      summaryTable2 <- summaryTable1 %>% filter(!(time == preIRS * T_YEAR & IRS %in% c(paste0("I-", 1:16))))
      CRs <- summaryTable2 %>% select(IRS, CR) %>% group_by(IRS) %>% summarise(meanCR = mean(CR))
      df2 <- df1 %>% left_join(CRs, by = "IRS")
      df2$IRS <- factor(df2$IRS, levels = paste0("I-", 1:16))
      dfspline <- as.data.frame(spline(df2$meanCR, df2$mP))
      epiCR<-ggplot(df2, aes(x=meanCR, y=mP, col = IRS)) + 
        geom_point(size = 5) + 
        geom_pointrange(aes(ymin=minP, ymax=maxP)) +
        # geom_line(data = dfspline, aes(x=x, y=y), color = "gray88", size = 1.9) + 
        # ggtitle(label) +
        theme_bw() +
        xlab("Effective Contact Rate \n(per Host, per Year)") + ylab("Prevalence \n(Short-term Rebound)") +
        theme(
          plot.title = element_text(color="black", size=sizeV, hjust = 0.5, vjust = 1),
          axis.title.x = element_text(color="black", size=sizeV),
          axis.title.y = element_text(color="black", size=sizeV),
          axis.text.x = element_text(color="black",
                                     size=sizeV, angle=0),
          axis.text.y = element_text(color="black",
                                     size=sizeV, angle=0),
          legend.text = element_text(color="black",
                                     size=sizeV, angle=0),
          legend.title = element_blank(),
          strip.text = element_text(color="black",
                                    size=sizeV, angle=0),
          plot.margin = unit(c(5,22,5,5), "points"),
          legend.position = c(0.9, 0.04)) + 
        scale_color_manual(values = Turbo(out.colors = 17)[2:(length(nums) + 1)]) +
        scale_fill_manual(values = Turbo(out.colors = 17)[2:(length(nums) + 1)]) +
        guides(col = "none", fill = "none") +
        coord_cartesian(ylim = c(0, 1))
      # scale_x_continuous(breaks = seq(2, 13, 2))) 
      epiCR
      index <- length(pAll2)
      pAll2[[index + 1]] <- epiCR
      
      p1_rmL <- p1 + guides(fill = "none", col = "none") + rremove("legend") 
      p2_rmL <- epiCR + guides(fill = "none", col = "none") + rremove("legend") 
      
      ggsave(paste0(saveDir7, "num-", min(nums), "-", max(nums), ".pdf"), p1_rmL, width = 6.45, height = 6)
      ggsave(paste0(saveDir6, "num-", min(nums), "-", max(nums), "-cr.pdf"), p2_rmL, width = 6, height = 6)
      
    }
  }
}

