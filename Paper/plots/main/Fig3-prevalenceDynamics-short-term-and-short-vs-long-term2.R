rm(list = ls())
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggpubr)
  library(gridExtra)
})
run <- 4
source("/home/qizhan/others/PhD/projects/intervention/writings/simulation4/plots/colorFunc.R")
saveDir0 <- "/project2/pascualmm/QZ/PhD/projects/intervention/natComRevision/plots/figures/main/Fig3/prevalenceDynamicsShortTerm2/"
if (!dir.exists(saveDir0)) {
  dir.create(saveDir0)
}
saveDir1 <- "/project2/pascualmm/QZ/PhD/projects/intervention/natComRevision/plots/figures/main/Fig3/prevalenceDynamicsShortTermVsLongTerm2/"
if (!dir.exists(saveDir1)) {
  dir.create(saveDir1)
}
readDir <- "/project2/pascualmm/QZ/PhD/projects/intervention/natComRevision/files/figures/main/Fig2/prev/"
numsPlotList <- list(
  "seasonal regionally-open" = list(201:209)
)
labelsList <- list(
  "seasonal regionally-open" = list("Seasonal \nRegionally-open \nBaseline Migration \nMedium Pool")
)
immsList <- list(
  "seasonal regionally-open" = list("Three-year")
)
T_YEAR <- 360
IRSDur <- 10
postIRS <- 10 - IRSDur
IRSType <- paste0(IRSDur, "yIRS")
seasonality <- c("seasonal")
openness <- c("regionally-open")
preIRSList <- list(
  "seasonal regionally-open" = 150
)
sizeV <- 20
pointSize <- 3.3
statsSummaryf <- function(x) {
  r <- quantile(x, probs = c(0.00, 0.25, 0.5, 0.75, 1.00))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}
for (i in 1:length(seasonality)) {
  s <- seasonality[i]
  saveDir2 <- paste0(saveDir0, s, "/")
  if (!dir.exists(saveDir2)) {
    dir.create(saveDir2)
  }
  saveDir3 <- paste0(saveDir1, s, "/")
  if (!dir.exists(saveDir3)) {
    dir.create(saveDir3)
  }
  for (j in 1:length(openness)) {
    o <- openness[j]
    saveDir4 <- paste0(saveDir2, o, "/")
    if (!dir.exists(saveDir4)) {
      dir.create(saveDir4)
    }
    saveDir5 <- paste0(saveDir3, o, "/")
    if (!dir.exists(saveDir5)) {
      dir.create(saveDir5)
    }
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
      if (!dir.exists(saveDir6)) {
        dir.create(saveDir6)
      }
      saveDir7 <- paste0(saveDir5, imm, "/")
      if (!dir.exists(saveDir7)) {
        dir.create(saveDir7)
      }
      
      summaryTable1 <- summaryTableCombined %>% filter(num %in% nums)
      summaryTable2 <- summaryTable1 %>% filter(time > (preIRS-3)*T_YEAR, time <= (preIRS + IRSDur + postIRS)*T_YEAR)
      summaryTable3 <- summaryTable2 %>% select(time, num, year, Prevalence, CR, rep, IRS)
      summaryTable3$IRS <- factor(summaryTable3$IRS, levels = c("preIRS", paste0("I-", 1:16)))
      
      if (s == "seasonal") {
        t1 <- (preIRS + 1)*T_YEAR + 300
        t2 <- (preIRS + 9)*T_YEAR + 300
      } else {
        t1 <- (preIRS + 2)*T_YEAR 
        t2 <- (preIRS + 10)*T_YEAR 
      }
      df1 <- summaryTable3 %>% filter(time == t1) 
      colnames(df1)[4] <- paste0("short_term_rebound_", colnames(df1)[4])
      df2 <- summaryTable3 %>% filter(time == t2) 
      colnames(df2)[4] <- paste0("longer_term_rebound_", colnames(df2)[4])
      
      df1_2 <- df1 %>% group_by(IRS) %>%
        summarise(x.min = min(short_term_rebound_Prevalence), 
                  x.lower = quantile(short_term_rebound_Prevalence, 0.25),
                  x.middle = quantile(short_term_rebound_Prevalence, 0.50),
                  x.upper = quantile(short_term_rebound_Prevalence, 0.75),
                  x.max = max(short_term_rebound_Prevalence)) %>% mutate(category = IRS) %>% select(-IRS)
      df2_2 <- df2 %>% group_by(IRS) %>%
        summarise(y.min = min(longer_term_rebound_Prevalence), 
                  y.lower = quantile(longer_term_rebound_Prevalence, 0.25),
                  y.middle = quantile(longer_term_rebound_Prevalence, 0.50),
                  y.upper = quantile(longer_term_rebound_Prevalence, 0.75),
                  y.max = max(longer_term_rebound_Prevalence)) %>% mutate(category = IRS) %>% select(-IRS)
      df <- df1_2 %>% left_join(df2_2, by = "category")
      
      p1 <- ggplot(df, aes(fill = category, color = category)) +
        # 2D box defined by the Q1 & Q3 values in each dimension, with outline
        geom_rect(aes(xmin = x.lower, xmax = x.upper, ymin = y.lower, ymax = y.upper), fill = NA, alpha = 0.3, size = pointSize/3) +
        geom_rect(aes(xmin = x.lower, xmax = x.upper, ymin = y.lower, ymax = y.upper), fill = NA, alpha = 0.3, size = pointSize/3) +
        
        # whiskers for x-axis dimension with ends
        geom_segment(aes(x = x.min, y = y.middle, xend = x.max, yend = y.middle), size = pointSize/3) + #whiskers
        geom_segment(aes(x = x.min, y = y.lower, xend = x.min, yend = y.upper), size = pointSize/3) + #lower end
        geom_segment(aes(x = x.max, y = y.lower, xend = x.max, yend = y.upper), size = pointSize/3) + #upper end
        
        # whiskers for y-axis dimension with ends
        geom_segment(aes(x = x.middle, y = y.min, xend = x.middle, yend = y.max), size = pointSize/3) + #whiskers
        geom_segment(aes(x = x.lower, y = y.min, xend = x.upper, yend = y.min), size = pointSize/3) + #lower end
        geom_segment(aes(x = x.lower, y = y.max, xend = x.upper, yend = y.max), size = pointSize/3) + #upper end
        
        # outliers
        # geom_point(data = df.outliers, aes(x = x.outliers, y = y.middle), size = 1, shape = 1) + # x-direction
        # geom_point(data = df.outliers, aes(x = x.middle, y = y.outliers), size = 1, shape = 1) + # y-direction
        
        xlab("Prevalence \n(Short-term Rebound)") + ylab("Prevalence \n(Long-term Rebound)") +
        theme_bw() + theme(
          plot.title = element_text(color="black", size=sizeV, hjust = 0.5),
          axis.title.x = element_text(color="black", size=sizeV),
          axis.title.y = element_text(color="black", size=sizeV),
          axis.text.x = element_text(color="black", size=sizeV, angle=0),
          axis.text.y = element_text(color="black", size=sizeV, angle=0),
          legend.text = element_text(color="black", size=sizeV, angle=0),
          legend.title = element_text(color="black", size=sizeV, angle=0),
          strip.text = element_text(color="black", size=sizeV, angle=0)) +
        scale_color_manual(values = Turbo(out.colors = 17)[2:(length(nums) + 1)]) +
        scale_fill_manual(values = Turbo(out.colors = 17)[2:(length(nums) + 1)]) +
        coord_cartesian(xlim = c(0,1), ylim = c(0,1)) + 
        guides(fill = "none", col = "none") +
        scale_x_continuous(expand = expansion(mult = 0.1)) + 
        scale_y_continuous(expand = expansion(mult = 0.1))
      print(p1)
      
      forwardYs <- 1
      if (s == "seasonal") {
        df1 <- summaryTable3 %>% filter(time == (preIRS + forwardYs)*T_YEAR + 300)
      } else {
        df1 <- summaryTable3 %>% filter(time == (preIRS + forwardYs + 1)*T_YEAR)
      }
      summaryTable2 <- summaryTable1 %>% filter(!(time == preIRS * T_YEAR & IRS %in% c(paste0("I-", 1:16))))
      CRs <- summaryTable2 %>% select(IRS, CR) %>% group_by(IRS) %>% summarise(meanCR = mean(CR))
      df2 <- df1 %>% left_join(CRs, by = "IRS")
      df2$IRS <- factor(df2$IRS, levels = paste0("I-", 1:16))
      epiCR <- ggplot(df2, aes(x = meanCR, y = Prevalence, col = IRS, alpha = 0.3)) + 
        stat_summary(fun.data = statsSummaryf, geom="boxplot", 
                     position=position_dodge(1)) +
        # geom_errorbar(aes(ymin = minP, ymax = maxP), width = pointSize/3, size = pointSize/4) +
        # geom_pointrange(aes(ymin = minP, ymax = maxP), size = pointSize) +
        theme_bw() +
        xlab("Effective Contact Rate \n(per Host, per Year)") + 
        ylab("Prevalence \n(Short-term Rebound)") +
        theme(
          plot.title = element_text(color="black", size = sizeV, hjust = 0.5, vjust = 1),
          axis.title.x = element_text(color="black", size = sizeV),
          axis.title.y = element_text(color="black", size = sizeV),
          axis.text.x = element_text(color="black", size = sizeV, angle=0),
          axis.text.y = element_text(color="black", size = sizeV, angle=0),
          legend.text = element_text(color="black", size = sizeV, angle=0),
          legend.title = element_blank(),
          strip.text = element_text(color="black", size = sizeV, angle=0),
          plot.margin = unit(c(5,22,5,5), "points"),
          legend.position = c(0.9, 0.04)) + 
        scale_color_manual(values = Turbo(out.colors = 17)[2:(length(nums) + 1)]) +
        scale_fill_manual(values = Turbo(out.colors = 17)[2:(length(nums) + 1)]) +
        guides(col = "none", fill = "none", alpha = "none") +
        coord_cartesian(ylim = c(0, 1))
      epiCR
      
      p1_rmL <- p1 + guides(fill = "none", col = "none") + rremove("legend") 
      p2_rmL <- epiCR + guides(fill = "none", col = "none") + rremove("legend") 
      
      ggsave(paste0(saveDir6, "num-", min(nums), "-", max(nums), "-cr.pdf"), p2_rmL, width = 4, height = 4)
      ggsave(paste0(saveDir7, "num-", min(nums), "-", max(nums), ".pdf"), p1_rmL, width = 6.45/1.5, height = 4)
    }
  }
}

