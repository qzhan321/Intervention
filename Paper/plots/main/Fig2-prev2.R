rm(list = ls())
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(ggpubr)
})
run <- 4
source(paste0("/home/qizhan/others/PhD/projects/intervention/writings/simulation", run, "/plots/colorFunc.R"))
readDir <- "/project2/pascualmm/QZ/PhD/projects/intervention/natComRevision/files/figures/main/Fig2/prev/"
openness <- c("regionally-open")
seasonality <- c("seasonal")
preIRSList <- list(
  "seasonal regionally-open" = 150
)
IRSDur <- 10
IRSType <- paste0(IRSDur, "yIRS")
postIRS <- 10-IRSDur

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

N <- 10000
sizeV <- 25

saveDir0 <- "/project2/pascualmm/QZ/PhD/projects/intervention/natComRevision/plots/figures/main/Fig2/prev2/"
if (!dir.exists(saveDir0)) {
  dir.create(saveDir0)
}
statsSummaryf <- function(x) {
  r <- quantile(x, probs = c(0.00, 0.25, 0.5, 0.75, 1.00))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}
pointSize <- 1.25
forwardYs <- 9
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
      
      saveDir3 <- paste0(saveDir2, imm, "/")
      if (!dir.exists(saveDir3)) {
        dir.create(saveDir3)
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
      p1 <- ggplot(summaryTable5, aes(timePlot, mP, col = IRS, fill = IRS))+ 
        geom_line(size = pointSize)+
        geom_ribbon(aes(y = mP, ymin = minP, ymax = maxP), alpha = .4, linetype = 0) +
        xlab("Years Since IRS Starts \n(Sustained IRS)") + ylab("Prevalence") +
        theme_bw() + theme(
          plot.title = element_text(color="black", size=sizeV, hjust = 0.5),
          axis.title.x = element_text(color="white", size=sizeV),
          axis.title.y = element_text(color="black", size=sizeV),
          axis.text.x = element_text(color="white", size=sizeV, angle=0),
          axis.text.y = element_text(color="black", size=sizeV, angle=0),
          legend.text = element_text(color="black", size=sizeV, angle=0),
          legend.title = element_text(color="black", size=sizeV, angle=0),
          strip.text = element_text(color="black", size=sizeV, angle=0),
          plot.margin = unit(c(5,10,5,5), "points")) + 
        coord_cartesian(ylim = c(0, 1)) +
        scale_color_manual(values = Turbo(out.colors = 17)[1:(length(nums)+1)], name = "") +
        scale_fill_manual(values = Turbo(out.colors = 17)[1:(length(nums)+1)], name = "") +
        scale_x_continuous(breaks = seq(-2, 10, 2)) 
      print(p1)
      
      summaryTable3 <- summaryTable2 %>% select(time, num, year, Prevalence, CR, rep, IRS)
      summaryTable3$IRS <- factor(summaryTable3$IRS, levels = c("preIRS", paste0("I-", 1:16)))
      if (s == "seasonal") {
        df1 <- summaryTable3 %>% filter(time == (preIRS + forwardYs)*T_YEAR + 300)
      } else {
        df1 <- summaryTable3 %>% filter(time == (preIRS + forwardYs + 1)*T_YEAR)
      }
      summaryTable2 <- summaryTable1 %>% filter(!(time == preIRS * T_YEAR & IRS %in% c(paste0("I-", 1:16))))
      CRs <- summaryTable2 %>% select(IRS, CR) %>% group_by(IRS) %>% summarise(meanCR = mean(CR))
      df2 <- df1 %>% left_join(CRs, by = "IRS")
      df2$IRS <- factor(df2$IRS, levels = paste0("I-", 1:16))
      epiCR <- ggplot(df2, aes(x = meanCR, y = Prevalence, col = IRS)) + 
        stat_summary(fun.data = statsSummaryf, geom="boxplot", position=position_dodge(1)) +
        theme_bw() +
        xlab("Effective Contact Rate \n(per Host, per Year)") + ylab("Prevalence") +
        theme(
          plot.title = element_text(color="black", size=sizeV, hjust = 0.5, vjust = 1),
          axis.title.x = element_text(color="white", size=sizeV),
          axis.title.y = element_text(color="white", size=sizeV),
          axis.text.x = element_text(color="white", size=sizeV, angle=0),
          axis.text.y = element_text(color="white", size=sizeV, angle=0),
          legend.text = element_text(color="black", size=sizeV, angle=0),
          legend.title = element_blank(),
          strip.text = element_text(color="black", size=sizeV, angle=0),
          plot.margin = unit(c(5,22,5,5), "points")) + 
        scale_color_manual(values = Turbo(out.colors = 17)[2:(length(nums) + 1)]) +
        scale_fill_manual(values = Turbo(out.colors = 17)[2:(length(nums) + 1)]) +
        guides(col = "none", fill = "none", alpha = "none") +
        coord_cartesian(ylim = c(0, 1))
      epiCR
      
      p1_rmL <- p1 + guides(fill = "none", col = "none") + rremove("legend") 
      p2_rmL <- epiCR + guides(fill = "none", col = "none") + rremove("legend") 
      
      ggsave(paste0(saveDir3, "num-", min(nums), "-", max(nums), "-time.pdf"), p1_rmL, width = 6, height = 6)
      ggsave(paste0(saveDir3, "num-", min(nums), "-", max(nums), "-cr.pdf"), p2_rmL, width = 6, height = 6)
      
      p1_legend <- get_legend(p1)
      ggsave(paste0(saveDir3, "lg.pdf"), p1_legend, width = 6, height = 6)
    }
  }
}


