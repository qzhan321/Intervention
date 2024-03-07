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
IRSDur <- 10
IRSType <- paste0(IRSDur, "yIRS")
postIRS <- 10-IRSDur

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

N <- 10000
pAll1 <- list()
pAll2 <- list()
sizeV <- 31

saveDir0 <- paste0("/project2/pascualmm/QZ/PhD/projects/intervention/writings/simulation", run, "/plots/")
dir.create(saveDir0)
saveDir1 <- paste0(saveDir0, "prevalenceDynamics/")
dir.create(saveDir1)

for (i in 1:length(seasonality)) {
  s <- seasonality[i]
  saveDir2 <- paste0(saveDir1, s, "/")
  dir.create(saveDir2)
  
  for (j in 1:length(openness)) {
    o <- openness[j]
    saveDir3 <- paste0(saveDir2, o, "/")
    dir.create(saveDir3)
    
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
      
      saveDir4 <- paste0(saveDir3, imm, "/")
      dir.create(saveDir4)
      
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
      summaryTable5$IRS <- factor(summaryTable5$IRS, levels = c("preIRS", paste0("I-", 1:16)))
      p1 <- ggplot(summaryTable5, aes(timePlot, mP, col = IRS, fill = IRS))+ 
        geom_line(size = 1.9)+
        geom_ribbon(aes(y = mP, ymin = minP, ymax = maxP), alpha = .2, linetype = 0) +
        # ggtitle(label) +
        xlab("Years Since IRS Starts \n(Sustained IRS)") + ylab("Prevalence") +
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
          plot.margin = unit(c(5,10,5,5), "points")) + 
        coord_cartesian(ylim = c(0, 1)) +
        scale_color_manual(values = Turbo(out.colors = 17)[1:(length(nums)+1)]) +
        scale_fill_manual(values = Turbo(out.colors = 17)[1:(length(nums)+1)]) +
        scale_x_continuous(breaks = seq(-2, 10, 2))
      p1 <- p1 + guides(fill = "none", col = "none")
      print(p1)
      index <- length(pAll1)
      pAll1[[index + 1]] <- p1
      
      
      forwardYs <- 9
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
      if (s == "seasonal" & o == "regionally-open" & k == 1) {
        yt <- "white"
      } else {
        yt <- "black"
      }
      epiCR<-ggplot(df2, aes(x=meanCR, y=mP, col = IRS)) + 
        geom_point(size = 5) + 
        geom_pointrange(aes(ymin=minP, ymax=maxP)) +
        # geom_line(data = dfspline, aes(x=x, y=y), color = "gray88", size = 1.9) + 
        # ggtitle(label) +
        theme_bw() +
        xlab("Effective Contact Rate \n(per Host, per Year)") + ylab("Prevalence") +
        theme(
          plot.title = element_text(color="black", size=sizeV, hjust = 0.5, vjust = 1),
          axis.title.x = element_text(color="black", size=sizeV),
          axis.title.y = element_text(color=yt, size=sizeV),
          axis.text.x = element_text(color="black",
                                     size=sizeV, angle=0),
          axis.text.y = element_text(color=yt,
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
      
      ggsave(paste0(saveDir4, "num-", min(nums), "-", max(nums), "-time.pdf"), p1_rmL, width = 6, height = 6)
      ggsave(paste0(saveDir4, "num-", min(nums), "-", max(nums), "-cr.pdf"), p2_rmL, width = 6, height = 6)
      
    }
  }
}


