rm(list = ls())
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(ggpubr)
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
  "seasonal semi-open" = list(201:211, 301:311),
  "seasonal regionally-open" = list(219:229, 237:245, 255:265),
  "non-seasonal closed" = list(101:111),
  "non-seasonal semi-open" = list(201:211, 301:311),
  "non-seasonal regionally-open" = list(201:210, 219:228, 237:246, 255:264)
)
labelsList <- list(
  "seasonal closed" = list("Seasonal \nClosed \n \n"),
  "seasonal semi-open" = list("Seasonal \nSemi-open \nBaseline Migration \n", 
                              "Seasonal \nSemi-open \nHigh Migration \n"),
  "seasonal regionally-open" = list("Seasonal \nRegionally-open \nBaseline Migration \nLarge Pool",
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
  "seasonal regionally-open" = list("Three-year", "Three-year", "Three-year"),
  "non-seasonal closed" = list("Three-year"),
  "non-seasonal semi-open" = list("Three-year", "Three-year"),
  "non-seasonal regionally-open" = list("Three-year", "Three-year", "Three-year", "Three-year")
)

T_YEAR <- 360

N <- 10000
sizeV <- 30

saveDir0 <- "/project2/pascualmm/QZ/PhD/projects/intervention/natComRevision/plots/figures/supp/fig1/prev/"
if (!dir.exists(saveDir0)) {
  dir.create(saveDir0)
}

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
      if (s == "seasonal" & o == "closed") {
        ytitleC <- "black"
        ytextC <- "black"
      } else if (s == "seasonal" & o == "regionally-open" & identical(nums, 255:265)) {
        ytitleC <- "black"
        ytextC <- "black"
      } else if (s == "non-seasonal" & o == "regionally-open" & identical(nums, 219:228)) {
        ytitleC <- "black"
        ytextC <- "black"
      } else {
        ytitleC <- "white"
        ytextC <- "white"
      }
      if (s == "seasonal" & o == "semi-open" & identical(nums, 301:311)) {
        xtitleC <- "black"
      } else if (s == "non-seasonal" & o == "semi-open" & identical(nums, 201:211)) {
        xtitleC <- "black"
      } else if (s == "non-seasonal" & o == "regionally-open" & identical(nums, 255:264)) {
        xtitleC <- "black"
      } else {
        xtitleC <- "white"
      }
      p1 <- ggplot(summaryTable5, aes(timePlot, mP, col = IRS, fill = IRS))+ 
        geom_line(size = 1.9)+
        geom_ribbon(aes(y = mP, ymin = minP, ymax = maxP), alpha = .4, linetype = 0) +
        ggtitle(label) +  
        xlab("Years Since IRS Starts \n(Sustained IRS)") + ylab("Prevalence") +
        theme_bw() + theme(
          plot.title = element_text(color="black", size=sizeV/1.25, hjust = 0.5),
          axis.title.x = element_text(color=xtitleC, size=sizeV),
          axis.title.y = element_text(color=ytitleC, size=sizeV),
          axis.text.x = element_text(color="black", size=sizeV, angle=0),
          axis.text.y = element_text(color=ytextC, size=sizeV, angle=0),
          legend.text = element_text(color="black", size=sizeV, angle=0),
          legend.title = element_text(color="black", size=sizeV, angle=0),
          strip.text = element_text(color="black", size=sizeV, angle=0),
          plot.margin = unit(c(5,10,5,5), "points")) + 
        coord_cartesian(ylim = c(0, 1)) +
        scale_color_manual(values = Turbo(out.colors = 17)[1:(length(nums)+1)], name = "") +
        scale_fill_manual(values = Turbo(out.colors = 17)[1:(length(nums)+1)], name = "") +
        scale_x_continuous(breaks = seq(-2, 10, 2)) 
      print(p1)
      
      
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
      epiCR <- ggplot(df2, aes(x=meanCR, y=mP, col = IRS)) + 
        geom_point(size = 5) + 
        geom_pointrange(aes(ymin=minP, ymax=maxP)) +
        ggtitle(label) +  
        xlab("Effective Contact Rate \n(per Host, per Year)") + ylab("Prevalence") +
        theme_bw() +
        theme(
          plot.title = element_text(color="black", size=sizeV/1.25, hjust = 0.5, vjust = 1),
          axis.title.x = element_text(color=xtitleC, size=sizeV),
          axis.title.y = element_text(color=ytitleC, size=sizeV),
          axis.text.x = element_text(color="black", size=sizeV, angle=0),
          axis.text.y = element_text(color=ytextC, size=sizeV, angle=0),
          legend.text = element_text(color="black", size=sizeV, angle=0),
          legend.title = element_blank(),
          strip.text = element_text(color="black", size=sizeV, angle=0),
          plot.margin = unit(c(5,22,5,5), "points")) + 
        scale_color_manual(values = Turbo(out.colors = 17)[2:(length(nums) + 1)]) +
        scale_fill_manual(values = Turbo(out.colors = 17)[2:(length(nums) + 1)]) +
        guides(col = "none", fill = "none") +
        coord_cartesian(ylim = c(0, 1))
      epiCR
      
      ggsave(paste0(saveDir3, "num-", min(nums), "-", max(nums), "-time.pdf"), p1 + guides(fill = "none", col = "none"), width = 6, height = 6)
      ggsave(paste0(saveDir3, "num-", min(nums), "-", max(nums), "-cr.pdf"), epiCR, width = 6.45, height = 7.45)
      
      if (s == "seasonal" & o == "regionally-open" & identical(nums, 219:229)) {
        p1_legend <- get_legend(p1)
        ggsave(paste0(saveDir3, "num-", min(nums), "-", max(nums), "-lg.pdf"), p1_legend, width = 6, height = 6)
      }
    }
  }
}


