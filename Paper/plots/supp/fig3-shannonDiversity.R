rm(list = ls())
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(ggpubr)
})
scaleFUN <- function(x) sprintf("%.2f", x)
run <- 4
source(paste0("/home/qizhan/others/PhD/projects/intervention/writings/simulation", run, "/plots/colorFunc.R"))
readDir <- paste0("/project2/pascualmm/QZ/PhD/projects/intervention/files", run, "/actualRuns/diversity/")
readDirCR <- paste0("/project2/pascualmm/QZ/PhD/projects/intervention/files", run, "/actualRuns/epi/")
state <- "true"
category <- "shannonDiversity"
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
IRSDur <- 10
IRSType <- paste0(IRSDur, "yIRS")
postIRS <- 10 - IRSDur

numsPlotList <- list(
  "seasonal closed" = list(101:110),
  "seasonal semi-open" = list(201:211, 301:311),
  "seasonal regionally-open" = list(219:229, 237:245, 255:265),
  "non-seasonal closed" = list(101:111),
  "non-seasonal semi-open" = list(201:211,301:311),
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
saveDir0 <- "/project2/pascualmm/QZ/PhD/projects/intervention/natComRevision/plots/figures/supp/fig3/shannonDiversity/"
if (!dir.exists(saveDir0)) {
  dir.create(saveDir0)
}

T_YEAR <- 360
state <- "true"
N <- 10000
sizeV <- 31
nums_w_reps_list <- list(
  "seasonal closed" = list(107:110),
  "seasonal semi-open" = list(208:211, 308:311),
  "seasonal regionally-open" = list(226:229, 242:245, 262:265),
  "non-seasonal closed" = list(108:110),
  "non-seasonal semi-open" = list(208:210, 307:310),
  "non-seasonal regionally-open" = list(207:209, 225:227, 242:245, 261:263)
)
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
    nums_w_reps <- nums_w_reps_list[[paste0(s, " ", o)]]
    
    for (k in 1:length(numsPlot)) {
      nums <- numsPlot[[k]]
      label <- labels[[k]]
      imm <- imms[[k]]
      saveDir3 <- paste0(saveDir2, imm, "/")
      if (!dir.exists(saveDir3)) {
        dir.create(saveDir3)
      }
      nums_w_reps_single <- nums_w_reps[[k]]
      shannon2VAll <- NULL
      for (m in 1:length(nums)) {
        num <- nums[m]
        if (num %in% nums_w_reps_single) {
          reps <- 0:2
        } else {
          reps <- 0
        }
        for (r in reps) {
          file <- paste0(readDir, category, "/", s, "/", o, "/", state, "/", category, "-", num, "-r", r, ".RData")
          load(file)
          df2 <- df %>% mutate(IRS = case_when(
            time<=preIRS*T_YEAR ~ "preIRS",
            time>preIRS*T_YEAR ~ paste0("I-", num - min(nums) + 1)
          ))
          df2App <- df2[df2$time==preIRS * T_YEAR | df2$time == (preIRS-1)*T_YEAR + 300,,drop=F]
          df2App$IRS = paste0("I-", num - min(nums) + 1)
          df3 <- rbind(df2, df2App)
          
          if (s == "seasonal") {
            tsdiff <- setdiff(seq(preIRS - 3, preIRS + IRSDur + postIRS - 1)*T_YEAR + 300, unique(df3$time))
          } else {
            tsdiff <- setdiff(seq(preIRS - 2, preIRS + IRSDur + postIRS)*T_YEAR, unique(df3$time))
          }
          
          if (length(tsdiff) > 0) {
            stopifnot(min(tsdiff) > preIRS * T_YEAR)
            df3App <- data.frame("num" = num, "shannon" = 0, "time" = tsdiff, "rep" = r, "IRS" = paste0("I-", num - min(nums) + 1))
            df4 <- rbind(df3, df3App)
          } else {
            df4 <- df3
          }
          shannon2VAll <- rbind(shannon2VAll, df4)
        }
      }
      
      file <- paste0(readDirCR, s, "/", o, "/", "summaryInfoTable-", IRSType, ".RData")
      load(file)
      summaryTable1 <- summaryTableCombined %>% filter(num %in% nums)
      forwardYs <- 9
      summaryTable2 <- summaryTable1 %>% filter(!(time == preIRS * T_YEAR & IRS %in% c(paste0("I-", 1:16))))
      CRs <- summaryTable2 %>% select(IRS, CR) %>% group_by(IRS) %>% summarise(meanCR = mean(CR))
      SD <- shannon2VAll %>% left_join(CRs, by = "IRS")
      SD2 <- SD %>% group_by(time, IRS, meanCR) %>% 
        summarise(mD=mean(shannon), maxD=max(shannon), minD=min(shannon))
      SD2 <- SD2 %>% mutate(timePlot = time/T_YEAR -preIRS)
      SD2$IRS <- factor(SD2$IRS, levels = c("preIRS", paste0("I-", 1:16)))
      if (s == "seasonal" & o == "closed") {
        ytitleC <- "black"
      } else if (s == "seasonal" & o == "regionally-open" & identical(nums, 255:265)) {
        ytitleC <- "black"
      } else if (s == "non-seasonal" & o == "regionally-open" & identical(nums, 219:228)) {
        ytitleC <- "black"
      } else {
        ytitleC <- "white"
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
      p1 <- ggplot(SD2, aes(timePlot, mD, col = IRS, fill = IRS)) + 
        geom_line(size = 1.9)+
        geom_ribbon(aes(y = mD, ymin = minD, ymax = maxD), alpha = .4, linetype = 0) +
        ggtitle(label) +
        xlab("Years Since IRS Starts \n(Sustained IRS)") + ylab("Shannon Diversity") +
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
        coord_cartesian(ylim = c(min(SD2$minD) - 1, max(SD2$maxD) + 0.5)) +
        scale_color_manual(values = Turbo(out.colors = 17)[1:(length(nums)+1)]) +
        scale_fill_manual(values = Turbo(out.colors = 17)[1:(length(nums)+1)]) +
        scale_x_continuous(breaks = seq(-2, 10, 2)) + guides(fill = "none", col = "none")
      p1
      
      if (s == "seasonal") {
        SD3 <- SD2 %>% filter(time %in% c((preIRS + forwardYs)*T_YEAR + 300))
      } else {
        SD3 <- SD2 %>% filter(time %in% c((preIRS + forwardYs + 1)*T_YEAR))
      }
      p2<-ggplot(SD3, aes(x = meanCR, y = mD, col = IRS)) + 
        geom_point(size = 5) +
        geom_pointrange(aes(ymin = minD, ymax = maxD)) + 
        ggtitle(label) +
        xlab("Effective Contact Rate \n(per Host, per Year)") + ylab("Shannon Diveristy") +
        theme_bw() +
        theme(
          plot.title = element_text(color="black", size=sizeV/1.25, hjust = 0.5, vjust = 1),
          axis.title.x = element_text(color=xtitleC, size=sizeV),
          axis.title.y = element_text(color=ytitleC, size=sizeV),
          axis.text.x = element_text(color="black", size=sizeV, angle=0),
          axis.text.y = element_text(color="black", size=sizeV, angle=0),
          legend.text = element_text(color="black", size=sizeV, angle=0),
          legend.title = element_blank(),
          strip.text = element_text(color="black", size=sizeV, angle=0),
          plot.margin = unit(c(5,22,5,5), "points")) + 
        scale_color_manual(values = Turbo(out.colors = 17)[2:(length(nums) + 1)]) +
        scale_fill_manual(values = Turbo(out.colors = 17)[2:(length(nums) + 1)]) +
        guides(col = "none", fill = "none") +
        coord_cartesian(ylim = c(min(SD2$minD) - 1, max(SD2$maxD) + 0.5)) +
        scale_y_continuous(labels = scaleFUN)
      print(p2)
      
      ggsave(paste0(saveDir3, "num-", min(nums), "-", max(nums), "-time.pdf"), p1, width = 6, height = 6)
      ggsave(paste0(saveDir3, "num-", min(nums), "-", max(nums), "-cr.pdf"), p2, width = 6.45, height = 7.45)
    }
  }
}
