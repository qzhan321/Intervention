rm(list = ls())
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(ggpubr)
})
scaleFUN <- function(x) sprintf("%.2f", x)
run <- 4
source(paste0("/home/qizhan/others/PhD/projects/intervention/writings/simulation", run, "/plots/colorFunc.R"))
readDir <- "/project2/pascualmm/QZ/PhD/projects/intervention/natComRevision/files/figures/main/Fig2/"
readDirCR <- "/project2/pascualmm/QZ/PhD/projects/intervention/natComRevision/files/figures/main/Fig2/prev/"
state <- "true"
category <- "shannonDiversity"
seasonality <- "seasonal"
openness <- "regionally-open"
preIRSList <- list(
  "seasonal regionally-open" = 150
)
IRSDur <- 10
IRSType <- paste0(IRSDur, "yIRS")
postIRS <- 10 - IRSDur

numsPlotList <- list(
  "seasonal regionally-open" = list(201:209)
)
labelsList <- list(
  "seasonal regionally-open" = list("Seasonal \nRegionally-open \nBaseline Migration \nMedium Pool")
)
immsList <- list(
  "seasonal regionally-open" = list("Three-year")
)

saveDir0 <- "/project2/pascualmm/QZ/PhD/projects/intervention/natComRevision/plots/figures/main/Fig2/"
if (!dir.exists(saveDir0)) {
  dir.create(saveDir0)
}
saveDir1 <- paste0(saveDir0, "shannonDiversity2/")
if (!dir.exists(saveDir1)) {
  dir.create(saveDir1)
}

forwardYs <- 9
T_YEAR <- 360
state <- "true"
N <- 10000
sizeV <- 25
nums_w_reps_list <- list(
  "seasonal regionally-open" = list(201:209)
)
statsSummaryf <- function(x) {
  r <- quantile(x, probs = c(0.00, 0.25, 0.5, 0.75, 1.00))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
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
    preIRS <- preIRSList[[paste0(s, " ", o)]]
    numsPlot <- numsPlotList[[paste0(s, " ", o)]]
    labels <- labelsList[[paste0(s, " ", o)]]
    imms <- immsList[[paste0(s, " ", o)]]
    nums_w_reps <- nums_w_reps_list[[paste0(s, " ", o)]]
    for (k in 1:length(numsPlot)) {
      nums <- numsPlot[[k]]
      label <- labels[[k]]
      imm <- imms[[k]]
      saveDir4 <- paste0(saveDir3, imm, "/")
      if (!dir.exists(saveDir4)) {
        dir.create(saveDir4)
      }
      nums_w_reps_single <- nums_w_reps[[k]]
      shannon2VAll <- NULL
      for (m in 1:length(nums)) {
        num <- nums[m]
        if (num %in% nums_w_reps_single) {
          reps <- 0:4
        } else {
          reps <- 0
        }
        for (r in reps) {
          file <- paste0(readDir, category, "/", s, "/", o, "/", state, "/", category, "-", num, "-r", r, ".RData")
          load(file)
          df2 <- df %>% mutate(IRS = case_when(
            time <= preIRS*T_YEAR ~ "preIRS",
            time > preIRS*T_YEAR ~ paste0("I-", num - min(nums) + 1)
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
      summaryTable2 <- summaryTable1 %>% filter(!(time == preIRS * T_YEAR & IRS %in% c(paste0("I-", 1:16))))
      CRs <- summaryTable2 %>% select(IRS, CR) %>% group_by(IRS) %>% summarise(meanCR = mean(CR))
      SD <- shannon2VAll %>% left_join(CRs, by = "IRS")
      SD2 <- SD %>% group_by(time, IRS, meanCR) %>% 
        summarise(mD = mean(shannon), maxD = max(shannon), minD = min(shannon))
      SD2 <- SD2 %>% mutate(timePlot = time/T_YEAR -preIRS)
      SD2$IRS <- factor(SD2$IRS, levels = c("preIRS", paste0("I-", 1:16)))
      
      p1 <- ggplot(SD2, aes(timePlot, mD, col = IRS, fill = IRS)) + 
        geom_line(size = pointSize) +
        geom_ribbon(aes(y = mD, ymin = minD, ymax = maxD), alpha = .4, linetype = 0) +
        xlab("Years Since IRS Starts \n(Sustained IRS)") + ylab("Shannon Diversity") +
        theme_bw() + theme(
          plot.title = element_text(color="black", size=sizeV, hjust = 0.5),
          axis.title.x = element_text(color="black", size=sizeV),
          axis.title.y = element_text(color="black", size=sizeV),
          axis.text.x = element_text(color="black", size=sizeV, angle=0),
          axis.text.y = element_text(color="black", size=sizeV, angle=0),
          legend.text = element_text(color="black", size=sizeV, angle=0),
          legend.title = element_text(color="black", size=sizeV, angle=0),
          plot.margin = unit(c(5,10,5,5), "points"),
          strip.text = element_text(color="black", size=sizeV, angle=0)) + 
        coord_cartesian(ylim = c(min(SD2$minD) - 1, max(SD2$maxD) + 0.5)) +
        scale_color_manual(values = Turbo(out.colors = 17)[1:(length(nums) + 1)]) +
        scale_fill_manual(values = Turbo(out.colors = 17)[1:(length(nums) + 1)]) +
        scale_x_continuous(breaks = seq(-2, 10, 2)) +
        scale_y_continuous(labels = scaleFUN) +
        guides(fill = "none", col = "none") 
      p1
      
      if (s == "seasonal") {
        SD2 <- SD %>% filter(time %in% c((preIRS + forwardYs)*T_YEAR + 300))
      } else {
        SD2 <- SD %>% filter(time %in% c((preIRS + forwardYs + 1)*T_YEAR))
      }
      p2 <- ggplot(SD2, aes(x = meanCR, y = shannon, col = IRS)) + 
        stat_summary(fun.data = statsSummaryf, geom="boxplot", position=position_dodge(1)) +
        theme_bw() +
        xlab("Effective Contact Rate \n(per Host, per Year)") + ylab("Shannon Diveristy") +
        theme(
          plot.title = element_text(color="black", size=sizeV, hjust = 0.5, vjust = 1),
          axis.title.x = element_text(color="black", size=sizeV),
          axis.title.y = element_text(color="white", size=sizeV),
          axis.text.x = element_text(color="black", size=sizeV, angle=0),
          axis.text.y = element_text(color="white", size=sizeV, angle=0),
          legend.text = element_text(color="black", size=sizeV, angle=0),
          legend.title = element_blank(),
          strip.text = element_text(color="black", size=sizeV, angle=0),
          plot.margin = unit(c(5,22,5,5), "points")) + 
        scale_color_manual(values = Turbo(out.colors = 17)[2:(length(nums) + 1)]) +
        scale_fill_manual(values = Turbo(out.colors = 17)[2:(length(nums) + 1)]) +
        scale_linetype(guide = "none") + 
        coord_cartesian(ylim = c(min(SD2$shannon) - 1, max(SD2$shannon) + 0.5)) +
        guides(col = "none", fill = "none") +
        scale_y_continuous(labels = scaleFUN)
      print(p2)
      
      p1_rmL <- p1 + guides(fill = "none", col = "none") + rremove("legend") 
      p2_rmL <- p2 + guides(fill = "none", col = "none") + rremove("legend") 

      ggsave(paste0(saveDir4, "num-", min(nums), "-", max(nums), "-time.pdf"), p1_rmL, width = 6, height = 6)
      ggsave(paste0(saveDir4, "num-", min(nums), "-", max(nums), "-cr.pdf"), p2_rmL, width = 6, height = 6)
    }
  }
}
