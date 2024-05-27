rm(list = ls())
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})
run <- 4
source(paste0("/home/qizhan/others/PhD/projects/intervention/writings/simulation", run, "/plots/colorFunc.R"))
saveDir0 <- paste0("/project2/pascualmm/QZ/PhD/projects/intervention/natComRevision/plots/figures/main/Fig3/extinctionProb/")
if (!dir.exists(saveDir0)) {
  dir.create(saveDir0)
}
readDir <- "/project2/pascualmm/QZ/PhD/projects/intervention/natComRevision/files/figures/main/Fig2/prev/" 

seasonality <- c("seasonal")
openness <- c("regionally-open")
nums <- 208:209
prefix <- "sim"
preIRS <- 150
T_YEAR <- 360
IRSType <- "10yIRS"
reps <- 0:4
ts <- seq(preIRS + 4, preIRS + 10, 2)
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
    
    dfAll <- NULL
    for (k in 1:length(nums)) {
      num <- nums[k]
      for (m in 1:length(ts)) {
        t <- ts[m]
        nExint <- 0
        for (n in 1:length(reps)) {
          r <- reps[n]
          summaryTable <- summaryTableCombined %>% filter(num == nums[k], time <= t*T_YEAR, time > preIRS*T_YEAR, rep == r)
          if (any(summaryTable$Prevalence == 0)) {
            minTime <- min(summaryTable %>% filter(Prevalence == 0) %>% .$time)
            stopifnot(all(summaryTable %>% filter(time >= minTime) %>% .$Prevalence == 0))
            nExint <- nExint + 1
          } 
        }
        df <- data.frame("p_extinc" = nExint/length(reps), 
                         "num" = num, "timeCheck" = t - preIRS, 
                         "IRS" = paste0("I-", num - 200))
        dfAll <- rbind(dfAll, df)
      }
    }
    
    dfAll$timeCheck <- factor(dfAll$timeCheck, levels = sort(unique(dfAll$timeCheck)))
    dfAll$IRS <- factor(dfAll$IRS, levels = unique(dfAll$IRS))
    sizeV <- 31
    p1 <- ggplot(dfAll, aes(IRS, p_extinc, alpha = timeCheck, col = IRS, fill = IRS))+ 
      geom_bar(stat = "identity", position = "dodge")+
      xlab("IRS Runs") + ylab("Probability of Extinction") +
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
      scale_color_manual(values = Turbo(out.colors = 17)[((min(nums)-200) + 1):((max(nums)-200) + 1)]) +
      scale_fill_manual(values = Turbo(out.colors = 17)[((min(nums)-200) + 1):((max(nums)-200) + 1)]) + 
      guides(fill = "none", col = "none") +
      scale_y_continuous(breaks = seq(0, 1, by = 0.2),
                         labels = seq(0, 1, by = 0.2)) + 
      scale_alpha_manual(name = "Years into \nSustained \nIRS", values = seq(0.4, 1, 0.2)) + 
      theme(legend.position = "right") + 
      coord_cartesian(ylim = c(0, 0.8)) + guides(alpha = guide_legend(ncol = 1))
    print(p1)
    ggsave(paste0(saveDir2, min(nums), "-", max(nums), "_extinctionProbability.pdf"), p1, width = 7.5, height = 6)
  }
}
