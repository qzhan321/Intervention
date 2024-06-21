rm(list = ls())
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(cowplot)
})
run <- 4
saveDir0 <- paste0("/project2/pascualmm/QZ/PhD/projects/intervention/natComRevision/plots/figures/supp/fig6/PTSCompleteOverlap/")
if (!dir.exists(saveDir0)) {
  dir.create(saveDir0)
}
readDir <- paste0("/project2/pascualmm/QZ/PhD/projects/intervention/files", run, "/actualRuns/PTSCompleteOverlap/")

numsPlotList <- list(
  "seasonal regionally-open" = list(255:265)
)
seasonality <- c("seasonal")
openness <- c("regionally-open")
preIRSList <- list(
  "seasonal regionally-open" = 150
)
labelsList <- list(
  "seasonal regionally-open" = list("Seasonal \nRegionally-open \nHigh Migration \nLarge Pool")
)
immsList <- list(
  "seasonal regionally-open" = list("Three-year")
)
T_YEAR <- 360
IRSDur <- 10
sizeV <- 25
pAll <- list()
type <- "dir"
state <- "true"
nums_w_reps <- NULL
colors <- c("cyan", "dark red")
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
    if (s == "seasonal") {
      layers <- c(preIRS-1, preIRS+1)*T_YEAR+300
    } else {
      layers <- c(preIRS, preIRS+2)*T_YEAR
    }
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
      
      saveDir4 <- paste0(saveDir3, type, "/")
      if (!dir.exists(saveDir4)) {
        dir.create(saveDir4)
      }
      
      saveDir5 <- paste0(saveDir4, state, "/")
      if (!dir.exists(saveDir5)) {
        dir.create(saveDir5)
      }
      
      for (m in 1:length(nums)) {
        num <- nums[m]
        
        if (num %in% nums_w_reps) {
          reps <- 0:2
        } else {
          reps <- 0
        }
        for (r in reps) {
          file <- paste0(readDir, s, "/", o, "/", type, "/", state, "/", "PTS_selectedTimes_", num, "_r", r, ".RData")
          load(file)
          temp1 <- dfAll %>% filter(time %in% layers) %>% mutate(state = state, type = type)
          PTSGroups <- seq(-0.025,1,0.025)
          breaksL <- rep("", length(PTSGroups)-1)
          index <- c(2,22,length(PTSGroups))
          breaksL[index-1] <- PTSGroups[index]
          temp2 <- temp1 %>% mutate(PTSGroup = cut(PTS, breaks = PTSGroups))
          temp3 <- temp2 %>% group_by(num, time, state, type, PTSGroup) %>% summarise(n=n())
          temp4 <- temp3 %>% group_by(num, time, state, type) %>% summarise(nTotal = sum(n))
          
          temp5 <- temp3 %>% left_join(temp4, by=c("num", "time", "state", "type")) %>% mutate(p = n/nTotal)
          temp5 %>% group_by(num, time, state, type) %>% summarise(sum(p))
          temp6 <- temp5 %>% mutate(IRS = paste0("I-", num - min(nums) + 1))
          temp6$PTSGroup <- factor(temp6$PTSGroup, levels = levels(cut(PTSGroups[-1], breaks = PTSGroups)))
          
          p1<-ggplot(temp6, aes(x=PTSGroup, y=p, col=as.factor(time), fill = as.factor(time)))+
            geom_bar(stat = "identity", position = "dodge")+ 
            ggtitle(paste0("I-", nums[m] - min(nums) + 1)) +
            theme_cowplot()+
            theme(plot.title = element_text(color="black", size=sizeV, hjust = 0.5),
              axis.title.x = element_text(color="black", size=sizeV),
              axis.title.y = element_text(color="black", size=sizeV),
              axis.text.x = element_text(color="black", size=sizeV, angle=0),
              axis.text.y = element_text(color="black", size=sizeV, angle=0),
              legend.text = element_text(color="black", size=sizeV, angle=0),
              legend.title = element_text(color="black", size=sizeV, angle=0),
              strip.text = element_text(color="black", size=sizeV, angle=0)) + 
            scale_fill_manual(name = "",values = colors, labels = c("Pre-IRS", "IRS")) +
            ylab("Frequency")+xlab("PTS")+
            scale_color_manual(name = "", labels = c("Pre-IRS", "IRS"), values = colors) +
            scale_x_discrete(labels = breaksL, expand = c(0.05,0)) + 
            theme(legend.position = c(0.5, 0.75)) +
            coord_cartesian(ylim = c(0, 0.8))
          print(p1)
          if (m > 1) {
            p1 <- p1 + theme(legend.position = "none")
          }
          ggsave(paste0(saveDir5, "num-", num, "-r", r, ".pdf"), p1, width = 7.00/1.24, height = 6.20/1.24)
        }
      }
    }
  }
}
