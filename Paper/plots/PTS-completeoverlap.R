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
source(paste0("/home/qizhan/others/PhD/projects/intervention/writings/simulation", run, "/plots/colorFunc.R"))
saveDir0 <- paste0("/project2/pascualmm/QZ/PhD/projects/intervention/writings/simulation", run, "/plots/PTSCompleteoverlap/")
dir.create(saveDir0)
readDir <- paste0("/project2/pascualmm/QZ/PhD/projects/intervention/files", run, "/actualRuns/PTSCompleteOverlap/")

numsPlotList <- list(
  "seasonal closed" = list(101:109),
  "seasonal semi-open" = list(201:208,
                              301:308),
  "seasonal regionally-open" = list(201:208,219:230,237:244,255:265),
  "non-seasonal closed" = list(101:110),
  "non-seasonal semi-open" = list(201:211,
                                  301:311),
  "non-seasonal regionally-open" = list(201:209,219:229,
                                        237:244,255:264)
)
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
labelsList <- list(
  "seasonal closed" = list("Seasonal \nClosed \n \n", "Seasonal \nClosed \n \n"),
  "seasonal semi-open" = list("Seasonal \nSemi-open \nBaseline Migration \n", "Seasonal \nSemi-open \nBaseline Migration \n", 
                              "Seasonal \nSemi-open \nHigh Migration \n", "Seasonal \nSemi-open \nHigh Migration \n"),
  "seasonal regionally-open" = list("Seasonal \nRegionally-open \nBaseline Migration \nMedium Pool",
                                    "Seasonal \nRegionally-open \nBaseline Migration \nLarge Pool",
                                    "Seasonal \nRegionally-open \nHigh Migration \nMedium Pool",
                                    "Seasonal \nRegionally-open \nHigh Migration \nLarge Pool",
                                    "Seasonal \nRegionally-open \nBaseline Migration \nMedium Pool",
                                    "Seasonal \nRegionally-open \nBaseline Migration \nLarge Pool",
                                    "Seasonal \nRegionally-open \nHigh Migration \nMedium Pool",
                                    "Seasonal \nRegionally-open \nHigh Migration \nLarge Pool"),
  "non-seasonal closed" = list("Non-seasonal \nClosed \n \n", "Non-seasonal \nClosed \n \n"),
  "non-seasonal semi-open" = list("Non-seasonal \nSemi-open \nBaseline Migration \n", "Non-seasonal \nSemi-open \nBaseline Migration \n",  
                                  "Non-seasonal \nSemi-open \nHigh Migration \n", "Non-seasonal \nSemi-open \nHigh Migration \n"),
  "non-seasonal regionally-open" = list("Non-seasonal \nRegionally-open \nBaseline Migration \nMedium Pool",
                                        "Non-seasonal \nRegionally-open \nBaseline Migration \nLarge Pool",
                                        "Non-seasonal \nRegionally-open \nHigh Migration \nMedium Pool",
                                        "Non-seasonal \nRegionally-open \nHigh Migration \nLarge Pool",
                                        "Non-seasonal \nRegionally-open \nBaseline Migration \nMedium Pool",
                                        "Non-seasonal \nRegionally-open \nBaseline Migration \nLarge Pool",
                                        "Non-seasonal \nRegionally-open \nHigh Migration \nMedium Pool",
                                        "Non-seasonal \nRegionally-open \nHigh Migration \nLarge Pool")
  
)
immsList <- list(
  "seasonal closed" = list("Three-year", "Four-year"),
  "seasonal semi-open" = list("Three-year", "Four-year",
                              "Three-year", "Four-year"),
  "seasonal regionally-open" = list("Three-year", "Three-year", "Three-year", "Three-year",
                                    "Four-year", "Four-year", "Four-year", "Four-year"),
  "non-seasonal closed" = list("Three-year", "Four-year"),
  "non-seasonal semi-open" = list("Three-year", "Four-year",
                                  "Three-year", "Four-year"),
  "non-seasonal regionally-open" = list("Three-year", "Three-year", "Three-year", "Three-year",
                                        "Four-year", "Four-year", "Four-year", "Four-year")
)
T_YEAR <- 360
IRSDur <- 10
sizeV <- 31
pAll <- list()
type <- "dir"
state <- "true"
nums_w_reps <- NULL
colors <- c("black", scales::hue_pal()(5))[c(3,6)]
colors <- c("cyan", "dark red")
for (i in 1:length(seasonality)) {
  s <- seasonality[i]
  saveDir1 <- paste0(saveDir0, s, "/")
  dir.create(saveDir1)
  
  for (j in 1:length(openness)) {
    o <- openness[j]
    saveDir2 <- paste0(saveDir1, o, "/")
    dir.create(saveDir2)
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
      dir.create(saveDir3)
      
      saveDir4 <- paste0(saveDir3, type, "/")
      dir.create(saveDir4)
      
      saveDir5 <- paste0(saveDir4, state, "/")
      dir.create(saveDir5)
      
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
          temp6$PTSGroup <- factor(temp6$PTSGroup, levels = unique(temp6$PTSGroup))
          
          p1<-ggplot(temp6, aes(x=PTSGroup, y=p, col=as.factor(time), fill = as.factor(time)))+
            geom_bar(stat = "identity", position = "dodge")+ 
            ggtitle(paste0("I-", nums[m] - min(nums) + 1)) +
            theme_cowplot()+
            theme(
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
            scale_fill_manual(name = "",values = colors, labels = c("Pre-IRS", "IRS")) +
            ylab("Frequency")+xlab("PTS")+
            scale_color_manual(name = "", labels = c("Pre-IRS", "IRS"), values = colors) +
            scale_x_discrete(labels = breaksL, expand = c(0.05,0)) + theme(legend.position = c(0.5, 0.75))
          print(p1)
          if (m > 1) {
            p1 <- p1 + theme(legend.position = "none")
          }
          ggsave(paste0(saveDir5, "num-", num, "-r", r, ".pdf"), p1, width = 7.50, height = 6.20)
        }
      }
    }
  }
}
