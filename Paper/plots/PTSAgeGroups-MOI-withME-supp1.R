rm(list = ls())
suppressPackageStartupMessages({
  library(RSQLite)
  library(reshape2)
  library(vegan)
  library(tidyverse)
  library(igraph)
  library(cowplot)
  library(scales)
  library(factoextra)
  library(ggplot2)
  library(ggpubr)
  library(gridExtra)
  library(patchwork)
  library(grid)
})
run <- 4
source(paste0("/home/qizhan/others/PhD/projects/intervention/writings/simulation", run, "/plots/colorFunc.R"))
readDir <- paste0("/project2/pascualmm/QZ/PhD/projects/intervention/files", run, "/actualRuns/PTSAgeDiff1/")
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

type <- "dir"
numsPlotList <- list(
  "seasonal closed" = list(101:107,115:121),
  "seasonal semi-open" = list(201:207, 216:222,  
                              301:307, 316:322),
  "seasonal regionally-open" = list(201:207,219:227,
                                    237:243,255:263,
                                    301:307,319:327,
                                    337:343,355:363),
  "non-seasonal closed" = list(101:107, 115:121),
  "non-seasonal semi-open" = list(201:207,
                                  216:222,301:307,316:322),
  "non-seasonal regionally-open" = list(201:207,219:227,
                                        237:243,255:263,
                                        301:307,319:327,
                                        337:343,355:363)
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
colors <- c("black", scales::hue_pal()(5))[c(3,6)]
N <- 10000
pAll1 <- list()
pAll2 <- list()
sizeV <- 1

saveDir0 <- paste0("/project2/pascualmm/QZ/PhD/projects/intervention/writings/simulation", run, "/plots/")
dir.create(saveDir0)
saveDir1 <- paste0(saveDir0, "PTSAgeGroupsAndMOI1/PTSZoomin/")
dir.create(saveDir1)
state <- "withME"
saveDir2 <- paste0(saveDir1, IRSType, "/")
dir.create(saveDir2)

MOI_chosen = 1:20
for (i in 1:length(seasonality)) {
  s <- seasonality[i]
  saveDir3 <- paste0(saveDir2, s, "/")
  dir.create(saveDir3)
  
  for (j in 1:length(openness)) {
    o <- openness[j]
    saveDir4 <- paste0(saveDir3, o, "/")
    dir.create(saveDir4)
    
    preIRS <- preIRSList[[paste0(s, " ", o)]]
    numsPlot <- numsPlotList[[paste0(s, " ", o)]]
    labels <- labelsList[[paste0(s, " ", o)]]
    imms <- immsList[[paste0(s, " ", o)]]
    
    for (k in 1:length(numsPlot)) {
      nums <- numsPlot[[k]]
      label <- labels[[k]]
      imm <- imms[[k]]
      
      if (k %in% c(2,4)) {
        xlimmax <- 0.05
      } else {
        xlimmax <- 0.1
      }
      
      saveDir5 <- paste0(saveDir4, imm, "/")
      dir.create(saveDir5)
      
      saveDir6 <- paste0(saveDir5, state, "/")
      dir.create(saveDir6)
      
      if (s == "seasonal") {
        layers <- c(preIRS - 1, preIRS + 1)*T_YEAR + 300
      } else {
        layers <- c(preIRS, preIRS + 2)*T_YEAR
      }
      for (m in 1:length(nums)) {
        num <- nums[m]
        file <- paste0(readDir, s, "/", o, "/", type, "/", state, "/", "PTSAcrossAgeGroups_", num, "_r0.RData")
        load(file)
        if (m == 1) {
          PTSAllSub <- PTSAll %>% filter(num %in% nums[m], time == layers[1], MOI %in% MOI_chosen)
          pts2x<-ggplot(PTSAllSub, aes(x=PTS, fill = ageGroup, col = ageGroup))+
            geom_density(alpha = 0.4, aes(group = ageGroup)) +
            ggtitle(paste0("Pre-IRS")) +
            theme_cowplot() +
            xlab("PTS") + ylab("Density") +
            theme_classic() + theme(
              plot.title = element_text(color="black", size=sizeV, hjust = 0.5),
              axis.title.x = element_text(color="white", size=sizeV),
              axis.title.y = element_text(color="white", size=sizeV),
              axis.text.x = element_text(color="white", 
                                         size=sizeV, angle=0),
              axis.text.y = element_text(color="white", 
                                         size=sizeV, angle=0),
              legend.text = element_text(color="black", 
                                         size=sizeV, angle=0),
              legend.title = element_text(color="black", 
                                          size=sizeV, angle=0),
              strip.text = element_text(color="black", 
                                        size=sizeV, angle=0),
              legend.position=c(0.55,0.6)) + 
            coord_cartesian(xlim = c(0, xlimmax)) +
            scale_fill_manual(name = "Age",values = c("magenta", "dark green")) +
            scale_x_continuous(expand = expansion(mult = 0.1)) + 
            scale_color_manual(name = "Age",values = c("magenta", "dark green")) +
            guides(fill = "none", col = "none")
          ggsave(paste0(saveDir6, "PTSAgeGroups-", num, "-Pre-IRS.pdf"), pts2x, width = 5, height = 5)
        }
        PTSAllSub <- PTSAll %>% filter(num %in% nums[m], time == layers[2], MOI %in% MOI_chosen)
        pts2x<-ggplot(PTSAllSub, aes(x=PTS, fill = ageGroup, col = ageGroup))+
          geom_density(alpha = 0.4, aes(group = ageGroup)) +
          ggtitle(paste0("I-", num - min(nums) + 1)) +
          theme_cowplot() + xlab("PTS") + ylab("Density") +
          theme_classic() + theme(
            plot.title = element_text(color="white", size=sizeV, hjust = 0.5),
            axis.title.x = element_text(color="white", size=sizeV),
            axis.title.y = element_text(color="white", size=sizeV),
            axis.text.x = element_text(color="white", 
                                       size=sizeV, angle=0),
            axis.text.y = element_text(color="white", 
                                       size=sizeV, angle=0),
            legend.text = element_text(color="black", 
                                       size=sizeV, angle=0),
            legend.title = element_text(color="black", 
                                        size=sizeV, angle=0),
            strip.text = element_text(color="black", 
                                      size=sizeV, angle=0),
            legend.position=c(0.55,0.6)) + 
          coord_cartesian(xlim = c(0, xlimmax)) +
          scale_fill_manual(name = "Age",values =  c("magenta", "dark green")) +
          scale_x_continuous(expand = expansion(mult = 0.1)) + 
          scale_color_manual(name = "Age",values = c("magenta", "dark green")) +
          guides(fill = "none", col = "none")
        pts2x
        ggsave(paste0(saveDir6, "PTSAgeGroups-", num, "-IRS.pdf"), pts2x, width = 5, height = 5)
        
      }
    }
  }
}

