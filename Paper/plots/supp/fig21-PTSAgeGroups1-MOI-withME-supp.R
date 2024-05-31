rm(list = ls())
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})
run <- 4
source(paste0("/home/qizhan/others/PhD/projects/intervention/writings/simulation", run, "/plots/colorFunc.R"))
readDir <- paste0("/project2/pascualmm/QZ/PhD/projects/intervention/files", run, "/actualRuns/PTSAgeDiff1/")
openness <- "regionally-open"
seasonality <- "non-seasonal"
preIRSList <- list(
  "non-seasonal regionally-open" = 150
)

IRSDur <- 10
IRSType <- paste0(IRSDur, "yIRS")
postIRS <- 10 - IRSDur

type <- "dir"
numsPlotList <- list(
  "non-seasonal regionally-open" = list(255:263)
)
labelsList <- list(
  "non-seasonal regionally-open" = list("Non-seasonal \nRegionally-open \nHigh Migration \nLarge Pool")
)
immsList <- list(
  "non-seasonal regionally-open" = list("Three-year")
)

T_YEAR <- 360
colors <- c("black", scales::hue_pal()(5))[c(3,6)]
sizeV <- 1

saveDir0 <- "/project2/pascualmm/QZ/PhD/projects/intervention/natComRevision/plots/figures/supp/fig21/"
if (!dir.exists(saveDir0)) {
  dir.create(saveDir0)
}
saveDir1 <- paste0(saveDir0, "PTSAgeGroups1Zoomin/")
if (!dir.exists(saveDir1)) {
  dir.create(saveDir1)
}
state <- "withME"
saveDir2 <- paste0(saveDir1, state, "/")
if (!dir.exists(saveDir2)) {
  dir.create(saveDir2)
}

MOI_chosen = 1:20
for (i in 1:length(seasonality)) {
  s <- seasonality[i]
  saveDir3 <- paste0(saveDir2, s, "/")
  if (!dir.exists(saveDir3)) {
    dir.create(saveDir3)
  }
  
  for (j in 1:length(openness)) {
    o <- openness[j]
    saveDir4 <- paste0(saveDir3, o, "/")
    if (!dir.exists(saveDir4)) {
      dir.create(saveDir4)
    }
    
    preIRS <- preIRSList[[paste0(s, " ", o)]]
    numsPlot <- numsPlotList[[paste0(s, " ", o)]]
    labels <- labelsList[[paste0(s, " ", o)]]
    imms <- immsList[[paste0(s, " ", o)]]
    
    for (k in 1:length(numsPlot)) {
      nums <- numsPlot[[k]]
      label <- labels[[k]]
      imm <- imms[[k]]
      
      xlimmax <- 0.05
      
      saveDir5 <- paste0(saveDir4, imm, "/")
      if (!dir.exists(saveDir5)) {
        dir.create(saveDir5)
      }
      
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
            xlab("PTS") + ylab("Density") +
            theme_classic() + theme(
              plot.title = element_text(color="white", size=sizeV, hjust = 0.5),
              axis.title.x = element_text(color="white", size=sizeV),
              axis.title.y = element_text(color="white", size=sizeV),
              axis.text.x = element_text(color="white", size=sizeV, angle=0),
              axis.text.y = element_text(color="white", size=sizeV, angle=0),
              legend.text = element_text(color="black", size=sizeV, angle=0),
              legend.title = element_text(color="black", size=sizeV, angle=0),
              strip.text = element_text(color="black", size=sizeV, angle=0),
              legend.position=c(0.55,0.6)) + 
            coord_cartesian(xlim = c(0, xlimmax)) +
            scale_fill_manual(name = "Age", values = c("magenta", "dark green")) +
            scale_x_continuous(expand = expansion(mult = 0.1)) + 
            scale_color_manual(name = "Age", values = c("magenta", "dark green")) +
            guides(fill = "none", col = "none")
          ggsave(paste0(saveDir5, "PTSAgeGroups-", num, "-Pre-IRS.pdf"), pts2x, width = 5, height = 5)
        }
        PTSAllSub <- PTSAll %>% filter(num %in% nums[m], time == layers[2], MOI %in% MOI_chosen)
        pts2x<-ggplot(PTSAllSub, aes(x = PTS, fill = ageGroup, col = ageGroup))+
          geom_density(alpha = 0.4, aes(group = ageGroup)) +
          ggtitle(paste0("I-", num - min(nums) + 1)) +
          xlab("PTS") + ylab("Density") +
          theme_classic() + theme(
            plot.title = element_text(color="white", size=sizeV, hjust = 0.5),
            axis.title.x = element_text(color="white", size=sizeV),
            axis.title.y = element_text(color="white", size=sizeV),
            axis.text.x = element_text(color="white", size=sizeV, angle=0),
            axis.text.y = element_text(color="white", size=sizeV, angle=0),
            legend.text = element_text(color="black", size=sizeV, angle=0),
            legend.title = element_text(color="black", size=sizeV, angle=0),
            strip.text = element_text(color="black", size=sizeV, angle=0),
            legend.position=c(0.55,0.6)) + 
          coord_cartesian(xlim = c(0, xlimmax)) +
          scale_fill_manual(name = "Age", values = c("magenta", "dark green")) +
          scale_x_continuous(expand = expansion(mult = 0.1)) + 
          scale_color_manual(name = "Age", values = c("magenta", "dark green")) +
          guides(fill = "none", col = "none")
        pts2x
        ggsave(paste0(saveDir5, "PTSAgeGroups-", num, "-IRS.pdf"), pts2x, width = 5, height = 5)
        
      }
    }
  }
}

