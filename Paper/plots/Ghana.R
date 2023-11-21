rm(list = ls())
suppressPackageStartupMessages({
  library(RSQLite)
  library(reshape2)
  library(vegan)
  library(tidyverse)
  library(igraph)
  library(cowplot)
  library(scales)
  #library(adegenet)
  library(factoextra)
  library(ggplot2)
  library(ggpubr)
})
sizeV <- 31
colors <- c("black", scales::hue_pal()(5))[c(3,6)]
saveDir <- "/project2/pascualmm/QZ/PhD/projects/intervention/writings/simulation4/plots/Ghana/MOI/"
pAll <- list()
# MOI
cols <- c("dark orange", "dark blue", "dark blue")
nums <- c(1,4,5)
labels <- c("Pre-IRS", "IRS", "Right Post-IRS")
for (i in 1:length(nums)) {
  num <- nums[i]
  col <- cols[i]
  label <- labels[i]
  load(paste0("~/others/projects/summer2022/round2/MOI/scripts/surveysMOI/survey_", num))
  pts2x<-ggplot(dfAll, 
                aes(MOI, col=as.factor(survey), fill =as.factor(survey)))+
    geom_histogram(bins = 40)+
    ggtitle(label) +
    theme_classic() + theme(
      plot.title = element_text(color="white", size=sizeV, hjust = 0.5),
      axis.title.x = element_text(color="black", size=sizeV),
      axis.title.y = element_text(color="black", size=sizeV),
      axis.text.x = element_text(color="black", 
                                 size=sizeV, angle=0),
      axis.text.y = element_text(color="black", 
                                 size=sizeV, angle=0),
      legend.text = element_text(face="bold", color="black", 
                                 size=10, angle=0),
      legend.title = element_text(face="bold", color="black", 
                                  size=11, angle=0),
      strip.text = element_blank()) +
    coord_cartesian(xlim = c(1,21)) + 
    scale_x_continuous(breaks = seq(1,21,4)) +
    scale_fill_manual(name = "",values = col) +
    ylab("Count") +
    xlab("MOI") +
    scale_color_manual(name = "", values = col) +
    guides("fill" = "none", "color" = "none")
  print(pts2x)
  pAll[[i]] <- pts2x
}
ggsave(paste0(saveDir, "Ghana-MOI-survey_1.pdf"), pAll[[1]] + rremove("xlab") + rremove("ylab"), width = 6, height = 6)
ggsave(paste0(saveDir, "Ghana-MOI-survey_4.pdf"), pAll[[2]] + rremove("xlab") + rremove("ylab"), width = 6, height = 6)
ggsave(paste0(saveDir, "Ghana-MOI-survey_5.pdf"), pAll[[3]] + rremove("xlab") + rremove("ylab"), width = 6, height = 6)




# PTS by age
rm(list = ls())
saveDir <- "/project2/pascualmm/QZ/PhD/projects/intervention/writings/simulation4/plots/Ghana/PTSAgeDiff/"
pAll <- list()
load(file = "/project2/pascualmm/QZ/papersOfficial/intervention/analysis/files/PTS/withME/Ghana/PTS-BC-byAge-separateMOI-1-to-20-combineYoungAgeGroups-add2015")
dfAll <- PTSdfAll %>% filter(state %in% c("pre IRS", "2y into IRS", "right post IRS (2015)")) %>% 
  mutate(state = case_when(state == "2y into IRS" ~ "IRS",
                           state == "pre IRS" ~ "Pre-IRS",
                           state == "right post IRS (2015)" ~ "Right Post-IRS"))
dfAll$state <- factor(dfAll$state, levels = c("Pre-IRS", "IRS", "Right Post-IRS"))
dfAll$age <- as.character(dfAll$age)
colnames(dfAll)[2] <- "Age"
ageSel <- c("1-10", ">20")
dfAll <- dfAll %>% mutate("MOI_label" = paste0("MOI = ", MOI))
sizeV <- 31
nums <- c(1,4,5)
labels <- c("Pre-IRS", "IRS", "Right Post-IRS")
for (i in 1:length(nums)) {
  num <- nums[i]
  label <- labels[i]
  p <- ggplot(dfAll %>% filter(Age %in% ageSel, survey %in% c(paste0("survey_", num))), aes(x = PTS, col = Age, fill = Age)) +
    geom_density(aes(group=Age, alpha = 0.1)) +
    ggtitle(label) +
    xlab("PTS") + ylab("Density") +
    theme_classic() + theme(
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
    scale_fill_manual(values=c("magenta", "dark green")) +
    scale_color_manual(values=c("magenta", "dark green")) +
    coord_cartesian(xlim = c(0,0.4)) + scale_alpha(guide = 'none') + 
    theme(legend.position = c(0.25, 0.65)) 
  print(p)
  if (num %in% c(4,5)) {
    p <- p + guides("fill" = "none", "color" = "none") 
  }
  pAll[[i]] <- p
}
ggsave(paste0(saveDir, "Ghana-PTS-survey_1.pdf"), pAll[[1]] + rremove("xlab") + rremove("ylab"), width = 6, height = 6)
ggsave(paste0(saveDir, "Ghana-PTS-survey_4.pdf"), pAll[[2]] + rremove("xlab") + rremove("ylab"), width = 6, height = 6)
ggsave(paste0(saveDir, "Ghana-PTS-survey_5.pdf"), pAll[[3]] + rremove("xlab") + rremove("ylab"), width = 6, height = 6)



# PTS by age zoom in
rm(list = ls())
saveDir <- "/project2/pascualmm/QZ/PhD/projects/intervention/writings/simulation4/plots/Ghana/PTSAgeDiffZoomin/"
pAll <- list()
load(file = "/project2/pascualmm/QZ/papersOfficial/intervention/analysis/files/PTS/withME/Ghana/PTS-BC-byAge-separateMOI-1-to-20-combineYoungAgeGroups-add2015")
dfAll <- PTSdfAll %>% filter(state %in% c("pre IRS", "2y into IRS", "right post IRS (2015)")) %>% 
  mutate(state = case_when(state == "2y into IRS" ~ "IRS",
                           state == "pre IRS" ~ "Pre-IRS",
                           state == "right post IRS (2015)" ~ "Right Post-IRS"))
dfAll$state <- factor(dfAll$state, levels = c("Pre-IRS", "IRS", "Right Post-IRS"))
dfAll$age <- as.character(dfAll$age)
colnames(dfAll)[2] <- "Age"
ageSel <- c("1-10", ">20")
dfAll <- dfAll %>% mutate("MOI_label" = paste0("MOI = ", MOI))
sizeV <- 31
nums <- c(1,4,5)
labels <- c("Pre-IRS", "IRS", "Right Post-IRS")
for (i in 1:length(nums)) {
  num <- nums[i]
  label <- labels[i]
  p <- ggplot(dfAll %>% filter(Age %in% ageSel, survey %in% c(paste0("survey_", num))), aes(x = PTS, col = Age, fill = Age)) +
    geom_density(aes(group=Age, alpha = 0.1)) +
    ggtitle(label) +
    xlab("PTS") + ylab("Density") +
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
                                size=sizeV, angle=0)) +
    scale_fill_manual(values=c("magenta", "dark green")) +
    scale_color_manual(values=c("magenta", "dark green")) +
    coord_cartesian(xlim = c(0,0.075)) + scale_alpha(guide = 'none') + 
    theme(legend.position = c(0.65, 0.65)) +
    guides(fill = "none", col = "none")
  print(p)
  pAll[[i]] <- p
}
ggsave(paste0(saveDir, "Ghana-PTS-survey_1.pdf"), pAll[[1]] + rremove("xlab") + rremove("ylab"), width = 6, height = 6)
ggsave(paste0(saveDir, "Ghana-PTS-survey_4.pdf"), pAll[[2]] + rremove("xlab") + rremove("ylab"), width = 6, height = 6)
ggsave(paste0(saveDir, "Ghana-PTS-survey_5.pdf"), pAll[[3]] + rremove("xlab") + rremove("ylab"), width = 6, height = 6)

