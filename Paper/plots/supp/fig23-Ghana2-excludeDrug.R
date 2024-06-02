rm(list = ls())
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyverse)
  library(ggplot2)
  library(ggpubr)
})
# MOI
sizeV <- 28
colors <- c("black", scales::hue_pal()(5))[c(3,6)]
saveDir <- "/project2/pascualmm/QZ/PhD/projects/intervention/natComRevision/plots/figures/supp/fig23/MOI2-excludeDrug/"
if (!dir.exists(saveDir)) {
  dir.create(saveDir)
}
readDir <- "/project2/pascualmm/QZ/PhD/projects/intervention/natComRevision/files/figures/main/Fig5/MOI/MOIEst/"
pAll <- list()
cols <- c("dark orange", "dark blue", "dark blue")
nums <- c(1,4,5)
prefix <- "survey"
labels <- c("Pre-IRS", "IRS", "Post-IRS")

readDirEpi <- "/home/qizhan/others/PhD/projects/intervention/natComRevision/utils/"
epi <- read.csv(paste0(readDirEpi, "Ghana_Survey_Merged_Epi_MOI_S1_S7_070721_UChicago_080822.csv"), header = T, row.names = 1)
epi$SeqID <- str_replace(epi$SeqID, "-", ".")

for (i in 1:length(nums)) {
  num <- nums[i]
  col <- cols[i]
  label <- labels[i]
  load(paste0(readDir, "survey_", num, ".RData"))
  MOIAllTemp <- outputList$indLevelMOI %>% left_join(epi %>% select("SeqID", "MalTreat2"), by = c("HostID" = "SeqID")) %>%
    filter(MalTreat2 == "No")
  MOIAll <- MOIAllTemp %>% group_by(MOI) %>% summarise(n = n()) %>%
    ungroup() %>% mutate(n_total = sum(n)) %>% mutate(p = n/n_total, survey = paste0(prefix, num))
  pts2x <- ggplot(MOIAll, aes(MOI, p, col = as.factor(survey), fill = as.factor(survey))) +
    geom_histogram(stat = "identity") +
    ggtitle(label) +
    theme_classic() + theme(
      plot.title = element_text(color="black", size=sizeV, hjust = 0.5),
      axis.title.x = element_text(color="black", size=sizeV),
      axis.title.y = element_text(color="black", size=sizeV),
      axis.text.x = element_text(color="black", size=sizeV, angle=0),
      axis.text.y = element_text(color="black", size=sizeV, angle=0),
      legend.text = element_text(color="black", size=sizeV, angle=0),
      legend.title = element_text(color="black", size=sizeV, angle=0),
      strip.text = element_blank()) +
    coord_cartesian(xlim = c(1,20)) + 
    scale_x_continuous(breaks = seq(1,20,4)) +
    ylab("Proportion") +
    xlab("MOI") +
    scale_fill_manual(name = "",values = col) +
    scale_color_manual(name = "", values = col) +
    guides("fill" = "none", "color" = "none")
  print(pts2x)
  pAll[[i]] <- pts2x
}
ggsave(paste0(saveDir, "Ghana-MOI-survey_1.pdf"), pAll[[1]], width = 6, height = 6)
ggsave(paste0(saveDir, "Ghana-MOI-survey_4.pdf"), pAll[[2]], width = 6, height = 6)
ggsave(paste0(saveDir, "Ghana-MOI-survey_5.pdf"), pAll[[3]], width = 6, height = 6)




# PTS by age
rm(list = ls())
saveDir <- "/project2/pascualmm/QZ/PhD/projects/intervention/natComRevision/plots/figures/supp/fig23/PTSAgeGroups1-2-excludeDrug/"
if (!dir.exists(saveDir)) {
  dir.create(saveDir)
}
pAll <- list()
readDir <- "/project2/pascualmm/QZ/PhD/projects/intervention/natComRevision/files/figures/main/Fig5/PTSAgeGroups1-excludeDrug/"
load(file = paste0(readDir, "PTS-BC-byAge-separateMOI-1-to-20-combineYoungAgeGroups-add2015-2"))
dfAll <- PTSdfAll %>% filter(state %in% c("pre IRS", "2y into IRS", "right post IRS (2015)")) %>% 
  mutate(state = case_when(state == "2y into IRS" ~ "IRS",
                           state == "pre IRS" ~ "Pre-IRS",
                           state == "right post IRS (2015)" ~ "Post-IRS"))
dfAll$state <- factor(dfAll$state, levels = c("Pre-IRS", "IRS", "Post-IRS"))
dfAll$age <- as.character(dfAll$age)
colnames(dfAll)[2] <- "Age"
ageSel <- c("1-10", ">20")
dfAll <- dfAll %>% mutate("MOI_label" = paste0("MOI = ", MOI))
sizeV <- 28
nums <- c(1,4,5)
labels <- c("Pre-IRS", "IRS", "Post-IRS")
for (i in 1:length(nums)) {
  num <- nums[i]
  label <- labels[i]
  p <- ggplot(dfAll %>% filter(Age %in% ageSel, survey %in% c(paste0("survey_", num))), aes(x = PTS, col = Age, fill = Age)) +
    geom_density(aes(group = Age, alpha = 0.1)) +
    ggtitle(label) +
    xlab("PTS") + ylab("Density") +
    theme_classic() + theme(
      plot.title = element_text(color="black", size=sizeV, hjust = 0.5),
      axis.title.x = element_text(color="black", size=sizeV),
      axis.title.y = element_text(color="black", size=sizeV),
      axis.text.x = element_text(color="black", size=sizeV, angle=0),
      axis.text.y = element_text(color="black", size=sizeV, angle=0),
      legend.text = element_text(color="black", size=sizeV, angle=0),
      legend.title = element_text(color="black", size=sizeV, angle=0),
      strip.text = element_text(color="black", size=sizeV, angle=0)) +
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
ggsave(paste0(saveDir, "Ghana-PTS-survey_1.pdf"), pAll[[1]], width = 6, height = 6)
ggsave(paste0(saveDir, "Ghana-PTS-survey_4.pdf"), pAll[[2]], width = 6, height = 6)
ggsave(paste0(saveDir, "Ghana-PTS-survey_5.pdf"), pAll[[3]], width = 6, height = 6)



# PTS zoom in
rm(list = ls())
saveDir <- "/project2/pascualmm/QZ/PhD/projects/intervention/natComRevision/plots/figures/supp/fig23/PTSAgeGroups1Zoomin-2-excludeDrug/"
if (!dir.exists(saveDir)) {
  dir.create(saveDir)
}
pAll <- list()
readDir <- "/project2/pascualmm/QZ/PhD/projects/intervention/natComRevision/files/figures/main/Fig5/PTSAgeGroups1-excludeDrug/"
load(file = paste0(readDir, "PTS-BC-byAge-separateMOI-1-to-20-combineYoungAgeGroups-add2015-2"))
dfAll <- PTSdfAll %>% filter(state %in% c("pre IRS", "2y into IRS", "right post IRS (2015)")) %>% 
  mutate(state = case_when(state == "2y into IRS" ~ "IRS",
                           state == "pre IRS" ~ "Pre-IRS",
                           state == "right post IRS (2015)" ~ "Post-IRS"))
dfAll$state <- factor(dfAll$state, levels = c("Pre-IRS", "IRS", "Post-IRS"))
dfAll$age <- as.character(dfAll$age)
colnames(dfAll)[2] <- "Age"
ageSel <- c("1-10", ">20")
dfAll <- dfAll %>% mutate("MOI_label" = paste0("MOI = ", MOI))
sizeV <- 28
nums <- c(1,4,5)
labels <- c("Pre-IRS", "IRS", "Post-IRS")
for (i in 1:length(nums)) {
  num <- nums[i]
  label <- labels[i]
  p <- ggplot(dfAll %>% filter(Age %in% ageSel, survey %in% c(paste0("survey_", num))), aes(x = PTS, col = Age, fill = Age)) +
    geom_density(aes(group = Age, alpha = 0.1)) +
    ggtitle(label) +
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
      legend.position = c(0.65, 0.65)) +
    scale_fill_manual(values=c("magenta", "dark green")) +
    scale_color_manual(values=c("magenta", "dark green")) +
    coord_cartesian(xlim = c(0,0.075)) + scale_alpha(guide = 'none') + 
    guides(fill = "none", col = "none")
  print(p)
  pAll[[i]] <- p
}
ggsave(paste0(saveDir, "Ghana-PTS-survey_1.pdf"), pAll[[1]] + rremove("xlab") + rremove("ylab"), width = 6, height = 6)
ggsave(paste0(saveDir, "Ghana-PTS-survey_4.pdf"), pAll[[2]] + rremove("xlab") + rremove("ylab"), width = 6, height = 6)
ggsave(paste0(saveDir, "Ghana-PTS-survey_5.pdf"), pAll[[3]] + rremove("xlab") + rremove("ylab"), width = 6, height = 6)

