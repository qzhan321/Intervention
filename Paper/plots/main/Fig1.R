rm(list = ls())
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(ggpubr)
})
saveDir0 <- paste0("/project2/pascualmm/QZ/PhD/projects/intervention/natComRevision/plots/figures/main/")
if (!dir.exists(saveDir0)) {
  dir.create(saveDir0)
}
saveDir1 <- paste0(saveDir0, "Fig1/")
if (!dir.exists(saveDir1)) {
  dir.create(saveDir1)
}

run <- 4
source(paste0("/home/qizhan/others/PhD/projects/intervention/writings/simulation", run, "/plots/colorFunc.R"))
readDir <- paste0("/project2/pascualmm/QZ/PhD/projects/intervention/files", run, "/actualRuns/epi/")
openness <- c("closed")
seasonality <- c("seasonal")

preIRS <- 200
T_YEAR <- 360

N <- 10000
T_YEAR <- 360
sizeV <- 35

IRSTypes <- c("10yIRS", "2yIRS")
pAll1 <- list()
pAll2 <- list()
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
    
    for (k in 1:length(IRSTypes)) {
      IRSType <- IRSTypes[k]
      if (IRSType == "2yIRS") {
        nums <- 1:14
      } else {
        nums <- 101:110
      }
      file <- paste0(readDir, s, "/", o, "/", "summaryInfoTable-", IRSType, ".RData")
      load(file)
      if (IRSType == "10yIRS") {
        summaryTable1 <- summaryTableCombined %>% filter(num %in% nums)
        summaryTable2 <- summaryTable1 %>% filter(!(time == preIRS * T_YEAR & IRS %in% paste0("I-", 1:14)))
        CRs <- summaryTable2 %>% select(IRS, CR) %>% group_by(IRS) %>% 
          summarise(meanCR = mean(CR)) %>%
          mutate(logmeanCR = log(meanCR))
        CRs <- CRs %>% mutate(IRS = ifelse(IRS == "preIRS", "Pre-IRS", IRS))
        CRs$IRS <- factor(CRs$IRS, levels = c("Pre-IRS", paste0("I-", 1:14)))
        p1 <- ggplot(CRs, aes(IRS, logmeanCR, col = IRS, fill = IRS))+ 
          geom_bar(stat = "identity")+
          xlab("") + ylab("Effective Contact Rate \nper Host per Year \n(log-scale)") +
          theme_bw() + theme(
            plot.title = element_text(color="black", size=sizeV, hjust = 0.5),
            axis.title.x = element_text(color="black", size=sizeV),
            axis.title.y = element_text(color="black", size=sizeV),
            axis.text.x = element_text(color="black", size=sizeV, angle=-30),
            axis.text.y = element_text(color="black", size=sizeV, angle=0),
            legend.text = element_text(color="black", size=sizeV, angle=0),
            legend.title = element_text(color="black", size=sizeV, angle=0),
            strip.text = element_text(color="black", size=sizeV, angle=0),
            plot.margin = unit(c(53,5,-24,5), "points")) +
          coord_cartesian(ylim = c(0, 3.1)) +
          scale_color_manual(values = Turbo(out.colors = 17)[1:(length(nums) + 1)]) +
          scale_fill_manual(values = Turbo(out.colors = 17)[1:(length(nums) + 1)]) +
          guides(fill = "none", col = "none")
        print(p1)
        index <- length(pAll1)
        pAll1[[index + 1]] <- p1
      }
      
      summaryTable1 <- summaryTableCombined %>% filter(num %in% nums)
      CRs2 <- summaryTable1 %>% group_by(time, IRS) %>%
        summarise(mCR = mean(n_total_bites), vCR = sd(n_total_bites),
                  maxCR = max(n_total_bites), minCR = min(n_total_bites)) 
      CRs2 <- CRs2 %>% mutate(IRS = ifelse(IRS == "preIRS", "Pre-IRS", IRS))
      CRs2$IRS <- factor(CRs2$IRS, levels = c("Pre-IRS", paste0("I-", 1:14)))
      p2 <- ggplot(CRs2, aes(x=time/T_YEAR - preIRS, y=mCR/N/30, col = IRS))+ 
        geom_line(size = 1.9)+
        geom_ribbon(aes(y = mCR/N/30, ymin = minCR/N/30, ymax = maxCR/N/30), alpha = .2, linetype = 0) +
        xlab("Years Since IRS Starts (Sustained IRS)") + ylab("Effective Contact Rate per Host per Day") +
        theme_bw() + theme(
          plot.title = element_text(color="black", size=sizeV, hjust = 0.5),
          axis.title.x = element_text(color="black", size=sizeV),
          axis.title.y = element_text(color="black", size=sizeV),
          axis.text.x = element_text(color="black", size=sizeV, angle=0),
          axis.text.y = element_text(color="black", size=sizeV, angle=0),
          legend.text = element_text(color="black", size=sizeV, angle=0),
          legend.title = element_text(color="black", size=sizeV, angle=0),
          strip.text = element_text(color="black", size=sizeV, angle=0)) + 
        theme(legend.position = c(0.6, 0.89)) +
        scale_color_manual(values = Turbo(out.colors = 17)[1:(length(nums)+1)]) +
        scale_fill_manual(values = Turbo(out.colors = 17)[1:(length(nums)+1)]) +
        guides(col = "none", fill = "none") + 
        scale_x_continuous(breaks = seq(-2,10,2))
      print(p2)
      index <- length(pAll2)
      pAll2[[index + 1]] <- p2
    }
  }
}
part1 <- pAll1[[1]] 
part2 <- pAll2[[1]] + rremove("xlab") + rremove("ylab") + theme(axis.text.x = element_text(color="white", size=sizeV, angle=0))
part3 <- pAll2[[2]] + rremove("xlab") + rremove("ylab")
ggsave(paste0(saveDir3, "IRSContactRate-1.pdf"), part1, width = 10.5, height = 6)
ggsave(paste0(saveDir3, "IRSContactRate-2.pdf"), part2, width = 6.4, height = 5.5)
ggsave(paste0(saveDir3, "IRSContactRate-3.pdf"), part3, width = 6.4, height = 5.5)

