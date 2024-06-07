rm(list = ls())
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(Hmisc)
})
T_YEAR <- 360
IRSDur <- 2
IRSType <- paste0(IRSDur, "yIRS")
postIRS <- 10 - IRSDur
run <- 4
source(paste0("/home/qizhan/others/PhD/projects/intervention/writings/simulation", run, "/plots/colorFunc.R"))
readDir <- paste0("/project2/pascualmm/QZ/PhD/projects/intervention/files", run, "/actualRuns/epi/")
openness <- "regionally-open"
seasonality <- c("seasonal", "non-seasonal")
file <- paste0("summaryInfoTable-", IRSType, "-pool.RData")
preIRS <- 150
dfAll <- NULL
for (i in 1:length(seasonality)) {
  s <- seasonality[i]
  load(paste0(readDir, s, "/", openness, "/", file))
  
  summaryTable1 <- summaryTableCombined %>% filter(num %in% c(1:99, 1001:1099))
  summaryTable2 <- summaryTable1 %>% mutate(migration = ifelse(num %in% c(1:18, 1001:1018), "Baseline Migration Rate", "High Migration Rate"))
  summaryTable3 <- summaryTable2 %>% filter(time > (preIRS-3)*T_YEAR, time <= (preIRS + IRSDur + postIRS)*T_YEAR)
  summaryTable4 <- summaryTable3 %>% group_by(time, state, migration) %>% summarise(
    mP = mean(Prevalence), vP = sd(Prevalence),
    maxP = max(Prevalence), minP = min(Prevalence),
    mD = mean(n_circulating_genes), vD = sd(n_circulating_genes),
    maxD = max(n_circulating_genes), minD = min(n_circulating_genes),
    mMOI = mean(MOI, na.rm = T), vMOI = sd(MOI, na.rm = T))
  summaryTable5 <- summaryTable4 %>% mutate(vP = ifelse(is.na(vP), 0, vP), vD = ifelse(is.na(vD), 0, vD), vMOI = ifelse(is.na(vMOI), 0, vMOI))
  summaryTable6 <- summaryTable5 %>% mutate(timePlot = time/T_YEAR - preIRS)
  dfAll <- rbind(dfAll, summaryTable6 %>% mutate(seasonality = capitalize(s)))
}
dfAll$state <- factor(dfAll$state, levels = c("Static", "Dynamic"))
dfAll$seasonality <- factor(dfAll$seasonality, levels = c("Seasonal", "Non-seasonal"))
dfAll$migration <- factor(dfAll$migration, levels = c("Baseline Migration Rate", "High Migration Rate"))
sizeV <- 25
pointSize <- 1.25
p1 <- ggplot(dfAll %>% filter(migration == "High Migration Rate"), 
             aes(timePlot, mP, linetype = state, fill = seasonality, col = seasonality)) +
  geom_line(size = pointSize)+
  geom_ribbon(aes(y = mP, ymin = minP, ymax = maxP), alpha = .4) +
  ggtitle("High Migration Rate") +
  xlab("Years Since IRS Starts \n(Transient IRS)") + ylab("Prevalence") +
  theme_bw() + theme(
    plot.title = element_text(color="black", size=sizeV, hjust = 0.5),
    axis.title.x = element_text(color="black", size=sizeV),
    axis.title.y = element_text(color="black", size=sizeV),
    axis.text.x = element_text(color="black", size=sizeV, angle=0),
    axis.text.y = element_text(color="black", size=sizeV, angle=0),
    legend.text = element_text(color="black", size=sizeV/1.25, angle=0),
    legend.title = element_text(color="black", size=sizeV/1.25, angle=0),
    strip.text = element_text(color="black", size=sizeV, angle=0),
    legend.box = "vertical",
    legend.direction = "vertical",
    legend.key = element_rect(colour = NA, fill = NA),
    legend.background = element_rect(fill='transparent', color = NA),
    legend.box.background = element_rect(fill='transparent', color = NA)) +
  coord_cartesian(ylim = c(0, 1)) +
  scale_x_continuous(breaks = seq(-2, 10, 2)) + 
  scale_color_manual(values = c("dark red", "dark green"), name = "Seasonality") + 
  scale_fill_manual(values = c("dark red", "dark green"), name  = "Seasonality") +
  scale_linetype_manual(values = c(22, "solid"), name = "Regional Pool") 

saveDir0 <- "/project2/pascualmm/QZ/PhD/projects/intervention/natComRevision/plots/figures/supp/fig26/prevalenceDynamicsTransientIRSPool/"
if (!dir.exists(saveDir0)) {
  dir.create(saveDir0)
}
ggsave(paste0(saveDir0, "transientIRS-pool-1.pdf"), p1, width = 9, height = 6)
