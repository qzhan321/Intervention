rm(list = ls())
suppressPackageStartupMessages({
  library(ggplot2)
  library(gridExtra)
  library(scales)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(ggpubr)
  library(grid)
  library(cowplot)
  library(Hmisc)
  library(cmocean)
})
run <- 4
source(paste0("/home/qizhan/others/PhD/projects/intervention/writings/simulation", run, "/plots/colorFunc.R"))
readDir <- paste0("/project2/pascualmm/QZ/PhD/projects/intervention/files", run, "/actualRuns/epi/")
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
IRSDur <- 2
IRSType <- paste0(IRSDur, "yIRS")
postIRS <- 10-IRSDur

numsPlotList <- list(
  "seasonal closed" = list(1:14),
  "seasonal semi-open" = list(1:14,101:115),
  "seasonal regionally-open" = list(1:12,19:34,37:48,55:70),
  "non-seasonal closed" = list(1:13),
  "non-seasonal semi-open" = list(1:15, 101:115),
  "non-seasonal regionally-open" = list(1:13,19:34,37:49,55:70)
)
labelsList <- list(
  "seasonal closed" = list("Seasonal \nClosed \n \n"),
  "seasonal semi-open" = list("Seasonal \nSemi-open \nBaseline Migration \n", "Seasonal \nSemi-open \nHigh Migration \n"),
  "seasonal regionally-open" = list("Seasonal \nRegionally-open \nBaseline Migration \nMedium Pool",
                                    "Seasonal \nRegionally-open \nBaseline Migration \nLarge Pool",
                                    "Seasonal \nRegionally-open \nHigh Migration \nMedium Pool",
                                    "Seasonal \nRegionally-open \nHigh Migration \nLarge Pool"),
  "non-seasonal closed" = list("Non-seasonal \nClosed \n \n"),
  "non-seasonal semi-open" = list("Non-seasonal \nSemi-open \nBaseline Migration \n", "Non-seasonal \nSemi-open \nHigh Migration \n"),
  "non-seasonal regionally-open" = list("Non-seasonal \nRegionally-open \nBaseline Migration \nMedium Pool",
                                        "Non-seasonal \nRegionally-open \nBaseline Migration \nLarge Pool",
                                        "Non-seasonal \nRegionally-open \nHigh Migration \nMedium Pool",
                                        "Non-seasonal \nRegionally-open \nHigh Migration \nLarge Pool")
)
immsList <- list(
  "seasonal closed" = list("Three-year"),
  "seasonal semi-open" = list("Three-year",
                              "Three-year"),
  "seasonal regionally-open" = list("Three-year", "Three-year", "Three-year", "Three-year"),
  "non-seasonal closed" = list("Three-year"),
  "non-seasonal semi-open" = list("Three-year",
                                  "Three-year"),
  "non-seasonal regionally-open" = list("Three-year", "Three-year", "Three-year", "Three-year")
)

T_YEAR <- 360

N <- 10000
pAll1 <- list()
pAll2 <- list()
sizeV <- 31

saveDir0 <- paste0("/project2/pascualmm/QZ/PhD/projects/intervention/writings/simulation", run, "/plots/")
dir.create(saveDir0)
saveDir1 <- paste0(saveDir0, "prevalenceDynamicsTransientIRS/")
dir.create(saveDir1)

for (i in 1:length(seasonality)) {
  s <- seasonality[i]
  saveDir2 <- paste0(saveDir1, s, "/")
  dir.create(saveDir2)
  
  for (j in 1:length(openness)) {
    o <- openness[j]
    saveDir3 <- paste0(saveDir2, o, "/")
    dir.create(saveDir3)
    
    file <- paste0(readDir, s, "/", o, "/", "summaryInfoTable-", IRSType, ".RData")
    load(file)
    
    preIRS <- preIRSList[[paste0(s, " ", o)]]
    numsPlot <- numsPlotList[[paste0(s, " ", o)]]
    labels <- labelsList[[paste0(s, " ", o)]]
    imms <- immsList[[paste0(s, " ", o)]]
    for (k in 1:length(numsPlot)) {
      nums <- numsPlot[[k]]
      label <- labels[[k]]
      imm <- imms[[k]]
      
      saveDir4 <- paste0(saveDir3, imm, "/")
      dir.create(saveDir4)
      
      summaryTable1 <- summaryTableCombined %>% filter(num %in% nums)
      summaryTable2 <- summaryTable1 %>% filter(time > (preIRS-3)*T_YEAR, time <= (preIRS + IRSDur + postIRS)*T_YEAR)
      summaryTable3 <- summaryTable2 %>% dplyr::group_by(time, IRS) %>%
        dplyr::summarise(mP = mean(Prevalence), vP = sd(Prevalence), 
                         maxP = max(Prevalence), minP = min(Prevalence),
                         mD = mean(n_circulating_genes), vD = sd(n_circulating_genes),
                         maxD = max(n_circulating_genes), minD = min(n_circulating_genes),
                         mMOI = mean(MOI, na.rm = T), vMOI = sd(MOI, na.rm = T))
      summaryTable4 <- summaryTable3 %>% mutate(vP = ifelse(is.na(vP), 0, vP), 
                                                vD = ifelse(is.na(vD), 0, vD),
                                                vMOI = ifelse(is.na(vMOI), 0, vMOI))
      summaryTable5 <- summaryTable4 %>% mutate(timePlot = time/T_YEAR - preIRS) 
      summaryTable5$IRS <- factor(summaryTable5$IRS, levels = c("preIRS", paste0("I-", 1:16)))
      p1 <- ggplot(summaryTable5, aes(timePlot, mP, col = IRS, fill = IRS))+ 
        geom_line(size = 1.9)+
        geom_ribbon(aes(y = mP, ymin = minP, ymax = maxP), alpha = .2, linetype = 0) +
        # ggtitle(label) +
        xlab("Years Since IRS Starts \n(Sustained IRS)") + ylab("Prevalence") +
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
        coord_cartesian(ylim = c(0, 1)) +
        scale_color_manual(values = Turbo(out.colors = 17)[1:(length(nums)+1)]) +
        scale_fill_manual(values = Turbo(out.colors = 17)[1:(length(nums)+1)]) +
        scale_x_continuous(breaks = seq(-2, 10, 2))
      p1 <- p1 + guides(fill = "none", col = "none")
      print(p1)
      index <- length(pAll1)
      pAll1[[index + 1]] <- p1
      
      
      forwardYs <- 9
      if (s == "seasonal") {
        df1 <- summaryTable5 %>% filter(time == (preIRS + forwardYs)*T_YEAR + 300)
      } else {
        df1 <- summaryTable5 %>% filter(time == (preIRS + forwardYs + 1)*T_YEAR)
      }
      summaryTable2 <- summaryTable1 %>% filter(!(time == preIRS * T_YEAR & IRS %in% c(paste0("I-", 1:16))))
      CRs <- summaryTable2 %>% select(IRS, CR) %>% group_by(IRS) %>% summarise(meanCR = mean(CR))
      df2 <- df1 %>% left_join(CRs, by = "IRS")
      df2$IRS <- factor(df2$IRS, levels = paste0("I-", 1:16))
      dfspline <- as.data.frame(spline(df2$meanCR, df2$mP))
      epiCR<-ggplot(df2, aes(x=meanCR, y=mP, col = IRS)) + 
        geom_point(size = 5) + 
        geom_pointrange(aes(ymin=minP, ymax=maxP)) +
        # geom_line(data = dfspline, aes(x=x, y=y), color = "gray88", size = 1.9) + 
        # ggtitle(label) +
        theme_bw() +
        xlab("Effective Contact Rate \n(per Host, per Year)") + ylab("Prevalence") +
        theme(
          plot.title = element_text(color="black", size=sizeV, hjust = 0.5, vjust = 1),
          axis.title.x = element_text(color="black", size=sizeV),
          axis.title.y = element_text(color="black", size=sizeV),
          axis.text.x = element_text(color="black",
                                     size=sizeV, angle=0),
          axis.text.y = element_text(color="black",
                                     size=sizeV, angle=0),
          legend.text = element_text(color="black",
                                     size=sizeV, angle=0),
          legend.title = element_blank(),
          strip.text = element_text(color="black",
                                    size=sizeV, angle=0),
          legend.position = c(0.9, 0.04)) + 
        scale_color_manual(values = Turbo(out.colors = 17)[2:(length(nums) + 1)]) +
        scale_fill_manual(values = Turbo(out.colors = 17)[2:(length(nums) + 1)]) +
        guides(col = "none", fill = "none") +
        coord_cartesian(ylim = c(0, 1))
      # scale_x_continuous(breaks = seq(2, 13, 2))) 
      epiCR
      index <- length(pAll2)
      pAll2[[index + 1]] <- epiCR
      
      p1_rmL <- p1 + guides(fill = "none", col = "none") + rremove("legend") 
      p2_rmL <- epiCR + guides(fill = "none", col = "none") + rremove("legend") 
      
      ggsave(paste0(saveDir4, "num-", min(nums), "-", max(nums), "-time.pdf"), p1_rmL, width = 6, height = 6)
      # ggsave(paste0(saveDir4, "num-", min(nums), "-", max(nums), "-cr.pdf"), p2_rmL, width = 6, height = 6)
      
    }
  }
}



rm(list=ls())
T_YEAR <- 360
IRSDur <- 2
IRSType <- paste0(IRSDur, "yIRS")
postIRS <- 10-IRSDur
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
  summaryTable1 <- summaryTableCombined
  summaryTable1 <- summaryTable1 %>% mutate(imm = ifelse(num %in% c(1:99, 1001:1099), "3y", "4y"))
  summaryTable1 <- summaryTable1 %>% mutate(migration = ifelse(num %in% c(1:18, 101:118, 1001:1018, 1101:1118), "Baseline Migration Rate", "High Migration Rate"))
  summaryTable2 <- summaryTable1 %>% filter(time > (preIRS-3)*T_YEAR, time <= (preIRS + IRSDur + postIRS)*T_YEAR)
  summaryTable3 <- summaryTable2 %>% dplyr::group_by(time, imm, state, migration) %>% dplyr::summarise(
    mP = mean(Prevalence), vP = sd(Prevalence), 
    maxP = max(Prevalence), minP = min(Prevalence),
    mD = mean(n_circulating_genes), vD = sd(n_circulating_genes),
    maxD = max(n_circulating_genes), minD = min(n_circulating_genes),
    mMOI = mean(MOI, na.rm = T), vMOI = sd(MOI, na.rm = T))
  summaryTable4 <- summaryTable3 %>% mutate(vP = ifelse(is.na(vP), 0, vP), 
                                            vD = ifelse(is.na(vD), 0, vD),
                                            vMOI = ifelse(is.na(vMOI), 0, vMOI))
  summaryTable5 <- summaryTable4 %>% mutate(timePlot = time/T_YEAR - preIRS) 
  dfAll <- rbind(dfAll, summaryTable5 %>% mutate(seasonality = capitalize(s)))
}
dfAll$state <- factor(dfAll$state, levels = c("Static", "Dynamic"))
dfAll$seasonality <- factor(dfAll$seasonality, levels = c("Seasonal", "Non-seasonal"))
dfAll$migration <- factor(dfAll$migration, levels = c("Baseline Migration Rate", "High Migration Rate"))
sizeV <- 31
p1 <- ggplot(dfAll %>% filter(migration == "High Migration Rate", imm == "3y"), aes(time/T_YEAR - preIRS, mP, linetype = state, fill = seasonality, col = seasonality))+ 
  geom_line(size = 1.9)+
  # facet_wrap(.~immRate, nrow=1)+
  geom_ribbon(aes(y = mP, ymin = minP, ymax = maxP), alpha = .2, linetype = 0) +
  ggtitle("High Migration Rate") +
  xlab("Years Since IRS Starts \n(Transient IRS)") + ylab("Prevalence") +
  theme_bw() + theme(
    plot.title = element_text(color="black", size=sizeV, hjust = 0.5),
    axis.title.x = element_text(color="black", size=sizeV),
    axis.title.y = element_text(color="white", size=sizeV),
    axis.text.x = element_text(color="black", 
                               size=sizeV, angle=0),
    axis.text.y = element_text(color="white", 
                               size=sizeV, angle=0),
    legend.text = element_text(color="black", 
                               size=sizeV, angle=0),
    legend.title = element_text(color="black", 
                                size=sizeV, angle=0),
    strip.text = element_text(color="black", 
                              size=sizeV, angle=0)) + coord_cartesian(ylim = c(0, 1)) +
  scale_x_continuous(breaks = seq(-2, 10, 2))
p1 <- p1 + guides(col = "none",
                  fill = "none",
                  linetype = "none")
# fill = guide_legend(title = NULL, order = 1, nrow = 2, keyheight = unit(0.1, 'inch')),
# linetype = guide_legend(title = NULL, order = 2, nrow = 2, keyheight = unit(0.1, 'inch'))) + 
# theme(legend.position = c(0.68, 0.075), legend.box = "horizontal")
p1 <- p1 + scale_color_manual(values = c("dark red", "dark red")) + scale_fill_manual(values = c("dark red", "dark red"))
print(p1)

p2 <- ggplot(dfAll %>% filter(migration == "Baseline Migration Rate", imm == "3y"), aes(time/T_YEAR - preIRS, mP, linetype = state, fill = seasonality, col = seasonality))+ 
  geom_line(size = 1.9)+
  # facet_wrap(.~immRate, nrow=1)+
  geom_ribbon(aes(y = mP, ymin = minP, ymax = maxP), alpha = .2, linetype = 0) +
  ggtitle("Baseline Migration Rate") +
  xlab("Years Since IRS Starts \n(Transient IRS)") + ylab("Prevalence") +
  theme_bw() + theme(
    plot.title = element_text(color="black", size=sizeV, hjust = 0.5),
    axis.title.x = element_text(color="black", size=sizeV),
    axis.title.y = element_text(color="black", size=sizeV),
    axis.text.x = element_text(color="black", 
                               size=sizeV, angle=0),
    axis.text.y = element_text(color="black", 
                               size=sizeV, angle=0),
    legend.text = element_text(color="black", 
                               size=sizeV/1.6, angle=0),
    legend.title = element_text(color="black", 
                                size=sizeV/1.6, angle=0),
    strip.text = element_text(color="black", 
                              size=sizeV, angle=0), 
    legend.box = "horizontal",
    legend.key = element_rect(colour = NA, fill = NA),
    legend.position = c(0.5, 0.998),
    legend.background = element_rect(fill='transparent', color = NA),
    legend.box.background = element_rect(fill='transparent', color = NA)) + 
  coord_cartesian(ylim = c(0, 1)) + 
  scale_x_continuous(breaks = seq(-2, 10, 2)) 
p2 <- p2 + guides(col = "none",
                  fill = "none",
                  linetype = "none")  
  # scale_linetype_manual(name = "", values = c("dashed", "solid"), label = c("Static Pool", "Dynamic Pool")) 
p2 <- p2 + scale_color_manual(values = c("dark red", "dark red")) + scale_fill_manual(values = c("dark red", "dark red"))
p1 <- p1 #+ rremove("xlab")
p2 <- p2 #+ guides(linetype = guide_legend(nrow = 1)) + scale_linetype_manual(name = "", values = c("dashed", "solid"), label = c("Static Pool", "Dynamic Pool"))
saveDir0 <- paste0("/project2/pascualmm/QZ/PhD/projects/intervention/writings/simulation", run, "/plots/")
dir.create(saveDir0)
saveDir1 <- paste0(saveDir0, "prevalenceDynamicsTransientIRSPool/")
dir.create(saveDir1)
ggsave(paste0(saveDir1, "transientIRS-pool-1.pdf"), p2, width = 6, height = 6)
ggsave(paste0(saveDir1, "transientIRS-pool-2.pdf"), p1, width = 6, height = 6)




rm(list = ls())
run <- 4
readDir <- paste0("/project2/pascualmm/QZ/PhD/projects/intervention/files", run, "/actualRuns/PTS/")
prefix<-"sim"
nums <- c(1,7,10)
preIRS_temp <- 200
seasonality <- "seasonal"
openness <- "closed"
quantiles <- round(c(seq(0.01,0.9, (0.9-0.01)/3)),2)
colfunc <- colorRampPalette(c("white", "gray75", "yellow3", "darksalmon", "dark red", "darkgreen", "dark blue", "black"))
colorsGray <- colfunc(100)
rows <- 1:1
cols <- 1:4
IRSDur <- 2
postIRS <- 8
T_YEAR <- 360

k = 100
colorsGray = rev(cmocean('deep')(k)[k:1])
type <- "dir"
state <- "true"
pAll3 <- list()
sizeV <- 31
for (q in 1:length(nums)) {
  num <- nums[q]
  print(num)
  load(file = paste0(readDir, seasonality, "/", openness, "/", type, "/", state, "/", prefix, "_", num, "_r0.RData"))
  
  rows_temp <- NULL
  cols_temp <- NULL
  for (k in 1:length(quantiles)) {
    a <- (k - 1) %% max(rows) 
    b <- (k - 1) %/% max(rows) 
    row <- rows[a + 1]
    col <- cols[b + 1]
    rows_temp <- c(rows_temp, row)
    cols_temp <- c(cols_temp, col)
  }
  
  quantiles_df <- data.frame("row" = rows_temp, "col" = cols_temp, 
                             "quantile" = quantiles)
  
  dfAll <- NULL
  for (m in 1:length(PTSDirMatrixList)) {
    PTStemp <- PTSDirMatrixList[[m]]
    if (length(PTStemp) > 1) {
      t12 <- names(PTSDirMatrixList)[m]
      t12 <- strsplit(t12, split = "-")
      t1 <- as.numeric(unlist(t12)[1])
      t2 <- as.numeric(unlist(t12)[2])
      if (t1 == t2) {
        PTStempV <- PTStemp[row(PTStemp)!=col(PTStemp)]
        PTSquans <- quantile(PTStempV, quantiles)
      } else {
        PTStempV <- as.vector(PTStemp)
        PTSquans <- quantile(PTStempV, quantiles)
      }
      
      df1 <- data.frame("PTS" = PTSquans, "t1" = t1, "t2" = t2, 
                        "quantile" = quantiles, "num" = num)
      df2 <- data.frame("PTS" = PTSquans, "t1" = t2, "t2" = t1, 
                        "quantile" = quantiles, "num" = num)
      df <- rbind(df1, df2)
      dfAll <- rbind(dfAll, df)
    }
  }
  
  dfAll <- left_join(dfAll, quantiles_df, by = "quantile")
  dfAll <- dfAll %>% filter(t1 <= (preIRS_temp + IRSDur + postIRS)*T_YEAR, 
                            t2 <= (preIRS_temp + IRSDur + postIRS)*T_YEAR)
  
  dat_text <- data.frame(
    row = rep(1:max(dfAll$row), each = max(dfAll$col)),
    col = rep(1:max(dfAll$col), max(dfAll$row)) 
  )
  dat_text$label <- format(quantiles, nsmall = 2)
  dat_text$PTS <- 0
  
  times <- seq(-3, IRSDur + postIRS - 1, 1) + round(300/360,2)
  
  dfAll <- dfAll %>% mutate(y1 = t1/T_YEAR - preIRS_temp, y2 = t2/T_YEAR - preIRS_temp)
  if (num == min(nums)) {
    p <- ggplot(dfAll, aes(x=y1, y=y2, fill = PTS)) + geom_tile() + 
      scale_fill_gradientn(colours = colorsGray, values = seq(0,1,(1-0)/99), 
                           limits = c(0,1)) + facet_grid(row~col) + 
      theme(
        plot.title = element_text(color="black", size=sizeV, hjust = 0.5),
        axis.title.x = element_text(color="black", size=sizeV),
        axis.title.y = element_text(color="black", size=sizeV),
        axis.text.x = element_text(color="black", 
                                   size=sizeV, angle=0),
        axis.text.y = element_text(color="black", 
                                   size=sizeV, angle=0),
        legend.text = element_text(color="black", 
                                   size=sizeV/1.25, angle=0),
        legend.title = element_text(color="black", 
                                    size=sizeV/1.25, angle=0),
        strip.text = element_text(color="black", 
                                  size=sizeV, angle=0)) + 
      xlab("Years since IRS starts (Transient IRS)") + ylab("Years since IRS starts") +
      # annotate('text', label = as.character(quantiles), x = min(dfAll$t1), y = max(dfAll$t2))
      geom_text(
        size    = sizeV/3,
        data    = dat_text,
        mapping = aes(x = Inf, y = Inf, label = label),
        hjust   = 2.25,
        vjust   = 2.5
        # fontface = "bold"
      ) + 
      # ggtitle(paste0("PTS quantiles")) +
      scale_x_continuous(breaks = round(times[seq(1, length(times), 2)])) +
      scale_y_continuous(breaks = round(times[seq(1, length(times), 2)])) +
      geom_hline(yintercept = c(0, IRSDur), col = "gray75", linetype = "dashed") +
      geom_vline(xintercept = c(0, IRSDur), col = "gray75", linetype = "dashed") +
      theme(strip.background = element_blank(), strip.text = element_blank())
    print(p)
  } else {
    p <- ggplot(dfAll, aes(x=y1, y=y2, fill = PTS)) + geom_tile() + 
      scale_fill_gradientn(colours = colorsGray, values = seq(0,1,(1-0)/99), 
                           limits = c(0,1)) + facet_grid(row~col) + 
      theme(
        plot.title = element_text(color="black", size=sizeV, hjust = 0.5),
        axis.title.x = element_text(color="black", size=sizeV),
        axis.title.y = element_text(color="black", size=sizeV),
        axis.text.x = element_text(color="black", 
                                   size=sizeV, angle=0),
        axis.text.y = element_text(color="black", 
                                   size=sizeV, angle=0),
        legend.text = element_text(color="black", 
                                   size=sizeV/1.25, angle=0),
        legend.title = element_text(color="black", 
                                    size=sizeV/1.25, angle=0),
        strip.text = element_text(color="black", 
                                  size=sizeV, angle=0)
      ) + 
      xlab("Years since IRS starts (Transient IRS)") + ylab("Years since IRS starts") +
      # ggtitle(paste0("PTS quantiles")) +
      scale_x_continuous(breaks = round(times[seq(1, length(times), 2)])) +
      scale_y_continuous(breaks = round(times[seq(1, length(times), 2)])) +
      geom_hline(yintercept = c(0, IRSDur), col = "gray75", linetype = "dashed") +
      geom_vline(xintercept = c(0, IRSDur), col = "gray75", linetype = "dashed") +
      theme(strip.background = element_blank(), strip.text = element_blank())
    print(p)
  }
  pAll3[[q]] <- p
}

p1_rmL <- pAll3[[1]] + rremove("xlab") + rremove("ylab") + rremove("legend") + ggtitle("Fast-rebound Regime")
p1_rmL
p2_rmL <- pAll3[[2]] + rremove("xlab") + rremove("ylab") + rremove("legend") + ggtitle("Transition Regime") 
p2_rmL
p3_rmL <- pAll3[[3]] + rremove("legend") + rremove("xlab") + rremove("ylab") + ggtitle("Slow-rebound Regime")
p3_rmL
p4_legend <- get_legend(pAll3[[1]])

saveDir0 <- paste0("/project2/pascualmm/QZ/PhD/projects/intervention/writings/simulation", run, "/plots/")
dir.create(saveDir0)
saveDir1 <- paste0(saveDir0, "/PTSHeatmap/")
dir.create(saveDir1)
ggsave(paste0(saveDir1, "transientIRS-heatmap-1.pdf"), p1_rmL, width = 12, height = 3.85)
ggsave(paste0(saveDir1, "transientIRS-heatmap-2.pdf"), p2_rmL, width = 12, height = 3.85)
ggsave(paste0(saveDir1, "transientIRS-heatmap-3.pdf"), p3_rmL, width = 12, height = 3.85)
ggsave(paste0(saveDir1, "transientIRS-heatmap-lg.pdf"), p4_legend, width = 6, height = 6)

