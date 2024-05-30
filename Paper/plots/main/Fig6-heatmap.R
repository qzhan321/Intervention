rm(list = ls())
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(cmocean)
  library(cowplot)
  library(ggpubr)
})
run <- 4
readDir <- paste0("/project2/pascualmm/QZ/PhD/projects/intervention/files", run, "/actualRuns/PTS/")
prefix<-"sim"
nums <- c(1,7,8,9,10)
preIRS <- 200
seasonality <- "seasonal"
openness <- "closed"
quantiles <- round(c(seq(0.01,0.9, (0.9-0.01)/3)),2)
rows <- 1:1
cols <- 1:4
IRSDur <- 2
postIRS <- 8
T_YEAR <- 360
times <- seq(-3, IRSDur + postIRS - 1, 1) + round(300/360,2)
k = 100
colorsGray = rev(cmocean('deep')(k)[k:1])
type <- "dir"
state <- "true"
sizeV <- 25
pAll <- list()
for (i in 1:length(nums)) {
  num <- nums[i]
  print(num)
  load(file = paste0(readDir, seasonality, "/", openness, "/", type, "/", state, "/", prefix, "_", num, "_r0.RData"))
  
  rows_temp <- NULL
  cols_temp <- NULL
  for (j in 1:length(quantiles)) {
    a <- (j - 1) %% max(rows) 
    b <- (j - 1) %/% max(rows) 
    row <- rows[a + 1]
    col <- cols[b + 1]
    rows_temp <- c(rows_temp, row)
    cols_temp <- c(cols_temp, col)
  }
  quantiles_df <- data.frame("row" = rows_temp, "col" = cols_temp, "quantile" = quantiles)
  
  dfAll <- NULL
  for (k in 1:length(PTSDirMatrixList)) {
    PTStemp <- PTSDirMatrixList[[k]]
    if (length(PTStemp) > 1) {
      t12 <- names(PTSDirMatrixList)[k]
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
      
      df1 <- data.frame("PTS" = unname(PTSquans), "t1" = t1, "t2" = t2, "quantile" = quantiles, "num" = num)
      df2 <- data.frame("PTS" = unname(PTSquans), "t1" = t2, "t2" = t1, "quantile" = quantiles, "num" = num)
      df <- rbind(df1, df2)
      dfAll <- rbind(dfAll, df)
    }
  }
  
  dfAll <- left_join(dfAll, quantiles_df, by = "quantile")
  dfAll <- dfAll %>% filter(t1 <= (preIRS + IRSDur + postIRS)*T_YEAR, 
                            t2 <= (preIRS + IRSDur + postIRS)*T_YEAR)
  
  dat_text <- data.frame(
    row = rep(1:max(dfAll$row), each = max(dfAll$col)),
    col = rep(1:max(dfAll$col), max(dfAll$row)) 
  )
  dat_text$label <- format(quantiles, nsmall = 2)
  dat_text$PTS <- 0
  
  dfAll <- dfAll %>% mutate(y1 = t1/T_YEAR - preIRS, y2 = t2/T_YEAR - preIRS)
  if (num == min(nums)) {
    p <- ggplot(dfAll, aes(x = y1, y = y2, fill = PTS)) + geom_tile() + 
      scale_fill_gradientn(colours = colorsGray, values = seq(0,1,(1-0)/99), limits = c(0,1)) + 
      facet_grid(row~col) + 
      theme(plot.title = element_text(color="black", size=sizeV, hjust = 0.5),
        axis.title.x = element_text(color="black", size=sizeV),
        axis.title.y = element_text(color="black", size=sizeV),
        axis.text.x = element_text(color="black", size=sizeV, angle=0),
        axis.text.y = element_text(color="black", size=sizeV, angle=0),
        legend.text = element_text(color="black", size=sizeV/1.25, angle=0),
        legend.title = element_text(color="black", size=sizeV/1.25, angle=0),
        strip.text = element_text(color="black", size=sizeV, angle=0)) + 
      xlab("Years since IRS starts (Transient IRS)") + ylab("Years since IRS starts") +
      geom_text(size = sizeV/3, data  = dat_text,
        mapping = aes(x = Inf, y = Inf, label = label),
        hjust   = 2.25, vjust   = 2.5
      ) + 
      scale_x_continuous(breaks = round(times[seq(1, length(times), 2)])) +
      scale_y_continuous(breaks = round(times[seq(1, length(times), 2)])) +
      geom_hline(yintercept = c(0, IRSDur), col = "gray75", linetype = "dashed") +
      geom_vline(xintercept = c(0, IRSDur), col = "gray75", linetype = "dashed") +
      theme(strip.background = element_blank(), strip.text = element_blank())
    print(p)
  } else {
    p <- ggplot(dfAll, aes(x = y1, y = y2, fill = PTS)) + geom_tile() + 
      scale_fill_gradientn(colours = colorsGray, values = seq(0,1,(1-0)/99), limits = c(0,1)) + 
      facet_grid(row~col) + 
      theme(plot.title = element_text(color="black", size=sizeV, hjust = 0.5),
        axis.title.x = element_text(color="black", size=sizeV),
        axis.title.y = element_text(color="black", size=sizeV),
        axis.text.x = element_text(color="black", size=sizeV, angle=0),
        axis.text.y = element_text(color="black", size=sizeV, angle=0),
        legend.text = element_text(color="black", size=sizeV/1.25, angle=0),
        legend.title = element_text(color="black", size=sizeV/1.25, angle=0),
        strip.text = element_text(color="black", size=sizeV, angle=0)
      ) + xlab("Years since IRS starts (Transient IRS)") + ylab("Years since IRS starts") +
      scale_x_continuous(breaks = round(times[seq(1, length(times), 2)])) +
      scale_y_continuous(breaks = round(times[seq(1, length(times), 2)])) +
      geom_hline(yintercept = c(0, IRSDur), col = "gray75", linetype = "dashed") +
      geom_vline(xintercept = c(0, IRSDur), col = "gray75", linetype = "dashed") +
      theme(strip.background = element_blank(), strip.text = element_blank())
    print(p)
  }
  pAll[[i]] <- p
}

p1_rmL <- pAll[[1]] + rremove("xlab") + rremove("ylab") + rremove("legend") + ggtitle("Fast-rebound Regime")
p2_rmL <- pAll[[2]] + rremove("xlab") + rremove("ylab") + rremove("legend") + ggtitle("Approaching Transition Regime") 
p3_rmL <- pAll[[3]] + rremove("xlab") + rremove("ylab") + rremove("legend") + ggtitle("Transition Regime") 
p4_rmL <- pAll[[4]] + rremove("xlab") + rremove("ylab") + rremove("legend") + ggtitle("Slow-rebound Regime") 
p5_rmL <- pAll[[5]] + rremove("legend") + rremove("xlab") + rremove("ylab") + ggtitle("Slow-rebound Regime")
p6_legend <- get_legend(pAll[[1]])

saveDir <- "/project2/pascualmm/QZ/PhD/projects/intervention/natComRevision/plots/figures/main/Fig6/PTSHeatmap/"
if (!dir.exists(saveDir)) {
  dir.create(saveDir)
}
ggsave(paste0(saveDir, "transientIRS-heatmap-1.pdf"), p1_rmL, width = 10, height = 3.208333)
ggsave(paste0(saveDir, "transientIRS-heatmap-2.pdf"), p2_rmL, width = 10, height = 3.208333)
ggsave(paste0(saveDir, "transientIRS-heatmap-3.pdf"), p3_rmL, width = 10, height = 3.208333)
ggsave(paste0(saveDir, "transientIRS-heatmap-4.pdf"), p4_rmL, width = 10, height = 3.208333)
ggsave(paste0(saveDir, "transientIRS-heatmap-5.pdf"), p5_rmL, width = 10, height = 3.208333)
ggsave(paste0(saveDir, "transientIRS-heatmap-lg.pdf"), p6_legend, width = 6, height = 6)

