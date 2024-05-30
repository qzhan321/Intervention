rm(list = ls())
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(Hmisc)
})
scaleFUN <- function(x) sprintf("%.2f", x)
color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
run <- 4
readDir <- paste0("/project2/pascualmm/QZ/PhD/projects/intervention/files", run, "/actualRuns/geneFreqTime/")

IRSTypes <- c("2yIRS")
seasonality <- c("seasonal")
openness <- c("closed")
N <- 10000
n_threshold <- 0
T_YEAR <- 360
sizeV <- 24
preIRS <- 200
numsPlot <- c(1,8,10)
state <- "true"
prefix <- "sim"
pAll <- list()
rep <- 0
for (i in 1:length(seasonality)) {
  s <- seasonality[i]
  for (j in 1:length(openness)) {
    o <- openness[j]
    for (k in 1:length(IRSTypes)) {
      IRSType <- IRSTypes[k]
      label <- paste0(capitalize(s), " ", capitalize(o))
      for (m in 1:length(numsPlot)) {
        numPlot <- numsPlot[m]
        file <- paste0(readDir, s, "/", o, "/", state, "/", prefix, "_", numPlot, "_rep_", rep, "_", IRSType)
        load(file)
        if (numPlot == 1) {
          cy <- "black"
          gt <- "Fast-rebound Regime"
        } else if (numPlot == 8) {
          cy <- "white"
          gt <- "Transition Regime"
        } else if (numPlot == 10) {
          cy <- "white"
          gt <- "Slow-rebound Regime"
        }
        geneFreqTimeTb1 <- df12 %>% filter(num %in% numPlot, time > (preIRS - 3)*T_YEAR)
        geneFreqTimeTb2 <- geneFreqTimeTb1 %>% mutate(gene_id_new = ifelse(n >= n_threshold, gene_id, "low-freq"))
        geneFreqTimeTb3 <- geneFreqTimeTb2 %>% group_by(time, num, r, gene_id_new) %>% summarise(p = sum(p))
        geneFreqTimeTb4 <- geneFreqTimeTb3 %>% mutate(year = time/T_YEAR - preIRS)
        
        geneFreqTimeTbSupp1 <- geneFreqTimeTb4 %>% group_by(gene_id_new) %>% summarise(sump = sum(p)) %>% arrange(desc(sump))
        geneFreqTimeTbSupp1$levelId <- 1:nrow(geneFreqTimeTbSupp1)
        geneFreqTimeTbSupp1$levelId <- factor(geneFreqTimeTbSupp1$levelId, levels = 1:nrow(geneFreqTimeTbSupp1))
        
        geneFreqTimeTb5 <- geneFreqTimeTb4 %>% left_join(geneFreqTimeTbSupp1, by = "gene_id_new")
        geneFreqTimeTb6 <- geneFreqTimeTb5 %>% arrange(desc(sump)) %>% select(-sump)
        
        geneFreqTimeTb7 <- NULL
        times <- unique(geneFreqTimeTb6$time)
        geneidsAll <- unique(geneFreqTimeTb6$gene_id_new)
        for (n in 1:length(times)) {
          t <- times[n]
          temp <- geneFreqTimeTb6 %>% filter(time == t)
          geneidsTemp <- unique(temp$gene_id_new)
          geneidsDiff <- setdiff(geneidsAll, geneidsTemp)
          tempSupp <- as_tibble(data.frame("time" = t, "num" = numPlot, "r" = rep, "gene_id_new" = geneidsDiff, "p" = 0, year = unique(temp$year)))
          tempSuppLevel <- geneFreqTimeTbSupp1 %>% filter(gene_id_new %in% geneidsDiff) %>% select(-sump)
          tempSupp <- tempSupp %>% left_join(tempSuppLevel, by = "gene_id_new")
          tempAll <- bind_rows(temp, tempSupp)
          geneFreqTimeTb7 <- bind_rows(geneFreqTimeTb7, tempAll)
        }
        
        colVec = c(sample(color, length(unique(geneFreqTimeTb7$gene_id_new))-1, replace=T), "grey")
        
        p1<-ggplot(geneFreqTimeTb7, aes(year, p, color = levelId, fill = levelId)) + 
          geom_area() +
          labs(x="Years Since IRS Starts \n(Transient IRS)", y = "Frequency") +
          ggtitle(gt) +
          theme_bw() +
          theme(plot.title = element_text(color="black", size=sizeV, hjust = 0.5),
            axis.title.x = element_text(color="black", size=sizeV),
            axis.title.y = element_text(color=cy, size=sizeV),
            axis.text.x = element_text(color="black", size=sizeV, angle=0),
            axis.text.y = element_text(color=cy, size=sizeV, angle=0),
            legend.text = element_text(color="black", size=sizeV, angle=0),
            legend.title = element_text(color="black", size=sizeV, angle=0),
            strip.text = element_text(color="black", size=sizeV, angle=0),
            legend.position = "none") + 
          scale_color_manual(values = colVec) +
          scale_fill_manual(values = colVec) +
          scale_y_continuous(labels = scaleFUN) +
          scale_x_continuous(breaks = seq(-2, 10, 2))
        print(p1)
        index <- length(pAll)
        pAll[[index + 1]] <- p1
      }
    }
  }
}
p1_rmL <- pAll[[1]] 
p2_rmL <- pAll[[2]] 
p3_rmL <- pAll[[3]] 
saveDir <- "/project2/pascualmm/QZ/PhD/projects/intervention/natComRevision/plots/figures/main/Fig6/geneFreqTime/"
if (!dir.exists(saveDir)) {
  dir.create(saveDir)
}
ggsave(paste0(saveDir, "geneFreqTime-1.pdf"), p1_rmL, width = 5, height = 5)
ggsave(paste0(saveDir, "geneFreqTime-2.pdf"), p2_rmL, width = 5, height = 5)
ggsave(paste0(saveDir, "geneFreqTime-3.pdf"), p3_rmL, width = 5, height = 5)
