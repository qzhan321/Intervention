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
  library(RColorBrewer)
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
sizeV <- 31
pAll <- list()
state <- "true"
prefix <- "sim"
for (i in 1:length(seasonality)) {
  s <- seasonality[i]
  for (j in 1:length(openness)) {
    o <- openness[j]
    for (k in 1:length(IRSTypes)) {
      IRSType <- IRSTypes[k]
      label <- paste0(capitalize(s), " ", capitalize(o))
      preIRS <- 200
      numsPlot <- c(1,8,10)
      rep <- 0
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
        } else {
          cy <- "white"
          gt <- "Slow-rebound Regime"
        }
        geneFreqTimeTb <- df12 %>% filter(num %in% numPlot, time > (preIRS - 3)*T_YEAR)
        geneFreqTimeTb2 <- geneFreqTimeTb %>% mutate(gene_id_new = ifelse(n >= n_threshold, gene_id, "low-freq"))
        
        geneFreqTimeTb3 <- geneFreqTimeTb2 
        geneFreqTimeTb4 <- geneFreqTimeTb3 %>% group_by(time, r, gene_id_new) %>% summarise(p = sum(p))
        geneFreqTimeTb5 <- geneFreqTimeTb4 %>% mutate(year = time/T_YEAR - preIRS)
        
        gfsupp <- geneFreqTimeTb5 %>% group_by(gene_id_new) %>% summarise(sump = sum(p)) %>%
          arrange(desc(sump))
        gfsupp$levelId <- 1:nrow(gfsupp)
        gfsupp$levelId <- factor(gfsupp$levelId, levels = 1:nrow(gfsupp))
        
        geneFreqTimeTb6 <- geneFreqTimeTb5 %>% left_join(gfsupp, by = "gene_id_new")
        geneFreqTimeTb7 <- geneFreqTimeTb6 %>% arrange(desc(sump)) %>% select(-sump)
        
        geneFreqTimeTb8 <- NULL
        times <- unique(geneFreqTimeTb7$time)
        geneidsAll <- unique(geneFreqTimeTb7$gene_id_new)
        for (n in 1:length(times)) {
          t <- times[n]
          gft <- geneFreqTimeTb7 %>% filter(time == t)
          geneidsTemp <- unique(gft$gene_id_new)
          geneidsDiff <- setdiff(geneidsAll, geneidsTemp)
          gft_supp <- as_tibble(data.frame("time" = t, "r" = rep, "gene_id_new" = geneidsDiff, "p" = 0, year = unique(gft$year)))
          gft_supp_l <- gfsupp %>% filter(gene_id_new %in% geneidsDiff) %>% select(-sump)
          gft_supp_all <- gft_supp %>% left_join(gft_supp_l, by = "gene_id_new")
          gft_all <- rbind(gft, gft_supp_all)
          geneFreqTimeTb8 <- rbind(geneFreqTimeTb8, gft_all)
        }
        
        colVec = c("grey", sample(color, length(unique(geneFreqTimeTb8$gene_id_new))-1, replace=T))
        
        p1<-ggplot(geneFreqTimeTb8, aes(year, p, color = levelId, fill = levelId)) + 
          geom_area() +
          labs(x="Years Since IRS Starts \n(Transient IRS)", y = "Frequency") +
          ggtitle(gt) +
          theme_bw() +
          theme(
            plot.title = element_text(color="black", size=sizeV, hjust = 0.5),
            axis.title.x = element_text(color="black", size=sizeV),
            axis.title.y = element_text(color=cy, size=sizeV),
            axis.text.x = element_text(color="black", 
                                       size=sizeV, angle=0),
            axis.text.y = element_text(color=cy, 
                                       size=sizeV, angle=0),
            legend.text = element_text(color="black", 
                                       size=sizeV, angle=0),
            legend.title = element_text(color="black", 
                                        size=sizeV, angle=0),
            strip.text = element_text(color="black", 
                                      size=sizeV, angle=0),
            # plot.margin = unit(c(0.35,0.035,0.035,0.2), "cm"),
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
saveDir <- paste0("/project2/pascualmm/QZ/PhD/projects/intervention/writings/simulation", run, "/plots/geneFreqTime/")
dir.create(saveDir)
ggsave(paste0(saveDir, "geneFreqTime-1.pdf"), p1_rmL, width = 6, height = 6)
ggsave(paste0(saveDir, "geneFreqTime-2.pdf"), p2_rmL, width = 6, height = 6)
ggsave(paste0(saveDir, "geneFreqTime-3.pdf"), p3_rmL, width = 6, height = 6)
