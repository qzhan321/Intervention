rm(list=ls())
suppressPackageStartupMessages({
  library(dplyr)
  library(optparse)
})

option_list <- list(
  make_option(c("-i", "--inputFile"), type = "character", default = NULL,
              help = "the path and the name of the input file containing host IDs and the number of their non-UPs A DBLa types"),
  make_option(c("-m", "--maxMOI"), type = "integer", default = 20,
              help = "the maximum possible value for MOI; can be some rough approximation based on the empirical dataset"),
  make_option(c("-r", "--repertoireSizeDist"), type = "character", default = NULL,
              help = "the path and the name of the repertoire size distribution file"),
  make_option(c("-t", "--typeOfRepertoireSizeDist"), type = "character", default = "count",
              help = "the type of the repertoire size distribution file, either the original count data of non-UPs A DBLa types for monoclonal infections, or its (smoothed) probability version"),
  make_option(c("-p", "--prior"), type = "character", default = "uniform",
              help = "user-specified prior type for MOI estimation [default= %default]; either a uniform or a negative binomial (negBinom)"),
  make_option(c("-s", "--params"), type = "character", default = "medium",
              help = "for a negative binomial prior, providing a parameter range [default= %default]; three options, low or medium or high, corresponding to a negative binomial with a low (~1), medium (~4), and high (~7) mean value"),
  make_option(c("-v", "--verbose"), type = "logical", default = TRUE,
              help = "whether output the prior distribution for MOI estimation or not [default= %default]"),
  make_option(c("-a", "--aggregate"), type = "character", default = "pool",
              help = "obtain the MOI distribution at the population level from individual MOI estimates; pool or mixtureDist: either pooling the maximum a posteriori MOI estimate of each individual host, or using mixture distribution, i.e., a weighted summation of all individual MOI distributions"),
  make_option(c("-o", "--output"), type = "character", default = "./out.RData",
              help = "the path and the name of the output file [default= %default]")
)
opt = parse_args(OptionParser(option_list = option_list))
print(opt)


# read input file
inputFile <- read.csv(opt$input, header = T, row.names = NULL)


# MOI priors
simPrior <- function(prior, size = 5, prob = 0.55, maxMOI = 20) {
  if (prior == "uniform") {
    MOI_priors <- rep(1/maxMOI, maxMOI)
    names(MOI_priors) <- as.character(1:maxMOI)
  } else if (prior == "negBinom") {
    set.seed(0)
    MOI_priors_temp <- rnbinom(10000000, size = size, prob = prob)
    MOI_priors_temp <- MOI_priors_temp[MOI_priors_temp >= 1 & MOI_priors_temp <= maxMOI]
    MOI_priors_temp_df <- data.frame("MOI" = MOI_priors_temp) %>% group_by(MOI) %>% 
      summarise(n = n()) %>% mutate(prop = n/sum(n))
    
    prior_baseline <- 0.00001 # a small non-zero baseline probability for these not likely MOI values given the specified parameters of the negative binomial
    MOI_priors <- rep(prior_baseline, maxMOI)
    names(MOI_priors) <- as.character(1:maxMOI)
    MOIValues <- as.character(MOI_priors_temp_df$MOI)
    MOIValuesProps <- MOI_priors_temp_df$prop
    MOI_priors[MOIValues] <- MOIValuesProps + prior_baseline
    MOI_priors <- MOI_priors/sum(MOI_priors)
  }
  return(MOI_priors)
} 
maxMOI <- opt$maxMOI
if (opt$prior == "uniform") {
  MOI_priors <- simPrior(prior = opt$prior, maxMOI = maxMOI)
} else if (opt$prior == "negBinom") {
  if (opt$params == "low") {
    set.seed(0)
    MOI_priors <- simPrior(prior = opt$prior, size = 2, prob = 0.75, maxMOI = maxMOI)
  } else if (opt$params == "medium") {
    set.seed(1)
    MOI_priors <- simPrior(prior = opt$prior, size = 5, prob = 0.55, maxMOI = maxMOI)
  } else if (opt$params == "high") {
    set.seed(2)
    MOI_priors <- simPrior(prior = opt$prior, size = 4, prob = 0.375, maxMOI = maxMOI)
  }
}


# serial convolution to derive the probability of sequencing and sampling s types given any possible value of MOI
repertoireSizeDistFile <- read.csv(opt$repertoireSizeDist, header = T, row.names = NULL)
repSizeLow <- min(repertoireSizeDistFile$DBLa_upsBC_rep_size) # the lowest number of non-upsA DBLa types sequenced in a single infection, from the user-supplied empirical repertoire size distribution
repSizeHigh <- max(repertoireSizeDistFile$DBLa_upsBC_rep_size) # the highest number of non-upsA DBLa types sequenced in a single infection, from the user-supplied empirical repertoire size distribution
if (opt$typeOfRepertoireSizeDist == "count") {
  repertoireSizeDistFile <- repertoireSizeDistFile %>% mutate(p = n/sum(n)) %>% select(-n)
} 
# for certain repertoire size, the probability may be 0; replace 0 with some small baseline probability
fillInMissingOrZero <- function(repertoireSizeDistFile) {
  if (!all(repertoireSizeDistFile$p > 0)) {
    repSizes_p_zero <- repertoireSizeDistFile$DBLa_upsBC_rep_size[repertoireSizeDistFile$p == 0]
    for (repSize_p_zero in repSizes_p_zero) {
      distances <- abs(repertoireSizeDistFile$DBLa_upsBC_rep_size-repSize_p_zero)
      neighbors_index_list <- sort(distances, index.return=TRUE)
      p_zero_index <- neighbors_index_list$ix[1]
      neighbors_index <- neighbors_index_list$ix[2:3]
      repertoireSizeDistFile$p[p_zero_index] <- mean(repertoireSizeDistFile$p[neighbors_index])
    }
  }  
  if (!all(seq(repSizeLow, repSizeHigh, 1) %in% repertoireSizeDistFile$DBLa_upsBC_rep_size)) {
    repSizes_missing <- setdiff(seq(repSizeLow, repSizeHigh, 1), repertoireSizeDistFile$DBLa_upsBC_rep_size)
    for (repSize_missing in repSizes_missing) {
      distances <- abs(repertoireSizeDistFile$DBLa_upsBC_rep_size-repSize_missing)
      neighbors_index <- sort(distances, index.return=TRUE)
      neighbors_index <- neighbors_index$ix[1:2]
      missing_neighbors <- repertoireSizeDistFile %>% filter(DBLa_upsBC_rep_size %in% repertoireSizeDistFile$DBLa_upsBC_rep_size[neighbors_index]) 
      repSize_missing_dist <- data.frame(DBLa_upsBC_rep_size = repSize_missing, p = mean(missing_neighbors$p))
      repertoireSizeDistFile <- rbind(repertoireSizeDistFile, repSize_missing_dist) 
    }
  }
  repertoireSizeDistFile <- repertoireSizeDistFile %>% arrange(DBLa_upsBC_rep_size) %>% mutate(p = p/sum(p))
  return(repertoireSizeDistFile)
}
repertoireSizeDistFile_no_missing <- fillInMissingOrZero(repertoireSizeDistFile)

prob_s_givenMOI <- function(repSizeLow, repSizeHigh, maxMOI, repSizeDist) {
  prob_s_givenMOI <- list()
  
  p <- repSizeDist$p
  names(p) <- as.character(repSizeDist$DBLa_upsBC_rep_size)
  prob_s_givenMOI[[1]] <- p
  
  for (MOI in 2:maxMOI) {
    convolution1 <- prob_s_givenMOI[[1]]
    convolution2 <- prob_s_givenMOI[[MOI-1]]
    
    prob_s_givenMOI_single <- rep(NA, repSizeHigh*MOI - repSizeLow*MOI + 1)
    names(prob_s_givenMOI_single) <- as.character(seq(repSizeLow*MOI, repSizeHigh*MOI))
    
    for (isoSize in (repSizeLow*MOI):(repSizeHigh*MOI)) {
      probAll <- 0
      for (convolution1Size in repSizeLow:repSizeHigh) {
        convolution1Prob = convolution1[as.character(convolution1Size)]
        convolution2Prob = convolution2[as.character(isoSize - convolution1Size)]
        prob <- convolution1Prob * convolution2Prob
        if(!is.na(prob)) {
          names(prob) <- NULL
          probAll <- probAll + prob
        }
      }
      prob_s_givenMOI_single[as.character(isoSize)] <- probAll
    }
    prob_s_givenMOI[[MOI]] <- prob_s_givenMOI_single
  }
  names(prob_s_givenMOI) <- as.character(1:maxMOI)
  return(prob_s_givenMOI)
}
prob_s_givenMOI <- prob_s_givenMOI(repSizeLow, repSizeHigh, maxMOI, repertoireSizeDistFile_no_missing)


# estimate MOI
MOIest <- function(prob_s_givenMOI, inputFile, maxMOI, repSizeLow, repSizeHigh) {
  prob_s_monoclonal <- prob_s_givenMOI[[1]]
  isoSizes <- inputFile$DBLa_upsBC_rep_size
  hosts <- inputFile$HostID
  
  MOIAll <- NULL
  for (i in 1:length(isoSizes)) {
    isoSize <- isoSizes[i]
    host <- hosts[i]
    
    if (isoSize < repSizeLow) {
      warning("The number of non-upsA DBLa type is fewer than the lower number presented in your repertoire size distribution. Automatically assign 1 to be the MOI value.")
      prob_MOI_given_isoSize <- c(1.0, rep(0.0, maxMOI - 1))
    } else if (isoSize > repSizeHigh*maxMOI) {
      warning("The number of non-upsA DBLa type is greater than the maximum possible one, i.e., the maximum repertoire size * maxMOI. Automatically assign maxMOI to be the MOI value.")
      prob_MOI_given_isoSize <- c(rep(0.0, maxMOI - 1), 1.0)
    } else if (isoSize >= repSizeLow & isoSize <= repSizeHigh*maxMOI) {
      denominator <- 0
      numerators <- rep(0, maxMOI)
      for (MOI in 1:maxMOI) {
        prob_s_givenMOI_single <- prob_s_givenMOI[[as.character(MOI)]][as.character(isoSize)]
        MOI_prior_single <- MOI_priors[as.character(MOI)]
        denominator_singleMOI <- prob_s_givenMOI_single * MOI_prior_single
        if (!is.na(denominator_singleMOI)) {
          names(denominator_singleMOI) <- NULL
          numerators[MOI] <- denominator_singleMOI
          denominator <- denominator + denominator_singleMOI
        }
      }
      stopifnot(denominator != 0)
      prob_MOI_given_isoSize <- numerators/denominator
    }
    if (opt$aggregate == "pool") {
      MOI_single <- data.frame("HostID" = host, "DBLa_upsBC_rep_size" = isoSize, "MOI" = c(1:maxMOI)[which.max(prob_MOI_given_isoSize)], "Prob" = max(prob_MOI_given_isoSize))
    } else if (opt$aggregate == "mixtureDist") {
      MOI_single <- data.frame("HostID" = host, "DBLa_upsBC_rep_size" = isoSize, "MOI" = 1:maxMOI, "Prob" = prob_MOI_given_isoSize)
    }
    MOIAll <- rbind(MOIAll, MOI_single)
  }
  return(MOIAll)
}

MOIestimation <- MOIest(prob_s_givenMOI, inputFile, maxMOI, repSizeLow, repSizeHigh)


# aggregate MOI estimates across individuals to get population level MOI distribution.
if (opt$aggregate == "pool") {
  MOIAllPop <- MOIestimation %>% group_by(MOI) %>% summarise("n" = n()) %>% mutate("n_total" = sum(n), "Prob" = n/n_total) %>% select(MOI, Prob)
} else if (opt$aggregate == "mixtureDist") {
  MOIAllPop <- MOIestimation %>% group_by(MOI) %>% summarise("n" = sum(Prob)) %>% mutate("n_total" = sum(n), "Prob" = n/n_total) %>% select(MOI, Prob)
} else if (opt$aggregate == "mixtureDist") {
}
if (opt$verbose == TRUE) {
  outputList <- list("prior" = MOI_priors,
                     "indLevelMOI" = MOIestimation,
                     "popLevelMOI" = MOIAllPop)
} else {
  outputList <- list("indLevelMOI" = MOIestimation,
                     "popLevelMOI" = MOIAllPop)
}


# save output 
save(outputList, file = opt$output)
