################################# 
## R script                    ##
## Project: TwinsRA            ##
## Beta inter intra.           ##
## Data: Shotgun metagenomics  ##
## Author: KB                  ##
## Date: 3/7/21                ##
#################################

### Load and save current R script ###

# Load R scripts
# load(file="/Users/KevinBu/Desktop/clemente_lab/Projects/twinsra/outputs/jobs20/beta_inter_intra_data.RData")
# load(file="/Users/KevinBu/Desktop/clemente_lab/Projects/twinsra/outputs/jobs20/beta_inter_intra_script.RData")

# Save R script
# Do this step prior to closing R
# save.image(file="/Users/KevinBu/Desktop/clemente_lab/Projects/twinsra/outputs/jobs20/beta_inter_intra_data.RData")

############################################################################
############################################################################
############################################################################

### Load libraries ###
library(reshape2)
library(phyloseq)
library(vegan)
library(ade4)
library(PMCMRplus)
library(ggplot2)
library(ggthemes)
library(ggrepel)
library(tidyverse)
library(dplyr)
library(tidyr)
library(scales)

############################################################################
############################################################################
############################################################################

### Statistics functions ###

# All plot statistics: mean, std deviation, median, min value, max value, 10%ile, 25%ile, 75%ile, 90%ile
stats.all = function(x) {
  mean <- mean(x)
  stddev <- sd(x)
  median <- median(x)
  val_min <- min(x)
  val_max <- max(x)
  per10 <- as.numeric(quantile(x, prob = c(0.10)))
  per25 <- as.numeric(quantile(x, prob = c(0.25)))
  per75 <- as.numeric(quantile(x, prob = c(0.75)))
  per90 <- as.numeric(quantile(x, prob = c(0.90)))
  return(c(mean = mean, sd = stddev, median = median, 
           val_min = val_min, val_max = val_max, 
           per10 = per10, per25 = per25, per75 = per75,  per90 = per90))
}

# Boxplot statistics: median, 25%ile, 75%ile
stats.boxplot <- function(x) {
  m <- median(x)
  per25 <- as.numeric(quantile(x, prob = c(0.25)))
  per75 <- as.numeric(quantile(x, prob = c(0.75)))
  return(c(y = m, ymin = per25, ymax = per75))
}

# Whiskers statistics: median, 10th percentile, 90th percentile
stats.whiskers = function(x) {
  m <- median(x)
  per10 <- as.numeric(quantile(x, prob = c(0.10)))
  per90 <- as.numeric(quantile(x, prob = c(0.90)))
  return(c(y = m, ymin = per10, ymax = per90))
}

# Outliers
min.outlier <- function(x) {
  subset(x, quantile(x, prob = c(0.10)) > x)
}

max.outlier <- function(x) {
  subset(x, quantile(x, prob = c(0.90)) < x)
}

############################################################################
############################################################################
############################################################################

# set working dir
dir = "/Users/KevinBu/Desktop/clemente_lab/Projects/twinsra/outputs/jobs02/"

# read in df
df = read.table('/Users/KevinBu/Desktop/clemente_lab/Projects/twinsra/outputs/jobs02/intra_beta_dist.tsv', 
                      sep = '\t', header = TRUE, row.names = 1, check.names = FALSE,
                      na.strings = "NA")

# rename factors
df$category[df$category == 'intra_UA_only'] <- "Intra-Unaffected"
df$category[df$category == 'intra_RA_only'] <- "Intra-RA"
df$category[df$category == 'intra_twin_pair'] <- "Intra-Twin"

# fix order of plotting
df$category <- factor(df$category, levels = c("Intra-Twin", "Intra-RA", "Intra-Unaffected")) 

# background theme
bkg <- theme_bw() +
  theme(axis.text.x = element_text(size = 8, face = "bold", color = "black")) +
  theme(axis.text.y = element_text(size = 10, color = "black")) +
  theme(axis.title.y = element_text(size = 12, color = "black", face = "bold")) +
  theme(axis.title.y = element_text(margin = unit(c(0, 8, 0, 0), "mm"))) +
  theme(legend.text = element_text(size = 8, color = "black")) +
  theme(legend.title = element_text(size = 10, face = "bold", color = "black"))

# choose colors
col1 <- c("#929aab", "#ce2525")
col2 <- c("#CD3414", "#0074e4", "#8f8787")
col2 <- c("#CD3414", "#8f8787", "#0074e4")
col2 <- c("#0074e4", "#CD3414", "#8f8787")

# for every plasma analyte in the table, calculate the change in abundance between siblings
for (i in 1:1) {
  
  # analyte of interest
  a <- colnames(df)[i+1]
  
  # create filenames
  filename_table = paste(a, "all_table.csv", sep = "_")
  filename_box.plot = paste("Beta", "all_box.plot.pdf", sep = "_")  
  
  # plot boxplot
  p <- ggplot(data = df, aes(x = category, y = distance, fill = category)) +
    stat_summary(fun.data = stats.whiskers, geom = "errorbar", 
                 color = "black", size = 0.8, width = 0.3) +
    stat_summary(fun.data = stats.boxplot, geom = "crossbar", 
                 color = "black", size = 0.5, width = 0.5) +
    geom_jitter(width = 0.1, size = 1.5) +
    scale_x_discrete(labels = c("Intra-Twin", "Intra-RA", "Intra-Unaffected")) +
    scale_fill_manual(values = col2) +      
    xlab(NULL) +
    ylab('Bray Curtis Difference') +
    bkg
  
  fpb = paste(dir, filename_box.plot, sep = "")
  pdf(file = fpb, height = 4, width = 6)
  plot(p)
  dev.off()
  
}
