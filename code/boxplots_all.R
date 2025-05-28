
### Load libraries ###
library(reshape2)
# library(phyloseq)
library(vegan)
# library(ade4)
# library(PMCMR)
# library(PMCMRplus)
library(ggplot2)
library(ggthemes)
library(ggrepel)
library(tidyverse)
library(dplyr)
library(tidyr)
library(scales)
library(ggpubr)

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

# Function to add asterisks for significance levels
add_significance <- function(p_value) {
  if (p_value < 0.0001) {
    return("****")
  } else if (p_value < 0.001) {
    return("***")
  } else if (p_value < 0.01) {
    return("**")
  } else if (p_value < 0.05) {
    return("*")
  } else {
    return("")
  }
}


############################################################################
############################################################################
############################################################################

### all data
paths = c('0_taxa')
        #  ('8_metabolites', 
        #  '10_acpa_fecal', '11_acpa_plasma', '12_olink', '13_metabolon', '15_taxa', '16_path', '17_ko',
        #  'L2_taxa', 'L3_taxa', 'L4_taxa', 'L5_taxa', 'L6_taxa', 'L7_taxa')
offsets = c(2)
      #   c(3, 3, 3, 3, 2, 2, 2, 2,
      #      2, 2, 2, 2, 2, 2)
nvars = c(1320)
# c(11, 116, 116, 86, 10, 690, 327, 429,
#          14, 119, 134, 160, 425, 690)
sizes = c(16)
# sizes = c(24, 16, 16, 24, 16, 24, 16, 16,
#          16, 16, 16, 16, 16, 16)
units = c(" abundance")
#units = c(" abundance (nmol/mg)", " abundance (MFI)", " abundance (MFI)", " level (NPX)", " abundance (ng/ml)", " abundance", " abundance",
#          " abundance", " abundance", " abundance", " abundance", " abundance", " abundance")
#data_fp = c("metabolites.txt", "acpa_fecal.txt", "acpa_plasma.txt", "olink.tsv", "metabolon.tsv", "taxa.tsv", "path.tsv", 'ko.tsv',
#            'L2_taxa.tsv','L3_taxa.tsv','L4_taxa.tsv', 'L5_taxa.tsv','L6_taxa.tsv','L7_taxa.tsv')
data_fp = c('taxa.tsv')

### statistics ###
# put a 7 here instead of 1 for pathways only
#for (k in 1:length(paths)){
for (k in length(paths)){
  path = paths[k]
  offset = offsets[k]
  nvar = nvars[k]
  unit = units[k]
  data = data_fp[k]
  size = sizes[k]
  
  # directory
  dir = paste0("/Users/KevinBu/Desktop/clemente_lab/Projects/twinspsd/outputs/jobs03/",path,"/")
  
  # read in data table
  d <- read.table(file = paste(dir, data, sep = ""),
                  header = TRUE, row.names = 1, sep = "\t", check.names = FALSE, # maybe set to true next time lol
                  na.strings = "NA")
  
  # create dataframes
  #d$Diagnosis <- sub("HC", "HC", d$Diagnosis)
  #d$Diagnosis <- sub("2_PsO", "PsO", d$Diagnosis)
  d.clean <- as.data.frame(d)
  d.all <- as.data.frame(d)
  
  #d <- as.data.frame(d)
  #if (path == '8_metabolites'){
  #  d.all <- d[!d$twinpair == "sib_04",] # remove 506 as no data for 507
  #}
  #d.clean <- d.all[!d.all$twinpair == "sib_10",] # remove 182, 183
  #d.clean <- d.clean[!d.clean$twinpair == "sib_09",] # remove 516, 517
  d.clean <- d.all
  
  # create tables for storing wilcoxon and ttest results
  stats.table.all <- matrix(data = NA, nrow = nvar, ncol = 5)
  colnames(stats.table.all) <- c("feature", "wilcoxon", "ttest", "wilcox stat", "ttest stat")
  
  stats.table.clean <- matrix(data = NA, nrow = nvar, ncol = 5)
  colnames(stats.table.clean) <- c("feature", "wilcoxon", "ttest", "wilcox stat", "ttest stat")
  
  # calculate wilcoxon and ttest between HC/PsO siblings for each metabolite
  for (i in 1:nvar) {
    stats.table.all[i,1] <- colnames(d.all)[i+offset]
    stats.table.clean[i,1] <- colnames(d.clean)[i+offset]
    
    #stats.table.all[i,2] <- wilcox.test(d.all[,i+offset] ~ Diagnosis, data = d.all, paired = TRUE)$p.value
    #stats.table.clean[i,2] <- wilcox.test(d.clean[,i+offset] ~ Diagnosis, data = d.clean, paired = TRUE)$p.value

    #stats.table.all[i,4] <- wilcox.test(d.all[,i+offset] ~ Diagnosis, data = d.all, paired = TRUE)$statistic
    #stats.table.clean[i,4] <- wilcox.test(d.clean[,i+offset] ~ Diagnosis, data = d.clean, paired = TRUE)$statistic
    
    #stats.table.all[i,3] <- t.test(d.all[,i+offset] ~ Diagnosis, data = d.all, paired = TRUE)$p.value
    #stats.table.clean[i,3] <- t.test(d.clean[,i+offset] ~ Diagnosis, data = d.clean, paired = TRUE)$p.value
    
    #stats.table.all[i,5] <- t.test(d.all[,i+offset] ~ Diagnosis, data = d.all, paired = TRUE)$statistic
    #stats.table.clean[i,5] <- t.test(d.clean[,i+offset] ~ Diagnosis, data = d.clean, paired = TRUE)$statistic
  
    stats.table.all[i,2] <- wilcox.test(x=d.all[d.all$Diagnosis == 'PsO',][,i+offset], 
                                        y=d.all[d.all$Diagnosis == 'HC',][,i+offset],
                                        paired=TRUE)$p.value
    stats.table.clean[i,2] <- wilcox.test(x=d.clean[d.clean$Diagnosis == 'PsO',][,i+offset], 
                                          y=d.clean[d.clean$Diagnosis == 'HC',][,i+offset],
                                          paired=TRUE)$p.value
    
    stats.table.all[i,4] <- wilcox.test(x=d.all[d.all$Diagnosis == 'PsO',][,i+offset], 
                                        y=d.all[d.all$Diagnosis == 'HC',][,i+offset],
                                        paired=TRUE)$statistic
    stats.table.clean[i,4] <- wilcox.test(x=d.clean[d.clean$Diagnosis == 'PsO',][,i+offset], 
                                          y=d.clean[d.clean$Diagnosis == 'HC',][,i+offset],
                                          paired=TRUE)$statistic
    
    stats.table.all[i,3] <- t.test(x=d.all[d.all$Diagnosis == 'PsO',][,i+offset], 
                                   y=d.all[d.all$Diagnosis == 'HC',][,i+offset],
                                   paired=TRUE)$p.value
    stats.table.clean[i,3] <- t.test(x=d.clean[d.clean$Diagnosis == 'PsO',][,i+offset], 
                                     y=d.clean[d.clean$Diagnosis == 'HC',][,i+offset],
                                     paired=TRUE)$p.value
    
    stats.table.all[i,5] <- t.test(x=d.all[d.all$Diagnosis == 'PsO',][,i+offset], 
                                   y=d.all[d.all$Diagnosis == 'HC',][,i+offset],
                                   paired=TRUE)$statistic
    stats.table.clean[i,5] <- t.test(x=d.clean[d.clean$Diagnosis == 'PsO',][,i+offset], 
                                     y=d.clean[d.clean$Diagnosis == 'HC',][,i+offset],
                                     paired=TRUE)$statistic
      
  }
  
  # Calculate FDR-adjusted p-values using p.adjust
  v = c('0_taxa')#, '16_path', '17_ko', 'L2_taxa', 'L3_taxa', 'L4_taxa', 'L5_taxa', 'L6_taxa')
  
  #   if (path == '15_taxa'){
  if (path %in% v){
      adjusted_p_values <- p.adjust(stats.table.clean[,2], method = "fdr")
  }
  else {
    adjusted_p_values <- p.adjust(stats.table.clean[,3], method = "fdr")
  }
  # Add adjusted p-values as a new column to the matrix
  stats.table.clean <- cbind(stats.table.clean, fdr_adjusted = adjusted_p_values)
  
  # save
  ft.all = paste(dir, path, "_stats_all.csv", sep = "")
  write.csv(file = ft.all, stats.table.all)
  
  ft.clean = paste(dir, path, "_stats_clean.csv", sep = "")
  write.csv(file = ft.clean, stats.table.clean)
  
  ############################################################################
  ############################################################################
  ############################################################################
  
  ### plots - clean ###
  
  # background theme
  bkg <- theme_bw() +
    theme(axis.text.x = element_text(size = 18, face = "bold", color = "black")) +
    theme(axis.text.y = element_text(size = 18, color = "black")) +
    theme(axis.title.y = element_text(size = size, color = "black", face = "bold")) +
    theme(axis.title.y = element_text(margin = unit(c(0, 8, 0, 0), "mm"))) +
    theme(legend.text = element_text(size = 18, color = "black")) +
    theme(legend.title = element_text(size = 24, face = "bold", color = "black"))
  
  
  # directory for storing files
  dir = paste0("/Users/KevinBu/Desktop/clemente_lab/Projects/twinspsd/outputs/jobs03/",path,"/clean/")
  
  # colors
  col1 <- c("#929aab", "#ce2525")
  col1 <- c("white","#5DB54B")#,"#FB9A99")
  
  col2 <- c("#f3a333", "#0074e4", "#8f8787")
  
  #col3 <- c("#0074e4")
  #col4 <- c("#f3a333", "#8f8787")
  
  #col3 <- c("white","#FB9A99") # PsA only
  #col4 <- c('white','purple') # PsD
  
  
  # line types
  line1 <- c("solid", "dashed", "dotted")
  line2 <- c("dashed")
  line3 <- c("solid", "dotted")  
  # for every feature in the table, calculate the change in Level between siblings
  for (i in 1:nvar) {
    
    # metabolite of interest
    m <- colnames(d.clean)[i+offset]
    
    # create filenames
    filename_table = paste(m, "clean_table.csv", sep = "_")
    filename_box.plot = paste(m, "clean_box.plot.pdf", sep = "_")  
    filename_line.plot = paste(m, "clean_line.plot.pdf", sep = "_")  
    
    # create categories by which to color lineplot: 1) Increase, 2) Decrease, 3) No change
    # subset and spread dataset into diagnosis columns
    # then calculate delta relative Level between PsO and HC siblings
    d.div <- d.clean %>%
      subset(select = c("twinpair", "Diagnosis", m)) %>%
      spread(key = "Diagnosis", value = m) %>%
      mutate(Level_sib.pair = (get("PsO") - get("HC"))) %>%
      mutate(
        Change.type_sib.pair = case_when(
          sign(Level_sib.pair) == 1 ~ "1_up",
          sign(Level_sib.pair) == -1 ~ "2_down",
          TRUE ~ "3_no.change"
        )
      ) 
    
    # merge with original dataset
    d.final <- gather(d.div, "HC", "PsO", key = "Diagnosis", value = "Abundance")
    d.final$Diagnosis <- factor(d.final$Diagnosis, levels = c("HC", "PsO")) # fix order of plotting
    
    # save dataset
    ft = paste(dir, filename_table, sep = "")
    write.csv(d.final, file = ft)
    
    # plot boxplot
    p <- ggplot(data = d.final, aes(x = Diagnosis, y = Abundance, fill = Diagnosis)) +
      stat_summary(fun.data = stats.whiskers, geom = "errorbar", 
                   color = "black", size = 0.8, width = 0.3) +
      stat_summary(fun.data = stats.boxplot, geom = "crossbar", 
                   color = "black", size = 0.5, width = 0.5) +
      geom_jitter(width = 0.1, size = 1.5) +
      scale_x_discrete(labels = c("HC", "PsO")) +
      scale_fill_manual(values = col1) +      
      xlab(NULL) +
      ylab(paste0(m, unit)) +
      bkg
    
    fpb = paste(dir, filename_box.plot, sep = "")
    pdf(file = fpb, height = 4.5, width = 5)
    plot(p)
    dev.off()
    
    # plot line plot
    if(m == "Octanoate" | m == "Propionate" | m == "Acetate" | m == "Butyrate") {
      p <- ggplot(data = d.final, aes(x = Diagnosis, y = Abundance, fill = Diagnosis)) +
        stat_summary(fun.data = stats.whiskers, geom = "errorbar", 
                     color = "black", size = 0.8, width = 0.3) +
        stat_summary(fun.data = stats.boxplot, geom = "crossbar",
                     color = "black", size = 0.5, width = 0.5) +
        geom_point(shape = 20, size = 3, color = "black") +
        geom_line(data = subset(d.final), aes(group = twinpair, color = Change.type_sib.pair, linetype = Change.type_sib.pair), size = 0.5) +
        geom_text_repel(data = subset(d.final, Diagnosis == "HC"), aes(label = twinpair), 
                        nudge_x = -0.5, size = 2, color = "black", segment.alpha = 0.3) +
        scale_x_discrete(labels = c("HC", "PsO")) +
        scale_fill_manual(values = col1) +
        scale_linetype_manual(values = line2, name = "Change", labels = c("Decrease")) +
        scale_color_manual(values = col3, name = "Change", labels = c("Decrease")) +      
        xlab(NULL) +
        ylab(paste0(m, unit)) +
        bkg
      } else if (m == 'Bilophila wadsworthia') {
          p <- ggplot(data = d.final, aes(x = Diagnosis, y = Abundance, fill = Diagnosis)) +
            stat_summary(fun.data = stats.whiskers, geom = "errorbar", 
                         color = "black", size = 0.8, width = 0.3) +
            stat_summary(fun.data = stats.boxplot, geom = "crossbar",
                         color = "black", size = 0.5, width = 0.5) +
            geom_point(shape = 20, size = 3, color = "black") +
            geom_line(data = subset(d.final), aes(group = twinpair, color = Change.type_sib.pair, linetype = Change.type_sib.pair), size = 0.5) +
            geom_text_repel(data = subset(d.final, Diagnosis == "HC"), aes(label = twinpair), 
                            nudge_x = -0.5, size = 2, color = "black", segment.alpha = 0.3) +
            scale_x_discrete(labels = c("HC", "PsO")) +
            scale_fill_manual(values = col1) +
            scale_linetype_manual(values = line3, name = "Change", labels = c("Increase", "No change")) +
            scale_color_manual(values = col4, name = "Change", labels = c("Increase", "No change")) +      
            xlab(NULL) +
            ylab(paste0(m, unit)) +
            bkg
        
      } else {
      p <- ggplot(data = d.final, aes(x = Diagnosis, y = Abundance, fill = Diagnosis)) +
      stat_summary(fun.data = stats.whiskers, geom = "errorbar", 
                   color = "black", size = 0.8, width = 0.3) +
      stat_summary(fun.data = stats.boxplot, geom = "crossbar",
                   color = "black", size = 0.5, width = 0.5) +
      geom_point(shape = 20, size = 3, color = "black") +
      geom_line(data = subset(d.final), aes(group = twinpair, color = Change.type_sib.pair, linetype = Change.type_sib.pair), size = 0.5) +
      geom_text_repel(data = subset(d.final, Diagnosis == "HC"), aes(label = twinpair), 
                      nudge_x = -0.5, size = 2, color = "black", segment.alpha = 0.3) +
      scale_x_discrete(labels = c("HC", "PsO")) +
      scale_fill_manual(values = col1) +
      scale_linetype_manual(values = line1, name = "Change", labels = c("Increase", "Decrease", "No change")) +
      scale_color_manual(values = col2, name = "Change", labels = c("Increase", "Decrease", "No change")) +      
      xlab(NULL) +
      ylab(paste0(m, unit)) +
      bkg
      }
    
    # add significance
    #pv = as.numeric(stats.table.clean[i,2]) # i,3 for ttest
    #if (pv < 0.05){
    #  p <- p + geom_signif(comparisons = list(c("HC", "PsO")), 
    #                       map_signif_level = TRUE,
    #                       textsize = 12, # 4,
    #                       step_increase = 0.05,
    #                       test = "t.test",
    #                       vjust = 0.5, # -0.5,
    #                       annotations = add_significance(pv))
    #}
    
    fpl = paste(dir, filename_line.plot, sep = "")
    pdf(file = fpl, height = 6, width = 8)
    plot(p)
    dev.off()
  }
}