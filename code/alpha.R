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
library(rstatix)
library(ggpubr)

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

# set working dir
dir = "/Users/KevinBu/Desktop/clemente_lab/Projects/twinspsd/outputs/jobs02/"

### Alpha Diversity Boxplots ###
df_alpha = read.table('/Users/KevinBu/Desktop/clemente_lab/Projects/twinspsd/inputs/df_merge_alpha.tsv', 
                      sep = '\t', header = TRUE, row.names = 1, check.names = FALSE,
                      na.strings = "NA")


# drop na 
df_alpha = df_alpha[!is.na(df_alpha$visit),]

# set factor level orders
df_alpha$Diagnosis <- factor(df_alpha$Diagnosis, levels = c("HC", "PsO", "PsA"))
df_alpha$DxPsD <- factor(df_alpha$DxPsD, levels = c("HC", "PsD"))

# create tables for storing wilcoxon and ttest results
stats.table.all <- matrix(data = NA, nrow = 1, ncol = 3)
# colnames(stats.table.all) <- c("alpha div", "wilcoxon", "ttest")
colnames(stats.table.all) <- c("alpha div", "kruskal", "anova")

# calculate adiv
# stats.table.all[1,1] <- colnames(df_alpha)[4]
stats.table.all[1,1] <- 'Alpha_Diversity'
#stats.table.all[1,2] <- wilcox.test(df_alpha[,4] ~ Diagnosis, data = df_alpha, paired = TRUE)$p.value
#stats.table.all[1,3] <- t.test(df_alpha[,4] ~ Diagnosis, data = df_alpha, paired = TRUE)$p.value
stats.table.all[1,2] <- kruskal.test(Alpha_Diversity ~ Diagnosis, data = df_alpha)$p.value
stats.table.all[1,3] <- anova_test(Alpha_Diversity ~ Diagnosis, data = df_alpha)$p

# save
ft.all = paste(dir, "alpha_stats.csv", sep = "")
write.csv(file = ft.all, stats.table.all)

# lms for rarefaction
df = read.table('/Users/KevinBu/Desktop/clemente_lab/Projects/twinspsd/inputs/df_merge_alpha_kd.tsv', 
                      sep = '\t', header = TRUE, row.names = 1, check.names = FALSE, stringsAsFactors = TRUE)

# df$Diagnosis <- as.character(df$Diagnosis)
df$Diagnosis <- factor(df$Diagnosis, levels=c("HC","PsO","PsA"))
contrasts(df$Diagnosis) <- contr.treatment(levels(df$Diagnosis))
model = glm(Alpha_Diversity ~ retained + Diagnosis, data=df)
summary(model)


############################################################################
############################################################################
############################################################################

### Analytes plots - all ###

# background theme
bkg <- theme_bw() +
  theme(axis.text.x = element_text(size = 24, face = "bold", color = "black")) +
  theme(axis.text.y = element_text(size = 18))+#, color = "black")) +
  theme(axis.title.y = element_text(size = 24, color = "black", face = "bold")) +
  theme(axis.title.y = element_text(margin = unit(c(0, 8, 0, 0), "mm"))) +
  theme(legend.text = element_text(size = 18, color = "black")) +
  theme(legend.title = element_text(size = 24, face = "bold", color = "black"))

# choose line types
line1 <- c("solid", "dashed", "dotted")

# choose colors
col1 <- c("white","#5DB54B","#FB9A99")
col4 <- c('white','purple') # PsD
col3 <- c("white","#FB9A99") # PsA only
col2 <- c("#f3a333", "#0074e4", "#8f8787")

# variable of interest
a <- 'Alpha_Diversity'

# create filenames
filename_table = paste(a, "all_table.csv", sep = "_")
filename_box.plot = paste(a, "all_box.plot.pdf", sep = "_")  
filename_line.plot = paste(a, "all_line.plot.pdf", sep = "_")  


# plot boxplot for all data
p <- ggplot(df_alpha, aes(x = Diagnosis, y = Alpha_Diversity, fill = Diagnosis)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position='jitter') + 
  # geom_jitter(width = 0.2, alpha = 0.7, size = 2) + 
  # theme(legend.position = "none") + 
  ylab("Alpha Diversity (Shannon)") +
  xlab(element_blank()) +
  geom_pwc(method = 'wilcox.test', 
           label = 'p.signif',  
           hide.ns = TRUE, 
           p.adjust.method = 'none',
           vjust = 0.5, # default 0, positive pushes it towards bar, negative further vertically away up; 1 is on the bar
           size = 0.8, # default 0.3
           label.size = 8, # default 3.88
           symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 0.1, Inf), 
                            symbols = c("****", "***", "**", "*", ".", "ns"))
  ) +
  scale_fill_manual(values = col1) +
  bkg 

fpb = paste(dir, filename_box.plot, sep = "")
pdf(file = fpb, height = 6, width = 8)
plot(p)
dev.off()


# paired data only merging PsD
df_alpha = read.table('/Users/KevinBu/Desktop/clemente_lab/Projects/twinspsd/inputs/df_merge_alpha_filt.tsv', 
                      sep = '\t', header = TRUE, row.names = 1, check.names = FALSE,
                      na.strings = "NA")

df_alpha$DxPsD <- factor(df_alpha$DxPsD, levels = c("HC", "PsD"))

# create categories by which to color lineplot: 1) Increase, 2) Decrease, 3) No change
# subset and spread dataset into diagnosis columns
# then calculate delta relative abundance between PsD and Unaffected siblings
d.div <- df_alpha %>%
  subset(select = c("twinpair", "DxPsD", a)) %>%
  spread(key = "DxPsD", value = a) %>%
  mutate(diff_sib.pair = (get("PsD") - get("HC"))) %>%
  mutate(
    Change.type_sib.pair = case_when(
      sign(diff_sib.pair) == 1 ~ "1_up",
      sign(diff_sib.pair) == -1 ~ "2_down",
      TRUE ~ "3_no.change"
    )
  ) 

# merge with original dataset
d.final <- gather(d.div, "HC", "PsD", key = "DxPsD", value = "Alpha_Diversity")

# save dataset
ft = paste(dir,  'PsD', filename_table,sep = "")
write.csv(d.final, file = ft)

# rewrite order of factors
d.final$Diagnosis <- factor(d.final$DxPsD, levels = c("HC", "PsD"))

# plot lineplot
p <- ggplot(data = d.final, aes(x = DxPsD, y = Alpha_Diversity, fill = DxPsD)) +
  stat_summary(fun.data = stats.whiskers, geom = "errorbar", 
               color = "black", size = 0.8, width = 0.3) +
  stat_summary(fun.data = stats.boxplot, geom = "crossbar",
               color = "black", size = 0.5, width = 0.5) +
  geom_point(shape = 20, size = 3, color = "black") +
  geom_line(data = subset(d.final), aes(group = twinpair, color = Change.type_sib.pair, linetype = Change.type_sib.pair), size = 0.5) +
  geom_text_repel(data = subset(d.final, Diagnosis == "HC"), aes(label = twinpair), 
                  nudge_x = -0.5, size = 2, color = "black", segment.alpha = 0.3) +
  scale_x_discrete(labels = c("HC", "PsD")) +
  scale_fill_manual(values = col4) +
  scale_linetype_manual(values = line1, name = "Entropy Difference", labels = c("Increase", "Decrease", "No change")) +
  scale_color_manual(values = col2, name = "Entropy Difference", labels = c("Increase", "Decrease", "No change")) +      
  xlab(NULL) +
  ylab("Shannon Entropy") +
  bkg

filename_line.plot = paste(a, "PsD_all_line.plot.pdf", sep = "_")  
fpl = paste(dir, filename_line.plot, sep = "")
pdf(file = fpl, height = 6, width = 8)
plot(p)
dev.off()


# paired data only for PsO
df_alpha = read.table('/Users/KevinBu/Desktop/clemente_lab/Projects/twinspsd/inputs/df_merge_alpha_filt.tsv', 
                      sep = '\t', header = TRUE, row.names = 1, check.names = FALSE,
                      na.strings = "NA")

df_alpha = df_alpha[df_alpha$Diagnosis %in% c('PsO','HC'),]
df_alpha <- df_alpha %>% group_by(twinpair) %>% mutate(n=n()) %>% filter(n>1) %>% ungroup() %>% select(-n)


df_alpha$Diagnosis <- factor(df_alpha$Diagnosis, levels = c("HC", "PsO"))

# create categories by which to color lineplot: 1) Increase, 2) Decrease, 3) No change
# subset and spread dataset into diagnosis columns
# then calculate delta relative abundance between PsD and Unaffected siblings
d.div <- df_alpha %>%
  subset(select = c("twinpair", "Diagnosis", a)) %>%
  spread(key = "Diagnosis", value = a) %>%
  mutate(diff_sib.pair = (get("PsO") - get("HC"))) %>%
  mutate(
    Change.type_sib.pair = case_when(
      sign(diff_sib.pair) == 1 ~ "1_up",
      sign(diff_sib.pair) == -1 ~ "2_down",
      TRUE ~ "3_no.change"
    )
  ) 

# merge with original dataset
d.final <- gather(d.div, "HC", "PsO", key = "Diagnosis", value = "Alpha_Diversity")

# save dataset
ft = paste(dir, 'PsO', filename_table, sep = "")
write.csv(d.final, file = ft)

# rewrite order of factors
d.final$Diagnosis <- factor(d.final$Diagnosis, levels = c("HC", "PsO"))

# plot lineplot
p <- ggplot(data = d.final, aes(x = Diagnosis, y = Alpha_Diversity, fill = Diagnosis)) +
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
  scale_linetype_manual(values = line1, name = "Entropy Difference", labels = c("Increase", "Decrease", "No change")) +
  scale_color_manual(values = col2, name = "Entropy Difference", labels = c("Increase", "Decrease", "No change")) +      
  xlab(NULL) +
  ylab("Shannon Entropy") +
  bkg

filename_line.plot = paste(a, "PsO_all_line.plot.pdf", sep = "_")  
fpl = paste(dir, filename_line.plot, sep = "")
pdf(file = fpl, height = 6, width = 8)
plot(p)
dev.off()



# paired data only for PsA
df_alpha = read.table('/Users/KevinBu/Desktop/clemente_lab/Projects/twinspsd/inputs/df_merge_alpha_filt.tsv', 
                      sep = '\t', header = TRUE, row.names = 1, check.names = FALSE,
                      na.strings = "NA")

df_alpha = df_alpha[df_alpha$Diagnosis %in% c('PsA','HC'),]
df_alpha <- df_alpha %>% group_by(twinpair) %>% mutate(n=n()) %>% filter(n>1) %>% ungroup() %>% select(-n)

df_alpha$Diagnosis <- factor(df_alpha$Diagnosis, levels = c("HC", "PsA"))

# create categories by which to color lineplot: 1) Increase, 2) Decrease, 3) No change
# subset and spread dataset into diagnosis columns
# then calculate delta relative abundance between PsD and Unaffected siblings
d.div <- df_alpha %>%
  subset(select = c("twinpair", "Diagnosis", a)) %>%
  spread(key = "Diagnosis", value = a) %>%
  mutate(diff_sib.pair = (get("PsA") - get("HC"))) %>%
  mutate(
    Change.type_sib.pair = case_when(
      sign(diff_sib.pair) == 1 ~ "1_up",
      sign(diff_sib.pair) == -1 ~ "2_down",
      TRUE ~ "3_no.change"
    )
  ) 

# merge with original dataset
d.final <- gather(d.div, "HC", "PsA", key = "Diagnosis", value = "Alpha_Diversity")

# save dataset
ft = paste(dir, 'PsA', filename_table, sep = "")
write.csv(d.final, file = ft)

# rewrite order of factors
d.final$Diagnosis <- factor(d.final$Diagnosis, levels = c("HC", "PsA"))

# plot lineplot
p <- ggplot(data = d.final, aes(x = Diagnosis, y = Alpha_Diversity, fill = Diagnosis)) +
  stat_summary(fun.data = stats.whiskers, geom = "errorbar", 
               color = "black", size = 0.8, width = 0.3) +
  stat_summary(fun.data = stats.boxplot, geom = "crossbar",
               color = "black", size = 0.5, width = 0.5) +
  geom_point(shape = 20, size = 3, color = "black") +
  geom_line(data = subset(d.final), aes(group = twinpair, color = Change.type_sib.pair, linetype = Change.type_sib.pair), size = 0.5) +
  geom_text_repel(data = subset(d.final, Diagnosis == "HC"), aes(label = twinpair), 
                  nudge_x = -0.5, size = 2, color = "black", segment.alpha = 0.3) +
  scale_x_discrete(labels = c("HC", "PsA")) +
  scale_fill_manual(values = col3) +
  scale_linetype_manual(values = line1, name = "Entropy Difference", labels = c("Increase", "Decrease", "No change")) +
  scale_color_manual(values = col2, name = "Entropy Difference", labels = c("Increase", "Decrease", "No change")) +      
  xlab(NULL) +
  ylab("Shannon Entropy") +
  bkg

filename_line.plot = paste(a, "PsA_all_line.plot.pdf", sep = "_")  
fpl = paste(dir, filename_line.plot, sep = "")
pdf(file = fpl, height = 6, width = 8)
plot(p)
dev.off()




