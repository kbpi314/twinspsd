### Load libraries ###

library(phyloseq)
library(vegan)
library(ade4)
library(PMCMR)
library(PMCMRplus)
library(ggplot2)
library(ggthemes)
library(ggrepel)

############################################################################
############################################################################
############################################################################

### Beta diversity pcoa ###

# background theme
bkg <- theme_bw() +
  theme(axis.text = element_text(size = 24, color = "black")) +
  theme(axis.title = element_text(size = 32, color = "black", face = "bold")) +
  theme(axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm"))) +
  theme(axis.title.y = element_text(margin = unit(c(0, 8, 0, 0), "mm"))) +
  theme(legend.text = element_text(size = 18, color = "black"))+ #, face = "bold")) +
  # theme(legend.title = element_blank()) +
  theme(legend.title = element_text(size = 24, face = "bold", color = "black"))
  theme(legend.justification = "right")# +
  #theme(panel.border = element_rect(color = "black", fill = NA, size = 1.5))

# function to specify that axis labels have 2 decimal places
f.dec <- function(x){
  format(round(x, 2), nsmall = 2)
}

# directory for storing files
# dir = "/Users/Lyusik/Dropbox/Julia_MiCRA/2019_Projects/RA_twins/16S/jobs/2_beta.div_pcoa/"
dir = "/Users/KevinBu/Desktop/clemente_lab/Projects/twinspsd/outputs/jobs02/"

# list of distance methods
# dists <- c("bray", "unifrac", "wunifrac")
dists <- c('bray')

# colors
col1 <- c("#929aab", "#ce2525")
col1 <- c("#929aab","#5DB54B")#,"#FB9A99")
# load data
df <- read.delim(file="/Users/KevinBu/Desktop/clemente_lab/Projects/twinspsd/outputs/jobs02/bray_curtis_pcoa.tsv",
              row.names=1)
df = subset(df, df$Diagnosis != "PsA")

# drop na 
df = df[!is.na(df$visit),]

# order factors for legend
df$Diagnosis <- factor(df$Diagnosis, levels=c('HC', 'PsO'))

# create filenames
filename_plot = paste("bdiv", 'unweighted_unifrac', "plot.pdf", sep = "_")

# plot beta diversity
p <- ggplot(df, aes(PC1, PC2)) + # data=df, aes(x = PC1, y = PC2, fill = Diagnosis)) +
  geom_point(data = df, aes(x = PC1, y = PC2, color = Diagnosis),size=4) +
  geom_line(aes(group = twinpair), lty = 2, colour = "black") + 
  scale_color_manual(values = c("HC" = "#929aab", "PsO" = "#5DB54B")) + #col1, labels = c("HC", "PsO")) +
  bkg + guides(shape = "none") +
  xlab("PC1 (20.2%)") +
  ylab("PC2 (9.4%)") +
  scale_x_continuous(labels = f.dec) + # 2 decimal places on x-axis
  scale_y_continuous(labels = f.dec)   # 2 decimal places on y-axis

# save plot
fp = paste(dir, filename_plot, sep = "")
pdf(file = fp, height = 6, width = 8)
plot(p)
dev.off()
