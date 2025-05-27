
library(tidyr)
retention_table <- read.csv("Anorexia_Nervosa_kneaddata_stats.tsv", sep='\t')
retention_table$other_filtering <- (retention_table$initial_reads - retention_table$trimmomatic_orphan - 
                                      retention_table$trimmomatic_paired - retention_table$trf_paired - 
                                      retention_table$trf_single - retention_table$human_paired) - retention_table$retained
retention_long <- gather(retention_table, step, count, initial_reads:other_filtering, factor_key = TRUE)
retention_long <- subset(retention_long, step != "initial_reads")

cols <- c("#bc7e9e", "#515d78","#870047", "#026fe9","#00a496","#9e91ca", "#b3cde0")

m <- max(retention_table$initial_reads)
breaks <- c(0, m/5, 2*m/5, 3*m/5, 4*m/5, m)
library(ggplot2)
ggplot(retention_long, aes(x=reorder(SampleID, -count), y=count, fill=factor(step, levels=c("trimmomatic_orphan", "trimmomatic_paired", "trf_single", "trf_paired", "human_paired", "other_filtering", "retained")))) +
  geom_bar(stat='identity') + theme_bw() +
  labs(x="Sample ID", y="Read Counts", fill="Removed By", title="Read Retention") +
  theme(text = element_text(size = 8.5, face = "bold"),
        element_line(size = 0.1),
        legend.position = "bottom",
        axis.text.x = element_text(angle = 290, hjust = 0),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        plot.margin = margin(t=20, b=20, r=20, l=20)) +
  guides(fill = guide_legend(nrow = 1)) +
  scale_fill_manual(values=cols, labels=c('Trimmomatic (Single)', 'Trimmomatic (Paired)', 'Tandem Repeats (Single)', 'Tandem Repeats (Paired)', 'Bowtie2 Human Reads', 'Other Filtering', 'Retained'))
ggsave('Anorexia_Nervosa_kneaddata_retention_plot.png', height=2200, width=3500, units='px')

subjs <- read.table("specimen_template.csv", sep=",", header=TRUE)$RawData.2
