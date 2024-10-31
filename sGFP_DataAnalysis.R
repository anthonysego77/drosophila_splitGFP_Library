# load libraries
library(tidyverse)
library(xlsx)
library(readxl)
library(ggpubr)

mRNAtoGFPExpressionData <- read_excel('~/Desktop/R-workingFolder/mRNACorrelations.xlsx')

#changed original xcel file to have only L CNS on sheet 3 to make it easier to work with
flyAtlas2CNSGenes <- read_excel('~/Desktop/R-workingFolder/FlyAtlas2_gene_data_v2.xlsx', sheet = "Sheet3")

#correlation plot
ggplot(mRNAtoGFPExpressionData, aes(x=normalized, y=mRNA)) +
  geom_point() +
  stat_cor(method = "pearson", label.x = -0, label.y = 80) +
  #stat_cor(method = "spearman", label.x = -0, label.y = 80) +
  
  geom_smooth(method = "lm", se=TRUE, aes(group=1), color = "black") +
  labs(x=" Average GFP Expression (normalized)", y="mRNA L3 CNS (FlyAtlas2)", title="Pearson Correlation")


# *3.2 correction value from flyatlas value to match flybase value
flyAtlas2CNSGenes <- flyAtlas2CNSGenes %>% mutate(flyBase_adjustment = `CNS L` * 3.2)

#remove all tRNA, snRNA, and other RNA related genes
flyAtlas2CNSGenes_filtered <- flyAtlas2CNSGenes %>%
  filter(!str_detect(Symbol, regex("RNA:", ignore_case = TRUE)))

log_data_CNS_Adj <- log10(flyAtlas2CNSGenes_filtered$flyBase_adjustment)  # Adjusted to extract the Adjustment column
clean_data <- log_data_CNS_Adj[is.finite(log_data_CNS_Adj)]  # Keep only finite values

#percentage of genes above log10(10) - 10 being the threshold for mRNA for detection in sGFP fluorescence in most cases
percent_above_1 <- df %>%
  summarize(percent_above_1 = mean(flyAtlas2CNSGenes_filtered$flyBase_adjustment > 1) * 100) %>%
  pull(percent_above_1)

percent_above_1

#head(clean_data)
#summary(clean_data)

#View(clean_data)

# Create the histogram plot
ggplot(data.frame(clean_data), aes(x = clean_data)) +
  geom_histogram(aes(y = ..count..), binwidth = 0.05, fill = "steelblue", color = "black") +  # Adjusted binwidth
  scale_x_continuous(trans = 'log10') + 
  labs(x = 'Log10 of 3.2 adjusted mRNA Expression in CNS', y = 'Frequency') +
  ggtitle('Histogram of mRNA expression by CNS') +
  #geom_vline(xintercept = log10(70), linetype = "dashed", color = "green") + 
  #annotate("text", label = "log10(70)", x=log10(100), y=100, angle = 90) 
  geom_vline(xintercept = log10(10), linetype = "dashed", color = "green") + 
  annotate("text", label = "", x=log10(8), y=100, angle = 90)+
  theme_minimal()

