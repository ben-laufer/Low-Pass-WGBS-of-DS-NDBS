# chr21 CNV plot
# Run CNV_Me first (https://github.com/hyeyeon-hwang/CNV_Me)
# By Hyeyeon Hwang

library(magrittr)
library(tidyr)
library(dplyr)
library(tibble)
library(ggplot2)

data_chr21 <- read.delim("DS_DBS_cnv_caller_output_chr21.txt", header = TRUE)
sample_names <- colnames(data_chr21)[1:3]
for (path in colnames(data_chr21)[-c(1:3)]) {
  name <- path %>% strsplit("[.]") %>% unlist()
  name <- name[8] %>% strsplit("_") %>% unlist()
  name <- name[1]
  sample_names <- c(sample_names, name)
}

colnames(data_chr21) <- sample_names
data_chr21 <- data_chr21 %>% unite("chr.start.end", chr:end, sep=".", remove=TRUE)

sample_info <- read.csv("sample_info_all.csv") 
colnames(sample_info) <- c("Name", "Diagnosis", "Sex", "Batch")



diag_groups <- character()
for (sample in colnames(data_chr21)[-1]){
  sample_info_idx <- match(sample, sample_info$Name)
  diag_groups <- c(diag_groups, as.character(sample_info$Diagnosis[sample_info_idx]))
}
#diag_groups[which(diag_groups == "Control")] <- "TD"
diag_groups[which(diag_groups == "DevelopmentalDelay")] <- "Developmental Delay"
diag_groups[which(diag_groups == "DownSyndrome")] <- "Down Syndrome"
diag_groups[which(diag_groups == "Control")] <- "Typical Development"
# 7 = 1st index sample, 7th idx in sample_info$Name, 3, 24

diag_groups <- factor(diag_groups, levels = c("Down Syndrome", 
                                              "Developmental Delay", 
                                              "Typical Development"))

data_chr21_diag <- data_chr21 %>% 
  pivot_longer(cols = colnames(data_chr21[-1]), names_to = "Sample_name") %>%
  add_column(Diagnosis = rep(diag_groups, nrow(data_chr21))) 


# to plot the mean values, use this dataframe
data_chr21_diag_mean <- data_chr21_diag %>% 
  group_by(chr.start.end, Diagnosis) %>%
  do(mutate(., Mean = mean(value)))




pdf("DS_DBS_cnv_plot.pdf", width = 10, height = 6)
# gives smoothed lines as many as # of sample, 3 colors for diagnosis group
ggplot(data_chr21_diag, #data_chr21_diag_mean, 
       aes(x = chr.start.end,
           y = value, 
           col = Diagnosis,
           group = Sample_name)) + 
  #geom_point(size = 1, alpha = 1/5) +
  geom_smooth(se=FALSE, size = 0.5) + 
  scale_y_continuous(limits = c(1, 4)) + 
  theme_classic() + 
  xlab("Chromosome 21 (5KB Bins)") +
  ylab("Copy Number") +
  theme(
    text = element_text(size = 20),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(color = "black", fill = NA),
    #legend.key.height = unit(3, "cm")
    legend.position = "bottom"
  )
dev.off()

