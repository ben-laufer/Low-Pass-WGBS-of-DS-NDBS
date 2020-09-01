# Analysis and plot of integration of WGBS data with RUNX1 ChIP-seq data
# Place Table S2 from Sanda et al. in a folder called runx1_chip in the working directory
# By Ben Laufer

# Load packages -----------------------------------------------------------

.libPaths("/share/lasallelab/programs/DMRichR/R_4.0")
AnnotationHub::setAnnotationHubOption("CACHE", "/share/lasallelab/programs/DMRichR/R_4.0")
ExperimentHub::setExperimentHubOption("CACHE", "/share/lasallelab/programs/DMRichR/R_4.0")
suppressPackageStartupMessages(library(DMRichR, attach.required = T))

setwd("/share/lasallelab/Ben/DS_DBS/DMRs/")

# Load data ---------------------------------------------------------------

load("cytosine_reports/consensus/consensus_bismark.RData")
load("cytosine_reports/consensus/consensus_bsseq.RData")

# Tidy diagnosis ----------------------------------------------------------

pData(bs.filtered.bsseq)$Diagnosis <- pData(bs.filtered.bsseq)$Diagnosis %>% 
  dplyr::recode_factor(.,
                       "DownSyndrome" = "Down Syndrome",
                       "DevelopmentalDelay" = "Developmental Delay",
                       "Control" = "Typical Development"
  )

# Process RUNX1 ChIP-seq data ---------------------------------------------
# Table S2 from Sanda et al.

RUNX1 <- readxl::read_xlsx("runx1_chip/mmc3.xlsx", sheet = 1, skip = 2, col_types = "text") %>%
  dplyr::filter(RUNX1 == 1) %>% # Select RUNX1 peaks
  dplyr::select(2:4) %>%
  magrittr::set_colnames(c("chr", "start", "end")) %>% 
  dplyr::mutate(chr = paste0("chr",chr)) %>% 
  GenomicRanges::makeGRangesFromDataFrame() %>%
  rtracklayer::liftOver(AnnotationHub::AnnotationHub()[["AH14221"]]) %>% # liftOver from hg18 to hg38
  unlist() %>%
  bsseq::getMeth(BSseq = bs.filtered.bsseq,
                   regions = .,
                   type = "smooth",
                   what = "perRegion") %>%
  na.omit() %>%
  dplyr::as_tibble() %>%
  dplyr::summarize_all(mean) %>%
  tidyr::gather(key = "sample",
                value = "CpG_Avg") %>%
  dplyr::left_join(pData(bs.filtered.bsseq) %>% 
                     as.data.frame() %>% 
                     tibble::rownames_to_column(var = "sample") %>%
                     dplyr::as_tibble())


# Stats -------------------------------------------------------------------

summary <- RUNX1 %>%
  group_by(Diagnosis, Sex) %>%
  summarise(
    count = n(),
    mean = mean(CpG_Avg, na.rm = TRUE),
    sd = sd(CpG_Avg, na.rm = TRUE)
  )

ANOVA <- RUNX1 %>%
  aov(CpG_Avg ~ Diagnosis + Sex, data = .)

Tukey <- ANOVA %>%
  TukeyHSD(which = "Diagnosis")

pairWise <- lm(CpG_Avg ~ Diagnosis + Sex, data = RUNX1) %>% 
  lsmeans::ref.grid(data = RUNX1) %>%
  lsmeans::lsmeans(~Diagnosis) %>%
  pairs(reverse = TRUE) %>%
  summary()

list("globalInput" = RUNX1,
     "Summary" = summary,
     "ANOVA" = broom::tidy(ANOVA),
     "Tukey" = broom::tidy(Tukey),
     "pairWise" = pairWise
)%>%
  openxlsx::write.xlsx("runx1_chip/runx1Stats.xlsx")


# Plot --------------------------------------------------------------------

library(ggsignif)

RUNXplot <- RUNX1 %>% 
  dplyr::mutate(CpG_Avg = CpG_Avg * 100) %>% 
  ggplot(aes(x = Diagnosis, y = CpG_Avg, fill = Diagnosis)) + 
  geom_boxplot() +
  labs(x = "", y = "Percent Methylation") +
  theme_classic(base_size = 30) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
  ylim(45, 60) +
  geom_signif(comparisons = list(
    c("Down Syndrome", "Developmental Delay"),
    c("Down Syndrome", "Typical Development")),
  map_signif_level = c("*" = 0.05), textsize = 8, step_increase = 0.15) +
  scale_x_discrete(labels=c("Down Syndrome" = "DS", "Developmental Delay" = "DD",
                            "Typical Development" = "TD"))

ggplot2::ggsave("runx1_chip/Smoothed RUNX1 peak methylation.pdf",
                plot = RUNXplot,
                device = NULL,
                width = 8.5,
                height = 8.5)
