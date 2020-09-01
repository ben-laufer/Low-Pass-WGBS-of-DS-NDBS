# LOLA chromHMM heatmap
# The analysis in this script requires the LOLA database to be installed when running DMRichR
# Ben Laufer

# Packages ----------------------------------------------------------------
rm(list = ls())

#.libPaths("/share/lasallelab/programs/DMRichR/R_3.6")
packages <-  c("tidyverse", "LOLA")
stopifnot(suppressMessages(sapply(packages, require, character.only = TRUE)))
options(readr.num_columns = 0)

# Load data ---------------------------------------------------------------

#setwd("/share/lasallelab/Ben/DS_DBS/DMRs")
setwd("/Users/benlaufer/Box/Down Syndrome/Dried Blood Spots/Results/DMRs/")

chromHMM <- map(c("DSvsTD/Extra/LOLA/Hypermethylated DMRs/ChromHMM/allEnrichments.tsv",
                  "DSvsTD/Extra/LOLA/Hypomethylated DMRs/ChromHMM/allEnrichments.tsv",
                  "DSvsDD/Extra/LOLA/Hypermethylated DMRs/ChromHMM/allEnrichments.tsv",
                  "DSvsDD/Extra/LOLA/Hypomethylated DMRs/ChromHMM/allEnrichments.tsv",
                  "DDvsTD/Extra/LOLA/Hypermethylated DMRs/ChromHMM/allEnrichments.tsv",
                  "DDvsTD/Extra/LOLA/Hypomethylated DMRs/ChromHMM/allEnrichments.tsv"
                  ),
            readr::read_tsv
            ) %>%
  `names<-` (c("DS vs TD",
               "DS vs TD",
               "DS vs DD",
               "DS vs DD",
               "DD vs TD",
               "DD vs TD"
               )
             ) %>% 
  data.table::rbindlist(idcol = "contrast") %>%
  as_tibble() %>%
  dplyr::select(contrast, userSet, qValue, oddsRatio, cellType, tissue, antibody) %>%
  dplyr::mutate(antibody = as.factor(antibody)) %>% 
  dplyr::mutate(antibody = dplyr::recode_factor(antibody,
                                                "01_TssA" = "Active TSS",
                                                "02_TssAFlnk" = "Flanking Active TSS",
                                                "03_TxFlnk" = "Transcription at Gene 5' and 3'",
                                                "04_Tx" = "Strong Transcription",
                                                "05_TxWk" = "Weak Transcription",
                                                "06_EnhG"= "Genic Enhancers",
                                                "07_Enh" = "Enhancers",
                                                "08_ZnfRpts" = "ZNF Genes & Repeats",
                                                "09_Het" = "Heterochromatin",
                                                "10_TssBiv" = "Bivalent/Poised TSS",
                                                "11_BivFlnk" = "Flanking Bivalent TSS/Enhancer",
                                                "12_EnhBiv" = "Bivalent Enhancer",
                                                "13_ReprPC" = "Repressed PolyComb",
                                                "14_ReprPCwk" = "Weak Repressed PolyComb",
                                                "15_Quies" = "Quiescent/Low"
                                                )
                ) %>%
  dplyr::mutate(userSet = dplyr::case_when(str_detect(userSet, "Hypermethylated DMRs") ~ "Hypermethylated",
                                           str_detect(userSet, "Hypomethylated DMRs") ~ "Hypomethylated"
                                           )
                ) %>% 
  dplyr::arrange(qValue) %>%
  dplyr::filter(tissue %in% (c("HSC & B-cell", "Blood & T-cell", "Brain"))) 

# Heatmap -----------------------------------------------------------------

heatData <- chromHMM %>%
  dplyr::group_by(contrast, userSet, tissue, antibody) %>%
  dplyr::summarise("Top log q value" = max(-log10(qValue))) %>%
  ungroup() 

(heatData %>%
    ggplot(aes(x = antibody,
            y = contrast,
            fill = `Top log q value`)) +
    geom_tile() +
    facet_grid(tissue~userSet) +
    labs(y = element_blank(),
         x = "Chromatin State"
    ) +
    theme_classic() + 
    theme(axis.text = element_text(size = 16),
          axis.title = element_text(size = 16),
          strip.text = element_text(size = 16),
          legend.text = element_text(size = 8),
          legend.title = element_text(size = 10),
          axis.text.x = element_text(angle = 45,
                                     hjust = 1)
    ) +
    guides(fill = guide_colourbar(barwidth = 0.75,
                                  barheight = 5)
    ) + 
    geom_text(data = heatData[(heatData$`Top log q value` > -log10(0.05)), ],
              label = "*",
              size = 8,
              colour = "white",
              show.legend = FALSE,
              nudge_y = -0.15) +
    scale_fill_gradient(
      name = expression("-log"[10](q)),
      low = "black",
      high = "red",
      limits = c(0,100)
    )
  ) %>% 
  ggsave("ChromHMM enrichments heatmap.pdf",
         plot = ., 
         width = 15,
         height = 8)

