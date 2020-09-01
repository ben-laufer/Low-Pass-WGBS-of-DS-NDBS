# Transcription factor family heatmaps for homer results
# Ref: https://github.com/robertamezquita/marge/blob/master/vignettes/marge-workflow.Rmd
# By Ben Laufer

# Packages ----------------------------------------------------------------
rm(list = ls())

#.libPaths("/share/lasallelab/programs/DMRichR/R_3.6")
BiocManager::install('robertamezquita/marge', ref = 'master')
packages <-  c("tidyverse", "marge")
stopifnot(suppressMessages(sapply(packages, require, character.only = TRUE)))

# Load data ---------------------------------------------------------------

#setwd("/share/lasallelab/Ben/DS_DBS/DMRs")
setwd("/Users/benlaufer/Box/Down Syndrome/Dried Blood Spots/Results/DMRs/")

homer <- map(c("DSvsTD/Extra/HOMER/hyper",
               "DSvsTD/Extra/HOMER/hypo",
               "DSvsDD/Extra/HOMER/hyper",
               "DSvsDD/Extra/HOMER/hypo",
               "DDvsTD/Extra/HOMER/hyper",
               "DDvsTD/Extra/HOMER/hypo"
               ),
             marge::read_known_results
             ) %>%
  `names<-` (c("DSvsTD_hyper",
               "DSvsTD_hypo",
               "DSvsDD_hyper",
               "DSvsDD_hypo",
               "DDvsTD_hyper",
               "DDvsTD_hypo"
               )
             ) %>% 
  data.table::rbindlist(idcol = "contrast") %>%
  as_tibble() %>%
  mutate(direction = case_when(str_detect(contrast, "hyper") ~ "Hypermethylated",
                               str_detect(contrast, "hypo") ~ "Hypomethylated"
                               )
         ) %>%
  mutate(contrast = case_when(str_detect(contrast, "DSvsTD") ~ "DS vs TD",
                              str_detect(contrast, "DSvsDD") ~ "DS vs DD",
                              str_detect(contrast, "DDvsTD") ~ "DD vs TD"
                              )
         )

families <- homer %>%
  arrange(desc(log_p_value)) %>%
  select(motif_family) %>%
  unique() %>%
  dplyr::slice(1:10) %>% 
  pull()

homer <- homer %>%
  group_by(contrast, direction, motif_family) %>%
  summarise("Top log p value" = max(log_p_value),
            "Top FDR" = min(fdr)) %>%
  ungroup() %>% 
  filter(motif_family %in% families)

(homer %>%
    ggplot(aes(x = motif_family,
               y = contrast,
               fill = `Top log p value`)) +
    geom_tile() +
    facet_grid(direction ~ .) +
    labs(y = element_blank(),
         x = "Motif Family"
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
    geom_text(data = homer[(homer$`Top FDR` < 0.05), ], # q-value (Benjamini)
              label = "*",
              size = 8,
              colour = "white",
              show.legend = FALSE,
              nudge_y = -0.15) +
    scale_fill_gradient(
      name = expression("-log"(p)), 
      low = "black",
      high = "red",
      limits = c(0,160)
      )
  ) %>% 
  ggsave("Transcription factor motif enrichments.pdf",
         plot = ., 
         width = 8.5,
         height = 6)
