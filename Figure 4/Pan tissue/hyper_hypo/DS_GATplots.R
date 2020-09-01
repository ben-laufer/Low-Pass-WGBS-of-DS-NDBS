#!/usr/bin/env Rscript

# DS GAT plots
# Ben Laufer

# Packages ----------------------------------------------------------------

if (!requireNamespace(c("tidyverse"), quietly = TRUE))
  install.packages("tidyverse", repos = "https://cloud.r-project.org")
suppressPackageStartupMessages(library(tidyverse))
options(readr.num_columns = 0)

# Tissue ------------------------------------------------------------------

cat("\n","Tidying DS tissue annotations...")
GAT <- readr::read_tsv("DS_tissues_hyper_hypo_results.tsv") %>%
  dplyr::mutate(annotation = dplyr::recode_factor(annotation,
                                                  "Fetal_Cortex" = "Fetal Frontal Cortex CpGs (450K)",
                                                  "Adult_Cerebellar_Folial_Cortex" = "Adult Cerebellar Folial Cortex CpGs (450K)",
                                                  "T_lymphocytes" = "Sorted Adult Peripheral T-lymphocytes (CD3+) CpGs (450K)",
                                                  "Neurons" = "Sorted Adult Frontal Cortex Neurons (NeuN+) CpGs (450K)",
                                                  "Adult_Frontal_Cortex" = "Adult Frontal Cortex CpGs (450K)",
                                                  "Glia" = "Sorted Adult Frontal Cortex Glia (NeuN-) CpGs (450K)",
                                                  "Mid-Gestation_Fetal_Cerebrum" = "Mid-Gestation Fetal Cerebrum CpGs (450K)",
                                                  "Whole-Blood" = "Whole-Blood DMR CpGs (450K)",
                                                  "PlacentaDMRs_sorted.bed" = "Placenta DMRs (RRBS)",
                                                  "Buccal_Epithelial_Cells" = "Buccal Epithelial Cell CpGs (450K)",
                                                  "neonatal_blood" = "Neonatal Blood CpGs (450K)",
                                                  "fc_wgbs_dmrs" = "Adult Frontal Cortex DMRs (WGBS)",
                                                  "neural_iPSC_derivatives" = "Neural iPSC Derivatives CpGs (450K)",
                                                  "HeartCpGs" = "Heart CpGs (450K)"
                                                  )
  ) %>%
  dplyr::mutate(fold = dplyr::case_when(fold < 1 ~ -1/fold,
                                        fold >= 1 ~ fold)
  ) %>%
  dplyr::mutate(signif = dplyr::case_when(qvalue <= 0.05 ~ 1,
                                          qvalue > 0.05 ~ 0)
  ) %>%
  dplyr::select(annotation,
                track,
                fold,
                signif)
cat("Done.")

cat("\n", "Plotting...")
(ggplot(data=GAT,
        aes(x = reorder(annotation,
                        fold),
                      y = fold,
                      fill = track)
        ) +
    geom_bar(stat = "identity") +
    coord_flip() +
    labs(y = "Fold Enrichment",
         x = element_blank()
         ) +
    theme_classic() + 
    theme(axis.text = element_text(size = 16),
          axis.title = element_text(size = 16),
          strip.text = element_text(size = 16),
          legend.position = "none"
          ) + 
    facet_grid(~track) +
    scale_y_continuous(limits = c(-8,8),
                       breaks = c(-8,-6,-4,-2,0,2,4,6,8)
    ) +
    geom_hline(yintercept = 0) +
    geom_text(data = GAT[(GAT$signif == 1 & GAT$fold > 0), ],
              label = "*",
              size = 8,
              show.legend = FALSE,
              nudge_y = 0.5,
              nudge_x = -0.05) +
    geom_text(data = GAT[(GAT$signif == 1 & GAT$fold < 0), ],
              label = "*",
              size = 8,
              show.legend = FALSE,
              nudge_y = -0.5,
              nudge_x = -0.05)
  ) %>%  
  ggsave("Gat_DS_tissues_hyper_hypo.pdf",
         plot = .,
         width = 11,
         height = 6)
cat("Done.")

