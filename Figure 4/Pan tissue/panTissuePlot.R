# Dwon syndrome pan-tissue enrichment plot
# This script requires data from run_DS_GAT.sh
# By Ben Laufer

# Packages ----------------------------------------------------------------

if (!requireNamespace(c("tidyverse"), quietly = TRUE))
  install.packages("tidyverse", repos = "https://cloud.r-project.org")
suppressPackageStartupMessages(library(tidyverse))
options(readr.num_columns = 0)

# Tissue ------------------------------------------------------------------

#setwd("/share/lasallelab/Ben/DS_DBS/DMRs")
setwd("/Users/benlaufer/Box/Down Syndrome/Dried Blood Spots/Results/DMRs/")

cat("\n","Tidying DS tissue annotations...")

tidyGAT <- . %>% 
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
  dplyr::mutate(track = dplyr::recode_factor(track,
                                             "Hypomethylated" = "Hypomethylated",
                                             "Hypermethylated" = "Hypermethylated"
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

DSvsTD <- readr::read_tsv("DSvsTD/Extra/GAT/Tissue/DS_tissues_hyper_hypo_results.tsv") %>%
  tidyGAT() %>%
  dplyr::mutate(contrast = "DS vs TD")

DSvsDD <- readr::read_tsv("DSvsDD/Extra/GAT/Tissue/DS_tissues_hyper_hypo_results.tsv") %>%
  tidyGAT() %>%
  dplyr::mutate(contrast = "DS vs DD")

DDvsTD <- readr::read_tsv("DDvsTD/Extra/GAT/Tissue/DS_tissues_hyper_hypo_results.tsv") %>%
  tidyGAT() %>%
  dplyr::mutate(contrast = "DD vs TD")

GAT <- dplyr::bind_rows(DSvsTD, DSvsDD, DDvsTD) %>%
  dplyr::mutate(contrast = forcats::as_factor(contrast))
  
cat("Done.")

cat("\n", "Plotting...")
(ggplot(data = GAT,
        aes(x = reorder(annotation,
                        fold),
            y = fold,
            fill = track)
        ) +
    geom_bar(stat = "identity",
             position = "dodge",
             width = 0.75,
             color = "Black") + # remove
    coord_flip() +
    labs(y = "Fold Enrichment",
         x = element_blank(),
         fill = "Direction"
    ) +
    theme_classic() + 
    theme(axis.text = element_text(size = 20),
          axis.title = element_text(size = 20),
          strip.text = element_text(size = 20),
          legend.title = element_text(size = 20),
          legend.text = element_text(size = 20),
          legend.position = "bottom" #,
          #panel.grid.major.y = element_line(size = 1)
    ) + 
    facet_grid(~contrast) + # track~contrast
    scale_y_continuous(limits = c(-7.5,7.5),
                       breaks = c(-6,-4,-2,0,2,4,6)
    ) +
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = seq_len(nlevels(GAT$annotation) -1) + 0.5,
               color = "grey") +
    scale_fill_manual(values = c("#00BFC4", "#F8766D"),
                      guide = guide_legend(reverse = TRUE)
                      ) +
    geom_text(data = GAT[(GAT$signif == 1 & GAT$fold > 0) & GAT$track == "Hypermethylated", ],
              label = "*",
              size = 8,
              show.legend = FALSE,
              nudge_y = 0.5,
              nudge_x = 0.075) +
    geom_text(data = GAT[(GAT$signif == 1 & GAT$fold > 0 & GAT$track == "Hypomethylated"), ],
              label = "*",
              size = 8,
              show.legend = FALSE,
              nudge_y = 0.5,
              nudge_x = -0.3) +
    geom_text(data = GAT[(GAT$signif == 1 & GAT$fold < 0), ],
              label = "*",
              size = 8,
              show.legend = FALSE,
              nudge_y = -0.5,
              nudge_x = -0.3)
  ) %>%  
  ggsave("Gat_DS_tissues_hyper_hypo_all.pdf",
         plot = .,
         width = 20,
         height = 10)
cat("Done.")

