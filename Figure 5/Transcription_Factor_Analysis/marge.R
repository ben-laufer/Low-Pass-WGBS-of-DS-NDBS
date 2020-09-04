# marge
# Parsing and summary heatmap for multiple DMR HOMER comparisons using MARGE
# By Ben Laufer

#' marge
#' @description Parse the HOMER results from multiple \code{DMRichR} runs.
#' @param contrast A character vector of the directory names for the different analyses
#'  i.e. c("DSvsTD", "DSvsDD", "DDvsTD").
#' @param direction A character vector of the DMR profiles to analyze
#'  c("both", hyper", "hypo").
#' @return A \code{tibble} of enrichment results.
#' @importFrom purrr map
#' @importFrom marge read_known_results
#' @importFrom dplyr filter mutate case_when select as_tibble
#' @importFrom tidyr separate
#' @importFrom magrittr %>%
#' @importFrom glue glue glue_collapse
#' @importFrom data.table rbindlist
#' @importFrom readxl read_xlsx
#' @references \url{https://github.com/robertamezquita/marge/blob/master/vignettes/marge-workflow.Rmd}
#' @export marge
#' 
marge <- function(contrast = contrast,
                  direction = c("both", "hyper", "hypo")
                  ){
  
  stopifnot(direction %in% c("both", "hyper", "hypo"))
  
  print(glue::glue("Tidying HOMER results for {tidyDirection} analyses from the {tidyContrasts} contasts",
                   tidyDirection = glue::glue_collapse({direction}, sep = ", ", last = " and "),
                   tidyContrasts = glue::glue_collapse({contrast}, sep = ", ", last = " and ")))
  
  outer(contrast,
        direction,
        function(contrast,
                 direction){
          glue::glue("{contrast}/Extra/HOMER/{direction}") 
        }) %>%
    as.vector() %>%
    purrr::map(function(file){
      marge::read_known_results(file)
    }) %>%
    `names<-` (outer(
      contrast,
      direction,
      function(contrast,
               direction){
        glue::glue("{contrast} {direction}") 
      }) %>%
        as.vector()
    ) %>% 
    data.table::rbindlist(idcol = "contrast") %>%
    dplyr::as_tibble() %>%
    tidyr::separate(contrast, c("contrast", "direction")) %>%
    dplyr::mutate(direction = dplyr::case_when(direction == "hyper" ~ "Hypermethylated",
                                               direction == "hypo" ~ "Hypomethylated")
    ) %>%
    return()
}

#' homerPlot
#' @description Plot HOMER transcription factor enrichment testing results
#'  from \code{DMRichR::marge()}.
#' @param data A \code{tibble} from \code{DMRichR::marge()}.
#' @return A \code{ggplot} object of enrichment results that can be viewed by calling it, 
#' saved with \code{ggplot2::ggsave()}, or further modified by adding \code{ggplot2} syntax.
#' @import ggplot2
#' @importFrom dplyr arrange desc distinct slice pull filter mutate select summarise group_by ungroup
#' @importFrom magrittr %>%
#' @importFrom glue glue glue_collapse
#' @references \url{https://github.com/robertamezquita/marge/blob/master/vignettes/marge-workflow.Rmd}
#' @export homerPlot 
#' 
homerPlot <- function(homer = homer){
  
  families <- homer %>%
    dplyr::arrange(dplyr::desc(log_p_value)) %>%
    dplyr::select(motif_family) %>%
    dplyr::distinct() %>%
    dplyr::slice(1:10) %>% 
    dplyr::pull()
  
  print(glue::glue("The top transcription factor families are {tidyFamilies}",
                   tidyFamilies = glue::glue_collapse({families}, sep = ", ", last = ", and ")))
  
  homer <- homer %>%
    dplyr::group_by(contrast, direction, motif_family) %>%
    dplyr::summarise("Top log p value" = max(log_p_value),
                     "Top FDR" = min(fdr)
                     ) %>%
    dplyr::ungroup() %>% 
    dplyr::filter(motif_family %in% families)
  
  homer %>%
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
      name = expression("-log"(p)), # expression("-log"[2](p))
      low = "black",
      high = "red" #,
      #limits = c(0,160)
    ) %>%
    return()
}

# Run ---------------------------------------------------------------------

if (!requireNamespace("marge", quietly = TRUE))
  BiocManager::install('robertamezquita/marge', ref = 'master')

packages <- c("tidyverse", "data.table", "marge")
stopifnot(suppressMessages(sapply(packages, require, character.only = TRUE)))
options(dplyr.summarise.inform = FALSE) 

setwd("/Users/benlaufer/Box/Down Syndrome/Dried Blood Spots/Results/DMRs/")

contrast <- c("DSvsTD", "DSvsDD", "DDvsTD")
direction <- c("hyper", "hypo")

marge(contrast = contrast,
      direction = direction) %>%
  dplyr::mutate(contrast = dplyr::case_when(contrast == "DSvsTD" ~ "DS vs TD",
                                            contrast == "DSvsDD" ~ "DS vs DD",
                                            contrast == "DDvsTD" ~ "DD vs TD")
                ) %>% 
  homerPlot() %>% 
  ggplot2::ggsave("Transcription factor motif enrichments.pdf",
                  plot = ., 
                  width = 8.5,
                  height = 6)
