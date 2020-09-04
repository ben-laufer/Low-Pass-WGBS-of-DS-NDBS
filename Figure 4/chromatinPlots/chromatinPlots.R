# chromatinPlots
# Parsing and summary heatmaps for multiple LOLA chromHMM and Roadmap enrichment analyses from multiple comparisons
# By Ben Laufer

#' chromatinReader
#' @description Parse LOLA chromHMM or Roadmap core histone modification enrichment
#'  results from multiple \code{DMRichR} runs.
#' @param contrast A character vector of the directory names for the different analyses
#'  i.e. c("DSvsTD", "DSvsDD", "DDvsTD") .
#' @param direction A character vector of the DMR profiles to analyze c("hyper", "hypo").
#' @param type A character vector of the type of results to parse c("chromHMM", "histone").
#' @return A \code{tibble} of enrichment results.
#' @importFrom purrr map
#' @importFrom readr read_tsv
#' @importFrom dplyr filter mutate case_when select as_tibble recode_factor arrange 
#' tally pull group_by left_join
#' @importFrom magrittr %>%
#' @importFrom stringr str_detect str_to_title
#' @importFrom glue glue glue_collapse
#' @importFrom data.table rbindlist
#' @importFrom forcats fct_unique
#' @export chromatinReader
#' 
chromatinReader <- function(contrast = contrast,
                            direction = direction,
                            type = c("chromHMM", "histone")
                            ){
  
  stopifnot(direction %in% c("All DMRs", "Hypermethylated DMRs", "Hypomethylated DMRs"))
  stopifnot(type %in% c("chromHMM", "histone"))
  
  print(glue::glue("Parsing {type} enrichment results for {tidyDirection} from {tidyContrasts}",
                   tidyDirection = glue::glue_collapse({direction}, sep = ", ", last = " and "),
                   tidyContrasts = glue::glue_collapse({contrast}, sep = ", ", last = ", and ")))
  
  directory <- dplyr::case_when(type == "chromHMM" ~ "ChromHMM",
                                type == "histone" ~ "RoadmapEpigenomics")
  
  data <- outer(contrast,
        direction,
        function(contrast,
                 direction){
          glue::glue("{contrast}/Extra/LOLA/{direction}/{directory}/allEnrichments.tsv") 
        }) %>%
    as.vector() %>%
    purrr::map(function(file){
      readr::read_tsv(file)
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
    dplyr::mutate(userSet = dplyr::case_when(stringr::str_detect(userSet, "All DMRs") ~ "All",
                                             stringr::str_detect(userSet, "Hypermethylated DMRs") ~ "Hypermethylated",
                                             stringr::str_detect(userSet, "Hypomethylated DMRs") ~ "Hypomethylated")
                  ) %>% 
    dplyr::arrange(qValue) 
  
  if(type == "chromHMM"){
    data <- data %>%
      dplyr::select(contrast, userSet, qValue, oddsRatio, cellType, tissue, antibody) %>%
      dplyr::mutate(antibody = forcats::as_factor(antibody)) %>%
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
                    ) 
    
    cat("\n", glue::glue("The following tissues are available for subsetting: {tissues}",
                         tissues = glue::glue_collapse({unique(data$tissue)},
                                                       sep = ", ", last = ", and ")), "\n")
    
  }else if(type == "histone"){
      data <- data %>%
        dplyr::left_join(readr::read_tsv("roadmap_epigenomics_index.txt"),
                         by = "filename") %>%
        dplyr::select(contrast, userSet, qValue, oddsRatio, donor, anatomy, antibody = antibody.x) %>%
        dplyr::mutate(antibody = tidyr::replace_na(antibody, "DNase")) %>% 
        dplyr::mutate(antibody = as.factor(antibody)) 
      
      # Get number of contrasts
      count <- outer(
        contrast,
        direction,
        function(contrast,
                 direction){
          glue::glue("{contrast} {direction}") 
        }) %>%
        length()
      
      # Subset for Core Histone Mods from same 127 samples: H3K4me1, H3K4me3, H3K27me3, H3K36me3, H3K9me3
      core <- data %>%
        dplyr::group_by(antibody) %>% 
        dplyr::tally() %>%
        dplyr::filter(n == 127*count) %>%
        dplyr::pull(antibody)
      
      data <- data %>% 
        dplyr::filter(antibody %in% core) %>%
        dplyr::mutate(anatomy = stringr::str_replace(anatomy, "_", " ")) %>% 
        dplyr::mutate(anatomy = stringr::str_to_title(anatomy)) %>%
        dplyr::mutate(anatomy = stringr::str_replace(anatomy, "Esc", "ESC")) %>%
        dplyr::mutate(anatomy = stringr::str_replace(anatomy, "Ipsc", "iPSC")) %>%
        dplyr::mutate(anatomy = stringr::str_replace(anatomy, "Gi", "GI")) %>%
        dplyr::rename(tissue = anatomy)
      
      print(glue::glue("Finshed parsing {mods} histone modifications for {tidyDirection} from {tidyContrasts}",
                       mods = glue::glue_collapse({data %>%
                           dplyr::select(tissue) %>%
                           levels()},
                           sep = ", ", last = " and "),
                       tidyDirection = glue::glue_collapse({direction}, sep = ", ", last = " and "),
                       tidyContrasts = glue::glue_collapse({contrast}, sep = ", ", last = ", and ")))
      
      print(glue::glue("The following tissues are available for subsetting: {tissues}",
                       tissues = glue::glue_collapse({unique(data$tissue)},
                                                     sep = ", ", last = ", and ")))
  }
  return(data)
}

#' chromatinMapper
#' @description Plot LOLA chromHMM and Roadmap enrichment testing results
#'  from \code{DMRichR::chromatinReader()}.
#' @param data A \code{tibble} from \code{DMRichR::chromatinReader()}.
#' @param type A character vector of the type of results to parse c("chromHMM", "histone").
#' @return A \code{ggplot} object of enrichment results that can be viewed by calling it, 
#' saved with \code{ggplot2::ggsave()}, or further modified by adding \code{ggplot2} syntax.
#' @import ggplot2
#' @importFrom dplyr group_by summarise ungroup
#' @importFrom magrittr %>%
#' @importFrom glue glue glue_collapse
#' @export chromatinMapper 
#' 
chromatinMapper <- function(data = data,
                            type = c("chromHMM", "histone")
                            ){
  
  stopifnot(type %in% c("chromHMM", "histone"))
  cat("\n", glue::glue("Plotting {type} annotation results"), "\n", "\n")
  
  data <- data %>%
    dplyr::group_by(contrast, userSet, tissue, antibody) %>%
    dplyr::summarise("Top log q value" = max(-log10(qValue))) %>%
    dplyr::ungroup() 
  
  p <- data %>%
      ggplot(aes(x = antibody,
                 y = contrast,
                 fill = `Top log q value`)) +
      geom_tile() +
      facet_grid(tissue~userSet) +
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
      geom_text(data = data[(data$`Top log q value` > -log10(0.05)), ],
                label = "*",
                size = 8,
                colour = "white",
                show.legend = FALSE,
                nudge_y = -0.15) +
      scale_fill_gradient(
        name = expression("-log"[10](q)),
        low = "black",
        high = "red"#,
        #limits = c(0,100)
        )
  if(type == "chromHMM"){
    p <- p +
      labs(y = element_blank(),
           x = "Chromatin State")
  }else if(type == "histone"){
    p <- p +
      labs(y = element_blank(),
           x = "Histone Modification")
  }
  return(p)  
}

# Run ---------------------------------------------------------------------

packages <- c("tidyverse", "data.table")
stopifnot(suppressMessages(sapply(packages, require, character.only = TRUE)))
options(readr.num_columns = 0)

setwd("/Users/benlaufer/Box/Down Syndrome/Dried Blood Spots/Results/DMRs/")

contrast <- c("DSvsTD", "DSvsDD", "DDvsTD")
direction <- c("Hypermethylated DMRs", "Hypomethylated DMRs")
type <- c("chromHMM", "histone")

purrr::walk(type,
            function(type){
             
             tissues <- dplyr::case_when(type == "chromHMM" ~ c("HSC & B-cell", "Blood & T-cell", "Brain"),
                                         type == "histone" ~ c("ESC", "Blood", "Brain"))
             
             chromatinReader(contrast = contrast,
                             direction = direction,
                             type = type) %>%
               dplyr::filter(tissue %in% tissues) %>%  
               dplyr::mutate(contrast = dplyr::case_when(stringr::str_detect(contrast, "DSvsTD") ~ "DS vs TD",
                                                         stringr::str_detect(contrast, "DSvsDD") ~ "DS vs DD",
                                                         stringr::str_detect(contrast, "DDvsTD") ~ "DD vs TD")
                             ) %>% 
               chromatinMapper(type = type) %>% 
               ggplot2::ggsave(glue::glue("{type}_heatmap.pdf"),
                               plot = ., 
                               width = 15,
                               height = 8)
           })
