# consensusPCA
# PCA of all samples using consensus DMRs
# Run consensusDMRs.sh first
# By Ben Laufer

#' PCA2
#' @description Performs and plots a PCA from individual smoothed methylation values
#' @param matrix Matrix of transposed individual methylation values
#' @param group Ordered character vector of sample groupings
#' @return A \code{ggplot} object that can be viewed by calling it,
#'  saved with \code{ggplot2::ggsave()}, or further modified by adding \code{ggplot2} syntax.
#' @import ggbiplot
#' @importFrom dplyr case_when
#' @importFrom forcats fct_rev
#' @importFrom glue glue
#' @export PCA2
#' 
PCA2 <- function(matrix = matrix,
                 group = NA){
  print(glue::glue("Performing PCA..."))
  data.pca <- prcomp(matrix, center = TRUE, scale. = TRUE)
  group <- factor(group, levels = unique(forcats::fct_rev(group)))
  
  cat("Plotting PCA...")
  PCA <- ggbiplot::ggbiplot(data.pca,
                            obs.scale = 1,
                            var.scale = 1,
                            groups = group,
                            ellipse = TRUE,
                            circle = FALSE,
                            var.axes = FALSE,
                            choices = 1:2) +
    scale_color_discrete(name = "Diagnosis") +
    theme_bw(base_size = 16) +
    geom_point(aes(colour = group), size = 4) +
    theme(legend.direction = "horizontal",
          legend.position = "bottom", # Change legend position
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 12),
          panel.grid.major = element_blank(),
          panel.border = element_rect(color = "black", size = 1.25),
          axis.ticks = element_line(size = 1.25),
          panel.grid.minor = element_blank())
  cat("Done", "\n")
  
  return(PCA)
}

# Run ---------------------------------------------------------------------

.libPaths("/share/lasallelab/programs/DMRichR/R_3.6")
packages <- c("DMRichR", "tidyverse")
stopifnot(suppressMessages(sapply(packages, require, character.only = TRUE)))
options(readr.num_columns = 0)

setwd("/share/lasallelab/Ben/DS_DBS/DMRs/cytosine_reports")

load("consensus/consensus_bismark.RData")
load("consensus/consensus_bsseq.RData")

# Tidy diagnosis

pData(bs.filtered.bsseq)$Diagnosis <- pData(bs.filtered.bsseq)$Diagnosis %>% 
  dplyr::recode_factor("DownSyndrome" = "Down Syndrome",
                       "DevelopmentalDelay" = "Developmental Delay",
                       "Control" = "Typical Development"
  )

# Read, getMeth, and plot

readr::read_tsv("consensus/ConsensusDMRs.bed",
                col_names = c("chr", "start", "end")
                ) %>% 
  GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE) %>%
  cbind(., data.frame(
    bsseq::getMeth(BSseq = bs.filtered.bsseq,
                   regions = .,
                   type = "smooth",
                   what = "perRegion"),
    check.names = FALSE)
    ) %>%
  dplyr::select(-seqnames, -start, -end, -width, -strand) %>% 
  na.omit() %>%
  as.matrix() %>%
  t() %>% 
  PCA2(group = bs.filtered.bsseq %>%
         pData() %>%
         dplyr::as_tibble() %>%
         dplyr::pull(!!testCovariate)
      ) %>%
  ggplot2::ggsave("consensus/conensusDMR_PCA.pdf",
                  plot = .,
                  device = NULL,
                  width = 10,
                  height = 6)
