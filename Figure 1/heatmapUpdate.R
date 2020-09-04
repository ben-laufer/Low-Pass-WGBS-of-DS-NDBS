# heatmapUpdate
# Update DMRichR heatmaps from multiple comparisons to have a consistent style
# By Ben Laufer

#' smoothPheatmap2
#' @description Plot a heatmap of normalized individual smoothed methylation value z scores for selected regions (i.e. significant DMRs)
#' @param bs.filtered.bsseq Smoothed \code{bsseq} object
#' @param sigRegions \code{GRanges} object of regions to plot a heatmap for
#' @param testCovariate The factor tested for differences between groups
#' @param filename Character specifying the name of the heatmap image file
#' @param ... Additional arguments passed onto \code{pheatmap()}
#' @return Saves a pdf image of the heatmap in the DMR folder
#' @import pheatmap pheatmap
#' @importFrom dplyr select_if
#' @importFrom RColorBrewer brewer.pal
#' @importFrom bsseq getMeth
#' @importClassesFrom bsseq BSseq 
#' @importMethodsFrom bsseq pData
#' @importFrom glue glue
#' @importFrom magrittr %>% set_colnames
#' @references \url{https://davetang.org/muse/2018/05/15/making-a-heatmap-in-r-with-the-pheatmap-package/}
#' @export smoothPheatmap2
#' 
smoothPheatmap2 <- function(bs.filtered.bsseq = bs.filtered.bsseq,
                            sigRegions = sigRegions,
                            testCovariate = testCovariate,
                            annotation_colors = annotation_colors,
                            filename = "DMRs.pdf",
                            ...){
  cat("\n[DMRichR] DMR heatmap \t\t\t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
  
  bsseq::getMeth(BSseq = bs.filtered.bsseq,
                 regions = sigRegions,
                 type = "smooth",
                 what = "perRegion") %>% 
    as.matrix() %>%
    pheatmap::pheatmap(.,
                       scale = "row",
                       annotation_col = pData(bs.filtered.bsseq) %>%
                         as.data.frame() %>%
                         dplyr::select_if(~ nlevels(.) > 1),
                       color = RColorBrewer::brewer.pal(11, name = "RdBu") %>%
                         rev(),
                       show_colnames = F,
                       border_color = NA,
                       main = glue::glue("{length(sigRegions)} DMRs"),
                       fontsize = 16,
                       filename = filename,
                       width = 11,
                       height = 8.5,
                       annotation_colors = annotation_colors,
                       ...
    ) %>%
    return()
}

# Run ---------------------------------------------------------------------

.libPaths("/share/lasallelab/programs/DMRichR/R_3.6")
suppressPackageStartupMessages(library(DMRichR, attach.required = T))

setwd("/share/lasallelab/Ben/DS_DBS/DMRs/")
testCovariate <- "Diagnosis"

contrasts <- c("DSvsTD", "DSvsDD", "DDvsTD")

mclapply(contrasts, function(contrast){
  
  load(glue::glue("{contrast}/Production/RData/bsseq.RData"))
  load(glue::glue("{contrast}/Production/RData/DMRs.RData"))
  filename <- glue::glue("./{contrast}_heatmap.pdf")
  
  pData(bs.filtered.bsseq)$Diagnosis <- pData(bs.filtered.bsseq)$Diagnosis %>% 
    dplyr::recode_factor("DownSyndrome" = "Down Syndrome",
                         "DevelopmentalDelay" = "Developmental Delay",
                         "Control" = "Typical Development")
  
  annotation_colors <- list(Diagnosis = c("Down Syndrome" = "#F8766D",
                                          "Developmental Delay" = "#00BA38",
                                          "Typical Development" = "#619CFF"),
                            Sex = c(Male = "#8DDBFF",
                                    Female = "#FEC0CB")
                            )
  
  sigRegions %>%
    smoothPheatmap2(bs.filtered.bsseq = bs.filtered.bsseq,
                    testCovariate = testCovariate,
                    annotation_colors = annotation_colors,
                    filename = filename)
  
},
mc.cores = length(contrasts),
mc.silent = TRUE)
