# Modify DMR heatmaps
# By Ben Laufer

# Load packages -----------------------------------------------------------

rm(list = ls())
.libPaths("/share/lasallelab/programs/DMRichR/R_3.6")
suppressPackageStartupMessages(library(DMRichR, attach.required = T))

setwd("/share/lasallelab/Ben/DS_DBS/DMRs/")
testCovariate <- "Diagnosis"

# Functions ---------------------------------------------------------------

smoothPheatmap <- function(regions = sigRegions,
                           bsseq = bs.filtered.bsseq,
                           testCovariate = testCovariate,
                           annotation_colors = annotation_colors,
                           filename = "DMRs.pdf",
                           ...){
  cat("\n[DMRichR] DMR heatmap \t\t\t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
  
  getMeth(BSseq = bsseq,
          regions = regions,
          type = "smooth",
          what = "perRegion") %>% 
    as.matrix() %>%
    pheatmap::pheatmap(.,
                       scale = "row",
                       annotation_col =  pData(bs.filtered.bsseq) %>%
                         as.data.frame() %>%
                         dplyr::select_if(~ nlevels(.) > 1),
                       color = RColorBrewer::brewer.pal(11,
                                                        name = "RdBu") %>%
                         rev(),
                       show_colnames = F,
                       #angle_col = 45,
                       border_color = "grey",
                       #breaks = c(-3,-2,-1,0,1,2,3),
                       main = glue::glue("Z-Scores of {length(regions)} Differentially Methylated Regions"),
                       fontsize = 16,
                       filename = filename,
                       width = 11,
                       height = 8.5,
                       annotation_colors = annotation_colors,
                       ...
    ) %>%
    return()
}

# DSvsTD ------------------------------------------------------------------

load("DSvsTD/Production/RData/bsseq.RData")
load("DSvsTD/Production/RData/DMRs.RData")
filename <- "./cytosine_reports/DSvsTD_heatmap.pdf"

pData(bs.filtered.bsseq)$Diagnosis <- pData(bs.filtered.bsseq)$Diagnosis %>% 
  dplyr::recode_factor(.,
                       "DownSyndrome" = "Down Syndrome",
                       "DevelopmentalDelay" = "Developmental Delay",
                       "Control" = "Typical Development"
  )

annotation_colors <- list(Diagnosis = c("Down Syndrome" = "#F8766D",
                                        "Typical Development" = "#619CFF"),
                          Sex = c(Male = "#8DDBFF",
                                  Female = "#FEC0CB")
)

sigRegions %>%
  smoothPheatmap(bsseq = bs.filtered.bsseq,
                 testCovariate = testCovariate,
                 annotation_colors = annotation_colors,
                 filename = filename)


rm(bs.filtered.bsseq, sigRegions, regions, filename, annotation_colors)

# DSvsDD ------------------------------------------------------------------

load("DSvsDD/Production/RData/bsseq.RData")
load("DSvsDD/Production/RData/DMRs.RData")
filename <- "./cytosine_reports/DSvsDD_heatmap.pdf"

annotation_colors <- list(Diagnosis = c("Down Syndrome" = "#F8766D",
                                        "Developmental Delay" = "#00BA38"),
                          Sex = c(Male = "#8DDBFF",
                                  Female = "#FEC0CB")
                          )

pData(bs.filtered.bsseq)$Diagnosis <- pData(bs.filtered.bsseq)$Diagnosis %>% 
  dplyr::recode_factor(.,
                       "DownSyndrome" = "Down Syndrome",
                       "DevelopmentalDelay" = "Developmental Delay",
                       "Control" = "Typical Development"
  )

sigRegions %>%
  smoothPheatmap(bsseq = bs.filtered.bsseq,
                 testCovariate = testCovariate,
                 annotation_colors = annotation_colors,
                 filename = filename)

rm(bs.filtered.bsseq, sigRegions, regions, filename, annotation_colors)

# DDvsTD ------------------------------------------------------------------

load("DDvsTD/Production/RData/bsseq.RData")
load("DDvsTD/Production/RData/DMRs.RData")
filename <- "./cytosine_reports/DDvsTD_heatmap.pdf"

annotation_colors <- list(Diagnosis = c("Developmental Delay" = "#00BA38",
                                        "Typical Development" = "#619CFF"),
                          Sex = c(Male = "#8DDBFF",
                                  Female = "#FEC0CB")
)

pData(bs.filtered.bsseq)$Diagnosis <- pData(bs.filtered.bsseq)$Diagnosis %>% 
  dplyr::recode_factor(.,
                       "DownSyndrome" = "Down Syndrome",
                       "DevelopmentalDelay" = "Developmental Delay",
                       "Control" = "Typical Development"
  )

sigRegions %>%
  smoothPheatmap(bsseq = bs.filtered.bsseq,
                 testCovariate = testCovariate,
                 annotation_colors = annotation_colors,
                 filename = filename)

rm(bs.filtered.bsseq, sigRegions, regions, filename, annotation_colors)
