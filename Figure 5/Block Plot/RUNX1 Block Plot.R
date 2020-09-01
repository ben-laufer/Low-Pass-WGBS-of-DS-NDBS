# RUNX1 block plot
# Run conensusBlocks.sh first
# By Ben Laufer

# Load packages -----------------------------------------------------------

.libPaths("/share/lasallelab/programs/DMRichR/R_3.6")
suppressPackageStartupMessages(library(DMRichR, attach.required = T))

setwd("/share/lasallelab/Ben/DS_DBS/DMRs")

# Load data ---------------------------------------------------------------

load("cytosine_reports/consensus/consensus_bismark.RData")
load("cytosine_reports/consensus/consensus_bsseq.RData")

packages <- c("BSgenome.Hsapiens.UCSC.hg38", "TxDb.Hsapiens.UCSC.hg38.knownGene", "org.Hs.eg.db")
glue::glue("Loading {packages}")
stopifnot(suppressMessages(sapply(packages, require, character.only = TRUE)))
goi <- BSgenome.Hsapiens.UCSC.hg38
TxDb <- TxDb.Hsapiens.UCSC.hg38.knownGene
annoDb <- "org.Hs.eg.db"

testCovariate <- "Diagnosis"
genome <- "hg38"

# Annotate pData ----------------------------------------------------------

pData <- bs.filtered.bsseq %>%
  pData() %>%
  dplyr::as_tibble() %>%
  dplyr::mutate(Diagnosis = dplyr::case_when(Diagnosis == "DownSyndrome" ~ "Down Syndrome",
                                             Diagnosis == "DevelopmentalDelay" ~ "Developmental Delay",
                                             Diagnosis == "Control" ~ "Typical Development"
                                                   )
                ) %>%
  dplyr::mutate(col = dplyr::case_when(Diagnosis == "Down Syndrome" ~ "firebrick3",
                                       Diagnosis == "Developmental Delay" ~ "forestgreen",
                                       Diagnosis == "Typical Development" ~ "mediumblue"
                                       )
                )

pData(bs.filtered) <- pData 
pData(bs.filtered.bsseq) <- pData 

# Plot --------------------------------------------------------------------

blocks <- readr::read_tsv("cytosine_reports/consensus/ConsensusBlocks.bed",
                            col_names = c("chr", "start", "end")
                          ) %>% 
  GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE) 

glue::glue("Plotting blocks...")

pdf("panTissue/panTissue_Blocks.pdf", height = 7.50, width = 11.50)

dmrseq::plotDMRs(bs.filtered.bsseq,
                 regions = blocks,
                 testCovariate = testCovariate,
                 annoTrack = getAnnot(genome),
                 regionCol = "#FF00001A",
                 qval = FALSE,
                 stat = FALSE)

dev.off()
