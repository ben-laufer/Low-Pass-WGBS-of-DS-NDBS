# venn
# Create an area-proportional Venn diagram based on overlaps of genomic coordinates from multiple comparisons
# By Ben Laufer

# Packages ----------------------------------------------------------------

.libPaths("/share/lasallelab/programs/DMRichR/R_3.6")

if (!requireNamespace("Vennerable", quietly = TRUE))
  BiocManager::install("js229/Vennerable")
  
packages <- c("DMRichR", "ChIPpeakAnno", "Vennerable", "svglite")
stopifnot(suppressMessages(sapply(packages, require, character.only = TRUE)))

# Load DMRs ---------------------------------------------------------------

setwd("/share/lasallelab/Ben/DS_DBS/DMRs/")

loadDMRs <- function(name){
  load(glue::glue("{name}/Production/RData/DMRs.RData"))
  assign(name, sigRegions, envir = .GlobalEnv)
  rm(sigRegions, regions)
}

contrasts <- c("DSvsTD", "DSvsDD", "DDvsTD")
  
purrr::walk(contrasts, loadDMRs)

# Plot Venn ---------------------------------------------------------------

dir.create("Overlaps")

res <- ChIPpeakAnno::makeVennDiagram(Peaks = list(DSvsTD,
                                                  DSvsDD,
                                                  DDvsTD
                                                  ), 
                                     NameOfPeaks = c("DSvsTD",
                                                     "DSvsDD",
                                                     "DDvsTD"
                                                     )
                                     )

# ref: https://support.bioconductor.org/p/67429/
venn_cnt2venn <- function(venn_cnt){
  n <- which(colnames(venn_cnt) == "Counts") - 1
  SetNames <- colnames(venn_cnt)[1:n]
  Weight <- venn_cnt[,"Counts"]
  names(Weight) <- apply(venn_cnt[,1:n], 1, paste, collapse="")
  Vennerable::Venn(SetNames = SetNames, Weight = Weight)
}
v <- venn_cnt2venn(res$vennCounts)

svg("Overlaps/Venn.svg",
    height = 8.5,
    width = 12)

plot(v)

dev.off()

# Manually edit in inkscape after
