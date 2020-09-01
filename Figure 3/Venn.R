# Create an area-proportional Venn plot based on genomic coordinates
# By Ben Laufer

# Packages ----------------------------------------------------------------

rm(list = ls())

.libPaths("/share/lasallelab/programs/DMRichR/R_3.6")
# install.packages("Vennerable", repos="http://R-Forge.R-project.org", type="source")
packages <-  c("DMRichR","ChIPpeakAnno", "Vennerable")
stopifnot(suppressMessages(sapply(packages, require, character.only = TRUE)))

setwd("/share/lasallelab/Ben/DS_DBS/DMRs/")

# Load DMRs ---------------------------------------------------------------

load("DSvsTD/Production/RData/DMRs.RData")
DSvsTD <- sigRegions
rm(sigRegions, regions)

load("DSvsDD/Production/RData/DMRs.RData")
DSvsDD <- sigRegions
rm(sigRegions, regions)

load("DDvsTD/Production/RData/DMRs.RData")
DDvsTD <- sigRegions
rm(sigRegions, regions)

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
  n <- which(colnames(venn_cnt)=="Counts") - 1
  SetNames=colnames(venn_cnt)[1:n]
  Weight=venn_cnt[,"Counts"]
  names(Weight) <- apply(venn_cnt[,1:n], 1, paste, collapse="")
  Vennerable::Venn(SetNames=SetNames, Weight=Weight)
}
v <- venn_cnt2venn(res$vennCounts)

svg("Overlaps/Venn.svg",
    height = 8.5,
    width = 12)

plot(v)

dev.off()

# Manually edit in inkscape after
