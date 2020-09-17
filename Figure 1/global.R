# global
# Process all samples for equal CpG coverage and smooth, test for global methylation differences, and create a single CpG density plot
# By Ben Laufer

# Load packages -----------------------------------------------------------

.libPaths("/share/lasallelab/programs/DMRichR/R_3.6")
packages <- c("DMRichR", "broom", "tidyverse")
stopifnot(suppressMessages(sapply(packages, require, character.only = TRUE)))

setwd("/share/lasallelab/Ben/DS_DBS/DMRs/cytosine_reports")

packages <- c("BSgenome.Hsapiens.UCSC.hg38", "TxDb.Hsapiens.UCSC.hg38.knownGene", "org.Hs.eg.db")
glue::glue("Loading {packages}")
stopifnot(suppressMessages(sapply(packages, require, character.only = TRUE)))

goi <- BSgenome.Hsapiens.UCSC.hg38
TxDb <- TxDb.Hsapiens.UCSC.hg38.knownGene
annoDb <- "org.Hs.eg.db"

testCovariate <- "Diagnosis"
genome <- "hg38"

# Load and process samples ------------------------------------------------

bs.filtered <- processBismark(files = list.files(path = getwd(), pattern = "*.txt.gz"),
                              meta = openxlsx::read.xlsx("sample_info_all.xlsx", colNames = TRUE) %>%
                                dplyr::mutate_if(is.character, as.factor),
                              testCovar = "Diagnosis",
                              adjustCovar = "Sex",
                              matchCovar = NULL,
                              Cov = 1,
                              mc.cores = 30,
                              per.Group = 0.75)

glue::glue("Saving Rdata...")
dir.create("cytosine_reports/consensus/")
save(bs.filtered, file = "consensus/consensus_bismark.RData")

# Individual smoothed values ----------------------------------------------

cat("\n[DMRichR] Smoothing individual methylation values \t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
start_time <- Sys.time()

bs.filtered.bsseq <- BSmooth(bs.filtered,
                             BPPARAM = MulticoreParam(workers = cores,
                                                      progressbar = TRUE)
)

bs.filtered.bsseq

glue::glue("Saving Rdata...")
save(bs.filtered.bsseq, file = "consensus/consensus_bsseq.RData")

glue::glue("Individual smoothing timing...")
end_time <- Sys.time()
end_time - start_time

# Tidy diagnosis ----------------------------------------------------------

pData(bs.filtered.bsseq)$Diagnosis <- pData(bs.filtered.bsseq)$Diagnosis %>% 
  dplyr::recode_factor("Control" = "Typical Development",
                       "DevelopmentalDelay" = "Developmental Delay",
                       "DownSyndrome" = "Down Syndrome"
  )

# ANOVA -------------------------------------------------------------------

cat("Testing for global methylation differences...")
global <- data.frame(DelayedMatrixStats::colMeans2(getMeth(BSseq = bs.filtered.bsseq, type = "smooth", what = "perBase"), na.rm = TRUE))
global$sample <- sampleNames(bs.filtered.bsseq)
names(global) <- c("CpG_Avg", "sample")
global <- dplyr::as_tibble(cbind(global, data.frame(pData(bs.filtered.bsseq))), rownames = NULL)

summary <- global %>%
  group_by(Diagnosis, Sex) %>%
  summarise(
    count = n(),
    mean = mean(CpG_Avg, na.rm = TRUE),
    sd = sd(CpG_Avg, na.rm = TRUE)
  )

ANOVA <- global %>%
  aov(CpG_Avg ~ Diagnosis + Sex, data = .)

Tukey <- ANOVA %>%
  TukeyHSD(which = "Diagnosis")

pairWise <- lm(CpG_Avg ~ Diagnosis + Sex, data = global) %>% 
  lsmeans::ref.grid(data = global) %>%
  lsmeans::lsmeans(~Diagnosis) %>%
  pairs(reverse = TRUE) %>%
  summary()

list("globalInput" = global,
     "Summary" = summary,
     "ANOVA" = broom::tidy(ANOVA),
     "Tukey" = broom::tidy(Tukey),
     "pairWise" = pairWise
     )%>%
  openxlsx::write.xlsx("Global/globalStats.xlsx")

cat("Done", "\n")

# Single CpG Density Plot -------------------------------------------------

densityPlot <- function(bsseq = bs.filtered.bsseq,
                        group = NA){
  print(glue::glue("[DMRichR] Density plot of {nrow(bsseq)} CpGs"))
  bs.filtered.bsseq %>% 
    bsseq::getMeth(BSseq = .,
                   type = "smooth",
                   what = "perBase") %>%
    dplyr::as_tibble() %>% 
    na.omit() %>%
    magrittr::set_colnames(paste(group, seq_along(1:length(group)))) %>%
    dplyr::transmute(Group1 = dplyr::select(., dplyr::contains(levels(group)[1])) %>% rowMeans()*100,
                     Group2 = dplyr::select(., dplyr::contains(levels(group)[2])) %>% rowMeans()*100,
                     Group3 = dplyr::select(., dplyr::contains(levels(group)[3])) %>% rowMeans()*100
    ) %>%
    magrittr::set_colnames(c(levels(group)[1], levels(group)[2], levels(group)[3])) %>% 
    tidyr::gather(key = "variable",
                  value = "value") %>%
    dplyr::mutate(variable = factor(.$variable)) %>% 
    dplyr::mutate(variable = factor(.$variable, levels = unique(forcats::fct_rev(group)))) %>% 
    ggplot(aes(value, color = variable)) +
    geom_density(size = 1.3) +
    labs(x = "Percent Methylation", y = "Density", color = "Diagnosis") +
    theme_classic() +
    scale_x_continuous(expand=c(0.05,0.05), breaks = c(0,25,50,75,100)) +
    scale_y_continuous(expand=c(0.0,0.001)) +
    theme(axis.text = element_text(size = 20), axis.title = element_text(size = 20),
          strip.text = element_text(size = 20), legend.text = element_text(size = 16),
          legend.position = "bottom", legend.title = element_text(size = 20)) %>%
    return()
}

bs.filtered.bsseq %>%
  densityPlot(group = bs.filtered.bsseq %>%
                pData() %>%
                dplyr::as_tibble() %>%
                dplyr::pull(!!testCovariate)
  ) %>% 
  ggplot2::ggsave("Global/Smoothed CpG Density Plot.pdf",
                  plot = .,
                  device = NULL,
                  width = 10,
                  height = 6)
