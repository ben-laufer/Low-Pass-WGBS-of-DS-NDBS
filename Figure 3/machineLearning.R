# WGBS machine learning analysis and heatmap
# By Hyeyeon Hwang

.libPaths("/share/lasallelab/programs/DMRichR/R_3.6")
suppressPackageStartupMessages(library(DMRichR, attach.required = T))
setwd("/share/lasallelab/Hyeyeon/projects/machine_learning_DS_DBS")
packages <- c("BSgenome.Hsapiens.UCSC.hg38", "TxDb.Hsapiens.UCSC.hg38.knownGene", "org.Hs.eg.db")
glue::glue("Loading {packages}")
stopifnot(suppressMessages(sapply(packages, require, character.only = TRUE)))

goi <- BSgenome.Hsapiens.UCSC.hg38
TxDb <- TxDb.Hsapiens.UCSC.hg38.knownGene
annoDb <- "org.Hs.eg.db"
testCovariate <- "Diagnosis"
genome <- "hg38"

# load RData
load("all_DMRs.RData")
regions <- sigRegions
load("all_bsseq.RData")
bsseq <- bs.filtered.bsseq

data <- getMeth(BSseq = bsseq, regions = regions, type = "smooth", what = "perRegion") %>% as.matrix() %>% t()
colnames(data) <- cbind(as.data.frame(seqnames(regions)), ranges(regions)) %>% 
  dplyr::as_tibble() %>% 
  dplyr::select(value, start, end) %>% 
  tidyr::unite("bed", c("value", "start", "end"), sep = ".") %>% 
  as.matrix()

# data_orig <- data

pData <- bsseq %>% 
  pData() %>% 
  dplyr::as_tibble() %>% 
  dplyr::mutate(Diagnosis = dplyr::case_when(Diagnosis == "DownSyndrome" ~ "Down Syndrome", 
                                             Diagnosis == "DevelopmentalDelay" ~ "Developmental Delay", 
                                             Diagnosis == "Control" ~ "Typical Development")) %>% 
  dplyr::mutate(col = dplyr::case_when(Diagnosis == "Down Syndrome" ~ "#F8766D", 
                                       Diagnosis == "Developmental Delay" ~ "#00BA38", 
                                       Diagnosis == "Typical Development" ~ "#619CFF")) %>% 
  dplyr::mutate(Diagnosis_2class = dplyr::case_when(Diagnosis == "Down Syndrome" ~ "Down Syndrome", 
                                                    Diagnosis == "Developmental Delay" ~ "Not Down Syndrome", 
                                                    Diagnosis == "Typical Development" ~ "Not Down Syndrome"))

pData <- pData %>% 
  dplyr::mutate_at(vars(Diagnosis, Diagnosis_2class), as.factor) %>% 
  tibble::add_column(Name = rownames(bsseq %>% pData()), .before = 1)

# reorder factor levels for legend order on heatmap
pData$Diagnosis <- factor(pData$Diagnosis, levels(pData$Diagnosis)[c(2, 1, 3)])

# for machine learning, use Diagnosis_2class (Down Syndrome, Not Down Syndrome)
data <- data %>% dplyr::as_tibble() %>% tibble::add_column(groups = pData$Diagnosis_2class, .before = 1) %>% 
  tibble::add_column(sampleID = pData$Name, .before = 1) %>% 
  tibble::add_column(Diagnosis_3class = pData$Diagnosis, .after = 2)


splitDmrs <- function(ranking) {
    DMRsplit <- ranking$DMR %>%
      strsplit(., split = "[.]") %>%
      as.data.frame() %>%
      t() %>%
      magrittr::set_colnames(c("chr", "start", "end")) %>%
      dplyr::as_tibble()
    
    ranking <- ranking %>%
      tibble::add_column(chr = DMRsplit$chr, .after = 2) %>%
      tibble::add_column(start = DMRsplit$start, .after = 3) %>%
      tibble::add_column(end = DMRsplit$end, .after = 4)
    
    return(ranking)
  }


getRfRanking <- function() {
    cat("\n", "Training random forest (RF) model for DMR ranking...")
    set.seed(5)
    borutaTrainObject <- Boruta::Boruta(groups ~ ., data = data %>% dplyr::select(-sampleID), doTrace = 0)  
    borutaTrainStats <- Boruta::attStats(borutaTrainObject)
    
    rfRanking <- tibble::tibble(DMR = rownames(borutaTrainStats), 
                                meanImp = borutaTrainStats$meanImp, 
                                decision = borutaTrainStats$decision) %>% 
      dplyr::arrange(dplyr::desc(meanImp)) %>% 
      tibble::add_column(Rank = 1:nrow(borutaTrainStats), .before = 1) 
    
    rfRanking <- rfRanking %>% splitDmrs()
    cat("Done.")
    return(rfRanking)
  }


getSvmRanking <- function() {
    cat("\n", "Training support vector machine (SVM) model for DMR ranking...")
    dataMatrix <- data %>% 
      dplyr::select(-c(groups, sampleID)) %>% 
      as.matrix()
    set.seed(5)
    sigfeatTrainObject <- sigFeature::sigFeature(dataMatrix, data$groups) 
    
    svmRanking <- tibble::tibble(Rank = 1:length(sigfeatTrainObject), 
                                 DMR = colnames(dataMatrix[, sigfeatTrainObject]))
    
    svmRanking <- svmRanking %>% splitDmrs()
    cat("Done.")
    return(svmRanking) 
  }  


getCommonDmrs <- function() {
    numPredictors <- ncol(data) - 2 
    numTopPercent <- ceiling(topPercent * .01 * numPredictors)
    allCommonFlag <- 0
    
    # Top percent of predictors is less than 10 AND number of predictors is at least 10
    if (numTopPercent < 10 && numPredictors >= 10) {
      commonDmrs <- intersect(rfRanking$DMR[1:10], svmRanking$DMR[1:10]) 
      cat("\n", glue::glue("{topPercent}% of {numPredictors} total DMRs is {numTopPercent} and less than 10."))
      cat("\n", glue::glue("Finding common DMRs in top 10 DMRs (instead of in top {topPercent}%) in RF and SVM predictor ranking lists."))
      case <- 10
      if(length(commonDmrs) == 10) { allCommonFlag <- 1 }
    } 
    # Top percent of predictors is less than 10 AND number of predictors is less than 10
    else if (numTopPercent < 10 && numPredictors < 10) {
      commonDmrs <- intersect(rfRanking$DMR[1:numPredictors], svmRanking$DMR[1:numPredictors]) 
      cat("\n", glue::glue("{topPercent}% of {numPredictors} total DMRs is {numTopPercent} and less than 10."))
      cat("\n", glue::glue("Finding common DMRs in top {numPredictors} DMRs (instead of in top {topPercent}%) in RF and SVM predictor ranking lists."))
      case <- numPredictors
      if(length(commonDmrs) == numPredictors) { allCommonFlag <- 1 }
    } 
    # Top percent of predictors is at least 10
    else {
      cat("\n", glue::glue("{topPercent}% of {numPredictors} total DMRs is {numTopPercent} and at least 10."))
      cat("\n", glue::glue("Finding common DMRs in top {topPercent}% of DMRs in RF and SVM predictor ranking lists."))
      commonDmrs <- intersect(rfRanking$DMR[1:numTopPercent], svmRanking$DMR[1:numTopPercent]) 
      case <- numTopPercent
      if(length(commonDmrs) == numTopPercent) { allCommonFlag <- 1 }
    }
    
    # No common DMRs
    if(length(commonDmrs) == 0) {
      cat("\n", glue::glue("There were 0 common DMRs. Rerun with a higher topPercent value for a greater number of common DMRs."))
    }
    
    # To find correct RF/SVM rank if all selected DMRs are the same in both RF and SVM lists
    if(allCommonFlag == 1) {
      commonDmrsRfRank <- numeric()
      commonDmrsSvmRank <- numeric()
      
      for(i in 1:length(commonDmrs)) {
        # rfRank 
        for (rfDmr in rfRanking$DMR) {
          if(commonDmrs[i] == rfDmr) {
            rfRank <- rfRanking$Rank[which(rfRanking$DMR == rfDmr)]
            commonDmrsRfRank <- commonDmrsRfRank %>% append(rfRank)
          }
        }
        # svmRank
        for (svmDmr in svmRanking$DMR) {
          if(commonDmrs[i] == svmDmr) {
            svmRank <- svmRanking$Rank[which(svmRanking$DMR == svmDmr)]
            commonDmrsSvmRank <- commonDmrsSvmRank %>% append(svmRank)
          }
        }
      }
    # RF / SVM rank otherwise
    } else {
      commonDmrsRfRank <- which(rfRanking$DMR %in% commonDmrs)
      commonDmrsSvmRank <- which(svmRanking$DMR %in% commonDmrs)
    }
    
    cat("\n")
    return(list(dmrs = commonDmrs, rfRank = commonDmrsRfRank , svmRank = commonDmrsSvmRank, case = case))
  }


annotateDmr <- function(dmrList, rfRank, svmRank, type) {
    cat(" Beginning DMR annotation...", "\n")
    
    annotatedDmrs <- dmrList %>% 
      # Set up DMR list before annotating 
      strsplit(., split = "[.]") %>%
      as.data.frame() %>%
      t() %>%
      magrittr::set_colnames(c("chr", "start", "end")) %>%
      dplyr::as_tibble() %>%
      # Annotate DMRs  
      GenomicRanges::makeGRangesFromDataFrame(ignore.strand = TRUE,
                                              seqnames.field = "chr",
                                              start.field = "start",
                                              end.field = "end") %>%
      ChIPseeker::annotatePeak(TxDb = TxDb,
                               annoDb = annoDb,
                               overlap = "all") %>%
      dplyr::as_tibble() %>%  
      # Add RF rank, SVM rank or both depending on "type"
      {if(type == "rf") tibble::add_column(., rank = rfRank, .before = 1) else .} %>%
      {if(type == "svm") tibble::add_column(., rank = svmRank, .before = 1) else .} %>%
      {if(type == "common") tibble::add_column(., RF.rank = rfRank, .before = 1) else .} %>%
      {if(type == "common") tibble::add_column(., SVM.rank = svmRank, .before = 2) else .} #%>%
      # Select only relevant columns
      # dplyr::select(-c(strand, geneChr, geneStart, geneEnd, geneStrand, geneId, transcriptId, distanceToTSS, ENSEMBL))
    
    return(annotatedDmrs)  
  }

# contains 3 class diagnosis labels
data_alldiag <- data
data <- data %>% dplyr::select(-Diagnosis_3class)

rfRanking <- getRfRanking()
svmRanking <- getSvmRanking()
topPercent = 1
commonDmrs <- getCommonDmrs()
annotatedCommonDmrs <- annotateDmr(dmrList = commonDmrs$dmrs, rfRank = commonDmrs$rfRank, svmRank = commonDmrs$svmRank, type = "common")

# Generate heatmap
heatmapData <- data[, which(colnames(data) %in% commonDmrs$dmrs)] %>% t()
colnames(heatmapData) <- data$sampleID
annot_col <-  data.frame(testCovariate = data_alldiag$Diagnosis_3class)
colnames(annot_col) <- testCovariate
rownames(annot_col) <- colnames(heatmapData)

annot_col2 <- annot_col %>% tibble::add_column(Sex = pData$Sex)
annot_colors2 = list(Diagnosis = c('Down Syndrome' = '#F8766D', 
                                   'Developmental Delay' = '#00BA38',
                                   'Typical Development' = '#619CFF'), 
                     Sex = c('Male' = '#8DDBFF', 'Female' = '#FEC0CB'))


heatmap_labs <- annotatedCommonDmrs[, c(3, 4, 5, 8, 18)]

heatmap_rownames_dmr <- rownames(heatmapData)
heatmap_dmr <- heatmap_labs %>% dplyr::as_tibble() %>% dplyr::select(seqnames, start, end) %>% tidyr::unite("dmr", c("seqnames", "start", "end"), sep = ".")
heatmap_labs <- heatmap_labs %>% tibble::add_column(dmr = heatmap_dmr, .before = 1)
heatmap_labs <- heatmap_labs[match(heatmap_rownames_dmr, heatmap_labs$dmr$dmr), ]

heatmap_labs <- heatmap_labs %>% 
  dplyr::mutate(labels = dplyr::case_when(stringr::str_detect(annotation, "Intron") ~ "Intron", 
                                          stringr::str_detect(annotation, "Intergenic") ~ "Intergenic", 
                                          stringr::str_detect(annotation, "Promoter") ~ "Promoter", 
                                          stringr::str_detect(annotation, "Exon") ~ "Exon", 
                                          stringr::str_detect(annotation, "3' UTR") ~ "3' UTR"))


labels_symbols_annot <- heatmap_labs %>% dplyr::select(SYMBOL, labels) %>% tidyr::unite("final_labels", c("SYMBOL", "labels"), sep = " ")

# set rownames of heatmapData as row labels to display on the heatmap
rownames(heatmapData) <- labels_symbols_annot$final_labels

# Create heatmap
pheatmap::pheatmap(mat = heatmapData, scale = "row", show_colnames = F, border_color = "black", 
                   main = glue::glue("Z-Scores of {nrow(heatmapData)} Differentially Methylated Regions"), 
                   annotation_col = annot_col2, fontsize = 16,
                   color = rev(RColorBrewer::brewer.pal(11, name = "RdBu")),
                   filename = "DS_DBS_machine_learning_heatmap.pdf", 
                   width = 30, height = 10, 
                   annotation_colors = annot_colors2)


