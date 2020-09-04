# GOplot
# Slim and compare GOfuncR terms from one comparison to other comparisons
# By Ben Laufer

#' REVIGO2
#' @description Slims top signficant Gene Ontology terms from enrichR, rGREAT, and GOfuncR.
#' The terms are ranked by dispensability.
#' @param GO A dataframe or list of dataframes returned
#' from \code{enrichR::enrichr()}, \code{rGREAT::getEnrichmentTables()}, or \code{GOfuncR::go_enrich()}.
#' @param tool A character vector of the name of the database (enrichR, rGREAT, or GOfuncR).
#' @return A \code{tibble} of top distinct and significant GO terms from an \code{enrichR} 
#' or \code{rGREAT} analysis.
#' @import enrichR
#' @import rGREAT
#' @import GOfuncR
#' @importFrom magrittr %>%
#' @importFrom dplyr filter as_tibble mutate select group_by slice ungroup
#' @importFrom rvest html_session html_form set_values submit_form
#' @importFrom data.table rbindlist
#' @importFrom glue glue
#' @references \url{https://github.com/hbc/revigoR/blob/master/rvest_revigo.R}
#' @references \url{http://revigo.irb.hr}
#' @export REVIGO2
REVIGO2 <- function(GO = GO,
                    tool = c("enrichR", "rGREAT", "GOfuncR")){
  
  print(glue::glue("Tidying results from {tool}..."))
  if(tool == "enrichR"){
    GO <- GO %>%
      data.table::rbindlist(idcol = "Database") %>%
      dplyr::filter(Database %in% c("GO_Biological_Process_2018",
                                    "GO_Cellular_Component_2018",
                                    "GO_Molecular_Function_2018"
      )
      ) %>% 
      dplyr::as_tibble() %>%
      dplyr::mutate(Term = stringr::str_extract(.$Term, "\\(GO.*")) %>%
      dplyr::mutate(Term = stringr::str_replace_all(.$Term, "[//(//)]",""), "") %>%
      dplyr::filter(P.value <= 0.05)
    
    goList <- paste(GO$Term, GO$P.value, collapse = "\n")
    
  }else if(tool == "rGREAT"){
    
    GO <-  GO %>%
      data.table::rbindlist(idcol = "Database") %>%
      dplyr::as_tibble() %>%
      dplyr::filter(Hyper_Raw_PValue <= 0.05)
    
    goList <- paste(GO$ID, GO$Hyper_Raw_PValue, collapse = "\n")
    
  }else if(tool == "GOfuncR"){
    
    GO <- GO$results %>%
      dplyr::filter(raw_p_overrep <= 0.05)
    
    goList <- paste(GO$node_id, GO$raw_p_overrep, collapse = "\n")
    
  }else{
    stop(glue("{tool} is not supported, please choose either enrichR, rGREAT, or GOfuncR [Case Sensitive]"))
  }
  
  print(glue::glue("Submiting results from {tool} to REVIGO..."))
  revigo_session <- rvest::html_session("http://revigo.irb.hr/")
  revigo_form <- rvest::html_form(revigo_session)[[1]]  
  filled_form <- rvest::set_values(revigo_form,
                                   'goList' = goList,
                                   'cutoff' = 0.4,
                                   'isPValue' = "yes",
                                   'measure' = "SIMREL")
  result_page <- rvest::submit_form(revigo_session,
                                    filled_form,
                                    submit = 'startRevigo')
  
  revigo_results <- list()
  for (i in 1:3){
    results_table <- rvest::html_table(result_page)[[i]]
    names(results_table) <- results_table[2,]
    revigo_results[[i]] <- results_table[3:nrow(results_table),]
  }
  names(revigo_results) <- c("Biological Process", "Cellular Component", "Molecular Function")
  
  revigo_results <- data.table::rbindlist(revigo_results, idcol = "Gene Ontology") %>%
    dplyr::filter(dispensability < 0.4) %>%
    dplyr::select("Gene Ontology",
                  Term = "description",
                  "term ID",
                  "-log10.p-value" = `log10 p-value`) %>% 
    dplyr::mutate("-log10.p-value" = -(as.numeric(`-log10.p-value`))) %>% 
    dplyr::mutate("Gene Ontology" = as.factor(`Gene Ontology`)) %>% 
    dplyr::group_by(`Gene Ontology`) %>%
    dplyr::slice(1:7) %>%
    dplyr::ungroup() %>% 
    return()
}

#' GOplot2
#' @description Slims and plots top signficant Gene Ontology terms from enrichR, rGREAT, and GOfuncR.
#' The terms are ranked by dispensability before being plotted, which then orders them by p-value. 
#' @param revigoResults A \code{tibble} from\code{DMRichR::REVIGO()}.
#' @return A \code{ggplot} object of top significant GO and pathway terms from an \code{enrichR} 
#' or \code{rGREAT} analysis that can be viewed by calling it, saved with \code{ggplot2::ggsave()}, 
#' or further modified by adding \code{ggplot2} syntax.
#' @import ggplot2
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate select
#' @importFrom forcats fct_rev
#' @importFrom ggsci scale_fill_d3
#' @importFrom Hmisc capitalize
#' @importFrom glue glue
#' @export GOplot2
#' 
GOplot2 <- function(revigoResults = revigoResults){
  
  print(glue::glue("Plotting slimmed gene ontology results"))
  
  revigoResults %>%
    dplyr::mutate(Term = stringr::str_trim(.$Term)) %>%
    dplyr::mutate(Term = Hmisc::capitalize(.$Term)) %>%
    dplyr::mutate(Term = stringr::str_wrap(.$Term, 45)) %>% 
    dplyr::mutate(Database = factor(.$Database)) %>% 
    dplyr::mutate(Term = factor(.$Term, levels = unique(.$Term[order(forcats::fct_rev(.$Database), .$`-log10.p-value`)]))) %>% 
    ggplot2::ggplot(aes(x = Term, y = `-log10.p-value`, fill = Database, group = Database)) +
    geom_bar(stat = "identity", position = position_dodge(), color = "Black") +
    facet_grid(~contrast) +
    geom_hline(yintercept = 0) +
    coord_flip() +
    scale_y_continuous(expand = c(0, 0)) +
    ggsci::scale_fill_d3(name = "Gene Ontology") +
    labs(y = expression("-log"[10](p))) +
    theme_classic() +
    theme(axis.text = element_text(size = 14),
          axis.title.y = element_blank(),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 12),
          legend.position = "bottom",
          strip.text = element_text(size = 14)
    )
}
 
# Run ---------------------------------------------------------------------

packages <- c("tidyverse", "DMRichR")
stopifnot(suppressMessages(sapply(packages, require, character.only = TRUE)))
options(readr.num_columns = 0)

setwd("/Users/benlaufer/Box/Down Syndrome/Dried Blood Spots/Results/DMRs")

# Get slimmed terms for DS vs TD

path <- "DSvsTD/Ontologies/All DMRs/GOfuncR.xlsx"

slim <- path %>%
  readxl::excel_sheets() %>%
  set_names() %>%
  map(readxl::read_excel, path = path) %>%
  REVIGO2(tool = "GOfuncR") %>%
  dplyr::pull("term ID") %>%
  as.character()

# Tidy GOfuncR analyses

tidy <- . %>% 
  dplyr::mutate(`-log10.p-value` = -log10(raw_p_overrep)) %>%
  dplyr::select(node_name, node_id, ontology, `-log10.p-value`) 

DSvsTD <- readxl::read_xlsx("DSvsTD/Ontologies/All DMRs/GOfuncR.xlsx") %>%
  tidy()

DSvsDD <- readxl::read_xlsx("DSvsDD/Ontologies/All DMRs/GOfuncR.xlsx") %>%
  tidy()

DDvsTD <- readxl::read_xlsx("DDvsTD/Ontologies/All DMRs/GOfuncR.xlsx") %>%
  tidy()

# Get top slimmed DS vs TD terms from DS vs DD and DD vs TD and plot results

GO <- DSvsTD %>%
  dplyr::filter(node_id %in% slim) %>% 
  dplyr::left_join(DSvsDD, by = c("node_name", "node_id", "ontology")) %>%
  dplyr::left_join(DDvsTD, by = c("node_name", "node_id", "ontology")) %>%
  dplyr::rename(DSvsTD = "-log10.p-value.x",
                DSvsDD = "-log10.p-value.y",
                DDvsTD = "-log10.p-value"
                ) %>%
  dplyr::mutate(ontology = dplyr::recode_factor(ontology,
                                                "biological_process" = "Biological Process",
                                                "cellular_component" = "Cellular Component",
                                                "molecular_function" = "Molecular Function"
                                                )
                ) %>% 
  dplyr::group_by(ontology) %>%
  dplyr::slice(1:7) %>%
  dplyr::ungroup() %>% 
  dplyr::rename(Term = node_name,
                Database = ontology,
                ID = node_id) %>% 
  tidyr::pivot_longer(cols = c(-Term, -Database, -ID),
                      names_to = "contrast",
                      values_to = "-log10.p-value"
                      ) %>%
  dplyr::mutate(contrast = dplyr::recode_factor(contrast,
                                                "DSvsTD" = "DS vs TD",
                                                "DSvsDD" = "DS vs DD",
                                                "DDvsTD" = "DD vs TD"
                                                )
                ) %>%
  dplyr::mutate(`-log10.p-value` = ifelse(is.na(`-log10.p-value`), 0, `-log10.p-value`)) %>%
  GOplot2() %>% 
  ggplot2::ggsave("GOfuncR_plot.pdf",
                  plot = .,
                  device = NULL,
                  height = 8.5,
                  width = 8.5)

