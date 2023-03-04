#' ---
#' title: "Create the tx2gene mapping"
#' author: "Emanuela Damieto"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#'    code_folding: hide
#' ---
#' # Setup
#' * Libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(here)
  library(readr)
  library(stringr)
  library(tidyr)
})

#' * Data
read_tsv(here("reference/gff3/agat_renamed_with_pos_short_removed_eggnog_added.gff.gz"),
         comment="#",col_names=FALSE,show_col_types=FALSE) %>% 
  filter(X3 %in% c("mRNA")) %>% 
  select(X9) %>% separate(X9,c("TXID","GENEID"),sep=";",extra="drop") %>% 
  mutate_all(str_replace_all,".*=","") %>% 
  write_tsv(here("reference/annotation/tx2gene_update.tsv.gz"))

#' # Session Info
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```