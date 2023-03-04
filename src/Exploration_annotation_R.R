#' ---
#' title: "Exploration annotation file"
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
  library(ggpubr)
  library(ggtree)
  library(here)
  library(LSD)
  library(stringr)
  library(tidyr)
})




#' # Import data
#' ```{r import annotation file, echo=TRUE}
#' ```

spruce_pine <- read.table(gzfile(here("reference/gff3/Picab02_codingAll.gff3.gz")), sep = "\t", quote="")

header <- c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes")
colnames(spruce_pine) <- header

genes <- subset(spruce_pine, type=="gene")
introns <- subset(spruce_pine, type=="intron")  
cds <- subset(spruce_pine, type=="CDS")
mrna <- subset(spruce_pine, type=="mRNA")


#' # Compare the cumulative intron length with the number of introns
#' ```{r comparison intron length vs # of intorns, echo=FALSE}
#' ```

parent <- sub(".*=","",introns$attributes)
table(parent)[which.max(table(parent))]
  
introns$parent <- parent
introns$length <- introns$end-introns$start+1
introns$counts <- 1
  

df_introns <- introns %>% group_by(parent)
  
df_introns_cum <- df_introns %>% summarise(
  counts= sum(counts),
  seqid=unique(seqid),
  length=sum(length),
  parent=unique(parent)
)
  
comparisonplot(df_introns_cum$counts, log10(df_introns_cum$length), xlab ="number of introns",
    ylab = "cumulative intron length (log scale)", main="Comparison plot of number of introns vs intron length", 
    cor=T, ylim=c(0,6))
