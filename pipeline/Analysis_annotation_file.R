#' ---
#' title: "Analysis of the annotation file"
#' author: "Emanuela Damieto"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    fig_width: 9
#'    fig_height: 6
#'    toc: true
#'    number_sections: true
#'    toc_depth: 4
#'    toc_float:
#'      collapsed: TRUE
#'      smooth_scroll: TRUE
#'    code_folding: hide
#'    theme: "flatly"
#'    highlight: pygments
#' ---
#' <p style="text-align: left;"><span style="color: #bd047c;"><em>emanuela.damieto@studenti.unitn.it</em></span></p>
#' 
#' # Import the file 
#' 
#' * Libraries
suppressPackageStartupMessages({
  library(gplots)
  library(here)
  library(emoji)
  library(GenomicRanges)
  library(dplyr)
})

#' ```{r import the annotation file,eval=FALSE,echo=FALSE}
#' Import the gff3 (general feature format) file and retrieve the subset of objects that we are interested in
#' ```
#spruce <- read.table(gzfile(here("reference/gff3/spruce_protein_coding_AGAT_renamed_positions.gff.gz")), sep = "\t")
spruce <- read.table(gzfile(here("reference/gff3/agat_renamed_with_pos_short_removed_eggnog_added.gff.gz")), sep = "\t", quote="")


header <- c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes")
colnames(spruce) <- header

introns <- subset(spruce, type=="intron")  
cds <- subset(spruce, type=="CDS")
mrna <- subset(spruce, type=="mRNA")



#' * Retrieve the transcript with the highest number of introns 
parent <- sub(".*=","",introns$attributes)
table_parents <- table(parent)
table_parents[which.max(table_parents)]
mrna[str_detect(mrna$attributes, "PICAB_chr02_M00000376376"),]
#' `r emoji("point_right")` **The transcript with the highest number of introns is PICAB_chr02_M00000386998 and it has 74 introns**

#' * Create a new dataset with cumulative length of introns per transcript 
introns$parent <- parent
introns$length <- introns$end-introns$start+1
introns$counts <- 1

ID <- sub(";.*","",mrna$attributes)
ID <- sub("ID=","",ID)
mrna$ID <- ID

df_introns_in_transcr <- introns %>% group_by(parent)

df_introns_in_transcr <- df_introns_in_transcr %>% summarise(
  counts= sum(counts),
  seqid=unique(seqid),
  length=sum(length),
  parent=unique(parent)
)

df_introns_in_transcr[which.max(df_introns_in_transcr$length),]
introns[introns$parent=="PICAB_chr04_M00001122551",]

#' * Some plots
ggplot(df_introns_in_transcr, aes(x=counts, y=length))+
  geom_point()+
  labs(x="# introns", y="Intron length", title="Scatterplot of intron length and occurrence")

ggplot(df_introns_in_transcr, aes(x=counts, y=length))+
  geom_point()+
  scale_y_log10()+
  labs(x="# introns", y="Intron length", title="Scatterplot of log intron length and occurrence")

comparisonplot(df_introns_in_transcr$counts, df_introns_in_transcr$length, xlab ="# introns", 
               ylab = "Intron length", main="Comparison plot of number of introns vs intron length", xlim=c(0,100))

comparisonplot(df_introns_in_transcr$counts, log10(df_introns_in_transcr$length), xlab ="# introns", 
               ylab = "log intron length", main="Comparison plot of number of introns vs log intron length", ylim=c(0,6))

comparisonplot(df_introns_in_transcr$counts, df_introns_in_transcr$length, xlab ="# introns", 
               ylab = "log intron length", log="y", main="Comparison plot of number of introns vs log intron length", xlim=c(0,100))


#' * Repeat with weighted intron length (longest/cumulative intron length)
introns_grouped <- introns %>% group_by(parent)
df_introns_in_transcr_weight <- introns_grouped %>% summarise(
  counts= sum(counts),
  seqid=unique(seqid),
  length=max(length)/sum(length),
  parent=unique(parent)
)

comparisonplot(df_introns_in_transcr_weight$counts, df_introns_in_transcr_weight$length, xlab ="# introns", 
               ylab = "weighted intron length", main="Comparison plot of number of introns vs weighted intron length", xlim=c(0,100), ylim=c(0,1))

ggplot(df_introns_in_transcr_weight, aes(x=counts, y=length))+
  geom_point()+
  labs(x="# introns", y="Intron length", title="Scatterplot of weighted intron length and occurrence")

#' * Split long and short introns based on the cumulative length 
long_introns <- df_introns_in_transcr %>% filter(log10(length)>4)
short_introns <- df_introns_in_transcr %>% filter(log10(length)<4)
long_introns_id <- long_introns$parent
short_introns_id <- short_introns$parent

#fai boxplot weighted intron length vs num of introns in long and short transcripts 
df_long_introns_weight <- df_introns_in_transcr_weight %>% filter(parent %in% long_introns_id,)
df_short_introns_weight <- df_introns_in_transcr_weight %>% filter(parent %in% short_introns_id,)

boxplot(list(long_introns=df_long_introns_weight$length, short_introns=df_short_introns_weight$length), ylab="Weighted intron length",main="Boxplot of weighted intron length", notch=TRUE, col="bisque")

#make the histogram
#fit the model??
#cumulative length vs weighted length for long and short introns??
#intron lenght vs CDS length 









#' * Analysis of intron length in long and short introns

#' * Analysis of CDS






