#' ---
#' title: "Drought root N. spruce Biological QA for expressed genes"
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
#' # Setup
#' This section and the next are relevant for reproducibility purposes. For results, please skip to section 3 (Quality Control)
#' 
#' * Libraries
suppressPackageStartupMessages({
  library(data.table)
  library(DESeq2)
  library(gplots)
  library(here)
  library(hyperSpec)
  library(magrittr)
  library(parallel)
  library(pheatmap)
  library(plotly)
  library(pvclust)
  library(tidyverse)
  library(tximport)
  library(vsn)
  library(emoji)
})

#' * Helper functions
source(here("UPSCb-common/src/R/featureSelection.R"))

#' * Graphics
hpal <- colorRampPalette(c("blue","white","red"))(100)

#' # Data
#' * Sample information
#' ```{r Instructions1,eval=FALSE,echo=FALSE}
#' # The csv file should contain the sample information, including the sequencing file name, 
#' # any relevant identifier, and the metadata of importance to the study design
#' # as columns, e.g. the SamplingTime for a time series experiment
#' ```

samples <- read_csv(here("data/drought_roots.csv"),
                      col_types=cols(.default=col_factor()))

#' * tx2gene translation table
#' ```{r Instructions2,eval=FALSE,echo=FALSE}
#' # This file is necessary if your species has more than one transcript per gene.
#' #
#' # It should then contain two columns, tab delimited, the first one with the transcript
#' # IDs and the second one the corresponding
#' #
#' # If your species has only one transcript per gene, e.g. Picea abies v1, then
#' # comment the next line
#' ```

tx2gene <- suppressMessages(read_delim(here("reference/annotation/Picab02_tx2gene.tsv.gz"), delim="\t", col_names=c("TXID","GENE"), skip=1))

#' * Raw data
filelist <- list.files(here("results/Salmon"), 
                          recursive = TRUE, 
                          pattern = "quant.sf",
                          full.names = TRUE)

#' * Sanity check to ensure that the data is sorted according to the sample info
#' ```{r Instructions3,eval=FALSE,echo=FALSE}
#' # This step is to validate that the salmon files are in the same order as 
#' # described in the samples object. If not, then they need to be sorted
#' ````
stopifnot(all(match(sub("[1-3]_151118_BC852HANXX_P2503_", "", sub("*_sortmerna_trimmomatic","",basename(dirname(filelist)))),
                    samples$SciLifeID) == 1:nrow(samples)))

samples_rep <-rbind(samples,samples)
samples_rep <- rbind(samples_rep,samples)

#' * Add filelist to samples as a new column
samples_rep %<>% mutate(Filenames = filelist)


#' * Export full rank samples
write_tsv(samples_rep,here("data/samples_full_rank.txt"))

#' * Read the expression at the gene level
#' ```{r Gene expression,eval=FALSE,echo=FALSE}
#' If the species has only one transcript per gene, or if you are conducting QA 
#' in transcript level replace with the following:
#' txi <- suppressMessages(tximport(files = samples$Filenames, type = "salmon",txOut=TRUE))
#' ```

txi <- suppressMessages(tximport(files = samples_rep$Filenames,
                                 type = "salmon",
                                 tx2gene=tx2gene))
counts <- txi$counts
samples_rep$Filenames <- sub("*_151118_BC852HANXX_P2503_", "_", sub("*_sortmerna_trimmomatic","",basename(dirname(samples_rep$Filenames))))
colnames(counts) <- samples_rep$Filenames
counts <- as.data.frame(counts)

#' * Import the list of the genes expressed in the control and C2d condition to avoid bias
expr_genes <- read.table(here("data/analysis/expressed_genes.txt"))
colnames(expr_genes) <- "Genes"
counts <- filter(counts, rownames(counts) %in% expr_genes$Genes)


#' 
#' <hr />
#' &nbsp;
#' 
#' # Quality Control
#' ## "Not expressed" genes
sel <- rowSums(counts) == 0
sprintf("%s%% percent (%s) of %s genes are not expressed",
        round(sum(sel) * 100/ nrow(counts),digits=1),
        sum(sel),
        nrow(counts))

#' ## Sequencing depth
#' * Let us take a look at the sequencing depth, colouring by Level 
#' ```{r Sequencing depth,eval=FALSE,echo=FALSE}
#' ```
dat <- tibble(x=colnames(counts),y=colSums(counts)) %>% 
  bind_cols(samples_rep)

ggplot(dat,aes(x,y,fill=Level)) + 
  geom_col() + 
  scale_y_continuous(name="reads") +
  facet_grid(~ factor(Level), scales = "free") +
  theme_bw() + 
  theme(axis.text.x=element_text(angle=90,size=4), axis.title.x=element_blank()) +
  labs(title ="Sample sequencing depth")

#' ```{r ggplot note,eval=FALSE,echo=FALSE}
#' # facet_grid divides the plot into subplots facet into columns, based on discrete variables 
#' ```

#' `r emoji("point_right")` **We observe some variation in the raw sequencing depth between conditions**
#' 

#' ## Per-gene mean expression
#' 
#' _i.e._ the mean raw count of every gene across samples is calculated
#' and displayed on a log10 scale.
#' 

ggplot(data.frame(value=log10(rowMeans(counts))),aes(x=value)) + 
  geom_density(na.rm = TRUE) +
  ggtitle("Gene mean raw counts distribution") +
  scale_x_continuous(name="mean raw counts (log10)") + 
  theme_bw()

#' ```{r mean raw counts note,eval=FALSE,echo=FALSE}
#' # Ideally gene mean raw counts around 100 -> log10 100 = 2
#' ```
#' 
#' `r emoji("point_right")` **The gene mean raw counts distribution is as expected since the highest peak is around 2**

#' ```{r Per sample expression,eval=FALSE,echo=FALSE}
#' # In the following, the second mutate also needs changing, I kept it 
#' # as an example to illustrate the first line. SampleID would be 
#' # a column in the samples object (the metadata) that uniquely identify
#' # the samples.
#' # If you have only a single metadata, then remove the second mutate call
#' # If you have more, add them as needed.
#' ```
#' 
#' ## Per-sample expression

dat <- as.data.frame(log10(counts)) %>% 
  utils::stack() %>% 
  mutate(Level=samples_rep$Level[match(ind,samples_rep$Filenames)]) 


ggplot(dat,aes(x=values,group=ind,col=Level)) + 
  geom_density(na.rm = TRUE) + 
  ggtitle("Sample raw counts distribution") +
  scale_x_continuous(name="per gene raw counts (log10)") + 
  theme_bw()

#' `r emoji("point_right")` **All samples have the same sequencing depth characteristics and there is no deviation when we look at one or the other condition**
#' 
#' * Export raw expression data
dir.create(here("data/analysis/salmon"),showWarnings=FALSE,recursive=TRUE)
write.csv(counts,file=here("data/analysis/salmon/raw-unormalised-gene-expression_data_expressed.csv"))
#' 
#' <hr />
#' &nbsp;
#' 
#' # Data normalisation 
#' ## Preparation
#' 
#' For visualization, the data is submitted to a variance stabilization
#' transformation using _DESeq2_. The dispersion is estimated independently
#' of the sample condition and replicate. 
#'  
#'  ```{r Data normalisation,eval=FALSE,echo=FALSE}
#'  # In the following, we provide the expected expression model, based on the study design.
#'  # It is technically irrelevant here, as we are only doing the quality assessment of the data, 
#'  # but it does not harm setting it correctly for the differential expression analyses that may follow.
#'  ```

dds <- DESeqDataSetFromTximport(
  txi=txi,
  colData = samples_rep,
  design = ~ Level)



dds <- dds[rownames(dds) %in% expr_genes$Genes,]

#' ## Size factors 
#' (_i.e._ the sequencing library size effect)
#' 

dds <- estimateSizeFactors(dds) %>%   
  suppressMessages()

#' boxplot of the sequencing libraries size factor:
boxplot(normalizationFactors(dds),
        main="Sequencing libraries size factor",
        las=2,log="y")
abline(h=1, col = "Red", lty = 3)


#' and without outliers:
boxplot(normalizationFactors(dds),
        main="Sequencing libraries size factor",
        las=2,log="y", outline=FALSE)
abline(h=1, col = "Red", lty = 3)

#' `r emoji("point_right")` **The sequencing libraries size factor is highly variable across samples but similar for the technical replicates**
#'
#' ```{r echo=FALSE,eval=FALSE}
#' # Assess whether there might be a difference in library size linked to a
#' # given metadata
#' # Developer: This would need to be ggplot2'ed

#' # ???????????? what is sizes???????????
#' # boxplot(split(sizes,dds$Level),las=2, main="Sequencing libraries size factor by Level")

#' # plot(sizes,log10(colSums(counts(dds))),ylab="log10 raw depth",xlab="scaling factor", col=rainbow(n=nlevels(dds$CHANGEME))[as.integer(dds$CHANGEME)],pch=19)
#' # legend("bottomright",fill=rainbow(n=nlevels(dds$CHANGEME)), legend=levels(dds$CHANGEME),cex=0.6)
#' ```


#' ## Validation
#' 
#' Let's look at standard deviations before and after VST normalization. 
#' This plot is to see whether there is a dependency of SD on the mean. 
#' 
#' Before:  
meanSdPlot(log1p(counts(dds))[rowSums(counts(dds))>0,])

#' After VST normalization, the red line should be almost horizontal which indicates
#' no dependency of variance on mean (homoscedastic).
#' 
#' ## Variance Stabilizing Transformation
vsd <- varianceStabilizingTransformation(dds, blind=TRUE)
vst <- assay(vsd)
vst <- vst - min(vst)

meanSdPlot(vst[rowSums(vst)>0,])

#' ```{r nice Sd plots, echo=FALSE,eval=FALSE}
#' meanSdPlot(log1p(counts(dds))[rowSums(counts(dds))>0,])$gg + ggtitle("Mean counts vs SD before VST normalization")
#' meanSdPlot(vst[rowSums(vst)>0,])$gg + ggtitle("Mean counts vs SD after VST normalization")
#' ```


#' `r emoji("point_right")` **We can conclude that the variance stabilization worked adequately even if the red line is not perfectly horizontal**
#' 
#' <hr />
#' &nbsp;
#' 
#' ## QC on the normalised data
#' ### PCA
pc <- prcomp(t(vst))
percent <- round(summary(pc)$importance[2,]*100)

#' ### Scree plot
#' 
#' We define the number of variables of the model (=1) `r nvar=1`


#' and the number of possible combinations (=8). 
#' ```{r PCA,eval=FALSE,echo=FALSE}
#' This needs to be adapted to your study design. Add or drop variables aas needed.
#' ```
#' `r nlevel=nlevels(dds$Level)`

#' We plot the percentage explained by different components
#' 
#' * the red line represents number of variables in the model  
#' * the orange line represents number of variable combinations.
#' 
nvar=1
nlevel=nlevels(dds$Level)
ggplot(tibble(x=1:length(percent),y=cumsum(percent)),aes(x=x,y=y)) +
  geom_line() + scale_y_continuous("variance explained (%)",limits=c(0,100)) +
  scale_x_continuous("Principal component") + 
  geom_vline(xintercept=nvar,colour="red",linetype="dashed",linewidth=0.5) + 
  geom_hline(yintercept=cumsum(percent)[nvar],colour="red",linetype="dashed",linewidth=0.5) +
  geom_vline(xintercept=nlevel,colour="orange",linetype="dashed",linewidth=0.5) + 
  geom_hline(yintercept=cumsum(percent)[nlevel],colour="orange",linetype="dashed",linewidth=0.5)+
  ggtitle("Percentage of variance explained by different components")
  
#' ### PCA plot
pc.dat <- bind_cols(PC1=pc$x[,1],
                    PC2=pc$x[,2],
                    as.data.frame(colData(dds)))

p <- ggplot(pc.dat,aes(x=PC1,y=PC2,col=Level,shape=Level,text=Filenames)) + 
  geom_point(size=2) + 
  ggtitle("Principal Component Analysis",subtitle="variance stabilized counts")

ggplotly(p) %>% 
  layout(xaxis=list(title=paste("PC1 (",percent[1],"%)",sep="")),
         yaxis=list(title=paste("PC2 (",percent[2],"%)",sep="")))

#' `r emoji("point_right")` **Some conditions clusterized better than others in the PCA plot**
#' 
#' ## Sequencing depth
#' Number of genes expressed per condition at different cutoffs:
conds <- factor(dds$Level)

dev.null <- rangeSamplesSummary(counts=vst,
                                conditions=conds,
                                nrep=3)
#' ## Heatmap
#' 
#' Filter for noise
#' 
sels <- rangeFeatureSelect(counts=vst,
                           conditions=conds,
                           nrep=3) %>% 
  suppressWarnings()
vst.cutoff <- 2
abline(h=10000, col="Red", lty=3)

#' ```{r heatmap,eval=FALSE,echo=FALSE}
#' * Heatmap of "all" genes
#' 
#' hm <- heatmap.2(t(scale(t(vst[sels[[vst.cutoff+1]],]))),
#'           distfun=pearson.dist,
#'           hclustfun=function(X){hclust(X,method="ward.D2")},
#'           labRow = NA,trace = "none",
#'           labCol = conds,
#'           col=hpal)
#' 
#' plot(as.hclust(hm$colDendrogram),xlab="",sub="", cex.axis=2)
#' ```

#' * Set the cut off to 5 in order to retrieve less than 10 000 genes
vst.cutoff <- 5

hm_reduced <- heatmap.2(t(scale(t(vst[sels[[vst.cutoff+1]],]))),
                distfun=pearson.dist,
                hclustfun=function(X){hclust(X,method="ward.D2")},
                labRow = NA,trace = "none",
                labCol = conds,
                col=hpal)

plot(as.hclust(hm_reduced$colDendrogram),xlab="",sub="", cex.axis=2)

#' * Using pheatmap 
mat <- t(scale(t(vst[sels[[vst.cutoff+1]],])))
hm <- pheatmap(mat,
               color = hpal,
               border_color = NA,
               clustering_distance_rows = "correlation",
               clustering_method = "ward.D2",
               show_rownames = F,
               labels_col = conds,
               angle_col = 90,
               legend = F)


#' `r emoji("point_right")` **The different conditions are not so different in gene expression level**
#' 

#' ## Clustering of samples
#' Done to assess the previous dendrogram's reproducibility
hm.pvclust <- pvclust(data = t(scale(t(vst[sels[[vst.cutoff+1]],]))),
                       method.hclust = "ward.D2", 
                       nboot = 1000, parallel = TRUE)

#' Plot the clustering with bp and au

plot(hm.pvclust, labels = conds, cex.axis=1.5)
pvrect(hm.pvclust)

#' `r emoji("point_right")` **Some conditions clusterize better than others**
#' 
#' <details><summary>bootstrapping results as a table</summary>
#' ```{r bootstrapping results as a table}
#' print(hm.pvclust, digits=3)
#' ```
#' </details>


#' # Combine technical replicates 
#' ```{tech rep, echo=FALSE, eval=FALSE}
#' # The block of code is meant to combine tech reps - as it is facoltative it is commented out
#' # First create a new variable in your sample object called BioID that identifies uniquely technical replicates, so one value for all tech rep of the same bio rep
#' samples_rep$Filenames <- filelist
#' samples_rep$BioID <- sub("[1-3]_151118_BC852HANXX_", "", sub("*_sortmerna_trimmomatic","",basename(dirname(samples_rep$Filenames))))
#' # Merging technical replicates
#' txi$counts <- sapply(split.data.frame(t(txi$counts),samples_rep$BioID),colSums)
#' txi$length <- sapply(split.data.frame(t(txi$length),samples_rep$BioID),colMaxs)
#' # Counts are now in alphabetic order, check and reorder if necessary
#' stopifnot(colnames(txi$counts) == samples_rep$BioID)
#' samples_rep <- samples_rep[match(colnames(txi$counts),samples_rep$BioID),]
#' # Recreate the dds
#' dds <- DESeqDataSetFromTximport(
#'        txi=txi,
#'        colData = samples_rep,
#'        design = ~ Level)
#'
#'dds <- dds[rownames(dds) %in% expr_genes$Genes,]
#' save(dds,file=here("data/analysis/salmon/dds_merge_lengthScaledTPM_expr_genes.rda"))
#' ```
#' 

samples_rep$Filenames <- filelist
samples_rep$BioID <- sub("[1-3]_151118_BC852HANXX_", "", sub("*_sortmerna_trimmomatic","",basename(dirname(samples_rep$Filenames))))
#' Merging technical replicates
txi$counts <- sapply(split.data.frame(t(txi$counts),samples_rep$BioID),colSums)
txi$length <- sapply(split.data.frame(t(txi$length),samples_rep$BioID),colMaxs)
#' Counts are now in alphabetic order, check and reorder if necessary
stopifnot(colnames(txi$counts) == samples_rep$BioID)
samples_rep <- samples_rep[match(colnames(txi$counts),samples_rep$BioID),]
#' Recreate the dds
dds <- DESeqDataSetFromTximport(
        txi=txi,
       colData = samples_rep,
       design = ~ Level)
dds <- dds[rownames(dds) %in% expr_genes$Genes,]

save(dds,file=here("data/analysis/salmon/dds_merge_expr_genes.rda"))

#' # Create the vst with counts normalized by the transcript length
#' ```{r Normalized counts,eval=FALSE,echo=FALSE}
#' It will needed for some of the downstream analysis. 
#' Here, we are generating counts from abundances, using the argument countsFromAbundance, 
#' scaled to library size and the average transcript length, averaged over samples. 
#' In this way we cannot compared different genes but just different transcripts of the same gene. 
#'
#' ```
txi <- suppressMessages(tximport(files = samples_rep$Filenames,
                                 type = "salmon",
                                 tx2gene=tx2gene,
                                 countsFromAbundance="lengthScaledTPM"))
counts <- txi$counts
samples_rep$Filenames <- sub("*_151118_BC852HANXX_P2503_", "_", sub("*_sortmerna_trimmomatic","",basename(dirname(samples_rep$Filenames))))
colnames(counts) <- samples_rep$Filenames
counts <- as.data.frame(counts)
counts <- filter(counts, rownames(counts) %in% expr_genes$Genes)


dds <- DESeqDataSetFromTximport(
  txi=txi,
  colData = samples_rep,
  design = ~ Level)

colnames(dds) <- samples_rep$Filenames

dds <- dds[rownames(dds) %in% expr_genes$Genes,]

vsd <- varianceStabilizingTransformation(dds,blind=FALSE)
vst <- assay(vsd)
vst <- vst - min(vst)

save(vst,file=here("data/analysis/DE/vst-aware-lengthScaledTPM_exprGenes.rda"))
write_delim(as.data.frame(vst) %>% rownames_to_column("ID"),
            here("data/analysis/DE/vst-aware-lengthScaledTPM_exprGenes.tsv"))

#' <hr />
#' &nbsp;
#' 
#' # Summary
#' `r emoji("star")` **The data are quite good so we can continue with our analysis**
#' 
#' 


#' # Session Info
#' <details><summary>Session Info</summary>
#' ```{r session info}
#' ```
sessionInfo()
#'
