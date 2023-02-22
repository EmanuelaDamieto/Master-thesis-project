#' ---
#' title: "Differential Expression analysis of drought roots of Norway spruce"
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
    library(data.table)
    library(DESeq2)
    library(gplots)
    library(here)
    library(hyperSpec)
    library(RColorBrewer)
    library(tidyverse)
    library(VennDiagram)
})

#' * Helper files
suppressMessages({
    source(here("UPSCb-common/Rtoolbox/src/plotEnrichedTreemap.R"))
    source(here("UPSCb-common/src/R/featureSelection.R"))
    source(here("UPSCb-common/src/R/volcanoPlot.R"))
    source(here("UPSCb-common/src/R/gopher.R"))
})

#' * Graphics
pal=brewer.pal(8,"Dark2")
hpal <- colorRampPalette(c("blue","white","red"))(100)
mar <- par("mar")

#' ```{r TODO, echo=FALSE,eval=FALSE}
#' TODO remember to add a function to check for correlation between logFC and 
#' transcript length - check for the effective length
#' ```
#' 

#' * Functions
#' 1. plot specific gene expression
"line_plot" <- function(dds=dds,vst=vst,gene_id=gene_id){
    message(paste("Plotting",gene_id))
    sel <- grepl(gene_id,rownames(vst))
    stopifnot(sum(sel)==1)

    p <- ggplot(bind_cols(as.data.frame(colData(dds)),
                          data.frame(value=vst[sel,])),
                aes(x=Level,y=value,col=Level,group=Level)) +
        geom_point() + geom_smooth() +
        scale_y_continuous(name="VST expression") + 
        ggtitle(label=paste("Expression for: ",gene_id))
    
    suppressMessages(suppressWarnings(plot(p)))
    return(NULL)
}

#' 2. extract the DE results. Default cutoffs are
#' from Schurch _et al._, RNA, 2016
"extract_results" <- function(dds,vst,contrast,
                              padj=0.01,lfc=0.5,
                              plot=TRUE,verbose=TRUE,
                              export=TRUE,default_dir=here("data/analysis/DE"),
                              default_prefix="DE-",
                              labels=colnames(dds),
                              sample_sel=1:ncol(dds),
                              expression_cutoff=0,
                              debug=FALSE,filter=c("median",NULL),...){
  
  # get the filter
  if(!is.null(match.arg(filter))){
    filter <- rowMedians(counts(dds,normalized=TRUE))
    message("Using the median normalized counts as default, set filter=NULL to revert to using the mean")
  }
  
  # validation
  if(length(contrast)==1){
    res <- results(dds,name=contrast,filter = filter)
  } else {
    res <- results(dds,contrast=contrast,filter = filter)
  }
  
  stopifnot(length(sample_sel)==ncol(vst))
  
  if(plot){
    par(mar=c(5,5,5,5))
    volcanoPlot(res)
    par(mar=mar)
  }
  
  # a look at independent filtering
  if(plot){
    plot(metadata(res)$filterNumRej,
         type="b", ylab="number of rejections",
         xlab="quantiles of filter")
    lines(metadata(res)$lo.fit, col="red")
    abline(v=metadata(res)$filterTheta)
  }
  
  if(verbose){
    message(sprintf("The independent filtering cutoff is %s, removing %s of the data",
                    round(metadata(res)$filterThreshold,digits=5),
                    names(metadata(res)$filterThreshold)))
    
    max.theta <- metadata(res)$filterNumRej[which.max(metadata(res)$filterNumRej$numRej),"theta"]
    message(sprintf("The independent filtering maximises for %s %% of the data, corresponding to a base mean expression of %s (library-size normalised read)",
                    round(max.theta*100,digits=5),
                    round(quantile(counts(dds,normalized=TRUE),probs=max.theta),digits=5)))
  }
  
  if(plot){
    qtl.exp=quantile(counts(dds,normalized=TRUE),probs=metadata(res)$filterNumRej$theta)
    dat <- data.frame(thetas=metadata(res)$filterNumRej$theta,
                      qtl.exp=qtl.exp,
                      number.degs=sapply(lapply(qtl.exp,function(qe){
                        res$padj <= padj & abs(res$log2FoldChange) >= lfc & 
                          ! is.na(res$padj) & res$baseMean >= qe
                      }),sum))
    if(debug){
      plot(ggplot(dat,aes(x=thetas,y=qtl.exp)) + 
             geom_line() + geom_point() +
             scale_x_continuous("quantiles of expression") + 
             scale_y_continuous("base mean expression") +
             geom_hline(yintercept=expression_cutoff,
                        linetype="dotted",col="red"))
      
      p <- ggplot(dat,aes(x=thetas,y=qtl.exp)) + 
        geom_line() + geom_point() +
        scale_x_continuous("quantiles of expression") + 
        scale_y_log10("base mean expression") + 
        geom_hline(yintercept=expression_cutoff,
                   linetype="dotted",col="red")
      suppressMessages(suppressWarnings(plot(p)))
      
      plot(ggplot(dat,aes(x=thetas,y=number.degs)) + 
             geom_line() + geom_point() +
             geom_hline(yintercept=dat$number.degs[1],linetype="dashed") +
             scale_x_continuous("quantiles of expression") + 
             scale_y_continuous("Number of DE genes"))
      
      plot(ggplot(dat,aes(x=thetas,y=number.degs[1] - number.degs),aes()) + 
             geom_line() + geom_point() +
             scale_x_continuous("quantiles of expression") + 
             scale_y_continuous("Cumulative number of DE genes"))
      
      plot(ggplot(data.frame(x=dat$thetas[-1],
                             y=diff(dat$number.degs[1] - dat$number.degs)),aes(x,y)) + 
             geom_line() + geom_point() +
             scale_x_continuous("quantiles of expression") + 
             scale_y_continuous("Number of DE genes per interval"))
      
      plot(ggplot(data.frame(x=dat$qtl.exp[-1],
                             y=diff(dat$number.degs[1] - dat$number.degs)),aes(x,y)) + 
             geom_line() + geom_point() +
             scale_x_continuous("base mean of expression") + 
             scale_y_continuous("Number of DE genes per interval"))
      
      p <- ggplot(data.frame(x=dat$qtl.exp[-1],
                             y=diff(dat$number.degs[1] - dat$number.degs)),aes(x,y)) + 
        geom_line() + geom_point() +
        scale_x_log10("base mean of expression") + 
        scale_y_continuous("Number of DE genes per interval") + 
        geom_vline(xintercept=expression_cutoff,
                   linetype="dotted",col="red")
      suppressMessages(suppressWarnings(plot(p)))
    }
  }
  
  sel <- res$padj <= padj & abs(res$log2FoldChange) >= lfc & ! is.na(res$padj) & 
    res$baseMean >= expression_cutoff
  
  if(verbose){
    message(sprintf(paste(
      ifelse(sum(sel)==1,
             "There is %s gene that is DE",
             "There are %s genes that are DE"),
      "with the following parameters: FDR <= %s, |log2FC| >= %s, base mean expression > %s"),
      sum(sel),padj,
      lfc,expression_cutoff))
  }
  
  # proceed only if there are DE genes
  if(sum(sel) > 0){
    val <- rowSums(vst[sel,sample_sel,drop=FALSE])==0
    if (sum(val) >0){
      warning(sprintf(paste(
        ifelse(sum(val)==1,
               "There is %s DE gene that has",
               "There are %s DE genes that have"),
        "no vst expression in the selected samples"),sum(val)))
      sel[sel][val] <- FALSE
    } 
    
    if(export){
      if(!dir.exists(default_dir)){
        dir.create(default_dir,showWarnings=FALSE,recursive=TRUE,mode="0771")
      }
      write.csv(res,file=file.path(default_dir,paste0(default_prefix,"-results.csv")))
      write.csv(res[sel,],file.path(default_dir,paste0(default_prefix,"-genes.csv")))
    }
    if(plot & sum(sel)>1){
      heatmap.2(t(scale(t(vst[sel,sample_sel]))),
                distfun = pearson.dist,
                hclustfun = function(X){hclust(X,method="ward.D2")},
                trace="none",col=hpal,labRow = FALSE,
                labCol=labels[sample_sel],...
      )
    }
  }
  return(list(all=rownames(res[sel,]),
              up=rownames(res[sel & res$log2FoldChange > 0,]),
              dn=rownames(res[sel & res$log2FoldChange < 0,])))
}

#' 3. extract and plot the enrichment results
extractEnrichmentResults <- function(enrichment,task="go",
                                     diff.exp=c("all","up","dn"),
                                     go.namespace=c("BP","CC","MF"),
                                     genes=NULL,export=TRUE,plot=TRUE,
                                     default_dir=here("data/analysis/DE"),
                                     default_prefix="DE",
                                     url="athaliana"){
    # process args
    diff.exp <- match.arg(diff.exp)
    de <- ifelse(diff.exp=="all","none",
                 ifelse(diff.exp=="dn","down",diff.exp))

    # sanity
    if( is.null(enrichment[[task]]) | length(enrichment[[task]]) == 0){
        message(paste("No enrichment for",task))
    } else {

        # write out
        if(export){
            write_tsv(enrichment[[task]],
                      file=here(default_dir,
                                paste0(default_prefix,"-genes_GO-enrichment.tsv")))
            if(!is.null(genes)){
                write_tsv(
                    enrichedTermToGenes(genes=genes,terms=enrichment[[task]]$id,url=url,mc.cores=16L),
                    file=here(default_dir,
                              paste0(default_prefix,"-enriched-term-to-genes.tsv"))
                )
            }
        }
        
        if(plot){
            sapply(go.namespace,function(ns){
                titles <- c(BP="Biological Process",
                            CC="Cellular Component",
                            MF="Molecular Function")
                suppressWarnings(tryCatch({plotEnrichedTreemap(enrichment,enrichment=task,
                                                               namespace=ns,
                                                               de=de,title=paste(default_prefix,titles[ns]))},
                                          error = function(e) {
                                              message(paste("Treemap plot failed for",ns, 
                                                            "because of:",e))
                                              return(NULL)
                                          }))
            })
        }
    }
}

#' * Data
#' ```{r load, echo=FALSE,eval=FALSE}
#' Here you are meant to load an RData object
#' that contains a DESeqDataSet object. If you ran the 
#' biological QA template, you need not change anything.
#' If you have technical replicates merge them before to perform the DE analysis
#' ```

load(here("data/analysis/salmon/dds_merge.rda"))

#' ## Normalisation for visualisation
#' ```{r normalisation, echo=FALSE, eval=FALSE, warning=FALSE, message=FALSE}
#' Here you have to normalize the dds object before to perform the DE analysis
#' ```        
vsd <- varianceStabilizingTransformation(dds,blind=FALSE)
vst <- assay(vsd)
vst <- vst - min(vst)
dir.create(here("data/analysis/DE"),showWarnings=FALSE)
save(vst,file=here("data/analysis/DE/vst-aware.rda"))
write_delim(as.data.frame(vst) %>% rownames_to_column("ID"),
            here("data/analysis/DE/vst-aware.tsv"))

#' # Differential Expression
#' ```{r import dds, echo=FALSE,eval=FALSE, warning=FALSE, message=FALSE}
#' Perform the DE analysis using the DESeq package
#' ```
dds <- DESeq(dds)


#' * Dispersion estimation
plotDispEsts(dds)

#' * Check the different contrasts
resultsNames(dds)

#' # Results
#' ```{r res, echo=FALSE,eval=FALSE}
#' Here you need to define the contrast you want to 
#' study 
#' 
#' The `contrast` can be given
#' by name, as a list (numerator/denominator) or as a vector of weight (e.g. c(0,1));
#' read the DESeq2 vignette for more info
#' 
#' The `label` argument is typically one (or a combination) of the metadata stored in colData
#' 
#' The function allows for looking at the independent filtering results using `debug=TRUE`
#' 
#' If you are not satisfied with the default from DESeq2, you can set your own cutoff using `expression_cutoff`
#' 
#' You can change the default output file prefix using `default_prefix`
#' 
#' You can select the set of samples to be added to the `heatmap`, using the `sample_sel` argument. It takes a logical vector.
#' 
#' ```

#' ```{r contrast, echo=FALSE,eval=FALSE}
#' Evaluate the different contrasts. The 80% water content in the soil is our control.
#' ```

#' ```{r different contrasts, echo=TRUE,eval=FALSE}
#' ### Contrast 60 vs 80
#' contrast_60_vs_80 <- extract_results(dds=dds,vst=vst,contrast="Level_60._vs_80.", labels = dds$Level, default_prefix="DE-60vs80")
#' 
#' ### Contrast 40 vs 80
#' contrast_40_vs_80 <- extract_results(dds=dds,vst=vst,contrast="Level_40._vs_80.", labels = dds$Level, default_prefix="DE-40vs80")

#' ### Contrast 30 vs 80
#' contrast_30_vs_80 <- extract_results(dds=dds,vst=vst,contrast="Level_30._vs_80.", labels = dds$Level, default_prefix="DE-30vs80")

#' ### Contrast 307d vs 80
#' contrast_307d_vs_80 <- extract_results(dds=dds,vst=vst,contrast="Level_30.7d_vs_80.", labels = dds$Level, default_prefix="DE-307dvs80")

#' ### Contrast collapse vs 80
#' contrast_Collapse_vs_80 <- extract_results(dds=dds,vst=vst,contrast="Level_Collapse_vs_80.", labels = dds$Level, default_prefix="DE-collapsevs80")

#' ### Contrast rehydrate vs 80
#' contrast_Rehydrate_vs_80 <- extract_results(dds=dds,vst=vst,contrast="Level_Rehydrate_vs_80.", labels = dds$Level, default_prefix="DE-Rehydratevs80")

#' ```

#' ## Contrast C2d vs control
contrast_C2d_vs_80 <- extract_results(dds=dds,vst=vst,contrast="Level_C2d_vs_80.", labels = dds$Level, default_prefix="DE-C2dvs80")

#' ### Show the heatmap just for the contrast we are interested in (C2d vs control)
contrast_C2d_vs_80 <- extract_results(dds=dds,vst=vst,contrast="Level_C2d_vs_80.", sample_sel = dds$Level %in% c("80%","C2d"), labels = dds$Level, default_prefix="DE-C2dvs80")

#' ### Change the value of the log2fc to remove bias due to the different amount of gene expressed in the two conditions 
contrast_C2d_vs_80 <- extract_results(dds=dds,vst=vst,contrast="Level_C2d_vs_80.", sample_sel = dds$Level %in% c("80%","C2d"), labels = dds$Level, default_prefix="DE-C2dvs80-lfc2", lfc=2)


#' ## Gene Ontology enrichment
#' ```{r go, echo=FALSE,eval=FALSE}
#' Once you have obtained a list of candidate genes, you most probably want
#' to annotate them.
#' 
#' We first identify the background; _i.e._ the
#' population of expressed genes. We select the genes expressed in a least
#' 2 replicate of one condition at a cutoff of `exp`.
#' 
#' Next we run the enrichment using 
#' the gofer3 REST API (interfaced through the gopher.R script loaded at the
#' beginning of this file).
#' 
#' Finally we export the go enrichment as a complete table.
#' We used to export another table consisting
#' of only the `id` and `padj` columns for using as input for _e.g._
#' REVIGO; but since flash is EOL and REVIGO not updated, we instead rely on 
#' the RtoolBox treemap.
#' 
#' In addition we now also export the list of genes that most likely resulted in
#' the corresponding go enrichment.
#' 
#' Make sure to change the `url` to match your species
#' 
#' 
#' res.list <- contrast_C2d_vs_80$all
#' background <- rownames(vst)[featureSelect(vst,dds$Level,exp=0.4)]
#' bla_enr <- gopher(contrast_C2d_vs_80$all, alpha = 0.05, task="go", background = background, url="picab02", endpoint = "enrichment")
#' 

# enr.list <- lapply(res.list,function(r){
#      lapply(r,gopher,background=background,task="go",url="picab02")
# })


#' dev.null <- lapply(names(enr.list),function(n){
#'     lapply(names(enr.list[[n]]),function(de){
#'         extractEnrichmentResults(enr.list[[n]][[de]],
#'                                 diff.exp=de,
#'                                 genes=res.list[[n]][[de]],
#'                                 default_prefix=paste(n,de,sep="-"),
#'                                 url="athaliana")
#'     })
#' })
#' ```
#' 

#' #### GO of genes expressed just under stress condition (C2d)
#' Use all the genes as background 
background <- rownames(vst)[featureSelect(vst,dds$Level,exp=0.0)]

res.list <- read.table(here("data/analysis/just_expr_C2d.txt"))
bla_enr <- gopher(res.list$V1, alpha = 0.05, task="go", background = background, url="picab02", endpoint = "enrichment")
dev.null <- extractEnrichmentResults(bla_enr,
                                    url="piceab02")


#' #### GO of genes expressed just in the control 
#' Use all the genes as background 
res.list_cnt <- read.table(here("data/analysis/just_expr_cnt.txt"))
bla_enr_cnt <- gopher(res.list_cnt$V1, alpha = 0.05, task="go", background = background, url="picab02", endpoint = "enrichment")
dev.null <- extractEnrichmentResults(bla_enr_cnt,
                                     url="piceab02")

#' # Session Info 
#' <details><summary>Session Info</summary>
#'  ```{r session info, echo=FALSE, message=FALSE}
#'  sessionInfo()
#'  ```



