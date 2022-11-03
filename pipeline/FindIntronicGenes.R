#The aim of this script is to analyze introni genes and introns in a gff3 file
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GenomicRanges")

# ANALYSIS OF NUMBER AND LENGTH OF INTRONS INSIDE GENES
## DATA PRE-PROCESSING

#Set the directory and import the data
setwd("/mnt/picea/home/edamieto/Git/Master-thesis-project/pipeline")
spruce_pine <- read.table("../data/gff3/spruce_pine-extended_master_complete.gff3", sep = "\t")

#gff = general feature format

#Add the header
header <- c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes")
colnames(spruce_pine) <- header

#installation of Genomic Ranges done :)
library("GenomicRanges")

#Subset the dataframe to retrieve the information that we want
genes <- subset(spruce_pine, type=="gene")
introns <- subset(spruce_pine, type=="intron")  
exons <- subset(spruce_pine, type=="exon")
cds <- subset(spruce_pine, type=="CDS")

#different types: gene, mRNA, exon, CDS, intron
#exon and CDS have scores
unique(spruce_pine["type"])

#check the dimension of each of them 
dim(spruce_pine[which(spruce_pine["type"]=="gene"),])[1]
#29226
dim(spruce_pine[which(spruce_pine["type"]=="exon"),])[1]
#143040
dim(spruce_pine[which(spruce_pine["type"]=="CDS"),])[1]
#112479
dim(spruce_pine[which(spruce_pine["type"]=="mRNA"),])[1]
#29226
dim(spruce_pine[which(spruce_pine["type"]=="intron"),])[1]
#113814

#TEST: The number of mRNA and genes was the same so I intersected them and they are the same
#gene is the parent of mRNA in gff3 files
mRNA <- subset(spruce_pine, type=="mRNA") 
gen_mRNA <- intersect(genes, mRNA) 
#all the mRNA and genes coincide

#Use genomic ranges package on genes and introns ;)
#Genomic Ranges works with S4 vectors
ggenes <- GRanges(genes)
gintrons <- GRanges(introns)

#Create the overlap
#This overlap identifies the introns within the genes
#from=query (introns), to=subject (genes)
overlap <- findOverlaps(gintrons, ggenes, type = "within")

head(overlap, 20)
#Count the number of introns in a gene
intr_in_genes <- as.data.frame(table(overlap@to))

intr_in_genes <- as.data.frame(t(apply(intr_in_genes,1, as.numeric)))
colnames(intr_in_genes) <- c("Gene", "Num_introns")

#Find the exon with the highest number of introns 
#max(intr_in_genes[2])

intr_in_genes[which.max(intr_in_genes$Num_introns),]
#     Gene  Num_introns
#2242 3388          74

ggenes[3388,]
#chr 2
#start 106246055, stop: 106732584 -> 486529 bp long
#JBrowse: transcript_126602

## DATA VISUALIZATION

#Try some plots
#Gene number vs #introns
plot(intr_in_genes, type= "l", lwd=2, main="Plot of gene number vs number of introns", xlab="Gene number", ylab="# introns")


par(mar=c(4,4,4,4))
boxplot(intr_in_genes[2], main="Boxplot of the number of introns", ylab = "Number of Introns",
        xlim = c(0.5, 1.4))
#unique(gcds@elementMetadata@listData$source)
hist(intr_in_genes[,2], main = "Histogram of the number of introns", 
     xlab = "Number of introns", xlim = c(0,80), ylim = c(0,14000))

#Num of introns vs genes with ggplot
library(ggplot2)
ggplot(intr_in_genes, aes(Gene, Num_introns)) +
  geom_point() +
  coord_cartesian(xlim=c(0,30000), ylim=c(0,80)) +
  labs(x="Genes", y="Number of introns", title="Scatterplot of the number of introns")

ggplot(intr_in_genes, aes(Gene, Num_introns)) +
  geom_point() +
  coord_cartesian(xlim=c(0,30000), ylim=c(20,80)) +
  labs(x="Genes", y="Number of introns", title="Scatterplot of the number of introns above 20")

ggplot(intr_in_genes, aes(Num_introns)) +
  geom_histogram(binwidth = 2, fill="light blue") +
  coord_cartesian(xlim=c(0,45), ylim=c(0,10000)) +
  labs(x="Number of introns", y="Counts", title="Histogram of the number of introns")

counts_intr <- as.data.frame(table(intr_in_genes$Num_introns))

counts_intr <- as.data.frame(t(apply(counts_intr,1, as.numeric)))
colnames(counts_intr) <- c("Num_introns", "Counts")

ggplot(counts_intr, aes(Num_introns, Counts)) +
  geom_line() +
  labs(title="Plot of the counts of the number of introns")

ggplot(counts_intr, aes(Counts,Num_introns)) +
  geom_point()+
  labs(title="Plot of the counts of the number of introns")


# ANALYSIS OF INTRON NUMBER IN GENES NOT PRESENT IN THE CHROMOSOMES 

#Select the seqid that are not in the chromosomes so they don't have chr in the seqid
seqid <- genes[1]


#There are genes that we know in which chromosome they are located but not exaclty where ex. PA_chr02_sUL006
#There are other genes that we don't know where they are located, 
# maybe they are small circular chromosome, e.g PA_cUP0116. In which we are interested in

#g1 <- grepl("^PA_chr[0-9]{2}$", seqid$seqid)
g2 <- grepl("PA_chr", seqid$seqid)

not_chr <- seqid[which(!g2),]
plot(table(not_chr), main="Counts of non chromosomic genes", xlab="Non chromosomic genes",
     ylab="Counts", type="h")
not_chr_counts <- as.data.frame(table(not_chr))


not_chr_counts[which.max(not_chr_counts$Freq),]
#     not_chr Freq
# 76 PA_sUP003  210

# ANALYSIS OF THE TRANSCRIPT WITH THE HIGHEST NUMBER OF INTRONS

#Interesting transcript: transcript_126602
#check how many exons
interesting_gene <- ggenes[3388,]
gexons <- GRanges(exons)
gcds <- GRanges(cds)
over_exons <- findOverlaps(gexons, interesting_gene, type="within")
#75 exons
over_inrons <- findOverlaps(gintrons, interesting_gene, type="within")
#74 introns 
over_CDS <- findOverlaps(gcds, interesting_gene, type="within")
#10 CDS

# ANALYSIS OF THE 10 GENES WITH THE HIGHEST NUMBER OF INTRONS

#Look at the first 10 genes with most introns
#sorted
intr_in_genes_sorted <- intr_in_genes[order(intr_in_genes$Num_introns, decreasing=TRUE),] 
hist(intr_in_genes_sorted$Num_introns)
plot(intr_in_genes_sorted)


ggplot(intr_in_genes_sorted, aes(Num_introns)) +
  geom_density()+
  labs(x="# introns", y="Density", title="Density plot of the number of introns within genes")

ggplot(intr_in_genes_sorted, aes(Num_introns)) +
  geom_boxplot() +
  labs(x="# introns", title="Boxplot of the number of introns within genes")

ggplot(intr_in_genes_sorted, aes(Num_introns)) +
  geom_freqpoly(binwidth=2) +
  labs(x="# introns", y="Counts", title="Boxplot of the number of introns within genes")

high_freqintr <- head(intr_in_genes_sorted, 10)
most_interesting_genes <- ggenes[high_freqintr$Gene,]

# DISTRIBUTION PLOT OF INTRON LENGTH AND  NUMBER OF INTRONS WITHIN GENES

#Distribution plot: intron length and freq of introns
#Intron length
#the length is the difference between the intron end and the intron start +1 (if not there is an intron with length 0)
intron_length <- introns$end-introns$start+1
boxplot(intron_length, main="Boxplot of intron length", ylab="Intron length")
max(intron_length)
# 602496

min(intron_length)
#1
introns[which.min(intron_length),]
#         seqid source   type   start     end score strand phase                     attributes
#228 PA_cUP0116      . intron 2488401 2488401     .      +     . Parent=transcript_138712.mrna1

hist(intron_length)

df <- data.frame(Gene=overlap@to, Introns=overlap@from,Intron_length=gintrons@ranges@width)
merged <- merge.data.frame(df, intr_in_genes, by="Gene", all.x = TRUE)

plot(merged$Num_introns, merged$Intron_length, main="Scatterplot of intron length and occurrence",
     xlab="# introns", ylab="Intron length")


ggplot(merged, aes(Num_introns, Intron_length))+
  geom_point()+
  labs(x="# introns", y="Intron length", title="Scatterplot of intron length and occurrence")

cor(merged$Num_introns, merged$Intron_length)
#-0.02274121

#Density plot
ggplot(merged, aes(Num_introns, Intron_length))+
  geom_density_2d()+
  labs(x="# introns", y="Intron length", title="Density plot of intron length and occurrence")

library(tidyverse)

#Group By gene
df_grouped <- merged %>% group_by(Gene)

#summarise by sum
df_new <- df_grouped %>% summarise(
  Intron_length=sum(Intron_length),
  Num_introns= unique(Num_introns)
)

#Num of introns vs intron length 
ggplot(df_new, aes(Num_introns, Intron_length))+
  geom_point()+
  labs(x="# introns", y="Intron length", title="Scatterplot of intron length and occurrence")


ggplot(df_new, aes(Num_introns, Intron_length))+
  geom_point()+
  geom_smooth(method="lm")+
  labs(x="# introns", y="Intron length", title="Scatterplot of intron length and occurrence")

ggplot(df_new, aes(Num_introns, Intron_length))+
  geom_point()+
  geom_smooth(method="gam")+
  labs(x="# introns", y="Intron length", title="Scatterplot of intron length and occurrence")

ggplot(df_new, aes(Num_introns, Intron_length))+
  geom_point()+
  geom_smooth(method="loess")+
  labs(x="# introns", y="Intron length", title="Scatterplot of intron length and occurrence")

#Density plot
ggplot(df_new, aes(Num_introns, Intron_length))+
  geom_density_2d()+
  labs(x="# introns", y="Intron length", title="Density plot of intron length and occurrence")

#try with the linear model 
linear_model <- lm(Intron_length~Num_introns, df_new)


ks.test(df_new$Num_introns, pnorm, mean(df_new$Num_introns), sd(df_new$Num_introns))
#p-value <0.05 (significant p-value). The null hypothesis is that the distribution is normal (comparison between data distribution and normal distribution)
#Since the p-value is smaller than the significance level, the null hypothesis should be rejected. 
#The data doesn't follow a normal distribution so it doesn't make sense run the anova!!

#To run ANOVA the data should be:
#1. Independence of the data 
#2. Normality of the distribution
#3. Homogeneity of variance 

#anova(linear_model)

summary(linear_model)
par(mfrow=c(2,2))
plot(linear_model)

#generalized additive models with integrated smoothness estimation
library(mgcv)
gen_add_mod <- gam(df_new$Intron_length~s(df_new$Num_introns, bs="cs"))
summary(gen_add_mod)
par(mfrow=c(1,1))
plot(gen_add_mod)


cor(df_new$Num_introns, df_new$Intron_length)
# 0.5964437

#Loess regression
loess_mod <- loess(Intron_length~Num_introns, df_new)
summary(loess_mod)
plot(loess_mod)

#PLOTTTTTTTTT :)
#Try with plot in log scale
ggplot(df_new, aes(Num_introns, log10(Intron_length)))+
  geom_point()+
  labs(x="# introns", y="log(Intron length)", title="Scatterplot of intron length and occurrence")


new_model <- lm(log10(Intron_length)~Num_introns, df_new)
summary(new_model)


ggplot(df_new, aes(Num_introns, log10(Intron_length)))+
  geom_point()+
  geom_smooth(method="gam")+
  labs(x="# introns", y="log(Intron length)", title="Scatterplot of intron length and occurrence")


ggplot(df_new, aes(x=Num_introns, y=log10(Intron_length)))+
  geom_point()+
  geom_smooth(method="lm", formula= y~log10(x))+
  labs(x="# introns", y="log(Intron length)", title="Scatterplot of intron length and occurrence")

new_model_2 <- lm(log10(Intron_length)~log10(Num_introns), df_new)
summary(new_model_2)

ggplot(df_new, aes(log10(Num_introns), log10(Intron_length)))+
  geom_point()+
  geom_smooth(method="lm", se=TRUE)+
  labs(x="log # introns", y="log(Intron length)", title="Scatterplot of intron length and occurrence")

#COMPARISON PLOT
library(LSD)

comparisonplot(df_new$Num_introns, df_new$Intron_length, xlab = "# introns", ylab = "Intron length",
               xlim=c(0,80), main="Comparison plot of intron length and occurrence")

comparisonplot(df_new$Num_introns, df_new$Intron_length, xlab = "# introns", ylab = "Intron length", cor=TRUE,
               xlim=c(0,80), main="Comparison plot of intron length and occurrence,")

comparisonplot(df_new$Num_introns, log(df_new$Intron_length), xlab = "# introns", ylab = "log Intron length", cor = TRUE,
               xlim=c(0,80), ylim=c(0,15), main="Comparison plot of intron length and occurrence,")

#zoom the x axis
comparisonplot(df_new$Num_introns, log(df_new$Intron_length), xlab = "# introns", ylab = "log Intron length", cor = TRUE,
               xlim=c(0,20), ylim=c(0,15), main="Comparison plot of intron length and occurrence,")

#there are two peaks, it could be interesting to analyze them if they are still present in the new file!


#Plot, plot and more plot :)
ggplot(df_new, aes(log10(Intron_length)))+
  geom_histogram()

ggplot(df_new, aes(log10(Intron_length)))+
  geom_density()+
  geom_vline(xintercept = 9.4, colour="red")

# ANALYSIS OF GENES WITHIN INTRONS

#Find genes inside introns
#In ovl object the query object is ggenes and the subject one is gintrons
ovl <- findOverlaps(ggenes, gintrons, type="within") #0 overlap 
ovl <- findOverlaps(ggenes, gintrons, type="within",ignore.strand=TRUE)
#ovl is not empty just if you ignore strand so this mean that every time the intron and the gene are in the opposite strand
# so opposite direction and this is interesting 

introns_with_genes <- gintrons[subjectHits(ovl),]
genes_within_introns<- ggenes[queryHits(ovl),]

#you could extract the ID for the genes and the parents for introns from the 9th column 

intersect(unique(genes[ovl@from,]$seqid),unique(genes_within_introns@seqnames@values))


#PLOT
#401 unique genes, 46 unique genes within introns, 10 unique not chromosomic genes within introns 

par(mar=c(5,5,2,5))
#barplot of chromosome 
barplot(sort(table(as.character(seqnames(ggenes[queryHits(ovl),])))),las=2, cex.names = 0.5,
        ylim = c(0,120), main="Barplot of the number of genes within introns", ylab = "# of genes within introns")

as.character(seqnames(ggenes[queryHits(ovl),]))


#CREA UN DF WITH GENE, INTRON TOTAL LENGTH, CHR
df_within <- data.frame(Gene=ovl@to, Intron=ovl@from, Intron_length=gintrons@ranges@width[ovl@from], chr= as.character(seqnames(ggenes[queryHits(ovl),])))

#Group By gene
df_grouped_within <- df_within %>% group_by(Gene)

#summarise by sum
df_new_within <- df_grouped_within %>% summarise(
  Intron_length=sum(Intron_length),
  chr=unique(chr)
)

#plot
ggplot(df_new_within, aes(chr, Intron_length))+
  geom_boxplot()+
  labs(title="Boxplot of cumulative intron length of the chromosomes", x="Chr", y="Cumulative intron length")+
  theme(axis.text.x = element_text(angle = 90)) 

ggplot(df_new_within, aes(chr, log10(Intron_length)))+
  geom_boxplot()+
  labs(title="Boxplot of cumulative intron length of the chromosomes", x="Chr", y="Cumulative intron length")+
  theme(axis.text.x = element_text(angle = 90)) 

## Consider just the genes within the chromosomes
df_chr_genes <- df_new_within %>% filter(str_detect(chr,  regex(paste0("^PA_chr[0-9]{2}$"))))

ggplot(df_chr_genes, aes(chr, Intron_length))+
  geom_boxplot()+
  labs(title="Boxplot of cumulative intron length of the chromosomes", x="Chr", y="Cumulative intron length")+
  theme(axis.text.x = element_text(angle = 90)) 

#Transform the y axis in log scale 
ggplot(df_chr_genes, aes(chr, log10(Intron_length)))+
  geom_boxplot()+
  labs(title="Boxplot of cumulative intron length of the chromosomes", x="Chr", y="Cumulative intron length")+
  theme(axis.text.x = element_text(angle = 90)) 

#Add standard deviation (notch)
ggplot(df_chr_genes, aes(chr, log10(Intron_length)))+
  geom_boxplot(notch=TRUE)+
  labs(title="Boxplot of cumulative intron length of the chromosomes", x="Chr", y="Cumulative intron length")+
  theme(axis.text.x = element_text(angle = 90)) 

#Notches are used in box plots to help visually assess whether the medians of distributions differ. 
#If the notches do not overlap, this is evidence that the medians are different.

#violin plot
ggplot(df_chr_genes, aes(chr, log10(Intron_length)))+
  geom_violin(fill= "light blue")+
  labs(title="Violin plot of cumulative intron length of the chromosomes", x="Chr", y="Cumulative intron length")+
  theme(axis.text.x = element_text(angle = 90)) 

## Consider the genes not in the chromosomes
df_not_chr_genes <- df_new_within %>% filter(str_detect(chr,  regex(paste0("PA_chr[0-9]")), negate = TRUE))
ggplot(df_not_chr_genes, aes(chr, Intron_length))+
  geom_boxplot()+
  labs(title="Boxplot of cumulative intron length of not chromosomic genes", x="Chr", y="Cumulative intron length")+
  theme(axis.text.x = element_text(angle = 90)) 

#Transform the y axis in log scale
ggplot(df_not_chr_genes, aes(chr, log10(Intron_length)))+
  geom_boxplot()+
  labs(title="Boxplot of cumulative intron length of not chromosomic genes", x="Chr", y="Cumulative intron length")+
  theme(axis.text.x = element_text(angle = 90)) 


#look at the complementary dataset -> genes not in introns, drop genes without introns
ovl_compl <- findOverlaps(ggenes, gintrons, type="any", ignore.strand= TRUE)
ovl_compl <- ovl_compl[!ovl_compl %in% ovl]

ovl_compl[ovl_compl@to==0,] #0 objects
ovl_compl[ovl_compl@from==0,] #0 objects

#create the complement dataframe
df_compl <- data.frame(Gene=ovl_compl@to, Intron=ovl_compl@from, Intron_length=gintrons@ranges@width[ovl_compl@from], chr= as.character(seqnames(ggenes[queryHits(ovl_compl),])))

#Group By gene
df_grouped_compl <- df_compl %>% group_by(Gene)

#summarise by sum
df_new_compl <- df_grouped_compl %>% summarise(
  Intron_length=sum(Intron_length),
  chr=unique(chr)
)

#barplot of chromosome 
par(mar=c(5,5,4,5))
barplot(sort(table(as.character(seqnames(ggenes[queryHits(ovl_compl),])))),las=2, cex.names = 0.5,
        ylim = c(0,10000), main="Barplot of the number of genes not within introns", ylab = "# of genes within introns")


#plot

## Consider just the genes within the chromosomes
df_chr_genes_compl <- df_new_compl %>% filter(str_detect(chr,  regex(paste0("^PA_chr[0-9]{2}$"))))

ggplot(df_chr_genes_compl, aes(chr, Intron_length))+
  geom_boxplot()+
  labs(title="Boxplot of cumulative intron length of the chromosomes", x="Chr", y="Cumulative intron length")+
  theme(axis.text.x = element_text(angle = 90)) 

#Transform the y axis in log scale 
ggplot(df_chr_genes_compl, aes(chr, log10(Intron_length)))+
  geom_boxplot()+
  labs(title="Boxplot of cumulative intron length of the chromosomes", x="Chr", y="Cumulative intron length")+
  theme(axis.text.x = element_text(angle = 90)) 

#Add standard deviation (notch)
ggplot(df_chr_genes_compl, aes(chr, log10(Intron_length)))+
  geom_boxplot(notch=TRUE)+
  labs(title="Boxplot of cumulative intron length of the chromosomes", x="Chr", y="Cumulative intron length")+
  theme(axis.text.x = element_text(angle = 90)) 

#violin plot
ggplot(df_chr_genes_compl, aes(chr, log10(Intron_length)))+
  geom_violin(fill= "light blue")+
  labs(title="Violin plot of cumulative intron length of the chromosomes", x="Chr", y="Cumulative intron length")+
  theme(axis.text.x = element_text(angle = 90)) 

















## GRanges object:
gr <- GRanges(
  seqnames=Rle(c("chr1", "chr2", "chr1", "chr3"), c(1, 3, 2, 4)),
  ranges=IRanges(1:10, width=10:1, names=head(letters,10)),
  strand=Rle(strand(c("-", "+", "*", "+", "-")), c(1, 2, 2, 3, 2)),
  score=1:10,
  GC=seq(1, 0, length=10)
)
gr

## GRangesList object:
gr1 <- GRanges(seqnames="chr2", ranges=IRanges(4:3, 6),
               strand="+", score=5:4, GC=0.45)
gr2 <- GRanges(seqnames=c("chr1", "chr1"),
               ranges=IRanges(c(7,13), width=3),
               strand=c("+", "-"), score=3:4, GC=c(0.3, 0.5))
gr3 <- GRanges(seqnames=c("chr1", "chr2"),
               ranges=IRanges(c(1, 4), c(3, 9)),
               strand=c("-", "-"), score=c(6L, 2L), GC=c(0.4, 0.1))
grl <- GRangesList("gr1"=gr1, "gr2"=gr2, "gr3"=gr3)

## Overlapping two GRanges objects:
table(!is.na(findOverlaps(gr, gr1, select="arbitrary")))
countOverlaps(gr, gr1)
findOverlaps(gr, gr1)
subsetByOverlaps(gr, gr1)

countOverlaps(gr, gr1, type="start")
findOverlaps(gr, gr1, type="start")
subsetByOverlaps(gr, gr1, type="start")

findOverlaps(gr, gr1, select="first")
findOverlaps(gr, gr1, select="last")

findOverlaps(gr1, gr)
findOverlaps(gr1, gr, type="start")
countOverlaps(gr1, gr, type="within")
findOverlaps(gr1, gr, type="equal")

## ---------------------------------------------------------------------
## MORE EXAMPLES
## ---------------------------------------------------------------------

table(!is.na(findOverlaps(gr, gr1, select="arbitrary")))
countOverlaps(gr, gr1)
findOverlaps(gr, gr1)
subsetByOverlaps(gr, gr1)

## Overlaps between a GRanges and a GRangesList object:

table(!is.na(findOverlaps(grl, gr, select="first")))
countOverlaps(grl, gr)
findOverlaps(grl, gr)
subsetByOverlaps(grl, gr)
countOverlaps(grl, gr, type="start")
findOverlaps(grl, gr, type="start")
subsetByOverlaps(grl, gr, type="start")
findOverlaps(grl, gr, select="first")

table(!is.na(findOverlaps(grl, gr1, select="first")))
countOverlaps(grl, gr1)
findOverlaps(grl, gr1)
subsetByOverlaps(grl, gr1)
countOverlaps(grl, gr1, type="start")
findOverlaps(grl, gr1, type="start")
subsetByOverlaps(grl, gr1, type="start")
findOverlaps(grl, gr1, select="first")

## Overlaps between two GRangesList objects:
countOverlaps(grl, rev(grl))
findOverlaps(grl, rev(grl))
subsetByOverlaps(grl, rev(grl))






