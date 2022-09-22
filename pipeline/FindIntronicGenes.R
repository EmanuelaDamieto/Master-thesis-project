#This script should be used to find intronic genes :)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GenomicRanges")


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
overlap <- findOverlaps(gintrons, ggenes, type = "within")

ovl <- findOverlaps(ggenes, gintrons, type="within")

head(overlap, 20)
intr_in_genes <- as.data.frame(table(overlap@to))

intr_in_genes <- as.data.frame(t(apply(intr_in_genes,1, as.numeric)))
colnames(intr_in_genes) <- c("Gene", "Num_introns")
#intr_in_genes <- sapply(intr_in_genes[1], as.integer)

#Find the exon with the highest number of introns 
#max(intr_in_genes[2])

intr_in_genes[which.max(intr_in_genes$Num_introns),]
#     Gene  Num_introns
#2242 3388          74

ggenes[3388,]
#chr 2
#start 106246055, stop: 106732584 -> 486529 bp long
#JBrowse: transcript_126602

#overlap between cds and introns
#findOverlaps(ggenes, gintrons, type = "within") #0 overlap
findOverlaps(gcds, gintrons, type = "within") #0 overlap


#Try some plots
#Gene number vs #introns
plot(intr_in_genes, type= "l", lwd=2, main="Plot of gene number vs number of introns", xlab="Gene number", ylab="# introns")


par(mar=c(4,4,4,4))
boxplot(intr_in_genes[2], main="Boxplot of the number of introns", ylab = "Number of Introns",
        xlim = c(0.5, 1.4))
#unique(gcds@elementMetadata@listData$source)
hist(intr_in_genes[,2], main = "Histogram of the number of introns", 
     xlab = "Number of introns", xlim = c(0,80), ylim = c(0,14000))


#PLOT
#NUm of introns vs genes
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



#Select the seqid that are not in the chromosomes so they don't have chr in the seqid
seqid <- genes[1]


#There are genes that we know in which chromosome they are located but not where ex. PA_chr02_sUL006
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

#PLOT AGAIN
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
ggplot(df_new, aes(Num_introns, log(Intron_length)))+
  geom_point()+
  labs(x="# introns", y="log(Intron length)", title="Scatterplot of intron length and occurrence")


new_model <- lm(log(Intron_length)~Num_introns, df_new)
summary(new_model)


ggplot(df_new, aes(Num_introns, log(Intron_length)))+
  geom_point()+
  geom_smooth(method="gam")+
  labs(x="# introns", y="log(Intron length)", title="Scatterplot of intron length and occurrence")


ggplot(df_new, aes(x=Num_introns, y=log(Intron_length)))+
  geom_point()+
  geom_smooth(method="lm", formula= y~log(x))+
  labs(x="# introns", y="log(Intron length)", title="Scatterplot of intron length and occurrence")

new_model_2 <- lm(log(Intron_length)~log(Num_introns), df_new)
summary(new_model_2)

ggplot(df_new, aes(log(Num_introns), log(Intron_length)))+
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
ggplot(df_new, aes(log(Intron_length)))+
  geom_histogram()+
  geom_title

ggplot(df_new, aes(log(Intron_length)))+
  geom_density()+
  geom_vline(xintercept = 9.4, colour="red")



#Find genes inside introns
overlap <- findOverlaps(gintrons, ggenes, type = "within")

ovl <- findOverlaps(ggenes, gintrons, type="within")

ggenes[subjectHits(overlap),]
gintrons[queryHits(overlap),]
#It should be the opposite if we manage to find the intronic genes :|

#you could extract the ID for the genes and the parents for introns from the 9th column 











