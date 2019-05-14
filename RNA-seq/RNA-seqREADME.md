# 3' RNA-seq (Quantseq)

Method for quantifying RNA expression levels.

## Getting Started

The following walk-through provides an example of the analysis pipeline for RNA-seq data. For the purposes of this example, we will be starting with demultiplexed fastq files. These files have been sequenced single ended 75bp at approximately 20 million paired reads per sample. 

The PeterMac bioinformatics core automatically applies Seqliner2.0 pipeline to the Fastq files to produce .bam and .bai files. if you wish to start with these skip to the **Counting Reads** section.

The packages used have been referenced at the end, the majority of this analysis is from the user manuals linked below, with only slight deviations from default settings. 

### Library preparation
 
RNA was extracted from cells using trizol and directzol RNA extraction kit according to manufacturers instructions. 
Quantseq library preparation kit was used 

### General overview of RNA-seq analysis

![GeneralOverview](https://github.com/madisonJK/ReferenceAnalysis/raw/master/RNA-seq/RNA-seq_Analysis_Overview.png)


### Software and Version Control

Hisat2 /2.1.0

samtools/1.4.1

cutadapt /1.14

FastQC /0.11.5

RNAseqQC /1.1.8

Subread /1.5.0-p3

R: 

Limma/

Deseq/ 

ggplot2/3.3.1

reshape2/1.4.3

pheatmap/1.0.12


## Aligning Fastq Files
 Fastq files are alined to the Mouse genome ( mm10) using Hisat. The resulting sam files are sorted, converted to bam and indexed with samtools. You will need a hisat indexed reference genome. 

```
module load hisat2
module load samtools

hisat2 -p 8 -x /data/reference/indexes/mouse/mm10/hisat2/grcm38/Mus_musculus.GRCm38.dna.toplevel -U $1 -S ${1}_hisat2.sam
samtools view -@ 4 -Sbo ${1}_hisat2.sam.bam ${1}_hisat2.sam
samtools sort -@ 4 -o ${1}_hisat2.sam.sorted.bam ${1}_hisat2.sam.bam
samtools index ${1}_hisat2.sam.sorted.bam ${1}_hisat2.sam.sorted.bam.bai
```
*Optional:* QC can be performed on the resulting Bam files using RNAseqc. You will need to input a sample file, which is a tab-delimited text file with 3 columns specifying ID, the filename of bam file, and comments. You will need to download the rRNA reference fileand a gtf reference file in addition to a fasta file of the reference genome. 

```
module load rna-seqc
module load java

RNA-SeQC -s "ID|filename.bam|comments"  -t ref/mm10genes.gtf -o m6a/Fastq/RNAseqc -r data/reference/indexes/mouse/mm10/fasta/Mus_musculus.GRCm38.dna.toplevel.fa -BWArRNA ref/mm10rRNA.fa ${1}

```

IGV tools can also be used to generate a .tdf for visualisation with IGV

```
module load igvtools

igvtools count -z 5 -w 25 -e 150 ${1}_hisat2.sam.sorted.bam ${1}_hisat2.sam.sorted.bam.tdf mm10

```

### Counting Reads

The resulting bam files can then be used as input into subreads *featureCounts* function:
*Unstranded is generally used, but can be altered by changing -s 1*

```
module load subread

featureCounts -t exon -g gene_id -F 'GTF/GFF' -O -M -s 0 -a data/reference/gtf/Mus_musculus.GRCm38.90.gtf  -o counts.txt *.sorted.bam
```

## Differential Methylation Analysis

### Limma/Voom

The remaining analysis and figure generation can be performed in R. This analysis is based on the [Combine RNA-seq tutorial](http://combine-australia.github.io/RNAseq-R/).



First, install and load libraries needed for analysis and figure generation
may need to use  ```BiocManager::install()``` instead of ```biocLite()``` if using Bioconducter 3.9

```
source("http://bioconductor.org/biocLite.R")
biocLite("edgeR")
biocLite("Glimma")
biocLite("edgeR")
biocLite("limma")
biocLite("gplots")
biocLite("org.Mm.eg.db")
biocLite("RColorBrewer")
biocLite("GenomicFeatures")
biocLite("topGO")
install.packages("dplyr")
biocLite("pheatmap")

library("GenomicFeatures")
library('limma')
library('edgeR')
library('Glimma')
library('gplots')
library('org.Mm.eg.db')
library('RColorBrewer')
library('biomaRt')
library(dplyr)
library(pheatmap)
library(ggfortify)
library(biomaRt)
library(reshape2)
library(ggrepel)

```

Read in the counts table generated above: 

```
seqdata <- read.delim("counts.txt", stringsAsFactors = FALSE)
#remove columns with extra data
countdata <- seqdata[,-(1:6)]
#rename rows with gene names 
rownames(countdata) <- seqdata[,1]
#rename columns with easier Identifiers - need to edit based on your sample names 
colnames(countdata)<-c("A1","A2","A3","B1","B2","B3")
```
If your counts are annotated with refseq_mrna IDs use the below code to covert to gene names. biomaRt can also convert ensembleIDs use listAttributes(mart) to see other options. 

```

mart <- useDataset("mmusculus_gene_ensembl", useMart("ensembl"))
genes <- rownames(countdata)
G_list <- getBM(filters= "refseq_mrna", attributes= c("refseq_mrna","external_gene_name"),values=genes,mart= mart)
seqdatamerge<-merge(countdata,G_list,by.x="row.names",by.y="refseq_mrna")
aggdata <- aggregate(seqdatamerge,by=list(seqdatamerge$external_gene_name),FUN = "mean")

```
Get rid of extra data columns (need to edit depending on number of samples), and write a table of final raw counts for easy reference later: 
```
allcounts<-aggdata[,c(-1,-2,-3,-9)]
rownames(allcounts)<-aggdata[,1]

write.table(allcounts,file= "GenenameCounts.txt",sep = "\t",quote = F)
```
Set a threshold and filter out genes with low/no counts. ```myCPM> 3``` refers to the number of samples required to surpass the count threshold (myCPM >( number of replicates - 1), if number of replicates is > 2, otherwise myCPM> 2 should be fine)

```
myCPM <- cpm(allcounts)
thresh <- myCPM > 3
keep <- rowSums(thresh) >= 2
counts.keep <- allcounts[keep,]

```
Make sample info based on your counts table and then create a contrast matrix. ncol = number of independent variables+1.
```

SI<-matrix(data = NA,nrow=ncol(counts.keep),ncol=2)
SI[,1]<-colnames(counts.keep)
SI[,2]<-c("A","A","A","B","B","B")
colnames(SI)<-c("Samplename","Variable")
SI<-as.data.frame(SI)

group <- paste(SI$Variable,sep=".")
group <- factor(group)
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)


cont.matrix <- makeContrasts(treatmentvscontrol =B-A,
                             levels=design)
```

Create DGE list, Voom normalise the data, and fit to a linear model
```
y <- DGEList(counts.keep)
y <- calcNormFactors(y)
logcounts <- cpm(y,log=TRUE)

plotMDS(y)


v <- voom(y,design,plot = TRUE, normalize.method = "quantile")

fit <- lmFit(v)

fit.cont <- contrasts.fit(fit, cont.matrix)
fit.cont <- eBayes(fit.cont)
summa.fit <- decideTests(fit.cont,coef=1)

```
Create a table of results sorted by pvalue and write to file

```
sigsorted<- topTable(fit.cont,coef=1,sort.by="p",n="Inf")
write.table(sigsorted, file= "results.txt", sep = "\t",quote = F,col.names = T) 
```
Also save normalised counts matrix for future reference

```
write.table(v$E,file = "normalisedCounts.txt",sep = "\t",quote = F,col.names = T) 
```

## Generating plots

MD plot
```
plotMD(fit.cont,coef=1,status=summa.fit)

```

### Heatmap


Make a heatmap of ony significant genes (p <0.05 and absolute log FC > 0.5) : 
```

filtered <- results[(results$P.Value <0.05)& ((results$logFC > 0.5) | (results$logFC < (-0.5))),]
filteredgenes<-rownames(filtered)
filteredheatmap<-v$E[filteredgenes,] 

pheatmap(filteredheatmap,scale="row",
         trace="none",main="Differentially expressed genes",
         border_color = NA,fontsize_row = 2)


```

Make a heatmap of the top 100 most significantly differentially expressed genes: 
```
heatmap<-v$E[order(fit.cont$p.value[,1]),]

pheatmap(heatmap[1:250,],scale="row",show_rownames = F,
         main="top 250 Differentially expressed genes", 
         color = colorRampPalette(rev(brewer.pal(n = 11, name ="RdYlBu")))(100))
  
```

Make a heatmap of specifc target genes of interest: 

```
pheatmap(v$E[c("Hoxa7","Hoxa9","Lmo2","Tal1","Fos","Zfp36"),],border_color = NA)
```

### Boxplot of target genes

Easy default boxplots
```
nice.col <- brewer.pal(4,name="Dark2")

Hoxa7<-boxplot(v$E["Hoxa7",]~group,vertical=TRUE,las=2,cex.axis=1,pch=16,cex=1.3,
                 col=nice.col,method="jitter",ylab="Normalised log2 expression",
                 main="Hoxa7")
Hoxa9<-boxplot(v$E["Hoxa9",]~group,vertical=TRUE,las=2,cex.axis=1,pch=16,cex=1.3,
               col=nice.col,method="jitter",ylab="Normalised log2 expression",
               main="Hoxa9")
```
Harder (but nicer looking) ggplot2 and reshape2 boxplots:

```
Targets1<- as.data.frame(t(v$E[c("Hoxa7","Hoxa9","Lmo2","Tal1","Fos","Zfp36"),]))
Targets1$group<-as.character(c("A","A","A","B","B","B"))
Targets1$group <- factor(Targets1$group,levels = c("A","B"),ordered = TRUE)
df.t <- melt(Targets1, id.var = "group")

ggplot(data = df.t, aes(x=variable, y=value)) + geom_boxplot(aes(fill=group)) + theme_light() +
facet_wrap( ~ variable, scales="free",nrow =1) +scale_fill_brewer(palette="Set1",direction = -1)+ 
  xlab("") +  ylab("Normalised CPM") + ggtitle(" RNA")


```

### Volcano plot 

easy default version: 
```
volcanoplot(fit.cont,coef=1,highlight=5,names=rownames(fit.cont))
```
ggplot2 version: 
```
res <- topTable(fit.cont,coef = 1,adjust="fdr",n=Inf)
res[,7] <- ifelse((res$P.Val < 0.05 & abs(res$logFC) > 1), "red", "grey34") 
res[,8]<-row.names(res)
size <- ifelse((res$P.Val < 0.05 & abs(res$logFC) > 0.5), 2, 1)

g <- ggplot(data=res, aes(x=logFC, y=-log10(P.Value))) +
  geom_point(size = size,colour=res[,7]) + ggtitle("RNA-seq") + labs(x="logFC",y="-log10(p-value)") +
  geom_text_repel(
    data = res[1:20,],
    aes(label = V8),
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines")
  ) +
  xlim(c(-6, 6)) + ylim(c(0,6)) +
  theme_minimal()
g
```


### Barcode plots

If you have corresponding genesets (up and down) you could make a barcode plot. genelist is a character string of genenames in the geneset of interest. 

```
genesetup<-read.table(genelist1)
genesetdn<-read.table(genelist2)


barcodeplot(sigsorted$logFC,index = sigsorted$genename %in% 
              genesetup, index2 = sigsorted$genename %in% 
              genesetdn,
            main="RNAseq+genesets")
```  
You can also perform statistical analysis on the association between the gene sets and your gene list:

```
index = rownames(v$E) %in% genesetup
index2 = rownames(v$E) %in% genesetdn
mroast(y=v$E,index2,design=design,contrast=c(1,-1),nrot=10000)
```

### making a ranked list for GSEA analysis

Make a ranked list sorted by logFC (could also sort by t statistic) for GSEA
```
results <- topTable(fit.cont,coef = 1,number = 20000)
resultssorted <- results[order(results$logFC),]
x <- cbind(rownames(resultssorted),resultssorted$logFC)
x<-as.data.frame(x)

write.table(x,"ranklist.rnk",sep="\t",quote = F,col.names = F,row.names = F)
```

If your data is mouse, you'll need to convert to human gene names to use human MSigDB sets: 
```
genes <- x$V1
G_list <- getBM(filters= "external_gene_name", attributes= c("hsapiens_homolog_associated_gene_name","external_gene_name"),values=genes,mart= mart)
ResultsMerge<-merge(G_list,x,by.x="external_gene_name",by.y="V1")
ResultsMerge<- ResultsMerge[-which(ResultsMerge$hsapiens_homolog_associated_gene_name==""), ]
ResultsMerge1<-ResultsMerge[,-1]
ResultsOrder<-ResultsMerge1[order(ResultsMerge1$V2),]
write.table(ResultsOrder,"Human.rnk",sep="\t",quote = F,col.names = F,row.names = F)
```



## Acknowledgments


* **Madison Kelly** - *Author of this workflow* [madisonJK](https://github.com/madisonJK)

* Stephin Vervoort - *Helping with initial RNAseq analysis*


## Further information and useful tutorials

[Combine RNA-seq tutorial](http://combine-australia.github.io/RNAseq-R/)



## References

Liao, Y., Smyth, G. K., & Shi, W. (2013). The Subread aligner: fast, accurate and scalable read mapping by seed-and-vote. Nucleic Acids Research, 41(10), e108.

Moll, P., Ante, M., Seitz, A., & Reda, T. (2014). QuantSeq 3′ mRNA sequencing for RNA quantification. Nature Methods, 11(12).

