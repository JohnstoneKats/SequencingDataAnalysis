# based on combine RNA seq workshop- see notes here: https://combine-australia.github.io/2018-09-26-RNAseq-Melbourne/
#RNA seq analysis from counts file generated from feature counts
#Generates following figures: PCA plot, MDS plot, heatmaps, volcano plots, individual gene counts plots
#Need to generate a sample and contrast matrix specific to your experiment

#load packages

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

#library("biomaRt", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.2")
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

#read in your counts table
setwd("")
seqdata <- read.delim("counts.txt", stringsAsFactors = FALSE)
#remove columns with extra data
countdata <- seqdata[,-(1:6)]
#rename rows with gene names 
rownames(countdata) <- seqdata[,1]
#rename columns with easier Identifiers - need to edit based on your sample names 
colnames(countdata)<-c("A1","A2","A3","B1","B2","B3")

#if your counts are annotated with refseq_mrna IDs use this code to covert to gene names 
mart <- useDataset("mmusculus_gene_ensembl", useMart("ensembl"))
genes <- rownames(countdata)
G_list <- getBM(filters= "refseq_mrna", attributes= c("refseq_mrna","external_gene_name"),values=genes,mart= mart)
seqdatamerge<-merge(countdata,G_list,by.x="row.names",by.y="refseq_mrna")

#will need to aggregate as multiple refseq IDs per gene
aggdata <- aggregate(seqdatamerge,by=list(seqdatamerge$external_gene_name),FUN = "mean")

#get rid of extra data columns (need to edit depending on number of samples)
allcounts<-aggdata[,c(-1,-2,-3,-9)]
rownames(allcounts)<-aggdata[,1]

#write table of final counts annotated with gene names for easy reference later
write.table(allcounts,file= "GenenameCounts.txt",sep = "\t",quote = F)

#set a threshold and filter out genes with low/no counts
myCPM <- cpm(allcounts)
thresh <- myCPM > 4
keep <- rowSums(thresh) >= 2
counts.keep <- allcounts[keep,]


#create DGE list
y <- DGEList(counts.keep)
y <- calcNormFactors(y)
logcounts <- cpm(y,log=TRUE)

plotMDS(y)


#make sample info ncol = n independent variables
SI<-matrix(data = NA,nrow=ncol(counts.keep),ncol=2)
SI[,1]<-colnames(counts.keep)
SI[,2]<-c("A","A","A","B","B","B")
colnames(SI)<-c("Samplename","Variable")
SI<-as.data.frame(SI)


#create concatination of group for design matrix
group <- paste(SI$Variable,sep=".")
group <- factor(group)
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)

#Voom transform the data (normalise)
v <- voom(y,design,plot = TRUE, normalize.method = "quantile")

#fit data to linear model
fit <- lmFit(v)


#make a contrast matrix& fit data #input comparison you want to make
cont.matrix <- makeContrasts(treatmentvscontrol =B-A,
                             levels=design)
fit.cont <- contrasts.fit(fit, cont.matrix)
fit.cont <- eBayes(fit.cont)
summa.fit <- decideTests(fit.cont,coef=1)


#create a table of results sorted by pvalue 
sigsorted<- topTable(fit.cont,coef=1,sort.by="p",n="Inf")


#make a MD and volcano plot
plotMD(fit.cont,coef=1,status=summa.fit)
volcanoplot(fit.cont,coef=1,highlight=5,names=rownames(fit.cont))

#plot normalised counts of individual genes
nice.col <- brewer.pal(4,name="Dark2")

Hoxa7<-boxplot(v$E["Hoxa7",]~group,vertical=TRUE,las=2,cex.axis=1,pch=16,cex=1.3,
                 col=nice.col,method="jitter",ylab="Normalised log2 expression",
                 main="Hoxa7")
Hoxa9<-boxplot(v$E["Hoxa9",]~group,vertical=TRUE,las=2,cex.axis=1,pch=16,cex=1.3,
               col=nice.col,method="jitter",ylab="Normalised log2 expression",
               main="Hoxa9")


#write a file of your results table
results <- topTable(fit.cont,coef = 1,number = 20000)
write.table(results,file = "results.txt",sep = "\t",quote = F,col.names = T)

#make a heatmap of most significanly differentially expressed genes
resultssorted <- results[order(results$logFC),]
filtered <- results[(results$P.Value <0.05)& ((results$logFC > 0.5) | (results$logFC < (-0.5))),]
filteredgenes<-rownames(filtered)
filteredheatmap<-v$E[filteredgenes,] 

pheatmap(filteredheatmap,scale="row",
         trace="none",main="Differentially expressed genes",
         border_color = NA,fontsize_row = 2)

heatmap<-v$E[order(fit.cont$p.value[,1]),]

pheatmap(heatmap[1:250,],scale="row",show_rownames = F,
         main="top 250 Differentially expressed genes", 
         color = colorRampPalette(rev(brewer.pal(n = 11, name ="RdYlBu")))(100))


#make boxplots of specific genes you're interested in
Targets1<- as.data.frame(t(v$E[c("Hoxa7","Hoxa9","Lmo2","Tal1","Fos","Zfp36"),]))
Targets1$group<-as.character(c("A","A","A","B","B","B"))
Targets1$group <- factor(Targets1$group,levels = c("A","B"),ordered = TRUE)
df.t <- melt(Targets1, id.var = "group")

ggplot(data = df.t, aes(x=variable, y=value)) + geom_boxplot(aes(fill=group)) + theme_light() +
facet_wrap( ~ variable, scales="free",nrow =1) +scale_fill_brewer(palette="Set1",direction = -1)+ 
  xlab("") +  ylab("Normalised CPM") + ggtitle(" RNA")

compare_means(c(Hoxa7,Hoxa9,Lmo2,Tal1,Fos,Zfp36)~group,Targets1,method = "t.test", paired = FALSE,
                             group.by = NULL, ref.group = NULL)

#make a heatmap of specific target genes 
pheatmap(v$E[c("Hoxa7","Hoxa9","Lmo2","Tal1","Fos","Zfp36"),],border_color = NA)


#make a ranked list sorted by logFC (could also sort by t statistic) for GSEA
results <- topTable(fit.cont,coef = 1,number = 20000)
resultssorted <- results[order(results$logFC),]
x <- cbind(rownames(resultssorted),resultssorted$logFC)
x<-as.data.frame(x)

write.table(x,"ranklist.rnk",sep="\t",quote = F,col.names = F,row.names = F)


# is mouse RNA seq, will need to convert to human to use MSigDB sets
genes <- x$V1
G_list <- getBM(filters= "external_gene_name", attributes= c("hsapiens_homolog_associated_gene_name","external_gene_name"),values=genes,mart= mart)
ResultsMerge<-merge(G_list,x,by.x="external_gene_name",by.y="V1")
ResultsMerge<- ResultsMerge[-which(ResultsMerge$hsapiens_homolog_associated_gene_name==""), ]
ResultsMerge1<-ResultsMerge[,-1]
ResultsOrder<-ResultsMerge1[order(ResultsMerge1$V2),]
write.table(ResultsOrder,"MK27Human.rnk",sep="\t",quote = F,col.names = F,row.names = F)


#if you have corresponding genesets (up and down) could make a barcode plot


genesetup<-read.table(#genesethere)
genesetdn<-read.table(#genesethere)


barcodeplot(results$logFC,index = results$genename %in% 
              genesetup, index2 = results$genename %in% 
              genesetdn,
            main="RNAseq+genesets")
  
#statistical analysis of above
index = rownames(v$E) %in% genesetup
index2 = rownames(v$E) %in% genesetdn
mroast(y=v$E,index2,design=design,contrast=c(1,-1),nrot=10000)



#volcanoPlot


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
