library(methylKit)
library(genomation)
library(genomation)
library(ggpubr)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(rgl)
library(biomaRt)
library(rowr)
library(dendextend)
library(BSgenome)
library("BSgenome.Mmusculus.UCSC.mm10")

##load in bam files and write to methylkit txt files for easy access later#####

sample.list<-list.files(path = "RRBS/", pattern = "*.sort.bam$", all.files = FALSE,
           full.names = FALSE, recursive = FALSE,
           ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)

PROCESS_BAM_FILES = TRUE


if( PROCESS_BAM_FILES ){
  for(sample_name in sample.list){ print(sample_name);
  my.methRaw = processBismarkAln( location = paste("RRBS/",sample_name,sep = ""),
                             sample.id=paste(sample_name), assembly="mm10", 
                            read.context="CpG", save.folder="/Users/kellymadison/Documents/R/MK26RNAseq/RRBS" )}}

##read in methylkit files 
methylkit.txt.list<-list.files( path = "RRBS/",pattern = "_CpG.txt$", all.files = FALSE,
           full.names = F, recursive = F,
           ignore.case = FALSE, include.dirs = F, no.. = T)

files = as.list( paste("RRBS/",methylkit.txt.list,sep = ""))

myobjall = methRead(files, sample.id=list("A1","A2","A3","A4","B4","B1","B2","B3"),
                    assembly="mm10", treatment=c(1,1,1,1,0,0,0,0)

#add chr for downstream analysis 
for (i in 1:8) {myobjall[[i]]$chr<-paste("chr",myobjall[[i]]$chr,sep = "") }

getMethylationStats(myobjall[[8]],plot=T,both.strands=FALSE)
getCoverageStats(myobjall[[7]],plot=TRUE,both.strands=FALSE)

filtered.myobj=filterByCoverage(myobjall,lo.count=10,lo.perc=NULL,hi.count=NULL,hi.perc=99.9)

methall = unite(filtered.myobj, destrand=TRUE, min.per.group=2L)
                    
#see correlation between samples 
getCorrelation(na.omit(methall),method = "pearson",plot = F)
corr<-getCorrelation(na.omit(methall),method = "pearson",plot = F)
corr<-read.table("correlation.txt",header = T, row.names = T)
pheatmap(corr)

#see clustering between samples
clusterSamples(methall, dist="correlation", method="ward", plot=TRUE)
PCASamples(methall, screeplot=TRUE)
PCASamples(na.omit(methall))



#####all cs annotation #####
gene.obj=readTranscriptFeatures("RRBS/mm10genes.bed")
allCann=annotateWithGeneParts(as(methall,"GRanges"),gene.obj)
allCTSS<-getAssociationWithTSS(allCann)
genomation::getTargetAnnotationStats(allCann,percentage=TRUE,precedence=TRUE)

genomation::plotTargetAnnotation(allCann,precedence=TRUE,
                     main="all c annotation")

genomation::getFeatsWithTargetsStats(allCann,percentage=TRUE)


cpg.obj=readFeatureFlank("RRBS/CpG_islands_mm10.bed", feature.flank.name=c("CpGi","shores")  )

AllCcpg=annotateWithFeatureFlank(as(methall,"GRanges"),
                                    cpg.obj$CpGi,cpg.obj$shores,
                                    feature.name="CpGi",flank.name="shores")
genomation::getTargetAnnotationStats(AllCcpg,percentage=TRUE,precedence=TRUE)
methall
genomation::plotTargetAnnotation(allCann,
                     main="")

genomation::getFeatsWithTargetsStats(AllCcpg,percentage=TRUE)

####calculate Diff C ######

myDiffuncorrected=calculateDiffMeth(methall,mc.cores=4)

myDiff20unc=getMethylDiff(myDiffuncorrected,difference=25,qvalue=0.05)

diffMethPerChr(myDiffuncorrected,plot=TRUE,qvalue.cutoff=0.01, meth.cutoff=25)

bedgraph(myDiff20unc, file.name = "DMC.bed", col.name = "meth.diff", unmeth = FALSE,
         log.transform = FALSE, negative = FALSE, add.on = "")

#####annotate diff Cs ####
annotateWithGeneParts(as(myDiff20unc,"GRanges"),gene.obj)

diffAnn=annotateWithGeneParts(as(myDiff20unc,"GRanges"),gene.obj)
TSSass<-getAssociationWithTSS(diffAnn)
getTargetAnnotationStats(diffAnn,percentage=TRUE,precedence=TRUE)

genomation::plotTargetAnnotation(diffAnn,precedence=TRUE,
                     main="differential methylation annotation")

genomation::getFeatsWithTargetsStats(diffAnn,percentage=TRUE)


hypo = getMethylDiff(myDiffuncorrected, difference = 25, qvalue = 0.05, 
                     type = "hypo",mc.cores = 4)
hyper = getMethylDiff(myDiffuncorrected, difference = 25, qvalue = 0.05, 
                      type = "hyper",mc.cores = 4)

hypoAnn<-annotateWithGeneParts(as(hypo,"GRanges"),gene.obj)
getTargetAnnotationStats(hypoAnn,percentage=TRUE,precedence=TRUE)
plotTargetAnnotation(hypoAnn,precedence=TRUE,
                     main="hypo c annotation")

hyperAnn<-annotateWithGeneParts(as(hyper,"GRanges"),gene.obj)
getTargetAnnotationStats(hyperAnn,percentage=TRUE,precedence=TRUE)
plotTargetAnnotation(hyperAnn,precedence=TRUE,
                     main="hyper c annotation")

######cpgi######


diffCpGann=annotateWithFeatureFlank(as(myDiff20unc,"GRanges"),
                                    cpg.obj$CpGi,cpg.obj$shores,
                                    feature.name="CpGi",flank.name="shores")

genomation::plotTargetAnnotation(diffCpGann,col=c("green","gray","white"),
                     main="")

getFeatsWithTargetsStats(diffCpGann,percentage=TRUE)

######-------####
CpGicounts=regionCounts(filtered.myobj,cpg.obj$CpGi)

getCoverageStats(CpGicounts[[2]],plot=TRUE,both.strands=FALSE)

methCpGi = unite(CpGicounts, destrand=TRUE, min.per.group=2L)

clusterSamples(methCpGi, dist="correlation", method="ward", plot=TRUE)
PCASamples(methCpGi)
mperc.meth.CpGi=percMethylation(methCpGi)
boxplot(mperc.meth.CpGi)
myDiffuncorrectedCpGi=calculateDiffMeth(methCpGi,mc.cores=4)
diffMethPerChr(myDiffuncorrectedCpGi,plot=TRUE,qvalue.cutoff=0.05, meth.cutoff=25)

cpgsig<- as.data.frame(mperc.meth.CpGi[myDiffuncorrectedCpGi$qvalue<0.05,])
boxplot(cpgsig)

myDiff25CpG=getMethylDiff(myDiffuncorrectedCpGi,difference=25,qvalue=0.05)

diffCpGAnn=annotateWithGeneParts(as(myDiff25CpG,"GRanges"),gene.obj)
TSSassdiffCpG<-getAssociationWithTSS(diffCpGAnn)

Diff25CpGAnn<-cbind(myDiff25CpG,TSSassdiffCpG)


mart <- useDataset("mmusculus_gene_ensembl", useMart("ensembl",host = "www.ensembl.org",
                                                    ensemblRedirect = FALSE))
genes <- Diff25CpGAnn$feature.name
G_list <- getBM(filters= "refseq_mrna", attributes= c("refseq_mrna","external_gene_name"),values=genes,mart= mart)

GenenamesDiff25CpG<-merge(Diff25CpGAnn,G_list,by.x="feature.name",by.y="refseq_mrna")



########summarise counts over promoters######
#can also summarise over exons, introns, or enhancers as well

promoters=regionCounts(filtered.myobj,gene.obj$promoters)

getCoverageStats(promoters[[2]],plot=TRUE,both.strands=FALSE)

methproms = unite(promoters, destrand=TRUE, min.per.group=2L)

getCorrelation(methproms,plot = F)
clusterSamples(methproms, dist="correlation", method="ward", plot=TRUE)
PCASamples(methproms, screeplot=TRUE)
PCASamples(methproms)
mperc.meth.proms=percMethylation(methproms)
boxplot(mperc.meth.proms)

myDiffuncorrectedproms=calculateDiffMeth(methproms,mc.cores=4)

promsig<- as.data.frame(mperc.meth.proms[myDiffuncorrectedproms$qvalue<0.05,])
boxplot(promsig)

boxplotshort<-as.data.frame(t(promsig))
boxplotshort$group<-c("A","A","A","A","B","B","B","B")
boxplotshort$group <- factor(boxplotshort$group,
                             levels = c("B","A"),ordered = TRUE)
boxplotmelt<- melt(boxplotshort, id.var = "group")

ggplot(data=subset(boxplotmelt, !is.na(value)), aes(x=group, y=value)) + geom_boxplot(aes(fill=group)) + 
  theme_light() + scale_fill_brewer(palette="Set1")+ xlab("") +
  ylab("percentage methylation")+ggtitle("Differentially methylated cytosines")
compare_means(value~group, boxplotmelt)

means <- aggregate(value~group, boxplotmelt, mean) 

pheatmap(promsig,border_color = NA,show_rownames = F,scale = "row")

myDiff25proms=getMethylDiff(myDiffuncorrectedproms,difference=25,qvalue=0.01)
diffMethPerChr(myDiffuncorrectedproms,plot=TRUE,qvalue.cutoff=0.05, meth.cutoff=25)

diffpromAnn=annotateWithGeneParts(as(myDiff25proms,"GRanges"),gene.obj)
TSSassdiffprom<-getAssociationWithTSS(diffpromAnn)

Diff25promAnn<-cbind(myDiff25proms,TSSassdiffprom)

mart <- useDataset("mmusculus_gene_ensembl", useMart("ensembl",host = "www.ensembl.org",
                                                     ensemblRedirect = FALSE))
genes <- Diff25promAnn$feature.name
G_list <- getBM(filters= "refseq_mrna", attributes= c("refseq_mrna","external_gene_name"),values=genes,mart= mart)

GenenamesDiff25prom<-merge(Diff25promAnn,G_list,by.x="feature.name",by.y="refseq_mrna")

#####average meth percentage accross diff gene feaatures #### 

mperc.meth.av<-cbind.fill(rowMeans(mperc.meth.proms[,5:8]),rowMeans(mperc.meth.exons[,5:8]),
           rowMeans(mperc.meth.CpGi[,5:8]),rowMeans(mperc.meth.introns[,5:8]),
           rowMeans(mperc.meth.LTen[,5:8])
)
colnames(mperc.meth.av)<-c("proms","exons","CpGi","introns","HSCenhancers")
boxplot(log(mperc.meth.av),pch = 20,cex = 0.4)

######Calculate percentage methylation#######

mperc.meth=percMethylation(methall)
allmeth<-getData(methall)
locperc.meth<-cbind(allmeth[,1:3],mperc.meth)


diffAnn=annotateWithGeneParts(as(locperc.meth,"GRanges"),gene.obj)
TSSass<-getAssociationWithTSS(diffAnn)
getTargetAnnotationStats(diffAnn,percentage=TRUE,precedence=TRUE)

anlocperc.meth<-merge(TSSass,locperc.meth,by.x = "target.row",by.y = "row.names")

ordered<-as.data.frame(mperc.meth[order(myDiffuncorrected$qvalue),])
allsig<- as.data.frame(mperc.meth[myDiffuncorrected$qvalue<0.001&abs(myDiffuncorrected$meth.diff)>25,])

boxplot(allsig[,c(5:8,1:4)])

kmeans(allsig,4)
pheatmap(ordered[1:100,],border_color = NA,show_rownames = F)
pheatmap(allsig, border_color = NA,show_rownames = F)

my_hclust_gene <- hclust(dist(allsig), method = "complete")
list <- as.data.frame(cutree(tree = as.dendrogram(my_hclust_gene), k = 6))

my_gene_col <- data.frame(cluster = ifelse(test = list == 1, yes = "cluster 1", no = "cluster 2"))
my_gene_col$cluster<-paste("cluster",list$`cutree(tree = as.dendrogram(my_hclust_gene), k = 6)`,sep = " ")
my_gene_col$cluster<-factor(my_gene_col$cluster)
my_gene_col<-as.data.frame(my_gene_col[,-1])

pheatmap(allsig, border_color = NA,show_rownames = F,annotation_row = my_gene_col)
as.dendrogram(my_hclust_gene) %>%
  plot(horiz = TRUE)


###ggplot boxplots####

#boxplot of percentage methylation of significant DMCs

boxplotshort<-as.data.frame(t(allsig))
boxplotshort$group<-c("A","A","A","A","B","B","B","B")
boxplotshort$group <- factor(boxplotshort$group,
                                 levels = c("UT","DOX"),ordered = TRUE)
boxplotmelt<- melt(boxplotshort, id.var = "group")

ggplot(data = boxplotmelt, aes(x=group, y=value)) + geom_boxplot(aes(fill=group)) + 
  theme_light() + scale_fill_brewer(palette="Set1")+ xlab("") +
  ylab("percentage methylation")+ggtitle("Differentially methylated cytosines")
compare_means(value~group, data=subset(boxplotmelt, !is.na(value)))
t.test 


#boxplot total methylation percentage overall##

boxplotshort<-as.data.frame(t(as.data.frame(mperc.meth)))
boxplotshort$group<-c("A","A","A","A","B","B","B","B")
boxplotshort$group <- factor(boxplotshort$group,
                             levels = c("A","B"),ordered = TRUE)
boxplotmelt<- melt(boxplotshort, id.var = "group")

ggplot(data=subset(boxplotmelt, !is.na(value)) , aes(x=group, y=value)) + geom_boxplot(aes(fill=group)) +
  theme_light() + scale_fill_brewer(palette="Set1")+ xlab("") +
  ylab("percentage methylation")+ggtitle("All measured cytosine methylation")
boxplot(mperc.meth)

compare_means(value~group, boxplotmelt,method = "wilcox.test")
  
######circos Plot ######

chrom.length = seqlengths(Mmusculus)  
chr.len = chrom.length[grep("_|chrM|chrX|chrY", names(chrom.length), invert = T)] 

ideoDMC(myDiffuncorrected, chrom.length = chr.len, difference = 20, qvalue = 0.05,circos = T, title = "test", hyper.col = "red", hypo.col = "blue")


ideoDMC <- function(myDiff, chrom.length, difference = 25, 
                    qvalue = 0.01, circos = FALSE, title = "test", hyper.col = "magenta", 
                    hypo.col = "green") {
  require(methylKit)
  require(GenomicRanges)
  require(ggbio)
  
  # chrom.length
  myIdeo <- GRanges(seqnames = names(chrom.length), ranges = IRanges(start = 1, width = chrom.length))
  seqlevels(myIdeo) = names(chrom.length)
  seqlengths(myIdeo) = (chrom.length)
  
  hypo = getMethylDiff(myDiff, difference = difference, qvalue = qvalue, 
                           type = "hypo",mc.cores = 4)
  hyper = getMethylDiff(myDiff, difference = difference, qvalue = qvalue, 
                            type = "hyper",mc.cores = 4)
  
  g.per = as(hyper, "GRanges")
  seqlevels(g.per,pruning.mode = "coarse") = seqlevels(myIdeo)
  seqlengths(g.per)=(chr.len)
  
  g.po = as(hypo, "GRanges")
  seqlevels(g.po,pruning.mode = "coarse") = seqlevels(myIdeo)
  seqlengths(g.po)=(chr.len)
  
  values(g.po)$id = "hypo"
  values(g.per)$id = "hyper"
  
  if (circos) {
    
    p <- ggplot() + layout_circle(myIdeo, geom = "ideo", fill = "gray70", radius = 39, trackWidth = 2)
    
    
    p <- p + layout_circle(c(g.po, g.per), geom = "point", size = 0.1, aes(x = midpoint, y = meth.diff, color = id), radius = 25, trackWidth = 30) + scale_colour_manual(values = c(hyper.col, hypo.col))
    p + layout_circle(myIdeo, geom = "text", aes(label = seqnames),vjust = 0, radius = 55, trackWidth = 7) + labs(title = "Dox vs UT")
    
  } else {
    
    p <- ggplot()+  layout_karyogram(myIdeo, cytoband = FALSE)
    p  + layout_karyogram(c(g.po, g.per), geom = "point", size = 1, 
                          aes(x = midpoint, 
                              y = meth.diff, color = id)) + scale_colour_manual(values = c(hyper.col, 
                                                                                           hypo.col)) + labs(title = title)
    # new alternative commented out
    #autoplot(c(g.po, g.per), layout = "karyogram", geom = "point", size = 0.65, 
    #aes(x = midpoint,y = meth.diff, color = id))   scale_colour_manual(values = c(hyper.col, 
    #                                                                                        hypo.col))  labs(title = title)
    
  }
}


