# Reduced Representation Bisulfite Sequencing

Method for the base pair resolution detection of methylated cytosines on a CGI-enriched sample of the genome. 

## Getting Started

The following walk-through provides an example of the analysis pipeline for RRBS data. For the purposes of this example, we will be starting with demultiplexed fastq files. These files have been sequenced paired end 75bp at approximately 20million paired reads per sample. 

### Library preparation
 
For a good explanation of library preparation click [here](https://www.epigenesys.eu/images/stories/protocols/pdf/20160127163832_p70.pdf)

In short, genomic DNA was digested with Msp1, Klenow (exo 5'-3'-) was used for end repair and poly-A tailing. Illumina truseq (methylated) DNA Adapters were ligated and samples were pooled. Zymo Methylgold kit was used for bisulfite conversion, Then library amplification was performed with PFU turbo CX. resulting libraries were sizeselected using a pippenprep for fragments between 250-500bp and sequenced on an Illumina Nextseq500. 

### General overview of analysis

![GeneralOverview](https://raw.githubusercontent.com/madisonJK/ReferenceAnalysis/RRBS/RRBS_Overview.png)


### Software and Version Control

Bismark/0.18.1

Bowtie2 /2.3.3

Trimgalore!/0.4.4

samtools/1.4.1

cutadapt /1.14

FastQC/0.11.5

R: Methylkit/1.11.0

ggplot2/3.3.1

reshape2/1.4.3

BSgenome/1.52.0

pheatmap/1.0.12


## Aligning Fastq Files

### Trimming fastq files

The first step is to trim the ends of the reads, trimming 2bp off the 3' end or the 5' end off the strand negated the bias introduced by end filling with unmethylated C's during the end repair step. 


```
module load trimgalore

s=$(basename $1)
t=$(dirname $1)
x=`echo $s | cut -d "_" -f 1-2`

trim_galore --paired --rrbs ${t}/${x}_R1_001.fastq.gz ${t}/${x}_R2_001.fastq.gz -o ${t}
```

### Aligning fastq files

The next step is to align to the genome, and assess the methylation percentage
see here: for additional information on how bismark works: 

The input for this step is the trimmed and validated files generated from trim_galore! above. 
The method outlined below first aligns in a paired end manner, then saves discordant reads as separate R1 and R2 fastq ffiles, then alignes these leftover reads using single end settings, ensuring that R2 is aligned to the reverse compliment

```
module load bismark

s=$(basename $1)
t=$(dirname $1)
x=`echo $s | cut -d "_" -f 1-2`

echo $s
echo $x

bismark --unmapped -X 1000 --multicore 4 ../../../data/reference/indexes/mouse/mm10/bowtie2/ -1 ${t}/${x}_R1_001_val_1.fq.gz -2 ${t}/${x}_R2_001_val_2.fq.gz

bismark --multicore 4 ../../../data/reference/indexes/mouse/mm10/bowtie2/ ${x}_R1_001_val_1.fq.gz_unmapped_reads_1.fq.gz

bismark --multicore 4 --pbat ../../../data/reference/indexes/mouse/mm10/bowtie2 ${x}_R2_001_val_2.fq.gz_unmapped_reads_2.fq.gz

```

This will leave you with 3 files (paired end alignments, and single end alignments. these are then pasted together, sorted and converted to bam files and indexed using samtools. 

```
module load samtools

samtools merge ${x}merged.bam ${x}_R1_001_val_1_bismark_bt2_pe.bam ${x}_R1_001_val_1.fq.gz_unmapped_reads_1_bismark_bt2.bam ${x}_R2_001_val_2.fq.gz_unmapped_reads_2_bismark_bt2.bam

samtools sort ${x}merged.bam > ${x}mergedsort.bam

samtools index ${x}mergedsort.bam
```

### generating bedgraphs for Visualisation

Input is the merged .sam files generated above. Output can be loaded into IGV or UCSC genome browser. 

```
bismark_methylation_extractor --bedgraph ${1}
```

## Differential Methylation Analysis

### DMC Analysis with MethylKit

The remaining analysis and plots can be performed in R 

first, load in bam files and write to methylkit.txt files for easy access later. 
.BAM files will need to be in the sam folder as indexed .BAI files

```
library(methylKit)
library(genomation)

sample.id<-list("A1","A2","A3","A4","B7","B8","B9","B10")

sample.list<-list.files(path = "RRBS/", pattern = "*.sort.bam$", all.files = FALSE,
           full.names = FALSE, recursive = FALSE,
           ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)

PROCESS_BAM_FILES = TRUE

if( PROCESS_BAM_FILES ){
  for(sample_name in sample.list){ print(sample_name);
  my.methRaw = processBismarkAln( location = paste("RRBS/",sample_name,sep = ""),
                             sample.id=paste(sample_name), assembly="mm10", 
                            read.context="CpG", save.folder="yourfolder/" )}}
```

read in methylkit files 

```
methylkit.txt.list<-list.files( path = "RRBS/",pattern = "_CpG.txt$", all.files = FALSE,
           full.names = F, recursive = F,
           ignore.case = FALSE, include.dirs = F, no.. = T)

files = as.list( paste("RRBS/",methylkit.txt.list,sep = ""))

myobjall = methRead(files, sample.id=list("A","A","A","A","B","B","B","B"),
                    assembly="mm10", treatment=c(1,1,1,1,0,0,0,0)

#optional: add chr for downstream analysis if your reference genome lacks "chr"
for (i in 1:8) {myobjall[[i]]$chr<-paste("chr",myobjall[[i]]$chr,sep = "") }
```

Perform some general QC 

```
getMethylationStats(myobjall[[1]],plot=T,both.strands=FALSE)
getCoverageStats(myobjall[[2]],plot=TRUE,both.strands=FALSE)
```

Filter out C's with low coverage (<10) and create a single object

```
filtered.myobj=filterByCoverage(myobjall,lo.count=10,lo.perc=NULL,hi.count=NULL,hi.perc=99.9)
methall = unite(filtered.myobj, destrand=TRUE, min.per.group=2L)

```
plot correlation and clustering

```
getCorrelation(na.omit(methall),method = "pearson",plot = F)
corr<-getCorrelation(na.omit(methall),method = "pearson",plot = F)
pheatmap(corr)

clusterSamples(methall, dist="correlation", method="ward", plot=TRUE)
PCASamples(methall, screeplot=TRUE)
```
annotate C's with genomic parts or CGIs

```
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

```
Calculate differential methylation of Cs, use 20% differential methylation and q value of 0.05 for cutoff

```

myDiffuncorrected=calculateDiffMeth(methall,mc.cores=4)


myDiff20unc=getMethylDiff(myDiffuncorrected,difference=25,qvalue=0.05)

```
plot differential methylation per chromosome 
```
diffMethPerChr(myDiffuncorrected,plot=TRUE,qvalue.cutoff=0.01, meth.cutoff=25)
```
save differential C's as a bedgraph
```
bedgraph(myDiff20unc, file.name = "MK26dmc.bed", col.name = "meth.diff", unmeth = FALSE,
         log.transform = FALSE, negative = FALSE, add.on = "")
```
find hypo/hypermethylated cytosines 

```
hypo = getMethylDiff(myDiffuncorrected, difference = 25, qvalue = 0.05, 
                     type = "hypo",mc.cores = 4)
hyper = getMethylDiff(myDiffuncorrected, difference = 25, qvalue = 0.05, 
                      type = "hyper",mc.cores = 4)
```

### DMR Analysis with MethylKit

For some analysis, looking at differentially methylated regions may make more sense. Count the methylated C's within a specificied window size (here it is a window sie of 1000bp) and combine all samples into a single table. 

```
methtiled = tileMethylCounts(filtered.myobj, cov.bases = 2, win.size = 1000, step.size = 1000)
methtiledunite<-unite(methtiled, destrand=TRUE, min.per.group=2L)
```

```methtiledunite``` can then be used the same as the ```methall``` object 


### Generating a CIRCOS plot

Create th ideoDMC function for plotting:

```
library(BSgenome)
library("BSgenome.Mmusculus.UCSC.mm10")

chrom.length = seqlengths(Mmusculus)  
chr.len = chrom.length[grep("_|chrM|chrX|chrY", names(chrom.length), invert = T)] 


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
    
  }
}

```

Use the ```myDiffuncorrected``` object generated above as input. Use circos = F for a linear chromosome plot.

```
ideoDMC(myDiffuncorrected, chrom.length = chr.len, difference = 20, qvalue = 0.05,circos = T, title = "test", hyper.col = "red", hypo.col = "blue")

```

### Generating a Heatmap

First, calculate percentage methylation:
```
mperc.meth=percMethylation(methall)

```

Then, plot the top 100 most significantly differentially methylated Cytosines (or regions) with pheatmap: 

```
ordered<-as.data.frame(mperc.meth[order(myDiffuncorrected$qvalue),])
pheatmap(ordered[1:100,],border_color = NA,show_rownames = F)
```

### Generating a boxplot

```
library(ggpubr)
library(reshape2)
library(ggplot2)

allsig<- as.data.frame(mperc.meth[myDiffuncorrected$qvalue<0.001&abs(myDiffuncorrected$meth.diff)>25,])
boxplotshort<-as.data.frame(t(allsig))
boxplotshort$group<-c("A","A","A","A","B","B","B","B")
boxplotshort$group <- factor(boxplotshort$group,
                                 levels = c("B","A"),ordered = TRUE)
boxplotmelt<- melt(boxplotshort, id.var = "group")

ggplot(data = boxplotmelt, aes(x=group, y=value)) + geom_boxplot(aes(fill=group)) + 
  theme_light() + scale_fill_brewer(palette="Set1")+ xlab("") +
  ylab("percentage methylation")+ggtitle("Differentially methylated cytosines")
compare_means(value~group, data=subset(boxplotmelt, !is.na(value)))
```


## Author of this workflow

* ** Madison Kelly ** - [madisonJK](https://github.com/madisonJK)


## Acknowledgments

* Stephin Vervoort - *establishing wetlab protocol and helping with initial RRBS analysis*


## Further information and useful tutorials

[Babraham RRBS guide](http://www.bioinformatics.babraham.ac.uk/projects/bismark/RRBS_Guide.pdf)

[Bismark User guide](https://www.bioinformatics.babraham.ac.uk/projects/bismark/Bismark_User_Guide.pdf)

[Trimming RRBS reads with TrimGalore!!](http://www.bioinformatics.babraham.ac.uk/projects/trim_galore/trim_galore_User_Guide_v0.3.7.pdf)

[QC, trimming and alignment of RRBS](https://www.epigenesys.eu/images/stories/protocols/pdf/20120720103700_p57.pdf)

[Preparation of RRBS libraries](https://www.epigenesys.eu/images/stories/protocols/pdf/20160127163832_p70.pdf)

[MethylKit User guide](https://bioconductor.org/packages/release/bioc/vignettes/methylKit/inst/doc/methylKit.html)



## References

