# 3' M6A-seq

M6A-seq is similar to meRIP-seq and is a method for quantifying methylation of RNA.
This walkthrough will outline the general steps of analysing a M6A-seq dataset from fastq files through to figure generation. 


## Getting Started

The following walk-through provides an example of the analysis pipeline for m6a-seq data. For the purposes of this example, we will be starting with demultiplexed fastq files. These files have been sequenced paired end 150bp at approximately 20 million paired reads per sample. 

The software used have been referenced at the end, the majority of this analysis is based on the user manuals linked below, with only slight deviations from default settings. 

### Library preparation
 
For  details of m6a-RIP seq  wetlab protocol see: 


### General overview of RNA-seq analysis

![GeneralOverview](https://github.com/madisonJK/ReferenceAnalysis/raw/master/RNA-seq/RNA-seq_Analysis_Overview.png)


### Software and Version Control

STAR

samtools/1.4.1

cutadapt /1.14

FastQC /0.11.5

RNAseqQC /1.1.8

Subread /1.5.0-p3

MACS


R: 

Limma/

Deseq/ 

exomePeak

MetPeak

ggplot2/3.3.1

reshape2/1.4.3

pheatmap/1.0.12

## Trimming adapter sequences
The fastq files need to be trimmed for illumina adapter sequences. This can be performed using cutadapt or trimgalore: 

```
module load cutadapt
s=$(basename $1)
x=`echo $s | cut -d "_" -f 1-2`
t=$(dirname $1)

echo ${1}

cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -G AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -q 30,30 -m 30 -o ${t}/${x}_R1.trimmed.fastq.gz -p ${t}/${x}_R2.trimmed.fastq.gz  ${t}/${x}_R1.fastq.gz ${t}/${x}_R2.fastq.gz
```

*Important note: ensure that R1 and R2 files are balanced and in order following trimming or alignment will fail*

## Aligning Fastq Files
Fastq files are aligned to the Mouse genome (mm10) using STAR. The resulting sam files are sorted, converted to bam and indexed with samtools. 

You will need to run star genomeGenerate first to generate the reference genome:

```
module load star
STAR --runMode genomeGenerate --runThreadN 8 --genomeDir ref/mm10STAR --genomeFastaFiles ../../../data/reference//indexes/mouse//mm10/fasta/Mus_musculus.GRCm38.dna.toplevel.fa --sjdbGTFfile ../../../data/reference/gtf/Mus_musculus.GRCm38.90.gtf

```

Next, align the trimmed Fastq files with STAR, then convert to Bam and index with samtools. These settings are are very low stringency for alignment, as the fragments appear small.

```
module load star
module load samtools
module load igvtools

s=$(basename $1)
t=$(dirname $1)
x=`echo $s | cut -d "_" -f 1-2`
echo $s
echo $x


STAR --outFileNamePrefix ${t}/${x} --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMismatchNmax 2 --outFilterMatchNmin 0 --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --genomeDir ref/mm10STAR --runThreadN 8 --readFilesIn ${t}/${x}_R1.trimmed.fastq.gz.sorted.fastq.gz ${t}/${x}_R2.trimmed.fastq.gz.sorted.fastq.gz


#samtools view -@ 8 -Sbo ${t}/${x}.sam.bam  ${t}/${x}.sam
#samtools sort -@ 8 -o .${t}/${x}sam.sorted.bam ${t}/${x}.sam.bam
samtools rmdup -s ${t}/${x}Aligned.sortedByCoord.out.bam ${t}/${x}.rmdup.bam
samtools index ${t}/${x}.rmdup.bam ${t}/${x}.rmdup.bam.bai
```

*Optional:* QC can be performed on the resulting Bam files using RNAseqc. You will need to input a sample file, which is a tab-delimited text file with 3 columns specifying ID, the filename of bam file, and comments. You will need to download the rRNA reference fileand a gtf reference file in addition to a fasta file of the reference genome. see RNA-SeQC --help for more info. 

```
module load rna-seqc
module load java

RNA-SeQC -s "ID|filename.bam|comments"  -t ref/mm10genes.gtf -o m6a/Fastq/RNAseqc -r data/reference/indexes/mouse/mm10/fasta/Mus_musculus.GRCm38.dna.toplevel.fa -BWArRNA ref/mm10rRNA.fa ${1}

```

IGV tools can also be used to generate a .tdf for visualisation with IGV

```
module load igvtools

igvtools count -z 5 -w 25 -e 225 ${t}/${x}.rmdup.bam ${t}/${x}.tdf mm10

```

Alternatively, you can also make wiggle file using STAR for visulisation with IGV:

```
STAR --runMode inputAlignmentsFromBAM --outFileNamePrefix ${t}/${x} --bamRemoveDuplicatesType UniqueIdentical --outWigType wiggle --outWigStrand Stranded ${t}/${x}Aligned.sortedByCoord.out.bam
```


### Calling Peaks with MACS2 

```
module load macs
macs2 callpeak -g mm -f BAMPE -t $1 -c $2 --nomodel --gsize 2.82e8 --extsize 200 --outdir ${s}_macs2 -n ${s}_normpeaks
```

### Counting Reads

The resulting bam files can then be used as input into subread's *featureCounts* function:
*The unstranded option is generally used to improve mapping rate, but this can be altered by changing -s 1*

```
module load subread

featureCounts -t exon -g gene_id -F 'GTF/GFF' -O -M -s 0 -a data/reference/gtf/Mus_musculus.GRCm38.90.gtf  -o counts.txt *.sorted.bam
```

## Differential RNA Methylation Analysis

### Calling peaks with exomePeak

The remaining analysis and figure generation can be performed in R. 


```
source("http://bioconductor.org/biocLite.R")
biocLite("edgeR")

library(ggrepel)

```

Read in the counts table generated above: 

```

```

## Generating plots

MD plot
```
plotMD(fit.cont,coef=1,status=summa.fit)

```



## Acknowledgments


* **Madison Kelly** - *Author of readme* [madisonJK](https://github.com/madisonJK)
* **Yogesh Kumar** - *Developed initial m6a analysis for the Das lab*




## Further information and useful tutorials

[Combine RNA-seq tutorial](http://combine-australia.github.io/RNAseq-R/)




## References
Cui, X., Zhang, L., Meng, J., Rao, M. K., Chen, Y., & Huang, Y. (2015). MeTDiff: a novel differential RNA methylation analysis for MeRIP-seq data. IEEE/ACM transactions on computational biology and bioinformatics, 15(2), 526-534.

Cui, X., Meng, J., Zhang, S., Chen, Y., & Huang, Y. (2016). A novel algorithm for calling mRNA m 6 A peaks by modeling biological variances in MeRIP-seq data. Bioinformatics, 32(12), i378-i385.


Dominissini, D., Moshitch-Moshkovitz, S., Schwartz, S., Salmon-Divon, M., Ungar, L., Osenberg, S., et al. (2013). Topology of the human and mouse m6A RNA methylomes revealed by m6A-seq. Nature, 485(7397), 201–206.

Liao, Y., Smyth, G. K., & Shi, W. (2013). The Subread aligner: fast, accurate and scalable read mapping by seed-and-vote. Nucleic Acids Research, 41(10), e108.

Zhang, M., Li, Q., & Xie, Y. (2018). A Bayesian hierarchical model for analyzing methylated RNA immunoprecipitation sequencing data. Quantitative Biology, 6(3), 275-286.
Meng, J., Lu, Z., Liu, H., Zhang, L., Zhang, S., Chen, Y., ... & Huang, Y. (2014). A protocol for RNA methylation differential analysis with MeRIP-Seq data and exomePeak R/Bioconductor package. Methods, 69(3), 274-281.




