# M6A-seq

M6A-seq is a method for quantifying methylation of RNA.
This walkthrough will outline the general steps of analysing a M6A-seq dataset from fastq files through to figure generation. 


## Getting Started

The following walk-through provides an example of the analysis pipeline for m6a-seq data. For the purposes of this example, we will be starting with demultiplexed fastq files. These files have been sequenced paired end 150bp at approximately 20 million paired reads per sample. 

Both the m6a IP and the total RNA-seq (input) for each sample are needed for analysis. 

The software used have been referenced at the end, the majority of this analysis is based on the user manuals linked below, with only slight deviations from default settings. 

### Library preparation
 
For  details of m6a-RIP seq  wetlab protocol see [Dominissi et al 2013](https://www.nature.com/articles/nature11112) 


### General overview of m6a-seq analysis

![GeneralOverview](m6aOverview.png)


### Software and Version Control

STAR/2.5.3a

samtools/1.4.1

cutadapt /1.14

FastQC /0.11.5

RNAseqQC /1.1.8

Subread /1.5.0-p3

macs/2.1.1


R: 

exomePeak/2.16.0

MetPeak/1.0.0

### Optional QC

You can run FastQC first to determine the quality of the raw data. fastqc will be run on the trimmed files again below anyway:
[1_fastQC.sbatch](1_fastQC.sbatch)
```
module load fastqc

fastqc -o fastqc -f fastq --noextract -t 8 *.fastq.gz

```

## Trimming adapter sequences
The fastq files need to be trimmed for illumina adapter sequences, and filter the resulting reads based on both quality and length. This can be performed using cutadapt or trimgalore(preferred method). Both m6a-IP and input samples are trimmed and mapped the same. 

[2.1_trimgalore.sbatch](2.1_trimgalore.sbatch)

```
module load trimgalore
module load fastqc

trim_galore --paired -q 30 --length 30 --fastqc R1.fastq.gz R2.fastq.gz -o trimmed

```

*Important note: ensure that R1 and R2 files are balanced and in order following trimming or alignment will fail*



## Aligning Fastq Files
Fastq files are then aligned to the mouse genome (mm10) using STAR. The resulting sam files are sorted, converted to bam, duplicates removed and indexed with samtools. 

You will need to run star genomeGenerate first to generate the reference genome, ensure ```--sjdbOverhang``` is set to fragment size - 1.

Reference files used here are the default reference files in /data/

```
module load star
STAR --runMode genomeGenerate --sjdbOverhang 149 --runThreadN 8 --genomeDir ref/mm10STAR --genomeFastaFiles ../../../data/reference//indexes/mouse//mm10/fasta/Mus_musculus.GRCm38.dna.toplevel.fa --sjdbGTFfile ../../../data/reference/gtf/Mus_musculus.GRCm38.90.gtf
```

Next, the trimmed Fastq files are aligned with STAR, which also converts to Bam and sorts. These settings are are very low stringency for alignment, as the fragments appear small.

[3_STARalignmm10PE.sbatch](3_STARalignmm10PE.sbatch)

```
module load star
module load samtools
module load igvtools

STAR --outFileNamePrefix Sample1 --outFilterScoreMinOverLread 0.2 --outFilterMatchNminOverLread 0.2 --outFilterMismatchNmax 2 --outFilterMatchNmin 0 --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --genomeDir ref/mm10STAR --runThreadN 16 --readFilesIn R1_val_1.fq.gz R2_val_2.fq.gz

```
Samtools is then used to remove dupliates and to index: 

```
module load samtools

samtools rmdup -s sample.combinedAligned.sortedByCoord.out.bam sample.sortedByCoord.rmdup.bam

samtools index sample.sortedByCoord.rmdup.bam sample.sortedByCoord.rmdup.bam.bai

```

an IGV viewable TDF can be created with igvtools

*ext size should be fragment length- read length, however it is difficult to calculate fraglength of RNA derived reads so ext is currently set to 0*

```
module load igvtools

igvtools count -z 5 -w 25 sample.sortedByCoord.rmdup.bam sample.sortedByCoord.rmdup.bam.tdf mm10

```

*Optional:* QC can be performed on the resulting Bam files using RNAseqc. You will need to input a sample file, which is a tab-delimited text file with 3 columns specifying ID, the filename of bam file, and comments. You will need to download the rRNA reference fileand a gtf reference file from UCSC table browser in addition to a fasta file of the reference genome. see RNA-SeQC --help for more info. 

*?may need to add readgroups using picard before RNAseqc will work?*

```
module load rna-seqc
module load java

RNA-SeQC -s "ID|Sample.Aligned.sortedByCoord.out.bam|comments"  -t ref/mm10genes.gtf -o m6a/Fastq/RNAseqc -r data/reference/indexes/mouse/mm10/fasta/Mus_musculus.GRCm38.dna.toplevel.fa -BWArRNA ref/mm10rRNA.fa 

```


### Calling Peaks with MACS2 

Peaks can then be called in the IP sample relative to the input RNA seq sample using macs2. Ensure theGenome size is set to the calculated transcriptome size for the genome used. in this case mm10 transcriptome size is used ( as in [Dominissi et al. 2013](https://www.ncbi.nlm.nih.gov/pubmed/22575960) )

[5_MACS2m6a.sbatch](5_MACS2m6a.sbatch)
```
module load macs
macs2 callpeak -g mm -f BAMPE -t IPsample.Aligned.sortedByCoord.out.bam -c Input.Aligned.sortedByCoord.out.bam --nomodel --gsize 2.82e8 --extsize 149 --outdir m6aIP_macs2 -n m6aIP_normpeaks
```
### Calling Peaks with exomePeak

As the samples are paired (INPUT and IP from same cell pool) it is suggested that peaks are called first on individual replicates, then merged.

The reference .gtf file needs to have the same chromosome annotation as the aligned bam files ( i.e. with or without "chr") 

See [exomePeak](https://github.com/ZhenWei10/exomePeak) for more info. 

*exomePeak.R*

```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("exomePeak")

library(exomePeak)

mm10<-"path/mm10mrna.gtf"

ip1 <- "path/IPsample.Aligned.sortedByCoord.out.bam"

input1<- "path/Input.Aligned.sortedByCoord.out.bam"

result=exomepeak(GENE_ANNO_GTF=mm10,IP_BAM = ip1,INPUT_BAM = input1,
                 EXPERIMENT_NAME="Sample")
```
or with replicates:

```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("exomePeak")

library(exomePeak)

mm10<-"path/mm10mrna.gtf"

ip1.rep1 <- "path/IPsample.Aligned.sortedByCoord.out.bam"
ip1.rep2 <- "path/IPsample.2.Aligned.sortedByCoord.out.bam"

input1.rep1<- "path/Input.Aligned.sortedByCoord.out.bam"
input1.rep2<- "path/Input.2.Aligned.sortedByCoord.out.bam"


ip1<-c(ip1.rep1,ip1.rep2)
input1<-c(input1.rep1,input1.rep2)

result=exomepeak(GENE_ANNO_GTF=mm10,IP_BAM = ip1,INPUT_BAM = input1,
                 EXPERIMENT_NAME="Sample")
```


### Calling Peaks with meTPeak

*Recommended exome-based Peak-calling Software*

MetPeak is run very similar to exomePeak and provides similar results. However, improvements in statistical modelling improves peak calling for low input samples. See [metPeak](https://github.com/compgenomics/MeTPeak) for more info. 

*MetPeak.R*

```
library("devtools")
install_github("compgenomics/MeTPeak",build_opts = c("--no-resave-data", "--no-manual"))

library(MeTPeak)

mm10<-"path/mm10mrnaNOCHR.gtf"

ip1 <- "path/IPsample.Aligned.sortedByCoord.out.bam"

input1<- "path/Input.Aligned.sortedByCoord.out.bam"

metpeak(GENE_ANNO_GTF=gtf,IP_BAM = ip1,INPUT_BAM = input1,
        EXPERIMENT_NAME="Sample.rep1")

```


#### Annotate peaks 

**With homer**

```
annotatePeaks.pl IPsample.peaks.bed mm10 -annStats IPsample.annstats.txt > IPsample.peaks.ann.bed

```

*Note: Depending on the reference genome used, you may have to add "chr" to the bed file chromosome notation first using
```awk '{print "chr"$0}' peaks.bed > peaks.chr.bed```*


**with bedtools closest** 

```
bedtools closest -D -a IPsample.peaks.bed -b mm10TSS.bed > IPpeaksmm10TSS.bed
```

#### Motif analysis

**with homer**

```
findMotifsGenome.pl peakfile.bed mm10 peakMotifanalysis
```
*note: If using on the cluster may need to set '''-preparsedDir dir/'''

- Can also use a total RNA background control. Create Peak files from the matched RNA input samples using macs2 or Csaw and add this as ``-bg RNApeaks.bed```

### Counting Reads

The resulting bam files can also then be used as input into subread's *featureCounts* function in order to count the reads in each gene. Instead of counting to the mm10 exome, could also provide the above calculated peak bed files. 
*The unstranded option is generally used but this can be altered by changing -s 1 if the stranded option was used in the mapping steps above*

*featurecounts.sbatch*

```
module load subread

featureCounts -t exon -g gene_id -F 'GTF/GFF' -O -M -s 0 -p -T 8 -a data/reference/gtf/Mus_musculus.GRCm38.90.gtf  -o counts.txt *.sorted.bam

```

## Differential m6a Analysis

### Using Bedtools Intersect 

Determnine which peaks are overlapping in two replicates 

```
module load bedtools
bedtools intersect -a IPsampleREP1.peaks.bed -b IPsampleREP2.peaks.bed > samplecombined.peaks.bed
```

Then, determine which peaks are overlapping/unique between 2 samples

```
module load bedtools
bedtools intersect -wao -a samplecombined.peaks.bed -b treated.samplecombined.peaks.bed > diffpeaks.bed
```


### Quantitative differential peak analysis with replicates with MeTPeak 

Can be performed with MeTPeak, see user guide for more details. 

## Generating plots

### Metagene plots

Can be created using: 
- Deeptools 
- guitR 
- metaplotR 

## Acknowledgments


* **Madison Kelly** - *Author of readme* [madisonJK](https://github.com/madisonJK)
* **Yogesh Kumar** - *Developed initial m6a analysis for the Das lab*


## Further information and useful tutorials

[Combine RNA-seq tutorial](http://combine-australia.github.io/RNAseq-R/)

[exomePeak Github](https://github.com/ZhenWei10/exomePeak)

[metPeak Github](https://github.com/compgenomics/MeTPeak)

[Trumpet Github](https://github.com/skyhorsetomoon/Trumpet)


## References
Cui, X., Zhang, L., Meng, J., Rao, M. K., Chen, Y., & Huang, Y. (2015). MeTDiff: a novel differential RNA methylation analysis for MeRIP-seq data. IEEE/ACM transactions on computational biology and bioinformatics, 15(2), 526-534.

Cui, X., Meng, J., Zhang, S., Chen, Y., & Huang, Y. (2016). A novel algorithm for calling mRNA m 6 A peaks by modeling biological variances in MeRIP-seq data. Bioinformatics, 32(12), i378-i385.

Dominissini, D., Moshitch-Moshkovitz, S., Schwartz, S., Salmon-Divon, M., Ungar, L., Osenberg, S., et al. (2013). Topology of the human and mouse m6A RNA methylomes revealed by m6A-seq. Nature, 485(7397), 201–206.

Liao, Y., Smyth, G. K., & Shi, W. (2013). The Subread aligner: fast, accurate and scalable read mapping by seed-and-vote. Nucleic Acids Research, 41(10), e108.

Zhang, M., Li, Q., & Xie, Y. (2018). A Bayesian hierarchical model for analyzing methylated RNA immunoprecipitation sequencing data. Quantitative Biology, 6(3), 275-286.
Meng, J., Lu, Z., Liu, H., Zhang, L., Zhang, S., Chen, Y., ... & Huang, Y. (2014). A protocol for RNA methylation differential analysis with MeRIP-Seq data and exomePeak R/Bioconductor package. Methods, 69(3), 274-281.




