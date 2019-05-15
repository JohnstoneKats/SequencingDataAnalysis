# ChIP-seq

This walkthrough will outline the general steps of analysing a ChIP-seq dataset from fastq files through to figure generation. 
see ```--help``` for more availavle options at each step. 

## Getting Started

The following walk-through provides an example of the analysis pipeline for ChIP-seq data. For the purposes of this example, we will be starting with demultiplexed fastq files. These files have been sequenced single ended 75bp at approximately 20 million paired reads per sample. An input control will be used below, but you could also use an IgG contorol as well. 

The PeterMac bioinformatics core automatically applies Seqliner2.0 pipeline to the Fastq files to produce .bam and .bai files. if you wish to start with these skip to the **Peak Calling** section.

The software used have been referenced at the end, the majority of this analysis is based on the user manuals linked below, with only slight deviations from default settings. 

*It is important to keep in mind there are MANY ways of analysing ChIP-seq data. Outlined below are just a few ways this analysis has been performed in the Johnstone/Kats lab. There are multiple useful links below in **Further information and useful tutorials** if you require alternative methods of analysis.*


### Library preparation
 
For more information about the ChIP-seq wetlab protocol please go [Here]()


### General overview of RNA-seq analysis

![GeneralOverview](https://github.com/madisonJK/ReferenceAnalysis/raw/master/ChIP-Seq/ChIP-seq_Overview.png)


### Software requirements and Version Control

Bowtie2 /2.3.3

Bedtools /2.26

samtools /1.4.1

FastQC /0.11.5

macs /2.1.1

Homer /4.8

Deeptools /2.5.3


R: 

Limma/

Deseq/ 

Csaw

ggplot2/3.3.1

reshape2/1.4.3

pheatmap/1.0.12



## Aligning Fastq Files

### QC
Perform QC on Fastq files using fastqc

*1_fastQC.sbatch*
```
module load fastqc

fastqc -o ${d} -f fastq --noextract -t 8 ${1}
```

### Align to reference genome with Bowtie2
Align Fastq files to the mouse genome (mm10) using Bowtie2. The resulting sam files are then sorted, converted to bam and indexed with samtools. IGV tools can then be used to visualise the TDF. You will need a Bowtie2 indexed reference genome. For paired ended data, use 2.2_ChIPseq_bowtie2_mouse_PE.sbatch. 
If you have an ERCC OR/AND an S2 spike-in, you will need to map to a merged reference genome e.g. mm10+ERCC OR mm10+DM3 OR mm10+DM3+ERCC.

*2_ChIPseq_bowtie2_mouse.sbatch*
```
module load bowtie2
module load samtools
module load igvtools

s=$(basename $1)
t=$(dirname $1)
x=`echo $s | cut -d "_" -f 1-2`
echo $s
echo $x

bowtie2 -p 32 -x /data/reference/indexes/mouse/mm10/bowtie2/Mus_musculus.GRCm38.dna.toplevel -U $1 -S ${1}.sam

samtools view -@ 8 -Sbo ${t}/${x}.sam.bam  ${t}/${x}.sam
samtools sort -@ 8 -o ${t}/${x}.sam.sorted.bam ${t}/${x}.sam.bam
samtools rmdup -s ${t}/${x}.sam.sorted.bam ${t}/${x}.sam.sorted.bam.rmdup.bam
samtools index ${t}/${x}.sam.sorted.bam.rmdup.bam ${t}/${x}.sam.sorted.bam.rmdup.bam.bai
igvtools count -z 5 -w 25 -e 225 ${t}/${x}.sam.sorted.bam.rmdup.bam ${t}/${x}.sam.sorted.bam.rmdup.bam.tdf mm10

```

## Peak Calling
This will produce a bed file of genomic coordinates of identified peaks or enriched regions. The resulting BED file can then be used as input for: peak annotation, motif analysis, differential peak analysis,input for deeptools compute matrix heatmap etc.

### MACS2
[HERE](https://hbctraining.github.io/Intro-to-ChIPseq/lessons/05_peak_calling_macs.html) is a great tutorial for peak calling with MACS. 

It is important to identify which settings will produce the best peaks

--broad : Use this option for broader-style enrichment.

--cutoffanalysis : Recommended to run this first to see how changing the p/q value impacts the number and width of the peaks.

see macs2 callpeak --help for more info

*3_MACS2.sbatch*

For broad peaks e.g. histone marks
```
module load macs

macs2 callpeak -g mm -f BAM -t IPSAMPLE.rmdup.bam -c INPUTSAMPLE.rmdup.bam --extsize 225 --broad --broad-cutoff 0.1 --cutoff-analysis --outdir IPSAMPLEmacs2 -n IPSAMPLEpeaks
```
For more narrow peaks e.g. transcription factor ChIP: 

```
module load macs

macs2 callpeak -g mm -f BAM -t $1 -c 181104_MK26RNA/SRA/SRR2132490.fastq.gz.sam.sorted.bam.rmdup.bam --extsize 225 --cutoff-analysis --outdir ${1}macs2 -n ${1}peaks
```


### Homer
Please see [HERE](http://homer.ucsd.edu/homer/ngs/peaks.html) for a more comprehensive tutorial on calling peaks with Homer.
First, a tag directory needs to be made. Make a tag directory for all samples, including input/igG control.
Homer uses .sam files, so you will need to open the final processed bam files as .sam files using samtools.

*3.1_HomerTagDir.sbatch*
```
module load homer
module load samtools

 
samtools view -h IP.bam > IP.sam
makeTagDirectory IPtagdir/ IP.sam -format sam

```
You can then use these tag directories to call peaks on your IP samples relative to your control.the peaks can then be converted to a bed file, which can then be annotated with homers annotatePeaks.pl. (You may need to add "chr" to the start if your reference genome does not include it - see line 3 below).

```
findPeaks IPtagdir/ -style factor -o auto -i InputControltagdir/ 
pos2bed.pl -o IPtagdir/IPpeaks.bed IPtagdir/IPpeaks.txt
awk '{print "chr"$0}' IPtagdir/IPpeaks.bed > IPtagdir/IPpeaks.chr.bed
annotatePeaks.pl IPtagdir/peaks.bed.chr.bed mm10 > IPtagdir/peaksannotated.txt
```
Homer can also be used for Superenhancer - style analysis: 

*3.1.1_HomerTagSuperE.sbatch*
 ``` 
 findPeaks IPtagdir/ -i InputControltagdir -style super -o IPtagdir/IPSuper.txt -superSlope -1000 -L 2 -typical IPnormalpeaks.txt

 ```

You can also annotate peaks using bedtools. You will need a reference genome in gtf or bed format. 

```
module load bedtools 
bedtools closest -D -a IPpeaks.chr.bed -b mm10reference.bed > IPpeaks.annotated.bed

```

### Finding enriched regions with Csaw
Csaw is an R package which employs window based analysis to calculate enriched regions in your IP sample compared to input control. The analysis is very similar to RNA-seq analysis.


## Motif analysis

Motif analysis can be performed on a peak bed file with Homer. 
You may sometimes have trouble with the homer config directories, try ```-preparsedDir IPsample/```

If you have the data, it may be a good idea to set the background as ATAC peaks from the same cell type to reduce bias for open regions of the genome. ```-bg ATACpeaks.bed```

*5_HomerFindMotifs.sbatch*

```
module load homer

findMotifsGenome.pl IPsamplePeaks.bed mm10 IPsample-motif 
```

## Generating Figures

### Heatmap with Deeptools 

you will need: bam files of ChIP, a reference bed file. this can include a reference genome with CDS or TSS, or it can be a peak summit bed file or a peak bed file generated by macs or homer above. 

You will first need to make bigwigs of your .bam files, including of the input control.
```
module load deeptools
bamCoverage -p 8 -e 225 --normalizeUsingRPKM -b IPsample.bam -o IPsample.bw
```

Then compute matrix for the regions you want to plot. 
scale-regions will this is good for enrichment accross gene features, e.g. metagene plots.

```
computeMatrix scale-regions -S IPsample.bw -R ref/mm10reference.bed -out IPsample-mm10genes.gz -b 1000 -a 1000
```

reference-point is better for plotting heatmap of enrichment relative to TSS, or peak summits. 
```
computeMatrix reference-point -S IPsample.bw  -R mm10TSS.bed -out IPsample-mm10TSS.gz -b 1000 -a 1000
```

From the created matrix you can then plot a heatmap:
*see plotHeatmap --help for more options*

```
plotHeatmap -m IPsample-mm10TSS.gz -out IPsample-mm10TSS.pdf   
``` 



## Acknowledgments


* **Madison Kelly** - *Author of this walkthrough* [madisonJK](https://github.com/madisonJK)

* Stephin Vervoort - *Helping with initial ChIP-seq analysis*


## Further information and useful tutorials

[Homer tutorials](http://homer.ucsd.edu/homer/index.html)

[Deeptools tutorials](https://deeptools.readthedocs.io/en/develop/)

[Csaw User Guide](https://www.bioconductor.org/packages/devel/workflows/vignettes/csawUsersGuide/inst/doc/csaw.pdf)

[Good ChIP-seq analysis tutorial](https://github.com/crazyhottommy/ChIP-seq-analysis) 


## References

Heinz, S., Benner, C., Spann, N., Bertolino, E., Lin, Y. C., Laslo, P., ... & Glass, C. K. (2010). Simple combinations of lineage-determining transcription factors prime cis-regulatory elements required for macrophage and B cell identities. Molecular cell, 38(4), 576-589.

Lun, A. T., & Smyth, G. K. (2015). csaw: a Bioconductor package for differential binding analysis of ChIP-seq data using sliding windows. Nucleic Acids Research, 44(5), e45.

Ramírez, F., Dündar, F., Diehl, S., Grüning, B. A., & Manke, T. (2014). deepTools: a flexible platform for exploring deep-sequencing data. Nucleic acids research, 42(W1), W187-W191.


Zhang, Y., Liu, T., Meyer, C. A., Eeckhoute, J., Johnson, D. S., Bernstein, B. E., ... & Liu, X. S. (2008). Model-based analysis of ChIP-Seq (MACS). Genome biology, 9(9), R137.



