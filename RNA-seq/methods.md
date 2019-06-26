# RNA-seq Methods

*breif general analysis methods to be used for protocol section in paper*


**Quant-seq**

1-5 x 105 sorted KL cells were lysed in Trizol reagent and RNA was extracted using Directzol RNA miniprep kit (Zymo research) according to manufacturer’s instructions. The Quant-seq 3’ mRNA-seq Library Prep Kit for Illumina (Lexogen) was used to generate libraries as per the manufacturer’s instructions. Libraries were pooled and sequenced with  75bp single end sequencing to a depth of 6-10 x 106 reads on a NextSeq (Illumina). Sequencing reads were demultiplexed using bcl2fastq (v2.17.1.14) and low-quality reads Q<30 were removed. The RNA sequencing reads were trimmed at the 5’ end using cutadapt (v1.14)67 to remove bias introduced by random primers used in the library preparation and 3’ end trimming was performed to eliminate poly-A-tail derived reads. Reads were mapped to the reference genome (mm10) using HISAT2. Reads were counted using subread software package (v1.5.0-p3)68. Differential gene expression analysis was performed using R package LIMMA (v3.32.4). R packages pheatmap (v1.0.8) and ggplot2 (v2.2.1) were used for figure generation. GSEA (v3.0) was used for analysing enrichment of gene sets and gene ontology (GO) analysis was performed using metascape [http://metascape.org]. 

from [Bcor loss perturbs myeloid differentiation and promotes leukaemogenesis](https://www.nature.com/articles/s41467-019-09250-6) 
