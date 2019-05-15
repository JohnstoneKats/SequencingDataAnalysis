# CRISPR screen analysis
genome wide screens and targeted screens 

# Trimming 

```
module load cutadapt 


cutadapt -g TGTGGAAAGGACGAAACACCG  -o ${1}.trimmed_P5.fastq.gz ${1}
cutadapt -a GTTTTAGAGCTAGAAATAGCAAG -o ${1}.trimmed_P7.fastq.gz ${1}.trimmed_P5.fastq.gz
cutadapt --maximum-length 20 -o ${1}.20bp.fastq.gz ${1}.trimmed_P7.fastq.gz  

```
