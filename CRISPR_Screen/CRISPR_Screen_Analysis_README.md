# CRISPR screen analysis
genome wide screens and targeted screens 

##Version control 

cutadapt /1.14

mageck /0.5.6

## Trimming 

```
module load cutadapt 

cutadapt -g TGTGGAAAGGACGAAACACCG  -o ${1}.trimmed_P5.fastq.gz ${1}
cutadapt -a GTTTTAGAGCTAGAAATAGCAAG -o ${1}.trimmed_P7.fastq.gz ${1}.trimmed_P5.fastq.gz
cutadapt --maximum-length 20 -o ${1}.20bp.fastq.gz ${1}.trimmed_P7.fastq.gz  

```

## Count guides with Mageck  

```
module load mageck
mageck count -l ListOfReferenceGuides.txt --fastq *trimmed_P7.fastq.gz --sample-label T0 Tend
```

## Statistical evaluation with Mageck

```
mageck test -k sample1.count.txt --norm-method total -t Tend -c T0 -n T0vsTend --pdf-report --remove-zero-threshold 0 --remove-zero control --additional-rra-parameters "--min-number-goodsgrna 2 --permutation 2000"
 
```
