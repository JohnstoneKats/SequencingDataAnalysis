# Downloading SRA datasets

Input is a "SRR..." number, this will download the files to a folder named SRA/. 

```
module load sratoolkit

fastq-dump -I --split-3 --skip-technical --readids --clip --dumpbase --gzip -O SRA/ -A SRR...
```
 
 For multiple SRAs: 
  
 ```
 for i in SRR1... SRR2... SRR3... ; do fastq-dump -I --split-3 --skip-technical --readids --clip --dumpbase --gzip -O SRA/ -A $i; done 
 ```
