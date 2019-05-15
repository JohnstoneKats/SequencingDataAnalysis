# Downloading SRA datasets

Input is a "SRR..." number
```
module load sratoolkit

fastq-dump -I --split-3 --skip-technical --readids --clip --dumpbase --gzip -O 1SRA -A $1
```
