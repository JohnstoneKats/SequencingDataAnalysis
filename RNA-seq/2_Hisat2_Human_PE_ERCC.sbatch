#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem-per-cpu=16G
#SBATCH --time=12:00:00   
#SBATCH --output=%j.stdout
#SBATCH --error=%j.error
#SBATCH --mail-type=ALL
#SBATCH --mail-user=madison.kelly@petermac.org
# SBATCH --workdir=$1
#SBATCH --partition=prod

### Run this command from the folder with the sbatch files
### for dir in `ls`; do sbatch trimmed_anbalysis.sbatch $dir; done
### indices files: /data/reference/indexes/mouse/mm10/hisat2

### If experimenting: use srun --pty /bin/bash

module load bowtie2
module load samtools
module load igvtools
module load hisat2

###hisat2 [options]* -x <ht2-idx> {-1 <m1> -2 <m2> | -U <r> | --sra-acc <SRA accession number>} [-S <sam>]

#file=$(ls | grep fastq.gz$)

s=$(basename $1)
t=$(dirname $1)
x=`echo $s | cut -d "_" -f 1-2`
echo $s
echo $x
###echo $file

hisat2 -p 16 -x /./researchers/stephin.vervoort/HG19_ERCC_hisat2/Homo_sapiens_GRCh37.73.dna.toplevel_ERCC92 -1 ${t}/${x}_R1_001.fastq.gz -2 ${t}/${x}_R2_001.fastq.gz -S stephinINTS6RNA/${x}_hisat2.sam


samtools view -@ 16 -Sbo stephinINTS6RNA/${x}_hisat2.sam.bam stephinINTS6RNA/${x}_hisat2.sam
samtools sort -@ 16 -o stephinINTS6RNA/${x}_hisat2.sam.sorted.bam stephinINTS6RNA/${x}_hisat2.sam.bam
#samtools rmdup -s ${1}_bowtie2.sam.sorted.bam ${1}_bowtie2.sam.sorted.bam.rmdup.bam
samtools index stephinINTS6RNA/${x}_hisat2.sam.sorted.bam stephinINTS6RNA/${x}_hisat2.sam.sorted.bam.bai
igvtools count -z 5 -w 25 stephinINTS6RNA/${x}_hisat2.sam.sorted.bam stephinINTS6RNA/${x}.tdf hg19


