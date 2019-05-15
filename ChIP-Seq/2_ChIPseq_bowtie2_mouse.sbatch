#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=6
#SBATCH --mem-per-cpu=32G
#SBATCH --time=01:00:00   
#SBATCH --output=%j.stdout 
#SBATCH --error=%j.stderr
#SBATCH --mail-type=ALL
#SBATCH --mail-user=madison.kelly@petermac.org
#SBATCH --partition=prod_short

### Run this command from the folder with the sbatch files
### for dir in `ls`; do sbatch trimmed_anbalysis.sbatch $dir; done
### indices files: /data/reference/indexes/mouse/mm10/hisat2

### If experimenting: use srun --pty /bin/bash

module load bowtie2

module load samtools
module load igvtools

###hisat2 [options]* -x <ht2-idx> {-1 <m1> -2 <m2> | -U <r> | --sra-acc <SRA accession number>} [-S <sam>]

#file=$(ls | grep fastq.gz$)

s=$(basename $1)
t=$(dirname $1)
x=`echo $s | cut -d "_" -f 1-2`
echo $s
echo $x
###echo $file
bowtie2 -p 32 -x /data/reference/indexes/mouse/mm10/bowtie2/Mus_musculus.GRCm38.dna.toplevel -U $1 -S ${1}.sam

samtools view -@ 8 -Sbo ${t}/${x}.sam.bam  ${t}/${x}.sam
samtools sort -@ 8 -o ${t}/${x}.sam.sorted.bam ${t}/${x}.sam.bam
samtools rmdup -s ${t}/${x}.sam.sorted.bam ${t}/${x}.sam.sorted.bam.rmdup.bam
samtools index ${t}/${x}.sam.sorted.bam.rmdup.bam ${t}/${x}.sam.sorted.bam.rmdup.bam.bai
igvtools count -z 5 -w 25 -e 225 ${t}/${x}.sam.sorted.bam.rmdup.bam ${t}/${x}.sam.sorted.bam.rmdup.bam.tdf mm10


