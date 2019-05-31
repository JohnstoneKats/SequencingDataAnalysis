#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --mem-per-cpu=128G
#SBATCH --time=8:00:00   
#SBATCH --output=%j.stdout 
#SBATCH --error=%j.stderr
#SBATCH --mail-type=ALL
#SBATCH --mail-user=madison.kelly@petermac.org
# SBATCH --workdir=$1
#SBATCH --partition=prod

### Run this command from the folder with the sbatch files
### for dir in `ls`; do sbatch trimmed_anbalysis.sbatch $dir; done
### indices files: /data/reference/indexes/mouse/mm10/hisat2

### If experimenting: use srun --pty /bin/bash

module load star
module load samtools
module load igvtools

s=$(basename $1)
t=$(dirname $1)
x=`echo $s | cut -d "_" -f 1-2`
echo $s
echo $x


STAR --outFileNamePrefix ${t}/${x} --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMismatchNmax 2 --outFilterMatchNmin 0 --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --genomeDir ref/mm10STAR --runThreadN 8 --readFilesIn ${t}/${x}_R1.trimmed.fastq.gz.sorted.fastq.gz ${t}/${x}_R2.trimmed.fastq.gz.sorted.fastq.gz

STAR --runMode inputAlignmentsFromBAM --outFileNamePrefix ${t}/${x} --bamRemoveDuplicatesType UniqueIdentical --outWigType wiggle --outWigStrand Stranded ${t}/${x}Aligned.sortedByCoord.out.bam
#samtools view -@ 8 -Sbo ${t}/${x}.sam.bam  ${t}/${x}.sam
#samtools sort -@ 8 -o .${t}/${x}sam.sorted.bam ${t}/${x}.sam.bam
samtools rmdup -s ${t}/${x}Aligned.sortedByCoord.out.bam ${t}/${x}.rmdup.bam
samtools index ${t}/${x}.rmdup.bam ${t}/${x}.rmdup.bam.bai
igvtools count -z 5 -w 25 -e 225 ${t}/${x}.rmdup.bam ${t}/${x}.tdf mm10



