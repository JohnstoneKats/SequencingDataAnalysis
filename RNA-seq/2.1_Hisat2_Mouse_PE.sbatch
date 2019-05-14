#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem-per-cpu=16G
#SBATCH --time=4:00:00   
#SBATCH --output=%j.stdout
#SBATCH --error=%j.error
#SBATCH --mail-type=ALL
#SBATCH --mail-user=user.name@petermac.org
#SBATCH --partition=prod

module load bowtie2
module load samtools
module load igvtools
module load hisat2

s=$(basename $1)
t=$(dirname $1)
x=`echo $s | cut -d "_" -f 1-2`
echo $s
echo $x

hisat2 -p 16 -x /./researchers/stephin.vervoort/HG19_ERCC_hisat2/Homo_sapiens_GRCh37.73.dna.toplevel_ERCC92 -1 ${t}/${x}_R1_001.fastq.gz -2 ${t}/${x}_R2_001.fastq.gz -S stephinINTS6RNA/${x}_hisat2.sam
samtools view -@ 16 -Sbo stephinINTS6RNA/${x}_hisat2.sam.bam stephinINTS6RNA/${x}_hisat2.sam
samtools sort -@ 16 -o stephinINTS6RNA/${x}_hisat2.sam.sorted.bam stephinINTS6RNA/${x}_hisat2.sam.bam
samtools index stephinINTS6RNA/${x}_hisat2.sam.sorted.bam stephinINTS6RNA/${x}_hisat2.sam.sorted.bam.bai
igvtools count -z 5 -w 25 stephinINTS6RNA/${x}_hisat2.sam.sorted.bam stephinINTS6RNA/${x}.tdf hg19


