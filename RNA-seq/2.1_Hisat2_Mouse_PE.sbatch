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

module load samtools
module load igvtools
module load hisat2

s=$(basename $1)
t=$(dirname $1)
x=`echo $s | cut -d "_" -f 1-2`
echo $s
echo $x

hisat2 -p 16 -x /data/reference/indexes/mouse/mm10/bowtie2/Mus_musculus.GRCm38.dna.toplevel -1 ${t}/${x}_R1_001.fastq.gz -2 ${t}/${x}_R2_001.fastq.gz -S ${t}/${x}_hisat2.sam
samtools view -@ 16 -Sbo ${t}/${x}_hisat2.sam.bam ${t}/${x}_hisat2.sam
samtools sort -@ 16 -o ${t}/${x}_hisat2.sorted.bam ${t}/${x}_hisat2.sam.bam
samtools index ${t}/${x}_hisat2.sorted.bam ${t}/${x}_hisat2.sorted.bam.bai
igvtools count -z 5 -w 25 ${t}/${x}_hisat2.sorted.bam ${t}/${x}.tdf mm10

rm ${1}_hisat2.sam
rm ${1}_hisat2.sam.bam

