#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem-per-cpu=32G
#SBATCH --time=4:00:00   
#SBATCH --output=%j.stdout 
#SBATCH --error=%j.stderr
#SBATCH --mail-type=ALL
#SBATCH --mail-user=user.name@petermac.org
#SBATCH --partition=prod


module load bowtie2
module load samtools
module load igvtools


s=$(basename $1)
t=$(dirname $1)
x=`echo $s | cut -d "_" -f 1`
echo $s
echo $x


bowtie2 -p 32 -x /data/reference/indexes/mouse/mm10/bowtie2/Mus_musculus.GRCm38.dna.toplevel -1 ${t}/${x}_1.fastq.gz -2 ${t}/${x}_2.fastq.gz -S ${t}/${x}.sam


samtools view -@ 8 -Sbo ${t}/${x}.sam.bam  ${t}/${x}.sam
samtools sort -@ 8 -o ${t}/${x}.sorted.bam ${t}/${x}.sam.bam
samtools rmdup -s ${t}/${x}.sorted.bam ${t}/${x}.sorted.rmdup.bam
samtools index ${t}/${x}.sorted.rmdup.bam ${t}/${x}.sorted.rmdup.bam.bai
igvtools count -z 5 -w 25 -e 225 ${t}/${x}.sorted.rmdup.bam ${t}/${x}.tdf mm10


rm ${t}/${x}.sam
rm ${t}/${x}.sam.bam
rm ${t}/${x}.sorted.bam

