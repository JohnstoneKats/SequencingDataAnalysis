#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mem-per-cpu=32G
#SBATCH --time=02:00:00   
#SBATCH --output=%j.stdout 
#SBATCH --error=%j.stderr
#SBATCH --mail-type=ALL
#SBATCH --mail-user=user.name@petermac.org
#SBATCH --partition=prod_short

module load bowtie2
module load samtools
module load igvtools

s=$(basename $1)
t=$(dirname $1)
x=`echo $s | cut -d "_" -f 1-2`
echo $s
echo $x

bowtie2 -p 32 -x /data/reference/indexes/mouse/mm10/bowtie2/Mus_musculus.GRCm38.dna.toplevel -U $1 -S ${1}.sam

samtools view -@ 8 -Sbo ${1}.sam.bam  ${1}.sam
samtools sort -@ 8 -o ${1}.sorted.bam ${t}/${x}.sam.bam
samtools rmdup -s ${1}.sorted.bam ${1}.sorted.rmdup.bam
samtools index ${1}.sorted.rmdup.bam ${1}.sorted.rmdup.bam.bai
igvtools count -z 5 -w 25 -e 225 ${1}.sorted.rmdup.bam ${1}.tdf mm10

rm ${1}.sam
rm ${1}.sam.bam
rm ${1}.sorted.bam
