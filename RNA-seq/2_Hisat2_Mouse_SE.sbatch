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
x=`echo $s | cut -d "_" -f 1-2`

echo $s
echo $x

hisat2 -p 8 -x /data/reference/indexes/mouse/mm10/hisat2/grcm38/Mus_musculus.GRCm38.dna.toplevel -U $1 -S ${1}_hisat2.sam
samtools view -@ 4 -Sbo ${1}_hisat2.sam.bam ${1}_hisat2.sam
samtools sort -@ 4 -o ${1}_hisat2.sorted.bam ${1}_hisat2.sam.bam
samtools index ${1}_hisat2.sorted.bam ${1}_hisat2.sorted.bam.bai
igvtools count -z 5 -w 25 -e 150 ${1}_hisat2.sorted.bam ${1}.tdf mm10

rm ${1}_hisat2.sam
rm ${1}_hisat2.sam.bam
