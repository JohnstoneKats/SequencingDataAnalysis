#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem-per-cpu=64G
#SBATCH --time=12:00:00   
#SBATCH --output=%j.stdout
#SBATCH --error=%j.error
#SBATCH --mail-type=ALL
#SBATCH --mail-user=madison.kelly@petermac.org
#SBATCH --partition=prod

### Run this command from the folder with the sbatch files
### for dir in `ls`; do sbatch trimmed_anbalysis.sbatch $dir; done
### indices files: /data/reference/indexes/mouse/mm10/hisat2

module load bowtie2
module load samtools
module load igvtools
module load bismark

s=$(basename $1)
t=$(dirname $1)
x=`echo $s | cut -d "_" -f 1-2`

echo $s
echo $x

bismark --unmapped -X 1000 --multicore 4 ../../../data/reference/indexes/mouse/mm10/bowtie2/ -1 ${t}/${x}_R1_001_val_1.fq.gz -2 ${t}/${x}_R2_001_val_2.fq.gz 

bismark --multicore 4 ../../../data/reference/indexes/mouse/mm10/bowtie2/ ${x}_R1_001_val_1.fq.gz_unmapped_reads_1.fq.gz 

bismark --multicore 4 --pbat ../../../data/reference/indexes/mouse/mm10/bowtie2 ${x}_R2_001_val_2.fq.gz_unmapped_reads_2.fq.gz

samtools merge ${x}merged.bam ${x}_R1_001_val_1_bismark_bt2_pe.bam ${x}_R1_001_val_1.fq.gz_unmapped_reads_1_bismark_bt2.bam ${x}_R2_001_val_2.fq.gz_unmapped_reads_2_bismark_bt2.bam

samtools sort ${x}merged.bam > ${x}mergedsort.bam

samtools index ${x}mergedsort.bam



