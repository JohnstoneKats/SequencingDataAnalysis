#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mem-per-cpu=64G
#SBATCH --time=2:00:00   
# SBATCH --output=%j.stdout
# SBATCH --error=%j.error
#SBATCH --mail-type=ALL
#SBATCH --mail-user=madison.kelly@petermac.org
#SBATCH --partition=prod_short
 
module load cutadapt

s=$(basename $1)
x=`echo $s | cut -d "_" -f 1-2`
t=$(dirname $1)

echo $1

cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -G AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -q 30,30 -m 30 -o ${t}/${x}_R1.trimmed.fastq.gz -p ${t}/${x}_R2.trimmed.fastq.gz  ${t}/${x}_R1.fastq.gz ${t}/${x}_R2.fastq.gz
 
