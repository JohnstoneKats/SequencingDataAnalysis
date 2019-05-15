#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem-per-cpu=32G
#SBATCH --time=1:00:00   
#SBATCH --output=%j.stdout
#SBATCH --error=%j.error
#SBATCH --mail-type=ALL
#SBATCH --mail-user=user.name@petermac.org
#SBATCH --partition=prod_short



module load macs

s=$(basename $1)

echo $s
echo $1 


macs2 callpeak -g mm -f BAM -t $1 -c 181104_MK26RNA/SRA/SRR2132490.fastq.gz.sam.sorted.bam.rmdup.bam --extsize 225 --broad -q 0.01 --broad-cutoff 0.1 --cutoff-analysis --outdir ${1}BROADq01q1macs2 -n ${s}BROADq01q1peaks


